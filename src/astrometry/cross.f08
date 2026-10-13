!
! Catalogue cross-match
!
! Copyright © 2013, 2018-20 F.Hroch (hroch@physics.muni.cz)
!
! This file is part of Munipack.
!
! Munipack is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Munipack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
!

program cross

  use iso_fortran_env
  use titsio

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)
  integer, parameter :: NFILE = 2

  integer :: eq
  character(len=FLEN_VALUE), dimension(NFILE) :: col_ra, col_dec, &
       col_pmra, col_pmdec, col_mag
  character(len=FLEN_KEYWORD) :: key_epoch
  character(len=4*FLEN_FILENAME) :: record, key, val, output
  character(len=4*FLEN_FILENAME), dimension(NFILE) :: file
  logical :: verbose = .false.
  real(dbl) :: tol,ftol
  integer :: n

  key_epoch = FITS_KEY_EPOCH
  col_ra = 'RA'
  col_dec = 'DEC'
  col_pmra = ''
  col_pmdec = ''
  col_mag = ''

  file=''
  output = ''
  tol = epsilon(tol)
  ftol = 1

  do

     read(*,'(a)',end=99) record

     eq = index(record,'=')
     if( eq == 0 ) stop 'Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'OUTPUT' ) then

        read(val,*) output

     else if( key == 'NFILES' ) then

        read(val,*) n
        if( NFILE /= n ) stop 'Cross-match requires exactly two files.'
        n = 0

     else if( key == 'COL_RA' ) then

        read(val,*) col_ra

     else if( key == 'COL_DEC' ) then

        read(val,*) col_dec

     else if( key == 'COL_PMRA' ) then

        read(val,*) col_pmra

     else if( key == 'COL_PMDEC' ) then

        read(val,*) col_pmdec

     else if( key == 'COL_MAG' ) then

        read(val,*) col_mag

     else if( key == 'FITS_KEY_EPOCH' ) then

        read(val,*) key_epoch

     else if( key == 'TOL' ) then

        read(val,*) tol

     else if( key == 'FTOL' ) then

        read(val,*) ftol

     else if( key == 'FILE' ) then

        n = n + 1
        read(val,*) file(n)

     end if

  end do
99 continue

  call crossing

  stop 0

contains


  subroutine crossing

    use crossmatch

    real(dbl), allocatable, dimension(:) :: ra1,dec1,flux1,ra2,dec2,flux2
    integer, dimension(:,:), allocatable :: idx
    character(len=FLEN_VALUE) :: catid1,catid2
    integer :: status

    status = 0
    call tabcoo(file(1),col_ra(1),col_dec(1),col_pmra(1),col_pmdec(1),&
         col_mag(1),key_epoch,ra1,dec1,flux1,catid1,status)
    call tabcoo(file(2),col_ra(2),col_dec(2),col_pmra(2),col_pmdec(2),&
         col_mag(2),key_epoch,ra2,dec2,flux2,catid2,status)
    if( status /= 0 .or. .not. (allocated(ra1) .and. allocated(ra2)) ) &
         stop 'Error: failed to read the coordinates of the tables.'

    allocate(idx(NFILE,max(size(ra1),size(ra2))))

    call join(ra1,dec1,ra2,dec2,tol,ftol,idx)

    call savecross(output,file,idx,status)

    if( allocated(ra1) ) deallocate(ra1,dec1,flux1)
    if( allocated(ra2) ) deallocate(ra2,dec2,flux2)
    deallocate(idx)

  end subroutine crossing


  subroutine tabcoo(filename,col_ra,col_dec,col_pmra,col_pmdec,&
         col_mag,key_epoch,ra,dec,flux,catid,status)

    use trajd

    character(len=*), intent(in) :: filename, col_ra, col_dec, col_pmra, &
         col_pmdec, col_mag, key_epoch
    real(dbl), allocatable, dimension(:), intent(out) :: ra, dec, flux
    character(len=*), intent(out) :: catid
    integer, intent(in out) :: status

    real(dbl), parameter :: nullval = real(0.0,dbl)
    integer, parameter :: ncols = 5
    integer, dimension(ncols) :: col, statuses
    integer :: nrows,i,l,frow,srows
    logical :: anyf
    character(len=FLEN_VALUE), dimension(ncols) :: label,colname
    real(dbl), allocatable, dimension(:) :: alpha, delta, pmalpha, pmdelta, mag
    real(dbl) :: t0,t,dt
    type(fitsfiles) :: fits

    if( status /= 0 ) return

    ! ICRS
    t0 = yearjd(2000.00_dbl)

    label = [ character(len=FLEN_VALUE) :: col_ra, col_dec, col_pmra, col_pmdec, &
         col_mag]

    if( label(1) == '' .or. label(2) == '' ) then
       write(error_unit,*) "Coordinate columns are not fully specified (empty)."
       return
    end if

    ! open and move to a table extension
    call fits_open_table(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read a table in the file `',trim(filename),"'."
       return
    end if

    call fits_get_num_rows(fits,nrows,status)
    if( .not. (nrows > 0 .and. status == 0) ) then
       write(error_unit,*) "Table `"//trim(filename)//"' looks empty."
       goto 666
    end if

    ! identification of catalogue
    call fits_read_key(fits,'EXTNAME',catid,status)
    call fits_read_key(fits,key_epoch,t,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       write(error_unit,*) "Epoch by `"//trim(key_epoch)// &
            "' keyword not found (assumed J2000.0)."
       status = 0
       t = t0
    end if

    ! find columns by labels
    statuses = 0
    col = 0
    do i = 1, ncols
       if( label(i) /= '' ) then
          call fits_get_colname(fits,.false.,label(i),colname(i),col(i),statuses(i))
       else
          col(i) = 0
       end if
    end do

    if( col(1) == 0 .or. col(2) == 0 ) then
       write(error_unit,*) "Error: coordinates column `"//trim(col_ra)//","// &
            trim(col_dec)//"' in `"//trim(filename)//"' not found."
       goto 666
    end if

    if( (label(3) /= '' .or. label(4) /= '') .and. &
         (statuses(3) /= 0 .or. statuses(4) /= 0 ) ) then
       write(error_unit,*) "Error: requested proper-motion columns `" // &
            trim(col_pmra)//"' or `"//trim(col_pmdec)//"' in `"// &
            trim(filename)//"' not found."
       goto 666
    end if

    if( label(5) /= '' .and. statuses(5) /= 0 ) then
       write(error_unit,*) "Error: a requested magnitude column `" // &
            trim(col_mag)//"' in `"//trim(filename)//"' not found."
       goto 666
    end if

    allocate(alpha(nrows),delta(nrows),pmalpha(nrows),pmdelta(nrows),mag(nrows))

    pmalpha = 0.0_dbl
    pmdelta = 0.0_dbl
    mag = 99.999

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,col(1),frow,nullval,alpha(i:l),anyf,status)
       call fits_read_col(fits,col(2),frow,nullval,delta(i:l),anyf,status)
       if( col(3) /= 0 .and. col(4) /= 0) then
          call fits_read_col(fits,col(3),frow,nullval,pmalpha(i:l),anyf,status)
          call fits_read_col(fits,col(4),frow,nullval,pmdelta(i:l),anyf,status)
       end if
       if( col(5) /= 0 ) then
          call fits_read_col(fits,col(5),frow,nullval,mag(i:l),anyf,status)
       end if
       if( status /= 0 ) goto 666
    end do

    pmalpha = pmalpha / 3.6e6_dbl  ! supposing mas/year
    pmdelta = pmdelta / 3.6e6_dbl

    ! compute coordinates for the epoch (the file is closed below,
    ! a second close failed with an invalid fitsfile pointer)
    allocate(ra(nrows),dec(nrows),flux(nrows))
    dt = (t - t0)/ 365.25

    ra  = alpha + dt*pmalpha
    dec = delta + dt*pmdelta
    flux = 10.0**((25.0 - mag)/2.5)

666 continue

    if( allocated(alpha) ) deallocate(alpha,delta,pmalpha,pmdelta,mag)

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine tabcoo

  subroutine savecross(output,filename,idx,status)

    character(len=*), intent(in) :: output
    character(len=*), dimension(:), intent(in) :: filename
    integer, dimension(:,:), intent(in) :: idx
    integer, intent(in out) :: status

    character(len=FLEN_VALUE), dimension(:), allocatable :: ttype, tform,tunit
    character(len=FLEN_KEYWORD) :: key
    integer :: nrows,i, l,frow, srows
    type(fitsfiles) :: fits

    status = 0
    if( fits_file_exist(output) ) call fits_file_delete(output)
    call fits_create_file(fits,output,status)

    allocate(ttype(NFILE),tform(NFILE),tunit(NFILE))

    tform = '1J'
    tunit = ''

    ttype = filename
    call fits_insert_btbl(fits,0,ttype,tform,tunit,'CROSSMATCH',status)

    do i = 1, NFILE
       write(key,'(a,i1)') 'FILE',i
       call fits_update_key(fits,key,filename(i),'FITS table',status)
    end do

    nrows = size(idx,2)
    call fits_get_rowsize(fits,srows,status)

    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       do n = 1, NFILE
          call fits_write_col(fits,n,frow,idx(n,i:l),status)
       end do
    end do

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    deallocate(ttype,tform,tunit)

  end subroutine savecross

end program cross
