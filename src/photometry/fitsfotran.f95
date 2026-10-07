!
!  FITS I/O for photometric calibration
!
!  Copyright © 2014-5, 2018 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Munipack.
!
!  Munipack is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Munipack is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with Munipack.  If not, see <http://www.gnu.org/licenses/>.

module fits_fotran

  use titsio

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)

  private :: idcol

contains

  subroutine traload(nametable,photosys_ref,photosys_instr, &
       tr,trerr,tr1,tr1err,status)

    character(len=*), intent(in) :: nametable
    character(len=*), intent(out) :: photosys_ref,photosys_instr
    real(dbl), dimension(:,:), allocatable, intent(out) :: tr,trerr,tr1,tr1err
    integer, intent(out) :: status

    character(len=FLEN_KEYWORD) :: key
    character(len=FLEN_VALUE) :: extname
    integer :: ncols,nrows,i
    integer, parameter :: frow = 1, extver = 0
    real(dbl), parameter :: nullval = 0.0_dbl
    logical :: anyf
    type(fitsfiles) :: fits

    status = 0

    ! open and move to first table extension
    call fits_open_file(fits,nametable,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read a table in the file `',trim(nametable),"'."
       return
    end if

    if( status /= 0 ) goto 666

    ! check reference frame and identification of this table
    call fits_make_keyn(FITS_KEY_PHOTSYS,1,key,status)
    call fits_read_key(fits,key,photosys_instr,status)
    call fits_make_keyn(FITS_KEY_PHOTSYS,2,key,status)
    call fits_read_key(fits,key,photosys_ref,status)

    extname = trim(FTHDUNAME) // '_FORWARD'
    call fits_movnam_hdu(fits,FITS_ANY_HDU, extname, extver,status)
    if( status /= 0 ) goto 666

    call fits_get_num_rows(fits,nrows,status)
    call fits_get_num_cols(fits,ncols,status)

    allocate(tr(nrows,ncols))
    do i = 1, ncols
       call fits_read_col(fits,i,frow,nullval,tr(:,i),anyf,status)
    end do

    extname = trim(FTHDUNAME) // '_BACKWARD'
    call fits_movnam_hdu(fits,FITS_ANY_HDU, extname, extver,status)
    if( status /= 0 ) goto 666
    call fits_get_num_rows(fits,nrows,status)
    call fits_get_num_cols(fits,ncols,status)
    allocate(tr1(nrows,ncols))
    do i = 1, ncols
       call fits_read_col(fits,i,frow,nullval,tr1(:,i),anyf,status)
    end do

    allocate(trerr(nrows,ncols),tr1err(nrows,ncols))
    trerr = -1
    tr1err = -1

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( status /= 0 ) then
       if( allocated(tr) ) deallocate(tr)
       if( allocated(tr1) ) deallocate(tr1)
       if( allocated(trerr) ) deallocate(trerr)
       if( allocated(tr1err) ) deallocate(tr1err)
    end if

  end subroutine traload


  subroutine trawrite(output,photosys_ref,photosys_instr,filters_ref, &
       filters_instr,ra,dec,ci0,ph,dph,cts,dcts, &
       tr,trerr,tr1,tr1err,filenames,airmass,ctph,extin)

    character(len=*), intent(in) :: output, photosys_ref, photosys_instr
    character(len=*), dimension(:), intent(in) :: filters_ref, filters_instr
    character(len=*), dimension(:), intent(in), optional :: filenames
    real(dbl), dimension(:), intent(in) :: ra,dec
    real(dbl), dimension(:), intent(in), optional :: airmass,ctph,extin
    real(dbl), intent(in) :: ci0
    real(dbl), dimension(:,:), intent(in) :: tr,trerr,tr1,tr1err, &
         ph,dph,cts,dcts

    integer, parameter :: frow = 1
    character(len=FLEN_VALUE), dimension(:), allocatable :: ttype, tform,tunit
    character(len=FLEN_VALUE) :: extname
    character(len=FLEN_KEYWORD) :: key
    character(len=FLEN_CARD) :: buf,a
    real(dbl) :: x
    integer :: n,nrows,ncols, nbands,status
    type(fitsfiles) :: fits

    status = 0
    if( fits_file_exist(output) ) call fits_file_delete(output)
    call fits_create_file(fits,output,status)
    call fits_insert_image(fits,8,0,[0],status)
    call fits_write_comment(fits,'Photometry system conversion table.',status)
    call fits_write_key(fits,FITS_KEY_CREATOR,FITS_VALUE_CREATOR,FITS_COM_CREATOR,status)
    call fits_write_comment(fits,MUNIPACK_VERSION,status)
    call fits_write_comment(fits,'Description: http://munipack.physics.muni.cz/dataform_phfotran.html',status)

    call fits_make_keyn(FITS_KEY_PHOTSYS,1,key,status)
    call fits_write_key(fits,key,photosys_instr,'instrumental photometry system',status)

    call fits_make_keyn(FITS_KEY_PHOTSYS,2,key,status)
    call fits_write_key(fits,key,photosys_ref,'reference photometry system (standard)',status)

    if( present(airmass) ) then
       x = sum(airmass)/size(airmass)
       call fits_write_key(fits,FITS_KEY_AIRMASS,x,-6,'average air-mass',status)
    end if

    ! mean colour index
    call fits_write_key(fits,'REFCI',ci0,-6,'reference colour index',status)

    if( present(filenames) ) then
       call fits_write_comment(fits,'Photometry transformation has been derived by processing..',status)
       call fits_write_comment(fits,' filename, filter, airmass, ctph, extinction:',status)
       do n = 1,size(filenames)
          if( size(airmass) > 0 ) then
             write(buf,'(f8.5)') airmass(n)
          else
             buf = "-"
          end if
          if( size(ctph) > 0 ) then
             write(a,'(1p,g0.5)') ctph(n)
          else
             a = "-"
          end if
          buf = trim(buf) // " " // trim(a)
          if( size(extin) > 0 ) then
             write(a,'(1p,g0.5)') extin(n)
          else
             a = "-"
          end if
          buf = trim(buf) // " " // trim(a)
          call fits_write_comment(fits,"'"//trim(filenames(n))//"' '"//  &
               trim(filters_instr(n))//"' "//trim(buf),status)
       end do
       call fits_write_comment(fits,'All elements are normalized for both unit area and time.',status)
    end if


    ! extension - FOTRAN
    nrows = size(tr,1)
    ncols = size(tr,2)
    allocate(ttype(ncols),tform(ncols),tunit(ncols))
    tform = '1D'
    tunit = '' !'an instrumental filter'
    ttype = filters_ref

    extname = trim(FTHDUNAME) // '_FORWARD'
    call fits_insert_btbl(fits,nrows,ttype,tform,tunit,extname,status)
    call fits_write_comment(fits,'Forward matrix (photon rate to counts rate).',status)
    do n = 1, ncols
       call fits_write_col(fits,n,frow,tr(:,n),status)
    end do

    extname = trim(FTHDUNAME) // '_BACKWARD'
    call fits_insert_btbl(fits,nrows,ttype,tform,tunit,extname,status)
    call fits_write_comment(fits,'Backward matrix (counts rate to photon rate).',status)
    do n = 1, ncols
       call fits_write_col(fits,n,frow,tr1(:,n),status)
    end do

    extname = trim(FTHDUNAME) // '_FORWARD_ERR'
    call fits_insert_btbl(fits,nrows,ttype,tform,tunit,extname,status)
    do n = 1, ncols
       call fits_write_col(fits,n,frow,trerr(:,n),status)
    end do

    extname = trim(FTHDUNAME) // '_BACKWARD_ERR'
    call fits_insert_btbl(fits,nrows,ttype,tform,tunit,extname,status)
    do n = 1, ncols
       call fits_write_col(fits,n,frow,tr1err(:,n),status)
    end do

    deallocate(ttype,tform,tunit)


    ! data info extension
    nrows = size(ra)
    nbands = size(ph,2)
    ncols = (2*2)*nbands + 2
    allocate(ttype(ncols),tform(ncols),tunit(ncols))
    tform = '1D'
    tunit(1:2) = 'deg'
    tunit(3:2+2*nbands) = 'ph/s/m2'
    tunit(3+2*nbands:2+4*nbands) = 'W/m2'
    tunit(3+4*nbands:) = 'cts/s/m2'
    ttype(1) = FITS_COL_RA
    ttype(2) = FITS_COL_DEC
    ttype(3:2+1*nbands) = 'PH_'//filters_ref
    ttype(3+1*nbands:2+2*nbands) = 'ePH_'//filters_ref
    ttype(3+2*nbands:2+3*nbands) = 'CTS_'//filters_instr
    ttype(3+3*nbands:) = 'eCTS_'//filters_instr

    call fits_insert_btbl(fits,nrows,ttype,tform,tunit,trim(FTHDUNAME)//'data',status)
    call fits_write_col(fits,1,frow,ra,status)
    call fits_write_col(fits,2,frow,dec,status)
    do n = 1,nbands
       call fits_write_col(fits,2+n,frow,ph(:,n),status)
    end do
    do n = 1,nbands
       call fits_write_col(fits,2+nbands+n,frow,dph(:,n),status)
    end do
    do n = 1,nbands
       call fits_write_col(fits,2+2*nbands+n,frow,cts(:,n),status)
    end do
    do n = 1,nbands
       call fits_write_col(fits,2+3*nbands+n,frow,dcts(:,n),status)
    end do

    deallocate(ttype,tform,tunit)

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine trawrite

  subroutine tradata(nametable,ra,dec,ph,dph,cts,dcts)

    character(len=*), intent(in) :: nametable
    real(dbl), dimension(:), allocatable, intent(out) :: ra,dec
    real(dbl), dimension(:,:), allocatable, intent(out) :: ph,dph,cts,dcts

    integer :: ncols,nrows,i,n,rcol,dcol,status
    integer, parameter :: frow = 1, felem = 1, extver = 0
    real(dbl), parameter :: nullval = 0.0_dbl
    integer, dimension(:), allocatable :: col_cts, col_ph, col_dcts, col_dph
    logical :: anyf
    type(fitsfiles) :: fits

    status = 0

    ! open and move to first table extension
    call fits_open_table(fits,nametable,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read a table in the file `',trim(nametable),"'."
       return
    end if

    call fits_movnam_hdu(fits,FITS_BINARY_TBL,trim(FTHDUNAME)//'data', extver,status)
    call fits_get_num_rows(fits,nrows,status)
    call fits_get_num_cols(fits,ncols,status)
    if( status /= 0 ) goto 666

    call fits_get_colnum(fits,.true.,FITS_COL_RA,rcol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_DEC,dcol,status)
    if( status /= 0 ) goto 666

    allocate(col_cts(ncols),col_ph(ncols),col_dph(ncols),col_dcts(ncols))

    call idcol(fits,'CTS_*',col_cts,n,status)
    call idcol(fits,'eCTS_*',col_dcts,n,status)
    call idcol(fits,'PH_*',col_ph,n,status)
    call idcol(fits,'ePH_*',col_dph,n,status)

    if( size(col_ph) /= size(col_cts) .or. n == 0 ) &
         stop 'CTS or PH columns badly identified.'

    allocate(cts(nrows,n),ph(nrows,n),dcts(nrows,n),dph(nrows,n),ra(nrows),dec(nrows))
    call fits_read_col(fits,rcol,frow,nullval,ra,anyf,status)
    call fits_read_col(fits,dcol,frow,nullval,dec,anyf,status)
    do i = 1,n
       call fits_read_col(fits,col_cts(i),frow,nullval,cts(:,i),anyf,status)
    end do
    do i = 1,n
       call fits_read_col(fits,col_dcts(i),frow,nullval,dcts(:,i),anyf,status)
    end do
    do i = 1,n
       call fits_read_col(fits,col_ph(i),frow,nullval,ph(:,i),anyf,status)
    end do
    do i = 1,n
       call fits_read_col(fits,col_dph(i),frow,nullval,dph(:,i),anyf,status)
    end do
    deallocate(col_cts,col_ph)

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine tradata

  subroutine idcol(fits,temp,cols,n,status)

    type(fitsfiles), intent(in) :: fits
    character(len=*), intent(in) :: temp
    integer, dimension(:), intent(out) :: cols
    integer, intent(out) :: n,status
    integer :: l

    n = 0
    do
       call fits_get_colnum(fits,.true.,temp,l,status)
       if( status == FITS_MULTIPLE_MATCH ) then
          n = n + 1
          cols(n) = l
       else
          exit
       end if
    end do
    status = 0

  end subroutine idcol

end module fits_fotran
