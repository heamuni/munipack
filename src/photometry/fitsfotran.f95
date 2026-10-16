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

  ! the model of the transformation: FLUX (the forward matrix of photon
  ! rates, assumed when the keyword is missing) or MAG (colour terms)
  character(len=*), parameter, private :: FITS_KEY_TRAMODEL = 'TRAMODEL'

  private :: idcol, trahead, tradatatab

contains

  subroutine traload(nametable,photosys_ref,photosys_instr, &
       tr,trerr,tr1,tr1err,status,filters)

    ! filters (optional) are the names of the columns of the forward
    ! matrix, the filters of the transformation

    character(len=*), intent(in) :: nametable
    character(len=*), intent(out) :: photosys_ref,photosys_instr
    real(dbl), dimension(:,:), allocatable, intent(out) :: tr,trerr,tr1,tr1err
    integer, intent(out) :: status
    character(len=FLEN_VALUE), dimension(:), allocatable, intent(out), &
         optional :: filters

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
    if( present(filters) ) then
       allocate(filters(ncols))
       do i = 1, ncols
          call fits_make_keyn('TTYPE',i,key,status)
          call fits_read_key(fits,key,filters(i),status)
       end do
    end if

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


  subroutine tramodel(nametable,model,status)

    ! The model of the transformation table: TRAMODEL of the primary
    ! header, FLUX (the linear mixing of photons, the forward matrix;
    ! the tables without the keyword) or MAG (the colour terms of
    ! the magnitudes)

    character(len=*), intent(in) :: nametable
    character(len=*), intent(out) :: model
    integer, intent(out) :: status

    type(fitsfiles) :: fits

    status = 0
    model = 'FLUX'
    call fits_open_file(fits,nametable,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read a table in the file `',trim(nametable),"'."
       return
    end if
    call fits_read_key(fits,FITS_KEY_TRAMODEL,model,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       model = 'FLUX'
       status = 0
    end if
    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)
    if( .not. (model == 'FLUX' .or. model == 'MAG') ) then
       write(error_unit,*) "Error: unknown model of the transformation `", &
            trim(model),"'."
       status = 1
    end if

  end subroutine tramodel


  subroutine traload_mag(nametable,photosys_ref,photosys_instr,filters, &
       cola,colb,cref,zp,dzp,k1,dk1,k2,dk2,status)

    ! The magnitude transformation (TRAMODEL = 'MAG'): for every filter
    ! m_std - m_inst = ZP + K1 (c - CREF) + K2 (c - CREF)**2,
    ! c = m_std(COLOR_A) - m_std(COLOR_B)

    character(len=*), intent(in) :: nametable
    character(len=*), intent(out) :: photosys_ref,photosys_instr
    character(len=FLEN_VALUE), dimension(:), allocatable, intent(out) :: &
         filters, cola, colb
    real(dbl), dimension(:), allocatable, intent(out) :: cref,zp,dzp,k1,dk1,k2,dk2
    integer, intent(out) :: status

    character(len=FLEN_KEYWORD) :: key
    integer, parameter :: frow = 1, extver = 0
    real(dbl), parameter :: nullval = 0.0_dbl
    character(len=FLEN_VALUE), dimension(3) :: scols
    character(len=FLEN_VALUE), dimension(7) :: dcols
    integer :: nrows, i, col
    logical :: anyf
    type(fitsfiles) :: fits

    status = 0
    call fits_open_file(fits,nametable,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read a table in the file `',trim(nametable),"'."
       return
    end if

    call fits_make_keyn(FITS_KEY_PHOTSYS,1,key,status)
    call fits_read_key(fits,key,photosys_instr,status)
    call fits_make_keyn(FITS_KEY_PHOTSYS,2,key,status)
    call fits_read_key(fits,key,photosys_ref,status)

    call fits_movnam_hdu(fits,FITS_BINARY_TBL,trim(FTHDUNAME)//'_MAG',extver,status)
    call fits_get_num_rows(fits,nrows,status)
    if( status /= 0 ) goto 666

    allocate(filters(nrows),cola(nrows),colb(nrows),cref(nrows),zp(nrows), &
         dzp(nrows),k1(nrows),dk1(nrows),k2(nrows),dk2(nrows))

    scols = [ character(len=FLEN_VALUE) :: 'FILTER', 'COLOR_A', 'COLOR_B' ]
    dcols = [ character(len=FLEN_VALUE) :: 'CREF', 'ZP', 'ZPERR', 'K1', &
         'K1ERR', 'K2', 'K2ERR' ]
    do i = 1, size(scols)
       call fits_get_colnum(fits,.true.,scols(i),col,status)
       select case(i)
       case(1)
          call fits_read_col(fits,col,frow,'',filters,anyf,status)
       case(2)
          call fits_read_col(fits,col,frow,'',cola,anyf,status)
       case(3)
          call fits_read_col(fits,col,frow,'',colb,anyf,status)
       end select
    end do
    do i = 1, size(dcols)
       call fits_get_colnum(fits,.true.,dcols(i),col,status)
       select case(i)
       case(1)
          call fits_read_col(fits,col,frow,nullval,cref,anyf,status)
       case(2)
          call fits_read_col(fits,col,frow,nullval,zp,anyf,status)
       case(3)
          call fits_read_col(fits,col,frow,nullval,dzp,anyf,status)
       case(4)
          call fits_read_col(fits,col,frow,nullval,k1,anyf,status)
       case(5)
          call fits_read_col(fits,col,frow,nullval,dk1,anyf,status)
       case(6)
          call fits_read_col(fits,col,frow,nullval,k2,anyf,status)
       case(7)
          call fits_read_col(fits,col,frow,nullval,dk2,anyf,status)
       end select
    end do

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine traload_mag


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
    integer :: n,nrows,ncols,status
    type(fitsfiles) :: fits

    status = 0
    call trahead(fits,output,'FLUX',photosys_ref,photosys_instr,filters_instr, &
         status,filenames,airmass,ctph,extin)

    ! mean colour index
    call fits_write_key(fits,'REFCI',ci0,-6,'reference colour index',status)

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

    call tradatatab(fits,filters_ref,filters_instr,ra,dec,ph,dph,cts,dcts,status)

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine trawrite


  subroutine trawrite_mag(output,photosys_ref,photosys_instr,filters_ref, &
       filters_instr,cola,colb,cref,zp,dzp,k1,dk1,k2,dk2,sys,nstars,rms, &
       ra,dec,ph,dph,cts,dcts,filenames,airmass,ctph,extin)

    ! the magnitude transformation (TRAMODEL = 'MAG', see traload_mag);
    ! sys the systematic errors, nstars the stars, rms the scatter of
    ! the residuals [mag] of the filters

    character(len=*), intent(in) :: output, photosys_ref, photosys_instr
    character(len=*), dimension(:), intent(in) :: filters_ref, filters_instr, &
         cola, colb
    character(len=*), dimension(:), intent(in), optional :: filenames
    real(dbl), dimension(:), intent(in) :: cref,zp,dzp,k1,dk1,k2,dk2,sys,rms
    integer, dimension(:), intent(in) :: nstars
    real(dbl), dimension(:), intent(in) :: ra,dec
    real(dbl), dimension(:), intent(in), optional :: airmass,ctph,extin
    real(dbl), dimension(:,:), intent(in) :: ph,dph,cts,dcts

    integer, parameter :: frow = 1, ncols = 13
    character(len=FLEN_VALUE), dimension(ncols) :: ttype, tform,tunit
    integer :: nrows,status
    type(fitsfiles) :: fits

    status = 0
    call trahead(fits,output,'MAG',photosys_ref,photosys_instr,filters_instr, &
         status,filenames,airmass,ctph,extin)

    nrows = size(filters_ref)
    ttype = [ character(len=FLEN_VALUE) :: 'FILTER', 'COLOR_A', 'COLOR_B', &
         'CREF', 'ZP', 'ZPERR', 'K1', 'K1ERR', 'K2', 'K2ERR', 'SYSERR', &
         'RMS', 'NSTARS' ]
    tform = '1D'
    tform(1:3) = '16A'
    tform(13) = '1J'
    tunit = 'mag'
    tunit(1:3) = ''
    tunit(9:10) = 'mag/mag**2'
    tunit(7:8) = 'mag/mag'
    tunit(13) = ''
    call fits_insert_btbl(fits,nrows,ttype,tform,tunit,trim(FTHDUNAME)//'_MAG',status)
    call fits_write_comment(fits,'Magnitude transformation of every filter:',status)
    call fits_write_comment(fits,' m_std - m_inst = ZP + K1*(c - CREF) + K2*(c - CREF)**2,',status)
    call fits_write_comment(fits,' c = m_std(COLOR_A) - m_std(COLOR_B), m_inst = -2.5 log10(counts/s/m2).',status)
    call fits_write_col(fits,1,frow,filters_ref,status)
    call fits_write_col(fits,2,frow,cola,status)
    call fits_write_col(fits,3,frow,colb,status)
    call fits_write_col(fits,4,frow,cref,status)
    call fits_write_col(fits,5,frow,zp,status)
    call fits_write_col(fits,6,frow,dzp,status)
    call fits_write_col(fits,7,frow,k1,status)
    call fits_write_col(fits,8,frow,dk1,status)
    call fits_write_col(fits,9,frow,k2,status)
    call fits_write_col(fits,10,frow,dk2,status)
    call fits_write_col(fits,11,frow,sys,status)
    call fits_write_col(fits,12,frow,rms,status)
    call fits_write_col(fits,13,frow,nstars,status)

    call tradatatab(fits,filters_ref,filters_instr,ra,dec,ph,dph,cts,dcts,status)

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine trawrite_mag


  subroutine trahead(fits,output,model,photosys_ref,photosys_instr, &
       filters_instr,status,filenames,airmass,ctph,extin)

    ! the primary HDU of a transformation table

    type(fitsfiles), intent(out) :: fits
    character(len=*), intent(in) :: output, model, photosys_ref, photosys_instr
    character(len=*), dimension(:), intent(in) :: filters_instr
    integer, intent(in out) :: status
    character(len=*), dimension(:), intent(in), optional :: filenames
    real(dbl), dimension(:), intent(in), optional :: airmass,ctph,extin

    character(len=FLEN_KEYWORD) :: key
    character(len=FLEN_CARD) :: buf,a
    real(dbl) :: x
    integer :: n

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

    if( model == 'MAG' ) then
       call fits_write_key(fits,FITS_KEY_TRAMODEL,model, &
            'transformation: colour terms of magnitudes',status)
    else
       call fits_write_key(fits,FITS_KEY_TRAMODEL,model, &
            'transformation: linear mixing of photon rates',status)
    end if

    if( present(airmass) ) then
       x = sum(airmass)/size(airmass)
       call fits_write_key(fits,FITS_KEY_AIRMASS,x,-6,'average air-mass',status)
    end if

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

  end subroutine trahead


  subroutine tradatatab(fits,filters_ref,filters_instr,ra,dec,ph,dph,cts,dcts,status)

    ! the data of the fit (read by tradata)

    type(fitsfiles), intent(in out) :: fits
    character(len=*), dimension(:), intent(in) :: filters_ref, filters_instr
    real(dbl), dimension(:), intent(in) :: ra,dec
    real(dbl), dimension(:,:), intent(in) :: ph,dph,cts,dcts
    integer, intent(in out) :: status

    integer, parameter :: frow = 1
    character(len=FLEN_VALUE), dimension(:), allocatable :: ttype, tform,tunit
    integer :: n,nrows,ncols,nbands

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

  end subroutine tradatatab

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

    ! status is intent(out): undefined on entry (it was not reset,
    ! so the search failed at random: 'CTS or PH columns badly identified')
    status = 0
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
