!
!  fitsgphot
!
!  Copyright © 2018 F.Hroch (hroch@physics.muni.cz)
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


module fitsgphot

  use titsio
  use iso_fortran_env

  implicit none

  integer, parameter, private :: rp = selected_real_kind(15)

contains

  subroutine fits_gphot_read(filename,hwhm,sep,raper,xcens,ycens, &
    sky,skyerr,skycorr,skyerrcorr,apcts,apcts_err,status, &
    ring,ecc,fwhm,naxes,errgain)

    ! ring: the sky annulus, ecc: the eccentricity of the apertures,
    ! fwhm: FWHM of FIND (-1 when unknown), naxes: the frame size,
    ! errgain: the gain of the photon noise in APCOUNTERR (APERRGN,
    ! -1 when unknown)

    character(len=*), intent(in) :: filename
    real, intent(out) :: hwhm, sep
    real, dimension(2), intent(out) :: ring
    real, intent(out) :: ecc, fwhm
    integer, dimension(2), intent(out) :: naxes
    real(rp), intent(out) :: errgain
    real, dimension(:), allocatable, intent(out) :: raper, xcens,ycens, &
         sky,skyerr,skycorr,skyerrcorr
    real, dimension(:,:), allocatable, intent(out) :: apcts,apcts_err
    integer, intent(in out) :: status

    integer :: nrows,xcol,ycol,ecol,scol,i,n,naper
    integer, parameter :: extver = 0, frow = 1
    real, parameter :: nullval = 0.0
    logical :: anyf
    character(len=FLEN_VALUE) :: key, label
    type(fitsfiles) :: fits

    status = 0
    call fits_open_file(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read the file `',trim(filename),"'."
       return
    end if

    ! the frame size
    naxes = 0
    call fits_read_key(fits,'NAXIS1',naxes(1),status)
    call fits_read_key(fits,'NAXIS2',naxes(2),status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       naxes = 0
       status = 0
    end if

    ! FWHM of the detection (the scale of the selection of the stars
    ! of the growth curve)
    fwhm = -1
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,FINDEXTNAME,extver,status)
    if( status == 0 ) then
       call fits_read_key(fits,FITS_KEY_FWHM,fwhm,status)
       if( status == FITS_KEYWORD_NOT_FOUND ) fwhm = -1
    end if
    status = 0

    call fits_movnam_hdu(fits,FITS_BINARY_TBL,APEREXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       write(error_unit,*) "Error: ",trim(APEREXTNAME)// &
            " extension not found. Has been aperture photometry by " // &
            "`munipack aphot "//trim(filename)//"' performed?"
       goto 666
    end if

    call fits_get_num_rows(fits,nrows,status)
    call fits_read_key(fits,FITS_KEY_NAPER,naper,status)
    call fits_read_key(fits,FITS_KEY_HWHM,hwhm,status)
    call fits_read_key(fits,trim(FITS_KEY_ANNULUS)//'2',sep,status)
    call fits_read_key(fits,trim(FITS_KEY_ANNULUS)//'1',ring(1),status)
    ring(2) = sep
    if( status /= 0 ) goto 666
    call fits_read_key(fits,FITS_KEY_ECCENTRICITY,ecc,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       ecc = 0
       status = 0
    end if
    ! the gain of the errors (aphot since BUGS P34), else unknown
    call fits_read_key(fits,'APERRGN',errgain,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       errgain = -1
       status = 0
    end if

    allocate(xcens(nrows),ycens(nrows),sky(nrows),skyerr(nrows), &
         skycorr(nrows),skyerrcorr(nrows),apcts(nrows,naper), &
         apcts_err(nrows,naper),raper(naper))

    do i = 1, naper
       call fits_make_keyn(FITS_KEY_APER,i,key,status)
       call fits_read_key(fits,key,raper(i),status)
    end do

    call fits_get_colnum(fits,.true.,FITS_COL_X,xcol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_Y,ycol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_SKY,scol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_SKYERR,ecol,status)
    call fits_read_col(fits,xcol,frow,nullval,xcens,anyf,status)
    call fits_read_col(fits,ycol,frow,nullval,ycens,anyf,status)
    call fits_read_col(fits,scol,frow,nullval,sky,anyf,status)
    call fits_read_col(fits,ecol,frow,nullval,skyerr,anyf,status)

    do i = 1, naper
       write(label,'(a,i0)') FITS_COL_APCOUNT,i
       call fits_get_colnum(fits,.true.,label,n,status)
       call fits_read_col(fits,n,frow,nullval,apcts(:,i),anyf,status)
       write(label,'(a,i0)') FITS_COL_APCOUNTERR,i
       call fits_get_colnum(fits,.true.,label,n,status)
       call fits_read_col(fits,n,frow,nullval,apcts_err(:,i),anyf,status)
    end do

    ! all the commands must be finished without eny error for
    ! correct head, any error indicates internal inconsistency

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( status /= 0 ) then
       if( allocated(xcens) ) deallocate(xcens,ycens,sky,skyerr, &
            skycorr,skyerrcorr,apcts,apcts_err,raper)
    end if


  end subroutine fits_gphot_read


  subroutine fits_find_save(filename,output,ghwhm,rflux90, &
       raper, xcens, ycens, sky,skyerr, gcount,gcount_err, growflag, &
       curve, curve_err, prof, status, method)

    character(len=*), intent(in) :: filename, output
    character(len=*), intent(in), optional :: method
    real, intent(in) :: ghwhm,rflux90
    real, dimension(:), intent(in) :: raper, xcens, ycens, sky,skyerr
    real(rp), dimension(:), intent(in) :: gcount,gcount_err
    integer, dimension(:), intent(in) :: growflag
    real(rp), dimension(:), intent(in) :: curve, curve_err, prof
    integer, intent(in out) :: status

    integer :: hdutype,n
    integer, parameter :: extver = 0, frow = 1, nbegin = 4
    real, parameter :: nullval = 0.0
    character(len=FLEN_VALUE), dimension(:), allocatable :: ttype, tform, tunit
    character(len=FLEN_COMMENT) :: com
    type(fitsfiles) :: fits


    if( status /= 0 ) return

    if( output == '' ) then
       call fits_open_file(fits,filename,FITS_READWRITE,status)
       if( status /= 0 ) then
          write(error_unit,*) &
            'Error: failed to open the file `',trim(filename),"' for a table update."
          return
       end if
    else
       call fits_precopy_file(fits,filename,output,FITS_READWRITE,.true.,status)
       if( status /= 0 ) then
          write(error_unit,*) &
            'Error: failed to create the file `',trim(output),"'."
          return
       end if
    end if

    ! grow photometry table
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,GROWEXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       status = 0
    else
       ! already presented ? remove it !
       call fits_delete_hdu(fits,hdutype,status)
       if( status /= 0 ) goto 666
    end if

    n = 7
    allocate(ttype(n), tform(n), tunit(n))

    tform = '1D'
    tform(7) = '1B'
    tunit = ''
    ttype(1) = FITS_COL_X
    ttype(2) = FITS_COL_Y
    ttype(3) = FITS_COL_SKY
    ttype(4) = FITS_COL_SKYERR
    ttype(5) = FITS_COL_GCOUNT
    ttype(6) = FITS_COL_GCOUNTERR
    ttype(7) = FITS_COL_GROWFLAG

    call fits_insert_btbl(fits,0,ttype,tform,tunit,GROWEXTNAME,status)
    call fits_write_key(fits,FITS_KEY_HWHM,ghwhm,-4, &
         '[pix] half width at half of maximum',status)
    call fits_write_key(fits,FITS_KEY_RF90,rflux90,-4, &
         '[pix] radius contains 90% of flux',status)
    if( present(method) ) &
         call fits_write_key(fits,'GROWMETH',method, &
         'growth curve: EE (enclosed energy), LEGACY (fit)',status)

    write(com,'(a,i0)') 'Count of stars used for curve construction: ', &
         count(growflag == 1)
    call fits_write_comment(fits,com,status)
    call fits_write_comment(fits,'GROWFLAG: 0 - star, not used,',status)
    call fits_write_comment(fits, &
         '          1 - star, used for growth curve construction,',status)
    call fits_write_comment(fits,'          2 - non-stellar object',status)

    call fits_write_col(fits,1,frow,xcens,status)
    call fits_write_col(fits,2,frow,ycens,status)
    call fits_write_col(fits,3,frow,sky,status)
    call fits_write_col(fits,4,frow,skyerr,status)
    call fits_write_col(fits,5,frow,gcount,status)
    call fits_write_col(fits,6,frow,gcount_err,status)
    call fits_write_col(fits,7,frow,growflag,status)
    deallocate(ttype,tform,tunit)

    ! store growth-curve
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,GROWCURVEXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       status = 0
    else
       ! already presented ? remove it !
       call fits_delete_hdu(fits,hdutype,status)
       if( status /= 0 ) goto 666
    end if

    n = 4
    allocate(ttype(n), tform(n), tunit(n))

    tform = '1D'
    tunit = ''
    ttype(1) = FITS_COL_R
    ttype(2) = FITS_COL_GROW
    ttype(3) = FITS_COL_GROWERR
    ttype(4) = FITS_COL_RPROF

    call fits_insert_btbl(fits,0,ttype,tform,tunit,GROWCURVEXTNAME,status)
    call fits_write_col(fits,1,frow,raper,status)
    call fits_write_col(fits,2,frow,curve,status)
    call fits_write_col(fits,3,frow,curve_err,status)
    call fits_write_col(fits,4,frow,prof,status)
    deallocate(ttype,tform,tunit)

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine fits_find_save

end module fitsgphot
