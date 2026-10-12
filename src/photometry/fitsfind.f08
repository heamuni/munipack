!
!  fitsfind
!
!  Copyright © 2013, 2018-22 F.Hroch (hroch@physics.muni.cz)
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


module fitsfind

  use titsio
  use iso_fortran_env

  implicit none

contains

  subroutine fits_find_read(filename,fkeys,data,readns,saturation,status)

    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: fkeys
    real, intent(out) :: readns,saturation
    real, dimension(:,:), allocatable, intent(out) :: data
    integer, intent(in out) :: status

    integer, parameter :: DIM = 2
    integer :: naxis, bitpix
    integer, dimension(DIM) :: naxes
    logical :: satkey, anyf
    type(fitsfiles) :: fits

    if( status /= 0 ) return

    call fits_open_image(fits,filename,FITS_READONLY,status)
        if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read an image in the file `',trim(filename),"'."
       return
    end if

    call fits_get_img_dim(fits,naxis,status)
    if( naxis /= DIM .and. status == 0 ) then
       write(error_unit,*) 'Error in FIND: Only 2D frames are supported.'
       goto 666
    end if
    call fits_get_img_size(fits,DIM,naxes,status)
    if( status /= 0 ) goto 666

    call fits_read_key(fits,fkeys(1),saturation,status)
    satkey = status == 0
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       call fits_get_img_type(fits,bitpix,status)
       if( status /= 0 ) goto 666
    end if

    call fits_read_key(fits,fkeys(2),readns,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       readns = 0
       status = 0
    end if

    if( status /= 0 ) goto 666

    allocate(data(naxes(1),naxes(2)))
    call fits_read_image(fits,0,0.0,data,anyf,status)
    if( status /= 0 ) goto 666

    if( .not. satkey ) then
       if( bitpix > 0 ) then
          saturation = 2.0**bitpix - 1
       else
          saturation = 0.99*maxval(data)
       end if
    end if

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( status /= 0 ) then
       if( allocated(data) ) deallocate(data)
    end if

  end subroutine fits_find_read



  subroutine fits_find_save(filename,output,fkeys,nstar, &
       fwhm,threshold,shrplo,shrphi,rndlo,rndhi, readns, &
       lothresh,  lobad, hibad, hmin, skymod, skyerr, skysig, maxsky, &
       stats, falserate, status)

    use oakleaf
    use findnoise, only: findstats

    ! results fills new FITS extension

    character(len=*), intent(in) :: filename, output
    character(len=*), dimension(:), intent(in) :: fkeys

    real, intent(in) :: fwhm, threshold, &
         shrplo,shrphi,rndlo,rndhi, lothresh, readns, lobad, hibad, hmin, &
         skymod, skyerr, skysig, falserate
    type(findstats), intent(in) :: stats
    integer, intent(in) :: maxsky, nstar
    integer, intent(in out) :: status

    integer, parameter :: extver = 0
    character(len=FLEN_CARD) :: buf
    character(len=FLEN_VALUE), dimension(5) :: ttype, tform, tunit
    real, dimension(:), allocatable :: xcen,ycen,hstar,round,sharp,ecc,incl
    real :: ecc_mean, incl_mean
    integer :: n, hdutype, frow, nrows, srows, i, l
    type(fitsfiles) :: fits

    if( status /= 0 ) return

    if( output == '' ) then
       call fits_open_file(fits,filename,FITS_READWRITE,status)
       if( status /= 0 ) then
          write(error_unit,*) 'Error: failed to open the file `', &
               trim(filename),"' for a table update."
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

    ! look for an older extension
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,FINDEXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       status = 0
    else if( status == 0 ) then
       ! already presented ? remove it !
       call fits_delete_hdu(fits,hdutype,status)
    end if
    if( status /= 0 ) goto 666

    ttype(1) = FITS_COL_X
    ttype(2) = FITS_COL_Y
    ttype(3) = FITS_COL_PEAKRATIO
    ttype(4) = FITS_COL_SHARP
    ttype(5) = FITS_COL_ROUND
    tform = '1E'
    tunit = ''

    call fits_insert_btbl(fits,nstar,ttype,tform,tunit,FINDEXTNAME,status)
    call fits_update_key(fits,fkeys(1),hibad,5,'[cts] saturation',status)
    call fits_update_key(fits,fkeys(2),readns,-7,'[ADU] read noise',status)
    call fits_update_key(fits,FITS_KEY_FWHM,fwhm,-2, &
         '[pix] standard FWHM of objects',status)
    call fits_update_key(fits,FITS_KEY_THRESHOLD,threshold,-4, &
         'threshold: Normal deviate of the false-alarm rate',status)
    call fits_update_key(fits,'HTHRESH',real(stats%x),-4, &
         '[sigma] the height threshold of the filtered sky',status)
    call fits_update_key(fits,'SKEWFILT',real(stats%skew),-4, &
         'skewness of the filtered sky (noise model)',status)
    call fits_update_key(fits,'KURTFILT',real(stats%kurt),-4, &
         'excess kurtosis of the filtered sky',status)
    if( falserate > 0 ) then
       call fits_update_key(fits,'FALSERAT',falserate,-4, &
            'expected noise peaks per frame (requested)',status)
    else
       call fits_update_key(fits,'FALSERAT',real(stats%nfalse),-4, &
            'expected noise peaks per frame (upper estimate)',status)
    end if
    call fits_update_key(fits,FITS_KEY_LOWBAD,lobad,-3, &
         '[cts] low good datum',status)
    call fits_update_key(fits,FITS_KEY_HIGHBAD,hibad,-3, &
         '[cts] high good datum',status)
    call fits_update_key(fits,FITS_KEY_RNDLO,rndlo,-3,'low round',status)
    call fits_update_key(fits,FITS_KEY_RNDHI,rndhi,-3,'high round',status)
    call fits_update_key(fits,FITS_KEY_SHRPLO,shrplo,-3,'low sharp',status)
    call fits_update_key(fits,FITS_KEY_SHRPHI,shrphi,-3,'high sharp',status)

    call fits_write_comment(fits,'Star detection parameters:',status)

    write(buf,*) 'Saturation (counts)=',hibad
    call fits_write_comment(fits,buf,status)

    write(buf,*) 'Read noise (ADU)=',readns
    call fits_write_comment(fits,buf,status)

    write(buf,*) 'Lower threshold (sigma)=',lothresh
    call fits_write_comment(fits,buf,status)

    write(buf,*) 'Levels range (counts) =',lobad, '..',hibad
    call fits_write_comment(fits,buf,status)

    write(buf,*) 'Hmin (counts) =',hmin
    call fits_write_comment(fits,buf,status)

    write(buf,*) 'Round range =',rndlo, '..',rndhi
    call fits_write_comment(fits,buf,status)

    write(buf,*) 'Sharp range =',shrplo, '..',shrphi
    call fits_write_comment(fits,buf,status)

    write(buf,*) 'Approximate sky value =',skymod,'+-',skyerr
    call fits_write_comment(fits,buf,status)

    write(buf,*) 'Estimated sky sigma =',skysig
    call fits_write_comment(fits,buf,status)

    write(buf,*) 'Pixels used for sky determination =',maxsky
    call fits_write_comment(fits,buf,status)

    allocate(xcen(nstar),ycen(nstar),sharp(nstar),round(nstar),hstar(nstar),&
         ecc(nstar),incl(nstar))
    do n = 1, nstar
       read(3) xcen(n),ycen(n),hstar(n),sharp(n),round(n),ecc(n),incl(n)
    end do

    ! sort arrays by height above lower threshold
    call sorter(xcen,ycen,hstar,sharp,round)

    ! mean eccentricity and inclination
    call rmean(ecc,ecc_mean)
    call rmean(incl,incl_mean)
    call fits_update_key(fits,FITS_KEY_ECCENTRICITY,ecc_mean,-2, &
         ' mean eccentricity',status)
    call fits_update_key(fits,FITS_KEY_INCLINATION,nint(incl_mean), &
         ' mean inclination',status)

    nrows = nstar
    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_write_col(fits,1,frow,xcen(i:l),status)
       call fits_write_col(fits,2,frow,ycen(i:l),status)
       call fits_write_col(fits,3,frow,hstar(i:l),status)
       call fits_write_col(fits,4,frow,sharp(i:l),status)
       call fits_write_col(fits,5,frow,round(i:l),status)
    end do

    deallocate(xcen,ycen,hstar,round,sharp,ecc,incl)

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine fits_find_save


  subroutine sorter(xcen,ycen,hstar,sharp,round)

    use quicksort

    real, dimension(:),intent(in out) :: xcen,ycen,hstar,sharp,round

    integer, parameter :: rp = selected_real_kind(15)
    real(rp), dimension(:), allocatable :: htmp
    real, dimension(:), allocatable :: tmp
    integer, dimension(:), allocatable :: id, idx
    integer :: i,n,m

    n = size(xcen)
    allocate(tmp(n),htmp(n),id(n),idx(n))
    id = [ (i, i = 1,n) ]
    htmp = hstar

    call qsort(htmp,id)
    ! hstar sorted into low to high order

    ! reverse sort
    m = n + 1
    forall( i = 1:n ) hstar(i) = real(htmp(m - i))
    forall( i = 1:n ) idx(i) = id(m - i)

    tmp = xcen
    xcen = tmp(idx)

    tmp = ycen
    ycen = tmp(idx)

    tmp = sharp
    sharp = tmp(idx)

    tmp = round
    round = tmp(idx)

    deallocate(tmp,htmp,id,idx)

  end subroutine sorter


end module fitsfind
