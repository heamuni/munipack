!
!  FITS tool for kombine
!
!  Copyright © 2018-26 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Munipack.
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

module fitskombi

  use astrotrafo
  use titsio
  use iso_fortran_env

  implicit none

  type :: KombiFits
     character(len=FLEN_FILENAME) :: filename = ''
     character(len=FLEN_VALUE) :: dateobs = '', filter = '', object = '', code = '',&
          radesys = ''
     real(REAL64), dimension(2) :: crval, crpix, crerr
     real(REAL64), dimension(2,2) :: cd
     real(REAL64) :: longitude, latitude, scale, rot, refl, epoch
     real, dimension(:,:), allocatable :: image
     integer, dimension(2) :: naxes = [ 0, 0 ]
     type(AstroWCS) :: wcs
     integer :: bitpix = -32
     integer :: naxis = 0
     real :: exptime = 1
     logical :: background = .true.
     logical :: geo = .false.
     logical :: status = .false.
   contains
     procedure :: Load,loadimg,reproject,midjd
  end type KombiFits

contains

  subroutine Load(fits,filename,fitskeys,background,verbose)

    ! opens a FITS file, only metadata are probed,
    ! %image components remains undefined

    class(KombiFits) :: fits
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: fitskeys
    logical, intent(in) :: background, verbose

    type(fitsfiles) :: fitsfile
    integer :: status

    fits%filename = filename
    fits%background = background

    status = 0
    call fits_open_image(fitsfile,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) 'Error: failed to read the file `',trim(filename),"'."
       fits%status = .false.
       return
    end if

    call fits_get_img_type(fitsfile,fits%bitpix,status)
    call fits_get_img_dim(fitsfile,fits%naxis,status)
    if( status /= 0 ) goto 666

    if( fits%naxis /= 2 ) then
       write(error_unit,*) &
            'Error: the assertion NAXIS /= 2 failed: `',trim(filename),"."
       goto 666
    end if

    call fits_get_img_size(fitsfile,2,fits%naxes,status)
    if( status /= 0 ) goto 666

    ! read date and time, dateobs
    call fits_write_errmark
    call fits_get_dateobs(fitsfile,fitskeys(1:2),fits%dateobs,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       if( verbose ) write(error_unit,*) "Warning: `",trim(filename), &
            "`: Date of observation by `",trim(fitskeys(1)),"',`",trim(fitskeys(2)),&
            "' keywords not found."
       status = 0
       fits%dateobs = ''
    end if

    ! exposure time
    call fits_write_errmark
    call fits_read_key(fitsfile,fitskeys(3),fits%exptime,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       if( verbose ) write(error_unit,*) "Warning: `",trim(filename), &
            "': Exposure time by `",trim(fitskeys(3)),"' keyword not found."
       status = 0
       fits%exptime = 1
    end if

    ! filter
    call fits_write_errmark
    call fits_read_key(fitsfile,fitskeys(4),fits%filter,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       if( verbose ) write(error_unit,*) "Warning: `",trim(filename), &
            "': Filter by `",trim(fitskeys(4)),"' keyword not found."
       status = 0
       fits%filter = ''
    end if

    ! object
    call fits_write_errmark
    call fits_read_key(fitsfile,fitskeys(5),fits%object,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       if( verbose ) write(error_unit,*) "Warning: `",trim(filename), &
            "': Object name by `",trim(fitskeys(5)),"' keyword not found."
       status = 0
       fits%object = ''
    end if

    ! Geographic coordinates
    call fits_write_errmark
    call fits_read_key(fitsfile,fitskeys(6),fits%longitude,status)
    call fits_read_key(fitsfile,fitskeys(7),fits%latitude,status)
    fits%geo = status == 0
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       if( verbose ) write(error_unit,*) &
            "Warning: `",trim(filename),"': Geographics coordinates by ", &
            trim(fitskeys(6)),' and ',trim(fitskeys(7)),' keyword not found.'
       status = 0
    end if

    ! get WCS calibration
    block
      character(len=:), allocatable :: header
      character(len=FLEN_VALUE), dimension(0) :: exclist
      integer :: nkeys

      call fits_hdr2str(fitsfile,1,exclist,nkeys,header,status)
      if( status == 0 ) then
         associate( wcs => fits % wcs )
           call wcs % init(header,nkeys,status)
           if( status == 0 ) then
              call wcs % get_code(fits % code)
              call wcs % get_affines(fits % scale, fits % rot, fits % refl)
              call wcs % get_crval(fits % crval)
              call wcs % get_crpix(fits % crpix)
              call wcs % get_radesys(fits % radesys)
              call wcs % get_epoch(fits % epoch)
           end if
         end associate
      end if
      deallocate(header)
    end block

    if( status /= 0 ) then
       if( verbose ) write(error_unit,*) &
            "Error: `",trim(filename),"': astrometry calibration not found."
       goto 666
    end if

666 continue

    call fits_close_file(fitsfile,status)
    call fits_report_error(error_unit,status)

    fits%status = status == 0

  end subroutine Load

  logical function loadimg(this)

    ! fills %image structure item

    class(KombiFits) :: this
    character(len=80) :: msg
    integer :: status
    logical :: anyf
    type(fitsfiles) :: fitsfile

    allocate(this%image(this%naxes(1),this%naxes(2)),stat=status,errmsg=msg)
    if( status /= 0 ) then
       write(error_unit,*) "Error: ",trim(this%filename),": ",trim(msg)
       loadimg = .false.
       return
    end if

    status = 0
    call fits_open_image(fitsfile,this%filename,FITS_READONLY,status)
    if( status == 0 ) then
       call fits_read_image(fitsfile,0,0.0,this%image,anyf,status)
       call fits_close_file(fitsfile,status)
       call fits_report_error(error_unit,status)
    end if

    loadimg = allocated(this%image) .and. status == 0

  end function loadimg

  subroutine reproject(this,method,jdref,crmov,alpha,delta,img,mask)

    use interpol
    use oakleaf

    class(KombiFits) :: this
    character(len=*), intent(in) :: method
    real(REAL64), intent(in) :: jdref
    real(REAL64), dimension(:), intent(in) :: crmov
    real(REAL64), dimension(:,:), intent(in) :: alpha, delta
    real, dimension(:,:), intent(out) :: img
    logical, dimension(:,:), intent(out) :: mask

    real, dimension(:,:), allocatable :: fx,fy,fxy
    integer :: i,j,n,m,i1,j1,nx,ny
    real :: x,y,back
    real(REAL64) :: xx,yy,dt,da,dd

    nx = size(this%image,1)
    ny = size(this%image,2)

    mask = .false.

    ! proper motion: the frame is sampled at positions shifted
    ! by the motion since the reference time jdref
    dt = this % midjd() - jdref
    da = dt * crmov(1)
    dd = dt * crmov(2)

    ! background
    if( this%background ) then
       m = max(minval(this%naxes) / 1000, 1)
       call rmean(pack(this%image(::m,::m),.true.),back)
    else
       back = 0.0
    end if

    if( method == 'BICUBIC' ) then
       allocate(fx(nx,ny),fy(nx,ny),fxy(nx,ny))
       call diff(this%image,fx,1,0)
       call diff(this%image,fy,0,1)
       call diff2(this%image,fxy)
    end if

    do i = 1, size(alpha,1)
       do j = 1, size(alpha,2)

          call this % wcs % cel2pix(alpha(i,j)+da,delta(i,j)+dd,xx,yy)

          ! points outside of the frame (and NaNs) are skipped,
          ! pixel (n,m) covers the area n-0.5 .. n+0.5, m-0.5 .. m+0.5
          if( .not. (0.5 <= xx .and. xx < nx + 0.5 .and. &
               0.5 <= yy .and. yy < ny + 0.5) ) cycle

          x = real(xx)
          y = real(yy)
          i1 = int(x)
          j1 = int(y)

          ! The stencil of the interpolation must lie inside of the frame,
          ! the bilinear and the nearest neighbourhood methods are used
          ! near the edges of the frame where the stencil is out of the frame.
          ! The derivatives (bicubic) are valid on 2 .. n-1.
          if( method == 'BICUBIC' .and. &
               2 <= i1 .and. i1 + 1 < nx .and. 2 <= j1 .and. j1 + 1 < ny ) then

             img(i,j) = bicubic(x-i1,y-j1,this%image(i1:i1+1,j1:j1+1), &
                  fx(i1:i1+1,j1:j1+1),fy(i1:i1+1,j1:j1+1),fxy(i1:i1+1,j1:j1+1))

          else if( method == 'BI3CONV' .and. &
               2 <= i1 .and. i1 + 2 <= nx .and. 2 <= j1 .and. j1 + 2 <= ny ) then

             ! the stencil -1 .. 2 around i1, j1
             img(i,j) = bi3conv(x-i1,y-j1,this%image(i1-1:i1+2,j1-1:j1+2))

          else if( method /= 'NEAR' .and. &
               1 <= i1 .and. i1 < nx .and. 1 <= j1 .and. j1 < ny ) then

             img(i,j) = bilinear(x-i1,y-j1,this%image(i1:i1+1,j1:j1+1))

          else

             ! nearest neighbourhood
             n = min(max(nint(x),1),nx)
             m = min(max(nint(y),1),ny)
             img(i,j) = this%image(n,m)

          end if
          img(i,j) = img(i,j) - back
          mask(i,j) = .true.

       end do
    end do

  end subroutine reproject

  function midjd(this) result(jd)

    ! returns JD at middle of this exposure

    class(KombiFits) :: this
    integer :: status
    real(REAL64) :: jd

    status = 0
    jd = fits_jd(this%dateobs,status)
    if( status /= 0 ) write(error_unit,'(a)') &
         "Error: Failed to convert `"//trim(this%dateobs)//"' to JD."

    jd = jd + this%exptime / (2*86400)

  end function midjd

end module fitskombi
