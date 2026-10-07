!
!  Colour FITS tool
!
!  Copyright © 2018-2022 F.Hroch (hroch@physics.muni.cz)
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

module fitscolour

  use titsio
  use iso_fortran_env

  implicit none

  type :: ColourFits
     character(len=FLEN_FILENAME) :: filename = ''
     real, dimension(:,:), allocatable :: image
     character(len=FLEN_VALUE) :: filter = ''
     character(len=FLEN_VALUE) :: dateobs = ''
     character(len=FLEN_VALUE) :: object = ''
     character(len=FLEN_VALUE) :: radesys = ''
     character(len=FLEN_VALUE), dimension(2) :: ctype
     real(REAL64), dimension(2) :: crval, crpix
     real(REAL64), dimension(2,2) :: cd
     integer, dimension(2) :: naxes = [ 0, 0 ]
     integer :: naxis = 0
     integer :: bitpix = 0
     real :: exptime = 1   ! [s]
     real :: area = 1      ! [m2]
     real :: leff = 1239.8e-9 ! [m] == 1eV
     real :: fref = 3.6e-9  ! [W/m2] reference flux
     real :: ctph = 1
     real :: scale = 1     ! degs per pixel
     real :: saturate = huge(1.0)
     real(REAL64) :: epoch
     logical :: astrometry = .false.
     logical :: calibrated = .false.
     logical :: status = .false.
   contains
     procedure :: Load, sky, aperture, apmean, apsum
  end type ColourFits

contains

  subroutine Load(fits,filename,fitskeys)

    ! opens a FITS file

    class(ColourFits) :: fits
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: fitskeys

    integer, parameter :: extver = 0
    character(len=80) :: msg
    real(REAL64), dimension(2) :: crder
    integer :: stat,status
    type(fitsfiles) :: fitsfile
    logical :: anyf

    fits%filename = filename
    fits%exptime = 1
    fits%area = 1

    status = 0
    call fits_open_image(fitsfile,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) 'Error: failed to read the file `',trim(filename),"'."
       fits%status = .false.
       return
    end if

    call fits_get_img_dim(fitsfile,fits%naxis,status)
    if( status /= 0 ) goto 666

    if( fits%naxis /= 2 .and. status == 0 ) then
       write(error_unit,*) "Error: `",trim(filename), &
            "': only 2D images are supported for colouring."
       goto 666
    end if

    call fits_get_img_type(fitsfile,fits%bitpix,status)
    call fits_get_img_size(fitsfile,2,fits%naxes,status)
    if( status /= 0 ) goto 666

    call fits_read_key(fitsfile,fitskeys(1),fits%filter,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       fits%filter = ''
    end if

    call fits_read_key(fitsfile,fitskeys(2),fits%exptime,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       fits%exptime = 1
    end if

    call fits_read_key(fitsfile,fitskeys(3),fits%area,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       fits%area = 1
    end if

    call fits_read_key(fitsfile,fitskeys(4),fits%object,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       fits%object = ''
    end if

    call fits_read_key(fitsfile,FITS_KEY_DATEOBS,fits%dateobs,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       fits%dateobs = ''
    end if

    ! add astrometry parameters
    call fits_read_wcs(fitsfile,fits%radesys,fits%epoch,fits%ctype, &
         fits%crval,fits%crpix,fits%cd,crder,status)
    fits%astrometry = status == 0
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
    else
       ! determine scale, degrees per pixel
       fits%scale = real(sqrt(fits%cd(1,1)**2 + fits%cd(1,2)**2))
    end if

    if( status /= 0 ) goto 666

    allocate(fits%image(fits%naxes(1),fits%naxes(2)),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) "Error: ",trim(msg)
       error stop 'ColourFits failed during memory allocation.'
    end if
    call fits_read_image(fitsfile,0,0.0,fits%image,anyf,status)

    ! photometry calibration
    call fits_write_errmark
    call fits_movnam_hdu(fitsfile,FITS_BINARY_TBL,PHOTOEXTNAME,extver,status)
    if( status == 0 ) then
       call fits_read_key(fitsfile,FITS_KEY_CTPH,fits%ctph,status)
       fits%calibrated = .true.
    else if( status == FITS_BAD_HDU_NUM ) then
       call fits_clear_errmark
       status = 0
    end if

    ! estimate saturate
    if( fits%bitpix > 0 ) then
       fits%saturate = 2.0**fits%bitpix - 1
    end if

666 continue

    if( status /= 0 .and. allocated(fits%image) ) deallocate(fits%image)

    call fits_close_file(fitsfile,status)
    call fits_report_error(error_unit,status)

    fits%status = status == 0

  end subroutine Load

  function sky(this)

    ! estimates background

    use oakleaf

    real :: sky
    class(ColourFits), intent(in) :: this
    integer :: dx,dy,n

    if( .not. this%status ) then
       sky = 0
       return
    end if

    n = int(sqrt(real(this%naxes(1)*this%naxes(2)) / 32768.0))
    dx = max(this%naxes(1)/n,1)
    dy = max(this%naxes(2)/n,1)
    call rmean(pack(this%image(::dx,::dy),.true.),sky)

  end function sky

  subroutine aperture(this,i0,j0,r0,back,data)

    class(ColourFits), intent(in) :: this
    integer, intent(in) :: i0, j0, r0
    real, intent(in) :: back
    real, dimension(:), allocatable, intent(out) :: data

    integer :: n,i,j,dx,dy,r2,width,height
    real, dimension(:), allocatable :: pixels

    ! init
    width = size(this%image,1)
    height = size(this%image,2)

    ! photometry in the aperture
    n = (2*r0 + 1)**2
    allocate(pixels(n))
    r2 = r0**2
    n = 0
    do j = max(j0-r0, 1), min(j0+r0, height)
       dy = (j - j0)**2
       do i = max(i0-r0, 1), min(i0+r0, width)
          dx = (i - i0)**2
          if( dx + dy <= r2 ) then
             n = n + 1
             pixels(n) = this%image(i,j) - back
          end if
       end do
    end do

    allocate(data(n))
    data = pixels(1:n)
    deallocate(pixels)

  end subroutine aperture

  function apmean(this,x,y,r,back)

    ! estimates mean flux in the aperture

    use oakleaf

    real :: apmean
    class(ColourFits), intent(in) :: this
    real, intent(in) :: x, y, r, back

    real, dimension(:), allocatable :: data

    call aperture(this,nint(x),nint(y),nint(r),back,data)
    call rmean(data,apmean)
    deallocate(data)

  end function apmean

  function apsum(this,x,y,r,back)

    ! estimates total flux in the aperture

    real :: apsum
    class(ColourFits), intent(in) :: this
    real, intent(in) :: x, y, r, back

    integer :: i0, j0, r0
    real, dimension(:), allocatable :: data

    r0 = nint(r)
    call locpeak(this,nint(x),nint(y),r0,i0,j0)
    call aperture(this,i0,j0,r0,back,data)
    apsum = sum(data)
    deallocate(data)

  end function apsum

  subroutine locpeak(this,x,y,r0,i0,j0)

    class(ColourFits), intent(in) :: this
    integer, intent(in) :: x, y, r0
    integer, intent(out) :: i0,j0

    integer :: i,j,width,height
    real :: qmax

    ! init
    width = size(this%image,1)
    height = size(this%image,2)
    i0 = x
    j0 = y

    ! locate the local maximum
    qmax = this%image(i0,j0)
    do j = max(y - r0, 1), min(y + r0, height)
       do i = max(x - r0, 1), min(x + r0, width)
          if( this%image(i,j) > qmax ) then
             qmax = this%image(i,j)
             i0 = i
             j0 = j
          end if
       end do
    end do

  end subroutine locpeak

end module fitscolour
