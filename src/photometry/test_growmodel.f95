
! gfortran -Wall -p -g -fcheck=all -I../lib  test_growmodel.f95 -L. -L../lib -ltitsio -lcfitsio -lm

program test_growmodel

  use titsio

  implicit none

  integer, parameter :: group = 1, extver = 0, frow = 1, felem = 1
  real, parameter :: nullval = 0.0
  integer, dimension(2) :: naxes, cens
  logical :: anyf
  real, dimension(:,:), allocatable :: psf
  real, dimension(12) :: raper, cts, grow
  integer :: zoom, bitpix, status, naxis, dim, i,j, m, n
  real :: f,d,b
  type(fitsfiles) :: fits

  status = 0
  call fits_open_file(fits,'art_01.fits',READONLY,status)
  if( status /= 0 ) stop 'failed to open'
  call fits_movnam_hdu(fits,IMAGE_HDU,'PSF',extver,status)
  call fits_get_img_param(fits,bitpix,naxis,naxes,status)
  call fits_get_key('ZOOM',zoom,status)
  dim = naxes(1) / 2
  allocate(psf(-dim:dim,-dim:dim))
  call fits_read_image(fits,psf,status)
  call fits_close_file(fits,status)
  call fits_report_error(error_unit,status)

  forall( i = 1:size(raper) ) raper(i) = exp(0.17338*((i-1)*1.570796327))

  cens = naxes / 2 + 1

  do m = 1, size(raper)
     n = nint(zoom * raper(m))
     f = 0
     do i = -n,n
        do j = -n,n
           if( ((i-zoom/10)**2 + j**2 <= n**2) .and. (-dim < i .and. i < dim) .and. (-dim < j .and. j < dim) ) then
              f = f + psf(i,j)
           end if
        end do
     end do
     cts(m) = f/zoom**2
     write(*,*) raper(m),real(f/zoom**2)
  end do

  deallocate(psf)

  grow(size(raper)) = 1
  b = -0.0000314*0.1
  b = 0
  do i = size(raper)-1,1,-1
     grow(i) = grow(i+1)*(cts(i)-b*raper(i)**2)/(cts(i+1)-b*raper(i+1)**2)
  end do
  do i = 1,size(raper)
     write(*,*) raper(i), cts(i), grow(i)
!     write(*,*) raper(i)/raper(size(raper)),cts(i),grow(i)
  end do

!  d = cts(1) / raper(1)
!  write(*,*) raper(1),cts(1) + 0*d*(0.1), d*0.1
!  do m = 2,size(raper)
!     write(*,*) raper(m),cts(m) + 0*d*(0.1), d
!     d = (cts(m) - cts(m-1)) / (raper(m) - raper(m-1))
!  end do


end program test_growmodel
