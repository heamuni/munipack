!
!  Interpolation visual test
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

! FCFLAGS="-Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace"
! gfortran $FCFLAGS -I ../lib testinterpol.f08 -L../lib -lfitsio -lcfitsio -lm


program testinterpol

  use interpol
  use titsio

  ! test different types of interpolation
  ! images like:
  ! https://matplotlib.org/gallery/images_contours_and_fields/interpolation_methods.html

  integer, parameter :: ncount = 70
  integer, parameter :: zoom = 10
  real, parameter :: room = zoom
  integer, parameter :: nsize = zoom*ncount
  real, dimension(0:ncount+1,0:ncount+1) :: q,fx,fy,fxy
  real, dimension(nsize,nsize) :: w
  type(fitsfiles) :: fits

  integer :: n,m,i,j,i1,i2,j1,j2,status
  real :: x,y
!  integer :: n
!  integer, dimension(:),allocatable :: mm

!  call random_seed(size=n)
!  allocate(mm(n))
!  call random_seed(get=mm)
!  write(*,*) n
!  write(*,*) mm
!  mm = 666
!  call random_seed(put=mm)


  do n = 1, ncount
     do m = 1, ncount
        call random_number(x)
        q(n,m) = x
     end do
  end do
  ! padding by mean
  q(0,:) = 0
  q(:,0) = 0
  q(ncount+1,:) = 0
  q(:,ncount+1) = 0

  do i = 1, ncount
     i1 = (i-1)*zoom + 1
     i2 = i*zoom
     do j = 1, ncount
        j1 = (j-1)*zoom + 1
        j2 = j*zoom
        w(i1:i2,j1:j2) = q(i,j)
     end do
  end do

  status = 0
  call fits_create_file(fits,'none.fits',status)
  call fits_insert_img(fits,-32,2,[nsize,nsize],status)
  call fits_write_image(fits,w,status)
  call fits_close_file(fits,status)

  ! nearest neighbourhood
  w = 1
  do i = 1, nsize
     do j = 1, nsize
        x = i / room + 0.5
        y = j / room + 0.5
        n = nint(x)
        m = nint(y)
        w(i,j) = q(n,m)
     end do
  end do

  status = 0
  call fits_create_file(fits,'nearest.fits',status)
  call fits_insert_img(fits,-32,2,[nsize,nsize],status)
  call fits_write_image(fits,w,status)
  call fits_close_file(fits,status)

  ! bi-linear interpolation
  w = 1
  do i = 1, nsize
     do j = 1, nsize
        x = i / room + 0.5
        y = j / room + 0.5
        i1 = int(x)
        i2 = i1 + 1
        j1 = int(y)
        j2 = j1 + 1
        w(i,j) = bilinear(x-i1,y-j1,q(i1:i2,j1:j2))
     end do
  end do

  call fits_create_file(fits,'bilinear.fits',status)
  call fits_insert_img(fits,-32,2,[nsize,nsize],status)
  call fits_write_2d(fits,w,status)
  call fits_close_file(fits,status)

  ! bi-cubic interpolation
  w = 1

  ! prepare first and second order difference
  call diff(q,fx,1,0)
  call diff(q,fy,0,1)
  call diff2(q,fxy)

  do i = 1, nsize
     do j = 1, nsize
        x = i / room + 0.5
        y = j / room + 0.5
        i1 = int(x)
        i2 = i1 + 1
        j1 = int(y)
        j2 = j1 + 1
        w(i,j) = bicubic(x-i1,y-j1,q(i1:i2,j1:j2),fx(i1:i2,j1:j2), &
             fy(i1:i2,j1:j2),fxy(i1:i2,j1:j2))
     end do
  end do

  call fits_create_file(fits,'bicubic.fits',status)
  call fits_insert_img(fits,-32,2,[nsize,nsize],status)
  call fits_write_2d(fits,w,status)
  call fits_close_file(fits,status)


  ! bi-cubic convolution interpolation
  w = 1
  do i = 1, nsize
     do j = 1, nsize
        x = i / room + 0.5
        y = j / room + 0.5
        i1 = int(x)
        i2 = i1 + 3
        j1 = int(y)
        j2 = j1 + 3
        if( 1 < i1 .and. i2 <= ncount .and. 1 < j1 .and. j2 <= ncount ) &
             w(i,j) = bi3conv(x-i1,y-j1,q(i1:i2,j1:j2))
     end do
  end do

  call fits_create_file(fits,'bi3conv.fits',status)
  call fits_insert_img(fits,-32,2,[nsize,nsize],status)
  call fits_write_2d(fits,w,status)
  call fits_close_file(fits,status)


end program testinterpol
