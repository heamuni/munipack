!
! trispace - create histogram from triangles of all possible permutations
!
! $ gfortran -Wall -g -p -fbounds-check -fimplicit-none  -I../lib match.f95  astrofits.f95 trispace.f95 -lcfitsio -lm
!
! Copyright © 2012 F.Hroch (hroch@physics.muni.cz)
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

program trispace

  use astrofits
  use matcher

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)
  integer, parameter :: uvdim = 100

  real(dbl), allocatable, dimension(:) :: x, y, x1, y1
  integer, dimension(:), allocatable :: id1
  real(dbl) :: w,h, u, v
  integer :: i,j,k,nmatch,iu,iv
  real(dbl),dimension(2) :: crpix
  integer, allocatable, dimension(:,:) :: uvhist

  call readfile('barnard_01R.fits',w,h,x,y,crpix,i)

  nmatch = size(x)
  allocate(id1(nmatch),x1(nmatch),y1(nmatch))
  x1 = x(1:nmatch)
  y1 = y(1:nmatch)

  allocate(uvhist(0:uvdim,0:uvdim))
  uvhist = 0

  do i = 1, nmatch
     do j = i+1, nmatch - 1
        do k = j+1, nmatch - 1
           call triangle(x1(i),y1(i),x1(j),y1(j),x1(k),y1(k),u,v)
!           write(*,*) i,j,k,u,v
           iu = int(uvdim*u)
           iv = int(uvdim*v)
!           write(*,*) u,v,iu,iv
           uvhist(iu,iv) = uvhist(iu,iv) + 1
        end do
     end do
  end do


  open(16,file='uvhist')
  do i = 0,uvdim
     do j = 0,uvdim
        write(16,'(2f6.3,1x,i0)') real(i)/real(uvdim),real(j)/real(uvdim),uvhist(i,j)
     end do
     write(16,*)
  end do
  close(16)

end program trispace

