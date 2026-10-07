!
!  Fourier transformation, very simple implementation for images
!
!
!  Copyright © 2016-2026 F.Hroch (hroch@physics.muni.cz)
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
!
!
!
! https://en.wikipedia.org/wiki/Discrete_Fourier_transform
!

module ftransform

  use iso_fortran_env
  implicit none

  real(REAL64), parameter, private :: pi = 3.1415926535897931_REAL64
  real(REAL64), parameter, private :: pi2 = 2*pi

contains

  pure subroutine ftra(x,z,m,w)

    ! for forward w = 1
    ! for backward w = -1

    complex(REAL64), dimension(:,:), intent(in) :: x
    complex(REAL64), dimension(:,:), intent(out) :: z
    integer, intent(in) :: m,w

    complex(REAL64), dimension(:,:), allocatable :: y
    integer :: i,j

    allocate(y(size(x,1),size(x,2)))

    do i = 1,size(x,1)
       y(i,:) = x(i,:)
       call fft(y(i,:),m,w)
    end do

    do j = 1,size(x,2)
       z(:,j) = y(:,j)
       call fft(z(:,j),m,w)
    end do

    deallocate(y)

  end subroutine ftra



  pure subroutine fft(a,m,d)

    ! This subroutine computes Fast Fourier transformation
    ! of the input/output array a with dimension a(2**m).
    ! d setup any of forward (d=1) and backward (d=-1) transformation.

    ! The code has been adopted from book in Czech (P.Prikryl: Numericke
    ! methody matematicke analyzy, SNTL, Prague 1985) without a reference.
    ! Numerical recipes introduces very similar procedure (in C++)
    ! with only reference to autor N.M. Brenner. See also:
    !    https://en.wikipedia.org/wiki/Fast_Fourier_transform

    complex(REAL64), dimension(:), intent(in out) :: a
    integer, intent(in) :: m,d

    complex(REAL64) :: u,w,t
    integer :: n,i,j,k,l,le,le1,ip
    real(REAL64) :: ang

    n = 2**m
    j = 1
    do i = 1,n-1
       if( i < j )then
          t = a(j)
          a(j) = a(i)
          a(i) = t
        endif
        k = n/2
        do while( k < j )
           j = j-k
           k = k/2
        end do
        j = j+k
     enddo

     do l = 1,m
        le = 2**l
        le1 = le/2
        u = cmplx(1.0_REAL64,0.0_REAL64,REAL64)
        ang = -d*pi/le1
        w = cmplx(cos(ang),-sin(ang),REAL64)
        do j = 1,le1
          do i = j,n,le
             ip = i+le1
             t = a(ip)*u
             a(ip) = a(i)-t
             a(i) = a(i)+t
          enddo
          u = u*w
        enddo
      enddo
    end subroutine fft


end module ftransform
