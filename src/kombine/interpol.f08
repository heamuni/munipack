!
! Interpolation routines suitable for kombine
!
! Copyright © 2018 F.Hroch (hroch@physics.muni.cz)
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
! All interpolation methods are implemented according to wiki pages:
!
!  https://en.wikipedia.org/wiki/Bicubic_interpolation
!  https://en.wikipedia.org/wiki/Bilinear_interpolation
!

module interpol

  implicit none


contains

  real function bilinear(x,y,f)

    ! https://en.wikipedia.org/wiki/Bilinear_interpolation

    real, intent(in) :: x,y
    real, dimension(:,:), intent(in) :: f
    real, dimension(2) :: dx, dy
    real, dimension(2,2) :: M

    dx = [1 - x, x]
    dy = [1 - y, y]
    M(1,:) = [ f(1,1), f(1,2) ]
    M(2,:) = [ f(2,1), f(2,2) ]

    dy = matmul(M,dy)
    bilinear = sum(dx*dy)

  end function bilinear


  real function bicubic(x,y,f,fx,fy,fxy)

    ! https://en.wikipedia.org/wiki/Bicubic_interpolation

    real, intent(in) :: x,y
    real, dimension(0:,0:), intent(in) :: f,fx,fy,fxy

    real, dimension(4,4) :: A,M,D
    real, dimension(4,4) :: B
    real, dimension(4) :: dx, dy

    ! interpolation matrix
    M(1,:) = [ 1, 0, 0, 0 ]
    M(2,:) = [ 0, 0, 1, 0 ]
    M(3,:) = [-3, 3,-2,-1 ]
    M(4,:) = [ 2,-2, 1, 1 ]

    ! functions and derivations matrix
    D(1,:) = [  f(0,0),  f(0,1),  fy(0,0),  fy(0,1) ]
    D(2,:) = [  f(1,0),  f(1,1),  fy(1,0),  fy(1,1) ]
    D(3,:) = [ fx(0,0), fx(0,1), fxy(0,0), fxy(0,1) ]
    D(4,:) = [ fx(1,0), fx(1,1), fxy(1,0), fxy(1,1) ]

    B = matmul(D,M)
    A = matmul(transpose(M),B)

    B = matmul(D,transpose(M))
    A = matmul(M,B)


    dx = [ 1.0, x, x**2, x**3 ]
    dy = [ 1.0, y, y**2, y**3 ]

    dy = matmul(A,dy)
    bicubic = sum(dx*dy)

  end function bicubic

  subroutine diff(f,f1,u,v)

    ! calculate first order differences along x,y-directions

    real, dimension(:,:), intent(in) :: f
    real, dimension(:,:), intent(out) :: f1
    integer, intent(in) :: u,v

    real, parameter :: h = 1 ! step
    integer :: i,j

    f1 = 0
    do i = 2, size(f,1) - 1
       do j = 2, size(f,2) - 1
          f1(i,j) = (f(i+u,j+v) - f(i-u,j-v)) / (2*h)
       end do
    end do

  end subroutine diff

  subroutine diff2(f,f2)

    ! second order mixed differences. both xy-directions

    real, dimension(:,:), intent(in) :: f
    real, dimension(:,:), intent(out) :: f2

    real, parameter :: h = 1
    integer :: i,j

    f2 = 0
    do i = 2, size(f,1) - 1
       do j = 2, size(f,2) - 1
          f2(i,j) = (f(i+1,j+1) - f(i+1,j-1) - f(i-1,j+1) + f(i-1,j-1)) / (4*h**2)
       end do
    end do


  end subroutine diff2


  real function bi3conv(x,y,f)

    real, intent(in) :: x,y
    real, dimension(-1:,-1:), intent(in) :: f

    ! bi-cubic interpolation by convolution
    ! https://en.wikipedia.org/wiki/Bicubic_interpolation
    ! Section: Bicubic convolution algorithm

    real, dimension(-1:2) :: b
    integer :: i

    do i = -1,2
       b(i) = bi3conv_p(x,f(-1:2,i))
    end do

    bi3conv = bi3conv_p(y,b)

  end function bi3conv

  real function bi3conv_p(x,f)

    real, intent(in) :: x
    real, dimension(-1:),intent(in) :: f

    ! bi-cubic interpolation by convolution
    ! https://en.wikipedia.org/wiki/Bicubic_interpolation
    ! Section: Bicubic convolution algorithm

    real, dimension(4,4) :: M
    real, dimension(4) :: dx, d

    ! interpolation matrix
    M(1,:) = [ 0, 2, 0, 0 ]
    M(2,:) = [-1, 0, 1, 0 ]
    M(3,:) = [ 2,-5, 4,-1 ]
    M(4,:) = [-1, 3,-3, 1 ]

    dx = [ 1.0, x, x**2, x**3 ]

    d = matmul(M,f(-1:2))

    bi3conv_p = 0.5 * sum(dx*d)

  end function bi3conv_p

end module interpol
