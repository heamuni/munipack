!
!  Test for BFSG method fitting the straight line of Rektorys
!
!  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Oakleaf.
!
!  Oakleaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Oakleaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.

module line_type

  use bfgs_method
  use iso_fortran_env
  implicit none

  type, extends(bfgs_function) :: line_function
     real(REAL64), dimension(:), pointer :: x,y
   contains
     procedure :: val => lvalue
     procedure :: grad => lgradient
     procedure :: hess => lhessian
  end type line_function

contains


  function lvalue(this,x) result(fval)
    class(line_function) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64) :: a,b,fval
    integer :: n

    a = x(1)
    b = x(2)
    associate( x => this%x, y => this%y )
      n = size(x)
      fval = sum((y - a*x - b)**2) !/ (2*n)
    end associate

  end function lvalue

  subroutine lgradient(this,x,fgrad)
    class(line_function) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), intent(out) :: fgrad
    real(REAL64) :: a,b
    integer :: n

    a = x(1)
    b = x(2)

    associate( x => this%x, y => this%y )
      n = size(x)
      fgrad(1) = -sum( (y - a*x - b)*x ) !/ n
      fgrad(2) = -sum( (y - a*x - b) ) !/ n
    end associate

  end subroutine lgradient

  subroutine lhessian(this,x,fhess)
    class(line_function) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:,:), intent(out) :: fhess
    integer :: n
    real(REAL64) :: a,b

    a = x(1)
    b = x(2)

    associate( x => this%x, y => this%y )
      n = size(x)
      fhess(1,1) = sum(x**2) !/ n
      fhess(2,1) = sum(x) !/ n
      fhess(1,2) = fhess(2,1)
      fhess(2,2) = size(x) !/ n
    end associate

  end subroutine lhessian

end module line_type

program what_about_line

  use line_type
  use bfgs_method
  use rektorys
  use QRHouse
  use iso_fortran_env
  implicit none

  real(REAL64), parameter :: range = 2
  real(REAL64), parameter :: amin = -30
  real(REAL64), parameter :: amax = 30
  real(REAL64), parameter :: bmin = -10
  real(REAL64), parameter :: bmax = 50

  type(line_function) :: fcn
  real(REAL64), dimension(2) :: p
  real(REAL64), dimension(2,2) :: B,Binv,H,Hinv
  real(REAL64) :: gtol, s0
  integer :: iflag, n, i

  fcn%x => x
  fcn%y => y

  p = 0
  n = size(x)
  gtol = sqrt(10*size(p)*epsilon(x))
  gtol = 1e3*size(p)*epsilon(x)

  call bfgsmin(fcn,p,H,gtol,666,.true.,'/var/tmp/line_step',iflag)
  write(*,*) 'p=',p,'    iflag=',iflag
  s0 = fcn % val(p)
  write(*,*) 'S0 = ',s0, ', s=',sqrt(s0/n)

  call qrinv(H,Hinv)
  do i = 1, 2
     write(*,'(a,2f25.15)') 'Hinv: ',Hinv(i,:)
  end do

  call fcn%hess(p,B)
  do i = 1, 2
     write(*,'(a,2f25.15))') 'B: ',B(i,:)
  end do

  do i = 1,2
     write(*,*) 'stderr',i,sqrt(s0*H(i,i)/(n - 2))
  end do
  call qrinv(B,Binv)
  do i = 1,2
     write(*,*) 'stderr',i,sqrt(s0*Binv(i,i)/(n - 2)),' (analytic)'
  end do

  write(*,*) 'Bp =',matmul(B,p), '  expected: 170, 1148.78'

  open(1,file='/var/tmp/line.par')
  write(1,*) 'XMIN = ',amin
  write(1,*) 'XMAX = ',amax
  write(1,*) 'YMIN = ',bmin
  write(1,*) 'YMAX = ',bmax
  close(1)

  call graph(fcn,amin,amax,bmin,bmax,'/var/tmp/line_fcn')

contains

    subroutine graph(fcn,x0,x1,y0,y1,filename)

    class(line_function), intent(in) :: fcn
    real(REAL64), intent(in) :: x0,x1,y0,y1
    character(len=*), intent(in) :: filename
    integer, parameter :: nx = 200
    integer, parameter :: ny = 200
    real(REAL64) :: x,dx,y,dy
    integer :: i,j

    dx = (x1 - x0) / (nx - 1)
    dy = (y1 - y0) / (ny - 1)

    open(1,file=filename)
    do j = 1, ny
       y = y0 + dy*j
       do i = 0, nx
          x = x0 + dx*i
          write(1,*) x,y,fcn % val([x,y])
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph


end program what_about_line
