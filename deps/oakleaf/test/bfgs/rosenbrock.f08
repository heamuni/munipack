!
!  Test for BFSG minimiser on Rosenbrock's function
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

module rosenbrock_line

  use bfgs_method
  use iso_fortran_env
  implicit none

  type, extends(bfgs_function) :: rosenbrock_function
   contains
     procedure :: val => rvalue
     procedure :: grad => rgradient
     procedure :: hess => rhessian
  end type rosenbrock_function

contains

  function rvalue(this,x) result(fval)
    class(rosenbrock_function) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64) :: fval

    fval = 100*(x(2) - x(1)**2)**2 + (1 - x(1))**2

  end function rvalue

  subroutine rgradient(this,x,fgrad)
    class(rosenbrock_function) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), intent(out) :: fgrad

    fgrad(1) =-400*(x(2) - x(1)**2)*x(1) - 2*(1 - x(1))
    fgrad(2) = 200*(x(2) - x(1)**2)

  end subroutine rgradient

  subroutine rhessian(this,x,fhess)
    class(rosenbrock_function) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:,:), intent(out) :: fhess

    fhess(1,1) = 800*x(1)**2 - 400*(x(2) - x(1)**2) + 2
    fhess(1,2) = -400*x(1)
    fhess(2,1) = -400*x(1)
    fhess(2,2) = 200

  end subroutine rhessian

end module rosenbrock_line


program rosenbrock

  use rosenbrock_line
  use bfgs_method
  use iso_fortran_env
  implicit none

  real(REAL64), parameter :: range = 2
  real(REAL64), parameter :: xmin = -2
  real(REAL64), parameter :: xmax = 2
  real(REAL64), parameter :: ymin = -1
  real(REAL64), parameter :: ymax = 2

  type(rosenbrock_function) :: fcn
  real(REAL64), dimension(2) :: x
  real(REAL64), dimension(2,2) :: H
  real(REAL64) :: gtol
  integer :: iflag

  x = [-1.2, 1.0 ]
!  x = [ -1 , 1 ]
!  x = 0
!  x = 666
  gtol = 0

  call bfgsmin(fcn,x,H,gtol,6666,.true.,'/var/tmp/rosenbrock_step',iflag)
  write(*,*) 'x=',x,'    iflag=',iflag
  write(*,*) 'H=',H

  open(1,file='/var/tmp/rosenbrock.par')
  write(1,*) 'XMIN = ',xmin
  write(1,*) 'XMAX = ',xmax
  write(1,*) 'YMIN = ',ymin
  write(1,*) 'YMAX = ',ymax
  close(1)

  call graph(fcn,xmin,xmax,ymin,ymax,'/var/tmp/rosenbrock_fcn')

contains

  subroutine graph(fcn,x0,x1,y0,y1,filename)

    class(rosenbrock_function), intent(in) :: fcn
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
          write(1,*) x,y,fcn%val([x,y])
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph


end program rosenbrock
