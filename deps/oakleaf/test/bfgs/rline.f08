!
!  Test for BFSG minimiser for robust regression
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


module the_function

  use huber_curvefit
  use iso_fortran_env
  implicit none

  type, extends(huberfit_function) :: line_function
   contains
     procedure :: val => lvalue
     procedure :: grad => lgradient
  end type line_function

contains

  function lvalue(this,x,p) result(fval)
    class(line_function) :: this
    real(REAL64), intent(in) :: x
    real(REAL64), dimension(:), intent(in) :: p
    real(REAL64) :: fval
    real(REAL64) :: a, b

    a = p(1)
    b = p(2)
    fval = a*x + b

  end function lvalue

  subroutine lgradient(this,x,p,fgrad)
    class(line_function) :: this
    real(REAL64), intent(in) :: x
    real(REAL64), dimension(:), intent(in) :: p
    real(REAL64), dimension(:), intent(out) :: fgrad

    fgrad(1) = x
    fgrad(2) = 1

  end subroutine lgradient

end module the_function

program what_about_rline

  use huber_curvefit
  use the_function
  use rektorys
  use iso_fortran_env
  implicit none

  real(REAL64), parameter :: range = 2
  real(REAL64), parameter :: amin = -10
  real(REAL64), parameter :: amax = 10
  real(REAL64), parameter :: bmin =-10
  real(REAL64), parameter :: bmax = 50
  real(REAL64), parameter :: smin = 0
  real(REAL64), parameter :: smax = 10

  type(line_function), target :: fcn
  real(REAL64), dimension(2) :: p, perr
  real(REAL64), dimension(2,2) :: H, cov
  real(REAL64) :: s, stdsig, gtol
  integer :: iflag,i,j,n

  p = [-2, 30]
  s = 2
  s = 1
  p = 0
  n = size(x)
  gtol = 1e3*size(p)*epsilon(gtol)

!  y = 3 + x / 2
  call huberfit(fcn,x,y,p,perr,stdsig,s,H,gtol,666, &
       .true.,'/var/tmp/rline_step',iflag)
  write(error_unit,*) 'rfit: iflag=',iflag
  write(error_unit,*) 'scale=',s,' stdsig=',stdsig

  do i = 1,2
     write(error_unit,*) i,': p=',p(i),'+-',perr(i)
  end do
!     ' s=',s,'    iflag=',iflag
!  write(*,*) 'B=',H

  cov = 0
  forall( i = 1:2, j = 1:2, H(i,j) > 0 ) cov(i,j) = H(i,j) / sqrt(H(i,i)*H(j,j))
  do i = 1,2
     write(error_unit,'(2(a,2f15.7))') 'H:',H(i,:),'   cov:',cov(i,:)
  end do

  open(1,file='/var/tmp/rline.par')
  write(1,*) 'XMIN = ',amin
  write(1,*) 'XMAX = ',amax
  write(1,*) 'YMIN = ',bmin
  write(1,*) 'YMAX = ',bmax
  close(1)

  call graph(fcn,s,amin,amax,bmin,bmax,x,y,'/var/tmp/rline_fcn')

contains

    subroutine graph(fcn,s,x0,x1,y0,y1,xdata,ydata,filename)

      use rfun

      class(line_function), intent(in) :: fcn
      real(REAL64), intent(in) :: x0,x1,y0,y1,s
      character(len=*), intent(in) :: filename
      real(REAL64), dimension(:), intent(in) :: xdata, ydata
      integer, parameter :: nx = 200
      integer, parameter :: ny = 200
      real(REAL64) :: x,dx,y,dy,f,r
      integer :: i,j,k

      dx = (x1 - x0) / (nx - 1)
      dy = (y1 - y0) / (ny - 1)

      open(1,file=filename)
      do j = 1, ny
         y = y0 + dy*j
         do i = 0, nx
            x = x0 + dx*i
            f = 0
            do k = 1, size(xdata)
               r = ydata(k) - fcn%val(xdata(k),[x,y])
               f = f + ihuber( r / s )
            end do
            write(1,*) x,y,f
         end do
         write(1,*)
      end do
      close(1)

    end subroutine graph

end program what_about_rline
