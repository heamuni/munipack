!
!  Copyright © 2023 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of OakLeaf.
!
!  OakLeaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  OakLeaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with OakLeaf.  If not, see <http://www.gnu.org/licenses/>.
!
!
!  zeroin in module version
!
!  This module provides modified version of zeroin procedure
!  by FMM library developed by Forsythe, Malcolm, and Moler
!  (http://www.netlib.org/fmm/index.html) and Brent's book
!  Algorithms for minimization without derivatives.
!
!  The changes by FH:
!   * limited counts of evaluation
!   * epsilon by the standard function
!   * Code for Fortran 2008+
!   * function pass replaced by class, suitable for multithread run

module zeroinmod

  use iso_fortran_env
  implicit none
  private
  public :: zeroin

  type, abstract, public :: zeroin_function
   contains
     procedure(f), public, deferred :: f
  end type zeroin_function

  abstract interface
     real(REAL64) function f(this,x)
       import :: REAL64, zeroin_function
       class(zeroin_function) :: this
       real(REAL64), intent(in) :: x
     end function f
  end interface

contains

  real(REAL64) function zeroin(f,ax,bx,tol,iflag)

    class(zeroin_function), intent(in) :: f
    real(REAL64), intent(in) :: ax,bx,tol
    integer, intent(out) :: iflag

!
!   a zero of the function  f(x)  is computed in the interval ax,bx.
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result ( .ge. 0.0d0)
!
!
!  output..
!
!  zeroin abcissa approximating a zero of  f  in the interval ax,bx
!
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  without  a  check.  zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
!  is the relative machine precision.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
!  Return flags:
!    iflag = 1, success
!    iflag = 2, solution is found at the interval edges ax,bx
!    iflag = 6, no convergence
!    iflag = 8, any improper setup ax > bx, or f(ax)*f(bx) > 0

    real(REAL64), parameter :: eps = epsilon(zeroin)
    real(REAL64) ::  a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s
    integer :: ieval, maxeval
    logical :: intx

    zeroin = 0
    if( .not. (bx > ax .and. tol > 0) ) then
       iflag = 8
       return
    end if

    ! maxeval ~= log2((bx - ax) / tol))
    maxeval = max(int(3*log((bx - ax) / tol)+1), 3*precision(zeroin))
    ieval = 0
    iflag = 6

    ! initialization
    a = ax
    b = bx
    fa = f%f(a)
    fb = f%f(b)
    intx = .true.

!    write(*,*) maxeval,fa,fb,a,b

    if( .not. (fa*fb < 0) ) then
       iflag = 8
       return
    end if

    do ieval = 1, maxeval

       if ( ieval == maxeval ) then
          iflag = 5
          zeroin = (a + b) / 2
          exit
       end if

       !
       ! begin step
       !
       if( intx ) then
          c = a
          fc = fa
          d = b - a
          e = d
       end if
       if (abs(fc) < abs(fb) ) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
       end if
       !
       ! convergence test
       !
       tol1 = 2*eps*abs(b) + 0.5*tol
       xm = 0.5_REAL64*(c - b)

!       write(*,*) '.',xm,tol1,a,b,fb,abs(xm) < tol1 .and. abs(fb) > epsilon(fb)
       if( abs(xm) < tol1 .and. abs(fb) > epsilon(fb) ) then
          iflag = 1
          zeroin = b
          exit
       end if

       !
       ! is bisection necessary
       !
       if (abs(e) < tol1 .or. abs(fa) <= abs(fb) ) then
          d = xm
          e = d
       else

          s = fb/fa
          !
          ! is quadratic interpolation possible
          !
          if ( abs(a-c) < epsilon(a) ) then
             !
             ! linear interpolation
             !
             p = 2*xm*s
             q = 1 - s
          else
             !
             ! inverse quadratic interpolation
             !
             q = fa/fc
             r = fb/fc
             p = s*(2*xm*q*(q - r) - (b - a)*(r - 1))
             q = (q - 1)*(r - 1)*(s - 1)
          end if
          !
          ! adjust signs
          !
          if (p > 0.0d0) then
             q = -q
          else
             p = -p
          end if
          s = e
          e = d
          !
          ! is interpolation acceptable
          !
          if ((2*p < 3*xm*q - abs(tol1*q)) .and. p < abs(0.5*s*q) ) then
             d = p / q
          else
             !
             ! bisection
             !
             d = xm
             e = d
          end if
       end if
       !
       ! complete step
       !
       a = b
       fa = fb
       if (abs(d) > tol1) then
          b = b + d
       else
          b = b + sign(tol1, xm)
       end if
       fb = f%f(b)
       intx = fb*(fc/abs(fc)) > real(0,REAL64)
    end do

!    write(*,*) a,b,fa,fb,tol,'.',iflag

    if( (abs(b - bx) < tol .or. abs(a - ax) < tol) .and. iflag == 1) iflag = 2

  end function zeroin

end module zeroinmod
