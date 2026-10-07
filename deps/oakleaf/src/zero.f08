!
!  zero
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


module zeroinmod

  use iso_fortran_env
  implicit none

  type, public :: zeroin_function
     logical :: break = .false.
   contains
     procedure :: f => fun
  end type zeroin_function

  abstract interface
     real(REAL64) function fun(this,x)
       import :: REAL64, zeroin_function
       class(zeroin_function) :: this
       real(REAL64), intent(in) :: x
     end function fun
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
  !
    real(REAL64), parameter :: eps = epsilon(zeroin)
    real(REAL64), parameter ::  a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s
!    double precision  dabs,dsign
    integer :: ieval, maxeval
    logical :: int
!
!  compute eps, the relative machine precision
!
!  FH: the loop is replaced by an internal function (to be faster)
!      eps = 1.0d0
!   10 eps = eps/2.0d0
!      tol1 = 1.0d0 + eps
!      if (tol1 .gt. 1.0d0) go to 10
  !



    ! maxeval ~= log2((bx - ax) / tol))
    maxeval = max(int(3*log((bx - ax) / tol)+1), precision(zeroin))
    ieval = 0
    iflag = 666

    ! initialization
    !
    a = ax
    b = bx
    fa = f%f(a)
    fb = f%f(b)
    int = .true.

    if( .not. (fa*fb < 0) ) then
       iflag = 2
       return
    end if

    do ieval = 1, maxeval

       if ( ieval == maxeval ) then
          flag = 1
          exit !!!!!!!!!!!!!!!
       end if


       !
       ! begin step
       !
!20     continue
       if( int ) then
          c = a
          fc = fa
          d = b - a
          e = d
       end if
!30  continue
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

       if( abs(xm) < tol1 .and. abs(fb) > epsilon(fb) ) then
          flag = 0
          exit
       end if

!    if (abs(xm) < tol1 .and. abs(fb) > epsilon(fb) ) then

    !
    ! is bisection necessary
    !
       if (abs(e) < tol1 .or. abs(fa) <= abs(fb) ) then
          d = xm
          e = d
       else

          s = fb/fa

!  if (dabs(fa) .le. dabs(fb)) go to 70
!
! is quadratic interpolation possible
!
          if ( abs(a-c) > epsilon(a) ) then
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
             s = fb/fa
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
          !  p = dabs(p)
          s = e
          e = d
          !
! is interpolation acceptable
!
          if ((2*p < 3*xm*q - abs(tol1*q)) .and. p < abs(0.5*s*q) ) then
            !  if (p .ge. dabs(0.5d0*e*q)) go to 70
             d = p / q
          else
             d = xm
             e = xm
             !            e = d
!            d = p/q
          end if
       end if
!  go to 80
!
! bisection
!
!70 d = xm
!  e = d
!
! complete step
!
80     a = b
       fa = fb
       if (abs(d) > tol1) then
          b = b + d
       else
          b = b + sign(tol1, xm)
       end if
       fb = f(b)
       ieval = ieval + 1
       int = fb*(fc/abs(fc)) > 0.0_REAL64
!       if ((fb*(fc/abs(fc))) > 0.0d0) then
!          go to 20
!       else
!          go to 30
!       end if
!
! done
   !
    end do
90  zeroin = b

    if( abs(b - bx) < tol .or. abs(b - ax) < tol ) iflag = 3

!   converge = ieval < maxeval .and.
!  return

  end function zeroin

end module zeroinmod
