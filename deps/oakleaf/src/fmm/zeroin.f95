!
!  zeroin
!
!  This procedure is adopted version of zeroin procedure
!  by FMM library developed by Forsythe, Malcolm, and Moler
!  (http://www.netlib.org/fmm/index.html) for Fortran 95.
!  Some changes to code are introduced by FH:
!   * limited counts of evaluation
!   * epsilon by standard function
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


double precision function zeroin(ax,bx,f,tol)

  implicit none

  double precision, intent(in) :: ax,bx,tol

  interface
     function f(x)
       double precision :: f
       double precision, intent(in) :: x
     end function f
  end interface

!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
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
  double precision :: eps = epsilon(zeroin)
  double precision  a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s
  double precision  dabs,dsign
  integer :: ieval, maxeval
!
!  compute eps, the relative machine precision
!
!  FH: the loop is replaced by an internal function (to be faster)
!      eps = 1.0d0
!   10 eps = eps/2.0d0
!      tol1 = 1.0d0 + eps
!      if (tol1 .gt. 1.0d0) go to 10
  !

  maxeval = max(int(3*log((bx - ax) / tol)), precision(zeroin))
  ieval = 0

! initialization
!
  a = ax
  b = bx
  fa = f(a)
  fb = f(b)
!
! begin step
!
20 c = a
  fc = fa
  d = b - a
  e = d
30 if (dabs(fc) .ge. dabs(fb)) go to 40
  a = b
  b = c
  c = a
  fa = fb
  fb = fc
  fc = fa
!
! convergence test
!
40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
  xm = .5*(c - b)
  if (dabs(xm) .le. tol1) go to 90
  if (fb .eq. 0.0d0) go to 90
  if (ieval > maxeval ) goto 90
!
! is bisection necessary
!
  if (dabs(e) .lt. tol1) go to 70
  if (dabs(fa) .le. dabs(fb)) go to 70
!
! is quadratic interpolation possible
!
  if (a .ne. c) go to 50
!
! linear interpolation
!
  s = fb/fa
  p = 2.0d0*xm*s
  q = 1.0d0 - s
  go to 60
!
! inverse quadratic interpolation
!
50 q = fa/fc
  r = fb/fc
  s = fb/fa
  p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
  q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
!
! adjust signs
!
60 if (p .gt. 0.0d0) q = -q
  p = dabs(p)
!
! is interpolation acceptable
!
  if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70
  if (p .ge. dabs(0.5d0*e*q)) go to 70
  e = d
  d = p/q
  go to 80
!
! bisection
!
70 d = xm
  e = d
!
! complete step
!
80 a = b
  fa = fb
  if (dabs(d) .gt. tol1) b = b + d
  if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
  fb = f(b)
  ieval = ieval + 1
  if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
  go to 30
!
! done
!
90 zeroin = b
  return

end function zeroin
