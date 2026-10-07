!
!  Copyright © 2024-2025 F.Hroch (hroch@physics.muni.cz)
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
!
!
!  fmin in module version
!
!  This module provides modified version of fmin procedure
!  by FMM library developed by Forsythe, Malcolm, and Moler
!  (http://www.netlib.org/fmm/index.html) and Brent's book
!  Algorithms for minimization without derivatives.
!
!  The changes by FH:
!   * limited counts of evaluation
!   * epsilon by the standard function
!   * Code for Fortran 2008+
!   * the optimised function is implemented as a class, suitable for multithreading
!
!   Additionally, the modified transcription of the original Brent's algorithm
!   in Algol60 has been added. The convergence limits, flags, etc. are supplemented.
!   There is an issue: if the result of multiplication of values is over
!   machine limit, the parabolic interpolation has overflow; the golden-search
!   is still possible, yet is left unimplemented.
!
!   iflag == 1 .. the minimum is found inside ax < x < bx
!   iflag == 2 .. the minimum is at the left endpoint x <= ax
!   iflag == 3 .. the minimum is at the ritg endpoint x >= bx
!   iflag == 5 .. improper input parameters
!   iflag == 6 .. no convergence

module fminim_REAL64

  use iso_fortran_env
  implicit none
  private
  public :: fmin_REAL64, goldfinger_REAL64

  type, abstract, public :: fmin_function_REAL64
   contains
     procedure(f), public, deferred :: f
  end type fmin_function_REAL64

  abstract interface
     real(REAL64) function f(this,x)
       import :: REAL64, fmin_function_REAL64
       implicit none
       class(fmin_function_REAL64) :: this
       real(REAL64), intent(in) :: x
     end function f
  end interface

  type, abstract, public :: goldfinger_function_REAL64
   contains
     procedure(gold), public, deferred :: f
  end type goldfinger_function_REAL64

  abstract interface
     real(REAL64) function gold(this,x)
       import :: REAL64, goldfinger_function_REAL64
       implicit none
       class(goldfinger_function_REAL64) :: this
       real(REAL64), intent(in) :: x
     end function gold
  end interface

  logical, parameter :: debug = .false.

contains

  function fmin_REAL64(f,ax,bx,tol) result(fmin)

    class(fmin_function_REAL64), intent(in) :: f
    real(REAL64) :: fmin
    real(REAL64), intent(in) :: ax,bx,tol

!
!      an approximation  x  to the point where  f  attains a minimum  on
!  the interval  (ax,bx)  is determined.
!
!
!  input..
!
!  ax    left endpoint of initial interval
!  bx    right endpoint of initial interval
!  f     function subprogram which evaluates  f(x)  for any  x
!        in the interval  (ax,bx)
!  tol   desired length of the interval of uncertainty of the final
!        result ( .ge. 0.0d0)
!
!
!  output..
!
!  fmin  abcissa approximating the point where  f  attains a minimum
!
!
!      the method used is a combination of  golden  section  search  and
!  successive parabolic interpolation.  convergence is never much slower
!  than  that  for  a  fibonacci search.  if  f  has a continuous second
!  derivative which is positive at the minimum (which is not  at  ax  or
!  bx),  then  convergence  is  superlinear, and usually of the order of
!  about  1.324....
!      the function  f  is never evaluated at two points closer together
!  than  eps*abs(fmin) + (tol/3), where eps is  approximately the square
!  root  of  the  relative  machine  precision.   if   f   is a unimodal
!  function and the computed values of   f   are  always  unimodal  when
!  separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
!  the abcissa of the global minimum of  f  on the interval  ax,bx  with
!  an error less than  3*eps*abs(fmin) + tol.  if   f   is not unimodal,
!  then fmin may approximate a local, but perhaps non-global, minimum to
!  the same accuracy.
!      this function subprogram is a slightly modified  version  of  the
!  algol  60 procedure  localmin  given in richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
!
    real(REAL64), parameter :: eps = sqrt(epsilon(1.0_REAL64))
    real(REAL64), parameter ::  c = (3 - sqrt(5.0_REAL64)) / 2
    real(REAL64) :: a,b,d,e,xm,p,q,r,tol1,tol2,u,v,w
    real(REAL64) :: fu,fv,fw,fx,x

!  FH: to prevent infinity loop, Fibonacci search is O(log(n))
!  see https://en.wikipedia.org/wiki/Fibonacci_search_technique
!  3x (?), precision(fmin) = 16, ... yes, we are realist....
    integer :: ieval, maxeval

    maxeval = max(int(3*log((bx - ax) / tol)), precision(1.0_REAL64))
    ieval = 0

!
!  c is the squared inverse of the golden ratio
!
!    c = (3 - sqrt(real(5,REAL64))) / 2
!
!  eps is approximately the square root of the relative machine
!  precision.
!
!  FH: the loop is replaced by an internal function (to be faster)
!      eps = 1.0d00
!   10 eps = eps/2.0d00
!      tol1 = 1.0d0 + eps
!      if (tol1 .gt. 1.0d00) go to 10
!  eps = epsilon(1.0d0)
!  eps = dsqrt(eps)
!
!  initialization
!
    a = ax
    b = bx
    v = a + c*(b - a)
    w = v
    x = v
    e = 0
    d = e
    fx = f%f(x)
    fv = fx
    fw = fx
!
!  main loop starts here
!
20  xm = (a + b) / 2
    tol1 = eps*abs(x) + tol/3
    tol2 = 2*tol1
!
! FH: number of evaluations
!
    ieval = ieval + 1
!
!  check stopping criterion
!
    if (abs(x - xm) <= (tol2 - (b - a)/2) ) go to 90
!
!  FH: maximum number of evaluations reached
    if( ieval > maxeval ) goto 90
!
! is golden-section necessary
!
    if (abs(e) <= tol1) go to 40
!
!  fit parabola
!
    r = (x - w)*(fx - fv)
    q = (x - v)*(fx - fw)
    p = (x - v)*q - (x - w)*r
    q = 2*(q - r)
    if (q > 0) p = -p
    q =  abs(q)
    r = e
    e = d
!
!  is parabola acceptable
!
    if (abs(p) >= abs(q*r/2)) go to 40
    if (p <= q*(a - x)) go to 40
    if (p >= q*(b - x)) go to 40
!
!  a parabolic interpolation step
!
    d = p/q
    u = x + d
!
!  f must not be evaluated too close to ax or bx
!
    if ((u - a) < tol2) d = sign(tol1, xm - x)
    if ((b - u) < tol2) d = sign(tol1, xm - x)
    go to 50
!
!  a golden-section step
!
40  if (x >= xm) e = a - x
    if (x < xm) e = b - x
    d = c*e
!
!  f must not be evaluated too close to x
!
50  if (abs(d) >= tol1) u = x + d
    if (abs(d) < tol1) u = x + sign(tol1, d)
    fu = f%f(u)
!
!  update  a, b, v, w, and x
!
    if (fu > fx) go to 60
    if (u >= x) a = x
    if (u < x) b = x
    v = w
    fv = fw
    w = x
    fw = fx
    x = u
    fx = fu
    go to 20
60  if (u < x) a = u
    if (u >= x) b = u
    if (fu <= fw) go to 70
    if (w == x) go to 70
    if (fu <= fv) go to 80
    if (v == x) go to 80
    if (v == w) go to 80
    go to 20
70  v = w
    fv = fw
    w = u
    fw = fu
    go to 20
80  v = u
    fv = fu
    go to 20
!
!  end of main loop
!
90  fmin = x

  end function fmin_REAL64


  subroutine goldfinger_REAL64(f,ax,bx,gtol,x,iflag)

    ! A minimalisation of a function: this is a particular implementation
    ! of locmin() procedure by Brent (1972). The minimum is searched
    ! on the interval [ax..bx].

    class(goldfinger_function_REAL64), intent(in) :: f
    real(REAL64), intent(in) :: ax, bx, gtol
    real(REAL64), intent(out) :: x
    integer, intent(out) :: iflag

    real(REAL64), parameter :: Kg = 1.44 ! Kg = 1/log2((1+sqrt(5))/2)
    real(REAL64), parameter :: log2 = 0.6931 ! = log(2)
    real(REAL64), parameter :: c = (3 - sqrt(5.0_REAL64))/2
    real(REAL64), parameter :: sqrteps = sqrt(epsilon(c))
    real(REAL64), parameter :: macheps = epsilon(c)
    real(REAL64) :: a,b,e,d,p,q,r,fw,fv,fx,fu,v,w,u,m,tol,t2
    integer :: iter, maxiter

    a = ax
    b = bx
    if( .not. (b > a .and. gtol > 0 ) ) then
       iflag = 5
       return
    end if

    iflag = 6
    maxiter = int(2*Kg*(log((b-a)/gtol)/log2)**2) + 1
    if( debug ) write(error_unit,'(a,i0,2(a,g0.5))') &
         'goldfinger 0: maxiter=',maxiter,', a=',a,', b=',b

    x = a + c*(b - a)
    w = x; v = x;
    e = 0; d = 0;
    fx = f%f(x)
    fv = fx; fw = fx;

    do iter = 1, maxiter

       tol = sqrteps*abs(x) + gtol
       t2 = 2*tol

       m = (a + b) / 2

       if( debug ) write(error_unit,'(a,i0,2(a,g0.5),a,es8.2)') &
            'goldfinger ',iter,': x=',x,', f(x)=',fx,', tol=',tol

       ! termination checks
       if( max(x - a, b - x) < t2 )  then
!          write(*,*) x-a,b-x,abs(b-bx) < epsilon(b),abs(a-ax) < epsilon(a)
          if( (bx - x) < t2 ) then
             iflag = 3
          else if( (x - ax) < t2 ) then
             iflag = 2
          else
             iflag = 1
          end if
          exit
       end if

       r = 0; q = r; p = q;
       if( abs(e) > tol ) then
          ! fit parabola
          r = (x - w)*(fx - fv)
          q = (x - v)*(fx - fw)
          p = (x - v)*q - (x - w)*r
          q = 2*(q - r)
          if( q > 0 ) then
             p = -p
          else
             q = -q
          end if
          r = e
          e = d
       end if

       if( abs(p) < abs(q*r/2) .and. p < q*(a-x) .and. p < q*(b-x) ) then
          ! A "parabolic interpolation" step
          d = p/q
          u = x + d
          if( (u - a) < t2 .or. (b - u) < t2 ) then
             if( x < m ) then
                d = tol
             else
                d = -tol
             end if
          end if
       else
          ! a "golden section" step
          if( x < m ) then
             e = b - x
          else
             e = a - x
          end if
          d = c*e
       end if

       ! f must not be evaluated too close x
       if( abs(d) > tol ) then
          u = x + d
       else
          u = x + sign(tol,d)
       end if
       fu = f%f(u)

       if( fu <= fx ) then
          if( u < x ) then
             b = x
          else
             a = x
          end if
          v = w; w = x; x = u;
          fv = fw; fw = fx; fx = fu
       else
          if( u < x ) then
             a = u
          else
             b = u
          end if
          if( fu <= fw .or. abs(w - x) < macheps ) then
             v = w; w = u;
             fv = fw; fw = fu
          else if( fu <= fv .or. abs(v - c) < macheps .or. abs(v - w) < macheps )then
             v = u
             fv = fu
          end if
       end if

    end do

  end subroutine goldfinger_REAL64


end module fminim_REAL64
