!
!  Generate noises with various statistical distributions
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


module noise

  use iso_fortran_env

  implicit none

contains

  function gnoise(mean,sig)

    ! Gets noise from N(mean,sig), eg. Normal distribution.

    real(REAL64), parameter :: sqrt2 = sqrt(2.0_REAL64)

    real(REAL64) :: gnoise,x
    real(REAL64), intent(in) :: mean, sig

    call random_number(x)
    gnoise = mean + sqrt2*sig*inverf(2*x-1)

  end function gnoise

  function gnoise2(t, s) result(coords)

    ! the polar method by Box-Muller, based on Brand: Data analysis, Chapter 4.

    real(REAL64), intent(in) :: t,s
    real(REAL64) :: u1,u2,v1,v2,r2,q
    real(REAL64), dimension(2) :: coords

    do
       call random_number(u1)
       call random_number(u2)
       v1 = 2*u1 - 1
       v2 = 2*u2 - 1
       r2 = v1**2 + v2**2
       if( r2 < 1 ) exit
    end do
    q = sqrt(-(2/r2)*log(r2))
    coords = q * [v1, v2] * s + t

  end function gnoise2

  function pnoise(lambda)

    ! Gets noise from Poisson(lambda) distribution.

    real(REAL64) :: pnoise
    real(REAL64), intent(in) :: lambda

    ! The approximation for low lambda uses Knuth's algorithm
    ! which is valid only for exp(-lambda) > epsilon(1.0_REAL64),
    ! eg. lambda < 740. Normal distribution is used for larger
    ! values. Approximation of N(lambda,sqrt(lambda)) will also
    ! valid with appropriate precision.

    if( lambda  < 500 ) then
       pnoise = pnoise_knuth(lambda)
    else
       pnoise = gnoise(lambda,sqrt(lambda))
    end if

  end function pnoise

  function pnoise_knuth(lam) result(k)

    ! Knuth's algorithm by http://en.wikipedia.org/wiki/Poisson_distribution
    ! This function is valid only for L=exp(-lam) > epsilon(L).

    real(REAL64), intent(in) :: lam
    real(REAL64) :: L,p,u
    integer :: k

    L = exp(-lam)
    k = 0
    p = 1
    do
       k = k + 1
       call random_number(u)
       p = p*u
       if( .not. (p > L) ) exit
    enddo
    k = k - 1

  end function pnoise_knuth


  ! The function inverf is adapted Numerical Recipes subroutine
  ! by chapter 6. Error functions. Its precision is better
  ! then 1e-6 on  -5 < x < 5.
  !
  ! Unfortunately, the erf(.) should be supported by compiler.
  ! It requires at least Fortran 2008.

    pure function inverf(p)

    real(REAL64), intent(in) :: p
    real(REAL64) :: inverf, pp,x,t,err
    integer :: j

    if( abs(p) < 1 ) then

       pp = 1 - abs(p)
       t = sqrt(-2.0_REAL64*log(pp/2))
       x = - 0.70711*((2.30753 + t*0.27061)/(1 + t*(0.99229 + t*0.04481)) - t)
       do j = 1,2
          err = 1 - erf(x) - pp
          x = x + err/(1.128379167099551257*exp(-x**2) - x*err)
       end do
       inverf = -sign(x,-p)

    else
       inverf = sign(huge(p),p)
    end if

  end function inverf


  pure function ierf(x)

    ! http://en.wikipedia.org/wiki/Error_function
    ! Unfortunatelly, this function has large numerical
    ! error (about 0.2!).

    real(REAL64), parameter :: pi = 3.1415926535897931_REAL64

    real(REAL64), intent(in) :: x
    real(REAL64) :: ierf, u, v, a, y, w

    a = 0.140012288686666_REAL64
    u = log(1 - x**2)
    v = 2/(pi*a)
    y = v + u/2

    w = sqrt(y**2 - u/a)
    ierf = sign(1.0_REAL64,x)*sqrt(w - y)

  end function ierf


end module noise
