!
!  Generate samples having noise of various statistical distributions
!
!  Copyright © 2016-2025 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of OakLeaf.
!
!  OakLeaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  OakLeaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with OakLeaf.  If not, see <http://www.gnu.org/licenses/>.
!

module noise_REAL64

  use iso_fortran_env
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter, private :: kind = REAL64

contains

  subroutine randomnoise

    integer :: n,istat
    integer, allocatable :: seed(:)

    ! code by https://gcc.gnu.org/onlinedocs/gcc-6.3.0/gfortran/
    call random_seed(size=n)
    allocate(seed(n))
    open(1, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(1) seed
       close(1)
       call random_seed(put=seed)
    end if

  end subroutine randomnoise

  ! Normal distribution
  ! https://en.wikipedia.org/wiki/Normal_distribution

  real(REAL64) function gnoise(t,s)

    ! the polar method by Box-Muller, based on Brand: Data analysis, Chapter 4.

    real(REAL64), intent(in) :: t,s
    real(REAL64) :: u1,u2,v1,v2,r2,q

    do
       call random_number(u1)
       call random_number(u2)
       v1 = 2*u1 - 1
       v2 = 2*u2 - 1
       r2 = v1**2 + v2**2
       if( r2 < 1 ) exit
    end do
    q = sqrt(-(2/r2)*log(r2))
    gnoise = v1*q*s + t

  end function gnoise

  subroutine gnoise_fill(t,s,x)

    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: x
    real(REAL64) :: u1,u2,v1,v2,r2,q
    integer :: i,n

    i = lbound(x,1)
    n = ubound(x,1)
    do
       do
          call random_number(u1)
          call random_number(u2)
          v1 = 2*u1 - 1
          v2 = 2*u2 - 1
          r2 = v1**2 + v2**2
          if( r2 < 1 ) exit
       end do
       q = sqrt(-(2/r2)*log(r2))
       x(i) = v1*q*s + t
       i = min(i+1,n)
       x(i) = v2*q*s + t
       if( i == n ) exit
       i = i + 1
    end do

  end subroutine gnoise_fill

  function gaussnoise(mean,sig) result(gnoise)

    ! Gets noise from N(mean,sig), eg. Normal distribution.

    real(kind), parameter :: sqrt2 = 1.4142135623730951_kind

    real(kind) :: gnoise,x
    real(kind), intent(in) :: mean, sig

    call random_number(x)
    gnoise = mean + sqrt2*sig*inverf(2*x - 1)

  end function gaussnoise

  subroutine gaussnoise_fill(t,sig,x)

    real(kind), intent(in) :: t,sig
    real(kind), dimension(:), intent(out) :: x
    integer :: i

    do i = 1, size(x)
       x(i) = gnoise(t,sig)
    end do

  end subroutine gaussnoise_fill

  ! Laplace distribution
  ! https://en.wikipedia.org/wiki/Laplace_distribution

  function lnoise(median,mad)

    real(kind) :: lnoise, u
    real(kind), intent(in) :: median, mad

    call random_number(u)
    u = u - 0.5_kind
    lnoise = median - mad*sign(1.0_kind,u)*log(1.0_kind - 2.0_kind*abs(u))

  end function lnoise

  ! Poisson distribution
  ! https://en.wikipedia.org/wiki/Poisson_distribution

  function pnoise(lambda)

    ! Gets noise from Poisson(lambda) distribution.

    real(kind) :: pnoise
    real(kind), intent(in) :: lambda

    ! The approximation for low lambda uses Knuth's algorithm
    ! which is valid only for exp(-lambda) > epsilon(1.0_kind),
    ! eg. lambda < 740. Normal distribution is used for larger
    ! values. Approximation of N(lambda,sqrt(lambda)) will also
    ! valid with appropriate precision.

    if( lambda  < 500 ) then
       pnoise = pnoise_knuth(lambda)
    else
       pnoise = gnoise(lambda,sqrt(lambda))
    end if

  end function pnoise

  subroutine pnoise_fill(lambda,x)

    real(kind), intent(in) :: lambda
    real(kind), dimension(:), intent(out) :: x
    integer :: i

    do i = 1, size(x)
       x(i) = pnoise(lambda)
    end do

  end subroutine pnoise_fill

  function pnoise_knuth(lam) result(k)

    ! Knuth's algorithm by http://en.wikipedia.org/wiki/Poisson_distribution
    ! This function is valid only for L=exp(-lam) > epsilon(L).

    real(kind), intent(in) :: lam
    real(kind) :: L,p,u
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

  elemental pure function inverf(p)

    real(kind), intent(in) :: p
    real(kind) :: inverf, pp,x,t,err
    integer :: j

    if( abs(p) < 1 ) then

       pp = 1 - abs(p)
       t = sqrt(-2.0_kind*log(pp/2))
       x = - 0.70711*((2.30753 + t*0.27061)/(1 + t*(0.99229 + t*0.04481)) - t)
       do j = 1,2
          err = 1 - erf(x) - pp
          x = x + err/(1.128379167099551257*exp(-x**2) - x*err)
       end do
       inverf = -sign(x,-p)

    else
       if( p > 0 ) then
          inverf = IEEE_VALUE(1.0_kind, IEEE_POSITIVE_INF)
       else
          inverf = IEEE_VALUE(1.0_kind, IEEE_NEGATIVE_INF)
       end if
    end if

  end function inverf


!!$  function ierf(x)
!!$
!!$    ! http://en.wikipedia.org/wiki/Error_function
!!$    ! Unfortunatelly, this function has large numerical
!!$    ! error (about 0.2!).
!!$
!!$    real(kind), parameter :: pi = 3.1415926535897931
!!$
!!$    real(kind), intent(in) :: x
!!$    real(kind) :: ierf, u, v, a, y, w
!!$
!!$    a = 0.140012288686666
!!$    u = log(1 - x**2)
!!$    v = 2/(pi*a)
!!$    y = v + u/2
!!$
!!$    w = sqrt(y**2 - u/a)
!!$    ierf = sign(1.0d0,x)*sqrt(w - y)
!!$
!!$  end function ierf


  ! Uniform distribution
  function unoise(a,b)

    ! Noise from uniform distribution U(a,b)

    real(kind) :: unoise, x
    real(kind), intent(in) :: a,b

    call random_number(x)
    unoise = 1.0/(b-a) * x

  end function unoise


end module noise_REAL64
