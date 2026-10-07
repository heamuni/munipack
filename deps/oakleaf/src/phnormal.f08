!
!  Robust estimation of mean photon rates
!
!  Copyright © 2020 F.Hroch (hroch@physics.muni.cz)
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
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!

module photonrate_normal

  use photonrate_common
  use iso_fortran_env

  implicit none

contains

  elemental pure function NormCDF(x)

    real(REAL64) :: NormCDF
    real(REAL64), intent(in) :: x
    real(REAL64), parameter :: s2 = sqrt(real(2.0,REAL64))

    NormCDF = (1 + erf(x/s2)) / 2

  end function NormCDF

  pure function mcdf(kmean,freq,sx,sy,tau,rho,phi,u) result(v)

    ! mean CDF for a set of kmeans

    real(REAL64), dimension(:), intent(in) :: u, kmean
    integer, dimension(:), intent(in) :: freq
    real(REAL64), intent(in) :: sx,sy,tau,rho,phi
    real(REAL64), dimension(:), allocatable :: v

    real(REAL64), dimension(:), allocatable :: sig
    real(REAL64) :: cosfi2
    integer :: k,l,n,m

    cosfi2 = cos(phi)**2

    k = size(kmean)
    n = 0
    v = [ (0, l = 1, size(u)) ]
    do l = 1, k
       m = freq(l)
       n = n + m
       sig = sqrt(kmean(l)*(tau/sy+u**2/sx)+rho**2*(1+((u/tau)**2-1)*cosfi2))
       v = v + m*NormCDF((kmean(l)/sig)*(u - tau))
    end do
    v = v / n

  end function mcdf

  real(REAL64) function invNormal(q)

    ! This function computes inverse of Normal function for q > 0.8.

    real(REAL64), intent(in) :: q

    real(REAL64), parameter :: pi2 = sqrt(2*PI)
    real(REAL64), parameter :: sqrt2 = sqrt(real(2.0,REAL64))
    real(REAL64), parameter :: tol = 1e-5
    real(REAL64) :: f, df, d, x, x0, q1, lnq, lnp
    integer :: iter
    logical :: converge

    if( (q < 0.5 .and. q > 0.2) .or. (q > 0.5 .and. q < 0.8 ) ) &
         stop 'invNormal is based on asymptotic expansion!'

    ! Initial estimate
    ! https://en.wikipedia.org/wiki/Normal_distribution#Quantile_function
    if( q > 0.5 ) then
       lnp = -2*log(1-q)
    else
       lnp = -2*log(q)
    end if
    x = -sqrt(lnp - log(lnp) - log(pi2))
    if( q > 0.5 ) x = -x
    x0 = x

    ! The asymptotic expansion 1-Phi(x) ~= exp(-x**2/2)/(sqrt(2*PI)*x)
    ! for 1 < x < 5 is used.
    if( q < 0.5 ) then
       lnq = log(q)
       q1 = 1 - q
    else
       lnq = log(1 - q)
       q1 = q
    end if

    converge = .false.
    do iter = 1, precision(x)
       f = -x**2/2 - log(pi2*x) - lnq
       df = -x - 1/x

       if( .not. (abs(df) > epsilon(df)) ) then
          if( debug ) write(error_unit,'(a)') &
               "Warning: fail in invNorm(): df == 0!"
          exit
       end if

       d = f /df
       x = x - d

       if( debug_iter ) write(error_unit,'(a,i2,2x,g0,2es12.2,es14.4)') &
            "z, incr., f, f': ",iter,real(x),d,f,df

       converge = abs(d) < tol
       if( converge ) exit

    end do

    ! Accurate computations, perhaps useless, if small precision is enought.
    if( converge ) then
       do iter = 1, precision(x)
          f = (1 + erf(x/sqrt2))/2 - q1
          df = exp(-x**2/2)/pi2

          if( .not. (abs(df) > epsilon(df)) ) exit

          d = f /df
          x = x - d

          if( debug_iter ) write(error_unit,'(a,i2,2x,g0,2es12.2,es14.4)') &
               "z, incr., f, f': ",iter,real(x),d,f,df

          converge = abs(d) < tol
          if( converge ) exit

       end do

    end if

    if( .not. converge ) then
       if( debug ) write(error_unit,'(a)') &
            "Warning: A divergence of invNormal occured."
       x = x0
    end if

    if( q < 0.5 ) then
       x = -x
    else
       x = x
    end if

    invNormal = x

  end function invNormal



end module photonrate_normal
