!
!  Huber's robust functions
!
!  Copyright © 2025-2026 F.Hroch (hroch@physics.muni.cz)
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

module rfun_huber_REAL64

  use iso_fortran_env
  implicit none
  private

  real(REAL64), parameter :: a = 1.345_REAL64
  real(REAL64), parameter :: a2 = a**2 / 2

  real(REAL64), parameter :: R = 0.05_REAL64
  real(REAL64), parameter :: u1 = a - R
  real(REAL64), parameter :: u2 = a + R
  real(REAL64), parameter :: za = 2*a - R/2
  real(REAL64), parameter :: r2 = 1 / (2 * R)
  real(REAL64), parameter :: da = (a-R)**2 + za*R - 2*R**2/3
  real(REAL64), parameter :: r6 = 1 / (6 * R)

  public :: huber_REAL64, dhuber_REAL64, ihuber_REAL64, &
       pshuber_REAL64, dpshuber_REAL64, ipshuber_REAL64, hubersqpf_REAL64

contains

  elemental pure function huber_REAL64(x) result(huber)

    real(REAL64) :: huber
    real(REAL64), intent(in) :: x

    huber = max(-a,min(a,x))
!    if( abs(x) < a )then
!       huber = x
!    else
!       huber = sign(a,x)
!    endif

  end function huber_REAL64


  elemental pure function dhuber_REAL64(x) result(dhuber)

    real(REAL64) :: dhuber
    real(REAL64), intent(in) :: x

    if( abs(x) < a )then
       dhuber = 1
    else
       dhuber = 0
    endif

  end function dhuber_REAL64

  elemental pure function ihuber_REAL64(x) result(ihuber)

    real(REAL64) :: ihuber
    real(REAL64), intent(in) :: x

    if( abs(x) < a ) then
       ihuber = x**2 / 2
    else
       ihuber = a*abs(x) - a2
    end if

  end function ihuber_REAL64

  elemental pure function pshuber_REAL64(x) result(huber)

    real(REAL64) :: huber
    real(REAL64), intent(in) :: x
    real(REAL64) :: u

    u = abs(x)
    if( u <= u1 ) then
       huber = x
    else if( u >= u2 ) then
       huber = sign(a,x)
    else !if( u1 < |x| < u2 ) then
       u = u - a
       huber = 0.5_REAL64 * sign(za + u*(1 - r2*u),x)
    end if

  end function pshuber_REAL64


  elemental pure function dpshuber_REAL64(x) result(dhuber)

    real(REAL64) :: dhuber
    real(REAL64), intent(in) :: x
    real(REAL64), parameter :: R1 = 1 / R
    real(REAL64) :: u

    u = abs(x)
    if( u <= u1 ) then
       dhuber = 1
    else if( u >= u2 ) then
       dhuber = 0
    else ! if( u1 < |x| < u2 ) then
       dhuber = 0.5_REAL64 * (1 - (u - a)*R1)
    endif

  end function dpshuber_REAL64

  elemental pure function ipshuber_REAL64(x) result(ihuber)

    real(REAL64) :: ihuber
    real(REAL64), intent(in) :: x

    real(REAL64) :: u

    u = abs(x)
    if( u < u1 ) then
       ihuber = x**2 / 2
    else if( u > u2 ) then
       ihuber = a*u - a2
    else ! if( u1 < |x| < u2 ) then
       u = u - a
       ihuber = 0.5_REAL64 * (da + u*(za + u*(0.5_REAL64 - r6*u)))
    end if

  end function ipshuber_REAL64


  elemental pure subroutine hubersqpf_REAL64(x,rho,psi)

    real(REAL64), intent(in) :: x
    real(REAL64), intent(out) :: rho,psi
    real(REAL64) :: u

    u = abs(x)
    if( u < u1 ) then
       psi = x
       rho = x**2 / 2
    else if( u > u2 ) then
       psi = sign(a,x)
       rho = a*u - a2
    else
       u = u - a
       psi = 0.5_REAL64 * sign(za + u*(1 - r2*u),x)
       rho = 0.5_REAL64 * (da + u*(za + u*(0.5_REAL64 - r6*u)))
    end if

  end subroutine hubersqpf_REAL64



  ! Pseudo-Huber (continous version by wiki, it is abandoned due numerical incompability)

!!$  elemental pure function pshuber_REAL64(x) result(pshuber)
!!$
!!$    real(REAL64) :: pshuber
!!$    real(REAL64), intent(in) :: x
!!$
!!$    pshuber = x/sqrt(1 + (x/delta)**2)
!!$
!!$  end function pshuber_REAL64
!!$
!!$
!!$  elemental pure function dpshuber_REAL64(x) result(dpshuber)
!!$
!!$    real(REAL64) :: dpshuber
!!$    real(REAL64), intent(in) :: x
!!$
!!$    dpshuber = 1/sqrt((1 + (x/delta)**2)**3)
!!$
!!$  end function dpshuber_REAL64
!!$
!!$  elemental pure function ipshuber_REAL64(x) result(ipshuber)
!!$
!!$    real(REAL64) :: ipshuber
!!$    real(REAL64), intent(in) :: x
!!$
!!$    ipshuber = delta**2*(sqrt(1 + (x/delta)**2) - 1)
!!$
!!$  end function ipshuber_REAL64

end module rfun_huber_REAL64
