!
!  Robust minimization functions
!
!  Copyright © 2011 - 2025 F.Hroch (hroch@physics.muni.cz)
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

module rfun_REAL64

  use iso_fortran_env
  implicit none
  private

  ! M-estimator functions
  ! Hogg in Launer, Wilkinson: Robustness in Statistics

  ! parameters of robust functions
  real(REAL64), parameter :: a = 1.345
  real(REAL64), parameter :: a2 = a**2 / 2
  real(REAL64), parameter :: c = 6
  real(REAL64), parameter :: c2 = c**2 / 2
  real(REAL64), parameter :: delta = 2
!  real(REAL64), parameter,public :: imaxtukey = c

  ! parameters for effective scale (see test/eparam.f08)
!  real(REAL64),parameter,public :: shuber = 0.8705
!  real(REAL64),parameter,public :: spshuber = 0.8652
!  real(REAL64),parameter,public :: stukey = 0.9008

  ! parameters for entropy determination (see test/eparam.f08)
!  real(REAL64),parameter,public :: ehuber = 2.0045
!  real(REAL64),parameter,public :: etukey = 2.1227

  public :: huber_REAL64, dhuber_REAL64, ihuber_REAL64, &
       pshuber_REAL64, dpshuber_REAL64, ipshuber_REAL64, &
       tukey_REAL64, dtukey_REAL64, itukey_REAL64, tukeys_REAL64, &
       square_REAL64, dsquare_REAL64, isquare_REAL64

contains

  !--------------------------------------------------------------
  !
  ! Robust functions
  !

  ! Huber

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


  ! Pseudo-Huber (continous version)

  elemental pure function pshuber_REAL64(x) result(pshuber)

    real(REAL64) :: pshuber
    real(REAL64), intent(in) :: x

    pshuber = x/sqrt(1 + (x/delta)**2)

  end function pshuber_REAL64


  elemental pure function dpshuber_REAL64(x) result(dpshuber)

    real(REAL64) :: dpshuber
    real(REAL64), intent(in) :: x

    dpshuber = 1/sqrt((1 + (x/delta)**2)**3)

  end function dpshuber_REAL64

  elemental pure function ipshuber_REAL64(x) result(ipshuber)

    real(REAL64) :: ipshuber
    real(REAL64), intent(in) :: x

    ipshuber = delta**2*(sqrt(1 + (x/delta)**2) - 1)

  end function ipshuber_REAL64

  ! Tukey

  elemental pure function tukey_REAL64(x) result(tukey)

    real(REAL64) :: tukey
    real(REAL64), intent(in) :: x

    if( abs(x) < c )then
       tukey = x*(1 - (x/c)**2)**2
    else
       tukey = 0
    endif

  end function tukey_REAL64

  elemental pure function dtukey_REAL64(x) result(dtukey)

    real(REAL64) :: dTukey
    real(REAL64), intent(in) :: x
    real(REAL64) :: t

    if( abs(x) < c )then
       t = (x/c)**2
       dtukey = 1 + t*(5*t - 6)
    else
       dtukey = 0
    endif

  end function dTukey_REAL64

  elemental pure function iTukey_REAL64(x) result(itukey)

    real(REAL64) :: itukey
    real(REAL64), intent(in) :: x
    real(REAL64) :: t

    if( abs(x) < c )then
       t = (x/c)**2
       iTukey = c2*t*(1 + t*(t/3 - 1))
    else
       iTukey = c
    endif

  end function iTukey_REAL64

  elemental pure subroutine tukeys_REAL64(x,tukey,dtukey,itukey)

    real(REAL64), intent(in) :: x
    real(REAL64), intent(out) :: tukey,dtukey
    real(REAL64), intent(out), optional :: itukey
    real(REAL64) :: t
    logical :: l

    l = present(itukey)

    if( abs(x) < c )then
       t = (x/c)**2
       tukey = x*(1 - t)**2
       dtukey = 1 + t*(5*t - 6)
       if(l) iTukey = c2*t*(1 + t*(t/3 - 1))
    else
       tukey = 0
       dtukey = 0
       if(l) iTukey = c
    endif

  end subroutine tukeys_REAL64

  !--------------------------------------------------------------
  !
  ! Non-Robust functions
  !

  elemental pure function square_REAL64(x) result(square)

    real(REAL64) :: square
    real(REAL64), intent(in) :: x

    square = x

  end function square_REAL64

  elemental pure function dsquare_REAL64(x) result(dsquare)

    real(REAL64) :: dsquare
    real(REAL64), intent(in) :: x

    dsquare = x ! to prevent compiler's warnings
    dsquare = 1

  end function dsquare_REAL64

  elemental pure function isquare_REAL64(x) result(isquare)

    real(REAL64) :: isquare
    real(REAL64), intent(in) :: x

    isquare = x**2 / 2

  end function isquare_REAL64

end module rfun_REAL64
