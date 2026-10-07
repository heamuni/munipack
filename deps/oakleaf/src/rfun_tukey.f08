!
!  Tukey's robust functions
!
!  Copyright © 2025 F.Hroch (hroch@physics.muni.cz)
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

module rfun_tukey_REAL64

  use iso_fortran_env
  implicit none
  private

  real(REAL64), parameter :: c = 6
  real(REAL64), parameter :: c2 = c**2 / 2

  public :: tukey_REAL64, dtukey_REAL64, itukey_REAL64, tukeys_REAL64, &
       tukeysqpf_REAL64

contains

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

  elemental pure subroutine tukeysqpf_REAL64(x,rho,psi)

    real(REAL64), intent(in) :: x
    real(REAL64), intent(out) :: rho,psi
    real(REAL64) :: t

    if( abs(x) < c )then
       t = (x/c)**2
       psi = x*(1 - t)**2
       rho = c2*t*(1 + t*(t/3 - 1))
    else
       psi = 0
       rho = c
    endif

  end subroutine tukeysqpf_REAL64


end module rfun_tukey_REAL64
