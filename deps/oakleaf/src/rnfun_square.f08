!
!  Square function in fashion of robust functions
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

module rnfun_square_REAL64

  use iso_fortran_env
  implicit none
  private

  public :: square_REAL64, dsquare_REAL64, isquare_REAL64

contains

    elemental pure function square_REAL64(x) result(square)

    real(REAL64) :: square
    real(REAL64), intent(in) :: x

    square = x

  end function square_REAL64

  elemental pure function dsquare_REAL64(x) result(dsquare)

    real(REAL64) :: dsquare
    real(REAL64), intent(in) :: x

    dsquare = 1

  end function dsquare_REAL64

  elemental pure function isquare_REAL64(x) result(isquare)

    real(REAL64) :: isquare
    real(REAL64), intent(in) :: x

    isquare = x**2 / 2

  end function isquare_REAL64

end module rnfun_square_REAL64
