!
!  The definition of the straight line of Rektorys:
!  Rektorys,K.(1963): Karel Rektorys a spolupracovníci: přehled užité matematiky.
!
!  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
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

module rektorys

  use iso_fortran_env
  implicit none

  real(REAL64), dimension(14), target :: &
       x = [ 3,4,6,6,6,7,8,9,11,11,12,12,14,16 ]
  real(REAL64), dimension(14), target :: &
       y = [ 24.82, 23.26, 14.77, 19.06, 14.79, 17.66, 11.83, &
       14.82, 5.16, 11.12, 8.04, 2.72, 0.74, 1.21 ]

end module rektorys
