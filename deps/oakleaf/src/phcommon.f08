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


module photonrate_common

  use iso_fortran_env

  implicit none

  real(REAL64), parameter :: PI = real(3.14159265358979323844,REAL64)
  real(REAL64), parameter :: RAD = 180 / PI

  ! The first and third quartiles: the deviation of the independent
  ! variable from N(x|0,1) having : 3/4 = N(Qdev|0,1), -Qdev for 1/4.
  real, parameter :: Qdev = 0.67448975

  ! debuging
  logical :: debug = .false.
  logical, parameter :: debug_iter = .true.

end module photonrate_common
