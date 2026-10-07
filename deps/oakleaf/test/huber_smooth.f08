!
!  Test for Huber, including smooth, function
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

! gfortran -Wall ../src/rfun_huber.f08 ../src/rfun.f08  huber_smooth.f08

program huber_smooth

  use rfun
  use iso_fortran_env
  implicit none

  real(REAL64),parameter :: dlt = sqrt(epsilon(1.0_REAL64))
  integer :: i
  real(REAL64) :: u


  do i = -1000,1000
     u = 0.01_REAL64*i
     write(*,*) u,huber(u),pshuber(u),dhuber(u),dpshuber(u),&
          ihuber(u),ipshuber(u),(ipshuber(u+dlt) - ipshuber(u-dlt))/(2*dlt),&
          (pshuber(u+dlt) - pshuber(u-dlt))/(2*dlt)
  end do

end program huber_smooth
