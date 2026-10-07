!
!  stat -  a library of basic statistical estimates
!
!  Copyright © 2011, 2017 F.Hroch (hroch@physics.muni.cz)
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
!

module means

  implicit none

  ! precision of real numbers
  integer, parameter, private :: rp = selected_real_kind(15)

  interface mean
     module procedure mean_double, mean_single
  end interface mean

contains

  subroutine mean_double(x,t,dt)

    real(rp), dimension(:), intent(in) :: x
    real(rp), intent(out) :: t,dt
    integer :: n

    n = size(x)
    if( n <= 0 ) then
       t = 0.0_rp
       dt = 0.0_rp
    else if( n == 1 ) then
       t = x(1)
       dt = 0.0_rp
    else
       t = sum(x)/n
       dt = sqrt(sum((x - t)**2)/(n - 1))
    endif

  end subroutine mean_double

  subroutine mean_single(x,t,dt)

    real, dimension(:), intent(in) :: x
    real, intent(out) :: t,dt
    integer :: n

    n = size(x)
    if( n <= 0 ) then
       t = 0.0
       dt = 0.0
    else if( n == 1 ) then
       t = x(1)
       dt = 0.0
    else
       t = sum(x)/n
       dt = sqrt(sum((x - t)**2)/(n - 1))
    endif

  end subroutine mean_single

end module means
