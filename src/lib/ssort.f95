!
!  selectsort - a simple sorting algorithm
!
!  by  Wirth,N: Algorithm + Data Structure = Programs, Prentice-Hall, 1975
!
!  Copyright © 2012 F.Hroch (hroch@physics.muni.cz)
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


module selectsort

  implicit none

  ! precision of real numbers
  integer, parameter, private :: rp = selected_real_kind(15)

contains

  subroutine ssort(a)

    real(rp), dimension(:), intent(in out) :: a
    real(rp) :: x
    integer :: n,i,j,k

    n = size(a)
    do i = 1, n - 1
       k = i
       x = a(i)
       do j = i + 1, n
          if( a(j) < x ) then
             k = j
             x = a(j)
          end if
       end do
       a(k) = a(i)
       a(i) = x
    end do

  end subroutine ssort

end module selectsort
