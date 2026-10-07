!
!  quick median - a quantile (median) algorithm with ~2*n
!
!  by  Wirth,N: Algorithm + Data Structure = Programs, Prentice-Hall, 1975
!
!  Copyright © 1997-2025 F.Hroch (hroch@physics.muni.cz)
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


module quickmedian_REAL64

  use iso_fortran_env
  implicit none
  private

  public :: qmed_REAL64

contains

  function qmed_REAL64(a,k) result(x)

    integer, intent(in) :: k
    real(REAL64), dimension(:), intent(in out) :: a
    real(REAL64) :: w,x
    integer :: l,r,i,j

    l = 1
    r = size(a)
    if( r == 0 ) then
       x = 0
       return
    end if

    x = a(1)
    do while( l < r )
       x = a(k)
       i = l
       j = r
       do
          do while( a(i) < x )
             i = i + 1
          enddo
          do while( x < a(j) )
             j = j - 1
          enddo
          if( i <= j ) then
             w = a(i)
             a(i) = a(j)
             a(j) = w
             i = i + 1
             j = j - 1
          endif
          if( i > j ) exit
       enddo
       if( j < k ) l = i
       if( k < i ) r = j
    enddo

  end function qmed_REAL64

end module quickmedian_REAL64
