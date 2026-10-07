!
!  QUICKSORT - recursive version of the QuickSort algorithm
!
!  by  Wirth,N: Algorithm + Data Structure = Programs, Prentice-Hall, 1975
!
!  Copyright © 1997-2011, 2017-8 F.Hroch (hroch@physics.muni.cz)
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


module quicksort

  ! precision of real numbers
  integer, parameter, private :: rp = selected_real_kind(15)

  private :: trid, tridi

contains

  subroutine qsort(a,idx)

    real(rp), dimension(:), intent(in out) :: a
    integer, dimension(:), optional, intent (in out) :: idx
    integer :: n

    n = size(a)
    if( n < 2 ) return

    if( present(idx) ) then
       if( n /= size(idx) ) stop 'qsort: n /= size(idx)'
       call tridi(1,n,a,idx)
    else
       call trid(1,n,a)
    end if

  end subroutine qsort


  recursive subroutine trid(l,r,a)

    integer,intent(in) :: l,r
    real(rp), dimension(:),intent(in out) :: a

    integer :: i,j
    real(rp) :: x,w     ! internal buffers

    i = l
    j = r
    x = a((l+r)/2)

    do
       do while ( a(i) < x )
          i = i + 1
       end do

       do while ( x < a(j) )
          j = j - 1
       end do

       if( i <= j )then
          w = a(i); a(i) = a(j); a(j) = w
          i = i + 1; j = j - 1
       endif

       if( i > j ) exit
    end do

    if( l < j ) call trid(l,j,a)
    if( i < r ) call trid(i,r,a)

  end subroutine trid

!--------------------------------------------------------------------------

  recursive subroutine tridi(l,r,a,idx)

    integer,intent(in) :: l, r
    integer, dimension(:), intent(in out) :: idx
    real(rp), dimension(:), intent(in out) :: a

    integer :: i,j,m
    real(rp) :: x,w     ! internal buffers

    i = l
    j = r
    x = a((l+r)/2)

    do
       do while ( a(i) < x )
          i = i + 1
       end do

       do while ( x < a(j) )
          j = j - 1
       end do

       if( i <= j )then
          w = a(i); a(i) = a(j); a(j) = w
          m = idx(i); idx(i) = idx(j); idx(j) = m
          i = i + 1; j = j - 1
       endif

       if( i > j ) exit
    end do

    if( l < j ) call tridi(l,j,a,idx)
    if( i < r ) call tridi(i,r,a,idx)

  end subroutine tridi


end module quicksort
