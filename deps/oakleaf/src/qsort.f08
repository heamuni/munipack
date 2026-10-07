!
!  QUICKSORT - interfaces for QuickSort algorithm
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
!  QUICKSORT - recursive version of the QuickSort algorithm
!
!  by  Wirth,N: Algorithm + Data Structure = Programs, Prentice-Hall, 1975
!
!  Copyright © 1997-2011, 2017, 2023-2024 F.Hroch (hroch@physics.muni.cz)
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
!  along with OakLeaf.  If not, see <http://www.gnu.org/licenses/>.


module quicksort_REAL64

  use iso_fortran_env
  implicit none

  private
  public :: qsort_REAL64
  !, qsort_lepage

contains

  subroutine qsort_REAL64(a)

    real(REAL64), dimension(:), intent(in out) :: a
    integer :: n

    n = size(a)
    if( n <= 1 ) return

    call sort(1,n,a)

  end subroutine qsort_REAL64

  recursive subroutine sort(l,r,a)

    integer, intent(in) :: l,r
    real(REAL64), dimension(:), intent(in out) :: a

    integer :: i,j
    real(REAL64) :: x,w

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
          w = a(i)
          a(i) = a(j)
          a(j) = w
          i = i + 1
          j = j - 1
       endif

       if( i > j ) exit
    end do

    if( l < j ) call sort(l,j,a)
    if( i < r ) call sort(i,r,a)

  end subroutine sort

  subroutine qsort_lepage(a,V)

    real(REAL64), dimension(:), intent(in out) :: a
    integer, dimension(:), intent(in out) :: V
    integer :: n

    n = size(a)
    if( n <= 1 ) return

    call trid_lepage(1,n,a,V)

  end subroutine qsort_lepage

  recursive subroutine trid_lepage(l,r,a,V)

    integer,intent(in) :: l,r
    real(REAL64), dimension(:),intent(in out) :: a
    integer, dimension(:),intent(in out) :: V

    integer :: i,j,k
    real(REAL64) :: x,w

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
          w = a(i)
          a(i) = a(j)
          a(j) = w
          k = V(i)
          V(i) = V(j)
          V(j) = k
          i = i + 1
          j = j - 1
       endif

       if( i > j ) exit
    end do

    if( l < j ) call trid_lepage(l,j,a,V)
    if( i < r ) call trid_lepage(i,r,a,V)

  end subroutine trid_lepage


end module quicksort_REAL64
