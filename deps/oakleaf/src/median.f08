!
!  Medians - common methods
!
!  Copyright © 2016-2025 F.Hroch (hroch@physics.muni.cz)
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
! Medians convenience functions:
! * qmedian - an interface to qmed(), ~2*n
! * median -  a general median on base of quicksort (fast ~n*log(n))

module medians_REAL64

  use iso_fortran_env
  implicit none
  private

  public :: median_REAL64, qmedian_REAL64

contains

  function qmedian_REAL64(b,q) result(qmedian)

    use quickmedian

    real(REAL64) :: qmedian
    real(REAL64),dimension(:),intent(in) :: b
    real(REAL64), intent(in), optional :: q
    real(REAL64), dimension(:), allocatable :: a

    integer :: n, nmed, status
    character(len=80) :: errmsg

    n = size(b)

    if( n == 0 ) then
       qmedian = 0
       return
    else if( n == 1 ) then
       qmedian = b(1)
       return
    end if

    if( present(q) ) then
       nmed = nint(q*n)
    else ! q = 0.5 (median)
       nmed = n/2 + 1
    end if

    ! this is right just only for odd-th elements of the sequence,
    ! we're ignore the right way, assuming huge dataset, n > 50

    allocate(a,source=b,stat=status,errmsg=errmsg)
    if( status /= 0 ) then
       write(error_unit,*) 'qmedian() allocate: ',trim(errmsg)
       error stop 'qmedian has no enough memory storage.'
    end if

    qmedian = qmed(a,nmed)

    deallocate(a)

  end function qmedian_REAL64


  ! Median on base of quicksort

  function median_REAL64(x) result(median)

    use quicksort

    real(REAL64) :: median
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), allocatable :: y
    integer :: n, status
    character(len=80) :: errmsg

    n = size(x)

    if( n == 0 ) then
       median = 0
       return
    else if( n == 1 ) then
       median = x(1)
       return
    end if

    allocate(y,source=x,stat=status,errmsg=errmsg)
    if( status /= 0 ) then
       write(error_unit,*) 'median() allocate: ',trim(errmsg)
       error stop 'median has no enough memory storage.'
    end if

    call qsort(y)
    median = (y(n/2+1) + y((n+1)/2)) / 2

    ! check
!    n2 = n / 2
!    write(*,*) n,n2,mod(n,2)
!    if( mod(n,2) == 0 ) then
!       ! even
!       median = (y(n2+1) + y(n2)) / 2.0_REAL64
!    else
!       ! odd
!       median = y(n2+1)
!    end if

    deallocate(y)

  end function median_REAL64

end module medians_REAL64
