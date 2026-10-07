!
!  Estimates on median base
!
!  Copyright © 2001 - 2018 F.Hroch (hroch@physics.muni.cz)
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
!
!

module madmean

  implicit none
  private

  ! print additional debuging results
  logical, parameter, private :: debug = .false.

  ! numerical precision of real numbers
  integer, parameter, private :: kind = selected_real_kind(15)

  ! 50% quantil of N(0,1)
  real(kind), parameter, private :: Q50 = 0.6745

  public :: madmed, madwmed

contains

  subroutine madmed(x,t,s)

    ! initial estimate of parameters

    use medians
    use quantilmean

    real(kind), dimension(:),intent(in) :: x
    real(kind), intent(out) :: t,s

    real(kind) :: mad,med
    integer :: n

    n = size(x)

    if( n > 50 ) then

       ! This choice of the threshold, used for median computation method,
       ! assumes that odd and even elements in sequence are the same (within
       ! requested precision). Note that qmed is the fastest known algorithm
       ! (~ 2*n) while median is simpler and slower (~ n*log(n)).

       med = qmedian(x)
       mad = qmedian(abs(x - med))

       ! correction for Normal distribution
       s = mad / q50
       t = med

    else if( n > 13 ) then

       ! correct way to compute median (two middle values)
       med = median(x)
       mad = median(abs(x - med))

       ! correction for Normal distribution
       s = mad / q50
       t = med

    else if( n > 2 ) then

       ! compute parameters by empirical CDF quantilies
       call qmean(x,t,stdsig=s)

    else if( n == 2 ) then

       t = (x(1) + x(2))/2
       s = abs(x(1) - x(2))/2

    else if( n == 1 ) then

       t = x(1)
       s = 0

    else

       t = 0
       s = 0

    end if

  end subroutine madmed

  subroutine madwmed(x,dx,t,s,d)

    ! This routine performs an initial estimate of both mean and scale
    ! using both median and absolute value.

    use medians
    use quantilmean

    real(kind), dimension(:),intent(in) :: x, dx
    real(kind), intent(out) :: t,s,d

    real(kind) :: mad,med
    integer :: n

    n = size(x)

    if( n > 50 ) then

       ! This choice of the threshold used for median computation method
       ! assumes that odd and even elements in sequence are the same (within
       ! requested precision). Note that qmed is the fastest known algorithm
       ! (~2*n) while median is simpler and slower (~n*log(n)).

       d = qmedian(dx)
       med = qmedian(x)
       mad = qmedian(abs(x - med)) / d

       ! correction for Normal distribution
       s = mad / q50
       t = med

    else if( n > 13 ) then

       ! correct way to compute median (two middle values)
       d = median(dx)
       med = median(x)
       mad = median(abs(x - med)) / d

       ! correction for Normal distribution
       s = mad / q50
       t = med

    else if( n > 2 ) then

       ! compute parameters by empirical CDF quantilies
       call qmean(dx,d,stdsig=s)
       call qmean(x,med,stdsig=s)
       call qmean(abs(x - med),mad,stdsig=s)
       mad = mad / d
       s = mad / q50
       t = med

    else if( n == 2 ) then

       d = (dx(1) + dx(2)) / 2
       t = (x(1) + x(2)) / 2
       s = abs(x(1) - x(2)) / d

    else if( n == 1 ) then

       d = dx(1)
       t = x(1)
       s = 1

    else

       d = 1
       t = 0
       s = 1

    end if

  end subroutine madwmed

end module madmean
