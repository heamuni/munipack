!
!  Quantile means
!
!  Copyright © 2016-8, 2023-2025 F.Hroch (hroch@physics.muni.cz)
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

module quantilmean_REAL64

  ! Subroutine qmean() determines statistical parameters by quantiles.
  !
  ! qmean should be called as:
  !
  !  call qmean(data,mean,stderr,stdsig,flag)
  !
  ! on input:
  !   data - array of data values to be estimated
  !
  ! on output are estimated:
  !   mean - robust mean
  !   stderr (optional) - standard error
  !   stdsig (optional) - standard deviation
  !   verbose (optional)
  !   flag (optional):
  !          0 == success, results are reliable
  !          1 == basic estimates, zero or a single point available
  !          6 == an unspecified error
  !          8 == failed to allocate memory
  !
  ! The results gives the robust estimate by median of the true value X
  ! of the sample x, with 68% probability, into the interval
  !
  !         mean - stderr  <  X  <  mean + stderr
  !
  ! Data distribution is approximated as Normal distribution
  !
  !          N(mean,stdsig).
  !
  ! Robust estimators has been prepared on base of
  !  * Hogg: Robustness in Statistics
  !  * Hubber: Robust Statistics
  !  * Davidon: Order statistics
  !  * my experiences

  use iso_fortran_env

  implicit none
  private

  real(REAL64), parameter :: q34 = 0.6745
  real(REAL64), parameter :: q58 = 0.3186

  public :: qmean_REAL64, byQ_REAL64, ecdf_REAL64, quantile_REAL64, &
       qbracket_REAL64, quantilefun_REAL64

contains

  subroutine qmean_REAL64(data,mean,stderr,stdsig,sortby,sigby,verbose,flag)

    use quicksort
    use quickmedian

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(out) :: mean
    real(REAL64), intent(out), optional :: stderr,stdsig
    character(len=*), intent(in), optional :: sortby,sigby
    logical, intent(in), optional :: verbose
    integer, intent(out), optional :: flag

    real(REAL64), dimension(:), allocatable :: x
    real(REAL64) :: t, s
    integer :: n, m, iflag
    character(len=88) :: tmeth, smeth, errmsg
    logical :: verb

    n = size(data)
    verb = .false.
    tmeth = "qmed"
    smeth = "MAD"
    if( present(flag) ) flag = 6
    if( present(verbose) ) verb = verbose
    if( present(sortby) ) tmeth = sortby
    if( present(sigby) ) smeth = sigby

    if( n == 0 ) then
       if( present(flag) ) flag = 5
       return
    else if( n == 1 ) then
       mean = data(1)
       if( present(stderr) )   stderr = 0
       if( present(stdsig) )   stdsig = 0
       if( present(flag) )     flag = 2
       return
    end if

    if( smeth == 'byQ' ) then
       m = (n*(n-1))/2
       if( .not. (m > n) ) then
          if( verb ) write(error_unit,*) 'qmean(): overflow for byQ'
          return
       end if
    else
       m = n
    end if
    allocate(x(m),stat=iflag,errmsg=errmsg)
    if( iflag /= 0 ) then
       if( verb ) write(error_unit,*) 'qmean(): ',trim(errmsg)
       if( present(flag) ) flag = 8
       return
    end if

    x(1:n) = data(1:n)
    if( tmeth == "qsort" ) then
       call qsort(x(1:n))
       t = (x((n+1)/2) + x(n/2+1)) / 2.0_REAL64
    else
       t = qmed(x(1:n),n/2)
    end if

    if( present(stderr) .or. present(stdsig) ) then

       if( smeth == "MAD" ) then
          x(1:n) = abs(data(1:n) - t)
          if( tmeth == "qsort" ) then
             call qsort(x(1:n))
             s = (x((n+1)/2) + x(n/2+1)) / (2*q34)
          else
             s = qmed(x(1:n),n/2) / q34
          end if
       else
          s =  byQ_REAL64(data(1:n),x)
       end if

    end if

    mean = t
    if( present(flag) ) flag = 0
    if( present(stdsig) ) stdsig = s
    if( present(stderr) ) stderr = s / sqrt(real(n))
    deallocate(x)

  end subroutine qmean_REAL64

  subroutine ecdf_REAL64(data,x,p)

    ! empirical cumulative distribution function

    use quicksort

    real(REAL64),dimension(:),intent(in) :: data
    real(REAL64),dimension(:),intent(out) :: x,p
    integer :: i,n
    real(REAL64) :: h

    n = size(data)
    h = 1.0_REAL64 / real(n + 1,REAL64)

    x = data
    call qsort(x)
    p = [ (i*h, i=1,n) ]

  end subroutine ecdf_REAL64

  function quantile_REAL64(q,x,y) result(t)

    ! Estimates q-quantile, this is inverse of CDF from linearly interpolated data
    ! IMPORTANT, the interpolation on tail data is not recommended at all.

    real(REAL64), intent(in) :: q
    real(REAL64), dimension(:), intent(in) :: x,y
    real(REAL64) :: t,dy,h,r
    integer :: n,m,low,high

    n = size(x)

    if( n == 0 ) then
       t = 0
       return
    else if( n == 1 ) then
       t = x(1)
       return
    end if

    if( q <= y(1) ) then
       t = x(1)
    else if( q >= y(n) ) then
       t = x(n)
    else

       ! we're assuming, the input CDF has equidistant steps in y-axis
       h = real(1,REAL64) / real(n + 1,REAL64)
       r = q / h
       m = nint(r)
!       write(*,*) h,r,m
       if( abs(m - r) < 10*epsilon(r) ) then
          t = x(m)
       else
          low = int(r)
          high = low + 1

          dy = y(high) - y(low)
          if( -abs(dy) > 10*epsilon(dy) ) then
             ! inverse by linear interpolation
             t = (x(high) - x(low))/dy*(q - y(low)) + x(low)
          else
             ! nearly singular
             t = (x(high) + x(low)) / 2
          end if
       end if
    end if
!    write(*,*) high,low,real(x(high)),real(x(low)),real(y(high)),real(y(low)),real(t)

  end function quantile_REAL64

  subroutine quantilefun_REAL64(data,q,x,flag)

    ! obsolete

    real(REAL64), dimension(:), intent(in) :: data, q
    real(REAL64), dimension(:), intent(out) :: x
    integer, intent(out) :: flag

    real(REAL64), dimension(:), allocatable :: xcdf,ycdf
    integer :: i,n,stat
    character(len=80) :: errmsg

    n = size(data)

    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'quantilefun: ',trim(errmsg)
       flag = 5
       return
    end if

    ! empirical cummulative distribution function of data
    call ecdf_REAL64(data,xcdf,ycdf)

!    write(*,*) real(xcdf)
!    write(*,*) real(ycdf)

    ! quantiles as inverse of ecdf
    do i = 1, size(q)
       x(i) = quantile_REAL64(q(i),xcdf,ycdf)
    end do

    deallocate(xcdf,ycdf)
    flag = 0

  end subroutine quantilefun_REAL64

  subroutine qbracket_REAL64(data,xmin,xmax,flag)

    ! obsolete

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(out) :: xmin, xmax
    integer, intent(out), optional :: flag

    real(REAL64), dimension(2), parameter :: q = [0.25, 0.75]
    real(REAL64), dimension(size(q)) :: x
    integer :: iflag

    call quantilefun_REAL64(data,q,x,iflag)
    xmin = x(1)
    xmax = x(2)
    if( present(flag) ) flag = iflag

  end subroutine qbracket_REAL64

  function byQ_REAL64(data,q) result(byQ)

    use quickmedian

    real(REAL64) :: byQ
    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), dimension(:), intent(out) :: q
    integer :: i,j,n,m

    n = size(data)
    m = size(q)
    if( (n*(n-1))/2 > m .or. n < 4 ) error stop 'byQ: insufficient array size'

    m = 0
    do i = 1, n
       do j = 1,i-1
          m = m + 1
          q(m) = abs(data(i) - data(j))
       end do
    end do

    byQ = qmed(q(1:m),m / 4) / (sqrt(2.0)*q58)

  end function byQ_REAL64

end module quantilmean_REAL64
