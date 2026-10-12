!
!  The excess scatter of residuals over their errors
!
!  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
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
!  A sample of residuals r(i) with the quoted errors e(i) scatters more
!  than the errors say when an unmodelled (systematic) effect is present.
!  The excess is described by one additional Gaussian scatter sigma, so
!  the errors of the residuals are sqrt(e**2 + sigma**2). The routine
!  excess_scatter estimates sigma by the maximum likelihood,
!
!      sum w(i)**2 * (r(i)**2/kappa - e(i)**2 - sigma**2) = 0,
!      w(i) = 1 / (e(i)**2 + sigma**2),
!
!  on the residuals within clip (default 3) of their errors (kappa is
!  the variance of the Normal distribution truncated at clip, so the
!  clipping does not bias sigma). The initial sigma is given by the
!  median absolute deviation, the clipping is iterated until the set of
!  residuals does not change. Sigma is zero when the residuals do not
!  scatter more than their errors. The residuals are supposed to have
!  zero mean (the residuals of a fit).
!

module exscatter

  use iso_fortran_env

  implicit none
  private

  integer, parameter, private :: dbl = REAL64

  public :: excess_scatter

contains

  subroutine excess_scatter(res,err,sigma,keep,clip)

    real(dbl), dimension(:), intent(in) :: res, err
    real(dbl), intent(out) :: sigma
    logical, dimension(:), intent(out), optional :: keep
    real(dbl), intent(in), optional :: clip

    real(dbl), parameter :: q50 = 0.6744897501960817_dbl
    real(dbl), dimension(:), allocatable :: u
    logical, dimension(:), allocatable :: mask, prev
    real(dbl) :: c, kappa, smax
    integer :: n, iter

    n = size(res)
    sigma = 0
    allocate(u(n),mask(n),prev(n))
    mask = err > 0
    if( present(keep) ) keep = mask
    if( count(mask) < 2 ) return

    c = 3
    if( present(clip) ) c = clip
    ! the variance of N(0,1) truncated at +-c
    kappa = 1 - 2*c*exp(-c**2/2)/sqrt(8*atan(1.0_dbl)) / erf(c/sqrt(2.0_dbl))

    smax = 10*maxval(abs(res),mask=mask) + tiny(smax)

    ! the initial estimate: the median of |res|/sqrt(err**2+sigma**2)
    ! equals the one of the Normal distribution
    if( madnorm(0.0_dbl) > q50 ) &
         sigma = bisect(madfun,0.0_dbl,smax)

    do iter = 1, 20
       prev = mask
       mask = err > 0
       where( mask ) mask = abs(res) <= c*sqrt(err**2 + sigma**2)
       if( count(mask) < 2 ) exit
       if( likefun(0.0_dbl) > 0 ) then
          sigma = bisect(likefun,0.0_dbl,smax)
       else
          sigma = 0
       end if
       if( iter > 1 .and. all(mask .eqv. prev) ) exit
    end do

    if( present(keep) ) keep = mask

  contains

    real(dbl) function madnorm(s)
      real(dbl), intent(in) :: s
      integer :: m
      real(dbl), dimension(:), allocatable :: a
      m = count(err > 0)
      a = pack(abs(res)/sqrt(err**2 + s**2), err > 0)
      call sort(a)
      if( mod(m,2) == 1 ) then
         madnorm = a(m/2+1)
      else
         madnorm = (a(m/2) + a(m/2+1)) / 2
      end if
    end function madnorm

    real(dbl) function madfun(s)
      real(dbl), intent(in) :: s
      madfun = madnorm(s) - q50
    end function madfun

    real(dbl) function likefun(s)
      real(dbl), intent(in) :: s
      likefun = sum((res**2/kappa - err**2 - s**2) / (err**2 + s**2)**2, &
           mask=mask)
    end function likefun

  end subroutine excess_scatter


  real(dbl) function bisect(f,a0,b0)

    ! the root of the decreasing function f in a0 .. b0, f(a0) > 0

    interface
       real(REAL64) function f(s)
         import :: REAL64
         real(REAL64), intent(in) :: s
       end function f
    end interface
    real(dbl), intent(in) :: a0, b0
    real(dbl) :: a, b, m
    integer :: i

    a = a0
    b = b0
    if( f(b) > 0 ) then
       bisect = b
       return
    end if
    do i = 1, 100
       m = (a + b) / 2
       if( f(m) > 0 ) then
          a = m
       else
          b = m
       end if
       if( b - a <= 1e-9_dbl*b ) exit
    end do
    bisect = (a + b) / 2

  end function bisect


  subroutine sort(a)

    ! heap sort

    real(dbl), dimension(:), intent(in out) :: a
    real(dbl) :: t
    integer :: n, i

    n = size(a)
    do i = n/2, 1, -1
       call sift(i,n)
    end do
    do i = n, 2, -1
       t = a(1)
       a(1) = a(i)
       a(i) = t
       call sift(1,i-1)
    end do

  contains

    subroutine sift(l,r)
      integer, intent(in) :: l, r
      integer :: i, j
      real(dbl) :: t
      i = l
      t = a(i)
      do
         j = 2*i
         if( j > r ) exit
         if( j < r ) then
            if( a(j+1) > a(j) ) j = j + 1
         end if
         if( a(j) <= t ) exit
         a(i) = a(j)
         i = j
      end do
      a(i) = t
    end subroutine sift

  end subroutine sort

end module exscatter
