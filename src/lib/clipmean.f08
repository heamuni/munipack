!
!  The mean clipped by the noise model of the detector
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
!  The mean of a sample of pixel values whose outliers (stars, cosmic
!  rays, hot pixels, satellites) are clipped at bounds given by the noise
!  model of the detector, so that the photons of a faint background are
!  never clipped. One routine for the means of Munipack that must be the
!  MEAN (not a mode, a median or a robust M-estimate) of a few
!  photo-electrons per pixel, whose skewed, discrete (Poisson)
!  distribution puts the robust means below the mean (BUGS P23, P30, P31):
!
!   * the sky of a star (aphot), of a frame (find), of an annulus
!     (psfphot): the pixels of a ring, one noise model (skymean.f08);
!   * a pixel of a stack of frames (kombine, dark): one value per frame,
!     each with its own noise model and the weight of its exposure.
!
!  The estimator is the mean (weighted by w) of the values v inside
!
!     t - k*s  <=  v  <=  t + q + sqrt(q**2 + (k*s)**2),
!
!  iterated from the median until the set of the values does not change;
!  t is the current mean, k = 5, s the larger of the observed scatter of
!  the used values and the one of a value by the noise model,
!
!     sqrt(max(t*w - bias, 0)/gain + rn2) / w     [the units of v],
!
!  gain [e-/ADU], rn2 the read noise squared [ADU**2] and bias [ADU] the
!  level of no photons (the pedestal of the frame less any subtracted
!  level); v*w is the value in the raw units, w is 1 by default. The upper
!  limit is the Bernstein bound of the Poisson tail: the probability that
!  a Poisson count exceeds its mean by more than q + sqrt(q**2 + (k*s)**2),
!  q = k**2/6 photo-electrons (= k**2/(6*gain*w) in the units of v), is
!  below exp(-k**2/2) for any mean, even when s is less than one electron
!  (the median absolute deviation of such data is zero, the clipping by
!  s alone would reject every photon). So the photons of the background
!  are never clipped (the bias of the truncation is below 0.001 e- per
!  pixel), while stars, cosmics and hot pixels, many s above, are. The
!  lower limit is k*s (the lower tail of the Poisson distribution is
!  shorter than Normal). The scatter by the data guards the clipping when
!  the model is wrong (no gain, processed frames), the model guards it
!  when the data have few distinct values. Without the gain (gain <= 0),
!  s is the observed scatter and q = k**2/6 (in the units of v). Fewer
!  than three values are not clipped.
!
!  Light within k*s of the background (the wings of stars) is not
!  rejected: such pixels should be masked before (aphot masks stars).
!
!  The outputs: the mean, its standard error (the observed scatter of
!  the used values / sqrt(their number); with the weights, the error of
!  the weighted mean sqrt(sum w (v - mean)**2 / (m - 1) / sum w)), the
!  scatter of the used values (unweighted, m - 1) and their number m;
!  the error and the scatter are zero for m < 2.

module clipmean

  use iso_fortran_env

  implicit none
  private

  integer, parameter :: dbl = REAL64

  ! the clipping [sigma]
  real(dbl), parameter, public :: clip_kappa = 5

  ! the work arrays of up to nsmall values are automatic (a pixel of
  ! a stack), allocated otherwise (the pixels of a ring or a frame)
  integer, parameter :: nsmall = 256

  public :: clipped_mean, clip_select

contains

  subroutine clipped_mean(v,gain,rn2,bias,mean,stderr,stdsig,nused,w)

    ! v .. the values
    ! gain, rn2, bias .. the noise model (see above): one element for all
    !      the values (the pixels of a frame), or one per value (a stack)
    ! w .. the weights (relative exposures of a stack; 1 when absent)
    ! mean, stderr, stdsig, nused .. see above

    real(dbl), dimension(:), intent(in) :: v, gain, rn2, bias
    real(dbl), intent(out) :: mean, stderr, stdsig
    integer, intent(out) :: nused
    real(dbl), dimension(:), intent(in), optional :: w

    integer :: n
    logical :: common, weighted

    n = size(v)
    weighted = present(w)
    common = size(gain) == 1 .and. size(rn2) == 1 .and. size(bias) == 1 &
         .and. .not. weighted
    if( .not. ((size(gain) == 1 .or. size(gain) == n) .and. &
         (size(rn2) == 1 .or. size(rn2) == n) .and. &
         (size(bias) == 1 .or. size(bias) == n)) ) &
         error stop 'clipped_mean: the noise model does not match the values.'
    if( weighted ) then
       if( size(w) /= n ) error stop 'clipped_mean: the weights do not match.'
    end if

    mean = 0
    stderr = 0
    stdsig = 0
    nused = 0
    if( n == 0 ) return

    if( n <= nsmall ) then
       block
         logical, dimension(n) :: keep, prev
         real(dbl), dimension(n) :: s, q, work
         call iterate(keep,prev,s,q,work)
       end block
    else
       block
         logical, dimension(:), allocatable :: keep, prev
         real(dbl), dimension(:), allocatable :: s, q, work
         if( common ) then
            allocate(keep(n),prev(n),s(1),q(1),work(n))
         else
            allocate(keep(n),prev(n),s(n),q(n),work(n))
         end if
         call iterate(keep,prev,s,q,work)
       end block
    end if

  contains

    subroutine iterate(keep,prev,s,q,work)

      ! keep, prev the sets of the used values; s, q the scatter and the
      ! quantum by the model (one element in the common case); work for
      ! the medians

      logical, dimension(:), intent(out) :: keep, prev
      real(dbl), dimension(:), intent(out) :: s, q, work

      integer, parameter :: maxit = 100
      real(dbl), parameter :: k = clip_kappa
      real(dbl) :: t, so, var, sw, c, lo, hi
      integer :: m, it, i

      keep = .true.
      if( n >= 3 ) then

         ! the median and the scatter by the median absolute deviation
         work = v
         t = median_of(work)
         work = abs(v - t)
         so = 1.4826_dbl*median_of(work)

         ! one quantum: a photo-electron in the units of v (1 without the
         ! model) times k**2/6
         if( common ) then
            if( gain(1) > 0 ) then
               q(1) = k**2 / (6*gain(1))
            else
               q(1) = k**2 / 6
            end if
         else
            do i = 1, n
               if( gain(ix(gain,i)) > 0 ) then
                  q(i) = k**2 / (6*gain(ix(gain,i))*wx(i))
               else
                  q(i) = k**2 / 6
               end if
            end do
         end if

         do it = 1, maxit
            prev = keep
            if( common ) then
               s(1) = max(so,model(t,1))
               lo = t - k*s(1)
               hi = t + q(1) + sqrt(q(1)**2 + (k*s(1))**2)
               keep = lo <= v .and. v <= hi
            else
               do i = 1, n
                  s(i) = max(so,model(t,i))
                  keep(i) = t - k*s(i) <= v(i) .and. &
                       v(i) <= t + q(i) + sqrt(q(i)**2 + (k*s(i))**2)
               end do
            end if
            m = count(keep)
            if( m < 2 ) then
               keep = prev
               exit
            end if
            if( weighted ) then
               t = sum(w*v,keep) / sum(w,keep)
               so = sqrt(sum(w*(v - t)**2,keep) / sum(w,keep) * m / (m - 1))
            else
               t = sum(v,keep) / m
               var = sum((v - t)**2,keep) / (m - 1)
               so = sqrt(var)
            end if
            if( it > 1 .and. all(keep .eqv. prev) ) exit
         end do

      end if

      ! the mean of the used values, its error and their scatter
      m = count(keep)
      nused = m
      if( weighted ) then
         sw = sum(w,keep)
         t = sum(w*v,keep) / sw
         if( m > 1 ) then
            c = sum(w*(v - t)**2,keep) / (m - 1)
            stderr = sqrt(c / sw)
            stdsig = sqrt(sum((v - t)**2,keep) / (m - 1))
         end if
      else
         t = sum(v,keep) / m
         if( m > 1 ) then
            stdsig = sqrt(sum((v - t)**2,keep) / (m - 1))
            stderr = stdsig / sqrt(real(m,dbl))
         end if
      end if
      mean = t

    end subroutine iterate

    real(dbl) function model(t,i)

      ! the scatter of the value i of the mean t by the noise model

      real(dbl), intent(in) :: t
      integer, intent(in) :: i
      integer :: ig

      ig = ix(gain,i)
      if( gain(ig) > 0 ) then
         model = sqrt(max(t*wx(i) - bias(ix(bias,i)),0.0_dbl)/gain(ig) + &
              max(rn2(ix(rn2,i)),0.0_dbl)) / wx(i)
      else
         model = 0
      end if

    end function model

    integer function ix(a,i)
      ! the element of a parameter of the value i (one for all, or one each)
      real(dbl), dimension(:), intent(in) :: a
      integer, intent(in) :: i
      if( size(a) == 1 ) then
         ix = 1
      else
         ix = i
      end if
    end function ix

    real(dbl) function wx(i)
      ! the weight of the value i
      integer, intent(in) :: i
      if( weighted ) then
         wx = w(i)
      else
         wx = 1
      end if
    end function wx

  end subroutine clipped_mean


  real(dbl) function median_of(a) result(x)

    ! the median of a, reordered (the mean of the two middle values of
    ! an even count)

    real(dbl), dimension(:), intent(in out) :: a
    integer :: n

    n = size(a)
    if( mod(n,2) == 1 ) then
       x = clip_select(a,n/2+1)
    else
       x = clip_select(a,n/2)
       ! the selection leaves the larger values above n/2
       x = (x + minval(a(n/2+1:n))) / 2
    end if

  end function median_of


  function clip_select(a,k) result(x)

    ! the k-th smallest element of a (Hoare's selection; a is reordered:
    ! a(1:k-1) <= a(k) <= a(k+1:))

    real(dbl), dimension(:), intent(in out) :: a
    integer, intent(in) :: k
    real(dbl) :: x, p, tmp
    integer :: l, r, i, j

    l = 1
    r = size(a)
    do while( l < r )
       p = a((l + r) / 2)
       i = l
       j = r
       do
          do while( a(i) < p )
             i = i + 1
          end do
          do while( p < a(j) )
             j = j - 1
          end do
          if( i <= j ) then
             tmp = a(i)
             a(i) = a(j)
             a(j) = tmp
             i = i + 1
             j = j - 1
          end if
          if( i > j ) exit
       end do
       if( j < k ) l = i
       if( k < i ) r = j
    end do
    x = a(k)

  end function clip_select

end module clipmean
