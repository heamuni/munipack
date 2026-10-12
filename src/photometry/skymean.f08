!
!  The mean level of sky by a clipped mean
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
!  The sky subtracted from a star must be the MEAN level of the sky:
!  the aperture sums the sky pixels, so any other location (the mode,
!  the median, a robust M-estimate) biases the counts by the difference
!  times the area. A sky of a few photo-electrons per pixel has a skewed,
!  discrete (Poisson) distribution, its robust means lie below its mean:
!  0.4 e-/pixel and the read noise 0.37 e- give -13% of the sky with the
!  Hampel estimator of rmean (BUGS P23), -30% at 0.1 e-/pixel.
!
!  The estimator is the mean of the pixels inside the interval
!
!     t - k*sig  <=  x  <=  t + q + sqrt(q**2 + (k*sig)**2),
!
!  iterated from the median; t is the current mean, k = 5, and sig is
!  the larger of the observed scatter of the pixels inside the interval
!  and the expected one by the noise model of the detector (when the gain
!  is given), sqrt((t - bias)/gain + readnoise**2) [ADU]. The upper limit
!  is the Bernstein bound of the Poisson tail: the probability that a
!  Poisson count exceeds its mean by more than q + sqrt(q**2 + (k*sig)**2),
!  q = k**2/6 photo-electrons (= k**2/(6*gain) ADU), is below exp(-k**2/2),
!  for any mean, even when sig is less than one electron (the median
!  absolute deviation of such data is zero, the clipping by sig alone
!  would reject every photon). So the sky photons are never clipped (the
!  bias of the truncation is below 0.001 e- per pixel), while stars,
!  cosmics and hot pixels, many sig above, are. The lower limit is k*sig
!  (the lower tail of the Poisson distribution is shorter than Normal).
!  The scatter by the data guards the clipping when the model is wrong
!  (no gain, frames processed); the model guards it when the data have
!  few distinct values.
!
!  Faint light within k*sig (the wings of neighbours) is not rejected:
!  its pixels should be masked before (aphot masks the detected stars).
!
!  The standard error is the observed scatter / sqrt(n); the scatter and
!  the number of used pixels are returned too.

module skymean

  use iso_fortran_env
  use titsio

  implicit none
  private

  integer, parameter :: dbl = REAL64

  ! the clipping [sigma]
  real(dbl), parameter, public :: sky_kclip = 5

  public :: sky_mean, sky_select, sky_noise_model

contains

  subroutine sky_mean(s,gain,rn2,bias,mean,stderr,stdsig,nused)

    ! s .. the sky pixels [ADU]
    ! gain [e-/ADU] (<= 0 when the noise model is unknown), rn2 the read
    !      noise squared [ADU**2], bias [ADU] the pedestal of the frame
    ! mean, stderr, stdsig [ADU]: the mean, its standard error and
    !      the scatter of pixels; nused: the number of the pixels used

    real(dbl), dimension(:), intent(in) :: s
    real(dbl), intent(in) :: gain, rn2, bias
    real(dbl), intent(out) :: mean, stderr, stdsig
    integer, intent(out) :: nused

    integer, parameter :: maxit = 100
    real(dbl), dimension(:), allocatable :: w
    real(dbl) :: t, tnew, sig, var, q, lo, hi, k
    integer :: n, m, mold, it

    n = size(s)
    k = sky_kclip
    nused = n
    if( n == 0 ) then
       mean = 0
       stderr = -1
       stdsig = -1
       return
    else if( n < 3 ) then
       mean = sum(s) / n
       stdsig = 0
       if( n > 1 ) stdsig = sqrt(sum((s - mean)**2) / (n - 1))
       stderr = stdsig / sqrt(real(n,dbl))
       return
    end if

    ! the initial estimate: the median and the scatter by the median
    ! absolute deviation (or by the model, when larger)
    allocate(w(n))
    w = s
    t = sky_select(w,(n+1)/2)
    w = abs(s - t)
    sig = max(1.4826_dbl*sky_select(w,(n+1)/2),model(t))

    ! one quantum: a photo-electron in ADU (1 ADU without the model)
    if( gain > 0 ) then
       q = k**2 / (6*gain)
    else
       q = k**2 / 6
    end if

    mold = -1
    var = 0
    do it = 1, maxit
       lo = t - k*sig
       hi = t + q + sqrt(q**2 + (k*sig)**2)
       m = count(lo <= s .and. s <= hi)
       if( m < 2 ) exit
       tnew = sum(s, lo <= s .and. s <= hi) / m
       var = sum((s - tnew)**2, lo <= s .and. s <= hi) / (m - 1)
       sig = max(sqrt(var),model(tnew))
       if( m == mold .and. abs(tnew - t) <= 1e-9_dbl*max(abs(t),sig) ) then
          t = tnew
          exit
       end if
       t = tnew
       mold = m
    end do

    mean = t
    nused = max(m,1)
    stdsig = sqrt(var)
    stderr = stdsig / sqrt(real(nused,dbl))

  contains

    real(dbl) function model(x)

      ! the scatter of a pixel of the mean x by the noise model

      real(dbl), intent(in) :: x

      if( gain > 0 ) then
         model = sqrt(max(x - bias,0.0_dbl)/gain + max(rn2,0.0_dbl))
      else
         model = 0
      end if

    end function model

  end subroutine sky_mean

  function sky_select(a,k) result(x)

    ! the k-th smallest element of a (Hoare's selection; a is reordered)

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

  end function sky_select

  subroutine sky_noise_model(filename,gain,rn2,bias,status)

    ! the noise model of the frame for the sky: the gain [e-/ADU] (GAIN;
    ! 1 for BUNIT = COUNTS; -1 when unknown), the read noise squared
    ! [ADU**2] (READNS [ADU] or RDNOISE [e-]; 0 when unknown) and the bias
    ! [ADU] (BIASLEV, the pedestal of artificial; 0 when unknown)

    character(len=*),intent(in) :: filename
    real(REAL64), intent(out) :: gain, rn2, bias
    integer, intent(in out) :: status

    type(fitsfiles) :: fits
    character(len=FLEN_VALUE) :: bunit
    real(REAL64) :: rn

    gain = -1
    rn2 = 0
    bias = 0
    if( status /= 0 ) return

    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) return

    call fits_write_errmark
    call fits_read_key(fits,FITS_KEY_GAIN,gain,status)
    if( status /= 0 .or. .not. (gain > 0) ) then
       call fits_clear_errmark
       status = 0
       gain = -1
    end if

    call fits_write_errmark
    call fits_read_key(fits,FITS_KEY_BUNIT,bunit,status)
    if( status == 0 ) then
       ! converted to electrons by phcorr
       if( bunit == 'COUNTS' ) gain = 1
    else
       call fits_clear_errmark
       status = 0
    end if

    call fits_write_errmark
    call fits_read_key(fits,FITS_KEY_READNS,rn,status)
    if( status == 0 ) then
       rn2 = rn**2
    else
       call fits_clear_errmark
       status = 0
       call fits_write_errmark
       call fits_read_key(fits,'RDNOISE',rn,status)
       if( status == 0 .and. gain > 0 ) then
          rn2 = (rn / gain)**2
       else
          call fits_clear_errmark
          status = 0
       end if
    end if

    call fits_write_errmark
    call fits_read_key(fits,'BIASLEV',bias,status)
    if( status /= 0 ) then
       call fits_clear_errmark
       status = 0
       bias = 0
    end if

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine sky_noise_model

end module skymean
