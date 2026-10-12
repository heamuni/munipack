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
!  The estimator is the mean of the pixels clipped by the noise model of
!  the detector (GAIN, READNS or RDNOISE, BIASLEV; lib/clipmean.f08,
!  shared with kombine and dark): the interval
!
!     t - k*sig  <=  x  <=  t + q + sqrt(q**2 + (k*sig)**2),
!
!  iterated from the median; t is the current mean, k = 5, and sig is
!  the larger of the observed scatter of the pixels inside the interval
!  and the expected one by the noise model of the detector (when the gain
!  is given), sqrt((t - bias)/gain + readnoise**2) [ADU]; q = k**2/6
!  photo-electrons (= k**2/(6*gain) ADU), the Bernstein bound of the
!  Poisson tail, so the sky photons are never clipped, while stars,
!  cosmics and hot pixels, many sig above, are.
!
!  Faint light within k*sig (the wings of neighbours) is not rejected:
!  its pixels should be masked before (aphot masks the detected stars).
!
!  The standard error is the observed scatter / sqrt(n); the scatter and
!  the number of used pixels are returned too.

module skymean

  use iso_fortran_env
  use titsio
  use clipmean, only: clip_kappa, clipped_mean

  implicit none
  private

  integer, parameter :: dbl = REAL64

  ! the clipping [sigma]
  real(dbl), parameter, public :: sky_kclip = clip_kappa

  public :: sky_mean, sky_noise_model

contains

  subroutine sky_mean(s,gain,rn2,bias,mean,stderr,stdsig,nused)

    ! s .. the sky pixels [ADU]
    ! gain [e-/ADU] (<= 0 when the noise model is unknown), rn2 the read
    !      noise squared [ADU**2], bias [ADU] the pedestal of the frame
    ! mean, stderr, stdsig [ADU]: the mean, its standard error and
    !      the scatter of pixels; nused: the number of the pixels used
    !      (stderr = stdsig = -1 without pixels)

    real(dbl), dimension(:), intent(in) :: s
    real(dbl), intent(in) :: gain, rn2, bias
    real(dbl), intent(out) :: mean, stderr, stdsig
    integer, intent(out) :: nused

    if( size(s) == 0 ) then
       mean = 0
       stderr = -1
       stdsig = -1
       nused = 0
       return
    end if

    call clipped_mean(s,[gain],[rn2],[bias],mean,stderr,stdsig,nused)

  end subroutine sky_mean

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
