!
!  The noise of the detection of stars: thresholds and shape limits
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
!  The threshold (BUGS P33)
!
!  find convolves the frame with a lowered Gaussian, h = sum w_k d_k
!  (sum w_k = 0, sum w_k G_k = 1), and looks for local maxima of h above
!  t sigma_h, sigma_h**2 = sigma**2 sum w_k**2 (sigma the scatter of a
!  sky pixel). The threshold t means a false-alarm probability of a pixel
!  Q(t) = erfc(t/sqrt(2))/2 only for Normal noise. A sky of a few
!  photo-electrons per pixel is Poisson: a pixel d = (N + R)/g + bias
!  [ADU] with N ~ Poisson(lam) electrons, the read noise R ~ N(0, rn**2)
!  and the gain g [e-/ADU]. The cumulant generating function of h in the
!  units of sigma_h, X = h/sigma_h, is (the convolution is a sum of
!  independent pixels)
!
!     K(u) = lam sum_k (exp(u y_k) - 1 - u y_k) + rhog u**2 / 2,
!
!  y_k = w_k / (g sigma_h), rhog = 1 - lam/(g sigma)**2 the fraction of
!  the variance which is not Poisson (read noise, quantization); its
!  cumulants are kappa_2 = 1, kappa_3 = lam sum y**3 (the skewness),
!  kappa_4 = lam sum y**4. On a Poisson sky of 0.53 e- per pixel with
!  the read noise 0.37 e- the skewness of the filtered sky is 0.105 for
!  FWHM 2.7 and 0.051 for 4.93, its excess kurtosis 0.22 and 0.08: the
!  tail is heavier than Normal, 5.7 (2.7) or 2.6 (4.93) times more
!  pixels exceed 4 sigma_h. The height threshold x is the one exceeded
!  with the probability Q(t): the tail by the saddle-point approximation
!  of Barndorff-Nielsen,
!
!     P(X > x) = Q(r*),   r* = r + ln(q/r) / r,
!     r = sqrt(2 (u x - K(u))),   q = u sqrt(K''(u)),   K'(u) = x,
!
!  is Q(t) when r*(u) = t. The approximation uses the exact cumulant
!  generating function and is accurate to the far tail: the probability
!  by an exponentially tilted Monte Carlo is within 1.3% of Q(t) at
!  t = 4 and 5 for skies of 0.13-200 e- (x = 4.59 and 6.04 for FWHM
!  2.7 on the sky above; the second-order Cornish-Fisher expansion
!  overcorrects below 1 e-: up to 5 times fewer pixels than Q(t)).
!  Without the gain the noise is taken Normal (x = t).
!
!  The expected number of noise peaks (--false-rate)
!
!  A peak of the noise is a local maximum of h above x. For a smooth
!  Gaussian random field the expected number of peaks above t in the
!  area A is the expected Euler characteristic of the excursion set
!  (Adler 1981; Worsley et al. 1996, Hum. Brain Mapp. 4, 58),
!
!     A lambda2 t exp(-t**2/2) / (2 pi)**(3/2),
!
!  lambda2 the variance of the derivative of the normalized field: by
!  the kernel, lambda2 = 2 (1 - rho1), rho1 the correlation of the
!  neighbouring pixels of h. The kernel of find is small and the field
!  is rough (rho1 = 0.08 for FWHM 2.7, 0.40 for 4.93): a peak is mostly
!  a single pixel above t, and the number of the pixels above t,
!  A Q(t) (Bonferroni), is the lower and better bound (Worsley et al.
!  1996: the lower of both). Normal noise, 2048**2: 120-135 peaks of the
!  sieve of find at t = 4 for FWHM 2-5 against A Q(4) = 132; 98 at
!  FWHM 10, where A Q is 1.3 times high; 17 per 1024**2 for FWHM 15
!  against 26 by the Euler characteristic (32 by A Q). The expected
!  count is an upper estimate (the shape tests reject some peaks).

module findnoise

  use iso_fortran_env

  implicit none
  private

  integer, parameter :: dbl = REAL64
  real(dbl), parameter :: pi = 3.14159265358979324_dbl

  public :: normal_tail, normal_deviate, height_threshold, noise_peaks, &
       rate_deviate, field_lambda2

  ! the statistics of a detection (reported in FIND)
  type, public :: findstats
     ! the threshold: the Normal deviate t, the height [sigma_h], the
     ! skewness and the excess kurtosis of the filtered sky, the expected
     ! number of noise peaks above the threshold in the frame
     real(dbl) :: t = 0, x = 0, skew = 0, kurt = 0, nfalse = 0
  end type findstats

contains

  real(dbl) function normal_tail(t)

    ! Q(t), the upper tail of the Normal distribution

    real(dbl), intent(in) :: t

    normal_tail = erfc(t / sqrt(2.0_dbl)) / 2

  end function normal_tail

  real(dbl) function normal_deviate(p)

    ! t of Q(t) = p (by bisection on the logarithm)

    real(dbl), intent(in) :: p
    real(dbl) :: a, b, t, lp
    integer :: i

    lp = log(p)
    a = -38
    b = 38
    do i = 1, 200
       t = (a + b) / 2
       if( log(normal_tail(t)) > lp ) then
          a = t
       else
          b = t
       end if
       if( b - a < 1e-12_dbl ) exit
    end do
    normal_deviate = (a + b) / 2

  end function normal_deviate

  subroutine height_threshold(y,lam,rhog,t,x,skew,kurt)

    ! The height x [sigma_h] of the convolution exceeded with the
    ! probability Q(t) by the noise of the sky: y the kernel in the units
    ! of sigma_h per electron, lam [e-] the Poisson mean of a pixel,
    ! rhog the fraction of the variance of a pixel which is Normal.
    ! skew, kurt are the skewness and the excess kurtosis of h.

    real(dbl), dimension(:), intent(in) :: y
    real(dbl), intent(in) :: lam, rhog, t
    real(dbl), intent(out) :: x, skew, kurt

    real(dbl) :: a, b, u, rs
    integer :: i

    skew = lam * sum(y**3)
    kurt = lam * sum(y**4)
    x = t
    if( .not. (lam > 0) .or. .not. (t > 0) .or. size(y) == 0 ) return

    ! bracket the saddle point (r* is increasing in u), then bisection
    a = 0
    b = t
    do i = 1, 60
       if( rstar(b) >= t ) exit
       a = b
       b = 2*b
    end do
    do i = 1, 100
       u = (a + b) / 2
       rs = rstar(u)
       if( rs < t ) then
          a = u
       else
          b = u
       end if
       if( b - a < 1e-10_dbl*b ) exit
    end do
    u = (a + b) / 2
    x = lam*sum(y*expm1(u*y)) + rhog*u

  contains

    real(dbl) function rstar(u)

      real(dbl), intent(in) :: u
      real(dbl) :: r, q, d2, e

      ! r**2/2 = u K'(u) - K(u), q = u sqrt(K''(u))
      e = lam*sum(exp(min(u*y,700.0_dbl))*(u*y - 1) + 1) + rhog*u**2/2
      d2 = lam*sum(y**2*exp(min(u*y,700.0_dbl))) + rhog
      r = sqrt(2*max(e,tiny(e)))
      q = u*sqrt(max(d2,tiny(d2)))
      if( u > 0 ) then
         rstar = r + log(q/r) / r
      else
         rstar = 0
      end if

    end function rstar

    elemental real(dbl) function expm1(z)

      ! exp(z) - 1 accurate for small z

      real(dbl), intent(in) :: z

      if( abs(z) < 1e-5_dbl ) then
         expm1 = z*(1 + z/2*(1 + z/3))
      else
         expm1 = exp(min(z,700.0_dbl)) - 1
      end if

    end function expm1

  end subroutine height_threshold

  real(dbl) function field_lambda2(w)

    ! The variance of the derivative of the convolution of white noise
    ! by the kernel w, normalized to the unit variance (finite
    ! differences of neighbouring pixels): 2 (1 - rho1).

    real(dbl), dimension(:,:), intent(in) :: w
    real(dbl) :: s2, cx, cy
    integer :: n1, n2

    n1 = size(w,1)
    n2 = size(w,2)
    s2 = sum(w**2)
    if( .not. (s2 > 0) ) then
       field_lambda2 = 2
       return
    end if
    cx = sum(w(1:n1-1,:)*w(2:n1,:)) / s2
    cy = sum(w(:,1:n2-1)*w(:,2:n2)) / s2
    field_lambda2 = 2*(1 - (cx + cy)/2)

  end function field_lambda2

  real(dbl) function noise_peaks(t,area,lambda2)

    ! The expected number of the noise peaks above t [sigma] in the area
    ! [pixels]: the lower of the number of pixels above t and the
    ! expected Euler characteristic of the excursion set.

    real(dbl), intent(in) :: t, area, lambda2

    noise_peaks = area * min(normal_tail(t), &
         lambda2*t*exp(-t**2/2) / (2*pi)**1.5_dbl)

  end function noise_peaks

  real(dbl) function rate_deviate(rate,area,lambda2)

    ! t at which the expected number of noise peaks is the rate
    ! (bisection; the number decreases with t above t = 1)

    real(dbl), intent(in) :: rate, area, lambda2
    real(dbl) :: a, b, t
    integer :: i

    a = 1
    b = 40
    do i = 1, 200
       t = (a + b) / 2
       if( log(noise_peaks(t,area,lambda2)) > log(rate) ) then
          a = t
       else
          b = t
       end if
       if( b - a < 1e-10_dbl ) exit
    end do
    rate_deviate = (a + b) / 2

  end function rate_deviate

end module findnoise
