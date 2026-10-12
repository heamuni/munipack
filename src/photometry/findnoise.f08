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
!
!  The shape of a point source (BUGS P32)
!
!  The sharpness of DAOFIND, sharp = (d_0 - mean of d_k, k /= 0)/h
!  (over the circle of the kernel, sky subtracted), and the roundness,
!  round = 2 (hx - hy)/(hx + hy) (the heights of the marginal fits), are
!  linear functions of the pixels divided by a height; to the first order
!
!     var(sharp) = sum_k (a_k - s w_k)**2 v_k / h**2,
!     var(round) = sum_k (cx_k - cy_k)**2 v_k / (h hm)**2,
!
!  a_k the coefficients of the numerator of sharp (1 in the centre,
!  -1/(n-1) else), cx, cy the coefficients of the marginal heights,
!  v_k = sigma**2 + h P_k / g the noise of the pixel by the sky and the
!  photon noise of the star of the height h, P the profile of the star
!  (unit height of the filter), hm = sum cx P, s the expected sharpness.
!  P is a Gaussian of the FWHM of the filter integrated over the pixels,
!  averaged over the sub-pixel position of the star (its sharpness is
!  0.396 for FWHM 4.93 and 0.526 for 2.7, the pixel phase changes them
!  by 0.001 and 0.005). On QUVIK-like frames (gain 0.101) the photon
!  noise of the star dominates: the scatter of the sharpness of real
!  stars at 4-5, 8-11, 25-50, 50-100 sigma_h is 0.45, 0.31, 0.16, 0.11
!  (MAD; 0.50, 0.33, 0.17, 0.12 by the model), so 39, 28, 13 and 4%
!  of them fell outside the fixed window 0.2 - 1.0.

module findnoise

  use iso_fortran_env

  implicit none
  private

  integer, parameter :: dbl = REAL64
  real(dbl), parameter :: pi = 3.14159265358979324_dbl

  public :: normal_tail, normal_deviate, height_threshold, noise_peaks, &
       rate_deviate, field_lambda2, shape_model, shape_noise

  ! the statistics of a detection (reported in FIND)
  type, public :: findstats
     ! the threshold: the Normal deviate t, the height [sigma_h], the
     ! skewness and the excess kurtosis of the filtered sky, the expected
     ! number of noise peaks above the threshold in the frame
     real(dbl) :: t = 0, x = 0, skew = 0, kurt = 0, nfalse = 0
     ! the shape tests: k [sigma] (0: the fixed windows), the reference
     ! sharpness and roundness (by the bright candidates or the model),
     ! the number of the bright candidates (0: by the model)
     real(dbl) :: shapesig = 0, sharpref = 0, roundref = 0
     integer :: nref = 0
     logical :: sharpnoise = .false., roundnoise = .false.
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

  subroutine shape_model(fwhm,mask,w,cx,cy,prof,sharp0,hm)

    ! The profile of a point source of the FWHM of the filter integrated
    ! over the pixels and averaged over its sub-pixel position, of the
    ! unit height of the filter (prof), its sharpness and the mean height
    ! of its marginal fits (hm); mask is the circle of the kernel w,
    ! cx, cy the coefficients of the heights of the marginal fits

    real(dbl), intent(in) :: fwhm
    logical, dimension(:,:), intent(in) :: mask
    real(dbl), dimension(:,:), intent(in) :: w, cx, cy
    real(dbl), dimension(:,:), intent(out) :: prof
    real(dbl), intent(out) :: sharp0, hm

    integer, parameter :: nph = 8
    real(dbl), dimension(size(w,1)) :: ex
    real(dbl), dimension(size(w,2)) :: ey
    real(dbl), dimension(size(w,1),size(w,2)) :: p
    real(dbl) :: s, dx, dy, h
    integer :: i, j, k, l, mx, my

    mx = size(w,1) / 2 + 1
    my = size(w,2) / 2 + 1
    s = sqrt(2.0_dbl) * fwhm / (2*sqrt(2*log(2.0_dbl)))
    prof = 0
    do k = 1, nph
       dx = (k - 0.5_dbl) / nph - 0.5_dbl
       do i = 1, size(w,1)
          ex(i) = pixel(i - mx - dx)
       end do
       do l = 1, nph
          dy = (l - 0.5_dbl) / nph - 0.5_dbl
          do j = 1, size(w,2)
             ey(j) = pixel(j - my - dy)
          end do
          do j = 1, size(w,2)
             p(:,j) = ex*ey(j)
          end do
          h = sum(w*p)
          if( h > 0 ) prof = prof + p / h
       end do
    end do
    prof = prof / nph**2
    ! the centre minus the mean of the others in the circle
    sharp0 = prof(mx,my) - &
         (sum(prof,mask) - prof(mx,my)) / max(count(mask) - 1,1)
    hm = (sum(cx*prof) + sum(cy*prof)) / 2

  contains

    real(dbl) function pixel(z)

      ! the Gaussian integrated over a pixel centred at z

      real(dbl), intent(in) :: z

      pixel = (erf((z + 0.5_dbl)/s) - erf((z - 0.5_dbl)/s)) / 2

    end function pixel

  end subroutine shape_model

  subroutine shape_noise(h,sig2,gain,sref,mask,w,cx,cy,prof,hm,ssharp,sround)

    ! The noise (the standard deviations) of sharp and round of a point
    ! source of the height h [ADU] of the filter: sig2 the variance of
    ! a sky pixel [ADU**2], gain [e-/ADU] (<= 0: the sky noise only),
    ! sref the expected sharpness, the profile by shape_model

    real(dbl), intent(in) :: h, sig2, gain, sref, hm
    logical, dimension(:,:), intent(in) :: mask
    real(dbl), dimension(:,:), intent(in) :: w, cx, cy, prof
    real(dbl), intent(out) :: ssharp, sround

    real(dbl), dimension(size(w,1),size(w,2)) :: v, a
    integer :: mx, my, n

    if( .not. (h > 0) ) then
       ssharp = huge(ssharp)
       sround = huge(sround)
       return
    end if

    mx = size(w,1) / 2 + 1
    my = size(w,2) / 2 + 1
    n = count(mask)
    v = sig2
    if( gain > 0 ) v = v + h*max(prof,0.0_dbl) / gain
    where( mask )
       a = -1.0_dbl / max(n - 1,1)
    elsewhere
       a = 0
    end where
    a(mx,my) = 1
    ssharp = sqrt(sum((a - sref*w)**2*v)) / h
    sround = sqrt(sum((cx - cy)**2*v)) / (h*max(hm,epsilon(hm)))

  end subroutine shape_noise

end module findnoise
