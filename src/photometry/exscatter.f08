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
!  The routine mixture_mean estimates the mean of values y(i) with
!  the errors e(i) together with the excess scatter sigma, robustly to
!  outliers, by the likelihood of a mixture (Hogg, Bovy & Lang 2010,
!  arXiv:1008.4686): a value is an inlier, N(mu, e(i)**2 + sigma**2), with
!  the probability 1 - f, or an outlier of a broad distribution,
!  N(mu, e(i)**2 + sigma**2 + B**2), with the probability f. The parameters
!  are found by the EM algorithm: the probabilities q(i) of the values to
!  be inliers, then mu as the mean weighted by q(i)/(e(i)**2 + sigma**2)
!  (and by (1 - q(i))/(e(i)**2 + sigma**2 + B**2)), sigma by the
!  restricted (REML) likelihood of the inliers,
!
!      sum u(i) w(i) r(i)**2 = sum u(i) (1 - h(i)),
!      u(i) = q(i) w(i), w(i) = 1/(e(i)**2 + sigma**2), h(i) = u(i)/sum u,
!
!  (r(i) = y(i) - mu, h(i) the leverage of the mean, which removes the
!  bias of few values), f as the mean of 1 - q(i), at most 1/2, and B**2
!  as the excess variance of the outliers, at least 9 (median(e**2) +
!  sigma**2) and bmin**2 (the outliers broader than three times the
!  inliers; B starts at 1 mag, bmin is 0.1 mag in the natural logarithm
!  of a flux by default). The start is robust: the median and sigma by
!  the median absolute deviation. Outliers then have a negligible weight
!  in both mu and sigma, there is no clipping (the values in the tails of
!  the inliers keep their weight), and sigma is zero when the inliers do
!  not scatter more than their errors.
!
!  The routine mixture_fit is the same likelihood for a polynomial,
!  y = beta(0) + beta(1) x + beta(2) x**2 (of the order 0, 1, 2; the colour
!  terms of the calibration, colterm.f08): the mean is replaced by
!  the weighted least squares with the weights q/v + (1-q)/(v + B**2),
!  the leverage by h(i) = u(i) x(i)' (X' U X)**(-1) x(i), the start of
!  the slope by the median of the slopes of pairs (Theil-Sen). The error
!  of x(i), ex(i), enters by the effective variance, e(i)**2 + g'**2 ex**2
!  - 2 g' cov(y,x) (g' the slope of the polynomial at x(i)). The covariance
!  of beta is (X' U X)**(-1). mixture_mean is mixture_fit of the order 0.
!
!  The routine mixture_fit_multi is the same for several colours (at most
!  ncolmax), y = beta(0) + sum_j beta(2j-1) x_j + beta(2j) x_j**2 (the order
!  of every colour 1 or 2; a passband depending on two colours, QUVIK
!  260-360 nm against SDSS u by u - g and g - r). The effective variance has
!  the gradient of the polynomial by all the colours, the covariances of
!  the colours of a star (u - g and g - r share g: -err(g)**2) and their
!  covariances with y; the start of the slopes is the least absolute
!  deviations (iteratively reweighted).
!
!  The observed colours scatter about the true ones by their errors, which
!  biases the slopes towards zero (the regression dilution, by the ratio
!  lambda of the variances of the true and the observed colours), and
!  the error of a colour containing the band of y (u - g against u) is in
!  both, a bias of cov/var(x). Two colours as correlated as u - g and g - r
!  of the stars (0.96) amplify the dilution to 4% and 10% (k = -0.3, +0.2),
!  more than the errors of 30-100 stars; one colour of faint stars (SDSS u
!  fainter than 20.5) by 2-4% of k, 1-2 errors of 100 stars. With errcorr,
!  the weighted least squares is corrected for the errors of the colours
!  (the moments, Fuller 1987, Measurement Error Models, 2.2; the terms of
!  the order 2 by the unbiased estimates of the powers for Gaussian errors;
!  see wls), one or several colours; capped when the colours are hardly
!  wider than their errors. It assumes the classical model of the errors,
!  the catalogue colour the true one plus its error.
!
!  Few stars (n - p + 1 < 16, p the number of the coefficients): f and B
!  cannot be estimated from few values, and the likelihood has two modes,
!  a deviating value either an outlier or a part of the scatter; EM from
!  the robust start, which prefers the first, biased sigma low by 15-25%
!  for 3-7 values and by 5-10% for 10-15 values (Monte Carlo). For few
!  values, f is fixed at 0.1 and B at 1 mag (at least 3 times
!  the inliers), EM runs from the robust start and from the plain
!  restricted likelihood of all values (q = 1), and the solution of
!  the larger posterior is taken: the restricted likelihood of the mixture,
!  sum ln((1-f) N(r,v) + f N(r,v+B**2)) - 1/2 ln det X'UX, with a weak
!  prior of sigma, half-normal of the scale 0.1 mag (a scatter of 0.2 mag
!  is unlikely, an outlier of 0.4 mag is not). Sigma is then nearly as
!  unbiased as by the plain REML (its square within 1-2% for 3-7 values)
!  and the outliers far from the others are still found; moderate ones
!  (a few sigma) cannot be told from the scatter. The results of 16 and
!  more values are those of the mixture above. With the prior, one degree
!  of freedom (n = p + 1) is enough for an estimate of sigma (3 stars and
!  a linear colour term); n = p, and n = p + 1 with a fixed sigma, is
!  the plain weighted fit.
!
!  The errors of beta (with sigma estimated): the covariance (X' U X)**(-1)
!  takes sigma as known, while few stars determine it poorly (from 3 stars
!  sigma = 0 is as likely as twice its true value), so the errors are too
!  small exactly when sigma dominates: the zero point within its error in
!  54-63% of sets of 3-10 stars for sigma >= 0.03 mag (Monte Carlo). The
!  covariance is the one of Knapp & Hartung (2003, Stat. Med. 22, 2693;
!  the random-effects meta-regression, the same model: known errors of
!  the values and an estimated excess variance) in the modification of
!  Roever, Knapp & Friede (2015, BMC Med. Res. Meth. 15, 99),
!
!      C = t**2 max(1, Q) (X' U X)**(-1),  Q = sum u r**2 / nu,
!
!  nu = sum q - p the degrees of freedom of the calibration stars, t
!  the quantile of Student's t of nu (at least 1) of the probability 68.27%
!  (1.84 for nu = 1, 1.32 for 2, 1.14 for 4, 1.05 for 10, 1.02 for 30). The
!  scale Q is the scatter of the calibration stars about the fit in
!  the units of their errors (with sigma): when sigma is underestimated
!  (from few stars, often zero) the stars scatter more and Q > 1; the floor
!  1 keeps the errors at least those by the errors of the stars. The t
!  quantile is the distribution of the mean divided by an error estimated
!  from nu degrees of freedom (exact for equal errors). The coverage of
!  the zero point is then 65-70% for 3-100 stars and sigma 0.01-0.06 mag
!  (54-69% before), of the slope of a colour term 66-76% for 5-100 stars;
!  for sigma = 0, 70-87% (few stars cannot tell sigma = 0 from 0.03 mag,
!  so the intervals that cover 0.03 cover more at 0). With 30 and more
!  stars the errors grow by 1-6% (100: 1-3%). A fixed sigma (fixsig >= 0) keeps
!  (X' U X)**(-1).
!
!  The interval of sigma (optional) is given by the profile of
!  the restricted likelihood of the calibration stars (q fixed),
!
!      -1/2 (sum q ln v + ln det X'UX + sum u r**2),
!
!  the values of sigma where it falls by 1/2 from its value at the estimate
!  (68%); the lower bound is zero when the likelihood at sigma = 0 is
!  within 1/2, so the interval is an upper limit for sigma = 0.
!

module exscatter

  use iso_fortran_env

  implicit none
  private

  integer, parameter, private :: dbl = REAL64

  public :: excess_scatter, mixture_mean, mixture_fit, mixture_fit_multi

  ! the colours of mixture_fit_multi (at most), the coefficients of
  ! the polynomial (the intercept, linear and quadratic terms)
  integer, parameter, public :: ncolmax = 3
  integer, parameter :: pmax = 1 + 2*ncolmax

  ! few stars, n - p + 1 < nfew (p the number of the coefficients):
  ! the fraction of outliers fixed at ffew, their width at 1 mag; the two
  ! modes of the likelihood compared with a prior of sigma, half-normal
  ! of the scale sfew [mag]
  integer, parameter :: nfew = 16
  real(dbl), parameter :: ffew = 0.1_dbl, sfew = 0.1_dbl

  ! the correction for the errors of the colours (wls): the corrected
  ! variance of every combination of the terms at least lammin of
  ! the observed one (the cap), the quadratic terms corrected (eivquad)
  real(dbl), parameter :: lammin = 0.5_dbl
  logical, parameter :: eivquad = .true.

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


  subroutine mixture_mean(y,err,mu,dmu,sigma,fout,member,bmin,sigint,nfit)

    ! The mean of y with the errors err and the excess scatter sigma,
    ! robust to outliers (see above); dmu is the error of mu, fout the
    ! fraction of outliers, member (optional) the probabilities of
    ! the values to be inliers (0 for err <= 0), bmin (optional) the
    ! smallest width B of the distribution of outliers, sigint (optional)
    ! the 68% interval of sigma (-1 when not estimated), nfit (optional)
    ! the number of the values used.

    real(dbl), dimension(:), intent(in) :: y, err
    real(dbl), intent(out) :: mu, dmu, sigma, fout
    real(dbl), dimension(:), intent(out), optional :: member
    real(dbl), intent(in), optional :: bmin
    real(dbl), dimension(2), intent(out), optional :: sigint
    integer, intent(out), optional :: nfit

    real(dbl), dimension(size(y),1) :: x0
    real(dbl), dimension(size(y),1,1) :: xx0
    real(dbl), dimension(0:2) :: beta
    real(dbl), dimension(3,3) :: cov

    x0 = 0
    xx0 = 0
    call mixfit(y,err,x0,x0,xx0,x0,[0],-1.0_dbl,beta,cov,sigma,fout,member, &
         bmin,sigint,nfit,dmu)
    mu = beta(0)

  end subroutine mixture_mean


  subroutine mixture_fit(y,e,x,ex,cxy,order,fixsig,beta,cov,sigma,fout,member, &
       bmin,sigint,nfit,errcorr,lambda,eivstat)

    ! The polynomial y = beta(0) + beta(1) x + beta(2) x**2
    ! (beta(order+1:) = 0) with the excess scatter sigma by the likelihood
    ! of a mixture (see above). y = ln(ph/ct), e its error, x the colour
    ! minus a reference (centred for the numerics), ex the error of
    ! the colour, cxy the covariance of y and x, fixsig < 0: sigma by
    ! the restricted likelihood, >= 0: sigma fixed at the value; cov
    ! the covariance of beta, fout the fraction of outliers, member
    ! (optional) the probabilities of the stars to be calibration stars
    ! (0 for e <= 0), bmin (optional) the smallest width B of
    ! the distribution of outliers, sigint (optional) the 68% interval of
    ! sigma (-1 when not estimated, sigma for a fixed one), nfit (optional)
    ! the number of the stars used, errcorr (optional, .false. by default)
    ! the correction for the errors of x (see wls), lambda (optional)
    ! the reliability ratio of x (1 without the errors), eivstat
    ! (optional) 0 when the correction is not applied, 1 applied, 2
    ! applied with the cap (see wls). A value is used when e > 0, ex >= 0
    ! and y, x finite.

    real(dbl), dimension(:), intent(in) :: y, e, x, ex, cxy
    integer, intent(in) :: order
    real(dbl), intent(in) :: fixsig
    real(dbl), dimension(0:2), intent(out) :: beta
    real(dbl), dimension(3,3), intent(out) :: cov
    real(dbl), intent(out) :: sigma, fout
    real(dbl), dimension(:), intent(out), optional :: member
    real(dbl), intent(in), optional :: bmin
    real(dbl), dimension(2), intent(out), optional :: sigint
    integer, intent(out), optional :: nfit
    logical, intent(in), optional :: errcorr
    real(dbl), intent(out), optional :: lambda
    integer, intent(out), optional :: eivstat

    real(dbl), dimension(size(y),1,1) :: xx0
    real(dbl), dimension(1) :: lam
    logical :: ec

    xx0 = 0
    ec = .false.
    if( present(errcorr) ) ec = errcorr
    call mixfit(y,e,reshape(x,[size(x),1]),reshape(ex,[size(ex),1]),xx0, &
         reshape(cxy,[size(cxy),1]),[order],fixsig,beta,cov,sigma,fout,member, &
         bmin,sigint,nfit,errcorr=ec,lambda=lam,eivstat=eivstat)
    if( present(lambda) ) lambda = lam(1)

  end subroutine mixture_fit


  subroutine mixture_fit_multi(y,e,x,ex,exx,cxy,order,fixsig,beta,cov,sigma, &
       fout,member,bmin,sigint,nfit,errcorr,lambda,eivstat)

    ! The polynomial of several colours (at most ncolmax),
    !
    !   y = beta(0) + sum_j beta(2j-1) x(:,j) + beta(2j) x(:,j)**2,
    !
    ! (beta(2j) = 0 for the order 1 of the colour j, order(j)) with
    ! the excess scatter sigma by the likelihood of a mixture, as
    ! mixture_fit: x(:,j) the colours minus a reference, ex(:,j) their
    ! errors, exx(:,j,k) the covariances of the colours j /= k of every
    ! star (the colours sharing a band; the diagonal is not used),
    ! cxy(:,j) the covariances of y and the colours; beta(0:2*ncol) and
    ! cov(1+2*ncol,1+2*ncol) in this order (the elements of the missing
    ! quadratic terms zero), lambda(ncol) (optional) the reliability
    ! ratios of the colours. A value is used when e > 0, all ex >= 0
    ! and y, x finite. The other arguments as mixture_fit.

    real(dbl), dimension(:), intent(in) :: y, e
    real(dbl), dimension(:,:), intent(in) :: x, ex, cxy
    real(dbl), dimension(:,:,:), intent(in) :: exx
    integer, dimension(:), intent(in) :: order
    real(dbl), intent(in) :: fixsig
    real(dbl), dimension(0:), intent(out) :: beta
    real(dbl), dimension(:,:), intent(out) :: cov
    real(dbl), intent(out) :: sigma, fout
    real(dbl), dimension(:), intent(out), optional :: member
    real(dbl), intent(in), optional :: bmin
    real(dbl), dimension(2), intent(out), optional :: sigint
    integer, intent(out), optional :: nfit
    logical, intent(in), optional :: errcorr
    real(dbl), dimension(:), intent(out), optional :: lambda
    integer, intent(out), optional :: eivstat

    real(dbl), dimension(size(x,2)) :: lam
    logical :: ec

    if( size(x,2) > ncolmax ) stop 'mixture_fit_multi: too many colours.'
    ec = .false.
    if( present(errcorr) ) ec = errcorr
    call mixfit(y,e,x,ex,exx,cxy,order,fixsig,beta,cov,sigma,fout,member, &
         bmin,sigint,nfit,errcorr=ec,lambda=lam,eivstat=eivstat)
    if( present(lambda) ) lambda = lam

  end subroutine mixture_fit_multi


  subroutine mixfit(y,e,x,ex,exx,cxy,order,fixsig,beta,cov,sigma,fout,member, &
       bmin,sigint,nfit,dbeta0,errcorr,lambda,eivstat)

    ! The implementation of mixture_fit, mixture_fit_multi and
    ! mixture_mean; dbeta0 (optional) is the error of beta(0). For
    ! the order 0 (the mean), the weighted mean and its error are computed
    ! directly (not by the inverse of a 1x1 matrix). The design matrix
    ! has the columns 1, x(:,1), [x(:,1)**2], x(:,2), ... (the column k
    ! the power cpow(k) of the colour cind(k)); one colour is computed
    ! exactly as before the colours were generalised (the robust start by
    ! Theil-Sen, the effective variance), several colours start by
    ! the least absolute deviations.

    real(dbl), dimension(:), intent(in) :: y, e
    real(dbl), dimension(:,:), intent(in) :: x, ex, cxy
    real(dbl), dimension(:,:,:), intent(in) :: exx
    integer, dimension(:), intent(in) :: order
    real(dbl), intent(in) :: fixsig
    real(dbl), dimension(0:), intent(out) :: beta
    real(dbl), dimension(:,:), intent(out) :: cov
    real(dbl), intent(out) :: sigma, fout
    real(dbl), dimension(:), intent(out), optional :: member
    real(dbl), intent(in), optional :: bmin
    real(dbl), dimension(2), intent(out), optional :: sigint
    integer, intent(out), optional :: nfit
    real(dbl), intent(out), optional :: dbeta0
    logical, intent(in), optional :: errcorr
    real(dbl), dimension(:), intent(out), optional :: lambda
    integer, intent(out), optional :: eivstat

    real(dbl), parameter :: q50 = 0.6744897501960817_dbl
    real(dbl), parameter :: fmax = 0.5_dbl, tol = 1e-10_dbl
    ! 1 mag in the natural logarithm of flux
    real(dbl), parameter :: bmag = 0.4_dbl*2.302585092994046_dbl
    real(dbl), dimension(:), allocatable :: yy, e0, xx, e2, r, q, v, w
    real(dbl), dimension(:,:), allocatable :: xm, xc, ex2, cc
    real(dbl), dimension(:,:,:), allocatable :: cxx
    logical, dimension(size(y)) :: mask
    real(dbl), dimension(pmax) :: th, th1, th2, colsc
    real(dbl), dimension(pmax,pmax) :: ainv
    real(dbl), dimension(ncolmax) :: xs, lamc
    integer, dimension(pmax) :: cind, cpow
    real(dbl) :: vb, vbmin, smax, f, f2, s2, l1, emin, e2med, lmax
    integer :: n, p, nc, j, k, eivst
    logical :: few, eiv, haserr

    beta = 0
    cov = 0
    sigma = max(fixsig,0.0_dbl)
    fout = 0
    lamc = 1
    eivst = 0
    if( present(member) ) member = 0
    if( present(dbeta0) ) dbeta0 = 0
    if( present(sigint) ) sigint = -1
    if( present(lambda) ) lambda = 1
    if( present(eivstat) ) eivstat = 0

    nc = size(x,2)
    mask = e > 0 .and. abs(y) < huge(y)
    do j = 1, nc
       mask = mask .and. ex(:,j) >= 0 .and. abs(x(:,j)) < huge(x)
    end do
    n = count(mask)
    if( present(nfit) ) nfit = n
    if( n == 0 ) return

    yy = pack(y,mask)
    e0 = pack(e,mask)
    allocate(xc(n,nc),ex2(n,nc),cc(n,nc),cxx(n,nc,nc))
    do j = 1, nc
       xc(:,j) = pack(x(:,j),mask)
       ex2(:,j) = pack(ex(:,j),mask)**2
       cc(:,j) = pack(cxy(:,j),mask)
       do k = 1, nc
          cxx(:,j,k) = pack(exx(:,j,k),mask)
       end do
    end do
    xx = xc(:,1)
    allocate(e2(n),r(n),q(n),v(n),w(n))

    ! the columns of the design matrix: the order is limited by the stars
    ! and by the range of the colours (a colour of a single value has no
    ! term)
    cind = 0
    cpow = 0
    p = 1
    do j = 1, nc
       xs(j) = maxval(abs(xc(:,j)))
       if( .not. (maxval(xc(:,j)) - minval(xc(:,j)) > &
            1e-6_dbl*max(xs(j),1.0_dbl)) ) cycle
       do k = 1, min(order(j),2)
          if( p < pmax ) then
             p = p + 1
             cind(p) = j
             cpow(p) = k
          end if
       end do
    end do
    p = max(min(p,n),1)
    if( p == 1 ) xs = 1
    xs = max(xs,epsilon(xs))
    colsc = 1
    do k = 2, p
       colsc(k) = xs(cind(k))**cpow(k)
    end do

    allocate(xm(n,p))
    xm(:,1) = 1
    do k = 2, p
       xm(:,k) = xc(:,cind(k))**cpow(k)
    end do
    ! the weighted least squares corrected for the errors of the colours
    ! (errcorr, see wls): with the errors of the colours (or their
    ! covariances with y) given; several colours as before without them
    ! (the correction by zero)
    haserr = any(cpow(2:p) == 1) .and. (any(ex2 > 0) .or. any(cc /= 0))
    eiv = .false.
    if( present(errcorr) ) eiv = errcorr .and. any(cpow(2:p) == 1) .and. &
         (nc > 1 .or. haserr)

    q = 1
    th = 0

    if( n < p + 1 .or. (fixsig >= 0 .and. n < p + 2) ) then
       ! too few stars for the scatter and outliers: the weighted fit
       call effvar(th)
       v = e2 + sigma**2
       w = 1 / v
       call wls(w,th,ainv)
       call effvar(th)
       v = e2 + sigma**2
       w = 1 / v
       call wls(w,th,ainv)
       call layout(th,ainv)
       if( present(dbeta0) ) dbeta0 = 1 / sqrt(sum(w))
       if( present(member) ) member = merge(1.0_dbl,0.0_dbl,mask)
       if( present(lambda) ) lambda = lamc(1:nc)
       if( present(eivstat) ) eivstat = eivst
       return
    end if

    vb = bmag**2
    vbmin = (0.1_dbl*bmag)**2
    if( present(bmin) ) vbmin = bmin**2

    ! few stars: the fraction and the width of the outliers are fixed
    few = fixsig < 0 .and. n - p + 1 < nfew

    ! the robust start: the slope by the median of the slopes of pairs
    ! (several colours: by the least absolute deviations), the intercept
    ! by the median, the scatter by the median absolute deviation of
    ! the normalised residuals
    th = 0
    if( nc == 1 ) then
       if( p >= 2 ) th(2) = theilsen(xx,yy)
       th(1) = median(yy - th(2)*xx)
    else
       if( p >= 2 ) call ladfit(th)
       th(1) = median(yy - matmul(xm(:,2:p),th(2:p)))
    end if
    call effvar(th)
    e2med = median(e2)
    emin = sqrt(minval(e2))
    r = yy - matmul(xm(:,1:p),th(1:p))
    smax = 10*maxval(abs(r)) + tiny(smax)
    if( fixsig < 0 ) then
       sigma = 0
       if( madnorm(0.0_dbl) > q50 ) sigma = bisect(madfun,0.0_dbl,smax)
    end if
    f = 0.1_dbl

    if( few ) then
       ! the likelihood of few stars has two modes (a deviating star as
       ! an outlier, or as a part of the scatter); EM from the robust start
       ! and from the plain restricted likelihood of all stars, the mode of
       ! the larger posterior (the likelihood with the prior of sigma)
       th2 = th
       s2 = sigma
       f = ffew
       call emloop(th,sigma,f)
       l1 = mixlike(th,sigma,f) - (sigma/(sfew*bmag))**2/2
       q = 1
       call remlfit(th2,s2)
       f2 = ffew
       call emloop(th2,s2,f2)
       if( mixlike(th2,s2,f2) - (s2/(sfew*bmag))**2/2 >= l1 ) then
          th = th2
          sigma = s2
       end if
       vb = max(bmag**2,9*(e2med + sigma**2))
    else
       call emloop(th,sigma,f)
    end if
    fout = f

    call effvar(th)
    call estep(th,sigma,fout)
    w = q / (e2 + sigma**2)
    if( .not. (sum(w) > 0) ) w = 1 / (e2 + sigma**2)
    call wls(w,th1,ainv)
    call layout(th,ainv)
    if( present(dbeta0) ) dbeta0 = 1 / sqrt(sum(w))
    if( present(lambda) ) lambda = lamc(1:nc)
    if( present(eivstat) ) eivstat = eivst
    ! the uncertainty of the estimated sigma in the errors of beta
    if( fixsig < 0 ) call smallsample(th)
    ! the fraction of outliers found
    if( few ) fout = sum(1 - q) / n

    if( present(member) ) member = unpack(q,mask,0.0_dbl)

    ! the interval of sigma by the profile of the restricted likelihood
    ! of the calibration stars (q and the effective variances fixed):
    ! the values where it decreases by 1/2 from its value at sigma
    if( present(sigint) ) then
       if( fixsig >= 0 ) then
          sigint = sigma
       else
          call reml_interval(sigint)
       end if
    end if

  contains

    subroutine layout(th,ainv)
      ! the coefficients and their covariance in the order of beta,
      ! 0 (the intercept), 2j-1 (x(:,j)), 2j (x(:,j)**2)
      real(dbl), dimension(:), intent(in) :: th
      real(dbl), dimension(:,:), intent(in) :: ainv
      integer, dimension(pmax) :: l
      integer :: i, k
      l(1) = 0
      do k = 2, p
         l(k) = 2*(cind(k) - 1) + cpow(k)
      end do
      do k = 1, p
         beta(l(k)) = th(k)
         do i = 1, p
            cov(l(i)+1,l(k)+1) = ainv(i,k)
         end do
      end do
    end subroutine layout

    subroutine smallsample(th)
      ! the covariance of beta (and dbeta0) with the uncertainty of
      ! the estimated sigma (see above): scaled by the scatter of
      ! the calibration stars about the fit th (weights w, at least 1)
      ! and by the square of the quantile of Student's t of nu degrees
      ! of freedom (at least 1)
      real(dbl), dimension(:), intent(in) :: th
      real(dbl) :: nu, qhk, fac
      nu = sum(q) - p
      if( .not. (nu > 0.5_dbl) ) return
      r = yy - matmul(xm(:,1:p),th(1:p))
      qhk = max(sum(w*r**2) / nu, 1.0_dbl)
      fac = tquant68(max(nu,1.0_dbl))**2 * qhk
      cov = fac * cov
      if( present(dbeta0) ) dbeta0 = sqrt(fac) * dbeta0
    end subroutine smallsample

    subroutine ladfit(th)
      ! the least absolute deviations of the polynomial (iteratively
      ! reweighted least squares), the start of several colours
      real(dbl), dimension(:), intent(out) :: th
      real(dbl), dimension(pmax) :: t1
      real(dbl), dimension(pmax,pmax) :: ai
      real(dbl), dimension(n) :: wl, rl
      real(dbl) :: dl
      integer :: iter
      logical :: eiv0
      eiv0 = eiv
      eiv = .false.
      wl = 1
      th = 0
      call wls(wl,th,ai)
      do iter = 1, 100
         rl = yy - matmul(xm(:,1:p),th(1:p))
         dl = max(1e-6_dbl*median(abs(rl)),1e-10_dbl)
         wl = 1 / max(abs(rl),dl)
         call wls(wl,t1,ai)
         if( all(abs(t1(1:p) - th(1:p))*colsc(1:p) <= 1e-9_dbl) ) then
            th = t1
            exit
         end if
         th = t1
      end do
      eiv = eiv0
    end subroutine ladfit

    subroutine emloop(th,sigma,f)
      ! the EM algorithm from th, sigma, f; the fraction f and the width B
      ! of the outliers fixed for few stars
      real(dbl), dimension(:), intent(in out) :: th
      real(dbl), intent(in out) :: sigma, f
      real(dbl), dimension(pmax) :: th1
      real(dbl) :: sig1, f1, sq
      integer :: iter, k
      logical :: conv

      if( few ) vb = max(bmag**2,9*(e2med + sigma**2))

      do iter = 1, 1000

         ! the effective variances by the current slope, E-step
         call effvar(th)
         call estep(th,sigma,f)

         ! M-step: beta, sigma (REML), f
         v = e2 + sigma**2
         w = q/v + (1 - q)/(v + vb)
         th1 = th
         call wls(w,th1,ainv)
         r = yy - matmul(xm(:,1:p),th1(1:p))
         smax = 10*maxval(abs(r)) + tiny(smax)
         if( fixsig >= 0 ) then
            sig1 = fixsig
         else if( remlfun(0.0_dbl) > 0 ) then
            sig1 = bisect(remlfun,0.0_dbl,smax)
         else
            sig1 = 0
         end if
         if( few ) then
            f1 = ffew
            vb = max(bmag**2,9*(e2med + sig1**2))
         else
            f1 = min(sum(1 - q)/n, fmax)
            sq = sum(1 - q)
            if( sq > 0.5_dbl ) &
                 vb = max(sum((1 - q)*(r**2 - e2 - sig1**2))/sq, &
                 9*(e2med + sig1**2), vbmin)
         end if

         conv = abs(sig1 - sigma) <= tol*emin .and. abs(f1 - f) <= tol
         do k = 1, p
            conv = conv .and. abs(th1(k) - th(k))*colsc(k) <= tol*emin
         end do
         th = th1
         sigma = sig1
         f = f1
         if( conv ) exit
      end do

    end subroutine emloop

    subroutine remlfit(th,sigma)
      ! the weighted fit and sigma by the restricted likelihood of all
      ! the stars (q = 1, no outliers), iterated from th, sigma
      real(dbl), dimension(:), intent(in out) :: th
      real(dbl), intent(in out) :: sigma
      real(dbl), dimension(pmax) :: th1
      real(dbl) :: sig1
      integer :: iter, k
      logical :: conv
      do iter = 1, 1000
         call effvar(th)
         v = e2 + sigma**2
         w = q / v
         th1 = th
         call wls(w,th1,ainv)
         r = yy - matmul(xm(:,1:p),th1(1:p))
         smax = 10*maxval(abs(r)) + tiny(smax)
         if( remlfun(0.0_dbl) > 0 ) then
            sig1 = bisect(remlfun,0.0_dbl,smax)
         else
            sig1 = 0
         end if
         conv = abs(sig1 - sigma) <= tol*emin
         do k = 1, p
            conv = conv .and. abs(th1(k) - th(k))*colsc(k) <= tol*emin
         end do
         th = th1
         sigma = sig1
         if( conv ) exit
      end do
    end subroutine remlfit

    real(dbl) function mixlike(th,s,f)
      ! the restricted log-likelihood of the mixture (q by the E-step):
      ! sum ln((1-f) N(r,v) + f N(r,v+B**2)) - 1/2 ln det X'UX
      real(dbl), dimension(:), intent(in) :: th
      real(dbl), intent(in) :: s, f
      real(dbl), dimension(pmax,pmax) :: am
      real(dbl), dimension(n) :: u
      integer :: i, j
      call effvar(th)
      call estep(th,s,f)
      mixlike = sum(log((1 - f)*exp(-r**2/(2*v))/sqrt(v) + &
           f*exp(-r**2/(2*(v + vb)))/sqrt(v + vb)))
      u = q / v
      am = 0
      do i = 1, p
         do j = 1, p
            am(i,j) = sum(u*xm(:,i)*xm(:,j))
         end do
      end do
      mixlike = mixlike - logdet(am,p) / 2
    end function mixlike

    subroutine reml_interval(sint)
      ! sint(1) = 0 when the likelihood at zero is within 1/2 (an upper
      ! limit sint(2) for sigma = 0); sint(2) = -1 when the likelihood
      ! does not fall by 1/2 below 10 mag (the calibration stars, sum q,
      ! hardly more than the coefficients)
      real(dbl), dimension(2), intent(out) :: sint
      real(dbl) :: d, a
      logical :: found
      integer :: i
      lmax = remlike(sigma)
      ! the upper bound: expanding steps, then bisection
      d = max(sigma,0.1_dbl*sqrt(e2med),tiny(d))
      a = sigma
      found = .false.
      do i = 1, 200
         if( remlike(sigma + d) <= lmax - 0.5_dbl ) then
            found = .true.
            exit
         end if
         a = sigma + d
         d = 2*d
         if( sigma + d > 10*bmag ) exit
      end do
      if( found ) then
         sint(2) = bisect(hifun,a,sigma + d)
      else
         sint(2) = -1
      end if
      sint(1) = 0
      if( sigma > 0 ) then
         if( remlike(0.0_dbl) < lmax - 0.5_dbl ) &
              sint(1) = sigma - bisect(lofun,0.0_dbl,sigma)
      end if
    end subroutine reml_interval

    real(dbl) function hifun(s)
      real(dbl), intent(in) :: s
      hifun = remlike(s) - (lmax - 0.5_dbl)
    end function hifun

    real(dbl) function lofun(t)
      real(dbl), intent(in) :: t
      lofun = remlike(sigma - t) - (lmax - 0.5_dbl)
    end function lofun

    real(dbl) function remlike(s)
      ! the restricted log-likelihood of the calibration stars,
      ! -1/2 (sum q ln v + ln det X'UX + sum u r**2), u = q/v, v = e2 + s**2,
      ! r the residuals of the weighted fit at s
      real(dbl), intent(in) :: s
      real(dbl), dimension(n) :: vs, u, rs
      real(dbl), dimension(pmax,pmax) :: am, ai
      real(dbl), dimension(pmax) :: b, t
      integer :: i, j
      vs = e2 + s**2
      u = q / vs
      am = 0
      b = 0
      do i = 1, p
         b(i) = sum(u*xm(:,i)*yy)
         do j = 1, p
            am(i,j) = sum(u*xm(:,i)*xm(:,j))
         end do
      end do
      call syminv(am,p,ai)
      t = 0
      t(1:p) = matmul(ai(1:p,1:p),b(1:p))
      rs = yy - matmul(xm(:,1:p),t(1:p))
      remlike = -(sum(q*log(vs)) + logdet(am,p) + sum(u*rs**2)) / 2
    end function remlike

    subroutine effvar(th)
      ! the variance of y - g(x) by the errors of y and x: g' the slopes
      ! of the polynomial by the colours, the covariances of the colours
      ! and of y with them
      real(dbl), dimension(:), intent(in) :: th
      real(dbl), dimension(n) :: gp
      real(dbl), dimension(n,nc) :: gpm
      integer :: j, k
      if( nc == 1 ) then
         if( p >= 2 ) then
            gp = th(2)
            if( p >= 3 ) gp = gp + 2*th(3)*xx
         else
            gp = 0
         end if
         e2 = e0**2 + gp**2*ex2(:,1) - 2*gp*cc(:,1)
      else
         gpm = 0
         do k = 2, p
            if( cpow(k) == 1 ) then
               gpm(:,cind(k)) = gpm(:,cind(k)) + th(k)
            else
               gpm(:,cind(k)) = gpm(:,cind(k)) + 2*th(k)*xc(:,cind(k))
            end if
         end do
         e2 = e0**2
         do j = 1, nc
            e2 = e2 + gpm(:,j)**2*ex2(:,j) - 2*gpm(:,j)*cc(:,j)
            do k = j + 1, nc
               e2 = e2 + 2*gpm(:,j)*gpm(:,k)*cxx(:,j,k)
            end do
         end do
      end if
      ! the covariance cannot make the variance smaller than the parts
      ! independent of the catalogue magnitude (a numerical guard)
      e2 = max(e2,0.01_dbl*e0**2)
    end subroutine effvar

    subroutine estep(th,sigma,f)
      real(dbl), dimension(:), intent(in) :: th
      real(dbl), intent(in) :: sigma, f
      real(dbl), dimension(n) :: lin, lout, den
      v = e2 + sigma**2
      r = yy - matmul(xm(:,1:p),th(1:p))
      lin = exp(-r**2/(2*v)) / sqrt(v)
      lout = exp(-r**2/(2*(v + vb))) / sqrt(v + vb)
      den = (1 - f)*lin + f*lout
      where( den > 0 )
         q = (1 - f)*lin / den
      elsewhere
         q = 0
      end where
    end subroutine estep

    subroutine wls(w,th,ainv)
      ! weighted least squares: th = (X'WX)**-1 X'W y, ainv = (X'WX)**-1;
      ! the mean directly. With the errors of the colours (eiv): the observed
      ! colours scatter about the true ones by their errors, which biases
      ! the slopes towards zero (the regression dilution), and the errors
      ! of the colours correlated with y (the colour contains the catalogue
      ! magnitude of the calibration) bias them too. The moments are
      ! corrected by the expected contribution of the errors (Fuller 1987,
      ! Measurement Error Models, 2.2): X'WX - h E, X'Wy - h e, E and e
      ! (errmom) the expected parts of the errors, ainv the covariance of
      ! the corrected estimate, A**-1 X'WX A**-1. The reliability ratios of
      ! the colours (lamc, the variance of the true colours over the observed
      ! one about their weighted mean) are computed whenever the errors are
      ! given. The correction is capped (h < 1, eivst = 2) so that
      ! the corrected variance of every combination of the terms is at
      ! least lammin = 0.5 of the observed one (the correction at most
      ! a factor of 2): when the errors are comparable to the range of
      ! the colours, the full correction (1/lambda) has a heavy tail (a small
      ! difference of two noisy moments), and a nearly singular matrix
      ! (the linear and quadratic terms, two correlated colours) has
      ! the same in the worst combination of the terms.
      real(dbl), dimension(:), intent(in) :: w
      real(dbl), dimension(:), intent(out) :: th
      real(dbl), dimension(:,:), intent(out) :: ainv
      real(dbl), dimension(pmax,pmax) :: am, ac, ai, em
      real(dbl), dimension(pmax) :: b, bc, eb
      real(dbl) :: sw, h, ha, hb
      integer :: i, j, iter
      th = 0
      ainv = 0
      if( p == 1 ) then
         sw = sum(w)
         th(1) = sum(w*yy) / sw
         ainv(1,1) = 1 / sw
         return
      end if
      am = 0
      b = 0
      do i = 1, p
         b(i) = sum(w*xm(:,i)*yy)
         do j = 1, p
            am(i,j) = sum(w*xm(:,i)*xm(:,j))
         end do
      end do
      call syminv(am,p,ainv)
      th(1:p) = matmul(ainv(1:p,1:p),b(1:p))
      if( .not. (eiv .or. haserr) ) return

      call errmom(w,em,eb)
      ! the reliability ratios of the colours (their linear terms)
      do i = 2, p
         if( cpow(i) == 1 ) lamc(cind(i)) = reliab(am,em,1.0_dbl,i)
      end do
      if( .not. eiv ) return

      eivst = 0
      h = 1
      if( .not. eivok(am,em,1.0_dbl) ) then
         ! the largest h of the correction within the cap, by bisection
         ha = 0
         hb = 1
         do iter = 1, 50
            h = (ha + hb) / 2
            if( eivok(am,em,h) ) then
               ha = h
            else
               hb = h
            end if
         end do
         h = ha
         if( .not. (h > 0) ) return
      end if
      ac(1:p,1:p) = am(1:p,1:p) - h*em(1:p,1:p)
      bc(1:p) = b(1:p) - h*eb(1:p)
      call syminv(ac,p,ai)
      if( .not. (ai(1,1) > 0) ) return
      th(1:p) = matmul(ai(1:p,1:p),bc(1:p))
      ainv(1:p,1:p) = matmul(ai(1:p,1:p),matmul(am(1:p,1:p),ai(1:p,1:p)))
      eivst = merge(1,2,h == 1)
    end subroutine wls

    logical function eivok(am,em,h)
      ! the correction by h within the cap: X'WX - h E at least lammin
      ! X'WX, v'(X'WX - h E)v >= lammin v'(X'WX)v for every v, i.e.
      ! (1 - lammin) X'WX - h E positive definite (the reliability ratio
      ! of every combination of the terms at least lammin; for one linear
      ! term, lambda >= lammin)
      real(dbl), dimension(:,:), intent(in) :: am, em
      real(dbl), intent(in) :: h
      real(dbl), dimension(pmax,pmax) :: a1, a2
      a1 = 0
      a1(1:p,1:p) = (1 - lammin)*am(1:p,1:p) - h*em(1:p,1:p)
      call syminv(a1,p,a2)
      eivok = a2(1,1) > 0
    end function eivok

    real(dbl) function reliab(am,em,h,k)
      ! the reliability ratio of the term k: the weighted variance about
      ! the mean of the column k of X corrected by h E over the observed one
      real(dbl), dimension(:,:), intent(in) :: am, em
      real(dbl), intent(in) :: h
      integer, intent(in) :: k
      real(dbl) :: vo, vc
      vo = am(k,k) - am(1,k)**2/am(1,1)
      vc = (am(k,k) - h*em(k,k)) - (am(1,k) - h*em(1,k))**2/am(1,1)
      if( vo > 0 ) then
         reliab = vc / vo
      else
         reliab = 1
      end if
    end function reliab

    subroutine errmom(w,em,eb)
      ! the expected parts of the errors of the colours in the moments of
      ! wls, E[X'WX] = Xi'W Xi + E, E[X'Wy] = Xi'W eta + e (Xi, eta
      ! the true values): the linear terms by the covariances of the colours
      ! of a star and with y (as before for several colours), the terms of
      ! the quadratic polynomials (eivquad) by the unbiased estimates of
      ! the products of the true colours for the Gaussian errors (the sums
      ! over the pairings of the factors, unbiased); without the errors
      ! of the colours all zero
      real(dbl), dimension(:), intent(in) :: w
      real(dbl), dimension(:,:), intent(out) :: em
      real(dbl), dimension(:), intent(out) :: eb
      integer, dimension(5) :: fl
      integer :: i, j, ni, nj
      em = 0
      eb = 0
      if( .not. haserr ) return
      do i = 2, p
         if( cpow(i) /= 1 ) cycle
         eb(i) = sum(w*cc(:,cind(i)))
         do j = 2, p
            if( cpow(j) /= 1 ) cycle
            if( cind(i) == cind(j) ) then
               em(i,j) = sum(w*ex2(:,cind(i)))
            else
               em(i,j) = sum(w*cxx(:,cind(i),cind(j)))
            end if
         end do
      end do
      if( .not. (eivquad .and. any(cpow(2:p) == 2)) ) return
      do i = 1, p
         call factors(i,fl,ni)
         do j = i, p
            if( cpow(i) /= 2 .and. cpow(j) /= 2 ) cycle
            call factors(j,fl(ni+1:),nj)
            em(i,j) = sum(w*ecorr(fl(1:ni+nj)))
            em(j,i) = em(i,j)
         end do
         if( cpow(i) == 2 ) then
            fl(ni+1) = 0
            eb(i) = sum(w*ecorr(fl(1:ni+1)))
         end if
      end do
    end subroutine errmom

    subroutine factors(k,fl,nf)
      ! the colours of the column k of X as a list of factors
      integer, intent(in) :: k
      integer, dimension(:), intent(out) :: fl
      integer, intent(out) :: nf
      nf = 0
      if( k > 1 ) then
         nf = cpow(k)
         fl(1:nf) = cind(k)
      end if
    end subroutine factors

    function serr(a,b) result(s)
      ! the covariance of the errors of the factors a, b (0 is y)
      integer, intent(in) :: a, b
      real(dbl), dimension(n) :: s
      if( a == 0 .and. b == 0 ) then
         s = 0
      else if( a == 0 ) then
         s = cc(:,b)
      else if( b == 0 ) then
         s = cc(:,a)
      else if( a == b ) then
         s = ex2(:,a)
      else
         s = cxx(:,a,b)
      end if
    end function serr

    function fval(a) result(v)
      ! the observed value of the factor a (0 is y)
      integer, intent(in) :: a
      real(dbl), dimension(n) :: v
      if( a == 0 ) then
         v = yy
      else
         v = xc(:,a)
      end if
    end function fval

    recursive function eprod(fl) result(v)
      ! the unbiased estimate of the product of the true values of
      ! the factors fl for the Gaussian errors: the sum over the pairings,
      ! (-1)**(pairs) prod cov(pair) prod x(unpaired),
      ! eprod(f) = x(f1) eprod(rest) - sum_j cov(f1,fj) eprod(rest without fj)
      integer, dimension(:), intent(in) :: fl
      real(dbl), dimension(n) :: v
      integer :: j, m
      m = size(fl)
      if( m == 0 ) then
         v = 1
         return
      end if
      v = fval(fl(1))*eprod(fl(2:m))
      do j = 2, m
         v = v - serr(fl(1),fl(j))*eprod([fl(2:j-1),fl(j+1:m)])
      end do
    end function eprod

    recursive function ecorr(fl) result(v)
      ! the observed product of the factors fl less its unbiased estimate
      ! (eprod), directly: ecorr(f) = x(f1) ecorr(rest) + sum_j cov(f1,fj)
      ! eprod(rest without fj)
      integer, dimension(:), intent(in) :: fl
      real(dbl), dimension(n) :: v
      integer :: j, m
      m = size(fl)
      v = 0
      if( m < 2 ) return
      v = fval(fl(1))*ecorr(fl(2:m))
      do j = 2, m
         v = v + serr(fl(1),fl(j))*eprod([fl(2:j-1),fl(j+1:m)])
      end do
    end function ecorr

    real(dbl) function madnorm(s)
      real(dbl), intent(in) :: s
      madnorm = median(abs(r)/sqrt(e2 + s**2))
    end function madnorm

    real(dbl) function madfun(s)
      real(dbl), intent(in) :: s
      madfun = madnorm(s) - q50
    end function madfun

    real(dbl) function remlfun(s)
      ! the REML equation of the calibration stars (decreasing in s),
      ! h the leverage of the polynomial (of the mean: u/sum(u))
      real(dbl), intent(in) :: s
      real(dbl), dimension(n) :: u, ws, h
      real(dbl), dimension(pmax,pmax) :: am, ai
      real(dbl) :: su
      integer :: i, j
      ws = 1 / (e2 + s**2)
      u = q * ws
      su = sum(u)
      if( .not. (su > 0) ) then
         remlfun = -1
         return
      end if
      if( p == 1 ) then
         h = u / su
      else
         am = 0
         do i = 1, p
            do j = 1, p
               am(i,j) = sum(u*xm(:,i)*xm(:,j))
            end do
         end do
         call syminv(am,p,ai)
         do i = 1, n
            h(i) = u(i)*dot_product(xm(i,1:p),matmul(ai(1:p,1:p),xm(i,1:p)))
         end do
      end if
      remlfun = sum(u*ws*r**2) - sum(u*(1 - h))
    end function remlfun

  end subroutine mixfit


  real(dbl) function tquant68(nu)

    ! the quantile t of Student's distribution of nu degrees of freedom
    ! (nu >= 1, real) of the probability of the standard error of the Normal
    ! distribution, P(|T| < t) = erf(1/sqrt(2)) = 0.6827, by bisection of
    ! P(|T| < t) = 1 - I_x(nu/2,1/2), x = nu/(nu + t**2); t = 1.8373 for
    ! nu = 1, 1.3213 for 2, 1.0 for nu -> infinity

    real(dbl), intent(in) :: nu

    real(dbl), parameter :: p68 = 0.6826894921370859_dbl
    real(dbl) :: a, b, m
    integer :: i

    if( nu > 1e7_dbl ) then
       tquant68 = 1
       return
    end if
    a = 1
    b = 2
    do i = 1, 60
       m = (a + b) / 2
       if( 1 - betai(nu/2,0.5_dbl,nu/(nu + m**2)) < p68 ) then
          a = m
       else
          b = m
       end if
       if( b - a < 1e-12_dbl ) exit
    end do
    tquant68 = (a + b) / 2

  end function tquant68


  real(dbl) function betai(a,b,x)

    ! the regularised incomplete beta function I_x(a,b) by the continued
    ! fraction (modified Lentz), Numerical Recipes 6.4

    real(dbl), intent(in) :: a, b, x

    real(dbl), parameter :: fpmin = 1e-300_dbl, eps = 1e-15_dbl
    real(dbl) :: bt, aa, c, d, del, h, qab, qam, qap, xx, aa1, bb
    integer :: m, m2
    logical :: swap

    if( x <= 0 ) then
       betai = 0
       return
    else if( x >= 1 ) then
       betai = 1
       return
    end if
    swap = x >= (a + 1) / (a + b + 2)
    if( swap ) then
       xx = 1 - x
       aa1 = b
       bb = a
    else
       xx = x
       aa1 = a
       bb = b
    end if
    bt = exp(log_gamma(aa1 + bb) - log_gamma(aa1) - log_gamma(bb) + &
         aa1*log(xx) + bb*log(1 - xx))
    qab = aa1 + bb
    qap = aa1 + 1
    qam = aa1 - 1
    c = 1
    d = 1 - qab*xx/qap
    if( abs(d) < fpmin ) d = fpmin
    d = 1 / d
    h = d
    do m = 1, 10000
       m2 = 2*m
       aa = m*(bb - m)*xx / ((qam + m2)*(aa1 + m2))
       d = 1 + aa*d
       if( abs(d) < fpmin ) d = fpmin
       c = 1 + aa/c
       if( abs(c) < fpmin ) c = fpmin
       d = 1 / d
       h = h*d*c
       aa = -(aa1 + m)*(qab + m)*xx / ((aa1 + m2)*(qap + m2))
       d = 1 + aa*d
       if( abs(d) < fpmin ) d = fpmin
       c = 1 + aa/c
       if( abs(c) < fpmin ) c = fpmin
       d = 1 / d
       del = d*c
       h = h*del
       if( abs(del - 1) < eps ) exit
    end do
    betai = bt*h/aa1
    if( swap ) betai = 1 - betai

  end function betai


  real(dbl) function logdet(a,p)

    ! the logarithm of the determinant of the symmetric positive definite
    ! matrix a(1:p,1:p) by the Cholesky decomposition (p <= pmax)

    real(dbl), dimension(pmax,pmax), intent(in) :: a
    integer, intent(in) :: p

    real(dbl), dimension(pmax,pmax) :: l
    real(dbl) :: s
    integer :: i, j

    logdet = 0
    l = 0
    do j = 1, p
       s = a(j,j) - sum(l(j,1:j-1)**2)
       if( .not. (s > 0) ) then
          logdet = -huge(s)
          return
       end if
       l(j,j) = sqrt(s)
       logdet = logdet + log(s)
       do i = j+1, p
          l(i,j) = (a(i,j) - sum(l(i,1:j-1)*l(j,1:j-1))) / l(j,j)
       end do
    end do

  end function logdet


  subroutine syminv(a,p,ai)

    ! the inverse of the symmetric positive definite matrix a(1:p,1:p)
    ! by the Cholesky decomposition (p <= pmax); ai = 0 for a singular one

    real(dbl), dimension(pmax,pmax), intent(in) :: a
    integer, intent(in) :: p
    real(dbl), dimension(pmax,pmax), intent(out) :: ai

    real(dbl), dimension(pmax,pmax) :: l, li
    integer :: i, j, k
    real(dbl) :: s

    ai = 0
    l = 0
    do j = 1, p
       s = a(j,j) - sum(l(j,1:j-1)**2)
       if( .not. (s > 0) ) return
       l(j,j) = sqrt(s)
       do i = j+1, p
          l(i,j) = (a(i,j) - sum(l(i,1:j-1)*l(j,1:j-1))) / l(j,j)
       end do
    end do
    ! the inverse of the lower triangular matrix
    li = 0
    do j = 1, p
       li(j,j) = 1 / l(j,j)
       do i = j+1, p
          s = 0
          do k = j, i-1
             s = s + l(i,k)*li(k,j)
          end do
          li(i,j) = -s / l(i,i)
       end do
    end do
    ai(1:p,1:p) = matmul(transpose(li(1:p,1:p)),li(1:p,1:p))

  end subroutine syminv


  real(dbl) function theilsen(x,y)

    ! the median of the slopes of the pairs of points (Theil-Sen);
    ! for many points, the pairs of every point with a regular subset
    ! of the others (at most about 2e6 pairs)

    real(dbl), dimension(:), intent(in) :: x, y
    real(dbl), dimension(:), allocatable :: s
    real(dbl) :: dx, dmin
    integer :: n, i, j, m, step
    integer(int64) :: npairs

    n = size(x)
    npairs = int(n,int64)*(n-1)/2
    step = int(max(1_int64,npairs/2000000_int64 + 1))
    allocate(s(int(min(npairs,2000000_int64)) + n))
    dmin = 1e-6_dbl*max(maxval(x) - minval(x),tiny(dmin))
    m = 0
    do i = 1, n - 1
       do j = i + 1, n, step
          dx = x(j) - x(i)
          if( abs(dx) > dmin .and. m < size(s) ) then
             m = m + 1
             s(m) = (y(j) - y(i)) / dx
          end if
       end do
    end do
    if( m > 0 ) then
       theilsen = median(s(1:m))
    else
       theilsen = 0
    end if

  end function theilsen


  real(dbl) function median(x)

    real(dbl), dimension(:), intent(in) :: x
    real(dbl), dimension(:), allocatable :: a
    integer :: n

    n = size(x)
    if( n == 0 ) then
       median = 0
       return
    end if
    allocate(a(n))
    a = x
    call sort(a)
    if( mod(n,2) == 1 ) then
       median = a(n/2+1)
    else
       median = (a(n/2) + a(n/2+1)) / 2
    end if

  end function median


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
