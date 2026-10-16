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

module exscatter

  use iso_fortran_env

  implicit none
  private

  integer, parameter, private :: dbl = REAL64

  public :: excess_scatter, mixture_mean, mixture_fit

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


  subroutine mixture_mean(y,err,mu,dmu,sigma,fout,member,bmin)

    ! The mean of y with the errors err and the excess scatter sigma,
    ! robust to outliers (see above); dmu is the error of mu, fout the
    ! fraction of outliers, member (optional) the probabilities of
    ! the values to be inliers (0 for err <= 0), bmin (optional) the
    ! smallest width B of the distribution of outliers.

    real(dbl), dimension(:), intent(in) :: y, err
    real(dbl), intent(out) :: mu, dmu, sigma, fout
    real(dbl), dimension(:), intent(out), optional :: member
    real(dbl), intent(in), optional :: bmin

    real(dbl), dimension(size(y)) :: x0
    real(dbl), dimension(0:2) :: beta
    real(dbl), dimension(3,3) :: cov

    x0 = 0
    call mixfit(y,err,x0,x0,x0,0,-1.0_dbl,beta,cov,sigma,fout,member,bmin,dmu)
    mu = beta(0)

  end subroutine mixture_mean


  subroutine mixture_fit(y,e,x,ex,cxy,order,fixsig,beta,cov,sigma,fout,member,bmin)

    ! The polynomial y = beta(0) + beta(1) x + beta(2) x**2
    ! (beta(order+1:) = 0) with the excess scatter sigma by the likelihood
    ! of a mixture (see above). y = ln(ph/ct), e its error, x the colour
    ! minus a reference (centred for the numerics), ex the error of
    ! the colour, cxy the covariance of y and x, fixsig < 0: sigma by
    ! the restricted likelihood, >= 0: sigma fixed at the value; cov
    ! the covariance of beta, fout the fraction of outliers, member
    ! (optional) the probabilities of the stars to be calibration stars
    ! (0 for e <= 0), bmin (optional) the smallest width B of
    ! the distribution of outliers. A value is used when e > 0, ex >= 0
    ! and y, x finite.

    real(dbl), dimension(:), intent(in) :: y, e, x, ex, cxy
    integer, intent(in) :: order
    real(dbl), intent(in) :: fixsig
    real(dbl), dimension(0:2), intent(out) :: beta
    real(dbl), dimension(3,3), intent(out) :: cov
    real(dbl), intent(out) :: sigma, fout
    real(dbl), dimension(:), intent(out), optional :: member
    real(dbl), intent(in), optional :: bmin

    call mixfit(y,e,x,ex,cxy,order,fixsig,beta,cov,sigma,fout,member,bmin)

  end subroutine mixture_fit


  subroutine mixfit(y,e,x,ex,cxy,order,fixsig,beta,cov,sigma,fout,member,bmin, &
       dbeta0)

    ! The implementation of mixture_fit and mixture_mean; dbeta0
    ! (optional) is the error of beta(0). For the order 0 (the mean),
    ! the weighted mean and its error are computed directly (not by
    ! the inverse of a 1x1 matrix).

    real(dbl), dimension(:), intent(in) :: y, e, x, ex, cxy
    integer, intent(in) :: order
    real(dbl), intent(in) :: fixsig
    real(dbl), dimension(0:2), intent(out) :: beta
    real(dbl), dimension(3,3), intent(out) :: cov
    real(dbl), intent(out) :: sigma, fout
    real(dbl), dimension(:), intent(out), optional :: member
    real(dbl), intent(in), optional :: bmin
    real(dbl), intent(out), optional :: dbeta0

    real(dbl), parameter :: q50 = 0.6744897501960817_dbl
    real(dbl), parameter :: fmax = 0.5_dbl, tol = 1e-10_dbl
    ! 1 mag in the natural logarithm of flux
    real(dbl), parameter :: bmag = 0.4_dbl*2.302585092994046_dbl
    real(dbl), dimension(:), allocatable :: yy, e0, xx, ex2, cc, e2, r, q, v, w
    real(dbl), dimension(:,:), allocatable :: xm
    logical, dimension(size(y)) :: mask
    real(dbl), dimension(3) :: th, th1
    real(dbl), dimension(3,3) :: ainv
    real(dbl) :: vb, vbmin, smax, f, sig1, f1, emin, e2med, sq, xs
    integer :: n, p, iter, k
    logical :: conv

    beta = 0
    cov = 0
    sigma = max(fixsig,0.0_dbl)
    fout = 0
    if( present(member) ) member = 0
    if( present(dbeta0) ) dbeta0 = 0

    mask = e > 0 .and. ex >= 0 .and. abs(y) < huge(y) .and. abs(x) < huge(x)
    n = count(mask)
    if( n == 0 ) return

    yy = pack(y,mask)
    e0 = pack(e,mask)
    xx = pack(x,mask)
    ex2 = pack(ex,mask)**2
    cc = pack(cxy,mask)
    allocate(e2(n),r(n),q(n),v(n),w(n))

    ! the order is limited by the stars and by the range of the colours
    p = min(order,2) + 1
    p = max(min(p,n),1)
    xs = maxval(abs(xx))
    if( .not. (maxval(xx) - minval(xx) > 1e-6_dbl*max(xs,1.0_dbl)) ) p = 1
    if( p == 1 ) xs = 1
    xs = max(xs,epsilon(xs))

    allocate(xm(n,3))
    xm(:,1) = 1
    xm(:,2) = xx
    xm(:,3) = xx**2

    q = 1
    th = 0

    if( n < p + 2 ) then
       ! too few stars for the scatter and outliers: the weighted fit
       call effvar(th)
       v = e2 + sigma**2
       w = 1 / v
       call wls(w,th,ainv)
       call effvar(th)
       v = e2 + sigma**2
       w = 1 / v
       call wls(w,th,ainv)
       beta(0:p-1) = th(1:p)
       cov(1:p,1:p) = ainv(1:p,1:p)
       if( present(dbeta0) ) dbeta0 = 1 / sqrt(sum(w))
       if( present(member) ) member = merge(1.0_dbl,0.0_dbl,mask)
       return
    end if

    vb = bmag**2
    vbmin = (0.1_dbl*bmag)**2
    if( present(bmin) ) vbmin = bmin**2

    ! the robust start: the slope by the median of the slopes of pairs,
    ! the intercept by the median, the scatter by the median absolute
    ! deviation of the normalised residuals
    th = 0
    if( p >= 2 ) th(2) = theilsen(xx,yy)
    th(1) = median(yy - th(2)*xx)
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
       f1 = min(sum(1 - q)/n, fmax)
       sq = sum(1 - q)
       if( sq > 0.5_dbl ) &
            vb = max(sum((1 - q)*(r**2 - e2 - sig1**2))/sq, &
            9*(e2med + sig1**2), vbmin)

       conv = abs(sig1 - sigma) <= tol*emin .and. abs(f1 - f) <= tol
       do k = 1, p
          conv = conv .and. abs(th1(k) - th(k))*xs**(k-1) <= tol*emin
       end do
       th = th1
       sigma = sig1
       f = f1
       if( conv ) exit
    end do
    fout = f

    call effvar(th)
    call estep(th,sigma,fout)
    w = q / (e2 + sigma**2)
    if( .not. (sum(w) > 0) ) w = 1 / (e2 + sigma**2)
    call wls(w,th1,ainv)
    beta(0:p-1) = th(1:p)
    cov(1:p,1:p) = ainv(1:p,1:p)
    if( present(dbeta0) ) dbeta0 = 1 / sqrt(sum(w))

    if( present(member) ) member = unpack(q,mask,0.0_dbl)

  contains

    subroutine effvar(th)
      ! the variance of y - g(x) by the errors of y and x
      real(dbl), dimension(3), intent(in) :: th
      real(dbl), dimension(n) :: gp
      if( p >= 2 ) then
         gp = th(2)
         if( p >= 3 ) gp = gp + 2*th(3)*xx
      else
         gp = 0
      end if
      e2 = e0**2 + gp**2*ex2 - 2*gp*cc
      ! the covariance cannot make the variance smaller than the parts
      ! independent of the catalogue magnitude (a numerical guard)
      e2 = max(e2,0.01_dbl*e0**2)
    end subroutine effvar

    subroutine estep(th,sigma,f)
      real(dbl), dimension(3), intent(in) :: th
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
      ! the mean directly
      real(dbl), dimension(:), intent(in) :: w
      real(dbl), dimension(3), intent(out) :: th
      real(dbl), dimension(3,3), intent(out) :: ainv
      real(dbl), dimension(3,3) :: am
      real(dbl), dimension(3) :: b
      real(dbl) :: sw
      integer :: i, j
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
    end subroutine wls

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
      real(dbl), dimension(3,3) :: am, ai
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


  subroutine syminv(a,p,ai)

    ! the inverse of the symmetric positive definite matrix a(1:p,1:p)
    ! by the Cholesky decomposition (p <= 3); ai = 0 for a singular one

    real(dbl), dimension(3,3), intent(in) :: a
    integer, intent(in) :: p
    real(dbl), dimension(3,3), intent(out) :: ai

    real(dbl), dimension(3,3) :: l, li
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
