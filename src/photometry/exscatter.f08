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

module exscatter

  use iso_fortran_env

  implicit none
  private

  integer, parameter, private :: dbl = REAL64

  public :: excess_scatter, mixture_mean

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

    real(dbl), parameter :: q50 = 0.6744897501960817_dbl
    real(dbl), parameter :: fmax = 0.5_dbl, tol = 1e-10_dbl
    ! 1 mag in the natural logarithm of flux
    real(dbl), parameter :: bmag = 0.4_dbl*2.302585092994046_dbl
    real(dbl), dimension(:), allocatable :: x, e2, r, q, v, w, a
    logical, dimension(size(y)) :: mask
    real(dbl) :: vb, vbmin, smax, f, mu1, sig1, f1, emin, e2med, sq
    integer :: n, iter, m

    mu = 0
    dmu = 0
    sigma = 0
    fout = 0
    if( present(member) ) member = 0

    mask = err > 0 .and. abs(y) < huge(y)
    n = count(mask)
    if( n == 0 ) return

    x = pack(y,mask)
    e2 = pack(err,mask)**2
    allocate(r(n),q(n),v(n),w(n))
    q = 1

    if( n < 3 ) then
       ! too few values for the scatter: the weighted mean
       w = 1 / e2
       mu = sum(w*x) / sum(w)
       dmu = 1 / sqrt(sum(w))
       if( present(member) ) member = merge(1.0_dbl,0.0_dbl,mask)
       return
    end if

    vb = bmag**2
    vbmin = (0.1_dbl*bmag)**2
    if( present(bmin) ) vbmin = bmin**2
    emin = sqrt(minval(e2))
    a = e2
    call sort(a)
    if( mod(n,2) == 1 ) then
       e2med = a(n/2+1)
    else
       e2med = (a(n/2) + a(n/2+1)) / 2
    end if

    ! the robust start: the median and the scatter by the median
    ! of |r|/sqrt(err**2 + sigma**2), equal to the one of N(0,1)
    a = x
    call sort(a)
    m = size(a)
    if( mod(m,2) == 1 ) then
       mu = a(m/2+1)
    else
       mu = (a(m/2) + a(m/2+1)) / 2
    end if
    r = x - mu
    smax = 10*maxval(abs(r)) + tiny(smax)
    if( madnorm(0.0_dbl) > q50 ) sigma = bisect(madfun,0.0_dbl,smax)
    f = 0.1_dbl

    do iter = 1, 1000

       ! E-step: the probabilities of inliers
       call estep(mu,sigma,f)

       ! M-step: mu, sigma (REML), f
       v = e2 + sigma**2
       w = q/v + (1 - q)/(v + vb)
       mu1 = sum(w*x) / sum(w)
       r = x - mu1
       smax = 10*maxval(abs(r)) + tiny(smax)
       if( remlfun(0.0_dbl) > 0 ) then
          sig1 = bisect(remlfun,0.0_dbl,smax)
       else
          sig1 = 0
       end if
       f1 = min(sum(1 - q)/n, fmax)
       sq = sum(1 - q)
       if( sq > 0.5_dbl ) &
            vb = max(sum((1 - q)*(r**2 - e2 - sig1**2))/sq, &
            9*(e2med + sig1**2), vbmin)

       if( abs(mu1 - mu) <= tol*emin .and. abs(sig1 - sigma) <= tol*emin &
            .and. abs(f1 - f) <= tol ) exit
       mu = mu1
       sigma = sig1
       f = f1
    end do
    mu = mu1
    sigma = sig1
    fout = f1

    call estep(mu,sigma,fout)
    w = q / (e2 + sigma**2)
    if( sum(w) > 0 ) then
       dmu = 1 / sqrt(sum(w))
    else
       dmu = 1 / sqrt(sum(1/(e2 + sigma**2)))
    end if

    if( present(member) ) member = unpack(q,mask,0.0_dbl)

  contains

    subroutine estep(mu,sigma,f)
      real(dbl), intent(in) :: mu, sigma, f
      real(dbl), dimension(n) :: lin, lout, den
      v = e2 + sigma**2
      r = x - mu
      lin = exp(-r**2/(2*v)) / sqrt(v)
      lout = exp(-r**2/(2*(v + vb))) / sqrt(v + vb)
      den = (1 - f)*lin + f*lout
      where( den > 0 )
         q = (1 - f)*lin / den
      elsewhere
         q = 0
      end where
    end subroutine estep

    real(dbl) function madnorm(s)
      real(dbl), intent(in) :: s
      real(dbl), dimension(n) :: b
      b = abs(r)/sqrt(e2 + s**2)
      call sort(b)
      if( mod(n,2) == 1 ) then
         madnorm = b(n/2+1)
      else
         madnorm = (b(n/2) + b(n/2+1)) / 2
      end if
    end function madnorm

    real(dbl) function madfun(s)
      real(dbl), intent(in) :: s
      madfun = madnorm(s) - q50
    end function madfun

    real(dbl) function remlfun(s)
      ! the REML equation of the inliers (decreasing in s)
      real(dbl), intent(in) :: s
      real(dbl), dimension(n) :: u, ws
      real(dbl) :: su
      ws = 1 / (e2 + s**2)
      u = q * ws
      su = sum(u)
      if( su > 0 ) then
         remlfun = sum(u*ws*r**2) - sum(u*(1 - u/su))
      else
         remlfun = -1
      end if
    end function remlfun

  end subroutine mixture_mean


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
