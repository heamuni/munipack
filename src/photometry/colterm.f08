!
!  Colour terms of the photometric calibration
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
!  A passband of an instrument differs from the passband of the catalogue
!  (or of the standard system), so the ratio of the expected photons and
!  the counts of a star depends on its colour. In magnitudes,
!
!      m_cat - m_inst = ZP + k1 (c - cref) + k2 (c - cref)**2,
!
!  c a colour index (a difference of two magnitudes), cref a reference
!  colour. Munipack calibrates the counts to photons by the ratio
!  t = ph/ct (CTPH); the same relation is
!
!      ln(ph/ct) = b0 + b1 (c - cref) + b2 (c - cref)**2,
!
!  t = exp(b0) the calibration of a star of the reference colour,
!  b1 = -k1/magln, b2 = -k2/magln (magln = 2.5/ln 10). The photons of
!  a star are ph = t exp(g(c)) ct, g(c) = b1 (c - cref) + b2 (c - cref)**2,
!  so every quantity derived from the photons (MAG, FLUX, ...) has
!  the colour term of the magnitudes.
!
!  The routine mixture_fit fits the polynomial (of the order 0, 1, 2)
!  together with the excess scatter sigma of the stars over their errors,
!  robustly to outliers, by the likelihood of a mixture (Hogg, Bovy &
!  Lang 2010, arXiv:1008.4686) as mixture_mean of exscatter.f08 does
!  for the mean: a star is a calibration star, N(X beta, v), with
!  the probability 1 - f, or an outlier, N(X beta, v + B**2), with
!  the probability f. The variance of a star v = e**2 + g'**2 ex**2
!  - 2 g' cov(y,c) + sigma**2 includes the error of the colour ex
!  propagated by the slope g' = dg/dc (the effective variance) and its
!  covariance with the catalogue magnitude of the calibration (the colour
!  B-V has the error of V). The EM algorithm: the probabilities q of
!  the stars to be calibration stars, beta by the weighted least squares
!  with the weights q/v + (1-q)/(v + B**2), sigma by the restricted (REML)
!  likelihood of the calibration stars,
!
!      sum u w r**2 = sum u (1 - h),   u = q w, w = 1/v,
!      h(i) = u(i) x(i)' (X' U X)**(-1) x(i)  (the leverage of the fit),
!
!  f the mean of 1 - q (at most 1/2), B**2 the excess variance of
!  the outliers (at least 9 (median e**2 + sigma**2) and (0.1 mag)**2).
!  The start is robust: the slope by the median of the slopes of pairs
!  of stars (Theil-Sen), the intercept by the median, sigma by the median
!  absolute deviation. Outliers have a negligible weight in the zero point,
!  in the colour coefficients and in sigma. The covariance of beta is
!  (X' U X)**(-1).
!

module colterm

  use iso_fortran_env

  implicit none
  private

  integer, parameter :: dbl = REAL64

  ! magnitudes per a unit of the natural logarithm of a flux, 2.5 / ln(10)
  real(dbl), parameter, public :: magln = 1.0857362047581294_dbl

  ! the reference colour is the mean colour of the calibration stars
  real(dbl), parameter, public :: cref_mean = huge(1.0_dbl)

  ! models of the colour dependence:
  ! COLT_NONE no colour term,
  ! COLT_CAT  the colour term of the calibration fitted on the catalogue
  !           colour (phcal --color-term), applied by the catalogue colour,
  ! COLT_TRA  the magnitude transformation of phfotran (--model mag),
  !           its coefficients fixed, applied by the colour of the
  !           calibrated magnitudes of the frames (phcal --tratab)
  integer, parameter, public :: COLT_NONE = 0, COLT_CAT = 1, COLT_TRA = 2

  type, public :: colour_term
     integer :: model = COLT_NONE
     ! the colour of every band: the names of the catalogue columns
     ! (COLT_CAT) or of the filters (COLT_TRA), c = A - B
     character(len=68), dimension(:), allocatable :: cola, colb
     ! COLT_TRA: the indices of the filters of the colour of every band
     integer, dimension(:), allocatable :: ia, ib
     ! the order of the polynomial, the number of the stars used
     integer, dimension(:), allocatable :: order, nstars
     ! ln(ph/ct) = beta(0) + beta(1) d + beta(2) d**2, d = c - cref
     ! (beta(0,:) is ln CTPH), the covariance of beta (in the order 0,1,2)
     real(dbl), dimension(:,:), allocatable :: beta
     real(dbl), dimension(:,:,:), allocatable :: cov
     ! the reference colour; the mean colour of the calibration stars,
     ! its scatter, the mean colour correction g of the calibration stars
     ! and its scatter (used for the stars without a colour)
     real(dbl), dimension(:), allocatable :: cref, cmean, cstd, gmean, gstd
     ! COLT_CAT: the catalogue (all its stars): positions [deg], colours,
     ! their errors (< 0 for an undefined colour), the search radius [deg]
     real(dbl), dimension(:), allocatable :: ra, dec, col, dcol
     real(dbl) :: tol = 0
  end type colour_term

  public :: colterm_init, mixture_fit, colour_corr

contains

  subroutine colterm_init(ctm,model,nbands)

    type(colour_term), intent(out) :: ctm
    integer, intent(in) :: model, nbands

    ctm%model = model
    allocate(ctm%cola(nbands),ctm%colb(nbands),ctm%ia(nbands),ctm%ib(nbands),&
         ctm%order(nbands),ctm%nstars(nbands),ctm%beta(0:2,nbands), &
         ctm%cov(3,3,nbands),ctm%cref(nbands),ctm%cmean(nbands), &
         ctm%cstd(nbands),ctm%gmean(nbands),ctm%gstd(nbands))
    ctm%cola = ''
    ctm%colb = ''
    ctm%ia = 0
    ctm%ib = 0
    ctm%order = 1
    ctm%nstars = 0
    ctm%beta = 0
    ctm%cov = 0
    ctm%cref = 0
    ctm%cmean = 0
    ctm%cstd = 0
    ctm%gmean = 0
    ctm%gstd = 0
    allocate(ctm%ra(0),ctm%dec(0),ctm%col(0),ctm%dcol(0))

  end subroutine colterm_init


  subroutine colour_corr(ctm,n,c,dc,g,dg)

    ! The colour correction g = ln(ph/ct) - ln CTPH of a star of the colour
    ! c (with the error dc, dc < 0: an undefined colour) in the band n,
    ! and its error dg: the error of the calibration at the colour (by
    ! the covariance of beta, the error of CTPH included) and the error of
    ! the colour. A star without a colour gets the mean correction of
    ! the calibration stars, its error inflated by their scatter.

    type(colour_term), intent(in) :: ctm
    integer, intent(in) :: n
    real(dbl), intent(in) :: c, dc
    real(dbl), intent(out) :: g, dg

    real(dbl), dimension(3) :: x
    real(dbl) :: d, gp, v

    if( dc >= 0 ) then
       d = c - ctm%cref(n)
       g = ctm%beta(1,n)*d + ctm%beta(2,n)*d**2
       gp = ctm%beta(1,n) + 2*ctm%beta(2,n)*d
       x = [ 1.0_dbl, d, d**2 ]
       v = dot_product(x,matmul(ctm%cov(:,:,n),x)) + (gp*dc)**2
    else
       g = ctm%gmean(n)
       v = ctm%cov(1,1,n) + ctm%gstd(n)**2
    end if
    dg = sqrt(max(v,0.0_dbl))

  end subroutine colour_corr


  subroutine mixture_fit(y,e,x,ex,cxy,order,fixsig,beta,cov,sigma,fout,member)

    ! The polynomial ln(ph/ct) = beta(0) + beta(1) x + beta(2) x**2
    ! (beta(order+1:) = 0) with the excess scatter sigma by the likelihood
    ! of a mixture (see above). y = ln(ph/ct), e its error, x the colour
    ! minus a reference (centred for the numerics), ex the error of
    ! the colour, cxy the covariance of y and x, fixsig < 0: sigma by
    ! the restricted likelihood, >= 0: sigma fixed at the value; cov
    ! the covariance of beta, fout the fraction of outliers, member
    ! (optional) the probabilities of the stars to be calibration stars
    ! (0 for e <= 0). A value is used when e > 0, ex >= 0 and y, x finite.

    real(dbl), dimension(:), intent(in) :: y, e, x, ex, cxy
    integer, intent(in) :: order
    real(dbl), intent(in) :: fixsig
    real(dbl), dimension(0:2), intent(out) :: beta
    real(dbl), dimension(3,3), intent(out) :: cov
    real(dbl), intent(out) :: sigma, fout
    real(dbl), dimension(:), intent(out), optional :: member

    real(dbl), parameter :: q50 = 0.6744897501960817_dbl
    real(dbl), parameter :: fmax = 0.5_dbl, tol = 1e-10_dbl
    ! 1 mag in the natural logarithm of flux
    real(dbl), parameter :: bmag = 1 / magln
    real(dbl), dimension(:), allocatable :: yy, e0, xx, ex2, cc, e2, r, q, v, w, a
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
       if( present(member) ) member = merge(1.0_dbl,0.0_dbl,mask)
       return
    end if

    vb = bmag**2
    vbmin = (0.1_dbl*bmag)**2

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
      ! weighted least squares: th = (X'WX)**-1 X'W y, ainv = (X'WX)**-1
      real(dbl), dimension(:), intent(in) :: w
      real(dbl), dimension(3), intent(out) :: th
      real(dbl), dimension(3,3), intent(out) :: ainv
      real(dbl), dimension(3,3) :: am
      real(dbl), dimension(3) :: b
      integer :: i, j
      am = 0
      b = 0
      do i = 1, p
         b(i) = sum(w*xm(:,i)*yy)
         do j = 1, p
            am(i,j) = sum(w*xm(:,i)*xm(:,j))
         end do
      end do
      call syminv(am,p,ainv)
      th = 0
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
      ! h the leverage of the polynomial
      real(dbl), intent(in) :: s
      real(dbl), dimension(n) :: u, ws, h
      real(dbl), dimension(3,3) :: am, ai
      integer :: i, j
      ws = 1 / (e2 + s**2)
      u = q * ws
      if( .not. (sum(u) > 0) ) then
         remlfun = -1
         return
      end if
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
      remlfun = sum(u*ws*r**2) - sum(u*(1 - h))
    end function remlfun

  end subroutine mixture_fit


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

end module colterm
