!
!  psffit - fitting of stars by a PSF
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
!  Every star is fitted separately by the weighted least squares
!  (Levenberg-Marquardt) for its flux, position (and the local sky).
!  The weights are the inverse of variances predicted by the model
!  (sky + all stars), not by data, so faint stars are not biased:
!
!     var = max(model - bias, 0) / gain + readnoise**2   [ADU**2]
!
!  (plus the quantisation noise of integer frames, 1/12 ADU**2). The
!  Poisson term is not floored: a sky of a fraction of an electron per
!  pixel has the variance of that fraction (a floor of one electron made
!  the errors of faint stars 1.2-1.3 times too large on a sky of 0.4 e-).
!
!  The fitted neighbours are subtracted: the sum of models of all
!  stars is kept in an image (mstar), the model of the star being
!  fitted is removed from it before its fit and added after it.
!

module psffit

  use iso_fortran_env
  use psfmodel

  implicit none
  private

  ! flags of the fit, combined by sum
  integer, parameter, public :: &
       FLAG_SATURATED = 1, & ! saturated pixel(s) in the core
       FLAG_EXCLUDED  = 2, & ! pixel(s) in the fitting radius excluded (bad, edge)
       FLAG_POSFIXED  = 4, & ! position not fitted (faint, or moved too far)
       FLAG_NOCONV    = 8, & ! not converged
       FLAG_FAILED    = 16,& ! fit failed, the flux is not valid
       FLAG_CROWDED   = 32,& ! a neighbour within the fitting radius
       FLAG_ADDED     = 64,& ! a saturated star added (not detected by find)
       FLAG_RESID     = 128  ! a star found in the residuals (not by find)

  type, public :: PsfType
     integer :: ptype = PSF_GAUSS
     integer :: zoom = PSF_DEFZOOM
     ! analytic PSF: sigma (GAUSS) or alpha (MOFFAT), beta, axis ratio,
     ! position angle [deg]
     real(REAL64) :: width = 1, beta = 2.5_REAL64, q = 1, pa = 0
     integer :: nsub = 0
     real(REAL64), dimension(:,:), allocatable :: sub
     ! the field-dependent PSF
     integer :: nbox = 0
     real(REAL64), dimension(:,:,:), allocatable :: stamps
     real(REAL64), dimension(:), allocatable :: rcoo, acoo
     real(REAL64) :: xcen = 0, ycen = 0, scale = 1
     ! the empirical PSF: the base (analytic or file PSF above) plus
     ! the residual tables res(-ntab:ntab,-ntab:ntab,npoly) sampled with
     ! the same zoom, the terms of the polynomial in the normalised field
     ! coordinates (x - pxc)/pxs, (y - pyc)/pys: 1, x, y, x**2, x*y, y**2
     logical :: empirical = .false.
     integer :: ntab = 0, npoly = 1
     real(REAL64), dimension(:,:,:), allocatable :: res
     real(REAL64) :: pxc = 0, pyc = 0, pxs = 1, pys = 1
   contains
     procedure :: sampled, sampled_base, prepare, fwhm, poly_terms
  end type PsfType

  type, public :: FrameType
     real(REAL64), dimension(:,:), allocatable :: data, mstar
     logical, dimension(:,:), allocatable :: valid, satur
     ! gain [e-/ADU], read noise squared [ADU**2], bias [ADU],
     ! quantisation noise [ADU**2]
     real(REAL64) :: gain = 1, rn2 = 0, bias = 0, qnoise = 0
   contains
     procedure :: variance
  end type FrameType

  type, public :: StarFit
     real(REAL64) :: x, y, flux, sky
     real(REAL64) :: xerr = -1, yerr = -1, fluxerr = -1, skyerr = -1
     real(REAL64) :: chi2 = -1, sharp = 0, round = 0
     real(REAL64) :: fitrad = 0, satrad = 0
     integer :: npix = 0, flag = 0, modrad = 0
     ! the group of the simultaneous fit (0 when fitted alone), its size
     integer :: group = 0, gsize = 1
     logical :: modelled = .false.
  end type StarFit

  public :: fit_star, fit_group, add_model, fit_shape, linear_flux

contains

  subroutine prepare(this,nmod)

    ! sampling of an analytic PSF up to the radius nmod pixels

    class(PsfType), intent(in out) :: this
    integer, intent(in) :: nmod

    if( this%ptype == PSF_FILE ) return
    if( allocated(this%sub) ) deallocate(this%sub)
    this%nsub = this%zoom*(nmod + 2)
    allocate(this%sub(-this%nsub:this%nsub,-this%nsub:this%nsub))
    call psf_analytic_sample(this%ptype,this%width,this%beta,this%q,this%pa, &
         this%zoom,this%nsub,this%sub)

  end subroutine prepare


  subroutine sampled(this,x,y,sub,nsub)

    ! the sampled PSF for a star at x,y (the empirical PSF: the base
    ! plus the residual tables)

    class(PsfType), intent(in) :: this
    real(REAL64), intent(in) :: x,y
    real(REAL64), dimension(:,:), allocatable, intent(in out) :: sub
    integer, intent(out) :: nsub

    real(REAL64), dimension(:), allocatable :: pt
    integer :: t, n

    call this%sampled_base(x,y,sub,nsub)
    if( this%empirical ) then
       n = this%ntab
       allocate(pt(this%npoly))
       call this%poly_terms(x,y,pt)
       do t = 1, this%npoly
          sub(-n:n,-n:n) = sub(-n:n,-n:n) + pt(t)*this%res(:,:,t)
       end do
    end if

  end subroutine sampled


  subroutine sampled_base(this,x,y,sub,nsub)

    ! the sampled PSF for a star at x,y without the residual tables
    ! of the empirical PSF (the size covers them)

    class(PsfType), intent(in) :: this
    real(REAL64), intent(in) :: x,y
    real(REAL64), dimension(:,:), allocatable, intent(in out) :: sub
    integer, intent(out) :: nsub

    real(REAL64), dimension(:,:), allocatable :: tmp
    integer :: n

    if( this%ptype == PSF_FILE ) then
       nsub = this%nbox
       if( this%empirical ) nsub = max(nsub,this%ntab)
       if( .not. allocated(sub) ) allocate(sub(-nsub:nsub,-nsub:nsub))
       if( nsub == this%nbox ) then
          call psf_field_interpolate(this%stamps,this%rcoo,this%acoo, &
               x - this%xcen,y - this%ycen,this%scale,sub)
       else
          n = this%nbox
          allocate(tmp(-n:n,-n:n))
          call psf_field_interpolate(this%stamps,this%rcoo,this%acoo, &
               x - this%xcen,y - this%ycen,this%scale,tmp)
          sub = 0
          sub(-n:n,-n:n) = tmp
       end if
    else
       nsub = this%nsub
       if( .not. allocated(sub) ) allocate(sub(-nsub:nsub,-nsub:nsub))
       sub = this%sub
    end if

  end subroutine sampled_base


  pure subroutine poly_terms(this,x,y,pt)

    ! the terms of the polynomial of the empirical PSF at x,y

    class(PsfType), intent(in) :: this
    real(REAL64), intent(in) :: x,y
    real(REAL64), dimension(:), intent(out) :: pt

    real(REAL64) :: u, v

    u = (x - this%pxc) / this%pxs
    v = (y - this%pyc) / this%pys
    pt(1) = 1
    if( size(pt) >= 3 ) then
       pt(2) = u
       pt(3) = v
    end if
    if( size(pt) >= 6 ) then
       pt(4) = u**2
       pt(5) = u*v
       pt(6) = v**2
    end if

  end subroutine poly_terms


  real(REAL64) function fwhm(this,x,y)

    ! FWHM at the position x,y [pix]

    class(PsfType), intent(in) :: this
    real(REAL64), intent(in) :: x,y
    real(REAL64), dimension(:,:), allocatable :: sub
    integer :: nsub

    if( this%ptype == PSF_FILE .or. this%empirical ) then
       call this%sampled(x,y,sub,nsub)
       fwhm = psf_sampled_fwhm(sub,this%zoom)
    else
       fwhm = psf_analytic_fwhm(this%ptype,this%width,this%beta,this%q)
    end if

  end function fwhm


  elemental real(REAL64) function variance(this,m)

    ! the variance [ADU**2] of a pixel with the mean m [ADU]

    class(FrameType), intent(in) :: this
    real(REAL64), intent(in) :: m

    ! the Poisson noise of the mean (not of one electron at least) and
    ! the read and the quantisation noise; a guard of (0.01 e-)**2 keeps
    ! the weight finite for a noiseless model (no sky, no read noise)
    variance = max(m - this%bias,0.0_REAL64) / this%gain + this%rn2 + this%qnoise
    variance = max(variance,1e-4_REAL64 / this%gain**2)

  end function variance


  subroutine add_model(frame,psf,star,sign)

    ! adds (sign=1) or subtracts (sign=-1) the model of the star to mstar

    type(FrameType), intent(in out) :: frame
    type(PsfType), intent(in) :: psf
    type(StarFit), intent(in) :: star
    integer, intent(in) :: sign

    real(REAL64), dimension(:,:), allocatable :: sub, p
    integer :: nsub, n, i0, j0, i1, i2, j1, j2

    if( star%flag >= FLAG_FAILED .and. mod(star%flag/FLAG_FAILED,2) == 1 ) &
         return
    n = star%modrad
    call psf%sampled(star%x,star%y,sub,nsub)
    allocate(p(-n:n,-n:n))
    i0 = nint(star%x)
    j0 = nint(star%y)
    call psf_pixelize(sub,nsub,psf%zoom,i0 - star%x,j0 - star%y,n,p)
    i1 = max(i0 - n,1)
    i2 = min(i0 + n,size(frame%mstar,1))
    j1 = max(j0 - n,1)
    j2 = min(j0 + n,size(frame%mstar,2))
    if( i1 > i2 .or. j1 > j2 ) return
    frame%mstar(i1:i2,j1:j2) = frame%mstar(i1:i2,j1:j2) + &
         sign*star%flux*p(i1-i0:i2-i0,j1-j0:j2-j0)

  end subroutine add_model


  subroutine linear_flux(frame,psf,star,localsky)

    ! The flux (and the sky) at the fixed position by the linear weighted
    ! least squares, used for the initial estimates. The model of the
    ! star must not be included in mstar.

    type(FrameType), intent(in) :: frame
    type(PsfType), intent(in) :: psf
    type(StarFit), intent(in out) :: star
    logical, intent(in) :: localsky

    real(REAL64), dimension(:), allocatable :: d, b, pk, w
    real(REAL64), dimension(:,:), allocatable :: sub, p
    integer, dimension(:), allocatable :: ii, jj
    real(REAL64) :: a11, a12, a22, g1, g2, det, f, s
    integer :: nsub, n, i0, j0, k, npix, iter

    call pixel_set(frame,star,ii,jj,npix,n)
    star%npix = npix
    if( npix < 4 ) then
       star%flux = 0
       star%fluxerr = -1
       return
    end if

    call psf%sampled(star%x,star%y,sub,nsub)
    allocate(p(-n:n,-n:n),d(npix),b(npix),pk(npix),w(npix))
    i0 = nint(star%x)
    j0 = nint(star%y)
    call psf_pixelize(sub,nsub,psf%zoom,i0 - star%x,j0 - star%y,n,p)
    do k = 1, npix
       d(k) = frame%data(ii(k),jj(k)) - frame%mstar(ii(k),jj(k))
       b(k) = frame%mstar(ii(k),jj(k))
       pk(k) = p(ii(k)-i0,jj(k)-j0)
    end do

    f = max(star%flux,0.0_REAL64)
    s = star%sky
    do iter = 1, 3
       w = 1 / frame%variance(b + f*pk + s)
       if( localsky ) then
          a11 = sum(w*pk**2)
          a12 = sum(w*pk)
          a22 = sum(w)
          g1 = sum(w*pk*d)
          g2 = sum(w*d)
          det = a11*a22 - a12**2
          if( .not. (det > 0) ) exit
          f = (a22*g1 - a12*g2) / det
          s = (a11*g2 - a12*g1) / det
          star%fluxerr = sqrt(a22 / det)
       else
          a11 = sum(w*pk**2)
          if( .not. (a11 > 0) ) exit
          f = sum(w*pk*(d - s)) / a11
          star%fluxerr = sqrt(1 / a11)
       end if
    end do
    star%flux = f
    star%sky = s

  end subroutine linear_flux


  subroutine pixel_set(frame,star,ii,jj,npix,n)

    ! valid pixels within the fitting radius, n is the box half-size

    type(FrameType), intent(in) :: frame
    type(StarFit), intent(in) :: star
    integer, dimension(:), allocatable, intent(out) :: ii, jj
    integer, intent(out) :: npix, n

    integer :: i, j, i0, j0, m
    real(REAL64) :: r2

    n = ceiling(star%fitrad) + 2
    i0 = nint(star%x)
    j0 = nint(star%y)
    m = (2*n + 1)**2
    allocate(ii(m),jj(m))
    r2 = star%fitrad**2
    npix = 0
    do j = j0 - n, j0 + n
       do i = i0 - n, i0 + n
          if( (i - star%x)**2 + (j - star%y)**2 > r2 ) cycle
          if( i < 1 .or. j < 1 .or. &
               i > size(frame%data,1) .or. j > size(frame%data,2) ) cycle
          if( .not. frame%valid(i,j) ) cycle
          npix = npix + 1
          ii(npix) = i
          jj(npix) = j
       end do
    end do

  end subroutine pixel_set


  subroutine fit_star(frame,psf,star,localsky,maxiter,snrpos,rcore)

    ! The fit of the star by the PSF for the flux, the position (and
    ! the sky). The model of the star must not be included in mstar,
    ! the star contains the initial estimates.

    type(FrameType), intent(in) :: frame
    type(PsfType), intent(in) :: psf
    type(StarFit), intent(in out) :: star
    logical, intent(in) :: localsky
    integer, intent(in) :: maxiter
    real(REAL64), intent(in) :: snrpos, rcore

    real(REAL64), dimension(:), allocatable :: d, b, w, m, mt
    real(REAL64), dimension(:,:), allocatable :: sub, jac, jact
    integer, dimension(:), allocatable :: ii, jj
    real(REAL64), dimension(4) :: par, ptry, g, dp, sb
    real(REAL64), dimension(4,4) :: a, ad, cov
    logical, dimension(4) :: free
    real(REAL64) :: chi2, chi2try, lambda, x0, y0, step, flim, sum1, sum2
    integer :: nsub, npix, n, k, iter, ntry, np, i, j
    logical :: converged, ok, posfixed

    star%flag = 0
    x0 = star%x
    y0 = star%y

    ! flags of the region: saturated core, excluded pixels
    call region_flags(frame,star,rcore)

    ! the maximal distance from the initial position
    flim = max(star%fitrad / 2,1.0_REAL64)

    call pixel_set(frame,star,ii,jj,npix,n)
    n = ceiling(star%fitrad + flim) + 2
    star%npix = npix
    np = 3
    if( localsky ) np = 4
    if( npix < np + 3 ) then
       star%flag = star%flag + FLAG_FAILED
       star%fluxerr = -1
       return
    end if

    call psf%sampled(star%x,star%y,sub,nsub)
    allocate(d(npix),b(npix),w(npix),m(npix),mt(npix),jac(npix,4),jact(npix,4))
    do k = 1, npix
       d(k) = frame%data(ii(k),jj(k)) - frame%mstar(ii(k),jj(k))
       b(k) = frame%mstar(ii(k),jj(k))
    end do

    ! the flux at the initial position decides on the position fit
    call linear_flux(frame,psf,star,localsky)
    posfixed = .not. (star%fluxerr > 0 .and. star%flux > snrpos*star%fluxerr)

    do   ! (once more with the fixed position when the fit runs away)

       free = [.true., .not. posfixed, .not. posfixed, localsky]
       par = [star%flux, star%x, star%y, star%sky]
       lambda = 1e-3_REAL64
       converged = .false.
       ok = .true.

       do iter = 1, maxiter

          call evaluate(par,m,jac)
          w = 1 / frame%variance(b + m + par(4))
          call normal(m,w,chi2,a,g)

          do ntry = 1, 12
             ad = a
             do i = 1, 4
                ad(i,i) = a(i,i)*(1 + lambda)
             end do
             call solve(ad,g,dp,ok)
             if( .not. ok ) exit
             ! limit the shift of the position
             step = maxval(abs(dp(2:3)))
             if( step > 0.5_REAL64 ) dp = dp * (0.5_REAL64 / step)
             ptry = par + dp
             call evaluate(ptry,mt,jact)
             chi2try = sum(w*(d - mt - ptry(4))**2)
             if( chi2try <= chi2 ) then
                lambda = max(lambda / 10,1e-12_REAL64)
                exit
             else
                lambda = lambda * 10
             end if
          end do
          if( .not. ok ) exit
          if( chi2try > chi2 ) then
             ! no improvement: the minimum within numerical precision
             converged = .true.
             exit
          end if
          par = ptry

          if( abs(dp(2)) < 1e-4_REAL64 .and. abs(dp(3)) < 1e-4_REAL64 .and. &
               abs(dp(1)) < 1e-6_REAL64*abs(par(1)) + 1e-3_REAL64*sqrt(1/a(1,1)) &
               .and. abs(dp(4)) < 1e-4_REAL64*sqrt(max(1/a(4,4),0.0_REAL64)) ) then
             converged = .true.
             exit
          end if
       end do

       if( .not. posfixed .and. (.not. ok .or. &
            hypot(par(2) - x0,par(3) - y0) > flim .or. .not. (par(1) > 0)) ) then
          ! the fit run away (or failed): repeat with the fixed position
          posfixed = .true.
          star%x = x0
          star%y = y0
          call linear_flux(frame,psf,star,localsky)
          cycle
       end if
       exit
    end do

    if( .not. ok ) then
       star%flag = star%flag + FLAG_FAILED
       star%fluxerr = -1
       return
    end if

    star%flux = par(1)
    star%x = par(2)
    star%y = par(3)
    star%sky = par(4)
    if( posfixed ) star%flag = star%flag + FLAG_POSFIXED
    if( .not. converged ) star%flag = star%flag + FLAG_NOCONV

    ! covariance by the weights of the final model
    call evaluate(par,m,jac)
    w = 1 / frame%variance(b + m + par(4))
    call normal(m,w,chi2,a,g)
    call invert(a,cov,ok)
    if( .not. ok ) then
       star%flag = star%flag + FLAG_FAILED
       star%fluxerr = -1
       return
    end if

    star%fluxerr = sqrt(max(cov(1,1),0.0_REAL64))
    if( .not. posfixed ) then
       star%xerr = sqrt(max(cov(2,2),0.0_REAL64))
       star%yerr = sqrt(max(cov(3,3),0.0_REAL64))
    else
       star%xerr = -1
       star%yerr = -1
    end if
    if( localsky ) then
       star%skyerr = sqrt(max(cov(4,4),0.0_REAL64))
    else if( star%skyerr > 0 ) then
       ! the error of the fixed sky propagates to the parameters
       do i = 1, 4
          sb(i) = sum(w*jac(:,i))
       end do
       sb = matmul(cov,sb)
       star%fluxerr = hypot(star%fluxerr,sb(1)*star%skyerr)
    end if
    star%chi2 = chi2 / max(npix - count(free),1)

    ! sharpness: the relative residual in the core (positive for objects
    ! sharper than PSF, negative for extended), roundness: the difference
    ! of ellipticities of data and model
    sum1 = 0
    sum2 = 0
    do k = 1, npix
       if( hypot(ii(k) - par(2),jj(k) - par(3)) <= rcore ) then
          sum1 = sum1 + (d(k) - m(k) - par(4))
          sum2 = sum2 + m(k)
       end if
    end do
    if( sum2 > 0 ) then
       star%sharp = sum1 / sum2
    else
       star%sharp = 0
    end if
    star%round = ellipticity(d - par(4)) - ellipticity(m)

  contains

    subroutine evaluate(par,m,jac)

      ! the model (without sky) and its derivatives by parameters

      real(REAL64), dimension(:), intent(in) :: par
      real(REAL64), dimension(:), intent(out) :: m
      real(REAL64), dimension(:,:), intent(out) :: jac

      real(REAL64), dimension(-n:n,-n:n) :: p, pu, pv
      integer :: i0, j0, k, i, j

      i0 = nint(par(2))
      j0 = nint(par(3))
      call psf_pixelize(sub,nsub,psf%zoom,i0 - par(2),j0 - par(3),n,p,pu,pv)
      jac = 0
      do k = 1, npix
         i = ii(k) - i0
         j = jj(k) - j0
         if( abs(i) <= n .and. abs(j) <= n ) then
            m(k) = par(1)*p(i,j)
            jac(k,1) = p(i,j)
            jac(k,2) = -par(1)*pu(i,j)
            jac(k,3) = -par(1)*pv(i,j)
         else
            m(k) = 0
         end if
         jac(k,4) = 1
      end do

    end subroutine evaluate

    subroutine normal(m,w,chi2,a,g)

      ! normal equations for the free parameters, unit diagonal for others

      real(REAL64), dimension(:), intent(in) :: m, w
      real(REAL64), intent(out) :: chi2
      real(REAL64), dimension(4,4), intent(out) :: a
      real(REAL64), dimension(4), intent(out) :: g

      real(REAL64), dimension(npix) :: r
      integer :: i, j

      r = d - m - par(4)
      chi2 = sum(w*r**2)
      a = 0
      g = 0
      do i = 1, 4
         if( .not. free(i) ) then
            a(i,i) = 1
            cycle
         end if
         g(i) = sum(w*jac(:,i)*r)
         do j = 1, i
            if( .not. free(j) ) cycle
            a(i,j) = sum(w*jac(:,i)*jac(:,j))
            a(j,i) = a(i,j)
         end do
      end do

    end subroutine normal

    real(REAL64) function ellipticity(f)

      real(REAL64), dimension(:), intent(in) :: f
      real(REAL64) :: sxx, syy, sxy, u, v, q
      integer :: k

      sxx = 0
      syy = 0
      sxy = 0
      do k = 1, npix
         q = max(f(k),0.0_REAL64)
         u = ii(k) - par(2)
         v = jj(k) - par(3)
         sxx = sxx + q*u**2
         syy = syy + q*v**2
         sxy = sxy + q*u*v
      end do
      if( sxx + syy > 0 ) then
         ellipticity = sqrt((sxx - syy)**2 + 4*sxy**2) / (sxx + syy)
      else
         ellipticity = 0
      end if

    end function ellipticity

  end subroutine fit_star


  subroutine fit_group(frame,psf,stars,localsky,maxiter,snrpos,rcore,ok)

    ! The simultaneous fit of a group of stars with overlapping fitting
    ! regions: the fluxes and positions of all members (and a common sky
    ! with localsky) in the union of their fitting circles. The models
    ! of the members must not be included in mstar. The faint members
    ! (S/N < snrpos at the initial positions) and the members which run
    ! away (more than a half of the fitting radius, or a non-positive
    ! flux) are fitted at their initial positions. ok is false when the
    ! normal equations are singular; the caller fits the members alone.

    type(FrameType), intent(in) :: frame
    type(PsfType), intent(in) :: psf
    type(StarFit), dimension(:), intent(in out) :: stars
    logical, intent(in) :: localsky
    integer, intent(in) :: maxiter
    real(REAL64), intent(in) :: snrpos
    real(REAL64), dimension(:), intent(in) :: rcore
    logical, intent(out) :: ok

    type MemberPsf
       real(REAL64), dimension(:,:), allocatable :: sub
       integer :: nsub, n
    end type MemberPsf

    type(MemberPsf), dimension(:), allocatable :: mp
    real(REAL64), dimension(:), allocatable :: d, b, w, m, mt, par, ptry, g, &
         dp, sb, x0, y0, flim, mpart
    real(REAL64), dimension(:,:), allocatable :: jac, jact, a, ad, cov, mm
    integer, dimension(:), allocatable :: ii, jj
    logical, dimension(:), allocatable :: free, posfixed, runaway
    logical, dimension(:,:), allocatable :: inside
    real(REAL64) :: chi2, chi2try, lambda, step, sum1, sum2, skyerr, r2
    integer :: nm, np, npix, mx, k, l, i, j, iter, ntry, ix1, ix2, jy1, jy2, &
         nfree, restart
    logical :: converged

    ok = .false.
    nm = size(stars)
    np = 3*nm + 1

    allocate(x0(nm),y0(nm),flim(nm),posfixed(nm),runaway(nm),mp(nm))
    do l = 1, nm
       stars(l)%flag = 0
       x0(l) = stars(l)%x
       y0(l) = stars(l)%y
       call region_flags(frame,stars(l),rcore(l))
       flim(l) = max(stars(l)%fitrad / 2,1.0_REAL64)
    end do

    ! the union of the fitting circles (valid pixels)
    ix1 = huge(ix1)
    ix2 = -huge(ix2)
    jy1 = huge(jy1)
    jy2 = -huge(jy2)
    do l = 1, nm
       ix1 = min(ix1,floor(stars(l)%x - stars(l)%fitrad))
       ix2 = max(ix2,ceiling(stars(l)%x + stars(l)%fitrad))
       jy1 = min(jy1,floor(stars(l)%y - stars(l)%fitrad))
       jy2 = max(jy2,ceiling(stars(l)%y + stars(l)%fitrad))
    end do
    ix1 = max(ix1,1)
    jy1 = max(jy1,1)
    ix2 = min(ix2,size(frame%data,1))
    jy2 = min(jy2,size(frame%data,2))
    k = max((ix2 - ix1 + 1)*(jy2 - jy1 + 1),1)
    allocate(ii(k),jj(k),inside(k,nm))
    npix = 0
    do j = jy1, jy2
       do i = ix1, ix2
          if( .not. frame%valid(i,j) ) cycle
          npix = npix + 1
          ii(npix) = i
          jj(npix) = j
          do l = 1, nm
             r2 = (i - stars(l)%x)**2 + (j - stars(l)%y)**2
             inside(npix,l) = r2 <= stars(l)%fitrad**2
          end do
          if( .not. any(inside(npix,:)) ) npix = npix - 1
       end do
    end do
    do l = 1, nm
       stars(l)%npix = count(inside(1:npix,l))
    end do
    if( npix < np + 3 .or. any(stars%npix < 4) ) return

    ! the PSFs of members, the pixelisation boxes covering all pixels
    mx = max(ix2 - ix1,jy2 - jy1)
    do l = 1, nm
       call psf%sampled(stars(l)%x,stars(l)%y,mp(l)%sub,mp(l)%nsub)
       mp(l)%n = mx + ceiling(flim(l)) + 2
    end do

    allocate(d(npix),b(npix),w(npix),m(npix),mt(npix),mpart(npix), &
         jac(npix,np),jact(npix,np),mm(npix,nm),par(np),ptry(np),g(np), &
         dp(np),sb(np),a(np,np),ad(np,np),cov(np,np),free(np))
    do k = 1, npix
       d(k) = frame%data(ii(k),jj(k)) - frame%mstar(ii(k),jj(k))
       b(k) = frame%mstar(ii(k),jj(k))
    end do

    do l = 1, nm
       par(3*l-2:3*l) = [max(stars(l)%flux,0.0_REAL64), stars(l)%x, stars(l)%y]
    end do
    if( localsky ) then
       par(np) = sum(stars%sky) / nm
       skyerr = -1
    else
       par(np) = sum(stars%sky) / nm
       if( all(stars%skyerr > 0) ) then
          skyerr = sum(stars%skyerr) / nm
       else
          skyerr = -1
       end if
    end if

    ! the linear fit of fluxes at the initial positions decides on the
    ! position fits
    posfixed = .true.
    call set_free
    do iter = 1, 3
       call evaluate(par,m,jac)
       w = 1 / frame%variance(b + m + par(np))
       call normal(par,m,w,chi2,a,g)
       call solve(a,g,dp,ok)
       if( .not. ok ) return
       par = par + dp
    end do
    call invert(a,cov,ok)
    if( .not. ok ) return
    do l = 1, nm
       posfixed(l) = .not. (cov(3*l-2,3*l-2) > 0 .and. &
            par(3*l-2) > snrpos*sqrt(cov(3*l-2,3*l-2)))
    end do

    do restart = 1, nm + 1

       call set_free
       lambda = 1e-3_REAL64
       converged = .false.
       ok = .true.

       do iter = 1, maxiter

          call evaluate(par,m,jac)
          w = 1 / frame%variance(b + m + par(np))
          call normal(par,m,w,chi2,a,g)

          do ntry = 1, 12
             ad = a
             do i = 1, np
                ad(i,i) = a(i,i)*(1 + lambda)
             end do
             call solve(ad,g,dp,ok)
             if( .not. ok ) exit
             ! limit the shifts of positions
             step = 0
             do l = 1, nm
                step = max(step,abs(dp(3*l-1)),abs(dp(3*l)))
             end do
             if( step > 0.5_REAL64 ) dp = dp * (0.5_REAL64 / step)
             ptry = par + dp
             call evaluate(ptry,mt,jact)
             chi2try = sum(w*(d - mt - ptry(np))**2)
             if( chi2try <= chi2 ) then
                lambda = max(lambda / 10,1e-12_REAL64)
                exit
             else
                lambda = lambda * 10
             end if
          end do
          if( .not. ok ) exit
          if( chi2try > chi2 ) then
             converged = .true.
             exit
          end if
          par = ptry

          converged = abs(dp(np)) < 1e-4_REAL64*sqrt(max(1/a(np,np),0.0_REAL64))
          do l = 1, nm
             converged = converged .and. abs(dp(3*l-1)) < 1e-4_REAL64 .and. &
                  abs(dp(3*l)) < 1e-4_REAL64 .and. abs(dp(3*l-2)) < &
                  1e-6_REAL64*abs(par(3*l-2)) + 1e-3_REAL64*sqrt(1/a(3*l-2,3*l-2))
          end do
          if( converged ) exit
       end do
       if( .not. ok ) return

       ! members which run away are fitted at the initial positions
       do l = 1, nm
          runaway(l) = .not. posfixed(l) .and. &
               (hypot(par(3*l-1) - x0(l),par(3*l) - y0(l)) > flim(l) .or. &
               .not. (par(3*l-2) > 0))
       end do
       if( .not. any(runaway) ) exit
       do l = 1, nm
          if( runaway(l) ) then
             posfixed(l) = .true.
             par(3*l-1) = x0(l)
             par(3*l) = y0(l)
             par(3*l-2) = max(stars(l)%flux,0.0_REAL64)
          end if
       end do
    end do

    ! covariance by the weights of the final model
    call evaluate(par,m,jac)
    w = 1 / frame%variance(b + m + par(np))
    call normal(par,m,w,chi2,a,g)
    call invert(a,cov,ok)
    if( .not. ok ) return
    if( .not. localsky .and. skyerr > 0 ) then
       do i = 1, np
          sb(i) = sum(w*jac(:,i))
       end do
       sb = matmul(cov,sb)
    else
       sb = 0
    end if

    nfree = count(free)
    do l = 1, nm
       stars(l)%flux = par(3*l-2)
       stars(l)%x = par(3*l-1)
       stars(l)%y = par(3*l)
       stars(l)%sky = par(np)
       if( localsky ) then
          stars(l)%skyerr = sqrt(max(cov(np,np),0.0_REAL64))
       end if
       if( posfixed(l) ) stars(l)%flag = stars(l)%flag + FLAG_POSFIXED
       if( .not. converged ) stars(l)%flag = stars(l)%flag + FLAG_NOCONV
       stars(l)%fluxerr = sqrt(max(cov(3*l-2,3*l-2),0.0_REAL64))
       if( skyerr > 0 ) stars(l)%fluxerr = hypot(stars(l)%fluxerr,sb(3*l-2)*skyerr)
       if( .not. posfixed(l) ) then
          stars(l)%xerr = sqrt(max(cov(3*l-1,3*l-1),0.0_REAL64))
          stars(l)%yerr = sqrt(max(cov(3*l,3*l),0.0_REAL64))
       else
          stars(l)%xerr = -1
          stars(l)%yerr = -1
       end if

       ! chi2, sharpness and roundness in the circle of the member,
       ! the other members subtracted
       mpart = par(3*l-2)*jac(:,3*l-2)
       k = stars(l)%npix
       chi2 = sum(w*(d - m - par(np))**2,inside(1:npix,l))
       stars(l)%chi2 = chi2 / max(k - count(free(3*l-2:3*l)),1)
       sum1 = 0
       sum2 = 0
       do k = 1, npix
          if( .not. inside(k,l) ) cycle
          if( hypot(ii(k) - par(3*l-1),jj(k) - par(3*l)) <= rcore(l) ) then
             sum1 = sum1 + (d(k) - m(k) - par(np))
             sum2 = sum2 + mpart(k)
          end if
       end do
       if( sum2 > 0 ) then
          stars(l)%sharp = sum1 / sum2
       else
          stars(l)%sharp = 0
       end if
       stars(l)%round = ellipticity(d - par(np) - (m - mpart),l) - &
            ellipticity(mpart,l)
    end do
    ok = .true.

  contains

    subroutine set_free

      integer :: l

      free = .true.
      do l = 1, nm
         free(3*l-1:3*l) = .not. posfixed(l)
      end do
      free(np) = localsky

    end subroutine set_free

    subroutine evaluate(par,m,jac)

      ! the model of all members (without sky) and its derivatives

      real(REAL64), dimension(:), intent(in) :: par
      real(REAL64), dimension(:), intent(out) :: m
      real(REAL64), dimension(:,:), intent(out) :: jac

      real(REAL64), dimension(:,:), allocatable :: p, pu, pv
      integer :: i0, j0, k, i, j, l, n

      m = 0
      jac = 0
      do l = 1, nm
         n = mp(l)%n
         allocate(p(-n:n,-n:n),pu(-n:n,-n:n),pv(-n:n,-n:n))
         i0 = nint(par(3*l-1))
         j0 = nint(par(3*l))
         call psf_pixelize(mp(l)%sub,mp(l)%nsub,psf%zoom,i0 - par(3*l-1), &
              j0 - par(3*l),n,p,pu,pv)
         do k = 1, npix
            i = ii(k) - i0
            j = jj(k) - j0
            if( abs(i) <= n .and. abs(j) <= n ) then
               m(k) = m(k) + par(3*l-2)*p(i,j)
               jac(k,3*l-2) = p(i,j)
               jac(k,3*l-1) = -par(3*l-2)*pu(i,j)
               jac(k,3*l) = -par(3*l-2)*pv(i,j)
            end if
         end do
         deallocate(p,pu,pv)
      end do
      jac(:,np) = 1

    end subroutine evaluate

    subroutine normal(par,m,w,chi2,a,g)

      real(REAL64), dimension(:), intent(in) :: par, m, w
      real(REAL64), intent(out) :: chi2
      real(REAL64), dimension(:,:), intent(out) :: a
      real(REAL64), dimension(:), intent(out) :: g

      real(REAL64), dimension(npix) :: r
      integer :: i, j

      r = d - m - par(np)
      chi2 = sum(w*r**2)
      a = 0
      g = 0
      do i = 1, np
         if( .not. free(i) ) then
            a(i,i) = 1
            cycle
         end if
         g(i) = sum(w*jac(:,i)*r)
         do j = 1, i
            if( .not. free(j) ) cycle
            a(i,j) = sum(w*jac(:,i)*jac(:,j))
            a(j,i) = a(i,j)
         end do
      end do

    end subroutine normal

    real(REAL64) function ellipticity(f,l)

      real(REAL64), dimension(:), intent(in) :: f
      integer, intent(in) :: l
      real(REAL64) :: sxx, syy, sxy, u, v, q
      integer :: k

      sxx = 0
      syy = 0
      sxy = 0
      do k = 1, npix
         if( .not. inside(k,l) ) cycle
         q = max(f(k),0.0_REAL64)
         u = ii(k) - par(3*l-1)
         v = jj(k) - par(3*l)
         sxx = sxx + q*u**2
         syy = syy + q*v**2
         sxy = sxy + q*u*v
      end do
      if( sxx + syy > 0 ) then
         ellipticity = sqrt((sxx - syy)**2 + 4*sxy**2) / (sxx + syy)
      else
         ellipticity = 0
      end if

    end function ellipticity

  end subroutine fit_group


  subroutine region_flags(frame,star,rcore)

    ! saturated pixels in the core, excluded pixels in the fitting radius

    type(FrameType), intent(in) :: frame
    type(StarFit), intent(in out) :: star
    real(REAL64), intent(in) :: rcore

    integer :: i, j, n, i0, j0
    real(REAL64) :: r
    logical :: sat, excl

    n = ceiling(star%fitrad) + 1
    i0 = nint(star%x)
    j0 = nint(star%y)
    sat = .false.
    excl = .false.
    do j = j0 - n, j0 + n
       do i = i0 - n, i0 + n
          r = hypot(i - star%x,j - star%y)
          if( r > star%fitrad ) cycle
          if( i < 1 .or. j < 1 .or. &
               i > size(frame%data,1) .or. j > size(frame%data,2) ) then
             excl = .true.
             cycle
          end if
          if( frame%satur(i,j) .and. r <= rcore ) sat = .true.
          if( .not. frame%valid(i,j) ) excl = .true.
       end do
    end do
    if( sat ) star%flag = star%flag + FLAG_SATURATED
    if( excl ) star%flag = star%flag + FLAG_EXCLUDED

  end subroutine region_flags


  subroutine fit_shape(frame,psf,star,fitbeta,elliptic,maxiter,pshape,ok,betamax)

    ! The fit of an isolated star by the analytic PSF with free shape:
    ! pshape = [width, beta, c1, c2], where the axis ratio and the angle
    ! are q = 1 - hypot(c1,c2), pa = atan2(c2,c1)/2. The sky is local.
    ! The width (and beta, c1, c2) on input are the initial values.
    ! The optional betamax is the largest beta: a Moffat profile with
    ! beta = 10 is Gaussian to 1e-6 of the flux beyond 4 FWHM, while beta
    ! converges slowly towards infinity for Gaussian stars.

    type(FrameType), intent(in) :: frame
    type(PsfType), intent(in) :: psf
    type(StarFit), intent(in out) :: star
    logical, intent(in) :: fitbeta, elliptic
    integer, intent(in) :: maxiter
    real(REAL64), dimension(4), intent(in out) :: pshape
    logical, intent(out) :: ok

    integer, parameter :: npar = 8
    real(REAL64), parameter :: rad = 45 / atan(1.0_REAL64)
    real(REAL64), intent(in), optional :: betamax
    real(REAL64), dimension(:), allocatable :: d, w, m, mt
    real(REAL64), dimension(:,:), allocatable :: jac, jact
    integer, dimension(:), allocatable :: ii, jj
    real(REAL64), dimension(npar) :: par, ptry, g, dp
    real(REAL64), dimension(npar,npar) :: a, ad
    logical, dimension(npar) :: free
    real(REAL64) :: chi2, chi2try, lambda, step
    integer :: npix, n, k, iter, ntry, i, j
    logical :: converged

    ok = .false.
    call pixel_set(frame,star,ii,jj,npix,n)
    if( npix < npar + 5 ) return
    n = n + 3

    allocate(d(npix),w(npix),m(npix),mt(npix),jac(npix,npar),jact(npix,npar))
    do k = 1, npix
       d(k) = frame%data(ii(k),jj(k)) - frame%mstar(ii(k),jj(k))
    end do

    free = [.true., .true., .true., .true., .true., fitbeta, elliptic, elliptic]
    par = [star%flux, star%x, star%y, star%sky, pshape]
    lambda = 1e-3_REAL64
    converged = .false.

    do iter = 1, maxiter
       call evaluate(par,m,jac)
       if( any(m /= m) ) return
       w = 1 / frame%variance(m)
       call normal(m,w,chi2,a,g)

       do ntry = 1, 12
          ad = a
          do i = 1, npar
             ad(i,i) = a(i,i)*(1 + lambda)
          end do
          call solve(ad,g,dp,ok)
          if( .not. ok ) return
          step = maxval(abs(dp(2:3)))
          if( step > 0.5_REAL64 ) dp = dp * (0.5_REAL64 / step)
          ! keep the shape parameters in their domains
          if( dp(5) < -0.5_REAL64*par(5) ) dp = dp * (-0.5_REAL64*par(5) / dp(5))
          ptry = par + dp
          if( ptry(6) < 1.05_REAL64 ) ptry(6) = 1.05_REAL64
          if( present(betamax) ) ptry(6) = min(ptry(6),betamax)
          if( hypot(ptry(7),ptry(8)) > 0.8_REAL64 ) &
               ptry(7:8) = ptry(7:8) * 0.8_REAL64 / hypot(ptry(7),ptry(8))
          call evaluate(ptry,mt,jact)
          chi2try = sum(w*(d - mt)**2)
          if( chi2try <= chi2 ) then
             lambda = max(lambda / 10,1e-12_REAL64)
             exit
          else
             lambda = lambda * 10
          end if
       end do
       if( chi2try > chi2 ) then
          converged = .true.
          exit
       end if
       par = ptry
       if( present(betamax) ) then
          ! a Gaussian profile: beta at the limit
          if( free(6) .and. par(6) >= betamax ) free(6) = .false.
       end if
       if( maxval(abs(dp(2:3))) < 1e-4_REAL64 .and. &
            abs(dp(5)) < 1e-5_REAL64*par(5) .and. abs(dp(6)) < 1e-4_REAL64 .and. &
            maxval(abs(dp(7:8))) < 1e-5_REAL64 ) then
          converged = .true.
          exit
       end if
    end do

    ok = converged .and. par(1) > 0 .and. par(5) > 0
    if( ok ) then
       star%flux = par(1)
       star%x = par(2)
       star%y = par(3)
       star%sky = par(4)
       pshape = par(5:8)
    end if

  contains

    subroutine evaluate(par,m,jac)

      real(REAL64), dimension(:), intent(in) :: par
      real(REAL64), dimension(:), intent(out) :: m
      real(REAL64), dimension(:,:), intent(out) :: jac

      real(REAL64), dimension(-n:n,-n:n) :: p, pu, pv, p1, p2
      real(REAL64), dimension(npar) :: q
      real(REAL64) :: h
      integer :: i0, j0, k, i, j, l

      i0 = nint(par(2))
      j0 = nint(par(3))
      call model(par,i0,j0,p,pu,pv)
      jac = 0
      do k = 1, npix
         i = max(min(ii(k) - i0,n),-n)
         j = max(min(jj(k) - j0,n),-n)
         m(k) = par(1)*p(i,j) + par(4)
         jac(k,1) = p(i,j)
         jac(k,2) = -par(1)*pu(i,j)
         jac(k,3) = -par(1)*pv(i,j)
         jac(k,4) = 1
      end do

      ! numerical derivatives by the shape parameters
      do l = 5, npar
         if( .not. free(l) ) cycle
         if( l == 5 ) then
            h = 1e-4_REAL64*par(5)
         else
            h = 1e-4_REAL64
         end if
         q = par
         q(l) = par(l) + h
         call model(q,i0,j0,p1)
         q(l) = par(l) - h
         call model(q,i0,j0,p2)
         do k = 1, npix
            i = max(min(ii(k) - i0,n),-n)
            j = max(min(jj(k) - j0,n),-n)
            jac(k,l) = par(1)*(p1(i,j) - p2(i,j)) / (2*h)
         end do
      end do

    end subroutine evaluate

    subroutine model(par,i0,j0,p,pu,pv)

      real(REAL64), dimension(:), intent(in) :: par
      integer, intent(in) :: i0, j0
      real(REAL64), dimension(-n:n,-n:n), intent(out) :: p
      real(REAL64), dimension(-n:n,-n:n), intent(out), optional :: pu, pv

      real(REAL64), dimension(:,:), allocatable :: sub
      real(REAL64) :: q, pa
      integer :: nsub

      nsub = psf%zoom*(n + 2)
      allocate(sub(-nsub:nsub,-nsub:nsub))
      q = 1 - hypot(par(7),par(8))
      pa = rad*atan2(par(8),par(7)) / 2
      call psf_analytic_sample(psf%ptype,par(5),par(6),q,pa,psf%zoom,nsub,sub)
      if( present(pu) ) then
         call psf_pixelize(sub,nsub,psf%zoom,i0 - par(2),j0 - par(3),n,p,pu,pv)
      else
         call psf_pixelize(sub,nsub,psf%zoom,i0 - par(2),j0 - par(3),n,p)
      end if

    end subroutine model

    subroutine normal(m,w,chi2,a,g)

      real(REAL64), dimension(:), intent(in) :: m, w
      real(REAL64), intent(out) :: chi2
      real(REAL64), dimension(npar,npar), intent(out) :: a
      real(REAL64), dimension(npar), intent(out) :: g

      real(REAL64), dimension(npix) :: r
      integer :: i, j

      r = d - m
      chi2 = sum(w*r**2)
      a = 0
      g = 0
      do i = 1, npar
         if( .not. free(i) ) then
            a(i,i) = 1
            cycle
         end if
         g(i) = sum(w*jac(:,i)*r)
         do j = 1, i
            if( .not. free(j) ) cycle
            a(i,j) = sum(w*jac(:,i)*jac(:,j))
            a(j,i) = a(i,j)
         end do
      end do

    end subroutine normal

  end subroutine fit_shape


  subroutine solve(a,b,x,ok)

    ! solution of the symmetric positive definite system by Cholesky

    real(REAL64), dimension(:,:), intent(in) :: a
    real(REAL64), dimension(:), intent(in) :: b
    real(REAL64), dimension(:), intent(out) :: x
    logical, intent(out) :: ok

    real(REAL64), dimension(size(b),size(b)) :: l
    real(REAL64) :: s
    integer :: n, i, j

    n = size(b)
    l = 0
    ok = .false.
    do j = 1, n
       s = a(j,j) - sum(l(j,1:j-1)**2)
       if( .not. (s > 0) ) return
       l(j,j) = sqrt(s)
       do i = j + 1, n
          l(i,j) = (a(i,j) - sum(l(i,1:j-1)*l(j,1:j-1))) / l(j,j)
       end do
    end do
    do i = 1, n
       x(i) = (b(i) - sum(l(i,1:i-1)*x(1:i-1))) / l(i,i)
    end do
    do i = n, 1, -1
       x(i) = (x(i) - sum(l(i+1:n,i)*x(i+1:n))) / l(i,i)
    end do
    ok = .true.

  end subroutine solve


  subroutine invert(a,c,ok)

    real(REAL64), dimension(:,:), intent(in) :: a
    real(REAL64), dimension(:,:), intent(out) :: c
    logical, intent(out) :: ok

    real(REAL64), dimension(size(a,1)) :: e
    integer :: i

    do i = 1, size(a,1)
       e = 0
       e(i) = 1
       call solve(a,e,c(:,i),ok)
       if( .not. ok ) return
    end do

  end subroutine invert

end module psffit
