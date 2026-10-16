!
!  psfempir - the empirical PSF from the stars of a frame
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
!  The empirical PSF is the base (an analytic core fitted on bright
!  stars, or the PSF of a file) plus residual tables sampled with ZOOM
!  samples per pixel, the same representation as a PSF file:
!
!     psf(x,y) = base(x,y) + sum_t P_t(x,y) R_t,
!
!  where P_t are the terms of a polynomial of the order 0, 1, 2 in the
!  field coordinates (DAOPHOT VARIABLE = 0, 1, 2) and the powers of the
!  field radius (psffit, poly_terms).
!
!  The tables are not built by depositing the residuals of stars on
!  the grid, but by the weighted linear least squares with the exact
!  forward model: a pixel of a star is the integral of the interpolant
!  of the samples over the pixel (psf_pixel_weights), so
!
!     data - neighbours - sky - flux*base = flux * sum_t P_t A R_t,
!
!  where A is the (linear) pixelisation at the sub-pixel position of
!  the star. Every star samples the grid at its own phase; the solution
!  combines all of them without a sub-pixel bias (nor the smoothing of
!  an interpolation). The weights are the inverse variances of the
!  model, outliers (cosmic rays, neighbours not subtracted) are clipped
!  at 5 sigma. A sky offset of every star is fitted too, on a ring of
!  pixels just outside of the tables (inside, a sky error would be a
!  pedestal of the tables). A regularisation by the squared Laplacian,
!  relative to the local diagonal of the normal matrix and chosen by
!  the cross-validation over the stars, keeps the tables smooth where
!  the stars sample the grid sparsely, and pulls them to zero at the
!  edge of the table.
!
!  The normal equations are solved by the conjugate gradients with
!  the Jacobi preconditioner; the system is never formed.
!
!  The residual tables do not move the centroid and carry the flux
!  tsum only: every table is corrected by a combination of the base at
!  the field centre and its derivatives so that its sum is tsum and its
!  first moments are zero. Without the normalisation (tsum = 0) the
!  flux of the model is that of the base (the analytic core is
!  normalised to the unit flux over the plane, the light beyond the
!  table is in its wings); psfphot sets tsum by the aperture fluxes of
!  the PSF stars (empirical_setsum). The centre of the PSF is the
!  centroid of the base.
!
!  The tables built on subsets of the stars give the noise of the
!  tables: the uncertainty of the fluxes by the PSF (psffit,
!  model_error) and the shrinkage of a correction of a PSF file.
!

module psfempir

  use iso_fortran_env
  use psfmodel
  use psffit

  implicit none
  private

  type :: EmpStar
     integer :: n, i0, j0
     real(REAL64) :: flux
     real(REAL64), dimension(:), allocatable :: pt
     integer, dimension(:), allocatable :: sx, lx, sy, ly
     real(REAL64), dimension(:,:), allocatable :: wx, wy, t, w
  end type EmpStar

  public :: empirical_build, empirical_setsum, empirical_fluxshape


contains

  subroutine empirical_build(frame,psf,stars,rt,clip,cv,final,smooth,chi2, &
       resrms,niter)

    ! Builds the residual tables psf%res on the stars (fitted, their
    ! models in frame%mstar when modelled), using the pixels in the
    ! radius rt [pix]. The current tables are the initial solution and,
    ! with clip, the residuals by them clip outliers. The regularisation
    ! smooth is chosen by the cross-validation with cv (and returned).
    ! The final build keeps the tables of subsets of the stars (psf%resk)
    ! for the uncertainty of the PSF. chi2 are the reduced chi2 of stars
    ! by the new PSF, resrms is the rms of the constant table relative to
    ! the peak of the base.

    type(FrameType), intent(in) :: frame
    type(PsfType), intent(in out) :: psf
    type(StarFit), dimension(:), intent(in) :: stars
    real(REAL64), intent(in) :: rt
    real(REAL64), intent(in out) :: smooth
    logical, intent(in) :: clip, cv, final
    real(REAL64), dimension(:), intent(out) :: chi2
    real(REAL64), intent(out) :: resrms
    integer, intent(out) :: niter

    integer, parameter :: maxiter = 1000
    real(REAL64), parameter :: tol = 1e-4_REAL64
    real(REAL64) :: ctol = tol

    type(EmpStar), dimension(:), allocatable :: es
    real(REAL64), dimension(:,:,:), allocatable :: x, r, z, p, q, b, diag, lamw
    real(REAL64), dimension(:), allocatable :: xs, rs, zz, ps, qs, bs, ds, xsa
    real(REAL64), dimension(:,:,:,:), allocatable :: xk
    logical, dimension(:), allocatable :: inuse
    real(REAL64), dimension(:,:), allocatable :: sub, pold, pbase, reff, pm, cb, &
         supp
    real(REAL64), dimension(:), allocatable :: eta, zs, tsum
    real(REAL64), dimension(3,3) :: cm
    real(REAL64), dimension(3) :: mom, coef
    real(REAL64) :: rz, rznew, alpha, bnorm, rnorm, s, h, sig, rsky, sm
    integer :: ns, nt, ntab, zoom, npoly, k, i, j, l, t, n, nsub, iter, i0, j0, nk
    logical :: ok

    ns = size(stars)
    ntab = psf%ntab
    zoom = psf%zoom
    npoly = psf%npoly

    ! the support of the tables: the circle of the radius rt, the tables
    ! are zero outside (no data constrain them there)
    allocate(supp(-ntab:ntab,-ntab:ntab))
    call support(ntab,zoom,rt,supp)

    ! the data of stars: the pixels within the table radius and a ring
    ! outside of it, where the tables are zero: the ring fixes the sky of
    ! the star, fitted together with the tables (a sky error would be
    ! absorbed by a pedestal of the tables otherwise)
    allocate(es(ns))
    rsky = rt + max(3.0_REAL64,rt / 4)
    nt = ceiling(rsky) + 1
    do k = 1, ns
       associate(e => es(k), st => stars(k))
         e%n = nt
         e%flux = st%flux
         e%i0 = nint(st%x)
         e%j0 = nint(st%y)
         allocate(e%pt(npoly),e%sx(-nt:nt),e%lx(-nt:nt),e%sy(-nt:nt), &
              e%ly(-nt:nt),e%wx(0:zoom+2,-nt:nt),e%wy(0:zoom+2,-nt:nt), &
              e%t(-nt:nt,-nt:nt),e%w(-nt:nt,-nt:nt))
         call psf%poly_terms(st%x,st%y,e%pt)
         allocate(pm(0:zoom+2,-nt:nt))
         call psf_pixel_weights(ntab,zoom,nt,e%i0 - st%x,e%sx,e%lx,e%wx,pm)
         call psf_pixel_weights(ntab,zoom,nt,e%j0 - st%y,e%sy,e%ly,e%wy,pm)
         deallocate(pm)

         ! the model of the star in mstar (as add_model)
         n = st%modrad
         allocate(pold(-n:n,-n:n))
         pold = 0
         if( st%modelled ) then
            call psf%sampled(st%x,st%y,sub,nsub)
            call psf_pixelize(sub,nsub,zoom,e%i0 - st%x,e%j0 - st%y,n,pold)
            deallocate(sub)
         end if
         ! the base
         allocate(pbase(-nt:nt,-nt:nt))
         call psf%sampled_base(st%x,st%y,sub,nsub)
         call psf_pixelize(sub,nsub,zoom,e%i0 - st%x,e%j0 - st%y,nt,pbase)
         deallocate(sub)

         e%t = 0
         e%w = 0
         do j = -nt, nt
            do i = -nt, nt
               i0 = e%i0 + i
               j0 = e%j0 + j
               if( (i0 - st%x)**2 + (j0 - st%y)**2 > rsky**2 ) cycle
               if( i0 < 1 .or. j0 < 1 .or. i0 > size(frame%data,1) .or. &
                    j0 > size(frame%data,2) ) cycle
               if( .not. frame%valid(i0,j0) ) cycle
               s = frame%data(i0,j0) - frame%mstar(i0,j0) - st%sky
               if( abs(i) <= n .and. abs(j) <= n ) s = s + st%flux*pold(i,j)
               e%t(i,j) = s - st%flux*pbase(i,j)
               ! the weights by the variance of the model, as the fits
               ! (the floor of one electron in the Poisson term is not
               ! needed with the regularisation relative to the local
               ! diagonal: with the mean of the diagonal, the exact
               ! variance of a sky of 0.4 e-/pixel gave the outer zone
               ! more weight and the empirical PSF of Gaussian stars came
               ! out 0.2% low in flux, BUGS P35)
               e%w(i,j) = 1 / frame%variance(frame%mstar(i0,j0) + st%sky)
            end do
         end do
         deallocate(pold,pbase)
       end associate
    end do

    allocate(x(-ntab:ntab,-ntab:ntab,npoly))
    x = psf%res
    do t = 1, npoly
       x(:,:,t) = supp*x(:,:,t)
    end do

    ! clipping of outliers by the current tables
    if( clip ) then
       allocate(reff(-ntab:ntab,-ntab:ntab),pm(-nt:nt,-nt:nt))
       do k = 1, ns
          associate(e => es(k))
            call combine(x,e%pt,reff)
            call forward(e,reff,pm)
            pm = (e%t - e%flux*pm)*sqrt(e%w)
            zs = pack(pm,e%w > 0)
            if( size(zs) < 10 ) cycle
            h = median(zs)
            sig = 1.4826_REAL64*median(abs(zs - h))
            sig = 5*max(sig,1.0_REAL64)
            where( abs(pm - h) > sig ) e%w = 0
          end associate
       end do
       deallocate(reff,pm)
    end if

    allocate(b(-ntab:ntab,-ntab:ntab,npoly),diag(-ntab:ntab,-ntab:ntab,npoly), &
         r(-ntab:ntab,-ntab:ntab,npoly),z(-ntab:ntab,-ntab:ntab,npoly), &
         p(-ntab:ntab,-ntab:ntab,npoly),q(-ntab:ntab,-ntab:ntab,npoly), &
         lamw(-ntab:ntab,-ntab:ntab,npoly), &
         reff(-ntab:ntab,-ntab:ntab),pm(-nt:nt,-nt:nt),eta(npoly),tsum(npoly), &
         xs(ns),rs(ns),zz(ns),ps(ns),qs(ns),bs(ns),ds(ns),inuse(ns))

    ! the flux of the tables (by the normalisation, zero without it)
    tsum = 0
    if( allocated(psf%tsum) ) then
       if( size(psf%tsum) == npoly ) tsum = psf%tsum
    end if

    ! the base at the field centre and its derivatives (the constraints)
    call psf%sampled_base(psf%pxc,psf%pyc,sub,nsub)
    allocate(cb(-ntab:ntab,-ntab:ntab),pold(-ntab:ntab,-ntab:ntab), &
         pbase(-ntab:ntab,-ntab:ntab))
    call constraint_basis(ntab,sub,nsub,supp,cb,pold,pbase,cm)

    ! the tables by all stars, the regularisation chosen by the
    ! cross-validation
    xs = 0
    inuse = .true.
    sm = smooth
    if( cv ) call cv_smooth(x,xs)
    call cgsolve(x,xs,niter)
    call project(x)

    ! The tables by nk subsets of the stars (interleaved by brightness)
    ! estimate the noise of the tables: the uncertainty of the fluxes by
    ! the PSF (psf%resk, model_error), and the shrinkage of a correction of
    ! a PSF file: in every radial zone and for every term, the table is
    ! kept where it is significant and shrinks towards zero where it is
    ! not (see shrink), so an exact base is not changed by the noise of
    ! the stars.
    nk = min(8,ns / 4)
    if( allocated(psf%resk) ) deallocate(psf%resk)
    if( nk >= 2 .and. (psf%ptype == PSF_FILE .or. final) ) then
       allocate(xk(-ntab:ntab,-ntab:ntab,npoly,nk),xsa(ns))
       do l = 1, nk
          xk(:,:,:,l) = x
          xsa = xs
          inuse = [(mod(k - 1,nk) == l - 1, k = 1, ns)]
          call cgsolve(xk(:,:,:,l),xsa,i)
          call project(xk(:,:,:,l))
       end do
       inuse = .true.
       if( psf%ptype == PSF_FILE ) then
          call shrink(x,xk)
          call project(x)
          do l = 1, nk
             call project(xk(:,:,:,l))
          end do
       end if
       if( final ) call move_alloc(xk,psf%resk)
    end if
    smooth = sm

    h = maxval(cb)
    resrms = sqrt(sum(x(:,:,1)**2,supp > 0) / count(supp > 0)) / h
    psf%res = x

    ! chi2 of stars by the new tables
    do k = 1, ns
       associate(e => es(k))
         call combine(x,e%pt,reff)
         call forward(e,reff,pm)
         l = count(e%w > 0)
         chi2(k) = sum(e%w*(e%t - e%flux*pm - xs(k))**2) / max(l - 4,1)
       end associate
    end do

  contains

    subroutine cgsolve(x,xs,niter)

      ! the tables x and the skies xs by the stars inuse (the initial
      ! values on input), by the preconditioned conjugate gradients

      real(REAL64), dimension(-ntab:,-ntab:,:), intent(in out) :: x
      real(REAL64), dimension(:), intent(in out) :: xs
      integer, intent(out) :: niter

      real(REAL64) :: rz, rznew, alpha, bnorm, rnorm, s
      integer :: k, t, iter

      ! the right side, the diagonal of the normal matrix
      b = 0
      diag = 0
      bs = 0
      ds = 1
      do k = 1, ns
         if( .not. inuse(k) ) cycle
         associate(e => es(k))
           bs(k) = sum(e%w*e%t)
           ds(k) = max(sum(e%w),tiny(1.0_REAL64))
           pm = e%flux*e%w*e%t
           call adjoint(e,pm,reff,.false.)
           do t = 1, npoly
              b(:,:,t) = b(:,:,t) + e%pt(t)*reff
           end do
           pm = e%flux**2*e%w
           call adjoint(e,pm,reff,.true.)
           do t = 1, npoly
              diag(:,:,t) = diag(:,:,t) + e%pt(t)**2*reff
           end do
         end associate
      end do
      where( .not. inuse ) xs = 0

      ! the regularisation relative to the local diagonal, the
      ! information of every sample (smoothed over the neighbours): the
      ! wing samples of bright stars carry orders of magnitude more
      ! weight than those of the core, a regularisation relative to the
      ! mean of the diagonal (as before) smoothed the core strongly
      do t = 1, npoly
         s = sum(diag(:,:,t),diag(:,:,t) > 0) / max(count(diag(:,:,t) > 0),1)
         call box3(diag(:,:,t),lamw(:,:,t))
         lamw(:,:,t) = sm*max(lamw(:,:,t),1e-6_REAL64*s)
         eta(t) = 1e-6_REAL64*s
         diag(:,:,t) = diag(:,:,t) + 20*lamw(:,:,t) + eta(t)
         b(:,:,t) = supp*b(:,:,t)
      end do

      call operator(x,xs,q,qs)
      r = b - q
      rs = bs - qs
      z = r / diag
      zz = rs / ds
      p = z
      ps = zz
      rz = sum(r*z) + sum(rs*zz)
      bnorm = sqrt(sum(b**2) + sum(bs**2))
      niter = 0
      do iter = 1, maxiter
         niter = iter
         rnorm = sqrt(sum(r**2) + sum(rs**2))
         if( .not. (rnorm > ctol*bnorm) ) exit
         call operator(p,ps,q,qs)
         s = sum(p*q) + sum(ps*qs)
         if( .not. (s > 0) ) exit
         alpha = rz / s
         x = x + alpha*p
         xs = xs + alpha*ps
         r = r - alpha*q
         rs = rs - alpha*qs
         z = r / diag
         zz = rs / ds
         rznew = sum(r*z) + sum(rs*zz)
         p = z + (rznew / rz)*p
         ps = zz + (rznew / rz)*ps
         rz = rznew
      end do

    end subroutine cgsolve

    subroutine cv_smooth(x,xs)

      ! The regularisation by the K-fold cross-validation (K = 4): for a
      ! value, the tables are built on the stars out of a fold and predict
      ! the pixels of the stars of the fold (their sky offsets fitted); the
      ! score is the sum of the squared normalised residuals of all folds.
      ! The values are the decades 0.001 .. 10 (relative to the local
      ! diagonal of the normal matrix); the search starts at 0.1 and walks
      ! to the lower score until the minimum is bracketed (the score is
      ! unimodal in practice). The solutions of the folds are less precise
      ! (1e-3) than the final one. Too few stars keep sm.

      real(REAL64), dimension(-ntab:,-ntab:,:), intent(in) :: x
      real(REAL64), dimension(:), intent(in) :: xs

      integer, parameter :: kf = 4, glo = -3, ghi = 1
      real(REAL64), dimension(:,:,:), allocatable :: xc
      real(REAL64), dimension(:), allocatable :: xsc
      real(REAL64), dimension(glo:ghi) :: score
      logical, dimension(glo:ghi) :: done
      integer :: g, step

      if( ns < 4*kf ) return
      allocate(xc(-ntab:ntab,-ntab:ntab,npoly),xsc(ns))
      done = .false.
      ctol = 1e-3_REAL64
      g = -1
      call cv_eval(x,xs,xc,xsc,g,score,done)
      call cv_eval(x,xs,xc,xsc,g + 1,score,done)
      if( score(g + 1) < score(g) ) then
         step = 1
         g = g + 1
      else
         step = -1
      end if
      do
         if( g + step < glo .or. g + step > ghi ) exit
         call cv_eval(x,xs,xc,xsc,g + step,score,done)
         if( .not. (score(g + step) < score(g)) ) exit
         g = g + step
      end do
      ctol = tol
      inuse = .true.
      sm = 10.0_REAL64**g


    end subroutine cv_smooth

    subroutine cv_eval(x,xs,xc,xsc,g,score,done)

      ! the score of the cross-validation of the regularisation 10**g

      real(REAL64), dimension(-ntab:,-ntab:,:), intent(in) :: x
      real(REAL64), dimension(:), intent(in) :: xs
      real(REAL64), dimension(-ntab:,-ntab:,:), intent(in out) :: xc
      real(REAL64), dimension(:), intent(in out) :: xsc
      integer, intent(in) :: g
      real(REAL64), dimension(-3:), intent(in out) :: score
      logical, dimension(-3:), intent(in out) :: done
      integer, parameter :: kf = 4
      real(REAL64) :: sw, sk
      integer :: f, k, it

      if( done(g) ) return
      sm = 10.0_REAL64**g
      score(g) = 0
      do f = 1, kf
         xc = x
         xsc = xs
         inuse = [(mod(k - 1,kf) /= f - 1, k = 1, ns)]
         call cgsolve(xc,xsc,it)
         call project(xc)
         do k = 1, ns
            if( inuse(k) ) cycle
            associate(e => es(k))
              call combine(xc,e%pt,reff)
              call forward(e,reff,pm)
              pm = e%t - e%flux*pm
              sw = sum(e%w)
              sk = 0
              if( sw > 0 ) sk = sum(e%w*pm) / sw
              score(g) = score(g) + sum(e%w*(pm - sk)**2)
            end associate
         end do
      end do
      done(g) = .true.

    end subroutine cv_eval

    subroutine box3(v,a)

      ! the mean of the 3x3 neighbourhood (inside of the table)

      real(REAL64), dimension(-ntab:,-ntab:), intent(in) :: v
      real(REAL64), dimension(-ntab:,-ntab:), intent(out) :: a
      real(REAL64), dimension(-ntab:ntab,-ntab:ntab) :: c, u
      integer :: i, j

      u = v
      c = 1
      a = u
      do j = -1, 1
         do i = -1, 1
            if( i == 0 .and. j == 0 ) cycle
            a(max(-ntab,-ntab-i):min(ntab,ntab-i),max(-ntab,-ntab-j):min(ntab,ntab-j)) = &
                 a(max(-ntab,-ntab-i):min(ntab,ntab-i),max(-ntab,-ntab-j):min(ntab,ntab-j)) + &
                 u(max(-ntab,-ntab-i)+i:min(ntab,ntab-i)+i,max(-ntab,-ntab-j)+j:min(ntab,ntab-j)+j)
            c(max(-ntab,-ntab-i):min(ntab,ntab-i),max(-ntab,-ntab-j):min(ntab,ntab-j)) = &
                 c(max(-ntab,-ntab-i):min(ntab,ntab-i),max(-ntab,-ntab-j):min(ntab,ntab-j)) + 1
         end do
      end do
      a = a / c

    end subroutine box3

    subroutine project(x)

      ! no shift, the flux of the tables tsum: the tables minus
      ! a combination of the base at the field centre and its derivatives
      ! have the sums tsum and zero first moments; the excess flux of the
      ! constant term is the error of the fluxes of the stars it was built
      ! by, so the tables are rescaled by it (the other terms by the same
      ! factor)

      real(REAL64), dimension(-ntab:,-ntab:,:), intent(in out) :: x

      call project_tables(ntab,x,tsum,supp,cb,pold,pbase,cm,sum(sub),.true.)

    end subroutine project

    subroutine shrink(x,xk)

      ! The shrinkage of the tables by radial zones: in every zone, the
      ! mean (the coherent part, it carries the flux) and the fluctuation
      ! about it are multiplied by a ramp of their power rho relative to
      ! the variance estimated by the scatter of the tables of the nk
      ! subsets: zero below the 90% quantile of Student's t of nk - 1
      ! degrees of freedom (rho = 2.0 for nk = 8), one above the 97.5%
      ! quantile (rho = 5.6), linear between; a zone without a significant
      ! correction is zero, a significant one is kept as it is. The tables
      ! of the subsets are shrunk by the same factors, so their scatter is
      ! the noise of the shrunk tables.

      real(REAL64), dimension(-ntab:,-ntab:,:), intent(in out) :: x
      real(REAL64), dimension(-ntab:,-ntab:,:,:), intent(in out) :: xk
      real(REAL64), dimension(*), parameter :: edges = [0.0_REAL64, 1.0_REAL64, &
           1.5_REAL64, 2.0_REAL64, 3.0_REAL64, 4.0_REAL64, 6.0_REAL64, 8.0_REAL64, &
           11.0_REAL64, 15.0_REAL64, 21.0_REAL64, 30.0_REAL64, 1e9_REAL64]
      ! Student's t quantiles 0.90 and 0.975 of 1..7 degrees of freedom
      real(REAL64), dimension(7), parameter :: t90 = [3.078_REAL64, &
           1.886_REAL64, 1.638_REAL64, 1.533_REAL64, 1.476_REAL64, 1.440_REAL64, &
           1.415_REAL64], t975 = [12.706_REAL64, 4.303_REAL64, 3.182_REAL64, &
           2.776_REAL64, 2.571_REAL64, 2.447_REAL64, 2.365_REAL64]
      logical, dimension(-ntab:ntab,-ntab:ntab) :: zone
      real(REAL64), dimension(-ntab:ntab,-ntab:ntab) :: rr, xbar, dev
      real(REAL64), dimension(size(xk,4)) :: mk
      real(REAL64) :: mx, mbar, vm, vf, fm, ff, rlo, rhi
      real(REAL64), parameter :: fluxfloor = 5e-4_REAL64
      integer :: i, j, t, m, n, kk, nk

      nk = size(xk,4)
      rlo = t90(min(nk - 1,7))**2
      rhi = t975(min(nk - 1,7))**2
      do j = -ntab, ntab
         do i = -ntab, ntab
            rr(i,j) = hypot(real(i,REAL64),real(j,REAL64)) / zoom
         end do
      end do
      do t = 1, npoly
         xbar = sum(xk(:,:,t,:),dim=3) / nk
         do m = 1, size(edges) - 1
            zone = rr >= edges(m) .and. rr < edges(m+1) .and. supp > 0
            n = count(zone)
            if( n == 0 ) cycle
            mx = sum(x(:,:,t),zone) / n
            do kk = 1, nk
               mk(kk) = sum(xk(:,:,t,kk),zone) / n
            end do
            mbar = sum(mk) / nk
            ! the variances of the mean of all stars
            vm = sum((mk - mbar)**2) / (nk*(nk - 1))
            ! the systematic floor: 0.05% of the flux in the zone (the
            ! background and the crowding of the frame are not known
            ! better, a smaller change of the profile is not trusted)
            vm = vm + (fluxfloor*zoom**2 / n)**2
            vf = 0
            do kk = 1, nk
               dev = xk(:,:,t,kk) - mk(kk) - (xbar - mbar)
               vf = vf + sum(dev**2,zone)
            end do
            vf = vf / (nk*(nk - 1))
            fm = ramp(mx**2,vm,rlo,rhi)
            ff = ramp(sum((x(:,:,t) - mx)**2,zone),vf,rlo,rhi)
            where( zone ) x(:,:,t) = fm*mx + ff*(x(:,:,t) - mx)
            do kk = 1, nk
               where( zone ) xk(:,:,t,kk) = fm*mk(kk) + ff*(xk(:,:,t,kk) - mk(kk))
            end do
         end do
      end do

    end subroutine shrink

    real(REAL64) function ramp(power,noise,rlo,rhi)

      ! zero for power/noise below rlo, one above rhi, linear between

      real(REAL64), intent(in) :: power, noise, rlo, rhi

      if( noise > 0 ) then
         ramp = min(1.0_REAL64,max(0.0_REAL64,(power / noise - rlo) / (rhi - rlo)))
      else
         ramp = 1
      end if

    end function ramp

    subroutine combine(x,pt,reff)

      real(REAL64), dimension(-ntab:,-ntab:,:), intent(in) :: x
      real(REAL64), dimension(:), intent(in) :: pt
      real(REAL64), dimension(-ntab:,-ntab:), intent(out) :: reff
      integer :: t

      reff = pt(1)*x(:,:,1)
      do t = 2, npoly
         reff = reff + pt(t)*x(:,:,t)
      end do

    end subroutine combine

    subroutine operator(v,vs,av,avs)

      ! the normal matrix times (v,vs): sum over stars of A'WA, the
      ! polynomial terms, the skies of stars vs, and the regularisation

      real(REAL64), dimension(-ntab:,-ntab:,:), intent(in) :: v
      real(REAL64), dimension(:), intent(in) :: vs
      real(REAL64), dimension(-ntab:,-ntab:,:), intent(out) :: av
      real(REAL64), dimension(:), intent(out) :: avs

      real(REAL64), dimension(-ntab:ntab,-ntab:ntab) :: reff, g, lap
      real(REAL64), dimension(-nt:nt,-nt:nt) :: pm
      integer :: k, t

      ! (the input is zero outside of the support)
      av = 0
      avs = 0
      do k = 1, ns
         if( .not. inuse(k) ) cycle
         associate(e => es(k))
           call combine(v,e%pt,reff)
           call forward(e,reff,pm)
           pm = e%w*(e%flux*pm + vs(k))
           avs(k) = sum(pm)
           pm = e%flux*pm
           call adjoint(e,pm,g,.false.)
           do t = 1, npoly
              av(:,:,t) = av(:,:,t) + e%pt(t)*g
           end do
         end associate
      end do
      do t = 1, npoly
         call laplace(v(:,:,t),lap)
         lap = lamw(:,:,t)*lap
         call laplace(lap,g)
         av(:,:,t) = supp*(av(:,:,t) + g + eta(t)*v(:,:,t))
      end do

    end subroutine operator

    subroutine forward(e,reff,pm)

      ! the pixels of the star by the samples reff (unit flux)

      type(EmpStar), intent(in) :: e
      real(REAL64), dimension(-ntab:,-ntab:), intent(in) :: reff
      real(REAL64), dimension(-nt:,-nt:), intent(out) :: pm

      real(REAL64), dimension(-ntab:ntab,-nt:nt) :: tmp
      integer :: i, j, l

      do j = -nt, nt
         tmp(:,j) = 0
         do l = 0, e%ly(j) - 1
            tmp(:,j) = tmp(:,j) + e%wy(l,j)*reff(:,e%sy(j)+l)
         end do
      end do
      pm = 0
      do j = -nt, nt
         if( e%ly(j) == 0 ) cycle
         do i = -nt, nt
            if( e%lx(i) == 0 ) cycle
            pm(i,j) = dot_product(e%wx(0:e%lx(i)-1,i), &
                 tmp(e%sx(i):e%sx(i)+e%lx(i)-1,j)) / zoom**2
         end do
      end do

    end subroutine forward

    subroutine adjoint(e,pm,g,squared)

      ! the transposed pixelisation of the pixels pm to the samples g
      ! (squared: of the squares of the weights, for the diagonal)

      type(EmpStar), intent(in) :: e
      real(REAL64), dimension(-nt:,-nt:), intent(in) :: pm
      real(REAL64), dimension(-ntab:,-ntab:), intent(out) :: g
      logical, intent(in) :: squared

      real(REAL64), dimension(-ntab:ntab,-nt:nt) :: tmp
      integer :: i, j, l, k1, k2

      tmp = 0
      do j = -nt, nt
         if( e%ly(j) == 0 ) cycle
         do i = -nt, nt
            if( e%lx(i) == 0 .or. pm(i,j) == 0 ) cycle
            k1 = e%sx(i)
            k2 = e%sx(i) + e%lx(i) - 1
            if( squared ) then
               tmp(k1:k2,j) = tmp(k1:k2,j) + e%wx(0:e%lx(i)-1,i)**2*pm(i,j) / zoom**4
            else
               tmp(k1:k2,j) = tmp(k1:k2,j) + e%wx(0:e%lx(i)-1,i)*pm(i,j) / zoom**2
            end if
         end do
      end do
      g = 0
      do j = -nt, nt
         do l = 0, e%ly(j) - 1
            if( squared ) then
               g(:,e%sy(j)+l) = g(:,e%sy(j)+l) + e%wy(l,j)**2*tmp(:,j)
            else
               g(:,e%sy(j)+l) = g(:,e%sy(j)+l) + e%wy(l,j)*tmp(:,j)
            end if
         end do
      end do

    end subroutine adjoint

    subroutine laplace(v,lap)

      ! the five-point Laplacian, zero outside of the table

      real(REAL64), dimension(-ntab:,-ntab:), intent(in) :: v
      real(REAL64), dimension(-ntab:,-ntab:), intent(out) :: lap

      lap = -4*v
      lap(-ntab+1:,:) = lap(-ntab+1:,:) + v(:ntab-1,:)
      lap(:ntab-1,:) = lap(:ntab-1,:) + v(-ntab+1:,:)
      lap(:,-ntab+1:) = lap(:,-ntab+1:) + v(:,:ntab-1)
      lap(:,:ntab-1) = lap(:,:ntab-1) + v(:,-ntab+1:)

    end subroutine laplace


  end subroutine empirical_build


  subroutine empirical_setsum(psf,rt)

    ! The tables (and those of the subsets) projected to the sums
    ! psf%tsum (the normalisation) with zero first moments, by the base
    ! at the field centre in the radius rt [pix].

    type(PsfType), intent(in out) :: psf
    real(REAL64), intent(in) :: rt

    real(REAL64), dimension(:,:), allocatable :: sub, cb, dbx, dby, supp
    real(REAL64), dimension(3,3) :: cm
    integer :: nsub, ntab, l

    if( .not. allocated(psf%tsum) ) return
    ntab = psf%ntab
    allocate(supp(-ntab:ntab,-ntab:ntab),cb(-ntab:ntab,-ntab:ntab), &
         dbx(-ntab:ntab,-ntab:ntab),dby(-ntab:ntab,-ntab:ntab))
    call support(ntab,psf%zoom,rt,supp)
    call psf%sampled_base(psf%pxc,psf%pyc,sub,nsub)
    call constraint_basis(ntab,sub,nsub,supp,cb,dbx,dby,cm)
    call project_tables(ntab,psf%res,psf%tsum,supp,cb,dbx,dby,cm,sum(sub),.false.)
    if( allocated(psf%resk) ) then
       do l = 1, size(psf%resk,4)
          call project_tables(ntab,psf%resk(:,:,:,l),psf%tsum,supp,cb,dbx,dby,cm, &
               sum(sub),.false.)
       end do
    end if

  end subroutine empirical_setsum


  subroutine empirical_fluxshape(psf,rt,cbn)

    ! the shape in which empirical_setsum adds flux to the tables: the
    ! base at the field centre in the radius rt [pix], unit sum

    type(PsfType), intent(in) :: psf
    real(REAL64), intent(in) :: rt
    real(REAL64), dimension(:,:), allocatable, intent(out) :: cbn

    real(REAL64), dimension(:,:), allocatable :: sub, supp
    integer :: nsub, ntab

    ntab = psf%ntab
    allocate(supp(-ntab:ntab,-ntab:ntab),cbn(-ntab:ntab,-ntab:ntab))
    call support(ntab,psf%zoom,rt,supp)
    call psf%sampled_base(psf%pxc,psf%pyc,sub,nsub)
    cbn = supp*sub(-ntab:ntab,-ntab:ntab)
    cbn = cbn / sum(cbn)

  end subroutine empirical_fluxshape


  subroutine support(ntab,zoom,rt,supp)

    ! the support of the tables: the circle of the radius rt [pix]

    integer, intent(in) :: ntab, zoom
    real(REAL64), intent(in) :: rt
    real(REAL64), dimension(-ntab:,-ntab:), intent(out) :: supp
    integer :: i, j

    do j = -ntab, ntab
       do i = -ntab, ntab
          if( i**2 + j**2 <= (zoom*rt)**2 ) then
             supp(i,j) = 1
          else
             supp(i,j) = 0
          end if
       end do
    end do

  end subroutine support


  subroutine constraint_basis(ntab,sub,nsub,supp,cb,dbx,dby,cm)

    ! the base sub at the field centre and its derivatives on the support
    ! (the directions removed from the tables), their sums and moments cm

    integer, intent(in) :: ntab, nsub
    real(REAL64), dimension(-nsub:,-nsub:), intent(in) :: sub
    real(REAL64), dimension(-ntab:,-ntab:), intent(in) :: supp
    real(REAL64), dimension(-ntab:,-ntab:), intent(out) :: cb, dbx, dby
    real(REAL64), dimension(3,3), intent(out) :: cm

    cb = sub(-ntab:ntab,-ntab:ntab)
    dbx = 0
    dby = 0
    dbx(-ntab+1:ntab-1,:) = (cb(-ntab+2:ntab,:) - cb(-ntab:ntab-2,:)) / 2
    dby(:,-ntab+1:ntab-1) = (cb(:,-ntab+2:ntab) - cb(:,-ntab:ntab-2)) / 2
    cb = supp*cb
    dbx = supp*dbx
    dby = supp*dby
    call moments(ntab,cb,cm(:,1))
    call moments(ntab,dbx,cm(:,2))
    call moments(ntab,dby,cm(:,3))

  end subroutine constraint_basis


  subroutine project_tables(ntab,x,tsum,supp,cb,dbx,dby,cm,flux,rescale)

    ! The tables minus a combination of the base cb and its derivatives
    ! dbx, dby have the sums tsum and zero first moments (no flux beyond
    ! the normalisation, no shift). With rescale, the excess flux of the
    ! constant table (the error of the fluxes of the stars it was built
    ! by; flux is the flux of the base) rescales all the tables.

    integer, intent(in) :: ntab
    real(REAL64), dimension(-ntab:,-ntab:,:), intent(in out) :: x
    real(REAL64), dimension(:), intent(in) :: tsum
    real(REAL64), dimension(-ntab:,-ntab:), intent(in) :: supp, cb, dbx, dby
    real(REAL64), dimension(3,3), intent(in) :: cm
    real(REAL64), intent(in) :: flux
    logical, intent(in) :: rescale

    real(REAL64), dimension(3) :: mom, coef
    real(REAL64) :: h
    integer :: t, pass
    logical :: ok

    h = 1
    if( rescale ) then
       call moments(ntab,x(:,:,1),mom)
       mom(1) = mom(1) - tsum(1)
       call solve3(cm,mom,coef,ok)
       if( ok ) h = 1 / (1 + coef(1)*sum(cb) / flux)
    end if
    do pass = 1, 2
       do t = 1, size(x,3)
          call moments(ntab,x(:,:,t),mom)
          mom(1) = mom(1) - tsum(t)
          call solve3(cm,mom,coef,ok)
          if( ok ) x(:,:,t) = h*(x(:,:,t) - coef(1)*cb - coef(2)*dbx - coef(3)*dby)
       end do
       ! the rescaled sums are restored
       if( h == 1 ) exit
       h = 1
    end do
    do t = 1, size(x,3)
       x(:,:,t) = supp*x(:,:,t)
    end do

  end subroutine project_tables


  subroutine moments(ntab,v,mom)

    ! the sum and the first moments of the table

    integer, intent(in) :: ntab
    real(REAL64), dimension(-ntab:,-ntab:), intent(in) :: v
    real(REAL64), dimension(3), intent(out) :: mom
    integer :: i

    mom(1) = sum(v)
    mom(2) = 0
    mom(3) = 0
    do i = -ntab, ntab
       mom(2) = mom(2) + i*sum(v(i,:))
       mom(3) = mom(3) + i*sum(v(:,i))
    end do

  end subroutine moments


  subroutine solve3(a,b,x,ok)

    ! a x = b by Gaussian elimination with pivoting

    real(REAL64), dimension(3,3), intent(in) :: a
    real(REAL64), dimension(3), intent(in) :: b
    real(REAL64), dimension(3), intent(out) :: x
    logical, intent(out) :: ok

    real(REAL64), dimension(3,4) :: m
    real(REAL64), dimension(4) :: row
    integer :: i, k, l

    m(:,1:3) = a
    m(:,4) = b
    ok = .false.
    do k = 1, 3
       l = k - 1 + maxloc(abs(m(k:3,k)),1)
       if( .not. (abs(m(l,k)) > 0) ) return
       row = m(k,:)
       m(k,:) = m(l,:)
       m(l,:) = row
       do i = k + 1, 3
          m(i,:) = m(i,:) - m(i,k) / m(k,k) * m(k,:)
       end do
    end do
    do i = 3, 1, -1
       x(i) = (m(i,4) - sum(m(i,i+1:3)*x(i+1:3))) / m(i,i)
    end do
    ok = .true.

  end subroutine solve3


  function median(x)

    real(REAL64) :: median
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), allocatable :: y
    real(REAL64) :: t
    integer :: n, i, j, k

    n = size(x)
    if( n == 0 ) then
       median = 0
       return
    end if
    allocate(y,source=x)
    j = 1
    do while( j < n / 3 )
       j = 3*j + 1
    end do
    do while( j >= 1 )
       do i = j + 1, n
          t = y(i)
          k = i
          do while( k > j )
             if( y(k-j) <= t ) exit
             y(k) = y(k-j)
             k = k - j
          end do
          y(k) = t
       end do
       j = j / 3
    end do
    if( mod(n,2) == 1 ) then
       median = y(n/2 + 1)
    else
       median = (y(n/2) + y(n/2 + 1)) / 2
    end if

  end function median

end module psfempir
