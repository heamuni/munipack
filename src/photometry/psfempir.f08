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
!  field coordinates (DAOPHOT VARIABLE = 0, 1, 2).
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
!  at 5 sigma. A light regularisation by the squared Laplacian keeps the
!  tables smooth where the stars sample the grid sparsely, and pulls
!  them to zero at the edge of the table.
!
!  The normal equations are solved by the conjugate gradients with
!  the Jacobi preconditioner; the system is never formed.
!
!  The residual tables carry no flux and do not move the centroid:
!  every table is corrected by a combination of the base at the field
!  centre and its derivatives so that its sum and its first moments are
!  zero. The flux of the model is that of the base (the analytic core
!  is normalised to the unit flux over the plane, the light beyond the
!  table is in its wings), and the centre of the PSF is the centroid
!  of the base.
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

  public :: empirical_build

contains

  subroutine empirical_build(frame,psf,stars,rt,clip,smooth,chi2,resrms,niter)

    ! Builds the residual tables psf%res on the stars (fitted, their
    ! models in frame%mstar when modelled), using the pixels in the
    ! radius rt [pix]. The current tables are the initial solution and,
    ! with clip, the residuals by them clip outliers. chi2 are the
    ! reduced chi2 of stars by the new PSF, resrms is the rms of the
    ! constant table relative to the peak of the base.

    type(FrameType), intent(in) :: frame
    type(PsfType), intent(in out) :: psf
    type(StarFit), dimension(:), intent(in) :: stars
    real(REAL64), intent(in) :: rt, smooth
    logical, intent(in) :: clip
    real(REAL64), dimension(:), intent(out) :: chi2
    real(REAL64), intent(out) :: resrms
    integer, intent(out) :: niter

    integer, parameter :: maxiter = 1000
    real(REAL64), parameter :: tol = 1e-4_REAL64

    type(EmpStar), dimension(:), allocatable :: es
    real(REAL64), dimension(:,:,:), allocatable :: x, r, z, p, q, b, diag
    real(REAL64), dimension(:,:), allocatable :: sub, pold, pbase, reff, pm, cb, &
         supp
    real(REAL64), dimension(:), allocatable :: lam, eta, zs
    real(REAL64), dimension(3,3) :: cm
    real(REAL64), dimension(3) :: mom, coef
    real(REAL64) :: rz, rznew, alpha, bnorm, rnorm, s, h, sig
    integer :: ns, nt, ntab, zoom, npoly, k, i, j, l, t, n, nsub, iter, i0, j0
    logical :: ok

    ns = size(stars)
    ntab = psf%ntab
    zoom = psf%zoom
    npoly = psf%npoly

    ! the support of the tables: the circle of the radius rt, the tables
    ! are zero outside (no data constrain them there)
    allocate(supp(-ntab:ntab,-ntab:ntab))
    do j = -ntab, ntab
       do i = -ntab, ntab
          if( i**2 + j**2 <= (zoom*rt)**2 ) then
             supp(i,j) = 1
          else
             supp(i,j) = 0
          end if
       end do
    end do

    ! the data of stars
    allocate(es(ns))
    nt = ceiling(rt) + 1
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
               if( (i0 - st%x)**2 + (j0 - st%y)**2 > rt**2 ) cycle
               if( i0 < 1 .or. j0 < 1 .or. i0 > size(frame%data,1) .or. &
                    j0 > size(frame%data,2) ) cycle
               if( .not. frame%valid(i0,j0) ) cycle
               s = frame%data(i0,j0) - frame%mstar(i0,j0) - st%sky
               if( abs(i) <= n .and. abs(j) <= n ) s = s + st%flux*pold(i,j)
               e%t(i,j) = s - st%flux*pbase(i,j)
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
            sig = 1.4826_REAL64*median(abs(zs - median(zs)))
            sig = 5*max(sig,1.0_REAL64)
            where( abs(pm) > sig ) e%w = 0
          end associate
       end do
       deallocate(reff,pm)
    end if

    ! the right side, the diagonal of the normal matrix
    allocate(b(-ntab:ntab,-ntab:ntab,npoly),diag(-ntab:ntab,-ntab:ntab,npoly), &
         r(-ntab:ntab,-ntab:ntab,npoly),z(-ntab:ntab,-ntab:ntab,npoly), &
         p(-ntab:ntab,-ntab:ntab,npoly),q(-ntab:ntab,-ntab:ntab,npoly), &
         reff(-ntab:ntab,-ntab:ntab),pm(-nt:nt,-nt:nt),lam(npoly),eta(npoly))
    b = 0
    diag = 0
    do k = 1, ns
       associate(e => es(k))
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

    ! the regularisation relative to the mean of the diagonal
    do t = 1, npoly
       s = sum(diag(:,:,t),diag(:,:,t) > 0) / max(count(diag(:,:,t) > 0),1)
       lam(t) = smooth*s
       eta(t) = 1e-6_REAL64*s
       diag(:,:,t) = diag(:,:,t) + 20*lam(t) + eta(t)
       b(:,:,t) = supp*b(:,:,t)
    end do

    ! the preconditioned conjugate gradients
    call operator(x,q)
    r = b - q
    z = r / diag
    p = z
    rz = sum(r*z)
    bnorm = sqrt(sum(b**2))
    niter = 0
    do iter = 1, maxiter
       niter = iter
       rnorm = sqrt(sum(r**2))
       if( .not. (rnorm > tol*bnorm) ) exit
       call operator(p,q)
       s = sum(p*q)
       if( .not. (s > 0) ) exit
       alpha = rz / s
       x = x + alpha*p
       r = r - alpha*q
       z = r / diag
       rznew = sum(r*z)
       p = z + (rznew / rz)*p
       rz = rznew
    end do

    ! no flux, no shift: the tables minus a combination of the base at
    ! the field centre and its derivatives have zero sum and moments
    call psf%sampled_base(psf%pxc,psf%pyc,sub,nsub)
    allocate(cb(-ntab:ntab,-ntab:ntab),pold(-ntab:ntab,-ntab:ntab), &
         pbase(-ntab:ntab,-ntab:ntab))
    cb = sub(-ntab:ntab,-ntab:ntab)
    pold = 0
    pbase = 0
    pold(-ntab+1:ntab-1,:) = (cb(-ntab+2:ntab,:) - cb(-ntab:ntab-2,:)) / 2
    pbase(:,-ntab+1:ntab-1) = (cb(:,-ntab+2:ntab) - cb(:,-ntab:ntab-2)) / 2
    cb = supp*cb
    pold = supp*pold
    pbase = supp*pbase
    call moments(cb,cm(:,1))
    call moments(pold,cm(:,2))
    call moments(pbase,cm(:,3))
    ! the model of the constant term is rescaled to the unit flux
    ! (and shifted), the other terms by the same factor
    call moments(x(:,:,1),mom)
    call solve3(cm,mom,coef,ok)
    if( ok ) then
       h = 1 / (1 + coef(1)*sum(cb) / sum(sub))
    else
       h = 1
    end if
    do t = 1, npoly
       call moments(x(:,:,t),mom)
       call solve3(cm,mom,coef,ok)
       if( ok ) x(:,:,t) = h*(x(:,:,t) - coef(1)*cb - coef(2)*pold - coef(3)*pbase)
    end do
    h = maxval(cb)
    resrms = sqrt(sum(x(:,:,1)**2,supp > 0) / count(supp > 0)) / h
    psf%res = x

    ! chi2 of stars by the new tables
    do k = 1, ns
       associate(e => es(k))
         call combine(x,e%pt,reff)
         call forward(e,reff,pm)
         l = count(e%w > 0)
         chi2(k) = sum(e%w*(e%t - e%flux*pm)**2) / max(l - 3,1)
       end associate
    end do

  contains

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

    subroutine operator(v,av)

      ! the normal matrix times v: sum over stars of A'WA, the polynomial
      ! terms, and the regularisation

      real(REAL64), dimension(-ntab:,-ntab:,:), intent(in) :: v
      real(REAL64), dimension(-ntab:,-ntab:,:), intent(out) :: av

      real(REAL64), dimension(-ntab:ntab,-ntab:ntab) :: reff, g, lap
      real(REAL64), dimension(-nt:nt,-nt:nt) :: pm
      integer :: k, t

      ! (the input is zero outside of the support)
      av = 0
      do k = 1, ns
         associate(e => es(k))
           call combine(v,e%pt,reff)
           call forward(e,reff,pm)
           pm = e%flux**2*e%w*pm
           call adjoint(e,pm,g,.false.)
           do t = 1, npoly
              av(:,:,t) = av(:,:,t) + e%pt(t)*g
           end do
         end associate
      end do
      do t = 1, npoly
         call laplace(v(:,:,t),lap)
         call laplace(lap,g)
         av(:,:,t) = supp*(av(:,:,t) + lam(t)*g + eta(t)*v(:,:,t))
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

    subroutine moments(v,mom)

      ! the sum and the first moments of the table

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

  end subroutine empirical_build


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
