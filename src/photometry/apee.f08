!
!  The enclosed energy of apertures (aperture corrections)
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
!  The enclosed energy EE of an aperture is the fraction of the total
!  flux of a star that the aperture contains; APCOUNT / EE is the total
!  flux. Two sources:
!
!   * apee_psf: a known PSF (the field-dependent PSF file of artificial
!     and psfphot, or an analytic Gauss or Moffat profile). The PSF at
!     the field position of the star is pixelised at its sub-pixel
!     position and summed with the partial-pixel weights of the aperture
!     photometry (daoapee), so EE is exactly the fraction of the model
!     counts the aperture sum contains. The error is zero.
!
!   * apee_auto: the frame's own stars. The apertures of every star are
!     its growth curve.
!     - Stars: all the apertures inside the sky ring valid (unsaturated,
!       inside the frame), no comparable neighbour (brighter than 1% in
!       the aperture of FWHM/2) within 5.5 FWHM and the plateau radius
!       + 1.5 FWHM (fainter detections within 2 FWHM are ignored: find
!       splits stars with extended wings; a brighter one excludes the
!       star), S/N > 30 in the aperture of 2 FWHM (at least the 30
!       brightest ones).
!     - The sky bias (a diagnostic, not removed): the median slope of
!       counts by area over the apertures larger than 4 FWHM (of stars
!       isolated by the largest aperture); zero within its noise for
!       the mean sky of daophotsb (BUGS P23).
!     - The total flux: the median of three apertures of the plateau of
!       the stacked growth curve, from the first aperture larger than
!       1.3 FWHM within 0.5% of the level of the curve (the median of
!       the apertures larger than 4 FWHM, or its maximum).
!     - EE = counts / total of every aperture is fitted by a polynomial
!       of the order 0, 1 or 2 in x,y, or in the distance from the field
!       centre (radial), weighted least squares, the excess scatter of
!       the stars over their errors added to the weights, outliers
!       clipped at 3 sigma. The order is lowered for too few stars
!       (3 per parameter) and while the higher terms are not significant
!       (chi2 at 99.9%); the model is not extrapolated beyond the stars.
!     - The error of EE of a star is the error of the polynomial at its
!       position (the covariance) and the excess scatter.
!

module aperture_ee

  use iso_fortran_env

  implicit none
  private

  integer, parameter, private :: dbl = REAL64

  ! the description of the enclosed energy, saved to the header
  type, public :: EEDescr
     character(len=80) :: source = ''      ! PSF or AUTO
     character(len=80) :: psftype = ''     ! FILE, GAUSS, MOFFAT
     character(len=4096) :: psffile = ''
     real(dbl) :: fwhm = -1, beta = -1, xcen = 0, ycen = 0, scale = -1
     integer :: order = -1, ordreq = -1, nstars = 0
     ! the variation by the distance from the field centre xcen, ycen
     logical :: radial = .false.
     real(dbl) :: rref1 = -1, rref2 = -1, skybias = 0
     real(dbl), dimension(:), allocatable :: sigma
     integer, dimension(:), allocatable :: korder
  end type EEDescr

  ! the isolation of stars for the growth curve [FWHM], their plateau
  ! starts at [FWHM], the plateau tolerance, the clipping [sigma]
  real(dbl), parameter, private :: iso_fwhm = 5.5_dbl, plat_fwhm = 1.3_dbl, &
       plat_tol = 0.005_dbl, clip = 3
  ! the sky bias is measured by apertures larger than [FWHM]
  real(dbl), parameter, private :: sky_fwhm = 4
  ! the bright stars: S/N above, at least
  real(dbl), parameter, private :: snr_min = 30
  integer, parameter, private :: nmin_bright = 30

  public :: apee_psf, apee_auto, apee_header

contains

  subroutine apee_psf(desc,stamps,rcoo,acoo,zoom,nbox,xcens,ycens,raper, &
       ecc,incl,ee,eeerr)

    ! EE of the stars at xcens,ycens by the PSF (desc%psftype FILE:
    ! stamps, rcoo, acoo, zoom, nbox and the field geometry desc%xcen,
    ! desc%ycen, desc%scale; GAUSS, MOFFAT: desc%fwhm, desc%beta)

    use psfmodel
    use daofotometr

    type(EEDescr), intent(in) :: desc
    real(dbl), dimension(:,:,:), intent(in) :: stamps
    real(dbl), dimension(:), intent(in) :: rcoo, acoo
    integer, intent(in) :: zoom, nbox
    real, dimension(:), intent(in) :: xcens, ycens, raper
    real, intent(in) :: ecc, incl
    real(dbl), dimension(:,:), intent(out) :: ee, eeerr

    real(dbl), dimension(:,:), allocatable :: sub, p
    real(dbl) :: u, v, width
    integer :: n, l, i0, j0, nsub, z, ptype

    ! the pixels to the largest aperture
    n = ceiling(maxval(raper)) + 2
    allocate(p(-n:n,-n:n))

    if( desc%psftype == 'FILE' ) then
       nsub = nbox
       z = zoom
       allocate(sub(-nsub:nsub,-nsub:nsub))
    else
       z = PSF_DEFZOOM
       nsub = z*(n + 2)
       allocate(sub(-nsub:nsub,-nsub:nsub))
       if( desc%psftype == 'GAUSS' ) then
          ptype = PSF_GAUSS
          width = desc%fwhm / (2*sqrt(2*log(2.0_dbl)))
       else
          ptype = PSF_MOFFAT
          width = psf_moffat_alpha(desc%fwhm,desc%beta)
       end if
       call psf_analytic_sample(ptype,width,desc%beta,1.0_dbl,0.0_dbl,z, &
            nsub,sub)
    end if

    do l = 1, size(xcens)
       if( desc%psftype == 'FILE' ) &
            call psf_field_interpolate(stamps,rcoo,acoo, &
            xcens(l) - desc%xcen,ycens(l) - desc%ycen,desc%scale,sub)
       i0 = nint(xcens(l))
       j0 = nint(ycens(l))
       u = i0 - real(xcens(l),dbl)
       v = j0 - real(ycens(l),dbl)
       call psf_pixelize(sub,nsub,z,u,v,n,p)
       call daoapee(n,p,i0,j0,xcens(l),ycens(l),raper,ecc,incl,ee(l,:))
    end do
    eeerr = 0

  end subroutine apee_psf


  subroutine apee_auto(desc,xcens,ycens,apcts,apcts_err,raper,ring,fwhm, &
       gain,naxes,verbose,ee,eeerr,ok)

    ! EE of all stars by the growth curves of the bright isolated ones,
    ! the requested order of the spatial variation is desc%ordreq

    use exscatter

    type(EEDescr), intent(in out) :: desc
    real, dimension(:), intent(in) :: xcens, ycens, raper, ring
    real, dimension(:,:), intent(in) :: apcts, apcts_err
    real, intent(in) :: fwhm
    real(dbl), intent(in) :: gain
    integer, dimension(:), intent(in) :: naxes
    logical, intent(in) :: verbose
    real(dbl), dimension(:,:), intent(out) :: ee, eeerr
    logical, intent(out) :: ok

    real(dbl), dimension(:,:), allocatable :: cts, var
    real(dbl), dimension(:), allocatable :: dnn, snr, g, tot, vtot, f, s, &
         res, w, xs, ys, cfit, xa, ya
    real(dbl), dimension(:,:), allocatable :: cov
    integer, dimension(:), allocatable :: idx, plat
    logical, dimension(:), allocatable :: cand, sel, keep, prev
    real(dbl) :: xc, yc, hx, hy, sig, gmax, chi2
    integer :: nst, nap, klim, kn, ks, k, l, m, nsel, order, npar, kmax, it, &
         ko, mo, ml, nl
    logical :: fitok, biased

    ok = .false.
    nst = size(xcens)
    nap = size(raper)
    ee = -1
    eeerr = -1
    allocate(desc%sigma(nap),desc%korder(nap))
    desc%sigma = -1
    desc%korder = -1
    desc%order = -1
    desc%nstars = 0
    if( nst < 3 .or. .not. (fwhm > 0) ) then
       write(error_unit,*) 'Warning: aperture correction: no stars or FWHM.'
       return
    end if

    ! counts and their variances (aphot's star term is counts/gain)
    allocate(cts(nst,nap),var(nst,nap))
    cts = apcts
    var = -1
    where( apcts_err > 0 )
       var = real(apcts_err,dbl)**2
    end where

    ! the apertures inside of the sky ring, the normalising aperture
    ! nearest to 2 FWHM
    klim = 0
    do k = 1, nap
       if( raper(k) <= ring(1) ) klim = k
    end do
    kn = minloc(abs(raper(1:max(klim,1)) - 2*fwhm),1)
    if( klim < 2 .or. raper(klim) < plat_fwhm*fwhm ) then
       write(error_unit,'(a,f0.1,a)') ' Warning: aperture correction: no aperture '//&
            'larger than ',plat_fwhm,' FWHM (the plateau of growth curves).'
       return
    end if

    ! The nearest comparable neighbours: brighter than 1% of the star
    ! (in the aperture nearest to FWHM/2). The fainter detections closer
    ! than 2 FWHM are not counted: they are mostly the secondary
    ! detections of the star itself (find splits stars with extended
    ! wings), the star is not used when the brighter one is so close.
    ! Saturated neighbours make the apertures of the star invalid.
    ks = minloc(abs(raper(1:klim) - fwhm/2),1)
    allocate(dnn(nst))
    dnn = huge(dnn)
    do l = 1, nst
       do m = l + 1, nst
          sig = hypot(real(xcens(l) - xcens(m),dbl),real(ycens(l) - ycens(m),dbl))
          if( sig > 40*fwhm ) cycle
          if( comparable(m,l,sig) ) dnn(l) = min(dnn(l),sig)
          if( comparable(l,m,sig) ) dnn(m) = min(dnn(m),sig)
       end do
    end do

    ! all the apertures inside the sky ring valid (unsaturated, inside
    ! the frame)
    allocate(cand(nst),sel(nst),snr(nst))
    do l = 1, nst
       cand(l) = all(cts(l,1:klim) > 0 .and. var(l,1:klim) > 0)
    end do
    where( cand )
       snr = cts(:,kn) / sqrt(var(:,kn))
    elsewhere
       snr = -1
    end where

    ! The sky bias (a diagnostic): the median, over the stars isolated by
    ! the largest aperture, of the slope of counts by area over the
    ! apertures larger than 4 FWHM (where the growth curve is flat). The
    ! sky of daophotsb is the mean of the ring, so the slope is zero
    ! within its noise (the robust mean used before was biased on a sky
    ! of a few photons per pixel, BUGS P23, and the bias was removed
    ! here). It is not removed: the slope includes the light of extended
    ! PSF wings too. At least three such apertures are needed (and the
    ! level of the growth curve is their median), the bias is zero
    ! otherwise.
    desc%skybias = 0
    biased = .false.
    block
      real(dbl), dimension(:), allocatable :: area, slopes, a
      integer, dimension(:), allocatable :: kout
      kout = pack([(k, k = 1, klim)],[(raper(k) >= sky_fwhm*fwhm, k = 1, klim)])
      if( size(kout) >= 3 ) then
         allocate(area(nap))
         area = 4*atan(1.0_dbl)*raper**2
         a = area(kout) - sum(area(kout))/size(kout)
         allocate(slopes(0))
         do l = 1, nst
            if( cand(l) .and. dnn(l) > raper(klim) + 2*fwhm ) &
                 slopes = [slopes, dot_product(cts(l,kout),a) / dot_product(a,a)]
         end do
         if( size(slopes) >= 5 ) then
            biased = .true.
            desc%skybias = median(slopes)
         end if
      end if
    end block

    ! the stacked growth curve of the bright isolated stars, its plateau
    call bright(cand .and. dnn > iso_fwhm*fwhm,sel,nsel)
    if( nsel < 3 ) then
       write(error_unit,'(a,i0,a)') ' Warning: aperture correction: ',nsel, &
            ' bright isolated stars (3 at least).'
       return
    end if
    allocate(g(klim))
    do k = 1, klim
       g(k) = median(pack(cts(:,k) / cts(:,kn),sel))
    end do
    ! the level of the curve: the median of the apertures larger than
    ! 4 FWHM (flat), or its maximum;
    ! the plateau: three apertures from the first one within 0.5%
    ! of the level
    if( biased ) then
       gmax = median(pack(g,raper(1:klim) >= sky_fwhm*fwhm))
    else
       gmax = maxval(g,mask=raper(1:klim) >= plat_fwhm*fwhm)
    end if
    kmax = klim
    do k = 1, klim
       if( raper(k) >= plat_fwhm*fwhm .and. g(k) >= gmax*(1 - plat_tol) ) then
          kmax = k
          exit
       end if
    end do
    plat = [(k, k = kmax, min(kmax+2,klim))]
    desc%rref1 = raper(plat(1))
    desc%rref2 = raper(plat(size(plat)))

    ! the stars isolated also by the plateau apertures
    call bright(cand .and. dnn > max(iso_fwhm*fwhm,desc%rref2 + 1.5*fwhm), &
         sel,nsel)
    if( nsel < 3 ) then
       write(error_unit,'(a,i0,a)') ' Warning: aperture correction: ',nsel, &
            ' bright isolated stars (3 at least).'
       return
    end if
    desc%nstars = nsel
    idx = pack([(l, l = 1, nst)],sel)

    ! the total flux: the median of the plateau apertures
    allocate(tot(nsel),vtot(nsel),xs(nsel),ys(nsel))
    do m = 1, nsel
       l = idx(m)
       tot(m) = median(cts(l,plat))
       vtot(m) = median(var(l,plat))
    end do

    ! the normalised field coordinates: x,y in -1 .. 1, or the distance
    ! from the field centre in 0 .. 1 (the half-diagonal)
    allocate(xa(nst),ya(nst))
    if( desc%radial ) then
       xc = desc%xcen
       yc = desc%ycen
       hx = max(hypot(naxes(1) / 2.0_dbl,naxes(2) / 2.0_dbl),1.0_dbl)
       xa = hypot(xcens - xc,ycens - yc) / hx
       ya = 0
    else
       xc = naxes(1) / 2.0_dbl
       yc = naxes(2) / 2.0_dbl
       hx = max(naxes(1) / 2.0_dbl,1.0_dbl)
       hy = max(naxes(2) / 2.0_dbl,1.0_dbl)
       xa = (xcens - xc) / hx
       ya = (ycens - yc) / hy
    end if
    xs = xa(idx)
    ys = ya(idx)
    ! the model is not extrapolated: the stars outside of the range
    ! of the fitted ones get the values at its border
    xa = min(max(xa,minval(xs)),maxval(xs))
    ya = min(max(ya,minval(ys)),maxval(ys))

    ! the order of the variation, lowered for too few stars
    order = max(min(desc%ordreq,2),0)
    do while( order > 0 .and. nsel < 3*nterms(model(order)) + 2 )
       order = order - 1
    end do
    if( order < desc%ordreq ) write(error_unit,'(a,i0,2a,i0,a,i0,a)') &
         ' Warning: aperture correction: ',nsel,' stars, the order of ', &
         'the variation lowered from ',desc%ordreq,' to ',order,'.'
    desc%order = order
    order = model(order)
    npar = nterms(order)

    allocate(f(nsel),s(nsel),res(nsel),w(nsel),keep(nsel),cfit(npar), &
         cov(npar,npar),prev(nsel))
    do k = 1, nap
       f = cts(idx,k) / tot
       s = sqrt(max(var(idx,k),0.0_dbl) / tot**2 + f**2*vtot / tot**2)
       where( .not. (cts(idx,k) > 0 .and. var(idx,k) > 0) ) s = -1
       keep = s > 0
       if( count(keep) < npar + 1 ) cycle

       sig = 0
       do it = 1, 10
          where( keep )
             w = 1 / (s**2 + sig**2)
          elsewhere
             w = 0
          end where
          call polyfit(order,xs,ys,f,w,cfit,cov,fitok)
          if( .not. fitok ) exit
          res = f - polyval(order,xs,ys,cfit)
          prev = keep
          call excess_scatter(res,s,sig,keep,clip)
          if( count(keep) < npar + 1 ) then
             fitok = .false.
             exit
          end if
          if( all(keep .eqv. prev) ) exit
       end do
       if( .not. fitok ) cycle
       where( keep )
          w = 1 / (s**2 + sig**2)
       elsewhere
          w = 0
       end where
       call polyfit(order,xs,ys,f,w,cfit,cov,fitok)
       if( .not. fitok ) cycle
       res = f - polyval(order,xs,ys,cfit)
       chi2 = sum(w*res**2)

       ! The variation must be significant: the order is lowered while
       ! the higher terms do not decrease chi2 at the 99.9% level, so
       ! the noise of the model is not added when EE does not vary.
       ko = desc%order
       mo = order
       do while( ko > 0 )
          ml = model(ko-1)
          nl = nterms(ml)
          call polyfit(ml,xs,ys,f,w,cfit(1:nl),cov(1:nl,1:nl),fitok)
          if( .not. fitok ) exit
          res = f - polyval(ml,xs,ys,cfit(1:nl))
          if( sum(w*res**2) - chi2 > chi2crit(nterms(mo) - nl) ) exit
          ko = ko - 1
          mo = ml
          chi2 = sum(w*res**2)
       end do
       ! the final fit of the chosen order, its excess scatter
       nl = nterms(mo)
       call polyfit(mo,xs,ys,f,w,cfit(1:nl),cov(1:nl,1:nl),fitok)
       if( .not. fitok ) cycle
       res = f - polyval(mo,xs,ys,cfit(1:nl))
       if( mo /= order ) then
          call excess_scatter(res,s,sig,prev,clip)
          where( keep )
             w = 1 / (s**2 + sig**2)
          elsewhere
             w = 0
          end where
          call polyfit(mo,xs,ys,f,w,cfit(1:nl),cov(1:nl,1:nl),fitok)
          if( .not. fitok ) cycle
          res = f - polyval(mo,xs,ys,cfit(1:nl))
       end if
       chi2 = sum(w*res**2)
       if( count(keep) > nl ) &
            cov(1:nl,1:nl) = cov(1:nl,1:nl) * max(1.0_dbl,chi2/(count(keep)-nl))
       desc%sigma(k) = sig
       desc%korder(k) = ko

       ee(:,k) = polyval(mo,xa,ya,cfit(1:nl))
       do l = 1, nst
          eeerr(l,k) = sqrt(max(polyvar(mo,xa(l),ya(l),cov(1:nl,1:nl)), &
               0.0_dbl) + sig**2)
       end do

       if( verbose ) write(error_unit,'(a,i0,a,f0.2,a,f0.4,a,i0,a,f0.4,a,i0,a,i0)') &
            ' Aperture correction: aperture ',k,' (',raper(k),' pix): EE(centre) = ', &
            cfit(1),', order ',ko,', excess scatter = ',sig,', stars ', &
            count(keep),'/',nsel
    end do

    ok = any(desc%sigma >= 0)
    if( verbose ) write(error_unit,'(a,i0,a,i0,a,f0.2,a,f0.2,a,g0.4,a)') &
         ' Aperture correction: ',nsel,' stars, order ',desc%order, &
         ', total flux by apertures ',desc%rref1,' .. ',desc%rref2, &
         ' pix, sky bias ',desc%skybias,' ADU/pix.'

  contains

    integer function model(order)
      ! the polynomial model of the order (see nterms)
      integer, intent(in) :: order
      if( desc%radial .and. order > 0 ) then
         model = -order
      else
         model = order
      end if
    end function model

    logical function comparable(m,l,d)

      ! is the star m a comparable neighbour of the star l at the distance d?

      integer, intent(in) :: m, l
      real(dbl), intent(in) :: d

      if( d < 2*fwhm ) then
         comparable = cts(m,ks) > cts(l,ks)
      else
         comparable = cts(m,ks) > 0.01*max(cts(l,ks),0.0_dbl)
      end if

    end function comparable

    subroutine bright(mask,sel,nsel)

      ! the stars in mask with S/N > 30 (in the aperture of 2 FWHM),
      ! or the brightest 30 ones

      logical, dimension(:), intent(in) :: mask
      logical, dimension(:), intent(out) :: sel
      integer, intent(out) :: nsel
      real(dbl), dimension(:), allocatable :: a
      real(dbl) :: lim
      integer :: n

      sel = mask .and. snr > snr_min
      nsel = count(sel)
      n = count(mask .and. snr > 0)
      if( nsel < nmin_bright .and. n > 0 ) then
         a = pack(snr,mask .and. snr > 0)
         call sortd(a)
         lim = a(max(n - nmin_bright + 1,1))
         sel = mask .and. snr >= lim .and. snr > 0
         nsel = count(sel)
      end if

    end subroutine bright

  end subroutine apee_auto


  real(dbl) function chi2crit(n)

    ! the 99.9% quantile of the chi2 distribution of n degrees of freedom

    integer, intent(in) :: n
    real(dbl), dimension(5), parameter :: q999 = &
         [10.828_dbl, 13.816_dbl, 16.266_dbl, 18.467_dbl, 20.515_dbl]

    chi2crit = q999(max(min(n,5),1))

  end function chi2crit

  ! The polynomial models (order below): 0, 1, 2 the polynomials
  ! of the order in x,y; -1, -2 the polynomials of the order 1, 2
  ! in the distance (x) from the field centre.

  integer function nterms(order)
    integer, intent(in) :: order
    if( order >= 0 ) then
       nterms = (order + 1)*(order + 2) / 2
    else
       nterms = 1 - order
    end if
  end function nterms

  pure subroutine terms(order,x,y,t)
    integer, intent(in) :: order
    real(dbl), intent(in) :: x, y
    real(dbl), dimension(:), intent(out) :: t
    t(1) = 1
    if( order < 0 ) then
       t(2) = x
       if( order <= -2 ) t(3) = x**2
       return
    end if
    if( order >= 1 ) then
       t(2) = x
       t(3) = y
    end if
    if( order >= 2 ) then
       t(4) = x**2
       t(5) = x*y
       t(6) = y**2
    end if
  end subroutine terms

  function polyval(order,x,y,c) result(p)
    integer, intent(in) :: order
    real(dbl), dimension(:), intent(in) :: x, y, c
    real(dbl), dimension(size(x)) :: p
    real(dbl), dimension(size(c)) :: t
    integer :: i
    do i = 1, size(x)
       call terms(order,x(i),y(i),t)
       p(i) = dot_product(c,t)
    end do
  end function polyval

  real(dbl) function polyvar(order,x,y,cov)
    integer, intent(in) :: order
    real(dbl), intent(in) :: x, y
    real(dbl), dimension(:,:), intent(in) :: cov
    real(dbl), dimension(size(cov,1)) :: t
    call terms(order,x,y,t)
    polyvar = dot_product(t,matmul(cov,t))
  end function polyvar

  subroutine polyfit(order,x,y,f,w,c,cov,ok)

    ! weighted least squares, the covariance is the inverse of the
    ! normal matrix

    integer, intent(in) :: order
    real(dbl), dimension(:), intent(in) :: x, y, f, w
    real(dbl), dimension(:), intent(out) :: c
    real(dbl), dimension(:,:), intent(out) :: cov
    logical, intent(out) :: ok

    real(dbl), dimension(size(c),size(c)) :: a
    real(dbl), dimension(size(c)) :: b, t
    integer :: i, j, n

    n = size(c)
    a = 0
    b = 0
    do i = 1, size(x)
       if( .not. (w(i) > 0) ) cycle
       call terms(order,x(i),y(i),t)
       do j = 1, n
          a(:,j) = a(:,j) + w(i)*t*t(j)
       end do
       b = b + w(i)*t*f(i)
    end do
    call invert(a,cov,ok)
    if( ok ) c = matmul(cov,b)

  end subroutine polyfit

  subroutine invert(a,ainv,ok)

    ! Gauss-Jordan inversion with partial pivoting (small matrices)

    real(dbl), dimension(:,:), intent(in) :: a
    real(dbl), dimension(:,:), intent(out) :: ainv
    logical, intent(out) :: ok

    real(dbl), dimension(size(a,1),2*size(a,1)) :: m
    real(dbl), dimension(2*size(a,1)) :: row
    real(dbl) :: scale
    integer :: n, i, p

    n = size(a,1)
    m = 0
    m(:,1:n) = a
    do i = 1, n
       m(i,n+i) = 1
    end do
    scale = maxval(abs(a))
    ok = scale > 0
    if( .not. ok ) return
    do i = 1, n
       p = maxloc(abs(m(i:n,i)),1) + i - 1
       if( abs(m(p,i)) < 1e3*epsilon(scale)*scale ) then
          ok = .false.
          return
       end if
       if( p /= i ) then
          row = m(i,:)
          m(i,:) = m(p,:)
          m(p,:) = row
       end if
       m(i,:) = m(i,:) / m(i,i)
       do p = 1, n
          if( p /= i ) m(p,:) = m(p,:) - m(p,i)*m(i,:)
       end do
    end do
    ainv = m(:,n+1:)

  end subroutine invert

  real(dbl) function median(a)
    real(dbl), dimension(:), intent(in) :: a
    real(dbl), dimension(:), allocatable :: b
    integer :: n
    n = size(a)
    if( n == 0 ) then
       median = 0
       return
    end if
    b = a
    call sortd(b)
    if( mod(n,2) == 1 ) then
       median = b(n/2+1)
    else
       median = (b(n/2) + b(n/2+1)) / 2
    end if
  end function median

  subroutine sortd(a)

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

  end subroutine sortd


  subroutine apee_header(filename,desc,status)

    ! the description of EE to the header of APERPHOT

    use titsio

    character(len=*), intent(in) :: filename
    type(EEDescr), intent(in) :: desc
    integer, intent(in out) :: status

    type(fitsfiles) :: fits
    character(len=FLEN_KEYWORD) :: key
    character(len=FLEN_VALUE) :: val
    integer :: k, l

    if( status /= 0 ) return

    call fits_open_file(fits,filename,FITS_READWRITE,status)
    if( status /= 0 ) return
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,APEREXTNAME,0,status)

    call fits_update_key(fits,'EESOURCE',desc%source, &
         'enclosed energy APEE by: PSF or AUTO (frame stars)',status)
    if( desc%source == 'PSF' ) then
       call fits_update_key(fits,'EEPSF',desc%psftype, &
            'PSF of EE: FILE, GAUSS, MOFFAT',status)
       if( desc%psftype == 'FILE' ) then
          ! the last characters of a long (path) name fit to the card
          l = len_trim(desc%psffile)
          if( l > 64 ) then
             val = '...'//desc%psffile(l-60:l)
          else
             val = desc%psffile
          end if
          call fits_update_key(fits,'EEPSFFIL',val,'PSF file of EE',status)
          call fits_update_key(fits,'EEPSFCX',desc%xcen,-7, &
               '[pix] field centre of the PSF file, x',status)
          call fits_update_key(fits,'EEPSFCY',desc%ycen,-7, &
               '[pix] field centre of the PSF file, y',status)
          call fits_update_key(fits,'EEPSFSCL',desc%scale,-10, &
               '[pix/deg] scale of RADCOO of the PSF file',status)
       else
          call fits_update_key(fits,'EEFWHM',desc%fwhm,-6, &
               '[pix] FWHM of the PSF of EE',status)
          if( desc%psftype == 'MOFFAT' ) &
               call fits_update_key(fits,'EEBETA',desc%beta,-5, &
               'Moffat beta of the PSF of EE',status)
       end if
    else
       call fits_update_key(fits,'EEVARREQ',desc%ordreq, &
            'requested order of the spatial variation of EE',status)
       call fits_update_key(fits,'EEVAR',desc%order, &
            'order of the spatial variation of EE (polynomial)',status)
       if( desc%radial ) then
          val = 'RADIAL'
       else
          val = 'XY'
       end if
       call fits_update_key(fits,'EEVARCOO',val, &
            'EE variation by: XY, RADIAL (from the field centre)',status)
       if( desc%radial ) then
          call fits_update_key(fits,'EEFCX',desc%xcen,-7, &
               '[pix] field centre for EE, x',status)
          call fits_update_key(fits,'EEFCY',desc%ycen,-7, &
               '[pix] field centre for EE, y',status)
       end if
       call fits_update_key(fits,'EENSTAR',desc%nstars, &
            'bright isolated stars for EE',status)
       call fits_update_key(fits,'EERTOT1',desc%rref1,-5, &
            '[pix] total flux: the smallest plateau aperture',status)
       call fits_update_key(fits,'EERTOT2',desc%rref2,-5, &
            '[pix] total flux: the largest plateau aperture',status)
       call fits_update_key(fits,'EESKYBIA',desc%skybias,-4, &
            '[ADU/pix] residual sky bias (diagnostic, not removed)',status)
       do k = 1, size(desc%sigma)
          call fits_make_keyn('EEORD',k,key,status)
          call fits_update_key(fits,key,desc%korder(k), &
               'order of the variation of APEE (significant)',status)
          call fits_make_keyn('EESIG',k,key,status)
          call fits_update_key(fits,key,desc%sigma(k),-3, &
               'scatter of EE of stars over the model',status)
       end do
    end if

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine apee_header

end module aperture_ee
