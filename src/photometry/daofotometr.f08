!
!     PHOTOMETRY, aperture photometry
!     Copyright (C) 1997-9,2010,2013-19 Filip Hroch, Masaryk University
!     Copyright (C) 1991 P.B. Stetson, Dominon Astrophysical Observatory
!
!
!  This file is part of Munipack.
!
!  Credits
!
!    Almost all this source is authored by P. B. Stetson.
!    I adapted it for Fortran 90 (allocatable arrays, precision,
!    array syntax), corrected errors, improve sky estimate by
!    variance stabilising estimator, elliptic aperture.
!
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
!===================================================================
!
! 2019: I did complete revision of the source. The various obsolete
!      parts has been removed: input by hand, output for every star
!      (but output for program to program has been added). Both
!      elliptic and testing estimators of edge fraction has been included.
!      The code has been improved by latest computer development
!      and my development of robust estimators and photometry in counts.
!      I tried to keep Stetson's programming spirit and conventions.

module daofotometr

  integer, parameter, private :: dbl = selected_real_kind(14)

  private :: circfrac, sqcover, pixgeom, pixfrac

contains

  subroutine  daophotsb (d,derr,xstar,ystar,raper,ring,ecc,incl,lobad,hibad,&
       verbose,plog,apcts,apcts_err,skystar,skystar_err,gain,rn2,bias,rmask)

    use iso_fortran_env
    use skymean

    implicit none

    real, dimension(:,:), intent(in) :: d,derr
    real, dimension(:), intent(in) :: xstar, ystar, raper,ring
    real, dimension(:,:), intent(out) :: apcts,apcts_err
    real, dimension(:), intent(out) :: skystar,skystar_err
    real, intent(in) :: lobad, hibad, ecc, incl
    logical, intent(in) :: verbose, plog
    ! the noise model of the frame for the sky (gain [e-/ADU], read noise
    ! squared [ADU**2], bias [ADU]); without it, the sky is clipped
    ! by the scatter of the data only
    real(dbl), intent(in), optional :: gain, rn2, bias
    ! the radii of the masks of the stars [pix]: the sky ring of a star
    ! leaves out the pixels inside the masks of the other stars
    real, dimension(:), intent(in), optional :: rmask

!
!=======================================================================
!
! This subroutine derives the concentric aperture photometry.  At
! present, this is the only place in all of DAOPHOT where sky values
! are derived for the individual stars.
!
!               OFFICIAL DAO VERSION:  1991 April 18
!
! Argument
!
!    WATCH (INPUT) governs whether information relating to the progress
!          of the reductions is to be typed on the terminal screen
!          during execution.
!
!=======================================================================
!

    integer, parameter  :: minsky = 7!, maxap = 12
    real, parameter :: pi = 3.141592653589793115997963
    real, parameter :: rad = 57.295779513082322865

!
! Parameters:
!
! MINSKY is the smallest number of pixels from which the sky may be
!        determined.  If for some star the number of sky pixels
!        is less than MINSKY, an error code will result and
!        control will return to the main program.
!
! MAXSKY the maximum number of pixels allowed in the sky annulus.
!        This and the user's requested inner sky radius will later
!        determine the maximum permitted outer sky radius.
!
! MAXAP  the maximum number of star apertures allowed.
!
    real(dbl), dimension(size(raper)) :: area, area2, apsum
    real(dbl), dimension(:), allocatable :: sky, skyall
    integer(int16), dimension(:,:), allocatable :: cover
    real, dimension(3) :: error
    real :: skymod, skysig, sigsq, skyvar, skyerr, datum, r, rsq, fractn, &
         edge, apmxsq, rout, dysq, xc, yc, e,c,s,sqrte,x,y,phi,q,rin
    real(dbl) :: g, rn, b, smean, serr, ssig, phpadu
    integer :: i, j, k, n, naper, nstar, lx, ly, nx, ny, mx, my, &
         nsky, nall, ncol, nrow, iflag, nused, ncover
    logical :: masking

    skystar = -1
    skystar_err = -1
    apcts = -1
    apcts_err = -1

    ncol = size(d,1)
    nrow = size(d,2)
    naper = size(raper)
    nstar = size(xstar)

!-----------------------------------------------------------------------
!
! SECTION 1
!
! Ascertain the name of the aperture photometry parameter table, and
! read it in.  Then set up all necessary variables for the forthcoming
! reductions. Finally, identify and open the input and output files.
!
    do i = 2, naper
       if( raper(i) < raper(i-1) ) &
            stop "Error in APHOT: some aperture radius is invalid."
    end do

    apmxsq = -1.0
    do i = 1, naper
       apmxsq = max(apmxsq, (raper(i)+0.5)**2)
    end do

!
! sky buffers allocation
!
    nsky = int(max(pi*((ring(2)+1)**2 - ring(1)**2),1.5))
    allocate(sky(nsky),skyall(nsky))

    ! the noise model
    g = -1
    rn = 0
    b = 0
    if( present(gain) ) g = gain
    if( present(rn2) ) rn = rn2
    if( present(bias) ) b = bias

    ! The photons (electrons) per ADU for the Poisson noise of the star:
    ! the gain of the frame (1 for frames in electrons, BUNIT = COUNTS),
    ! 1 when unknown (the caller warns). The data are in ADU, so the
    ! variance of the counts of a star is counts/gain [ADU**2].
    if( g > 0 ) then
       phpadu = g
    else
       phpadu = 1
    end if

    ! the masks: the number of stars covering a pixel
    masking = present(rmask)
    if( masking ) call starmask(xstar,ystar,rmask,ncol,nrow,cover)

!
! NAPER   is the number of apertures, whose radii are stored in
!         elements 1 through NAPER of the array PAR.
!
! APMXSQ  is the outermost edge of the largest aperture-- if the
!         distance squared of the center of a pixel from the centroid of
!         the star is greater than APMXSQ, then we know that no part
!         of the pixel is to be included in any aperture.
!
! Now define the other variables whose values are in the table.
!

    if( ring(1) >= ring(2) ) &
         stop "Error in APHOT: invalid sky ring radius."

    rin  = ring(1)
    rout = ring(2)

! If progress is being monitored, type out column headers.
!
!  if (watch > 0.5) &
!  if( verbose ) &
!       write (*,"(/13X, 'STAR', 5X, 'X', 7X, 'Y', 9X, 'MAG.(1)', 8X, 'SKY')")
!
!-----------------------------------------------------------------------
!
! SECTION 2
!
! Derive aperture photometry object by object.
!
! Get the coordinates of next object to be measured.
!
    lx = 1
    ly = 1
    nx = ncol
    ny = nrow

    e = ecc
    c = cos(incl/rad)
    s = sin(incl/rad)
    sqrte = sqrt(1 - e**2)

    do n = 1, nstar

       xc = xstar(n)
       yc = ystar(n)
!       istar = n
!
! Compute the limits of the submatrix.
!
       lx = max(1, int(xc-rout)+1)
       mx = min(ncol, int(xc+rout))
       ly = max(1, int(yc-rout)+1)
       my = min(nrow, int(yc+rout))
       edge = min(xc-0.5, (ncol+0.5)-xc, yc-0.5, (nrow+0.5)-yc)
!
! EDGE is the distance of the star's centroid from the outermost
! extremum of the array.
!

!
! Initialize star counts and aperture area.
!
       do  i = 1, naper
          apsum(i) = real(0.0,dbl)
!
! If this star aperture extends outside the array, the magnitude
! in this aperture will be no good.
!
          if (edge < raper(i)) apsum(i) = -huge(apsum)  ! Null magnitude
          area(i) = real(0.0,dbl)
          area2(i) = real(0.0,dbl)
       enddo
!
! Now read through the submatrix, picking out the data we want.
!
       nsky = 0
       nall = 0
       do  j = ly, my
          dysq = (j - yc)**2
          do i = lx,mx
             rsq = dysq + (i - xc)**2
             datum = d(i,j)

             ! elliptic apertures
             call pixgeom(i,j,xc,yc,e,c,s,sqrte,x,y,r,phi,q)

!
! Is this pixel within the sky annulus?
!
!           write(*,*) r,rin,rout,datum,lobad,hibad
             if( rin <= r .and. r <= rout .and. &
                  lobad <= datum .and. datum < hibad ) then
                nall = nall + 1
                skyall(nall) = datum
                ncover = 0
                if( masking ) then
                   ! covered by another star? (the star's own mask is
                   ! subtracted by the same test as in starmask)
                   ncover = cover(i,j)
                   if( (i - xc)**2 + (j - yc)**2 < rmask(n)**2 ) &
                        ncover = ncover - 1
                end if
                if( ncover == 0 ) then
                   nsky = nsky + 1
                   sky(nsky) = datum
                end if
             endif
!
! The inclusion of partial pixels inside the aperture is done as
! follows:  if the distance of the center of the current pixel from the
! centroid of the star [radius vector r(i,j)] is exactly equal to the
! radius of the aperture [R(k)], then one-half of the counts in the
! pixel are included.  If r(i,j) < R(k)-0.5, then the entire pixel is
! included, while if r(i,j) > R(k)+0.5, the pixel is wholly excluded.
! In between, viz. for  R(k)-0.5 < r(i,j) < R(k)+0.5, the fraction of
! the counts included varies linearly.  Therefore a circular aperture
! is approximated by an irregular (not even convex) polygon.
!
! If this pixel falls completely outside the LARGEST aperture, go on
! to the next pixel.  Notice that APMXSQ has actually been defined
! as (R(k)+0.5)**2 for the largest value of R(k), in accordance with
! the formula used for the partial pixels.
!
             if (rsq <= apmxsq) then

                ! DAOPHOT original code adds to 0.5 to r by the commented
                ! formula; Munipack uses unmodified r value, the one-half
                ! is added in the fractn determination.
                ! r = sqrt(rsq) - 0.5
!
                do k = 1, naper
!
! if this pixel falls completely outside THIS aperture, go on to the
! next aperture.
!
                   if ( r <= raper(k)+1 ) then

                    ! determination of fraction of pixel in circular
                    ! apertures uses traditional estimator, while
                    ! the elliptic the improved one
                      fractn = pixfrac(raper(k),e,sqrte,x,y,r,phi,q)

!
! fractn is the fraction of the pixel that falls inside the
! (irregular) aperture.
!
! If the pixel is bad, set the total counts in this aperture to a number
! so negative that it will never be positive again.
!                                                          ! Null magnitude
                      if (datum < lobad .or. datum >= hibad ) apsum(k)=-huge(1.0)
                      apsum(k) = apsum(k) + fractn*datum
                      area(k) = area(k) + fractn
                      area2(k) = area2(k) + fractn**2
                   endif
                enddo
             endif
          enddo ! i
       enddo ! j
!
! We have accumulated the brightnesses of individual sky pixels in the
! one-dimensional array SKY.  Pixels falling above or below the BAD
! limits have already been eliminated.  Now sort SKY to place the
! pixels in order of increasing brightness.
!
       ! too much of the ring masked (a crowded field): all the ring
       if( masking .and. nsky < max(minsky,nall/4) ) then
          nsky = nall
          sky(1:nsky) = skyall(1:nall)
       end if

       if (nsky < minsky)  then
          write(error_unit,*) "There aren't enough pixels in the sky annulus."
          write(error_unit,*) 'Object at coordinates:',xc,yc
          write(error_unit,*) ' Are you sure your bad pixel thresholds are all right?'
          write(error_unit,*) ' If so, then you need a larger outer sky radius.'
          write(error_unit,*) nsky,minsky,size(sky),lobad,hibad
          goto 3333
       end if
!
! The mean of the sky (not the mode: the aperture sums the sky pixels),
! clipped against stars, cosmics and bad pixels (see skymean).
!
       call sky_mean(sky(1:nsky),g,rn,b,smean,serr,ssig,nused)
       skymod = real(smean)
       skyerr = real(serr)
       skysig = real(ssig)
       iflag = 0
       if( nused < minsky ) iflag = 5

       skyvar = skysig**2
       sigsq = skyerr**2    ! equivalent of: sigsq = skyvar/nsky

!
! SKYMOD has units of (ADU/pixel), and SKYSIG is the pixel-to-pixel
! scatter of SKYMOD, in units of (ADU/pixel).  SKYVAR is the
! variance (square of the standard deviation) of the sky brightness,
! (ADU/pixel)**2, and SIGSQ is the square of the standard error of the
! mean sky brightness.
!
! Subtract the sky from the integrated brightnesses in the apertures,
! convert the results to magnitudes, and compute standard errors.
!
       do i = 1, naper
!
! If the modal sky value could not be determined, set the magnitude
! to 99.999: and total count sum to -1.
!
          apsum(i) = apsum(i) - skymod*area(i)
          if( apsum(i) > 0 .and. iflag < 5 ) then
!
! If the star + sky is fainter than the sky, or if the star aperture
! extends beyond the limits of the picture, or if there is a bad pixel
! in the star aperture, set the magnitude to 99.999.
!
             ! the noise of the sky pixels in the aperture: the sum of
             ! the squares of their weights (a partial pixel of the
             ! fraction f adds f**2 of its variance)
             error(1) = real(area2(i)*skyvar)
             error(2) = real(apsum(i)/phpadu)
             error(3) = real(sigsq*area(i)**2)

! For Munipack, the counts are in ADU (the units of the frame); the
! Poisson variance of the star is counts/gain [ADU**2].

             apcts(n,i) = real(apsum(i))
             apcts_err(n,i) = real(sqrt(sum(error)))
!             if( i == 6 ) write(*,*) area(i),apsum(i),sqrt(error)
          else
             apcts(n,i) = -1
             apcts_err(n,i) = -1
          end if


!
! These variables ERRORn are the respective variances (squares of the
! mean errors) for: (1) random noise inside the star aperture, including
! readout noise and the degree of contamination by other stars in the
! neighborhood, as estimated by the scatter in the sky values (this
! standard error increases as the square root of the area of the
! aperture); (2) the Poisson statistics of the observed star brightness
! (counts/gain); (3) the uncertainty of the mean sky brightness (this
! standard error increases directly with the area of the aperture).
!
       enddo
!
! Write out the answers.
!
!      if (watch > 0.5) then
!      if( verbose ) then
!         write (*,"(10X, I5, 2F8.1, F9.3, ' +-', F6.3, 3x, g0.3)") &
!              istar, xc, yc, apmag, magerr, skymod
!         write (*,"(/1X, I5, 14F9.3)") istar, xc, yc, (apmag(i), i=1,naper)
!         write (*,"(4X, F9.3, 2F6.2, F8.3, 11F9.3)") skymod, min(99.99,skysig),&
!              min(999.99, max(-99.99,skyskw)), (magerr(i), i=1,naper)
!      endif

       if( plog ) write(*,'(a,2(f0.3,1x),es15.5)') '=APHOT> ',xc,yc,apcts(n,1)

       skystar(n) = skymod
       skystar_err(n) = skyerr

3333 continue
    enddo ! over stars
!
!-----------------------------------------------------------------------
!
! Normal return.
!
! Estimate magnitude limit, close files, and return.
!

    if( verbose ) write(*,"(a,i0,a)") " Found ",nstar," star(s)."
    deallocate(sky,skyall)

  end subroutine daophotsb

  subroutine daocontam(xstar,ystar,raper,rmask,apcts,apcts_err,ncont)

    ! The apertures contaminated by the light of a brighter neighbour
    ! (by more than three times their error) are invalidated (-1): a
    ! faint star a few pixels from a bright one keeps its light in the
    ! aperture. The light of the neighbour m in the aperture of the
    ! radius a of the star n at the distance d is the light of m in the
    ! annulus d-a .. d+a times the fraction of the annulus the aperture
    ! covers, at least the light of m inside a-d (the aperture over the
    ! centre of m).
    !
    ! The light of m in the annulus is its flux times the mean profile of
    ! the frame: the median growth curve of isolated bright stars (no
    ! detection within the last aperture plus three reference radii),
    ! normalised by their last aperture; the flux of m is its reference
    ! aperture (nearest to 2 pix) divided by the profile there. The
    ! growth curve of m itself cannot be used at the distance of n: its
    ! apertures there contain n, so n's own light was taken for m's and
    ! invalidated the stars beside any brighter neighbour closer than the
    ! largest aperture (28 of 150 on a Gaussian field, pairs 15-20 px).
    ! With fewer than three isolated stars the growth curve of m is used,
    ! the light of n taken out of it (n's growth curve within d/2, spread
    ! over d - r90 .. d + r90 around m, r90 enclosing 90% of it).
    !
    ! The estimate must be significant too (three times its error by the
    ! errors of m's apertures): the growth curves of faint detections are
    ! noise. A neighbour without valid apertures (saturated, at the
    ! border) contaminates the apertures reaching its mask. ncont is the
    ! number of invalidated apertures.

    use oakleaf, only: median

    real, dimension(:), intent(in) :: xstar, ystar, raper, rmask
    real, dimension(:,:), intent(in out) :: apcts, apcts_err
    integer, intent(out) :: ncont

    real, parameter :: kappa = 3
    real, parameter :: snmin = 30
    real, dimension(:,:), allocatable :: g, ge, q
    real, dimension(:), allocatable :: bright, rcurve, gm, ee
    integer, dimension(:), allocatable :: nval
    logical :: profile
    integer :: n, m, k, nstar, naper, kref, np, kt
    real :: d, a, rlo, rhi, c, e, f, reach, tn, r90, tm, tmerr, risol

    nstar = size(xstar)
    naper = size(raper)
    ncont = 0
    if( nstar < 2 .or. naper < 1 ) return

    ! the growth curves (a copy: the flags do not change them), valid
    ! from the first aperture to the last valid one in a row
    allocate(g(nstar,naper),ge(nstar,naper),nval(nstar),bright(nstar), &
         rcurve(nstar),gm(naper),ee(naper))
    g = apcts
    ge = apcts_err
    do n = 1, nstar
       nval(n) = 0
       do k = 1, naper
          if( .not. (g(n,k) > 0 .and. apcts_err(n,k) > 0) ) exit
          nval(n) = k
       end do
    end do

    ! the brightness: the counts in the aperture nearest to 2 pix
    kref = minloc(abs(raper - 2),1)
    do n = 1, nstar
       if( nval(n) >= kref ) then
          bright(n) = g(n,kref)
       else if( nval(n) > 0 ) then
          bright(n) = g(n,nval(n))
       else
          ! no valid aperture: saturated (or at the border)
          bright(n) = huge(1.0)
       end if
       if( nval(n) > 0 ) then
          rcurve(n) = raper(nval(n))
       else
          rcurve(n) = rmask(n)
       end if
    end do
    reach = raper(naper) + max(maxval(rcurve),maxval(rmask))

    ! the mean profile by isolated bright stars
    risol = raper(naper) + 3*raper(kref)
    allocate(q(nstar,naper))
    np = 0
    do n = 1, nstar
       if( nval(n) < naper ) cycle
       if( .not. (g(n,kref) > snmin*ge(n,kref)) ) cycle
       if( any(hypot(xstar - xstar(n),ystar - ystar(n)) < risol .and. &
            (xstar /= xstar(n) .or. ystar /= ystar(n))) ) cycle
       np = np + 1
       q(np,:) = g(n,:) / g(n,naper)
    end do
    profile = np >= 3
    if( profile ) then
       do k = 1, naper
          ee(k) = median(q(1:np,k))
          if( k > 1 ) ee(k) = max(ee(k),ee(k-1))
       end do
       ee = ee / ee(naper)
    end if

    do n = 1, nstar
       if( nval(n) == 0 ) cycle
       do m = 1, nstar
          if( m == n .or. .not. (bright(m) > bright(n)) ) cycle
          if( abs(xstar(m) - xstar(n)) > reach .or. &
               abs(ystar(m) - ystar(n)) > reach ) cycle
          d = hypot(xstar(m) - xstar(n),ystar(m) - ystar(n))
          if( nval(m) > 0 ) then
             if( profile ) then
                ! the flux of m by its reference aperture
                kt = min(kref,nval(m))
                tm = g(m,kt) / ee(kt)
                tmerr = ge(m,kt) / ee(kt)
             else
                ! the growth curve of m without the light of n
                tn = max(curve(g(n,:),nval(n),min(d/2,rcurve(n))),0.0)
                r90 = max(rfrac(g(n,:),nval(n),0.9*tn),raper(1))
                do k = 1, nval(m)
                   gm(k) = g(m,k) - tn * &
                        min(max((raper(k) - d) / (2*r90) + 0.5,0.0),1.0)
                end do
             end if
          end if
          do k = 1, naper
             if( .not. (apcts(n,k) > 0 .and. apcts_err(n,k) > 0) ) cycle
             a = raper(k)
             if( nval(m) == 0 ) then
                ! a saturated neighbour: its mask
                if( d - a < rmask(m) ) then
                   apcts(n,k) = -1
                   apcts_err(n,k) = -1
                   ncont = ncont + 1
                end if
                cycle
             end if
             rlo = max(d - a,0.0)
             rhi = d + a
             f = min(a**2 / (rhi**2 - rlo**2),1.0)
             if( profile ) then
                if( rlo >= raper(naper) ) cycle
                c = tm * (curve(ee,naper,rhi) - curve(ee,naper,rlo)) * f
                if( a > d ) c = max(c,tm * curve(ee,naper,a - d))
                e = tmerr * c / tm
             else
                if( rlo >= rcurve(m) ) cycle
                c = (curve(gm,nval(m),rhi) - curve(gm,nval(m),rlo)) * f
                e = curverr(m,rhi) * f
                if( a > d ) then
                   c = max(c,curve(gm,nval(m),a - d))
                   e = max(e,curverr(m,a - d))
                end if
             end if
             if( c > kappa*apcts_err(n,k) .and. c > kappa*e ) then
                apcts(n,k) = -1
                apcts_err(n,k) = -1
                ncont = ncont + 1
             end if
          end do
       end do
    end do

  contains

    real function curve(gc,nv,r)

      ! the growth curve gc (its first nv apertures) at the radius r:
      ! linear in r**2 between the apertures, from zero at the centre,
      ! constant beyond the last valid aperture

      real, dimension(:), intent(in) :: gc
      integer, intent(in) :: nv
      real, intent(in) :: r
      integer :: k
      real :: r1, g1

      curve = gc(nv)
      r1 = 0
      g1 = 0
      do k = 1, nv
         if( r <= raper(k) ) then
            curve = g1 + (gc(k) - g1)*(r**2 - r1**2) / (raper(k)**2 - r1**2)
            return
         end if
         r1 = raper(k)
         g1 = gc(k)
      end do

    end function curve

    real function rfrac(gc,nv,t)

      ! the radius where the growth curve gc reaches t (the inverse of
      ! curve; the last valid aperture when it does not)

      real, dimension(:), intent(in) :: gc
      integer, intent(in) :: nv
      real, intent(in) :: t
      integer :: k
      real :: r1, g1

      rfrac = raper(nv)
      r1 = 0
      g1 = 0
      do k = 1, nv
         if( t <= gc(k) ) then
            if( gc(k) > g1 ) then
               rfrac = sqrt(r1**2 + (t - g1)*(raper(k)**2 - r1**2) / (gc(k) - g1))
            else
               rfrac = raper(k)
            end if
            return
         end if
         r1 = raper(k)
         g1 = gc(k)
      end do

    end function rfrac

    real function curverr(m,r)

      ! the error of the growth curve of the star m at the radius r:
      ! the error of the aperture enclosing r (the last valid one beyond)

      integer, intent(in) :: m
      real, intent(in) :: r
      integer :: k

      curverr = ge(m,nval(m))
      do k = 1, nval(m)
         if( r <= raper(k) ) then
            curverr = ge(m,k)
            return
         end if
      end do

    end function curverr

  end subroutine daocontam

  subroutine starmask(xstar,ystar,rmask,ncol,nrow,cover)

    ! the number of the stars whose masks (discs of the radius rmask
    ! centred at the stars) cover a pixel

    use iso_fortran_env, only: int16

    real, dimension(:), intent(in) :: xstar, ystar, rmask
    integer, intent(in) :: ncol, nrow
    integer(int16), dimension(:,:), allocatable, intent(out) :: cover
    integer :: i, j, n

    allocate(cover(ncol,nrow))
    cover = 0
    do n = 1, size(xstar)
       do j = max(1,int(ystar(n) - rmask(n)) - 1), &
            min(nrow,int(ystar(n) + rmask(n)) + 1)
          do i = max(1,int(xstar(n) - rmask(n)) - 1), &
               min(ncol,int(xstar(n) + rmask(n)) + 1)
             if( (i - xstar(n))**2 + (j - ystar(n))**2 < rmask(n)**2 ) &
                  cover(i,j) = min(cover(i,j) + 1_int16,huge(cover))
          end do
       end do
    end do

  end subroutine starmask

  subroutine daomaskrad(d,xstar,ystar,lobad,hibad,gain,rn2,bias,rmax,rmask)

    ! The radii of the masks of stars: the light of a star is followed
    ! outwards by the mean of annuli 1 pixel wide while it is above the sky
    ! by 2 standard errors of the mean (by the scatter of the sky); the
    ! mask is one pixel larger. The sky and its scatter are of the whole
    ! frame. Saturated pixels count as bright. So the radius follows the
    ! actual profile and brightness (faint stars get small masks, bright,
    ! saturated stars and extended halos large ones), at most rmax.

    use skymean

    real, dimension(:,:), intent(in) :: d
    real, dimension(:), intent(in) :: xstar, ystar
    real, intent(in) :: lobad, hibad, rmax
    real(dbl), intent(in) :: gain, rn2, bias
    real, dimension(:), intent(out) :: rmask

    integer, parameter :: nsamp = 250000
    real(dbl), dimension(:), allocatable :: s
    real(dbl), dimension(:), allocatable :: asum
    integer, dimension(:), allocatable :: acnt
    real(dbl) :: sky, skyerr, skysig, rr
    integer :: ncol, nrow, nstep, n, i, j, k, nused, kmax

    ncol = size(d,1)
    nrow = size(d,2)

    ! the sky of the frame by a sample of pixels
    nstep = max(int(sqrt(real(ncol)*real(nrow)/nsamp)),1)
    s = pack(real(d(1:ncol:nstep,1:nrow:nstep),dbl), &
         lobad <= d(1:ncol:nstep,1:nrow:nstep) .and. &
         d(1:ncol:nstep,1:nrow:nstep) < hibad)
    call sky_mean(s,gain,rn2,bias,sky,skyerr,skysig,nused)

    kmax = max(int(rmax),1)
    allocate(asum(0:kmax+1),acnt(0:kmax+1))
    do n = 1, size(xstar)
       asum = 0
       acnt = 0
       do j = max(1,int(ystar(n)) - kmax - 1), min(nrow,int(ystar(n)) + kmax + 1)
          do i = max(1,int(xstar(n)) - kmax - 1), min(ncol,int(xstar(n)) + kmax + 1)
             rr = hypot(real(i - xstar(n),dbl),real(j - ystar(n),dbl))
             k = nint(rr)
             if( k > kmax ) cycle
             if( d(i,j) >= hibad ) then
                asum(k) = asum(k) + huge(1.0)
                acnt(k) = acnt(k) + 1
             else if( d(i,j) >= lobad ) then
                asum(k) = asum(k) + (d(i,j) - sky)
                acnt(k) = acnt(k) + 1
             end if
          end do
       end do
       rmask(n) = rmax
       do k = 1, kmax
          if( acnt(k) > 0 ) then
             if( asum(k) / acnt(k) < 2*skysig/sqrt(real(acnt(k),dbl)) ) then
                rmask(n) = min(k + 1.0,rmax)
                exit
             end if
          end if
       end do
    end do

  end subroutine daomaskrad

  subroutine daoapee(n,psf,i0,j0,xc,yc,raper,ecc,incl,ee)

    ! The enclosed energy: the fraction of the flux of a model star
    ! that the apertures of daophotsb contain. The pixelised model
    ! psf(-n:n,-n:n) (the whole flux is 1) has its element (0,0) at
    ! the pixel (i0,j0) of the frame, the star centre is xc,yc.
    ! The pixels are weighted exactly as in daophotsb (the same partial
    ! pixels, the same limits), so an aperture sum divided by ee is
    ! the total flux of the star.

    integer, intent(in) :: n, i0, j0
    real(dbl), dimension(-n:n,-n:n), intent(in) :: psf
    real, intent(in) :: xc, yc, ecc, incl
    real, dimension(:), intent(in) :: raper
    real(dbl), dimension(:), intent(out) :: ee

    real, parameter :: rad = 57.295779513082322865
    real :: apmxsq, rsq, dysq, x, y, r, phi, q, c, s, sqrte, fractn
    integer :: i, j, k

    apmxsq = -1.0
    do k = 1, size(raper)
       apmxsq = max(apmxsq, (raper(k)+0.5)**2)
    end do

    c = cos(incl/rad)
    s = sin(incl/rad)
    sqrte = sqrt(1 - ecc**2)

    ee = 0
    do j = j0 - n, j0 + n
       dysq = (j - yc)**2
       do i = i0 - n, i0 + n
          rsq = dysq + (i - xc)**2
          if( rsq <= apmxsq ) then
             call pixgeom(i,j,xc,yc,ecc,c,s,sqrte,x,y,r,phi,q)
             do k = 1, size(raper)
                if( r <= raper(k)+1 ) then
                   fractn = pixfrac(raper(k),ecc,sqrte,x,y,r,phi,q)
                   ee(k) = ee(k) + fractn*psf(i-i0,j-j0)
                end if
             end do
          end if
       end do
    end do

  end subroutine daoapee

  pure subroutine pixgeom(i,j,xc,yc,e,c,s,sqrte,x,y,r,phi,q)

    ! the pixel (i,j) in the frame of the (elliptic) aperture
    ! centred at xc,yc: the rotated offsets, the radius, the angle
    ! and the radius factor of the ellipse

    integer, intent(in) :: i,j
    real, intent(in) :: xc,yc,e,c,s,sqrte
    real, intent(out) :: x,y,r,phi,q
    real :: u,v

    u = i - xc
    v = j - yc
    x = c*u + s*v
    y =-s*u + c*v
    if( abs(x) > 0 .or. abs(y) > 0 ) then
       r = sqrt(x**2 + y**2)
       phi = atan2(y,x*sqrte)
    else
       r = 0
       phi = 0
    end if
    q = sqrt(1 - e**2*sin(phi)**2)

  end subroutine pixgeom

  pure real function pixfrac(a,e,sqrte,x,y,r,phi,q)

    ! the fraction of the pixel (by pixgeom) inside the aperture
    ! of the radius (major semiaxis) a

    real, intent(in) :: a,e,sqrte,x,y,r,phi,q
    real :: b,xt,yt

    if( abs(e) < epsilon(e) ) then
       pixfrac = max(0.0, min(1.0,a - r + 0.5))
    else
       b = a*sqrte
       xt = a*cos(phi)
       yt = b*sin(phi)
       pixfrac = circfrac(a,b,phi,xt-x,yt-y,a*q>r)
       ! pixfrac = sqcover(x,y,a)
    end if

  end function pixfrac

  pure real function circfrac(a,b,phi,dx,dy,outer)

    ! Estimates relative fraction of pixel intesrected with aperture circle edge
    ! The pixel is approximated by circle, the edge by a line.

    real, parameter :: pi = 3.14159265358979323844
    real, parameter :: rc = 0.564 ! = 1/sqrt(pi) radius of circumcircle
    real, parameter :: Area = pi*rc**2      ! circle area, should be == 1

    real, intent(in) :: a,b,phi,dx,dy
    logical, intent(in) :: outer
    real :: k,q,D,w,s,x1,x2,y1,y2,t,alpha,dA

    ! singular points
    if( abs(a) < epsilon(a) .or. abs(phi) < epsilon(phi) ) then
       if( outer ) then
          circfrac = 1
       else
          circfrac = 0
       end if
       return
    end if

    ! the intersection
    ! https://en.wikipedia.org/wiki/Intersection_(Euclidean_geometry)

    ! tangent: y = kx + q
    k = - (b/a) / tan(phi)
    q = dy - k*dx

    ! discriminant
    w = 1 + k**2
    D = rc**2*w - q**2

    if( D > 0 ) then
       s = sqrt(D)
       x1 = (-q*k + s) / w
       x2 = (-q*k - s) / w
       y1 = (q + k*s) / w
       y2 = (q - k*s) / w
       t = sqrt((y2 - y1)**2 + (x2 - x1)**2)
       alpha = 2*asin(min(t/rc/2,1.0))
       ! https://en.wikipedia.org/wiki/Circular_segment
       dA = rc**2*(alpha - sin(alpha)) / 2
       if( outer ) then
          circfrac = (Area - dA) / Area
       else
          circfrac = dA / Area
       end if
    else
       if( outer ) then
          circfrac = 1
       else
          circfrac = 0
       end if
    end if

  end function circfrac

  real function sqcover(xpix,ypix,raper)

    ! This routine computes coveradge of pixel on the edge by numerical way.
    ! It's dumy, slow, designed for testing purposes.

    real, intent(in) :: xpix,ypix,raper

    real :: x,y
    integer :: i,j,n,k

    n = 0
    k = 0
    do i = -50, 50
       x = xpix + i / 100.0
       do j = -50, 50
          y = ypix + j / 100.0
          if( sqrt(x**2 + y**2) < raper ) k = k + 1
          n = n + 1
       end do
    end do

    sqcover = real(k) / real(n)

  end function sqcover

end module daofotometr
