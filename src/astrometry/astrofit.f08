!
! astrofit - the astrometry solution
!
! Copyright © 2011-2026 F.Hroch (hroch@physics.muni.cz)
!
! This file is part of Munipack.
!
! Munipack is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Munipack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
!

module astrofitting

  use iso_fortran_env
  implicit none
  private

  public :: astrofit

contains

  subroutine astrofit(fit,ptype,pcode,a,d,x,y,sig,xsig,xcen,ycen, &
       acen,dacen,dcen,ddcen,sc,dsc,refl,rms,rot,drot,pv,ares,dres, &
       id1_init,id2_init,id1_final,id2_final,rmscheck,verbose,plog,reliable, &
       scale0,scaletol,frame,ntrial)

    ! scale0, scaletol (optional): the expected scale [pix/deg] and its
    ! relative tolerance; a solution with a scale out of the tolerance is rejected;
    ! frame (optional): the size of the frame [pix] (NAXIS1, NAXIS2),
    ! ntrial (optional, 1 by default): the number of the transformations
    ! tried by the match (the starting pairs of the sequences), both for
    ! the test of the significance of the solution

    use robfit
    use deform
    use astrotrafo
    use estimator
    use nearmatch

    character(len=*), intent(in) :: fit,ptype,pcode
    logical, intent(in) :: verbose, plog, rmscheck
    logical, intent(out) :: reliable
    real(REAL64), dimension(:), intent(in) :: a, d, x, y
    real(REAL64), intent(in) :: xcen,ycen,sig,xsig
    real(REAL64), intent(out) :: acen,dcen,dacen,ddcen,sc,dsc,rms,rot,drot,refl
    real(REAL64), dimension(:,0:), intent(out) :: pv
    real(REAL64), dimension(:), allocatable, intent(out) :: ares, dres
    integer, dimension(:), intent(in) :: id1_init,id2_init
    integer, dimension(:), allocatable, intent(out) :: id1_final,id2_final
    real(REAL64), intent(in), optional :: scale0, scaletol
    integer, dimension(:), intent(in), optional :: frame
    real(REAL64), intent(in), optional :: ntrial

    real(REAL64), dimension(:), allocatable :: alpha,delta,xstar,ystar
    integer, dimension(:), allocatable :: id1, id2
    type(AstroWCS) :: wcs
    character(len=80) :: fitmethod
    real(REAL64) :: dlim, s
    logical :: status
    integer :: n,nra,ndec,ndat

    reliable = .false.
    fitmethod = fit

    ! check dimensions
    if( size(a) /= size(d) ) &
         stop 'Astrofit: alpha and delta sizes does not corresponds.'
    if( size(x) /= size(y) ) &
         stop 'Astrofit: xstar and ystar sizes does not corresponds.'

    ! check stars > 1
    if( .not. (size(id1_init) > 1) ) then
       rms = huge(rms)
       return
    end if

    id1 = id1_init
    id2 = id2_init

    ! select matched stars and compute centre base image coordinates
    n = size(id1)
    allocate(alpha(n),delta(n),xstar(n),ystar(n))
    alpha = a(id1)
    delta = d(id1)
    xstar = x(id2) - xcen
    ystar = y(id2) - ycen

    ! intial estimates of all scale, rotation, reflection, and centre of projection
    call inestim(pcode,alpha,delta,xstar,ystar,acen,dcen,sc,dsc, &
         rot,drot,refl,rms,verbose)

    if( .not. (rms > 0) ) then
       ! rms is the statistical scale for robust method, it should be > 0
       fitmethod = 'SQUARES'
       if( verbose .and. fitmethod /= fit ) write(error_unit,*) &
            'Warning in astrofit(): the fitting method switched to "squares"', &
            ' as the statistical scale is zero.'
    end if

    ! the initial robust estimation
    if( plog ) &
         write(*,'(a,en12.3,en13.3,f8.1,2f14.5)') "=AFIT> ",rms,sc,rot,acen,dcen

    ! estimate preliminary solution based on a small set of matched stars
    s = rms
    call robmin(pcode,'ROBUST',alpha,delta,refl,xstar,ystar,s,&
         acen,dacen,dcen,ddcen,sc,dsc,rot,drot,rms,verbose,reliable)
    if( .not. reliable ) return
    if( .not. scale_prior() ) return

    ! select such much stars as possible for final fit
    deallocate(alpha,delta,xstar,ystar)
    deallocate(id1,id2)

    ! In fact, the limit xsig is limit of clipping of data, sig [pix], s[deg]
    dlim = xsig*max(sig / sc,rms) ! [deg]
!    write(error_unit,*) 'll',dlim,sc,rms*sc,sig,sig/sc
    if( verbose ) &
         write(error_unit,'(a,f0.2)') '# Cross match limit for stars [pix]: ',dlim*sc
    call listmatch(pcode,a,d,x,y,acen,dcen,xcen,ycen,sc,rot,refl,dlim,id1,id2,status)
    n = size(id1)
    if( .not. n > 1 ) then
       reliable = .false.
       return
    end if

    if( size(id1) /= size(id2) ) stop 'Error in astrofit(): unconformed sequences.'
    allocate(alpha(n),delta(n),xstar(n),ystar(n))
    alpha = a(id1)
    delta = d(id1)
    xstar = x(id2) - xcen
    ystar = y(id2) - ycen

    if( verbose ) &
         write(error_unit,'(a,i0)') '# cross-matched stars to fit: ',n

    fitmethod = fit
    if( .not. (s > 0) ) then
       fitmethod = 'SQUARES'
       if( verbose .and. fitmethod /= fit ) write(error_unit,*) &
            'Warning in astrofit(): the fitting method switched to "squares"', &
            ' as the statistical scale is zero.'
    end if

    if( n > 5 ) then

       call robmin(pcode,fitmethod,alpha,delta,refl,xstar,ystar,s,&
            acen,dacen,dcen,ddcen,sc,dsc,rot,drot,rms,verbose,reliable)
       if( .not. reliable ) return

       if( ptype == 'TPV' ) then

          if( fitmethod == 'ROBUST' .and. n > size(pv) + 1 ) then

             call tpvstaged
             if( .not. reliable ) return
          else if( fitmethod == 'ROBUST' ) then
             if( verbose ) write(error_unit,*) &
                  'Warning in astrofit(): too few stars for deformations.'
             pv = 0
             pv(1,1) = 1
             pv(2,1) = 1
          else
             if( verbose ) write(error_unit,*) &
                  'Warning in astrofit(): deformations needs robust methods.'
             pv = 0
             pv(1,1) = 1
             pv(2,1) = 1
          end if
       end if

       if(plog) write(*,'(a,en12.3,en13.3,f8.1,2f14.5)') &
            "=RFIT> ",rms,sc,rot,acen,dcen
       if( .not. scale_prior() ) return
       if( rmscheck ) then
          if( .not. significant() ) return
       end if

    else

       write(error_unit,*) &
            "Warning in astrofit: low-precision solution having a few points only."
       reliable = n > 2
       return

    end if

    if( (rmscheck .and. rms*sc < xsig*sig ) .or. .not. rmscheck ) then

       block
         real(REAL64), parameter :: rad = 45 / atan(1.0_REAL64)
         real(REAL64), dimension(n) :: a,d

         if( ptype == 'TPV' ) then
            call wcs % init(ptype,[acen,dcen],[xcen,ycen],sc,rot,refl,pv)
         else
            call wcs % init(ptype,[acen,dcen],[xcen,ycen],sc,rot,refl)
         end if

         ! calculate residuals
         call wcs % pix2cel(x(id2),y(id2),a,d)
         if( ptype /= '' ) then
            ! RA residuals as arcs on the sky (scaled by cos(delta)),
            ! the difference is reduced into [-180,180) to handle 0/360 wrap
            ares = (modulo(alpha - a + 180, 360.0_REAL64) - 180) * cos(delta / rad) ! [deg]
         else
            ! relative astrometry: rectangular coordinates [pix]
            ares = alpha - a
         end if
         dres = delta - d  ! [deg]

       end block

       if( verbose ) then

          ! the sign test, var = n*p*(1-p), where p = 0.5
          ndat = n
          nra = count(ares > 0)
          ndec = count(dres > 0)

          write(error_unit,'(a,i7,"/",f0.1,2(i5,"+-",f0.1))') &
               '# sign test (total/expected, RA+, Dec+): ',&
               ndat,ndat/2.0,nra,sqrt(nra*0.25),ndec,sqrt(ndec*0.25)

          call wcs % print_header
       end if

       id1_final = id1
       id2_final = id2
    else
       ! the solution is rejected, no residuals are available
       reliable = .false.
       allocate(id1_final(0),id2_final(0))
       if( verbose ) then
          write(error_unit,*) 'rms*scale < xsig*sig: ',rms*sc,xsig*sig
          write(error_unit,*)

          write(error_unit,*) 'Most common occurrences:'
          write(error_unit,*) &
               '  * crowded field mesh: increase --minmatch 7 (or more)'
          write(error_unit,*) &
               '  * wide-field deformation, use TPV as the projection'
       end if
    end if

    deallocate(alpha,delta,xstar,ystar)
    deallocate(id1,id2)

  contains

    logical function significant()

      ! The number of the cross-matched stars must be significantly larger
      ! than the number of chance coincidences of the frame stars with
      ! the catalogue (a false match gives a few stars only, comparable
      ! with the chance). The test counts the stars whose residuals of the
      ! solution are within rtest = 3 times their scatter, at least --sig,
      ! at most the limit of the cross-match. Two stars fix the solution
      ! (they match by construction), the others are tested against the
      ! Poisson distribution of chance coincidences by the density of the
      ! catalogue stars in the overlap of the catalogue and the frame;
      ! the probability is multiplied by the number of transformations
      ! tried by the match (trials, the Bonferroni correction).

      real(REAL64), parameter :: pmax = 1e-3
      integer, parameter :: nfix = 2
      type(AstroWCS) :: twcs
      real(REAL64), dimension(n) :: xc, yc
      real(REAL64) :: lambda, pval, rtest, trials
      integer :: nt

      if( ptype == 'TPV' ) then
         call twcs % init(ptype,[acen,dcen],[xcen,ycen],sc,rot,refl,pv)
      else
         call twcs % init(ptype,[acen,dcen],[xcen,ycen],sc,rot,refl)
      end if
      call twcs % cel2pix(alpha,delta,xc,yc)
      rtest = min(dlim*sc,max(3*rms*sc,sig))
      nt = count(hypot(xc - (xstar + xcen),yc - (ystar + ycen)) < rtest)

      trials = 1
      if( present(ntrial) ) trials = max(ntrial,1.0_REAL64)
      call chance(pcode,a,d,x,y,acen,dcen,xcen,ycen,sc,rot,refl,rtest, &
           nt,nfix,trials,lambda,pval,frame)
      significant = pval <= pmax
      if( verbose ) write(error_unit,'(a,i0,a,f0.2,a,i0,a,es9.2e3,a,es9.2e3,a,es9.2e3,a,l1)') &
           '# significance: ',nt,' stars within ',rtest,' pix (',nfix, &
           ' fix the solution), expected by chance: ',lambda,', trials: ', &
           trials,', probability by chance: ',pval,', significant: ',significant
      if( .not. significant ) then
         if( verbose ) write(error_unit,'(a)') &
              '# the match is not significant against chance coincidences, rejected.'
         reliable = .false.
         allocate(id1_final(0),id2_final(0))
      end if

    end function significant

    logical function scale_prior()

      ! the scale of the solution against the scale prior (if any),
      ! a solution out of the tolerance is rejected

      scale_prior = .true.
      if( .not. (present(scale0) .and. present(scaletol)) ) return
      if( .not. (scale0 > 0 .and. scaletol > 0) ) return
      scale_prior = abs(sc / scale0 - 1) <= scaletol
      if( .not. scale_prior ) then
         if( verbose ) write(error_unit,'(a,f0.1,a,f0.1,a,f0.1,a)') &
              '# the scale of the solution ',sc,' [pix/deg] is out of the expected ', &
              scale0,' +- ',100*scaletol,'%, rejected.'
         reliable = .false.
         if( .not. allocated(id1_final) ) allocate(id1_final(0),id2_final(0))
      end if

    end function scale_prior

    subroutine tpvstaged

      ! The staged fit of the TPV deformation. The stars cross-matched by
      ! the TAN solution miss the edges of the field where the deformation
      ! moves them by more than the match limit, so a polynomial fitted on
      ! them is extrapolated there. The deformation of the order 2 is
      ! fitted first, all stars are cross-matched again by that solution
      ! (within the limit by the scatter of the fit, xsig*max(sig,rms)),
      ! and refitted; when the matched set does not change, the order is
      ! raised to 3 and the iterations continue until the set is stable.

      integer, parameter :: maxstage = 12
      real(REAL64), dimension(size(pv,1),0:size(pv,2)-1) :: pvi
      integer, dimension(:), allocatable :: j1, j2
      real(REAL64) :: tlim
      integer :: order, stage, nset
      logical :: stable

      pvi = 0
      order = 2
      stable = .false.
      do stage = 1, maxstage
         call defmin(pcode,alpha,delta,refl,xstar,ystar,s,&
              acen,dcen,sc,rot,rms,pv,verbose,reliable,order,pvi)
         if( .not. reliable ) return

         ! the cross-match by the solution
         tlim = xsig*max(sig / sc,rms)
         call defmatch(pcode,a,d,x,y,acen,dcen,xcen,ycen,sc,rot,refl,pvi, &
              tlim,j1,j2)
         nset = size(j1)
         stable = nset == n
         if( stable ) stable = all(j1 == id1) .and. all(j2 == id2)
         if( verbose ) write(error_unit,'(a,i0,a,i0,a,i0,a,f0.2,a,i0,a,l1)') &
              '# TPV stage ',stage,', order ',order,': ',n,' stars fitted, rms ', &
              rms*sc,' pix, cross-matched by the solution: ',nset,', stable: ',stable
         if( stable .and. order == 3 ) exit

         if( .not. stable ) then
            if( nset <= size(pv) + 1 ) then
               ! too few stars by the solution: the last fit is kept
               if( order == 3 ) exit
            else
               n = nset
               call move_alloc(j1,id1)
               call move_alloc(j2,id2)
               deallocate(alpha,delta,xstar,ystar)
               allocate(alpha(n),delta(n),xstar(n),ystar(n))
               alpha = a(id1)
               delta = d(id1)
               xstar = x(id2) - xcen
               ystar = y(id2) - ycen
            end if
         end if
         if( stable .or. stage >= maxstage/2 ) order = 3
      end do

      ! the final solution of the order 3 on the last set
      if( .not. (stable .and. order == 3) ) then
         call defmin(pcode,alpha,delta,refl,xstar,ystar,s,&
              acen,dcen,sc,rot,rms,pv,verbose,reliable,3,pvi)
      end if

    end subroutine tpvstaged

  end subroutine astrofit

  subroutine chance(type,a,d,x,y,acen,dcen,xcen,ycen,sc,rot,refl,r, &
       n,nfix,ntrial,lambda,pval,frame)

    ! The expected number lambda of chance coincidences of frame stars
    ! (x,y) with catalogue stars (a,d) within the radius r [pix],
    ! and the probability pval to get n or more such coincidences
    ! (Poisson), nfix of them fixing the transformation (they are not
    ! tested), multiplied by the number of trials (the transformations
    ! tried). The stars are supposed to be distributed uniformly over
    ! the overlap of the catalogue (its bounding box) and the frame
    ! (frame = NAXIS1, NAXIS2; the bounding box of the frame stars
    ! when it is not given).

    use astrotrafo

    character(len=*), intent(in) :: type
    real(REAL64), dimension(:), intent(in) :: a,d,x,y
    real(REAL64), intent(in) :: acen,dcen,xcen,ycen,sc,rot,refl,r,ntrial
    integer, intent(in) :: n, nfix
    real(REAL64), intent(out) :: lambda, pval
    integer, dimension(:), intent(in), optional :: frame

    real(REAL64), parameter :: pi = 3.14159265358979323846_REAL64
    real(REAL64), dimension(size(a)) :: u,v,xc,yc
    real(REAL64) :: xmin,xmax,ymin,ymax,area,p,t
    type(AstrometryProjection) :: tp
    type(AstrometryGeometry) :: tg
    integer :: ncat,nfr,k,m

    call tp % init(type,acen,dcen)
    call tp % proj(a,d,u,v)
    call tg % init(sc,rot,refl)
    call tg % uv2pix(u,v,xc,yc)
    xc = xc + xcen
    yc = yc + ycen

    ! the overlap of the catalogue and the frame
    if( present(frame) ) then
       xmin = 0.5_REAL64
       xmax = frame(1) + 0.5_REAL64
       ymin = 0.5_REAL64
       ymax = frame(2) + 0.5_REAL64
    else
       xmin = minval(x)
       xmax = maxval(x)
       ymin = minval(y)
       ymax = maxval(y)
    end if
    xmin = max(xmin,minval(xc))
    xmax = min(xmax,maxval(xc))
    ymin = max(ymin,minval(yc))
    ymax = min(ymax,maxval(yc))
    area = (xmax - xmin)*(ymax - ymin)
    lambda = 0
    pval = 1
    if( .not. (area > 0) ) return

    ncat = count(xmin <= xc .and. xc <= xmax .and. ymin <= yc .and. yc <= ymax)
    nfr = count(xmin <= x .and. x <= xmax .and. ymin <= y .and. y <= ymax)

    ! probability of a chance coincidence for a frame star
    p = min(max(ncat - nfix,0)*pi*r**2 / area, 1.0_REAL64)
    lambda = max(nfr - nfix,0)*p
    m = n - nfix

    ! the upper tail of Poisson distribution, P(N >= m)
    if( m <= 0 .or. m <= lambda ) then
       pval = 1
    else if( lambda > 0 ) then
       pval = 0
       do k = m, m + 1000
          t = exp(k*log(lambda) - lambda - log_gamma(k + 1.0_REAL64))
          pval = pval + t
          if( t < epsilon(t)*pval ) exit
       end do
    else
       pval = 0
    end if
    pval = min(ntrial*pval,1.0_REAL64)

  end subroutine chance

end module astrofitting
