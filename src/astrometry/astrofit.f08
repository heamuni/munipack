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
       id1_init,id2_init,id1_final,id2_final,rmscheck,verbose,plog,reliable)

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
    if( .not. n > 1 ) return
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

          if( fitmethod == 'ROBUST' ) then

             call defmin(pcode,alpha,delta,refl,xstar,ystar,s,&
                  acen,dcen,sc,rot,rms,pv,verbose,reliable)
             if( .not. reliable ) return
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
         ares = (alpha - a) !* cos(delta / rad) ! [deg]
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

  end subroutine astrofit

end module astrofitting
