!
!  FITS I/O for photometry calibration
!
!  Copyright © 2014-2026 F.Hroch (hroch@physics.muni.cz)
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

module sfits

  use titsio
  use iso_fortran_env

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: rad = 57.295779513082322865_dbl

  type qFITS
     character(len=FLEN_FILENAME) :: filename
     real, dimension(:,:), allocatable :: frame
     integer :: naxis, nrows
     integer, dimension(2) :: naxes
     real(dbl), dimension(:),allocatable :: ra,dec,cts,dcts,sky,dsky,&
          ph,dph,phsky,dphsky,mag,dmag
     real(dbl) :: exptime,coo_err,aper,area,scale,rsky,rskydev,skymag,&
          radf90, airmass, jd, hjd, longitude, latitude
     real(dbl), dimension(2) :: annuls
     character(len=FLEN_VALUE) :: filter, phtype
     ! the enclosed energy (APEE) of the aperture applied, its source
     logical :: apee = .false.
     character(len=FLEN_VALUE) :: eesource = ''
     logical :: geodefined = .false.
     logical :: status
  end type qFITS

  private :: readfits, qfits_init, qfits_deallocate, qfits_setph, colkeys, &
       colcomments

contains



  subroutine writecal(filenames,outputs,keys,advanced,phsystable, &
       filters,catid,photsys,area,init_area,saper,tra,ctph,dctph,quantities,phres,&
       verbose,syserr,apcorr,ctm,sysint,nsys)

    ! syserr (optional) are the systematic errors [mag] of the calibration
    ! of the frames (by ctph), added in quadrature to the errors of all
    ! objects; apcorr (optional) is the constant aperture correction
    ! (for the header only, ctph includes it); ctm (optional) the colour
    ! term (colterm): COLT_CAT applies the colour term by the catalogue
    ! colour of every star (the stars without it by the mean correction
    ! of the calibration stars), COLT_TRA the magnitude transformation by
    ! the colours of the calibrated magnitudes (tra is the identity then);
    ! sysint (optional) the intervals of the estimated systematic errors
    ! [mag] and nsys the numbers of their stars (header keys PHSYSERE,
    ! PHSYSNST; nsys < 0 for a fixed systematic error, no keys)

    use phsysfits
    use oakleaf
    use jamming
    use calibre
    use fotran
    use mfits
    use colterm

    real, parameter :: pi = 3.14159

    character(len=*), dimension(:), intent(in) :: filenames,outputs,keys, &
         quantities, filters
    character(len=*), intent(in) :: photsys, catid, phsystable
    real(dbl), dimension(:,:), intent(in) :: tra
    real(dbl), dimension(:), intent(in) :: ctph,dctph
    type(photores), dimension(:), intent(in), optional :: phres
    real(dbl), intent(in) :: area
    integer, intent(in) :: saper
    logical, intent(in) :: init_area, advanced, verbose
    real(dbl), dimension(:), intent(in), optional :: syserr
    real(dbl), intent(in), optional :: apcorr
    type(colour_term), intent(in), optional :: ctm
    real(dbl), dimension(:,:), intent(in), optional :: sysint
    integer, dimension(:), intent(in), optional :: nsys

    real(dbl), parameter :: nullval = 0
    integer, parameter :: dim2 = 2
    integer, dimension(:), allocatable :: ndat
    integer, dimension(:,:), allocatable :: idx,pairs
    type(qFITS), dimension(:), allocatable :: fitses
    type(type_phsys) :: phsyscal
    real(dbl), dimension(:), allocatable :: atol
    real(dbl), dimension(:,:), allocatable :: ras,decs,ct,dct,sky,dsky,ph,dph, &
         phsky,dphsky,mag,dmag,phr,dphr
    real, dimension(:,:,:), allocatable :: cube,tcube
    real(dbl), dimension(size(ctph)) :: xctph,sqarcsec,airmass,aring,nsky, &
         leff,lfwhm,rskies,rskydevs,srel
    real(dbl) :: apc
    real(dbl), dimension(1,size(ctph)) :: skymag,dskymag,skyrate,dskyrate
    integer, allocatable, dimension(:) :: order
    integer :: i,j,k,n,status,nrows_max, width, height, nq, ncat, nfiles, nrows, npht, nt
    real(dbl) :: rsky,rskydev,k0,r,lref
    ! the colour term: the colours of the stars (i,n,j), their errors (< 0
    ! undefined), the source of the colour (COLT_CAT: the bit j-1 set for
    ! the colour j of the catalogue, 0 the mean colours of the calibration
    ! stars; COLT_TRA: 2 the calibrated magnitudes), the covariances of
    ! the colours of a star (i,j,k)
    real(dbl), dimension(:,:,:), allocatable :: scol,sdcol,scv
    integer, dimension(:,:), allocatable :: csrc
    integer, dimension(:), allocatable :: cidx
    real(dbl), dimension(size(ctph)) :: m0,dm0,mm,dmm,cc,dcc
    real(dbl), dimension(3,3) :: cvp, ap
    real(dbl), dimension(3) :: cpp
    real(dbl), dimension(1) :: x1,dx1,y1,dy1
    real(dbl) :: g,dg,f
    real(dbl), dimension(2,size(ctph)) :: sint
    integer, dimension(size(ctph)) :: nst
    logical :: lcol

    nfiles = size(filenames)
    nq = size(quantities)
    npht = size(tra,1)

    ! the relative systematic errors (2.5/ln 10 = 1.0857 mag)
    srel = 0
    if( present(syserr) ) srel = max(syserr,0.0_dbl) / 1.0857362047581294_dbl
    apc = 1
    if( present(apcorr) ) apc = apcorr
    sint = -1
    nst = -1
    if( present(sysint) ) sint = sysint
    if( present(nsys) ) nst = nsys

    status = 0
    allocate(fitses(nfiles))

    ! init photometry parameters
    call phselect(phsystable,photsys,phsyscal)
    call phsyspairs(phsyscal,filters,pairs)
    call file_order(filenames,keys,phsyscal%filter,order)

    nrows_max = 0
    ! read frames
    do n = 1, nfiles
       call qfits_init(fitses(n))
       call readfits(filenames(order(n)),keys,saper,fitses(n),verbose,status)
       if( status /= 0 ) stop 'FITS read failed.'
       if( fitses(n)%nrows > nrows_max ) nrows_max = fitses(n)%nrows
       if( init_area ) fitses(n)%area = area
       airmass(n) = fitses(n)%airmass
       ! passband parameters of the frame
       leff(n) = -1
       do i = 1, phsyscal%nbands
          if( phsyscal%filter(i) == fitses(n)%filter ) then
             leff(n) = phsyscal%lam_eff(i)
             lfwhm(n) = phsyscal%lam_fwhm(i)
          end if
       end do
    end do

    ! convert frames from counts to photons
    if( nfiles == 1 ) then
       fitses(1)%frame = real(ctph(1)*fitses(1)%frame)
    else
       ! uncorrect handling for multiple frames with small mutual offsets !
       ! this code supposes all the fames with no mutual offset and the same dimensions
       width = fitses(1)%naxes(1)
       height = fitses(1)%naxes(2)
       allocate(cube(width,height,nfiles),tcube(width,height,nfiles))
       do n = 1,nfiles
          cube(:,:,n) = fitses(n)%frame
       end do

       do i = 1,width
          do j = 1,height
             !             call tra_frame(tra,pairs,ctph,cube(i,j,:),tcube(i,j,:))
             forall(n=1:npht)
                tcube(i,j,n) = real(ctph(n)*cube(i,j,n)*tra(n,n))
             end forall
          end do
       end do

       do n = 1,nfiles
          fitses(n)%frame = tcube(:,:,n)
       end do

       deallocate(cube,tcube)
    end if


    ! match both frames to be able to determine of colour indexes
    if( nfiles > 1 ) then

       allocate(ras(nrows_max,nfiles),decs(nrows_max,nfiles),atol(nfiles), &
            ndat(nfiles),idx(nrows_max,nfiles))
       idx = 0
       ras = 0
       decs = 0
       do n = 1,nfiles
          nrows = fitses(n)%nrows
          ndat(n) = nrows
          ras(1:nrows,n) = fitses(n)%ra
          decs(1:nrows,n) = fitses(n)%dec
       end do

       atol = 5.0*fitses%coo_err

       call jamframes(atol, ndat, ras, decs, idx)

       deallocate(ras,decs,atol,ndat)

    else ! nfiles == 1
       allocate(idx(nrows_max,nfiles))
       forall(i=1:size(idx,1))
          idx(i,1) = i
       end forall
    end if

    ncat = 0
    do i = 1,size(idx,1)
       if( all(idx(i,:) > 0) ) ncat = ncat + 1
    end do

    allocate(ct(ncat,nfiles),dct(ncat,nfiles),sky(ncat,nfiles), &
         dsky(ncat,nfiles),ras(ncat,nfiles),decs(ncat,nfiles))


    ncat = 0
    do i = 1,size(idx,1)
       if( all(idx(i,:) > 0) )then
          ncat = ncat  + 1
          do j = 1,nfiles
             k = idx(i,j)
             ras(ncat,j) = fitses(j)%ra(k)
             decs(ncat,j) = fitses(j)%dec(k)
             ct(ncat,j) = fitses(j)%cts(k)
             dct(ncat,j) = fitses(j)%dcts(k)
             sky(ncat,j) = fitses(j)%sky(k)
             dsky(ncat,j) = fitses(j)%dsky(k)
          end do
       end if
    end do
    deallocate(idx)

    do n = 1,nfiles
       deallocate(fitses(n)%ra,fitses(n)%dec)
       allocate(fitses(n)%ra(ncat),fitses(n)%dec(ncat))
       fitses(n)%ra = ras(:,n)
       fitses(n)%dec= decs(:,n)
    end do
    deallocate(ras,decs)

    ! the colours of the stars by the catalogue (the positions on
    ! the first frame)
    lcol = .false.
    if( present(ctm) ) lcol = ctm%model /= COLT_NONE
    if( lcol ) then
       nt = ctm%nterm
       allocate(scol(ncat,nfiles,nt),sdcol(ncat,nfiles,nt),csrc(ncat,nfiles), &
            scv(ncat,nt,nt))
       scol = 0
       sdcol = -1
       csrc = 0
       scv = 0
       if( ctm%model == COLT_CAT ) then
          allocate(cidx(ncat))
          call cooid(fitses(1)%ra,fitses(1)%dec,ctm%ra,ctm%dec,ctm%tol,cidx)
          do i = 1, ncat
             k = cidx(i)
             if( k > 0 ) then
                do j = 1, nt
                   if( ctm%dcol(k,j) >= 0 ) then
                      scol(i,:,j) = ctm%col(k,j)
                      sdcol(i,:,j) = ctm%dcol(k,j)
                      csrc(i,:) = ibset(csrc(i,:),j-1)
                   end if
                end do
                scv(i,:,:) = ctm%ccov(k,:,:)
             end if
          end do
          deallocate(cidx)
          if( verbose ) then
             if( nt == 1 ) then
                write(*,'(a,i0,a,i0,a)') ' Colour term: ',count(csrc(:,1) == 1), &
                     ' of ',ncat,' stars with a catalogue colour.'
             else
                write(*,'(a,i0,a,i0,a)') ' Colour terms: ', &
                     count(csrc(:,1) == 2**nt - 1),' of ',ncat, &
                     ' stars with all the catalogue colours,'
                do j = 1, nt
                   write(*,'(a,i0,a,i0,a)') '   colour ',j,': ', &
                        count(btest(csrc(:,1),j-1)),' stars'
                end do
             end if
          end if
       end if
    end if


    allocate(ph(ncat,nfiles),dph(ncat,nfiles),phsky(ncat,nfiles), &
         dphsky(ncat,nfiles))
    ph = -1
    dph = -1
    phsky = -1
    dphsky = -1

    do n = 1,size(fitses)
       xctph(n) = fitses(n)%area * fitses(n)%exptime
       !xctph(n) = 1
       ! scale is in [pix/deg], annuls and aper in [deg]
       ! square arcsec per pixel
       sqarcsec(n) = (3600 / fitses(n)%scale)**2
       ! unused, areas in pixels:
       aring(n) = pi*(fitses(n)%annuls(2)**2 - fitses(n)%annuls(1)**2) *  &
            fitses(n)%scale**2
       nsky(n) = pi*fitses(n)%aper**2*fitses(n)%scale**2
    end do

!    write(*,*) aring
!    write(*,*) nsky

    allocate(mag(ncat,nfiles),dmag(ncat,nfiles),phr(ncat,nfiles), &
         dphr(ncat,nfiles))

    ! No correction for the colour dependence of the extinction is applied
    ! (r < 0: no extinction model in the header). The correction used for
    ! more filters (colorex) estimated the extinction from CTPH, which is
    ! the efficiency of the whole instrument, and corrected the photons of
    ! the outputs only, not of the calibration: the magnitudes had a colour
    ! term against the calibration stars themselves (BUGS P43). The colour
    ! terms of an instrument, the second order extinction among them, are
    ! the transformation of phfotran (--tratab), fitted by the stars.
    if( nfiles > 1 .and. any(leff < 0) ) &
         stop 'Filter of a frame not found in photometry system.'
    k0 = 0
    r = -1
    lref = 0

    do i = 1,ncat
       if( nfiles == 1 )then
!          write(*,*) sqrt(ctph(1)),sqarcsec(1), dsky(i,1)
          where( sky(i,:) > 0 .and. dsky(i,:) > 0 )
             phsky(i,:) = ctph(1) * sky(i,:) !* sqarcsec(1)
             dphsky(i,:) = phsky(i,:) * (dsky(i,:)/sky(i,:))
!             dphsky(i,:) = sqrt(ctph(1)) * dsky(i,:) !* sqarcsec(1)
          elsewhere
             phsky(i,:) = -1
             dphsky(i,:) = -1
          end where

          where( ct(i,:) > 0 .and. dct(i,:) > 0 )
             ph(i,:) = ctph(1)*ct(i,:)
             dph(i,:) = sqrt(ctph(1)) * dct(i,:)  ! pure Poisson
             dph(i,:) = sqrt(ctph(1)**2*dct(i,:)**2 + ct(i,:)**2*dctph(1)**2 &
                  + (srel(1)*ph(i,:))**2)

             ! "gain" update
!             dph(i,:) = dph(i,:) / sqrt(ctph(1))

!             dph(i,:) = sqrt(dct(i,:)**2 + ct(i,:)**2*dctph(1)**2 + ctph(1)**2*dph(i,:)**2)


             ! Errors are estimated under Poisson statistics
             ! assumption, it can be underestimated for faint stars.

             ! fully optimal estimate: independent scaling Po and N parts
!             dph(i,:) = max(dct(i,:)**2 - ct(i,:),0.0)
!             dph(i,:) = sqrt(ph(i,:) + dph(i,:)*ctph(1)**2)


!             dph(i,:) = sqrt(ct(i,:)**2 * dctph(1)**2 + &
!                  ctph(1)**2 * dct(i,:)**2)
!             dph(i,:) = ph(i,:)*(dct(i,:)/ct(i,:))
             !dph(i,:) = sqrt(ph(i,:) + dphsky(i,:)**2*nsky(1)**2 + &
             !     dphsky(i,:)**2*pi*nsky(1)**2)
          elsewhere
             ph(i,:) = -1
             dph(i,:) = -1
          end where

          ! the colour term by the colour of the star: the photons
          ! ctph exp(g) ct, the error of the calibration at the colour
          ! (CTPH included) and of the colour
          if( lcol ) then
             if( ph(i,1) > 0 ) then
                call colour_corr(ctm,1,scol(i,1,:),sdcol(i,1,:),g,dg,scv(i,:,:))
                f = exp(g)
                ph(i,1) = ctph(1)*f*ct(i,1)
                dph(i,1) = sqrt((ctph(1)*f*dct(i,1))**2 + (ph(i,1)*dg)**2 &
                     + (srel(1)*ph(i,1))**2)
             end if
          end if

!          write(*,*) sqrt(abs(ph(i,:))),sqrt(dphsky(i,:)**2*nsky(1)**2), &
          !               sqrt(dphsky(i,:)**2*pi*nsky(1)**2)
       else
          ! the errors of the counts and of the calibration as for one frame
          if( lcol .and. ctm%model == COLT_CAT ) then
             ! the colour term of every filter by the catalogue colour
             do n = 1, nfiles
                call colour_corr(ctm,n,scol(i,n,:),sdcol(i,n,:),g,dg,scv(i,:,:))
                f = exp(g)
                if( ct(i,n) > 0 .and. dct(i,n) > 0 ) then
                   dct(i,n) = sqrt((ctph(n)*f*dct(i,n))**2 + &
                        (ctph(n)*f*ct(i,n)*dg)**2) / xctph(n)
                else
                   dct(i,n) = -1
                end if
                ct(i,n) = ctph(n) * f * ct(i,n) / xctph(n)
             end do
          else
             where( ct(i,:) > 0 .and. dct(i,:) > 0 )
                dct(i,:) = sqrt(ctph**2*dct(i,:)**2 + ct(i,:)**2*dctph**2) / xctph
             elsewhere
                dct(i,:) = -1
             end where
             ct(i,:) = ctph * ct(i,:) / xctph
          end if
          sky(i,:) = ctph * sky(i,:) / xctph
          dsky(i,:) = ctph * dsky(i,:) / xctph
          call fotra(tra,pairs,ct(i,:),dct(i,:),phr(i,:),dphr(i,:))
          call fotra(tra,pairs,sky(i,:),dsky(i,:),phsky(i,:),dphsky(i,:))
          where( phr(i,:) > 0 )
             dphr(i,:) = sqrt(dphr(i,:)**2 + (srel*phr(i,:))**2)
          end where
          if( lcol .and. ctm%model == COLT_TRA ) then
             ! the magnitude transformation: the magnitudes calibrated
             ! by CTPH, the colour terms by the colours of the standard
             ! magnitudes (solved), back to the photons
             if( all(phr(i,:) > 0 .and. dphr(i,:) > 0) ) then
                do n = 1, nfiles
                   x1 = phr(i,n)
                   dx1 = dphr(i,n)
                   call phsysphmag1(phsyscal,filters(n),x1,dx1,y1,dy1)
                   m0(n) = y1(1)
                   dm0(n) = magln * dphr(i,n) / phr(i,n)
                end do
                call colour_solve(ctm,m0,dm0,mm,dmm,cc,dcc,status)
                if( status == 0 ) then
                   do n = 1, nfiles
                      y1 = mm(n)
                      dy1 = 0
                      call phsysmagph1(phsyscal,filters(n),y1,dy1,x1,dx1)
                      phr(i,n) = x1(1)
                      dphr(i,n) = x1(1) * dmm(n) / magln
                   end do
                   scol(i,:,1) = cc
                   sdcol(i,:,1) = dcc
                   csrc(i,:) = 2
                else
                   phr(i,:) = -1
                end if
                status = 0
             else
                phr(i,:) = -1
             end if
          end if
          where( phr(i,:) > 0 )
             ph(i,:) = phr(i,:) * xctph
             dph(i,:) = dphr(i,:) * xctph !! correct !!!
          elsewhere
             ph(i,:) = - 1
             dph(i,:) = -1
          end where
          where( phsky(i,:) > 0 )
             phsky(i,:) = phsky(i,:) * xctph
             dphsky(i,:) = dphsky(i,:) * xctph !! correct needs!!!
          elsewhere
             phsky(i,:) = -1
             dphsky(i,:) = -1
          end where
       end if

    end do
!    write(*,*) ct(4,1),ph(4,1)


!    write(*,*) ctph,dctph
!    write(*,*) ct
!    write(*,*) dct
!    write(*,*) ph
!    write(*,*) dph


    do i = 1, nfiles
       ! a mean sky by the sky in rings
       call rmean(phsky(:,i),rsky,stdsig=rskydev)
       rskies(i) = rsky
       rskydevs(i) = rskydev
       skyrate(1,i) = (rsky / xctph(i)) / sqarcsec(i)
       dskyrate(1,i)= (rskydev / (xctph(i)) / sqarcsec(i))
!       write(*,*) rsky,skyrate(1,i),xctph(i),sqarcsec(i)
    end do

    if( nfiles == 1 ) then
       call phsysphmag1(phsyscal,filters(1),ph(:,1)/xctph(1), &
            dph(:,1)/xctph(1),mag(:,1),dmag(:,1))
!            dph(:,1)/sqrt(xctph(1)),mag(:,1),dmag(:,1))
       call phsysphmag1(phsyscal,filters(1),skyrate(1,:),dskyrate(1,:), &
            skymag(1,:),dskymag(1,:))
!       write(*,*) trim(filenames(1)),xctph(1),ph(4,1),mag(4,1)
    else
       call phsysphmag(phsyscal,filters,pairs,phr,dphr,mag,dmag)
       call phsysphmag(phsyscal,filters,pairs,skyrate,dskyrate,skymag,dskymag)
    end if

!    write(*,*) xctph
!    write(*,*) mag
!    write(*,*) dmag

    do n = 1,size(fitses)
       call qfits_setph(fitses(n),ph(:,n),dph(:,n),phsky(:,n),dphsky(:,n))

       fitses(n)%mag = mag(:,n)
       fitses(n)%dmag = dmag(:,n)
       ! meaning of rsky,rskydev has been changed during April 2026
       fitses(n)%rsky = rskies(n) ! skyrate(1,n)
       fitses(n)%rskydev = rskydevs(n) ! dskyrate(1,n)
       fitses(n)%skymag = skymag(1,n)

    end do

    ! the stars without a colour: the mean colour of the calibration stars
    ! and its scatter (as applied by colour_corr)
    if( lcol ) then
       do n = 1, nfiles
          if( ctm%model == COLT_CAT .and. nt == 1 ) then
             where( csrc(:,n) == 0 )
                scol(:,n,1) = ctm%cmean(1,n)
                sdcol(:,n,1) = ctm%cstd(1,n)
             end where
          else if( ctm%model == COLT_CAT ) then
             ! several colours: the missing ones predicted by the known
             ! ones (as applied by colour_corr), their conditional scatter
             do i = 1, ncat
                if( csrc(i,n) == 2**nt - 1 ) cycle
                call colour_predict(ctm,n,scol(i,n,:),sdcol(i,n,:),cpp(1:nt), &
                     cvp(1:nt,1:nt),ap(1:nt,1:nt))
                do j = 1, nt
                   if( .not. btest(csrc(i,n),j-1) ) then
                      scol(i,n,j) = cpp(j)
                      sdcol(i,n,j) = sqrt(max(cvp(j,j),0.0_dbl))
                   end if
                end do
             end do
          else
             where( csrc(:,n) == 0 )
                scol(:,n,1) = ctm%cmean(1,n)
                sdcol(:,n,1) = ctm%cstd(1,n)
             end where
          end if
       end do
    end if

    do n = 1,size(fitses)
       i = order(n)
       if( lcol ) then
          call savefits(filenames(i),outputs(i),fitses(n),advanced, &
               photsys,phsyscal,area,init_area,catid,quantities, &
               ctph(n),dctph(n),status,k0,r,lref,srel(n),apc,phres(n), &
               ctm,n,scol(:,n,:),sdcol(:,n,:),csrc(:,n),sysint=sint(:,n),nsys=nst(n))
       else if( present(phres) ) then
          call savefits(filenames(i),outputs(i),fitses(n),advanced, &
               photsys,phsyscal,area,init_area,catid,quantities, &
               ctph(n),dctph(n),status,k0,r,lref,srel(n),apc,phres(n), &
               sysint=sint(:,n),nsys=nst(n))
       else
          call savefits(filenames(i),outputs(i),fitses(n),advanced, &
               photsys,phsyscal,area,init_area,catid,quantities, &
               ctph(n),dctph(n),status,k0,r,lref,srel(n),apc, &
               sysint=sint(:,n),nsys=nst(n))
       end if
       call qfits_deallocate(fitses(n))
    end do
    if( lcol ) deallocate(scol,sdcol,csrc,scv)

    deallocate(ct,dct,sky,dsky,ph,dph,mag,dmag,phr,dphr, &
         phsky,dphsky,pairs,order)

  end subroutine writecal

  subroutine readfits(filename,keys,saper,fits,verbose,status)

    ! this is candidate for merge with the analogical subroutine in mfits.f95

    use mfits, only: phmethod, apee_use, apee_apply
    use astrotrafo
    use astrosphere
    use trajd

    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: keys
    integer, intent(in) :: saper
    type(qFITS), intent(in out) :: fits
    logical, intent(in) :: verbose
    integer, intent(in out) :: status

    real, parameter :: nullval = -huge(1.0)
    real(REAL64), parameter :: dullval = -huge(real(0.0,REAL64))
    integer, parameter :: group = 1, extver = 0
    integer :: frow
    real(dbl), dimension(:), allocatable :: xcen, ycen
    logical :: anyf
    character(len=FLEN_CARD) :: key, dateobs
    character(len=FLEN_CARD), dimension(0) :: exclist
    character(len=:), allocatable :: header
    integer :: w,h,i,l,ccol,dcol,xcol,ycol,scol,ecol,nrows,maper, naper, srows, nkeys
    real(dbl), dimension(2) :: crval, crder
    real(dbl) :: ls, jd0, scale,rot,refl
    real :: hwhm, aper
    logical :: have_grow
    type(fitsfiles) :: fitsfile
    type(AstroWCS) :: wcs

    fits%filename = filename

    call fits_open_image(fitsfile,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read an image in the file `',trim(filename),"'."
       return
    end if

    call fits_get_img_size(fitsfile,2,fits%naxes,status)
    if( status /= 0 ) goto 666

    ! WCS
    call fits_hdr2str(fitsfile,1,exclist,nkeys,header,status)
    if( status == 0 ) then
       call wcs % init(header,nkeys,status)
       deallocate(header)
    end if
    if( status /= 0 ) then
       write(error_unit,*) 'Error: Astrometry keywords not found in header.'
       goto 666
    end if
    call wcs % get_affines(scale,rot,refl)
    call wcs % get_crval(crval)
    call wcs % get_crder(crder)
    fits%scale = scale

    ! Citizen date
    call fits_get_dateobs(fitsfile,(/keys(3),keys(8)/),dateobs,status)
    if( status == 0 ) then
       fits%jd = fits_jd(dateobs,status)

       block
         real(REAL64), parameter :: one = real(1.0,REAL64)
         integer :: year, month, day, hour, minute
         real(dbl) :: second
         call fits_str2date(dateobs,year,month,day,hour,minute,second,status)

         ! Julian date of 1. january
         jd0 = datjd(real(year,REAL64),one,one)

       end block

    else
       write(error_unit,*) "Warning: Failed to get date or time."
       fits%jd = 0
       jd0 = 0
       status = 0
    end if

    call fits_read_key(fitsfile,keys(1),fits%exptime,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       write(error_unit,*) "Warning: Exposure time identified by FITS keyword `", &
            trim(keys(1)),"' not found (default to 1 sec)."
       fits%exptime = 1
       status = 0
    end if

    call fits_read_key(fitsfile,keys(2),fits%area,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
!       write(error_unit,*) "Warning: Area identified by FITS keyword `",trim(keys(2)), &
!            "' not found (default to 1 [m2])."
       fits%area = 1
       status = 0
    end if

!    call ftgkyd(25,keys(3),fits%airmass,buf,status)
!    if( status == KEYWORD_NOT_FOUND ) then
!       fits%airmass = 1
!       status = 0
!    end if

    call fits_read_key(fitsfile,keys(5),fits%filter,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       write(error_unit,*) "Warning: Filter identified by FITS keyword `",trim(keys(5)), &
            "' not found (default to empty `')."
       fits%filter = ''
       status = 0
    end if

    call fits_read_key(fitsfile,keys(6),fits%longitude,status)
    call fits_read_key(fitsfile,keys(7),fits%latitude,status)
    fits%geodefined = status == 0
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       if( verbose ) write(error_unit,*) &
            "Warning: Geographical longitude or latitude identified by FITS keywords `",&
            trim(keys(6)),"',`",trim(keys(7)),"' not found. Airmass is undefined."
       status = 0
    end if

    ! Julian date at midpoint of exposure
    fits%jd = fits%jd + fits%exptime / 2.0_dbl / 86400.0_dbl

    ! Heliocentric JD
    ls = longsun(fits%jd - jd0)
    fits%hjd = fits%jd + helcor(crval(1),crval(2),ls)

    ! airmass at the time
    if( fits%geodefined ) &
         fits%airmass = xairmass(fits%jd,fits%longitude,fits%latitude, &
         crval(1),crval(2))

    ! frame
    w = fits%naxes(1)
    h = fits%naxes(2)
    allocate(fits%frame(w,h))
    call fits_read_image(fitsfile,group,nullval,fits%frame,anyf,status)

    ! move to aperture photometry extension
    call fits_movnam_hdu(fitsfile,FITS_BINARY_TBL,APEREXTNAME,0,status)
    if( status == FITS_BAD_HDU_NUM .and. phmethod == 'PSF' ) then
       ! PSF photometry without aperture photometry: no annulus
       status = 0
       fits%annuls = 0
       fits%aper = -1
    else if( status == FITS_BAD_HDU_NUM ) then
       write(error_unit,*) &
            "Failed to find an aperture photometry extension in ",trim(filename)
       goto 666
    else

    ! sky annulus
    do i = 1, 2
       call fits_make_keyn(FITS_KEY_ANNULUS,i,key,status)
       call fits_read_key(fitsfile,key,fits%annuls(i),status)
    end do
    fits%annuls = fits%annuls / scale

    ! aperture
    call fits_read_key(fitsfile,FITS_KEY_HWHM,hwhm,status)
    if( saper == 0 ) then
       call fits_read_key(fitsfile,FITS_KEY_NAPER,naper,status)
       do i = 1, naper
          call fits_make_keyn(FITS_KEY_APER,i,key,status)
          call fits_read_key(fitsfile,key,aper,status)
          ! the first aperture larger than 2*sigma of Gaussian profile
          if( aper > 2*hwhm / sqrt(2*log(2.0)) ) then
             maper = i
             exit
          end if
       end do
    else
       maper = saper
    end if
    call fits_make_keyn(FITS_KEY_APER,maper,key,status)
    call fits_read_key(fitsfile,key,fits%aper,status)
    fits%aper = fits%aper / scale !* tastr%scale

    end if

    ! move to growth-curve photometry extension
    call fits_movnam_hdu(fitsfile,FITS_BINARY_TBL,GROWEXTNAME,0,status)
    have_grow = .not. (status == FITS_BAD_HDU_NUM)
    if( status == FITS_BAD_HDU_NUM ) then
       ! no problem, we're continuing with aperture extension
       status = 0
    end if

    if( phmethod == 'PSF' ) then
       call fits_movnam_hdu(fitsfile,FITS_BINARY_TBL,PSFEXTNAME,0,status)
       if( status == FITS_BAD_HDU_NUM ) then
          write(error_unit,*) &
               "Failed to find a PSF photometry extension in ",trim(filename)
          goto 666
       end if
       call fits_read_key(fitsfile,FITS_KEY_HWHM,hwhm,status)
       call fits_get_colnum(fitsfile,.true.,FITS_COL_PSFCOUNT,ccol,status)
       call fits_get_colnum(fitsfile,.true.,FITS_COL_PSFCOUNTERR,dcol,status)
       fits%radf90 = -1
       fits%aper = -1
       fits%phtype = 'PSFPHOT'
    else if( phmethod == 'GROW' .and. .not. have_grow ) then
       write(error_unit,*) &
            "Failed to find a growth-curve photometry extension in ",trim(filename)
       goto 666
    else if( .not. have_grow .or. (saper > 0 .and. phmethod /= 'GROW') &
         .or. phmethod == 'APER' ) then
       call fits_movnam_hdu(fitsfile,FITS_BINARY_TBL,APEREXTNAME,0,status)
       write(key,'(a,i0)') FITS_COL_APCOUNT,maper
       call fits_get_colnum(fitsfile,.true.,key,ccol,status)
       write(key,'(a,i0)') FITS_COL_APCOUNTERR,maper
       call fits_get_colnum(fitsfile,.true.,key,dcol,status)
       fits%radf90 = -1
       fits%phtype = 'APHOT'
    else
       call fits_read_key(fitsfile,FITS_KEY_RF90,fits%radf90,status)
       fits%radf90 = fits%radf90 / fits%scale
       fits%aper = -1
       call fits_get_colnum(fitsfile,.true.,FITS_COL_GCOUNT,ccol,status)
       call fits_get_colnum(fitsfile,.true.,FITS_COL_GCOUNTERR,dcol,status)
       fits%phtype = 'GPHOT'
    end if

    call fits_get_colnum(fitsfile,.true.,FITS_COL_X,xcol,status)
    call fits_get_colnum(fitsfile,.true.,FITS_COL_Y,ycol,status)
    call fits_get_colnum(fitsfile,.true.,FITS_COL_SKY,scol,status)
    call fits_get_colnum(fitsfile,.true.,FITS_COL_SKYERR,ecol,status)
    if( status /= 0 ) goto 666

    call fits_get_num_rows(fitsfile,nrows,status)
    allocate(xcen(nrows),ycen(nrows))
    allocate(fits%sky(nrows),fits%dsky(nrows), &
         fits%cts(nrows),fits%dcts(nrows))

    call fits_get_rowsize(fitsfile,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fitsfile,xcol,frow,dullval,xcen(i:l),anyf,status)
       call fits_read_col(fitsfile,ycol,frow,dullval,ycen(i:l),anyf,status)
       call fits_read_col(fitsfile,scol,frow,dullval,fits%sky(i:l),anyf,status)
       call fits_read_col(fitsfile,ecol,frow,dullval,fits%dsky(i:l),anyf,status)
       call fits_read_col(fitsfile,ccol,frow,dullval,fits%cts(i:l),anyf,status)
       call fits_read_col(fitsfile,dcol,frow,dullval,fits%dcts(i:l),anyf,status)
       if( status /= 0 ) goto 666
    end do

    ! the total flux by the enclosed energy of the aperture
    if( fits%phtype == 'APHOT' .and. apee_use ) then
       call apee_apply(fitsfile,maper,fits%cts,fits%dcts,status,fits%apee, &
            fits%eesource)
       if( status /= 0 ) goto 666
    end if

    call fits_close_file(fitsfile,status)

    !  coordinates conversions
    allocate(fits%ra(nrows),fits%dec(nrows))
    call wcs % pix2cel(xcen,ycen,fits%ra,fits%dec)

    fits%nrows = nrows
    ! coordinate errors, when CRDER is undefined, HWHM is used as the matching
    ! tolerance (5*coo_err) like in mfits
    if( all(crder > 0) ) then
       fits%coo_err = maxval(crder)
    else
       fits%coo_err = (hwhm / scale) / 5
    end if

    deallocate(xcen,ycen)
    return

666 continue

    call fits_close_file(fitsfile,status)
    call fits_report_error(error_unit,status)

  end subroutine readfits


  subroutine savefits(filename,output,fits,advanced,photsys,phsyscal,&
       area,init_area,catid, qlabels,ctph, dctph, status, k0, r, lref, &
       srel, apcorr, phres, ctm, nband, scol, sdcol, csrc, sysint, nsys)

    ! ctm (optional): the colour term of the band nband, the colours
    ! of the stars scol(:,j), their errors sdcol and sources csrc (columns
    ! COLOR, COLORERR, COLORSRC are added to the table, COLORj, COLORjERR
    ! for the colours j >= 2); sysint, nsys
    ! (optional): the interval [mag] of the estimated systematic error
    ! and the number of its stars (nsys < 0: a fixed value, no keys)

    use phsysfits
    use photoconv
    use calibre
    use colterm

    character(len=*), intent(in) :: filename,output,photsys,catid
    character(len=*), dimension(:), intent(in) :: qlabels
    type(qFITS), intent(in) :: fits
    type(type_phsys), intent(in) :: phsyscal
    real(dbl), intent(in) :: area, ctph, dctph, k0, r, lref, srel, apcorr
    logical, intent(in) :: init_area, advanced
    type(photores), intent(in), optional :: phres
    integer, intent(in out) :: status
    type(colour_term), intent(in), optional :: ctm
    integer, intent(in), optional :: nband
    real(dbl), dimension(:,:), intent(in), optional :: scol, sdcol
    integer, dimension(:), intent(in), optional :: csrc
    real(dbl), dimension(2), intent(in), optional :: sysint
    integer, intent(in), optional :: nsys

    integer, parameter :: group = 1, extver = 0
    character(len=*), parameter :: bunit = units(3)
    character(len=FLEN_VALUE), dimension(:), allocatable :: ttype, tform,tunit
    character(len=FLEN_CARD), dimension(2) :: com
    real(dbl), dimension(:,:), allocatable :: q,dq
    character(len=FLEN_CARD) :: key, buf
    real(dbl) :: photflam,photzpt,photplam,photbw
    integer :: n, m, i, l, width, height, hdutype, nq, nrows, ncols, chdu, nhdu, &
         srows, frow, nt, j
    logical :: lsys
    type(fitsfiles) :: fitsfile, fitscal

    width = fits%naxes(1)
    height = fits%naxes(2)
    nq = size(qlabels)

    call fits_open_image(fitsfile,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read an image in the file `',trim(filename),"'."
       return
    end if

    call fits_create_scratch(fitscal,status)
    call fits_copy_header(fitsfile,fitscal,status)
    call fits_update_key(fitscal,'BITPIX',-32,'',status)
    !call fits_delete_key(fitscal,'BSCALE,status)
    !call fits_delete_key(fitscal,'BZERO,status)

    if( status /= 0 ) goto 666

    call fits_get_hdu_num(fitsfile,chdu)

    call fits_update_key(fitscal,FITS_KEY_PHOTSYS,photsys,'photometry filter system',status)
    call fits_update_key(fitscal,FITS_KEY_BUNIT,bunit,'Physical units of array values',status)

    if( init_area ) then
       call fits_update_key(fitscal,FITS_KEY_AREA,area,-5, &
            '[m2] Area of input aperture',status)
    end if

    if( photsys /= '' ) then
       ! compatibility with HST calibration
       call phstkeys(phsyscal,fits%filter,photflam,photzpt,photplam,photbw)

       call fits_update_key(fitscal,FITS_KEY_PHOTFLAM,photflam,5, &
            '[erg/s/cm2/A] flux for 1 photon/s/cm2',status)
       call fits_update_key(fitscal,FITS_KEY_PHOTZPT,photzpt,-6, &
            'mag zero-point for fluxes in [erg/s/cm2/A]',status)
       call fits_update_key(fitscal,FITS_KEY_PHOTPLAM,photplam,-6, &
            '[A] effective wavelength',status)
       call fits_update_key(fitscal,FITS_KEY_PHOTBW,photbw,-6,'[A] passband FWHM',status)
    end if

    call fits_write_image(fitscal,group,fits%frame,status)
    if( status == FITS_NUMERICAL_OVERFLOW ) then
       call fits_report_error(error_unit,status)
       status = 0
    end if

    if( advanced ) then
       ! copy all others extensions
       call fits_get_num_hdus(fitsfile,nhdu,status)
       do n = 1,nhdu
          if( n /= chdu ) then
             call fits_movabs_hdu(fitsfile,n,hdutype,status)
             call fits_copy_hdu(fitsfile,fitscal,0,status)
             if( status /= 0 ) goto 666
          end if
       end do
    end if
    call fits_close_file(fitsfile,status)
    call fits_report_error(error_unit,status)

    ! prepare HDU with photometry
    ncols = 6 + 2*nq
    nt = 0
    if( present(ctm) ) nt = ctm%nterm
    if( present(ctm) ) ncols = ncols + 3 + 2*(nt - 1)
    allocate(ttype(ncols),tform(ncols),tunit(ncols))

    tform = '1D'
    tunit(1:2) = 'deg'
    tunit(3:) = bunit

    ttype(1:6) = [ character(len=FLEN_VALUE) :: FITS_COL_RA, FITS_COL_DEC, &
         FITS_COL_SKY, FITS_COL_SKYERR, FITS_COL_PHOTON, FITS_COL_PHOTONERR ]

    do n = 1, nq
       i = 2*n + 5
       ttype(i) = qlabels(n)
       ttype(i+1) = trim(qlabels(n))//'ERR'
       call quantity(qlabels(n),tunit(i))
       tunit(i+1) = tunit(i)
    end do
    if( present(ctm) ) then
       i = 6 + 2*nq
       ttype(i+1:i+3) = [ character(len=FLEN_VALUE) :: 'COLOR', 'COLORERR', &
            'COLORSRC' ]
       tunit(i+1:i+2) = 'mag'
       tunit(i+3) = ''
       tform(i+3) = '1J'
       do j = 2, nt
          m = i + 3 + 2*(j - 2)
          write(ttype(m+1),'(a,i0)') 'COLOR',j
          ttype(m+2) = trim(ttype(m+1))//'ERR'
          tunit(m+1:m+2) = 'mag'
       end do
    end if

    call fits_insert_btbl(fitscal,ncols,ttype,tform,tunit,PHOTOEXTNAME,status)
    call fits_update_key(fitscal,FITS_KEY_JD,fits%jd,-18, &
         'Julian date at midpoint of exposure',status)
    call fits_update_key(fitscal,FITS_KEY_HJD,fits%hjd,-12,'Heliocentric JD',status)

    if( fits%geodefined ) then
       call fits_update_key(fitscal,FITS_KEY_AIRMASS,fits%airmass,-5,&
            'The airmass of centre at JD',status)
    else
       call fits_update_key(fitscal,FITS_KEY_AIRMASS,fits%airmass,-1,&
            'Consider the airmass as UNDEFINED',status)
    end if
    call fits_update_key(fitscal,FITS_KEY_PHOTOTYP,fits%phtype,'photometry method',status)
    call fits_update_key(fitscal,FITS_KEY_CTPH,ctph,13,'counts per photons, calibration',status)
    call fits_update_key(fitscal,FITS_KEY_CTPHERR,dctph,1,'std. error of CTPH ',status)
    call fits_update_key(fitscal,FITS_KEY_PHSYSERR,1.0857362047581294_dbl*srel,-5, &
         '[mag] systematic error added to errors',status)
    lsys = .false.
    if( present(nsys) .and. present(sysint) ) lsys = nsys >= 0
    if( lsys ) then
       if( sysint(2) >= 0 ) &
            call fits_update_key(fitscal,FITS_KEY_PHSYSERE, &
            max(sysint(2) - 1.0857362047581294_dbl*srel,0.0_dbl),-5, &
            '[mag] upper error of PHSYSERR (68%, limit if 0)',status)
       call fits_update_key(fitscal,FITS_KEY_PHSYSNST,nsys, &
            'calibration stars of PHSYSERR',status)
    end if
    if( present(ctm) ) call colkeys(fitscal,ctm,nband,status)
    if( fits%apee ) then
       call fits_update_key(fitscal,FITS_KEY_APCORTYP,'EE', &
            'aperture correction: APEE (enclosed energy)',status)
       if( fits%eesource /= '' ) &
            call fits_update_key(fitscal,'EESOURCE',fits%eesource, &
            'APEE by: PSF or AUTO (frame stars)',status)
    else if( abs(apcorr - 1) > epsilon(apcorr) ) then
       call fits_update_key(fitscal,FITS_KEY_APCORTYP,'CONST', &
            'aperture correction: the constant APCORR',status)
       call fits_update_key(fitscal,FITS_KEY_APCORR,apcorr,-6, &
            'aperture correction (total / aperture flux)',status)
    else if( fits%phtype == 'APHOT' ) then
       call fits_update_key(fitscal,FITS_KEY_APCORTYP,'NONE', &
            'aperture correction: none (aperture flux)',status)
    end if

    if( fits%aper > 0 ) &
         call fits_update_key(fitscal,FITS_KEY_APER,fits%aper,3, &
         '[deg] calibration aperture radius',status)
    if( fits%radf90 > 0 ) &
         call fits_update_key(fitscal,FITS_KEY_RF90,fits%radf90,3, &
         '[deg] radius contains 90% of flux',status)

    com(1) = '[deg] inner sky annulus radius'
    com(2) = '[deg] outer sky annulus radius'
    do i = 1, 2
       call fits_make_keyn(FITS_KEY_ANNULUS,i,key,status)
       call fits_update_key(fitscal,key,fits%annuls(i),3,com(i),status)
    end do
    call fits_update_key(fitscal,FITS_KEY_PHOTSYS,photsys,'photometry filter system',status)
    call fits_update_key(fitscal,FITS_KEY_FILTER,fits%filter,'filter in the system',status)
    if( status /= 0 ) goto 666

    if( r > 1 ) then
       call fits_write_comment(fitscal,'Extinction model k = k0*(lam/lref)**-r with parameters:',status)
       call fits_write_key(fitscal,FITS_KEY_EXTINK0,k0,-3,'mean extinction k0',status)
       call fits_write_key(fitscal,FITS_KEY_EXTINR,r,-3,'extinction power r',status)
       call fits_write_key(fitscal,FITS_KEY_EXTINREF,lref,-3,'[m] extinction reference wavelenght lref',status)
    end if

    call fits_update_key(fitscal,FITS_KEY_SKYMAG,fits%skymag,-5, &
         '[mag/arcsec2] mean sky in magnitude',status)
    call fits_update_key(fitscal,FITS_KEY_SKYMEAN,fits%rsky,-5, &
         '[ph] mean sky',status)
    call fits_update_key(fitscal,FITS_KEY_SKYSTD,fits%rskydev,-2, &
         '[ph] standard deviation of the mean sky',status)

    call fits_update_key(fitscal,FITS_KEY_CREATOR,FITS_VALUE_CREATOR,FITS_COM_CREATOR,status)
    if( status /= 0 ) goto 666

    call fits_write_comment(fitscal,BEGIN_PHOTOCAL,status)

    if( catid /= '' ) then
       call fits_write_comment(fitscal,"Reference photometric sequence: "//trim(catid),status)
    end if
    write(buf,'(a,g0.5,a,es8.1)') "Counts per photons = ",ctph," +- ",dctph
    call fits_write_comment(fitscal,buf,status)
    write(buf,'(a,f7.5,a)') "Systematic error = ",1.0857362047581294_dbl*srel, &
         " mag (added in quadrature)"
    call fits_write_comment(fitscal,buf,status)
    if( lsys ) then
       if( sysint(2) >= 0 ) then
          write(buf,'(a,f7.5,a,f8.5,a,i0,a)') "Systematic error 68% interval = ", &
               sysint(1)," .. ",sysint(2)," mag (",nsys," stars)"
       else
          write(buf,'(a,i0,a)') "Systematic error interval not determined (",nsys, &
               " stars)"
       end if
       call fits_write_comment(fitscal,buf,status)
       if( nsys < nsyswarn ) then
          if( sysint(2) >= 0 ) then
             write(buf,'(a,i0,a)') "Warning: the systematic error from ",nsys, &
                  " stars only is uncertain;"
          else
             write(buf,'(a)') "Warning: too few stars for the systematic error;"
          end if
          call fits_write_comment(fitscal,buf,status)
          call fits_write_comment(fitscal, &
               "  a fixed --syserr determined from many frames is preferable.",status)
       end if
    end if
    if( fits%apee ) then
       call fits_write_comment(fitscal, &
            "Counts divided by the enclosed energy of the aperture (APEE)",status)
    end if
    if( present(ctm) ) call colcomments(fitscal,ctm,nband,status)

    if( present(phres) ) then

       write(buf,'(a,i0)') "Objects used = ",phres%ndat
       call fits_write_comment(fitscal,buf,status)
       call fits_write_comment(fitscal,'   Catalogue RA,DEC [deg]    Photons/s/m2    Counts/s/m2       rel.dev.',status)
       do i = 1, phres%ndat
          write(buf,'(2f13.8,2en15.3,f14.4)') &
               phres%ra(i),phres%dec(i),phres%pht(i),phres%cts(i),phres%res(i)
          call fits_write_comment(fitscal,buf,status)
       end do
    end if

    call fits_write_comment(fitscal,MUNIPACK_VERSION,status)
    call fits_write_comment(fitscal,'Description: http://munipack.physics.muni.cz/dataform_photometry.html',status)
    call fits_write_comment(fitscal,END_PHOTOCAL,status)

    nrows = size(fits%ra)

    allocate(q(nrows,nq),dq(nrows,nq))
    do n = 1, nq
       call phsysconv(qlabels(n),fits%filter,phsyscal,fits%area,fits%exptime,&
            1.0_dbl,fits%ph,fits%dph,q(:,n),dq(:,n))
    end do

    call fits_get_rowsize(fitscal,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)
       call fits_write_col(fitscal,1,frow,fits%ra(i:l),status)
       call fits_write_col(fitscal,2,frow,fits%dec(i:l),status)
       call fits_write_col(fitscal,3,frow,fits%phsky(i:l),status)
       call fits_write_col(fitscal,4,frow,fits%dphsky(i:l),status)
       call fits_write_col(fitscal,5,frow,fits%ph(i:l),status)
       call fits_write_col(fitscal,6,frow,fits%dph(i:l),status)
       do n = 1,nq
          m = 2*n + 5
          call fits_write_col(fitscal,m,frow,q(i:l,n),status)
          call fits_write_col(fitscal,m+1,frow,dq(i:l,n),status)
       end do
       if( present(ctm) ) then
          m = 2*nq + 7
          call fits_write_col(fitscal,m,frow,scol(i:l,1),status)
          call fits_write_col(fitscal,m+1,frow,sdcol(i:l,1),status)
          call fits_write_col(fitscal,m+2,frow,csrc(i:l),status)
          do j = 2, nt
             m = 2*nq + 7 + 3 + 2*(j - 2)
             call fits_write_col(fitscal,m,frow,scol(i:l,j),status)
             call fits_write_col(fitscal,m+1,frow,sdcol(i:l,j),status)
          end do
       end if
    end do

    deallocate(ttype,tform,tunit,q,dq)

    ! add residuals to a next extension
    if( advanced .and. present(phres) ) then

       ncols = 5
       allocate(ttype(ncols),tform(ncols),tunit(ncols))

       tform = '1D'
       tunit(1:2) = 'deg'
       tunit(3:) = ''

       ttype(1:5) = [ character(len=FLEN_VALUE) :: FITS_COL_RA, FITS_COL_DEC, &
            FITS_COL_PHOTON, FITS_COL_COUNT, 'RESIDUALS' ]

       call fits_insert_btbl(fitscal,0,ttype,tform,tunit,EXT_PHRES,status)
       call fits_update_key(fitscal,FITS_KEY_CTPH,ctph,13,'counts per photons, calibration',status)
       call fits_update_key(fitscal,FITS_KEY_CTPHERR,dctph,1,'std. error of CTPH ',status)

       call fits_write_comment(fitscal,'This table contains relative residuals of photometry',status)
       call fits_write_comment(fitscal,'calibration in '//trim(PHOTOEXTNAME)//' extension.',status)
       write(buf,'(a,g0.5,a,es8.1)') "Counts per photons = ",ctph," +- ",dctph
       call fits_write_comment(fitscal,buf,status)

       call fits_get_rowsize(fitscal,srows,status)
       do i = 1, nrows, srows
          frow = i
          l = min(i+srows,nrows)

          call fits_write_col(fitscal,1,frow,phres%ra(i:l),status)
          call fits_write_col(fitscal,2,frow,phres%dec(i:l),status)
          call fits_write_col(fitscal,3,frow,phres%pht(i:l),status)
          call fits_write_col(fitscal,4,frow,phres%cts(i:l),status)
          call fits_write_col(fitscal,5,frow,phres%res(i:l),status)
       end do

       deallocate(ttype,tform,tunit)

    end if

    if( status == 0 ) then
       if( fits_file_exist(output) ) call fits_file_delete(output)
       call fits_file_duplicate(fitscal,output,status)
    end if

666 continue

    call fits_delete_file(fitscal,status)
    call fits_report_error(error_unit,status)

  end subroutine savefits


  subroutine colkeys(fitscal,ctm,n,status)

    ! the keywords of the colour term of the band n; the keys of
    ! the colours j >= 2 (several --color-term): COLTNTRM, COLTERMj,
    ! COLTREFj, COLTjK1, COLTjK1E, COLTjK2, COLTjK2E, COLTjCOV, COLTMEAj,
    ! COLTSTDj and COLTCjk, the correlations of the linear coefficients

    use colterm

    type(fitsfiles), intent(in out) :: fitscal
    type(colour_term), intent(in) :: ctm
    integer, intent(in) :: n
    integer, intent(in out) :: status

    character(len=FLEN_VALUE) :: c
    character(len=FLEN_KEYWORD) :: key
    character(len=FLEN_COMMENT) :: com
    integer :: j, k, nt, i1, i2

    nt = ctm%nterm
    c = trim(ctm%cola(1,n))//'-'//trim(ctm%colb(1,n))
    if( ctm%model == COLT_CAT ) then
       call fits_update_key(fitscal,'COLTMODE','CATALOGUE', &
            'colour term fitted by the catalogue colour',status)
    else
       call fits_update_key(fitscal,'COLTMODE','TRANSFORM', &
            'colour term of the transformation (tratab)',status)
    end if
    call fits_update_key(fitscal,'COLTERM',c, &
         'colour index c of the colour term',status)
    call fits_update_key(fitscal,'COLTREF',ctm%cref(1,n),-6, &
         '[mag] reference colour: CTPH of c = COLTREF',status)
    call fits_update_key(fitscal,'COLTK1',-magln*ctm%beta(1,n),-6, &
         '[mag/mag] k1 of (c - COLTREF)',status)
    call fits_update_key(fitscal,'COLTK1E',magln*sqrt(max(ctm%cov(2,2,n),0.0_dbl)),&
         -2,'[mag/mag] std. error of COLTK1',status)
    if( ctm%order(1,n) >= 2 ) then
       call fits_update_key(fitscal,'COLTK2',-magln*ctm%beta(2,n),-6, &
            '[mag/mag2] k2 of (c - COLTREF)**2',status)
       call fits_update_key(fitscal,'COLTK2E', &
            magln*sqrt(max(ctm%cov(3,3,n),0.0_dbl)),-2, &
            '[mag/mag2] std. error of COLTK2',status)
    end if
    if( ctm%model == COLT_CAT ) then
       call fits_update_key(fitscal,'COLTCOV', &
            -ctm%cov(1,2,n)/sqrt(max(ctm%cov(1,1,n)*ctm%cov(2,2,n), &
            tiny(1.0_dbl))),-4,'correlation of CTPH and COLTK1',status)
       call fits_update_key(fitscal,'COLTMEAN',ctm%cmean(1,n),-5, &
            '[mag] mean colour of the calibration stars',status)
       call fits_update_key(fitscal,'COLTSTD',ctm%cstd(1,n),-4, &
            '[mag] std. deviation of their colours',status)
       if( nt == 1 ) then
          call fits_update_key(fitscal,'COLTNSTA',ctm%nstars(n), &
               'calibration stars with a colour',status)
       else
          call fits_update_key(fitscal,'COLTNSTA',ctm%nstars(n), &
               'calibration stars with all the colours',status)
       end if
    end if
    if( nt == 1 ) return

    ! several colours: the colour j of the keys with j
    call fits_update_key(fitscal,'COLTNTRM',nt,'number of colour terms',status)
    do j = 2, nt
       c = trim(ctm%cola(j,n))//'-'//trim(ctm%colb(j,n))
       write(key,'(a,i0)') 'COLTERM',j
       write(com,'(a,i0,a,i0)') 'colour index c',j,' of the colour term ',j
       call fits_update_key(fitscal,key,c,com,status)
       write(key,'(a,i0)') 'COLTREF',j
       write(com,'(a,i0,a,i0)') '[mag] reference colour of c',j,': COLTREF',j
       call fits_update_key(fitscal,key,ctm%cref(j,n),-6,com,status)
       write(key,'(a,i0,a)') 'COLT',j,'K1'
       write(com,'(a,i0,a,i0,a)') '[mag/mag] k1 of (c',j,' - COLTREF',j,')'
       call fits_update_key(fitscal,key,-magln*ctm%beta(2*j-1,n),-6,com,status)
       write(key,'(a,i0,a)') 'COLT',j,'K1E'
       write(com,'(a,i0,a)') '[mag/mag] std. error of COLT',j,'K1'
       call fits_update_key(fitscal,key, &
            magln*sqrt(max(ctm%cov(2*j,2*j,n),0.0_dbl)),-2,com,status)
       if( ctm%order(j,n) >= 2 ) then
          write(key,'(a,i0,a)') 'COLT',j,'K2'
          write(com,'(a,i0,a,i0,a)') '[mag/mag2] k2 of (c',j,' - COLTREF',j,')**2'
          call fits_update_key(fitscal,key,-magln*ctm%beta(2*j,n),-6,com,status)
          write(key,'(a,i0,a)') 'COLT',j,'K2E'
          write(com,'(a,i0,a)') '[mag/mag2] std. error of COLT',j,'K2'
          call fits_update_key(fitscal,key, &
               magln*sqrt(max(ctm%cov(2*j+1,2*j+1,n),0.0_dbl)),-2,com,status)
       end if
       write(key,'(a,i0,a)') 'COLT',j,'COV'
       write(com,'(a,i0,a)') 'correlation of CTPH and COLT',j,'K1'
       call fits_update_key(fitscal,key, &
            -ctm%cov(1,2*j,n)/sqrt(max(ctm%cov(1,1,n)*ctm%cov(2*j,2*j,n), &
            tiny(1.0_dbl))),-4,com,status)
       write(key,'(a,i0)') 'COLTMEA',j
       write(com,'(a,i0,a)') '[mag] mean colour c',j,' of the calibration stars'
       call fits_update_key(fitscal,key,ctm%cmean(j,n),-5,com,status)
       write(key,'(a,i0)') 'COLTSTD',j
       write(com,'(a,i0,a)') '[mag] std. deviation of their colours c',j
       call fits_update_key(fitscal,key,ctm%cstd(j,n),-4,com,status)
    end do
    ! the correlations of the linear coefficients (u - g and g - r of
    ! the stars are correlated, so are their coefficients)
    do j = 1, nt - 1
       do k = j + 1, nt
          i1 = 2*j
          i2 = 2*k
          write(key,'(a,2i0)') 'COLTC',j,k
          write(com,'(a,i0,a,i0)') 'correlation of k1 of the colours ',j,' and ',k
          call fits_update_key(fitscal,key, &
               ctm%cov(i1,i2,n)/sqrt(max(ctm%cov(i1,i1,n)*ctm%cov(i2,i2,n), &
               tiny(1.0_dbl))),-4,com,status)
       end do
    end do

  end subroutine colkeys


  subroutine colcomments(fitscal,ctm,n,status)

    use colterm

    type(fitsfiles), intent(in out) :: fitscal
    type(colour_term), intent(in) :: ctm
    integer, intent(in) :: n
    integer, intent(in out) :: status

    character(len=FLEN_CARD) :: buf
    integer :: j

    if( ctm%nterm > 1 ) then
       call fits_write_comment(fitscal,"Colour terms: m_cat - m_inst = ZP + " // &
            "sum_j k1_j*(c_j - ref_j)",status)
       if( any(ctm%order(:,n) >= 2) ) &
            call fits_write_comment(fitscal,"  + k2_j*(c_j - ref_j)**2",status)
       do j = 1, ctm%nterm
          write(buf,'(a,i0,3a,f0.4,a,f0.5,a,f0.5)') "  c",j," = ", &
               trim(ctm%cola(j,n))//'-'//trim(ctm%colb(j,n)),", ref = ", &
               ctm%cref(j,n),", k1 = ",-magln*ctm%beta(2*j-1,n)," +- ", &
               magln*sqrt(max(ctm%cov(2*j,2*j,n),0.0_dbl))
          if( ctm%order(j,n) >= 2 ) then
             write(buf,'(a,a,f0.5,a,f0.5)') trim(buf),", k2 = ", &
                  -magln*ctm%beta(2*j,n)," +- ", &
                  magln*sqrt(max(ctm%cov(2*j+1,2*j+1,n),0.0_dbl))
          end if
          call fits_write_comment(fitscal,buf,status)
       end do
       call fits_write_comment(fitscal,"  applied by the catalogue colours (COLORSRC: " // &
            "bit j-1 set", status)
       call fits_write_comment(fitscal,"  for c_j); a missing colour is predicted by " // &
            "the others", status)
       call fits_write_comment(fitscal,"  (the colours of the calibration stars), " // &
            "error by its scatter.",status)
       return
    end if

    write(buf,'(a,f0.4,a)') "Colour term: m_cat - m_inst = ZP + k1*(c - ", &
         ctm%cref(1,n),")"
    if( ctm%order(1,n) >= 2 ) buf = trim(buf) // " + k2*(c - ref)**2"
    call fits_write_comment(fitscal,buf,status)
    write(buf,'(3a,f0.5,a,f0.5)') "  c = ",trim(ctm%cola(1,n))//'-'// &
         trim(ctm%colb(1,n)),", k1 = ",-magln*ctm%beta(1,n)," +- ", &
         magln*sqrt(max(ctm%cov(2,2,n),0.0_dbl))
    if( ctm%order(1,n) >= 2 ) then
       write(buf,'(a,a,f0.5,a,f0.5)') trim(buf),", k2 = ",-magln*ctm%beta(2,n), &
            " +- ",magln*sqrt(max(ctm%cov(3,3,n),0.0_dbl))
    end if
    call fits_write_comment(fitscal,buf,status)
    if( ctm%model == COLT_CAT ) then
       call fits_write_comment(fitscal,"  applied by the catalogue colour (COLORSRC = 1);", &
            status)
       call fits_write_comment(fitscal,"  the stars without it get the mean correction", &
            status)
       call fits_write_comment(fitscal,"  of the calibration stars, error by its scatter (0).",&
            status)
    else
       call fits_write_comment(fitscal,"  applied by the colours of the standard magnitudes", &
            status)
       call fits_write_comment(fitscal,"  solved from all filters (COLORSRC = 2).",status)
    end if

  end subroutine colcomments


  subroutine tra_frame(tra,pairs,ctph,xcts,xpht)

    use fotran

    integer, dimension(:,:), intent(in) :: pairs
    real(dbl), dimension(:,:), intent(in) :: tra
    real(dbl), dimension(:), intent(in) :: ctph
    real, dimension(:), intent(in) :: xcts
    real, dimension(:), intent(out) :: xpht

    real(dbl), dimension(:), allocatable :: cts,dcts,pht,dpht
    integer :: n

    n = size(xcts)
    allocate(cts(n),dcts(n),pht(n),dpht(n))
    cts = xcts * ctph
    dcts = 0
    dpht = 0
    call fotra(tra,pairs,cts,dcts,pht,dpht)
    xpht = real(pht)
    deallocate(cts,dcts,pht,dpht)

  end subroutine tra_frame


  subroutine qfits_init(fits)

    type(qFITS), intent(out) :: fits
    fits%status = .false.
    fits%filename = ''
    fits%nrows = 0
    fits%naxis = 0
    fits%naxes = 0
    fits%aper = -1
    fits%area = 1
    fits%exptime = 1
    fits%coo_err = -1
    fits%scale = 0
    fits%annuls = 0
    fits%filter = ''
    fits%rsky = 0
    fits%rskydev = 0
    fits%skymag = 99.999
    fits%radf90 = -1
    fits%airmass = 1
    fits%jd = 0.0
    fits%latitude = 0.0
    fits%longitude = 0.0

  end subroutine qfits_init

  subroutine qfits_deallocate(fits)

    type(qFITS), intent(in out) :: fits

    if( allocated(fits%frame) ) deallocate(fits%frame)
    if( allocated(fits%ra) ) deallocate(fits%ra)
    if( allocated(fits%dec) ) deallocate(fits%dec)
    if( allocated(fits%cts) ) deallocate(fits%cts)
    if( allocated(fits%dcts) ) deallocate(fits%dcts)
    if( allocated(fits%sky) ) deallocate(fits%sky)
    if( allocated(fits%dsky) ) deallocate(fits%dsky)
    if( allocated(fits%ph) ) deallocate(fits%ph)
    if( allocated(fits%dph) ) deallocate(fits%dph)
    if( allocated(fits%phsky) ) deallocate(fits%phsky)
    if( allocated(fits%dphsky) ) deallocate(fits%dphsky)

    if( allocated(fits%dmag) ) deallocate(fits%dmag)
    if( allocated(fits%mag) ) deallocate(fits%mag)

  end subroutine qfits_deallocate

  subroutine qfits_setph(fits,ph,dph,phsky,dphsky)

    type(qFITS), intent(in out) :: fits
    real(dbl), dimension(:), intent(in) :: ph,dph,phsky,dphsky
    integer :: ncat

    ncat = size(ph)
    allocate(fits%ph(ncat),fits%dph(ncat),fits%phsky(ncat),fits%dphsky(ncat))
    allocate(fits%mag(ncat),fits%dmag(ncat))
    fits%ph = ph
    fits%dph = dph
    fits%phsky = phsky
    fits%dphsky = dphsky

  end subroutine qfits_setph

end module sfits
