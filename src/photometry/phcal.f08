!
!  photometric calibration
!
!  Copyright © 2012-7, 2019-2024 F.Hroch (hroch@physics.muni.cz)
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
!  Problems:
!    * calibration frame is not referenced when a frame is referenced
!      to a frame (done)
!    * standard deviation has bad value (done)
!
!  Adds:
!    * Re-design of command line parameters: better defaults to --col-*,
!      rename all mags to be VO-compatible (Vmag,e_Vmag). User must specify
!      only filter/photometry system (or nothing when filter/phsystem is
!      presented in FITS header).
!

program phcal

  use titsio
  use phsysfits
  use mfits
  use iso_fortran_env

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)

  character(len=4*FLEN_FILENAME) :: line, key, val
  character(len=FLEN_FILENAME), dimension(:), allocatable :: filename, &
       output
  character(len=FLEN_FILENAME) :: ref = '', cat = '', tratable = '', &
       phsystable = 'photosystems.fits'
  character(len=FLEN_VALUE), dimension(:), allocatable :: qlabels
  character(len=FLEN_VALUE) :: col_ra = FITS_COL_RA, col_dec = FITS_COL_DEC, &
       photsys_ref = '', photsys_instr = ''
  character(len=FLEN_VALUE), dimension(:), allocatable :: col_mag, col_magerr, &
       filters
  logical :: verbose = .false., plog = .false., list = .false.
  logical :: init_area = .false. , init_photsys = .false.
  logical :: cal_manual = .false.
  logical :: advanced = .false.
  real(dbl) :: utol = -1.0/3600.0
  real(dbl) :: area
  real(dbl) :: threshold = 0.1
  real(dbl) :: maxerr = 0.1
  real(dbl) :: apcorr = 1.0_dbl
  ! the systematic error [mag] added to the calibrated errors,
  ! < 0: estimated by the residuals of the calibration stars
  real(dbl) :: syserr = -1
  real(dbl), dimension(:), allocatable :: ctph_manual,dctph_manual
  ! the colour terms (--color-term, repeatable, at most ncoltmax): the
  ! catalogue columns A,B of every colour A - B (and of their errors), their
  ! number (2 or 4), the orders of the polynomials and the reference
  ! colours (huge: the mean colour of the calibration stars), one value
  ! for all the colours or one for every colour
  integer, parameter :: ncoltmax = 3
  character(len=FLEN_VALUE), dimension(4,ncoltmax) :: col_color
  integer, dimension(ncoltmax) :: ncolor
  integer :: ncterm = 0
  integer, dimension(:), allocatable :: color_order
  real(dbl), dimension(:), allocatable :: color_ref
  ! the fit corrected for the errors of the catalogue colours (with
  ! the error columns; --color-err-correct)
  logical :: color_errcorr = .true.
  integer :: eq,nfile,n,nmag,nmagerr,nfilters,nctph,saper,ncol1
  character(len=FLEN_KEYWORD), dimension(8) :: keys
  integer :: nqlabels

  keys = [ character(len=FLEN_KEYWORD) :: FITS_KEY_EXPTIME, FITS_KEY_AREA, &
       FITS_KEY_DATEOBS, FITS_KEY_PHOTSYS, FITS_KEY_FILTER, &
       FITS_KEY_LONGITUDE, FITS_KEY_LATITUDE, FITS_KEY_TIMEOBS ]

!  nqlabels = 7
!  qlabels = (/'PHOTRATE','FLUX    ','FNU     ','FLAM    ','MAG     ', &
!       'ABMAG   ','STMAG   '/)

  nqlabels = 2
  allocate(qlabels(nqlabels))
  qlabels = [ character(len=FLEN_VALUE) :: 'FLUX', 'MAG' ]

  allocate(col_mag(0), col_magerr(0), filters(0))
  col_color = ''
  ncolor = 0
  allocate(color_order(1),color_ref(1))
  color_order = 1
  color_ref = huge(1.0_dbl)
  allocate(ctph_manual(0), dctph_manual(0))
  allocate(filename(0), output(0))
  nfilters = 0
  nmag = 0
  nmagerr = 0
  nfile = 0
  area = 1
  saper = 0

  do
     read(*,'(a)',end=20) line

     eq = index(line,'=')
     if( eq == 0 ) stop 'Malformed input record.'

     key = line(:eq-1)
     val = line(eq+1:)

     if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'PIPELOG' ) then

        read(val,*) plog

     else if( key == 'ADVANCED' ) then

        read(val,*) advanced

     else if( key == 'FITS_KEY_EXPTIME' ) then

        read(val,*) keys(1)

     else if( key == 'FITS_KEY_AREA' ) then

        read(val,*) keys(2)

     else if( key == 'FITS_KEY_DATEOBS' ) then

        read(val,*) keys(3)

     else if( key == 'FITS_KEY_PHOTOSYS' ) then ! remove?

        read(val,*) keys(4)

     else if( key == 'FITS_KEY_FILTER' ) then

        read(val,*) keys(5)

     else if( key == 'FITS_KEY_LONGITUDE' ) then

        read(val,*) keys(6)

     else if( key == 'FITS_KEY_LATITUDE' ) then

        read(val,*) keys(7)

     else if( key == 'FITS_KEY_TIMEOBS' ) then

        read(val,*) keys(8)

     else if( key == 'NFILTERS' ) then

        read(val,*) nfilters

        deallocate(filters)
        allocate(filters(nfilters))
        filters = ''

     else if( key == 'FILTERS' ) then

        ! a set of filters
        read(val,*) filters

     else if( key == 'COL_RA' ) then

        read(val,*) col_ra

     else if( key == 'COL_DEC' ) then

        read(val,*) col_dec

     else if( key == 'COL_NMAG' ) then

        read(val,*) nmag

        deallocate(col_mag)
        allocate(col_mag(nmag))
        col_mag = ''

     else if( key == 'COL_MAG' ) then

        ! a set columns
        read(val,*) col_mag

     else if( key == 'COL_NMAGERR' ) then

        read(val,*) nmagerr

        deallocate(col_magerr)
        allocate(col_magerr(nmagerr))
        col_magerr= ''

     else if( key == 'COL_MAGERR' ) then

        read(val,*) col_magerr

     else if( key == 'COL_NCOLOR' ) then

        ! the next colour term (every --color-term)
        read(val,*) ncol1
        if( .not. (ncol1 == 2 .or. ncol1 == 4) ) &
             stop 'Error: --color-term takes two columns (A,B) or four (A,B,errA,errB).'
        ncterm = ncterm + 1
        if( ncterm > ncoltmax ) stop 'Error: at most three colour terms (--color-term).'
        ncolor(ncterm) = ncol1

     else if( key == 'COL_COLOR' ) then

        if( ncterm == 0 ) stop 'Malformed input record (COL_COLOR).'
        read(val,*) col_color(1:ncolor(ncterm),ncterm)

     else if( key == 'COLOR_NORDER' ) then

        read(val,*) ncol1
        deallocate(color_order)
        allocate(color_order(ncol1))
        color_order = 1

     else if( key == 'COLOR_ORDER' ) then

        read(val,*) color_order
        if( any(color_order < 1 .or. color_order > 2) ) &
             stop 'Error: --color-order is 1 or 2.'

     else if( key == 'COLOR_NREF' ) then

        read(val,*) ncol1
        deallocate(color_ref)
        allocate(color_ref(ncol1))
        color_ref = huge(1.0_dbl)

     else if( key == 'COLOR_REF' ) then

        read(val,*) color_ref

     else if( key == 'COLOR_ERRCORR' ) then

        read(val,*) color_errcorr

     else if( key == 'TOL' ) then

        read(val,*) utol

     else if( key == 'SAPER' ) then

        read(val,*) saper

     else if( key == 'PHOTOMETRY' ) then

        read(val,*) phmethod
        if( .not. (phmethod == 'APER' .or. phmethod == 'GROW' .or. &
             phmethod == 'PSF') ) stop 'Error: photometry is APER, GROW or PSF.'

     else if( key == 'AREA' ) then

        read(val,*) area
        init_area = .true.

     else if( key == 'THRESHOLD' ) then

        read(val,*) threshold

     else if( key == 'MAXERR' ) then

        read(val,*) maxerr

     else if( key == 'APCORR' ) then

        read(val,*) apcorr

     else if( key == 'APEE' ) then

        ! the enclosed energy of apertures (APEE of APERPHOT) applied
        read(val,*) apee_use

     else if( key == 'SYSERR' ) then

        read(val,*) syserr

     else if( key == 'NQUANTITIES' ) then

        read(val,*) nqlabels

        deallocate(qlabels)
        allocate(qlabels(nqlabels))
        qlabels = ''

     else if( key == 'QUANTITIES' ) then

        read(val,*) qlabels

     else if( key == 'PHOTSYS_REF' ) then

        read(val,*) photsys_ref

     else if( key == 'PHOTSYS_INSTR' ) then

        read(val,*) photsys_instr
        init_photsys = .true.

     else if( key == 'PHSYSTABLE' ) then

        read(val,*) phsystable

     else if( key == 'TRATABLE' ) then

        read(val,*) tratable

     else if( key == 'LIST' ) then

        read(val,*) list

     else if( key == 'NCTPH' ) then

        read(val,*) nctph
        deallocate(ctph_manual,dctph_manual)
        allocate(ctph_manual(nctph),dctph_manual(nctph))
        ctph_manual = 1
        dctph_manual= 0

     else if( key == 'CTPH' ) then

        ! manual calibration
        read(val,*) ctph_manual
        dctph_manual = 0
        cal_manual = .true.

     else if( key == 'CAT' ) then

        ! calibration agains to the catalogue
        read(val,*) cat

     else if( key == 'REF' ) then

        ! calibration agains to the already calibrated frame
        read(val,*) ref

     else if( key == 'NFILES' ) then

        read(val,*) nfile

        deallocate(filename,output)
        allocate(filename(nfile),output(nfile))
        filename = ''
        output = ''
        n = 0

     else if( key == 'FILE' ) then

        n = n + 1
        if( n > size(filename) ) stop 'Too many files.'

        read(line(eq+1:),*) filename(n), output(n)

     end if

  end do

20 continue

  ! a constant aperture correction replaces the enclosed energy
  if( abs(apcorr - 1) > epsilon(apcorr) ) apee_use = .false.

  ! check of requested quantities before any output is created,
  ! PHRATE is an alias of PHOTRATE
  do n = 1, size(qlabels)
     if( qlabels(n) == 'PHRATE' ) qlabels(n) = 'PHOTRATE'
     if( .not. any(qlabels(n) == [ character(len=8) :: 'PHOTON', &
          'PHOTRATE', 'PHOTNU', 'PHOTLAM', 'FLUX', 'FNU', 'FLAM', 'MAG', &
          'ABMAG', 'STMAG' ]) ) then
        write(error_unit,*) "Error: unknown quantity `",trim(qlabels(n)),"'."
        stop 'An unknown quantity for calibrated data requested.'
     end if
  end do

  if( list ) then

     call listphsys(phsystable)

  else

     if( size(filename) == 0 ) stop 'No frames to process.'

     if( .not. (area > 0) ) &
          stop 'Error: The specified area is negative or zero, very funny!'

     if( ncterm > 0 .and. cat == '' ) &
          stop 'Error: --color-term needs a reference catalogue (-c).'
     if( ncterm > 0 .and. .not. (size(color_order) == 1 .or. &
          size(color_order) == ncterm) ) &
          stop 'Error: --color-order needs one value, or one for every --color-term.'
     if( ncterm > 0 .and. .not. (size(color_ref) == 1 .or. &
          size(color_ref) == ncterm) ) &
          stop 'Error: --color-ref needs one value, or one for every --color-term.'

     if( cat /= '' ) then
        call calibrate   ! catcal ?
     else if( ref /= '' ) then
        call framecal
     else if( cal_manual ) then
        call mancal
     else
        stop 'Calibration type unknown.'
     end if

  end if

  deallocate(filename,output,col_mag,col_magerr,filters,ctph_manual, &
       dctph_manual,color_order,color_ref)

  stop 0
!  write(error_unit,*) 'STOP 0'


contains

  ! calibration
  subroutine calibrate

    use phsysfits
    use calibre
    use minpacks
    use photoconv
    use jamming
    use sfits
    use fits_fotran
    use colterm

    character(len=FLEN_VALUE) :: catid,photosys_instr,photsys1,photsys2,catid1, &
         tmodel
    character(len=FLEN_VALUE), dimension(:), allocatable :: filter_frames, &
         filters_tab, colmag, tcola, tcolb
    real(dbl), allocatable, dimension(:) :: refra, refdec, ctph,dctph, &
         exptime, areas, tol, aper, sys, cra, cdec, &
         tcref, tzp, tdzp, tk1, tdk1, tk2, tdk2
    real(dbl), allocatable, dimension(:,:) :: tr,tr1,trerr,tr1err,refph,drefph,&
         cts,dcts,dn,ddn,ra,dec,mag,dmag,cmag,cdmag,magx,dmagx,sysint,ccol,dccol
    real(dbl), allocatable, dimension(:,:,:) :: ccl,dccl,ccy,ccv,cerr
    real(dbl) :: q, d
    type(type_phsys) :: phsyscal
    type(photores), dimension(:), allocatable :: phres
    type(colour_term) :: ctm
    integer, dimension(:,:), allocatable :: pairs
    integer, dimension(:), allocatable :: forder, perm, nsys
    integer :: status,nfiles,ncat,i,nf,ia,ib,j,k,l

    if( photsys_ref == '' ) &
         stop 'No identificator for reference photometry system.'

    if( nmag == 0 ) &
         stop 'Magnitude column(s) are not specified.'

    if( nmagerr > 0 .and. nmag /= nmagerr ) &
         stop 'Magnitude error column(s) does not match magnitudes itself.'

    if( nmagerr == 0 ) write(error_unit,*) &
         'Warning: Magnitude error column(s) are missing (consider include).'

    nfiles = size(filename)

    call readcat(cat,(/col_ra,col_dec/),col_mag,col_magerr, &
         refra,refdec,mag,dmag,catid,status)
    ! filters as results of readcat (?!)
    if( status /= 0 ) stop 'Failed to read a catalogue.'
    if( .not. allocated(mag) ) stop 'An empty catalogue.'

    ! the colour terms: the colours A - B of the catalogue (and their
    ! errors by the error columns, zero without them), undefined (an error
    ! < 0) when a magnitude is missing; the covariances of the colours
    ! sharing a band (u - g and g - r: -err(g)**2)
    if( ncterm > 0 ) then
       if( tratable /= '' ) stop &
            'Error: --color-term and --tratab exclude each other (the table has the colours).'
       call colour_check
       do j = 1, ncterm
          call readcat(cat,(/col_ra,col_dec/),col_color(1:2,j), &
               col_color(3:ncolor(j),j),cra,cdec,cmag,cdmag,catid1,status)
          if( status /= 0 ) stop 'Failed to read the colour columns of the catalogue.'
          if( j == 1 ) then
             allocate(ccol(size(cmag,1),ncterm),dccol(size(cmag,1),ncterm), &
                  cerr(size(cmag,1),2,ncterm),ccv(size(cmag,1),ncterm,ncterm))
             cerr = 0
          end if
          if( size(cmag,1) /= size(ccol,1) ) stop 'Colour columns of different lengths.'
          where( cmag(:,1) < 99 .and. cmag(:,2) < 99 )
             ccol(:,j) = cmag(:,1) - cmag(:,2)
             dccol(:,j) = 0
          elsewhere
             ccol(:,j) = 0
             dccol(:,j) = -1
          end where
          if( ncolor(j) == 4 ) then
             where( dccol(:,j) >= 0 .and. cdmag(:,1) < 9 .and. cdmag(:,2) < 9 )
                dccol(:,j) = sqrt(cdmag(:,1)**2 + cdmag(:,2)**2)
                cerr(:,1,j) = cdmag(:,1)
                cerr(:,2,j) = cdmag(:,2)
             end where
          end if
          deallocate(cra,cdec,cmag,cdmag)
       end do
       ccv = 0
       do j = 1, ncterm
          do k = 1, ncterm
             if( j == k ) cycle
             where( dccol(:,j) >= 0 .and. dccol(:,k) >= 0 )
                ccv(:,j,k) = shared(1,j,1,k)*cerr(:,1,j)**2 &
                     - shared(1,j,2,k)*cerr(:,1,j)**2 &
                     - shared(2,j,1,k)*cerr(:,2,j)**2 &
                     + shared(2,j,2,k)*cerr(:,2,j)**2
             end where
          end do
       end do
    end if

    ! reference system
    call phselect(phsystable,photsys_ref,phsyscal)

    call readframes(filename,keys,phsyscal%filter,photosys_instr,threshold,&
         maxerr,ra,dec,dn,ddn,utol,exptime,areas,filter_frames,tol,init_area,&
         saper,aper,verbose,status,forder)
    if( status /= 0 ) stop 'Failed to read frames.'
    if( size(filter_frames) /= size(filters) ) &
         stop 'Filters and frames dimensions differs.'
    if( init_area ) areas = area

    ! The frames are arranged by the filters of the photometry system,
    ! the catalogue columns (--col-mag, --col-magerr) follow the filters
    ! in the order given by -f: the filters and the columns are arranged
    ! as the frames, by the filter names of the frames.
    allocate(perm(nfiles))
    perm = 1
    if( nfiles > 1 ) then
       if( size(mag,2) /= size(filters) ) &
            stop 'Magnitude columns and filters dimensions differs.'
       call filters_arrange(filter_frames,filters,perm)
       filters = filters(perm)
       mag = mag(:,perm)
       dmag = dmag(:,perm)
    end if
    call phsyspairs(phsyscal,filters,pairs)
    ! the catalogue columns of the frames (arranged)
    allocate(colmag(nfiles))
    if( nfiles > 1 ) then
       colmag = col_mag(perm)
    else
       colmag = col_mag(1)
    end if

    if( verbose ) then
       write(*,*) 'Filename, filter, exptime [s], area [m2]:'
       do i = 1,size(filename)
          write(*,'(2(a,3x),2(f0.3,1x))') trim(filename(forder(i))), &
               trim(filter_frames(i)),exptime(i), areas(i)
       end do
       if( nfiles > 1 ) then
          write(*,*) 'Filter, catalogue magnitude column:'
          do i = 1,size(filters)
             write(*,'(1x,a,3x,a)') trim(filters(i)),trim(col_mag(perm(i)))
          end do
       end if
    end if
    deallocate(perm,forder)

!    write(*,*) size(mag), tol
!    write(*,*) size(cts,1)
    if( allocated(ccol) ) then
       ! the colours of all stars of the catalogue for the outputs,
       ! the colours of the calibration stars joint with the magnitudes
       call colterm_init(ctm,COLT_CAT,nfiles,ncterm)
       do j = 1, ncterm
          ctm%cola(j,:) = col_color(1,j)
          ctm%colb(j,:) = col_color(2,j)
          ctm%order(j,:) = color_order(min(j,size(color_order)))
          ctm%cref(j,:) = color_ref(min(j,size(color_ref)))
          ctm%colerr(j) = ncolor(j) == 4
       end do
       ctm%errcorr = color_errcorr
       ctm%ra = refra
       ctm%dec = refdec
       ctm%col = ccol
       ctm%dcol = dccol
       ctm%ccov = ccv
       ctm%tol = tol(1)
       nf = size(mag,2)
       l = nf + 2*ncterm + ncterm**2
       allocate(magx(size(mag,1),l),dmagx(size(mag,1),l))
       magx(:,1:nf) = mag
       do j = 1, ncterm
          magx(:,nf+2*j-1) = ccol(:,j)
          magx(:,nf+2*j) = dccol(:,j)
          do k = 1, ncterm
             magx(:,nf+2*ncterm+(j-1)*ncterm+k) = ccv(:,j,k)
          end do
       end do
       dmagx(:,1:nf) = dmag
       dmagx(:,nf+1:) = 0
       call jamcat(tol, refra, refdec, magx, dmagx, ra,dec,dn,ddn,cts,dcts)
       mag = magx(:,1:nf)
       dmag = dmagx(:,1:nf)
       deallocate(ccol,dccol,ccv)
       allocate(ccol(size(magx,1),ncterm),dccol(size(magx,1),ncterm), &
            ccv(size(magx,1),ncterm,ncterm))
       do j = 1, ncterm
          ccol(:,j) = magx(:,nf+2*j-1)
          dccol(:,j) = magx(:,nf+2*j)
          do k = 1, ncterm
             ccv(:,j,k) = magx(:,nf+2*ncterm+(j-1)*ncterm+k)
          end do
       end do
       deallocate(magx,dmagx)
    else
       call jamcat(tol, refra, refdec, mag, dmag, ra,dec,dn,ddn,cts,dcts)
    end if

    ncat = size(cts,1)
    if( ncat == 0 ) stop 'Catalogue to frame joint stars not found.'
    allocate(refph(ncat,size(mag,2)),drefph(ncat,size(mag,2)))

    if( photsys_ref /= '' ) then
       call phsysmagph(phsyscal,filters,pairs,mag(1:ncat,:),dmag(1:ncat,:), &
            refph,drefph)

       ! photon rates [ph/s/m2] to expected photons on n-th frame;
       ! columns of refph corresponds to frames (both ordered by filters)
       do n = 1,nfiles
          q = exptime(n)*areas(n)
          refph(:,n) = q*refph(:,n)
          drefph(:,n) = q*drefph(:,n)
       end do

    else
       ! undefined photosystem, relative rates
       call relmagph(mag(1:ncat,:),dmag(1:ncat,:),refph,drefph)
    end if

!    do i = 1, ncat
!       write(*,*) mag(i,1),dmag(i,1),1.086*drefph(i,1)/refph(i,1),refph(i,1)
!    end do

    tmodel = ''
    if( tratable /= '' ) then
       call tramodel(tratable,tmodel,status)
       if( status /= 0 ) &
            stop 'Instrumental to standard conversion table not found.'
    end if

    if( tmodel == 'MAG' ) then
       ! the magnitude transformation: the colour terms of every filter
       ! (fixed) by the colours of the catalogue, the counts converted
       ! by the colours of the calibrated magnitudes (writecal)
       call traload_mag(tratable,photsys2,photsys1,filters_tab,tcola,tcolb, &
            tcref,tzp,tdzp,tk1,tdk1,tk2,tdk2,status)
       if( status /= 0 ) &
            stop 'Instrumental to standard conversion table not found.'
       if( photsys2 /= photsys_ref ) &
            stop 'Reference photometry system does not match the conversion table.'
       if( photsys_instr /= '' ) photosys_instr = photsys_instr
       if( photsys1 /= photosys_instr ) then
          write(error_unit,*) "Instrumental photometry systems: `", &
               trim(photsys1),"' (table), `",trim(photosys_instr),"' (frames)."
          stop 'Instrumental photometry system does not match the conversion table.'
       end if
       if( size(filters_tab) /= size(filters) ) &
            stop "Amount of filters doesn't correspons to transformation table."
       if( any(filters_tab /= filters) ) then
          write(error_unit,*) 'Filters of the table: ', &
               (trim(filters_tab(i))//' ',i=1,size(filters_tab))
          stop 'Filters do not match the conversion table.'
       end if
       call colterm_init(ctm,COLT_TRA,nfiles)
       do n = 1, nfiles
          ia = findloc(filters,tcola(n),1)
          ib = findloc(filters,tcolb(n),1)
          if( ia == 0 .or. ib == 0 ) &
               stop 'Filters of the colours of the conversion table not found.'
          ctm%ia(n) = ia
          ctm%ib(n) = ib
          ctm%cola(1,n) = tcola(n)
          ctm%colb(1,n) = tcolb(n)
          ctm%cref(1,n) = tcref(n)
          ctm%beta(1,n) = -tk1(n) / magln
          ctm%beta(2,n) = -tk2(n) / magln
          ctm%cov(2,2,n) = (tdk1(n) / magln)**2
          ctm%cov(3,3,n) = (tdk2(n) / magln)**2
          ctm%order(1,n) = merge(2,1,tk2(n) /= 0)
       end do
       ctm%tol = tol(1)
       deallocate(filters_tab,tcola,tcolb,tcref,tzp,tdzp,tk1,tdk1,tk2,tdk2)
       allocate(tr(size(refph,2),size(refph,2)),tr1(size(refph,2),size(refph,2)), &
            trerr(size(refph,2),size(refph,2)),tr1err(size(refph,2),size(refph,2)))
       tr = 0.0_dbl
       forall(n=1:size(refph,2))
          tr(n,n) = 1.0_dbl
       end forall
       tr1 = tr
    else if( tratable/= '' ) then
       call traload(tratable,photsys2,photsys1,tr,trerr,tr1,tr1err,status, &
            filters_tab)
       if( status /= 0 ) &
            stop 'Instrumental to standard conversion table not found.'
       if( photsys2 /= photsys_ref ) &
            stop 'Reference photometry system does not match the conversion table.'
       ! the instrumental system: --photsys-instr, or PHOTSYS of the frames
       if( photsys_instr /= '' ) photosys_instr = photsys_instr
       if( photsys1 /= photosys_instr ) then
          write(error_unit,*) "Instrumental photometry systems: `", &
               trim(photsys1),"' (table), `",trim(photosys_instr),"' (frames)."
          stop 'Instrumental photometry system does not match the conversion table.'
       end if
       if( size(tr,1) /= size(filters) .or. size(tr,2) /= size(filters) ) &
            stop "Amount of filters doesn't correspons to transformation table."
       ! both the filters and the table in the order of the photometry system
       if( any(filters_tab /= filters) ) then
          write(error_unit,*) 'Filters of the table: ', &
               (trim(filters_tab(i))//' ',i=1,size(filters_tab))
          stop 'Filters do not match the conversion table.'
       end if
       deallocate(filters_tab)
       ! The counts are converted to photons by the inverse of the forward
       ! matrix (photons to counts) that is used for the calibration, so
       ! a calibration star gets back its photons. The backward matrix of
       ! the table is fitted independently with the same (tridiagonal)
       ! structure, while the inverse of a tridiagonal matrix is full:
       ! 0.2-0.3% colour-dependent errors with off-diagonal terms of 5-13%.
       call qrinv(tr,tr1)
    else
       allocate(tr(size(refph,2),size(refph,2)),tr1(size(refph,2),size(refph,2)), &
            trerr(size(refph,2),size(refph,2)),tr1err(size(refph,2),size(refph,2)))
       tr = 0.0_dbl
       forall(n=1:size(refph,2))
          tr(n,n) = 1.0_dbl
       end forall
       tr1 = tr
    end if

    if( verbose ) then
       write(*,*) 'Photons to counts transformation:'
       do i = 1,size(filters)
          write(*,'(1x,a,6f8.3)') trim(filters(i)),tr(i,:)
       end do
       write(*,*) 'Counts to photons transformation:'
       do i = 1,size(filters)
          write(*,'(1x,a,6f8.3)') trim(filters(i)),tr1(i,:)
       end do
!       write(*,*) 'Photon and count rates are in [*/s/m2].'
    end if

!    do i = 1, ncat
!       write(*,*) real(refph(i,1)),real(cts(i,1))
!    end do

!    do i = 1, size(cts,1)
!       write(*,*) dcts(i,1), cts(i,1), dcts(i,1)**2-cts(i,1)
!    end do


    ! The reference photon flux is computed by magnitudes,
    ! or by flux normalised per both second and square meter.
    ! Statistical errors of flux can not be easy scaled
    ! on arbitraty exposure time or telescope diameter,
    ! if ones are mix of Poisson and Normal distributions
    ! of deviations.
    !
    ! On the other side, our photons, carefully selected
    ! for bright stars, has nearly Poisson's rates and
    ! errors can be properly scaled.

    do n = 1,0*nfiles
       q = exptime(n)*areas(n)
       dcts(1:ncat,n) = max(dcts(1:ncat,n)**2 - cts(1:ncat,n),0.0)
       dcts(1:ncat,n) = sqrt(cts(1:ncat,n) / q + dcts(1:ncat,n) / q**2)
       cts(1:ncat,n) = cts(1:ncat,n) / q
       dcts(1:ncat,n) = dcts(1:ncat,n) / sqrt(q)
    end do


!    block
!      real :: s
!      s = 0
!      open(1,file='n')
!      do n = 1, size(dn,1)
!         if( dn(n,1) > 0 ) then
!            s = s + sqrt(dn(n,1))/dn(n,1)
!            write(1,*) n, sqrt(dn(n,1))/sqrt(dn(1,1)),s/(sqrt(dn(1,1))/dn(1,1))
!         end if
!      end do
!      close(1)
!    end block

    do i = 1, size(cts,1)

 !      write(*,*) (dcts(i,1)**2 - cts(i,1))/cts(i,1), (drefph(i,1)**2 - refph(i,1))/refph(i,1)/30
    end do

!    if( index(cat,'cone') > 0 ) then
!       do i = 1, size(cts,1)
!          !       dcts(i,1) = sqrt(cts(i,1))
!          drefph(i,1) = sqrt(refph(i,1) + drefph(i,1) / 10)
!       end do
!    end if

    allocate(phres(nfiles),ctph(nfiles),dctph(nfiles),sys(nfiles), &
         sysint(2,nfiles),nsys(nfiles))
    do n = 1, nfiles
       call photores_init(phres(n),ncat)
       phres(n)%ra = refra(1:ncat)
       phres(n)%dec = refdec(1:ncat)
    enddo

    if( ctm%model /= COLT_NONE ) then
       ! the colours of the calibration stars, their errors, the covariance
       ! of ln(ph) and the colour (the colour of a band containing
       ! the catalogue magnitude of the band: -dmag**2/magln for A,
       ! +dmag**2/magln for B, with the errors of the colour only)
       allocate(ccl(ncat,nfiles,ctm%nterm),dccl(ncat,nfiles,ctm%nterm), &
            ccy(ncat,nfiles,ctm%nterm))
       ccy = 0
       do n = 1, nfiles
          if( ctm%model == COLT_CAT ) then
             do j = 1, ncterm
                ccl(:,n,j) = ccol(1:ncat,j)
                dccl(:,n,j) = dccol(1:ncat,j)
                if( ncolor(j) == 4 ) then
                   d = 0
                   if( col_color(1,j) == colmag(n) ) d = d - 1
                   if( col_color(2,j) == colmag(n) ) d = d + 1
                   where( dccl(:,n,j) >= 0 .and. dmag(1:ncat,n) < 9 )
                      ccy(:,n,j) = d * dmag(1:ncat,n)**2 / magln
                   end where
                end if
             end do
          else
             ia = ctm%ia(n)
             ib = ctm%ib(n)
             ccl(:,n,1) = mag(1:ncat,ia) - mag(1:ncat,ib)
             where( dmag(1:ncat,ia) < 9 .and. dmag(1:ncat,ib) < 9 )
                dccl(:,n,1) = sqrt(dmag(1:ncat,ia)**2 + dmag(1:ncat,ib)**2)
             elsewhere
                dccl(:,n,1) = 0
             end where
             d = 0
             if( ia == n ) d = d - 1
             if( ib == n ) d = d + 1
             where( dmag(1:ncat,ia) < 9 .and. dmag(1:ncat,ib) < 9 )
                ccy(:,n,1) = d * dmag(1:ncat,n)**2 / magln
             end where
          end if
       end do
       if( ctm%model == COLT_CAT ) then
          call calibr_colour(refph,drefph,cts(1:ncat,:),dcts(1:ncat,:),ccl,dccl, &
               ccy,ctm,ctph,dctph,phres,verbose,syserr,sys,sysint,nsys, &
               ccv=ccv(1:ncat,:,:))
       else
          call calibr_colour(refph,drefph,cts(1:ncat,:),dcts(1:ncat,:),ccl,dccl, &
               ccy,ctm,ctph,dctph,phres,verbose,syserr,sys,sysint,nsys)
       end if
       deallocate(ccl,dccl,ccy)
    else if( nfiles == 1 ) then
       call calibr(refph,drefph,cts(1:ncat,:),dcts(1:ncat,:),ctph,dctph, &
            phres,verbose,syserr,sys,sysint,nsys)
    else
       call caliber(pairs,filters,tr,tr1,refph,drefph, &
            cts(1:ncat,:),dcts(1:ncat,:),ctph,dctph,phres,verbose,syserr,sys, &
            sysint,nsys)
    end if
    call syscheck(filters,sys,sysint,nsys)

    ! apperture correction
    ctph = apcorr * ctph
    dctph = apcorr * dctph

    if( verbose ) then
       if( nfiles > 1 ) then
          write(*,*) 'Relative residuals:'
          write(*,'(5x,tr5,5a10)') filters
          do i = 1, ncat
             write(*,'(i5,5f10.4)') i,(phres(n)%res(i),n=1,size(phres))
          end do
       end if
       write(*,*) 'Final solutions:'
       do i = 1,size(ctph)
          write(*,'(3a,g15.5,a,1pg8.1)') 'ctph(',trim(filters(i)),') =', &
               ctph(i),'+-',dctph(i)
       end do
    end if

    ! per second and square meter
    q = exptime(1)*areas(1)
    d = dctph(1) / ctph(1)
!    ctph = q / ctph
!    dctph = ctph * d

    if( ctm%model /= COLT_NONE ) then
       call writecal(filename,output,keys,advanced,phsystable,filters, &
            catid,photsys_ref,area,init_area,saper,tr1,ctph,dctph,qlabels,phres, &
            verbose,sys,apcorr,ctm,sysint=sysint,nsys=nsys)
    else
       call writecal(filename,output,keys,advanced,phsystable,filters, &
            catid,photsys_ref,area,init_area,saper,tr1,ctph,dctph,qlabels,phres, &
            verbose,sys,apcorr,sysint=sysint,nsys=nsys)
    end if

    do i = 1,size(phres)
       call photores_destroy(phres(i))
    end do
    deallocate(phres)

    if( allocated(pairs) ) deallocate(pairs)

    deallocate(refra,refdec,refph,drefph,filter_frames,exptime,areas,aper, &
         tol,tr,tr1,trerr,tr1err,cts,dcts,dn,ddn,ra,dec,mag,dmag,ctph,dctph,sys, &
         sysint,nsys, &
         colmag)

  end subroutine calibrate

  subroutine syscheck(filters,sys,sysint,nsys)

    ! The systematic error estimated from few stars is uncertain: a warning
    ! (the header gets it as a comment). A fixed --syserr: nsys = -1
    ! (no interval, no number of stars in the header).

    use calibre, only: nsyswarn

    character(len=*), dimension(:), intent(in) :: filters
    real(dbl), dimension(:), intent(in) :: sys
    real(dbl), dimension(:,:), intent(in) :: sysint
    integer, dimension(:), intent(in out) :: nsys

    integer :: i

    if( syserr >= 0 ) then
       nsys = -1
       return
    end if

    do i = 1, size(nsys)
       if( nsys(i) < nsyswarn ) then
          if( sysint(2,i) >= 0 ) then
             write(error_unit,'(3a,i0,a,f6.4,a,f6.4,a,f6.4,a)',advance='no') &
                  ' Warning: the systematic error of `',trim(filters(i)), &
                  "' is estimated from ",nsys(i),' stars only: ',sys(i), &
                  ' mag (68% interval ',sysint(1,i),' .. ',sysint(2,i),' mag)'
          else
             write(error_unit,'(3a,i0,a,f6.4,a)',advance='no') &
                  ' Warning: the systematic error of `',trim(filters(i)), &
                  "' cannot be determined from ",nsys(i),' stars (',sys(i),' mag)'
          end if
          write(error_unit,'(a)') &
               '; consider a fixed --syserr determined from many frames.'
       end if
    end do

  end subroutine syscheck

  subroutine relmagph(mag,dmag,ph,dph)

    use photoconv

    real(dbl), dimension(:,:), intent(in) :: mag,dmag
    real(dbl), dimension(:,:), intent(out) :: ph,dph

    do n = 1, size(mag,2)
       call mag2rate(mag(:,n),dmag(:,n),ph(:,n),dph(:,n))
!       dph(:,n) = sqrt(ph(:,n))
    end do


  end subroutine relmagph


  subroutine framecal

    ! framecal suppose calibration in single filter only

    use jamming
    use calibre
    use sfits

    character(len=FLEN_VALUE) :: catid,photosys_frames
    character(len=FLEN_VALUE), dimension(:), allocatable :: filter_ref, &
         filters_fram
    real(dbl), allocatable, dimension(:) :: refra, refdec, ctph,dctph, &
         exptime, areas, tol, aper, sys
    real(dbl), allocatable, dimension(:,:) :: refph,drefph, &
         cts,dcts,ra,dec,sysint
    integer, allocatable, dimension(:) :: nsys
    real(dbl), dimension(1,1), parameter :: tratab = real(1.0,dbl)
    type(photores), dimension(:), allocatable :: phres
    type(type_phsys) :: phsyscal
    integer :: status,nfiles,ncat,i,n
    real(dbl) :: reftime, refarea, refaper, q

    nfiles = size(filename)
    status = 0

    nfilters = 1
    allocate(filters_fram(1),filter_ref(1))

    ! reference system
    call phselect(phsystable,photsys_ref,phsyscal)

    call readref(ref,keys,maxerr,threshold,refra,refdec,refph,drefph, &
         filter_ref(1),reftime,refarea,refaper,status)
    if( status /= 0 ) stop 'Failed to read a reference frame.'
    catid = trim(ref)

    call readframes(filename,keys,phsyscal%filter,photosys_frames,threshold,&
         maxerr,ra,dec,cts,dcts,utol,exptime,areas,filters_fram,tol,init_area,&
         saper,aper,verbose,status)
    if( status /= 0 ) stop 'Failed to read frames.'
    if( init_area ) areas = area
    if( any(abs(refaper - aper) > 0.1) ) &
         write(error_unit,*) 'Warning: synthetic apertures looks different.'

    ! Rates of the reference frame are rescaled on calibrated frames;
    ! the reference frames are supposed as the best frame, or the frame
    ! with its exposure longer than the exposure of calibrated ones.
    q = (reftime * refarea) / (exptime(1)*areas(1))
    !    write(*,*) q,sqrt(q)
!    q = reftime * refarea
!    q = exptime(1) * areas(1)
    !    q = 1
    refph = refph / q
!    drefph = drefph / sqrt(q) ! pure Poisson
    drefph = drefph / q    ! Normal-like

!    drefph = max(drefph**2 - refph,0.0)
!    drefph = sqrt(refph / q + drefph/q**2)
!    refph = refph / q

!    do i = 1, size(refph,1)
!       write(*,*) sqrt(refph(i,1)),drefph(i,1)
!    end do

    if( verbose ) then
       write(*,*) 'Filename, filter, exptime [s], area [m2]:'
       write(*,'(2(a,3x),2(f0.3,1x),a)') trim(ref), &
               trim(filter_ref(1)),reftime,refarea,'   (reference)'
       do i = 1,size(filename)
          write(*,'(2(a,3x),2(f0.3,1x))') trim(filename(i)), &
               trim(filters_fram(i)),exptime(i),areas(i)
       end do
    end if


!!$    block
!!$      real(dbl), allocatable, dimension(:) :: refra0, refdec0
!!$      real(dbl), allocatable, dimension(:,:) :: refph0,drefph0, &
!!$           cts0,dcts0,ra0,dec0
!!$      integer :: m
!!$
!!$      allocate(ctph(1),dctph(1),phres(1))
!!$      allocate(refra0,source=refra)
!!$      allocate(refdec0,source=refdec)
!!$      allocate(refph0,source=refph)
!!$      allocate(drefph0,source=drefph)
!!$      allocate(cts0,source=cts)
!!$      allocate(dcts0,source=dcts)
!!$      allocate(ra0,source=ra)
!!$      allocate(dec0,source=dec)
!!$
!!$      open(1,file='n')
!!$
!!$      do n = 1, size(refra0)
!!$
!!$         deallocate(refra,refdec,refph,drefph)
!!$         allocate(refra(n),refdec(n),refph(n,1),drefph(n,1))
!!$
!!$         refra = refra0(:n)
!!$         refdec = refdec0(:n)
!!$         refph(:,1) = refph0(:n,1)
!!$         drefph(:,1) = drefph0(:n,1)
!!$         ra = ra0
!!$         dec = dec0
!!$         cts = cts0
!!$         dcts = dcts0
!!$         call jamref(tol, refra, refdec, refph, drefph, ra,dec,cts,dcts)
!!$         m = size(cts,1)
!!$         call calibr(refph(1:m,:),drefph(1:m,:),cts(1:m,:),dcts(1:m,:),&
!!$              ctph,dctph,phres,verbose)
!!$         write(1,*) n,ctph,dctph,(cts(1,1)/ctph(1)-refph(1,1))/refph(1,1)
!!$
!!$      end do
!!$      close(1)
!!$    end block
!!$
!!$


    call jamref(tol, refra, refdec, refph, drefph, ra,dec,cts,dcts)

    ncat = size(cts,1)
    if( ncat == 0 ) stop 'Frame to frame joint stars not found.'

!    do n = 1,nfiles
!       q = (exptime(n) * areas(n)) / (reftime * refarea)
!       q = 1
!       cts(1:ncat,n) = cts(1:ncat,n) / q
!       dcts(1:ncat,n) = dcts(1:ncat,n) / sqrt(q)
!    end do

    allocate(phres(nfiles),ctph(nfiles),dctph(nfiles),sys(nfiles), &
         sysint(2,nfiles),nsys(nfiles))
    do n = 1, nfiles
       call photores_init(phres(n),ncat)
       phres(n)%ra = refra(1:ncat)
       phres(n)%dec = refdec(1:ncat)
    enddo

!    write(*,*) refph(1:ncat,:)
!    write(*,*) cts(1:ncat,:)

    call calibr(refph(1:ncat,:),drefph(1:ncat,:),cts(1:ncat,:),dcts(1:ncat,:),&
         ctph,dctph,phres,verbose,syserr,sys,sysint,nsys)
    call syscheck(filters_fram,sys,sysint,nsys)

    ! pull back of original rates
!    q = reftime * refarea
!    q = (reftime * refarea) / (exptime(1)*areas(1))
!    ctph = q * ctph
!    dctph = q * dctph

    ! apperture correction
    ctph = apcorr * ctph
    dctph = apcorr * dctph

    ! update rates
    q = (reftime * refarea) / (exptime(1)*areas(1))
    !    q = exptime(1)*areas(1)
    q = 1
    do i = 1, size(ctph)
       phres(i)%pht = phres(i)%pht / q
       phres(i)%cts = phres(i)%cts / q
!       phres(i)%pht = phres(i)%pht / (reftime * refarea)
!       phres(i)%cts = phres(i)%cts / (exptime(1)*areas(1))
    end do

    ! per second and square meter
!    q = (reftime * refarea) / (exptime(1)*areas(1))
    !q = exptime(1)*areas(1)
!    d = dctph(1) / ctph(1)
!    ctph = q * ctph
!    dctph = d * ctph

    call writecal(filename,output,keys,advanced,phsystable, &
         filters_fram,catid,photsys_ref,area,init_area,saper,tratab, &
         ctph,dctph,qlabels,phres,verbose,sys,apcorr,sysint=sysint,nsys=nsys)

    deallocate(refra,refdec,refph,drefph,filter_ref,filters_fram,exptime, &
         areas,tol,cts,dcts,ra,dec,ctph,dctph,sys,sysint,nsys)

    do i = 1,size(phres)
       call photores_destroy(phres(i))
    end do

    deallocate(phres)

  end subroutine framecal

  ! manual calibration
  subroutine mancal

    use sfits
    use fits_fotran
    use minpacks

    character(len=FLEN_VALUE) :: catid,photsys1,photsys2
    character(len=FLEN_VALUE), dimension(:), allocatable :: filter_frames
    real(dbl), allocatable, dimension(:,:) :: tr,tr1,trerr,tr1err
    type(type_phsys) :: phsyscal
    integer, dimension(:), allocatable :: order, perm
    integer :: n,status

    ! the calibrations (-C) follow the filters (-f) in the order given
    ! by the user, writecal arranges the frames by the filters of
    ! the photometry system: both are arranged by the filter names
    if( size(filename) > 1 .and. size(filters) == size(filename) .and. &
         size(ctph_manual) == size(filters) ) then
       call phselect(phsystable,photsys_ref,phsyscal)
       call file_order(filename,keys,phsyscal%filter,order,filters=filter_frames)
       allocate(perm(size(filters)))
       call filters_arrange(filter_frames,filters,perm)
       filters = filters(perm)
       ctph_manual = ctph_manual(perm)
       dctph_manual = dctph_manual(perm)
       deallocate(order,perm,filter_frames)
    end if

    if( tratable/= '' ) then
       call tramodel(tratable,photsys1,status)
       if( status == 0 .and. photsys1 == 'MAG' ) stop &
            'Error: a magnitude transformation (MAG) needs the catalogue calibration (-c).'
       call traload(tratable,photsys2,photsys1,tr,trerr,tr1,tr1err,status)
       if( status /= 0 ) &
            stop 'Instrumental to standard conversion table not found.'
       if( size(tr,1) /= size(filters) ) &
            stop "Dimensions of transformation table and filters doesn't corresponds."
       if( photsys2 /= photsys_ref ) &
            stop 'Reference photometry system does not match the conversion table.'
       if( photsys1 /= photsys_instr ) &
            stop 'Instrumental photometry system does not match the conversion table.'
       ! the counts to photons by the inverse of the forward matrix
       ! (as in calibrate)
       call qrinv(tr,tr1)

    else

       allocate(tr(size(filename),size(filename)),tr1(size(filename),size(filename)), &
            trerr(0,0),tr1err(0,0))
       tr = 0.0_dbl
       forall(n=1:size(tr,1))
          tr(n,n) = 1.0_dbl
       end forall
       tr1 = tr
    end if

    catid = ''

    call writecal(filename,output,keys,advanced,phsystable,filters, &
         catid,photsys_ref,area,init_area,saper,tr1,apcorr*ctph_manual, &
         apcorr*dctph_manual,qlabels,verbose=verbose, &
         syserr=spread(max(syserr,0.0_dbl),1,size(ctph_manual)),apcorr=apcorr)

    deallocate(tr,tr1,trerr,tr1err)

  end subroutine mancal

  subroutine colour_check

    ! the colours of the terms must be independent: two colours of
    ! the same bands, or three of three bands (u - g, g - r, u - r) make
    ! the fit singular

    character(len=FLEN_VALUE), dimension(2*ncoltmax) :: bands
    real(dbl), dimension(ncoltmax,2*ncoltmax) :: a
    real(dbl), dimension(ncoltmax,ncoltmax) :: gm
    real(dbl) :: det
    integer :: nb, j, k

    nb = 0
    do j = 1, ncterm
       if( col_color(1,j) == col_color(2,j) ) &
            stop 'Error: --color-term A,B with A = B.'
       do k = 1, 2
          if( .not. any(bands(1:nb) == col_color(k,j)) ) then
             nb = nb + 1
             bands(nb) = col_color(k,j)
          end if
       end do
    end do
    a = 0
    do j = 1, ncterm
       a(j,findloc(bands(1:nb),col_color(1,j),1)) = 1
       a(j,findloc(bands(1:nb),col_color(2,j),1)) = -1
    end do
    ! the Gram determinant (an integer) is zero for dependent colours
    gm = 0
    gm(1:ncterm,1:ncterm) = matmul(a(1:ncterm,1:nb),transpose(a(1:ncterm,1:nb)))
    select case (ncterm)
    case (1)
       det = gm(1,1)
    case (2)
       det = gm(1,1)*gm(2,2) - gm(1,2)*gm(2,1)
    case default
       det = gm(1,1)*(gm(2,2)*gm(3,3) - gm(2,3)*gm(3,2)) &
            - gm(1,2)*(gm(2,1)*gm(3,3) - gm(2,3)*gm(3,1)) &
            + gm(1,3)*(gm(2,1)*gm(3,2) - gm(2,2)*gm(3,1))
    end select
    if( .not. (abs(det) > 0.5_dbl) ) &
         stop 'Error: the colours of --color-term are not independent.'

  end subroutine colour_check

  real(dbl) function shared(a,j,b,k)
    ! 1 when the band a (1: A, 2: B) of the colour j is the band b of
    ! the colour k, else 0
    integer, intent(in) :: a, j, b, k
    shared = 0
    if( col_color(a,j) == col_color(b,k) ) shared = 1
  end function shared

end program phcal
