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
  real(dbl), dimension(:), allocatable :: ctph_manual,dctph_manual
  integer :: eq,nfile,n,nmag,nmagerr,nfilters,nctph,saper
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

     else if( key == 'TOL' ) then

        read(val,*) utol

     else if( key == 'SAPER' ) then

        read(val,*) saper

     else if( key == 'AREA' ) then

        read(val,*) area
        init_area = .true.

     else if( key == 'THRESHOLD' ) then

        read(val,*) threshold

     else if( key == 'MAXERR' ) then

        read(val,*) maxerr

     else if( key == 'APCORR' ) then

        read(val,*) apcorr

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
       dctph_manual)

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

    character(len=FLEN_VALUE) :: catid,photosys_instr,photsys1,photsys2
    character(len=FLEN_VALUE), dimension(:), allocatable :: filter_frames
    real(dbl), allocatable, dimension(:) :: refra, refdec, ctph,dctph, &
         exptime, areas, tol, aper
    real(dbl), allocatable, dimension(:,:) :: tr,tr1,trerr,tr1err,refph,drefph,&
         cts,dcts,dn,ddn,ra,dec,mag,dmag
    real(dbl) :: q, d
    type(type_phsys) :: phsyscal
    type(photores), dimension(:), allocatable :: phres
    integer, dimension(:,:), allocatable :: pairs
    integer :: status,nfiles,ncat,i

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

    ! reference system
    call phselect(phsystable,photsys_ref,phsyscal)
    call phsyspairs(phsyscal,filters,pairs)

    call readframes(filename,keys,phsyscal%filter,photosys_instr,threshold,&
         maxerr,ra,dec,dn,ddn,utol,exptime,areas,filter_frames,tol,init_area,&
         saper,aper,verbose,status)
    if( status /= 0 ) stop 'Failed to read frames.'
    if( size(filter_frames) /= size(filters) ) &
         stop 'Filters and frames dimensions differs.'
    if( init_area ) areas = area

    if( verbose ) then
       write(*,*) 'Filename, filter, exptime [s], area [m2]:'
       do i = 1,size(filename)
          write(*,'(2(a,3x),2(f0.3,1x))') trim(filename(i)), &
               trim(filter_frames(i)),exptime(i), areas(i)
       end do
    end if

!    write(*,*) size(mag), tol
!    write(*,*) size(cts,1)
    call jamcat(tol, refra, refdec, mag, dmag, ra,dec,dn,ddn,cts,dcts)

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

    if( tratable/= '' ) then
       call traload(tratable,photsys2,photsys1,tr,trerr,tr1,tr1err,status)
       if( status /= 0 ) &
            stop 'Instrumental to standard conversion table not found.'
       if( photsys2 /= photsys_ref ) &
            stop 'Reference photometry system does not match the conversion table.'
       if( photsys1 /= photsys_instr ) &
            stop 'Instrumental photometry system does not match the conversion table.'
       if( size(tr,1) /= size(filters) ) &
            stop "Amount of filters doesn't correspons to transformation table."
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

    allocate(phres(nfiles),ctph(nfiles),dctph(nfiles))
    do n = 1, nfiles
       call photores_init(phres(n),ncat)
       phres(n)%ra = refra(1:ncat)
       phres(n)%dec = refdec(1:ncat)
    enddo

    if( nfiles == 1 ) then
       call calibr(refph,drefph,cts(1:ncat,:),dcts(1:ncat,:),ctph,dctph, &
            phres,verbose)
    else
       call caliber(pairs,filters,tr,tr1,refph,drefph, &
            cts(1:ncat,:),dcts(1:ncat,:),ctph,dctph,phres,verbose)
    end if

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

    call writecal(filename,output,keys,advanced,phsystable,filters, &
         catid,photsys_ref,area,init_area,saper,tr1,ctph,dctph,qlabels,phres,verbose)

    do i = 1,size(phres)
       call photores_destroy(phres(i))
    end do
    deallocate(phres)

    if( allocated(pairs) ) deallocate(pairs)

    deallocate(refra,refdec,refph,drefph,filter_frames,exptime,areas,aper, &
         tol,tr,tr1,trerr,tr1err,cts,dcts,dn,ddn,ra,dec,mag,dmag,ctph,dctph)

  end subroutine calibrate

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
         exptime, areas, tol, aper
    real(dbl), allocatable, dimension(:,:) :: refph,drefph, &
         cts,dcts,ra,dec
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

    allocate(phres(nfiles),ctph(nfiles),dctph(nfiles))
    do n = 1, nfiles
       call photores_init(phres(n),ncat)
       phres(n)%ra = refra(1:ncat)
       phres(n)%dec = refdec(1:ncat)
    enddo

!    write(*,*) refph(1:ncat,:)
!    write(*,*) cts(1:ncat,:)

    call calibr(refph(1:ncat,:),drefph(1:ncat,:),cts(1:ncat,:),dcts(1:ncat,:),&
         ctph,dctph,phres,verbose)

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
         ctph,dctph,qlabels,phres,verbose)

    deallocate(refra,refdec,refph,drefph,filter_ref,filters_fram,exptime, &
         areas,tol,cts,dcts,ra,dec,ctph,dctph)

    do i = 1,size(phres)
       call photores_destroy(phres(i))
    end do

    deallocate(phres)

  end subroutine framecal

  ! manual calibration
  subroutine mancal

    use sfits
    use fits_fotran

    character(len=FLEN_VALUE) :: catid,photsys1,photsys2
    real(dbl), allocatable, dimension(:,:) :: tr,tr1,trerr,tr1err
    integer :: n,status

    if( tratable/= '' ) then
       call traload(tratable,photsys2,photsys1,tr,trerr,tr1,tr1err,status)
       if( status /= 0 ) &
            stop 'Instrumental to standard conversion table not found.'
       if( size(tr,1) /= size(filters) ) &
            stop "Dimensions of transformation table and filters doesn't corresponds."
       if( photsys2 /= photsys_ref ) &
            stop 'Reference photometry system does not match the conversion table.'
       if( photsys1 /= photsys_instr ) &
            stop 'Instrumental photometry system does not match the conversion table.'

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
         apcorr*dctph_manual,qlabels,verbose=verbose)

    deallocate(tr,tr1,trerr,tr1err)

  end subroutine mancal

end program phcal
