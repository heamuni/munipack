!
! Kombine     Combine a set of frames.
!
! Copyright © 1998-2026 F.Hroch (hroch@physics.muni.cz)
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

program kombine

  use fitskombi
  use oakleaf
  use astrotrafo
  use titsio
  use iso_fortran_env
  use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN, IEEE_IS_NAN,&
       IEEE_POSITIVE_INF, IEEE_NEGATIVE_INF

  implicit none

  character(len=FLEN_FILENAME) :: filename
  character(len=4*FLEN_FILENAME) :: key,val, record
  character(len=80) :: msg

  character(len=FLEN_KEYWORD), dimension(7) :: fitskeys
  character(len=FLEN_FILENAME) :: output = 'kombine.fits'

  ! spherical sprojection
  logical :: init_sproj = .false.  ! initialised by user
  logical :: init_crval = .false.  ! initialised by user
  logical :: init_crpix = .false.
  logical :: init_scale=.false.
  logical :: init_angle=.false.
  logical :: init_reflex = .false.
  logical :: init_jdref = .false.
  character(len=FLEN_VALUE) :: ptype = "TAN"
  real(REAL64) :: scale = 3600      ! pix per deg
  real(REAL64) :: angle = 0         ! rotation in degrees
  logical :: reflex = .false.
  real(REAL64) :: jdref = 0         ! proper motion reference JD
  real(REAL64), dimension(2) :: crval, crpix
  real(REAL64), dimension(2) :: crmov = 0
  real(REAL64) :: longitude, latitude
  logical :: geodefined = .false.

  ! information keywords
  character(len=FLEN_VALUE) :: object = ''
  character(len=FLEN_VALUE) :: filter = ''
  character(len=FLEN_VALUE) :: dateobs = ''
  integer :: bitpix = -32        ! default output bitpix
  integer, dimension(2) :: naxes = -1 ! the code is restricted on 2D frames
  integer :: nfiles = 0          ! input files count
  real :: exptime = 1            ! mean exposure time
  logical :: robust = .true.     ! output image is computed by robust mean
  ! the robust mean: CLIP (the mean clipped by the noise model, default)
  ! or HAMPEL (the M-estimator of OakLeaf, biased low on a few photons)
  character(len=FLEN_VALUE) :: meanmode = 'CLIP'
  logical :: background = .true. ! subtract background
  logical :: fluxcons = .false.  ! scale by pixel areas to conserve total flux
  logical :: verbose = .false.   ! verbose output
  logical :: pipelog = .false.

  ! selection of method of interpolation:
  !
  !  NEAR(default)  BILINEAR    BICUBIC   BI3CONV
  !  fast           slow        slower    the slowest
  !  pixelized      good        perfect   for undersampled
  !
  ! warning: all BI-* methods can give undesired fringing patterns on the background
  character(len=10) :: method = 'NEAR'

  integer :: eq, stat, status
  type(KombiFits), dimension(:), allocatable :: files
  type(AstroWCS) :: wcs
  real, dimension(:,:), allocatable :: IMAGE, EXPMASK, IMGERR, IMGDEV

  ! mandatory keys
  fitskeys = [ character(len=FLEN_KEYWORD) :: FITS_KEY_DATEOBS, FITS_KEY_TIMEOBS, &
       FITS_KEY_EXPTIME, FITS_KEY_FILTER, FITS_KEY_OBJECT, FITS_KEY_LONGITUDE, &
       FITS_KEY_LATITUDE ]

  do
     read(*,'(a)',iostat=stat,iomsg=msg) record
     if( stat == IOSTAT_END ) exit
     if( stat > 0 ) then
        write(error_unit,*) trim(msg)
        error stop 'An input error.'
     end if

     eq = index(record,'=')
     if( eq == 0 ) stop 'An input record malformed.'

     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'OUTPUT' ) then

        read(val,*) output

     else if( key == 'FITS_KEY_OBJECT' ) then

        read(val,*) fitskeys(5)

     else if( key == 'FITS_KEY_FILTER' ) then

        read(val,*) fitskeys(4)

     else if( key == 'FITS_KEY_EXPTIME' ) then

        read(val,*) fitskeys(3)

     else if( key == 'FITS_KEY_DATEOBS' ) then

        read(val,*) fitskeys(1)

     else if( key == 'FITS_KEY_TIMEOBS' ) then

        read(val,*) fitskeys(2)

     else if( key == 'FITS_KEY_LONGITUDE' ) then

        read(val,*) fitskeys(6)

     else if( key == 'FITS_KEY_LATITUDE' ) then

        read(val,*) fitskeys(7)

     else if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'PIPELOG' ) then

        read(val,*) pipelog

     else if( key == 'BITPIX' ) then

        read(val,*) bitpix

     else if( key == 'SCALE' ) then

        read(val,*) scale
        init_scale = .true.

     else if( key == 'ANGLE' ) then

        read(val,*) angle
        init_angle = .true.

     else if( key == 'REFLEX' ) then

        read(val,*) reflex
        init_reflex = .true.

     else if( key == 'CRVAL' ) then

        read(val,*) crval
        init_crval = .true.

     else if( key == 'CRMOV' ) then

        read(val,*) crmov

     else if( key == 'JDREF' ) then

        read(val,*) jdref
        init_jdref = .true.

     else if( key == 'CRPIX' ) then

        read(val,*) crpix
        init_crpix = .true.

     else if( key == 'NAXES' ) then

        read(val,*) naxes

     else if( key == 'PROJECTION' ) then

        read(val,*) ptype
        init_sproj = .true.

     else if( key == 'INTERPOL' ) then

        read(val,*) method
        if( .not. (method == 'NEAR' .or. method == 'BILINEAR' .or. &
             method == 'BICUBIC' .or. method == 'BI3CONV') ) &
             stop 'Specified interpolation method is not implemented yet.'

     else if( key == 'ROBUST' ) then

        read(val,*) robust

     else if( key == 'MEAN' ) then

        read(val,*) meanmode
        if( .not. (meanmode == 'CLIP' .or. meanmode == 'HAMPEL') ) &
             stop 'The robust mean is CLIP or HAMPEL.'

     else if( key == 'BACKGROUND' ) then

        read(val,*) background

     else if( key == 'FLUXCONSERVE' ) then

        read(val,*) fluxcons

     else if( key == 'NFILES' ) then

        read(val,*) nfiles
        allocate(files(nfiles))
        nfiles = 0

     else if( key == 'FILE' ) then

        read(val,*) filename

        nfiles = nfiles + 1
        call files(nfiles) % Load(filename,fitskeys,background,verbose)
        if( .not. files(nfiles) % status ) then
           nfiles = nfiles - 1
           write(error_unit,*) &
                "Error: file `",trim(filename),"' read failed. Ignored."
        end if

     end if

  end do ! read

  if( nfiles == 0 ) stop 'No input images.'

  ! Section: initialisation -----------------------------------
  ! Unspecified parameters are derived from the first frame
  block
    real(REAL64) :: jd, jdmin
    integer :: n

    exptime = 1

    if( dateobs == '' ) then
       jdmin = huge(jdmin)
       do n = 1, nfiles
          jd = files(n)%midjd()
          if( jd < jdmin ) then
             jdmin = jd
             dateobs = files(n)%dateobs
          end if
       end do
    end if

    if( filter == '' ) then
       filter = files(1)%filter
       do n = 1, nfiles
          if( files(n)%filter /= filter ) then
             write(error_unit,*) "Warning: the filter `",trim(files(n)%filter), &
                  "' of `",trim(files(n)%filename),"' does not match --",&
                  "`",trim(filter),"' due `",trim(files(1)%filename),"' is assumed."
          end if
       end do
    end if

    if( object == '' ) then
       object = files(1)%object
       do n = 1, nfiles
          if( files(n)%object /= object ) then
             write(error_unit,*) "Warning: the object `",trim(files(n)%object), &
                  "'does not match for `",trim(files(n)%filename), &
                  "' (the object is `",trim(object),"' by `", &
                  trim(files(1)%filename),"')"
          end if
       end do
    end if

    if( files(1)%geo ) then
       geodefined = .true.
       longitude = files(1)%longitude
       latitude = files(1)%latitude
    end if

    if( any(naxes <= 0) ) &
       naxes = files(1)%naxes

    block
      character(len=FLEN_VALUE) :: xcode
      real(REAL64), dimension(2) :: xcrval, xcrpix
      real(REAL64) :: xscale, xrot, xrefl

      ! the centre of the output frame is the centre of the first frame
      ! (its CRVAL belongs to its CRPIX, which needn't be the centre)
      xcode = files(1)%code
      xcrpix = (naxes + 1) / 2.0_REAL64
      call files(1)%wcs%pix2cel((files(1)%naxes(1) + 1) / 2.0_REAL64, &
           (files(1)%naxes(2) + 1) / 2.0_REAL64, xcrval(1), xcrval(2))
      xscale = files(1)%scale
      xrot = files(1)%rot
      xrefl = files(1)%refl

      if( init_sproj ) xcode = ptype
      if( init_crval ) xcrval = crval
      if( init_crpix ) xcrpix = crpix
      if( init_scale ) xscale = scale
      if( init_angle ) xrot = angle
      if( init_reflex )then
         if( reflex ) then
            xrefl = -1
         else
            xrefl = 1
         end if
      end if

      call wcs % init(xcode,xcrval,xcrpix,xscale,xrot,xrefl)
    end block

    ! reference time for proper motion
    if( .not. init_jdref ) jdref = files(1)%midjd()

  end block ! init

  ! time to create the output image
  allocate(IMAGE(naxes(1),naxes(2)),IMGERR(naxes(1),naxes(2)), &
       EXPMASK(naxes(1),naxes(2)),stat=stat,errmsg=msg)
  if( stat /= 0 ) then
     write(error_unit,*) "Error: ",trim(msg)
     stop 'Failed to allocate memory.'
  end if
  IMAGE = 0
  EXPMASK = 0

  ! Section: processing ---------------------------------------------
  block

    real, dimension(:,:,:), allocatable :: bigimg
    logical, dimension(:,:), allocatable :: mask, fill
    real(REAL64), dimension(:,:), allocatable :: alpha, delta
    real(REAL64), dimension(:), allocatable :: xgrid, ygrid
    real, dimension(:,:), allocatable :: img
    real, dimension(:), allocatable :: fbuf, wexp, xgain, xrn2, xoff, wbuf
    real(REAL64), dimension(:,:), allocatable :: s0, s1, s2
    integer, dimension(:,:), allocatable :: ncov
    logical, dimension(:), allocatable :: isntnan
    real :: amean, stddev
    integer :: i,j,n,m, idev

    allocate(alpha(naxes(1),naxes(2)),delta(naxes(1),naxes(2)),fbuf(nfiles), &
         img(naxes(1),naxes(2)),mask(naxes(1),naxes(2)),isntnan(nfiles), &
         fill(naxes(1),naxes(2)),xgrid(naxes(1)),ygrid(naxes(1)),&
         stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Failed to allocate an image.'
    end if

    ! determine nodes of sky coordinate grid,
    ! the nodes are centres of pixels: (i,j) in FITS (WCSLIB) convention
    do i = 1, naxes(1)
       xgrid(i) = real(i,REAL64)
    end do
    do j = 1, naxes(2)
       ygrid = real(j,REAL64)
       call wcs % pix2cel(xgrid,ygrid,alpha(:,j),delta(:,j))
    end do

    ! compute mean exposure time
    fbuf = files(:)%exptime
    call rmean(fbuf,exptime)

    ! mask of touched pixels, value .true. is for pixels which needs to be
    ! additionaly set with image mean
    fill = .true.

    ! The mean of the exposure-normalised values v = img / w, w the
    ! exposure relative to the mean one, weighted by w (the arithmetic
    ! mean of the rates), and its standard error by the scatter of the
    ! values: c = sum w (v - mean)**2 / (m - 1), err**2 = c / sum w (the
    ! error of the mean for m equal exposures), NaN for m < 2 (mosaics:
    ! m is the coverage of the pixel).
    if( robust ) then
       ! compute the average by robust mean

       ! allocate memory
       allocate(bigimg(naxes(1),naxes(2),nfiles),IMGDEV(naxes(1),naxes(2)), &
            stat=stat,errmsg=msg)
       if( stat /= 0 ) then
          write(error_unit,*) trim(msg)
          error stop 'Failed to allocate so huge image.'
       end if

       IMGERR = IEEE_VALUE(IMGERR, IEEE_QUIET_NAN) !needs appropriate handling
       IMGDEV = 0

       do n = 1, nfiles
          if( files(n)%loadimg() ) then
             call files(n)%reproject(method,jdref,crmov,alpha,delta,img,mask,fluxcons)
             where( mask )
                bigimg(:,:,n) = img / (files(n)%exptime / exptime) ! area?
             elsewhere
                bigimg(:,:,n) = IEEE_VALUE(bigimg(:,:,n), IEEE_QUIET_NAN)  ! =nan
             end where
             fill = fill .and. .not. mask
             deallocate(files(n)%image)
          end if
       end do

       ! the exposure weights and the noise model of the frames
       allocate(wexp(nfiles),xgain(nfiles),xrn2(nfiles),xoff(nfiles),wbuf(nfiles))
       do n = 1, nfiles
          wexp(n) = files(n)%exptime / exptime
          xgain(n) = files(n)%gain
          xrn2(n) = files(n)%rn2
          ! the raw level (before the background subtraction) above the bias
          xoff(n) = files(n)%back - files(n)%bias
       end do

       do i = 1, naxes(1)
          do j = 1, naxes(2)
             if( .not. fill(i,j) ) then
                isntnan = .not. IEEE_IS_NAN(bigimg(i,j,:))
                m = count(isntnan)
                fbuf(1:m) = pack(bigimg(i,j,:),isntnan)
                if( meanmode == 'HAMPEL' ) then
                   call rmean(fbuf(1:m),IMAGE(i,j),IMGERR(i,j),stddev)
                   if( m < 2 ) IMGERR(i,j) = IEEE_VALUE(IMGERR(i,j),IEEE_QUIET_NAN)
                else
                   wbuf(1:m) = pack(wexp,isntnan)
                   call clipped_mean(fbuf(1:m),wbuf(1:m),pack(xgain,isntnan), &
                        pack(xrn2,isntnan),pack(xoff,isntnan),IMAGE(i,j), &
                        IMGERR(i,j),stddev)
                end if
                fbuf(1:m) = fbuf(1:m) - IMAGE(i,j)
                idev = maxloc(abs(fbuf(1:m)),1)
                if( stddev > 0 ) then
                   IMGDEV(i,j) = fbuf(idev) / stddev
                else
                   if( fbuf(idev) > 0 ) then
                      IMGDEV(i,j) = IEEE_VALUE(IMGDEV(i,j),IEEE_POSITIVE_INF)
                   else
                      IMGDEV(i,j) = IEEE_VALUE(IMGDEV(i,j),IEEE_NEGATIVE_INF)
                   end if
                end if
             end if
          end do
       end do

    else

       ! compute by the arithmetic mean: the sums of the weights w, of
       ! img = w v, of img**2 / w = w v**2 and the count of the values
       allocate(s0(naxes(1),naxes(2)),s1(naxes(1),naxes(2)), &
            s2(naxes(1),naxes(2)),ncov(naxes(1),naxes(2)),stat=stat,errmsg=msg)
       if( stat /= 0 ) then
          write(error_unit,*) trim(msg)
          error stop 'Failed to allocate so huge image.'
       end if
       s0 = 0
       s1 = 0
       s2 = 0
       ncov = 0
       do n = 1, nfiles
          if( files(n)%loadimg() ) then
             call files(n)%reproject(method,jdref,crmov,alpha,delta,img,mask,fluxcons)
             where( mask )
                s0 = s0 + (files(n)%exptime / exptime)
                s1 = s1 + img
                s2 = s2 + real(img,REAL64)**2 / (files(n)%exptime / exptime)
                ncov = ncov + 1
             end where
             fill = fill .and. .not. mask
             deallocate(files(n)%image)
          end if
       end do

       IMAGE = 0
       IMGERR = IEEE_VALUE(IMGERR, IEEE_QUIET_NAN)
       where( s0 > 0 ) IMAGE = real(s1 / s0)
       where( ncov > 1 .and. s0 > 0 ) &
            IMGERR = real(sqrt(max(s2 - s1**2/s0,0.0_REAL64) / (ncov - 1) / s0))
       EXPMASK = real(s0)
       deallocate(s0,s1,s2,ncov)
    end if

    ! untouched pixels are filled by a mean of the touched ones
    n = max(minval(naxes) / 1000, 1)
    if( any(.not. fill(::n,::n)) ) then
       call rmean(pack(IMAGE(::n,::n),.not. fill(::n,::n)),amean)
    else
       amean = 0
    end if
    where( fill ) IMAGE = amean

    ! scale on proper values: the sum of the frames (the mean times their
    ! number, for every pixel, also of mosaics); the error too
    IMAGE = nfiles * IMAGE
    IMGERR = nfiles * IMGERR

  end block

  ! Section: FITS save ----------------------------------------------
  block
    integer :: n,naxis
    character(len=FLEN_COMMENT) :: buf, code, radesys
    character(len=FLEN_VALUE), dimension(2) :: ctype
    real(REAL64), dimension(2) :: crval,crpix,crder
    real(REAL64), dimension(2,2) :: cd
    real(REAL64) :: totexp, scale, rot, refl, epoch
    real, parameter :: minvalue = 0
    real :: maxvalue
    type(fitsfiles) :: fits
    type(AstrometryGeometry) :: tg

    status = 0
    call fits_create_scratch(fits,status)
    naxis = size(naxes)
    call fits_insert_image(fits,bitpix,naxis,naxes,status)

    if( object /= '' ) then
       call fits_update_key(fits,fitskeys(5),object,'Object',status)
    end if
    if( filter /= '' ) then
       call fits_update_key(fits,fitskeys(4),filter,'Filter',status)
    end if

    totexp = sum(files(:)%exptime)
    call fits_update_key(fits,fitskeys(3),totexp,-3,'Total exposure time',status)

    if( dateobs /= '' ) then
       call fits_update_key(fits,fitskeys(1),dateobs,'The start of whole set',status)
    end if

    if( geodefined ) then
       call fits_update_key(fits,fitskeys(6),longitude,-5, &
            '[deg] geographic longitude (-east)',status)
       call fits_update_key(fits,fitskeys(7),latitude,-5, &
            '[deg] geographic latitude (+north)',status)
    end if

    ! print output info
    if( verbose ) then
       write(error_unit,*) 'Output image: ',trim(output)
       write(error_unit,*) 'Dimension: ',naxes(1),'x',naxes(2)
       write(error_unit,*) 'Object: ',trim(object)
       write(error_unit,*) 'Filter: ',trim(filter)
       write(error_unit,*) 'Date of observation: ',trim(dateobs)
       write(error_unit,*) 'Total exp. time: ',totexp
    end if

    call wcs % get_code(code)
    ctype = [ 'RA---'//code(1:3), 'DEC--'//code(1:3) ]
    call wcs % get_crval(crval)
    call wcs % get_crpix(crpix)
    call wcs % get_affines(scale,rot,refl)
    call wcs % get_radesys(radesys)
    call wcs % get_epoch(epoch)
    call tg % init(scale,rot,refl)
    call tg % get_fcd(cd)

    call fits_update_wcs(fits,radesys,epoch,ctype,crval,crpix,cd,crder,status)

    if( ptype == '' ) then
       call fits_write_errmark
       call fits_delete_key(fits,'CTYPE1',status)
       call fits_delete_key(fits,'CTYPE2',status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
       end if
    end if

    write(buf,'(a,i0,a)') &
         'The image is a collection of ',nfiles,' exposure(s):'
    call fits_write_comment(fits,buf,status)
    do n = 1, nfiles
       call fits_write_comment(fits,"'"//trim(files(n)%filename)//"'",status)
    enddo

    if( robust .and. meanmode == 'HAMPEL' ) then
       buf = 'robust (Hampel)'
    else if( robust ) then
       buf = 'clipped'
    else
       buf = 'arithmetic'
    end if
    call fits_write_comment(fits, &
         'The averadge by '//trim(buf)//' mean.',status)
    call fits_write_comment(fits,'STDERR: the error of IMAGE (the sum).',status)

    if( method == 'NEAR' ) then
       buf = 'Nearest-neighbour'
    else if( method == 'BILINEAR' ) then
       buf = 'Bi-linear'
    else if( method == 'BICUBIC' ) then
       buf = 'Bicubic'
    else if( method == 'BI3CONV' ) then
       buf = 'Bicubic convolution '
    else
       buf = 'Unspecified'
    end if
    call fits_write_comment(fits,trim(buf) // &
         ' method of interpolation has been used.',status)

    ! range cut-off
    if( bitpix > 0 )then
       maxvalue = 2.0**bitpix - 1.0
       IMAGE = max(minvalue,min(maxvalue,IMAGE))
    end if

    call fits_write_image(fits,0,IMAGE,status)

    ! stderr
    call fits_insert_image(fits,bitpix,naxis,naxes,status)
    call fits_update_key(fits,'EXTNAME',EXT_STDERR,'',status)
    call fits_write_comment(fits,'Standard errors',status)
    call fits_write_image(fits,0,IMGERR,status)

    ! maximal deviations for the robust mean
    if( robust ) then
       call fits_insert_image(fits,bitpix,naxis,naxes,status)
       call fits_update_key(fits,'EXTNAME',EXT_MAXDEV,'',status)
       call fits_write_comment(fits,'Maximal deviations in sigmas',status)
       call fits_write_image(fits,0,IMGDEV,status)
    end if

    if( status == 0 ) then
       if( fits_file_exist(output) ) call fits_file_delete(output)
       call fits_file_duplicate(fits,output,status)
    end if
    call fits_delete_file(fits,status)
    call fits_report_error(error_unit,status)

  end block

  if( status == 0 )then
     stop 0
  else
     stop 'Kombine failed.'
  end if

contains

  subroutine clipped_mean(v,w,gain,rn2,off,mean,err,sdev)

    ! The mean of the values v of a pixel (exposure-normalised, weights
    ! w the relative exposures) clipped by the noise model, as the sky
    ! of aphot (skymean.f08): the interval
    !
    !    t - k*s  <=  v  <=  t + q + sqrt(q**2 + (k*s)**2),
    !
    ! iterated from the median, t the current mean, k = 5, s the larger
    ! of the observed scatter of the used values and the one of a value
    ! by the noise model, sqrt((t*w + off)/gain + rn2) / w [ADU] (off
    ! the raw level of the frame above the bias: the subtracted background
    ! less BIASLEV), q = k**2/(6*gain*w) one photo-electron times k**2/6
    ! (the Bernstein bound of the Poisson tail: no photon of a sky of a
    ! fraction of an electron is clipped, while the robust means of such
    ! values are biased low). Cosmic rays, hot pixels and satellites, many
    ! s above, are rejected. Without the gain, s is the scatter, q = k**2/6.
    ! err is the standard error of the weighted mean (as the arithmetic
    ! mean), sdev the scatter of the used values.

    real, dimension(:), intent(in) :: v, w, gain, rn2, off
    real, intent(out) :: mean, err, sdev

    integer, parameter :: maxit = 50
    real(REAL64), parameter :: k = 5
    real(REAL64), dimension(size(v)) :: s, q, lo, hi, d
    logical, dimension(size(v)) :: keep, prev
    real(REAL64) :: t, sw, so, c
    integer :: m, n, it

    n = size(v)
    if( n == 0 ) then
       mean = 0
       err = IEEE_VALUE(err,IEEE_QUIET_NAN)
       sdev = 0
       return
    end if

    keep = .true.
    if( n >= 3 ) then
       ! the median and the scatter by the median absolute deviation
       t = vmedian(real(v,REAL64))
       d = abs(v - t)
       so = 1.4826_REAL64*vmedian(d)
       where( gain > 0 )
          q = k**2 / (6*gain*w)
       elsewhere
          q = k**2 / 6
       end where
       do it = 1, maxit
          ! the scatter of a value by the noise model
          where( gain > 0 )
             s = sqrt(max(t*w + off,0.0_REAL64)/gain + rn2) / w
          elsewhere
             s = 0
          end where
          s = max(so,s)
          lo = t - k*s
          hi = t + q + sqrt(q**2 + (k*s)**2)
          prev = keep
          keep = lo <= v .and. v <= hi
          m = count(keep)
          if( m < 2 ) then
             keep = prev
             exit
          end if
          t = sum(w*v,keep) / sum(w,keep)
          so = sqrt(sum(w*(v - t)**2,keep) / sum(w,keep) * m / (m - 1))
          if( it > 1 .and. all(keep .eqv. prev) ) exit
       end do
    end if

    m = count(keep)
    sw = sum(w,keep)
    t = sum(w*v,keep) / sw
    mean = real(t)
    if( m > 1 ) then
       c = sum(w*(v - t)**2,keep) / (m - 1)
       err = real(sqrt(c / sw))
       sdev = real(sqrt(sum((v - t)**2,keep) / (m - 1)))
    else
       err = IEEE_VALUE(err,IEEE_QUIET_NAN)
       sdev = 0
    end if

  end subroutine clipped_mean

  real(REAL64) function vmedian(a)

    ! the median of a small array (insertion sort of a copy)

    real(REAL64), dimension(:), intent(in) :: a
    real(REAL64), dimension(size(a)) :: b
    real(REAL64) :: x
    integer :: i, j, n

    n = size(a)
    b = a
    do i = 2, n
       x = b(i)
       j = i - 1
       do while( j >= 1 )
          if( b(j) <= x ) exit
          b(j+1) = b(j)
          j = j - 1
       end do
       b(j+1) = x
    end do
    if( mod(n,2) == 1 ) then
       vmedian = b(n/2+1)
    else
       vmedian = (b(n/2) + b(n/2+1)) / 2
    end if

  end function vmedian

end program kombine
