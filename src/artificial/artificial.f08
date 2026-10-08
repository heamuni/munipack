!
!  artificial frames
!
!  Copyright © 2016-2019, 2022-2026 F.Hroch (hroch@physics.muni.cz)
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
!  * Checking proper combination of parameters during computations (to prevent
!    frames below horizon with atmosphere, negative values of radii, etc.)?
!  * Add galaxies (profiles as NFW, Hubble, King, spirals) or globular clusters
!    with Plummer? Fractals for nebulae?
!

program artificial

  use spray
  use generator
  use artgalaxy
  use flatmodule
  use fitskeys
  use lcio
  use phsysfits
  use trajd
  use titsio
  use iso_fortran_env

  implicit none

  real(REAL64), parameter :: pi = 4*atan(1.0_REAL64)
  real(REAL64), parameter :: J2000 = 2451544.5_REAL64
  character(len=*), parameter :: artificial_fits = 'artificial_???.fits'

  character(len=4*FLEN_FILENAME) :: record,key,val
  character(len=FLEN_FILENAME) :: cat = '',  phsystable = '', fitskeyfile = ''

  character(len=FLEN_FILENAME), codimension[*] :: psffile = ''
  character(len=FLEN_FILENAME), codimension[*] :: galcat = ''
  character(len=FLEN_FILENAME), codimension[*] :: flatmodelfile = ''
  character(len=3), codimension[*] :: prjcode = 'TAN'
  logical :: verbose = .false.
  integer, dimension(2) :: naxes = [ 666, 666 ]
  real(REAL64), dimension(2) :: crpix = [ 333, 333 ]
  integer :: bitpix = -32
  real(REAL64) :: gain = 1
  real(REAL64) :: saturate = 10.0_REAL32**int(log10(huge(1.0_REAL32)))
  real(REAL64) :: skymag = 21  ! mag per 1 sq. arcsec (at zenit)
  real(REAL64), dimension(2) :: skygrad = [ 0,0 ]  ! mag per 1 sq. arcsec per pix
  integer :: maglim = 13
  real(REAL64) :: exptime = 1 ! in seconds
  real(REAL64) :: area = 1  ! 1 sec per 1 sq.meter, equiv. r=0.56m
  real(REAL64) :: hwhm = 1  ! in pixels
  real(REAL64) :: eccentricity = 0 ! 0 (circle) .. 1 (line)
  real(REAL64) :: inclination = 0  ! inclination of major semiaxis, degs
  real(REAL64) :: airy = 0.1   ! Airy disk radius in pix for r=0.56m
  real(REAL64) :: beta = 2     ! Moffat's exponent
  real(REAL64) :: scale = 3600.0  ! 1 arcsex = 1 pix
  real(REAL64) :: fov = 666. / 3600.0  ! field of view, used for scale setup
  real(REAL64), dimension(2) :: crval = 0 ! coordinates of centre [deg]
  real(REAL64) :: angle = 0  ! a rotation around the center
  real(REAL64) :: qefficiency = 1  ! quantum efficiency
  real(REAL64) :: gnoise = 0       ! Gaussian noise instead of skymag
  logical      :: noise = .true.   ! data has added noise componets
  real(REAL64) :: extink = 0.1     ! extinction coefficient
  real(REAL64) :: longitude = 0    ! geographic longitude
  real(REAL64) :: latitude = 0     ! geographic latitude
  logical      :: atmosphere = .false.
  real(REAL64) :: cosmics_freq = 0                 ! cosmics ray frequency in [s]
  real(REAL64) :: cosmics_cmax = 1e3  ! counts in maximum of cosmics ray
  real(REAL64) :: cosmics_len = 10    ! cosmics ray lengths
  integer      :: cosmics_nevents = 0 ! Po(freq*exptime)
  real(REAL64) :: lcmag = 15                ! lc magnitude
  real(REAL64) :: lcamp = 1                 ! lc amplitude
  real(REAL64) :: lcjd0 = J2000             ! lc reference time
  real(REAL64) :: lcper = 1                 ! lc period
  real(REAL64), dimension(2) :: lcoo = 0    ! lc coordinates
  character(len=FLEN_VALUE) :: lcmodel = '', lctable = '', lcfourier = ''
  character(len=FLEN_VALUE), codimension[*] :: filter = 'V'
  character(len=FLEN_VALUE), codimension[*] :: profile = 'SEEING'
  character(len=FLEN_VALUE), codimension[*] :: spread = 'AUTO'
  ! 'FFT' for FFT convolution, 'RANDOM' for random offset, 'AUTO' for 'RANDOM'
  ! of short exposures and 'TEST' for direct convolution
  character(len=FLEN_KEYWORD), dimension(:), allocatable :: keys_include
  character(len=FLEN_CARD), dimension(:), allocatable :: cards_include
  integer :: mcount_keys = 0
  integer :: nstars_coo
  real(REAL64), dimension(:), allocatable :: xcoo,ycoo,flux
  real, dimension(:,:), allocatable :: flatmodel
  integer, dimension(2) :: flatmodel_naxes = 0
  type(type_phsys) :: phsys

  character(len=FLEN_FILENAME), codimension[*] :: filemask = ''
  integer, codimension[*] :: nfiles = 1
  real(REAL64), codimension[*] :: jd0 = J2000         ! 2000-01-01T00:00:00
  real(REAL64), codimension[*] :: timestep = 60       ! time step in seconds
  type(event_type), codimension[*] :: events
  integer, codimension[*] :: mcaf = 0

  logical :: lcmode = .false.
  logical :: datetime_init = .false.
  logical :: lccoo_init = .false.
  logical :: fov_defined = .false.
  logical :: area_defined = .false.
  logical :: airy_defined = .false.
  logical :: diameter_defined = .false.

  integer :: eq,n,m,i,istat
  character(len=81) :: errmsg
  integer :: year = 2000, month = 1, day = 1, hour = 0, minute = 0
  real(REAL64) :: sec = 0
  real(REAL64) :: radii = 0.56              ! radii
  real(REAL64) :: diameter = 2*0.56         ! diameter
  character(len=FLEN_VALUE) :: phsystem = 'Johnson'
  character(len=FLEN_VALUE), dimension(5) :: labels = ''

  type(ArtGalaxies) :: galaxies
  integer :: ngal = 0

  real(REAL64), dimension(:,:,:), allocatable :: mpsf
  real(REAL64), dimension(:), allocatable :: rcoo, acoo
  real(REAL64) :: rgrid, agrid
  integer :: nbox, zoom
  integer, dimension(3) :: mpsf_shape = 0

  real(REAL64), dimension(:), allocatable :: lctime, lclux, lcdux
  complex(REAL64), dimension(:), allocatable :: lcfcoeff
  integer :: nlc_four = 0, nlc_tab = 0

  if( this_image() == 1 ) then

     do
        read(*,'(a)',iostat=istat,iomsg=errmsg) record
        if( is_iostat_end(istat) ) exit
        if( istat /= 0 ) then
           write(error_unit,*) 'Error: ',trim(errmsg)
           stop 'An std. input error in artificial.'
        end if
        eq = index(record,'=')
        if( eq == 0 ) stop 'Malformed input record.'
        key = record(:eq-1)
        val = record(eq+1:)

        if( key == 'VERBOSE' ) then

           read(val,*) verbose

!        else if( key == 'PIPELOG' ) then

!           read(val,*) plog

        else if( key == 'PSFFILE' ) then

           read(val,*) psffile

        else if( key == 'FITSKEYFILE' ) then

           read(val,*) fitskeyfile

        else if( key == 'FLATMODELFILE' ) then

           read(val,*) flatmodelfile

        else if( key == 'PROFILE' ) then

           read(val,*) profile
           if( .not. any(profile == ['SEEING','MOFFAT','GAUSS ']) ) then
              write(error_unit,*) "Error: unknown PSF `",trim(profile), &
                   "' (use SEEING, MOFFAT or GAUSS)."
              stop 'Error: unknown PSF.'
           end if

        else if( key == 'SPREAD' ) then

           read(val,*) spread
           if( .not. any(spread == ['AUTO  ','FFT   ','RANDOM','TEST  ']) ) then
              write(error_unit,*) "Error: unknown spread method `", &
                   trim(spread),"' (use AUTO, FFT or RANDOM)."
              stop 'Error: unknown spread method.'
           end if

        else if( key == 'NAXES' ) then

           read(val,*) naxes

        else if( key == 'BITPIX' ) then

           read(val,*) bitpix

        else if( key == 'GAIN' ) then

           read(val,*) gain

        else if( key == 'SATURATE' ) then

           read(val,*) saturate

        else if( key == 'SKYMAG' ) then

           read(val,*) skymag

        else if( key == 'SKYGRADX' ) then

           read(val,*) skygrad(1)

        else if( key == 'SKYGRADY' ) then

           read(val,*) skygrad(2)

        else if( key == 'GNOISE' ) then

           read(val,*) gnoise

        else if( key == 'MAGLIM' ) then

           read(val,*) maglim

        else if( key == 'HWHM' ) then

           read(val,*) hwhm

        else if( key == 'ECCENTRICITY' ) then

           read(val,*) eccentricity
           if( .not. abs(eccentricity) < 1 ) &
                stop 'Error: eccentricity < 1 should be.'

        else if( key == 'INCLINATION' ) then

           read(val,*) inclination

        else if( key == 'AIRY' ) then

           read(val,*) airy
           airy_defined = .true.

        else if( key == 'BETA' ) then

           read(val,*) beta

        else if( key == 'FILTER' ) then

           read(val,*) filter

        else if( key == 'EXPTIME' ) then

           read(val,*) exptime

        else if( key == 'AREA' ) then

           read(val,*) area
           if( .not. (area > 0) ) stop 'Error: a non-positive area.'
           area_defined = .true.

        else if( key == 'QEFF' ) then

           read(val,*) qefficiency

        else if( key == 'EXTINK' ) then

           read(val,*) extink

        else if( key == 'COSMICS_FREQ' ) then

           read(val,*) cosmics_freq

        else if( key == 'COSMICS_CMAX' ) then

           read(val,*) cosmics_cmax

        else if( key == 'COSMICS_LEN' ) then

           read(val,*) cosmics_len

        else if( key == 'FOV' ) then

           read(val,*) fov
           fov_defined = .true.

        else if( key == 'DIAMETER' ) then

           read(val,*) diameter
           if( .not. (diameter > 0) ) stop 'Error: a non-positive diameter.'
           diameter_defined = .true.

        elseif( key == 'PROJECTION' ) then

           read(val,*) prjcode

        else if( key == 'ANGLE' ) then

           read(val,*) angle

        else if( key == 'SCALE' ) then

           read(val,*) scale

        else if( key == 'CRVAL' ) then

           read(val,*) crval

        else if( key == 'CAT' ) then

           read(val,*) cat

        else if( key == 'GALCAT' ) then

           read(val,*) galcat

        else if( key == 'COL_RA' ) then

           read(val,*) labels(1)

        else if( key == 'COL_DEC' ) then

           read(val,*) labels(3)

        else if( key == 'COL_PMRA' ) then

           read(val,*) labels(2)

        else if( key == 'COL_PMDEC' ) then

           read(val,*) labels(4)

        else if( key == 'COL_MAG' ) then

           read(val,*) labels(5)

        else if( key == 'PHSYSTABLE' ) then

           read(val,*) phsystable

        else if( key == 'PHOTSYS' ) then

           read(val,*) phsystem

        else if( key == 'GEOGRAPHIC' ) then

           read(val,*) longitude, latitude

        else if( key == 'TIMESTEP' ) then

           read(val,*) timestep

        else if( key == 'DATE' ) then

           read(val,*) year,month,day
           datetime_init = .true.

        else if( key == 'TIME' ) then

           read(val,*) hour,minute,sec
           datetime_init = .true.

        else if( key == 'ATMOSPHERE' ) then

           read(val,*) atmosphere

        else if( key == 'LCMODEL' ) then

           read(val,*) lcmodel

        else if( key == 'LCTABLE' ) then

           read(val,*) lctable

        else if( key == 'LCFOURIER' ) then

           read(val,*) lcfourier

        else if( key == 'LCMAG' ) then

           read(val,*) lcmag

        else if( key == 'LCAMP' ) then

           read(val,*) lcamp

        else if( key == 'LCJD0' ) then

           read(val,*) lcjd0

        else if( key == 'LCPER' ) then

           read(val,*) lcper
           if( .not. (lcper > 0) ) stop 'Error: a non-positive period.'

        else if( key == 'LCCOO' ) then

           read(val,*) lcoo
           lccoo_init = .true.

        else if( key == 'NOISE' ) then

           read(val,*) noise

        else if( key == 'NFILES' ) then

           read(val,*) nfiles

        else if( key == 'OUTPUT' ) then

           read(val,*) filemask

        end if

     end do

     crpix = naxes / 2
     if( fov_defined ) then
        scale = naxes(2) / fov
     else
        fov = naxes(2) / scale
     end if
     ! The diameter defines the radius (and both the area and Airy disk
     ! unless they are given explicitly); the area defines the radius
     ! when the diameter is not given.
     if( area_defined .and. .not. diameter_defined ) then
        radii = sqrt(area / pi)
     else
        radii = diameter / 2
     end if
     if( .not. area_defined ) area = pi*radii**2
     if( .not. airy_defined ) airy = 7e-2 / radii

     if( bitpix > 0 ) then
        if( bitpix > digits(bitpix) ) then
           saturate = min(saturate,real(huge(bitpix)-1,kind(saturate)))
        else
           saturate = min(saturate,2.0_REAL64**bitpix-1)
        end if
     end if
     if( .not. (saturate > 0) ) stop 'saturate > 0 unsatisfied.'

     if( spread == 'AUTO' ) then
        if( exptime < 3 ) then
           spread = 'RANDOM'
        else
           spread = 'FFT'
        end if
     end if

     if( .not. lccoo_init ) lcoo = crval

     ! photometry system, filter
     call phselect(phsystable,phsystem,phsys)

     block
       use noise
       cosmics_nevents = int(pnoise(cosmics_freq*exptime))
     end block

     ! stars
     if( cat /= '' ) then
        call catstars(cat,labels,prjcode,scale,angle,crval,crpix,phsys, &
             filter,xcoo,ycoo,flux)
     else
        call genstars(maglim,naxes,xcoo,ycoo,flux)
     end if
     nstars_coo = size(xcoo)

     if( psffile /= '' ) then
        call LoadSprayFile(psffile,mpsf,rcoo,acoo,rgrid,agrid,zoom,nbox)
        mpsf_shape = shape(mpsf)
     end if

     if( galcat /= '' ) then
        call galaxies%Load(galcat)
        ngal = size(galaxies%gs)
     end if

     lcmode = lcmodel /= '' .or. lctable /= '' .or. lcfourier /= ''
     if( lctable /= '' ) then
        call LoadLCTable(lctable,lctime,lclux,lcdux)
        nlc_tab = size(lctime)
     end if

     if( lcfourier /= ''    )then
        call LoadFourier(lcfourier,lcfcoeff)
        nlc_four = size(lcfcoeff)
     end if

     ! the flat-field model
     if( flatmodelfile /= '' ) then

        allocate(flatmodel(naxes(1),naxes(2)))
        call flatmodel_load(flatmodelfile,flatmodel,istat)

        if( istat == 0 ) then
           flatmodel_naxes = shape(flatmodel)
        else
           flatmodelfile = ''
        end if

     end if

     ! additional fits keys
     if( fitskeyfile /= '' ) then
        call include_keywords(fitskeyfile,keys_include,cards_include)
        mcount_keys = size(keys_include)
     end if

     if( verbose ) then
        write(error_unit,*)
        write(error_unit,*) "LIST OF PARAMETERS"
        write(error_unit,*)
        if( psffile == '' ) then
           write(error_unit,*) "Point spread function (PSF): ",trim(profile)
           write(error_unit,*) "Seeing spread method: ",trim(spread)
           write(error_unit,'(a,f0.2)') " HWHM = ",hwhm
           write(error_unit,'(a,f0.3)') " Eccentricity = ",eccentricity
           write(error_unit,'(a,f0.1)') " Inclination = ",inclination
           if( profile == 'MOFFAT' ) then
              write(error_unit,'(a,f0.2)') " Moffat exponent (beta) = ",beta
           else if( profile == 'SEEING' ) then
              write(error_unit,'(a,f0.2)') " Airy disk radius (airy) = ",airy
           end if
        else
           write(error_unit,*) 'PSF by file: ',trim(psffile)
        end if
        write(error_unit,'(a,f0.2)') " Quantum efficiency = ",qefficiency
        write(error_unit,*) "Photometry system: ",trim(phsystem)
        write(error_unit,*) "Filter: ",trim(filter)
        write(error_unit,*) "Atmosphere modelling: ",atmosphere
        if( atmosphere ) then
           write(error_unit,'(a,f0.2)') " Extinction coefficient = ",extink
        end if
        write(error_unit,'(a,f0.1,a)') " Sky magnitude = ",skymag," mag/arcsec2"
        write(error_unit,'(a,2(g0.2,2x),a)') &
             " Sky gradient = ",skygrad," mag/arcsec2/pix"
        write(error_unit,'(a,f0.2,a)') " Gaussian noise added = ",gnoise," [cts]"
        write(error_unit,'(a,f0.3,a)') " Area = ",area," [m2]"
        write(error_unit,'(a,f0.3,a)') " Exposure time = ",exptime," [s]"
        write(error_unit,'(a,f0.3,a)') " Gain = ",gain," [e-/adu]"
        write(error_unit,'(a,es0.3,a)') " Saturation = ",saturate," [cts]"
        write(error_unit,'(a,f0.3,a)') " Time span = ",timestep," [s]"
        write(error_unit,'(a,i4.4,2("-",i2.2))') " Exposure start date: ", &
             year,month,day
        write(error_unit,'(a,2(i2.2,":"),f0.3)') " Exposure start time: ", &
             hour,minute,sec
        write(error_unit,'(a,i0," x ",i0)') " Width and height of image(s): ",naxes
        write(error_unit,'(a,2f8.3,a)') " Equatorial coordinates: ",crval," [deg]"
        write(error_unit,'(a,f0.1,a,g0.3,a)') ' Scale: ',scale," [pix/deg],  ", &
             3600/scale," [arcsec/pix]"
        write(error_unit,'(a,f0.1,a)') ' Field rotation = ',angle," [deg]"
        write(error_unit,'(a,f0.2,a)') ' Field of view = ',fov," [deg]"
        write(error_unit,*) 'Celestial projection = ',prjcode
        if( cat == '' ) then
           write(error_unit,*) "Star coordinates and fluxes: random generated"
           write(error_unit,'(a,i0)') " Magnitude limit = ",maglim
        else
           write(error_unit,*) "Star coordinates and fluxes: by `",trim(cat),"'"
        end if
        write(error_unit,'(a,2f8.3,a)') " Observatory coordinates: ", &
             longitude, latitude," [deg]"
        write(error_unit,'(a,i0)') " Cosmic ray events (max): #",cosmics_nevents
        write(error_unit,*) "Galaxy catalogue: `",trim(galcat),"'"
        write(error_unit,*) "Flat model file: `",trim(flatmodelfile),"'"
        write(error_unit,*) "Additional FITS keywords by: `",trim(fitskeyfile),"'"
        write(error_unit,*) "Light curve modelling: ",lcmode
        if( lcmode ) then
           if( lcmodel /= '' ) then
              write(error_unit,*) "Light curve model: ",trim(lcmodel)
           else if( lctable /= '' ) then
              write(error_unit,*) "Light curve in table: `",trim(lctable),"'"
           else if( lcfourier /= '' ) then
              write(error_unit,*) "Light curve by Fourier coefficients in: `", &
                   trim(lcfourier),"':"
           end if
           write(error_unit,'(a,f0.3)') " Light curve mean magnitude = ",lcmag
           write(error_unit,'(a,f0.3)') " Light curve amplitude = ",lcamp
           write(error_unit,'(a,f0.6)') " Light curve JD0 = ",lcjd0
           write(error_unit,'(a,f0.6,a)') " Light curve period = ",lcper, &
                " [day(s)]"
           write(error_unit,'(a,2f7.3,a)') " Coordinates of the object: ",lcoo, &
                " [deg]"
        end if
        write(error_unit,*)
     end if

     if( datetime_init ) jd0 = datetimejd(year,month,day,hour,minute,sec)
     timestep = timestep / 86400.0_REAL64

     do i = 1, num_images()
        nfiles[i] = nfiles
        filemask[i] = filemask
        timestep[i] = timestep
        jd0[i] = jd0
        prjcode[i] = prjcode
        filter[i] = filter
        psffile[i] = psffile
        profile[i] = profile
        spread[i] = spread
        galcat[i] = galcat
        flatmodelfile[i] = flatmodelfile
     end do

  end if ! this_image() == 1

  call co_broadcast(verbose,1)
  call co_broadcast(naxes,1)
  call co_broadcast(crpix,1)
  call co_broadcast(bitpix,1)
  call co_broadcast(gain,1)
  call co_broadcast(saturate,1)
  call co_broadcast(skymag,1)
  call co_broadcast(skygrad,1)
  call co_broadcast(maglim,1)
  call co_broadcast(exptime,1)
  call co_broadcast(area,1)
  call co_broadcast(hwhm,1)
  call co_broadcast(eccentricity,1)
  call co_broadcast(inclination,1)
  call co_broadcast(airy,1)
  call co_broadcast(beta,1)
  call co_broadcast(scale,1)
  call co_broadcast(crval,1)
  call co_broadcast(angle,1)
  call co_broadcast(qefficiency,1)
  call co_broadcast(gnoise,1)
  call co_broadcast(noise,1)
  call co_broadcast(extink,1)
  call co_broadcast(longitude,1)
  call co_broadcast(latitude,1)
  call co_broadcast(atmosphere,1)
  call co_broadcast(cosmics_cmax,1)
  call co_broadcast(cosmics_len,1)
  call co_broadcast(cosmics_nevents,1)
  call co_broadcast(phsys,1)
  call co_broadcast(lcmode,1)
  call co_broadcast(lcmag,1)
  call co_broadcast(lcamp,1)
  call co_broadcast(lcjd0,1)
  call co_broadcast(lcper,1)
  call co_broadcast(lcoo,1)
  call co_broadcast(lcmodel,1)
  call co_broadcast(lctable,1)
  call co_broadcast(lcfourier,1)
  call co_broadcast(nlc_four,1)
  call co_broadcast(nlc_tab,1)
  call co_broadcast(nbox,1)
  call co_broadcast(zoom,1)
  call co_broadcast(rgrid,1)
  call co_broadcast(agrid,1)
  call co_broadcast(mpsf_shape,1)
  call co_broadcast(ngal,1)
  call co_broadcast(nstars_coo,1)
  call co_broadcast(mcount_keys,1)
  call co_broadcast(flatmodel_naxes,1)

  sync all

  ! allocate shared variables
  if( this_image() > 1 ) then
     n = nstars_coo
     allocate(xcoo(n),ycoo(n),flux(n))

     if( all(mpsf_shape > 0) ) then
        m = mpsf_shape(1)
        n = nbox
        allocate(mpsf(m,-n:n,-n:n),rcoo(m),acoo(m))
     end if
     if( ngal > 0 ) allocate(galaxies%gs(ngal))
     if( nlc_tab > 0 ) &
          allocate(lctime(nlc_tab),lclux(nlc_tab),lcdux(nlc_tab))
     if( nlc_four > 0 ) allocate(lcfcoeff(nlc_four))
     if( all(flatmodel_naxes > 0) ) &
          allocate(flatmodel(flatmodel_naxes(1),flatmodel_naxes(2)))
     if( mcount_keys > 0 ) &
          allocate(keys_include(mcount_keys),cards_include(mcount_keys))
  end if

  sync all

  ! setup the variables
  call co_broadcast(xcoo,1)
  call co_broadcast(ycoo,1)
  call co_broadcast(flux,1)

  if( all(flatmodel_naxes > 0) ) call co_broadcast(flatmodel,1)

  if( mcount_keys > 0 ) then
     call co_broadcast(keys_include,1)
     call co_broadcast(cards_include,1)
  end if

  if( all(mpsf_shape > 0) ) then
     call co_broadcast(mpsf,1)
     call co_broadcast(rcoo,1)
     call co_broadcast(acoo,1)
  end if
  if( ngal > 0 ) then
     do i = 1, ngal
        call co_broadcast(galaxies%gs(i),1)
     end do
  end if
  if( nlc_tab > 0 ) then
     call co_broadcast(lctime,1)
     call co_broadcast(lclux,1)
     call co_broadcast(lcdux,1)
  end if
  if( nlc_four > 0 ) call co_broadcast(lcfcoeff,1)

  sync all
  if( verbose .and. (this_image() == 1 .and. num_images() > 1) ) then
     write(error_unit,'(a,i0,a)') ' Launching of the artificial sky in ',&
          num_images(),' parallel processes ...'
     write(error_unit,*)
     write(error_unit,*) "Filename, JD(begin), airmass, seeing radius:"
  end if
  ! setup finished and synced

  ! lauching processes
  block
    character(len=FLEN_FILENAME) :: file, fmt, mask
    integer :: n,l1,l2,mp,m

    mp = this_image()
    mask = filemask

    if( mask == '' ) then
       if( nfiles == 1 ) then
          mask = artificial_fits(:10) // artificial_fits(15:)
       else
          mask = artificial_fits
       end if
    end if

    ! determine '?' in filemask
    l1 = index(mask,'?')
    if( l1 > 0 ) then
       l1 = l1 - 1
       l2 = index(mask,'?',back=.true.) + 1
       n = l2 - l1 - 1
       write(fmt,'(a,i0,a)') '(a,i0.',n,',a)'
    end if
    do i = mp, nfiles, num_images()
       if( l1 > 0 ) then
          write(file,fmt) mask(:l1),i,mask(l2:)
       else
          file = mask
       end if
       call artfits(file,jd0 + (i-1)*timestep)
       mcaf[mp] = mcaf[mp] + 1
    end do
    if( verbose ) write(error_unit,'(2(a,i0))') 'finished process #',this_image(),&
         ', total images generated: ',mcaf[mp]

    ! the last finished process send an event to all images
    critical
      m = 0
      do i = 1, num_images()
         m = m + mcaf[i]
      end do
    end critical
    if( m == nfiles ) then
       do i = 1, num_images()
          event post(events[i])
       end do
    end if
  end block

  event wait(events)

  !  stop 0
  write(error_unit,*) 'STOP 0'

contains

  subroutine artfits(filename,jd)

    use marker
    use astrosphere
    use lcurves
    use cosmic_ray_events

    real(REAL64), parameter :: rad = 57.295779513082322865_REAL64

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: jd

    integer :: status
    character(len=FLEN_COMMENT) :: com
    character(len=FLEN_VALUE) :: datetime
    character(len=3) :: pcode
    real(REAL64), dimension(:,:), allocatable :: sky
    real(REAL64) :: c,s,x,qe,see,extin,hwhms,jd2, sec, back, skyint, ccmag
    integer :: year, month, day, hour, minute, m, i
    real(REAL64), dimension(2) :: bgrad
    real(REAL64), dimension(:), allocatable :: xstar,ystar,fluxes,xvar, yvar, fvar
    type(fitsfiles) :: fits
    type(SprayType), target :: spsf
    type(FileSprayType), target :: fpsf
    type(LightCurve), target :: lcwave
    type(LightCurveDeltaCep), target :: lcdcep
    type(LightCurveAlgol), target :: lcalgol
    type(LightCurveTabular), target :: lctab
    type(LightCurveFourier), target :: lcfour
    class(LightCurve), pointer :: lc => null()
    class(SprayType), pointer :: psf => null()

    ! background at zenit
    skyint = magiconv(phsys,filter, [skymag])
    back = sqrt(skyint**2 + gnoise**2)
    back = back * exptime * area * qefficiency * (3600/scale)**2
    bgrad = 1.086 * skygrad * back

    if( atmosphere ) then

       ! compute airmass
       x = xairmass(jd,longitude,latitude,crval(1),crval(2))
       if( x < 0 ) then
          call jdatetime(jd,year,month,day,hour,minute,sec)
          write(error_unit,'(a,i4.4,2(a,i2.2),a,2(i2.2,a),f6.3,a)') &
               "Warning: Airmass undefined (below horizon) for date  ", &
               year,"-",month,"-",day,"  and time ",hour,":",minute,":",sec," UT."
          return
       end if

       ! simulate seeing and extinction
       extin = exp(-extink*(x-1))          ! relative extinction against X=1
       qe = extin * qefficiency

       ! http://www.astro.auth.gr/~seeing-gr/seeing_gr_files/theory/node17.html
       ! Fries parameter for spread
       see = x**0.6
       hwhms = hwhm*see

       ! sky increases by airmass
       back  = back * extin
       bgrad = 1.086 * skygrad * back

       if( verbose ) write(error_unit,'(a,f15.5,f8.3,2f6.2)') &
            trim(filename),jd,x,extin,see
    else
       hwhms = hwhm
       qe = qefficiency

       if( verbose ) &
            write(error_unit,'(a,f15.5)') trim(filename),jd
    endif

    ! simulate profile
    if( psffile == '' ) then
       call spsf%SprayMe(profile,spread,hwhms,eccentricity,inclination,airy,beta)
       psf => spsf
    else
       call fpsf%SprayCan(mpsf,rcoo,acoo,rgrid,agrid,scale,zoom,nbox)
       psf => fpsf
    end if

    ! whatever LC is generated
    if( lcmode ) then
       ! variable stars

       if( lcmodel == 'WAVE' ) then
          lcwave = LightCurve(lcmag, lcamp, lcper, lcjd0, lcoo)
          lc => lcwave
       else if( lcmodel == 'CEPHEID' ) then
          lcdcep = LightCurveDeltaCep(lcmag, lcamp, lcper, lcjd0, lcoo)
          lc => lcdcep
       else if( lcmodel == 'ALGOL' ) then
          lcalgol = LightCurveAlgol(lcmag, lcamp, lcper, lcjd0, lcoo)
          lc => lcalgol
       else if( lcmodel /= '' ) then
          stop "Error: a light curve model unrecognised."
       end if

       if( lctable /= '' ) then
          lctab = LightCurveTabular(lcmag, lcamp, lcper, lcjd0, lcoo)
          call lctab%SetLCTable(lctime,lclux,lcdux)
          lc => lctab
        end if

        if( lcfourier /= ''    )then
           lcfour = LightCurveFourier(lcmag, lcamp, lcper, lcjd0, lcoo)
           call lcfour%SetFourCoefficients(lcfcoeff)
           lc => lcfour
        end if

        ccmag = lc%magvar(jd)
        allocate(xvar(1),yvar(1),fvar(1))
        call rain(prjcode,[lc%coo(1)],[lc%coo(2)],[ccmag],scale,angle,crval,crpix, &
             phsys,filter,xvar,yvar,fvar)
    else
       allocate(xvar(0),yvar(0),fvar(0))
    end if

    ! Frame generator
    allocate(sky(naxes(1),naxes(2)))
    sky = 0

    call background(back,bgrad,noise,sky)

    if( ngal > 0 ) then
       call galaxies%Render(prjcode,crval,crpix,scale,angle,exptime,area,qe, &
            phsys,filter,sky)
    end if

    ! merge both lists: variables and ordinary stars
    m = size(xvar)
    n = m + size(xcoo)
    allocate(xstar(n),ystar(n),fluxes(n))
    xstar(1:m) = xvar;  xstar(m+1:) = xcoo
    ystar(1:m) = yvar;  ystar(m+1:) = ycoo
    fluxes(1:m) = fvar; fluxes(m+1:) = flux

    ! update fluxes
    fluxes = fluxes * exptime * area * qe

    ! create the artificial sky
    call spreadstars(psf,xstar,ystar,fluxes,crpix,noise,sky)

    ! flat model
    if( all(flatmodel_naxes > 0) ) sky = sky * flatmodel

    ! cosmics
    call cosmics(cosmics_nevents,cosmics_cmax,cosmics_len,sky)

    ! make it realistic
    sky = sky / gain
    sky = min(sky,saturate)

    ! watermark
    call watermark(sky,5*sqrt(back))

    ! save the generated artificial sky to FITS file
    status = 0
    if( fits_file_exist(filename) ) call fits_file_delete(filename)
    call fits_create_file(fits,filename,status)
    if( status /= 0 ) then
       write(error_unit,*) "Error: `",trim(filename),"' failed to create."
       goto 666
    end if

    ! convert JD to date
    call jdatetime(jd,year,month,day,hour,minute,sec)
    jd2 = jd + exptime / 2.0_REAL64 / 86400.0_REAL64

    ! the artifical frame
    call fits_insert_image(fits,bitpix,size(naxes),naxes,status)
    write(datetime,&
         "(i4.4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',f6.3,tl6,i2.2)") &
         year, month, day, hour, minute, sec, nint(sec)
    ! There is the trick as Brno: f6.3 edition does not leave leading zeros
    ! unlike to i2.2. I found on slack-overflow the advice --- go back
    ! and re-write the integer part. This is .., The Dude said.
    call fits_write_key(fits,FITS_KEY_DATEOBS,datetime,'date and time in UT',&
         status)
    call fits_write_key(fits,FITS_KEY_FILTER,filter,'filter',status)
    call fits_write_key(fits,FITS_KEY_GAIN,gain,-3,'gain',status)
    call fits_write_key(fits,FITS_KEY_SATURATE,saturate,-3,'saturate value',status)
    call fits_write_key(fits,FITS_KEY_JD,jd2,-16, &
         '[day] Julian date at half of exposure',status)
    call fits_write_key(fits,FITS_KEY_EXPTIME,exptime,-3,'[s] exposure time', &
         status)
    call fits_write_key(fits,FITS_KEY_AREA,area,-2,'[m2] detection area',status)
    call fits_write_key(fits,FITS_KEY_HWHM,hwhms,-3, &
         '[pix] half with at half maximum',status)
    call fits_write_key(fits,FITS_KEY_ECCENTRICITY,eccentricity,-3, &
         'eccentricity of ellipsis (0=circle)',status)
    call fits_write_key(fits,FITS_KEY_INCLINATION,inclination,-1, &
         '[deg] inclination of major semiaxis',status)
    if( atmosphere ) &
         call fits_write_key(fits,FITS_KEY_AIRMASS,x,-4,'airmass',status)
    call fits_write_key(fits,'PSF',profile,'PSF',status)
    if( profile == 'MOFFAT' ) then
       call fits_write_key(fits,'MOFFATB',beta,-3,'Moffat exponent',status)
    else if ( profile == 'SEEING' ) then
       call fits_write_key(fits,'AIRYDISK',airy,-3,'[pix] Airy disk radius', &
            status)
       call fits_write_key(fits,'SPREAD',spread,'seeing spread method ',status)
    end if
    pcode = prjcode
    call fits_write_key(fits,'CTYPE1','RA---'//pcode, &
         'the coordinate type for the first axis',status)
    call fits_write_key(fits,'CTYPE2','DEC--'//pcode, &
         'the coordinate type for the second axis',status)
    call fits_write_key(fits,'CRPIX1',real(crpix(1)),-6,'[pix]',status)
    call fits_write_key(fits,'CRPIX2',real(crpix(2)),-6,'[pix]',status)
    call fits_write_key(fits,'CRVAL1',crval(1),-17,'[deg]',status)
    call fits_write_key(fits,'CRVAL2',crval(2),-17,'[deg]',status)
    s = sin(angle/rad)
    c = cos(angle/rad)
    call fits_write_key(fits,'CD1_1',-c/scale,13,'',status)
    call fits_write_key(fits,'CD1_2',s/scale,13,'',status)
    call fits_write_key(fits,'CD2_1',s/scale,13,'',status)
    call fits_write_key(fits,'CD2_2',c/scale,13,'',status)
    call fits_write_key(fits,FITS_KEY_LONGITUDE,longitude,-5, &
         '[deg] geographic longitude (-east)',status)
    call fits_write_key(fits,FITS_KEY_LATITUDE,latitude,-5, &
         '[deg] geographic latitude (+north)',status)

    ! include additional records
    do i = 1, mcount_keys
       call fits_update_card(fits,keys_include(i),cards_include(i),status)
    end do

    call fits_write_comment(fits,'Artificial sky by Munipack.',status)
    call fits_write_comment(fits,'Parameters:',status)
    write(com,'(a,es10.3,1x,2(g0.2,1x))') &
         'Background: level, gradient = ',back,bgrad
    call fits_write_comment(fits,com,status)
    write(com,'(a,l2)') 'Noise: ',noise
    call fits_write_comment(fits,com,status)
    if( psffile == '' ) then
       write(com,'(2a)') 'PSF = ',trim(profile)
       call fits_write_comment(fits,com,status)
       if( profile == 'SEEING' ) then
          write(com,'(2a)') 'Seeing spread method = ',trim(spread)
          call fits_write_comment(fits,com,status)
       end if
    else
       write(com,'(2a)') 'PSF file = ',trim(psffile)
       call fits_write_comment(fits,com,status)
    end if
    if( cat == '' ) then
       write(com,'(a,i0)') 'Magnitude limit  = ',maglim
       call fits_write_comment(fits,com,status)
    end if
    write(com,'(a,g0.3)') 'The efficiency  = ', qe
    call fits_write_comment(fits,com,status)
    if( atmosphere ) then
       write(com,'(a,g0.3)') 'Extinction coefficient  = ', extink
       call fits_write_comment(fits,com,status)
       write(com,'(a,g0.3)') 'Seeing factor  = ', see
       call fits_write_comment(fits,com,status)
    end if
    if( lcmode ) then
       if( lcmodel /= '' ) then
          write(com,'(2a)') 'Light curve model = ',trim(lcmodel)
          call fits_write_comment(fits,com,status)
       end if
       if( lcfourier /= '' ) then
          write(com,'(2a)') 'Fourier coefficients table = ', trim(lcfourier)
          call fits_write_comment(fits,com,status)
       end if
       if( lctable /= '' ) then
          write(com,'(2a)') 'Spline interpolation table = ', trim(lctable)
          call fits_write_comment(fits,com,status)
       end if
       write(com,'(a,f0.3)') "Light curve current magnitude = ",ccmag
       call fits_write_comment(fits,com,status)
       write(com,'(a,2(f0.3,1x))') &
            "Light curve mean magnitude, amplitude = ",lcmag,lcamp
       call fits_write_comment(fits,com,status)
       write(com,'(a,2(f0.6,1x),a)') &
            "Light curve elements (JD0, period) = ",lcjd0,lcper," [day(s)]"
       call fits_write_comment(fits,com,status)
       write(com,'(a,2(f0.5,1x),a)') &
            "Light curve object coordinates = ",lcoo," [deg]"
       call fits_write_comment(fits,com,status)
    end if
    if( galcat /= '' ) then
       write(com,'(2a)') 'Galaxy file = ',trim(galcat)
       call fits_write_comment(fits,com,status)
    end if
    if( flatmodelfile /= '' .and. allocated(flatmodel) ) then
       write(com,'(2a)') 'Flat-model file = ',trim(flatmodelfile)
       call fits_write_comment(fits,com,status)
    end if
    call fits_write_comment(fits,MUNIPACK_VERSION,status)
    call fits_write_image(fits,0,real(sky),status)

    ! PSF extension
    if( psffile == '' ) then
       call fits_insert_img(fits,-32,2,shape(psf%psf),status)
       call fits_update_key(fits,'EXTNAME','PSF','',status)
       call fits_update_key(fits,'ZOOM',psf%zoom,'zoom',status)
       call fits_write_image(fits,0,real(psf%psf),status)
    end if

    ! Star table extension
    block
      character(len=FLEN_VALUE), dimension(3) :: ttype, tform, tunit
      integer :: nrows, srows, i, l, frow

      nrows = size(xstar)
      ttype(1) = FITS_COL_X
      ttype(2) = FITS_COL_Y
      ttype(3) = FITS_COL_COUNT
      tform  = '1E'
      tunit = ''
      tunit(1:2) = 'pix'
      tunit(3) = 'electron'
      call fits_insert_btbl(fits,nrows,ttype,tform,tunit,'CATALOGUE',status)
      call fits_write_comment(fits,'COUNT is the mean number of electrons ' // &
           'of star (before flat, gain),',status)
      call fits_write_comment(fits,'the image is in ADU: ' // &
           'divide by gain (and apply the flat model).',status)

      call fits_get_rowsize(fits,srows,status)
      do i = 1, nrows, srows
         frow = i
         l = min(i+srows,nrows)

         call fits_write_col(fits,1,frow,xstar(i:l),status)
         call fits_write_col(fits,2,frow,ystar(i:l),status)
         call fits_write_col(fits,3,frow,fluxes(i:l),status)
      end do
    end block

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine artfits


  function magiconv(phsys,filter,mag)

    use photoconv

    real(REAL64) :: magiconv
    character(len=*), intent(in) :: filter
    real(REAL64), dimension(:), intent(in) :: mag
    type(type_phsys), intent(in) :: phsys
    real(REAL64), dimension(size(mag)) :: ph,dph,dmag

    dmag = 9.999
    call phsysmagph1(phsys,filter,mag,dmag,ph,dph)
    magiconv = ph(1)

  end function magiconv

end program artificial
