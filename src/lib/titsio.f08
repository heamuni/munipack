!
!  FITSIO definitions for Munipack
!
!  Copyright © 2020-2026 F.Hroch (hroch@physics.muni.cz)
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

module titsio

  use fitsio_mmviii
  use iso_fortran_env

  implicit none

  ! version identification like:
  ! character(len=*), parameter, private :: VERSION = "0"
  include 'fversion.inc'

  character(len=*), parameter :: MUNIPACK_VERSION = &
       "Munipack "//VERSION//", (C)1997-2026 F. Hroch <hroch@physics.muni.cz>"

  character(len=*), parameter :: FITS_VALUE_CREATOR = &
       "Munipack "//VERSION
  character(len=*), parameter :: FITS_COM_CREATOR = &
       "https://munipack.physics.muni.cz"

  character(len=*), parameter :: FINDEXTNAME = 'FIND'
  character(len=*), parameter :: APEREXTNAME = 'APERPHOT'
  character(len=*), parameter :: GROWEXTNAME = 'GROWPHOT'
  character(len=*), parameter :: GROWCURVEXTNAME = 'GROWCURVE'
  character(len=*), parameter :: GROWFUNCEXTNAME = 'GROWFUNC'
  character(len=*), parameter :: GROWDATEXTNAME = 'GROWDATA'
  character(len=*), parameter :: PSFEXTNAME  = 'PSFPHOT'
  character(len=*), parameter :: PHOTOEXTNAME= 'PHOTOMETRY'
  character(len=*), parameter :: MEXTNAMETS =  'TIMESERIES'
  character(len=*), parameter :: MEXTNAMETSC = 'CATALOGUE'
  character(len=*), parameter :: FHDUNAME =    'PHOTOSYS'
  character(len=*), parameter :: FTHDUNAME =   'FOTRAN'
  character(len=*), parameter :: EXT_STDERR =  'STDERR'
  character(len=*), parameter :: EXT_MAXDEV =   'MAXDEV'
  character(len=*), parameter :: EXT_PHRES =   'PHRES'
  character(len=*), parameter :: ASTROREPORT = 'ASTROREPORT'

  character(len=*), parameter :: BEGIN_ASTROMETRY = &
       '=== Astrometric Solution by Munipack ==='
  character(len=*), parameter :: END_ASTROMETRY = &
       '=== End of Astrometric Solution by Munipack ==='
  character(len=*), parameter :: BEGIN_PHOTOCAL = &
       '=== Photometric Calibration by Munipack ==='
  character(len=*), parameter :: END_PHOTOCAL = &
       '=== End of Photometric Calibration by Munipack ==='

  character(len=*), parameter :: FITS_KEY_NAPER = 'NAPER'
  character(len=*), parameter :: FITS_KEY_APER = 'APER'
  character(len=*), parameter :: FITS_KEY_SAPER = 'SAPER'
  character(len=*), parameter :: FITS_KEY_ANNULUS = 'ANNULUS'
  character(len=*), parameter :: FITS_KEY_FWHM = 'FWHM'
  character(len=*), parameter :: FITS_KEY_HWHM = 'HWHM'
  character(len=*), parameter :: FITS_KEY_ECCENTRICITY = 'ECCENTR'
  character(len=*), parameter :: FITS_KEY_INCLINATION = 'INCL'
  character(len=*), parameter :: FITS_KEY_RF90 = 'RADFLX90'
  character(len=*), parameter :: FITS_KEY_EXTINK0 = 'EXTIN_K0'
  character(len=*), parameter :: FITS_KEY_EXTINR = 'EXTIN_R'
  character(len=*), parameter :: FITS_KEY_EXTINREF = 'EXTINREF'
  character(len=*), parameter :: FITS_KEY_THRESHOLD = 'THRESH'
  character(len=*), parameter :: FITS_KEY_LOWBAD = 'LOWBAD'
  character(len=*), parameter :: FITS_KEY_HIGHBAD = 'HIGHBAD'
  character(len=*), parameter :: FITS_KEY_RNDLO = 'RNDLO'
  character(len=*), parameter :: FITS_KEY_RNDHI = 'RNDHI'
  character(len=*), parameter :: FITS_KEY_SHRPLO = 'SHRPLO'
  character(len=*), parameter :: FITS_KEY_SHRPHI = 'SHRPHI'
  character(len=*), parameter :: FITS_KEY_PHOTPLAM = 'PHOTPLAM'
  character(len=*), parameter :: FITS_KEY_PHOTZPT = 'PHOTZPT'
  character(len=*), parameter :: FITS_KEY_PHOTFLAM = 'PHOTFLAM'
  character(len=*), parameter :: FITS_KEY_PHOTBW = 'PHOTBW'
  character(len=*), parameter :: FITS_KEY_CTPH = 'CTPH'
  character(len=*), parameter :: FITS_KEY_CTPHERR = 'CTPHERR'
  character(len=*), parameter :: FITS_KEY_PHSYSERR = 'PHSYSERR'
  character(len=*), parameter :: FITS_KEY_PHSYSERE = 'PHSYSERE'
  character(len=*), parameter :: FITS_KEY_PHSYSNST = 'PHSYSNST'
  character(len=*), parameter :: FITS_KEY_APCORTYP = 'APCORTYP'
  character(len=*), parameter :: FITS_KEY_APCORR = 'APCORR'
  character(len=*), parameter :: FITS_KEY_CSPACE = 'CSPACE'
  character(len=*), parameter :: FITS_KEY_REFRAME = 'REFRAME'
  character(len=*), parameter :: FITS_KEY_SATURATE = 'SATURATE'
  character(len=*), parameter :: FITS_KEY_GAIN = 'GAIN'
  character(len=*), parameter :: FITS_KEY_READNS = 'READNS'
  character(len=*), parameter :: FITS_KEY_AREA = 'AREA'
  character(len=*), parameter :: FITS_KEY_EXPTIME = 'EXPTIME'
  character(len=*), parameter :: FITS_KEY_PHOTSYS = 'PHOTSYS'
  character(len=*), parameter :: FITS_KEY_FILTER = 'FILTER'
  character(len=*), parameter :: FITS_KEY_FILTREF = 'FILTREF'
  character(len=*), parameter :: FITS_KEY_OBJECT = 'OBJECT'
  character(len=*), parameter :: FITS_KEY_OBSERVER = 'OBSERVER'
  character(len=*), parameter :: FITS_KEY_ORIGIN = 'ORIGIN'
  character(len=*), parameter :: FITS_KEY_AUTHOR = 'AUTHOR'
  character(len=*), parameter :: FITS_KEY_INSTRUME = 'INSTRUME'
  character(len=*), parameter :: FITS_KEY_TELESCOP = 'TELESCOP'
  character(len=*), parameter :: FITS_KEY_BIBREF = 'BIBREF'
  character(len=*), parameter :: FITS_KEY_LONGITUDE = 'LONGITUD'
  character(len=*), parameter :: FITS_KEY_LATITUDE =  'LATITUDE'
  character(len=*), parameter :: FITS_KEY_EPOCH = 'JEPOCH'
  character(len=*), parameter :: FITS_KEY_RADESYS = 'RADESYS'
  character(len=*), parameter :: FITS_KEY_CREATOR = 'CREATOR'
  character(len=*), parameter :: FITS_KEY_ORIGHDU = 'ORIGHDU'
  character(len=*), parameter :: FITS_KEY_BUNIT = 'BUNIT'
  character(len=*), parameter :: FITS_KEY_SKYMAG = 'SKYMAG'
  character(len=*), parameter :: FITS_KEY_SKYMEAN = 'SKYMEAN'
  character(len=*), parameter :: FITS_KEY_SKYSIG = 'SKYSIG'
  character(len=*), parameter :: FITS_KEY_SKYSTD = 'SKYSTD'
  character(len=*), parameter :: FITS_KEY_IMAGETYP = 'IMAGETYP'
  character(len=*), parameter :: FITS_KEY_DATEOBS = 'DATE-OBS'
  character(len=*), parameter :: FITS_KEY_TIMEOBS = 'TIME-OBS'
  character(len=*), parameter :: FITS_KEY_TEMPERATURE = 'TEMPERAT'
  character(len=*), parameter :: FITS_KEY_AIRMASS = 'AIRMASS'
  character(len=*), parameter :: FITS_KEY_JD = 'JD'
  character(len=*), parameter :: FITS_KEY_HJD = 'HJD'
  character(len=*), parameter :: FITS_KEY_FILENAME = 'FILENAME'
  character(len=*), parameter :: FITS_KEY_PHOTOTYP = 'PHOTOTYP'

  ! definitions of column labels, common
  character(len=*), parameter :: FITS_COL_FILENAME = 'FILENAME'
  character(len=*), parameter :: FITS_COL_TIME = 'TIME'
  character(len=*), parameter :: FITS_COL_X = 'X'
  character(len=*), parameter :: FITS_COL_Y = 'Y'
  character(len=*), parameter :: FITS_COL_RA = 'RAJ2000'
  character(len=*), parameter :: FITS_COL_DEC = 'DEJ2000'
  character(len=*), parameter :: FITS_COL_PMRA = 'pmRA'
  character(len=*), parameter :: FITS_COL_PMDEC = 'pmDE'
  character(len=*), parameter :: FITS_COL_SKY = 'SKY'
  character(len=*), parameter :: FITS_COL_SKYERR = 'SKYERR'
  character(len=*), parameter :: FITS_COL_AZIMUTH = 'AZIMUTH'
  character(len=*), parameter :: FITS_COL_ZENITD = 'ZENITD'
  character(len=*), parameter :: FITS_COL_AIRMASS = 'AIRMASS'
  character(len=*), parameter :: FITS_COL_R = 'R'
  character(len=*), parameter :: FITS_COL_GROW = 'GROWCURVE'
  character(len=*), parameter :: FITS_COL_GROWERR = 'GROWCURVEERR'
  character(len=*), parameter :: FITS_COL_RPROF = 'RADIALPROFILE'
  character(len=*), parameter :: FITS_COL_RESGROW = 'RESGROW'
  character(len=*), parameter :: FITS_COL_GROWFLAG = 'GROWFLAG'

  ! definitions of column labels, find
  character(len=*), parameter :: FITS_COL_PEAKRATIO = 'PEAKRATIO'
  character(len=*), parameter :: FITS_COL_SHARP = 'SHARP'
  character(len=*), parameter :: FITS_COL_ROUND = 'ROUND'

  ! definitions of column labels, general photometry
  character(len=*), parameter :: FITS_COL_COUNT =     'COUNT'
  character(len=*), parameter :: FITS_COL_COUNTERR =  'COUNTERR'
  character(len=*), parameter :: FITS_COL_PHOTON =    'PHOTON'
  character(len=*), parameter :: FITS_COL_PHOTONERR = 'PHOTONERR'
  character(len=*), parameter :: FITS_COL_PHOTRATE =    'PHOTRATE'
  character(len=*), parameter :: FITS_COL_PHOTRATEERR = 'PHOTRATEERR'

  ! definitions of column labels, aperture photometry
  character(len=*), parameter :: FITS_COL_APCOUNT =    'APCOUNT'
  character(len=*), parameter :: FITS_COL_APCOUNTERR = 'APCOUNTERR'
  character(len=*), parameter :: FITS_COL_APEE =       'APEE'
  character(len=*), parameter :: FITS_COL_APEEERR =    'APEEERR'

  ! definitions of column labels, growth-curve photometry
  character(len=*), parameter :: FITS_COL_GCOUNT =    'GCOUNT'
  character(len=*), parameter :: FITS_COL_GCOUNTERR = 'GCOUNTERR'

  ! definitions of column labels, L-photometry
  character(len=*), parameter :: FITS_COL_LBCOUNT =    'LBCOUNT'
  character(len=*), parameter :: FITS_COL_LBCOUNTERR = 'LBCOUNTERR'

  ! definitions of column labels, PSF photometry
  character(len=*), parameter :: FITS_COL_PSFCOUNT =   'PSFCOUNT'
  character(len=*), parameter :: FITS_COL_PSFCOUNTERR ='PSFCOUNTERR'
  character(len=*), parameter :: FITS_COL_PSFPEAK =    'PSFPEAK'
  character(len=*), parameter :: FITS_COL_PSFPEAKERR = 'PSFPEAKERR'

  ! definitions of column labels, output photometry quantities
  character(len=*), parameter :: FITS_COL_PHOTNU =     'PHOTNU'
  character(len=*), parameter :: FITS_COL_PHOTNUERR =  'PHOTNUERR'
  character(len=*), parameter :: FITS_COL_PHOTLAM =    'PHOTLAM'
  character(len=*), parameter :: FITS_COL_PHOTLAMERR = 'PHOTLAMERR'
  character(len=*), parameter :: FITS_COL_FLUX =       'FLUX'
  character(len=*), parameter :: FITS_COL_FLUXERR =    'FLUXERR'
  character(len=*), parameter :: FITS_COL_FNU =        'FNU'
  character(len=*), parameter :: FITS_COL_FNUERR =     'FNUERR'
  character(len=*), parameter :: FITS_COL_FLAM =       'FLAM'
  character(len=*), parameter :: FITS_COL_FLAMERR=     'FLAMERR'
  character(len=*), parameter :: FITS_COL_MAG =        'MAG'
  character(len=*), parameter :: FITS_COL_MAGERR =     'MAGERR'
  character(len=*), parameter :: FITS_COL_ABMAG =      'ABMAG'
  character(len=*), parameter :: FITS_COL_ABMAGERR =   'ABMAGERR'
  character(len=*), parameter :: FITS_COL_STMAG =      'STMAG'
  character(len=*), parameter :: FITS_COL_STMAGERR =   'STMAGERR'
  character(len=*), parameter :: FITS_COL_RATE =       'RATE'
  character(len=*), parameter :: FITS_COL_RATEERR =    'RATEERR'

  character(len=*), parameter :: FITS_COL_FILTER =     'FILTER'

  character(len=*), parameter :: FITS_COL_LAMEFF =     'LAM_EFF'
  character(len=*), parameter :: FITS_COL_LAMFWHM =    'LAM_FWHM'
  character(len=*), parameter :: FITS_COL_NUEFF =      'NU_EFF'
  character(len=*), parameter :: FITS_COL_NUFWHM =     'NU_FWHM'
  character(len=*), parameter :: FITS_COL_FNUREF =     'FNU_REF'
  character(len=*), parameter :: FITS_COL_FLAMREF =    'FLAM_REF'


contains

  subroutine fits_get_dateobs(fitsfile,keys,dateobs,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), dimension(:), intent(in) :: keys
    character(len=*), intent(out) :: dateobs
    integer, intent(in out) :: status

    character(len=FLEN_VALUE) :: date, time
    character(len=FLEN_COMMENT) :: com
    integer :: status1,status2

    dateobs = ''
    if( status /= 0 ) return

    status1 = 0
    status2 = 0

    com = ''
    date = ''
    time = ''
    call fits_read_key(fitsfile,keys(1),date,com,status1)
    call fits_read_key(fitsfile,keys(2),time,com,status2)

    if( status1 == 0 .and. status2 == 0 ) then
       ! full date is in dateobs, full time is in timeobs
       ! this time specification is obsolete now (since 2000)
       dateobs = trim(date) // "T" // trim(time)
       status = 0
    else if ( status1 == 0 ) then
       ! both date and time are in one record, the standard way
       dateobs = date
       status = 0
    else
       status = status1
    end if

  end subroutine fits_get_dateobs

  function fits_jd(dateobs,status) result(jd)

    use trajd

    real(REAL64) :: jd
    character(len=*), intent(in) :: dateobs
    integer ::year, month, day, hour, minute
    real(REAL64) :: second, d
    integer :: status

    jd = 0
    if( .not. (status == 0 .and. dateobs /= '') ) return

    ! decode date-time string
    call fits_str2date(dateobs,year,month,day,hour,minute,second,status)
    if( status == 0 ) then

       ! Evaluate Julian date
       d = day + (hour + (minute + second/60.0_REAL64) / 60.0_REAL64) / 24.0_REAL64
       jd = datjd(real(year,REAL64),real(month,REAL64),d)
    end if

  end function fits_jd

  subroutine fits_read_wcs(fitsfile,radesys,epoch,ctype,crval,crpix,cd,crder,status)
    ! this procedure is used by colouring and lists
    ! it might be replaced by WCSLIB

    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(out) :: radesys
    real(REAL64), intent(out) :: epoch
    character(len=*), dimension(:), intent(out) :: ctype
    real(REAL64), dimension(:), intent(out) :: crval,crpix,crder
    real(REAL64), dimension(:,:), intent(out) :: cd
    integer, intent(in out) :: status

    character(len=FLEN_CARD) :: buf
    character(len=FLEN_KEYWORD) :: keyword, key
    integer :: n,m

    buf = ''
    keyword = ''
    key = ''
    ctype = ''
    radesys = ''

    ! read mandatory parameters
    call fits_read_key(fitsfile,FITS_KEY_RADESYS,radesys,buf,status)
    if( status /= 0 ) then
       status = 0
       radesys = ''
    end if

    call fits_read_key(fitsfile,FITS_KEY_EPOCH,epoch,buf,status)
    if( status /= 0 ) then
       status = 0
       epoch = IEEE_VALUE(epoch, IEEE_QUIET_NAN)
    end if

    ! read astrometric calibration
    do n = 1, size(ctype)
       call fits_make_keyn('CTYPE',n,keyword,status)
       call fits_read_key(fitsfile,keyword,ctype(n),buf,status)
    end do
    do n = 1, size(crval)
       call fits_make_keyn('CRVAL',n,keyword,status)
       call fits_read_key(fitsfile,keyword,crval(n),buf,status)
    end do
    do n = 1, size(crpix)
       call fits_make_keyn('CRPIX',n,keyword,status)
       call fits_read_key(fitsfile,keyword,crpix(n),buf,status)
    end do
    do n = 1, size(cd,1)
       call fits_make_keyn('CD',n,keyword,status)
       do m = 1, size(cd,2)
          call fits_make_keyn(trim(keyword)//"_",m,key,status)
          call fits_read_key(fitsfile,key,cd(n,m),buf,status)
       end do
    end do

    if( status /= 0 ) return

    ! optional keywords
    call fits_write_errmark
    do n = 1, size(crder)
       call fits_make_keyn('CRDER',n,keyword,status)
       call fits_read_key(fitsfile,keyword,crder(n),buf,status)
    end do
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       status = 0
       crder = 0
    end if

  end subroutine fits_read_wcs

  subroutine fits_update_wcs(fitsfile,radesys,epoch,ctype,crval,crpix,cd,crder,status)

    use, intrinsic :: ieee_arithmetic, only: IEEE_IS_NORMAL

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: radesys
    real(REAL64), intent(in) :: epoch
    character(len=*), dimension(:), intent(in) :: ctype
    real(REAL64), dimension(:), intent(in) :: crval,crpix
    real(REAL64), dimension(:,:), intent(in) :: cd
    real(REAL64), dimension(:), intent(in) :: crder
    integer, intent(in out) :: status

    integer, parameter :: digits = 15
    character(len=FLEN_KEYWORD) :: keyword, key
    integer :: n,m

    keyword = ''
    key = ''

    do n = 1, size(ctype)
       call fits_make_keyn('CTYPE',n,keyword,status)
       call fits_update_key(fitsfile,keyword,ctype(n), &
            'coordinate projection type',status)
    end do
    do n = 1, size(crval)
       call fits_make_keyn('CRVAL',n,keyword,status)
       call fits_update_key(fitsfile,keyword,crval(n),digits, &
            '[deg] spherical coordinates of center of projection',status)
    end do
    do n = 1, size(crpix)
       call fits_make_keyn('CRPIX',n,keyword,status)
       call fits_update_key(fitsfile,keyword,crpix(n),-6, &
            '[pix] reference pixel in focal plane',status)
    end do
    do n = 1, size(cd,1)
       call fits_make_keyn('CD',n,keyword,status)
       do m = 1, size(cd,2)
          call fits_make_keyn(trim(keyword)//"_",m,key,status)
          call fits_update_key(fitsfile,key,cd(n,m),digits, &
               '[deg/pix] scaled rotation matrix',status)
       end do
    end do

    ! optional
    if( all(crder > 0) ) then
       do n = 1, size(crder)
          call fits_make_keyn('CRDER',n,keyword,status)
          call fits_update_key(fitsfile,keyword,crder(n),1, &
               '[deg] standard errors',status)
       end do
    end if

    do n = 1, size(ctype)
       call fits_make_keyn('CUNIT',n,keyword,status)
       call fits_update_key(fitsfile,keyword,'deg','units of the CRVALx axis',status)
    end do

    ! mandatory keywords
    call fits_update_key(fitsfile,FITS_KEY_RADESYS,radesys, &
         'reference coordinate system',status)
    if( IEEE_IS_NORMAL(epoch) ) &
       call fits_update_key(fitsfile,FITS_KEY_EPOCH,epoch,-6, &
            'epoch of equatorial coordinates',status)

  end subroutine fits_update_wcs


end module titsio
