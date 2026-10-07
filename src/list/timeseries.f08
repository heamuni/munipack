!
!  time series of photometric quantities (light curves) extracted
!              from calibrated images
!
!  Copyright © 2012-4, 2018-20 F.Hroch (hroch@physics.muni.cz)
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


program timeseries

  use savelcfits
  use fitscat
  use fitsts
  use titsio
  use iso_fortran_env

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)

  character(len=FLEN_FILENAME) :: output = 'timeserie.fits'
  character(len=FLEN_FILENAME) :: cat = ''
  character(len=FLEN_VALUE) :: timetype = 'JD' ! time type: JD,MJD,HJD,PHASE
  character(len=FLEN_VALUE) :: timestamp = 'MID' ! stamp: BEGIN, MID, END
  character(len=FLEN_VALUE) :: cootype = 'EQT'   ! coordinates EQT, RECT
  character(len=FLEN_FILENAME), dimension(:), allocatable :: files
  character(len=FLEN_KEYWORD), dimension(:), allocatable :: keys
  character(len=FLEN_VALUE), dimension(:), allocatable :: cols, coocols
  real(dbl), dimension(:), allocatable :: xcoo, ycoo
  character(len=FLEN_KEYWORD), dimension(5) :: fitskeys
  character(len=FLEN_VALUE) :: extname = PHOTOEXTNAME ! default data extension
  logical :: stdout = .false.    ! results to std.out.
  logical :: printname = .false. ! print filenames
  logical :: printtime = .true.  ! print time-like
  logical :: printhorizon = .false.  ! print horizonal coordinates
  logical :: printairmass = .false.  ! print air-mass
  logical :: verbose = .false.   ! be verbose
  real(dbl) :: tolerance = -1    ! cross-match tolerance
  real(dbl) :: lc_epoch = 0      ! reference time for phase
  real(dbl) :: lc_period = 1     ! period
  logical :: lc_set = .false.    ! light curve paramaters defined?
  logical :: status = .false.    ! output file state
  integer :: nvalid = 0          ! count of valid fitses
  logical :: cmatch = .true.     ! full star match for alls?
  logical :: keys_found = .true. ! full star match for alls?
  integer :: n,eq,ncoo,nfiles,ncols,ncoocols,nkeys,stat
  character(len=4*FLEN_FILENAME) :: record,key,val
  character(len=80) :: msg

  type(TSFits), dimension(:), allocatable :: fitses
  type(CatFits), allocatable :: cfits

  nfiles = 0
  ncoo = 0
  ncols = 0
  nkeys = 0
  allocate(cols(0),keys(0),xcoo(0),ycoo(0))

  ! default coordinates column labels
  ncoocols = 2
  allocate(coocols(2))
  coocols = [ FITS_COL_RA, FITS_COL_DEC ]

  ! mandatory keys
  fitskeys = [ FITS_KEY_DATEOBS, FITS_KEY_TIMEOBS, FITS_KEY_EXPTIME, &
       FITS_KEY_LONGITUDE, FITS_KEY_LATITUDE ]

  do
     read(*,'(a)',iostat=stat,iomsg=msg) record
     if( stat == IOSTAT_END ) exit
     if( stat > 0 ) then
        write(error_unit,*) trim(msg)
        error stop 'An input error.'
     end if

     eq = index(record,'=')
     if( eq == 0 ) error stop 'Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'OUTPUT' ) then

        read(val,*) output

     else if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'EXTNAME' ) then

        read(val,*) extname

     else if( key == 'NCOOCOL' ) then

        read(val,*) ncoocols
        deallocate(coocols)
        allocate(coocols(ncoocols))

     else if( key == 'COOCOL' ) then

        read(val,*) coocols

     else if( key == 'NCOL' ) then

        read(val,*) ncols
        deallocate(cols)
        allocate(cols(ncols))

     else if( key == 'COL' ) then

        read(val,*) cols

     else if( key == 'NKEY' ) then

        read(val,*) nkeys
        deallocate(keys)
        allocate(keys(nkeys))

     else if( key == 'KEY' ) then

        read(val,*) keys

     else if( key == 'COOTYPE' ) then

        read(val,*) cootype

     else if( key == 'TIMETYPE' ) then

        read(val,*) timetype

     else if( key == 'TIMESTAMP' ) then

        read(val,*) timestamp

     else if( key == 'NCOO' ) then

        read(val,*) ncoo

        deallocate(xcoo,ycoo)
        allocate(xcoo(ncoo),ycoo(ncoo))
        ncoo = 0

     else if( key == 'COO' ) then

        ncoo =  ncoo + 1
        read(val,*) xcoo(ncoo), ycoo(ncoo)

     else if( key == 'TOL' ) then

        read(val,*) tolerance

     else if( key == 'LC_EPOCH' ) then

        read(val,*) lc_epoch
        lc_set = .true.

     else if( key == 'LC_PERIOD' ) then

        read(val,*) lc_period
        lc_set = .true.

     else if( key == 'FITS_KEY_DATEOBS' ) then

        read(val,*) fitskeys(1)

     else if( key == 'FITS_KEY_TIMEOBS' ) then

        read(val,*) fitskeys(2)

     else if( key == 'FITS_KEY_EXPTIME'  ) then

        read(val,*) fitskeys(3)

     else if( key == 'FITS_KEY_LONGITUDE'  ) then

        read(val,*) fitskeys(4)

     else if( key == 'FITS_KEY_LATITUDE'  ) then

        read(val,*) fitskeys(5)

     else if( key == 'STDOUT' ) then

        read(val,*) stdout

     else if( key == 'PRINTNAME' ) then

        read(val,*) printname

     else if( key == 'PRINTTIME' ) then

        read(val,*) printtime

     else if( key == 'PRINTHORIZON' ) then

        read(val,*) printhorizon

     else if( key == 'PRINTAIRMASS' ) then

        read(val,*) printairmass

     else if( key == 'CATALOGUE' ) then

        read(val,*) cat

     else if( key == 'NFILES' ) then

        read(val,*) nfiles

        allocate(files(nfiles))
        nfiles = 0

     else if( key == 'FILE' ) then

        nfiles = nfiles + 1
        read(val,*) files(nfiles)

     end if

  end do

  if( nfiles == 0 ) stop 'No files on input.'

  if( .not. lc_set .and. timetype == 'PHASE' .and. verbose) &
       write(error_unit,*) &
       'Warning: Timetype PHASE invoked but light curve elements was not given.'

  ! read coordinates
  allocate(cfits)
  if( ncoo == 0 .and. cat /= '' ) then
     call cfits%Catalogue(cat)
     if( .not. cfits%status ) stop 'Failed to read a catalogue.'
  else if( ncoo > 0 .and. cat == '' )then
     call cfits%Coords(xcoo,ycoo)
  else !left it with defaults
     call cfits%ZeroInit()
  end if

  ! processsing
  allocate(fitses(nfiles))
  nvalid = 0
  do n = 1, nfiles
     call fitses(n)%Load(files(n),extname,fitskeys,keys,coocols,cootype, &
          cols,tolerance,cfits,verbose)
     if( fitses(n)%status ) then
        nvalid = nvalid + 1
        cmatch = fitses(n)%cmatch .and. cmatch
        keys_found = fitses(n)%key_found .and. keys_found
        if ( stdout ) then
           call fitses(n)%report(timetype,timestamp,lc_epoch,lc_period, &
                printname,printtime,printhorizon,printairmass)
        end if
     end if
  end do

  if( nvalid == 0 ) stop 'Error: There is no valid FITS file.'

  if( .not. cmatch .and. .not. verbose) write(error_unit,*) &
       "Warning: One or more stars not found. Use --verbose for details."

  if( .not. keys_found .and. .not. verbose ) write(error_unit,*) &
       "Warning: One or more keywords not found. Use --verbose for details."

  ! save data
  call savelc(output,fitses,timetype,timestamp,cfits,cols, &
       lc_epoch,lc_period,printname,printtime,printhorizon,printairmass,status)

  deallocate(fitses,cfits,files,keys,cols,coocols,xcoo,ycoo)

  if( status ) then
     stop 0
  else
     stop 'An error occurred during timeseries session.'
  end if

end program timeseries
