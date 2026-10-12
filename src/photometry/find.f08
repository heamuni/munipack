!
!  find stars
!
!  Copyright © 2013-4, 2016-25 F.Hroch (hroch@physics.muni.cz)
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


program find

  use titsio
  use iso_fortran_env

  implicit none

  logical, parameter :: debug = .false.
  integer, parameter :: maxsky = 4e6

  character(len=4*FLEN_FILENAME) :: record,key,val
  character(len=FLEN_FILENAME) :: filename, output
  character(len=FLEN_KEYWORD), dimension(2) :: fkeys
  character(len=80) :: msg
  real :: lothresh = 5.0
  real :: threshold = 7.0
  real :: fwhm = 3.0
  real :: readns_init = -1
  real :: saturation = -1
  real :: shrplo = 0.2, shrphi = 1.0
  real :: rndlo = -1.0, rndhi = 1.0
  ! the expected number of noise peaks per frame setting the threshold
  ! (<= 0: by the threshold), the shape tests by the noise [sigma]
  ! (<= 0: the fixed windows only), the windows given by the user
  real :: falserate = -1
  real :: shapesig = 3
  logical :: fixsharp = .false., fixround = .false.
  logical :: verbose = .false., plog = .false.
  logical :: ex, exitus = .true.
  integer :: eq, stat

  fkeys(1) = FITS_KEY_SATURATE
  fkeys(2) = FITS_KEY_READNS

  do
     read(*,'(a)',iostat=stat,iomsg=msg) record
     if( stat == IOSTAT_END ) exit
     if( stat > 0 ) then
        write(error_unit,*) trim(msg)
        error stop 'An input error.'
     end if

     eq = index(record,'=')
     if( eq == 0 ) stop 'Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'VERBOSE' ) then

        read(val,*) verbose
        if( verbose ) &
             write(error_unit,*) 'Filename, sky value [cts], No. of stars:'

     else if( key == 'PIPELOG' ) then

        read(val,*) plog

     else if( key == 'FWHM' ) then

        read(val,*) fwhm

     else if( key == 'READNOISE' ) then

        read(val,*) readns_init

     else if( key == 'SATURATE' ) then

        read(val,*) saturation

     else if( key == 'THRESHOLD' ) then

        read(val,*) threshold

     else if( key == 'FALSE_RATE' ) then

        read(val,*) falserate
        if( .not. (falserate > 0) ) stop 'Error in find: --false-rate > 0 expected.'

     else if( key == 'SHAPE_SIGMA' ) then

        read(val,*) shapesig

     else if( key == 'LOWER_THRESHOLD' ) then

        read(val,*) lothresh

     else if( key == 'ROUND_LOWER' ) then

        read(val,*) rndlo
        fixround = .true.

     else if( key == 'ROUND_HIGHER' ) then

        read(val,*) rndhi
        fixround = .true.

     else if( key == 'SHARP_LOWER' ) then

        read(val,*) shrplo
        fixsharp = .true.

     else if( key == 'SHARP_HIGHER' ) then

        read(val,*) shrphi
        fixsharp = .true.

     else if( key == 'FITS_KEY_SATURATE' ) then

           read(val,*) fkeys(1)

     else if( key == 'FITS_KEY_READNOISE' ) then

           read(val,*) fkeys(2)

     else if( key == 'FILE' ) then

        read(val,*) filename, output

        if( verbose ) &
             write(error_unit,'(a)',advance="no") trim(filename)//": "

        call the_finder(ex)
        exitus = exitus .and. ex

     end if

  end do

  if( exitus ) then
     stop 0
  else
     stop 'Some error(s) occurred during this star FIND.'
  end if

contains

  subroutine the_finder(exitus)

    use mdaosky
    use mdaofind
    use fitsfind
    use findnoise, only: findstats
    use skymean, only: sky_noise_model

    logical, intent(out) :: exitus

    real, dimension(:,:), allocatable :: data
    real :: lobad, hibad, hmin, readns,skymod, skyerr, skysig, satur
    real(REAL64) :: gain, rn2, bias
    type(findstats) :: stats
    integer :: nstep,nstar,status

    if( fwhm <= 0.0 ) stop 'Error in find: Assumption FWHM > 0 unsatisfied.'

    status = 0
    call fits_find_read(filename,fkeys,data,readns,satur,status)
    if( status /= 0 ) goto 666

    ! by default, saturation is assigned:
    !  * if presented in the file and not specified by user: by the file
    !  * if specified: by user
    if( saturation > 0 ) then
       hibad = saturation
    else
       hibad = satur
    end if
    if( readns_init > 0.0 ) readns = readns_init

    ! Sky estimate: the mean, clipped by the noise model of the frame
    call sky_noise_model(filename,gain,rn2,bias,status)
    if( status /= 0 ) goto 666
    if( readns_init > 0.0 ) rn2 = real(readns_init,REAL64)**2
    nstep = max(int(log(float(size(data))/maxsky) / log(2.0)),1)
    call daosky(data,nstep,verbose,hibad,skymod,skyerr,skysig,gain,rn2,bias)
    if( debug ) write(*,*) 'Sky:',skymod,'+-',skyerr,' sig=',skysig

    ! Detected stars are temporary stored in a scratch file under unit 3.
    ! The file is written by daofind and re-read by fits_find_save().
    open(3,status='scratch',form='unformatted')

    call daofind(data,debug,plog,hibad,skymod,skysig,fwhm,lothresh,threshold,&
         shrplo,shrphi,rndlo,rndhi,readns,falserate,shapesig,fixsharp,fixround, &
         nstar,lobad,hmin,stats,gain,bias)

    rewind(3)
    call fits_find_save(filename,output,fkeys,nstar,fwhm,real(stats%t),&
         shrplo,shrphi,rndlo,rndhi,readns,lothresh, lobad, hibad, hmin, &
         skymod, skyerr, skysig, maxsky, stats, falserate, status)

    if( verbose ) then
       write(error_unit,*) skymod, nstar
       write(error_unit,'(a,f0.3,a,f0.3,a,f0.4,a,g0.4,a)') &
            '  threshold ',stats%t,' (Normal deviate) = ',stats%x, &
            ' sigma of the filtered sky (skewness ',stats%skew, &
            '); expected noise peaks per frame ',stats%nfalse,'.'
       if( stats%sharpnoise .or. stats%roundnoise ) &
            write(error_unit,'(a,f0.1,a,f0.3,a,f0.3,a,i0,a)') &
            '  shape tests within ',stats%shapesig,' sigma of sharp ', &
            stats%sharpref,', round ',stats%roundref,' (',stats%nref, &
            ' bright stars; 0: the model).'
    end if
    close(3)

666 continue

    exitus = status == 0

    if( allocated(data) ) deallocate(data)

  end subroutine the_finder

end program find
