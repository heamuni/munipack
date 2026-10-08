!
!  grow-curve photometry
!
!  Copyright © 2016-8 F.Hroch (hroch@physics.muni.cz)
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


program gphot

  use titsio
  use iso_fortran_env

  implicit none

  character(len=4*FLEN_FILENAME) :: record,key,val
  character(len=FLEN_FILENAME) :: filename,output
  logical :: verbose = .false., plog = .false.
  integer :: eq
  real :: threshold = 0.1
  logical :: ex, exitus = .true.

  do
     read(*,'(a)',end=20) record

     eq = index(record,'=')
     if( eq == 0 ) stop 'Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'PIPELOG' ) then

        read(val,*) plog

     else if( key == 'THRESHOLD' ) then

        read(val,*) threshold

     else if( key == 'FILE' ) then

        read(val,*) filename, output

        call growfits(filename,output,ex)
        exitus = exitus .and. ex

     end if

  end do

20 continue

  if( exitus ) then
     stop 0
  else
     stop 'An error during growth-curve photometry occurred.'
  end if

contains

  subroutine growfits(filename,output,exitus)

    use grow_curve
    use grow_report
    use fitsgphot

    character(len=*),intent(in) :: filename,output
    logical, intent(out) :: exitus

    integer, parameter :: rp = selected_real_kind(15)

    integer :: status,nrows,naper
    real, dimension(:), allocatable :: xcens,ycens,sky,skyerr,skycorr,skyerrcorr
    real, dimension(:,:), allocatable :: apcts,apcts_err
    real(rp), dimension(:), allocatable :: gcount,gcount_err
    integer, dimension(:), allocatable :: growflag
    real(rp), dimension(:), allocatable :: curve, curve_err, prof
    real, dimension(:), allocatable :: raper
    real :: hwhm, rflux90, ghwhm, sep
    type(grow_reporter) :: reporter
    logical :: report_init

    report_init = .false.
    status = 0

    call fits_gphot_read(filename,hwhm,sep,raper,xcens,ycens, &
         sky,skyerr,skycorr,skyerrcorr,apcts,apcts_err,status)
    if( status /= 0 ) goto 666

    naper = size(raper)
    nrows = size(xcens)
    allocate(gcount(nrows),gcount_err(nrows), growflag(nrows), &
         curve(naper), curve_err(naper), prof(naper) )

    if( verbose ) then
       report_init = .true.
       call grow_report_init(reporter,naper,nrows,raper)
    end if

    ! angular separation of stars intended as base for growth-curve
    ! construction, the stars can share only the ring of sky
    !  sep = ring(2)

    ! growth-curve photometry
    call growphot(xcens,ycens,apcts,apcts_err,sky,skyerr,raper,hwhm,sep,&
         threshold, &
         gcount,gcount_err,skycorr,skyerrcorr,curve,curve_err,growflag, &
         prof,rflux90,ghwhm,verbose,reporter)

    ! update sky: the growth-curve fit models the aperture counts
    ! as f*curve + b*area, so b (skycorr) is the residual background
    ! per pixel left after subtraction of the sky; it is added back.
    ! The uncertainty of the sky estimate and of b are combined.
    sky = sky + skycorr
    where( skyerr > 0 .and. skyerrcorr > 0 )
       skyerr = sqrt(skyerr**2 + skyerrcorr**2)
!    elsewhere
!       skyerr = 0
    end where

    call fits_find_save(filename,output,ghwhm,rflux90, &
       raper, xcens, ycens, sky,skyerr, gcount,gcount_err, growflag, &
       curve, curve_err, prof, status)


!    if( verbose ) call grow_report_dump(15,reporter,status)

666 continue

    if( allocated(apcts) ) deallocate(raper,xcens,ycens,sky,skyerr,skycorr, &
         skyerrcorr,apcts,apcts_err)

    if( allocated(gcount) ) deallocate(gcount,gcount_err,growflag, &
         curve,curve_err,prof)

    if( verbose .and. report_init ) call grow_report_terminate(reporter)

    exitus = status == 0

  end subroutine growfits


end program gphot
