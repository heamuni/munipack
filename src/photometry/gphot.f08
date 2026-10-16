!
!  grow-curve photometry
!
!  Copyright © 2016-8, 2026 F.Hroch (hroch@physics.muni.cz)
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
!  The growth curve is the enclosed energy of the apertures by the
!  bright isolated stars of the frame (the machinery of aphot
!  --apcorr auto, apee.f08), varying over the field; the total flux of
!  every star is the estimate of the least variance by its apertures
!  and the growth curve (growee.f08). METHOD = 'LEGACY' runs the
!  growth-curve fit of the versions before (growcurve.f95).

program gphot

  use titsio
  use iso_fortran_env

  implicit none

  character(len=4*FLEN_FILENAME) :: record,key,val
  character(len=FLEN_FILENAME) :: filename,output
  character(len=FLEN_VALUE) :: method = 'EE'
  logical :: verbose = .false., plog = .false.
  integer :: eq
  real :: threshold = 0.1
  ! the order of the spatial variation of the growth curve, by the
  ! distance from the field centre (eecen when set)
  integer :: eevar = 2
  logical :: eeradial = .false., eecen_set = .false.
  real(REAL64), dimension(2) :: eecen = 0
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

     else if( key == 'METHOD' ) then

        read(val,*) method
        if( .not. (method == 'EE' .or. method == 'LEGACY') ) &
             stop 'Unknown METHOD (EE, LEGACY).'

     else if( key == 'EEVAR' ) then

        read(val,*) eevar
        if( eevar < 0 .or. eevar > 2 ) stop 'EEVAR is 0, 1 or 2.'

     else if( key == 'EERADIAL' ) then

        read(val,*) eeradial

     else if( key == 'EECENTRE' ) then

        read(val,*) eecen
        eecen_set = .true.

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
    use grow_ee
    use aperture_ee
    use fitsgphot

    character(len=*),intent(in) :: filename,output
    logical, intent(out) :: exitus

    integer, parameter :: rp = selected_real_kind(15)

    integer :: status,nrows,naper,kend,k
    real, dimension(:), allocatable :: xcens,ycens,sky,skyerr,skycorr,skyerrcorr
    real, dimension(:,:), allocatable :: apcts,apcts_err
    real(rp), dimension(:), allocatable :: gcount,gcount_err,area
    real(rp), dimension(:,:), allocatable :: ee,eeerr
    integer, dimension(:), allocatable :: growflag, nused
    logical, dimension(:), allocatable :: used
    real(rp), dimension(:), allocatable :: curve, curve_err, prof
    real, dimension(:), allocatable :: raper
    real, dimension(2) :: ring
    integer, dimension(2) :: naxes
    real :: hwhm, rflux90, ghwhm, sep, ecc, fwhm, fwhm_ee
    real(rp) :: errgain
    type(grow_reporter) :: reporter
    type(EEDescr) :: desc
    logical :: report_init, ok

    report_init = .false.
    status = 0

    call fits_gphot_read(filename,hwhm,sep,raper,xcens,ycens, &
         sky,skyerr,skycorr,skyerrcorr,apcts,apcts_err,status, &
         ring,ecc,fwhm,naxes,errgain)
    if( status /= 0 ) goto 666

    naper = size(raper)
    nrows = size(xcens)
    allocate(gcount(nrows),gcount_err(nrows), growflag(nrows), &
         curve(naper), curve_err(naper), prof(naper) )

    if( method == 'LEGACY' ) then

       ! The growth-curve routines use the width as the Gaussian sigma
       ! (their `hwhm' argument), while the HWHM keyword is a true HWHM.
       hwhm = hwhm / sqrt(2*log(2.0))

       if( verbose ) then
          report_init = .true.
          call grow_report_init(reporter,naper,nrows,raper)
       end if

       ! growth-curve photometry
       call growphot(xcens,ycens,apcts,apcts_err,sky,skyerr,raper,hwhm,sep,&
            threshold, &
            gcount,gcount_err,skycorr,skyerrcorr,curve,curve_err,growflag, &
            prof,rflux90,ghwhm,verbose,reporter)

       ! update sky: the growth-curve fit models the aperture counts
       ! as f*curve + b*area, so b (skycorr) is the residual background
       ! per pixel left after subtraction of the sky; it is added back.
       ! The error of the mean sky of the ring cancels in sky + b (b absorbs
       ! it), so the error of the corrected sky is the error of b by the
       ! noise of the pixels of the apertures (growfit); it was combined with
       ! the error of the ring in quadrature (24-41% too large).
       sky = sky + skycorr
       where( skyerr > 0 .and. skyerrcorr > 0 )
          skyerr = skyerrcorr
       end where

    else

       ! the growth curve: the enclosed energy of the apertures by the
       ! bright isolated stars (the plateau, the variation over the field)
       desc%source = 'AUTO'
       desc%predict = .true.
       desc%ordreq = eevar
       desc%radial = eeradial
       if( eecen_set ) then
          desc%xcen = eecen(1)
          desc%ycen = eecen(2)
       else
          desc%xcen = naxes(1) / 2
          desc%ycen = naxes(2) / 2
       end if
       if( fwhm > 0 ) then
          fwhm_ee = fwhm
       else
          fwhm_ee = 2*hwhm
       end if
       allocate(ee(nrows,naper),eeerr(nrows,naper),used(nrows),nused(nrows), &
            area(naper))
       call apee_auto(desc,xcens,ycens,apcts,apcts_err,skyerr,raper,ring, &
            ecc,fwhm_ee,naxes,verbose,ee,eeerr,ok,used)
       if( .not. ok ) then
          write(error_unit,*) 'Error: ',trim(filename), &
               ': the growth curve failed (too few bright isolated stars).'
          status = 1
          goto 666
       end if

       ! the apertures to the end of the plateau
       kend = minloc(abs(raper - real(desc%rref2)),1)
       area = 4*atan(1.0_rp)*real(raper,rp)**2*sqrt(1 - real(ecc,rp)**2)

       call grow_ee_phot(area,kend,errgain,apcts,apcts_err,skyerr,ee,eeerr, &
            gcount,gcount_err,skycorr,skyerrcorr,nused)

       ! the sky of the star by the ring and the apertures
       where( skyerrcorr > 0 )
          sky = sky + skycorr
          skyerr = skyerrcorr
       end where

       where( used )
          growflag = 1
       elsewhere
          growflag = 0
       end where

       ! the growth curve of the frame, the profile, HWHM and the radius
       ! of 90% of the total flux
       call grow_ee_curve(raper,area,ee,eeerr,curve,curve_err,prof)
       block
         real(rp), dimension(naper) :: c
         c = curve
         do k = 1, naper
            if( c(k) < 0 ) c(k) = 1
         end do
         rflux90 = estim_rflux90(raper,c)
       end block
       ghwhm = estim_hwhm(raper,prof)
       if( .not. (ghwhm > 0) ) ghwhm = hwhm

       if( verbose ) then
          write(*,'(a,i0,a,f0.2,a,f0.2,a)') ' Growth curve: ',count(used), &
               ' stars, total flux by the apertures ',desc%rref1,' .. ', &
               desc%rref2,' pix.'
          write(*,'(a,f0.2,a,f0.2,a)') ' Growth curve: HWHM = ',ghwhm, &
               ' pix, 90% of the flux within ',rflux90,' pix.'
          write(*,'(a,i0,a,i0,a)') ' GCOUNT of ',count(gcount > 0),' of ', &
               nrows,' stars; the stars by the largest aperture of GCOUNT:'
          do k = 1, kend
             write(*,'(a,i3,f8.2,a,i6)') '  aperture ',k,raper(k),' pix:', &
                  count(nused == k)
          end do
       end if

    end if

    call fits_find_save(filename,output,ghwhm,rflux90, &
       raper, xcens, ycens, sky,skyerr, gcount,gcount_err, growflag, &
       curve, curve_err, prof, status, method)

    if( method /= 'LEGACY' .and. status == 0 ) then
       if( output == '' ) then
          call apee_header(filename,desc,status,GROWEXTNAME)
       else
          call apee_header(output,desc,status,GROWEXTNAME)
       end if
    end if

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
