!
!  aperture photometry
!
!  Copyright © 2010-19 F.Hroch (hroch@physics.muni.cz)
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
!  Ideas:
!    * The aperture radii distribution by [1 + erf(x/FWHM)]/2.
!

program aphot

  use titsio
  use fitsaphot
  use daofotometr

  implicit none

  real, parameter :: pi = 3.141592653

  character(len=4*FLEN_FILENAME) :: record,key,val
  character(len=FLEN_FILENAME) :: filename,output
  character(len=FLEN_KEYWORD), dimension(1) :: fkeys
  character(len=80) :: msg
  real, dimension(:), allocatable :: raper, xstar, ystar
  real, dimension(2) :: ring
  real ::  ecc_par = 0           ! eccentricity: 0 (circle) .. 1 (line)
  real :: incl_par = 0           ! inclination of major semiaxis, degs
  logical :: ellipticity = .false. ! enable elliptic apertures
  logical :: snap = .false.      ! locate a local maximum brightness (by hand)
  logical :: ecc_defined = .false., incl_defined = .false.
  logical :: verbose = .false., plog = .false.
  logical :: ex, exitus = .true.
  integer :: i,eq,stat,naper,nstars
  real :: theta, pitch

  ! Aperture radii are defined by the logarithmic spiral
  ! (https://en.wikipedia.org/wiki/Logarithmic_spiral); the angles
  ! are in the interval 0 .. pi with step pi / (naper - 1). The pitch
  ! parameter is ln(20)/pi; it is computed from the inner ring radius.
  ! pitch 9.84 deg, raper(1) = 1, .. raper(12) = 20
  naper = 12
  allocate(raper(naper))
  pitch = log(20.0) / pi
  do i = 1, naper
     theta = (i - 1.0) / (naper - 1.0) * pi
     raper(i) = exp(pitch*theta)
  end do
  ring = [20.0, 30.0]

  nstars = 0
  allocate(xstar(0),ystar(0))

  fkeys(1) = FITS_KEY_SATURATE

  do
     read(*,'(a)',iostat=stat,iomsg=msg) record
     if( stat == IOSTAT_END ) exit
     if( stat > 0 ) then
        write(error_unit,*) trim(msg)
        error stop 'Some input error.'
     end if

     eq = index(record,'=')
     if( eq == 0 ) stop 'Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'PIPELOG' ) then

        read(val,*) plog

     else if( key == 'FITS_KEY_SATURATE' ) then

        read(val,*) fkeys(1)

     else if( key == 'NAPER' ) then

        read(val,*) naper
        deallocate(raper)
        allocate(raper(naper))

     else if( key == 'APER' ) then

        read(val,*) raper

     else if( key == 'RING' ) then

        read(val,*) ring

     else if( key == 'ELLIPTICITY' ) then

        read(val,*) ellipticity

     else if( key == 'ECCENTRICITY' ) then

        read(val,*) ecc_par
        ecc_defined = .true.
        if( .not. abs(ecc_par) < 1 ) stop 'Eccentricity >= 1.'

     else if( key == 'INCLINATION' ) then

        read(val,*) incl_par
        incl_defined = .true.

     else if( key == 'NSTARS' ) then

        read(val,*) nstars
        deallocate(xstar,ystar)
        allocate(xstar(nstars),ystar(nstars))
        nstars = 0

     else if( key == 'STAR' ) then

        nstars = nstars + 1
        read(val,*) xstar(nstars),ystar(nstars)

     else if( key == 'SNAP' ) then

        read(val,*) snap

     else if( key == 'FILE' ) then

        read(val,*) filename, output

        if( verbose ) write(*,*) trim(filename)//": "

        if( nstars > 0 ) then
           call byhand(exitus)
        else
           call the_aphot(ex)
           exitus = exitus .and. ex
        end if

     end if

  end do

  deallocate(raper,xstar,ystar)

  if( exitus ) then
     stop 0
  else
     stop 'Some error(s) occurred during this APHOT run.'
  end if

contains

  subroutine the_aphot(exitus)

    logical, intent(out) :: exitus

    real, dimension(:,:), allocatable :: data,stderr
    real, dimension(:), allocatable :: xcens,ycens,sky,sky_err
    real, dimension(:,:), allocatable :: apcts,apcts_err
    real :: lobad,hibad,fwhm,hwhm,ecc,incl
    integer :: nrows, status

    status = 0
    call fits_aphot_read(filename,data,stderr,xcens,ycens, &
         lobad,hibad,fwhm,ecc,incl,status)
    if( status /= 0 ) goto 666

    if(  ecc_defined )  ecc =  ecc_par
    if( incl_defined ) incl = incl_par

    if( .not. ellipticity ) then
       ecc = 0
       incl = 0
    end if

    nrows = size(xcens)
    allocate(apcts(nrows,naper),apcts_err(nrows,naper), &
         sky(nrows),sky_err(nrows))

    call daophotsb(data,stderr,xcens,ycens,raper,ring,ecc,incl,lobad,hibad, &
         verbose,plog,apcts,apcts_err,sky,sky_err)

    ! estimate of width parameter for Gaussian
    call estim_hwhm(data,xcens,ycens,sky,fwhm,lobad,hibad,hwhm)
    if( hwhm < 0 ) hwhm = fwhm / 2

    call fits_aphot_save(filename, output, hwhm, ecc, incl, raper, ring, &
         xcens, ycens, apcts,apcts_err,sky,sky_err, status)

666 continue

    exitus = status == 0

    if( allocated(data) ) deallocate(data,stderr)
    if( allocated(apcts) ) deallocate(xcens,ycens,apcts,apcts_err,sky,sky_err)

  end subroutine the_aphot

  subroutine byhand(exitus)

    use mdaosky

    logical, intent(out) :: exitus

    real, dimension(:,:), allocatable :: data, stderr
    real, dimension(:), allocatable :: sky, sky_err
    real, dimension(:,:), allocatable :: apcts,apcts_err
    real :: lobad, hibad, saturate, skyavg,skyerr,skysig
    integer :: i, nrow, nstep, status

    status = 0
    call fits_aphot_image(filename,fkeys,data,stderr,saturate,status)
    if( status /= 0 ) goto 666

    hibad = saturate
    nstep = max(int(log(float(size(data))/4e4) / log(2.0)),1)
    call daosky(data,nstep,verbose,hibad,skyavg,skyerr,skysig)
    lobad = max(0.0,skyavg - 3*skysig)

    nrow = size(xstar)
    allocate(apcts(nrow,naper),apcts_err(nrow,naper),sky(nrow),sky_err(nrow))

    if( snap ) then
       block
         integer :: i,j,n,i0,j0,imin,imax,jmin,jmax
         real :: fmax
         logical :: found

         do n = 1, nrow

            imin = max(int(xstar(n)-raper(1)-1),1)
            jmin = max(int(ystar(n)-raper(1)-1),1)
            imax = min(int(xstar(n)+raper(1)+1),size(data,1))
            jmax = min(int(ystar(n)+raper(1)+1),size(data,2))
            found = .false.
            fmax = lobad
            do i = imin, imax
               do j = jmin, jmax
                  if( data(i,j) > fmax ) then
                     fmax = data(i,j)
                     i0 = i
                     j0 = j
                     found = .true.
                  end if
               end do
            end do

            if( found ) then
               xstar(n) = i0
               ystar(n) = j0
            end if

         end do
       end block
    end if

    call daophotsb(data,stderr,xstar,ystar,raper,ring,ecc_par,incl_par, &
         lobad,hibad,verbose,plog,apcts,apcts_err,sky,sky_err)

    do i = 1, nrow
       write(*,*) '=BYHAND> ',i,xstar(i),ystar(i),sky(i),sky_err(i), &
            apcts(i,:),apcts_err(i,:)
    end do

666 continue

    exitus = status == 0

    if( allocated(data) ) deallocate(data,stderr)
    if( allocated(apcts) ) deallocate(apcts,apcts_err,sky,sky_err)

  end subroutine byhand


end program aphot
