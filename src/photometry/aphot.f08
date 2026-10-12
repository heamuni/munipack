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
  use aperture_ee
  use skymean, only: sky_noise_model
  use psfmodel, only: psf_file_load

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

  ! the enclosed energy of apertures (aperture corrections): by a PSF
  ! (eemode = 'PSF', the PSF file or an analytic one) or by the frame's
  ! stars (eemode = 'AUTO')
  character(len=FLEN_VALUE) :: eemode = '', psftype = ''
  character(len=FLEN_FILENAME) :: psffile = ''
  real(REAL64) :: fwhm_par = -1, beta_par = 2.5, scale_par = -1
  real(REAL64), dimension(2) :: psfcen_par = -1
  logical :: psfcen_set = .false.
  integer :: eevar = 2
  logical :: eeradial = .false.
  real(REAL64), dimension(:,:,:), allocatable :: stamps
  real(REAL64), dimension(:), allocatable :: rcoo, acoo
  integer :: zoom = 1, nbox = 0

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

     else if( key == 'PSF' ) then

        read(val,*) psftype
        if( .not. (psftype == 'GAUSS' .or. psftype == 'MOFFAT' .or. &
             psftype == 'FILE') ) stop 'Error: PSF is GAUSS, MOFFAT or FILE.'
        eemode = 'PSF'

     else if( key == 'PSFFILE' ) then

        read(val,*) psffile
        psftype = 'FILE'
        eemode = 'PSF'

     else if( key == 'FWHM' ) then

        read(val,*) fwhm_par

     else if( key == 'BETA' ) then

        read(val,*) beta_par

     else if( key == 'SCALE' ) then

        read(val,*) scale_par

     else if( key == 'PSFCENTRE' ) then

        read(val,*) psfcen_par
        psfcen_set = .true.

     else if( key == 'APCORR' ) then

        read(val,*) eemode
        if( eemode /= 'AUTO' ) stop 'Error: APCORR is AUTO.'

     else if( key == 'APCORRRADIAL' ) then

        read(val,*) eeradial

     else if( key == 'APCORRVAR' ) then

        read(val,*) eevar
        if( eevar < 0 .or. eevar > 2 ) stop 'Error: APCORRVAR is 0, 1 or 2.'

     else if( key == 'FILE' ) then

        read(val,*) filename, output

        if( verbose ) write(*,*) trim(filename)//": "

        if( psftype == 'FILE' .and. .not. allocated(stamps) ) then
           if( psffile == '' ) stop 'Error: PSF=FILE requires a PSF file.'
           stat = 0
           call psf_file_load(psffile,stamps,rcoo,acoo,zoom,nbox,verbose,stat)
           if( stat /= 0 ) stop 'Error: failed to read the PSF file.'
        end if

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
    real :: lobad,hibad,fwhm,hwhm,ecc,incl,fwhm_ee
    integer :: nrows, status
    real(REAL64), dimension(:,:), allocatable :: ee, eeerr
    real(REAL64) :: scale, gain, ngain, nrn2, nbias
    real, dimension(:), allocatable :: rmask
    integer, dimension(2) :: naxes
    type(EEDescr) :: desc
    logical :: eeok

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
         sky(nrows),sky_err(nrows),rmask(nrows))

    ! the sky: the clipped mean by the noise model, the detected stars
    ! masked
    call sky_noise_model(filename,ngain,nrn2,nbias,status)
    if( status /= 0 ) goto 666
    call daomaskrad(data,xcens,ycens,lobad,hibad,ngain,nrn2,nbias,ring(2)/2, &
         rmask)
    call daophotsb(data,stderr,xcens,ycens,raper,ring,ecc,incl,lobad,hibad, &
         verbose,plog,apcts,apcts_err,sky,sky_err,ngain,nrn2,nbias,rmask)

    ! estimate of width parameter for Gaussian
    call estim_hwhm(data,xcens,ycens,sky,fwhm,lobad,hibad,hwhm)
    if( hwhm < 0 ) hwhm = fwhm / 2

    ! the enclosed energy of the apertures
    eeok = .false.
    if( eemode /= '' ) then
       allocate(ee(nrows,naper),eeerr(nrows,naper))
       call fits_aphot_info(filename,naxes,scale,gain,status)
       if( status /= 0 ) goto 666
       if( eemode == 'PSF' ) then
          desc%source = 'PSF'
          desc%psftype = psftype
          if( psftype == 'FILE' ) then
             desc%psffile = psffile
             if( psfcen_set ) then
                desc%xcen = psfcen_par(1)
                desc%ycen = psfcen_par(2)
             else
                ! the same reference pixel as artificial and psfphot
                desc%xcen = naxes(1) / 2
                desc%ycen = naxes(2) / 2
             end if
             if( scale_par > 0 ) then
                desc%scale = scale_par
             else
                desc%scale = scale
             end if
             if( .not. (desc%scale > 0) .and. size(rcoo) > 1 ) then
                write(error_unit,*) 'Error: the scale [pix/deg] is undefined ', &
                     '(no WCS in the header): use --scale.'
                status = 1
                goto 666
             end if
          else
             if( fwhm_par > 0 ) then
                desc%fwhm = fwhm_par
             else if( fwhm > 0 ) then
                desc%fwhm = fwhm
             else
                desc%fwhm = 2*hwhm
             end if
             desc%beta = beta_par
             if( psftype == 'MOFFAT' .and. .not. (beta_par > 1) ) &
                  stop 'Error: Moffat beta > 1 required.'
             allocate(stamps(0,0,0),rcoo(0),acoo(0))
          end if
          call apee_psf(desc,stamps,rcoo,acoo,zoom,nbox,xcens,ycens,raper, &
               ecc,incl,ee,eeerr)
          if( psftype /= 'FILE' ) deallocate(stamps,rcoo,acoo)
          eeok = .true.
       else
          desc%source = 'AUTO'
          desc%ordreq = eevar
          desc%radial = eeradial
          if( psfcen_set ) then
             desc%xcen = psfcen_par(1)
             desc%ycen = psfcen_par(2)
          else
             desc%xcen = naxes(1) / 2
             desc%ycen = naxes(2) / 2
          end if
          if( fwhm > 0 ) then
             fwhm_ee = fwhm
          else
             fwhm_ee = 2*hwhm
          end if
          call apee_auto(desc,xcens,ycens,apcts,apcts_err,raper,ring,fwhm_ee, &
               gain,naxes,verbose,ee,eeerr,eeok)
          if( .not. eeok ) write(error_unit,*) 'Warning: ',trim(filename), &
               ': the aperture correction failed, APEE not saved.'
       end if
    end if

    if( eeok ) then
       call fits_aphot_save(filename, output, hwhm, ecc, incl, raper, ring, &
            xcens, ycens, apcts,apcts_err,sky,sky_err, status, ee, eeerr)
       if( output == '' ) then
          call apee_header(filename,desc,status)
       else
          call apee_header(output,desc,status)
       end if
    else
       call fits_aphot_save(filename, output, hwhm, ecc, incl, raper, ring, &
            xcens, ycens, apcts,apcts_err,sky,sky_err, status)
    end if

666 continue

    exitus = status == 0

    if( allocated(data) ) deallocate(data,stderr)
    if( allocated(apcts) ) deallocate(xcens,ycens,apcts,apcts_err,sky,sky_err)
    if( allocated(rmask) ) deallocate(rmask)

  end subroutine the_aphot

  subroutine byhand(exitus)

    use mdaosky

    logical, intent(out) :: exitus

    real, dimension(:,:), allocatable :: data, stderr
    real, dimension(:), allocatable :: sky, sky_err
    real, dimension(:,:), allocatable :: apcts,apcts_err
    real, dimension(:), allocatable :: rmask
    real :: lobad, hibad, saturate, skyavg,skyerr,skysig
    real(REAL64) :: ngain, nrn2, nbias
    integer :: i, nrow, nstep, status

    status = 0
    call fits_aphot_image(filename,fkeys,data,stderr,saturate,status)
    call sky_noise_model(filename,ngain,nrn2,nbias,status)
    if( status /= 0 ) goto 666

    hibad = saturate
    nstep = max(int(log(float(size(data))/4e4) / log(2.0)),1)
    call daosky(data,nstep,verbose,hibad,skyavg,skyerr,skysig,ngain,nrn2,nbias)
    ! the same low-bad threshold as in find (lothresh = 5 sigma), not clamped
    ! to zero (it would cut half of the sky on sky-subtracted frames)
    lobad = skyavg - 5*skysig

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

    allocate(rmask(nrow))
    call daomaskrad(data,xstar,ystar,lobad,hibad,ngain,nrn2,nbias,ring(2)/2, &
         rmask)
    call daophotsb(data,stderr,xstar,ystar,raper,ring,ecc_par,incl_par, &
         lobad,hibad,verbose,plog,apcts,apcts_err,sky,sky_err, &
         ngain,nrn2,nbias,rmask)

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
