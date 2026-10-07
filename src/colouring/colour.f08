!
! the colouring book -- the engine
!
!
! Copyright © 2010-2026 F.Hroch (hroch@physics.muni.cz)
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

module colouring_book

  use iso_fortran_env

  implicit none
  private

  ! Planck's constant * speed of light [in electron-Volts]
  real, parameter :: hcevolt = 1.2398e-06 ! = h c / 1eV

  ! D65 white point for XYZ values for 2 deg observer
  ! https://en.wikipedia.org/wiki/Illuminant_D65
!  real, dimension(3), parameter :: white_point = [ 1.07304, 1.0, 0.948115 ]
  real, dimension(3), parameter :: white_point = [ 1.08883, 1.0, 0.95047 ]

  integer, parameter, public :: &
       WHITE_SPOT = 5, &         ! use a spot to calibrate white
       WHITE_STAR = 4, &         ! use a star to calibrate white
       WHITE_CTPH = 3, &         ! use quantum efficiency, aka CTPH
       WHITE_CAL = 2, &          ! use calibrated frames
       WHITE_WEIGHTS = 1         ! use weights by user to calibrate white

  public :: colour

contains

  subroutine colour(cspace,outspace,bands,bitpix,cmatrix,x,y,r, &
       white,estimbacks,weights,ctphs,backgrounds,output,verbose,status)

    use fitscolour
!    use phio
    use oakleaf
    use titsio

    character(len=*), intent(in) :: cspace, outspace, output
    type(ColourFits) , dimension(:), intent(in) :: bands
    integer, intent(in) :: bitpix, white
    real, dimension(:,:), intent(in) :: cmatrix
    real, intent(in) :: x, y, r
    real, dimension(:), allocatable, intent(in) :: weights, ctphs, backgrounds
    logical, intent(in) :: estimbacks, verbose
    integer, intent(out) :: status

    integer :: nbands, width, height
    integer :: i,j,k,stat

    character(len=FLEN_CARD) :: com
    character(len=80) :: msg, squnit
    integer, dimension(3) :: naxes
    real, dimension(:,:,:), allocatable :: cube
    real, dimension(size(bands)) :: e, w, cts, sq, backs, fvega
    logical :: astrometry
    type(fitsfiles) :: fits

    status = 1
    nbands = size(bands)

    width = bands(1)%naxes(1)
    height = bands(1)%naxes(2)

    astrometry = all(bands(:)%astrometry)
    e = hcevolt / bands(:)%leff           ! effective photon energy in eV
    if( astrometry ) then
       sq = (3600*bands(:)%scale)**2     ! scale arcsec2 per pixels
    else
       sq = 1
    end if

    ! backgrounds
    if( allocated(backgrounds) ) then
       backs = backgrounds
    else if( estimbacks ) then
       do k = 1, size(bands)
          backs(k) = bands(k)%sky()
          if( verbose ) write(error_unit,*) "Info: background for `", &
               trim(bands(k)%filename),"' estimated on: ",backs(k)
       end do
    else ! backgrounds are set to zero
       backs = 0
    end if

    if( white == WHITE_SPOT ) then
       ! determine weights by averadging of white area
       do k = 1, nbands
          cts(k) = bands(k)%apmean(x,y,r,backs(k))
          if( verbose ) then
             write (error_unit,*) 'Info: white spot mean in ',&
                  trim(bands(k)%filter),' is: ',cts(k)
          end if
       end do
       if( .not. all(cts > 0) ) &
            stop 'Error: white spot mean determination failed.'

       call weights_setup(cts,cmatrix,w)
!       w = w / w(2)

    else if ( white == WHITE_STAR ) then
       ! determine weights by equalising of flux of a white star

       do k = 1, nbands
          cts(k) = bands(k)%apsum(x,y,r,backs(k))
          if( verbose ) then
             write (error_unit,*) 'Info: white star in ',&
                  trim(bands(k)%filter),' has total counts in aperture: ',cts(k)
          end if
       end do
       if( .not. all(cts > 0) ) &
            stop 'Error: white star determination failed.'

       ! NEEDS TO BE ADJUSTED PROPERLY !!!!!!!!!!!!!!!!!!!!!!!!!!1
       goto 33
       block
         real, dimension(3), parameter :: mags = [ 12.57, 11.95, 11.58 ]
         real, dimension(3) :: flx = 1
         real, parameter :: bv = mags(1) - mags(2)
         real, parameter :: vr = mags(2) - mags(3)
         character(len=666) :: develop_the_code_PLEASE

         flx(1) = 10**(0.4*bv)
         flx(3) = 10**(0.4*vr)
         cts = cts / flx
       end block

       fvega = [1.162, 1.041, 1.052 ]
       cts = cts / fvega
!       write(error_unit,*) cts / cts(2)

       call weights_setup(cts*e,cmatrix,w)
!       call weights_setup(cts/cts(2)*e,cmatrix,w)
!       w = w / w(2)
33     continue

       call weights_setup(cts*e,cmatrix,w)
!       w = w / w(2)

    else if( white == WHITE_WEIGHTS ) then
       ! by user

       w = weights

    else if( white == WHITE_CTPH .or. white == WHITE_CAL ) then
       ! passed CTPH or photon rates are already calibrated

       ! NEEDS REVISION !!!!

!       write(*,*) bands(:)%fref / bands(2)%fref
!       write(*,*) bands(:)%exptime / bands(2)%exptime
!       write(*,*) bands(:)%ctph / bands(2)%ctph
!       write(*,*) bands(:)%ctph * bands(:)%fref / (bands(2)%ctph*bands(2)%fref)
!       write(*,*) e / e(2)

!       f = 1 * (bands(:)%exptime * bands(:)%area)
!       f = bands(2)%fref / bands(:)%fref
!       f = bands(:)%exptime  / bands(2)%exptime
       cts = (bands(:)%ctph * bands(:)%fref) / (bands(2)%fref * bands(2)%ctph)
!       f = bands(:)%fref !/ bands(:)%exptime
!       f = e * (bands(:)%exptime * bands(:)%area) !/ (bands(2)%exptime*bands(2)%area)
!       e = 1
       !       f = f / f(2)
!       f = 1 / [ 3.166, 1.0, 2.0]
!       f = [ 0.28, 1.0, 0.89 ]
       !       write(*,*) f
       cts = 1 / [ 1.4, 1.0, 0.5 ]
       cts = [ 1.43E+10,  8.60E+09, 1.32E+10 ]  ! A0, photon flux
!       cts = [ 7.09E+05,  7.03E+05,  1.59E+06 ]
       ! mean density flux in the filter
!       cts = [ 1.64, 9.26, 7.09 ]  ! A0
!       cts = [  7.74, 7.26, 8.45 ] ! F5
       cts = [ 5.9, 6.393, 10.10 ] ! D65
       !cts = [ 1.15, 1.24, 1.95 ]
!       cts = [ 7.74, 7.26, 8.45 ]
!       cts = [ 3.18E-13, 2.53E-13, 4.62E-13 ]
!       cts = [ 6.50E-09, 3.12E-09, 3.90E-09 ]
!       cts = 1
!       cts = [ 1.43E+10,  8.60E+09, 1.32E+10 ] / [ 7.09E+05,  7.03E+05,  1.59E+06 ]
!       f = 1 / [ 1.6, 1.0, 1.53 ]
       !       f = 1 / [ 3.3, 3.0, 2.3 ]
       call weights_setup(cts,cmatrix,w)
       w = w / (bands(:)%exptime * bands(:)%area / (bands(2)%exptime*bands(2)%area))
!       w = w * (bands(:)%exptime * bands(:)%area)
!       if( white == WHITE_CTPH ) w = w * ctphs
!       e = hcevolt / bands(:)%leff
!       w = w * (bands(:)%exptime * bands(:)%area)
!       w = w / w(2)
    else
       ! guesstimate, no data for adjusting available

       cts = 1
       call weights_setup(cts*e,cmatrix,w)
!       w = w / w(2)
       write(error_unit,*) "Warning: weights are undefined, colours becomes false."

    end if

    if( verbose ) then
       write(error_unit,'(a)') &
            ' Info: filter, wavelength[nm], energy[eV], background[cts], weight, exptime[s]:'
       do k = 1, nbands
          write(error_unit,'(2a,3x,f7.1,1x,f5.2,3x,en12.2,f12.2,f10.3)') " Info: ", &
               trim(bands(k)%filter),1e9*bands(k)%leff,e(k),backs(k),w(k),bands(k)%exptime
       end do
    end if

    if( .not. any(w > 0) ) &
         stop 'Error: all channel weights, and thinking, must be positive.'

    ! memory
    allocate(cube(width,height,nbands),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) "Error: ",trim(msg)
       error stop 'Failed to allocate the colour cube.'
    end if

    ! scale inputs, convert to fluxes
    do k = 1, nbands
!       cube(:,:,k) = w(k) * (bands(k)%image - backs(k)) * e(k) / sq(k)
       cube(:,:,k) = w(k) * (bands(k)%image - backs(k))
    end do

    ! the rotation in the colour space, the transformation
    forall( i = 1:width, j = 1:height )
        cube(i,j,:) = matmul(cmatrix,cube(i,j,:))
    end forall

    ! cut to the appropriate range
    block
      real :: cmax

      if( bitpix > 0 ) then
         cmax = 2.0**bitpix - 1
         cube = max(0.0,min(cmax,cube))
      else
         cube = max(0.0,cube)
      end if

    end block

    ! write out Colour FITS
    status = 0
    if( fits_file_exist(output) ) call fits_file_delete(output)
    call fits_create_file(fits,output,status)
    if( status /= 0 ) then
       write(error_unit,*) 'Error: failed to create the file `',trim(output),"'."
       return
    end if

    naxes = [ width, height, nbands ]
    call fits_insert_img(fits,bitpix,3,naxes,status)

    ! by V filter
    do k = 1, nbands
       if( bands(k)%filter == 'V' ) then

          if( bands(k)%object /= '' ) &
               call fits_update_key(fits,FITS_KEY_OBJECT,bands(k)%object, &
               'object by V filter',status)

          call fits_update_key(fits,FITS_KEY_DATEOBS,bands(k)%dateobs, &
               'date by V filter',status)

          call fits_update_key(fits,FITS_KEY_EXPTIME,bands(k)%exptime,-3,&
               '[s] original exposure by V filter',status)

          call fits_update_key(fits,FITS_KEY_AREA,bands(k)%area,-3, &
               '[m2] original area',status)

          if( bands(k)%astrometry ) then
             call fits_update_wcs(fits,bands(k)%radesys,bands(k)%epoch, &
                  bands(k)%ctype,bands(k)%crval,bands(k)%crpix,bands(k)%cd, &
                  real([-1,-1],REAL64),status)

             ! adds spectral components
             call fits_update_key(fits,'CTYPE3','WAVE','spectral coordinate type',status)
             call fits_update_key(fits,'CD3_3',150,'',status)
             call fits_update_key(fits,'CRVAL3',380,'[nm] Z component',status)
             call fits_update_key(fits,'CRPIX3',1,'[pix] Z component',status)

          end if
       end if
    end do

    call fits_update_key(fits,FITS_KEY_CSPACE,outspace, &
         'the colour space of stored data',status)
    call fits_write_comment(fits, &
         "Original colour-space: '"//trim(cspace)//"'",status)
    call fits_write_comment(fits,'Original file, weight, background:',status)
    do k = 1, nbands
       write(com,'(a,2x,g0.3,1x,g0.5)') trim(bands(k)%filename),w(k),backs(k)
       call fits_write_comment(fits,com,status)
    end do

    if( white == WHITE_CAL .or. white == WHITE_CTPH ) then
       if( astrometry ) then
          squnit = 'arcsec2'
       else
          squnit = 'pixel2'
       end if
       call fits_update_key(fits,FITS_KEY_BUNIT,'eV/s/m2/'//squnit, &
            'image data represents intensity (energy)',status)
    end if

    call fits_write_comment(fits,MUNIPACK_VERSION,status)
    call fits_update_key(fits,'CREATOR','Munipack', &
         'Created by colouring utility of Munipack',status)
    call fits_write_cube(fits,0,cube,status)
    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    deallocate(cube)

  end subroutine colour

  subroutine weights_setup(flx,cmatrix,w)

    !    use minpacks
    use oakleaf

    real, dimension(:), intent(in) :: flx
    real, dimension(:,:), intent(in) :: cmatrix
    real, dimension(:), intent(out) :: w

    real(REAL64), dimension(size(flx)) :: x

!    call qrsolve(real(cmatrix,REAL64),real(white_point,REAL64),x)
    call qrsol(real(cmatrix,REAL64),real(white_point,REAL64),x)
!    write(error_unit,*) 'x:',x / x(2)
!    write(error_unit,*) 'flx:',flx / (cts(2)*e(2))
!    write(error_unit,*) 'flx:',flx / (flx(2))

    w = real( x / flx )
!    w = real(x / (cts*e))
!    w = real( x / (cts/e))
!    w = x

!    write(*,*) 'cts:',cts
!    write(*,*) 'e:',e

    w = w / w(2)

  end subroutine weights_setup

end module colouring_book
