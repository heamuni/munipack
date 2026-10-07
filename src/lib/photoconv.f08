!
!  photoconv - photometry conversions
!
!  Copyright © 2013 - 15, 2017-2021 F.Hroch (hroch@physics.muni.cz)
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


module photoconv

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)

  ! http://en.wikipedia.org/wiki/Planck_constant
  real(dbl), parameter, private :: planck = 6.62606957e-34_dbl    ! [J.s] in SI

  ! http://en.wikipedia.org/wiki/Speed_of_light
  real(dbl), parameter, private :: c = 299792458.0_dbl            ! [m/s] in SI

  ! reciprocal flux for star of magnitude 25
  real(dbl), parameter, private :: flux25 = 1e10

  ! spectral flux density for star of magnitude 0 in V (545nm) band for AB system
  real(dbl), parameter, private :: ABspflux = 3.631e-23 ! [W/m2/Hz]

  ! spectral flux density for star of magnitude 0 for ST system
  ! (3.6335e-9 erg/s/cm2/A, STMAG = -2.5 log10(FLAM) - 21.10 in cgs)
  real(dbl), parameter, private :: STspflux = 3.6335e-11 ! [W/m2/nm]

  ! gauss-hermite coefficients
  real(dbl), parameter, private :: h0 = 1.1816359006_dbl, &
       h12 = 0.295408975151_dbl, &
       a = 1.22474487139_dbl, &
       a1 = 0.707106781187_dbl
  ! h*c (planck * light speed in [J]
  real(dbl), parameter, private :: hc = 1.98644568e-25

  ! electron-Volt
  real(dbl), parameter, private :: evolt = 1.602176634e-19

  ! sqrt(pi)
  real(dbl), parameter, private :: sqrtpi = 1.772453850905515882_dbl

  ! sqrt(pi/2)
  real(dbl), parameter, private :: sqrtpi2 = 1.2533141373155001_dbl

  ! 2.5 / ln(10), the factor converting relative errors to magnitudes
  real(dbl), parameter, private :: magerr = 1.0857362047581294_dbl

  ! Passband width convention
  ! -------------------------
  ! The photometric system table provides LAM_EFF and LAM_FWHM.
  ! The passband is approximated by a Gaussian with the standard deviation
  ! sigma = LAM_FWHM / 2 (LAM_FWHM is treated as 2*sigma, not as the true
  ! FWHM = 2.355*sigma of the Gaussian). Its effective (equivalent) width,
  ! the integral of the normalised transmission, is
  !
  !      W = sqrt(2*pi) * sigma = sqrt(pi/2) * LAM_FWHM = 1.2533 * LAM_FWHM
  !
  ! (a Gaussian with the true FWHM would give 1.0645 * FWHM). The calibration
  ! (mag2ph0, ph2mag0) relates the photon rate ph [ph/s/m2] and the zero
  ! magnitude spectral flux FLAM_REF [W/m2/nm] by
  !
  !      ph = FLAM_REF * 10**(-0.4*mag) * (W/1nm) * LAM_EFF / (h*c)
  !
  ! All derived quantities use the same width W and only the wavelength
  ! columns (LAM_EFF, LAM_FWHM); the frequency width is derived
  ! as W_nu = c * W / LAM_EFF**2, NU_EFF = c / LAM_EFF:
  !
  !      PHOTLAM = ph / (W/1nm)                  [ph/s/m2/nm]
  !      PHOTNU  = ph / W_nu                     [ph/s/m2/Hz]
  !      FLUX    = ph * h*c / LAM_EFF            [W/m2] (written in eV/s/m2)
  !      FLAM    = ph * h*c / LAM_EFF / (W/1nm)  [W/m2/nm]
  !      FNU     = ph * h*c / LAM_EFF / W_nu     [W/m2/Hz]
  !      STMAG   = -2.5*log10(FLAM / 3.6335e-11 W/m2/nm)
  !      ABMAG   = -2.5*log10(FNU / 3.631e-23 W/m2/Hz)


  ! known quantities
  character(len=8), dimension(12), parameter :: quantities = [ &
       'COUNT   ', &
       'RATE    ', &
       'PHOTON  ', &
       'PHOTRATE', &
       'PHOTNU  ', &
       'PHOTLAM ', &
       'FLUX    ', &
       'FNU     ', &
       'FLAM    ', &
       'MAG     ', &
       'ABMAG   ', &
       'STMAG   ' ]

  ! units
  character(len=14), dimension(12), parameter :: units = [ &
       'count         ', &
       'count/s/m2    ', &
       'photon        ', &
       'photon/s/m2   ', &
       'photon/s/m2/Hz', &
       'photon/s/m2/nm', &
       'eV/s/m2       ', &
       'W/m2/Hz       ', &
       'W/m2/nm       ', &
       'mag           ', &
       'abmag         ', &
       'stmag         ']


contains

  subroutine quantity(key,unit)

    ! simulates searching in map stucture: ("key","unit")

    character(len=*), intent(in) :: key
    character(len=*), intent(out) :: unit
    integer :: i

    unit = ''
    do i = 1,size(quantities)
       if( key == quantities(i) ) then
          unit = units(i)
          return
       end if
    end do

  end subroutine quantity

  subroutine squantity(key,unit)

    ! quantity for surface photometry

    character(len=*), intent(in) :: key
    character(len=*), intent(out) :: unit

    call quantity(key,unit)
    unit = trim(unit)//'/arcsec2'

  end subroutine squantity

  subroutine ctrate(gain,area,exptime,arcscale,cts,dcts,rate,drate)

    real(dbl), intent(in) :: gain,area,exptime,arcscale
    real(dbl), dimension(:), intent(in) :: cts,dcts
    real(dbl), dimension(:), intent(out) :: rate,drate
    real(dbl) :: f

    f = gain/(area*exptime*arcscale**2)

    where( cts >= 0.0_dbl)
       rate = f * cts
       drate = f * dcts
    elsewhere
       rate = -1.0_dbl
       drate = -1.0_dbl
    end where

  end subroutine ctrate

  subroutine phrate(area,exptime,arcscale,ph,dph,rate,drate)

    ! photons (in CCD array) to photon rates

    real(dbl), intent(in) :: area,exptime,arcscale
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: rate,drate
    real(dbl) :: f

    f = area * exptime * arcscale**2

    where( ph >= 0.0_dbl)
       rate = ph / f
       drate = dph / f
    elsewhere
       rate = -1.0_dbl
       drate = -1.0_dbl
    end where

  end subroutine phrate

  subroutine ct2ph(gain,area,exptime,arcscale,ctph,dctph,cts,dcts,ph,dph)

    ! convert counts rate to photon rate [ph/s/m2]

    real(dbl), intent(in) :: gain,area,exptime,arcscale,ctph,dctph
    real(dbl), dimension(:), intent(in) :: cts,dcts
    real(dbl), dimension(:), intent(out) :: ph,dph
    real(dbl) :: f

    f = gain/(area*exptime)/arcscale**2/ctph

    where( cts >= 0.0_dbl)
       ph = f * cts
       dph = f * dcts
       !dph = sqrt(f**2*dctph**2 + f**2*dcts**2)
    elsewhere
       ph = -1.0_dbl
       dph = -1.0_dbl
    end where

  end subroutine ct2ph

  subroutine ph2photlam(lfwhm,ph,dph,q,dq)

    ! photon rate per unit of wavelegth [ph/s/m2/nm],
    ! ph / W, W = sqrt(pi/2)*lfwhm in nm (see the passband width convention)

    real(dbl), intent(in) :: lfwhm
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: q,dq
    real(dbl) :: f

    f = 1.0_dbl / (sqrtpi2 * lfwhm / 1e-9_dbl)

    where( ph >= 0.0_dbl )
       q = f * ph
       dq = f * dph
    elsewhere
       q = -1.0_dbl
       dq = -1.0_dbl
    end where

  end subroutine ph2photlam

  subroutine ph2photnu(leff,lfwhm,ph,dph,q,dq)

    ! photon rate per unit of frequency [ph/s/m2/Hz],
    ! ph / W_nu, W_nu = c * W / leff**2, W = sqrt(pi/2)*lfwhm

    real(dbl), intent(in) :: leff,lfwhm
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: q,dq
    real(dbl) :: f

    f = leff**2 / (c * sqrtpi2 * lfwhm)

    where( ph >= 0.0_dbl )
       q = f * ph
       dq = f * dph
    elsewhere
       q = -1.0_dbl
       dq = -1.0_dbl
    end where

  end subroutine ph2photnu


  subroutine ph2flux(leff,ph,dph,q,dq)

    ! (light intensity, energy) flux integrated over the passband [W/m2],
    ! in convenient electron-Volts [eV/s/m2]: ph * h*c/leff

    real(dbl), intent(in) :: leff
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: q,dq
    real(dbl) :: f

    f = hc / leff / evolt

    where( ph >= 0.0_dbl )
       q = f * ph
       dq = f * dph
    elsewhere
       q = -1.0_dbl
       dq = -1.0_dbl
    end where

  end subroutine ph2flux

    subroutine ph2flux0(leff,ph,dph,flux,dflux)

    ! (light intensity, energy) flux integrated over frequency band [W/m2]

    real(dbl), intent(in) :: leff
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: flux,dflux
    real(dbl) :: f

    f = planck * c / leff

    where( ph >= 0.0_dbl )
       flux = f * ph
       dflux = f * dph
    elsewhere
       flux = -1.0_dbl
       dflux = -1.0_dbl
    end where

  end subroutine ph2flux0

  subroutine cts2flux0(leff,cts,dcts,flux,dflux)
    ! replaced by ph2flux0

    ! (light intensity, energy) flux integrated over frequency band [W/m2]

    real(dbl), intent(in) :: leff
    real(dbl), dimension(:), intent(in) :: cts,dcts
    real(dbl), dimension(:), intent(out) :: flux,dflux
    real(dbl) :: f

    f = planck * c / leff

    where( cts >= 0.0_dbl )
       flux = f * cts
       dflux = f * dcts
    elsewhere
       flux = -1.0_dbl
       dflux = -1.0_dbl
    end where

  end subroutine cts2flux0

  subroutine ph2fllux(leff,ph,dph,flux,dflux)

    ! (light intensity, energy) flux integrated over wavelength band [W/m2]

    real(dbl), intent(in) :: leff
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: flux,dflux
    real(dbl) :: f

    f = planck * c / leff

    where( ph >= 0.0_dbl )
       flux = f * ph
       dflux = f * dph
    elsewhere
       flux = -1.0_dbl
       dflux = -1.0_dbl
    end where

  end subroutine ph2fllux


  subroutine ph2fnu(leff,lfwhm,ph,dph,q,dq)

    ! spectral flux density per frequency [W/m2/Hz],
    ! ph * (h*c/leff) / W_nu = ph * h * leff / W, W = sqrt(pi/2)*lfwhm

    real(dbl), intent(in) :: leff,lfwhm
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: q,dq
    real(dbl) :: f

    f = planck * leff / (sqrtpi2 * lfwhm)

    where( ph >= 0.0_dbl )
       q = f * ph
       dq = f * dph
    elsewhere
       q = -1.0_dbl
       dq = -1.0_dbl
    end where

  end subroutine ph2fnu

  subroutine ph2flam0(leff,lfwhm,ph,dph,q,dq)

    ! spectral flux density per wavelength [W/m2/nm],
    ! ph * (h*c/leff) / W, W = sqrt(pi/2)*lfwhm in nm

    real(dbl), intent(in) :: leff,lfwhm
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: q,dq
    real(dbl) :: f

    f = (hc / leff) / (sqrtpi2 * lfwhm / 1e-9_dbl)

    where( ph >= 0.0_dbl )
       q = f * ph
       dq = f * dph
    elsewhere
       q = -1.0_dbl
       dq = -1.0_dbl
    end where

  end subroutine ph2flam0

  subroutine mag2flux0(lwidth,zeroflam,mag,dmag,flux,dflux)

    ! magnitudes to flux in the passband [W/m2]

    real(dbl), dimension(:), intent(in) :: lwidth,zeroflam
    real(dbl), dimension(:), intent(in) :: mag,dmag
    real(dbl), dimension(:), intent(out) :: flux,dflux

    where( mag < 99 )
       flux = sqrtpi2 * 1e9 * zeroflam * lwidth * 10.0_dbl**(-0.4_dbl * mag) ! [nm->m]
    elsewhere
       flux = -1
    end where

    where( dmag < 9 )
       dflux = flux * dmag / magerr
    elsewhere
       dflux = 0
    end where

  end subroutine mag2flux0

  subroutine mag2ph(leff,lwidth,zeroflam,mag,dmag,ph,dph)

    ! magnitudes to photon flux in a passband [ph/s/m2]
    ! 3-th order aproximation

    real(dbl), dimension(:), intent(in) :: leff,lwidth,zeroflam
    real(dbl), dimension(:), intent(in) :: mag,dmag
    real(dbl), dimension(:), intent(out) :: ph,dph

!    real(dbl), dimension(size(mag)) :: flux,flam,e,mag0,n0,w,lw
    real(dbl), dimension(size(mag)) :: flux,flam,w,lw

    integer :: i

!    e = (/0.65,0.81,0.88,0.92,0.96/)
!    mag0 = 0
!    mag0 = (/-0.23611,-20.533,-21.226,-20.984,-0.35776/)
!    mag0 = (/-0.23611,-0.693,0.04,-20.984,-0.35776/)
!    mag0 = (/-0.23611,-0.77305+0.04,0.04,-0.24246+0.04,-0.35776/)
!    n0 = (/4.34413517E+09,1.34928722E+10,8.57578752E+09,1.31331963E+10,9.31594240E+09/)
    w = (/0.572981,0.683449,0.76545,0.645865,0.634657/)

!    flux = zeroflux * 10.0_dbl**(-0.4_dbl * (mag - mag0))

!    flux = w * flux
!    flux = e *flux
!    call flux2flam(leff,lwidth/2,flux,flam)
!    write(*,'(a,5es13.5)') 'flux:',flux
!    write(*,'(a,5es13.5)') 'flam',flam
!    call flam2plux(1e9*leff,1e9*lwidth/2,flam,ph)
!    write(*,'(a,5es13.5)') 'ph',ph

!    stop
!    ph = n0 * 10.0_dbl**(-0.4_dbl * (mag - mag0))

!    dph = dmag*ph/1.0857

    !
!    flux = e * zeroflux * 10.0_dbl**(-0.4_dbl * (mag - mag0))

!    write(*,*) real(zeroflam)
!    write(*,*) real(lwidth)
!    write(*,*) real(mag)

    flux = zeroflam * 10.0_dbl**(-0.4_dbl * mag) ! [nm -> m]
    flam = zeroflam * 10.0_dbl**(-0.4_dbl * mag) / lwidth ! [nm -> m]
!    flux = 1e9 * e * zeroflam * lwidth * 10.0_dbl**(-0.4_dbl * mag) ! [nm -> m]
 !    flux = 10.0_dbl**(-0.4_dbl * (mag - mag0))
!    flux = flux * w
!    ph = flux*leff

    if( size(mag) == 1 ) then

       ph = flux*leff

    else

       lw = lwidth / 2
       ph(1) = flux(1)*leff(1)*(1 + a1**2/2*lw(1)**2/leff(1)/(leff(2)-leff(1)) &
            *( flux(2)/flux(1) * lw(1)/lw(2) - 1))
!       write(*,*) ph(1)
!       ph(1) = flam(1)*leff(1)*(1 + a1**2/2*lw(1)**2/leff(1)/(leff(2)-leff(1)) &
!            *( flam(2)/flam(1) - 1))
!       write(*,*) ph(1)

!       write(*,*) flam(1),flam(2),flam(2)/flam(1)
!       write(*,*) flux(1),flux(2),lw(1),lw(2),flux(2)/flux(1) * lw(2)/lw(1)
!       stop

       do i = 2,size(mag)
!       ph(i) = 0.886227*flux(i)*leff(i)*(1 + &
!       ph(i) = flux(i)*leff(i)*(1 - &
!            0.707107**2/2*lw(i)**2/leff(i)/(leff(i+1)-leff(i)) + &
!            0.707107**2/2*lw(i)*lw(i+1)/leff(i)/(leff(i+1)-leff(i)) * flux(i+1)/flux(i))
          ph(i) = flux(i)*leff(i)*(1 + a1**2/2*lw(i)**2/leff(i)/(leff(i)-leff(i-1)) &
               *( flux(i-1)/flux(i) * lw(i)/lw(i-1) - 1))
!          ph(i) = flam(i)*leff(i)*(1 + a1**2/2*lw(i)**2/leff(i)/(leff(i)-leff(i-1)) &
!               *( flam(i-1)/flam(i) - 1))
!       write(*,*) 0.707107**2/2*lw(i)**2/leff(i)/(leff(i+1)-leff(i)), &
!            0.707107**2/2*lw(i)*lw(i+1)/leff(i)/(leff(i+1)-leff(i)) * flux(i+1)/flux(i)
!       write(*,*) flux(i)*leff(i)/hc,ph(i)/hc,real((flux(i)*leff(i)-ph(i))/ph(i))
       end do
!    do i = 1,size(mag)
!       ph(i) = flux(i)*leff(i)
!    end do

    end if

    ph = ph / hc
    dph = dmag*ph/magerr

!    write(*,'(a,5es13.5)') 'ph',ph


  end subroutine mag2ph

  subroutine mag2ph3(leff,lwidth,flamref,mag,dmag,ph,dph)

    ! magnitudes to photon flux in a passband [ph/s/m2]
    ! 3-th order aproximation

    real(dbl), dimension(:), intent(in) :: leff,lwidth,flamref
    real(dbl), dimension(:), intent(in) :: mag,dmag
    real(dbl), dimension(:), intent(out) :: ph,dph

!    real(dbl), dimension(size(mag)) :: flux,flam,e,mag0,n0,w,lw
    real(dbl), dimension(size(mag)) :: flux,flam,w

!    integer :: i

!    e = (/0.65,0.81,0.88,0.92,0.96/)
!    mag0 = 0
!    mag0 = (/-0.23611,-20.533,-21.226,-20.984,-0.35776/)
!    mag0 = (/-0.23611,-0.693,0.04,-20.984,-0.35776/)
!    mag0 = (/-0.23611,-0.77305+0.04,0.04,-0.24246+0.04,-0.35776/)
!    n0 = (/4.34413517E+09,1.34928722E+10,8.57578752E+09,1.31331963E+10,9.31594240E+09/)
    w = (/0.572981,0.683449,0.76545,0.645865,0.634657/)

    flux = flamref * 1e9 * 10.0_dbl**(-0.4_dbl * (mag - 0.04))

    flux(3) = 1e9 * flamref(3) * 10.0_dbl**(-0.4_dbl * mag(3)) * sqrtpi * lwidth(3)/2 * w(3)
    flux(2) = flux(3)*10**(-0.4*(mag(2)-mag(3)))*1.9644 !*  (w(3)/w(2))
    flux(1) = flux(2)*10**(-0.4*(mag(1)-mag(2)))*0.39476
    flux(4) = flux(3)*10**(-0.4*(mag(3)-mag(4)))/0.79986
    flux(5) = flux(4)*10**(-0.4*(mag(4)-mag(5)))/1.8034

!    write(*,'(a,5es13.5)') 'flux:',flux

!    flux = flamref * 1e9 * 10.0_dbl**(-0.4_dbl * mag) * sqrtpi * lwidth/2 / w
!    flux = flux / (w/w(3))
!    flux = flamref * 1e9 * 10.0_dbl**(-0.4_dbl * mag) * sqrtpi * lwidth/2 / w

!    flux = w * flux
!    flux = e *flux
    call flux2flam(leff,lwidth/2,flux,flam)
    flam = 1e-9*flam
    write(*,'(a,5es13.5)') 'flux:',flux
    write(*,'(a,5es13.5)') 'flam',flam
!    call flam2plux(1e9*leff,1e9*lwidth/2,flam,ph)


!    flam = 1e9 * flamref * 10.0_dbl**(-0.4_dbl * mag)
!    write(*,'(a,5es13.5)') 'flam',1e9 * flamref * 10.0_dbl**(-0.4_dbl * mag) / (lwidth/2) / sqrtpi

!    flam(3) = 1e9 * flamref(3) * 10.0_dbl**(-0.4_dbl * mag(3))
!    flam(2) = flam(3)*10**(-0.4*(mag(2)-mag(3)))
!    write(*,'(a,5es13.5)') 'flam',flam

    flam = 1e9*flam
!    flam = flux / lwidth
    call flam2plux(leff,lwidth/2,flam,ph)
    write(*,'(a,5es13.5)') 'ph',ph

    dph = ph * dmag / magerr

  end subroutine mag2ph3


  subroutine ph2mag3(leff,lwidth,zeroflux,ph,dph,mag,dmag)

    ! converts photon flux to magnitudes in a passband [ph/s/m2]
    ! 3-th order aproximation

    real(dbl), dimension(:), intent(in) :: leff,lwidth,zeroflux
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: mag,dmag

!    real(dbl), dimension(size(mag)) :: flux,flam,e,mag0,n0,w,lw
    real(dbl), dimension(size(mag)) :: flux,flam,w

!    integer :: i

!    e = (/0.65,0.81,0.88,0.92,0.96/)
!    mag0 = 0
!    mag0 = (/-0.23611,-20.533,-21.226,-20.984,-0.35776/)
!    mag0 = (/-0.23611,-0.693,0.04,-20.984,-0.35776/)
!    mag0 = (/-0.23611,-0.77305+0.04,0.04,-0.24246+0.04,-0.35776/)
!    n0 = (/4.34413517E+09,1.34928722E+10,8.57578752E+09,1.31331963E+10,9.31594240E+09/)
    w = (/0.572981,0.683449,0.76545,0.645865,0.634657/)

!    write(*,*)
!    flux = zeroflux * 10.0_dbl**(-0.4_dbl * mag)

!    flux = w * flux
!    flux = e *flux
    call ph2flam(leff,lwidth/2,ph,flam)
!    write(*,'(a,5es13.5)') 'ph:',ph
!    write(*,'(a,5es13.5)') 'flam',flam
    call flam2flux(leff,lwidth/2,flam,flux)
!    write(*,'(a,5es13.5)') 'flux',flux
    mag = -2.5*log10(flux/zeroflux)
    dmag = (dph/ph) * mag * magerr !?

  end subroutine ph2mag3

!!$  subroutine mag2pht(leff,lwidth,zeroflux,exptime,area,mag,dmag,ph,dph)
!!$
!!$    ! magnitudes to photon flux in a passband on area and time [ph]
!!$
!!$    real(dbl), dimension(:), intent(in) :: leff,lwidth,zeroflux,exptime,area
!!$    real(dbl), dimension(:), intent(in) :: mag,dmag
!!$    real(dbl), dimension(:), intent(out) :: ph,dph
!!$    real(dbl), dimension(size(mag)) :: c
!!$
!!$    call mag2ph(leff,lwidth,zeroflux,mag,dmag,ph,dph)
!!$
!!$    c = exptime * area
!!$    ph = c * ph
!!$    dph = c *dph
!!$
!!$  end subroutine mag2pht

  subroutine flux2flam(leff,lw,flux,flam)

    ! flam deconvolution from flux

    real(dbl), dimension(:), intent(in) :: leff,lw,flux
    real(dbl), dimension(:), intent(out) :: flam

    real(dbl), dimension(size(flux)) :: w
    real(dbl), dimension(size(flux)-2) :: b,d,s,p,x
    real(dbl) :: hh
    integer :: i,j,n,m,iter

    n = size(flux)
    m = n - 2

    w = 1
    w = (/0.572981,0.683449,0.76545,0.645865,0.634657/)

    ! initial estimate of flux density (used 1 and n-th points)
    flam = w*flux/lw/sqrtpi

    do iter = 1,2

       b = w(2:n-1) * flux(2:n-1)
       p(1) = 0.0_dbl
       s(m) = 0.0_dbl

       do i = 2,n-1
          j = i - 1
          d(j) = (h0+2*h12*(1+(lw(i)*a)**2/(leff(i)-leff(i-1))/(leff(i)-leff(i+1))))*lw(i)
          hh = 2*h12*lw(i)**3*a**2/(leff(i+1)-leff(i-1))
          if( j > 1 ) p(j) = hh/(leff(i) - leff(i-1))
          if( j < m ) s(j) = hh/(leff(i+1) - leff(i))
          if( j == 1 ) b(j) = b(j) - hh*flam(i-1)/(leff(i) - leff(i-1))
          if( j == m ) b(j) = b(j) - hh*flam(i+1)/(leff(i+1) - leff(i))
       end do

       x = b
       call tridig(d,p,s,x)
       flam(2:n-1) = x

       ! more accurate estimate of points on edges
       hh = h12*lw(1)**3*a**2/(leff(3)-leff(1))
       p(1) = hh*(flam(3) - flam(2))/(leff(3) - leff(2))
       s(1) = hh*flam(2)/(leff(2) - leff(1))
       flam(1) = (w(1)*flux(1) - p(1) + s(1))/((h0 + 2*h12)*lw(1)+hh/(leff(2) - leff(1)))

       hh = h12*lw(n)**3*a**2/(leff(n-2)-leff(n))
       p(1) = hh*(flam(n-2) - flam(n-1))/(leff(n-2) - leff(n-1))
       s(1) = hh*flam(n-1)/(leff(n-1) - leff(n))
       flam(n) = (w(n)*flux(n) - p(1) + s(1))/((h0 + 2*h12)*lw(n)+hh/(leff(n-1)-leff(n)))

    end do


!!$    ! initial estimate of flux density (used 1 and n-th points)
!!$    w = 1
!!$    w = (/0.572981,0.683449,0.76545,0.645865,0.634657/)
!!$!    flam = w*flux/lw/sqrtpi*1e-9
!!$    write(*,*) flux
!!$    flam = w * flux / lw /sqrtpi * 1e-9
!!$    write(*,*) real(flam)
!!$
!!$!    ph = zeroflux * leff * 10.0_dbl**(-0.4_dbl*mag) / hc
!!$
!!$    do iter = 1,2
!!$
!!$       b = w(2:n-1) * flux(2:n-1)
!!$       p(1) = 0.0_dbl
!!$       s(m) = 0.0_dbl
!!$
!!$       do i = 2,n-1
!!$          j = i - 1
!!$          d(j) = (h0+2*h12*(1-(lw(i)*a)**2/(leff(i)-leff(i-1))/(leff(i+1)-leff(i))))*lw(i)
!!$          hh = 2*h12*lw(i)**3*a**2/(leff(i+1)-leff(i-1))
!!$          if( j > 1 ) p(j) = hh/(leff(i) - leff(i-1))
!!$          if( j < m ) s(j) = hh/(leff(i+1) - leff(i))
!!$          if( j == 1 ) b(j) = b(j) - hh*flam(i-1)/(leff(i) - leff(i-1))
!!$          if( j == m ) b(j) = b(j) - hh*flam(i+1)/(leff(i+1) - leff(i))
!!$       end do
!!$
!!$       x = b
!!$       call tridig(d,p,s,x)
!!$       flam(2:n-1) = x
!!$
!!$       ! more accurate estimate of points on edges
!!$       hh = h12*lw(1)**3*a**2/(leff(3)-leff(1))
!!$       p(1) = hh*(flam(3) - flam(2))/(leff(3) - leff(2))
!!$       s(1) = hh*flam(2)/(leff(2) - leff(1))
!!$       flam(1) = (w(1)*flux(1) - p(1) + s(1))/((h0 + 2*h12)*lw(1)+hh/(leff(2) - leff(1)))
!!$
!!$       hh = h12*lw(n)**3*a**2/(leff(n-2)-leff(n))
!!$       p(1) = hh*(flam(n-2) - flam(n-1))/(leff(n-2) - leff(n-1))
!!$       s(1) = hh*flam(n-1)/(leff(n-1) - leff(n))
!!$       flam(n) = (w(n)*flux(n) - p(1) + s(1))/((h0 + 2*h12)*lw(n)+hh/(leff(n-1)-leff(n)))
!!$
!!$!       write(*,*) real(flam)
!!$
!!$    end do
!!$
!!$!    flam = w*flux/lw/sqrtpi
!!$!    flam = 1d-9 * flam ! conversion [W/m3] -> [W/m2/nm]
!!$
!!$    write(*,*) real(flam)

!    stop

!!$
!!$
!!$    w = (/0.572981,0.683449,0.76545,0.645865,0.634657/)
!!$
!!$    ! initial estimate of flux density (used 1 and n-th points)
!!$    flam = w*flux/lw/sqrtpi
!!$
!!$    do iter = 1,2
!!$
!!$       b = w(2:n-1) * flux(2:n-1)
!!$       p(1) = 0
!!$       s(m) = 0
!!$
!!$       do i = 2,n-1
!!$          j = i - 1
!!$          d(j) = (h0+2*h12*(1+(lw(i)*a)**2/(leff(i)-leff(i-1))/(leff(i)-leff(i+1))))*lw(i)
!!$          hh = 2*h12*lw(i)**3*a**2/(leff(i+1)-leff(i-1))
!!$          if( j > 1 ) p(j) = hh/(leff(i) - leff(i-1))
!!$          if( j < m ) s(j) = hh/(leff(i+1) - leff(i))
!!$          if( j == 1 ) b(j) = b(j) - hh*flam(i-1)/(leff(i) - leff(i-1))
!!$          if( j == m ) b(j) = b(j) - hh*flam(i+1)/(leff(i+1) - leff(i))
!!$       end do
!!$
!!$       x = b
!!$       call tridig(d,p,s,x)
!!$       flam(2:n-1) = x
!!$
!!$       ! more accurate estimate of points on edges
!!$       hh = h12*lw(1)**3*a**2/(leff(3)-leff(1))
!!$       p(1) = hh*(flam(3) - flam(2))/(leff(3) - leff(2))
!!$       s(1) = hh*flam(2)/(leff(2) - leff(1))
!!$       flam(1) = (w(1)*flux(1) - p(1) + s(1))/((h0 + 2*h12)*lw(1)+hh/(leff(2) - leff(1)))
!!$
!!$       hh = h12*lw(n)**3*a**2/(leff(n-2)-leff(n))
!!$       p(1) = hh*(flam(n-2) - flam(n-1))/(leff(n-2) - leff(n-1))
!!$       s(1) = hh*flam(n-1)/(leff(n-1) - leff(n))
!!$       flam(n) = (w(n)*flux(n) - p(1) + s(1))/((h0 + 2*h12)*lw(n)+hh/(leff(n-1)-leff(n)))
!!$
!!$    end do


  end subroutine flux2flam

  subroutine flam2flux(leff,lw,flam,flux)

    ! flam deconvolution from flux

    real(dbl), dimension(:), intent(in) :: leff,lw,flam
    real(dbl), dimension(:), intent(out) :: flux

    real(dbl), dimension(size(flam)) :: w
    real(dbl) :: hh,ww,dd,u1,u2,df,df1,t
    integer :: i,n

    w = (/0.572981,0.683449,0.76545,0.645865,0.634657/)
    w = 1

    n = size(flux)

    do i = 2,n-1
       ww = h0 + 2*h12*(1 + (lw(i)*a)**2/(leff(i)-leff(i-1))/(leff(i)-leff(i+1)))
       dd = 2*h12*lw(i)**3*a**2
       u1 = dd/(leff(i-1)-leff(i))/(leff(i-1)-leff(i+1))
       u2 = dd/(leff(i+1)-leff(i))/(leff(i+1)-leff(i-1))
       flux(i) = ww*lw(i)*flam(i) + u1*flam(i-1) + u2*flam(i+1)
    end do

    hh = h0 + 2*h12
    i = 1
    t = (lw(i)*a)**2
    df = (flam(i+1) - flam(i)) / (leff(i+1) - leff(i))
    df1 = (flam(i+2) - flam(i+1)) / (leff(i+2) - leff(i+1))
    flux(i) = hh*flam(i)*lw(i) - h12*t*lw(i)/(leff(i+2) - leff(i))*(df - df1)

    i = n
    t = (lw(i)*a)**2
    df = (flam(i-1) - flam(i)) / (leff(i-1) - leff(i))
    df1 = (flam(i-2) - flam(i-1)) / (leff(i-2) - leff(i-1))
    flux(i) = hh*flam(i)*lw(i) - h12*t*lw(i)/(leff(i-2) - leff(i))*(df - df1)

    flux = flux / w

!!$    do i = 2,n-1
!!$       ww = h0 + 2*h12*(1 + (lw(i)*a)**2/(leff(i)-leff(i-1))/(leff(i)-leff(i+1)))
!!$       dd = 2*h12*lw(i)**3*a**2
!!$       u1 = dd/(leff(i-1)-leff(i))/(leff(i-1)-leff(i+1))
!!$       u2 = dd/(leff(i+1)-leff(i))/(leff(i+1)-leff(i-1))
!!$       flux(i) = ww*lw(i)*flam(i) + u1*flam(i-1) + u2*flam(i+1)
!!$    end do
!!$
!!$    hh = h0 + 2*h12
!!$    i = 1
!!$    t = (lw(i)*a)**2
!!$    df = (flam(i+1) - flam(i)) / (leff(i+1) - leff(i))
!!$    df1 = (flam(i+2) - flam(i+1)) / (leff(i+2) - leff(i+1))
!!$    flux(i) = hh*flam(i)*lw(i) + h12*t*lw(i)/(leff(i+2) - leff(i))*(df - df1)
!!$
!!$    i = n
!!$    t = (lw(i)*a)**2
!!$    df = (flam(i-1) - flam(i)) / (leff(i-1) - leff(i))
!!$    df1 = (flam(i-2) - flam(i-1)) / (leff(i-2) - leff(i-1))
!!$    flux(i) = hh*flam(i)*lw(i) + h12*t*lw(i)/(leff(i-2) - leff(i))*(df - df1)
!!$
!!$    flux = flux / w


  end subroutine flam2flux

  subroutine ph2flam(leff,lw,ph,flam)

    ! flam deconvolution from photons

    real(dbl), dimension(:), intent(in) :: leff,lw,ph
    real(dbl), dimension(:), intent(out) :: flam

    real(dbl), dimension(size(ph)) :: w,wplux
    real(dbl), dimension(size(ph)-2) :: b,d,s,p,x
    real(dbl) :: hh,ww,dd,u,df1,t,t2,uu,u0,u1,u2
    integer :: i,j,n,m,iter

    w = (/0.572981,0.683449,0.76545,0.645865,0.634657/)

    w = 1
    n = size(ph)
    m = n - 2


    wplux = w * ph

    ! initial estimate of flux density (used 1-th and n-th points)
    flam = hc*wplux/lw/sqrtpi/leff

!    write(*,*) real(flam)

    do iter = 1,2

       b = wplux(2:n-1)
       p(1) = 0
       s(m) = 0

       do i = 2,n-1

          j = i - 1

          t2 = (lw(i)*a)**2
          u1 = 2*t2*(2*leff(i) - leff(i+1))/(leff(i-1)-leff(i))/(leff(i-1)-leff(i+1))
          u2 = 2*t2*(2*leff(i) - leff(i-1))/(leff(i+1)-leff(i-1))/(leff(i+1)-leff(i))
          uu = leff(i)**2 + leff(i+1)*leff(i-1) - leff(i)*(leff(i+1)+leff(i-1)) + t2
          u0 = 2*leff(i)*uu - 2*t2*(2*leff(i)-(leff(i+1)+leff(i-1)))
          u0 = u0 / (leff(i)-leff(i-1))/(leff(i)-leff(i+1))

!          write(*,*) u0,u1,u2

          uu = lw(i) / hc
          d(j) = uu * (h0*leff(i) + h12*u0)
          if( j > 1 ) p(j) = h12 * u1 * uu
          if( j < m ) s(j) = h12 * u2 * uu
          if( j == 1 ) b(j) = b(j) - flam(i-1)*u1* uu * h12
          if( j == m ) b(j) = b(j) - flam(i+1)*u2* uu * h12

!          write(*,*) d(j),p(j),s(j)

       end do

!       write(*,*) real(d)
!       write(*,*) real(p)
!       write(*,*) real(s)
!       write(*,*) real(b)

       x = b
       call tridig(d,p,s,x)
       flam(2:n-1) = x

!       write(*,*) real(flam)

       hh = h0 + 2*h12
       i = 1
       t = (a*lw(i))**2
       df1 = (flam(i+2) - flam(i+1)) / (leff(i+2) - leff(i+1))
       u = 1 - leff(i)/(leff(i+2) - leff(i))
       ww = hh*leff(i) - h12*t*u/(leff(i+1) - leff(i))
       dd = h12*t*u/(leff(i+1) - leff(i))*flam(i+1) + &
            h12*t*leff(i)/(leff(i+2) - leff(i))*df1
       flam(i) = (wplux(i)*hc/lw(i) - dd)/ww

       i = n
       t = (a*lw(i))**2
       df1 = (flam(i-2) - flam(i-1)) / (leff(i-2) - leff(i-1))
       u = 1 - leff(i)/(leff(i-2) - leff(i))
       ww = hh*leff(i) - h12*t*u/(leff(i-1) - leff(i))
       dd = h12*t*u/(leff(i-1) - leff(i))*flam(i-1) + &
            h12*t*leff(i)/(leff(i-2) - leff(i))*df1
       flam(i) = (wplux(i)*hc/lw(i) - dd)/ww

!       write(*,*) real(flam)

!       stop

    end do



!!$    ! initial estimate of flux density, flam (used 1-th and n-th points)
!!$    flam = hc*ph/lw/sqrtpi/leff * w
!!$!    flam(2:n-1) = hc*plux(2:n-1)/lw(2:n-1)/sqrtpi/leff(2:n-1) * w(2:n-1)
!!$!    flam(1) = 1.61426515E-12
!!$!    flam(n) = 7.64483623E-14
!!$
!!$    write(*,*) real(leff)
!!$    write(*,*) real(lw)
!!$    write(*,*) real(flam)
!!$
!!$    do iter = 1,2
!!$
!!$       b = w(2:n-1) * ph(2:n-1) !?
!!$!       b = ph(2:n-1)
!!$       p(1) = 0
!!$       s(m) = 0
!!$
!!$       do i = 2,n-1
!!$
!!$          j = i - 1
!!$          t = lw(i)*a
!!$          u = (leff(i) - leff(i-1))*(leff(i+1) - leff(i))
!!$          write(*,*) u,leff(i) - leff(i-1),leff(i+1) - leff(i)
!!$          d(j) = h0*leff(i) + 2*h12*leff(i)*(1-t**2/u) &
!!$               + 2*h12*t*(leff(i+1) + leff(i-1) - 2*leff(i))/u
!!$          d(j) = d(j)*lw(i) / hc
!!$
!!$          write(*,*) h0*leff(i),h12*leff(i)*(1-t**2/u),2*h12*t*(leff(i+1) + leff(i-1) - 2*leff(i))/u,u
!!$
!!$          hh = 2*h12*lw(i)**3*a**2/(leff(i+1)-leff(i-1))
!!$          u1 = hh*(2*leff(i)-leff(i+1))/(leff(i) - leff(i-1))
!!$          u2 = hh*(2*leff(i)-leff(i-1))/(leff(i+1) - leff(i))
!!$          if( j > 1 ) p(j) = u1 / hc
!!$          if( j < m ) s(j) = u2 / hc
!!$          if( j == 1 ) b(j) = b(j) - flam(i-1)*u1 / hc / w(i-1)
!!$          if( j == m ) b(j) = b(j) - flam(i+1)*u2 / hc / w(i+1)
!!$
!!$          write(*,*) d(j),p(j),s(j)
!!$
!!$       end do
!!$
!!$!       write(*,*) real(d)
!!$!       write(*,*) real(p)
!!$!       write(*,*) real(s)
!!$!       write(*,*) real(b)
!!$
!!$       x = b
!!$       call tridig(d,p,s,x)
!!$       flam(2:n-1) = x
!!$
!!$       write(*,*) real(flam)
!!$
!!$       hh = h0 + 2*h12
!!$       i = 1
!!$       t = (a*lw(i))**2
!!$       df1 = (flam(i+2) - flam(i+1)) / (leff(i+2) - leff(i+1))
!!$       u = 1 - leff(i)/(leff(i+2) - leff(i))
!!$       ww = hh*leff(i) - h12*t*u/(leff(i+1) - leff(i))
!!$       dd = h12*t*u/(leff(i+1) - leff(i))*flam(i+1) + &
!!$            h12*t*leff(i)/(leff(i+2) - leff(i))*df1
!!$       flam(i) = (ph(i)*hc/lw(i) - dd)/ww
!!$
!!$       i = n
!!$       t = (a*lw(i))**2
!!$       df1 = (flam(i-2) - flam(i-1)) / (leff(i-2) - leff(i-1))
!!$       u = 1 - leff(i)/(leff(i-2) - leff(i))
!!$       ww = hh*leff(i) - h12*t*u/(leff(i-1) - leff(i))
!!$       dd = h12*t*u/(leff(i-1) - leff(i))*flam(i-1) + &
!!$            h12*t*leff(i)/(leff(i-2) - leff(i))*df1
!!$       flam(i) = (ph(i)*hc/lw(i) - dd)/ww
!!$
!!$       write(*,*) real(flam)
!!$
!!$!       stop
!!$
!!$    end do

  end subroutine ph2flam


  subroutine flam2plux(leff,lw,flam,plux3)

    ! photon flux from flam

    real(dbl), dimension(:), intent(in) :: leff,lw,flam
    real(dbl), dimension(:), intent(out) :: plux3

    real(dbl), dimension(size(flam)) :: w

    real(dbl) :: hh,df,df1,t,t2,uu,u0,u1,u2
    integer :: i,n

    w = 1
    w = (/0.572981,0.683449,0.76545,0.645865,0.634657/)
    n = size(flam)

    do i = 2,n-1

       t2 = (lw(i)*a)**2
       u1 = 2*t2*(2*leff(i) - leff(i+1))/(leff(i-1)-leff(i))/(leff(i-1)-leff(i+1))
       u2 = 2*t2*(2*leff(i) - leff(i-1))/(leff(i+1)-leff(i-1))/(leff(i+1)-leff(i))
       uu = leff(i)**2 + leff(i+1)*leff(i-1) - leff(i)*(leff(i+1)+leff(i-1)) + t2
       u0 = 2*leff(i)*uu - 2*t2*(2*leff(i)-(leff(i+1)+leff(i-1)))
       u0 = u0 / (leff(i)-leff(i-1))/(leff(i)-leff(i+1))

       plux3(i) = lw(i)*((h0*leff(i) + h12*u0)*flam(i)+h12*(u1*flam(i-1)+u2*flam(i+1)))

    end do

    i = 1
    hh = h0 + 2*h12
    t = (a*lw(i))**2
    df = (flam(i+1) - flam(i)) / (leff(i+1) - leff(i))
    df1 = (flam(i+2) - flam(i)) / (leff(i+2) - leff(i+1))
    plux3(i) = lw(i)*(hh*leff(i)*flam(i) + h12*df*t + &
         h12*t*leff(i)/(leff(i+2) - leff(i))*(df1 - df))

    i = n
    t = (a*lw(i))**2
    df = (flam(i-1) - flam(i)) / (leff(i-1) - leff(i))
    df1 = (flam(i-2) - flam(i)) / (leff(i-2) - leff(i-1))
    plux3(i) = lw(i)*(hh*leff(i)*flam(i) + h12*df*t + &
         h12*t*leff(i)/(leff(i-2) - leff(i))*(df1 - df))

    plux3 = plux3 / hc / w


!!$
!!$    do i = 2,n-1
!!$       plux3(i) = (h0*leff(i) &
!!$            + 2*h12*(1-(lw(i)*a)**2/(leff(i) - leff(i-1))/(leff(i+1) - leff(i)))*leff(i) &
!!$            + 2*h12*lw(i)*(leff(i+1) + leff(i-1) - 2*leff(i))/(leff(i) - leff(i-1))/(leff(i+1) - leff(i)) )*flam(i)*lw(i) + &
!!$            2*h12*lw(i)**3*a**2/(leff(i+1)-leff(i-1))* &
!!$            (flam(i-1)*(2*leff(i) - leff(i+1))/(leff(i) - leff(i-1)) + flam(i+1)*(2*leff(i) - leff(i-1))/(leff(i+1) - leff(i)))
!!$    end do
!!$
!!$    dd = lw(1)*a
!!$    df = (flam(2) - flam(1)) / (leff(2) - leff(1))
!!$    df1 = (flam(3) - flam(2)) / (leff(3) - leff(2))
!!$    hh = h12*lw(1)**3*a**2/(leff(3)-leff(1))
!!$    plux3(1) = lw(1)*(h0*flam(1)*leff(1) + &
!!$         h12*(leff(1)*(2*flam(1) + (df1 - df)*dd**2/(leff(3)-leff(1))) +&
!!$         2*df*dd**2))
!!$
!!$    dd = lw(n)*a
!!$    df = (flam(n-1) - flam(n)) / (leff(n-1) - leff(n))
!!$    df1 = (flam(n-2) - flam(n-1)) / (leff(n-2) - leff(n-1))
!!$    hh = h12*lw(n)**3*a**2/(leff(n-2)-leff(n))
!!$    plux3(n) = lw(n)*(h0*flam(n)*leff(n) + &
!!$         h12*(leff(n)*(2*flam(n) + (df1 - df)*dd**2/(leff(n-2)-leff(n))) +&
!!$         2*df*dd**2))
!!$
!!$    plux3 = plux3 / hc / w * 1d-18 ! due flam [W/m2/nm]
!!$


!!$    do i = 2,n-1
!!$
!!$       plux3(i) = (h0*leff(i) &
!!$            + 2*h12*(1-(lw(i)*a)**2/(leff(i) - leff(i-1))/(leff(i+1) - leff(i)))*leff(i) &
!!$            + 2*h12*lw(i)*a*(leff(i+1) + leff(i-1) - 2*leff(i))/(leff(i) - leff(i-1))/(leff(i+1) - leff(i)) )*flam(i)*lw(i) + &
!!$            2*h12*lw(i)**3*a**2/(leff(i+1)-leff(i-1))* &
!!$            (flam(i-1)*(2*leff(i) - leff(i+1))/(leff(i) - leff(i-1)) + flam(i+1)*(2*leff(i) - leff(i-1))/(leff(i+1) - leff(i)))
!!$    end do
!!$
!!$    i = 1
!!$    hh = h0 + 2*h12
!!$    t = (a*lw(i))**2
!!$    df = (flam(i+1) - flam(i)) / (leff(i+1) - leff(i))
!!$    df1 = (flam(i+2) - flam(i)) / (leff(i+2) - leff(i+1))
!!$    plux3(i) = lw(i)*(hh*leff(i)*flam(i) + h12*df*t + &
!!$         h12*t*leff(i)/(leff(i+2) - leff(i))*(df1 - df))
!!$
!!$    i = n
!!$    t = (a*lw(i))**2
!!$    df = (flam(i-1) - flam(i)) / (leff(i-1) - leff(i))
!!$    df1 = (flam(i-2) - flam(2)) / (leff(i-2) - leff(i-1))
!!$    plux3(i) = lw(i)*(hh*leff(i)*flam(i) + h12*df*t + &
!!$         h12*t*leff(i)/(leff(i-2) - leff(i))*(df1 - df))
!!$    plux3 = plux3 / hc / w
!!$

  end subroutine flam2plux


  subroutine mag2ph0(leff,lfwhm,flamref,mag,dmag,ph,dph)

    ! magnitudes to photon flux in a passband [ph/s/m2]
    ! 0-th order Gauss-Hermite integration

    real(dbl), intent(in) :: leff,lfwhm,flamref
    real(dbl), dimension(:), intent(in) :: mag,dmag
    real(dbl), dimension(:), intent(out) :: ph,dph
    real(dbl) :: ph0

    ph0 = 1e9 * sqrtpi2 * flamref * lfwhm *  leff / hc

    where( mag < 99 )
       ph = ph0 * 10.0_dbl**(-0.4_dbl*mag)
    elsewhere
       ph = -1
    end where

    where( dmag < 9 .and. ph > 0 )
       dph = ph * dmag / magerr
!    elsewhere( dmag >= 9 .and. mag < 99 )
!       dph = sqrt(ph)
    elsewhere
       dph = -1
    end where

  end subroutine mag2ph0

  subroutine ph2mag0(leff,lfwhm,flamref,ph,dph,mag,dmag)

    ! photons [ph/s/m2] to flux [W/m2] (and to magnitudes)

    real(dbl), intent(in) :: flamref,leff,lfwhm
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: mag,dmag
    real(dbl) :: ph0

    ph0 = 1e9 * sqrtpi2 * flamref * lfwhm * leff / hc

    where( ph > 0.0_dbl )
       mag = -2.5_dbl*log10(ph / ph0)
    elsewhere
       mag = 99.99999_dbl
    end where

    where( dph > 0 .and. ph > 0.0_dbl )
       dmag = magerr*dph/ph
    elsewhere
       dmag = 9.99999_dbl
    end where

  end subroutine ph2mag0

  subroutine flux2mag0(lfwhm,flamref,flux,dflux,mag,dmag)

    ! flux [W/m2] to magnitudes

    real(dbl), intent(in) :: flamref,lfwhm
    real(dbl), dimension(:), intent(in) :: flux,dflux
    real(dbl), dimension(:), intent(out) :: mag,dmag
    real(dbl) :: f0

    f0 = sqrtpi2 * lfwhm * flamref * 1e9

    where( flux > 0.0_dbl )
       mag = -2.5_dbl*log10(flux / f0)
       dmag = magerr*dflux/flux
    elsewhere
       mag = 99.99999_dbl
       dmag = 9.99999_dbl
    end where

  end subroutine

!!$  subroutine mag2pht0(leff,lwidth,zeroflux,exptime,area,mag,dmag,ph,dph)
!!$
!!$    ! magnitudes to photon flux in a passband on area and time [ph]
!!$
!!$    real(dbl), dimension(:), intent(in) :: leff,lwidth,zeroflux,exptime,area
!!$    real(dbl), dimension(:), intent(in) :: mag,dmag
!!$    real(dbl), dimension(:), intent(out) :: ph,dph
!!$    real(dbl), dimension(size(mag)) :: c
!!$
!!$    call mag2ph(leff,lwidth,zeroflux,mag,dmag,ph,dph)
!!$
!!$    c = exptime * area
!!$    ph = c * ph
!!$    dph = c *dph
!!$
!!$  end subroutine mag2pht0


  subroutine ph2mag(leff,lfwhm,flamref,ph,dph,mag,dmag)

    ! photons [ph/s/m2] to flux [W/m2] (and to magnitudes)

    real(dbl), dimension(:), intent(in) :: flamref,leff,lfwhm
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: mag,dmag

    real(dbl), dimension(size(mag)) :: flux,lw
    real(dbl) :: u,v,det
    integer :: i

    if( size(ph) == 1 ) then

       flux = hc / leff * ph
       mag = -2.5_dbl * log10(flux / flamref)

    else

       lw = lfwhm / 2

       u = a**2/2/(leff(2) - leff(1))
       det = 1 - u*(lw(2)**2/leff(2) + lw(1)**2/leff(1))
       det = 1
       v = 1 + u*lw(2)**2/leff(2)*(1-(leff(1)/leff(2))**3*(ph(2)/ph(1)))
       flux(1) = ph(1)*v/leff(1)/det

       do i = 2,size(ph)

          u = a**2/2/(leff(i) - leff(i-1))
          det = 1 - u*(lw(i)**2/leff(i) + lw(i-1)**2/leff(i-1))
          det = 1
          v = 1 + u*lw(i)**2/leff(i)*(1-(leff(i-1)/leff(i))**3*(ph(i-1)/ph(1)))
          flux(i) = ph(i)*v/leff(i)/det

       end do

       flux = hc * flux
       mag = -2.5_dbl * log10(flux / flamref)
       dmag = magerr * (dph/ph) * mag

    end if

  end subroutine ph2mag

  subroutine mag2rate(mag,dmag,r,dr)

    ! relative rates from magnitudes

    real(dbl), dimension(:), intent(in) :: mag,dmag
    real(dbl), dimension(:), intent(out) :: r,dr

    r = flux25*10.0_dbl**(-0.4_dbl*mag)
    dr = dmag*r/magerr

  end subroutine mag2rate

  subroutine x2mag(f0,flux,dflux,mag,dmag)

    ! (relative) flux (or anythink) to magnitudes

    real(dbl), intent(in) :: f0
    real(dbl), dimension(:), intent(in) :: flux,dflux
    real(dbl), dimension(:), intent(out) :: mag,dmag

    where( flux > 0.0_dbl .and. f0 > 0.0_dbl )
       mag = -2.5_dbl*log10(flux/f0)
       dmag = magerr*dflux/flux
    elsewhere
       mag = 99.99999_dbl
       dmag = 9.99999_dbl
    end where

  end subroutine x2mag

  subroutine fnu2abmag(fnu,dfnu,mag,dmag)

    ! fnu to AB-magnitudes

    real(dbl), dimension(:), intent(in) :: fnu,dfnu
    real(dbl), dimension(:), intent(out) :: mag,dmag

    where( fnu > 0.0_dbl )
       mag = -2.5_dbl*log10(fnu/ABspflux)
       dmag = magerr*dfnu/fnu
    elsewhere
       mag = 99.99999_dbl
       dmag = 9.99999_dbl
    end where

  end subroutine fnu2abmag

  subroutine flam2stmag(flam,dflam,mag,dmag)

    ! flam to ST-magnitudes

    real(dbl), dimension(:), intent(in) :: flam,dflam
    real(dbl), dimension(:), intent(out) :: mag,dmag

    where( flam > 0.0_dbl )
       mag = -2.5_dbl*log10(flam/STspflux)
       dmag = magerr*dflam/flam
    elsewhere
       mag = 99.99999_dbl
       dmag = 9.99999_dbl
    end where

  end subroutine flam2stmag

  subroutine ph2abmag(leff,lfwhm,ph,dph,mag,dmag)

    ! photons to AB-magnitudes, -2.5*log10(FNU / ABspflux), FNU as in ph2fnu

    real(dbl), intent(in) :: leff,lfwhm
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: mag,dmag
    real(dbl) :: f

    f = planck * leff / (sqrtpi2 * lfwhm) / ABspflux

    where( ph > 0.0_dbl )
       mag = -2.5_dbl*log10(f*ph)
       dmag = magerr*dph/ph
    elsewhere
       mag = 99.99999_dbl
       dmag = 9.99999_dbl
    end where

  end subroutine ph2abmag

  subroutine ph2stmag(leff,lfwhm,ph,dph,mag,dmag)

    ! photons to ST-magnitudes, -2.5*log10(FLAM / STspflux), FLAM as in ph2flam0

    real(dbl), intent(in) :: leff,lfwhm
    real(dbl), dimension(:), intent(in) :: ph,dph
    real(dbl), dimension(:), intent(out) :: mag,dmag
    real(dbl) :: f

    f = (hc / leff) / (sqrtpi2 * lfwhm / 1e-9_dbl) / STspflux

    where( ph > 0.0_dbl )
       mag = -2.5_dbl*log10(f*ph)
       dmag = magerr*dph/ph
    elsewhere
       mag = 99.99999_dbl
       dmag = 9.99999_dbl
    end where

  end subroutine ph2stmag

  subroutine hstparam(lam_eff,lam_fwhm,photflam,photzpt,photplam,photbw)

    ! HST-like keywords: PHOTFLAM is the spectral flux density [erg/s/cm2/A]
    ! corresponding to the photon rate 1 photon/s/cm2 in the passband,
    ! (h*c/lam_eff) / W, W = sqrt(pi/2)*lam_fwhm; PHOTZPT is the ST
    ! zero point, STMAG = -2.5*log10(FLAM) + PHOTZPT for FLAM in cgs.

    real(dbl), intent(in) :: lam_eff, lam_fwhm
    real(dbl), intent(out) :: photflam,photzpt,photplam,photbw

    ! [J/photon] -> [erg] (1e7), per [A] (W / 1e-10 m)
    photflam = 1e7_dbl * (hc / lam_eff) / (sqrtpi2 * lam_fwhm / 1e-10_dbl)
    ! [W/m2/nm] -> [erg/s/cm2/A] is factor 100
    photzpt = 2.5_dbl * log10(100 * STspflux)
    photplam = lam_eff / 1e-10
    photbw = lam_fwhm / 1e-10

  end subroutine hstparam


  ! solution of tri-diagonal system of linear equations
  subroutine tridig(d,p,s,b)

    real(dbl), dimension(:), intent(in out) :: d,p,s,b
    real(dbl) :: r
    integer :: i,k,n

    n = size(d)

    do i = 2,n
       r = p(i) / d(i-1)
       d(i) = d(i) - r* s(i-1)
       b(i) = b(i) - r* b(i-1)
    end do

    b(n) = b(n) / d(n)
    do i = 1,n-1
       k = n - i
       b(k) = (b(k) -s(k)*b(k+1)) / d(k)
    end do

  end subroutine tridig

end module photoconv
