!
!  Spherical astronomy module
!
!  Copyright © 1996 - 2026 F.Hroch (hroch@physics.muni.cz)
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



module astrosphere

  use iso_fortran_env
  implicit none

  real(REAL64), parameter, private :: rad =  57.29577951308232286464772_REAL64

contains

  pure function gmst(jd)

    real(REAL64) :: gmst
    real(REAL64), intent(in) :: jd


!  Greenwich sidereal time in hours
!
!  jd is a full Julian date
!
!  The precision is better than 1 second.
!  According to Astronomical Almanac 2000.


    real(REAL64) :: tu,t

    tu = (jd - 2451545.0_REAL64)/36525.0_REAL64
    t = 24110.54841_REAL64 + tu*(8640184.812866_REAL64 + tu*(0.093104_REAL64-6.2e-6_REAL64*tu))
    gmst = mod(t/3600.0_REAL64 + 24.0_REAL64*(jd - aint(jd)) + 12.0_REAL64,24.0_REAL64)

  end function gmst


  pure function lmst(jd,longitude)

    real(REAL64) :: lmst
    real(REAL64), intent(in) :: jd, longitude


!  local sidereal time in hours
!
!  jd is a full Julian date
!  lambda is a longitude in degrees: -west ... +east
!
!  The precision is better than 1 second.
!  According to Astronomical Almanac 2000.

    lmst = mod(gmst(jd) + longitude/15.0_REAL64,24.0_REAL64)

  end function lmst



  pure function hangle(lmst,ra)

!   hour angle in degrees

    real(REAL64) :: hangle
    real(REAL64), intent(in) :: lmst, ra

    hangle = mod(lmst - ra,360.0_REAL64)

  end function hangle


  pure subroutine eq2hor(ha, dec, latitude, az, elev)

    real(REAL64), intent(in) :: ha,dec,latitude
    real(REAL64), intent(out) :: az, elev


!
! equatorial to horizontal coordinates
!
!  all arguments in degrees
!

    real(REAL64) :: sinh, cosh, sind, cosd, sinl, cosl, x,y,z,r

    sinh = sin(ha/RAD)
    cosh = cos(ha/RAD)
    sind = sin(dec/RAD)
    cosd = cos(dec/RAD)
    sinl = sin(latitude/RAD)
    cosl = cos(latitude/RAD)

    x = -cosh*cosd*sinl + sind*cosl
    y = -sinh*cosd
    z = cosh*cosd*cosl + sind*sinl

    r = hypot(x,y)
    if( abs(r) > epsilon(r) )then
       az = RAD*atan2(y,x)
    else
       az = 0.0_REAL64
    end if
    if( az < 0_REAL64 ) az = az + 360.0_REAL64
    elev = RAD*atan2(z,r)

  end subroutine eq2hor


  pure subroutine hor2eq(az, elev, latitude, ha,dec)

    real(REAL64), intent(in) :: az,elev,latitude
    real(REAL64), intent(out) :: ha, dec

!
!  horizontal to equatorial coordinates
!
!  all arguments in degrees
!

    real(REAL64) :: sina, cosa, sine, cose, sinl, cosl, x, y, z, r

    sina = sin(az/RAD)
    cosa = cos(az/RAD)
    sine = sin(elev/RAD)
    cose = cos(elev/RAD)
    sinl = sin(latitude/RAD)
    cosl = cos(latitude/RAD)

    x = -cosa*cose*sinl + sine*cosl
    y = -sina*cose
    z = cosa*cose*cosl + sine*sinl

    r = hypot(x,y)
    if( abs(r) > epsilon(r) )then
       ha = RAD*atan2(y,x)
    else
       ha = 0.0_REAL64
    endif
    dec = RAD*atan2(z,r)

  end subroutine hor2eq


  pure function refract(z)

    real(REAL64) :: refract
    real(REAL64), intent(in) :: z

!
!     compute refraction angle in degrees
!
!    Smart: Textbook on spherical astronomy
!
!    constants for pressure 760mmHg, 10deg C with
!    suffucient accuracy for z < 75 deg
!

    real(REAL64) :: tanz

    tanz = tan(z/RAD)
    refract = (58.16_REAL64*tanz - 0.067_REAL64*tanz*tanz*tanz)/3600.0_REAL64

  end function refract


  pure function airmass(z)

    real(REAL64) :: airmass
    real(REAL64), intent(in) :: z

    !
    !     compute airmass,
    !
    !      young&irvine: aj,72,945,(1967)
    !
    ! the airmass is limited on the given range of zenit distances

    real(REAL64) :: secz

    if( 0 <= z .and. z < 86.5 ) then
       secz = 1.0_REAL64/cos(z/RAD)
       airmass = secz*(1.0_REAL64 - 1.2e-3_REAL64*(secz**2 - 1.0_REAL64))
    else
       airmass = -1
    end if

  end function airmass

  pure function xairmass(jd,long,lat,ra,dec)

    real(REAL64) :: xairmass
    real(REAL64), intent(in) :: jd,long,lat,ra,dec
    real(REAL64) :: t,h,ha,a

    t = lmst(jd,long)
    ha = hangle(15.0_REAL64*t,ra)
    call eq2hor(ha,dec,lat,a,h)
    xairmass = airmass(90.0_REAL64 - h)

  end function xairmass


  pure function longsun(d)

    real(REAL64), intent(in) :: d   ! days since 1. january
    real(REAL64) :: longsun

   !  approx (!!!!) of length of the Sun
    longsun = mod(279.465 + 0.985647*d,360.0_REAL64)

  end function longsun


  pure function helcor(alpha,delta,ls)

    real(REAL64), intent(in) ::  alpha,delta,ls
    real(REAL64) ::  helcor

    ! heliocentric correction in days, angles in degrees

    helcor = 0.9174077_REAL64*sin(alpha/rad)*cos(delta/rad) + &
         0.3979486_REAL64*sin(delta/rad)
    helcor = helcor*sin(ls/rad) + cos(ls/rad)*cos(alpha/rad)*cos(delta/rad)
    helcor = -0.0057755_REAL64*helcor

  end function helcor

  pure function phase(jd,min0,per)

    real(REAL64), intent(in) ::  jd,min0,per
    real(REAL64) ::  phase

    phase = mod(jd - min0,per) / per

    ! Phase is negative for jd < min0.

  end function phase

end module astrosphere
