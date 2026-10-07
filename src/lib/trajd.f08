!
!  Julian date conversion module
!
!  Copyright © 1996-2012, 2026 F.Hroch (hroch@physics.muni.cz)
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
! 1996-12-23 F.Hroch
! 1996 to f90
! 2000 some change according to Numerical Recipes (removed)
! 2026 another update of Julian dates conversions
!
! --------------------------------------------------------------------

module trajd

  use iso_fortran_env

  implicit none

contains

  !----------------------------------------------------------------------

  pure integer function jdn (y,m,d)

    ! Julian date number by https://en.wikipedia.org/wiki/Julian_day

    integer, intent(in) :: y,m,d
    integer :: c,e,f,g

    c = (m - 14) / 12
    e = 1461 * (y + 4800 + c)
    f = 367 * (m - 2 - 12 * c)
    g = 3 * ((y + 4900 + c) / 100)
    jdn = e / 4 + f / 12 - g / 4 + d - 32075

  end function jdn

  pure function datjd (year,month,day) result(jd)

    real(REAL64), intent(in) :: year,month,day
    integer :: intd
    real(REAL64) :: d, dt, jd

    d = day - 0.5
    intd = int(d)
    dt = d - intd
    jd = jdn(nint(year),nint(month),intd) + dt

  end function datjd

  pure function datetimejd(year,month,day,hour,min,sec) result(jd)

    integer, intent(in) :: year,month,day,hour,min
    real(REAL64), intent(in) :: sec
    real(REAL64) :: y,m,d, jd

    y = year
    m = month
    d = day + (hour + (min + sec/60.0_REAL64) / 60.0_REAL64) / 24.0_REAL64
    jd = datjd(y,m,d)

  end function datetimejd


  !----------------------------------------------------------------------

  pure function yearjd(year)

    ! Compute Julian date from input year with fraction

    real(REAL64), intent(in) :: year
    real(REAL64) :: yearjd,month,day,f
    integer :: y

    y = int(year)
    f = year - y
    month = int(12*f)
    day = 365.25*(year - (y + month/12.0_REAL64)) + 1
    month = month + 1

    yearjd = datjd(aint(year),month,day)

  end function yearjd


  !----------------------------------------------------------------------

  pure function mjd (jd)

    real(REAL64), intent(in) ::  jd
    real(REAL64) :: mjd

    mjd = jd - 2400000.5_REAL64

  end function mjd

  !----------------------------------------------------------------------

  pure function hjd (jd,helcor)

    real(REAL64), intent(in) ::  jd,helcor
    real(REAL64) :: hjd

    hjd = jd + helcor

  end function hjd

  !----------------------------------------------------------------------

  pure subroutine jdn2date(jdn,year,month,day)

    ! The algorithm due Edward Graham Richards published in
    ! Explanatory Supplement to the Astronomical Almanac, (pp. 76–104), 2013
    ! (https://en.wikipedia.org/wiki/Julian_day)

    integer, intent(in) :: jdn
    integer, intent(out) :: year,month,day

    integer, parameter :: y = 4716
    integer, parameter :: j = 1401
    integer, parameter :: m = 2
    integer, parameter :: n = 12
    integer, parameter :: r = 4
    integer, parameter :: p = 1461
    integer, parameter :: v = 3
    integer, parameter :: u = 5
    integer, parameter :: s = 153
    integer, parameter :: w = 2
    integer, parameter :: BB = 274277
    integer, parameter :: C = -38
    integer :: e,f,g,h

    if( jdn < 0 ) then
       f = jdn + j
    else
       f = jdn + j + (((4*jdn + BB) / 146097) * 3) / 4 + C
    end if

    e = r * f + v
    g = mod(e,p) / r
    h = u * g + w
    day = mod(h,s) / u + 1
    month = mod(h / s + m, n) + 1
    year = e / p - y + (n + m - month) / n

  end subroutine jdn2date


  pure subroutine jdat(jd,y,m,d)

    real(REAL64), intent(in) ::  jd
    real(REAL64), intent(out) ::  y,m,d

    integer ::  year,month,day, jdn
    real(REAL64) :: dt

    jdn = nint(jd)
    dt = jd - jdn
    call jdn2date(jdn,year,month,day)
    y = year
    m = month
    d = day + dt + 0.5_REAL64

  end subroutine jdat

! ---------------

  pure subroutine jdat2(jd,year,month,day)

    ! by Meeus, Jan, Astronomical algorithms

    real(REAL64), parameter :: break = 2299163.0_REAL64

    real(REAL64), intent(in) ::  jd
    real(REAL64), intent(out) ::  year,month,day
    real(REAL64) :: f,jd1
    integer :: z, alpha, a,b,c,d,e

    jd1 = jd + 0.5_REAL64
    z = int(jd1)
    f = jd1 - z

    if( z < 2299161 ) then
       a = z
    else
       alpha = int((z - 1867216.25_REAL64)/36524.25_REAL64)
       a = z + 1 + alpha - int(alpha/4.0_REAL64)
    end if

    b = a + 1524
    c = int((b - 122.1)/365.25)
    d = int(325.25 * c)
    e = int((b - d)/30.6001)

    day = b - d - int(30.6001*e) + f
    if( e < 14 ) then
       month = e - 1
    else
       month = e - 13
    end if
    if( month > 2 ) then
       year = c - 4716
    else
       year = c - 4715
    end if

  end subroutine jdat2

  pure subroutine jdatetime(jd,year,month,day,hour,min,sec)

    real(REAL64), intent(in) ::  jd
    integer, intent(out) ::  year,month,day,hour,min
    real(REAL64), intent(out) ::  sec

    real(REAL64) ::  y,m,d

    call jdat(jd,y,m,d)
    year = nint(y)
    month = nint(m)
    day = int(d)
    hour = int(24*(d - day))
    min = int(1440*(d - (day + hour / 24.0_REAL64)))
    sec = 86400*(d - (day + (hour + min/60.0_REAL64)/24.0_REAL64))

  end subroutine jdatetime


  ! the epoch
  pure function jdepoch(jd)

    real(REAL64), intent(in) :: jd
    real(REAL64) :: y,m,d,jd0,jdepoch

    call jdat(jd,y,m,d)
    jd0 = datjd(y,1.0_REAL64,1.0_REAL64)
    jdepoch = y + (jd - jd0)/365.25_REAL64

  end function jdepoch

end module trajd
