!
! Celestial projections
!
! Copyright © 2011-2026 F.Hroch (hroch@physics.muni.cz)
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
! All angles are in degrees!

module projections

  use iso_fortran_env
  implicit none

  real(REAL64), private, parameter :: rad = 180.0_REAL64/3.1415926535897931160_REAL64

contains

  pure elemental subroutine carree(alpha, delta, alpha0, delta0, x,y)

    real(REAL64), intent(in) :: alpha, delta, alpha0, delta0
    real(REAL64), intent(out) :: x,y

    x = alpha - alpha0
    y = delta - delta0

  end subroutine carree

  pure elemental subroutine invcarree(x,y,alpha0,delta0,alpha,delta)

    real(REAL64), intent(in) :: alpha0, delta0, x, y
    real(REAL64), intent(out) :: alpha, delta

    alpha = x + alpha0
    delta = y + delta0

  end subroutine invcarree


  pure elemental subroutine gnomon(alpha, delta, alpha0, delta0, x,y)

    real(REAL64), intent(in) :: alpha, delta, alpha0, delta0
    real(REAL64), intent(out) :: x,y
    real(REAL64) :: c,p,q,r,v,w,s

    c = cos(delta / rad)
    p = sin(delta / rad)
    q = c*sin((alpha - alpha0) / rad)
    r = c*cos((alpha - alpha0) / rad)
    v = sin(delta0 / rad)
    w = cos(delta0 / rad)
    s = p*v + r*w

    if( abs(s) > 0 ) then
       x = - q / s * rad
       y = (w*p - v*r) / s * rad
    else
       x = 0
       y = 0
    end if

  end subroutine gnomon

  pure elemental subroutine invgnomon(xx,yy,alpha0,delta0,alpha,delta)

    real(REAL64), intent(in) :: alpha0, delta0, xx, yy
    real(REAL64), intent(out) :: alpha, delta
    real(REAL64) :: p,q,r,v,w,t,x,y

    x = xx / rad
    y = yy / rad
    v = sin(delta0 / rad)
    w = cos(delta0 / rad)
    t = sqrt(1.0_REAL64 + x**2 + y**2)
    p = (v + w*y)/t
    q = -x/t
    r = (w - v*y)/t

    delta = rad * asin(p)
    if( abs(q) > 0 .and. abs(r) > 0 ) then
       alpha = rad * atan2(q,r) + alpha0
    else
       alpha = alpha0
    end if

  end subroutine invgnomon

  subroutine gnomond(alpha, delta, alpha0, delta0, x,y)

    real(REAL64), intent(in) :: alpha, delta, alpha0, delta0
    real(REAL64), intent(out) :: x,y

    call gnomon(alpha/rad, delta/rad,alpha0/rad,delta0/rad,x,y)
    x = rad*x
    y = rad*y
    call gnomon(alpha,delta,alpha0,delta0,x,y)

  end subroutine gnomond

  subroutine invgnomond(x,y,alpha0,delta0,alpha,delta)

    real(REAL64), intent(in) :: alpha0, delta0, x, y
    real(REAL64), intent(out) :: alpha, delta

    call invgnomon(x/rad,y/rad,alpha0/rad,delta0/rad,alpha,delta)
    alpha = rad*alpha
    delta = rad*delta

    call invgnomon(x,y,alpha0,delta0,alpha,delta)

  end subroutine invgnomond

  subroutine identity(alpha, delta, alpha0, delta0, x,y)

    real(REAL64), intent(in) :: alpha, delta, alpha0, delta0
    real(REAL64), intent(out) :: x,y

    x = alpha - alpha0
    y = delta - delta0

  end subroutine identity

  subroutine invidentity(x,y,alpha0,delta0,alpha,delta)

    real(REAL64), intent(in) :: alpha0, delta0, x, y
    real(REAL64), intent(out) :: alpha, delta

    alpha = x + alpha0
    delta = y + delta0

  end subroutine invidentity


end module projections
