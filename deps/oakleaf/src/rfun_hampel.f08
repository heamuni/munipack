!
!  Hampel's robust functions
!
!  Copyright © 2022-2026 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Oakleaf.
!
!  Oakleaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Oakleaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!
!  There are two implementations of the smooth Hampel's function.
!  The implementation based on splines is more clear, but about 25% slower.
!  The if-then-else constructions looks to be implemented faster:
!  https://en.wikipedia.org/wiki/Branch_predictor

module rfun_hampel_REAL64

  use iso_fortran_env
  implicit none
  private

  real(REAL64), parameter :: a = 1.65_REAL64
  real(REAL64), parameter :: b = 3.3_REAL64
  real(REAL64), parameter :: c = 8.525_REAL64

  real(REAL64), parameter :: R = 0.275_REAL64
  real(REAL64), parameter :: q = a/(c - b)
  real(REAL64), parameter :: half = 0.5_REAL64
  real(REAL64), parameter :: halfq = q / 2
  real(REAL64), parameter :: rr = 1 / (2*R)
  real(REAL64), parameter :: za = 2*a - R/2
  real(REAL64), parameter :: zb = -2*a/q + R/2
  real(REAL64), parameter :: zc = -R/2
  real(REAL64), parameter :: a1 = a - R, a2 = a + R
  real(REAL64), parameter :: b1 = b - R, b2 = b + R
  real(REAL64), parameter :: c1 = c - R, c2 = c + R
  real(REAL64), parameter :: da = (a-R)**2 + za*R - 2*R**2/3
  real(REAL64), parameter :: d1 = (za*R + R**2/3 + da)/2 - a*(a + R)
  real(REAL64), parameter :: db = zb*R - R**2/3 - 2*(a*(b-R) + d1)/q
  real(REAL64), parameter :: &
         d2 = -(db + zb*R + 2*R**2/3)/2 + (b+R)**2/2 - c*(b+R)
  real(REAL64), parameter :: &
       dc = zc*R - 2*R**2/3 - 2*(c*(c-R) - (c-R)**2/2 + d2)
  real(REAL64), parameter :: d3 = -q*(zc*R + R**2/3 + dc)/2
  real(REAL64), parameter :: d = a*b - a**2/2 - q*(b*c - b**2/2)
  real(REAL64), parameter :: six = q*c**2/2 + d

  integer, parameter :: sdim = 32
  real(REAL64), parameter :: h = 0.275_REAL64
  real(REAL64), dimension(0:sdim), parameter :: spline_b = [ &
       0.0_REAL64, h, 2*h, 3*h, 4*h, 5*h, &
       za/2, (za+h/2)/2, a, a, a, a, &
       -q*zb/2, q*(c-(b+1*h)), q*(c-(b+2*h)), q*(c-(b+3*h)), &
       q*(c-(b+4*h)), q*(c-(b+5*h)), q*(c-(b+6*h)), q*(c-(b+7*h)), &
       q*(c-(b+8*h)), q*(c-(b+9*h)), q*(c-(b+10*h)), q*(c-(b+11*h)), &
       q*(c-(b+12*h)), q*(c-(b+13*h)), q*(c-(b+14*h)), q*(c-(b+15*h)), &
       q*(c-(b+16*h)), q*(c-(b+17*h)), q*(c-(b+18*h)), &
       -q*zc/2, 0.0_REAL64 ]

  real(REAL64), dimension(0:sdim), parameter :: spline_c = [ &
       half, half, half, half, half, half, half/2, &
       0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, -q/4, &
       -q/2, -q/2, -q/2, -q/2, -q/2, -q/2, -q/2, -q/2, -q/2, -q/2, -q/2, &
       -q/2, -q/2, -q/2, -q/2, -q/2, -q/2, -q/2, -q/4, 0.0_REAL64 ]

  real(REAL64), dimension(0:sdim), parameter :: spline_d = [ &
       0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, &
       -1/b, -1/b, &
       0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, &
       -q/b, -q/b, &
       0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, &
       0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, &
       0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, 0.0_REAL64, &
       q/b, q/b, 0.0_REAL64 ]

  real(REAL64), dimension(0:sdim), parameter :: spline_y = [ &
       0.0_REAL64, h**2/2, (2*h)**2/2, (3*h)**2/2, (4*h)**2/2, &
       (da + (-h)*(za + (-h)*(half - (-h)/(6*R))))/2, da/2, &
       (da + h*(za + h*(half - h/(6*R))))/2, &
       a*(8*h) + d1, a*(9*h) + d1, a*(10*h) + d1, &
       -q*(db + (-h)*(zb + (-h)*(half + (-h)/(6*R))))/2,&
       -q*(db + (0*h)*(zb + (0*h)*(half + (0*h)/(6*R))))/2,&
       -q*(db + (1*h)*(zb + (1*h)*(half + (1*h)/(6*R))))/2,&
       q*((14*h)*(c - (14*h)/2) + d2), q*((15*h)*(c - (15*h)/2) + d2), &
       q*((16*h)*(c - (16*h)/2) + d2), q*((17*h)*(c - (17*h)/2) + d2), &
       q*((18*h)*(c - (18*h)/2) + d2), q*((19*h)*(c - (19*h)/2) + d2), &
       q*((20*h)*(c - (20*h)/2) + d2), q*((21*h)*(c - (21*h)/2) + d2), &
       q*((22*h)*(c - (22*h)/2) + d2), q*((23*h)*(c - (23*h)/2) + d2), &
       q*((24*h)*(c - (24*h)/2) + d2), q*((25*h)*(c - (25*h)/2) + d2), &
       q*((26*h)*(c - (26*h)/2) + d2), q*((27*h)*(c - (27*h)/2) + d2), &
       q*((28*h)*(c - (28*h)/2) + d2), q*((29*h)*(c - (29*h)/2) + d2), &
       -q*(dc + (-h)*(zc + (-h)*(half - (-h)/(6*R))))/2, -q*dc/2, d3 ]
  public :: spline_y, spline_b, spline_c, spline_d

  public :: ipshampel_REAL64, ishampel_REAL64, ihampel_REAL64, &
       pshampel_REAL64, shampel_REAL64, hampel_REAL64, &
       dpshampel_REAL64, dshampel_REAL64, dhampel_REAL64, &
       hampelsqpf_REAL64

contains

  elemental pure function pshampel_REAL64(x) result(pshampel)

    real(REAL64):: pshampel
    real(REAL64), intent(in) :: x
    real(REAL64) :: u

    u = abs(x)
    if( u < a1 ) then
       pshampel = x
    else if( u < a2 ) then
       u = u - a
       pshampel = half*sign(za + u*(1 - rr*u),x)
    else if( u < b1 ) then
       pshampel = sign(a,x)
    else if( u < b2 ) then
       u = u - b
       pshampel = halfq*sign(zb + u*(1 + rr*u),x)
    else if( u < c1 ) then
       pshampel = sign(q*(c - u),x)
    else if( u < c2 ) then
       u = u - c
       pshampel = halfq*sign(zc + u*(1 - rr*u),x)
    else
       pshampel = 0
    endif

  end function pshampel_REAL64

  elemental pure function dpshampel_REAL64(x) result(dpshampel)

    real(REAL64) :: dpshampel
    real(REAL64), intent(in) :: x
    real(REAL64), parameter :: R1 = 1 / R
    real(REAL64) :: u

    u = abs(x)
    if( u < a1 ) then
       dpshampel = 1
    else if( u < a2 ) then
       dpshampel = half*(1 - (u - a)*R1)
    else if( u < b1 ) then
       dpshampel = 0
    else if( u < b2 ) then
       dpshampel = -halfq*(1 + (u - b)*R1)
    else if( u < c1 ) then
       dpshampel = -q
    else if( u < c2 ) then
       dpshampel = -halfq*(1 - (u - c)*R1)
    else
       dpshampel = 0
    endif

  end function dpshampel_REAL64

  elemental pure function ipshampel_REAL64(x) result(ipshampel)

    real(REAL64):: ipshampel
    real(REAL64), intent(in) :: x
    real(REAL64), parameter :: r6 = 1 / (6*R)
    real(REAL64) :: u

    u = abs(x)
    if( u < a1 )then
       ipshampel = half*x**2
    else if( u < a2 ) then
       u = u - a
       ipshampel = half*(da + u*(za + u*(half - r6*u)))
    else if( u < b1 ) then
       ipshampel = a*u + d1
    else if( u < b2 ) then
       u = u - b
       ipshampel = -halfq*(db + u*(zb + u*(half + r6*u)))
    else if( u < c1 ) then
       ipshampel = q*(u*(c - half*u) + d2)
    else if( u < c2 ) then
       u = u - c
       ipshampel = -halfq*(dc + u*(zc + u*(half - r6*u)))
    else ! if ( u >= c ) then
       ipshampel = d3
    end if

  end function ipshampel_REAL64


  elemental pure function shampel_REAL64(x) result(shampel)

    real(REAL64):: shampel
    real(REAL64), intent(in) :: x
    real(REAL64), dimension(0:sdim), parameter :: b = spline_b
    real(REAL64), dimension(0:sdim), parameter :: c2 = 2*spline_c
    real(REAL64), dimension(0:sdim), parameter :: d3 = 3*spline_d
    real(REAL64), parameter :: h1 = 1 / h
    real(REAL64) :: u, t
    integer :: k

    u = abs(x)
    k = int(u*h1)

    select case(k)
    case(0:sdim)
       t = u - k*h
       shampel = sign(b(k) + t*(c2(k) + d3(k)*t),x)
    case default
       shampel = 0
    end select

  end function shampel_REAL64

  elemental pure function dshampel_REAL64(x) result(dshampel)

    real(REAL64):: dshampel
    real(REAL64), intent(in) :: x
    real(REAL64), dimension(0:sdim), parameter :: c2 = 2*spline_c
    real(REAL64), dimension(0:sdim), parameter :: d6 = 6*spline_d
    real(REAL64), parameter :: h1 = 1 / h
    real(REAL64) :: u, t
    integer :: k

    u = abs(x)
    k = int(u*h1)

    select case(k)
    case(0:sdim)
       t = u - k*h
       dshampel = c2(k) + d6(k)*t
    case default
       dshampel = 0
    end select

   end function dshampel_REAL64

  elemental pure function ishampel_REAL64(x) result(ishampel)

    real(REAL64):: ishampel
    real(REAL64), intent(in) :: x
    real(REAL64), dimension(0:sdim), parameter :: y = spline_y
    real(REAL64), dimension(0:sdim), parameter :: b = spline_b
    real(REAL64), dimension(0:sdim), parameter :: e = spline_c
    real(REAL64), dimension(0:sdim), parameter :: d = spline_d
    real(REAL64), parameter :: h1 = 1 / h
    real(REAL64) :: u,t
    integer :: k

    u = abs(x)
    k = int(u*h1)

    select case(k)
    case(0:sdim)
       t = u - k*h
       ishampel = y(k) + t*(b(k) + t*(e(k) + t*d(k)))
    case default
       ishampel = d3
    end select

  end function ishampel_REAL64

  elemental pure subroutine hampelsqpf_REAL64(x,rho,psi)

    real(REAL64), intent(in) :: x
    real(REAL64), intent(out) :: rho, psi
    real(REAL64), parameter :: r6 = 1 / (6*R)
    real(REAL64) :: u

    u = abs(x)

    if( u < a1 ) then
       psi = x
       rho = half*x**2
    else if( u < a2 ) then
       u = u - a
       psi = half*sign(za + u*(1 - rr*u),x)
       rho = half*(da + u*(za + u*(half - r6*u)))
    else if( u < b1 ) then
       psi = sign(a,x)
       rho = a*u + d1
    else if( u < b2 ) then
       u = u - b
       psi = halfq*sign(zb + u*(1 + rr*u),x)
       rho =-halfq*(db + u*(zb + u*(half + r6*u)))
    else if( u < c1 ) then
       psi = sign(q*(c - u),x)
       rho = q*(u*(c - half*u) + d2)
    else if( u < c2 ) then
       u = u - c
       psi = halfq*sign(zc + u*(1 - rr*u),x)
       rho =-halfq*(dc + u*(zc + u*(half - r6*u)))
    else
       psi = 0
       rho = d3
    endif

  end subroutine hampelsqpf_REAL64

  elemental pure function hampel_REAL64(x) result(hampel)

    real(REAL64):: hampel
    real(REAL64), intent(in) :: x
    real(REAL64) :: u

    u = abs(x)
    if( u < a )then
       hampel = x
    else if( a <= u .and. u < b )then
       hampel = sign(a,x)
    else if( b <= u .and. u < c )then
       hampel = sign(q*(c - u),x)
    else ! if( u >= c )then
       hampel = 0
    endif

  end function hampel_REAL64

 elemental pure function dhampel_REAL64(x) result(dhampel)

    real(REAL64) :: dhampel
    real(REAL64),intent(in) :: x
    real(REAL64) :: u

    u = abs(x)
    if( u < a )then
       dhampel = 1
    else if( b <= u .and. u < c )then
       dhampel = -q
    else
       dhampel = 0
    endif

  end function dhampel_REAL64

  elemental pure function ihampel_REAL64(x) result(ihampel)

    real(REAL64):: ihampel
    real(REAL64), intent(in) :: x
    real(REAL64) :: u

    u = abs(x)
    if( u < a )then
       ihampel = x**2 / 2
    else if( a <= u .and. u < b )then
       ihampel = a*(u - a/2)
    else if( b <= u .and. u < c )then
       ihampel = q*u*(c - u/2) + d
    else !if( u >= c )then
       ihampel = six
    endif

  end function ihampel_REAL64

end module rfun_hampel_REAL64
