!
!  Colour extinction
!
!  Copyright © 2015, 2017-8 F.Hroch (hroch@physics.muni.cz)
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


module colorex

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)

  ! sqrt(pi/2) - obsolete
  real(dbl), parameter, private :: sqrtpi2 = 1.2533141373155001_dbl
  ! sqrt(2*pi)
  real(dbl), parameter, private :: sqrt2pi = 2.50662721600161359133_dbl

  ! sqrt(2)
  real(dbl), parameter, private :: sqrt2 = 1.41421356237309504878_dbl

  ! Constants for Gauss-Hermite Integration of second order
  real(dbl), parameter, private :: &
       h1 = 0.886226925453_dbl, &
       a1 = 0.707106781187_dbl

  private :: extin, attenuation

contains

  subroutine correx(leff,dlam,airmass,k0,r,lref,flux)

    real(dbl), dimension(:), intent(in) :: leff,dlam,airmass
    real(dbl), dimension(:), intent(in out) :: flux
    real(dbl), intent(in) :: k0,r,lref

    real(dbl), dimension(size(leff)) :: flam
    real(dbl) :: a0,ax,f1,f2,dx
    integer :: i,n

    n = size(leff)

!    write(*,'(a,5f15.3)') "in:",flux

    ! sanity check for undefined fluxes, correction is omitted
    if( any(flux < 0) ) return

    ! compute flux density from fluxes (replace by flam2flux from photoconv module?)
!    flam = flux / (1e9 * sqrtpi2 * dlam)
    flam = flux / (sqrtpi2 * dlam)

!    write(*,'(a,5f15.3)') 'flam:',flam
!    write(*,'(a,5f15.3)') 'ctph:',ctph

    ! estimate effective extinction
!    call extin(leff,ctph,airmass,lref,k0,r)

!    write(*,*) 'e',k0,r

!    do i = 1,n
!       flam(i) = flam(i) / attenuation(k0,r,leff(i),airmass(i))
!    end do

    ! approximation of derivations near edges
    call attenuation(k0,r,lref,leff(1),airmass(1),a0,ax)
!    a0 = attenuation(k0,r,leff(1),airmass(1))
!    ax = attenuation1(k0,r,leff(1),airmass(1))
!    flux(1) = 2.0_dbl*(flam(1)*a0/a0 + &
!         0*2*1e9*(flam(2) - flam(1))/(leff(2) - leff(1))*(dlam(1)*a1)**2*ax)
!    write(*,*) 'l',flam(1)*a0, 1e9*(flam(2) - flam(1))/(leff(2) - leff(1))*(dlam(1)*a1)**2*ax
!    write(*,*) (flam(2) - flam(1))/(leff(2) - leff(1))*(sqrt2pi*dlam(1)*a1)**2*(ax/a0)/flam(1)*1e9
!    write(*,*) (flam(2) - flam(1))/flam(1), dlam(1)/(leff(2) - leff(1)),ax/a0

!    write(*,*) (flam(2) - flam(1))/flam(1)* dlam(1)/(leff(2) - leff(1)) * ax/a0, &
!         (flam(2) - flam(1))/flam(1)* dlam(1)/(leff(2) - leff(1)) * ax/a0 * 2* dlam(1)*a1**2*1e9
!    write(*,*) flam(1)*a0,2*1e9*(flam(2) - flam(1))/(leff(2) - leff(1))*(dlam(1)*a1)**2*ax
!    flux(1) = flux(1)*(1 + &
!         2*1e9*(flam(2) - flam(1))/(leff(2) - leff(1))*(dlam(1)*a1)**2*ax/a0/flam(1))
    dx = sqrt2*dlam(1)*a1
    f1 = (flam(2) - flam(1))/(leff(2) - leff(1))
!    flux(1) = flux(1)*(1 - 1e9*f1*ax/a0*dx**2/flam(1))
    flux(1) = flux(1)*(1 - f1*ax/a0*dx**2/flam(1))

!    a0 = attenuation(k0,r,leff(n),airmass(n))
!    ax = attenuation1(k0,r,leff(n),airmass(n))
!    flux(n) = 2.0_dbl*(flam(n)*a0/a0 - &
!         0*2*1e9*(flam(n) - flam(n-1))/(leff(n) - leff(n-1))*(dlam(n)*a1)**2*ax)

!    write(*,*) flam(n)*a0,2*1e9*(flam(n) - flam(n-1))/(leff(n) - leff(n-1))*(dlam(n)*a1)**2*ax
    do i = 2,n-1
       dx = sqrt2*dlam(i)*a1
       call attenuation(k0,r,lref,leff(i),airmass(i),a0,ax)
!       a0 = attenuation(k0,r,leff(i),airmass(i))
!       ax = attenuation1(k0,r,leff(i),airmass(i))
       f1 = (flam(i-1) - flam(i))/(leff(i-1) - leff(i))
       f2 = (flam(i+1) - flam(i))/(leff(i+1) - leff(i))
!       flux(i) = flux(i)*(1 + &
       !            2*1e9*(f2-f1)/2*(dlam(i)*a1)**2*(ax/a0)/flam(i))
!       flux(i) = flux(i)*(1 - ((f2-f1)*0 + 1e9*(f1+f2)*ax/a0*dx)* dx/2/flam(i))
       flux(i) = flux(i)*(1 - (f1+f2)/2*ax/a0*dx**2/flam(i))
!       write(*,*) (f2-f1)* dx/2/flam(i),(f1+f2)*1e9*ax/a0*dx* dx/2/flam(i),((f2-f1) + (f1+f2)*1e9*ax/a0*dx)* dx/2/flam(i)

    end do
!    stop
!    a0 = attenuation(k0,r,leff(n),airmass(n))
!    ax = attenuation1(k0,r,leff(n),airmass(n))
    call attenuation(k0,r,lref,leff(n),airmass(n),a0,ax)
    dx = sqrt2*dlam(n)*a1
    f1 = (flam(n) - flam(n-1))/(leff(n) - leff(n-1))
    flux(n) = flux(n)*(1 - f1*dx**2*ax/a0/flam(n))
!    flux(n) = flux(n)*(1 - 1e9*f1*dx**2*ax/a0/flam(n))


    ! approximation at internal points
!    do i = 2,n-1
!       dx = sqrt2*dlam(i)*a1
!       f1 = flam(i) + (flam(i+1) - flam(i))/(leff(i+1) - leff(i))*dx
!       f2 = flam(i) - (flam(i) - flam(i-1))/(leff(i) - leff(i-1))*dx
!       flux(i) = f1*attenuation(k0,r,leff(i)+dx,airmass(i)) + &
!            f2*attenuation(k0,r,leff(i)-dx,airmass(i))
!       flux(i) = sqrt2*dlam(i)*h1*flux(i) * 1e9
!    end do

!    flux = sqrt2*dlam*h1*flux * 1e9

!    flux = 2*h1*flux

!    write(*,*) flux
!    write(*,'(a,5f15.3)') "out:",flux
!    write(*,*) airmass
!    write(*,*) (2.5*log10(flux(i-1)/flux(i)),i=2,n)

  end subroutine correx

!!$  function attenuation(k0,r,lam,x)
!!$
!!$    real(dbl) :: attenuation
!!$    real(dbl), intent(in) :: k0,r,lam,x
!!$    real(dbl) :: k
!!$
!!$    k = (1e9*lam / k0)**(-r)
!!$    attenuation = exp(-k*x)
!!$
!!$  end function attenuation
!!$
!!$  function attenuation1(k0,r,lam,x)
!!$
!!$    real(dbl) :: attenuation1
!!$    real(dbl), intent(in) :: k0,r,lam,x
!!$    real(dbl) :: k,k1
!!$
!!$    k = (1e9*lam / k0)**(-r)
!!$    k1 = -r*(1e9*lam / k0)**(-r-1) / k0
!!$    attenuation1 = -k1*x*exp(-k*x)
!!$
!!$  end function attenuation1

  subroutine attenuation(k0,r,lref,lam,x,a0,a1)

    real(dbl), intent(in) :: k0,r,lam,x,lref
    real(dbl), intent(out) :: a0,a1
    real(dbl) :: k,k1

    k = k0*(lam /lref)**(-r)
    k1 = -r*k/lam
    a0 = exp(-k*x)
    a1 = -k1*x*a0

  end subroutine attenuation

  subroutine extin(leff,a,x,lref,k0,r)

    use oakleaf

    real(dbl), dimension(:), intent(in) :: leff,a,x
    real(dbl), intent(in) :: lref
    real(dbl), intent(out) :: k0,r

    real(dbl), dimension(size(leff)) :: k,t,s,dt,ds
    real(dbl) :: sig,dk,dr,amin,a0

    amin = minval(a)
    if( amin < 1 ) then
       a0 = 1.11/amin
    else
       a0 = 1
    end if
    k = log(a0*a) / x   ! (-) sign is removed because ctph is reciprocal to a

    t = log(k)
    s = log(leff/lref)
!    dt = log(0.1*k)
!    ds = log(0.01*1e9*leff)
    dt = 1e-2
    ds = 5e-2

!    write(*,*) a0
!    write(*,*) k
!    write(*,*) s
!    write(*,*) t
!    write(*,*) ds
!    write(*,*) dt

    call rline(s,t,k0,r,dk,dr,ds,dt,sig)

!    write(*,*) sig,k0,r
!    write(*,*) 'r=',r
    r = -r

    ! check range, from Mie r=1 to Rayleigh r=4 scaterrings,
    ! this is fall-back way when some points are uncorrect,
    ! r is very sensitive to bad points and I determined
    ! a mean value r=3 for La Silla which means scattering
    ! on fluctations rather than on the dust
    if( .not. (1 < r .and. r < 4) ) then
!    if( .true. ) then
       r = 3
       t = log(k) + r*log(leff/lref)
!       write(*,*) t
       !write(*,*) k
       call rmean(t,k0,dk,sig)
!       write(*,*) sig,k0
    end if

    k0 = exp(k0)

  end subroutine extin

  subroutine scaterr(leff,a,x,k0,r,lref)

    real(dbl), dimension(:), intent(in) :: leff,a,x
    real(dbl), intent(out) :: k0,r,lref

    ! reference effective wavelenght
    lref = leff(max(size(leff)/2+1,1))

    call extin(leff,a,x,lref,k0,r)

  end subroutine scaterr

end module colorex
