!
!  An initial construction of an individual growth-curve of a star.
!
!  Copyright © 2016 F.Hroch (hroch@physics.muni.cz)
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

module grow_init

  implicit none

  ! numerical precision of real numbers
  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: pi = 3.14159265358979312_dbl

  ! print debug informations ?
  logical, parameter, private :: verbose = .false.

  real(dbl), dimension(:), allocatable, private :: cts,dcts,aper,raper
  real(dbl), private :: maxback, hwhm, par_a, par_b
  integer, private :: nhwhm, naper
  logical, private :: second

  private :: minall, the_constructor

contains




  subroutine growinit(xraper,xcts,xdcts,skyerr,rhwhm,mhwhm,grow,dgrow,valid,a,b)

    use NelderMead

    real, dimension(:), intent(in) :: xraper
    real(dbl), dimension(:), intent(in) :: xcts,xdcts
    real, intent(in) :: rhwhm
    real(dbl), intent(in) :: skyerr
    integer, intent(in) :: mhwhm
    real(dbl), dimension(:), intent(out) :: grow,dgrow
    logical, intent(out) :: valid
    real(dbl), intent(in), optional :: a,b
    real(dbl), dimension(2) :: p,dp
    real(dbl) :: pmin
    integer :: ifault

    if( present(a) .and. present(b) ) then
       second = .true.
       par_a = a
       par_b = b
!       write(*,*) a,b
    else
       second = .false.
       par_a = 0
       par_b = 0
    end if

    naper = size(xraper)
    allocate(raper(naper),aper(naper),cts(naper),dcts(naper))

    raper = xraper
    aper = pi*raper**2
    cts = xcts
    dcts = xdcts
    hwhm = rhwhm
    nhwhm = mhwhm

    ! estimate parameters
    p(1) = 1
    p(2) = 0
    dp(1) = 0.01
    dp(2) = 0.1*skyerr
    maxback = 10*skyerr
    call nelmin1(minall,p,dp,pmin,ifault)
    valid = ifault == 0 .and. p(2) < 0.99*maxback

    ! update grow curve
    if( valid ) call the_constructor(p(1),p(2),grow,dgrow,valid)

    deallocate(raper,aper,cts,dcts)

  end subroutine growinit


  function minall(p)

    real(dbl) :: minall
    real(dbl), dimension(:), intent(in) :: p

    real(dbl), dimension(:), allocatable :: grow,dgrow
    real(dbl) :: s,t,b,d,a,x,g,x1,x2,gx
    integer :: i,n
    logical :: valid

    t = p(1)
    b = p(2)

    ! range check of parameters
    if( .not. (0.5 < t .and. t < 1.0 ) .or.  abs(b) > maxback) then
       minall = 1e5
       return
    end if

    ! determine grow curve with actual parameters
    allocate(grow(naper),dgrow(naper))
    call the_constructor(t,b,grow,dgrow,valid)
    if( .not. valid ) then
       minall = 1e6
       return
    end if

    ! asymptotic parameter, 1-order
    n = (naper + nhwhm) / 2
    a = (raper(nhwhm) / hwhm) * (1 - grow(nhwhm))
    a = (raper(n) / hwhm) * (1 - grow(n))
!    a = (raper(nhwhm) / hwhm)**2 * (1 - grow(nhwhm))
!    write(*,*) a


!    g = (grow(naper)*raper(naper) - grow(n)*raper(n)) / (raper(naper) - raper(n))
!    a = (grow(naper) - grow(n)) / (1/raper(n) - 1/raper(naper))

!    x1 = raper(naper) / hwhm
!    x2 = raper(n) / hwhm
!    g = (grow(naper)*x2 - grow(n)*x1) / (x2 - x1)
!    a = (grow(naper) - grow(n)) / (1/x1 - 1/x2)

!    write(*,*) a,g

    if( second ) then
       gx = 1 - 1/ (par_a*raper(naper) - par_b)
    else
       gx = 1
    end if

    ! asympotic estimate for larger apertures
    s = 0
    do i = (nhwhm+naper)/2,naper
       d = 1-0.00177*0 - a / (raper(i) / hwhm)
       d = 1 - a / (raper(i) / hwhm)**2
       x = raper(i) / hwhm
!       d = 1 - exp(-x**2)/1.77245*(1/x-0.5/x**2)
       !       d = 1-exp(-x)/(x)
       if( second ) then
          d = (1 - 1/ (par_a*raper(i) - par_b))/gx
       else
          d = 1 - a / x
       end if
       s = s + abs(grow(i) - d)
    end do

    minall = s

    deallocate(grow,dgrow)

  end function minall

  subroutine the_constructor(t,b,grow,dgrow,valid)

    real(dbl), intent(in) :: b,t
    real(dbl), dimension(:), intent(out) :: grow,dgrow
    logical, intent(out) :: valid

    real(dbl), dimension(size(grow)) :: flux,dflux
    integer :: i

    valid = .false.

    flux = cts - b*aper
    dflux = dcts

    if( any(flux < epsilon(flux)) ) return

    ! set up initial point of grow curve
    grow(nhwhm) = t
    dgrow(nhwhm) = t * dcts(nhwhm)/cts(nhwhm)

    ! compute grow curve
    do i = nhwhm-1,1,-1
       grow(i) = grow(i+1)*(flux(i) / flux(i+1))
       dgrow(i) = grow(i)*sqrt((dcts(i)/cts(i))**2 + (dcts(i+1)/cts(i+1))**2)/1.41
    end do

    do i = nhwhm+1,size(grow)
       grow(i) = grow(i-1)*(flux(i) / flux(i-1))
       dgrow(i) = grow(i-1)*sqrt((dcts(i)/cts(i))**2 + (dcts(i-1)/cts(i-1))**2)/1.41
    end do

    valid = .true.

  end subroutine the_constructor

end module grow_init
