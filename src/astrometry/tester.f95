!
! tester - astrometry
!
!
! $ gfortran -Wall -g -p -fbounds-check -fimplicit-none projections.f95 atrafo.f95 \
!            lsfit.f95 tester.f95 -L../lib -lminpack -lm
!
! Copyright © 2011-2 F.Hroch (hroch@physics.muni.cz)
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

program tester

  use astrotrafo
  use lsfit
  use absfit
  use robfit
  use robustmean

  implicit none

  character(len=*), parameter :: type = 'GNOMONIC'
  integer, parameter :: ns = 9999
  integer, parameter :: dbl = selected_real_kind(15)
  real(dbl), dimension(ns) :: x,y,a,d,a2,d2,x0,y0,dx,dy
  real(dbl), parameter :: a0 = 330.0, d0 = 3.0, r = 0.6
  real(dbl), parameter :: scale0 = 1.e3, pa0 = 0.0
  real(dbl), parameter :: w = 1000, h = 800
  real, parameter :: ge = (1.0/3600.0)*scale0

  real(dbl) :: q,acen,dacen,dcen,ddcen,sc,dsc,pa,dpa,rms,s0,sres,mad,u,v
  integer :: i
  type(TypeProjection) :: t



!  do i = 1,100
!     write(*,*) i/100.0,invdist(i/100.0),0.5*erf(x/sqrt(2.0))+0.5
!  end do
!  stop

  do i = 1,ns
     x(i) = gnoise(0.0,1.0)
!     write(*,*) x(i)
  end do

  call histogram(x,'ds')
  
  call rmean(x,u,v)
  write(*,'(3f10.3)') sqrt(sum(x**2)/ns),sum(abs(x))/ns,v*sqrt(real(ns))
  stop


  do i = 1, ns
     call random_number(q)
     a(i) = a0 + r*(q - 0.5) 
     call random_number(q)
     d(i) = d0 + r*(q - 0.5)
  end do

  t%type = type
  t%acen = a0
  t%dcen = d0
  t%xcen = w/2.0
  t%ycen = h/2.0
  t%scale = scale0
  t%pa = pa0

  call trafo(t,a,d,x0,y0)

  do i = 1, ns
     x(i) = x0(i) + gnoise(0.0,ge)
     y(i) = y0(i) + gnoise(0.0,ge)
!     write(*,*) gnoise(0.14),gnoise(0.14)
!     write(*,'(a,e25.15)') '>', (1.0-ee)*gnoise(0.0,ge)+ee*gnoise(0.0,ge*3.0)
  end do

!  stop

  do i = 1,ns
     write(*,'(2f10.5,2f20.3,2f10.5)') a(i),d(i),x(i),y(i)
  end do

  acen = a0+0.1
  dcen = d0-0.1
  sc = 1.0
  pa = 0
  call lsmin(type,a,d,w/2.0,h/2.0,w,h,x,y,acen,dacen,dcen,ddcen,sc,dsc,pa,dpa,rms,s0,.true.)
!  call absmin(type,a,d,w/2.0,h/2.0,w,h,x,y,acen,dcen,sc,pa,rms,s0,.true.)
  mad = rms
!  call robmin(type,a,d,w/2.0,h/2.0,w,h,x,y,mad,acen,dacen,dcen,ddcen,sc,dsc,pa,dpa,rms,s0,.true.)

  write(*,99) acen,dacen,dacen/acen
  write(*,99) dcen,ddcen,ddcen/dcen
  write(*,99) sc,dsc,dsc/sc
  write(*,99) pa,dpa,dpa/pa
  write(*,*) rms,rms*3600
  write(*,*) s0,s0*3600
99 format(f20.10,'+-',f20.10,'  ',1es20.5)

  ! residuals
  t%acen = acen
  t%dcen = dcen
  t%scale = sc
  t%pa = pa
  sres = 0.0
  call invtrafo(t,x,y,a2,d2)
  do i = 1,ns
     write(*,'(2f10.5,2f10.5,2f8.3)') a(i),d(i),a2(i),d2(i), &
          3600*(a(i)-a2(i))/cos(dcen/57.3),3600*(d(i)-d2(i))
     sres = sres + ((a(i) - a2(i))/cos(dcen/57.3))**2 + (d(i) - d2(i))**2
     dx(i) = 3600*(a(i)-a2(i))/cos(dcen/57.3)
     dy(i) = 3600*(d(i)-d2(i))
  end do

  write(*,*) sqrt(sres/ns/2.0)*3600.0

  call histogram(dx,'histx')
  call histpgram(dy,'histy')


contains

   function gnoise(t,s)

    real :: gnoise
    real :: x,t,s
    real :: e,ee = 0.1

    call random_number(e)
    call random_number(x)
    if( e < ee ) then
       gnoise = 3.0*s*invdist(x)
    else
       gnoise = s*invdist(x)
    end if
!    gnoise = 0.5*erf(x/sqrt(2.0)/s)+0.5
    

  end function gnoise


  function invdist(xx)

    real :: invdist,xx

    ! rational approximation of an inverse to a cumulative function
    ! of Gaussian distribution  with precision better than 0.00045
    ! J.Andel: Statistical methods, Matfyz Press, Prague 1991
    
    real :: w,f,x
    logical :: interval

    x = xx
    if( x < 0.0 ) then
       invdist = 0.0
    elseif( x > 1.0 )then
       invdist = 1.0
    else
       interval = x < 0.5
       if( .not. interval ) x = 1.0 - x + epsilon(1.0)
       w = sqrt(-2.0*log(x));
       f = -w + (2.515517 + w*(0.802853 + w*0.010328))/ &
            (1.0 + w*(1.432788 + w*(0.189269 + w*0.001308)));
       if( interval ) then
          invdist = f
       else
          invdist = -f
       endif
    endif
    
  end function invdist

  
  subroutine histogram(x,name)

    character(len=*) :: name
    real(dbl), dimension(:) :: x

    integer, parameter :: nhist = 99
    integer, dimension(nhist) :: hist
    real(dbl), dimension(nhist) :: bins
    real(dbl), parameter :: h = 0.2
    integer :: i,j
    real(dbl) :: ss

    do i = 1, nhist
       bins(i) = h*(i - nhist/2)
!       write(*,*) i,bins(i)
    end do

  
    hist = 0
    do i = 1,size(x)
       do j = 2,nhist
          if( bins(j-1) < x(i) .and. x(i) <= bins(j)) hist(j) = hist(j) + 1
       end do
    end do

    open(1,file=name)
    ss = sum(hist)*h
    do i = 1,nhist
       write(1,*) bins(i)-h/2.0,hist(i)/ss
    end do
    close(1)

  end subroutine histogram


end program tester
