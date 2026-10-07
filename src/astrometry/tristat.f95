!
! triangle statistics
!
! fortran -Wall -g -p -fbounds-check -fimplicit-none -I../lib astrofits.f95 tristat.f95 -L../lib -lastrometry -lcfitsio -lm
!
!
! Copyright © 2012 F.Hroch (hroch@physics.muni.cz)
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

program tristat

  use astrotrafo
  use astrofits

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)
  logical , parameter :: monte = .false.
  integer, parameter :: nr = 100
  real(dbl), parameter :: qe = 0.1

  real(dbl), dimension(:), allocatable :: x,y,xe,ye,a,d,xhist
  integer, dimension(-nr:nr) :: hist
  real(dbl), dimension(2) :: crpix
  character(len=80) :: idcat,rframe

  integer :: ns,i,j,k,t,status,ihist
  real(dbl) :: u,v,ue,ve,q,w,h,s

  type(TypeProjection) :: at

  
  hist = 0

  if( monte ) then
     ns = 333
     allocate(x(ns),y(ns),xe(ns),ye(ns),xhist(ns*ns*ns))

     do i = 1, ns
        call random_number(q)
        x(i) = q - 0.5
        call random_number(q)
        y(i) = q - 0.5
     end do

     do i = 1, ns
        xe(i) = x(i) + gnoise(qe)
        ye(i) = y(i) + gnoise(qe)
!        x(i) = x(i) + gnoise(qe)
!        y(i) = y(i) + gnoise(qe)
     end do

  else
     
     status = 0
     call readfile('barnard_01R.fits',w,h,x,y,crpix,status)
     ns = size(x)

     call readtable("bcat.fits","RAJ2000","DEJ2000",a,d,rframe,idcat,status)
     at%type = "GNOMONIC"
     at%acen = 269.44
     at%dcen = 4.74
     at%xcen = 0.0
     at%ycen = 0.0
     at%scale = 1.0
     at%pa = 0.0

     allocate(xe(ns),ye(ns),xhist(ns*ns*ns))
     call trafo(at,a(1:ns),d(1:ns),xe,ye)

  end if


  xhist = 0.0
  ihist = 0

  do i = 1, ns
     do j = i+1, ns - 1
        do k = j+1, ns - 1
           call triangle(x(i),y(i),x(j),y(j),x(k),y(k),u,v)
           call triangle(xe(i),ye(i),xe(j),ye(j),xe(k),ye(k),ue,ve)
!           t = nr*(u - ue)
!           write(*,*) (u-ue)**2+(v-ve)**2
!           write(*,*) abs(u-ue)+abs(v-ve)
!           t = nr*((x(i) - x(j))-(xe(i) - xe(j)))**2 + ((y(i) - y(j))-(ye(i) - ye(j)))**2
!           write(*,*) i,j,t
!           t = nr*t
           t = nr*sqrt((u-ue)**2+(v-ve)**2)
!           t = nr*(u-ue)
!           t = nr*(abs(u-ue)+abs(v-ve))
!           t = nr*(u-ue)
           if( -nr < t .and. t < nr ) hist(t) = hist(t) + 1
           ihist = ihist + 1
!           xhist(ihist) = sqrt((u-ue)**2)
           xhist(ihist) = sqrt((u-ue)**2+(v-ve)**2)
!           xhist(ihist) = abs(u-ue) + abs(v-ve)
!           xhist(ihist) = v-ve
           write(*,*) xhist(ihist)
        end do
     end do
  end do


!  s = sum(hist)
!  do i = -nr,nr
!     if( i /= 0 ) write(*,*) i/100.0,hist(i)/s
!  end do

  call histogram(xhist(1:ihist),'uveuclid')

  deallocate(x,y,xe,ye,xhist)

contains

subroutine triangle(x1,y1,x2,y2,x3,y3,u,v)

  integer, parameter :: rp = selected_real_kind(15)
  real(rp), parameter :: eps = epsilon(1.0_rp)

  real(rp), intent(in) :: x1,y1,x2,y2,x3,y3
  real(rp), intent(out) :: u,v

  real(rp), dimension(3) :: d
  real(rp) :: x
  integer :: i,j,jj

  ! compute sides of triangle
  d(1) = (x1 - x2)**2 + (y1 - y2)**2
  d(2) = (x1 - x3)**2 + (y1 - y3)**2
  d(3) = (x2 - x3)**2 + (y2 - y3)**2

  ! sorting
  do i = 2, 3
     do j = 3, i, -1

        jj = j - 1
        if( d(jj) > d(j) )then

           x = d(jj)
           d(jj) = d(j)
           d(j) = x

        endif
     enddo
  enddo

  ! compute transformed coordinates
  if( d(3) > eps )then
     u = sqrt(d(1)/d(3))
     v = sqrt(d(2)/d(3))
  else
     ! it would be better put u,v to the definition domain
     u = 0.0
     v = 0.0
  endif

end subroutine triangle

   function gnoise(noise)

    real(dbl) :: gnoise, noise
    real(dbl) :: x

    call random_number(x)
    gnoise = noise*invdist(x)
    

  end function gnoise


  function invdist(xx)

    real(dbl) :: invdist,xx

    ! rational approximation of an inverse to a cumulative function
    ! of Gaussian distribution  with precision better than 0.00045
    ! J.Andel: Statistical methods, Matfyz Press, Prague 1991
    
    real(dbl) :: w,f,x
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
    real(dbl), parameter :: h = 0.02
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

end program tristat
