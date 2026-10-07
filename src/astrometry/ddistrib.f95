!
! compute 2-D distribution of gaussian distributions with offset
!
! fortran -Wall -g -p -fbounds-check -fimplicit-none ddistrib.f95
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



program ddistrib

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)
  real(dbl), dimension(-200:600) :: y,f,g
  real(dbl), dimension(-500:500,-500:500) :: ff,gg,yy
  real(dbl), dimension(-500:500) :: ww
  integer :: i,j
  real(dbl) :: u,v,s

  g = 0.0
  g(0:300) = 2.0
  f = 0.0
  do i = -200,400
     f(i) = 0.5*exp(i/100.0-4.0)
  end do
  
  do i = -200,600
     y(i) = dconv1d(400-i,f,g)
!     write(*,*) i/100.0,g(i),f(i),y(i)/50
  end do
  
  do i = -500,500
     do j = -500,500
        u = i/100.0
        v = j/100.0
        ff(i,j) = exp(-(u**2 + v**2)/2.0)
        gg(i,j) = exp(-((u-1.0)**2 + (v-0.0)**2)/2.0)
!        gg(i,j) = 1.0
     end do
  end do

  do i = -500,500,10
     do j = -500,500,10
        yy(i,j) = dconv2d(i,j,ff,gg)
!        write(*,*) i/100.0,j/100.0,ff(i,j),gg(i,j),yy(i,j)
     end do
  end do

  do i = -500,500,10
     s = 0.0
     do j = -500,500,10
        s = s + yy(i,j)
     end do
     ww(i) = s
     write(*,*) i,s
  end do
        


!  do i = 1,100    
!     d = i/10.0
!     write(*,*) d,dconv1d(d)
!  end do
  
contains

  function dconv1d(j,f,g) result(s)

    real(dbl) :: s
    real(dbl), dimension(:), intent(in) :: f,g
    integer, intent(in) :: j
    integer :: i,l,nmin,nmax
    integer, dimension(1) :: ub,lb

    lb = lbound(f)
    ub = ubound(f)
    nmin = lb(1)
    nmax = ub(1)

    s = 0.0
    do i = nmin,nmax
       l = i - j
       if( nmin < l .and. l < nmax ) &
            s = s + f(i)*g(l)
    end do


  end function dconv1d

  function dconv2d(k,l,f,g) result(s)

    real(dbl) :: s
    real(dbl), dimension(:,:), intent(in) :: f,g
    integer, intent(in) :: k,l
    integer :: i,j,u,v,nxmin,nxmax,nymin,nymax
    integer, dimension(2) :: ub,lb

    lb = lbound(f)
    ub = ubound(f)
    nxmin = lb(1)
    nxmax = ub(1)
    nymin = lb(2)
    nymax = ub(2)

    s = 0.0
    do i = nxmin,nxmax
       u = i - k
       do j = nymin, nymax
          v = j - l
          if( (nxmin < u .and. u < nxmax) .and. (nymin < v .and. v < nymax) ) &
               s = s + f(i,j)*g(u,v)
       end do
    end do


  end function dconv2d

end program ddistrib
