!
! The initial estimation of astrometrical parameters
!
! Copyright © 2015-2026 F.Hroch (hroch@physics.muni.cz)
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

module estimator

  use iso_fortran_env
  implicit none
  private

  real(REAL64), parameter :: rad = 57.295779513082322865_REAL64

  public :: inestim

contains

  subroutine inestim(type,a,d,xx,yy,acen,dcen,sc,dsc,rot,drot,refl,rms,verbose)

    use astrotrafo

    character(len=*),intent(in) :: type
    real(REAL64),dimension(:),intent(in) :: a,d,xx,yy
    real(REAL64),intent(in out) :: acen,dcen
    real(REAL64),intent(out) :: sc,rot,dsc,drot,refl,rms
    logical, intent(in) :: verbose
    real(REAL64), dimension(:), allocatable :: x,y,u,v
    real(REAL64) :: xoff,yoff,dxoff,dyoff,uoff,voff
    integer :: ndat
    type(AstrometryProjection) :: tt
    type(AstrometryGeometry) :: tg

    if( size(a) < 1 ) stop 'InEstim: no points available.'
    ! asi spis < 2, test pro jeden, dva, tri body

    if( verbose ) write(error_unit,*) "=== Parameters estimation ==="

    ndat = size(a)
    allocate(x(ndat),y(ndat),u(ndat),v(ndat))

    ! projection of spherical coordinates to a tangent plane
    call tt % init(type,acen,dcen)
    call tt % proj(a,d,u,v)

    ! the scale
    call escaler(u,v,xx,yy,sc,dsc)
    if( verbose ) write(error_unit,'(a,f10.1,g10.2)') &
         '# estim init scale [pix/deg]: ',sc,dsc

    ! detect the reflexion
    call ereflex(u,v,xx,yy,refl)
    if( verbose ) write(error_unit,'(a,f5.0)') '# estim init reflexion: ',refl

    ! the angle of rotation
    call tg % init(sc,0.0_REAL64,refl)
    call tg % uv2pix(xx,yy,x,y)
    call erot(u,v,x,y,rot,drot)
    if( verbose ) write(error_unit,'(a,f8.3,g10.2)') &
         '# estim init angle [deg]: ',rot,drot

    ! offsets
    ! origin of rectangular coordinates is moved to center of frame
    ! because that's point the frame is rotated around
    call tg % init(sc,rot,refl)
    call tg % uv2pix(u,v,x,y)
    call eoff(x,y,xx,yy,xoff,dxoff,yoff,dyoff)
    if( verbose ) write(error_unit,'(a,2f12.1,2g10.2)') &
         '# estim init offset [pix]: ',xoff,yoff,dxoff,dyoff

    ! update the center of projection
    call tg % pix2uv(xoff,yoff,uoff,voff)
    call tt % invproj(uoff,voff,acen,dcen)

    if( verbose ) write(error_unit,'(a,2f12.3)') &
         '# estim init centre of projection [deg]: ',acen,dcen

    rms = hypot(dxoff,dyoff) / sc ! [deg]

  end subroutine inestim


  subroutine escaler(u,v,x,y,c,dc)

    use oakleaf

    real(REAL64), dimension(:),intent(in) :: x,y,u,v
    real(REAL64), intent(out) :: c,dc

    integer :: i,j,n
    real(REAL64),dimension(:),allocatable :: cc
    real(REAL64) :: d1,d2

    allocate(cc(size(x)**2))

    n = 0
    do i = 2, size(x)
       do j = 1,i-1
          d1 = hypot(x(j) - x(i),y(j) - y(i))
          d2 = hypot(u(j) - u(i),v(j) - v(i))
          if( d1 > 0 .and. d2 > 0 ) then
             n = n + 1
             cc(n) = d1 / d2
          end if
!         if(debug) write(*,'(a,2es15.8,3f15.5)') '#',x(i),y(i),u(i),v(i),d2/d1
       end do
    enddo

    if( n > 1 ) then
       call rmean(cc(1:n),c,stdsig=dc)
    else
       c = 0
       dc = -1
       stop 'Error in escaler: the scale is undefined due lack of points.'
    end if

  end subroutine escaler

  subroutine erot(u,v,x,y,f,df)

    use oakleaf

    real(REAL64), dimension(:),intent(in) :: x,y,u,v
    real(REAL64), intent(out) :: f,df

    integer :: i,j,n
    real(REAL64),dimension(:),allocatable :: fs
    real(REAL64) :: d1,d2,x1,y1,x2,y2,c,s,d12

    n = size(x) - 1
    allocate(fs(n))

    n = 0
    do i = 1, size(x) - 1
       j = i + 1
       x1 = x(j) - x(i)
       y1 = y(j) - y(i)
       d1 = hypot(x1,y1)
       x2 = u(j) - u(i)
       y2 = v(j) - v(i)
       d2 = hypot(x2,y2)
       d12 = d1*d2
       if( d12 > 0 ) then
          n = n + 1
          ! cos
          c = (x1*x2 + y1*y2) / d12
          ! sin
          s = (x1*y2 - x2*y1) / d12
          !if(debug)write(*,'(a,6f12.5)') '# erot ',x(i),y(i),u(i),v(i),cc(i),cs(i)
          fs(n) = rad * sign(atan2(abs(s),c),s)
!          fs(n) = rad * atan2(s,c)
       endif
    enddo

    if( n > 1 ) then
       call rmean(fs(1:n),f,stdsig=df)
    else if( n == 1 ) then
       f = fs(1)
       df = 0
    else
       f = 0
       df = 0
       write(error_unit,*) 'Warning in erot: the angle is undefined.'
    end if

  end subroutine erot


  pure function angle(x1,y1,x2,y2,x3,y3)

    real(REAL64), intent(in) :: x1,y1,x2,y2,x3,y3
    real(REAL64) :: angle,u1,u2,v1,v2,c,s,r1,r2,r12

    u1 = x1 - x2
    u2 = x3 - x2
    v1 = y1 - y2
    v2 = y3 - y2
    r1 = hypot(u1, v1)
    r2 = hypot(u2, v2)
    r12 = r1*r2
    if( r12 > 0 ) then
       c = (u1*u2 + v1*v2) / r12
       s = (u1*v2 - u2*v1) / r12
       angle = sign(atan2(abs(s),c),s)
    else
       angle = 0
    end if

  end function angle


  pure subroutine ereflex(u,v,x,y,refl)

    real(REAL64), dimension(:),intent(in) :: x,y,u,v
    real(REAL64), intent(out) :: refl

    integer :: i,j,k
    integer, dimension(:), allocatable :: r
    real(REAL64) :: f1,f2

    allocate(r(size(x)-2))

    do i = 1, size(x)-2
       j = i + 1
       k = i + 2
       f1 = angle(x(i),y(i),x(j),y(j),x(k),y(k))
       f2 = angle(u(i),v(i),u(j),v(j),u(k),v(k))

       r(i) = nint(sign(1.0_REAL64,f1*f2))
!       write(error_unit,'(a,3f7.2,i3)') '# ereflex ',f1,f2,f1*f2,r(i)
    enddo

    refl = sign(1,sum(r))

  end subroutine ereflex

  subroutine eoff(u,v,x,y,xoff,dxoff,yoff,dyoff)

    use oakleaf

    real(REAL64), dimension(:), intent(in) :: x,y,u,v
    real(REAL64), intent(out) :: xoff,dxoff,yoff,dyoff

    call rmean(u-x,xoff,stdsig=dxoff)
    call rmean(v-y,yoff,stdsig=dyoff)

  end subroutine eoff

end module estimator
