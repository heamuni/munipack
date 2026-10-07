!
!  Visualise traces of cosmic ray events
!
!  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
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


module cosmic_ray_events

  use iso_fortran_env

  implicit none
  private

  public :: cosmics

contains

  subroutine cosmics(nevents,fmax,dmax,sky)

    integer, intent(in) :: nevents
    real(REAL64), intent(in) :: fmax, dmax
    real(REAL64), dimension(:,:), intent(in out) :: sky

    real(REAL64), parameter :: pi2 = 8*atan(1.0_REAL64)
    real(REAL64), dimension(2) :: r1, r2, r3, r
    real(REAL64) :: x1,y1,rr,aa,t,f,dt,t0,ff,x3,y3,phi,dr
    integer :: n,k,m,i,j,width,height

    width = size(sky,1)
    height = size(sky,2)

    do n = 1, nevents

       call random_number(x1)
       call random_number(y1)
       call random_number(rr)
       call random_number(aa)
       call random_number(x3)
       call random_number(y3)

!       x1 = 2*x1 - 0.5
!       y1 = 2*y1 - 0.5
       r1 = [ x1*width, y1*height ]
       dr = rr*dmax
       phi = pi2*aa
       r2 = r1 + dr*[ cos(phi), sin(phi) ]
!       x2 = 2*x2 - 0.5
!       y2 = 2*y2 - 0.5
!       r2 = [ x2*width, y2*height]
       r = r2 - r1
       r3 = r1 + x3*r
       r3 = r3 + sign(0.07*y3*[r(2),r(1)],y3-0.5)

!       write(*,*) 'r1=',r1
!       write(*,*) 'r2=',r2
!       write(*,*) 'r3=',r3,y3*[-r(2),r(1)],dot_product(r,[-r(2),r(1)])

       call random_number(t0)
       call random_number(dt)
       call random_number(ff)

       m = 2*int(maxval(abs(r)))
       do k = 1, m
          t = real(k) / real(m)
          ! linear Bezier curve
          ! r = (1 - t)*r1 + t*r2
          ! quadratic Bezier curve
          r = (1 - t)**2*r1 + 2*(1 - t)*t*r3 + t**2*r2
          f = fmax*ff*exp(-((t-t0)/dt)**2)
          i = nint(r(1))
          j = nint(r(2))
          if( 1 <= i .and. i <= width .and. 1 <= j .and. j <= height ) then
             sky(i,j) = sky(i,j) + f
          end if
       end do
    end do

  end subroutine cosmics


end module cosmic_ray_events
