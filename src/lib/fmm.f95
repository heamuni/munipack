!
!  Fortran 95+ module interface to selected routines from the book
!
!       "Computer Methods for Mathematical Computations",
!           by Forsythe, Malcolm, and Moler (1977)
!
!  known as FMM (http://www.netlib.org/fmm/index.html).
!
!
!  Copyright © 2015-6 F.Hroch (hroch@physics.muni.cz)
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

module fmm

    implicit none

    interface

       subroutine spline (n, x, y, b, c, d)

         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) ::n
         real(dbl), dimension(n), intent(in) :: x,y
         real(dbl), dimension(n), intent(out) :: b,c,d

       end subroutine spline

       double precision function seval(n, u, x, y, b, c, d)

         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n
         real(dbl), intent(in) :: u
         real(dbl), dimension(n), intent(in) :: x,y,b,c,d

       end function seval


       double precision function fmin(ax,bx,f,tol)

         double precision, intent(in) :: ax,bx,tol

         interface
            function f(x)
              double precision :: f
              double precision, intent(in) :: x
            end function f
         end interface

       end function fmin

       double precision function zeroin(ax,bx,f,tol)

         double precision, intent(in) :: ax,bx,tol

         interface
            function f(x)
              double precision :: f
              double precision, intent(in) :: x
            end function f
         end interface

       end function zeroin

    end interface

end module fmm
