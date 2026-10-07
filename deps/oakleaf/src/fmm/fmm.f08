!
!  Fortran 95+ module interface to selected routines from the book
!
!       "Computer Methods for Mathematical Computations",
!           by Forsythe, Malcolm, and Moler (1977)
!
!  known as FMM (http://www.netlib.org/fmm/index.html).
!
!
!  Copyright © 2015-8 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of OakLeaf.
!
!  OakLeaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  OakLeaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with OakLeaf.  If not, see <http://www.gnu.org/licenses/>.
!
!

module fmm

  implicit none
  public

  interface

     subroutine spline (n, x, y, b, c, d)

       use iso_fortran_env

       integer, intent(in) ::n
       real(REAL64), dimension(n), intent(in) :: x,y
       real(REAL64), dimension(n), intent(out) :: b,c,d

     end subroutine spline

     real(REAL64) function seval(n, u, x, y, b, c, d)

       use iso_fortran_env

       integer, intent(in) :: n
       real(REAL64), intent(in) :: u
       real(REAL64), dimension(n), intent(in) :: x,y,b,c,d

     end function seval

     real(REAL64) function fmin(ax,bx,f,tol)

       use iso_fortran_env

       real(REAL64), intent(in) :: ax,bx,tol

       interface
          function f(x)
            use iso_fortran_env
            real(REAL64) :: f
            real(REAL64), intent(in) :: x
          end function f
       end interface

     end function fmin

     real(REAL64) function zeroin(ax,bx,f,tol)

       use iso_fortran_env

       real(REAL64), intent(in) :: ax,bx,tol

       interface
          function f(x)
            use iso_fortran_env
            real(REAL64) :: f
            real(REAL64), intent(in) :: x
          end function f
       end interface

     end function zeroin

  end interface

end module fmm
