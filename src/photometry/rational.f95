!
!  Robust estimate of a rational function
!
!  Copyright © 2015 F.Hroch (hroch@physics.muni.cz)
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

module robustrational

  implicit none

  ! numerical precision of real numbers
  integer, parameter, private :: dbl = selected_real_kind(15)

  ! print debug informations ?
  logical, parameter, private :: verbose = .false.

  ! estimate jacobian by differences or derivations ?
  !logical, parameter, private :: analytic = .false.
  ! currently, only differences are implemented

  ! limit to detect non-convergent series, reccomended values: >100
  real(dbl), parameter, private :: siglim = 1000.0

  real(dbl), dimension(:), allocatable, private :: xdata, ydata
  real(dbl), private :: par1

  private :: ratfun

contains

  subroutine robrat(x,y,p,dp,ab)

    use NelderMead

    real(dbl), dimension(:), intent(in) :: x,y
    real(dbl), intent(in) :: ab
    real(dbl), dimension(:), intent(in out) :: p,dp
    integer :: icount, numres, ifault, ndata
    real(dbl),dimension(size(p)) :: p0
    real(dbl) :: reqmin,rms

    ndata = size(x)
    p0 = p
    par1 = ab
    reqmin = epsilon(reqmin)
    allocate(xdata(ndata),ydata(ndata))
    xdata = x
    ydata = y
    call nelmin(ratfun,2,p0,p,rms,reqmin,dp,1,1000000,icount,numres,ifault)

    deallocate(xdata,ydata)

  end subroutine robrat

  function ratfun(p) result(s)

    real(dbl), dimension(:), intent(in) :: p
    real(dbl) :: s

!    s = sum(abs(p(1)*(xdata + p(2))/(1.0_dbl + p(3)*xdata) - ydata))
    s = sum(abs(p(1)*(xdata + p(2))/(1.0_dbl + par1*xdata) - ydata))

  end function ratfun


end module robustrational
