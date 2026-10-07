!
!  Generate artifical sky frame
!
!
!  Copyright © 2016-2026 F.Hroch (hroch@physics.muni.cz)
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


module marker

  use iso_fortran_env
  implicit none

contains

  subroutine watermark(sky,gray)

    use artpicture

    real(REAL64), dimension(:,:), intent(in out) :: sky
    real(REAL64), intent(in) :: gray
    real, dimension(:,:), allocatable :: mask
    integer :: i,j,x,y

    allocate(mask(xdimask,ydimask))
    call genmask(mask)

    x = size(sky,1) - xdimask - 10
    y = size(sky,2) - ydimask - 10
    do i = 1, xdimask
       do j = 1, ydimask
          sky(i+x,j+y) = sky(i+x,j+y) + gray*mask(i,ydimask-j+1)
       end do
    end do

    deallocate(mask)

  end subroutine watermark


end module marker
