!
! Cross-match
!
! Copyright © 2013 F.Hroch (hroch@physics.muni.cz)
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

module crossmatch

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)

contains

  subroutine join(ra1,dec1,ra2,dec2,tol,ftol,idx)

    real(dbl), dimension(:), intent(in) :: ra1, dec1, ra2, dec2
    real(dbl), intent(in) :: tol, ftol
    integer, dimension(:,:), intent(out) :: idx
    integer :: i,j,jmin
    real(dbl) :: r,rmin

    idx = 0
    
    do i = 1,size(ra2)
       rmin = tol
       jmin = 0
       do j = 1,size(ra1)
          r = spmetr(ra1(j),dec1(j),ra2(i),dec2(i))
          if( r < rmin )then
             rmin = r
             jmin = j
          end if
       end do

       if( jmin > 0 ) then
          idx(1,i) = jmin
          idx(2,jmin) = i
       end if

    end do

  end subroutine join

  function spmetr(a1,d1,a2,d2)

    ! angular distance in degrees

    real(dbl), parameter :: rad = 57.295779513082322865_dbl
    real(dbl) :: spmetr
    real(dbl), intent(in) :: a1,d1,a2,d2
    real(dbl) :: da,dd,sind
    
    da = a1 - a2
    dd = d1 - d2
    if( abs(da) > 10.0 .and. abs(dd) > 10.0 ) then

       sind = sin(d1/rad)*sin(d2/rad) + cos(d1/rad)*cos(d2/rad)*cos((a1 - a2)/rad)
       spmetr = rad*acos(min(sind,1.0_dbl))

    else
       ! for small angles, haversine formula is used
       ! http://en.wikipedia.org/wiki/Great-circle_distance
       da = da / rad
       dd = dd / rad
       sind = sqrt(sin(dd/2)**2 + cos(d1/rad)*cos(d2/rad)*sin(da/2)**2)
       spmetr = 2*rad*asin(min(sind,1.0_dbl))

    end if

  end function spmetr



end module crossmatch

