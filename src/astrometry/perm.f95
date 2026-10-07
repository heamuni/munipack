
! perm - creates all possible permutations
!
! $ gfortran -Wall -g -p -fbounds-check -fimplicit-none  -I../lib match.f95  perm.f95
!
! Copyright © 2012-4 F.Hroch (hroch@physics.muni.cz)
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

program perm

  implicit none

  integer, dimension(:), allocatable :: id
  integer :: n,nmatch,nperm

  nperm = 0

  nmatch = 3
  allocate(id(nmatch))

  n = 0
  call tperm

  write(*,*) 'nperm=',nperm

contains

  
  recursive subroutine tperm

    use matcher

    integer :: i


    if( n == size(id) ) then

       nperm = nperm + 1
       write(*,*) id

    else

       do i = 1, size(id)

          if( .not. locate(i,id(1:n)) ) then

             n = n + 1
             id(n) = i
             call tperm
             n = n  - 1
          end if
          
       end do
    end if
    
  end subroutine tperm


function locate(d,id)

  logical :: locate
  integer, intent(in) :: d
  integer, dimension(:), intent(in) :: id
  integer :: i

  locate = .false.
  do i = 1, size(id)
     if( id(i) == d ) then
        locate = .true.
        return
     end if
  end do
 
end function locate


end program perm

