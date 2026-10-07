!
! Matching on base of looking for nearest stars
!
! Copyright © 2017-2026 F.Hroch (hroch@physics.muni.cz)
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

module nearmatch

  use iso_fortran_env
  implicit none

contains

  subroutine listmatch(type,alpha,delta,xpix,ypix,acen,dcen,xcen,ycen, &
       scale,rot,refl,tol,i1,i2,status)

    use astrotrafo

    character(len=*), intent(in) :: type
    real(REAL64), dimension(:), intent(in) :: alpha,delta,xpix,ypix
    real(REAL64), intent(in) :: acen,dcen,xcen,ycen,scale,rot,refl,tol
    integer, dimension(:), allocatable, intent(out) :: i1,i2
    logical, intent(out) :: status
    real(REAL64), dimension(:), allocatable :: x,y,u,v
    integer, dimension(:), allocatable :: id1,id2
    type(AstrometryProjection) :: tp
    type(AstrometryGeometry) :: tg
    real(REAL64) :: smin,smin0,s
    integer :: nsize,msize,n,i,j,idx
    logical :: found

    status = .false.

    if( size(alpha) /= size(delta) .or. size(xpix) /= size(ypix) ) &
         stop 'Listmach arrays does not conforms.'

    ! select maximum amount of stars by list
    nsize = size(alpha)
    msize = size(xpix)
    allocate(u(nsize),v(nsize),id1(nsize),id2(nsize),x(msize),y(msize))

    call tp % init(type,acen,dcen)
    call tp % proj(alpha,delta,u,v)

    call tg % init(scale,rot,refl)
    call tg % pix2uv(xpix-xcen,ypix-ycen,x,y)

    smin0 = tol

    n = 0
    do i = 1, nsize
       smin = smin0
!       write(*,'(i5,5f15.5)') i,u(i),v(i),x(i),y(i),smin
!       write(*,*) alpha(i),delta(i),smin
       idx = 0
       do j = 1, msize
          s = hypot(u(i) - x(j),v(i) - y(j))
          if( s < smin ) then
             ! the nearest and the brighest star (a flux-sorted list is assumed)
             smin = s
             idx = j
          end if
       end do
       if( idx > 0 ) then

          ! check whatever the star is already presented,
          ! it will be ignored, because brighter ones are prefered
          ! (and we are supposing flux-sorted lists).
          found = .false.
          do j = 1, n
             if( idx == id2(j) ) then
                found = .true.
                exit
             end if
          end do

          if( .not. found ) then
             n = n + 1
             id1(n) = i
             id2(n) = idx
          end if

          ! write(*,'(f10.3,2i5,2f15.2)') smin,i,idx,u(i)-x(idx),v(i)-y(idx)

          if( n == size(id1) ) exit
       end if
    end do

    allocate(i1(n),i2(n))
    i1 = id1(1:n)
    i2 = id2(1:n)

    status = n > 0

  end subroutine listmatch


end module nearmatch
