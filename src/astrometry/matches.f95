
! tester - matching
!
! $ gfortran -Wall -g -p -fbounds-check -fimplicit-none match.f95 matches.f95
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

! dummy
module match_print

contains

  subroutine matchprint(id1,id2,p,d)
    integer, parameter :: dbl = selected_real_kind(15)
    integer, dimension(:), intent(in) :: id1,id2
    real(dbl), intent(in) :: p,d
    integer :: n
    character(len=80) :: fmt

    if( size(id1) /= size(id2) ) stop 'astromatch size(id1) /= size(id2)'

    n = size(id2)

    write(fmt,"(a,i0,a,i0,a)") "(",n,"i3,a,",n,"i3,a,f5.3,f8.5)"
!    write(*,fmt) id1,' -> ',id2,' : ',p,d

!    write(fmt,"(a,i0,a)") "(a,i0,",2*n,"(1x,f0.3),f6.3,f8.5)"
!    write(*,fmt) "=M> ",n,(xstar(id2(i)),ystar(id2(i)),i=1,n),p,d
    
  end subroutine matchprint

end module match_print



program matches

  use matcher
  use match_print

  implicit none

  integer, parameter :: dp = selected_real_kind(15)
  integer :: n, np, nmatch, i, j
  integer, parameter :: npmax = 10000
  integer, allocatable, dimension(:)  :: p,pn,id1,id2
  real(dp), allocatable, dimension(:) :: x1,y1,x2,y2
  real(dp) :: q, pmax


  np = 1
  n = 20
  nmatch = 7

  ! random stars
  allocate(x1(n),y1(n),x2(n),y2(n))!,id1(n),id2(n))
  do i = 1, n
     call random_number(q)
     x1(i) = 100.0*q
     call random_number(q)
     y1(i) = 100.0*q
!     id1(i) = i
  end do

  ! permutations
  write(*,*) 'number of permutations (factorial) = ',fact(n)
  allocate(p(n),pn(n))
  do i = 1, n
     p(i) = i
  enddo
!  pn = p
  call per(n)
  
  write(*,*) pn


  ! random transformation
  do i = 1, n
     j = pn(i)
!     id2(i) = i
     x2(i) = x1(j) + gnoise(0.01)
     y2(i) = y1(j) + gnoise(0.01)
     write(*,*) i,j
  end do

  
  call match(x1,y1,x2,y2,nmatch,id1,id2,pmax,matchprint)

  write(*,*) pmax
  write(*,*) id1
  write(*,*) id2

  ! residuals
  do i = 1,size(id1)
!     j = id2(i)
!     write(*,'(4f10.3)') x1(id1(i)),y1(id1(i)),x2(id2(i)),y2(id2(i))
!     write(*,'(4f10.3)') x1(i),y1(i),x2(i),y2(i)
  end do


  deallocate(p,pn)
  deallocate(x1,y1,x2,y2,id1,id2)

contains

   function gnoise(noise)

    real :: gnoise, noise
    real :: x

    call random_number(x)
    gnoise = noise*invdist(x)
    

  end function gnoise


  function invdist(xx)

    real :: invdist,xx

    ! rational approximation of an inverse to a cumulative function
    ! of Gaussian distribution  with precision better than 0.00045
    ! J.Andel: Statistical methods, Matfyz Press, Prague 1991
    
    real :: w,f,x
    logical :: interval

    x = xx
    if( x < 0.0 ) then
       invdist = 0.0
    elseif( x > 1.0 )then
       invdist = 1.0
    else
       interval = x < 0.5
       if( .not. interval ) x = 1.0 - x + epsilon(1.0)
       w = sqrt(-2.0*log(x));
       f = -w + (2.515517 + w*(0.802853 + w*0.010328))/ &
            (1.0 + w*(1.432788 + w*(0.189269 + w*0.001308)));
       if( interval ) then
          invdist = f
       else
          invdist = -f
       endif
    endif
    
  end function invdist



  ! factorial
  recursive function fact(n)  result(ff)

     implicit none

     integer :: ff
     integer, intent(in) :: n

     if( n > 0 ) then
        ff = n*fact(n-1)
     else
        ff = 1
     endif
!     write(*,*) n,ff

   end function fact


   ! all permutations
   recursive subroutine per(level)

     implicit none

     integer, intent(in) :: level
     integer :: i,pp

     if( np > npmax ) return

     if( level > 1 )then
        do i = 1,n
!           if( i /= n ) then
           pp = p(i)
           p(i) = p(1)
           p(1) = pp
           call per(level-1)
           pp = p(i)
           p(i) = p(1)
           p(1) = pp
!        end if
        enddo
     else
!        write(*,'(a,100I3)') 'p:',p
        continue

        np = np + 1
!        if( np == n ) then
        if( np == npmax ) then
           pn = p
!           write(*,'(a,100I3)') 'pn:',pn
        end if
           
     end if

!     write(*,*) 'p:',p
   end subroutine per


end program matches


