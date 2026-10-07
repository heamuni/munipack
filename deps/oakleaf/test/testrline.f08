!
!  Copyright © 2025 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of OakLeaf.
!
!  OakLeaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  OakLeaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with OakLeaf.  If not, see <http://www.gnu.org/licenses/>.
!


! gfortran -fcheck=all -g -p  -Wall testrline.f95 -L. -L../minpack -lnoise -lrstat -lsort -llmin -lminpacks -lminpack -lfmm -lm

program testrline

  use oakleaf
  use noise_REAL64
  use toolbox
  use testsuite_REAL64

  implicit none
  integer, parameter :: dbl = selected_real_kind(15)
  integer :: voltage

  call get_environment_variable("VOLTAGE",status=voltage)

  if( .not. (voltage == 0) ) then
     call tester(rektoris)
  else
     call randomnoise

     call tester(rektoris)
     call tester(normal)
     call tester(reghuber)
  end if

contains

  subroutine normal

    integer, parameter :: nmax = 70
    real(dbl), dimension(:),allocatable :: x,y,dx,dy
    real(dbl) :: a,b,da,db,sig,t
    integer :: i
    logical :: reliable, cond

    allocate(x(nmax),y(nmax),dx(nmax),dy(nmax))
    do i = 1,size(x)

       call random_number(t)
       x(i) = t
       if( mod(i,10) > 0 ) then
          y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,0.1_dbl)
       else
          y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,1.0_dbl)
       end if

    end do
    i = size(x)/3
    !  y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,0.1_dbl)
    i = size(x)/2
    !  y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,0.1_dbl)
    i = size(x)/4
    !  y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,0.1_dbl)

    dx = 0.1 !* 0.707
    dy = 0.1 !* 0.707
    dx = 0
    dy = 1

    call rline(x,y,a,b,da,db,dx,dy,sig,reliable,verbose=.true.)

    write(*,*) 'a=',a,'+-',da
    write(*,*) 'b=',b,'+-',db
    write(*,*) 'sig=',sig,' reliable:',reliable
    cond = reliable .and. abs(a-1) < 5*0.1 .and. abs(b-1) < 5*0.1
    write(*,'(a,l2)') 'TEST normal PASS:',cond

!    open(1,file='/tmp/rline2')
!    do i = 1, size(x)
!       write(1,*) x(i),y(i),y(i) - (a + b*x(i))
!    end do
!    close(1)

    deallocate(x,y,dx,dy)


  end subroutine normal


  subroutine rektoris

    integer, parameter :: ndat = 14
    real(dbl), dimension(ndat) :: x,y
    real(dbl) :: a,b,da,db,sig
    logical :: reliable, cond

    x = [ 3,4,6,6,6,7,8,9,11,11,12,12,14,16 ]
    y = [ 24.82, 23.26, 14.77, 19.06, 14.79, 17.66, 11.83, &
          14.82, 5.16, 11.12, 8.04, 2.72, 0.74, 1.21 ]

    call rline(x,y,a,b,da,db,sigma=sig,reliable=reliable,verbose=.true.)

    write(*,*) "rline:", reliable
    write(*,*) "sig:",sig
    write(*,*) "a:",a,"+-",da
    write(*,*) "b:",b,"+-",db
    cond = reliable .and. abs(a-29.3) < 0.1 .and. abs(b+1.9) < 0.1
    write(*,'(a,l2)') 'TEST rektoris PASS:',cond

  end subroutine rektoris

  subroutine reghuber

    ! The example data widely discussed in Regression chapter
    ! of Huber's monograph. The data has been generated as
    ! the line: f(x) = -2 - x + N(0,0.6) for first five points.
    ! The last point has value intentionally set to zero to simulate
    ! a gross error.

    integer, parameter :: ndat = 6
    real(dbl), dimension(ndat) :: x,y
    real(dbl) :: a,b,da,db,sig, tol
    logical :: reliable, cond

    x = [ -4, -3, -2, -1, 0, 10 ]
    y = [ 2.48, 0.73, -0.04, -1.44, -1.32, 0.0 ]

    a = -2
    b = -1
    call rline(x,y,a,b,da,db,sigma=sig,reliable=reliable,verbose=.true.)

    write(*,*) "rline:", reliable
    write(*,*) "sig:",sig
    write(*,*) "a:",a,"+-",da
    write(*,*) "b:",b,"+-",db
    tol = 0.6 / 2 ! = 0.6/ sqrt(6 - 2)
    cond = reliable .and. abs(a+2) < tol*abs(a) .and. abs(b+1) < tol
    write(*,'(a,l2)') 'TEST huber PASS:',cond

  end subroutine reghuber


end program testrline
