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

!
! gfortran -Wall  -fcheck=all  teststat.f95 -L. -L../minpack -lrstat -lsort -llmin -lminpacks -lminpack -lm
!

program testmean

  use oakleaf
  use noise_REAL64
!  use toolbox
  use testsuite_REAL64
  use iso_fortran_env

  implicit none

  integer, parameter :: rp = REAL64
  character(len=15), parameter :: fmt = '(a,l3)'
  integer :: voltage

  call randomnoise

  call get_environment_variable("VOLTAGE",status=voltage)

  if( .not. (voltage == 0) ) then

     call tester(build)

  else

     call tester(identity)
     call tester(patolog1)
     call tester(patolog2)
     call tester(patolog3)
     call tester(patolog7)
     call tester(dark1)
     call tester(dark2)
     call tester(normal01)
     call tester(mixed09)
     call tester(poisson1)
     call tester(poisson2)

  end if


!  call histogram(t,-5.0_rp,5.0_rp,50,'/tmp/h')

!  stop 0
!50 continue




!!$  ! double Normal
!!$  write(*,*)
!!$  do i = 1, size(t)
!!$     call random_number(u)
!!$     if( u > 0.5 ) then
!!$        t(i) = gnoise(1.0_rp,0.1_rp)
!!$     else
!!$        t(i) = gnoise(0.0_rp,0.1_rp)
!!$     end if
!!$  end do
!!$  u = sum(t)/size(t)
!!$  s = sqrt(sum((t-u)**2)/(size(t)-1.0))
!!$  v = s / sqrt(size(t)-0.0)
!!$  write(*,*) "Amean 0.5: N(0,1) + 0.5: N(0,1): ",u,v,s
!!$  call rmean(t,u,v,s)
!!$  write(*,*) "Rmean 0.5: N(0,1) + 0.5: N(0,1): ",u, v,s
!!$
!!$45 continue
!!$  write(*,*) '--- N(0,1) ---'
!!$  do i = 1, size(t)
!!$     call random_number(u)
!!$     t(i) = gnoise(0.0_rp,1.0_rp)
!!$  end do
!!$  u = sum(t)/size(t)
!!$  s = sqrt(sum((t-u)**2)/(size(t)-1.0))
!!$  v = s / sqrt(size(t)-0.0)
!!$  write(*,*) "Amean: N(0,1): ",u,v,s
!!$  call rmean(t,u,v,s)
!!$  write(*,*) "Rmean: N(0,1): ",u,v,s



contains


  subroutine build

    real(rp), dimension(:), allocatable :: t
    type(suite) ::rob
    logical :: cond

    integer :: n

    n = 10**4
    allocate(t(n))
    call gnoise_fill(0.0_rp,1.0_rp,t)
    write(*,*) 'N(0,1) n=',n
    call rob%robstat(t)
    cond = abs(rob%t - 0) < 5/sqrt(real(N))
    write(*,'(a,l2)') 'TEST build PASS:',cond
    write(*,*)

  end subroutine build

  subroutine identity

    real(rp), dimension(100) :: x
    type(suite) :: std, rob, qt
    logical :: cond

    x = 0

    write(*,*) 'Identity..'
    call rob%robstat(x)
    call std%stdstat(x)
    call qt%qtstat(x)
    cond = abs(rob%t-0) < epsilon(x) .and. abs(rob%s) < epsilon(x)
    write(*,fmt) 'TEST identity PASS:',cond

  end subroutine identity

  subroutine patolog1

    real(rp), dimension(16) :: w
    type(suite) :: std, rob, qt
    logical :: cond

    w = [ 10109.000000000000 ,  10109.000000000000 ,   10109.000000000000 , &
       10109.000000000000 ,   10109.000000000000  ,      10109.000000000000  , &
       10109.000000000000 ,   10109.000000000000   ,     10109.000000000000  , &
       729415.00000000000 ,   10109.000000000000  ,      10109.000000000000 , &
       10109.000000000000 ,   10109.000000000000  ,      10109.000000000000 , &
       10109.000000000000 ]

    write(*,*) 'Patolog1..'
    call rob%robstat(w)
    call std%stdstat(w)
    call qt%qtstat(w)
    cond = abs(rob%t-10109) < epsilon(w) .and. abs(rob%s) < epsilon(w)
    write(*,fmt) 'TEST patolog1 PASS:',cond

  end subroutine patolog1

  subroutine patolog2

    real(rp), dimension(9) :: w
    type(suite) :: std, rob, qt
    logical :: cond

    w = [ 0.00000000,       0.00000000,       0.00000000,       0.00000000, &
       0.00000000,       0.00000000,       109.288963,       18.7713432, &
       0.00000000 ]

    write(*,*) 'Patolog2 zeros..'
    call rob%robstat(w)
    call std%stdstat(w)
    call qt%qtstat(w)
    cond = abs(rob%t) < 10*epsilon(rob%dt)
    write(*,fmt) 'TEST patolog2 PASS:',cond

  end subroutine patolog2

  subroutine patolog3

    real(rp), dimension(4) :: w
    type(suite) :: std, rob, qt
    logical :: cond

    w = [2.25176118E-05, 4.98100162E-05, 2.31736030E-05,   9.98998730E-05]

    write(*,*) 'Patolog 3'
    call rob%robstat(w)
    call std%stdstat(w)
    call qt%qtstat(w)
    cond = abs(rob%t-std%t) < 5*sqrt(std%dt**2 + rob%dt**2)
    write(*,fmt) 'TEST patolog3 PASS:',cond

  end subroutine patolog3


  subroutine patolog7

    real(rp), dimension(14) :: t
    type(suite) :: std, rob
    logical :: cond

    t = [1.01148117, 1.00847173, 1.00453901,1.00900149,1.00946915,1.00475299,&
       1.01501632, 1.01952410, 1.02912927, 1.01704729, 1.01167893, 1.00378621,&
       1.00812364,       1.01006317 ]

    write(*,*) 'Patolog7, scaling'
    call rob%robstat(t)
    call std%stdstat(t)
    cond = abs(rob%t - std%t) < 5*sqrt(std%dt**2 + rob%dt**2)
    write(*,fmt) 'TEST patolog7 PASS:',cond

  end subroutine patolog7


  subroutine dark1

    real(rp), dimension(7) :: t = [16, 12, 99,95,18,87,10 ]
    type(suite) :: std, rob
    logical :: cond

    write(*,*) 'dark1'
    call rob%robstat(t)
    call std%stdstat(t)
    cond = abs(rob%t-std%t) < 5*sqrt(std%dt**2 + rob%dt**2)
    write(*,fmt) 'TEST dark1 PASS:',cond

  end subroutine dark1


  subroutine dark2

    real(rp), dimension(7) :: t = [ 142, 141, 149, 149, 142, 148,149 ]
    type(suite) :: std, rob, qt
    logical :: cond

    write(*,*) 'dark2 favorite'
    call rob%robstat(t)
    call std%stdstat(t)
    call qt%qtstat(t)
    cond = abs(rob%t-std%t) < 5*sqrt(std%dt**2 + rob%dt**2)
    write(*,fmt) 'TEST dark2 PASS:',cond

  end subroutine dark2


  subroutine normal01

    real(rp), dimension(:), allocatable :: t
    type(suite) :: std, rob, qt
    logical :: cond

    integer :: m,n

    do m = 1,6
       n = 10**m
       allocate(t(n))
       call gnoise_fill(0.0_rp,1.0_rp,t)
       write(*,*) 'N(0,1) n=',n
       call rob%robstat(t)
       call std%stdstat(t)
       call qt%qtstat(t)
       cond = abs(rob%t - std%t) < 5*sqrt(std%dt**2 + rob%dt**2) .and. &
            abs(rob%t - 0) < 5*rob%dt
       write(*,'(a,i0,a,l2)') 'TEST normal01(',n,') PASS:',cond
       write(*,*)
       deallocate(t)
    end do

  end subroutine normal01

  subroutine mixed09

    real(rp), dimension(:), allocatable :: t
    type(suite) :: std, rob, qt
    logical :: cond

    integer :: m,n,k,l

    do m = 1,6
       n = 10**m
       allocate(t(n))
       call gnoise_fill(0.0_rp,1.0_rp,t)

       k = int(0.9*n)
       do l = 1, k
          t(l) = gnoise(0.0_rp,1.0_rp)
       end do
       do l = k+1, n
          t(l) = gnoise(1.0_rp,10.0_rp)
       end do

       write(*,*) '0.9*N(0,1)+0.1*N(0,10)  n=',n
       call rob%robstat(t)
       call std%stdstat(t)
       call qt%qtstat(t)
       cond = abs(rob%t - 0) < 5/sqrt(real(N))
       write(*,'(a,i0,a,l2)') 'TEST mixed09(',n,') PASS:',cond
       write(*,*)
       deallocate(t)
    end do

  end subroutine mixed09

  subroutine poisson1

    real(rp), dimension(:), allocatable :: t,dt
    type(suite) :: std, rob, qt
    logical :: cond

    integer :: m,n,l
    real(rp) :: mean

    mean = 10**6

    do m = 1,6
       n = 10**m
       allocate(t(n),dt(n))

       do l = 1, n
          t(l) = pnoise(mean)
       end do

       write(*,*) 'Po(10**6)  n=',n
       call rob%robstat(t)
       call std%stdstat(t)
       call qt%qtstat(t)

       cond = abs(rob%t - mean) < 5*rob%dt
       write(*,'(a,i0,a,l2)') 'TEST poisson1(',n,') PASS:',cond

       write(*,*)
       deallocate(t,dt)
    end do

  end subroutine poisson1

  subroutine poisson2

    real(rp), dimension(:), allocatable :: t,dt
    logical :: cond
    type(suite) :: std, rob, qt

    integer :: m,n,k,l
    real(rp) :: mean

    mean = 10**6

    do m = 1,6
       n = 10**m
       allocate(t(n),dt(n))

       k = int(0.9*n)
       do l = 1, k
          t(l) = pnoise(mean)
       end do
       do l = k+1, n
          t(l) = pnoise(10*mean)
       end do

       write(*,*) '0.9*Po(10000)+0.1*Po(10**5)  n=',n
       call rob%robstat(t)
       call std%stdstat(t)
       call qt%qtstat(t)

       cond = abs(rob%t-mean) < 5*rob%dt
       write(*,'(a,i0,a,l2)') 'TEST poisson2(',n,') PASS:',cond

       write(*,*)
       deallocate(t,dt)
    end do

  end subroutine poisson2


end program testmean
