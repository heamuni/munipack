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

program testrfactor

  use oakleaf
  use noise_REAL64
!  use toolbox
  use testsuite_REAL64

  implicit none

  integer, parameter :: rp = selected_real_kind(15)
  character(len=15), parameter :: fmt = '(a,l3)'
  integer :: voltage

  call randomnoise

  call get_environment_variable("VOLTAGE",status=voltage)

  if( .not. (voltage == 0) ) then

     call tester(build)

  else
!     call tester(bmuma)
!     call tester(poisson_mix)
!     call tester(poisson_mix_contamined)

!     call tester(poisson_normal1)
!     call tester(poisson_normal)
!     stop 0
!     call tester(poisson3)
     call single
     call identity
!     call tester(bmuma)
     call tester(poisson1)
     call tester(poisson2)

  end if


contains


  subroutine build

    real(rp), dimension(:), allocatable :: x,dx,y,dy
    type(suite) :: poi
    logical :: cond

    integer :: n,l
    real(rp) :: mean

    mean = 10**6

    n = 10**4
    allocate(x(n),dx(n),y(n),dy(n))

    do l = 1, n
       x(l) = pnoise(mean)
       y(l) = pnoise(mean)
    end do

    write(*,*) 'Po(10**6)  n=',n
    dx = sqrt(x)
    dy = sqrt(y)
    call poi%poistat(x,dx,y,dy)
    cond = abs(poi%t - 1) < 5*poi%dt
    write(*,'(a,l2)') 'TEST build PASS:',cond
    write(*,*)
    deallocate(x,dx,y,dy)

  end subroutine build


  subroutine poisson1

    real(rp), dimension(:), allocatable :: x,dx,y,dy
    type(suite) :: std, rob, qt, poi
    logical :: cond

    integer :: m,n,l
    real(rp) :: mean, fac

    mean = 10**4
    fac = 4

    do m = 1,4
       n = 10**m
       allocate(x(n),dx(n),y(n),dy(n))

       do l = 1, n
          x(l) = pnoise(mean)
          y(l) = pnoise(fac*mean)
       end do

       write(*,*) 'Po(10**6)  n=',n
       call std%stdstat(y/x)
       call rob%robstat(y/x)
       call qt%qtstat(y/x)

       dx = sqrt(x)
       dy = sqrt(y)
!       dx = 1
!       dy = 1
       call poi%poistat(x,dx,y,dy)
       cond = abs(1/poi%t - fac) < 0.1 !5*poi%dt
       write(*,'(a,i0,a,l2)') 'TEST poisson1(',n,') PASS:',cond
       write(*,*)
       deallocate(x,dx,y,dy)
    end do

  end subroutine poisson1

  subroutine poisson2

    real(rp), dimension(:), allocatable :: x,dx,y,dy
    type(suite) :: std, rob, qt, poi
    logical :: cond

    integer :: m,n,l,k
    real(rp) :: mean, fac

    mean = 10**4
    fac = 4

    do m = 1,4
       n = 10**m
       allocate(x(n),dx(n),y(n),dy(n))

       k = int(0.9*n)
       do l = 1, k
          x(l) = pnoise(mean)
          y(l) = pnoise(fac*mean)
       end do
       do l = k+1, n
          x(l) = pnoise(mean)
          y(l) = pnoise(10*mean)
       end do

       write(*,*) '0.9*Po(10000)+0.1*Po(10**5)  n=',n
       call std%stdstat(y/x)
       call rob%robstat(y/x)
       call qt%qtstat(y/x)

       dx = sqrt(x)
       dy = sqrt(y)
       call poi%poistat(x,dx,y,dy)
       cond = abs(1/poi%t - fac) < 0.1 !5*poi%dt
       write(*,'(a,i0,a,l2)') 'TEST poisson2(',n,') PASS:',cond
       write(*,*)
       deallocate(x,dx,y,dy)
    end do

  end subroutine poisson2

  subroutine poisson3

    real(rp), dimension(:), allocatable :: x,dx,y,dy
    type(suite) :: poi
    logical :: cond

    integer :: m,n,l
    real(rp) :: mean, fac

    mean = 8
    fac = 4

    do m = 1,4
       n = 10**m
       allocate(x(n),dx(n),y(n),dy(n))

       do l = 1, n
          x(l) = pnoise(mean)
          y(l) = pnoise(mean/fac)
       end do

       write(*,*) 'Po(8)  n=',n
       dx = 1
       dy = 1
       call poi%poistat(x,dx,y,dy)
       cond = abs(poi%t - fac) < 5*poi%dt
       write(*,'(a,i0,a,l2)') 'TEST poisson3(',n,') PASS:',cond
       write(*,*)
       deallocate(x,dx,y,dy)
    end do

  end subroutine poisson3

!!$  subroutine bmuma
!!$
!!$    real(rp), dimension(6) :: x,dx,y,dy
!!$    type(suite) :: std, rob, qt, poi!, stdw, robw,
!!$    logical :: cond
!!$
!!$    x = [1.58682E+07, 4.53635E+05, 4.41066E+05, 4.14017E+05, &
!!$         4.55078E+05, 3.11341E+05]
!!$    dx = [ 3.990E+03, 712., 704., 682., 713., 608. ]
!!$    y = [ 1.51535E+06, 44100., 42143., 39892., 42909., 30737. ]
!!$    dy = [ 1.233E+03, 221., 216., 209., 218., 189. ]
!!$
!!$    write(*,*) 'BM UMa'
!!$    call std%stdstat(x/y)
!!$    call rob%robstat(x/y)
!!$    call qt%qtstat(x/y)
!!$
!!$!    call stdw%stdwstat(x/y,sqrt(dx**2 + dy**2))
!!$!    call robw%robwstat(x/y,sqrt(dx**2 + dy**2))
!!$    call poi%poistat(y,dy,x,dx)
!!$    cond = abs(poi%t - 10.463) < 5*poi%dt
!!$    write(*,'(a,l2)') 'TEST bmuma PASS:',cond
!!$    write(*,*)
!!$
!!$    block
!!$      real(REAL64) :: t
!!$      call phrate(x,y,t,verbose=.true.)
!!$    end block
!!$
!!$  end subroutine bmuma

  subroutine single

    real, dimension(6) :: x,dx,y,dy
    real :: t,dt,s,c
    integer :: flag

    x = [1.58682E+07, 4.53635E+05, 4.41066E+05, 4.14017E+05, &
         4.55078E+05, 3.11341E+05]
    dx = [ 3.990E+03, 712., 704., 682., 713., 608. ]
    y = [ 1.51535E+06, 44100., 42143., 39892., 42909., 30737. ]
    dy = [ 1.233E+03, 221., 216., 209., 218., 189. ]

    call fmean(x,dx,y,dy,t,dt,s,c,flag=flag)
    write(*,*) 'Single precision: ',t,dt,s,c,flag

  end subroutine single

  subroutine identity

    real(rp), dimension(:), allocatable :: x,dx,y,dy
    type(suite) :: poi
    logical :: cond

    integer :: n,l
    real(rp) :: mean

    mean = 10**2

    n = 10
    allocate(x(n),dx(n),y(n),dy(n))

    do l = 1, n
       x(l) = max(pnoise(mean),1d-3)
    end do
    dx = sqrt(x)
    y = x
    dy = dx

    write(*,*) 'Po(10**2)  n=',n
    call poi%poistat(x,dx,y,dy)
    cond = abs(poi%t - 1) < 5*epsilon(poi%t)
    write(*,'(a,l2)') 'TEST identity PASS:',cond
    write(*,*)
    deallocate(x,dx,y,dy)

  end subroutine identity

  subroutine poisson_normal

    real(rp), dimension(:), allocatable :: x,dx,y,dy
    type(suite) :: std, rob, qt, poi
    logical :: cond

    integer :: m,n,l
    real(rp) :: mean, fac, sig, zero

    mean = 10**4
    fac = 4
    zero = 0.0
    sig = 1e3

    do m = 16,16,2
       n = 2**m
       allocate(x(n),dx(n),y(n),dy(n))

       do l = 1, n
          x(l) = pnoise(mean) + gnoise(zero,sig)
          y(l) = pnoise(fac*mean) + gnoise(zero,sig*fac)
       end do

!       call histogram(y,fac*(mean-5*sig),fac*(mean+5*sig),1000,'/tmp/g')
!       stop 'histo'

       write(*,*) 'Po(1e4)+N(0,1e3)  n=',n
       call std%stdstat(y/x)
       call rob%robstat(y/x)
       call qt%qtstat(y/x)

       call gnoise_fill(zero,sig,dx)
       dx = sqrt(x + dx**2)
!       dx = sqrt(x + gnoise(zero,sig)**2)
       call gnoise_fill(zero,sig,dy)
       !       dy = sqrt(y + fac**2*gnoise(zero,sig)**2)
       dy = sqrt(y + fac**2*dy**2)
       call poi%poistat(x,dx,y,dy)
!       call poi%poistat(y,dy,x,dx)
       cond = abs(poi%t - fac) < 5*poi%dt
       write(*,'(a,i0,a,l2)') 'TEST poisson_normal(',n,') PASS:',cond
       write(*,*)
       deallocate(x,dx,y,dy)

    end do

  end subroutine poisson_normal

  subroutine poisson_normal1

    ! simulates a few stars with magnitudes distributed by exponential law

    real(rp), dimension(:), allocatable :: x,dx,y,dy
    type(suite) :: std, rob, qt, poi
    logical :: cond

    integer :: m,n,l,mg
    real(rp) :: mean, mag, fac, sig, zero, e !, s, r, e, t, r1

    mean = 10**4
    fac = 4
    zero = 0.0
    sig = 2e4

    n = 1000000
    allocate(x(n),dx(n),y(n),dy(n))

    n = 0
    do mg = 10,12
       mag = mg / 1.0
       mean = 10**((25 - mag)/2.5)
!       mag = 25 - 2.5*log10(mean)
       m = max(int(exp(mag*1.1)/1e5),1)
!       m = 1000
!       m = min(int(exp(mag*1.1)/5e4),2)
!       m = 5
       write(*,*) mean,mag,m
       do l = 1, m
          if( n == size(x) ) goto 42
          n = n  + 1
          x(n) = pnoise(mean) + gnoise(zero,sig)
          y(n) = pnoise(fac*mean) + gnoise(zero,sig*fac)
          e = x(n) + gnoise(zero,sig)**2
          if( e > 0 ) then
             dx(n) = sqrt(e)
          else
             dx(n) = 1e-7
          end if
          e = y(n) + fac**2*gnoise(zero,sig)**2
          if( e > 0 ) then
             dy(n) = sqrt(e)
          else
             dy(n) = 1e-7
          end if
          if( .not. (x(n) > 0 .and. y(n) > 0) ) n = n - 1
!               dy(n) = sqrt(y(n) + fac**2*max(gnoise(zero,sig)**2,1d-7))
!          write(*,*) x(n),y(n)
       end do
    end do
42 continue

!    call histogram(y,0.5d6,4.5d6,1000,'/tmp/q')
!    stop 'histo'

!    open(2,file='/tmp/r')
!    s = 0
!    do l = 250, 550
!       t = l / 100.0
!       s = 0
!       r = 0
!       r1 = 0
!       do m = 1, n
!          e = (y(m) - t*x(m)) / sqrt(dy(m)**2 + t**2*dx(m)**2)
!          s = s + itukey(e)
!          r = r + ihuber(e)
!          r1 = r1 + ihuber(e) + log(sqrt(dy(m)**2 + t**2*dx(m)**2))
!!          s = s + ((y(m) - t*x(m)) / sqrt(y(m) + t*x(m)**2))**2
!       end do
!       write(2,*) t, s / n, r / n, r1 / n**2
!    end do
!    close(2)

    write(*,*) 'Po(..by mag...)+N(0,1e3)  n=',n
    call std%stdstat(y(:n)/x(:n))
    call rob%robstat(y(:n)/x(:n))
    call qt%qtstat(y(:n)/x(:n))
    call poi%poistat(x(:n),dx(:n),y(:n),dy(:n))
    cond = abs(poi%t - fac) < 5*poi%dt
    write(*,'(a,i0,a,l2)') 'TEST poisson_normal(',n,') PASS:',cond
    write(*,*)
    deallocate(x,dx,y,dy)

  end subroutine poisson_normal1


end program testrfactor
