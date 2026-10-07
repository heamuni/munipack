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

module testsuite_REAL64

  use oakleaf
  use iso_fortran_env

  implicit none

  ! testing suite

  integer, parameter, private :: rp = REAL64

  type :: suite
     integer :: n
     real(rp) :: t,dt,s,sig
   contains
     procedure :: stdstat, robstat, qtstat,poistat!,phstat
  end type suite

contains

  subroutine stdstat(this,x)

    class(suite) :: this
    real(rp), dimension(:), intent(in) :: x

    associate( t => this%t, dt => this%dt, s => this%s, sig => this%sig, &
         n => this%n)

      n = size(x)
      t = sum(x) / n
      s = sqrt(sum( (x-t)**2 / (n-1) ))
      sig = s
      dt = s / sqrt(real(n))
      write(*,*) "Amean:",t,dt,sig,s
    end associate

  end subroutine stdstat

  subroutine robstat(this,x)

    real(rp), dimension(:), intent(in) :: x
    class(suite) :: this

    integer :: iflag
    real :: t1,t2

    associate( t => this%t, dt => this%dt, s => this%s, sig => this%sig, &
         n => this%n)
      call cpu_time(t1)
      call rmean(x,t,dt,sig,s,verbose=.true.,flag=iflag)
      call cpu_time(t2)

      write(*,*) "Rmean:",t,dt,sig,s,iflag,nint(1e3*(t2-t1)),"ms"!
    end associate

  end subroutine robstat

  subroutine qtstat(this,x)

    real(rp), dimension(:), intent(in) :: x
    class(suite) :: this

    associate( t => this%t, dt => this%dt, s => this%s, n => this%n)
      call qmean(x,t,s)
      write(*,*) "Qmean:",t,s
    end associate

  end subroutine qtstat

  subroutine poistat(this,x,dx,y,dy)

    real(rp), dimension(:), intent(in) :: x,dx,y,dy
    class(suite) :: this

    logical :: reli
    real :: t1,t2

    associate( t => this%t, dt => this%dt, s => this%s, n => this%n)
      call cpu_time(t1)
      block
        ! just temporary solution
        integer, parameter :: rpkind = selected_real_kind(13)
        real(rpkind) :: tt,dtt,ss
        real(rpkind), dimension(:), allocatable :: xx, dxx, yy, dyy
        integer :: N

        N = size(x)
        allocate(xx(N),yy(N),dxx(N),dyy(N))
        xx = real(x,rpkind)
        dxx = real(dx,rpkind)
        yy = real(y,rpkind)
        dyy = real(dy,rpkind)
        call fmean(xx,dxx,yy,dyy,tt,dtt,ss,reliable=reli,verbose=.true.)
        t = real(tt,rp)
        dt = real(dtt,rp)
        s = real(ss,rp)
      end block
!      call fmean(x,dx,y,dy,t,dt,s,reliable=reli,verbose=.true.)
      call cpu_time(t2)

      write(*,*) "fmean:",t,dt,s,reli,nint(1e3*(t2-t1)),"ms"
    end associate

  end subroutine poistat

!!$  subroutine phstat(this,x,y)
!!$
!!$    real(rp), dimension(:), intent(in) :: x,y
!!$    class(suite) :: this
!!$
!!$    logical :: reli
!!$    real :: t1,t2
!!$
!!$    associate( t => this%t, dt => this%dt, n => this%n)
!!$      call cpu_time(t1)
!!$      call phrate(x,y,t,dt,reliable=reli,verbose=.true.)
!!$      call cpu_time(t2)
!!$
!!$      write(*,*) "phrate:",t,dt,reli,nint(1e3*(t2-t1)),"ms"
!!$    end associate
!!$
!!$  end subroutine phstat


  subroutine tester(proc)

    interface
       subroutine proc
       end subroutine proc
    end interface

    write(*,*) "Starting .."
    call proc
    write(*,*)
    write(*,*)

  end subroutine tester


end module testsuite_REAL64
