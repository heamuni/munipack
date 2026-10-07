!
!  Robust estimate of a plane
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

module robustplane

  implicit none

  ! numerical precision of real numbers
  integer, parameter, private :: dbl = selected_real_kind(15)

  ! print debug informations ?
  logical, parameter, private :: verbose = .false.

  real(dbl), dimension(:), allocatable, private :: xdata, ydata, zdata, &
       xsig, ysig, zsig
  real(dbl), private :: rsig

  private :: planefun,planedif,planelog,res

contains

  subroutine rplane(x,dx,y,dy,z,dz,q,dq,sig,info)

    use minpacks
    use NelderMead
    use oakleaf, only: median, huber, dhuber

    real(dbl), parameter :: macheps = epsilon(1.0_dbl)
    integer, parameter :: maxeval = 1000000
    real(dbl), dimension(:), intent(in) :: x,y,z,dx,dy,dz
    real(dbl), dimension(:), intent(in out) :: q,dq
    real(dbl), intent(out) :: sig
    integer, intent(out) :: info
    integer :: icount, numres, ifault, ndata, nprint
    real(dbl), dimension(size(q)+1) :: p,dp,p0
    real(dbl), dimension(size(p),size(p)) :: jac,hess
    real(dbl),dimension(:), allocatable :: f,df,r
    real(dbl) :: s,d,mad,sum2,sum3
    integer :: i

    if( size(q) /= 3 ) stop 'Plane needs to have tree parameters.'

    if( verbose ) then
       nprint = 1
    else
       nprint = 0
    end if

    info = 3 ! no initialisation
    ndata = size(x)
    allocate(xdata(ndata),ydata(ndata),zdata(ndata),xsig(ndata),ysig(ndata),zsig(ndata),&
         f(ndata),df(ndata),r(ndata))
    xdata = x
    ydata = y
    zdata = z
    xsig = dx
    ysig = dy
    zsig = dz

    ! init by absolute deviations
    p(1:3) = q
    dp(1:3) = dq
    p0 = p
    call nelmin(planefun,3,p0(1:3),p(1:3),mad,macheps,dp(1:3),1,maxeval,icount,numres,ifault)
    if( verbose ) write(*,'(a,3f15.3,i5)') 'rplane init:',p(1:3),ifault
    if( ifault /= 0 ) then
       info = 2 ! no convergence
       if( verbose ) &
            write(*,*) "rplane init finished prematurely without convergence."
       goto 666
    end if
    q = p(1:3)

    r = abs(zdata - (p(1) + p(2)*xdata + p(3)*ydata))
    sig = median(r) / 0.6745
    if( verbose ) write(*,'(a,f0.5)') 'sig0 = ',sig
    rsig = sig

    if( verbose ) write(*,*) 'Winsorisation:'
    do i = 1,size(zdata)
       r(i) = (zdata(i) - (p(1) + p(2)*xdata(i) + p(3)*ydata(i))) / sig
       if( abs(r(i)) > 3 ) then
          d = p(1) + p(2)*xdata(i) + p(3)*ydata(i) + sign(3*sig,r(i))
          if( verbose ) write(*,'(4g15.5)') i,d,zdata(i),r(i)
          zdata(i) = d
       end if
    end do

    ! locate proper minimum
    p(4) = 1
    p0 = p
    dp = abs(0.1*p + 0.01*sig + 1e-3)
    dp(4) = 0.1
    dp = 0.1*dp
    call nelmin(planelog,size(p),p0,p,s,macheps,dp,1,maxeval,icount,numres,ifault)
    if( verbose ) write(*,'(a,i3,4g13.5)') 'Approximate solution: ',ifault,p
    if( ifault /= 0 ) then
       info = 2 ! no convergence
       if( verbose ) &
            write(*,*) "rplane finished prematurely without likelihood convergence."
       goto 666
    end if
    q = p(1:3)

    ! robust estimate
    call lmdif2(planedif,p,epsilon(p),nprint,info)
    if( verbose ) write(*,'(a,4f10.5,i5)') 'robust:',p,info
    if( info == 5 ) then
       info = 2 ! no convergence
       if( verbose ) &
            write(*,*) "rplane finished prematurely without lmdif convergence."
       goto 666
    end if
    q = p(1:3)

    ! estimation of uncertainties
    call res(p(1),p(2),p(3),r)
    f = huber(r)
    df = dhuber(r)
    sum2 = sum(df)
    sum3 = sum(f**2)

    if( sum2 > epsilon(sum2) ) then

       ! Huber (6.6)
       ! compute jac!!
       call qrinv(jac,hess)
       sig = sig*sqrt(sum3/sum2*ndata/(ndata-1.0))
       do i = 1, size(q)
          if( abs(hess(i,i)) > 0 ) then
             dq(i) = sig * sqrt(abs(hess(i,i)))
          else
             dq(i) = sig / sqrt(ndata - 1.0)
          end if
       end do

       if( verbose ) then
          write(*,*) 'jac:',real(jac(1,:))
          write(*,*) 'jac:',real(jac(2,:))
          write(*,*) 'jac:',real(jac(3,:))
          write(*,*) 'jac:',real(jac(4,:))
          write(*,*) 'hess:',real(hess(1,:))
          write(*,*) 'hess:',real(hess(2,:))
          write(*,*) 'hess:',real(hess(3,:))
          write(*,*) 'hess:',real(hess(4,:))
       end if

    else
       dq = 0
    end if
    info = 0

666 continue

    deallocate(xdata,ydata,zdata,xsig,ysig,zsig,f,df,r)

  end subroutine rplane


  function planefun(p) result(s)

    real(dbl), dimension(:), intent(in) :: p
    real(dbl) :: s

    s = sum(abs(p(1) + p(2)*xdata + p(3)*ydata - zdata)) / size(xdata)

  end function planefun


  subroutine res(a,b,c,r)

    real(dbl), intent(in) :: a,b,c
    real(dbl), dimension(:), intent(out) :: r

    r = (zdata - (a + b*xdata + c*ydata)) / &
         sqrt(zsig**2 + b**2*xsig**2 + c**2*ysig**2 + rsig**2)

  end subroutine res


  function planelog(p)

    use oakleaf

    real(dbl), dimension(:), intent(in) :: p
    real(dbl) :: planelog
    real(dbl), dimension(:), allocatable :: r,f,ds
    real(dbl) :: a,b,c,s
    integer :: n

    n = size(xdata)
    a = p(1)
    b = p(2)
    c = p(3)
    s = p(4)

    if( s < epsilon(s) ) then
       planelog = 100*n
       return
    end if

    allocate(r(n),f(n),ds(n))

    call res(a,b,c,r)
    r = r / s
    f = ihuber(r)

    planelog = sum(f) + n*log(s)
!    write(*,'(7f10.5)') p,planelog,sum(f),n*log(s)

    deallocate(r,f,ds)

  end function planelog



  subroutine planedif(m,np,p,fvec,iflag)

    use oakleaf

    integer, intent(in) :: m,np
    integer, intent(in out) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(out) :: fvec
    real(dbl), dimension(:), allocatable :: r,f,ds
    integer :: n
    real(dbl) :: a,b,c,s

    if( iflag == 0 ) then
!       write(*,'(6g13.5)') p
!       write(*,'(6g13.5)') fvec
       return
    end if

    n = size(xdata)
    allocate(r(n),f(n),ds(n))

    a = p(1)
    b = p(2)
    c = p(3)
    s = p(4)

    ds = sqrt(zsig**2 + b**2*xsig**2 + c**2*ysig**2 + rsig**2)
    call res(a,b,c,r)
    r = r / s
    f = huber(r)

    fvec(1) = sum(f/ds)
    fvec(2) = sum(f*(xdata + r*b*xsig**2/ds)/ds)
    fvec(3) = sum(f*(ydata + r*c*ysig**2/ds)/ds)
    fvec(4) = sum(f*r)/s - n

    fvec = - fvec / s

    deallocate(r,f,ds)

  end subroutine planedif


end module robustplane
