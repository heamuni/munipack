!
!  Robust estimation of mean photon rates
!
!  Copyright © 2020 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of OakLeaf.
!
!  OakLeaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  OakLeaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!


module photonrate_phiopt

  use photonrate_common
  use iso_fortran_env

  implicit none
  private

  public :: phiopt

  logical :: debug_graph = .false.

  real(REAL64), dimension(:), allocatable :: xpdf, ypdf, res
  real(REAL64), dimension(:), pointer :: y, x, mju
!  real(REAL64), dimension(3) :: xlast, ylast
!  integer :: nlast
  real(REAL64) :: sx, sy, tau, rho, funmin, funmax, scale
!  logical :: second


contains

  subroutine phiopt(sxs,sys,taus,rhos,mjus,cts,pht,s,phi,dphi,singular)

    use fmm
!    use minpacks

    real(REAL64), intent(in) :: sxs, sys, taus, rhos
    real(REAL64), intent(out) :: s, phi, dphi
    real(REAL64), dimension(:), intent(in), target :: mjus, cts,pht
    logical, intent(out) :: singular

    character(len=80) :: errmsg
!    real(REAL64), dimension(3,3) :: mat
!    real(REAL64), dimension(3) :: koff
    real(REAL64) :: phimin, phimax, tol, xphi, yphi
    integer :: n,stat
    logical :: reli

    sx = sxs
    sy = sys
    tau = taus
    rho = rhos
    mju => mjus
    y => cts
    x => pht
    dphi = -1

!    call rphi(tau,rho,cts,pht,mju,sx,sy,s,phi,dphi,reli)
!    if( debug ) write(error_unit,*) 'rphi:',real(RAD*phi),real(RAD*dphi),reli

!    if( reli ) return


    n = size(cts)
    allocate(xpdf(n),ypdf(n),res(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,*) 'Error `',trim(errmsg),"' (phiopt)."
       phi = PI / 4
       return
    end if



    phimin = 0
    phimax = PI / 2
    tol = 1 / RAD
    funmin = huge(funmin)
    funmax = -huge(funmin)
!    xlast = 0
!    ylast = 0
!    nlast = 0

!    second = .false.
!    xphi = fmin(phimin,phimax,cdfmax,tol)
!    second = .true.
!    yphi = fmin(phimin,phimax,cdfmax,tol)

    !    phiopt = (xphi + yphi)  / 2

!    xphi = fmin(phimin,phimax,xsymmetry,tol)
!    yphi = fmin(phimin,phimax,ysymmetry,tol)
!    if( debug ) write(error_unit,*) 'Final phi:', real(RAD*xphi), real(RAD*yphi)
    phi = fmin(phimin,phimax,symmetry,tol)
    write(*,*) funmax, funmin
    singular = funmax / funmin < 1.1
    s = scale
    if( debug ) write(error_unit,*) 'Final phi:', real(RAD*phi),singular,real(scale)

!!$    write(*,*) real(xlast)
!!$    write(*,*) real(ylast)
!!$
!!$    if( nlast >= 3 ) then
!!$       mat(:,1) = xlast**2
!!$       mat(:,2) = xlast
!!$       mat(:,3) = 1
!!$       write(*,*) real(mat(1,:))
!!$       write(*,*) real(mat(2,:))
!!$       write(*,*) real(mat(3,:))
!!$!       call qrinv(mat,mat1)
!!$!       write(*,*) mat1
!!$!       write(*,*) matmul(mat1,ylast)
!!$       call qrsolve(mat,ylast,koff)
!!$       write(*,*) real(koff),minval(ylast)
!!$
!!$       if( koff(1) > 0 ) then
!!$          s0 = minval(ylast)
!!$          dphi = sqrt(s0/koff(1)/(n-1))
!!$          write(*,*) dphi*RAD
!!$       end if
!!$    end if


    if( debug_graph ) then
       block
         use quantilmean

         real(REAL64) :: sigu, sigv
         integer :: i
         open(1,file='/tmp/symmetry')
         do i = 0, 90
            write(1,*) i,symmetry(i/RAD)
         end do
         close(1)


         sigu = rho*sin(phi)
         sigv = rho*cos(phi)/tau

         res = (x - mju) / sqrt(mju/sx + sigv**2)
!         res = x - mju
         call ecdf(res,xpdf,ypdf)
         open(1,file='/tmp/xsymmetry')
         do i = 1, size(xpdf)
            write(1,*) xpdf(i),ypdf(i)
         end do
         close(1)

         res = (y - tau*mju) / sqrt(mju*tau/sy + sigu**2)
!         res = y - tau*mju
         call ecdf(res,xpdf,ypdf)
         open(1,file='/tmp/ysymmetry')
         do i = 1, size(xpdf)
            write(1,*) xpdf(i),ypdf(i)
         end do
         close(1)

       end block
    end if

    deallocate(xpdf, ypdf,res)

  end subroutine phiopt

!!$  real(REAL64) function cdfmax(phi)
!!$
!!$    use quantilmean
!!$
!!$    real(REAL64), intent(in) :: phi
!!$
!!$    real(REAL64), parameter :: q25 = 1 / real(4.0,REAL64)
!!$    real(REAL64), parameter :: q75 = 3 / real(4.0,REAL64)
!!$    real(REAL64) :: sigu, sigv, xmax, ymax, x25, x75
!!$
!!$    sigu = rho*sin(phi)
!!$    sigv = rho*cos(phi)/tau
!!$
!!$    ! pht
!!$    if( .not. second ) then
!!$       sig = sqrt(mju/sx + sigv**2)
!!$       res = (x - mju) / sig
!!$
!!$       call ecdf(res,xpdf,ypdf)
!!$       x25 = quantile(q25,xpdf,ypdf)
!!$       x75 = quantile(q75,xpdf,ypdf)
!!$
!!$       cdfmax = max(Qdev + x25, x75 - Qdev)
!!$    end if
!!$
!!$    ! cts
!!$    if( second ) then
!!$       sig = sqrt(mju*tau/sy + sigu**2)
!!$       res = (y - tau*mju) / sig
!!$
!!$       call ecdf(res,xpdf,ypdf)
!!$       x25 = quantile(q25,xpdf,ypdf)
!!$       x75 = quantile(q75,xpdf,ypdf)
!!$
!!$       cdfmax = max(Qdev + x25, x75 - Qdev)
!!$    end if
!!$
!!$!    cdfmax = max(xmax, ymax)
!!$
!!$    write(*,*) real(RAD*phi), real(cdfmax), real(x25), real(x75)
!!$
!!$  end function cdfmax

  real(REAL64) function symmetry(phi)

    use quantilmean

    real(REAL64), intent(in) :: phi

    real(REAL64), parameter :: q25 = 1 / real(4.0,REAL64)
    real(REAL64), parameter :: q75 = 3 / real(4.0,REAL64)
    real(REAL64) :: sigu, sigv, x25, x75, y25, y75

    sigu = rho*sin(phi)
    sigv = rho*cos(phi)/tau

    res = (x - mju) / sqrt(mju/sx + sigv**2)
!    res = x - mju
    call ecdf(res,xpdf,ypdf)
    x25 = quantile(q25,xpdf,ypdf)
    x75 = quantile(q75,xpdf,ypdf)

    res = (y - tau*mju) / sqrt(mju*tau/sy + sigu**2)
!    res = y - tau*mju
    call ecdf(res,xpdf,ypdf)
    y25 = quantile(q25,xpdf,ypdf)
    y75 = quantile(q75,xpdf,ypdf)

!    symetry = max(abs(y75 - x75),abs(y25 - x25))
!    symmetry = (Qdev + x25)**2 + (x75 - Qdev)**2  &
!             + (Qdev + y25)**2 + (y75 - Qdev)**2
!    symmetry = (abs(Qdev + x25) + abs(x75 - Qdev))/2  &
!             - (abs(Qdev + y25) + abs(y75 - Qdev))/2
    symmetry = abs((x75 - x25)/2 - (y75 - y25)/2)

!    block
!      integer :: i
!      do i = 2, size(ypdf)
!         if( ypdf(i-1) <= 0.25 .and. 0.25 <= ypdf(i) ) then
!            write(*,*) xpdf(i),xpdf(i-1)
!         end if
!      end do
!    end block

    !    write(*,*) real(x25),real(x75),real(y25),real(y75)
!    if( debug_iter ) write(error_unit,*) real(RAD*phi), real(symmetry), &
!         real(abs(x75-y75)), real(abs(x25 - y25))
    if( debug_iter ) write(error_unit,*) real(RAD*phi), real(symmetry), &
         real((x75-x25)/2), real((y75 - y25)/2)

    if( symmetry < funmin ) funmin = symmetry
    if( symmetry > funmax ) funmax = symmetry
    scale = ((x75 - x25) + (y75 - y25))/(2*Qdev)

!    nlast = nlast + 1
!    xlast(1:2) = xlast(2:3)
!    xlast(3) = phi
!    ylast(1:2) = ylast(2:3)
!    ylast(3) = symmetry


!    if( nlast < 3 ) then
!       nlast = nlast + 1
!       last(nlast) = symmetry
!    else
!       last(2:3) = last(1:2)
!       last(1) = symmetry
!    end if

  end function symmetry


  real(REAL64) function xsymmetry(phi)

    use quantilmean

    real(REAL64), intent(in) :: phi

    real(REAL64), parameter :: q25 = 1 / real(4.0,REAL64)
    real(REAL64), parameter :: q75 = 3 / real(4.0,REAL64)
    real(REAL64) :: sigu, sigv, x25, x75, y25, y75

    sigu = rho*sin(phi)
    sigv = rho*cos(phi)/tau

    res = (x - mju) / sqrt(mju/sx + sigv**2)
    call ecdf(res,xpdf,ypdf)
    x25 = quantile(q25,xpdf,ypdf)
    x75 = quantile(q75,xpdf,ypdf)

!    xsymmetry = abs((x75-x25)/(2*Qdev) - 1)
    xsymmetry = abs(x75/Qdev - 1) + abs(x25/Qdev + 1)

    write(*,*) real(RAD*phi),real(xsymmetry),real(x75/Qdev),real(x25/Qdev)

  end function xsymmetry

  real(REAL64) function ysymmetry(phi)

    use quantilmean

    real(REAL64), intent(in) :: phi

    real(REAL64), parameter :: q25 = 1 / real(4.0,REAL64)
    real(REAL64), parameter :: q75 = 3 / real(4.0,REAL64)
    real(REAL64) :: sigu, sigv, x25, x75, y25, y75

    sigu = rho*sin(phi)
    sigv = rho*cos(phi)/tau

    res = (y - tau*mju) / sqrt(mju*tau/sy + sigu**2)
    call ecdf(res,xpdf,ypdf)
    y25 = quantile(q25,xpdf,ypdf)
    y75 = quantile(q75,xpdf,ypdf)
!    ysymmetry = abs((y75-y25)/(2*Qdev) - 1)
    ysymmetry = abs(y75/Qdev - 1) + abs(y25/Qdev + 1)

    write(*,*) real(RAD*phi),real(ysymmetry),real(y75/Qdev),real(y25/Qdev)
  end function ysymmetry


  subroutine rphi(tau,rho,cts,pht,mju,sx,sy,s,phi,dphi,converge)

    use rfun

    real(REAL64), intent(in) :: tau,rho, sx, sy, s
    real(REAL64), dimension(:), intent(in) :: cts,pht,mju
    real(REAL64), intent(in out) :: phi
    real(REAL64), intent(out) :: dphi
    logical, intent(out) :: converge

    real(REAL64) :: tol,d,f,df,cosfi,sinfi,sigv,sigu,rx,ry,sigx,sigy, &
         sigy_f, sigx_f, sigy_ff, sigx_ff, psix, psiy, dpsix, dpsiy, &
         psi,dpsi
    integer :: iter, i, n


    converge = .false.

    n = size(cts)

    if( n == 0 ) stop 'in rphi'

!    phi = PI / 4
    tol = 1 / RAD / sqrt(real(n))
    if( debug_iter ) write(error_unit,*) 'rphi tol:',real(RAD*tol)
    do iter = 1, precision(phi)

       f = 0
       df = 0
       do concurrent (i = 1:n)

          cosfi = cos(phi)
          sinfi = sin(phi)
          sigv = (rho/tau)*cosfi
          sigu = rho*sinfi

          rx = pht(i) - mju(i)
          ry = cts(i) - mju(i)*tau
          sigx = sqrt(mju(i) + sigv**2)
          sigy = sqrt(mju(i)*tau + sigu**2)
          sigy_f = rho*sigu*cosfi/sigy
          sigx_f = -(rho/tau)*sigv*sinfi/sigx

          psix = (rx / sigx) / s
          psiy = (ry / sigy) / s
!         psix = tukey((rx / sigx) / s)
!         psiy = tukey((ry / sigy) / s)
!          psix = huber(rx / sigx / s)
!          psiy = huber(ry / sigy / s)


          psi = (psiy*(ry/sigy)/s - 1)*(sigy_f/sigy) &
              + (psix*(rx/sigx)/s - 1)*(sigx_f/sigx)

          sigy_ff = rho*(rho/sigy)*(cosfi**2 - sinfi**2) &
               - (sigu*cosfi/sigy)*(sigy_f/sigy)
          sigx_ff = (rho/tau)*((rho/sigx)*(sinfi**2 - cosfi**2)) &
               + (rho/sigx)*(sigx_f/sigx)*cosfi*sinfi

          dpsix = 1
          dpsiy = 1
!         dpsix = dtukey((rx / sigx) / s)
!         dpsiy = dtukey((ry / sigy) / s)
!          dpsix = dhuber(rx / sigx / s)
!          dpsiy = dhuber(ry / sigy / s)


          dpsi = -(dpsiy*(ry/sigy)**2*(sigy_f/sigy)/s**2 &
               + psiy*(ry/sigy)*(sigy_f/sigy))/s &
               *(sigy_f/sigy) + (psiy*(ry/sigy)/s - 1)*(sigy_ff/sigy - (sigy_f/sigy)**2) &
               -(dpsix*(rx/sigx)**2*(sigx_f/sigx)/s**2 &
               + psix*(rx/sigx)*(sigx_f/sigx))/s &
               *(sigx_f/sigx) + (psix*(rx/sigx)/s - 1)*(sigx_ff/sigx - (sigx_f/sigx)**2)


          f = f + psi
          df = df + dpsi

       end do
       f = f / n
       df = df / n

       goto 33
       block
         real(REAL64) :: f1, f2, fd, eps
         eps = epsilon(phi)**(1.0/3.0)

         f1 = eloglike(mju,cts,pht,tau,rho,phi+eps,sx,sy)
         f2 = eloglike(mju,cts,pht,tau,rho,phi-eps,sx,sy)
         fd = (f1 - f2)/(2*eps)

         write(*,*) f,fd,f-fd
         write(*,*) f,f+df*eps
!         stop 666

       end block
33     continue

       if( .not. (df < 0) ) exit ! when a local minimum is encountered

       d = f / df
       phi = phi - d

       if( debug_iter ) write(error_unit,'(a,i2,2x,g0,2es12.2,es14.4)') &
            "phi, incr., f, f': ",iter,real(phi),d,f,df

       converge = abs(d) < tol .and. df < 0

       ! exit of iterations: the absolute error must be at least |d| < tol
       if( converge ) exit

       if( phi < 0 .or. phi > PI/2 ) exit

    end do

    if( converge ) then
       dphi = 1 / sqrt(-df*(2*n - 1))
    end if

  end subroutine rphi

  real(REAL64) function eloglike(mju,cts,pht,tau,rho,phi,sx,sy)

    real(REAL64), dimension(:), intent(in) :: mju,cts,pht
    real(REAL64), intent(in) :: tau,rho,phi,sx,sy

    real(REAL64) :: rx,ry,sigu,sigv,sigx,sigy,s
    integer :: n,i

    n = size(cts)
    sigu = rho*sin(phi)
    sigv = (rho/tau)*cos(phi)

    s = 0
    do concurrent (i = 1:n)
       rx = pht(i) - mju(i)
       ry = cts(i) - mju(i)*tau
       sigx = sqrt(mju(i)/sx + sigv**2)
       sigy = sqrt(mju(i)*tau/sy + sigu**2)
!       s = s + itukey(ry/sigy) + itukey(rx/sigx) + log(sigy) + log(sigx)
       s = s + (ry/sigy)**2/2 + (rx/sigx)**2/2 + log(sigy) + log(sigx)
    end do

    eloglike = -s / n

  end function eloglike


end module photonrate_phiopt
