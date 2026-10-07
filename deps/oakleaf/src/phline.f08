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

! NOTE: neslo by pakove body resit obecneji pomoci ortogonalnich polynomu?



module photonrate_line

  use photonrate_common
  use iso_fortran_env

  implicit none
  private

  public :: phline

  logical, parameter :: debug_graph = .true.
  logical, parameter :: debug_test = .false.

 ! Does to use Jacobian analytic, or numerical?
  logical, parameter :: analytic = .true.

  ! see phopt.f08
  integer, parameter :: nprint = 1

  ! to winsorise, limits
  real(REAL64) :: rmin =-1.618
  real(REAL64) :: rmax = 1.618
  real(REAL64) :: zmin =-1.618
  real(REAL64) :: zmax = 1.618

  real(REAL64), parameter :: sqrt2 = sqrt(real(2.0,REAL64))

  real(REAL64), dimension(:), pointer :: kmean
  real(REAL64), dimension(:), allocatable :: cts, pht, mju, xold
  real(REAL64), dimension(:), allocatable, target :: backd
  integer, dimension(:), pointer :: iset
!  real(REAL64), dimension(3) :: lambda
  real(REAL64) :: nu, lambda, scale, rhomax, tmin, tmax, backsx, backsy, backtau, &
       backrho, beta, &
       tauback, backsigu, backsigv, siguback,sigvback ! to be removed
!       tauback, siguback,sigvback, backs, backsx, backsy, &
  !       backsigu, backsigv
  integer :: ndata = 0
  integer :: nfev

  logical :: robust = .false.

contains

  subroutine phline(y,x,nset,isets,sxs,sys, &
       tau,dtau,rho,phi,sigu,dsigu,sigv,dsigv,kmeans,mjus,robustness)

    use rfun
    use medians

    real(REAL64), dimension(:), target, intent(in) :: x,y, kmeans
    real(REAL64), dimension(:), intent(out) :: mjus
    integer, intent(in) :: nset
    integer, dimension(:), target, intent(in) :: isets
    real(REAL64), intent(in) :: sxs,sys
    real(REAL64), intent(in out) :: tau, dtau, sigu, sigv, rho, phi
    real(REAL64), intent(out) :: dsigu, dsigv
    logical, intent(in) :: robustness

    logical, dimension(:), allocatable :: valid
    real(REAL64), dimension(:), allocatable :: delta
    integer, dimension(:), allocatable :: idx
!    real(REAL64), dimension(3) :: c
    real(REAL64), dimension(2) :: sc
    real(REAL64) :: c,stol,s,d,sx,sy,dphi,drho,rho2,tau0,rho0, &
         eta,omega,eta_tol, omega_tol,sig2

    integer :: iter,i,n
    logical :: converge, reli, xreli, yreli

    if( size(x) /= size(y) .or. size(x) < 2 ) &
         stop  'An improper input in phline: size(x) /= size(y) .or. size < 2'

    if( size(x) == 2 ) then
       stop 'Please, implement me.'
    end if


    ndata = size(y)
!    cts => y
!    pht => x
!    mju => mjus
    kmean => kmeans
    iset => isets
!    sx = sxs
!    sy = sys
!    robust = robustness
!    rhomax = 10*sqrt(sigu**2 + tau**2*sigv**2 + 1)
    rhomax = 10*(rho + 1)
    if( rhomax < 1e4 ) rhomax = 5e4
    tmin = 0.1
    tmin = max(tau-10*dtau,tmin)
    tmax = tau + 10*dtau
    tau0 = tau
    rho0 = rho

    robust = robustness


!    write(*,*) tau,tmin,tmax,rhomax
!    stop

    if( robust ) then
!       scale = stukey
       rmin = -2.5
       rmax = 2.5
       !       scale = 1
       ! winsorization setup
       rmin = -1.618
       rmax =  1.618
       scale = 0.9068012
       scale = sqrt(shuber)
       beta = sqrt(shuber)
       beta = shuber
       !       beta = sqrt(stukey)
!       beta = stukey
!       beta = 0.75

       zmin = -1.618
       zmax =  1.618

       ! filter out the out of scale range values
       block
         real(REAL64) :: d,e,mu,sigx, sigy

         allocate(valid(size(y))) ! check
         allocate(idx(size(y)))

         n = 0
         do i = 1, size(y)
            mu =  (y(i) / tau + x(i)) / 2
            d = mu - x(i)
            e = y(i) - tau*(x(i) + d)
            !sigy = sqrt(y(i)*(1 + tau)*tau + rho**2)
            !sigx = sigy / tau
            sigx = sqrt(x(i)*(1 + 1/tau) + (rho/tau)**2) / 2
            sigy = sqrt(y(i)*(1 + tau) + rho**2) / 2
            valid(i) = abs(d/sigx) < 3*zmax .or. abs(e/sigy) < 3*zmax
            valid(i) = .true.
            if( valid(i) ) then
               n = n + 1
               idx(n) = i
            end if
            if( .not. valid(i) .and. .true. ) then
               write(error_unit,'(a,i6,1x,2(g0,2x),2x,2(f10.3,1x))') &
                    'Out:',i,real(x(i)),real(y(i)),d/sigx,e/sigy
            end if
         end do

!         n = count(valid)
         allocate(cts(n),pht(n),mju(n),delta(n)) ! check
         n = 0
         do i = 1, size(valid)
            if( valid(i) ) then
               n = n + 1
               cts(n) = y(i)
               pht(n) = x(i)
            end if
         end do
!         cts = pack(y,valid)
!         pht = pack(x,valid)
!         cts = y(idx(1:n))
!         pht = x(idx(1:n))
       end block
    else
       rmin = -huge(rmin)
       rmax = huge(rmax)
       zmin = -huge(zmin)
       zmax =  huge(zmax)
       scale = 1
       beta = 1
       allocate(cts,source=y) ! check!!!!
       allocate(pht,source=x)
       allocate(valid(size(mjus)))
       valid = .true.
       n = size(cts)
       allocate(mju(n),delta(n),idx(n)) ! check
       idx = [ (i,i=1,n) ]
    end if
    mju = (cts/tau + pht) / 2
    delta = mju - pht
    !    scale = 1

    if( debug_test ) then
       robust = .false.
       tau = 2.02
       sigu = 500
       sigv = 1000
       rho = sqrt(sigu**2 + tau**2*sigv**2)
!       phi = atan2(sigu,tau*sigv)
       rho = 2e3
       sx = 1
       sy = 1
       s = 1.0
       lambda = 1
       nu = 1
       delta = 0
!       delta(1:5) = [ 0.01*x(1), -0.02*x(2), -0.01*x(3), -0.0*x(4), 0.02*x(5) ]
!       delta(1:3) = [ 0.01*sqrt(x(1)), -0.02*sqrt(x(2)), -0.01*sqrt(x(3)) ]
       delta = 1*sqrt(pht)
!       delta = [ 0.01*sqrt(x(1)) ]
       call sigtest(tau,rho)
       call linetest(tau,rho,lambda,mju)
!       call disptest(tau,rho,delta)
!       call linertest(tau,rho)
!       call linetest(tau,rho,delta)
!       call cscaletest(tau,rho,delta,sx,sy)
!       call line1test(tau,rho,s,delta)
!       call lineftest(tau,rho,phi,sx,sy,delta)
!       call graph_1d(tau,rho,s,delta,2,1d2)
!       call xscaletest(pht,mju,tau,sigv,sx,lambda(1),nu)
!       call yscaletest(cts,mju,tau,sigu,sy,lambda(2),nu)
       stop 'All tests.'
    end if

!    su = scale+0.01
!    sv = scale-0.01

!    rho = max(rho,real(1.0,REAL64))

    ! initial estimates
    goto 60
    block
      use fmm

      real(REAL64) :: smin, smax, tol
      smin = 0
      smax = 10*rho
      smax = 1e6
      tol = 1+0.01*rho
!      smax = PI/2
      !      tol = 0.01
      lambda = 0
      nu = 0
      tauback = tau

      if( debug_iter ) write(error_unit,*) &
           'scale on the interval: ',real(smin),real(smax)
      rho = fmin(smin,smax,xscalemin,tol)
      reli = sigv < smax - 2*tol
      write(*,*) 'min MAD(x)=',real(rho),reli

!      sigv = rho/tau/sqrt(2.0)
!      if( debug_iter ) then
!         sx = median(abs(pht - mju)/sqrt(mju + sigv**2)) / Qdev
!!         sx = median(abs(pht - cts/tau)/sqrt(mju*(1+1/tau) + (sigv/tau)**2)) / Qdev
!         write(error_unit,'(a,g0,4x,a,f0.4,l2)') &
!              'MAD(sigv)=',real(sigv),'sx=',sx,reli
!      end if

!    if( debug_iter ) write(error_unit,*) &
!         'xscale on the interval: ',real(smin),real(smax)
!      tol = 1+0.01*sigv
      rho = fmin(smin,smax,yscalemin,tol)
      reli = sigu < smax - 2*tol
      write(*,*) 'min MAD(y)=',real(rho),reli

!      sigu = rho/sqrt(2.0)
!      if( debug_iter ) then
!         sy = median(abs(cts - mju*tau)/sqrt(mju*tau + sigu**2)) / Qdev
!         write(error_unit,'(a,g0,4x,a,f0.4,l2)') &
!              'MAD(sigu)=',real(sigu),'sy=',sy,reli
!      end if

      tauback = tau
!      rho = 5e3
!      delta = mju - x
!      call graph_xscale('/tmp/sv',sigv)
!      call graph_yscale('/tmp/su',sigu)
!      call graph_ecdf(tau,rho,phi,cts,pht,delta)

    end block
60  continue


!    do i = 1, size(cts)
!       write(*,*) i, &
!            pht(i)/sqrt(pht(i)*(1+1/tau)+(rho/tau)**2), &
!            (cts(i)-tau*pht(i))/sqrt(cts(i)*(tau+1) + rho**2)
!            (pht(i) - mju(i))/sqrt(mju(i)*(1+1/tau) + (rho/tau)**2), &
!            (cts(i) - mju(i)*tau)/sqrt(mju(i)*tau*(tau+1) + rho**2)
!    end do

!    stop


!    su = median(abs(cts - mju*tau)/sqrt(mju*tau + sigu**2))
!    sv = median(abs(pht - mju)/sqrt(mju + sigv**2))

!    write(error_unit,*) 'Initial scales:',real(su),real(sv)


!    call xscale1(sigv,dsigv,sv,reli)
!    call xscale2(sigv,dsigv,sv,reli)
!    call yscale1(tau,sigu,dsigu,su,reli)

!    tauback = tau
!    call xscale(sigv,dsigv,sv,xreli)
!    call yscale(sigu,dsigu,su,yreli)


!    lambda = 0
!    nu = 0
!    stop 0


!!$    tauback = tau
!!$    su = median(abs(cts - mju*tau)/sqrt(mju*tau + sigu**2)) /  Qdev
!!$    sv = median(abs(pht - mju)/sqrt(mju + sigv**2)) /  Qdev
!!$    write(*,*) '1111111111',su,sv
!!$    call xscale(sigv,dsigv,sv,reli)
!!$    call yscale(sigu,dsigu,su,reli)
!!$    write(*,*) '222222222222',su,sv,sigu,sigv
!!$    su = median(abs(cts - mju*tau)/sqrt(mju*tau + sigu**2)) /  Qdev
!!$    sv = median(abs(pht - mju)/sqrt(mju + sigv**2)) /  Qdev
!!$    write(*,*) '33333333333',su,sv,sigu,sigv
!!$    call xscale(sigv,dsigv,sv,reli)
!!$    call yscale(sigu,dsigu,su,reli)
!!$    write(*,*) '44444444',su,sv,sigu,sigv
!!$    su = median(abs(cts - mju*tau)/sqrt(mju*tau + sigu**2)) /  Qdev
!!$    sv = median(abs(pht - mju)/sqrt(mju + sigv**2)) /  Qdev
!!$    write(*,*) '5555555',su,sv,sigu,sigv
!!$    call xscale(sigv,dsigv,sv,reli)
!!$    call yscale(sigu,dsigu,su,reli)
!!$    write(*,*) '666666',su,sv,sigu,sigv
!!$    stop


!    sigv = 1e3
!    sigu = 5e2
!    sv = 1
    !    su = 1
!!$    block
!!$      real(REAL64) :: alpha
!!$      alpha = atan(tau)
!!$      write(*,*) rad*alpha,sigu/cos(alpha),sigv/sin(alpha)
!!$      write(*,*) sigv*tau
!!$    end block
!!$    rho = sqrt(sigu**2 + tau**2*sigv**2)
!!$    phi = atan2(sigu,tau*sigv)
!!$    !    phi = atan2(tau*sigv,sigu)
!!$!    rho = 2e3
!!$    !   phi = 14/RAD
!!$    rho = sqrt(sigu**2 + sigv**2)
!!$    rho = sigv!*tau
!!$    phi = atan2(sigv,sigu)
!!$    write(error_unit,*) "rho,phi",real(rho),real(RAD*phi)
!!$    delta = (pht - cts/tau)/2 / sqrt(mju*tau*(1+tau) + rho**2)
!!$    delta = mju - x
!    delta = 0

!    delta = 0
!    call graph_rhos(tau,mju,cts,pht,rho,delta)


!    sv = median(abs(pht - cts/tau)/sqrt(mju*(1+1/tau)+(rho/tau)**2))/Qdev
!    sv = 1
!    su = median(abs(cts - tau*mju)/sqrt(mju*tau*(1+tau)+rho**2))/Qdev
!    su = 1
!    sv = 1

!    sx = median(abs(pht - cts/tau)/sqrt((pht)*(1+1/tau)+(rho/tau)**2))/Qdev
!    sy = median((1+tau)*abs(cts - tau*(pht+delta))/sqrt((pht+delta)*tau*(1+tau)+rho**2))/Qdev
!    write(*,*) 'su,sv:',su,sv

!    rho = sigv
!    sv = 1

!    sx = 1
!    sy = 1

!    rho2 = sigu**2 + tau**2*sigv**2
!    if( rho2 > 0 ) then
!       rho = sqrt(rho2)
!    else
!       rho = 0
!    end if
    if( sigu > 0 .and. sigv > 0 ) then
       phi = atan2(sigu,tau*sigv)
    else
       phi = 0 ! it also implicates rho==0, so we can leave it undefined
    end if
!    rho = 2e3
!    phi = 1/tau
!    sigv = rho/tau*cos(phi)
    !    sigu = rho*sin(phi)


!!$    do i = 1, -size(cts)
!!$
!!$       d = abs(cts(i) - tau*pht(i))/sqrt(pht(i)*tau*(tau+1) + rho**2)
!!$       if( d < rmax ) then
!!$          mju(i) = (cts(i)/tau + pht(i))/2
!!$          delta(i) = mju(i) - pht(i)
!!$       else
!!$          mju(i) = pht(i)
!!$          delta(i) = 0
!!$       end if
!!$!       write(*,*) real(cts(i)),real(pht(i)),real(d),real(mju(i)),nint(delta(i))
!!$    end do
!    delta = mju - pht
!    delta = 0

!!$    do i = 1, size(cts)
!!$       d = abs(cts(i) - tau*pht(i))/sqrt(pht(i)*tau*(tau+1) + rho**2)
!!$       write(*,*) real(cts(i)),real(pht(i)),real(d),real(mju(i)),nint(delta(i))
!!$    end do


!    sx = 2 * median(abs(pht - mju)/sqrt(mju*(1+1/tau) + (rho/tau)**2)) / Qdev
!    sy = 2 * median(abs(cts - mju*tau)/sqrt(mju*tau*(tau+1) + rho**2)) / Qdev

!    delta = 0
!    call thescale(tau,delta,sx,sy,rho,reli)
!    write(*,*) 'thescale()',real(rho),reli
!    stop



!    sx = xdisp(tau,delta,sigv,sx)
!    sy = ydisp(tau,delta,sigu,sy)
!    c(2:3) = scale - [ sx, sy ]

!    call thescale(tau,delta,sx,sy,rho,reli)
!    c(1) = cscale(tau,rho,delta,sx,sy)

!    c = linedisp(tau,rho,phi,delta,sx,sy)
!    write(*,*) 'Constrains init:',real(c),real(sx),real(sy)!,real(sigu),real(sigv)

!    stop

!    tau = 2
!    rho = 4123
!    phi = 14/57.3
    !    write(*,*) tau,rho,phi,sx,sy
!    block
!      real(REAL64), dimension(2) :: s

!      s = linescale(tau,rho,delta,sx,sy)
!      s2 = s2**2
!      write(*,*) 'linescale:',real(s)!,sum(s**2), sum(c**2)/2
!      call thexscale(tau,rho,delta,sx,reli)
!      call theyscale(tau,rho,delta,sy,reli)
!      stop
!    end block


!    block
!      real(REAL64), dimension(size(mju)+2) :: z, fvec
!      z(1:2) = [ tau, rho ]
!      z(3:) = delta
!      fvec = gradf(tau,rho,delta,sx,sy)
!      write(*,*) loglike(z,sx,sy), sqrt(sum(fvec**2)/(size(mju)+2))
!    end block
!    stop


!    block
!      call graph_rhos(tau,rho,cts,pht,delta,sx,sy)
!      stop
!    end block




!    c = 0
    stol = 0.01
!    nu = 1
!    lambda = nu*c

    write(*,*) 'Initial values:',real(tau),real(rho)!,real(delta(1))

    ! winsorising of initial values for dispersion esimates
!    do i = 1, size(cts)
!    end do


!    sx = 1
!    sy = 1
!    delta = mju - pht

!    block
!      call graph_scale(tau,mju,cts,pht,rho)
!      stop
!    end block

!    stop

!    su = median(abs(cts - mju*tau)/sqrt(mju*tau + sigu**2))
!    sv = median(abs(pht - mju)/sqrt(mju + sigv**2))
!    write(*,*) 'su,sv:',su,sv,sigu,sigv

!    call graph_disp(tau,rho,delta)
!    stop

    block ! plot functions

      real(REAL64) :: tau1, rho1, dt, dr

      tau1 = 2
      rho1 = 4123
      dt = 0.01
      dr = 5e2
      call graph_fun('/tmp/phl',tau1,dt,rho1,dr,mju)

    end block


    goto 20
!    mjus = (cts/tau + pht) / 2
    block
      real(REAL64) :: tau1, rho1, dt, dr
      tau1 = tau
      rho1 = rho
      nfev = 0
      call therline(tau1,dtau,rho1,drho,delta,reli)
      write(error_unit,*) 'Nfev=',nfev

      dt = 0.01
      dr = 1e3
      call graph_rlike('/tmp/phline',tau1,dt,rho1,dr,mju)
      !    mjus = (cts/tau + pht) / 2
!      if( reli ) then
!         mjus = (cts/tau + pht) / 2
!      end if
    end block
20 continue

    block
      real(REAL64) :: tau1, rho1, dt, dr
      tau1 = tau
      rho1 = rho

      mju = (cts/tau + pht) / 2
      call fitline(tau1,dtau,rho1,drho,delta,reli)

    end block


!    call thecline(tau,dtau,rho,drho,delta,reli)
    goto 22
    block
      real(REAL64) :: d,e,mu,sigx, sigy

      do i = 1, size(cts)
         d = delta(i)
         e = cts(i) - tau*(pht(i) + d)
         sigy = sqrt(cts(i)*(1 + tau)*tau + rho**2)
         sigx = sqrt(pht(i)*(1 + 1/tau) + (rho/tau)**2)
         if( abs(d/sigx) <= zmax .and. abs(e/sigy) <= zmax ) then
            mjus(i) = (cts(i)/tau + pht(i)) / 2
         else if( abs(d/sigx) <= zmax .and. abs(e/sigy) > zmax ) then
            mjus(i) = pht(i)
         else if( abs(d/sigx) > zmax .and. abs(e/sigy) <= zmax ) then
            mjus(i) = cts(i) / tau
         else
            ! an arbitrary choice
            mjus(i) = (cts(i)/tau + pht(i)) / 2
         end if
!            write(error_unit,'(a,i6,1x,2(g0,2x),2x,2(f10.3,1x),3x,g0)') &
!                 'Outlier!:',i,real(x(i)),real(y(i)),d/sigx,e/sigy,real(mjus(i))
!         write(error_unit,'(a,i6,1x,2(g0,2x),2x,2(f10.3,1x),3x,g0)') &
!              'Updates:',i,real(x(i)),real(y(i)),d/sigx,e/sigy,real(mjus(i))
      end do
    end block
22  continue


!    write(*,*) real(cts)
!    write(*,*) real(pht)


!    goto 44

    converge = .false.

    !c = logvar(tau,rho) - scale
!    nu = 1
    lambda = 0.5 !(nu/2)*c
    eta_tol = 1e-3 ! tolerance for c
    eta_tol = sqrt(2.0/ndata)      ! tolerance for c by variance of chi2
    omega_tol = 10*sqrt(epsilon(tau)) ! tolerance for fvec**2
    nu = 10
    omega = 1/nu
    eta = 1/nu**0.1
    nfev = 0

    write(error_unit,*) 'Tolerances(eta,omega):',real(eta_tol), real(omega_tol)

    do iter = 1, -3*precision(tau)

!       call xscale(tau,delta,sx,sigv,reli)
!       call yscale(tau,delta,sy,sigu,reli)

!       tauback = tau
!       call xscale(sigv,dsigv,sv,xreli)
!       call yscale(sigu,dsigu,su,yreli)
!       stop

!       reli = xreli .and. yreli
!       tauback = tau
!       write(*,*) '-----',sigu,sigv,su,sv

!       if( reli ) then
!          siguback = sigu
!          sigvback = sigv
!       else
!          call graph_xscale('/tmp/sv',sigv)
!          call graph_yscale('/tmp/su',sigu)
!          exit
!       end if


!       backrho = rho
!       call theline1(tau,dtau,rho,drho,phi,dphi,delta,reli)
!       call theline(tau,dtau,rho,drho,sv,su,delta,reli)
!       call theline(tau,dtau,sigu,dsigu,sigv,dsigv,su,sv,delta,reli)
       !       call thefline(tau,dtau,rho,drho,phi,dphi,sx,sy,delta,reli)
!       sx = 1
       !       sy = 1

       ! to pass omega???
       call theline(tau,dtau,rho,drho,omega,reli)

       if( reli ) then
          c = cvar(tau,rho) !/ ndata
!          sig2 = c + (ndata - 1)*beta*(
!          c = sig2 - scale
       else
!          sig2 = -1
          c = huge(c)
       end if


       ! update scale
!       sx = 1
!       sy = 1
!       call thescale(tau,delta,sx,sy,rho,reli)
!       sc = linescale(tau,rho,phi,delta,sx,sy)
!       sx = sc(1)
!       sy = sc(2)


!       tauback = tau
!       mju = pht + delta
!       call thescale(rho,drho,s,reli)

!       sigv = rho/tau*cos(phi)
!       sigu = rho*sin(phi)

!       sx = xdisp(tau,delta,sigv,sx)
!       sy = ydisp(tau,delta,sigu,sy)
!       c = [sx, sy] - scale
       !       c = [ xdisp(tau,delta,sigv,sx), ydisp(tau,delta,sigu,sy) ]
       !       c = disp(tau,rho,phi,delta,sx,sy)

!       c = scale - sc

       if( debug_iter ) then
          block
            character(len=666) :: filename
            real(REAL64) :: tau0, rho0, dt, dr
!           real(REAL64), dimension(size(delta)+1) :: z
!            z(1:2) = [ tau, rho ]
!           z(1) = tau
!           z(2:) = delta
!            write(error_unit,*) 'Loglike:',loglike(z,rho,sx,sy), &
            !                 '    s=',real(sqrt(sum(delta**2)/size(delta)))
!            if( reli ) then
!            sig2 =
!            c = logvar(tau,rho) - scale
!            else
!               sig2 = -1
!               c = 666
!            end if
!            write(error_unit,*) 'c:',real(c)
!            write(error_unit,*) 'lambda:',real(lambda),'      c:',real(c)

!            write(error_unit,*) 's:',real(sx),real(sy)

!103         format(a,i0,(1x,a,g0,a,es12.3),2x,(a,g0,a,es8.1,2x),l2)
!5         (1x,a,f0.4),1x,1(1x,a,g0),1x,l1)
!            write(error_unit,103) 'phline #',iter,'tau=',real(tau),'+-',dtau, &
!                 'rho=',real(rho),'+-',drho,reli

104         format(a,i0,1x,(a,g0),1x,(a,i0),2x,1(a,f0.3,1x),(a,es7.1,1x),(a,f0.1),l3)
            write(error_unit,104) 'phline #',iter,'tau=',real(tau),'rho=', &
                 nint(rho),'λ=',lambda,'|c|=',abs(c),'log(ν)=',log10(nu),reli

            tau0 = tau
            rho0 = rho
            dt = 0.01
            dr = 1e3
            write(filename,'(a,i0)') '/tmp/phline_',iter
            call graph_like(filename,tau0,dt,rho0,dr,mju)
          end block

       end if

       if( .not. reli ) exit

!       write(error_unit,*) real(abs(c)),real(eta),real(omega)
       if( abs(c) < eta ) then

          converge = abs(c) < eta_tol .and. omega < omega_tol
!          converge = abs(c) < 1/sqrt(2.0*ndata) .and. omega < omega_tol
          if( converge ) exit

          lambda = lambda - nu*c
          eta = eta / nu**0.9
          omega = omega / nu

       else

          nu = 100*nu
          eta = 1 / nu**0.1
          omega = 1 / nu

       end if
       omega = max(omega, omega_tol/2)


!       mju = pht + delta
!       sc = linescale(tau,rho,delta,sx,sy)
!       write(*,*) 'sc:',real(sc)
!       sx = sc(1)
!       sy = sc(2)
!       call thexscale(tau,rho,delta,sx,reli)
!       call theyscale(tau,rho,delta,sy,reli)
       !       c = linedisp(tau,rho,delta,sx,sy)
!       c(2:3) = scale - [sx, sy]

!       sx = 1
!       sy = 1
!       c(1) = cscale(tau,rho,delta,sx,sy)
       !       c = rdisp(tau,rho) - scale

!       c = [sx, sy] - scale
!       c = s - scale
!       converge = abs(c) < stol
!       if( converge .or. .not. reli ) exit

       !       mju = pht + delta

       if ( .false. ) then
          lambda = lambda - nu*c
          nu = 2*nu
       end if

    end do
    write(error_unit,*) 'Nfev=',nfev

44  continue


    !    tau = 2
    delta = mju - pht
    call graph_ecdf(tau,rho,phi,cts,pht,delta)
    call graph_muline(tau,rho,delta,cts,pht)



!    call graph_rhos(tau,rho,cts,pht,delta,sx,sy)
!    call graph_mutau(tau,mju,cts,pht,sigu,sigv)

!    stop 'phline'

    n = 0
    do i = 1, -size(mjus)
       if( valid(i) ) then
          n = n + 1
          mjus(i) = pht(n) + delta(n)
       else
          mjus(i) = (y(i)/tau + x(i)) / 2
       end if
    end do

    deallocate(cts,pht,mju,idx)

    if( .not. reli ) then
       tau = tau0
       rho = rho0
       mjus = (y/tau + x) / 2
    end if

  end subroutine phline


  subroutine fitline(tau,dtau,rho,drho,delta,reli)

    use sqps

    real(REAL64), intent(in out) :: tau, rho
    real(REAL64), intent(out) :: dtau, drho
    real(REAL64), dimension(:), intent(out) :: delta
    logical, intent(out) :: reli

    integer :: maxfev = 666
    real(REAL64), dimension(2) :: x
    real(REAL64), dimension(2,2) :: cov
    character(len=666) :: msg
    character(len=6) :: label
    real(REAL64) :: tol, lam
    integer :: npar, info, n, i, K, nprint

    n = size(mju)
    npar = 2
    nprint = 1

    x(1) = tau
    x(2) = rho
    lam = 1

!    call locsqp(slike,x,lam,quasihess=.true.,nprint=nprint,info=info,msg=msg)

    reli = tmin < x(1) .and. x(1) < tmax &
         .and. abs(x(2)) < rhomax &
         .and. info == 1

    if( debug_iter ) then
       write(error_unit,'(a,i0,3a,l2)') ' fitline(): info=', &
            info,': ',trim(msg),'  reliable:',reli
       write(error_unit,*) 'solution:',real(x),' lambda=',real(lam)
    end if

    if( reli ) then
       tau = x(1)
       rho = x(2)
!       call locsqp(slike,x,lam,cov,quasihess=.false.,nprint=nprint,info=info,msg=msg)
       if( info == 1 ) then
          tau = x(1)
          rho = x(2)
          if( cov(1,1) > 0 ) then
             dtau = sqrt(cov(1,1))
          end if
          if( cov(2,2) > 0 ) then
             drho = sqrt(cov(2,2))
          end if
       end if
    end if

    if( debug_iter .and. (0 < info .and. info < 5) ) then
       block
         integer, parameter :: maxpar = 7
         real(REAL64), dimension(npar,npar) :: cor
!         real(REAL64), dimension(3) :: xvec,fvec
         real(REAL64) :: q
         integer :: i,j,n,iflag
         character(len=3) :: dots
         character(len=80) :: fmt

         cor = 1
         do i = 1, npar
            do j = i+1, npar
               if( cov(i,i) > 0 .and. cov(j,j) > 0 &
                    .and. cov(i,i) < 1e10 .and. cov(j,j) < 1e10) then
                  q = sqrt(cov(i,i)*cov(j,j))
!                  write(*,*) cov(i,i),cov(j,j),q,cov(i,j)
                  if( 1e-10 < q .and. q < 1e10 .and. abs(cov(i,j)) < 1e10 ) then
                     cor(i,j) = cov(i,j) / q
                  else
                     cor(i,j) = huge(cor)
                  end if
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         n = min(maxpar,npar)
         if( n < npar ) then
            dots = '...'
         else
            dots = ''
         end if

         write(error_unit,*) 'Covariance (Hessian in minimum):'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(es9.2,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cov(i,1:n),dots
         end do

         write(error_unit,*) 'Correlation:'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(f6.3,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cor(i,1:n),dots
         end do

       end block
    end if

!    call graph_line(x(1),x(2),cts,pht,x(3:npar-1),x(npar))
!    call graph_rho(x(1),x(2),cts,pht,x(3:npar-1),x(npar))

  end subroutine fitline

  subroutine slike(x,f,fvec,fjac,c,cvec,cjac,iflag)

    use rfun

    real(REAL64), intent(out) :: f,c
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), intent(out) :: fvec, cvec
    real(REAL64), dimension(:,:), intent(out) :: fjac, cjac
    integer, intent(in out) :: iflag

    character(len=666) :: msg
    real(REAL64), dimension(:), allocatable :: sigx, sigy, zx, zy, psix, psiy, &
         dpsix, dpsiy, sigx_t, sigx_r, sigy_t, sigy_r
    real(REAL64), dimension(2,2) :: s2jac, dcc
    real(REAL64), dimension(2) :: dc, s2grad
    real(REAL64) :: tau, rho, ds
    integer :: i, j, stat

    tau = x(1)
    rho = x(2)

    f = loglike(tau,rho,mju)

    allocate(sigx(ndata),sigy(ndata),zx(ndata),zy(ndata),psix(ndata),psiy(ndata), &
         dpsix(ndata),dpsiy(ndata),sigx_t(ndata),sigx_r(ndata),sigy_t(ndata), &
         sigy_r(ndata),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       if( debug_iter ) write(error_unit,*) 'phline(slike): ',trim(msg)
       iflag = -1
       return
    end if

    sigx = sigmax(tau,rho,pht)
    sigy = sigmay(tau,rho,cts)
    zx = (mju - pht) / sigx
    zy = (cts - tau*mju) / sigy
    if( robust ) then
       psix = huber(zx)
       psiy = huber(zy)
       dpsix = dhuber(zx)
       dpsiy = dhuber(zy)
!       call hubers(zx,psix,dpsix) update !
       !       call hubers(zy,psiy,dpsiy) update !
!      call tukeys(zx,psix,dpsix)
       !      call tukeys(zy,psiy,dpsiy)

    else
       psix = zx
       psiy = zy
       dpsix = 1
       dpsiy = 1
    end if



    ds = sum(dpsix + dpsiy)
!    if( .not. (ds > 0) ) then
!       if( debug_iter ) write(error_unit,*) 'phline(slike): assumption ds > 0'
!       iflag = -1
!       return
!    end if

    c = sum(psix**2 + psiy**2) - beta*ds
    c = sum(psix**2 + psiy**2) - 2*ndata
    c = sum(psix**2 + psiy**2) - beta*(2*ndata)
!    c = sum(psix**2 + psiy**2) - ds
    write(*,*) sum(psix**2 + psiy**2),ds,c,beta
!    c = c / ndata

    do i = 1, ndata
       call grad_sigx(sigx(i),pht(i),tau,rho,sigx_t(i),sigx_r(i))
       call grad_sigy(sigy(i),cts(i),tau,rho,sigy_t(i),sigy_r(i))
    end do

    fvec = grad_ressum(tau,rho,mju,zx,zy,psix,psiy,dpsix,dpsiy,sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r) !/ ndata
    cvec = grad_sigvar(mju,zx,zy,psix,psiy,dpsix,dpsiy,sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r)  !/ ndata


    if( iflag == 2 ) then

       do i = 1, 2
          do j = 1, 2
             cjac(i,j) = dc(i)*dc(j)
          end do
       end do

       block
         real(REAL64), dimension(:), allocatable :: sigx_tt,sigx_tr,sigx_rr, &
              sigy_tt,sigy_tr,sigy_rr

         allocate(sigx_tt(ndata),sigx_tr(ndata),sigx_rr(ndata), &
              sigy_tt(ndata),sigy_tr(ndata),sigy_rr(ndata),stat=stat,errmsg=msg)
         if( stat /= 0 ) then
            if( debug_iter ) write(error_unit,*) 'phline(thelike): ',trim(msg)
            iflag = -1
            return
         end if

         do i = 1, ndata
            call jac_sigx(sigx(i),pht(i),tau,rho,sigx_t(i),sigx_r(i), &
                 sigx_tt(i),sigx_tr(i),sigx_rr(i))
            call jac_sigy(sigy(i),cts(i),tau,rho,sigy_t(i),sigy_r(i), &
                 sigy_tt(i),sigy_tr(i),sigy_rr(i))
         end do

         fjac = jac_ressum(tau,rho,mju,sigx,sigy,zx,zy,psix,psiy,dpsix,dpsiy, &
              sigx_t,sigx_r,sigy_t,sigy_r, &
              sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr) !/ ndata
         cjac = jac_cvar(tau,rho,mju,sigx,sigy,zx,zy,psix,psiy,dpsix,dpsiy, &
              sigx_t,sigx_r,sigy_t,sigy_r, &
              sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr) !/ ndata

       end block

    end if

    if( iflag == 0 ) then
       if( debug_iter ) then
          block
            real(REAL64) :: logl, var!, lag
            logl = loglike(tau,rho,mju) / ndata
            var = sum(psix**2 + psiy**2) / (2*ndata*beta)
            var = sum(psix**2 + psiy**2) / ds
!            lag = f - 0.5*c
!            var = cvar(tau,rho)
            write(error_unit,*) ' loglike=',real(logl),' disp=',real(var)!, &
!                 'lag=',real(lag)
            write(error_unit,*) 'x = ',real(x)
            write(error_unit,*) 'fvec = ',real(fvec)
            write(error_unit,*) 'cvec = ',real(cvec)
          end block
       end if

!       if( x(1) < tmin .or. x(1) > tmax .or. abs(x(2)) > rhomax ) iflag = -1
       return
    end if



  end subroutine slike


  subroutine therline(tau,dtau,rho,drho,delta,reli)

    use minpacks
!    use medians

    real(REAL64), intent(in out) :: tau, rho
    real(REAL64), intent(out) :: dtau, drho
    real(REAL64), dimension(:), intent(out) :: delta
    logical, intent(out) :: reli

    integer :: maxfev = 666
    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(3) :: x, diag
    real(REAL64), dimension(3,3) :: cov, fjac
    character(len=666) :: msg
    character(len=6) :: label
    real(REAL64) :: tol
    integer :: npar, info, n, i, K

    n = size(mju)
    !    K = size(kmean)
    K = size(mju)
    npar = 3
    maxfev = 666*npar
    tol = sqrt(epsilon(x))
    allocate(xold(npar))

    x(1) = tau
    x(2) = rho
    x(3) = -0.5 !-0.99
    diag(1) = 1/tau
    diag(2) = 1/(rho + 1)
    diag(npar) = 1

!    diag(1) = 1
!    diag(2) = 1e-3
!    diag(npar) = 1

    write(*,*) real(diag(1:3))

    if( analytic ) then
!       tol = 1e2*epsilon(x)!*(1 + abs(maxval(x)))
       label = 'hybrj3'
       call hybrj3(rlike,x,tol,maxfev,factor=factor,diag=diag,nprint=nprint,info=info,msg=msg)
    else
!       tol = 1e-7
       label = 'hybrd3'
!       call hybrd3(cliked,x,tol,maxfev,factor=factor,diag=diag,nprint=nprint,info=info,msg=msg)
       call hybrd3(rliked,x,tol,maxfev,factor=factor,nprint=nprint,info=info,msg=msg)
    end if

    reli = tmin < x(1) .and. x(1) < tmax &
         .and. abs(x(2)) < rhomax &
         .and. info == 1

    if( debug_iter ) then
       write(error_unit,'(3a,i0,3a,l2)') ' theline() ',label,' info=', &
            info,': ',trim(msg),'  reliable:',reli
       write(error_unit,*) 'Lambda=',real(x(npar))
    end if

!       jac = linerjac(tau,rho)
!       call qrinv(jac,cov)
!       do i = 1, 2
!          write(*,*) i,real(cov(i,i))
!       end do

!       do i = 1, size(delta)
!          write(*,*) i,real(delta(i)-x(2+i))
!       end do

    if( reli ) then
       ! Sometimes, the lambda parameter of LM-method is incidentally tunneled
       ! into the mirrored part of the parameter space. Fortunately, the like
       ! function is pretty symmetric, so I do to ignore the "tunnelling".
       ! The correct, in principle, solution should to keep the condition:
       ! rho > 0.
       tau = x(1)
       rho = abs(x(2))
       lambda = x(3)
!       delta = x(3:npar-1)



       ! covariance matrix
       ! jac =
       ! cov =


!       if( all(diag > 0) ) then
!          dtau =  sqrt(cov(1,1)/(2*n - npar))
!          drho = sqrt(cov(2,2)/(2*n - npar))
!          dmu = sqrt(cov(3:,3:)/(2*n - npar))
!          dphi = sqrt(cov(3,3)/(2*n - npar))
          !       dsigv = sqrt(cov(3,3)/(2*n - npar))
!       end if
    end if

    if( debug_iter .and. (0 < info .and. info < 5) ) then
       block
         integer, parameter :: maxpar = 7
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64), dimension(3) :: xvec,fvec
         real(REAL64) :: q
         integer :: i,j,n,iflag
         character(len=3) :: dots
         character(len=80) :: fmt

         iflag = 2
         xvec = x
         xvec(3) = 0
         call rlike(npar,xvec,fvec,fjac,npar,iflag)
         fjac = matmul(transpose(fjac),fjac)
!         call qrinv(fjac,cov)
         cov = fjac

         cor = 1
         do i = 1, npar
            do j = i+1, npar
               if( cov(i,i) > 0 .and. cov(j,j) > 0 &
                    .and. cov(i,i) < 1e10 .and. cov(j,j) < 1e10) then
                  q = sqrt(cov(i,i)*cov(j,j))
!                  write(*,*) cov(i,i),cov(j,j),q,cov(i,j)
                  if( 1e-10 < q .and. q < 1e10 .and. abs(cov(i,j)) < 1e10 ) then
                     cor(i,j) = cov(i,j) / q
                  else
                     cor(i,j) = huge(cor)
                  end if
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         n = min(maxpar,npar)
         if( n < npar ) then
            dots = '...'
         else
            dots = ''
         end if

         write(error_unit,*) 'Covariance (Hessian in minimum):'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(es9.2,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cov(i,1:n),dots
         end do

         write(error_unit,*) 'Correlation:'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(f6.3,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cor(i,1:n),dots
         end do

       end block
    end if

!    call graph_line(x(1),x(2),cts,pht,x(3:npar-1),x(npar))
!    call graph_rho(x(1),x(2),cts,pht,x(3:npar-1),x(npar))

    deallocate(xold)

  end subroutine therline

  subroutine rlike(npar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(npar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

    real(REAL64) :: s0, logl, d, lambda
    integer :: i, k

    if( iflag == 0 ) then
       if( debug_iter ) then
          if( xold(1) > 0 ) then
             s0 = sum((x - xold)**2)
          else
             s0 = -1
          end if
          logl = rloglike(x(1),x(2),mju)
          d = rdisp(x(1),x(2)) / scale
          write(error_unit,*) 's0 = ',real(s0),' loglike=',real(logl),&
               ' disp=',real(d)
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
          xold = x
       end if

       if( x(1) < tmin .or. x(1) > tmax .or. abs(x(2)) > rhomax ) iflag = -1

       ! convergence check ?
!       if( x(1) < tmin .or. x(1) > tmin .or. abs(x(2)) > rhomax ) &
!            iflag = -1
!       if( x(1) < tmin .or. x(1) > tmin .or. x(2) < 0 .or. &
!            x(3) < 0.1 .or. x(3) > 10 .or. x(4) < 0.1 .or. x(4) > 10 ) then
!          fvec = x(1)*maxval(cts)
!       end if
       return

    else if( iflag == 1 ) then

       lambda = x(npar)
       fvec(1:2) = gradrline(x(1),x(2))*(1 + lambda)
       fvec(npar) = rdisp(x(1),x(2)) - scale
       nfev = nfev + 1

    else if( iflag == 2 ) then

       lambda = x(npar)
       fjac(1:2,1:2) = linerjac(x(1),x(2))*(1 + lambda)
       fjac(npar,1:npar-1) = gradrline(x(1),x(2))
       fjac(1:npar-1,npar) = fjac(npar,1:npar-1)
       fjac(npar,npar) = 0

    end if

  end subroutine rlike


  subroutine rliked(npar,x,fvec,iflag)

    integer, intent(in) :: npar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(npar), intent(out) :: fvec
!    real(REAL64), dimension(3) :: c
!    real(REAL64), dimension(:,:), allocatable :: dc
    !    real(REAL64), dimension(:), allocatable :: sc
!    real(REAL64), parameter :: sx = real(1.0,REAL64), sy = real(1.0,REAL64)
    real(REAL64) :: logl, d, lambda
    integer :: i, k
    logical :: r

    if( iflag == 0 ) then
       if( debug_iter ) then
          logl = rloglike(x(1),x(2),mju)
          r = robust
          robust = .false.
          d = rdisp(x(1),x(2))
          robust = r
          d = rdisp(x(1),x(2))
          write(error_unit,*) 'rliked: loglike=',real(logl),' disp=',real(d)
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
       end if

       if( x(1) < tmin .or. x(1) > tmax .or. abs(x(2)) > rhomax ) iflag = -1

       return
    end if

    nfev = nfev + 1
    lambda = x(npar)
    fvec(1:2) = gradrline(x(1),x(2))*(1 + lambda)
    fvec(npar) = rdisp(x(1),x(2)) - scale

  end subroutine rliked

  real(REAL64) function rdisp(tau,rho)

    use rfun

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), parameter :: a = 1.345
    real(REAL64) :: s, d, mu, e, zx, zy, sigx, sigy, psi, dpsi, ds
    integer :: i,n

    n = size(cts)

    s = 0
    ds = 0
    do i = 1, n
!       mu = (cts(i)/tau + pht(i)) / 2
       mu = mju(i)
       d = mu - pht(i)
!       d = delta(i)
!       mu = pht(i) + d
       e = cts(i) - tau*mu

       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)
       zx = d / sigx
       zy = e / sigy
       if( robust ) then
!          zx = max(-a,min(zx,a))
!          zy = max(-a,min(zy,a))
!          ipsi = (zx**2 + zy**2)/2
!          ipsi = ihuber(zx) + ihuber(zy)
          psi = huber(zx)**2 + huber(zy)**2
          dpsi = dhuber(zx) + dhuber(zy)
       else
          psi = zx**2 + zy**2
          dpsi = 2
       end if

       s = s + psi
       ds = ds + dpsi

    end do

    if( ds > 0 ) then
       rdisp = s / ds
    else
       stop '(ds > 0) assumption failed.'
    end if

  end function rdisp


  function gradrline(tau,rho) result(vec)

    use rfun

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(2) :: vec

    real(REAL64) :: mu,d,e,zx,zy,psix,psiy,px,py,sigx,sigy, &
         sigx_t, sigx_r, sigy_t, sigy_r, ds, dpsi
    integer :: i,n

    n = size(cts)

    vec = 0
    ds = 0
    do i = 1, n
       !       mu = (cts(i)/tau + pht(i)) / 2
       mu = mju(i)
       d = mu - pht(i)
!       d = pht(i) - mu
       e = cts(i) - tau*mu
!       d = (tau*pht(i) - cts(i)) / 2 / tau
!       e = (tau*pht(i) - cts(i)) / 2

       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)
       call grad_sigx(sigx,pht(i),tau,rho,sigx_t,sigx_r)
       call grad_sigy(sigy,cts(i),tau,rho,sigy_t,sigy_r)

       zx = d / sigx
       zy = e / sigy

       if( robust ) then
          psix = huber(zx)
          psiy = huber(zy)
          dpsi = dhuber(zx) + dhuber(zy)
       else
          psix = zx
          psiy = zy
          dpsi = 2
       end if

       px = psix*zx
       py = psiy*zy

       vec(1) = vec(1) - (px*(sigx_t/sigx) + py*(sigy_t/sigy) + psiy*(mu/sigy))
       vec(2) = vec(2) - (px*(sigx_r/sigx) + py*(sigy_r/sigy))
!       vec(2+i) = psix/sigx - tau*psiy/sigy
       ds = ds + dpsi

    end do

    vec = vec / ds
!    vec = vec / (2*n)

  end function gradrline

  function linerjac(tau,rho) result(jac)

    use rfun

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(2,2) :: jac
    real(REAL64) :: d,e,ds,mu,zx,zy,psix,psiy,dpsix,dpsiy, &
         sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r, &
         sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr,&
         chi2_tt, chi2_tr, chi2_td, chi2_rr, chi2_rd, chi2_dd, &
         px, py, zx2, zy2, sxt2, syt2, sxtr, sytr, sxr2, syr2, dzy_t, dpx, dpy
    integer :: i,j,n

    n = size(cts)
    jac = 0
    ds = 0
    do i = 1, n
       !       mu = (cts(i)/tau + pht(i)) / 2
       mu = mju(i)
       d = mu - pht(i)
!       mu = pht(i) + d
       e = cts(i) - tau*mu
       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)
       call grad_sigx(sigx,pht(i),tau,rho,sigx_t,sigx_r)
       call grad_sigy(sigy,cts(i),tau,rho,sigy_t,sigy_r)
       call jac_sigx(sigx,pht(i),tau,rho,sigx_t,sigx_r,sigx_tt,sigx_tr,sigx_rr)
       call jac_sigy(sigy,cts(i),tau,rho,sigy_t,sigy_r,sigy_tt,sigy_tr,sigy_rr)
       zx = d / sigx
       zy = e / sigy

       if( robust ) then
          psix = huber(zx)
          psiy = huber(zy)
          dpsix = dhuber(zx)
          dpsiy = dhuber(zy)
       else
          psix = zx
          psiy = zy
          dpsix = 1
          dpsiy = 1
       end if

       px = psix*zx
       py = psiy*zy

       zx2 = zx**2
       zy2 = zy**2
       sxt2 = (sigx_t/sigx)**2
       syt2 = (sigy_t/sigy)**2
       sxtr = (sigx_t/sigx)*(sigx_r/sigx)
       sytr = (sigy_t/sigy)*(sigy_r/sigy)
       sxr2 = (sigx_r/sigx)**2
       syr2 = (sigy_r/sigy)**2
       dzy_t = mu/sigy + zy*(sigy_t/sigy)
       dpx = dpsix*zx + psix
       dpy = dpsiy*zy + psiy
       chi2_tt = -dpsix*zx2*sxt2 - 2*px*sxt2 + px*(sigx_tt/sigx) &
            -dpsiy*dzy_t*zy*(sigy_t/sigy) &
            - psiy*(mu/sigy + 2*zy*(sigy_t/sigy))*(sigy_t/sigy) &
            + py*(sigy_tt/sigy) &
            -(dpsiy*dzy_t + psiy*(sigy_t/sigy))*(mu/sigy)
!            -dpsiy*dzy_t*(mu/sigy) &
!            -psiy*(mu/sigy)*(sigy_t/sigy)
!               + dpsiy*zy2*syt2 + 2*py*syt2 - py*(sigy_tt/sigy) &
!            - dpsiy*(mu/sigy)**2 - (dpsiy*zy + psiy)*(sigy_t/sigy)*(mu/sigy)
!               + (dpsiy*(1-zy*(sigy_t/sigy))-psiy)*(mu/sigy)**2
!       chi2_tt = dpsix*zx2*sxt2 + dpsiy*zy2*syt2 + dpsiy*(mu/sigy)**2
!       chi2_tr = -dpsix*zx2*sxtr + 2*px*sxtr - px*(sigx_tr/sigx) &
       chi2_tr =-dpx*zx*(sigx_t/sigx)*(sigx_r/sigx) &
               - dpy*dzy_t*(sigy_r/sigy) &
               - px*(sigx_t/sigx)*(sigx_r/sigx) + px*(sigx_tr/sigx) &
               - py*(sigy_t/sigy)*(sigy_r/sigy) + py*(sigy_tr/sigy)

!               + dpsiy*zy2*sytr + 2*py*sytr - py*(sigy_tr/sigy) &
!               - psiy*(mu/sigy)*(sigy_r/sigy)
!       chi2_td =-px*(sigx_t/sigx) + py*tau*(sigy_t/sigy) &
!            + (dpsiy - psiy)*(mu/sigy)*(tau/sigy)
       chi2_td =-dpx*(sigx_t/sigx)/sigx &
               + dpsiy*dzy_t*tau/sigy + psiy*(sigy_t/sigy)*tau/sigy - psiy/sigy
       chi2_rr = dpsix*zx2*sxr2 + (2*sxr2 - sigx_rr/sigx)*px &
               + dpsiy*zy2*syr2 + (2*syr2 - sigy_rr/sigy)*py
       chi2_rd =-dpx*(sigx_r/sigx)/sigx + dpy*(sigy_r/sigy)*tau/sigy
       chi2_dd = dpsix/sigx**2 + dpsiy*(tau/sigy)**2

       jac(1,1:2) = jac(1,1:2) - [ chi2_tt, chi2_tr]
       jac(2,2) = jac(2,2) + chi2_rr
!       j = i + 2
!       jac(1,j) = chi2_td
!       jac(2,j) = chi2_rd
!       jac(j,j) = chi2_dd

       ds = ds + dpsix + dpsiy

    end do

    jac = jac / ds
!    jac = jac / (2*n)

    do i = 1, size(jac,1)
       do j = 1, i - 1
          jac(i,j) = jac(j,i)
       end do
    end do

  end function linerjac


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine thecline(tau,dtau,rho,drho,delta,reli)

    use minpacks
!    use medians

    real(REAL64), intent(in out) :: tau, rho
    real(REAL64), intent(out) :: dtau, drho
    real(REAL64), dimension(:), intent(out) :: delta
    logical, intent(out) :: reli

    integer :: maxfev = 666
    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(:), allocatable :: x, diag
    real(REAL64), dimension(:,:), allocatable :: cov, jac
    character(len=666) :: msg
    character(len=6) :: label
    real(REAL64) :: tol
    integer :: npar, info, n, i, K

    n = size(mju)
    !    K = size(kmean)
    K = size(mju)
    npar = K + 3
    maxfev = 666*npar
    allocate(x(npar),diag(npar),xold(npar))
    xold = -1

    x(1) = tau
    x(2) = rho + 1
    x(3:npar-1) = delta
    x(npar) = 0 ! lambda
    diag(1) = 1/tau
    diag(2) = 1/(rho + 1)
    diag(3:npar-1) = 1/sqrt(min(pht,cts/tau)*(1/tau+1) + (rho/tau)**2)
    diag(npar) = 1

    diag(1) = 1
    diag(2) = 1e-3
    diag(3:npar-1) = 1e-4
    diag(npar) = 1

    write(*,*) real(diag(1:3))

    if( analytic ) then
       tol = 1e2*epsilon(x)!*(1 + abs(maxval(x)))
       tol = sqrt(epsilon(x))
       label = 'hybrj3'
       x(npar) = -0.9 ! lambda
!       jac = linejac(tau,rho,delta)
!       do i = 1, npar-1
!          diag(i) = sqrt(abs(jac(i,i)))
!       end do
!    write(*,*) 'diag:',real(diag)
       call hybrj3(clike,x,tol,maxfev,factor=factor,diag=diag,nprint=nprint,info=info,msg=msg)
    else
       tol = 1e-7
       label = 'hybrd3'
!       call hybrd3(cliked,x,tol,maxfev,factor=factor,diag=diag,nprint=nprint,info=info,msg=msg)
       call hybrd3(cliked,x,tol,maxfev,factor=factor,nprint=nprint,info=info,msg=msg)
    end if

    reli = tmin < x(1) .and. x(1) < tmax &
         .and. abs(x(2)) < rhomax &
         .and. info == 1

    if( debug_iter ) then
       write(error_unit,'(3a,i0,3a,l2)') ' theline() ',label,' info=', &
            info,': ',trim(msg),'  reliable:',reli
       write(error_unit,*) 'Lambda=',real(x(npar))
    end if

       jac = linejac(tau,rho,delta)
       do i = 1, npar-1
          write(*,*) i,real(jac(i,i))
       end do

       do i = 1, size(delta)
          write(*,*) i,real(delta(i)-x(2+i))
       end do

    if( reli ) then
       ! Sometimes, the lambda parameter of LM-method is incidentally tunneled
       ! into the mirrored part of the parameter space. Fortunately, the like
       ! function is pretty symmetric, so I do to ignore the "tunnelling".
       ! The correct, in principle, solution should to keep the condition:
       ! rho > 0.
       tau = x(1)
       rho = abs(x(2))
       delta = x(3:npar-1)



       ! covariance matrix
       ! jac =
       ! cov =


!       if( all(diag > 0) ) then
!          dtau =  sqrt(cov(1,1)/(2*n - npar))
!          drho = sqrt(cov(2,2)/(2*n - npar))
!          dmu = sqrt(cov(3:,3:)/(2*n - npar))
!          dphi = sqrt(cov(3,3)/(2*n - npar))
          !       dsigv = sqrt(cov(3,3)/(2*n - npar))
!       end if
    end if

!!$    if( debug_iter .and. (0 < info .and. info < 5) ) then
!!$       block
!!$         integer, parameter :: maxpar = 7
!!$         real(REAL64), dimension(npar,npar) :: cor
!!$         real(REAL64) :: q
!!$         integer :: i,j,n
!!$         character(len=3) :: dots
!!$         character(len=80) :: fmt
!!$
!!$         cor = 1
!!$         do i = 1, npar
!!$            do j = i+1, npar
!!$               if( cov(i,i) > 0 .and. cov(j,j) > 0 &
!!$                    .and. cov(i,i) < 1e10 .and. cov(j,j) < 1e10) then
!!$                  q = sqrt(cov(i,i)*cov(j,j))
!!$!                  write(*,*) cov(i,i),cov(j,j),q,cov(i,j)
!!$                  if( 1e-10 < q .and. q < 1e10 .and. abs(cov(i,j)) < 1e10 ) then
!!$                     cor(i,j) = cov(i,j) / q
!!$                  else
!!$                     cor(i,j) = huge(cor)
!!$                  end if
!!$               else
!!$                  cor(i,j) = huge(cor)
!!$               end if
!!$               cor(j,i) = cor(i,j)
!!$            end do
!!$         end do
!!$
!!$         n = min(maxpar,npar)
!!$         if( n < npar ) then
!!$            dots = '...'
!!$         else
!!$            dots = ''
!!$         end if
!!$
!!$         write(error_unit,*) 'Covariance (Hessian in minimum):'
!!$         write(fmt,'(a,i0,a)') '(i2,a,',n,'(es9.2,1x),a)'
!!$         do i = 1, n
!!$            write(error_unit,fmt) i,': ',cov(i,1:n),dots
!!$         end do
!!$
!!$         write(error_unit,*) 'Correlation:'
!!$         write(fmt,'(a,i0,a)') '(i2,a,',n,'(f6.3,1x),a)'
!!$         do i = 1, n
!!$            write(error_unit,fmt) i,': ',cor(i,1:n),dots
!!$         end do
!!$
!!$       end block
!!$    end if

!    call graph_line(x(1),x(2),cts,pht,x(3:npar-1),x(npar))
    call graph_rho(x(1),x(2),cts,pht,x(3:npar-1),x(npar))

  end subroutine thecline

  subroutine clike(npar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(npar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

    real(REAL64) :: s0, logl, d, lambda
    character(len=3) :: dots
    integer :: i, k, m

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          m = min(size(x),7)
          if( size(x) > 7 ) dots = '...'
          s0 = sqrt(sum(x(3:npar-1)**2)/(size(x)-3))
          if( xold(1) > 0 ) then
             s0 = sum((x - xold)**2)
          end if
          logl = cloglike(x(1),x(2),x(3:npar-1))
          d = disp(x(1),x(2),x(3:npar-1))
          write(error_unit,*) 's0 = ',real(s0),' loglike=',real(logl),&
               ' disp=',real(d)
          write(error_unit,*) 'x = ',real(x(1:m)),trim(dots),real(x(npar))
          write(error_unit,*) 'fvec = ',real(fvec(1:m)),trim(dots),real(fvec(npar))
          xold = x
       end if

       if( x(1) < tmin .or. x(1) > tmax .or. abs(x(2)) > rhomax ) iflag = -1

       ! convergence check ?
!       if( x(1) < tmin .or. x(1) > tmin .or. abs(x(2)) > rhomax ) &
!            iflag = -1
!       if( x(1) < tmin .or. x(1) > tmin .or. x(2) < 0 .or. &
!            x(3) < 0.1 .or. x(3) > 10 .or. x(4) < 0.1 .or. x(4) > 10 ) then
!          fvec = x(1)*maxval(cts)
!       end if
       return

    else if( iflag == 1 ) then

       lambda = x(npar)
       fvec(1:npar-1) = gradline(x(1),x(2),x(3:npar-1))*(1 + lambda)
       fvec(npar) = disp(x(1),x(2),x(3:npar-1)) - scale

    else if( iflag == 2 ) then

       lambda = x(npar)
       fjac(1:npar-1,1:npar-1) = linejac(x(1),x(2),x(3:npar-1))*(1 + lambda)
       fjac(npar,1:npar-1) = gradline(x(1),x(2),x(3:npar-1))
       fjac(1:npar-1,npar) = fjac(npar,1:npar-1)
       fjac(npar,npar) = 0

    end if

  end subroutine clike


  subroutine cliked(npar,x,fvec,iflag)

    integer, intent(in) :: npar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(npar), intent(out) :: fvec
!    real(REAL64), dimension(3) :: c
!    real(REAL64), dimension(:,:), allocatable :: dc
    !    real(REAL64), dimension(:), allocatable :: sc
!    real(REAL64), parameter :: sx = real(1.0,REAL64), sy = real(1.0,REAL64)
    real(REAL64) :: s0, logl, d, lambda
    character(len=3) :: dots
    integer :: i, k, m

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          m = min(size(x),7)
          if( size(x) > 7 ) dots = '...'
          s0 = sqrt(sum(x(3:npar-1)**2)/(size(x)-3))
          logl = cloglike(x(1),x(2),x(3:npar-1))
          d = disp(x(1),x(2),x(3:npar-1))
          write(error_unit,*) 's0 = ',real(s0),' loglike=',real(logl),&
               ' disp=',real(d)
          write(error_unit,*) 'x = ',real(x(1:m)),trim(dots),real(x(npar))
          write(error_unit,*) 'fvec = ',real(fvec(1:m)),trim(dots),real(fvec(npar))
       end if

       if( x(1) < tmin .or. x(1) > tmax .or. abs(x(2)) > rhomax ) iflag = -1

       ! convergence check ?
!       if( x(1) < tmin .or. x(1) > tmin .or. abs(x(2)) > rhomax ) &
!            iflag = -1
!       if( x(1) < tmin .or. x(1) > tmin .or. x(2) < 0 .or. &
!            x(3) < 0.1 .or. x(3) > 10 .or. x(4) < 0.1 .or. x(4) > 10 ) then
!          fvec = x(1)*maxval(cts)
!       end if
       return
    end if

    lambda = x(npar)
!    fvec(1:npar-1) = linegrad2(x(1),x(2),x(3:npar-1)) &
    !         + lambda*grad_disp(x(1),x(2),x(3:npar-1))
    fvec(1:npar-1) = linegrad2(x(1),x(2),x(3:npar-1))*(1 - lambda)
    fvec(npar) = disp(x(1),x(2),x(3:npar-1)) - scale

  end subroutine cliked

  function gradline(tau,rho,delta) result(vec)

    use rfun

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(size(delta)+2) :: vec

    real(REAL64) :: mu,d,e,zx,zy,psix,psiy,px,py,sigx,sigy, &
         sigx_t, sigx_r, sigy_t, sigy_r
    integer :: i,n

    n = size(delta)

    vec = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       e = cts(i) - tau*mu

       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)
       call grad_sigx(sigx,pht(i),tau,rho,sigx_t,sigx_r)
       call grad_sigy(sigy,cts(i),tau,rho,sigy_t,sigy_r)

       zx = d / sigx
       zy = e / sigy

       if( robust ) then
          psix = huber(zx)
          psiy = huber(zy)
       else
          psix = zx
          psiy = zy
       end if

       px = psix*zx
       py = psiy*zy

       vec(1) = vec(1) - (px*(sigx_t/sigx) + py*(sigy_t/sigy) + psiy*(mu/sigy))
       vec(2) = vec(2) - (px*(sigx_r/sigx) + py*(sigy_r/sigy))
       vec(2+i) = psix/sigx - tau*psiy/sigy

    end do

    vec = vec / (2*n)

  end function gradline

  function linegrad2(tau,rho,delta) result(vec)

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(size(delta)+2) :: vec
    real(REAL64), dimension(3) :: gradx, grady

    real(REAL64) :: mu,d,e,r,s,zx,zy,sigx,sigy,rx,ry,rho2,tt, sig2, &
         grad_t, grad_r, grad_d, sigx_t, sigx_r, sigy_t, sigy_r, &
         gradx_t, gradx_r, gradx_d, grady_t, grady_r, grady_d
    integer :: i,n

    n = size(delta)

    tt = tau*(tau + 1)
    rho2 = rho**2

    vec = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d

       sigx = sqrt(pht(i)*(1 + 1/tau) + (rho/tau)**2) / 2
       zx = d / sigx

       if( zx < zmin .or. zx > zmax ) then
          zx = max(zmin,min(zx,zmax))
          gradx_t = 0
          gradx_r = 0
          gradx_d = 0
       else
          sigx_t = -(pht(i) + 2*rho**2/tau)/(8*sigx*tau**2)
          sigx_r = (rho/tau**2)/(4*sigx)
          gradx_t =-(d/sigx)**2*(sigx_t/sigx)
          gradx_r =-(d/sigx)**2*(sigx_r/sigx)
          gradx_d = (d/sigx)/sigx
       end if


       sigy = sqrt(cts(i)*(1 + tau) + rho**2) / 2
       e = cts(i) - tau*mu
       zy = e / sigy

       if( zy < zmin .or. zy > zmax ) then
          zy = max(zmin,min(zy,zmax))
          grady_t = 0
          grady_r = 0
          grady_d = 0
       else
          sigy_t = cts(i)/(8*sigy)
          sigy_r = (rho/sigy)/4
          grady_t =-(e/sigy)**2*(sigy_t/sigy) - (e/sigy)*(mu/sigy)
          grady_r =-(e/sigy)**2*(sigy_r/sigy)
          grady_d =-(e/sigy)*(tau/sigy)
       end if

!       grad_t = -((e/sigy)**2*(sigy_t/sigy) + (d/sigx)**2*(sigx_t/sigx)) &
!            - (e/sigy)*(mu/sigy)
!       grad_r = -((e/sigy)**2*(sigy_r/sigy) + (d/sigx)**2*(sigx_r/sigx))
!       grad_d = (e/sigy)*(-tau/sigy) + (d/sigx)/sigx

       vec(1:2) = vec(1:2) + [gradx_t + grady_t, gradx_r + grady_r ]
       vec(i+2) = gradx_d + grady_d

    end do
    vec = vec / (2*n)

  end function linegrad2

  function linejac(tau,rho,delta) result(jac)

    use rfun

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(size(delta)+2,size(delta)+2) :: jac
    real(REAL64) :: d,e,mu,zx,zy,psix,psiy,dpsix,dpsiy, &
         sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r, &
         sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr,&
         chi2_tt, chi2_tr, chi2_td, chi2_rr, chi2_rd, chi2_dd, &
         px, py, zx2, zy2, sxt2, syt2, sxtr, sytr, sxr2, syr2, dzy_t, dpx, dpy
    integer :: i,j,n

    n = size(delta)
    jac = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       e = cts(i) - tau*mu
       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)
       call grad_sigx(sigx,pht(i),tau,rho,sigx_t,sigx_r)
       call grad_sigy(sigy,cts(i),tau,rho,sigy_t,sigy_r)
       call jac_sigx(sigx,pht(i),tau,rho,sigx_t,sigx_r,sigx_tt,sigx_tr,sigx_rr)
       call jac_sigy(sigy,cts(i),tau,rho,sigy_t,sigy_r,sigy_tt,sigy_tr,sigy_rr)
       zx = d / sigx
       zy = e / sigy

       if( robust ) then
          psix = huber(zx)
          psiy = huber(zy)
          dpsix = dhuber(zx)
          dpsiy = dhuber(zy)
       else
          psix = zx
          psiy = zy
          dpsix = 1
          dpsiy = 1
       end if

       px = psix*zx
       py = psiy*zy

       zx2 = zx**2
       zy2 = zy**2
       sxt2 = (sigx_t/sigx)**2
       syt2 = (sigy_t/sigy)**2
       sxtr = (sigx_t/sigx)*(sigx_r/sigx)
       sytr = (sigy_t/sigy)*(sigy_r/sigy)
       sxr2 = (sigx_r/sigx)**2
       syr2 = (sigy_r/sigy)**2
       dzy_t = mu/sigy + zy*(sigy_t/sigy)
       dpx = dpsix*zx + psix
       dpy = dpsiy*zy + psiy
       chi2_tt = -dpsix*zx2*sxt2 - 2*px*sxt2 + px*(sigx_tt/sigx) &
            -dpsiy*dzy_t*zy*(sigy_t/sigy) &
            - psiy*(mu/sigy + 2*zy*(sigy_t/sigy))*(sigy_t/sigy) &
            + py*(sigy_tt/sigy) &
            -(dpsiy*dzy_t + psiy*(sigy_t/sigy))*(mu/sigy)
!            -dpsiy*dzy_t*(mu/sigy) &
!            -psiy*(mu/sigy)*(sigy_t/sigy)
!               + dpsiy*zy2*syt2 + 2*py*syt2 - py*(sigy_tt/sigy) &
!            - dpsiy*(mu/sigy)**2 - (dpsiy*zy + psiy)*(sigy_t/sigy)*(mu/sigy)
!               + (dpsiy*(1-zy*(sigy_t/sigy))-psiy)*(mu/sigy)**2
!       chi2_tt = dpsix*zx2*sxt2 + dpsiy*zy2*syt2 + dpsiy*(mu/sigy)**2
!       chi2_tr = -dpsix*zx2*sxtr + 2*px*sxtr - px*(sigx_tr/sigx) &
       chi2_tr =-dpx*zx*(sigx_t/sigx)*(sigx_r/sigx) &
               - dpy*dzy_t*(sigy_r/sigy) &
               - px*(sigx_t/sigx)*(sigx_r/sigx) + px*(sigx_tr/sigx) &
               - py*(sigy_t/sigy)*(sigy_r/sigy) + py*(sigy_tr/sigy)

!               + dpsiy*zy2*sytr + 2*py*sytr - py*(sigy_tr/sigy) &
!               - psiy*(mu/sigy)*(sigy_r/sigy)
!       chi2_td =-px*(sigx_t/sigx) + py*tau*(sigy_t/sigy) &
!            + (dpsiy - psiy)*(mu/sigy)*(tau/sigy)
       chi2_td =-dpx*(sigx_t/sigx)/sigx &
               + dpsiy*dzy_t*tau/sigy + psiy*(sigy_t/sigy)*tau/sigy - psiy/sigy
       chi2_rr = dpsix*zx2*sxr2 + (2*sxr2 - sigx_rr/sigx)*px &
               + dpsiy*zy2*syr2 + (2*syr2 - sigy_rr/sigy)*py
       chi2_rd =-dpx*(sigx_r/sigx)/sigx + dpy*(sigy_r/sigy)*tau/sigy
       chi2_dd = dpsix/sigx**2 + dpsiy*(tau/sigy)**2

       jac(1,1:2) = jac(1,1:2) - [ chi2_tt, chi2_tr]
       jac(2,2) = jac(2,2) + chi2_rr
       j = i + 2
       jac(1,j) = chi2_td
       jac(2,j) = chi2_rd
       jac(j,j) = chi2_dd

    end do

    jac = jac / (2*n)

    do i = 1, size(jac,1)
       do j = 1, i - 1
          jac(i,j) = jac(j,i)
       end do
    end do

  end function linejac

  subroutine sigmas(tau,rho,x,y,sigx,sigy)

    real(REAL64), intent(in) :: tau,rho,x,y
    real(REAL64), intent(out) :: sigx,sigy

    sigx = sqrt(x*(1+1/tau)+(rho/tau)**2) / 2
    sigy = sqrt(y*(1+tau)+rho**2) / 2

  end subroutine sigmas

  subroutine grad_sigx(sigx,x,tau,rho,sigx_t,sigx_r)

    real(REAL64), intent(in) :: sigx,x,tau,rho
    real(REAL64), intent(out) :: sigx_t,sigx_r
    real(REAL64) :: tau2

    tau2 = tau**2
    sigx_t = -(x*tau + 2*rho**2)/sigx / (8*tau*tau2)
    sigx_r = (rho/sigx) / (4*tau2)

  end subroutine grad_sigx

  subroutine grad_sigy(sigy,y,tau,rho,sigy_t,sigy_r)

    real(REAL64), intent(in) :: sigy,y,tau,rho
    real(REAL64), intent(out) :: sigy_t,sigy_r

    sigy_t = (y / sigy) / 8
    sigy_r = (rho / sigy) / 4

  end subroutine grad_sigy

  subroutine jac_sigx(sigx,x,tau,rho,sigx_t,sigx_r,sigx_tt,sigx_tr,sigx_rr)

    real(REAL64), intent(in) :: sigx,x,tau,rho,sigx_t,sigx_r
    real(REAL64), intent(out) :: sigx_tt,sigx_tr,sigx_rr
    real(REAL64) :: w, rs, ts, t2, t3

    t2 = tau**2
    t3 = t2*tau
    w = x*tau + 2*rho**2
    rs = rho / sigx
    ts = sigx_t/sigx
    sigx_tt = (-x/sigx + (w/sigx)*(3/tau + ts))/(8*t3)
    sigx_tr =-rs*(1 + (tau/2)*ts)/(2*t3)
    sigx_rr = (1/sigx - rs*(sigx_r/sigx))/(4*t2)

  end subroutine jac_sigx

  subroutine jac_sigy(sigy,y,tau,rho,sigy_t,sigy_r,sigy_tt,sigy_tr,sigy_rr)

    real(REAL64), intent(in) :: tau,rho,sigy,y,sigy_t,sigy_r
    real(REAL64), intent(out) :: sigy_tt,sigy_tr,sigy_rr
    real(REAL64) :: ys, rs

    ys = (y / sigy) / 8
    rs = sigy_r / sigy
    sigy_tt = -ys*(sigy_t/sigy)
    sigy_tr = -ys*rs
    sigy_rr = (1/sigy - (rho/sigy)*rs)/4

  end subroutine jac_sigy


  function linegrad1(tau,rho,delta) result(vec)

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(size(delta)+2) :: vec
    real(REAL64), dimension(3) :: gradx, grady

    real(REAL64) :: mu,d,r,s,sx,sy,sigx,sigy,rx,ry,rho2,tt, sig2
    integer :: i,n

    n = size(delta)

    tt = tau*(tau + 1)
    rho2 = rho**2

    vec = 0
    do i = 1, n

       sx = sqrt(pht(i)*(1+1/tau) + (rho/tau)**2)
       sy = sqrt(cts(i)*(1+tau) + rho2)
       s = sx
!       if( delta(i)/s < rmin ) then
!          d = rmin*s
!       else if( delta(i)/s > rmax ) then
!          d = rmax*s
!       else
!          d = delta(i)
!       end if
       d = delta(i)

       mu = pht(i) + d
       sigy = sqrt(mu*tt + rho2) / 2
       sigx = sigy / tau

       r = cts(i) - mu*tau
!       s = sqrt(cts(i)*(1+tau) + rho2)
       s = sigy
       if( r/s < rmin ) then
          r = rmin*s
       else if( r/s > rmax ) then
          r = rmax*s
       end if
!       r = cts(i) - mu*tau

       rx = d / sigx
       ry = r / sigy
       gradx = xline_grad(rx,tau,rho,d,mu,sigx)
       grady = yline_grad(ry,tau,rho,r,mu,sigy)

       vec(1:2) = vec(1:2) + (gradx(1:2) + grady(1:2))
       vec(i+2) = gradx(3) + grady(3)

    end do
    vec = vec / n

  end function linegrad1

  function linegrad(tau,rho,delta) result(vec)

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(size(delta)+2) :: vec
    real(REAL64), dimension(3) :: gradx, grady

    real(REAL64) :: mu,d,r,sigx,sigy,rx,ry,rho2,tt, sig2
    integer :: i,n

    n = size(delta)

    tt = tau*(tau + 1)
    rho2 = rho**2

    vec = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - mu*tau

       sig2 = mu*tt + rho2

       if( sig2 > 0 ) then
!       if( .true. ) then

          sigy = sqrt(sig2) / 2
          sigx = sigy / tau
          rx = d / sigx
          ry = r / sigy

          if( robust ) then

!             if( rx < rmin ) then
             rx = max(rmin,min(rx,rmax))

!                tt1 = 1 + 1/tau
!                call quadra(1/rmin**2,-tt1,-pht(i)*tt1-(rho/tau)**2,d1,d2)
!                d = d1
!                mu = pht(i) + d
!                sigy = sqrt(sig2) / 2
!                sigx = sigy / tau


!             rx = max(rmin,rx)
             gradx = xline_grad(rx,tau,rho,d,mu,sigx)

             ry = max(rmin,min(ry,rmax))
             !             ry = max(rmin,ry)
             if( sqrt(ry**2 + rx**2) > rmax ) then
                gradx = xline_grad(-rmax,tau,rho,-rmax*sigx,mu,sigx)
                grady = yline_grad(-rmax,tau,rho,-rmax*sigy,mu,sigy)
             else
                gradx = xline_grad(rx,tau,rho,d,mu,sigx)
                grady = yline_grad(ry,tau,rho,r,mu,sigy)
             end if

          else

             !rx = max(min(rx,rmax),rmin)
             gradx = xline_grad(rx,tau,rho,d,mu,sigx)

             !          ry = max(min(ry,rmax),rmin)
             grady = yline_grad(ry,tau,rho,r,mu,sigy)

          end if

       else ! sig2 <= 0
!          sig2 = rho2
!          gradx = grad_xlike(rmin,tau,rho,d,mu,sx,sigx)
          !          grady = grad_ylike(rmin,tau,rho,r,mu,sy,sigy)
          gradx = 0
          grady = 0
       end if

       vec(1:2) = vec(1:2) + (gradx(1:2) + grady(1:2))
       vec(i+2) = gradx(3) + grady(3)

    end do
    vec = vec / n

  end function linegrad

  function xline_grad(rx,tau,rho,d,mu,sigx) result(grad)

    real(REAL64), intent(in) :: tau,rho,d,mu,rx,sigx
    real(REAL64), dimension(3) :: grad
    real(REAL64) :: sigx_t, sigx_r, sigx_d, q

!    call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)

    q = 1 - rx * (d / sigx)
    grad(1) = q*(sigx_t/sigx)
    grad(2) = q*(sigx_r/sigx)
    grad(3) = q*(sigx_d/sigx) + rx/sigx

  end function xline_grad

  function yline_grad(ry,tau,rho,r,mu,sigy) result(grad)

    real(REAL64), intent(in) :: tau,rho,r,mu,ry,sigy
    real(REAL64), dimension(3) :: grad
    real(REAL64) :: sigy_t, sigy_r, sigy_d, rs, q

!    call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)

    q = 1 - ry * (r/sigy)
    grad(1) = q*(sigy_t/sigy) - ry*(mu/sigy)
    grad(2) = q*(sigy_r/sigy)
    grad(3) = q*(sigy_d/sigy) - ry*tau/sigy

  end function yline_grad

  function disp(tau,rho,delta)

    use rfun

    real(REAL64) :: disp
    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64) :: s, d, mu, r, zx, zy, sigx, sigy, ipsi!, rho2, tt
    integer :: i,n

    n = size(delta)
!    tt = tau*(tau + 1)
!    rho2 = rho**2

    s = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
!       sigy = sqrt(mu*tt + rho2) / 2
!       sigx = sigy / tau
!       sigx = sqrt(pht(i)*(1 + 1/tau) + (rho/tau)**2) / 2
       !       sigy = sqrt(cts(i)*(1 + tau) + rho**2) / 2

       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)
       zx = d / sigx
       zy = r / sigy
       if( robust ) then

!!$          if( zx < rmin ) then
!!$             block
!!$               real(REAL64) :: a,b,c,det,d1,d2,d12
!!$
!!$               a = (tau/rmin)**2
!!$               b = -tt
!!$               c = -(pht(i)*tt + rho2)
!!$               det = b**2 - 4*a*c
!!$               if( det > 0 ) then
!!$                  d12 = sqrt(det)
!!$                  d1 = -(b + d12)/(2*a)
!!$                  d2 = -(b - d12)/(2*a)
!!$                  write(*,*) i,real(zx),real(d),real(d1),real(d2),&
!!$                       real(rmin*sqrt(pht(i)*tt+rho2)/tau)
!!$               else if( abs(det) < epsilon(det) ) then
!!$                  d1 = -b/(2*a)
!!$                  d2 = d1
!!$                  write(*,*) d,d1,d2
!!$               else
!!$                  write(*,*) d,'imag'
!!$               end if
!!$               stop
!!$             end block
!!$          end if

!          write(*,*) real(sigx),real(sigy)
!          write(*,*) real(sigx),real(sigy)
!          zx = d / sigx
!          zy = (cts(i) - tau*mu) / sigy

!          write(*,*) real(zx),real(zy)
          zx = max(zmin,min(zx,zmax))
          zy = max(zmin,min(zy,zmax))
          ipsi = (zx**2 + zy**2)/2


!          ipsi = itukey(zx) + itukey(zy)

!          ipsi = 0
!          if( abs(tukey(zx)) > epsilon(zx) .and. abs(tukey(zy)) > epsilon(zy) ) then

!          if( abs(tukey(zx)) > epsilon(zx) ) ipsi = ipsi + itukey(zx)
!          if( abs(tukey(zy)) > epsilon(zy) ) ipsi = ipsi + itukey(zy)

!             ipsi = itukey(zx) + itukey(zy)
!          end if
!          rho = max(min(z,rmax),rmin)**2/2
       else
          ipsi = (zx**2 + zy**2)/2
       end if

       s = s + ipsi !+ log(sigx*scale) + log(sigy*scale)

!       if(

!       s = s + cfun(max(min(zx,rmax),rmin)) +

    end do

    disp = s / n

  end function disp

  function grad_disp(tau,rho,delta) result(grad)

    use rfun

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(size(delta)+2) :: grad
    real(REAL64), dimension(3) :: xgrad, ygrad
    real(REAL64) :: s, d, e, mu, r, zx, zy, sigx, sigy, rho2, tt, &
         sigx_t,sigx_r,sigx_d,sigy_t,sigy_r,sigy_d, psix, psiy, grad_t, &
         grad_r, grad_d, gradx_t, gradx_r, gradx_d, grady_t, grady_r, grady_d
    integer :: i,n

    n = size(delta)
    tt = tau*(tau + 1)
    rho2 = rho**2

    grad = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
!       sigx = sqrt(pht(i)*(1 + 1/tau) + (rho/tau)**2) / 2
!       sigy = sqrt(cts(i)*(1 + tau) + rho**2) / 2
!       sigy = sqrt(mu*tt + rho2) / 2
!       sigx = sigy / tau

!       zx = d / sigx
!       zy = r / sigy

       if( robust ) then

!!$          if( rmin < zx .and. zx < rmax ) then
!!$!             call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!!$          else
!!$             zx = max(rmin,min(zx,rmax))
!!$             sigx_t = 0
!!$             sigx_r = 0
!!$             sigx_d = 0
!!$          end if
!!$!          call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!!$
!!$          if( rmin < zy .and. zy < rmax ) then
!!$!             call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)
!!$          else
!!$             zy = max(rmin,min(zy,rmax))
!!$             sigy_t = 0
!!$             sigy_r = 0
!!$             sigy_d = 0
!!$          end if
!!$!          call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)
!!$
!!$          xgrad = xgradlike(zx,sigx,sigx_t,sigx_r,sigx_d)
!!$          ygrad = ygradlike(tau,mu,zy,sigy,sigy_t,sigy_r,sigy_d)

          sigx = sqrt(pht(i)*(1 + 1/tau) + (rho/tau)**2) / 2
          zx = d / sigx
          if( zx < zmin .or. zx > zmax ) then
             zx = max(zmin,min(zx,zmax))
             gradx_t = 0
             gradx_r = 0
             gradx_d = 0
          else
             sigx_t = -(pht(i) + 2*rho**2/tau)/(8*sigx*tau**2)
             sigx_r = (rho/tau**2)/(4*sigx)
             gradx_t =-(d/sigx)**2*(sigx_t/sigx)
             gradx_r =-(d/sigx)**2*(sigx_r/sigx)
             gradx_d = (d/sigx)/sigx
          end if

          sigy = sqrt(cts(i)*(1 + tau) + rho**2) / 2
          e = cts(i) - tau*mu
          zy = e / sigy
          if( zy < zmin .or. zy > zmax ) then
             zy = max(zmin,min(zy,zmax))
             grady_t = 0
             grady_r = 0
             grady_d = 0
          else
             sigy_t = cts(i)/(8*sigy)
             sigy_r = rho/sigy/4
             grady_t = (e/sigy)*(-mu/sigy) - (e/sigy)**2*(sigy_t/sigy)
             grady_r =-(e/sigy)**2*(sigy_r/sigy)
             grady_d =-(e/sigy)*(tau/sigy)
          end if

          grad_t = (e/sigy)*(-mu/sigy) - (e/sigy)**2*(sigy_t/sigy) &
               - (d/sigx)**2*(sigx_t/sigx)
          grad_r = -((e/sigy)**2*(sigy_r/sigy) + (d/sigx)**2*(sigx_r/sigx))
          grad_d = (e/sigy)*(-tau/sigy) + (d/sigx)/sigx

          grad(1) = grad(1) + (gradx_t + grady_t)
          grad(2) = grad(2) + (gradx_r + grady_r)
          grad(i+2) = gradx_d + grady_d

!          grad(1) = grad(1) + xgrad(1) + ygrad(1)
!          grad(2) = grad(2) + xgrad(2) + ygrad(2)
!          grad(i+2) = xgrad(3) + ygrad(3)


!!$          psix = tukey(d/sigx)
!!$          psiy = tukey(r/sigy)
!!$
!!$          if( abs(psix) > epsilon(psix) .and. abs(psiy) > epsilon(psiy) ) then
!!$
!!$             grad(1) = grad(1) + psix*zx*(sigx_t/sigx)+ psiy*zy*(sigy_t/sigy) &
!!$                  + psiy*(mu/sigy)
!!$             grad(2) = grad(2) + psix*zx*(sigx_r/sigx)+ psiy*zy*(sigy_r/sigy)
!!$             grad(i+2) = psix*zx*(sigx_d/sigx)+ psiy*zy*(sigy_d/sigy) &
!!$                  -psix/sigx +psiy*(tau/sigy)
!!$          end if
       else

!          call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!          call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)

!          xgrad = xgradlike(zx,sigx,sigx_t,sigx_r,sigx_d)
!          ygrad = ygradlike(tau,mu,zy,sigy,sigy_t,sigy_r,sigy_d)

!          grad(1) = grad(1) + xgrad(1) + ygrad(1)
!          grad(2) = grad(2) + xgrad(2) + ygrad(2)
!          grad(i+2) = xgrad(3) + ygrad(3)

!          grad(1) = grad(1) + (d/sigx)**2*(sigx_t/sigx) + (r/sigy)**2*(sigy_t/sigy) &
!               + (mu/sigy)*(r/sigy)
!          grad(2) = grad(2) + (d/sigx)**2*(sigx_r/sigx) + (r/sigy)**2*(sigy_r/sigy)
!          grad(i+2) = (d/sigx)**2*(sigx_d/sigx) + (r/sigy)**2*(sigy_d/sigy) &
          !               - (d/sigx)/sigx + (tau/sigy)*(r/sigy)

          sigx = sqrt(pht(i)*(1 + 1/tau) + (rho/tau)**2) / 2
          zx = d / sigx
          sigx_t = -(pht(i) + 2*rho**2/tau)/(8*sigx*tau**2)
          sigx_r = (rho/tau**2)/(4*sigx)
          gradx_t =-(d/sigx)**2*(sigx_t/sigx)
          gradx_r =-(d/sigx)**2*(sigx_r/sigx)
          gradx_d = (d/sigx)/sigx

          sigy = sqrt(cts(i)*(1 + tau) + rho**2) / 2
          e = cts(i) - tau*mu
          zy = e / sigy
          sigy_t = cts(i)/(8*sigy)
          sigy_r = rho/sigy/4
          grady_t = (e/sigy)*(-mu/sigy) - (e/sigy)**2*(sigy_t/sigy)
          grady_r =-(e/sigy)**2*(sigy_r/sigy)
          grady_d =-(e/sigy)*(tau/sigy)

          grad_t = (e/sigy)*(-mu/sigy) - (e/sigy)**2*(sigy_t/sigy) &
               - (d/sigx)**2*(sigx_t/sigx)
          grad_r = -((e/sigy)**2*(sigy_r/sigy) + (d/sigx)**2*(sigx_r/sigx))
          grad_d = (e/sigy)*(-tau/sigy) + (d/sigx)/sigx

          grad(1) = grad(1) + (gradx_t + grady_t)
          grad(2) = grad(2) + (gradx_r + grady_r)
          grad(i+2) = gradx_d + grady_d

       end if
    end do

!    if( robust ) then
!       grad = -grad / (2*n)
!    else
!       grad = -grad / n
!    end if

    grad = -grad / (2*n)

  end function grad_disp

  function xgradlike(zx,sigx,sigx_t,sigx_r,sigx_d) result(grad)

    real(REAL64), dimension(3) :: grad
    real(REAL64), intent(in) :: zx,sigx,sigx_t,sigx_r,sigx_d
    real(REAL64) :: ds,ds2

    ds = zx
    ds2 = ds**2
    grad(1) = ds2*(sigx_t/sigx)
    grad(2) = ds2*(sigx_r/sigx)
    grad(3) = ds2*(sigx_d/sigx) - ds/sigx

  end function xgradlike

  function ygradlike(tau,mu,zy,sigy,sigy_t,sigy_r,sigy_d) result(grad)

    real(REAL64), dimension(3) :: grad
    real(REAL64), intent(in) :: tau,mu,zy,sigy,sigy_t,sigy_r,sigy_d
    real(REAL64) :: rs,rs2

    rs = zy
    rs2 = rs**2
    grad(1) = rs2*(sigy_t/sigy) + rs*(mu/sigy)
    grad(2) = rs2*(sigy_r/sigy)
    grad(3) = rs2*(sigy_d/sigy) + rs*(tau/sigy)

  end function ygradlike


  subroutine disptest(tau,rho,delta)

    use minpacks

    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in) :: tau,rho

    real(REAL64), dimension(size(delta)+2) :: grad,eps,z,dc
    real(REAL64), dimension(-1:1) :: c
    integer :: i, j, l

    write(error_unit,*) 'disptest():'

    eps = epsilon(eps)**(1.0/3.0)

    grad = grad_disp(tau,rho,delta)
    do i = 1, size(grad)
       do l = -1,1,2
          z(1) = tau
          z(2) = rho
          z(3:) = delta
          z(i) = z(i) + l*eps(i)
          c(l) = disp(z(1),z(2),z(3:))
       end do
       dc(i) = (c(1) - c(-1)) / (2*eps(i))
       write(error_unit,*) i,real(grad(i)),real(dc(i)),real(dc(i)-grad(i))
    end do

  end subroutine disptest


!!$  real(REAL64) function xdisp(tau,delta,sigv,sx)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: tau, sigv, sx
!!$    real(REAL64), dimension(:), intent(in) :: delta
!!$    real(REAL64) :: d, s, mu, sigv2, rx, sigx, f
!!$    integer :: i,n
!!$
!!$    n = size(pht)
!!$    sigv2 = sigv**2
!!$
!!$    s = 0
!!$    do i = 1, n
!!$       d = delta(i)
!!$       mu = pht(i) + d
!!$       sigx = sqrt(mu + sigv2)
!!$       rx = (d/sigx)/sx
!!$       if( robust ) then
!!$          f = 2*itukey(rx)
!!$       else
!!$          f = rx**2
!!$       end if
!!$       s = s + f
!!$    end do
!!$    xdisp = scale - s/n
!!$
!!$  end function xdisp
!!$
!!$  real(REAL64) function ydisp(tau,delta,sigu,sy)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: tau, sigu, sy
!!$    real(REAL64), dimension(:), intent(in) :: delta
!!$    real(REAL64) :: r, s, mu, sigu2, ry, sigy, f
!!$    integer :: i,n
!!$
!!$    n = size(cts)
!!$    sigu2 = sigu**2
!!$
!!$    s = 0
!!$    do i = 1, n
!!$       mu = pht(i) + delta(i)
!!$       r = cts(i) - tau*mu
!!$       sigy = sqrt(mu*tau + sigu2)
!!$       ry = (r/sigy)/sy
!!$       if( robust ) then
!!$          f = 2*itukey(ry)
!!$       else
!!$          f = ry**2
!!$       end if
!!$       s = s + f
!!$    end do
!!$    ydisp = scale - s/n
!!$
!!$  end function ydisp
!!$
!!$
!!$  function disp(tau,rho,phi,delta,sx,sy) result(c)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: tau, rho, phi, sx,sy
!!$    real(REAL64), dimension(:), intent(in) :: delta
!!$    real(REAL64), dimension(2) :: c
!!$    real(REAL64) :: r, d, mu, sigv, sigu, sigu2, sigv2, rx, ry, sigx, sigy, &
!!$         fx, fy
!!$    integer :: i,n
!!$
!!$    n = size(cts)
!!$    sigu = rho*sin(phi)
!!$    sigv = rho*(cos(phi)/tau)
!!$    sigu2 = sigu**2
!!$    sigv2 = sigv**2
!!$
!!$    c = 0
!!$    do i = 1, n
!!$       d = delta(i)
!!$       mu = pht(i) + d
!!$       r = cts(i) - tau*mu
!!$       sigx = sqrt(mu + sigv2)
!!$       sigy = sqrt(mu*tau + sigu2)
!!$       rx = (d/sigx)/sx
!!$       ry = (r/sigy)/sy
!!$       if( robust ) then
!!$          fx = 2*itukey(rx)
!!$          fy = 2*itukey(ry)
!!$       else
!!$          fx = rx**2
!!$          fy = ry**2
!!$       end if
!!$       c = c + [ fx, fy ]
!!$    end do
!!$    c = scale - c/n
!!$
!!$  end function disp


  function gradc(tau,rho,phi,delta,sx,sy) result(dc)

    use rfun

    real(REAL64), intent(in) :: tau, rho, phi, sx,sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(:,:), allocatable :: dc
    real(REAL64) :: r, d, mu, sinfi, cosfi, sigu, sigv, sigu2, sigv2, rx, ry, &
         sigx, sigy, psix, psiy, dx, dy, sigx_t, sigy_t, sigx_r, sigy_r, &
         sigx_f, sigy_f, sigx_d, sigy_d
    integer :: i,n

    n = size(cts)
    allocate(dc(2,3+n))

    sinfi = sin(phi)
    cosfi = cos(phi)
    sigu = rho*sinfi
    sigv = rho*(cosfi/tau)
    sigu2 = sigu**2
    sigv2 = sigv**2

    dc = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
       sigx = sqrt(mu + sigv2)
       sigy = sqrt(mu*tau + sigu2)
       rx = (d/sigx)/sx
       ry = (r/sigy)/sy
       if( robust ) then
          psix = tukey(rx)
          psiy = tukey(ry)
       else
          psix = rx
          psiy = ry
       end if
       dx = (psix/sx)*(d/sigx)
       dy = (psiy/sy)*(r/sigy)

       sigy_t = (mu / sigy) / 2
       sigx_t = -(cosfi**2/tau**3)*(rho/sigx) * rho
       sigy_r = (rho/sigy)*sinfi**2
       sigx_r = (rho/sigx)*(cosfi/tau)**2
       sigy_f = (sinfi*cosfi) * (rho/sigy) * rho
       sigx_f = -(sinfi*cosfi)/tau**2 * (rho/sigx) * rho
!       dt = 2*(dx*(sigx_t/sigx) + dy*(sigy_t/sigy))
!       dr = 2*(dx*(sigx_r/sigx) + dy*(sigy_r/sigy))
!       df = 2*(dx*(sigx_f/sigx) + dy*(sigy_f/sigy))
       dc(1,1:3) = dc(1,1:3) + 2*dx*[ sigx_t/sigx, sigx_r/sigx, sigx_f/sigx ]
       dc(2,1:3) = dc(2,1:3) + 2*dy*[ sigy_t/sigy, sigy_r/sigy, sigy_f/sigy ]

       sigx_d = 1/(2*sigx)
       sigy_d = (tau/2) / sigy
!       dd = 2*(dx*(sigx_d/sigx) + dy*(sigy_d/sigy))
       dc(1,3+i) = 2*dx*(sigx_d/sigx)
       dc(2,3+i) = 2*dy*(sigy_d/sigy)
    end do
    dc = dc / n

  end function gradc


  function linescale(tau,rho,delta,sx,sy) result(s)

    use rfun

    real(REAL64), intent(in) :: tau, rho, sx,sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(2) :: s

    real(REAL64) :: r, d, mu, rx, ry, sigx, sigy, fx, fy, sigv2, rho2, tt
    integer :: i,n

    n = size(delta)
    rho2 = rho**2
    tt = tau*(tau + 1)
!    sigv2 = ((rho/tau)*cos(phi))**2

    s = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
!       sigx = sqrt(mu + sigv2) / 2
       sigy = sqrt(mu*tt + rho2) / 2
       sigx = sigy / tau
       rx = (d/sigx) / sx
       ry = (r/sigy) / sy
       if( robust ) then
          fx = itukey(rx)
          fy = itukey(ry)
!          fx = max(min(rx,rmax),rmin)**2/2
!          fy = max(min(ry,rmax),rmin)**2/2
       else
          fx = rx**2/2 + log(sx*sigx)
          fy = ry**2/2 + log(sy*sigy)
       end if
!       fx = rx**2/2
       s = s + [ fx, fy ]
    end do
!10  continue

    s = 2*s/n

  end function linescale


  function linedisp(tau,rho,delta,sx,sy) result(c)

    use rfun

    real(REAL64), intent(in) :: tau, rho, sx,sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(2) :: c, s2

    real(REAL64) :: r, d, mu, rx, ry, sigx, sigy, fx, fy, rho2
    integer :: i,n

    c = scale - linescale(tau,rho,delta,sx,sy)
    return

    n = size(cts)
    rho2 = rho**2

!    sigu = rho*sin(phi)
!    sigv = (rho/tau)*cos(phi)
!    sigu2 = sigu**2
!    sigv2 = sigv**2
!    qx = sqrt2 / sx
!    qy = sqrt2 / sy

!!$    c = 0
!!$    do i = 1, n
!!$       d = delta(i)
!!$       mu = pht(i) + d
!!$       r = cts(i) - tau*mu
!!$!       sigx = sqrt(mu*(1 + 1/tau) + (rho/tau)**2)/2
!!$       !       sigy = sqrt(mu*tau*(1 + tau) + rho**2)/2
!!$       sigx = sqrt(mu + sigv2)
!!$       sigy = sqrt(mu*tau + sigu2)
!!$       rx = (d/sigx)/sx
!!$       ry = (r/sigy)/sy
!!$       if( robust ) then
!!$          fx = 2*itukey(rx)
!!$          fy = 2*itukey(ry)
!!$       else
!!$          fx = rx**2
!!$          fy = ry**2
!!$       end if
!!$       c = c + [ fx, fy ]
!!$    end do

!    goto 10
    s2 = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
!       sigx = sqrt(mu + sigv2) / 2
       sigy = sqrt(mu*tau*(tau+1) + rho2) / 2
       sigx = sigy / tau
!       rx = qx*(d/sigx)
!       ry = qy*(r/sigy)!/sqrt(tau)
       rx = (d/sigx) / sx
       ry = (r/sigy) / sy
       if( robust ) then
          fx = itukey(rx)
          fy = itukey(ry)
       else
          fx = rx**2/2
          fy = ry**2/2
       end if
!       fx = rx**2/2
       s2 = s2 + [ fx, fy ]
    end do
!10  continue

!    c = sqrt(s2/(n-1))
    c = scale - 2*s2/n


  end function linedisp


  function linegradc(tau,rho,delta,sx,sy) result(dc)

    use rfun

    real(REAL64), intent(in) :: tau, rho, sx,sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(:,:), allocatable :: dc

    real(REAL64) :: r, d, mu, rx, ry, sigx, sigy, psix, psiy, dx, dy, &
         sigx_t, sigy_t, sigx_r, sigy_r, sigx_d, sigy_d!, &
 !        cosfi, sinfi, sigu, sigv, sigu2, sigv2, q, qx, qy, c2, s2, cs
    integer :: i,n

    n = size(cts)
    allocate(dc(2,n+1))

!    sinfi = sin(phi)
!    cosfi = cos(phi)
!    c2 = cosfi**2
!    s2 = sinfi**2
!    cs = cosfi*sinfi
!    sigu = rho*sinfi
!    sigv = (rho/tau)*cosfi
!    sigu2 = sigu**2
!    sigv2 = sigv**2
!    qx = sqrt2 / sx
!    qy = sqrt2 / sy

!!$    dc = 0
!!$    do i = 1, n
!!$       d = delta(i)
!!$       mu = pht(i) + d
!!$       r = cts(i) - tau*mu
!!$!       sigx = sqrt(mu*(1 + 1/tau) + (rho/tau)**2)/2
!!$!       sigy = sqrt(mu*tau*(1 + tau) + rho**2)/2
!!$       sigx = sqrt(mu + sigv2)
!!$       sigy = sqrt(mu*tau + sigu2)
!!$       rx = (d/sigx)/sx
!!$       ry = (r/sigy)/sy
!!$       if( robust ) then
!!$          psix = tukey(rx)
!!$          psiy = tukey(ry)
!!$       else
!!$          psix = rx
!!$          psiy = ry
!!$       end if
!!$       dx = (psix/sx)*(d/sigx)
!!$       dy = (psiy/sy)*(r/sigy)
!!$
!       sigy_t = (mu / sigy) / 2
!       sigx_t = -(cosfi**2/tau**3)*(rho/sigx) * rho
!       sigy_r = (rho/sigy)*sinfi**2
!       sigx_r = (rho/sigx)*(cosfi/tau)**2
!       sigy_f = (sinfi*cosfi) * (rho/sigy) * rho
!       sigx_f = -(sinfi*cosfi)/tau**2 * (rho/sigx) * rho
!       dt = 2*(dx*(sigx_t/sigx) + dy*(sigy_t/sigy))
!       dr = 2*(dx*(sigx_r/sigx) + dy*(sigy_r/sigy))
!       df = 2*(dx*(sigx_f/sigx) + dy*(sigy_f/sigy))
!!$       sigx_t = -(mu/sigx + (2*rho/tau)*(rho/sigx))/(2*tau**2)/4
!!$       sigy_t = (2*tau + 1)/2 * (mu/sigy)/4
!!$       sigx_r = (rho/sigx)/tau**2/4
!!$       sigy_r = rho / sigy/4
!!$       sigx_d = (1 + 1/tau)/2 / sigx/2
!!$       sigy_d = tau*(1 + tau)/2/sigy/2
!!$
!!$       dc(1,1:2) = dc(1,1:2) + 2*dx*[ sigx_t/sigx, sigx_r/sigx ]
!!$       !       dc(2,1:2) = dc(2,1:2) + 2*dy*[ -mu/sigy-(r/sigy)*(sigy_t/sigy), sigy_r/sigy ]
!!$       dc(2,2) = dc(2,2) + 2*dy*sigy_r/sigy
!!$       dc(2,1) = dc(2,1) + 2*(psiy/sy)*(mu/sigy + (r/sigy)*(sigy_t/sigy))
!!$
!!$!       sigx_d = 1/(2*sigx)/2
!!$!       sigy_d = (tau/2) / sigy/2
!!$
!!$!       sigx_d = (1 + 1/tau)/2 / sigx/2
!!$!       sigy_d = tau*(1 + tau)/2/sigy/2
!!$
!!$!       dd = 2*(dx*(sigx_d/sigx) + dy*(sigy_d/sigy))
!!$       dc(1,2+i) = -2*(psix/sx)*(1/sigx - (d/sigx)*(sigx_d/sigx))
!!$       dc(2,2+i) = 2*(psiy/sy)*(tau/sigy + (r/sigy)*(sigy_d/sigy))
!!$
!!$
!!$       sigy_t = (mu / sigy) / 2
!!$       sigx_t = -(cosfi**2/tau**3)*(rho/sigx) * rho
!!$       sigy_r = (rho/sigy)*sinfi**2
!!$       sigx_r = (rho/sigx)*(cosfi/tau)**2
!!$       sigy_f = (sinfi*cosfi) * (rho/sigy) * rho
!!$       sigx_f = -(sinfi*cosfi)/tau**2 * (rho/sigx) * rho
!!$       sigx_d = 1/(2*sigx)
!!$       sigy_d = (tau/2) / sigy
!!$
!!$       dc(1,1:3) = dc(1,1:3) + 2*dx*[ sigx_t/sigx, sigx_r/sigx, sigx_f/sigx ]
!!$       dc(1,3+i) = -2*(psix/sx)*(1/sigx - (d/sigx)*(sigx_d/sigx))
!!$       dc(2,1) = dc(2,1) + 2*(psiy/sy)*(mu/sigy + (r/sigy)*(sigy_t/sigy))
!!$       dc(2,2:3) = dc(2,2:3) + 2*dy*[ sigy_r/sigy, sigy_f/sigy ]
!!$       dc(2,3+i) = 2*(psiy/sy)*(tau/sigy + (r/sigy)*(sigy_d/sigy))
!!$
!!$    end do


!    goto 20
!    q = sqrt2
!    q = 2
    dc = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
!       sigx = sqrt(mu + sigv2) / 2
!       sigy = sqrt(mu*tau + sigu2)
       sigy = sqrt(mu*tau*(tau+1) + rho**2) / 2
       sigx = sigy / tau
       rx = (d/sigx) / sx
       ry = (r/sigy) / sy
       if( robust ) then
          psix = tukey(rx)
          psiy = tukey(ry)
!          psix = max(min(rx,rmax),rmin)
!          psiy = max(min(ry,rmax),rmin)
       else
          psix = rx
          psiy = ry
       end if
!       psix = rx
       dx = (psix/sx)*(d/sigx)
       dy = (psiy/sy)*(r/sigy)

!       sigy_t = (mu / sigy) / 2
!       sigx_t = -c2*(rho/sigx) * (rho / tau**3)
!       sigy_r = (rho/sigy)*s2
!       sigx_r = (rho/sigx)*(cosfi/tau)**2
!       sigy_f = cs * (rho/sigy) * rho
!       sigx_f = -cs * (rho/sigx) * (rho / tau**2)
!       sigx_d = 1/(2*sigx)
!       sigy_d = (tau/2) / sigy

!       call grad_sigx(sigx,mu,tau,rho,cosfi,sinfi,sigx_t,sigx_r,sigx_f,sigx_d)
!       call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!       sigx_f = 0
!       call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)
!       sigy_f = 0

!       dc(1,1:2) = dc(1,1:2) + dx*[ sigx_t/sigx, sigx_r/sigx ]
       dc(1,1) = dc(1,1) + dx*(sigx_t/sigx)
       dc(1,1+i) = -(psix/sx)*(1/sigx - (d/sigx)*(sigx_d/sigx))
       dc(2,1) = dc(2,1) + (psiy/sy)*(mu/sigy + (r/sigy)*(sigy_t/sigy))
!       dc(2,2) = dc(2,2) + dy* sigy_r/sigy
       dc(2,1+i) = (psiy/sy)*(tau/sigy + (r/sigy)*(sigy_d/sigy))



!!$       dc(1,1:3) = dc(1,1:3) + q*dx*[ sigx_t/sigx, sigx_r/sigx, sigx_f/sigx ]
!!$       dc(1,3+i) = -q*(psix/sx)*(1/sigx - (d/sigx)*(sigx_d/sigx))
!!$       dc(2,1) = dc(2,1) + q*(psiy/sy)*(mu/sigy + (r/sigy)*(sigy_t/sigy))
!!$!       dc(2,1) = dc(2,1) &
!!$!            + q*(psiy/sy)*(mu/sigy + (r/sigy)*(sigy_t/sigy) + (r/sigy)/2/tau)/sqrt(tau)
!!$       dc(2,2:3) = dc(2,2:3) + q*dy*[ sigy_r/sigy, sigy_f/sigy ]!/sqrt(tau)
!!$       dc(2,3+i) = q*(psiy/sy)*(tau/sigy + (r/sigy)*(sigy_d/sigy))!/sqrt(tau)

    end do
!20  continue

    dc = 2 * dc / n

  end function linegradc


  function linejacc(tau,rho,delta,sx,sy) result(dd)

    use rfun

    real(REAL64), intent(in) :: tau, rho, sx,sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(:,:,:), allocatable :: dd

    real(REAL64) :: r, d, mu, rho2, rx, ry, sigx, sigy, psix, psiy, dpsix, dpsiy, &
         dx, dy, px, py, dpx, dpy, ddx, ddy, &
         sigx_t, sigy_t, sigx_r, sigy_r, sigx_d, sigy_d, &
         sigx_tt,sigy_tt, sigx_tr,sigy_tr, sigx_td,sigy_td, &
         sigx_rr,sigy_rr, sigx_rd,sigy_rd, sigx_fd,sigy_fd, sigx_dd,sigy_dd, &
         st2, str, std, sr2, srd, sd, sdd, st, sr
!    real(REAL64) :: r, d, mu, rx, ry, sigx, sigy, psix, psiy, dx, dy, &
!         sigx_t, sigy_t, sigx_r, sigy_r, sigx_f, sigy_f, sigx_d, sigy_d, &
!         cosfi, sinfi, sigu, sigv, sigu2, sigv2, q, qx, qy, c2, s2, cs, &
!         sigx_tt,sigy_tt, sigx_tr,sigy_tr, sigx_tf,sigy_tf, sigx_td,sigy_td, &
!         sigx_rr,sigy_rr, sigx_rf,sigy_rf, sigx_rd,sigy_rd, &
!         sigx_ff,sigy_ff, sigx_fd,sigy_fd, sigx_dd,sigy_dd

    integer :: i,j,l,n

    n = size(delta)
    allocate(dd(2,n+1,n+1))
    rho2 = rho**2
!!$    sinfi = sin(phi)
!!$    cosfi = cos(phi)
!!$    c2 = cosfi**2
!!$    s2 = sinfi**2
!!$    cs = cosfi*sinfi
!!$    sigu = rho*sinfi
!!$    sigv = (rho/tau)*cosfi
!!$    sigu2 = sigu**2
!!$    sigv2 = sigv**2
!    qx = sqrt2 / sx
!    qy = sqrt2 / sy

    dd = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
       sigy = sqrt(mu*tau*(tau+1) + rho2) / 2
       sigx = sigy / tau

       rx = (d/sigx) / sx
       ry = (r/sigy) / sy
       if( robust ) then
          psix = tukey(rx)
          psiy = tukey(ry)
          dpsix = dtukey(rx)
          dpsiy = dtukey(ry)
       else
          psix = rx
          psiy = ry
          dpsix = 1
          dpsiy = 1
       end if
       dx = (psix/sx)*(d/sigx)
       dy = (psiy/sy)*(r/sigy)
       ddx = dpsix/sx**2*(d/sigx)**2
       ddy = dpsiy/sy**2*(r/sigy)**2
       px = psix/sx
       py = psiy/sy
       dpx = dpsix/sx**2
       dpy = dpsiy/sy**2

!       call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!       call jac_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d, &
!            sigx_tt,sigx_tr,sigx_td,sigx_rr,sigx_rd,sigx_dd)

       st2 = (sigx_t/sigx)**2
       dd(1,1,1) = dd(1,1,1) + ddx*st2 + dx*(2*st2 - sigx_tt/sigx)

!       write(*,*) ddx*st2, dx*2*st2, -dx* sigx_tt/sigx

!       dd(1,1,1) = dd(1,1,1) + (d/sx)**2*(-6*sigx_t**2/sigx**4+2*sigx_tt/sigx**3)

!       dd(1,1,1) = dd(1,1,1) + (dpsix/sx**2)*((d/sigx)*(sigx_t/sigx))**2 &
!            +(psix/sx)*2*(d/sigx)*(sigx_t/sigx)**2 &
!            - (psix/sx)*(d/sigx)*(sigx_tt/sigx)
       str = (sigx_t/sigx)*(sigx_r/sigx)
!       dd(1,1,2) = dd(1,1,2) + ddx*str + dx*(2*str - sigx_tr/sigx)

       !       std = (sigx_t/sigx)*(sigx_d/sigx)
!       st = sigx_t/sigx
       sd = d / sigx
       sdd = sigx_d/sigx
       std = sigx_t/sigx
       dd(1,1,i+1) = -dpx*sd*std*(1/sigx - sd*sdd) &
            + dx*(2*std*sdd - sigx_td/sigx) - (px/sigx)*std

       sr2 = (sigx_r/sigx)**2
!       dd(1,2,2) = dd(1,2,2) + ddx*sr2 + dx*(2*sr2 - sigx_rr/sigx)
       srd = (sigx_r/sigx)*(sigx_d/sigx)
       sr = sigx_r/sigx
       dd(1,2,i+1) = -dpx*sd*sr*(1/sigx - sd*sdd) &
            + dx*(2*srd - sigx_rd/sigx) - (px/sigx)*sr

       dd(1,i+1,i+1) = dpx*(1/sigx - sd*sdd)**2 &
            +px*(-(2/sigx)*sdd +2*sd*sdd**2 - sd*(sigx_dd/sigx))


!       call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)
!       call jac_sigy(sigy,mu,sigy_t,sigy_r,sigy_d, &
!            sigy_tt,sigy_tr,sigy_td,sigy_rr,sigy_rd,sigy_dd)

       st = sigy_t / sigy
       sr = r / sigy
       dd(2,1,1) = dd(2,1,1) + dpy*(mu/sigy + sr*st)**2 &
            - py*(-2*(mu/sigy)*st - 2*sr*st**2 + sr*(sigy_tt/sigy))
       str = (sigy_t/sigy)*(sigy_r/sigy)
       srd = sigy_r/sigy
!       dd(2,1,2) = dd(2,1,2) + dpy*(mu/sigy + sr*st)*sr*srd &
!            - py*(- (mu/sigy)*srd -2*sr*str + sr*(sigy_tr/sigy) )
!       dpy*((mu/sigy)*sr*(sigy_r/sigy) + sr**2*str) &


       sd = sigy_d / sigy
       dd(2,1,i+1) = dpy*(mu/sigy + sr*st)*(tau/sigy + sr*sd) &
            -py*(1/sigy - (mu/sigy)*sd - (tau/sigy)*st - 2*sr*st*sd &
            + sr*(sigy_td/sigy))

       sr2 = (sigy_r/sigy)**2
!       dd(2,2,2) = dd(2,2,2) + dpy*sr**2*sr2 &
!            +py*(2*sr*sr2 - sr*(sigy_rr/sigy))

       dd(2,2,i+1) = dpy*(tau/sigy + sr*sd)*sr*srd &
            -py*(-(tau/sigy + 2*sr*sd)*srd + sr*(sigy_rd/sigy))

       dd(2,i+1,i+1) = dpy*(tau/sigy + sr*sd)**2 &
            -py*(-2*(tau/sigy)*sd - 2*sr*sd**2 + sr*(sigy_dd/sigy))

    end do

    dd = -2 * dd / n

    do concurrent ( i = 1:n )
       do concurrent( j = 1:i-1 )
          do concurrent ( l = 1:2 )
             dd(l,i,j) = dd(l,j,i)
          end do
       end do
    end do

  end function linejacc


  subroutine theline(tau,dtau,rho,drho,tol,reli)

    use minpacks

    real(REAL64), intent(in out) :: tau, rho
    real(REAL64), intent(out) :: dtau, drho
    real(REAL64), intent(in) :: tol
    logical, intent(out) :: reli

    integer :: maxfev = 666
    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(2) :: x, diag
    real(REAL64), dimension(2,2) :: cov
    character(len=666) :: msg
    character(len=5) :: label
!    real(REAL64) :: tol
    integer :: npar, info, n, i, K

    n = size(mju)
    K = size(kmean)
    npar = 2
    maxfev = 666*npar
!    tol = sqrt(epsilon(x))
!    maxfev = 200*npar
!    allocate(x(npar),diag(npar),cov(npar,npar))


!    x(1) = tau
!    x(2) = rho

    x = [ tau, rho ]
!    diag = [ 1/tau, 1/(rho + 1) ]
!    x(2) = rho
!    x(3) = sx
!    x(4) = sy
!    x(3:) = delta
!    diag(1:3) = [ 0.001, 1e2, 0.1 ]
!    diag(1:3) = diag(1:3)**2
!    diag(4:) = sqrt((pht+delta)*(1+1/tau) + (rho/tau)**2)
!    diag = 1 / diag**2
!    diag(2) = -diag(2)
!    diag = 1

!    backsx = sx
!    backsy = sy
!    backrho = rho

    if( analytic ) then
!       tol = 10*n*epsilon(x)*(1 + abs(maxval(x)))
       label = 'lmder'
       call lmder3(thelike,x,tol=tol,maxev=maxfev,factor=factor,nprint=nprint,info=info,msg=msg)
!       if( 0 < info .and. info < 5 ) then
!          block
!            real(REAL64), dimension(:,:), allocatable :: fjac
!            fjac = jacf(x(1),x(2),x(3:),sx,sy)
!            call qrinv(fjac,cov)
!          end block
!       end if
    else
!       tol = 1e-5
       label = 'lmdif'
       call lmdif3(theliked,x,tol=tol,maxev=maxfev,factor=factor,nprint=nprint,info=info,msg=msg)
!       call lmdif3(theliked,x)
       if( 0 < info .and. info < 5 ) then
          block

            real(REAL64), dimension(npar) :: fvec,z
            real(REAL64), dimension(npar,npar) :: djac
            real(REAL64), dimension(-1:1,npar) :: dg
            real(REAL64), parameter :: eps = epsilon(djac)**(1.0/3.0)
            integer :: i,j,l
!            integer :: iflag
!            iflag = 1
!            call fdjac2(liked,npar,npar,x,fvec,fjac,npar,iflag,epsilon(x),wa4)
!            fjac = jacf(x(1),x(2),x(3:),sx,sy)
!            call qrinv(fjac,cov)

!!$            do i = 1, npar
!!$               do j = i, npar
!!$                  do l = -1,1,2
!!$                     z = x
!!$                     z(i) = z(i) + l*eps
!!$!                     c = linedisp(z(1),z(2),z(3:),sx,sy)
!!$                     !dc = linegradc(z(1),z(2),z(3:),sx,sy)
!!$!                     dg(l,:) = gradlike(z(1),rho,z(2:),sx,sy) !&
!!$!                  + dc(1,:)*(nu*c(1)-lambda(1)) + dc(2,:)*(nu*c(2)-lambda(2))
!!$                  end do
!!$                  djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps)
!!$               end do
!!$            end do
!!$            call qrinv(djac,cov)
!            do i = 1, npar
!               write(*,*) i,real(djac(i,i)),real(cov(i,i))
!            end do
          end block
       end if
    end if
    diag = [ (cov(i,i),i=1,npar) ]
!    diag = [ (cov(i,i),i=1,2) ]
!    if( diag(3) < 0 ) diag(3) = 1
!    if( cov(3,3) < 0 ) then
!       cov(3,:) = 0
!       cov(:,3) = 0
!       cov(3,3) = 1
!    end if

!    goto 22
!    diag = [ (cov(i,i),i=1,npar) ]
!22 continue

    reli = tmin < x(1) .and. x(1) < tmax &
         .and. abs(x(2)) < rhomax &
!         .and. 0.1 < x(3) .and. x(3) < 10 &
!         .and. 0.1 < x(4) .and. x(4) < 10 &
!         .and. all(diag > 0) &
         .and. (0 < info .and. info < 5)

!    do i = 1, npar
!       write(*,*) i,diag(i)
!    end do

!!$    if( any(diag < 0) .and. info > 0 .and. debug_iter ) then
!!$       block
!!$         integer :: i
!!$         real(REAL64), dimension(npar,npar) :: fjac
!!$!         fjac = jacf1(x(1),siguback,sigvback,sx,sy,x(2:))
!!$
!!$         do i = 1, npar
!!$!            if( .not. (diag(i) > 0) .or. i == 1 ) &
!!$            write(*,*) i,real(diag(i)),real(cts(i)),real(pht(i)),&
!!$                 real(x(i1))!,real(fjac(i,i))
!!$         end do
!!$       end block
!!$    end if

    if( debug_iter ) then
       write(error_unit,'(3a,i0,3a,es8.1,a,l2)') ' theline() ',label,' info=', &
            info,': ',trim(msg),' tol:',tol,'  reliable:',reli
!       write(error_unit,*) real(x)
    end if

    if( reli ) then
       ! Sometimes, the lambda parameter of LM-method is incidentally tunneled
       ! into the mirrored part of the parameter space. Fortunately, the like
       ! function is pretty symmetric, so I do to ignore the "tunnelling".
       ! The correct, in principle, solution should to keep the condition:
       ! rho > 0.
       tau = x(1)
       rho = abs(x(2))
!       phi = x(3)
 !      sx = x(3)
!       sy = x(4)
!       delta = x(2:)

!!$       if( all(diag > 0) ) then
!!$          dtau =  sqrt(cov(1,1)/(2*n - npar))
!!$          drho = sqrt(cov(2,2)/(2*n - npar))
!!$!          dmu = sqrt(cov(3:,3:)/(2*n - npar))
!!$!          dphi = sqrt(cov(3,3)/(2*n - npar))
!!$          !       dsigv = sqrt(cov(3,3)/(2*n - npar))
!!$       end if
    end if

    if( debug_iter .and. (0 < info .and. info < 5) ) then
!    if( .false. ) then
       block
!         integer, parameter :: maxpar = 7
         real(REAL64), dimension(npar,npar) :: cor,cov,fjac
         real(REAL64), dimension(npar) :: fvec
         real(REAL64) :: q
         integer :: i,j,n,iflag
!         character(len=3) :: dots
         character(len=80) :: fmt

         iflag = 2
         call thelike(npar,npar,x,fvec,fjac,npar,iflag)
         fjac = matmul(transpose(fjac),fjac)
         call qrinv(fjac,cov)
!         cov = fjac


         cor = 1
         do i = 1, npar
            do j = i+1, npar
               if( cov(i,i) > 0 .and. cov(j,j) > 0 &
                    .and. cov(i,i) < 1e10 .and. cov(j,j) < 1e10) then
                  q = sqrt(cov(i,i)*cov(j,j))
!                  write(*,*) cov(i,i),cov(j,j),q,cov(i,j)
                  if( 1e-10 < q .and. q < 1e10 .and. abs(cov(i,j)) < 1e10 ) then
                     cor(i,j) = cov(i,j) / q
                  else
                     cor(i,j) = huge(cor)
                  end if
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

!!$         n = min(maxpar,npar)
!!$         if( n < npar ) then
!!$            dots = '...'
!!$         else
!!$            dots = ''
!!$         end if

         n = npar
         write(error_unit,*) 'Covariance (Hessian in minimum):         Correlation: '
         write(fmt,'(a,i0,a,i0,a)') '(i2,a,',n,'(es9.2,1x),17x,',n,'(f6.3,1x))'
!         write(*,*) ":",trim(fmt),":"
         do i = 1, n
            write(error_unit,fmt) i,': ',cov(i,1:n),cor(i,1:n)
         end do

!         write(error_unit,*) 'Correlation:'
!         write(fmt,'(a,i0,a)') '(i2,a,',n,'(f6.3,1x),a)'
!         do i = 1, n
!            write(error_unit,fmt) i,': ',cor(i,1:n),dots
!         end do

       end block
    end if


  end subroutine theline

  subroutine thelike(npar,mpar,x,fvec,fjac,ldfjac,iflag)

    use rfun

    integer, intent(in) :: npar, mpar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(in out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(in out) :: fjac

    character(len=80) :: msg
    real(REAL64), dimension(:), allocatable :: sigx, sigy, zx, zy, psix, psiy, &
         dpsix, dpsiy, sigx_t, sigx_r, sigy_t, sigy_r
    real(REAL64), dimension(2,2) :: s2jac, cjac, dcc
    real(REAL64), dimension(2) :: dc, s2grad
    real(REAL64) :: tau, rho, c, ds
    integer :: i, j, stat

    tau = x(1)
    rho = x(2)

    if( iflag == 0 ) then

       if( debug_iter ) then
          block
            real(REAL64) :: logl
            logl = loglike(tau,rho,mju)
            c = cvar(tau,rho) !/ ndata
            write(error_unit,*) 'thelike: loglike=',real(logl),' cvar=',real(c)
            write(error_unit,*) 'x = ',real(x)
            write(error_unit,*) 'fvec = ',real(fvec)
          end block
       end if

       ! sanity checks
       if( tau < tmin .or. tau > tmax ) then
          write(error_unit,*) 'phline(thelike): tau < tmin .or. tau > tmax'
          iflag = -1
       end if
       if( abs(rho) > rhomax ) then
          write(error_unit,*) 'phline(thelike): abs(rho) > rhomax'
          iflag = -1
       end if
       if( nu > 1e16 ) then
          write(error_unit,*) 'phline(thelike): nu > 1e16'
          iflag = -1
       end if

       return
    end if

    allocate(sigx(ndata),sigy(ndata),zx(ndata),zy(ndata),psix(ndata),psiy(ndata), &
         dpsix(ndata),dpsiy(ndata),sigx_t(ndata),sigx_r(ndata),sigy_t(ndata), &
         sigy_r(ndata),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       if( debug_iter ) write(error_unit,*) 'phline(thelike): ',trim(msg)
       iflag = -1
       return
    end if

    sigx = sigmax(tau,rho,pht)
    sigy = sigmay(tau,rho,cts)
    zx = (mju - pht) / sigx
    zy = (cts - tau*mju) / sigy
    if( robust ) then
       psix = huber(zx)
       psiy = huber(zy)
       dpsix = dhuber(zx)
       dpsiy = dhuber(zy)
!       call hubers(zx,psix,dpsix) update !
       !       call hubers(zy,psiy,dpsiy) update !
      call tukeys(zx,psix,dpsix)
      call tukeys(zy,psiy,dpsiy)
    else
       psix = zx
       psiy = zy
       dpsix = 1
       dpsiy = 1
    end if

    ds = sum(dpsix + dpsiy)
    if( .not. (ds > 0) ) then
       if( debug_iter ) write(error_unit,*) 'phline(thelike): assumption ds > 0'
       iflag = -1
       return
    end if

    !    c = sigvar(psix,psiy,dpsix,dpsiy,ds) - (ndata - 1)*beta*ds
    c = sum(psix**2 + psiy**2) - beta*ds
!    c = sum(psix**2 + psiy**2)/ds - beta
!    c = c / ndata

    do i = 1, ndata
       call grad_sigx(sigx(i),pht(i),tau,rho,sigx_t(i),sigx_r(i))
       call grad_sigy(sigy(i),cts(i),tau,rho,sigy_t(i),sigy_r(i))
    end do

    s2grad = grad_ressum(tau,rho,mju,zx,zy,psix,psiy,dpsix,dpsiy,sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r)
    dc =     grad_sigvar(mju,zx,zy,psix,psiy,dpsix,dpsiy,sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r)
!    dc = dc / ds

!    dc = dc / ndata

!    write(error_unit,*) s2grad/ndata, dc, c

    if( iflag == 1 ) then

       fvec =  s2grad / ndata &
!            + grad_sigma(sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r) / ndata &
            + (nu*c - lambda)*dc
!       fvec =  ds + grad_sigma(sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r,dpsix,dpsiy)

       nfev = nfev + 1
!       write(*,*) fvec

    else if( iflag == 2 ) then

       do i = 1, 2
          do j = 1, 2
             cjac(i,j) = dc(i)*dc(j)
          end do
       end do

       block
         real(REAL64), dimension(:), allocatable :: sigx_tt,sigx_tr,sigx_rr, &
              sigy_tt,sigy_tr,sigy_rr

         allocate(sigx_tt(ndata),sigx_tr(ndata),sigx_rr(ndata), &
              sigy_tt(ndata),sigy_tr(ndata),sigy_rr(ndata),stat=stat,errmsg=msg)
         if( stat /= 0 ) then
            if( debug_iter ) write(error_unit,*) 'phline(thelike): ',trim(msg)
            iflag = -1
            return
         end if

         do i = 1, ndata
            call jac_sigx(sigx(i),pht(i),tau,rho,sigx_t(i),sigx_r(i), &
                 sigx_tt(i),sigx_tr(i),sigx_rr(i))
            call jac_sigy(sigy(i),cts(i),tau,rho,sigy_t(i),sigy_r(i), &
                 sigy_tt(i),sigy_tr(i),sigy_rr(i))
         end do

         s2jac = jac_ressum(tau,rho,mju,sigx,sigy,zx,zy,psix,psiy,dpsix,dpsiy, &
              sigx_t,sigx_r,sigy_t,sigy_r, &
              sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr)
         dcc = jac_cvar(tau,rho,mju,sigx,sigy,zx,zy,psix,psiy,dpsix,dpsiy, &
              sigx_t,sigx_r,sigy_t,sigy_r, &
              sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr)
!         dcc = dcc / ndata

         fjac = s2jac / ndata &
!              + jac_sigma(sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r, &
!              sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr) / ndata &
              + nu*cjac + (nu*c - lambda)*dcc

       end block

    end if

  end subroutine thelike

  subroutine theliked(npar,mpar,x,fvec,iflag)

    integer, intent(in) :: npar, mpar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(in out) :: fvec

    real(REAL64), dimension(size(x),size(fvec)) :: fjac
    integer :: nflag

    if( iflag == 0 ) then
       nflag = 0
    else if( iflag == 1 ) then
       nflag = 1
    end if
    call thelike(npar,mpar,x,fvec,fjac,npar,nflag)
    if( nflag < 0 ) iflag = nflag
!    write(*,*) fvec

  end subroutine theliked



  function grad_ressum(tau,rho,mju,zx,zy,psix,psiy,dpsix,dpsiy,sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r) &
       result(vec)

    use rfun
    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: mju,sigx,sigy,zx,zy,psix,psiy, &
         dpsix,dpsiy,sigx_t,sigx_r,sigy_t,sigy_r
    real(REAL64), dimension(2) :: vec

    real(REAL64), dimension(size(zx)) :: px, py, sx, sx_t, sx_r, sy, sy_t, sy_r, psi
    integer :: i

!!$    where( abs(zx) > 2 )
!!$       sigx = sigmax(tau,rho,mju)
!!$       zx = (mju - pht) / sigx
!!$    end where
!!$    where( abs(zy) > 2 )
!!$       sigy = sigmay(tau,rho,tau*mju)
!!$       zy = (cts - tau*mju) / sigy
!!$    end where

    sx = sigx
    sx_t = sigx_t
    sx_r = sigx_r

    sy = sigy
    sy_t = sigy_t
    sy_r = sigy_r

    px = psix*zx - real(1,REAL64)
    py = psiy*zy - real(1,REAL64)
    psi = psiy

    if( .false. .and. robust ) then
       do i = 1, ndata
          if( abs(zx(i)) > 3 ) then
             sx(i) = sigmax(tau,rho,mju(i))
             call grad_sigx(sx(i),mju(i),tau,rho,sx_t(i),sx_r(i))
             px(i) = (pht(i) - mju(i)) / sx(i)
             px(i) = px(i)*huber(px(i)) - 1
             px(i) = 0
          end if
          if( abs(zy(i)) > 3 ) then
             sy(i) = sigmay(tau,rho,tau*mju(i))
             call grad_sigy(sy(i),tau*mju(i),tau,rho,sy_t(i),sy_r(i))
             py(i) = (cts(i) - tau*mju(i)) / sy(i)
             psi(i) = huber(py(i))
             py(i) = py(i)*psi(i) - 1
!             psi(i) = huber(100*py(i))
             py(i) = 0
          end if
       end do
    end if

!    px = psix**2/2
!    py = psiy**2/2
    !    mu = huber(mju/sigy)

    vec(1) = -sum(px*(sx_t/sx) + py*(sy_t/sy) + psi*(mju/sy))
    vec(2) = -sum(px*(sx_r/sx) + py*(sy_r/sy))
    return

    vec(1) = -sum(px*(sigx_t/sigx) + py*(sigy_t/sigy) + psiy*(mju/sigy) - &
         sigx_t/sigx - sigy_t/sigy)
!    vec(1) = -sum(px*(sigx_t/sigx) + py*(sigy_t/sigy) + psiy*mu)
    vec(2) = -sum(px*(sigx_r/sigx) + py*(sigy_r/sigy) - &
         sigx_r/sigx - sigy_r/sigy)

    write(*,*) 'ressum 1:',real(sum(px*(sigx_t/sigx))), &
         real(sum(py*(sigy_t/sigy))),real(sum(psiy*(mju/sigy)))
    write(*,*) 'ressum 2:',real(sum(px*(sigx_r/sigx))), &
         real(sum(py*(sigy_r/sigy)))

  end function grad_ressum

  function jac_ressum(tau,rho,mju,sigx,sigy,zx,zy,psix,psiy,dpsix,dpsiy, &
       sigx_t,sigx_r,sigy_t,sigy_r,sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr) &
       result(jac)

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: mju,sigx,sigy,zx,zy,psix,psiy, &
         dpsix,dpsiy,sigx_t,sigx_r,sigy_t,sigy_r,&
         sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr
    real(REAL64), dimension(2,2) :: jac

    real(REAL64) :: &
!         sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr,&
         chi2_tt, chi2_tr, chi2_td, chi2_rr, chi2_rd, chi2_dd, &
         px, py, zx2, zy2, sxt2, syt2, sxtr, sytr, sxr2, syr2, dzy_t, dpx, dpy
    integer :: i,j

    jac = 0
    do i = 1, ndata

!       call jac_sigx(sigx(i),pht(i),tau,rho,sigx_t(i),sigx_r(i), &
!            sigx_tt,sigx_tr,sigx_rr)
!       call jac_sigy(sigy(i),cts(i),tau,rho,sigy_t(i),sigy_r(i), &
!            sigy_tt,sigy_tr,sigy_rr)

       px = psix(i)*zx(i)
       py = psiy(i)*zy(i)

       zx2 = zx(i)**2
       zy2 = zy(i)**2
       sxt2 = (sigx_t(i)/sigx(i))**2
       syt2 = (sigy_t(i)/sigy(i))**2
       sxtr = (sigx_t(i)/sigx(i))*(sigx_r(i)/sigx(i))
       sytr = (sigy_t(i)/sigy(i))*(sigy_r(i)/sigy(i))
       sxr2 = (sigx_r(i)/sigx(i))**2
       syr2 = (sigy_r(i)/sigy(i))**2
       dzy_t = mju(i)/sigy(i) + zy(i)*(sigy_t(i)/sigy(i))
       dpx = dpsix(i)*zx(i) + psix(i)
       dpy = dpsiy(i)*zy(i) + psiy(i)
       chi2_tt = -dpsix(i)*zx2*sxt2 - 2*px*sxt2 + px*(sigx_tt(i)/sigx(i)) &
            -dpsiy(i)*dzy_t*zy(i)*(sigy_t(i)/sigy(i)) &
            - psiy(i)*(mju(i)/sigy(i) + 2*zy(i)*(sigy_t(i)/sigy(i)))*(sigy_t(i)/sigy(i)) &
            + py*(sigy_tt(i)/sigy(i)) &
            -(dpsiy(i)*dzy_t + psiy(i)*(sigy_t(i)/sigy(i)))*(mju(i)/sigy(i))
!            -dpsiy*dzy_t*(mu/sigy) &
!            -psiy*(mu/sigy)*(sigy_t/sigy)
!               + dpsiy*zy2*syt2 + 2*py*syt2 - py*(sigy_tt/sigy) &
!            - dpsiy*(mu/sigy)**2 - (dpsiy*zy + psiy)*(sigy_t/sigy)*(mu/sigy)
!               + (dpsiy*(1-zy*(sigy_t/sigy))-psiy)*(mu/sigy)**2
!       chi2_tt = dpsix*zx2*sxt2 + dpsiy*zy2*syt2 + dpsiy*(mu/sigy)**2
!       chi2_tr = -dpsix*zx2*sxtr + 2*px*sxtr - px*(sigx_tr/sigx) &
       chi2_tr =-dpx*zx(i)*(sigx_t(i)/sigx(i))*(sigx_r(i)/sigx(i)) &
               - dpy*dzy_t*(sigy_r(i)/sigy(i)) &
               - px*(sigx_t(i)/sigx(i))*(sigx_r(i)/sigx(i)) + px*(sigx_tr(i)/sigx(i)) &
               - py*(sigy_t(i)/sigy(i))*(sigy_r(i)/sigy(i)) + py*(sigy_tr(i)/sigy(i))

!               + dpsiy*zy2*sytr + 2*py*sytr - py*(sigy_tr/sigy) &
!               - psiy*(mu/sigy)*(sigy_r/sigy)
!       chi2_td =-px*(sigx_t/sigx) + py*tau*(sigy_t/sigy) &
!            + (dpsiy - psiy)*(mu/sigy)*(tau/sigy)
!       chi2_td =-dpx*(sigx_t(i)/sigx(i))/sigx(i) &
!               + dpsiy(i)*dzy_t*tau/sigy(i) + psiy(i)*(sigy_t(i)/sigy(i))*tau/sigy(i) - psiy(i)/sigy(i)
       chi2_rr = dpsix(i)*zx2*sxr2 + (2*sxr2 - sigx_rr(i)/sigx(i))*px &
               + dpsiy(i)*zy2*syr2 + (2*syr2 - sigy_rr(i)/sigy(i))*py
!       chi2_rd =-dpx*(sigx_r(i)/sigx(i))/sigx(i) + dpy*(sigy_r(i)/sigy(i))*tau/sigy(i)
!       chi2_dd = dpsix(i)/sigx(i)**2 + dpsiy(i)*(tau/sigy(i))**2

       jac(1,1:2) = jac(1,1:2) - [ chi2_tt, chi2_tr]
       jac(2,2) = jac(2,2) + chi2_rr

    end do

    jac(2,1) = jac(1,2)

!    jac = jac / (2*ndata)

!    do i = 1, size(jac,1)
!       do j = 1, i - 1
!          jac(i,j) = jac(j,i)
!       end do
!    end do

  end function jac_ressum

  real(REAL64) function sigvar(psix,psiy,dpsix,dpsiy,ds)

    real(REAL64), dimension(:), intent(in) :: psix,psiy,dpsix,dpsiy
    real(REAL64), intent(in) :: ds
!    real(REAL64) :: s2
!    integer :: n

!    ds = sum(dpsix + dpsiy)

!    if( ds > 0 ) then
!    n = ndata
!    s2 = sum(psix**2 + psiy**2)
!    sigvar = (s2 / ds) * (real(n) / real(n-1))
!    sigvar = s2 !- (n - 1)*scale*ds
    sigvar = sum(psix**2 + psiy**2)
!    else
!       stop '(ds > 0) assumption failed.'
!    end if

  end function sigvar

  function grad_sigvar(mju,zx,zy,psix,psiy,dpsix,dpsiy,sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r) result(vec)

!    real(REAL64), intent(in) :: ds
    real(REAL64), dimension(:), intent(in) :: mju,zx,zy,psix,psiy,dpsix,dpsiy,&
         sigx,sigy,sigx_t,sigy_t,sigx_r,sigy_r
    real(REAL64), dimension(2) :: vec
    real(REAL64), dimension(size(zx)) :: px, py
!    real(REAL64) :: ds

    px = psix*dpsix*zx
    py = psiy*dpsiy*zy
!    px = psix*zx
!    py = psiy*zy
!    integer :: n

!    n = ndata
!    ds = sum(dpsix + dpsiy)
!    if( ds > 0 ) then
       !       vec = (2*n*ds2)/ds * (n/ds) * (real(n) / real(n-1))
!    vec = 2*(s2grad / ds) * (real(n) / real(n-1))

    vec(1) = -2*sum(px*(sigx_t/sigx) + py*(sigy_t/sigy) + psiy*dpsiy*(mju/sigy))
!    vec(1) = -2*sum(px*(sigx_t/sigx) + py*(sigy_t/sigy) + psiy*(mju/sigy))
    vec(2) = -2*sum(px*(sigx_r/sigx) + py*(sigy_r/sigy))

!    write(*,*) 'sigvar 1:',real(sum(px*(sigx_t/sigx))),&
!         real(sum(py*(sigy_t/sigy))),real(sum(psiy*dpsiy*(mju/sigy)))
!    write(*,*) 'sigvar 2:',real(sum(px*(sigx_r/sigx))), &
!         real(sum(py*(sigy_r/sigy)))

!    vec = vec / ds

!    else
!       stop '(ds > 0) assumption failed.'
!    end if

  end function grad_sigvar

  function jac_cvar(tau,rho,mju,sigx,sigy,zx,zy,psix,psiy,dpsix,dpsiy, &
       sigx_t,sigx_r,sigy_t,sigy_r,sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr) &
       result(jac)

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: mju,sigx,sigy,zx,zy,psix,psiy, &
         dpsix,dpsiy,sigx_t,sigx_r,sigy_t,sigy_r, &
         sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr
    real(REAL64), dimension(size(zx)) :: sx_t, sx_r, sy_t, sy_r, &
         sx_tt, sx_tr, sx_rr, sy_tt, sy_tr, sy_rr
    real(REAL64), dimension(2,2) :: jac, jacx, jacy
    real(REAL64) :: ds

!    real(REAL64) :: &
!         sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr,&
!         chi2_tt, chi2_tr, chi2_td, chi2_rr, chi2_rd, chi2_dd, &
!         px, py, zx2, zy2, sxt2, syt2, sxtr, sytr, sxr2, syr2, dzy_t, dpx, dpy
!    integer :: i,j

    sx_t = -zx*(sigx_t/sigx)
    sx_r = -zx*(sigx_r/sigx)
    sx_tt = (2*(sigx_t/sigx)**2 - sigx_tt/sigx)*zx
    sx_tr = (2*(sigx_t/sigx)*(sigx_r/sigx) - sigx_tr/sigx)*zx
    sx_rr = (2*(sigx_r/sigx)**2 - sigx_rr/sigx)*zx

    sy_t = -zy*(sigy_t/sigy) - mju/sigy
    sy_r = -zy*(sigy_r/sigy)
    sy_tt = (2*(sigy_t/sigy)**2 - sigy_tt/sigy)*zy + 2*(mju/sigy)*(sigy_t/sigy)
    sy_tr = (2*(sigy_t/sigy)*(sigy_r/sigy) - sigy_tr/sigy)*zy + &
         (mju/sigy)*(sigy_r/sigy)
    sy_rr = (2*(sigy_r/sigy)**2 - sigy_rr/sigy)*zy

    jacx(1,1) = sum((dpsix*sx_t)**2 + psix*dpsix*sx_tt)
    jacy(1,1) = sum((dpsiy*sy_t)**2 + psiy*dpsiy*sy_tt)
    jacx(1,2) = sum(dpsix*sx_t*sx_r + psix*dpsix*sx_tr)
    jacy(1,2) = sum(dpsiy*sy_t*sy_r + psiy*dpsiy*sy_tr)
    jacx(2,2) = sum((dpsix*sx_r)**2 + psix*dpsix*sx_rr)
    jacy(2,2) = sum((dpsiy*sy_r)**2 + psiy*dpsiy*sy_rr)

    jacx(2,1) = jacx(1,2)
    jacy(2,1) = jacy(1,2)

    ds = sum(dpsix + dpsiy)
    jac = 2*(jacx + jacy) / ds



  end function jac_cvar

  function grad_sigma(sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r) result(svec)

    real(REAL64), dimension(:), intent(in) :: sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r
    real(REAL64), dimension(2) :: svec

    svec(1) = sum(sigx_t/sigx + sigy_t/sigy)
    svec(2) = sum(sigx_r/sigx + sigy_r/sigy)

!    svec = svec / (1*ndata)
    !svec = svec / sum(dpsix + dpsiy)

  end function grad_sigma

  function jac_sigma(sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r, &
       sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr) result(jac)

    real(REAL64), dimension(:), intent(in) :: sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r,&
         sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr
    real(REAL64), dimension(2,2) :: jac

    jac(1,1) = sum(sigx_tt/sigx - (sigx_t/sigx)**2 + sigy_tt/sigy - (sigy_t/sigy)**2)
    jac(1,2) = sum(sigx_tr/sigx - (sigx_t/sigx)*(sigx_r/sigx) &
         + sigy_tr/sigy - (sigy_t/sigy)*(sigy_r/sigy))
    jac(2,2) = sum(sigx_rr/sigx - (sigx_r/sigx)**2 + sigy_rr/sigy - (sigy_r/sigy)**2)

    jac(2,1) = jac(1,2)

  end function jac_sigma

  real(REAL64) elemental pure function sigmax(tau,rho,x)

    real(REAL64), intent(in) :: tau,rho,x

    sigmax = sqrt(x*(1+1/tau)+(rho/tau)**2) / 2

  end function sigmax

  real(REAL64) elemental pure function sigmay(tau,rho,y)

    real(REAL64), intent(in) :: tau,rho, y

    sigmay = sqrt(y*(1+tau)+rho**2) / 2

  end function sigmay

  real(REAL64) function cvar(tau,rho)

    use rfun

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(:), allocatable :: sigx, sigy, zx, zy, psix, psiy, &
         dpsix, dpsiy
    !    real(REAL64) :: ds
    integer :: n

    n = ndata
    allocate(sigx(n),sigy(n),zx(n),zy(n),psix(n),psiy(n),dpsix(n),dpsiy(n))

    sigx = sigmax(tau,rho,pht)
    sigy = sigmay(tau,rho,cts)
    zx = (mju - pht) / sigx
    zy = (cts - tau*mju) / sigy
    if( robust ) then
       psix = huber(zx)
       psiy = huber(zy)
       dpsix = dhuber(zx)
       dpsiy = dhuber(zy)
!      call tukeys(zx,psix,dpsix)
!      call tukeys(zy,psiy,dpsiy)
    else
       psix = zx
       psiy = zy
       dpsix = 1
       dpsiy = 1
    end if
!    ds = sum(dpsix + dpsiy)

!    logvar = sigvar(psix,psiy,dpsix,dpsiy,ds)
!    cvar = sum(psix**2 + psiy**2) - (n - 1)*beta*sum(dpsix + dpsiy)
    cvar = sum(psix**2 + psiy**2) - beta*(2*ndata) !sum(dpsix + dpsiy)
!    cvar = sum(psix**2 + psiy**2)/sum(dpsix + dpsiy) - beta
!    write(*,*) sum(psix**2 + psiy**2),beta,sum(dpsix + dpsiy)

  end function cvar

  subroutine funky(tau,rho,lambda,f,fvec,fjac,c,cvec,cjac)

    use rfun

    real(REAL64), intent(in) :: tau,rho,lambda
    real(REAL64), intent(out) :: f,c
    real(REAL64), dimension(:), intent(out) :: fvec,cvec
    real(REAL64), dimension(:,:), intent(out) :: fjac,cjac

    real(REAL64), dimension(:), allocatable :: sigx, sigy, zx, zy, psix, psiy, &
         dpsix, dpsiy, sigx_t, sigx_r, sigy_t, sigy_r
    real(REAL64) :: ds
    integer :: i

    allocate(sigx(ndata),sigy(ndata),zx(ndata),zy(ndata),psix(ndata),psiy(ndata), &
         dpsix(ndata),dpsiy(ndata),sigx_t(ndata),sigx_r(ndata),sigy_t(ndata), &
         sigy_r(ndata))

    sigx = sigmax(tau,rho,pht)
    sigy = sigmay(tau,rho,cts)
    zx = (mju - pht) / sigx
    zy = (cts - tau*mju) / sigy
    if( robust ) then
       psix = huber(zx)
       psiy = huber(zy)
       dpsix = dhuber(zx)
       dpsiy = dhuber(zy)
!       psix = tanh(zx)
!       psiy = tanh(zy)
!       dpsix = 1-tanh(zx)**2
!       dpsiy = 1-tanh(zy)**2
!       call hubers(zx,psix,dpsix) update !
!       call hubers(zy,psiy,dpsiy) update !
!      call tukeys(zx,psix,dpsix)
!      call tukeys(zy,psiy,dpsiy)
    else
       psix = zx
       psiy = zy
       dpsix = 1
       dpsiy = 1
    end if

    ds = sum(dpsix + dpsiy)

    f = loglike(tau,rho,mju)
!    c = sum(psix**2 + psiy**2) - ds
!    c = sum(psix**2 + psiy**2)/ds - beta
    c = sum(psix**2 + psiy**2) - ds
!    c = sum(psix**2 + psiy**2) - beta*(2*ndata)

    do i = 1, ndata
       call grad_sigx(sigx(i),pht(i),tau,rho,sigx_t(i),sigx_r(i))
       call grad_sigy(sigy(i),cts(i),tau,rho,sigy_t(i),sigy_r(i))
    end do

    fvec = grad_ressum(tau,rho,mju,zx,zy,psix,psiy,dpsix,dpsiy,sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r)
    cvec = grad_sigvar(mju,zx,zy,psix,psiy,dpsix,dpsiy,sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r)
!    cvec = cvec / ds

    fjac = 0
    cjac = 0

  end subroutine funky



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine like(npar,mpar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, mpar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

!    real(REAL64), dimension(4) :: vec
!    real(REAL64), dimension(4,4) :: jac
    !    real(REAL64) :: tau, sigu, sigv
!    real(REAL64), dimension(2) :: c
    real(REAL64), dimension(3) :: c
    real(REAL64), dimension(:,:), allocatable :: dc
    real(REAL64), dimension(:,:,:), allocatable :: dd
!    real(REAL64), dimension(:), allocatable :: su_d, sv_d
!    real(REAL64) :: su_t, su_u, sv_v
    character(len=3) :: dots
    integer :: i,j,k,m

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          m = min(size(x),8)
          if( size(x) > 8 ) dots = '...'
          write(error_unit,*) 's0 = ',real(sqrt(sum(x(2:)**2)/(size(x)-1)))!, &
!               ' loglike=',real(loglike(x,backrho,backsx,backsy))
          write(error_unit,*) 'x = ',real(x(1:m)),trim(dots)
          write(error_unit,*) 'fvec = ',real(fvec(1:m)),trim(dots)
          block
            integer :: i,n
            character(len=80) :: fmt
            n = min(6,npar)
            write(fmt,'(a,i0,a)') '(a,i0,a,',n,'(1x,es9.2),a)'
            do i = 1, -n
               write(error_unit,fmt) &
                    'fjac(',i,',:)=',fjac(i,1:n),trim(dots)
            end do
          end block
       end if

       ! convergence check
       if( x(1) < tmin .or. x(1) > tmax ) iflag = -1
!       if( x(1) < tmin .or. x(1) > tmax .or. abs(x(2)) > rhomax ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

!       if( any(x < 0) ) then
!          fvec = 1e10
!          return
!       end if


       c(2:3) = linedisp(x(1),backrho,x(2:),backsx,backsy)
       dc(2:3,:) = linegradc(x(1),backrho,x(2:),backsx,backsy)

!       fvec = gradf(x(1),x(2),x(3:),backsx,backsy) &
       !            + [ (dc(k,:)*(nu*c(k) - lambda(k)), k=1,2) ]
       fvec = gradf(x(1),backrho,x(2:),backsx,backsy)
!       fvec = fvec + dc(1,:)*(nu*c(1) - lambda(1)) + dc(2,:)*(nu*c(2) - lambda(2))



    else if( iflag == 2 ) then

       return
!       sv = 1
!       su = 1
!       fjac = jacf(tau,sigu,sigv,s)
       c(2:3) = linedisp(x(1),backrho,x(2:),backsx,backsy)
       dc(2:3,:) = linegradc(x(1),backrho,x(2:),backsx,backsy)
       dd = linejacc(x(1),backrho,x(2:),backsx,backsy)
       fjac = jacf(x(1),backrho,x(2:),backsx,backsy)

       do i = 1, size(fjac,1)
          do j = i, size(fjac,2)
             fjac(i,j) = fjac(i,j) - sum(lambda*dd(:,i,j)) &
                  +nu*(sum(dc(:,i)*dc(:,j)) + sum(c(:)*dd(:,i,j)))
          end do
       end do

       do concurrent ( i = 1:size(fjac,1) )
          do concurrent ( j = i+1:size(fjac,2) )
             fjac(j,i) = fjac(i,j)
          end do
       end do

!       fjac(1,1) = fjac(1,1) - lambda(1)*su_tt + nu*(su_t**2 - c(1)*su_tt)
!       do i = 1, n
!          fjac(i,i) = fjac(i,i) - lambda(1)*su_dd(i) -lambda(2)*sv_dd(i) &
!          + nu*(su_d(i)**2 - c(1)*su_dd(i)) + sv_d(i)**2 - c(2)*sv_dd(i))
!       end do
!       fjac(1,2:) = fjac(1,2:) &
!            - lambda(1)*su_td + nu*(su_t*su_d - c(1)*su_td)

    end if

  end subroutine like

  subroutine liked(npar,mpar,x,fvec,iflag)

    integer, intent(in) :: npar, mpar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec
    real(REAL64), dimension(3) :: c
    real(REAL64), dimension(:,:), allocatable :: dc
!    real(REAL64), dimension(:), allocatable :: sc
    character(len=3) :: dots
    integer :: i, k, m

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          m = min(size(x),8)
          if( size(x) > 8 ) dots = '...'
          write(error_unit,*) 's0 = ',real(sqrt(sum(x(2:)**2)/(size(x)-1)))!, &
!               ' loglike=',real(loglike(x,backrho,backsx,backsy))
          write(error_unit,*) 'x = ',real(x(1:m)),trim(dots)
          write(error_unit,*) 'fvec = ',real(fvec(1:m)),trim(dots)
       end if

!       if( x(1) < tmin .or. x(1) > tmax .or. abs(x(2)) > rhomax &
!            .or. d() iflag = -1

       ! convergence check ?
!       if( x(1) < tmin .or. x(1) > tmin .or. x(2) < 0 .or. &
!            x(3) < 0.1 .or. x(3) > 10 .or. x(4) < 0.1 .or. x(4) > 10 ) then
!          fvec = x(1)*maxval(cts)
!       end if
       return
    end if

    allocate(dc(3,npar))

    fvec = gradlike(x(1),backrho,x(2:),backsx,backsy)
    c(1) = cscale(x(1),backrho,x(2:),backsx,backsy)
    dc(1,:) = cscale_grad(x(1),backrho,x(2:),backsx,backsy)
!    c(2:3) = linedisp(x(1),backrho,x(2:),backsx,backsy)
!    dc(2:3,:) = linegradc(x(1),backrho,x(2:),backsx,backsy)
    do i = 1, npar
!       fvec(i) = fvec(i) + dc(1,i)*(nu*c(1) - lambda(1))    !+ sum(dc(:,i)*(nu*c - lambda))
    end do
!    return
!    fvec = gradf(x(1),x(2),x(3:),backsx,backsy)
!    c = linedisp(x(1),x(2),x(3:),backsx,backsy)
!    dc = linegradc(x(1),x(2),x(3:),backsx,backsy)
!    fvec = fvec + dc(1,:)*(nu*c(1) - lambda(1)) + dc(2,:)*(nu*c(2) - lambda(2))

!    fvec = gradf(x(1),x(2),x(3),x(4:),backsx,backsy) &
!         + [ (dc(k,:)*(nu*c(k) - lambda(k)), k=1,2) ]


!    fvec(3) = fvec(3) - lambda(1) + nu*c(1)
!    fvec(4) = fvec(4) - lambda(2) + nu*c(2)

  end subroutine liked

  function gradlike(tau,rho,delta,sx,sy) result(vec)

    real(REAL64), intent(in) :: tau,rho,sx,sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(size(delta)+2) :: vec
    real(REAL64), dimension(3) :: gradx, grady

    real(REAL64) :: mu,d,r,sigx,sigy,rx,ry,rho2,tt, sig2
    integer :: i,n

    n = size(delta)

    tt = tau*(tau + 1)
    rho2 = rho**2

    vec = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - mu*tau

       sig2 = mu*tt + rho2

       if( sig2 > 0 ) then
!       if( .true. ) then

          sigy = sqrt(sig2) / 2
          sigx = sigy / tau
          rx = (d / sigx)/sx
          ry = (r / sigy)/sy

          if( robust ) then
             if( rx < rmin ) then
                gradx = grad_xlike(rmin,tau,rho,d,mu,sx,sigx)
             else if( rx > rmax ) then
                gradx = grad_xlike(rmax,tau,rho,d,mu,sx,sigx)
             else
                gradx = grad_xlike(rx,tau,rho,d,mu,sx,sigx)
             endif

             if( ry < rmin ) then
                grady = grad_ylike(rmin,tau,rho,r,mu,sy,sigy)
             else if( ry > rmax ) then
                grady = grad_ylike(rmax,tau,rho,r,mu,sy,sigy)
             else
                grady = grad_ylike(ry,tau,rho,r,mu,sy,sigy)
             endif

          else

             !rx = max(min(rx,rmax),rmin)
             gradx = grad_xlike(rx,tau,rho,d,mu,sx,sigx)

             !          ry = max(min(ry,rmax),rmin)
             grady = grad_ylike(ry,tau,rho,r,mu,sy,sigy)

          end if

!       if( ry < rmin ) then
!          y = ylike(rmin)
!       else if( ry > rmax ) then
!          y = ylike(rmax)
!       else
!          y = ylike(ry)
!       endif

!       rx = max(min((d / sigx)/sx,rmax),rmin)
!       ry = max(min((r / sigy)/sy,rmax),rmin)

       else ! sig2 <= 0
!          sig2 = rho2
!          gradx = grad_xlike(rmin,tau,rho,d,mu,sx,sigx)
          !          grady = grad_ylike(rmin,tau,rho,r,mu,sy,sigy)
          gradx = 0
          grady = 0
       end if

       vec(1:2) = vec(1:2) + (gradx(1:2) + grady(1:2))
       vec(i+2) = gradx(3) + grady(3)
       !vec(1) = vec(1) + (gradx(1) + grady(1))
       !vec(i+1) = gradx(3) + grady(3)

    end do
    vec = vec / n

  end function gradlike

  function grad_xlike(rx,tau,rho,d,mu,sx,sigx) result(grad)

    real(REAL64), intent(in) :: tau,rho,d,mu,rx,sx,sigx
    real(REAL64), dimension(3) :: grad
    real(REAL64) :: sigx_t, sigx_r, sigx_d, rs, q

!    call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)

    rs = rx / sx
    q = 1 - rs * (d / sigx)
    grad(1) = q*(sigx_t/sigx)
    grad(2) = q*(sigx_r/sigx)
    grad(3) = q*(sigx_d/sigx) + rs/sigx

  end function grad_xlike

  function grad_ylike(ry,tau,rho,r,mu,sy,sigy) result(grad)

    real(REAL64), intent(in) :: tau,rho,r,mu,ry,sy,sigy
    real(REAL64), dimension(3) :: grad
    real(REAL64) :: sigy_t, sigy_r, sigy_d, rs, q

!    call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)

    rs = ry / sy
    q = 1 - rs * (r/sigy)
    grad(1) = q*(sigy_t/sigy) - rs*(mu/sigy)
    grad(2) = q*(sigy_r/sigy)
    grad(3) = q*(sigy_d/sigy) - rs*tau/sigy

  end function grad_ylike

  function gradf(tau,rho,delta,sx,sy) result(vec)

    use rfun

    real(REAL64), intent(in) :: tau,rho,sx,sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(2+size(delta)) :: vec
    real(REAL64) :: mu,d,r,sigx,sigy,psix,psiy,rx,ry,r_t,r_d, rho2,&
         sigx_t,sigy_t,sigy_r,sigx_r,sigy_d, sigx_d,dt,dr,dd
    integer :: i, n

    n = size(pht)
    rho2 = rho**2

    vec = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - mu*tau
       r_t = -mu
       r_d = -tau
!       sigx = sqrt(mu*(1 + 1/tau) + (rho/tau)**2)/2
!       sigy = sqrt(mu*tau*(1+tau) + rho**2)/2
!       sigx_t = -(mu/sigx + (2*rho/tau)*(rho/sigx))/(2*tau**2)/4
!       sigy_t = (2*tau + 1)/2 * (mu/sigy)/4
!       sigx_r = (rho/sigx)/tau**2/4
!       sigy_r = rho / sigy/4
!       sigx_d = (1 + 1/tau)/2 / sigx/2
!       sigy_d = tau*(1 + tau)/2/sigy/2

!       sigx = sqrt(mu + sigv2) / 2
!       sigy = sqrt(mu*tau + sigu2)
       sigy = sqrt(mu*tau*(tau+1) + rho2) / 2
       sigx = sigy / tau

!       sigy_t = (mu / sigy) / 2
!       sigx_t = -c2*(rho/sigx) * (rho/tau**3)
!       sigy_r = (rho/sigy)*s2
!       sigx_r = (rho/sigx)*(cosfi/tau)**2
!       sigy_f = cs * (rho/sigy) * rho
!       sigx_f = -cs * (rho/sigx) * (rho/tau**2)
!       sigx_d = 1/(2*sigx)
!       sigy_d = tau / sigy / 2

!       call grad_sigx(sigx,mu,tau,rho,cosfi,sinfi,sigx_t,sigx_r,sigx_f,sigx_d)
!       call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!       sigx_f = 0
!       call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)
!       sigy_f = 0

       rx = (d / sigx)/sx
       ry = (r / sigy)/sy
       if( robust ) then
          psix = tukey(rx)
          psiy = tukey(ry)
       else
          psix = rx
          psiy = ry
       end if

       dt = -(psiy/sy)*(mu/sigy + (r/sigy)*(sigy_t/sigy)) + sigy_t/sigy &
            -(psix/sx)*(d/sigx)*(sigx_t/sigx) + sigx_t /sigx
       dr = -(psix/sx)*(d/sigx)*(sigx_r/sigx) + sigx_r/sigx &
            -(psiy/sy)*(r/sigy)*(sigy_r/sigy) + sigy_r/sigy
!       df = -(psix/sx)*(d/sigx)*(sigx_f/sigx) + sigx_f/sigx &
!            -(psiy/sy)*(r/sigy)*(sigy_f/sigy) + sigy_f/sigy
       dd =  (psix/sx)*(1/sigx - (d/sigx)*(sigx_d/sigx)) &
            -(psiy/sy)*(tau/sigy + (r/sigy)*(sigy_d/sigy)) &
            + sigx_d/sigx + sigy_d/sigy

       vec(1:2) = vec(1:2) + [ dt, dr ]
       vec(2+i) = dd


!!$       dt = -(psix/sx)*(d/sigx)*(sigx_t/sigx) + sigx_t/sigx &
!!$            +(psiy/sy)*(r_t/sigy - (r/sigy)*(sigy_t/sigy)) + sigy_t/sigy
!!$       dr = -(psix/sx)*(d/sigx)*(sigx_r/sigx) + sigx_r/sigx &
!!$            -(psiy/sy)*(r/sigy)*(sigy_r/sigy) + sigy_r/sigy
!!$       dd = (psix/sx)*(1/sigx - (d/sigx)*(sigx_d/sigx)) + sigx_d/sigx &
!!$           +(psiy/sy)*(r_d/sigy - (r/sigy)*(sigy_d/sigy)) + sigy_d/sigy
!!$       vec(1:2) = vec(1:2) + [ dt, dr ]
!!$       vec(2+i) = dd

    end do
    vec = vec / n

  end function gradf

  function jacf(tau,rho,delta,sx,sy) result(jac)

    use rfun

    real(REAL64), intent(in) :: tau,rho,sx,sy
    real(REAL64), dimension(:), intent(in) :: delta

    real(REAL64), parameter :: z = real(0.0,REAL64)
    real(REAL64), dimension(1+size(delta),1+size(delta)) :: jac
    real(REAL64) :: d,r,mu,sigx,sigy,psix,psiy,dpsix,dpsiy,rx,ry,r_t,r_d, &
         sigx_t,sigy_t,sigy_r,sigx_r,sigx_d,sigy_d, &
         sigx_tt,sigy_tt, sigx_tr,sigy_tr, sigx_td,sigy_td, &
         sigx_rr,sigy_rr, sigx_rd,sigy_rd, &
         sigy_fd, sigx_dd,sigy_dd, &
         dtt, dtr, dtf, dtd, drr, drf, drd, dff, dfd, ddd
!    real(REAL64) :: d,r,mu,sigx,sigy,psix,psiy,dpsix,dpsiy,rx,ry,r_t,r_d, &
!         sigx_t,sigy_t,sigy_r,sigx_r,sigx_f,sigy_f,sigx_d,sigy_d, &
!         sigx_tt,sigy_tt, sigx_tr,sigy_tr, sigx_tf,sigy_tf, sigx_td,sigy_td, &
!         sigx_rr,sigy_rr, sigx_rf,sigy_rf, sigx_rd,sigy_rd, &
!         sigx_ff,sigy_ff, sigx_fd,sigy_fd, sigx_dd,sigy_dd, &
!         dtt, dtr, dtf, dtd, drr, drf, drd, dff, dfd, ddd, &
!         sinfi, cosfi, c2,s2,cs,sigu, sigv, sigu2, sigv2
    integer :: i, n

    n = size(delta)
!    sinfi = sin(phi)
!    cosfi = cos(phi)
!    c2 = cosfi**2
!    s2 = sinfi**2
!    cs = sinfi*cosfi
!    sigu = rho*sinfi
!    sigv = rho*(cosfi/tau)
!    sigu2 = sigu**2
!    sigv2 = sigv**2

    jac = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - mu*tau
!       sigx = sqrt(mu + sigv2)
!       sigy = sqrt(mu*tau + sigu2)
       sigy = sqrt(mu*tau*(tau+1) + rho**2) / 2
       sigx = sigy / tau

       r_t = -mu
       r_d = -tau

!!$       sigy_t = (mu / sigy) / 2
!!$       sigx_t = -c2*(rho/sigx) * (rho/tau**3)
!!$       sigy_r = (rho/sigy)*s2
!!$       sigx_r = (rho/sigx)*(cosfi/tau)**2
!!$       sigy_f = cs * (rho/sigy) * rho
!!$       sigx_f = -cs * (rho/sigx) * (rho/tau**2)
!!$       sigx_d = 1/(2*sigx)
!!$       sigy_d = tau / sigy / 2
!!$
!!$       sigx_tt = c2/tau**3*(3/tau*(rho/sigx) + (rho/sigx)*(sigx_t/sigx))*rho
!!$       sigx_tr = c2/tau**3*((rho/sigx)**2*sigx_r - 2*(rho/sigx))
!!$       sigx_tf = (2*cs*rho + c2*(rho/sigx)*sigx_f)*(rho/sigx)/tau**3
!!$       sigx_td = c2/tau**3*(rho/sigx)**2*sigx_d
!!$       sigx_rr = c2/tau**2*(1 - (rho/sigx)*sigx_r)/sigx
!!$       sigx_rf = -(2*cs + c2*(sigx_f/sigx))*(rho/sigx)/tau**2
!!$       sigx_rd = -c2*(rho/sigx)*(sigx_d/sigx)/tau**2
!!$       sigx_ff = (-(2*c2-1)*rho + cs*(rho/sigx)*sigx_f)*(rho/sigx)/tau**2
!!$       sigx_fd = cs*sigx_d*(rho/sigx)**2/tau**2
!!$       sigx_dd = -(sigx_d/sigx)/2/sigx
!!$
!!$       sigy_tt = -(mu/sigy)*(sigy_t/sigy)/2
!!$       sigy_tr = -(mu/sigy)*(sigy_r/sigy)/2
!!$       sigy_tf = -(mu/sigy)*(sigy_f/sigy)/2
!!$       sigy_td = 1/(2*sigy) - (mu/sigy)*(sigy_d/sigy)/2
!!$       sigy_rr = s2*(1/sigy - (rho/sigy)*(sigy_r/sigy))
!!$       sigy_rf = (2*cs - s2*(sigy_f/sigy))*(rho/sigy)
!!$       sigy_rd = -s2*(rho/sigy)*(sigy_d/sigy)
!!$       sigy_ff = rho*(c2 - s2)*(rho/sigy) - cs*(rho/sigy)**2*sigy_f
!!$       sigy_fd = -cs*sigy_d*(rho/sigy)**2
!!$       sigy_dd = -(tau/2)*(sigy_d/sigy**2)


!       call grad_sigx(sigx,mu,tau,rho,cosfi,sinfi,sigx_t,sigx_r,sigx_f,sigx_d)
!       call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!       sigx_f = 0
!       call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)
!       sigy_f = 0

!       call jac_sigx(sigx,mu,tau,rho,c2,cs,sigx_t,sigx_r,sigx_f,sigx_d, &
!            sigx_tt,sigx_tr,sigx_tf,sigx_td,sigx_rr,sigx_rf,sigx_rd, &
!            sigx_ff,sigx_fd,sigx_dd)
!       call jac_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d, &
!            sigx_tt,sigx_tr,sigx_td,sigx_rr,sigx_rd,sigx_dd)
!       sigx_tf = 0
!       sigx_rf = 0
!       sigx_ff = 0
!       sigx_fd = 0

!       call jac_sigy(sigy,mu,sigy_t,sigy_r,sigy_d, &
!            sigy_tt,sigy_tr,sigy_td,sigy_rr,sigy_rd,sigy_dd)
!       sigy_tf = 0
!       sigy_rf = 0
!       sigy_ff = 0
!       sigy_fd = 0

       rx = (d / sigx)/sx
       ry = (r / sigy)/sy
       if( robust ) then
          psix = tukey(rx)
          psiy = tukey(ry)
          dpsix = dtukey(rx)
          dpsiy = dtukey(ry)
       else
          psix = rx
          psiy = ry
          dpsix = 1
          dpsiy = 1
       end if

       dtt = (dpsix/sx**2)*((d/sigx)*(sigx_t/sigx))**2 &
            -(psix/sx)*(d/sigx)*(-2*(sigx_t/sigx)**2+sigx_tt/sigx) &
            -(sigx_t/sigx)**2 + sigx_tt/sigx &
            +(dpsiy/sy**2)*(mu/sigy + (r/sigy)*(sigy_t/sigy))**2 &
            +(psiy/sy)*(2*(mu/sigy)*(sigy_t/sigy)+2*(r/sigy)*(sigy_t/sigy)**2 &
            -(r/sigy)*(sigy_tt/sigy)) -(sigy_t/sigy)**2 + sigy_tt/sigy

       dtr = (dpsix/sx**2)*(d/sigx)**2*(sigx_t/sigx)*(sigx_r/sigx) &
            +(psix/sx)*(d/sigx)*(2*(sigx_t/sigx)*(sigx_r/sigx) &
            -sigx_tr/sigx) - (sigx_t/sigx)*(sigx_r/sigx) + sigx_tr/sigx &
            +(dpsiy/sy**2)*(mu/sigy + (r/sigy)*(sigy_t/sigy)) &
            *(r/sigy)*(sigy_r/sigy) &
            +(psiy/sy)*((sigy_r/sigy)*(mu/sigy)-(r/sigy)*(sigy_tr/sigy) &
            +2*(r/sigy)*(sigy_t/sigy)*(sigy_r/sigy)) &
            -(sigy_t/sigy)*(sigy_r/sigy) + sigy_tr/sigy

!!$       dtf = (dpsix/sx**2)*(d/sigx)**2*(sigx_t/sigx)*(sigx_f/sigx) &
!!$            +(psix/sx)*(d/sigx)*(2*(sigx_t/sigx)*(sigx_f/sigx) &
!!$            -sigx_tf/sigx) - (sigx_t/sigx)*(sigx_f/sigx) + sigx_tf/sigx &
!!$            +(dpsiy/sy**2)*(mu/sigy + (r/sigy)*(sigy_t/sigy)) &
!!$            *(r/sigy)*(sigy_f/sigy) &
!!$            +(psiy/sy)*((sigy_f/sigy)*(mu/sigy)-(r/sigy)*(sigy_tf/sigy) &
!!$            +2*(r/sigy)*(sigy_t/sigy)*(sigy_f/sigy)) &
!!$            -(sigy_t/sigy)*(sigy_f/sigy) + sigy_tf/sigy

       dtd = (dpsix/sx**2)*(-d/sigx)*(sigx_t/sigx) &
            *(1/sigx-(d/sigx)*(sigx_d/sigx)) &
            +(psix/sx)*(d/sigx)*(2*(sigx_t/sigx)*(sigx_d/sigx)-sigx_td/sigx)&
           -(psix/sx)*(sigx_t/sigx**2)-(sigx_t/sigx)*(sigx_d/sigx)+sigx_td/sigx&
           +(dpsiy/sy**2)*(mu/sigy+(r/sigy)*(sigy_t/sigy)) &
           *(tau/sigy+(r/sigy)*(sigy_d/sigy)) &
           -(psiy/sy)*(1/sigy-(mu/sigy)*(sigy_d/sigy)-(tau/sigy)*(sigy_t/sigy)&
           -2*(r/sigy)*(sigy_t/sigy)*(sigy_d/sigy)+(r/sigy)*(sigy_td/sigy))&
           -(sigy_t/sigy)*(sigy_d/sigy) + sigy_td/sigy


       drr = (dpsix/sx**2)*((d/sigx)*(sigx_r/sigx))**2 &
            +(psix/sx)*(d/sigx)*(2*(sigx_r/sigx)**2-sigx_rr/sigx) &
            -(sigx_r/sigx)**2 + sigx_rr/sigx &
            +(dpsiy/sy**2)*((r/sigy)*(sigy_r/sigy))**2 &
            +(psiy/sy)*(r/sigy)*(2*(sigy_r/sigy)**2-(sigy_rr/sigy)) &
            -(sigy_r/sigy)**2 + sigy_rr/sigy

!!$       drf = (dpsix/sx**2)*(d/sigx)**2*(sigx_r/sigx)*(sigx_f/sigx) &
!!$            +(psix/sx)*(d/sigx)*(2*(sigx_r/sigx)*(sigx_f/sigx) - sigx_rf/sigx) &
!!$            -(sigx_r/sigx)*(sigx_f/sigx) + sigx_rf/sigx &
!!$            +(dpsiy/sy**2)*(r/sigy)**2*(sigy_r/sigy)*(sigy_f/sigy) &
!!$            +(psiy/sy)*(r/sigy)*(2*(sigy_r/sigy)*(sigy_f/sigy) - sigy_rf/sigy) &
!!$            -(sigy_r/sigy)*(sigy_f/sigy) + sigy_rf/sigy

       drd = (dpsix/sx**2)*(1/sigx-(d/sigx)*(sigx_d/sigx)) &
            *(-d/sigx)*(sigx_r/sigx) &
            +(psix/sx)*(-sigx_r/sigx**2+2*(d/sigx)*(sigx_r/sigx)*(sigx_d/sigx) &
            -(d/sigx)*(sigx_rd/sigx)) &
            - (sigx_r/sigx)*(sigx_d/sigx) + sigx_rd/sigx &
            +(dpsiy/sy**2)*(tau/sigy + (r/sigy)*(sigy_d/sigy)) &
            *(r/sigy)*(sigy_r/sigy) &
            +(psiy/sy)*((tau/sigy)*(sigy_r/sigy)-(r/sigy)*(sigy_rd/sigy) &
            +2*(r/sigy)*(sigy_r/sigy)*(sigy_d/sigy)) &
            -(sigy_r/sigy)*(sigy_d/sigy) + sigy_rd/sigy

!!$       dff = (dpsix/sx**2)*((d/sigx)*(sigx_f/sigx))**2 &
!!$            +(psix/sx)*(d/sigx)*(2*(sigx_f/sigx)**2 - sigx_ff/sigx) &
!!$            -(sigx_f/sigx)**2 + sigx_ff/sigx &
!!$            +(dpsiy/sy**2)*((r/sigy)*(sigy_f/sigy))**2 &
!!$            +(psiy/sy)*(r/sigy)*(2*(sigy_f/sigy)**2 - (sigy_ff/sigy)) &
!!$            -(sigy_f/sigy)**2 + sigy_ff/sigy
!!$
!!$
!!$       dfd = (dpsix/sx**2)*(1/sigx-(d/sigx)*(sigx_d/sigx)) &
!!$            *(-d/sigx)*(sigx_f/sigx) &
!!$            +(psix/sx)*(-sigx_f/sigx**2+2*(d/sigx)*(sigx_f/sigx)*(sigx_d/sigx) &
!!$            -(d/sigx)*(sigx_fd/sigx)) &
!!$            - (sigx_f/sigx)*(sigx_d/sigx) + sigx_fd/sigx &
!!$            +(dpsiy/sy**2)*(tau/sigy + (r/sigy)*(sigy_d/sigy)) &
!!$            *(r/sigy)*(sigy_f/sigy) &
!!$            +(psiy/sy)*((tau/sigy)*(sigy_f/sigy) &
!!$            +2*(r/sigy)*(sigy_f/sigy)*(sigy_d/sigy) &
!!$            -(r/sigy)*(sigy_fd/sigy)) &
!!$            -(sigy_f/sigy)*(sigy_d/sigy) + sigy_fd/sigy

       ddd = (dpsix/sx**2)*(1/sigx-(d/sigx)*(sigx_d/sigx))**2 &
            -(psix/sx)*(2*sigx_d/sigx**2 - 2*(d/sigx)*(sigx_d/sigx)**2 &
            +(d/sigx)*(sigx_dd/sigx)) &
            -(sigx_d/sigx)**2 + sigx_dd/sigx &
            +(dpsiy/sy**2)*(tau/sigy + (r/sigy)*(sigy_d/sigy))**2 &
            +(psiy/sy)*(2*(tau/sigy)*(sigy_d/sigy)+2*(r/sigy)*(sigy_d/sigy)**2 &
            -(r/sigy)*(sigy_dd/sigy)) &
            -(sigy_d/sigy)**2 + sigy_dd/sigy


       jac(1,1:2) = jac(1,1:2) + [dtt, dtr ]
       jac(2,1:2) = jac(2,1:2) + [dtr, drr ]
!       jac(3,:) = jac(3,:) + [dtf, drf, dff ]
       jac(1,1+i) = dtd
       jac(1+i,1) = dtd
       jac(2,1+i) = drd
       jac(1+i,2) = drd
!       jac(3,2+i) = dfd
!       jac(2+i,2) = dfd
       jac(i+1,i+1) = ddd

    end do
    jac = jac / n

  end function jacf


  subroutine thefline(tau,dtau,rho,drho,phi,dphi,sx,sy,delta,reli)

    use minpacks

    real(REAL64), intent(in out) :: tau, rho, phi, sx, sy
    real(REAL64), intent(out) :: dtau, drho, dphi
    real(REAL64), dimension(:), intent(out) :: delta
    logical, intent(out) :: reli

    integer :: maxfev = 666
    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(:), allocatable :: x, diag
    real(REAL64), dimension(:,:), allocatable :: cov
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: tol
    integer :: npar, info, n, i, K

    n = size(mju)
    K = size(kmean)
    npar = 3 + K
    maxfev = 666*npar
    allocate(x(npar),diag(npar),cov(npar,npar))

    backsx = sx
    backsy = sy

    x(1) = tau
    x(2) = rho
    x(3) = phi
    x(4:) = delta
    tol = 10*n*epsilon(x)*(1 + abs(maxval(x)))
    diag(1:3) = [ 0.001,1e3, 0.01 ]
    diag(4:) = sqrt((pht+delta)*(1+1/tau) + (rho/tau)**2)
    diag = 1 / diag**2


    if( analytic ) then
       label = 'lmder'
!       call lmder3(flike,x,cov,eps,maxfev,factor,nprint,info,msg)
    else
       label = 'lmdif'
       call lmdif3(fliked,x,cov,tol,maxfev,factor=factor,nprint=nprint, &
            info=info,msg=msg)
    end if
    diag = [ (cov(i,i),i=1,npar) ]
    diag(3) = 1

    reli = tmin < x(1) .and. x(1) < tmax &
         .and. abs(x(2)) < rhomax & !.and. abs(x(3)) < rhomax &
!         .and. 0.1 < x(4) .and. x(4) < 10 &
!         .and. 0.1 < x(5) .and. x(5) < 10 &
         .and. all(diag > 0) .and. (0 < info .and. info < 5)

    if( .not. all(diag > 0) .and. .false.) then
       block
         integer :: i
         real(REAL64), dimension(npar,npar) :: fjac
!         fjac = jacf1(x(1),siguback,sigvback,su,sv,x(2:))

         do i = 1, K
            if( .not. (diag(i) > 0) .or. i == 1 ) &
                 write(*,*) i,real(diag(i)),real(cts(i)),real(pht(i)),real(mju(i)),real(fjac(i,i))
         end do
       end block
    end if

    if( debug_iter ) then
       write(error_unit,'(3a,i0,3a,l2)') ' theline() ',label,' info=', &
            info,': ',trim(msg),'  reliable:',reli
    end if

    if( reli ) then
       ! Sometimes, the lambda parameter of LM-method is incidentally found
       ! in the mirrored part of the parameter space. Fortunately, the function
       ! is pretty symmetric, so I do to ignore the "tunnelling". The correct,
       ! in principle, solution should to keep the condition: sigv > 0.
       tau = x(1)
       rho = x(2)
       phi = x(3)
!       su = x(4)
!       sv = x(5)
       delta = x(4:)

       dtau =  sqrt(cov(1,1)/(2*n - npar))
       drho = sqrt(cov(2,2)/(2*n - npar))
!       dphi = sqrt(cov(3,3)/(2*n - npar))
    end if

    if( debug_iter ) then
       block
         integer, parameter :: maxpar = 7
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64) :: q
         integer :: i,j,n
         character(len=3) :: dots
         character(len=80) :: fmt

         cor = 1
         do i = 1, npar
            do j = i+1, npar
               if( cov(i,i) > 0 .and. cov(j,j) > 0 ) then
                  q = sqrt(cov(i,i)*cov(j,j))
                  if( abs(q) > 0 ) then
                     cor(i,j) = cov(i,j) / q
                  else
                     cor(i,j) = huge(cor)
                  end if
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         n = min(maxpar,npar)
         if( n < npar ) then
            dots = '...'
         else
            dots = ''
         end if

         write(error_unit,*) 'Covariance (Hessian in minimum):'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(es9.2,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cov(i,1:n),dots
         end do

         write(error_unit,*) 'Correlation:'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(f6.3,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cor(i,1:n),dots
         end do

       end block
    end if


  end subroutine thefline

  subroutine flike(npar,mpar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, mpar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

!    real(REAL64), dimension(4) :: vec
!    real(REAL64), dimension(4,4) :: jac
    !    real(REAL64) :: tau, sigu, sigv
    real(REAL64), dimension(2) :: c
    real(REAL64), dimension(:), allocatable :: su_d, sv_d
    real(REAL64) :: su_t, su_u, sv_v
    character(len=3) :: dots

!    n = size(delta)
!    tau = x(1)
!    sigu = x(2)
!    sigv = x(3)
!    s = x(4)
!    delta = x(2:)

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          if( size(x) > 8 ) dots = '...'
          write(error_unit,*) 's0 = ',real(sum(x(2:)**2)/(size(x)-1)),&
               real(sqrt(sum(fvec(2:)**2)/(size(fvec)-1)))
          write(error_unit,*) 'x = ',x(1:8),trim(dots)
          write(error_unit,*) 'fvec = ',fvec(1:8),trim(dots)
          block
            integer :: i,n
            character(len=80) :: fmt
            n = min(7,npar)
            write(fmt,'(a,i0,a)') '(a,i0,a,',n,'(1x,es9.2),a)'
            do i = 1, n
               write(error_unit,fmt) &
                    'fjac(',i,',:)=',fjac(i,1:n),trim(dots)
            end do
          end block
       end if

       ! convergence check
!       if( tau < tmin .or. tau > tmax ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

!       sv = 1
!       su = 1
!       fvec = fgradf(x(1),x(2),x(3),x(4),x(5),x(6:))

!       su_t = - fvec(4) / fvec(1)
!       su_u = - fvec(4) / fvec(2)
!       sv_v = - fvec(5) / fvec(3)
!       su_d = - fvec(4) / fvec(6:)
!       sv_d = - fvec(5) / fvec(6:)
       c = x(4:5) - scale

!       fvec(1) = fvec(1) - lambda(1)*su_t + nu*c(1)*su_t
!       fvec(2) = fvec(2) - lambda(1)*su_u + nu*c(1)*su_u
!       fvec(3) = fvec(3) - lambda(2)*sv_v + nu*c(2)*sv_v
!       fvec(4) = fvec(4) - lambda(1) + nu*c(1)
!       fvec(5) = fvec(5) - lambda(2) + nu*c(2)
!       fvec(6:) = fvec(6:) - lambda(1)*su_d - lambda(2)*sv_d  &
!            + nu*(c(1)*su_d + c(2)*sv_d)


!       su = ydisp(x(1),x(2:))
!       c = [ su, sv ] - scale

 !      fvec(1) = fvec(1) - lambda(1)*su_t + nu*c(1)*su_t
 !      fvec(2:) = fvec(2:) - lambda(1)*su_d - lambda(2)*sv_d &
 !           +(nu/2)*(c(1)*su_d + c(2)*sv_d)

!       fvec = gradf(tau,sigu,sigv,s)
!       fvec(4) = fvec(4) - lambda + nu*(s - scale)

    else if( iflag == 2 ) then

!       sv = 1
!       su = 1
!       fjac = jacf(tau,sigu,sigv,s)
!       fjac = jacf(x(1),x(2),x(3),x(4),x(5),x(6:))

!       fjac(1,1) = fjac(1,1) - lambda(1)*su_tt + nu*(su_t**2 - c(1)*su_tt)
!       do i = 1, n
!          fjac(i,i) = fjac(i,i) - lambda(1)*su_dd(i) -lambda(2)*sv_dd(i) &
!          + nu*(su_d(i)**2 - c(1)*su_dd(i)) + sv_d(i)**2 - c(2)*sv_dd(i))
!       end do
!       fjac(1,2:) = fjac(1,2:) &
!            - lambda(1)*su_td + nu*(su_t*su_d - c(1)*su_td)

    end if

  end subroutine flike

  subroutine fliked(npar,mpar,x,fvec,iflag)

    integer, intent(in) :: npar, mpar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(in out) :: fvec
!    real(REAL64) :: su,sv
    real(REAL64), dimension(2) :: c
    real(REAL64), dimension(:,:), allocatable :: dc
!    real(REAL64) :: su_t, su_u, sv_v
    character(len=3) :: dots

!    su = 1
!    sv = 1

!    real(REAL64) :: tau, sigu, sigv, s

!    tau = x(1)
!    sigu = x(2)
!    sigv = x(3)
!    s = x(4)

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          if( size(x) > 8 ) dots = '...'
          write(error_unit,*) 's0 = ',real(sqrt(sum(x(4:)**2)/(size(x)-1))), &
               ' loglike=',real(floglike(x,backsx,backsy))
          write(error_unit,*) 'x = ',real(x(1:8)),trim(dots)
          write(error_unit,*) 'fvec = ',real(fvec(1:8)),trim(dots)
       end if

       ! convergence check ?
       return
    end if

!    write(error_unit,*) 'x = ',real(x),trim(dots)


    fvec = fgradf(x(1),x(2),x(3),x(4:),backsx,backsy)

!    c = disp(x(1),x(2),x(3),x(4:),backsx,backsy)
    dc = gradc(x(1),x(2),x(3),x(4:),backsx,backsy)

!    fvec = fvec + dc(1,:)*(nu*c(1) - lambda(1)) + dc(2,:)*(nu*c(2) - lambda(2))
!    fvec = fvec + sum(dc*(nu*c - lambda),1)


!!$    block
!!$      integer :: i, l
!!$      real(REAL64), dimension(size(x)) :: fd,eps,z
!!$      real(REAL64), dimension(-1:1) :: f
!!$
!!$      eps = epsilon(eps)**(1.0/3.0)
!!$
!!$      do i = 1, size(fvec)
!!$         do l = -1,1,2
!!$            z = x
!!$            z(i) = z(i) + l*eps(i)
!!$            c = z(4:5) - scale
!!$            f(l) = loglike(z) !- sum(lambda*c) + (nu/2)*sum(c**2)
!!$         end do
!!$         fd(i) = (f(1) - f(-1)) / (2*eps(i))
!!$         write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
!!$              real(fd(i)-fvec(i)),real(eps(i))
!!$      end do
!!$      fvec = fd
!!$
!!$    end block


!    write(error_unit,*) 'fvec = ',real(fvec),trim(dots)


!    write(*,*) lambda, nu
!    write(*,*)
!    stop


!    su_t = - fvec(1) / fvec(4)
!    su_u = - fvec(2) / fvec(4)
!    sv_v = - fvec(3) / fvec(5)
!    su_d = - fvec(6:) / fvec(4)
!    sv_d = - fvec(6:) / fvec(5)
!    c = x(4:5) - scale

!    fvec(1) = fvec(1) - lambda(1)*su_t + nu*c(1)*su_t
!    fvec(2) = fvec(2) - lambda(1)*su_u + nu*c(1)*su_u
!    fvec(3) = fvec(3) - lambda(2)*sv_v + nu*c(2)*sv_v
!   fvec(4) = fvec(4) - lambda(1) + nu*c(1)
!   fvec(5) = fvec(5) - lambda(2) + nu*c(2)
!    fvec(6:) = fvec(6:) - lambda(1)*su_d - lambda(2)*sv_d  &
!         + nu*(c(1)*su_d + c(2)*sv_d)


!    fvec = gradf(x(1),x(2),x(3),x(4),x(5),x(6:))
!    fvec(4) = fvec(4) - lambda(1) + nu*(x(4) - scale)
!    fvec(5) = fvec(5) - lambda(2) + nu*(x(5) - scale)
!    fvec = gradf(x(1),x(2),x(3),x(4),x(5),x(6:))
!    fvec(4) = fvec(4) - sum(lambda) + nu*(s - scale)

  end subroutine fliked

  function fgradf(tau,rho,phi,delta,sx,sy) result(vec)

    use rfun

    real(REAL64), intent(in) :: tau,rho,phi,sx,sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(3 + size(delta)) :: vec
    real(REAL64) :: sigu,sigv,sigu2, sigv2, mu,d, t,r,sigx,sigy,psix,psiy, &
         sigy_t,sigx_t,sigy_r,sigx_r, sigy_f,sigx_f, sigy_d, sigx_d, &
         sinfi, cosfi, dt,dr,df,ds,dd,bu,bv
    integer :: i, n

    n = size(pht)

    sinfi = sin(phi)
    cosfi = cos(phi)
    sigu = rho*sinfi
    sigv = rho*(cosfi/tau)
    sigu2 = sigu**2
    sigv2 = sigv**2

    vec = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - mu*tau
       sigx = sqrt(mu + sigv2)
       sigy = sqrt(mu*tau + sigu2)
       sigy_t = (mu / sigy) / 2
       sigx_t = -(cosfi**2/tau**3)*(rho/sigx) * rho
       sigy_r = (rho/sigy)*sinfi**2
       sigx_r = (rho/sigx)*(cosfi/tau)**2
       sigy_f = (sinfi*cosfi) * (rho/sigy) * rho
       sigx_f = -(sinfi*cosfi)/tau**2 * (rho/sigx) * rho

       if( robust ) then
          psix = tukey((d / sigx)/sx)
          psiy = tukey((r / sigy)/sy)
       else
          psix = (d / sigx) / sx
          psiy = (r / sigy) / sy
       end if

       dt = -(psiy/sy)*(mu/sigy + (r/sigy)*(sigy_t/sigy)) + sigy_t/sigy &
            -(psix/sx)*(d/sigx)*(sigx_t/sigx) + sigx_t /sigx
       dr = -(psix/sx)*(d/sigx)*(sigx_r/sigx) + sigx_r/sigx &
            -(psiy/sy)*(r/sigy)*(sigy_r/sigy) + sigy_r/sigy
       df = -(psix/sx)*(d/sigx)*(sigx_f/sigx) + sigx_f/sigx &
            -(psiy/sy)*(r/sigy)*(sigy_f/sigy) + sigy_f/sigy
!       bv = (1 - (psix/sv)*(d/sigx))/sv
!       bu = (1 - (psiy/su)*(r/sigy))/su

       sigx_d = 1/(2*sigx)
       sigy_d = tau / sigy / 2
       dd =  (psix/sx)*(1/sigx - (d/sigx)*(sigx_d/sigx)) &
            -(psiy/sy)*(tau/sigy + (r/sigy)*(sigy_d/sigy)) &
            + sigx_d/sigx + sigy_d/sigy

       vec(1:3) = vec(1:3) + [ dt, dr, df ]!, bu, bv ]
       vec(3+i) = dd

    end do
    vec = vec / n

  end function fgradf



  subroutine theline0(tau,dtau,sigu,dsigu,sigv,dsigv,su,sv,delta,reli)

    use minpacks

    real(REAL64), intent(in out) :: tau, sigu, sigv, su, sv
    real(REAL64), intent(out) :: dtau, dsigu, dsigv
    real(REAL64), dimension(:), intent(out) :: delta
    logical, intent(out) :: reli

    integer :: maxfev = 666
    real(REAL64), parameter :: factor = 1
    real(REAL64), dimension(:), allocatable :: x, diag
    real(REAL64), dimension(:,:), allocatable :: cov
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: eps
    integer :: npar, info, n, i, K

    n = size(mju)
    K = size(kmean)
    npar = 5 + K
!    maxfev = 666*npar
    allocate(x(npar),diag(npar),cov(npar,npar))

    siguback = sigu
    sigvback = sigv


    x(1) = tau
    x(2) = sigu
    x(3) = sigv
    x(4) = su
    x(5) = sv
    x(6:) = mju - pht
    eps = 10*n*epsilon(x)*(1 + abs(maxval(x)))

    if( analytic ) then
       label = 'lmder'
!       call lmder3(like,x,cov,eps,maxfev,factor,nprint,info,msg)
    else
       label = 'lmdif'
!       call lmdif3(liked,x,cov,eps,maxfev,factor,nprint,info,msg)
    end if
    diag = [ (cov(i,i),i=1,npar) ]

    reli = tmin < x(1) .and. x(1) < tmax &
         .and. abs(x(2)) < rhomax .and. abs(x(3)) < rhomax &
         .and. 0.1 < x(4) .and. x(4) < 10 &
         .and. 0.1 < x(5) .and. x(5) < 10 &
         .and. all(diag > 0) .and. (0 < info .and. info < 5)

    if( .not. all(diag > 0) .and. .false. ) then
       block
         integer :: i
         real(REAL64), dimension(npar,npar) :: fjac
         fjac = jacf1(x(1),siguback,sigvback,su,sv,x(2:))

         do i = 1, K
            if( .not. (diag(i) > 0) .or. i == 1 ) &
                 write(*,*) i,real(diag(i)),real(cts(i)),real(pht(i)),real(mju(i)),real(fjac(i,i))
         end do
       end block
    end if

    if( debug_iter ) then
       write(error_unit,'(3a,i0,3a,l2)') ' theline() ',label,' info=', &
            info,': ',trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       ! Sometimes, the lambda parameter of LM-method is incidentally found
       ! in the mirrored part of the parameter space. Fortunately, the function
       ! is pretty symmetric, so I do to ignore the "tunnelling". The correct,
       ! in principle, solution should to keep the condition: sigv > 0.
       tau = x(1)
       sigu = abs(x(2))
       sigv = abs(x(3))
       su = x(4)
       sv = x(5)
       delta = x(6:)

       dtau =  sqrt(cov(1,1)/(2*n - npar))
       dsigu = sqrt(cov(2,2)/(2*n - npar))
       dsigv = sqrt(cov(3,3)/(2*n - npar))
    end if

    if( debug_iter ) then
       block
         integer, parameter :: maxpar = 7
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64) :: q
         integer :: i,j,n
         character(len=3) :: dots
         character(len=80) :: fmt

         cor = 1
         do i = 1, npar
            do j = i+1, npar
               if( cov(i,i) > 0 .and. cov(j,j) > 0 ) then
                  q = sqrt(cov(i,i)*cov(j,j))
                  if( abs(q) > 0 ) then
                     cor(i,j) = cov(i,j) / q
                  else
                     cor(i,j) = huge(cor)
                  end if
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         n = min(maxpar,npar)
         if( n < npar ) then
            dots = '...'
         else
            dots = ''
         end if

         write(error_unit,*) 'Covariance (Hessian in minimum):'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(es9.2,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cov(i,1:n),dots
         end do

         write(error_unit,*) 'Correlation:'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(f6.3,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cor(i,1:n),dots
         end do

       end block
    end if


  end subroutine theline0

  subroutine like0(npar,mpar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, mpar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

!    real(REAL64), dimension(4) :: vec
!    real(REAL64), dimension(4,4) :: jac
    !    real(REAL64) :: tau, sigu, sigv
    real(REAL64), dimension(2) :: c
    real(REAL64), dimension(:), allocatable :: su_d, sv_d
    real(REAL64) :: su_t, su_u, sv_v
    character(len=3) :: dots

!    n = size(delta)
!    tau = x(1)
!    sigu = x(2)
!    sigv = x(3)
!    s = x(4)
!    delta = x(2:)

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          if( size(x) > 8 ) dots = '...'
          write(error_unit,*) 's0 = ',real(sum(x(2:)**2)/(size(x)-1)),&
               real(sqrt(sum(fvec(2:)**2)/(size(fvec)-1)))
          write(error_unit,*) 'x = ',x(1:8),trim(dots)
          write(error_unit,*) 'fvec = ',fvec(1:8),trim(dots)
          block
            integer :: i,n
            character(len=80) :: fmt
            n = min(7,npar)
            write(fmt,'(a,i0,a)') '(a,i0,a,',n,'(1x,es9.2),a)'
            do i = 1, n
               write(error_unit,fmt) &
                    'fjac(',i,',:)=',fjac(i,1:n),trim(dots)
            end do
          end block
       end if

       ! convergence check
!       if( tau < tmin .or. tau > tmax ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

!       sv = 1
!       su = 1
       fvec = gradf0(x(1),x(2),x(3),x(4),x(5),x(6:))

       su_t = - fvec(4) / fvec(1)
       su_u = - fvec(4) / fvec(2)
       sv_v = - fvec(5) / fvec(3)
       su_d = - fvec(4) / fvec(6:)
       sv_d = - fvec(5) / fvec(6:)
       c = x(4:5) - scale

!       fvec(1) = fvec(1) - lambda(1)*su_t + nu*c(1)*su_t
!       fvec(2) = fvec(2) - lambda(1)*su_u + nu*c(1)*su_u
!       fvec(3) = fvec(3) - lambda(2)*sv_v + nu*c(2)*sv_v
!       fvec(4) = fvec(4) - lambda(1) + nu*c(1)
!       fvec(5) = fvec(5) - lambda(2) + nu*c(2)
!       fvec(6:) = fvec(6:) - lambda(1)*su_d - lambda(2)*sv_d  &
!            + nu*(c(1)*su_d + c(2)*sv_d)


!       su = ydisp(x(1),x(2:))
!       c = [ su, sv ] - scale

 !      fvec(1) = fvec(1) - lambda(1)*su_t + nu*c(1)*su_t
 !      fvec(2:) = fvec(2:) - lambda(1)*su_d - lambda(2)*sv_d &
 !           +(nu/2)*(c(1)*su_d + c(2)*sv_d)

!       fvec = gradf(tau,sigu,sigv,s)
!       fvec(4) = fvec(4) - lambda + nu*(s - scale)

    else if( iflag == 2 ) then

!       sv = 1
!       su = 1
!       fjac = jacf(tau,sigu,sigv,s)
       fjac = jacf0(x(1),x(2),x(3),x(4),x(5),x(6:))

!       fjac(1,1) = fjac(1,1) - lambda(1)*su_tt + nu*(su_t**2 - c(1)*su_tt)
!       do i = 1, n
!          fjac(i,i) = fjac(i,i) - lambda(1)*su_dd(i) -lambda(2)*sv_dd(i) &
!          + nu*(su_d(i)**2 - c(1)*su_dd(i)) + sv_d(i)**2 - c(2)*sv_dd(i))
!       end do
!       fjac(1,2:) = fjac(1,2:) &
!            - lambda(1)*su_td + nu*(su_t*su_d - c(1)*su_td)

    end if

  end subroutine like0

  subroutine liked0(npar,mpar,x,fvec,iflag)

    integer, intent(in) :: npar, mpar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec
!    real(REAL64) :: su,sv
    real(REAL64), dimension(2) :: c
    real(REAL64), dimension(:), allocatable :: su_d, sv_d
    real(REAL64) :: su_t, su_u, sv_v
    character(len=3) :: dots

!    su = 1
!    sv = 1

!    real(REAL64) :: tau, sigu, sigv, s

!    tau = x(1)
!    sigu = x(2)
!    sigv = x(3)
!    s = x(4)

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          if( size(x) > 8 ) dots = '...'
          write(error_unit,*) 's0 = ',real(sum(x(2:)**2)/(size(x)-1))!, &
!               ' loglike=',real(loglike(x,))
          write(error_unit,*) 'x = ',real(x(1:8)),trim(dots)
          write(error_unit,*) 'fvec = ',real(fvec(1:8)),trim(dots)
       end if

       ! convergence check ?
       return
    end if

!    write(error_unit,*) 'x = ',real(x),trim(dots)


    fvec = gradf0(x(1),x(2),x(3),x(4),x(5),x(6:))

!!$    block
!!$      integer :: i, l
!!$      real(REAL64), dimension(size(x)) :: fd,eps,z
!!$      real(REAL64), dimension(-1:1) :: f
!!$
!!$      eps = epsilon(eps)**(1.0/3.0)
!!$
!!$      do i = 1, size(fvec)
!!$         do l = -1,1,2
!!$            z = x
!!$            z(i) = z(i) + l*eps(i)
!!$            c = z(4:5) - scale
!!$            f(l) = loglike(z) !- sum(lambda*c) + (nu/2)*sum(c**2)
!!$         end do
!!$         fd(i) = (f(1) - f(-1)) / (2*eps(i))
!!$         write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
!!$              real(fd(i)-fvec(i)),real(eps(i))
!!$      end do
!!$      fvec = fd
!!$
!!$    end block


!    write(error_unit,*) 'fvec = ',real(fvec),trim(dots)


!    write(*,*) lambda, nu
!    write(*,*)
!    stop


    su_t = - fvec(1) / fvec(4)
    su_u = - fvec(2) / fvec(4)
    sv_v = - fvec(3) / fvec(5)
    su_d = - fvec(6:) / fvec(4)
    sv_d = - fvec(6:) / fvec(5)
    c = x(4:5) - scale

!    fvec(1) = fvec(1) - lambda(1)*su_t + nu*c(1)*su_t
!    fvec(2) = fvec(2) - lambda(1)*su_u + nu*c(1)*su_u
!    fvec(3) = fvec(3) - lambda(2)*sv_v + nu*c(2)*sv_v
!   fvec(4) = fvec(4) - lambda(1) + nu*c(1)
!   fvec(5) = fvec(5) - lambda(2) + nu*c(2)
!    fvec(6:) = fvec(6:) - lambda(1)*su_d - lambda(2)*sv_d  &
!         + nu*(c(1)*su_d + c(2)*sv_d)


!    fvec = gradf(x(1),x(2),x(3),x(4),x(5),x(6:))
!    fvec(4) = fvec(4) - lambda(1) + nu*(x(4) - scale)
!    fvec(5) = fvec(5) - lambda(2) + nu*(x(5) - scale)
!    fvec = gradf(x(1),x(2),x(3),x(4),x(5),x(6:))
!    fvec(4) = fvec(4) - sum(lambda) + nu*(s - scale)

  end subroutine liked0

  function gradf0(tau,sigu,sigv,su,sv,delta) result(vec)

    use rfun

    real(REAL64), intent(in) :: tau,sigu,sigv,su,sv
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(5 + size(delta)) :: vec
    real(REAL64) :: sigu2, sigv2, mu,d, t,r,sigx,sigy,psix,psiy, &
         sigy_t,sigy_u,sigx_v,sigy_d, sigx_d, dt,du,dv,ds,dd,bu,bv
    integer :: i, n

    n = size(pht)

    sigu2 = sigu**2
    sigv2 = sigv**2

    vec = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - mu*tau
       sigx = sqrt(mu + sigv2)
       sigy = sqrt(mu*tau + sigu2)
       sigy_t = (mu / sigy) / 2
       sigy_u = sigu / sigy
       sigx_v = sigv / sigx

       if( robust ) then
          psix = tukey((d / sigx)/sv)
          psiy = tukey((r / sigy)/su)
       else
          psix = (d / sigx) / sv
          psiy = (r / sigy) / su
       end if

       dt = -(psiy/su)*(mu/sigy + (r/sigy)*(sigy_t/sigy)) + sigy_t/sigy
       du = sigy_u/sigy*(1 - (psiy/su)*(r/sigy))
       dv = sigx_v/sigx*(1 - (psix/sv)*(d/sigx))
       bv = (1 - (psix/sv)*(d/sigx))/sv
       bu = (1 - (psiy/su)*(r/sigy))/su

       sigx_d = 1/(2*sigx)
       sigy_d = tau / sigy / 2
       dd =  (psix/sv)*(1/sigx - (d/sigx)*(sigx_d/sigx)) &
            -(psiy/su)*(tau/sigy + (r/sigy)*(sigy_d/sigy)) &
            + sigx_d/sigx + sigy_d/sigy

       vec(1:5) = vec(1:5) + [ dt, du, dv, bu, bv ]
       vec(5+i) = dd

    end do
    vec = vec / n

  end function gradf0

  function jacf0(tau,sigu,sigv,su,sv,delta) result(jac)

    use rfun

    real(REAL64), intent(in) :: tau,sigu,sigv,su,sv
    real(REAL64), dimension(:), intent(in) :: delta

    real(REAL64), parameter :: z = real(0.0,REAL64)
    real(REAL64), dimension(5+size(delta),5+size(delta)) :: jac
    real(REAL64) :: sigu2, sigv2,d,r,mu,sigx,sigy,psix,psiy,dpsix,dpsiy, &
         sigy_t,sigy_u,sigx_v,sigy_tt,dtt,duu,dvv,dss,dtu,dts,dsu,dsv,&
         sigx_vv, sigy_uu, dbu, dbv, ddd, dtd
    integer :: i, n

    n = size(delta)
    sigu2 = sigu**2
    sigv2 = sigv**2

    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - mu*tau
       sigx = sqrt(mu + sigv2)
       sigy = sqrt(mu*tau + sigu2)
       sigy_t = mu / (2*sigy)
       sigy_u = sigu / sigy
       sigx_v = sigv / sigx

       sigy_tt = - (mu / (2*sigy)) * (sigy_t/sigy)

       if( robust ) then
          psix = tukey(d/(sv*sigx))
          psiy = tukey(r/(su*sigy))
          dpsix = dtukey(d/(sv*sigx))
          dpsiy = dtukey(r/(su*sigy))
       else
          psix = d / (sv*sigx)
          psiy = r / (su*sigy)
          dpsix = 1
          psiy = 1
       end if

       dtt = dpsiy*(mu/sigy + (r/sigy)*(sigy_t/sigy))**2 &
            +psiy*((mu/sigy)*(sigy_t/sigy) + &
            2*(mu/sigy**2 + (r/sigy)*(sigy_t))*(sigy_t/sigy)) &
            - (sigy_t/sigy)**2 + sigy_tt/sigy
       duu = dpsiy*((r/sigy)*(sigy_u))**2 &
            -psiy*(r/sigy)*(2*(sigy_u/sigy)**2 + sigy_u/sigy) &
            - (sigy_u/sigy)**2 + (sigy_uu/sigy)
       dvv = dpsix*((d/sigx)*(sigx_v))**2 &
            -psix*(d/sigx)*(2*(sigx_v/sigx)**2 + sigx_v/sigx) &
            - (sigx_v/sigx)**2 + (sigx_vv/sigx)
       dss = (dpsix*(d/sigx)**2 + dpsiy*(r/sigy)**2)/sv**4 &
            -2*(psix*(d/sigx) + psiy*(r/sigy))/sv**3 - 2/sv**2
       dtu = 1
       dts = 1
       dsu = 1
       dsv = 1
       dbu = 1
       dbv = 1

       jac(1,:) = jac(1,:) + [dtt, dtu, z,   dts, dtu ]
       jac(2,:) = jac(2,:) + [dtu, duu, z,   dsu, dtu ]
       jac(3,:) = jac(3,:) + [  z,   z, dvv, dsv, dtu ]
       jac(4,:) = jac(4,:) + [dts, dsu, dsv, dss, dtu ]
       jac(5,:) = jac(5,:) + [dts, dsu, dsv, dss, dtu ]
       jac(5,5+i) = dtd
       jac(5+i,5) = dtd
       jac(i+5,i+5) = ddd

    end do
    jac = jac / n

  end function jacf0



  subroutine theline1(tau,dtau,delta,reli)

    use minpacks

    real(REAL64), intent(in out) :: tau
    real(REAL64), dimension(:), intent(in out) :: delta
    real(REAL64), intent(out) :: dtau!, ddelta
    logical, intent(out) :: reli

    integer :: maxfev = 666
    real(REAL64), parameter :: factor = 1
    real(REAL64), dimension(:), allocatable :: x, diag
    real(REAL64), dimension(:,:), allocatable :: cov
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: tol
    integer :: npar, info, n, i, K

    n = size(mju)
    K = size(kmean)
    npar = 1 + K
    maxfev = 666*npar
    allocate(x(npar),diag(npar),cov(npar,npar))

!    siguback = sigu
!    sigvback = sigv


    x(1) = tau
    x(2:) = delta
!    x(2) = sigu
!    x(3) = sigv
!    x(4) = su
!    x(5) = sv
!    do concurrent (i = 1:K)
!       x(1+i) = (-1)**i*0.1*sqrt(pht(i))
!    end do
!    x(2:) = 0
    tol = 10*n*epsilon(x)*(1 + abs(maxval(x)))

    if( analytic ) then
       label = 'lmder'
!       call lmder3(like1,x,cov,eps,maxfev,factor,nprint,info,msg)
    else
       label = 'lmdif'
       call lmdif3(liked1,x,cov,tol,maxfev,factor=factor,nprint=nprint,&
            info=info,msg=msg)
    end if
    diag = [ (cov(i,i),i=1,npar) ]

    reli = tmin < x(1) .and. x(1) < tmax &
!         .and. abs(x(2)) < rhomax .and. abs(x(3)) < rhomax &
!         .and. 0.1 < x(4) .and. x(4) < 10 &
!         .and. 0.1 < x(5) .and. x(5) < 10 &
         .and. all(abs(diag) > tol) .and. (0 < info .and. info < 5)

    if( .not. all(diag > 0) .and. .false.) then
       block
         integer :: i
         real(REAL64), dimension(npar,npar) :: fjac
!         fjac = jacf1(x(1),siguback,sigvback,su,sv,x(2:))

         do i = 1, K
            if( .not. (diag(i) > 0) .or. i == 1 ) &
                 write(*,*) i,real(diag(i)),real(cts(i)),real(pht(i)),real(mju(i)),real(fjac(i,i))
         end do
       end block
    end if

!    write(*,*) real(diag)

    if( debug_iter ) then
!!$       msg = ''
!!$       if( 0 < info .and. info < 5 ) then
!!$          msg = ' (success)'
!!$       else if( info == -1 ) then
!!$          msg = ' (out of range)'
!!$       else if( info == 5 ) then
!!$          msg = ' (no convergence)'
!!$       else if( .not. all(diag > 0) ) then
!!$          msg = ' (no minimum)'
!!$       end if
       write(error_unit,'(3a,i0,a,a,l2)') ' theline() ',label,' info: ', &
            info,trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       ! Sometimes, the lambda parameter of LM-method is incidentally found
       ! in the mirrored part of the parameter space. Fortunately, the function
       ! is pretty symmetric, so I do to ignore the "tunnelling". The correct,
       ! in principle, solution should to keep the condition: sigv > 0.
       tau = x(1)
!       sigu = abs(x(2))
!       sigv = abs(x(3))
!       su = x(4)
!       sv = x(5)
       delta = x(2:)

       dtau = sqrt(cov(1,1)/(2*n - npar))
!       do concurrent i = 1:size(delta)
!          ddelta = sqrt(cov(i,i)/(2*n - npar))
!       end do
    end if

    if( debug_iter ) then
       block
         integer, parameter :: maxpar = 7
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64) :: q
         integer :: i,j,n
         character(len=3) :: dots
         character(len=80) :: fmt

         cor = 1
         do i = 1, npar
            do j = i+1, npar
               if( cov(i,i) > 0 .and. cov(j,j) > 0 ) then
                  q = sqrt(cov(i,i)*cov(j,j))
                  if( abs(q) > 0 ) then
                     cor(i,j) = cov(i,j) / q
                  else
                     cor(i,j) = huge(cor)
                  end if
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         n = min(maxpar,npar)
         if( n < npar ) then
            dots = '...'
         else
            dots = ''
         end if

         write(error_unit,*) 'Covariance (Hessian in minimum):'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(es9.2,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cov(i,1:n),dots
         end do

         write(error_unit,*) 'Correlation:'
         write(fmt,'(a,i0,a)') '(i2,a,',n,'(f6.3,1x),a)'
         do i = 1, n
            write(error_unit,fmt) i,': ',cor(i,1:n),dots
         end do

       end block
    end if


  end subroutine theline1


  subroutine like1(npar,mpar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, mpar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(in out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(in out) :: fjac

!    real(REAL64), dimension(4) :: vec
!    real(REAL64), dimension(4,4) :: jac
!    real(REAL64) :: tau, sigu, sigv
    real(REAL64) :: su,sv
    character(len=3) :: dots

!    n = size(delta)
!    tau = x(1)
!    sigu = x(2)
!    sigv = x(3)
!    s = x(4)
!    delta = x(2:)

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          if( size(x) > 8 ) dots = '...'
          write(error_unit,*) 's0 = ',real(sum(x(2:)**2)/(size(x)-1)),&
               real(sqrt(sum(fvec(2:)**2)/(size(fvec)-1)))
          write(error_unit,*) 'x = ',x(1:8),trim(dots)
          write(error_unit,*) 'fvec = ',fvec(1:8),trim(dots)
          block
            integer :: i,n
            character(len=80) :: fmt
            n = min(7,npar)
            write(fmt,'(a,i0,a)') '(a,i0,a,',n,'(1x,es9.2),a)'
            do i = 1, n
               write(error_unit,fmt) &
                    'fjac(',i,',:)=',fjac(i,1:n),trim(dots)
            end do
          end block
       end if

       ! convergence check
!       if( tau < tmin .or. tau > tmax ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

       sv = 1
       su = 1
!       fvec = gradf1(x(1),siguback,sigvback,su,sv,x(2:))

!       su = ydisp(x(1),x(2:))
!       c = [ su, sv ] - scale

!       fvec(1) = fvec(1) - lambda(1)*su_t + nu*c(1)*su_t
!       fvec(2:) = fvec(2:) - lambda(1)*su_d - lambda(2)*sv_d &
!            +(nu/2)*(c(1)*su_d + c(2)*sv_d)

!       fvec = gradf(tau,sigu,sigv,s)
!       fvec(4) = fvec(4) - lambda + nu*(s - scale)

    else if( iflag == 2 ) then

       sv = 1
       su = 1
!       fjac = jacf(tau,sigu,sigv,s)
 !      fjac = jacf1(x(1),siguback,sigvback,su,sv,x(2:))

!       fjac(1,1) = fjac(1,1) - lambda(1)*su_tt + nu*(su_t**2 - c(1)*su_tt)
!       do i = 1, n
!          fjac(i,i) = fjac(i,i) - lambda(1)*su_dd(i) -lambda(2)*sv_dd(i) &
!          + nu*(su_d(i)**2 - c(1)*su_dd(i)) + sv_d(i)**2 - c(2)*sv_dd(i))
!       end do
!       fjac(1,2:) = fjac(1,2:) &
!            - lambda(1)*su_td + nu*(su_t*su_d - c(1)*su_td)

    end if

  end subroutine like1

  subroutine liked1(npar,mpar,x,fvec,iflag)

    integer, intent(in) :: npar, mpar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(in out) :: fvec
 !   real(REAL64) :: su,sv
    character(len=3) :: dots

 !   su = 1
 !   sv = 1

!    real(REAL64) :: tau, sigu, sigv, s

!    tau = x(1)
!    sigu = x(2)
!    sigv = x(3)
!    s = x(4)

    if( iflag == 0 ) then
       if( debug_iter ) then
          dots = ''
          if( size(x) > 8 ) dots = '...'
          write(error_unit,*) 's0 = ',real(sum(x(2:)**2)/(size(x)-1))
          write(error_unit,*) 'x = ',real(x(1:8)),trim(dots)
          write(error_unit,*) 'fvec = ',real(fvec(1:8)),trim(dots)
       end if

       ! convergence check ?
       return
    end if

!    fvec = gradf(x(1),x(2),x(3),x(4),x(5),x(6:))
!    fvec(4) = fvec(4) - lambda(1) + nu*(x(4) - scale)
!    fvec(5) = fvec(5) - lambda(2) + nu*(x(5) - scale)
!    fvec = gradf1(x(1),backrho,backs,x(2:))
    fvec = gradf1a(x(1),backsigv,backsigu,backsx,backsy,x(2:))
!    fvec(4) = fvec(4) - sum(lambda) + nu*(s - scale)

  end subroutine liked1


  function gradf1a(tau,sigv,sigu,sx,sy,delta) result(vec)

    use rfun

    real(REAL64), intent(in) :: tau,sigv,sigu,sx,sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(1 + size(delta)) :: vec
!    real(REAL64), dimension(size(delta)) :: dd
    real(REAL64) :: mu, r,d,sigx,sigy,psix,psiy, r_t, r_d, &
         sigx_t,sigy_t,sigy_d, sigx_d, dt, dd, rho, phi
    integer :: i, n

    n = size(pht)

!    sx = backsx
!    sy = backsy
!    sigu2 = sigu**2
!    sigv2 = sigv**2

    rho = sqrt(sigu**2 + (tau*sigv)**2)
    phi = atan2(sigu,tau*sigv)

    vec = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = (cts(i) - mu*tau)*(1 + tau)
       sigx = sqrt(mu*(1 + 1/tau) + (rho/tau)**2)
       sigy = sqrt(mu*tau*(1+tau) + rho**2)

       r_t = -mu*(1 + tau) + (cts(i) - mu*tau)
       r_d = -tau*(1 + tau)

       r = cts(i) - mu*tau
       sigx = sqrt(mu + sigv**2)
       sigy = sqrt(mu*tau + sigu**2)

       sigx_t = -(mu/sigx + (2*rho/tau)*(rho/sigx))/(2*tau**2)
       sigy_t = (2*tau + 1)/2 * (mu/sigy)
       sigx_d = (1 + 1/tau)/2 / sigx
       sigy_d = tau*(1 + tau)/2/sigy

       if( robust ) then
          psix = tukey((d / sigx)/sx)
          psiy = tukey((r / sigy)/sy)
       else
          psix = (d / sigx) / sx
          psiy = (r / sigy) / sy
       end if

       dt = -(psix/sx)*(d/sigx)*(sigx_t/sigx) + sigx_t/sigx &
            +(psiy/sy)*(r_t/sigy - (r/sigy)*(sigy_t/sigy)) + sigy_t/sigy

       dd = (psix/sx)*(1/sigx - (d/sigx)*(sigx_d/sigx)) + sigx_d/sigx &
            +(psiy/sy)*(r_d/sigy - (r/sigy)*(sigy_d/sigy)) + sigy_d/sigy

       vec(1) = vec(1) + dt
       vec(1+i) = dd


!       dt = -(psiy/su)*(mu/sigy + (r/sigy)*(sigy_t/sigy)) + sigy_t/sigy
!       du = sigy_u/sigy*(1 - (psiy/su)*(r/sigy))
!       dv = sigx_v/sigx*(1 - (psix/sv)*(delta(i)/sigx))
!       ds = -psix*(d/sigx)/sv**2 - psiy*(r/sigy)/su**2 !+ 2/s

!       sigx_d = 1/(2*sigx)
!       sigy_d = tau / sigy / 2
!       dd =  (psix/sv)*(1/sigx - (d/sigx)*(sigx_d/sigx)) &
!            -(psiy/su)*(tau/sigy + (r/sigy)*(sigy_d/sigy)) &
!            + sigx_d/sigx + sigy_d/sigy

!       vec(1:4) = vec(1:4) + [ dt, du, dv, ds ]
!       vec(4+i) = dd
!       vec(1) = vec(1) + dt
!       vec(1+i) = dd
    end do
    vec = vec / n
!    vec(1:4) = vec(1:4) / n
!    vec(4) =  vec(4) + 2/s

  end function gradf1a


  function gradf1(tau,rho,s,delta) result(vec)

    use rfun

    real(REAL64), intent(in) :: tau,rho,s
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(1 + size(delta)) :: vec
!    real(REAL64), dimension(size(delta)) :: dd
    real(REAL64) :: mu, r,d,sigx,sigy,psix,psiy, r_t, r_d, &
         sigx_t,sigy_t,sigy_d, sigx_d, dt, dd,sx,sy
    integer :: i, n

    n = size(pht)

!    sx = backsx
!    sy = backsy
!    sigu2 = sigu**2
!    sigv2 = sigv**2

    vec = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = (cts(i) - mu*tau)*(1 + tau)
       sigx = sqrt(mu*(1 + 1/tau) + (rho/tau)**2)
       sigy = sqrt(mu*tau*(1+tau) + rho**2)

       r_t = -mu*(1 + tau) + (cts(i) - mu*tau)
       r_d = -tau*(1 + tau)

       sigx_t = -(mu/sigx + (2*rho/tau)*(rho/sigx))/(2*tau**2)
       sigy_t = (2*tau + 1)/2 * (mu/sigy)
       sigx_d = (1 + 1/tau)/2 / sigx
       sigy_d = tau*(1 + tau)/2/sigy

       if( robust ) then
          psix = tukey((d / sigx)/sx)
          psiy = tukey((r / sigy)/sy)
       else
          psix = (d / sigx) / sx
          psiy = (r / sigy) / sy
       end if

       dt = -(psix/sx)*(d/sigx)*(sigx_t/sigx) + sigx_t/sigx &
            +(psiy/sy)*(r_t/sigy - (r/sigy)*(sigy_t/sigy)) + sigy_t/sigy

       dd = (psix/sx)*(1/sigx - (d/sigx)*(sigx_d/sigx)) + sigx_d/sigx &
            +(psiy/sy)*(r_d/sigy - (r/sigy)*(sigy_d/sigy)) + sigy_d/sigy

       vec(1) = vec(1) + dt
       vec(1+i) = dd


!       dt = -(psiy/su)*(mu/sigy + (r/sigy)*(sigy_t/sigy)) + sigy_t/sigy
!       du = sigy_u/sigy*(1 - (psiy/su)*(r/sigy))
!       dv = sigx_v/sigx*(1 - (psix/sv)*(delta(i)/sigx))
!       ds = -psix*(d/sigx)/sv**2 - psiy*(r/sigy)/su**2 !+ 2/s

!       sigx_d = 1/(2*sigx)
!       sigy_d = tau / sigy / 2
!       dd =  (psix/sv)*(1/sigx - (d/sigx)*(sigx_d/sigx)) &
!            -(psiy/su)*(tau/sigy + (r/sigy)*(sigy_d/sigy)) &
!            + sigx_d/sigx + sigy_d/sigy

!       vec(1:4) = vec(1:4) + [ dt, du, dv, ds ]
!       vec(4+i) = dd
!       vec(1) = vec(1) + dt
!       vec(1+i) = dd
    end do
    vec = vec / n
!    vec(1:4) = vec(1:4) / n
!    vec(4) =  vec(4) + 2/s

  end function gradf1


  function jacf1(tau,sigu,sigv,su,sv,delta) result(jac)

    use rfun

    real(REAL64), intent(in) :: tau,sigu,sigv,su,sv
    real(REAL64), dimension(:), intent(in) :: delta

    real(REAL64), parameter :: z = real(0.0,REAL64)
    real(REAL64), dimension(1+size(delta),1+size(delta)) :: jac
    real(REAL64) :: sigu2, sigv2, s, d,r,mu,sigx,sigy,psix,psiy,dpsix,dpsiy, &
         sigy_t,sigy_tt,dtt,&
         sigx_dd,sigy_dd,sigy_td,ddd,dtd,sigx_d,sigy_d
    integer :: i, n

    n = size(mju)
    s = scale

!    s = sv
    sigu2 = sigu**2
    sigv2 = sigv**2

    jac = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - mu*tau
       sigx = sqrt(mu + sigv2)
       sigy = sqrt(mu*tau + sigu2)
       sigy_t = mu / (2*sigy)
!       sigy_u = sigu / sigy
!       sigx_v = sigv / sigx
       sigx_d = 1/(2*sigx)
       sigy_d = tau / sigy / 2

       sigy_tt = - (mu / (2*sigy)) * (sigy_t/sigy)
!       sigy_tt = - (mu/sigy)**2/sigy/4
       sigx_dd = -(sigx_d/sigx)/(2*sigx)
       sigy_dd = -tau*(sigy_d/sigy)/sigy
       sigy_td = 1/(2*sigy) - (mu/sigy)*(sigy_d/sigy)/2

       if( robust ) then
          psix = tukey(d/(s*sigx))
          psiy = tukey(r/(s*sigy))
          dpsix = dtukey(d/(s*sigx))
          dpsiy = dtukey(r/(s*sigy))
       else
          psix = d / (s*sigx)
          psiy = r / (s*sigy)
          dpsix = 1
          dpsiy = 1
       end if

       dtt = dpsiy*(mu/sigy + (r/sigy)*(sigy_t/sigy))**2 &
            +psiy*(2*(mu/sigy)*(sigy_t/sigy) &
            - (r/sigy)*(sigy_tt/sigy)+ 2*(r/sigy)*(sigy_t/sigy)**2) &
            - (sigy_t/sigy)**2 + sigy_tt/sigy

       ddd = (dpsix/s**2)*(1/sigx - (d/sigx)*(sigx_d/sigx))**2 &
            +(dpsiy/s**2)*(tau/sigy + (r/sigy)*(sigy_d/sigy))**2 &
            -(psix/s)*(-2*(sigx_d/sigx)/sigx + (d/sigx)*(sigx_dd/sigx) &
            + 2*(d/sigx)*(sigx_d/sigx)**2) &
            -(psiy/s)*(-2*tau*(sigy_d/sigy)/sigy + (r/sigy)*(sigy_dd/sigy) &
            - 2*(r/sigy)*(sigy_d/sigy)**2) &
            -(sigx_d/sigx)**2 + sigx_dd/sigx -(sigy_d/sigy)**2 + sigy_dd/sigy

       dtd = (dpsiy/s**2)*(mu/sigy + (r/sigy)*(sigy_t/sigy)) &
            *(tau/sigy + (r/sigy)*(sigy_d/sigy)) &
            -(psiy/s)*(1/sigy - (mu/sigy)*(sigy_d/sigy) &
            -tau*(sigy_t/sigy)/sigy + (r/sigy)*(sigy_td/sigy) &
            -2*(r/sigy)*(sigy_t/sigy)*(sigy_d/sigy)) &
            -(sigy_d/sigy)*(sigy_t/sigy) + sigy_td/sigy

       jac(1,1) = jac(1,1) + dtt
       jac(1,1+i) = dtd
       jac(1+i,1) = dtd
       jac(i+1,i+1) = ddd

    end do
    jac = jac / n

  end function jacf1


  subroutine lineftest(tau,rho,phi,sx,sy,delta)

    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in) :: tau,rho,phi,sx,sy

    real(REAL64), dimension(3+size(delta)) :: fd,fvec,eps,z
    real(REAL64), dimension(3+size(delta),3+size(delta)) :: fjac, djac
    real(REAL64), dimension(-1:1,3+size(delta)) :: dg
    real(REAL64), dimension(-1:1) :: f
    real(REAL64), dimension(2) :: c
    integer :: i, j, l, iflag

    write(error_unit,*) 'lineftest():'

    eps = epsilon(eps)**(1.0/3.0)
    fvec = fgradf(tau,rho,phi,delta,sx,sy)
    z(1:3) = [ tau, rho, phi ]
    z(4:) = delta
    iflag = 1
    backsx = sx
    backsy = sy
    call fliked(size(fvec),size(fvec),z,fvec,iflag)
    do i = 1, size(fvec)
       do l = -1,1,2
          z(1:3) = [ tau, rho, phi ]
          z(4:) = delta
          z(i) = z(i) + l*eps(i)
!          c = z(4:5) - scale
          f(l) = floglike(z,sx,sy) !- sum(lambda*c) + (nu/2)*sum(c**2)
       end do
       fd(i) = (f(1) - f(-1)) / (2*eps(i))
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
            real(fd(i)-fvec(i)),real(eps(i))
    end do
    return

!!$    fjac = jacf(tau,sigu,sigv,su,sv,delta)
!!$    do i = 1,size(fjac,1)
!!$       do j = i, size(fjac,2)
!!$          do l = -1,1,2
!!$             z = [ tau, sigu, sigv, su, sv ]
!!$!             z(1) = tau
!!$!             z(2:) = delta
!!$             z(i) = z(i) + l*eps(i)
!!$             dg(l,:) = gradf(z(1),z(2),z(3),z(4),z(5),z(6:))
!!$          end do
!!$          djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps(i))
!!$
!!$          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
!!$               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
!!$       end do
!!$    end do

  end subroutine lineftest

  real(REAL64) function floglike(z,sx,sy)

    use rfun

    real(REAL64), dimension(:), intent(in) :: z
    real(REAL64), intent(in) :: sx, sy

    real(REAL64), dimension(size(z)-3) :: delta
    real(REAL64) :: tau,sigu2, sigv2,d,r,mu,sigx,sigy,rhox,rhoy,sigu,sigv, &
         su,sv,rho,phi
    integer :: i, n

    tau = z(1)
    rho = z(2)
    phi = z(3)
!    su = z(4)
!    sv = z(5)
    delta = z(4:)

    n = size(pht)
    sigu = rho*sin(phi)
    sigv = rho*cos(phi)/tau
    sigu2 = sigu**2
    sigv2 = sigv**2

    floglike = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
       sigx = sqrt(mu  + sigv2) /2
       !       sigy = sqrt(tau*mu + sigu2)
       sigy = sqrt(mu*tau*(tau+1) + rho**2) / 2

       if( robust ) then
          rhox = itukey(d/sigx/sx)
          rhoy = itukey(r/sigy/sy)
       else
          rhox = (d/sigx/sx)**2/2
          rhoy = (r/sigy/sy)**2/2
       end if

       floglike = floglike + rhox + rhoy + log(sigx) + log(sigy)

    end do
    floglike = floglike / n + log(sx) + log(sy)

  end function floglike

  subroutine linertest(tau,rho)

    use minpacks

    real(REAL64), intent(in) :: tau,rho

    real(REAL64), dimension(2) :: fd,fvec,z
    real(REAL64), dimension(2,2) :: fjac, djac, qinv
    real(REAL64), dimension(-1:1,2) :: fv
    real(REAL64), dimension(-1:1) :: f
    real(REAL64), dimension(2) :: c, cd
    real(REAL64), dimension(2,-1:1) :: cc
    real(REAL64), dimension(2,2,-1:1) :: dcc
    real(REAL64), dimension(2,2) :: dcd
    real(REAL64), dimension(:,:,:), allocatable :: dd
    real(REAL64), dimension(:,:), allocatable :: dc, da
    real(REAL64) :: h
    integer :: i, j, l, iflag

    write(error_unit,*) 'linertest():'

    h = epsilon(h)**(1.0/3.0)
    mju = (cts/tau + pht) / 2
!    fvec = gradlike(tau,rho,delta,sx,sy)
!    fvec = linegrad2(tau,rho,delta)
    fvec = gradrline(tau,rho)
    z(1) = tau
    z(2) = rho
    iflag = 1
!    call liked(size(fvec),size(fvec),z,fvec,iflag)
    do i = 1, size(fvec)
       do l = -1,1,2
          z(1) = tau
          z(2) = rho
          z(i) = z(i) + l*h
          f(l) = rloglike(z(1),z(2),mju)
       end do
       fd(i) = (f(1) - f(-1)) / (2*h)
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),real(fd(i)-fvec(i))
    end do
!    stop

    ! second derivatives
    write(error_unit,*) 'linetest() -- Jacobian:'
    !  fjac = jacf(tau,rho,delta,sx,sy)
    fjac = linerjac(tau,rho)
!    iflag = 2
!    call like(size(fjac,1),size(fjac,2),z,fvec,fjac,size(fjac,1),iflag)
    do i = 1, size(fjac,1)
       do l = -1,1,2
          z(1) = tau
          z(2) = rho
          z(i) = z(i) + l*h
          fv(l,:) = gradrline(z(1),z(2))
       end do
       djac(i,:) = (fv(1,:) - fv(-1,:)) / (2*h)

       do j = i, size(fjac,2)
          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do

  end subroutine linertest

  real(REAL64) function rloglike(tau,rho,mju)

    use rfun

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: mju

    real(REAL64) :: d,e,s,ds,mu,sigx,sigy,rhox,rhoy,zx, zy, dpsi
    integer :: i, n

    n = size(pht)

    s = 0
    ds = 0
    do i = 1, n

!       d = delta(i)
!       mu = pht(i) + d
       !       mu = (cts(i)/tau + pht(i)) / 2
       mu = mju(i)
       d = mu - pht(i)
!       d = pht(i) - mu
       e = cts(i) - tau*mu
!       d = (tau*pht(i) - cts(i)) / 2 / tau
!       e = (tau*pht(i) - cts(i)) / 2
!       e = -e

       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)

       zx = d / sigx
       zy = e / sigy
       if( robust ) then
!          rhox = ihuber(zx)
!          rhoy = ihuber(zy)
          rhox = huber(zx)**2
          rhoy = huber(zy)**2
          dpsi = dhuber(zx) + dhuber(zy)
       else
          rhox = zx**2
          rhoy = zy**2
          dpsi = 2
       end if

       s = s + rhox + rhoy !+ log(sigx) + log(sigy)
       ds = ds + dpsi

    end do
    rloglike = s / (2*ds)
!    rloglike = s / (2*n)

  end function rloglike



  subroutine linetest3(tau,rho,delta)

    use minpacks

    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in) :: tau,rho

    real(REAL64), dimension(size(delta)+2) :: fd,fvec,z
    real(REAL64), dimension(size(delta)+2,size(delta)+2) :: fjac, djac, qinv
    real(REAL64), dimension(-1:1,size(delta)+2) :: fv
    real(REAL64), dimension(-1:1) :: f
    real(REAL64), dimension(2) :: c, cd
    real(REAL64), dimension(2,-1:1) :: cc
    real(REAL64), dimension(2,size(delta)+2,-1:1) :: dcc
    real(REAL64), dimension(2,size(delta)+2) :: dcd
    real(REAL64), dimension(:,:,:), allocatable :: dd
    real(REAL64), dimension(:,:), allocatable :: dc, da
    real(REAL64) :: h
    integer :: i, j, l, iflag

    write(error_unit,*) 'linetest():'

    h = epsilon(h)**(1.0/3.0)
!    fvec = gradlike(tau,rho,delta,sx,sy)
!    fvec = linegrad2(tau,rho,delta)
    fvec = gradline(tau,rho,delta)
    z(1) = tau
    z(2) = rho
    z(3:) = delta
!    eps = eps*abs(loglike(z,sx,sy))
    iflag = 1
!    backsx = sx
!    backsy = sy
!    call liked(size(fvec),size(fvec),z,fvec,iflag)
    do i = 1, size(fvec)
       do l = -1,1,2
          z(1) = tau
          z(2) = rho
          z(3:) = delta
          z(i) = z(i) + l*h
!          c = z(3:4) - scale
          !          f(l) = loglike(z) - sum(lambda*c) + (nu/2)*sum(c**2)
 !         c = linedisp(z(1),z(2),z(3:),sx,sy)
          f(l) = loglike(z(1),z(2),z(3:)) !- sum(lambda*c) + (nu/2)*sum(c**2)
       end do
       fd(i) = (f(1) - f(-1)) / (2*h)
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),real(fd(i)-fvec(i))
    end do
!    stop

    ! test for derivatives of linedisp
!!$    write(error_unit,*) 'linetest() -- linedisp:'
!!$    da = linegradc(tau,rho,delta,sx,sy)
!!$    do i = 1, size(delta)+2
!!$       do l = -1,1,2
!!$          z(1) = tau
!!$          z(2) = rho
!!$          z(3:) = delta
!!$          z(i) = z(i) + l*eps(i)
!!$          cc(:,l) = linedisp(z(1),z(2),z(3:),sx,sy)
!!$       end do
!!$       cd = (cc(:,1) - cc(:,-1)) / (2*eps(i))
!!$       do l = 1,2
!!$          write(error_unit,*) i,l,real(da(l,i)),real(cd(l)),real(cd(l)-da(l,i))
!!$       end do
!!$    end do

!!$!    write(*,*) eps(1)
!!$    write(error_unit,*) 'linetest() -- jacc:'
!!$!    dd = linejacc(tau,rho,delta,sx,sy)
!!$    do i = 1, size(dd,2)
!!$       do j = i, size(dd,3)
!!$          do l = -1,1,2
!!$             z(1) = tau
!!$             z(2) = rho
!!$             z(3:) = delta
!!$             z(i) = z(i) + l*h
!!$ !            da = linegradc(z(1),z(2),z(3:),sx,sy)
!!$             dcc(:,:,l) = da
!!$!             write(*,*) real(z),'>',real(da)
!!$          end do
!!$!          write(*,*) real(dcc(1,:,1))
!!$!          write(*,*) real(dcc(1,:,-1))
!!$!          stop
!!$!          da = linegradc(tau+eps(i),rho,delta,sx,sy)
!!$!          dcc(:,:,1) = da
!!$!          da = linegradc(tau-eps(i),rho,delta,sx,sy)
!!$!          dcc(:,:,-1) = da
!!$          dcd(:,:) = (dcc(:,:,1) - dcc(:,:,-1)) / (2*h)
!!$!          write(*,*) real(dcd(1,:))
!!$!          stop
!!$          do l = 1,2
!!$             write(error_unit,'(3i5,3(2x,g0))') i,j,l,real(dd(l,i,j)), &
!!$                  real(dcd(l,j)), real(dd(l,i,j)-dcd(l,j))
!!$          end do
!!$       end do
!!$    end do
!    return


    ! second derivatives
    write(error_unit,*) 'linetest() -- Jacobian:'
    !  fjac = jacf(tau,rho,delta,sx,sy)
    fjac = linejac(tau,rho,delta)
!    iflag = 2
!    call like(size(fjac,1),size(fjac,2),z,fvec,fjac,size(fjac,1),iflag)
    do i = 1, size(fjac,1)
       do l = -1,1,2
          z(1) = tau
          z(2) = rho
          z(3:) = delta
          z(i) = z(i) + l*h
          fv(l,:) = gradline(z(1),z(2),z(3:))
!             c = z(3:4) - scale
!             f(l) = loglike(z) - lambda(1)*c(1) + (nu/2)*c(1)**2
!             c = linedisp(z(1),z(2),z(3:),sx,sy)
!             dc = linegradc(z(1),z(2),z(3:),sx,sy)
!             dg(l,:) = gradf(z(1),z(2),z(3:),sx,sy) &
!                  + dc(1,:)*(nu*c(1)-lambda(1)) + dc(2,:)*(nu*c(2)-lambda(2))
!             call liked(size(fvec),size(fvec),z,dg(l,:),iflag)
       end do
       djac(i,:) = (fv(1,:) - fv(-1,:)) / (2*h)

       do j = i, size(fjac,2)
          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do




!    call qrinv(djac,qinv)
!    z(1) = 1e-3
!    z(2) = 0.1*rho
!    z(3) = 0.1
!    z(4:) = sqrt((pht+delta)*(1+1/tau) + (rho/tau)**2)
!    do i = 1, size(qinv,1)
!       write(*,*) i,real(qinv(i,i)),real(z(i)**2/qinv(i,i))
!    end do

  end subroutine linetest3


  subroutine linetest(tau,rho,lambda,mju)

    use minpacks

    real(REAL64), intent(in) :: tau,rho,lambda
    real(REAL64), dimension(:), intent(in) :: mju

    real(REAL64), dimension(2) :: fd,fvec,z,cvec
    real(REAL64), dimension(2,2) :: fjac, djac, qinv, jjac, cjac
    real(REAL64), dimension(2,-1:1) :: fv
    real(REAL64), dimension(-1:1) :: f
!    real(REAL64), dimension(2) :: c, cd
    real(REAL64), dimension(2,-1:1) :: cc
    real(REAL64), dimension(2,2,-1:1) :: dcc
    real(REAL64), dimension(2,2) :: dcd
    real(REAL64), dimension(:,:,:), allocatable :: dd
    real(REAL64), dimension(:,:), allocatable :: dc, da
    real(REAL64) :: h,c
    integer :: i, j, l, iflag

    write(error_unit,*) 'linetest():'

    h = epsilon(h)**(1.0/3.0)
    !    fvec = grad_line(tau,rho,lambda,mju)
    iflag = 1
    z = [ tau, rho ]
    !call thelike(2,2,z,fvec,fjac,2,iflag)
    call slike(z,f(1),fvec,fjac,c,cvec,cjac,iflag)
    fvec = fvec - lambda*cvec
    do i = 1, size(fvec)
       do l = -1,1,2
          z = [ tau, rho ]
          z(i) = z(i) + l*h
!          c = cvar(z(1),z(2)) !/ ndata  !logvar(z(1),z(2)) - scale
!          f(l) = loglike(z(1),z(2),mju) - lambda*c + (nu/2)*c**2
          f(l) = loglike(z(1),z(2),mju) - lambda*cvar(z(1),z(2))
       end do
       fd(i) = (f(1) - f(-1)) / (2*h)
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),real(fd(i)-fvec(i))
    end do

    write(error_unit,*) 'linetest() -- constraints:'
    cvec = cvar_sig(tau,rho) !/ ndata
    do i = 1, size(fvec)
       do l = -1,1,2
          z = [ tau, rho ]
          z(i) = z(i) + l*h
          f(l) = cvar(z(1),z(2)) !/ ndata !logvar(z(1),z(2)) - scale
       end do
       fd(i) = (f(1) - f(-1)) / (2*h)
       write(error_unit,*) i,real(cvec(i)),real(fd(i)),real(fd(i)-cvec(i))
    end do



    ! second derivatives
    write(error_unit,*) 'linetest() -- Jacobian:'
    !    fjac = jac_line(tau,rho,lambda,mju)
    iflag = 2
    z = [ tau, rho ]
    call thelike(2,2,z,fvec,fjac,2,iflag)
    do i = 1, size(fjac,1)
       do l = -1,1,2
          z = [ tau, rho ]
          z(i) = z(i) + l*h
          iflag = 1
          call thelike(2,2,z,fv(:,l),jjac,2,iflag)
!          fv(l,:) = grad_line(z(1),z(2),lambda,mju)
       end do
       djac(i,:) = (fv(:,1) - fv(:,-1)) / (2*h)

       do j = i, size(fjac,2)
          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real(fjac(i,j)-djac(i,j))
       end do
    end do

  end subroutine linetest


  real(REAL64) function loglike(tau,rho,mju)

    use rfun

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: mju

    real(REAL64) :: d,e,s,ds,sigx,sigy,rhox,rhoy,zx, zy, dpsi
    integer :: i

    s = 0
    ds = 0
    do i = 1, ndata

       d = mju(i) - pht(i)
       e = cts(i) - tau*mju(i)

       sigx = sigmax(tau,rho,pht(i))
       sigy = sigmay(tau,rho,cts(i))

       zx = d / sigx
       zy = e / sigy
       if( robust ) then
!          rhox = huber(zx)**2/2
!          rhoy = huber(zy)**2/2
          rhox = ihuber(zx)
          rhoy = ihuber(zy)
          dpsi = dhuber(zx) + dhuber(zy)
!          dpsi = 2
       else
          rhox = zx**2/2
          rhoy = zy**2/2
          dpsi = 2
       end if

       s = s + rhox + rhoy + log(sigx) + log(sigy)
       ds = ds + dpsi

    end do
!    loglike = s / ds
    loglike = s !/ ndata

  end function loglike


  function cvar_sig(tau,rho) result(dc)

    use rfun

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(2) :: dc!, s2grad
    real(REAL64), dimension(:), allocatable :: sigx, sigy, zx, zy, psix, psiy, &
         dpsix, dpsiy, sigx_t, sigx_r, sigy_t, sigy_r
    real(REAL64) :: ds
    integer :: i

    allocate(sigx(ndata),sigy(ndata),zx(ndata),zy(ndata),psix(ndata),psiy(ndata), &
         dpsix(ndata),dpsiy(ndata),sigx_t(ndata),sigx_r(ndata),sigy_t(ndata), &
         sigy_r(ndata))

    sigx = sigmax(tau,rho,pht)
    sigy = sigmay(tau,rho,cts)
    zx = (mju - pht) / sigx
    zy = (cts - tau*mju) / sigy
    if( robust ) then
       psix = huber(zx)
       psiy = huber(zy)
       dpsix = dhuber(zx)
       dpsiy = dhuber(zy)
    else
       psix = zx
       psiy = zy
       dpsix = 1
       dpsiy = 1
    end if

    do i = 1, ndata
       call grad_sigx(sigx(i),pht(i),tau,rho,sigx_t(i),sigx_r(i))
       call grad_sigy(sigy(i),cts(i),tau,rho,sigy_t(i),sigy_r(i))
    end do
    ds = sum(dpsix + dpsiy)

    !    logvar = sigvar(psix,psiy,dpsix,dpsiy)
    !c = sigvar(psix,psiy,dpsix,dpsiy) - scale
!    s2grad = grad_ressum(mju,sigx,sigy,zx,zy,psix,psiy,dpsix,dpsiy,sigx_t,sigx_r,sigy_t,sigy_r)
    dc = grad_sigvar(mju,zx,zy,psix,psiy,dpsix,dpsiy,sigx,sigy,sigx_t,sigx_r,sigy_t,sigy_r)! / ds

  end function cvar_sig


  real(REAL64) function loglike3(tau,rho,delta)

    use rfun

    real(REAL64), parameter :: a = real(1.345,REAL64)
    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: delta
!    real(REAL64), intent(in) :: sx, sy

!    real(REAL64), dimension(size(z)-2) :: delta
    real(REAL64) :: d,e,s,mu,sigx,sigy,rhox,rhoy,rho2, zx, zy, s2
    integer :: i, n

!    tau = z(1)
!    rho = z(2)
!    phi = z(3)
!    sx = z(3)
!    sy = sx
!    sy = z(4)
!    delta = z(3:)

    n = size(pht)
    rho2 = rho**2
!    sigu = rho*sin(phi)
!    sigv = rho*cos(phi)/tau
!    sigu2 = sigu**2
!    sigv2 = sigv**2

    s = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       e = cts(i) - tau*mu
!       write(*,*) '.',i,d,mu,rho,tau
!       sigx = sqrt(mu*(1 + 1/tau)  + (rho/tau)**2)/2
!       sigy = sqrt(mu*tau*(1 + tau) + rho**2)/2
       !       sigx = sqrt(mu  + sigv2) / 2
!       s2 = mu*tau*(tau+1) + rho2
!       if( .not. (s2 > 0) ) exit
!       sigy = sqrt(s2) / 2
!       sigx = sigy / tau

!       sigx = sqrt(pht(i)*(1 + 1/tau) + (rho/tau)**2) / 2
!       sigy = sqrt(cts(i)*(1 + tau) + rho**2) / 2
       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)


       zx = d / sigx
       zy = e / sigy
       if( robust ) then
!          rhox = itukey(zx)
!          rhoy = itukey(zy)

          rhox = max(-a,min(zx**2/2,a))
          rhoy = max(-a,min(zy**2/2,a))

          rhox = ihuber(zx)
          rhoy = ihuber(zy)

!          rhox = max(min(rx,rmax),rmin)**2/2
!          rhoy = max(min(ry,rmax),rmin)**2/2
       else
          rhox = zx**2/2
          rhoy = zy**2/2
       end if

!       rx = max(zmin,min(rx,zmax))
!       ry = max(zmin,min(ry,zmax))

!       rhox = zx**2/2
!       rhoy = zy**2/2
       s = s + rhox + rhoy !+ log(sigx) + log(sigy)

    end do
!    if( sx > 0 .and. sy > 0 ) then
    loglike3 = s / (2*n) !+ log(sx) + log(sy)
!    else
!       loglike = 666
!    end if

  end function loglike3


  real(REAL64) function cloglike(tau,rho,delta)

    use rfun

    real(REAL64) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: delta

!    real(REAL64), dimension(size(z)-1) :: delta
!    real(REAL64), dimension(:), pointer :: delta
    real(REAL64) :: d,r,s,mu,sigx,sigy,rhox,rhoy,rho2, rx, ry
    integer :: i, n

!    tau = z(1)
!    rho = z(2)
    !    delta = z(3:)
!    delta => z(3:)

    n = size(delta)
    rho2 = rho**2

    s = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
!       write(*,*) '.',i,d,mu,rho,tau
!       sigx = sqrt(mu*(1 + 1/tau)  + (rho/tau)**2)/2
!       sigy = sqrt(mu*tau*(1 + tau) + rho**2)/2
       !       sigx = sqrt(mu  + sigv2) / 2
!       s2 = mu*tau*(tau+1) + rho2
!       if( .not. (s2 > 0) ) exit
!       sigy = sqrt(s2) / 2
       !       sigx = sigy / tau

       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)
       rx = d / sigx
       ry = r / sigy
       if( robust ) then
          rhox = itukey(rx)
          rhoy = itukey(ry)
!          rhox = max(min(rx,rmax),rmin)**2/2
!          rhoy = max(min(ry,rmax),rmin)**2/2
       else
          rhox = rx**2/2
          rhoy = ry**2/2
       end if

!       write(*,*) rx,ry,sigx,sigy
       s = s + rhox + rhoy + log(sigx) + log(sigy)

    end do
!    write(*,*) s/n
!    stop
    cloglike = s / n

  end function cloglike



  subroutine sigtest(tau,rho)

!    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in) :: tau,rho

    real(REAL64), dimension(2) :: z
    real(REAL64), dimension(2,2) :: djacx, djacy
    real(REAL64), dimension(-1:1,2) :: dx, dy
    real(REAL64), dimension(-1:1) :: sigx, sigy
    real(REAL64), dimension(2) :: dsigx, dsigy
    real(REAL64) :: h,sigx_t,sigx_r,sigy_t,sigy_r, &
         sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr
    integer :: i, l, k

    write(error_unit,*) 'sigtest():'

    h = epsilon(h)**(1.0/3.0)
    do k = 1, 1
       call sigmas(tau,rho,pht(k),cts(k),sigx(0),sigy(0))
       call grad_sigx(sigx(0),pht(k),tau,rho,sigx_t,sigx_r)
       call grad_sigy(sigy(0),cts(k),tau,rho,sigy_t,sigy_r)
       do i = 1, 2
          do l = -1,1,2
             z(1:2) = [ tau, rho ]
             z(i) = z(i) + l*h
             call sigmas(z(1),z(2),pht(k),cts(k),sigx(l),sigy(l))
          end do
          dsigx(i) = (sigx(1) - sigx(-1)) / (2*h)
          dsigy(i) = (sigy(1) - sigy(-1)) / (2*h)
       end do
       write(error_unit,*) k,'x_t',real(sigx_t),real(dsigx(1)),real(sigx_t-dsigx(1))
       write(error_unit,*) k,'x_r',real(sigx_r),real(dsigx(2)),real(sigx_r-dsigx(2))
       write(error_unit,*) k,'y_t',real(sigy_t),real(dsigy(1)),real(sigy_t-dsigy(1))
       write(error_unit,*) k,'y_r',real(sigy_r),real(dsigy(2)),real(sigy_r-dsigy(2))
    end do

    do k = 1, 1
       call sigmas(tau,rho,pht(k),cts(k),sigx(0),sigy(0))
       call grad_sigx(sigx(0),pht(k),tau,rho,sigx_t,sigx_r)
       call grad_sigy(sigy(0),cts(k),tau,rho,sigy_t,sigy_r)
       call jac_sigx(sigx(0),pht(k),tau,rho,sigx_t,sigx_r,sigx_tt,sigx_tr,sigx_rr)
       call jac_sigy(sigy(0),cts(k),tau,rho,sigy_t,sigy_r,sigy_tt,sigy_tr,sigy_rr)
       do i = 1,2
          do l = -1,1,2
             z(1:2) = [ tau, rho ]
             z(i) = z(i) + l*h
             call sigmas(z(1),z(2),pht(k),cts(k),sigx(l),sigy(l))
             call grad_sigx(sigx(l),pht(k),z(1),z(2),dx(l,1),dx(l,2))
             call grad_sigy(sigy(l),cts(k),z(1),z(2),dy(l,1),dy(l,2))
          end do
          djacx(i,:) = (dx(1,:) - dx(-1,:)) / (2*h)
          djacy(i,:) = (dy(1,:) - dy(-1,:)) / (2*h)
       end do

       write(error_unit,*) k,'x_tt',real(sigx_tt),real(djacx(1,1)),&
            real(sigx_tt - djacx(1,1))
       write(error_unit,*) k,'x_tr',real(sigx_tr),real(djacx(1,2)),&
            real(sigx_tr - djacx(1,2))
       write(error_unit,*) k,'x_rr',real(sigx_rr),real(djacx(2,2)),&
            real(sigx_rr - djacx(2,2))
       write(error_unit,*) k,'y_tt',real(sigy_tt),real(djacy(1,1)),&
            real(sigy_tt - djacy(1,1))
       write(error_unit,*) k,'y_tr',real(sigy_tr),real(djacy(1,2)),&
            real(sigy_tr - djacy(1,2))
       write(error_unit,*) k,'y_rr',real(sigy_rr),real(djacy(2,2)),&
            real(sigy_rr - djacy(2,2))
!200    format(i3,2x,a,3(2x,g0))
    end do

  end subroutine sigtest


!!$  subroutine testgradsig(z,gradx,grady)
!!$
!!$    real(REAL64), dimension(:), intent(in) :: z
!!$    real(REAL64), dimension(:), intent(out) :: gradx, grady
!!$    real(REAL64) :: tau,rho,sigx, sigy, sigx_t,sigy_t,sigy_r,sigx_r
!!$    integer :: i, n
!!$
!!$    tau = z(1)
!!$    rho = z(2)
!!$
!!$    n = size(pht)
!!$    gradx = 0
!!$    grady = 0
!!$    do i = 1, n
!!$
!!$       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)
!!$       call grad_sigx(sigx,pht(i),tau,rho,sigx_t,sigx_r)
!!$       call grad_sigy(sigy,cts(i),tau,rho,sigy_t,sigy_r)
!!$
!!$       gradx(1:2) = gradx(1:2) + [sigx_t, sigx_r ]
!!$       grady(1:2) = grady(1:2) + [sigy_t, sigy_r ]
!!$       gradx(2+i) = 0
!!$       grady(2+i) = 0
!!$
!!$    end do
!!$
!!$    gradx = gradx / n
!!$    grady = grady / n
!!$
!!$  end subroutine testgradsig
!!$
!!$  subroutine testjacsig(z,jacx,jacy)
!!$
!!$    real(REAL64), dimension(:), intent(in) :: z
!!$    real(REAL64), dimension(:,:), intent(out) :: jacx, jacy
!!$    real(REAL64) :: tau,rho,sigx, sigy, sigx_t,sigy_t,sigy_r,sigx_r,&
!!$         sigx_tt,sigx_tr,sigx_rr,sigy_tt,sigy_tr,sigy_rr
!!$    integer :: i, j, n
!!$
!!$    tau = z(1)
!!$    rho = z(2)
!!$
!!$    n = size(pht)
!!$    jacx = 0
!!$    jacy = 0
!!$    do i = 1, n
!!$
!!$       call sigmas(tau,rho,pht(i),cts(i),sigx,sigy)
!!$       call grad_sigx(sigx,pht(i),tau,rho,sigx_t,sigx_r)
!!$       call grad_sigy(sigy,cts(i),tau,rho,sigy_t,sigy_r)
!!$
!!$       call jac_sigx(sigx,pht(i),tau,rho,sigx_t,sigx_r,sigx_tt,sigx_tr,sigx_rr)
!!$       call jac_sigy(sigy,cts(i),tau,rho,sigy_t,sigy_r,sigy_tt,sigy_tr,sigy_rr)
!!$
!!$       jacx(1,1:2) = jacx(1,1:2) + [ sigx_tt, sigx_tr ]
!!$       jacx(2,1:2) = jacx(2,1:2) + [ sigx_tr, sigx_rr ]
!!$       jacx(2+i,2+i) = 0
!!$       jacx(1,2+i) = 0
!!$       jacx(2,2+i) = 0
!!$
!!$       jacy(1,1:2) = jacy(1,1:2) + [ sigy_tt, sigy_tr ]
!!$       jacy(2,1:2) = jacy(2,1:2) + [ sigy_tr, sigy_rr ]
!!$       jacy(2+i,2+i) = 0
!!$       jacy(1,2+i) = 0
!!$       jacy(2,2+i) = 0
!!$
!!$    end do
!!$
!!$    jacx = jacx / n
!!$    jacy = jacy / n
!!$
!!$    do i = 1, size(jacx,1)
!!$       do j = 1, i - 1
!!$          jacx(i,j) = jacx(j,i)
!!$          jacy(i,j) = jacy(j,i)
!!$       end do
!!$    end do
!!$
!!$  end subroutine testjacsig
!!$
!!$  subroutine gradsig(z,fx,fy,gradx,grady,jacx,jacy)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: z
!!$    real(REAL64), intent(out) :: fx, fy
!!$    real(REAL64), dimension(:), intent(out) :: gradx, grady
!!$    real(REAL64), dimension(:,:) :: jacx, jacy
!!$    real(REAL64) :: tau,rho,sigx, sigy, mu,d,r,psix,psiy,rx,ry,r_t,r_d, &
!!$         sigx_t,sigy_t,sigy_r,sigx_r,sigy_d, sigx_d, &
!!$         sigx_tt,sigy_tt, sigx_tr,sigy_tr, sigx_td,sigy_td, &
!!$         sigx_rr,sigy_rr, sigx_rd,sigy_rd, sigx_dd,sigy_dd
!!$
!!$    integer :: i, j, n
!!$
!!$    tau = z(1)
!!$    rho = z(2)
!!$
!!$    n = size(pht)
!!$
!!$    fx = 0
!!$    fy = 0
!!$    gradx = 0
!!$    grady = 0
!!$    jacx = 0
!!$    jacy = 0
!!$    do i = 1, n
!!$
!!$       d = z(2+i)
!!$       mu = pht(i) + d
!!$       r = cts(i) - mu*tau
!!$       r_t = -mu
!!$       r_d = -tau
!!$
!!$!       sigx = sqrt(mu + sigv2) / 2
!!$!       sigy = sqrt(mu*tau + sigu2)
!!$       sigy = sqrt(mu*(tau + 1)*tau + rho**2) / 2
!!$       sigx = sigy / tau
!!$
!!$!       sigy_t = (mu / sigy) / 2
!!$!       sigx_t = -c2*(rho/sigx) * (rho/tau**3)
!!$!       sigy_r = (rho/sigy)*s2
!!$!       sigx_r = (rho/sigx)*(cosfi/tau)**2
!!$!       sigy_f = cs * (rho/sigy) * rho
!!$!       sigx_f = -cs * (rho/sigx) * (rho/tau**2)
!!$!       sigx_d = 1/(2*sigx)
!!$!       sigy_d = tau / sigy / 2
!!$
!!$ !      call grad_sigx(sigx,mu,tau,rho,cosfi,sinfi,sigx_t,sigx_r,sigx_f,sigx_d)
!!$!       call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!!$!       sigx_f = 0
!!$!       call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)
!!$!       sigy_f = 0
!!$
!!$
!!$       sigx_tt = c2/tau**3*(3/tau*(rho/sigx) + (rho/sigx)*(sigx_t/sigx))*rho
!!$       sigx_tr = c2/tau**3*((rho/sigx)**2*sigx_r - 2*(rho/sigx))
!!$       sigx_tf = (2*cs*rho + c2*(rho/sigx)*sigx_f)*(rho/sigx)/tau**3
!!$       sigx_td = c2/tau**3*(rho/sigx)**2*sigx_d
!!$       sigx_rr = c2/tau**2*(1 - (rho/sigx)*sigx_r)/sigx
!!$       sigx_rf = -(2*cs + c2*(sigx_f/sigx))*(rho/sigx)/tau**2
!!$       sigx_rd = -c2*(rho/sigx)*(sigx_d/sigx)/tau**2
!!$       sigx_ff = (-(2*c2-1)*rho + cs*(rho/sigx)*sigx_f)*(rho/sigx)/tau**2
!!$       sigx_fd = cs*sigx_d*(rho/sigx)**2/tau**2
!!$       sigx_dd = -(sigx_d/sigx)/2/sigx
!!$
!!$       sigy_tt = -(mu/sigy)*(sigy_t/sigy)/2
!!$       sigy_tr = -(mu/sigy)*(sigy_r/sigy)/2
!!$       sigy_tf = -(mu/sigy)*(sigy_f/sigy)/2
!!$       sigy_td = 1/(2*sigy) - (mu/sigy)*(sigy_d/sigy)/2
!!$       sigy_rr = s2*(1/sigy - (rho/sigy)*(sigy_r/sigy))
!!$       sigy_rf = (2*cs - s2*(sigy_f/sigy))*(rho/sigy)
!!$       sigy_rd = -s2*(rho/sigy)*(sigy_d/sigy)
!!$       sigy_ff = rho*(c2 - s2)*(rho/sigy) - cs*(rho/sigy)**2*sigy_f
!!$       sigy_fd = -cs*sigy_d*(rho/sigy)**2
!!$       sigy_dd = -(tau/2)*(sigy_d/sigy**2)
!!$
!!$!       call jac_sigx(sigx,mu,tau,rho,c2,cs,sigx_t,sigx_r,sigx_f,sigx_d, &
!!$!            sigx_tt,sigx_tr,sigx_tf,sigx_td,sigx_rr,sigx_rf,sigx_rd, &
!!$!            sigx_ff,sigx_fd,sigx_dd)
!!$!       call jac_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d, &
!!$!            sigx_tt,sigx_tr,sigx_td,sigx_rr,sigx_rd,sigx_dd)
!!$!       sigx_tf = 0
!!$!       sigx_rf = 0
!!$!       sigx_ff = 0
!!$!       sigx_fd = 0
!!$
!!$!       call jac_sigy(sigy,mu,sigy_t,sigy_r,sigy_d, &
!!$!            sigy_tt,sigy_tr,sigy_td,sigy_rr,sigy_rd,sigy_dd)
!!$!       sigy_tf = 0
!!$!       sigy_rf = 0
!!$!       sigy_ff = 0
!!$!       sigy_fd = 0
!!$
!!$       fx = fx + sigx
!!$       fy = fy + sigy
!!$
!!$       gradx(1:2) = gradx(1:2) + [sigx_t, sigx_r ]
!!$       grady(1:2) = grady(1:2) + [sigy_t, sigy_r ]
!!$       gradx(2+i) = sigx_d
!!$       grady(2+i) = sigy_d
!!$
!!$       jacx(1,1:2) = jacx(1,1:2) + [ sigx_tt, sigx_tr ]
!!$       jacx(2,1:2) = jacx(2,1:2) + [ sigx_tr, sigx_rr ]
!!$!       jacx(3,1:2) = jacx(3,1:2) + [ sigx_tf, sigx_rf ]
!!$       jacx(2+i,2+i) = sigx_dd
!!$       jacx(1,2+i) = sigx_td
!!$       jacx(2,2+i) = sigx_rd
!!$!       jacx(3,3+i) = sigx_fd
!!$
!!$       jacy(1,1:2) = jacy(1,1:2) + [ sigy_tt, sigy_tr ]
!!$       jacy(2,1:2) = jacy(2,1:2) + [ sigy_tr, sigy_rr ]
!!$!       jacy(3,1:3) = jacy(3,1:3) + [ sigy_tf, sigy_rf, sigy_ff ]
!!$       jacy(2+i,2+i) = sigy_dd
!!$       jacy(1,2+i) = sigy_td
!!$       jacy(2,2+i) = sigy_rd
!!$!       jacy(3,3+i) = sigy_fd
!!$
!!$    end do
!!$    fx = fx / n
!!$    fy = fy / n
!!$    gradx = gradx / n
!!$    grady = grady / n
!!$    jacx = jacx / n
!!$    jacy = jacy / n
!!$
!!$    do i = 1, size(jacx,1)
!!$       do j = 1, i - 1
!!$          jacx(i,j) = jacx(j,i)
!!$          jacy(i,j) = jacy(j,i)
!!$       end do
!!$    end do
!!$
!!$  end subroutine gradsig

!!$  subroutine grad_sigx(sigx,mu,tau,rho,cosfi,sinfi,sigx_t,sigx_r,sigx_f,sigx_d)
!!$
!!$    real(REAL64), intent(in) :: sigx,mu,tau,rho,cosfi,sinfi
!!$    real(REAL64), intent(out) :: sigx_t,sigx_r,sigx_f,sigx_d
!!$    real(REAL64) :: ct2, rs
!!$
!!$    ct2 = (cosfi/tau)**2
!!$    rs = rho / sigx / 4
!!$    sigx_t = -(rs/tau) * (rho*ct2)
!!$    sigx_r = rs*ct2
!!$    sigx_f = -rs* (cosfi*sinfi) * (rho/tau**2)
!!$    sigx_d = 1/(8*sigx)
!!$
!!$  end subroutine grad_sigx

!!$  subroutine grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!!$
!!$    real(REAL64), intent(in) :: sigx,mu,tau,rho
!!$    real(REAL64), intent(out) :: sigx_t,sigx_r,sigx_d
!!$
!!$    sigx_t = -(mu + 2*rho**2/tau)/sigx / (8*tau**2)
!!$    sigx_r = (rho/sigx) / (4*tau**2)
!!$    sigx_d = (1 + 1/tau)/sigx / 8
!!$
!!$  end subroutine grad_sigx
!!$
!!$
!!$  subroutine grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)
!!$
!!$    real(REAL64), intent(in) :: sigy,mu,tau,rho
!!$    real(REAL64), intent(out) :: sigy_t,sigy_r,sigy_d
!!$
!!$    real(REAL64), parameter :: half = real(1.0/2.0,REAL64)
!!$
!!$    sigy_t = (tau + half) * (mu/sigy) / 4
!!$    sigy_r = rho / sigy / 4
!!$    sigy_d = tau*(tau + 1)/8 / sigy
!!$
!!$  end subroutine grad_sigy


!!$  subroutine jac_sigx(sigx,mu,tau,rho,c2,cs,sigx_t,sigx_r,sigx_f,sigx_d, &
!!$       sigx_tt,sigx_tr,sigx_tf,sigx_td,sigx_rr,sigx_rf,sigx_rd, &
!!$       sigx_ff,sigx_fd,sigx_dd)
!!$
!!$    real(REAL64), intent(in) :: sigx,mu,tau,rho,c2,cs,sigx_t,sigx_r, &
!!$         sigx_f,sigx_d
!!$    real(REAL64), intent(out) :: sigx_tt,sigx_tr,sigx_tf,sigx_td, &
!!$         sigx_rr,sigx_rf,sigx_rd, sigx_ff,sigx_fd,sigx_dd
!!$
!!$    sigx_tt = c2/tau**3*(3/tau*(rho/sigx) + (rho/sigx)*(sigx_t/sigx))*rho / 4
!!$    sigx_tr = c2/tau**3*((rho/sigx)**2*sigx_r - 2*(rho/sigx)) / 4
!!$    sigx_tf = (2*cs*rho + c2*(rho/sigx)*sigx_f)*(rho/sigx)/tau**3 / 4
!!$    sigx_td = c2/tau**3*(rho/sigx)**2*sigx_d / 4
!!$    sigx_rr = c2/tau**2*(1 - (rho/sigx)*sigx_r)/sigx / 4
!!$    sigx_rf = -(2*cs + c2*(sigx_f/sigx))*(rho/sigx)/tau**2 / 4
!!$    sigx_rd = -c2*(rho/sigx)*(sigx_d/sigx)/tau**2 /4
!!$    sigx_ff = (-(2*c2-1)*rho + cs*(rho/sigx)*sigx_f)*(rho/sigx)/tau**2 /4
!!$    sigx_fd = cs*sigx_d*(rho/sigx)**2/tau**2 /4
!!$    sigx_dd = -(sigx_d/sigx)/2/sigx /4
!!$
!!$  end subroutine jac_sigx

!!$  subroutine jac_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d, &
!!$       sigx_tt,sigx_tr,sigx_td,sigx_rr,sigx_rd,sigx_dd)
!!$
!!$    real(REAL64), intent(in) :: sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d
!!$    real(REAL64), intent(out) :: sigx_tt,sigx_tr,sigx_td, &
!!$         sigx_rr,sigx_rd, sigx_dd
!!$    real(REAL64) :: w, t2, t3
!!$
!!$    w = mu + 2*rho**2/tau
!!$    t2 = 4*tau**2
!!$    t3 = 8*tau**3
!!$    sigx_tt = (3*w/sigx - mu/sigx + tau*(w/sigx)*(sigx_t/sigx))/t3
!!$    sigx_tr = -(4*rho/sigx - tau*(w/sigx)*(sigx_r/sigx))/t3
!!$    sigx_td = tau*(w/sigx)*(sigx_d/sigx)/t3 - tau/sigx/t3
!!$    sigx_rr = (1/sigx - (rho/sigx)*(sigx_r/sigx))/t2
!!$    sigx_rd = -(rho/sigx)*(sigx_d/sigx)/t2
!!$    sigx_dd = -((1+1/tau)/sigx)*(sigx_d/sigx)/8
!!$
!!$  end subroutine jac_sigx
!!$
!!$
!!$  subroutine jac_sigy(sigy,mu,sigy_t,sigy_r,sigy_d, &
!!$       sigy_tt,sigy_tr,sigy_td,sigy_rr,sigy_rd,sigy_dd)
!!$
!!$    real(REAL64), intent(in) :: sigy,mu,sigy_t,sigy_r,sigy_d
!!$    real(REAL64), intent(out) :: sigy_tt,sigy_tr,sigy_td,sigy_rr,sigy_rd,sigy_dd
!!$
!!$    real(REAL64), parameter :: quarter = real(1.0/4.0,REAL64)
!!$
!!$
!!$    sigy_tt = (mu/4 - sigy_t**2) / sigy
!!$    sigy_tr = - sigy_t*sigy_r / sigy
!!$    sigy_td =  (1/mu - sigy_d/sigy) * sigy_t
!!$    sigy_rr = (quarter - sigy_r**2) / sigy
!!$    sigy_rd = - sigy_r*sigy_d / sigy
!!$    sigy_dd = - sigy_d**2 / sigy
!!$
!!$  end subroutine jac_sigy


  subroutine cscaletest(tau,rho,delta,sx,sy)

    use minpacks

    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in) :: tau,rho,sx,sy

    real(REAL64), dimension(size(delta)+1) :: dc,cgrad,eps,z
    real(REAL64), dimension(-1:1) :: c
    integer :: i, j, l

    write(error_unit,*) 'cscaletest():'

    eps = epsilon(eps)**(1.0/3.0)

    cgrad = cscale_grad(tau,rho,delta,sx,sy)
    do i = 1, size(cgrad)
       do l = -1,1,2
          z(1) = tau
          z(2:) = delta
          z(i) = z(i) + l*eps(i)
          c(l) = cscale(z(1),rho,z(2:),sx,sy)
       end do
       dc(i) = (c(1) - c(-1)) / (2*eps(i))
       write(error_unit,*) i,real(cgrad(i)),real(dc(i)),real(dc(i)-cgrad(i))
    end do


  end subroutine cscaletest

  subroutine linetest0(tau,sigu,sigv,su,sv,delta)

    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in) :: tau,sigu,sigv,su,sv

    real(REAL64), dimension(5+size(delta)) :: fd,fvec,eps,z
    real(REAL64), dimension(5+size(delta),5+size(delta)) :: fjac, djac
    real(REAL64), dimension(-1:1,5+size(delta)) :: dg
    real(REAL64), dimension(-1:1) :: f
    real(REAL64), dimension(2) :: c
    integer :: i, j, l, iflag

    write(error_unit,*) 'linetest():'

    eps = epsilon(eps)**(1.0/3.0)
    !    fvec = gradf(tau,sigu,sigv,su,sv,delta)
    z(1:5) = [ tau, sigu, sigv, su, sv ]
    z(6:) = delta
    iflag = 1
    call liked0(size(fvec),size(fvec),z,fvec,iflag)
    do i = 1, size(fvec)
       do l = -1,1,2
          z(1:5) = [ tau, sigu, sigv, su, sv ]
          z(6:) = delta
          z(i) = z(i) + l*eps(i)
          c = z(4:5) - scale
!          f(l) = loglike0(z) - sum(lambda*c) + (nu/2)*sum(c**2)
       end do
       fd(i) = (f(1) - f(-1)) / (2*eps(i))
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
            real(fd(i)-fvec(i)),real(eps(i))
    end do
    return

    fjac = jacf0(tau,sigu,sigv,su,sv,delta)
    do i = 1,size(fjac,1)
       do j = i, size(fjac,2)
          do l = -1,1,2
             z = [ tau, sigu, sigv, su, sv ]
!             z(1) = tau
!             z(2:) = delta
             z(i) = z(i) + l*eps(i)
             dg(l,:) = gradf0(z(1),z(2),z(3),z(4),z(5),z(6:))
          end do
          djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps(i))

          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do

  end subroutine linetest0



  real(REAL64) function loglike0(z)

    use rfun

    real(REAL64), dimension(:), intent(in) :: z

    real(REAL64), dimension(size(z)-5) :: delta
    real(REAL64) :: tau,sigu2, sigv2,d,r,mu,sigx,sigy,rhox,rhoy,sigu,sigv,su,sv
    integer :: i, n

    tau = z(1)
    sigu = z(2)
    sigv = z(3)
    su = z(4)
    sv = z(5)
    delta = z(6:)

    n = size(pht)
    sigu2 = sigu**2
    sigv2 = sigv**2

    loglike0 = 0
    do i = 1, n

       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
       sigx = sqrt(mu  + sigv2)
       sigy = sqrt(tau*mu + sigu2)

       if( robust ) then
          rhox = itukey(d/sigx/sv)
          rhoy = itukey(r/sigy/su)
       else
          rhox = (d/sigx/sv)**2/2
          rhoy = (r/sigy/su)**2/2
       end if

       loglike0 = loglike0 + rhox + rhoy + log(sigx) + log(sigy)

    end do
    loglike0 = loglike0 / n + log(su) + log(sv)

  end function loglike0




  subroutine line1test(tau,rho,s,delta)

    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in) :: tau,rho,s

    real(REAL64), dimension(1+size(delta)) :: fd,fvec,eps,z
    real(REAL64), dimension(1+size(delta),1+size(delta)) :: fjac, djac
    real(REAL64), dimension(-1:1,1+size(delta)) :: dg
    real(REAL64), dimension(-1:1) :: f
    real(REAL64), dimension(2) :: c
    integer :: i, j, l

    write(error_unit,*) 'line1test():'

    eps = epsilon(eps)**(1.0/3.0)
    fvec = gradf1(tau,rho,s,delta)
    do i = 1, size(fvec)
       do l = -1,1,2
!          z(1:5) = [ tau, sigu, sigv, su, sv ]
!          z(6:) = delta
          z(1) = tau
          z(2:) = delta
          z(i) = z(i) + l*eps(i)
          c = z(4:5) - scale
          f(l) = loglike1(z,rho,s) !- sum(lambda*c) + (nu/2)*sum(c**2)
       end do
       fd(i) = (f(1) - f(-1)) / (2*eps(i))
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
            real(fd(i)-fvec(i)),real(eps(i))
    end do
    return

!    fjac = jacf1(tau,sigu,sigv,su,sv,delta)
    do i = 1,size(fjac,1)
       do j = i, size(fjac,2)
          do l = -1,1,2
!             z = [ tau, sigu, sigv, s ]
             z(1) = tau
             z(2:) = delta
             z(i) = z(i) + l*eps(i)
 !            dg(l,:) = gradf1(z(1),sigu,sigv,su,sv,z(2:))
 !            dg(l,:) = loglike(z(1),z(2),z(3),z(4))
          end do
          djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps(i))

          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do

  end subroutine line1test



  function loglike1(z,rho,s) result(loglike)

    use rfun

    real(REAL64), dimension(:), intent(in) :: z
    real(REAL64), intent(in) :: rho,s

    real(REAL64), dimension(size(z)-1) :: delta
    real(REAL64) :: tau,r,mu,sigx,sigy,rhox,rhoy,loglike
    integer :: i, n

    tau = z(1)
!    sigu = z(2)
!    sigv = z(3)
!    su = z(4)
    !    sv = z(5)
    delta = z(2:)


    n = size(pht)

    loglike = 0
    do i = 1, n

!       rx = pht(i) - mju(i)
       !       ry = cts(i) - mju(i)*tau
       mu = pht(i) + delta(i)
       r = (cts(i) - tau*mu)*(1 + tau)
       sigx = sqrt(mu*(1 + 1/tau) + (rho/tau)**2)
       sigy = sqrt(mu*tau*(1 + tau) + rho**2)

       if( robust ) then
          rhox = itukey(delta(i)/sigx/s)
          rhoy = itukey(r/sigy/s)
       else
          rhox = (delta(i)/sigx/s)**2/2
          rhoy = (r/sigy/s)**2/2
       end if

       loglike = loglike + rhox + rhoy + log(sigx) + log(sigy)

    end do
    loglike = loglike / n + log(s)

  end function loglike1


  subroutine graph_1d(tau,rho,s,delta,idx,dt)

    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in) :: tau,rho,s,dt
    integer :: idx


    real(REAL64), dimension(1+size(delta)) :: z, z0
    integer :: i

!    z(1:5) = [ tau, sigu, sigv, su, sv ]
!    z(6:) = delta
    z(1) = tau
    z(2:) = delta
    z0 = z

    open(1,file='/tmp/1d')
    do i = -100, 100, 5
       z = z0
       z(idx) = z0(idx) + i * dt
       write(1,*) z(idx),loglike1(z,rho, s)
    end do
    close(1)

  end subroutine graph_1d

  subroutine thexscale(tau,rho,delta,sx,reli)

    use fmm

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in out) :: sx
    logical, intent(out) :: reli

    real(REAL64) :: smin, smax, tol

    smin = 0.3
    smax = 3
    tol = 1e-3
    write(*,*) 'thexscale(): ',real(smin),' < rho < ',real(smax),real(tol)

    backtau = tau
    backrho = rho

    allocate(backd,source=delta)
    sx = fmin(smin,smax,sxfun,tol)
    reli = smin + tol < sx .and. sx < smax - tol

    deallocate(backd)

    write(*,*) 'thexscale() estimate=',real(sx),reli

  end subroutine thexscale

  function sxfun(sx)

    use rfun

    real(REAL64) :: sxfun
    real(REAL64), intent(in) :: sx
    real(REAL64), dimension(:), pointer :: delta
    real(REAL64) :: s, d, mu, r, z, sig, ipsi, tau, rho, rhot2, tt
    integer :: i,n

    delta => backd
    n = size(delta)
    tau = backtau
    rho = backrho
    tt = tau*(tau + 1)
    rhot2 = (rho/tau)**2

    s = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       sig = sqrt(mu*(1/tau + 1) + rhot2) / 2

       z = (d/sig)/sx
       if( robust ) then
          ipsi = itukey(z)
       else
          ipsi = z**2/2
       end if

       s = s + ipsi
    end do

    sxfun = abs(2*s/n - scale)**2

  end function sxfun

  subroutine theyscale(tau,rho,delta,sy,reli)

    use fmm

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in out) :: sy
    logical, intent(out) :: reli

    real(REAL64) :: smin, smax, tol

    smin = 0.3
    smax = 3
    tol = 1e-3
    write(*,*) 'thexscale(): ',real(smin),' < rho < ',real(smax),real(tol)

    backtau = tau
    backrho = rho

    allocate(backd,source=delta)
    sy = fmin(smin,smax,syfun,tol)
    reli = smin + tol < sy .and. sy < smax - tol

    deallocate(backd)

    write(*,*) 'theyscale() estimate=',real(sy),reli

  end subroutine theyscale

  function syfun(sy)

    use rfun

    real(REAL64) :: syfun
    real(REAL64), intent(in) :: sy
    real(REAL64), dimension(:), pointer :: delta
    real(REAL64) :: s, d, mu, r, z, sig, ipsi, tau, rho, rho2, tt
    integer :: i,n

    delta => backd
    n = size(delta)
    tau = backtau
    rho = backrho
    tt = tau*(tau + 1)
    rho2 = rho**2

    s = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
       sig = sqrt(mu*tau*(tau + 1) + rho2) / 2

       z = (r/sig)/sy
       if( robust ) then
          ipsi = itukey(z)
       else
          ipsi = z**2/2
       end if

       s = s + ipsi
    end do

    syfun = abs(2*s/n - scale)**2

  end function syfun


  subroutine thescale(tau,delta,sx,sy,rho,reli)

    use fmm

    real(REAL64), intent(in) :: tau, sx, sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), intent(in out) :: rho
    logical, intent(out) :: reli

    real(REAL64) :: rhomin, rhomax, tol

    rhomin = 0
    rhomax = 3*rho
    tol = rho / 1e3
    write(*,*) 'thescale(): 0 < rho < ',real(rhomax),'  tol:',real(tol)

    backtau = tau
    backsx = sx
    backsy = sy

    allocate(backd,source=delta)
    rho = fmin(rhomin,rhomax,scalefun1,tol)
    reli = rho < rhomax - tol

    block
      integer :: i
      real(REAL64) :: rho, e, df


      e = epsilon(e)**(1.0/3.0)
      open(1,file='/tmp/rho')
      do i = 1, 100
         rho = i*(rhomax/100.0)
         df = (scalefun1(rho+e) - scalefun1(rho-e))/(2*e)
         write(1,*) rho, scalefun1(rho), df,cscale(tau,rho,delta,sx,sy)
      end do
      close(1)
    end block

    deallocate(backd)

    write(*,*) 'thescale() estimate=',real(rho),reli


  end subroutine thescale

  function scalefun(rho)

    use rfun

    real(REAL64) :: scalefun
    real(REAL64), intent(in) :: rho
    real(REAL64), dimension(:), pointer :: delta
    real(REAL64) :: s, d, mu, r, z, sig, ipsi, tau, rho2, tt
    integer :: i,n

    delta => backd
    n = size(delta)
    tau = backtau
    tt = tau*(tau + 1)
    rho2 = rho**2

    s = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
       sig = sqrt(mu*tau*(tau + 1) + rho2) / 2
!!$       sigy = sqrt(mu*tt + rho2) / 2
!!$       sigx = sigy / tau
!!$
!!$       rx = (d/sigx)/sx
!!$       ry = (r/sigy)/sy
!!$       if( robust ) then
!!$          rhox = itukey(rx)
!!$          rhoy = itukey(ry)
!!$!          rhox = max(min(rx,rmax),rmin)**2
!!$!          rhoy = max(min(ry,rmax),rmin)**2
!!$       else
!!$          rhox = rx**2/2
!!$          rhoy = ry**2/2
!!$       end if

!       s = s + rhox + rhoy

!       sig = sqrt(mu*tau*(tau + 1) + rho2) / 2
!       s = s + ((cts(i) - tau*mu)/sig)**2/2 + log(sig)


       z = r / sig
       if( robust ) then
          ipsi = itukey(z)
!          rho = max(min(z,rmax),rmin)**2/2
       else
          ipsi = z**2/2
       end if

       s = s + ipsi + log(sig)

    end do

!    scalefun = abs(s/(2*n) - 2*scale)
!    scalefun = abs(s/n - scale)
    scalefun = s / n

!    write(*,*) rho,scalefun

  end function scalefun

  function scalefun1(rho)

    use rfun

    real(REAL64) :: scalefun1
    real(REAL64), intent(in) :: rho
    real(REAL64), dimension(:), pointer :: delta
    real(REAL64) :: s, d, mu, r, zx, zy, sigx, sigy, ipsi, tau, rho2, tt, sx, sy
    integer :: i,n

    delta => backd
    n = size(delta)
    tau = backtau
    sx = backsx
    sy = backsy
    tt = tau*(tau + 1)
    rho2 = rho**2

    s = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       r = cts(i) - tau*mu
       sigy = sqrt(mu*tau*(tau + 1) + rho2) / 2
       sigx = sigy / tau

       zx = (d / sigx) / sx
       zy = (r / sigy) / sy
       if( robust ) then
          ipsi = itukey(zx) + itukey(zy)
!          rho = max(min(z,rmax),rmin)**2/2
       else
          ipsi = (zx**2 + zy**2)/2
       end if

       s = s + ipsi + log(sigx*scale) + log(sigy*scale)

!       if(

!       s = s + cfun(max(min(zx,rmax),rmin)) +

    end do

    scalefun1 = s / n + log(sx) + log(sy)

  end function scalefun1

  function cfun(z,sig)

    real(REAL64) :: cfun
    real(REAL64), intent(in) :: z,sig

    cfun = z**2/2 + log(sig)

  end function cfun


  function liscale(tau,rho,delta)

    use rfun

    real(REAL64) :: liscale
    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64) :: s, mu, z, d, sig, tt, rho2
    integer :: i, n

    n = size(delta)
    tt = tau*(tau + 1)
    rho2 = rho**2

    s = 0
    do i = 1, n
       mu = pht(i) + delta(i)
       sig = sqrt(mu*tt + rho2) / 2

       z = (cts(i) - tau*mu) / sig
       if( robust ) then
          d = itukey(z)
!          d = max(min(z,rmax),rmin)**2/2
       else
          d = z**2/2
       end if

       s = s + (d + log(sig))

    end do

    liscale = s / n

  end function liscale

  function liscale_grad(tau,rho,delta) result(dc)

    use rfun

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(:), allocatable :: dc

    real(REAL64) :: mu, psi, sig, sig_t, sig_r, sig_d, z, z_t, z_r, z_d, &
         tt, rho2
    integer :: i,n

    n = size(delta)
    tt = tau*(tau + 1)
    rho2 = rho**2
    allocate(dc(n+1))

    dc = 0
    do i = 1, n
       mu = pht(i) + delta(i)
       sig = sqrt(mu*tt + rho2) / 2

       z = (cts(i) - tau*mu) / sig
       if( robust ) then
          psi = tukey(z)
!          d = max(min(z,rmax),rmin)**2/2
       else
          psi = z**2/2
       end if

!       call grad_sigy(sig,mu,tau,rho,sig_t,sig_r,sig_d)

       z_t = -2*(mu/sig + z*(sig_t/sig))
       z_r = -2*z*(sig_r/sig)
       z_d = -2*(tau/sig + z*(sig_d/sig))

!       dc(1:2) = dc(1:2) + [psi*z_t + sig_t/sig, psi*z_r + sig_r/sig ]
       dc(1) = dc(1) + (psi*z_t + sig_t/sig)
       dc(i+1) = psi*z_d + sig_d/sig

    end do

    dc = dc / n

  end function liscale_grad

  function cscale(tau,rho,delta,sx,sy)

    use rfun

    real(REAL64) :: cscale
    real(REAL64), intent(in) :: tau,rho, sx, sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64) :: s, mu, r, d, zx, zy, psix, psiy, sigx, sigy, tt, rho2, &
         zx_r, zy_r, sigx_r, sigy_r, dx, dy
    integer :: i, n

    n = size(delta)
    tt = tau*(tau + 1)
    rho2 = rho**2

    s = 0
    do i = 1, n
       mu = pht(i) + delta(i)

       if( mu > 0 ) then
          d = delta(i)
          r = cts(i) - tau*mu
          sigy = sqrt(mu*tt + rho2) / 2
          sigx = sigy / tau

          zx = (d / sigx) / sx
          zy = (r / sigy) / sy
          if( robust ) then
             psix = tukey(zx)
             psiy = tukey(zy)
             !          d = max(min(z,rmax),rmin)**2/2
          else
             psix = zx
             psiy = zy
          end if

!          if( robust ) then
!             zx = max(min(zx,rmax),rmin)
!             zy = max(min(zy,rmax),rmin)
!          else


          sigy_r = rho / sigy / 4
          sigx_r = sigy_r / tau

          zx_r = -(d/sigx)*(sigx_r/sigx)
          zy_r = -(r/sigy)*(sigy_r/sigy)

          dx = psix*zx_r + sigx_r/sigx
          dy = psiy*zy_r + sigy_r/sigy

!          dx = max(min(dx,rmax),rmin)
!          dy = max(min(dy,rmax),rmin)

          s = s + dx + dy
!          s = s + psix*zx_r + psiy*zy_r + sigx_r/sigx + sigy_r/sigy
!          s = s + (psi*z_r + sig_r/sig)

       else
          s = s + 1e10
       end if

    end do

    cscale = s / n

  end function cscale

  function cscale_grad(tau,rho,delta,sx,sy) result(dc)

    use rfun

    real(REAL64), intent(in) :: tau, rho, sx, sy
    real(REAL64), dimension(:), intent(in) :: delta
    real(REAL64), dimension(:), allocatable :: dc

    real(REAL64) :: r, d, mu, psix, psiy, dpsix, dpsiy, sigx, sigy, &
         sigx_t, sigx_r, sigx_d, sigy_t, sigy_r, sigy_d, zx, zy, &
         zx_t, zx_r, zx_d, zy_t, zy_r, zy_d, zx_rt, zx_rd, zy_rt, zy_rd, &
         tt, rho2, sigx_tr, sigx_rd, sigy_tr, sigy_rd
    integer :: i,n

    n = size(delta)
    tt = tau*(tau + 1)
    rho2 = rho**2
    allocate(dc(n+1))

    dc = 0
    do i = 1, n
       d = delta(i)
       mu = pht(i) + d
       if( .not. ((mu*tt + rho2) > 0) ) then
!          write(*,*) i,mu,rho,delta(i)
          mu = pht(i)
       end if
       sigy = sqrt(mu*tt + rho2) / 2
       sigx = sigy / tau

       r = cts(i) - tau*mu
       zx = (d / sigx) / sx
       zy = (r / sigy) / sy
       if( robust ) then
          psix = tukey(zx)
          psiy = tukey(zy)
          dpsix = dtukey(zx)
          dpsiy = dtukey(zy)
!          d = max(min(z,rmax),rmin)**2/2
       else
          psix = zx
          psiy = zy
          dpsix = 1
          dpsiy = 1
       end if

!       call grad_sigx(sigx,mu,tau,rho,sigx_t,sigx_r,sigx_d)
!       call grad_sigy(sigy,mu,tau,rho,sigy_t,sigy_r,sigy_d)

       sigx_tr = -(rho/sigx)*(sigx_t/sigx)/4
       sigx_rd = -(rho/sigx)*(sigx_d/sigx)/4

       sigx_tr = -(rho/sigx)*(sigx_t/sigx)/4
       sigx_rd = -(rho/sigx)*(sigx_d/sigx)/4

       sigy_tr = -(rho/sigy)*(sigy_t/sigy)/4
       sigy_rd = -(rho/sigy)*(sigy_d/sigy)/4

       zx_t = -(d/sigx)*(sigx_t/sigx)
       zx_r = -zx*(sigx_r/sigx)
       zx_d = -(d/sigx)*(sigx_d/sigx)

       zy_t = -(mu/sigy + (r/sigy)*(sigy_t/sigy))
       zy_r = -zy*(sigy_r/sigy)
       zy_d = -(tau/sigy + (r/sigy)*(sigy_d/sigy))

       zy_rt = -(-(mu/sigy + 2*(r/sigy)*(sigy_t/sigy))*(sigy_r/sigy) &
            + (r/sigy)*(sigy_tr/sigy))
       zy_rd = -(-(tau/sigy + 2*(r/sigy)*(sigy_d/sigy))*(sigy_r/sigy) &
            + (r/sigy)*(sigy_rd/sigy))

!       dc(1:2) = dc(1:2) + [psi*z_t + sig_t/sig, psi*z_r + sig_r/sig ]
       dc(1) = dc(1) + &
            (dpsix*zx_t*zx_r + psix*zx_rt - (sigx_r/sigx)*(sigx_t/sigx) &
            + sigx_tr/sigx) &
            + (dpsiy*zy_t*zy_r + psiy*zy_rt - (sigy_r/sigy)*(sigy_t/sigy) &
            + sigy_tr/sigy)
       dc(i+1) = &
            (dpsix*zx_d*zx_r + psix*zx_rd - (sigx_r/sigx)*(sigx_d/sigx) + sigx_rd/sigx) &
            +(dpsiy*zy_d*zy_r + psiy*zy_rd - (sigy_r/sigy)*(sigy_d/sigy) + sigy_rd/sigy)


    end do

    dc = dc / n

  end function cscale_grad

  subroutine thescale1(rho,drho,s,reli)

    use minpacks
    use medians
    use fmm

    real(REAL64), intent(in out) :: rho,s
    real(REAL64), intent(out) :: drho
    logical, intent(out) :: reli

    integer, parameter :: npar = 2
    integer, parameter :: maxfev = 666
    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(npar) :: x
    real(REAL64), dimension(2,2) :: cov, fjac
    character(len=2880) :: msg
    character(len=5) :: label
    real(REAL64) :: tol, smin, smax
    integer :: info, n
    logical :: w

!    write(*,*) real(pht)
!    write(*,*) real(mju)
    write(*,*) real(rho),real(s), real(tauback)


    n = size(pht)
    x = [ rho, s ]
    tol = 10*n*epsilon(x)*(1 + abs(rho))

!    tol = 1e-3
!    write(error_unit,*) 'tol',tol

    if( analytic ) then
       label = 'hybrj'
!       call hybrj3(scalelike,x,tol,maxfev,factor,nprint,info,msg)
    else
       label = 'hybrd'
!       call hybrd3(scaleliked,x,tol,maxfev,factor,nprint,info,msg)
    end if
    reli = abs(x(1)) < rhomax .and. 0.1 < x(2) .and. x(2) < 10 &
         .and. (info == 1 .or. (3 <= info .and. info <= 5))


    if( debug_iter ) then
       write(error_unit,'(3a,i0,3a,l2)') ' thescale() ',label,' info=', &
            info,': ',trim(msg),'  reliable:',reli
    end if

    if( reli ) then
       rho = abs(x(1))
       s = x(2)

!       fjac = xdpsisum(mju,pht,sigv,s)
!       call qrinv(fjac,cov)

       if( cov(1,1) > 0 ) then
          drho = sqrt(cov(1,1)/(n - npar))
       else
          drho = -1
       end if

    if( debug_iter ) then
       block
         real(REAL64), dimension(2,2) :: cor
         real(REAL64) :: q,r
         integer :: i,j

         cor = 1
         do i = 1, size(cor,1)
            do j = i+1, size(cor,2)
               if( cov(i,i) > tol .and. cov(j,j) > tol ) then
                  q = sqrt(cov(i,i)*cov(j,j))
                  write(*,*) cov(i,i),cov(j,j),q
                  cor(i,j) = cov(i,j) / q
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         write(error_unit,*) 'Covariance (Hessian in minimum):'
         do i = 1,2
            write(error_unit,'(i2,a,*(es10.3,1x))') i,': ',cov(i,:)
         end do

         write(error_unit,*) 'Correlation:'
         do i = 1,2
            write(error_unit,'(i2,a,*(f6.3,1x))') i,': ',cor(i,:)
         end do

!         do i = 1, -size(pht)
!            r = (pht(i) - mju(i))/sqrt(mju(i) + sigv**2)
!            w = abs(r) < 1.5
!            write(error_unit,'(i3,2(2x,g0),2x,f7.3,l3)') &
!                 i,real(pht(i)),real(mju(i)),real(r),w
!         end do

       end block
    end if




    else
       !       dsigv = -1

       smin = 0
       smax = rhomax
       tol = 1+0.01*rho

       if( debug_iter ) write(error_unit,*) &
            'thescale failback on the interval: ',real(smin),real(smax)
       rho = fmin(smin,smax,thescalemin,tol)
       reli = rho < smax - 2*tol

       if( debug_iter ) then
          s = median(abs(pht - mju)/sqrt(mju + rho**2)) / Qdev
          write(error_unit,'(a,g0,4x,a,f0.4,l2)') &
               'rho mad=',real(rho),'s=',s,reli
       end if

    end if

    stop 'thescale'
  end subroutine thescale1

  subroutine scaleliked(npar,x,fvec,iflag)

    integer, intent(in) :: npar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(npar), intent(out) :: fvec

    real(REAL64) :: rho, s

    rho = x(1)
    s = x(2)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
       end if

       if( abs(rho) > rhomax .or. s < 0.1 .or. s > 10 ) iflag = -1
       return
    end if

    fvec = grads(tauback,mju,rho,s)
   ! fvec(2) = fvec(2) - lambda(1) + nu*(s - scale)

  end subroutine scaleliked

  function grads(tau,mju,rho,s) result(vec)

    use rfun

    real(REAL64), intent(in) :: tau,rho,s
    real(REAL64), dimension(:), intent(in) :: mju
    real(REAL64), dimension(2) :: vec
    real(REAL64) :: mu,d,r,sigx,sigy,psix,psiy, sigy_r,sigx_r, dr,ds
    integer :: i, n

    n = size(pht)

    vec = 0
    do i = 1, n

!       d = delta(i)
       !       mu = pht(i) + d
       d = mju(i) - pht(i)
       mu = mju(i)
       r = (cts(i) - mu*tau)*(1 + tau)
       sigx = sqrt(mu*(1 + 1/tau) + (rho/tau)**2)
       sigy = sqrt(mu*tau*(1+tau) + rho**2)
       sigx_r = (rho/sigx)/tau**2
       sigy_r = rho / sigy

       if( robust ) then
          psix = tukey((d / sigx)/s)
          psiy = tukey((r / sigy)/s)
       else
          psix = (d / sigx) / s
          psiy = (r / sigy) / s
       end if

       dr = -(psix/s)*(d/sigx)*(sigx_r/sigx) + sigx_r/sigx &
            -(psiy/s)*(r/sigy)*(sigy_r/sigy) + sigy_r/sigy
       ds = (1 - ((psix/s)*(d/sigx) + (psiy/s)*(r/sigy))) / s
       vec = vec + [ dr, ds ]

    end do
    vec = vec / n

  end function grads

  real(REAL64) function thescalemin(rho)

    use medians
    use quantilmean

    real(REAL64), intent(in) :: rho

    real(REAL64), dimension(:), allocatable :: xcdf, ycdf
    real(REAL64), dimension(2), parameter :: q = [1, 3] / real(4.0,REAL64)
    real(REAL64), dimension(2) :: x
    character(len=80) :: errmsg
    real(REAL64) :: xmad
    integer :: i, n, stat

!    xmad = median(abs(pht - mju)/sqrt(mju + sigv**2))
!    xscalemin = abs(xmad/Qdev - scale)

    n = size(pht)
    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,'(2a)') 'Error in thescalemin: ',trim(errmsg)
       return
    end if

    ! empirical CDF
!    call ecdf((pht - mju)/sqrt(mju + sigv**2),xcdf,ycdf)
!    call ecdf((pht - cts/tauback)/sqrt(mju*(1 + 1/tauback) + sigv**2),xcdf,ycdf)
    call ecdf((pht - cts/tauback)/sqrt(mju*(1 + 1/tauback) + rho**2/tauback**2),xcdf,ycdf)

    do i = 1, 2
       x(i) = quantile(q(i),xcdf,ycdf)
    end do
    xmad = (x(2) - x(1)) / 2
    thescalemin = abs(xmad/Qdev - scale)

  end function thescalemin


  subroutine scaletest(tau,rho,s,lambda,nu)

!    real(REAL64), dimension(:), intent(in) :: x,mju
    real(REAL64), intent(in) :: tau,rho,s,lambda,nu

    real(REAL64), dimension(2) :: fd,fvec,eps,z
    real(REAL64), dimension(2,2) :: fjac, djac
    real(REAL64), dimension(-1:1,2) :: dg
    real(REAL64), dimension(-1:1) :: f
    integer :: i, j, l

    write(error_unit,*) 'scaletest():'

    eps = epsilon(eps)**(1.0/3.0)
    fvec = grads(tau,mju,rho,s)
    fvec(2) = fvec(2) - lambda + nu*(s - scale)
    do i = 1, size(fvec)
       do l = -1,1,2
          z = [ rho, s ]
          z(i) = z(i) + l*eps(i)
          f(l) = sloglike(mju,tau,z(1),z(2),lambda,nu)
       end do
       fd(i) = (f(1) - f(-1)) / (2*eps(i))
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
            real(fd(i)-fvec(i)),real(eps(i))
    end do

    return
!    fjac = xdpsisum(mju,x,sigv,s)
    fjac(2,2) = fjac(2,2) + nu
    do i = 1,size(fjac,1)
       do j = i, size(fjac,2)
          do l = -1,1,2
             z = [ rho, s ]
             z(i) = z(i) + l*eps(i)
!             dg(l,:) = xpsisum(mju,x,z(1),z(2))
             dg(l,2) = dg(l,2) - lambda + nu*(z(2) - scale)
          end do
          djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps(i))

          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do

  end subroutine scaletest

  real(REAL64) function sloglike(mju,tau,rho,s,lambda,nu)

    use rfun

    real(REAL64), dimension(:), intent(in) :: mju
    real(REAL64), intent(in) :: tau,rho,s,lambda,nu

    real(REAL64) :: d,r,c,mu,sigx,sigy,rhox,rhoy,rx,ry
    integer :: i, n

    n = size(pht)

    sloglike = 0
    do i = 1, n

       d = mju(i) - pht(i)
       mu = mju(i)
       r = (cts(i) - tau*mu)*(1 + tau)
       sigx = sqrt(mu*(1 + 1/tau)  + (rho/tau)**2)
       sigy = sqrt(mu*tau*(1 + tau) + rho**2)

       r = (cts(i) - tau*mu)
       sigx = sqrt(mu  + (rho*cos(14/RAD)/tau)**2)
       sigy = sqrt(mu*tau + (rho*sin(14/RAD))**2)

       rx = (d/sigx)/s
       ry = (r/sigy)/s
       if( robust ) then
          rhox = itukey(rx)
          rhoy = itukey(ry)
       else
          rhox = (rx)**2/2
          rhoy = (ry)**2/2
       end if

       sloglike = sloglike + rhox + rhoy + log(sigx) + log(sigy)
!       sloglike = sloglike + rhox + 0*rhoy + log(sigx) + log(sigy)*0

    end do
    c = s - scale
    sloglike = sloglike / n + log(s) - lambda*c + (nu/2)*c**2

  end function sloglike


  subroutine xscale1(sigv,dsigv,s,reli)

    use minpacks
    use medians

    real(REAL64), intent(in out) :: sigv, s
    real(REAL64), intent(out) :: dsigv
    logical, intent(out) :: reli

    real(REAL64), dimension(3,3) :: A
    real(REAL64), dimension(3) :: b,d,x
    real(REAL64) :: l
    integer :: iter

    reli = .false.

    l = 0
    x = -0.5

    l = 0
    !    l = lambda(2)

    s = scale
    sigv = median(abs(mju - pht))
    write(*,*) sigv,s

    x = [ sigv, s, l]
!    x = [ 1e3, 1.0, -1.0 ]
    do iter = 1, 17

       A = 0
       b = 0
       A(1:2,1:2) = xdpsisum(mju,pht,sigv,s)
       b(1:2) = xpsisum(mju,pht,sigv,s)

       A(2,2) = A(2,2) + 2*l
       A(2,3) = 2*(s - scale)
       A(3,2) = A(2,3)
!       b(2) = b(2) + 2*l*(s - scale)
       b(3) = (s - scale)**2

!       l = x(3)
!       A(1,:) = [ 2*real(l), 0.0, 2*real(x(1))]
!       A(2,:) = [    0.0, 2*real(l), 2*real(x(2))]
!       A(3,:) = [ 2*real(x(1)),2*real(x(2)),0.0]
!       b = [1 + 2*x(1)*l, 1 + 2*x(2)*l, sum(x(1:2)**2)-2]

       call qrsolve(A,-b,d)
       x = x + d
       sigv = x(1)
       s = x(2)
       l = x(3)
       l = d(3)

       write(*,*) real(x(1:2)),real(l)

!       reli = abs(d(2)) < 0.01
       if( reli ) exit

    end do
    reli = .true.

    write(error_unit,*) 'xscale:',real(sigv),reli

!    stop 'x'

  end subroutine xscale1

  subroutine xscale2(sigv,dsigv,s,reli)

    use minpacks
    use medians

    real(REAL64), intent(in out) :: sigv, s
    real(REAL64), intent(out) :: dsigv
    logical, intent(out) :: reli

    real(REAL64), dimension(3,3) :: A
    real(REAL64), dimension(3) :: b,d,x
    real(REAL64) :: l
    integer :: iter

    reli = .false.

    l = 0
    x = -0.5

    l = 0
    !    l = lambda(2)

    s = scale
    sigv = median(abs(mju - pht))
    write(*,*) sigv,s

    x = [ sigv, s, l]
!    x = [ 1e3, 1.0, -1.0 ]
    do iter = 1, 17

       A = 0
       b = 0
       A(1:2,1:2) = xdpsisum(mju,pht,sigv,s)
       b(1:2) = xpsisum(mju,pht,sigv,s)

!       A(2,2) = A(2,2) + 2*l
       A(2,3) = -1
       A(3,2) = A(2,3)
       b(2) = b(2) - l
       b(3) = -(s - scale)

!       l = x(3)
!       A(1,:) = [ 2*real(l), 0.0, 2*real(x(1))]
!       A(2,:) = [    0.0, 2*real(l), 2*real(x(2))]
!       A(3,:) = [ 2*real(x(1)),2*real(x(2)),0.0]
!       b = [1 + 2*x(1)*l, 1 + 2*x(2)*l, sum(x(1:2)**2)-2]

       call qrsolve(A,-b,d)
       x = x + d
       sigv = x(1)
       s = x(2)
       l = x(3)
!       l = d(3)

       write(*,*) real(x)

!       reli = abs(d(2)) < 0.01
       if( reli ) exit

    end do
    reli = .true.

    write(error_unit,*) 'xscale:',real(sigv),reli

!    stop 'x'

  end subroutine xscale2


  subroutine yscale1(tau,sigu,dsigu,s,reli)

    use minpacks
    use medians

    real(REAL64), intent(in) :: tau
    real(REAL64), intent(in out) :: sigu, s
    real(REAL64), intent(out) :: dsigu
    logical, intent(out) :: reli

    real(REAL64), dimension(3,3) :: A
    real(REAL64), dimension(3) :: b,d,x
    real(REAL64) :: l
    integer :: iter

    reli = .false.

    s = scale
    sigu = median(abs(cts - tau*mju))
    write(*,*) sigu,s

    l = 0
   ! l = lambda(2)
    x = [ sigu, s, l]
    do iter = 1, 17

       A = 0
       b = 0
       A(1:2,1:2) = ydpsisum(mju,cts,tau,sigu,s)
       b(1:2) = ypsisum(mju,cts,tau,sigu,s)

       A(2,2) = A(2,2) + 2*l
       A(2,3) = 2*(s - scale)
       A(3,2) = A(2,3)
       b(2) = b(2) + 2*l*(s - scale)
       b(3) = (s - scale)**2

       call qrsolve(A,-b,d)
       x = x + d
       sigu = x(1)
       s = x(2)
       l = x(3)
!       l = d(3)

       write(*,*) real(x(1:2)),real(l)

!       reli = abs(d(2)) < 0.01
       if( reli ) exit

    end do
    reli = .true.

    write(error_unit,*) 'yscale:',real(sigu),reli

  end subroutine yscale1


  subroutine xscale(sigv,dsigv,s,reli)

    use minpacks
    use medians
    use fmm

    real(REAL64), intent(in out) :: sigv, s
    real(REAL64), intent(out) :: dsigv
    logical, intent(out) :: reli

    integer, parameter :: npar = 3
    integer, parameter :: maxfev = 666
    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(npar) :: x
    real(REAL64), dimension(2,2) :: cov, fjac
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: tol, smin, smax
    integer :: info, n
    logical :: w

!    write(*,*) real(pht)
!    write(*,*) real(mju)
    write(*,*) real(sigv),real(s), real(tauback)


    n = size(pht)
    x(1:2) = [ sigv, s ]
    x(3) = 0.1
    tol = 10*n*epsilon(x)*(1 + abs(sigv))

!    tol = 1e-3
!    write(error_unit,*) 'tol',tol

    if( analytic ) then
       label = 'hybrj'
!       call hybrj3(xlikesig,x,tol,maxfev,factor,nprint,info,msg)
    else
       label = 'hybrd'
!       call hybrd3(xlikedsig,x,tol,maxfev,factor,nprint,info,msg)
    end if
    reli = abs(x(1)) < rhomax .and. 0.1 < x(2) .and. x(2) < 10 &
         .and. (info == 1 .or. (3 <= info .and. info <= 5))


    if( debug_iter ) then
!!$       msg = ''
!!$       if( info == 1 .or. (3 <= info .and. info <= 5) ) then
!!$          msg = ' (success)'
!!$       else if( info == -1 ) then
!!$          msg = ' (out of range)'
!!$       else if( info == 1 ) then
!!$          msg = ' (no convergence)'
!!$       end if
       write(error_unit,'(3a,i0,3a,l2)') ' xscale() ',label,' info=', &
            info,': ',trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       ! Sometimes, the lambda parameter of LM-method is incidentally found
       ! in the mirrored part of the parameter space. Fortunately, the function
       ! is pretty symmetric, so I do to ignore the "tunnelling". The correct,
       ! in principle, solution should to keep the condition: sigv > 0.
       sigv = abs(x(1))
       s = x(2)

       fjac = xdpsisum(mju,pht,sigv,s)
       call qrinv(fjac,cov)

       if( cov(1,1) > 0 ) then
          dsigv = sqrt(cov(1,1)/(n - npar))
       else
          dsigv = -1
       end if

    if( debug_iter ) then
       block
         real(REAL64), dimension(2,2) :: cor
         real(REAL64) :: q,r
         integer :: i,j

         cor = 1
         do i = 1, size(cor,1)
            do j = i+1, size(cor,2)
               if( cov(i,i) > 0 .and. cov(j,j) > 0 ) then
                  q = sqrt(cov(i,i)*cov(j,j))
                  write(*,*) cov(i,i),cov(j,j),q
                  cor(i,j) = cov(i,j) / q
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         write(error_unit,*) 'Covariance (Hessian in minimum):'
         do i = 1,2
            write(error_unit,'(i2,a,*(es10.3,1x))') i,': ',cov(i,:)
         end do

         write(error_unit,*) 'Correlation:'
         do i = 1,2
            write(error_unit,'(i2,a,*(f6.3,1x))') i,': ',cor(i,:)
         end do

         do i = 1, -size(pht)
            r = (pht(i) - mju(i))/sqrt(mju(i) + sigv**2)
            w = abs(r) < 1.5
            write(error_unit,'(i3,2(2x,g0),2x,f7.3,l3)') &
                 i,real(pht(i)),real(mju(i)),real(r),w
         end do

       end block
    end if




    else
       !       dsigv = -1

       smin = 0
       smax = rhomax
       tol = 1+0.01*sigv

       if( debug_iter ) write(error_unit,*) &
            'xscale failback on the interval: ',real(smin),real(smax)
       sigv = fmin(smin,smax,xscalemin,tol)
       reli = sigv < smax - 2*tol

       if( debug_iter ) then
          s = median(abs(pht - mju)/sqrt(mju + sigv**2)) / Qdev
          write(error_unit,'(a,g0,4x,a,f0.4,l2)') &
               'sigv mad=',real(sigv),'s=',s,reli
       end if

    end if




!    stop 'xscale'
  end subroutine xscale

  function xpsilike(mu,x,sigv,s) result(psi)

    use rfun

    real(REAL64), intent(in) :: mu,x,sigv,s
    real(REAL64), dimension(2) :: psi
    real(REAL64) :: rx, sigx, psix, sigx_v

    rx = x - mu
    sigx = sqrt(mu + sigv**2)
    sigx_v = sigv/sigx

    if( robust ) then
       psix = tukey((rx / sigx) / s)
    else
       psix = (rx / sigx) / s
    end if

    psi(1) = (1 - (psix/s)*(rx/sigx))*(sigx_v/sigx)
    psi(2) = -(psix/s**2)*(rx/sigx) + 1/s

  end function xpsilike


  function xdpsilike(mu,x,sigv,s) result(dpsi)

    use rfun

    real(REAL64), intent(in) :: mu,x,sigv,s
    real(REAL64), dimension(2,2) :: dpsi
    real(REAL64) :: rx, sigx, psix, dpsix, sigx_v, sigx_vv

    rx = x - mu
    sigx = sqrt(mu + sigv**2)
    sigx_v = sigv/sigx
    sigx_vv = -(sigv/sigx)*(sigx_v/sigx) + 1/sigx

    if( robust ) then
       psix = tukey((rx / sigx) / s)
       dpsix = dtukey((rx / sigx) / s)
    else
       psix = (rx / sigx) / s
       dpsix = 1
    end if

    dpsi(1,1) = (dpsix/s**2)*((rx/sigx)*(sigx_v/sigx))**2 &
         - (psix/s)*(rx/sigx)*((sigx_vv/sigx) - 2*(sigx_v/sigx)**2) &
         - (sigx_v/sigx)**2 + sigx_vv/sigx

    dpsi(1,2) = (dpsix/s**3)*(rx/sigx)**2*(sigx_v/sigx) &
         + (psix/s**2)*(rx/sigx)*(sigx_v/sigx)
    dpsi(2,2) = (dpsix/s**4)*(rx/sigx)**2 + 2*(psix/s**3)*(rx/sigx) &
         - 1/s**2

    dpsi(2,1) = dpsi(1,2)

  end function xdpsilike


  function xpsisum(mju,x,sigv,s) result(sums)

    real(REAL64), dimension(:), intent(in) :: mju,x
    real(REAL64), intent(in) :: sigv,s
    real(REAL64), dimension(2) :: sums

    real(REAL64), dimension(2) :: psi
    integer :: i, n

    n = size(x)

    sums = 0
    !    do concurrent (i = 1:n)
    do i = 1,n
       psi = xpsilike(mju(i),x(i),sigv,s)
!       write(*,*) real(psi),real((mju(i)-x(i))/sqrt(mju(i) + sigv**2))
       sums = sums + psi
    end do
    sums = sums / n
!    stop

  end function xpsisum

  function xdpsisum(mju,x,sigv,s) result(dsums)

    real(REAL64), dimension(:), intent(in) :: mju,x
    real(REAL64), intent(in) :: sigv,s
    real(REAL64), dimension(2,2) :: dsums

    real(REAL64), dimension(2,2) :: dpsi
    integer :: i, n

    n = size(x)

    dsums = 0
!    do concurrent (i = 1:n)
    do i = 1,n
       dpsi = xdpsilike(mju(i),x(i),sigv,s)
       dsums = dsums + dpsi
    end do
    dsums = dsums / n

  end function xdpsisum


  subroutine xscaletest(x,mju,sigv,s,lambda,nu)

    real(REAL64), dimension(:), intent(in) :: x,mju
    real(REAL64), intent(in) :: sigv,s,lambda,nu

    real(REAL64), dimension(2) :: fd,fvec,eps,z
    real(REAL64), dimension(2,2) :: fjac, djac
    real(REAL64), dimension(-1:1,2) :: dg
    real(REAL64), dimension(-1:1) :: f
    integer :: i, j, l

    write(error_unit,*) 'xscaletest():'

    eps = epsilon(eps)**(1.0/3.0)
    fvec = xpsisum(mju,x,sigv,s)
    fvec(2) = fvec(2) - lambda + nu*(s - scale)
    do i = 1, size(fvec)
       do l = -1,1,2
          z = [ sigv, s ]
          z(i) = z(i) + l*eps(i)
!          f(l) = xeloglikes(mju,x,z(1),z(2),lambda,nu)
       end do
       fd(i) = (f(1) - f(-1)) / (2*eps(i))
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
            real(fd(i)-fvec(i)),real(eps(i))
    end do

    fjac = xdpsisum(mju,x,sigv,s)
    fjac(2,2) = fjac(2,2) + nu
    do i = 1,size(fjac,1)
       do j = i, size(fjac,2)
          do l = -1,1,2
             z = [ sigv, s ]
             z(i) = z(i) + l*eps(i)
             dg(l,:) = xpsisum(mju,x,z(1),z(2))
             dg(l,2) = dg(l,2) - lambda + nu*(z(2) - scale)
          end do
          djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps(i))

          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do

  end subroutine xscaletest

!  subroutine xlikesig(npar,mpar,x,fvec,fjac,ldfjac,iflag)
  subroutine xlikesig(npar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(npar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

!    real(REAL64), dimension(2) :: psi
!    real(REAL64), dimension(2,2) :: dpsi
    real(REAL64) :: sigv, s, lambda
    integer :: i,n

    n = size(cts)
    sigv = x(1)
    s = x(2)
    lambda = x(3)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
          do i = 1, npar
             write(error_unit,'(a,i0,a,*(1x,es10.3))') &
                  'fjac(',i,',:) = ',real(fjac(i,:))
          end do
       end if

       ! convergence check
       if( abs(sigv) > rhomax .or. s < 0.1 .or. s > 10 ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

       fvec(1:2) = xpsisum(mju,pht,sigv,s)
!       fvec(2) = fvec(2) - lambda(1) + nu*(s - scale)
       fvec(2) = fvec(2) - lambda
       fvec(3) = -(s - scale)

    else if( iflag == 2 ) then

       fjac(1:2,1:2) = xdpsisum(mju,pht,sigv,s)
       !       fjac(2,2) = fjac(2,2) + nu
       fjac(3,:) = [ 0, -1, 0 ]
       fjac(:,3) = [ 0, -1, 0 ]

!       do concurrent (i = 2:npar)
!          do concurrent (j = 1:i-1)
!             fjac(i,j) = fjac(j,i)
!          end do
!       end do

    end if

  end subroutine xlikesig

!  subroutine xlikedsig(npar,mpar,x,fvec,iflag)
  subroutine xlikedsig(npar,x,fvec,iflag)

    integer, intent(in) :: npar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(npar), intent(out) :: fvec

!    real(REAL64), dimension(2) :: psi
    real(REAL64) :: sigv, s, lambda

    sigv = x(1)
    s = x(2)
    lambda = x(3)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
       end if

       if( abs(sigv) > rhomax .or. s < 0.1 .or. s > 10 ) iflag = -1
       return
    end if

    fvec(1:2) = xpsisum(mju,pht,sigv,s)
    fvec(2) = fvec(2) - lambda
    fvec(3) = -(s - scale)
!    fvec(2) = fvec(2) - lambda(1) + nu*(s - scale)

  end subroutine xlikedsig


!  real(REAL64) function xscalemin(sigv)
  real(REAL64) function xscalemin(rho)

    use medians
    use quantilmean
    use photonrate_normal

!    real(REAL64), intent(in) :: sigv
    real(REAL64), intent(in) :: rho

    real(REAL64), dimension(:), allocatable :: xcdf, ycdf
    real(REAL64), dimension(2), parameter :: q = [1, 3] / real(4.0,REAL64)
    real(REAL64), dimension(2) :: x
    character(len=80) :: errmsg
    real(REAL64) :: xmad, tau
    integer :: i, n, stat

!!$    block
!!$      real(REAL64) :: tau,phi,rho
!!$
!!$      rho = 2e3
!!$      phi = sigv
!!$      tau = tauback
!!$
!!$      xmad = median(abs(pht - mju)/sqrt(mju + (rho/tau*cos(phi))**2))
!!$      xscalemin = abs(xmad/Qdev - scale)**2
!!$      write(*,*) phi,xmad/Qdev
!!$      return
!!$    end block

    tau = tauback
!    xmad = median(abs(pht - mju)/sqrt(mju + (sigv/tauback)**2))
    xmad = median(abs(pht - mju)/sqrt(mju*(1 + 1/tau) + (rho/tau)**2))*2
    xscalemin = abs(xmad/Qdev - scale)
    write(*,*) rho,xmad/Qdev,xscalemin
    return

    n = size(pht)
    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,'(2a)') 'Error in xscalemin: ',trim(errmsg)
       return
    end if

    ! empirical CDF
!    call ecdf((pht - mju)/sqrt(mju + sigv**2),xcdf,ycdf)
!    call ecdf((pht - cts/tauback)/sqrt(mju*(1 + 1/tauback) + sigv**2),xcdf,ycdf)
!    call ecdf((pht - cts/tauback)/sqrt(mju*(1 + 1/tauback) + sigv**2/tauback**2),xcdf,ycdf)

    do i = 1, 2
       x(i) = quantile(q(i),xcdf,ycdf)
    end do
    xmad = (x(2) - x(1)) / 2
    xscalemin = abs(xmad/Qdev - scale)

!    smax = 0
!    do i = 1, size(pht)
!       smax = max(smax,abs(ycdf(i) - Phi(xcdf(i))))
!    end do
    xscalemin = maxval(abs(ycdf - NormCDF(xcdf)), abs(xcdf) < 1)
!    write(*,*) sigv,xscalemin

  end function xscalemin


!  real(REAL64) function yscalemin(sigu)
  real(REAL64) function yscalemin(rho)

    use medians
    use quantilmean
    use photonrate_normal

    real(REAL64), intent(in) :: rho!,sigu

    real(REAL64), dimension(:), allocatable :: xcdf, ycdf
    real(REAL64), dimension(2), parameter :: q = [1, 3] / real(4.0,REAL64)
    real(REAL64), dimension(2) :: x
    character(len=80) :: errmsg
    real(REAL64) :: ymad, tau
    integer :: i, n, stat

    tau = tauback
!    ymad = median(abs(cts - tau*mju)/sqrt(mju*tau + sigu**2))
    ymad = median(abs(cts - tau*mju)/sqrt(mju*tau*(tau + 1) + rho**2))*2
    ymad = median(abs(cts - tau*pht)/sqrt((cts + tau*pht)/2*(tau + 1) + rho**2))
    yscalemin = abs(ymad/Qdev - scale)
    yscalemin = abs(ymad/Qdev - 1)
    write(*,*) rho,ymad/Qdev,yscalemin
    return

    n = size(pht)
    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,'(2a)') 'Error in yscalemin: ',trim(errmsg)
       return
    end if

    ! empirical CDF
    !    call ecdf((cts - mju*tauback)/sqrt(mju*tauback + sigu**2),xcdf,ycdf)
    call ecdf(2*(cts - tau*mju)/sqrt(mju*tau*(tau + 1) + rho**2),xcdf,ycdf)

    do i = 1, 2
       x(i) = quantile(q(i),xcdf,ycdf)
    end do
    ymad = (x(2) - x(1)) / 2
    yscalemin = abs(ymad/Qdev - scale)

!    yscalemin = maxval(abs(ycdf - NormCDF(xcdf)), abs(xcdf) < 1)
    write(*,*) rho,yscalemin


  end function yscalemin



  subroutine yscale(sigu,dsigu,s,reli)

    use minpacks

    real(REAL64), intent(in out) :: sigu, s
    real(REAL64), intent(out) :: dsigu
    logical, intent(out) :: reli

    integer, parameter :: npar = 3
    integer, parameter :: maxfev = 666
    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(npar) :: par
    real(REAL64), dimension(2,2) :: cov, fjac
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: tol
    integer :: info, n

!    write(*,*) real(cts)
!    write(*,*) real(mju)
!    write(*,*) real(sigu),real(s), real(tauback)


    n = size(cts)
    par(1:2) = [ sigu, s ]
    par(3) = 0.1
    tol = 10*n*epsilon(par)*(1 + abs(sigu))

    if( analytic ) then
       label = 'hybrj'
!       call hybrj3(ylikesig,par,tol,maxfev,factor,nprint,info,msg)
    else
       label = 'hybrd'
!       call hybrd3(ylikedsig,par,tol,maxfev,factor,nprint,info,msg)
    end if

    reli = abs(par(1)) < rhomax .and. 0.1 < par(2) .and. par(2) < 10 &
         .and. (info == 1 .or. (3 <= info .and. info <= 5))
!         .and. (0 < info .and. info < 5)

    if( debug_iter ) then
!!$       msg = ''
!!$       if ( info == 1 .or. (3 <= info .and. info <= 5) ) then
!!$          msg = ' (success)'
!!$       else if( info == -1 ) then
!!$          msg = ' (out of range)'
!!$       else if( info == 1 ) then
!!$          msg = ' (no convergence)'
!!$       end if
       write(error_unit,'(3a,i0,3a,l2)') ' yscale() ',label,' info=', &
            info,': ',trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       fjac = ydpsisum(mju,cts,tauback,par(1),par(2))
       call qrinv(fjac,cov(1:2,1:2))

       sigu = abs(par(1))
       s = par(2)

       if( cov(1,1) > 0 ) then
          dsigu = sqrt(cov(1,1)/(n - npar))
       else
          dsigu = -1
       end if

    if( debug_iter ) then
       block
         real(REAL64), dimension(2,2) :: cor
         real(REAL64) :: q
         integer :: i,j

         cor = 1
         do i = 1, size(cor,1)
            do j = i+1, size(cor,2)
               if( cov(i,i) > 0 .and. cov(j,j) > 0 ) then
                  q = sqrt(cov(i,i)*cov(j,j))
                  cor(i,j) = cov(i,j) / q
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         write(error_unit,*) 'Covariance (Hessian in minimum):'
         do i = 1, 2
            write(error_unit,'(i2,a,*(es10.3,1x))') i,': ',cov(i,:)
         end do

         write(error_unit,*) 'Correlation:'
         do i = 1,2
            write(error_unit,'(i2,a,*(f6.3,1x))') i,': ',cor(i,:)
         end do

       end block
    end if


    else
       dsigu = -1
    end if


!    stop 'yscale'


  end subroutine yscale


  function ypsilike(mu,y,tau,sigu,s) result(psi)

    use rfun

    real(REAL64), intent(in) :: mu,y,tau,sigu,s
    real(REAL64), dimension(2) :: psi
    real(REAL64) :: ry, sigy, psiy, sigy_u

    ry = y - mu*tau
    sigy = sqrt(mu*tau + sigu**2)
    sigy_u = sigu/sigy

    if( robust ) then
       psiy = tukey((ry / sigy) / s)
    else
       psiy = (ry / sigy) / s
    end if

    psi(1) = (1 - (psiy/s)*(ry/sigy))*(sigy_u/sigy)
    psi(2) = -(psiy/s**2)*(ry/sigy) + 1/s

  end function ypsilike


  function ydpsilike(mu,y,tau,sigu,s) result(dpsi)

    use rfun

    real(REAL64), intent(in) :: mu,y,tau,sigu,s
    real(REAL64), dimension(2,2) :: dpsi
    real(REAL64) :: ry, sigy, psiy, dpsiy, sigy_u, sigy_uu

    ry = y - mu*tau
    sigy = sqrt(mu*tau + sigu**2)
    sigy_u = sigu/sigy
    sigy_uu = 1/sigy - (sigu/sigy)*(sigy_u/sigy)

    if( robust ) then
       psiy = tukey((ry / sigy) / s)
       dpsiy = dtukey((ry / sigy) / s)
    else
       psiy = (ry / sigy) / s
       dpsiy = 1
    end if

    dpsi(1,1) = (dpsiy/s**2)*((ry/sigy)*(sigy_u/sigy))**2 &
         - (psiy/s)*(ry/sigy)*((sigy_uu/sigy) - 2*(sigy_u/sigy)**2) &
         - (sigy_u/sigy)**2 + sigy_uu/sigy

    dpsi(1,2) = (dpsiy/s**3)*(ry/sigy)**2*(sigy_u/sigy) &
         + (psiy/s**2)*(ry/sigy)*(sigy_u/sigy)
    dpsi(2,2) = (dpsiy/s**4)*(ry/sigy)**2 + 2*(psiy/s**3)*(ry/sigy) &
         - 1/s**2

    dpsi(2,1) = dpsi(1,2)

  end function ydpsilike


  function ypsisum(mju,y,tau,sigu,s) result(sums)

    real(REAL64), dimension(:), intent(in) :: mju,y
    real(REAL64), intent(in) :: tau,sigu,s
    real(REAL64), dimension(2) :: sums

    real(REAL64), dimension(2) :: psi
    integer :: i, n

    n = size(y)

    sums = 0
    !    do concurrent (i = 1:n)
    do i = 1,n
       psi = ypsilike(mju(i),y(i),tau,sigu,s)
       sums = sums + psi
    end do
    sums = sums / n

  end function ypsisum

  function ydpsisum(mju,y,tau,sigu,s) result(dsums)

    real(REAL64), dimension(:), intent(in) :: mju,y
    real(REAL64), intent(in) :: tau,sigu,s
    real(REAL64), dimension(2,2) :: dsums

    real(REAL64), dimension(2,2) :: dpsi
    integer :: i, n

    n = size(y)

    dsums = 0
!    do concurrent (i = 1:n)
    do i = 1,n
       dpsi = ydpsilike(mju(i),y(i),tau,sigu,s)
       dsums = dsums + dpsi
    end do
    dsums = dsums / n

  end function ydpsisum


!  subroutine ylikesig(npar,mpar,x,fvec,fjac,ldfjac,iflag)
  subroutine ylikesig(npar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(npar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

!    real(REAL64), dimension(2) :: psi
!    real(REAL64), dimension(2,2) :: dpsi
    real(REAL64) :: sigu, s, lambda
    integer :: i

!    n = size(cts)
    sigu = x(1)
    s = x(2)
    lambda = x(3)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec),':',real(sum(fvec**2))
          do i = 1, npar
             write(error_unit,'(a,i0,a,*(1x,es10.3))') &
                  'fjac(',i,',:) = ',real(fjac(i,:))
          end do
       end if

       ! convergence check
       if( abs(sigu) > rhomax .or. s < 0.1 .or. s > 10 ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

       fvec(1:2) = ypsisum(mju,cts,tauback,sigu,s)
       fvec(2) = fvec(2) - lambda
       fvec(3) = -(s - scale)
!       fvec(2) = fvec(2) - lambda(2) + nu*(s - scale)

    else if( iflag == 2 ) then

       fjac(1:2,1:2) = ydpsisum(mju,cts,tauback,sigu,s)
       fjac(3,:) = [ 0, -1, 0 ]
       fjac(:,3) = [ 0, -1, 0 ]

!       fjac(2,2) = fjac(2,2) + nu

!       do concurrent (i = 2:npar)
!          do concurrent (j = 1:i-1)
!             fjac(i,j) = fjac(j,i)
!          end do
!       end do

    end if

  end subroutine ylikesig

  subroutine ylikedsig(npar,x,fvec,iflag)

    integer, intent(in) :: npar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(npar), intent(out) :: fvec

!    real(REAL64), dimension(2) :: psi
    real(REAL64) :: sigu, s, lambda

    sigu = x(1)
    s = x(2)
    lambda = x(3)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x),':',real(sum(fvec**2))
          write(error_unit,*) 'fvec = ',real(fvec)
       end if

       if( abs(sigu) > rhomax .or. s < 0.1 .or. s > 10 ) iflag = -1
       return
    end if

    fvec(1:2) = ypsisum(mju,cts,tauback,sigu,s)
    fvec(2) = fvec(2) - lambda
    fvec(3) = -(s - scale)

    !fvec(2) = fvec(2) - lambda(2) + nu*(s - scale)

  end subroutine ylikedsig


  subroutine yscaletest(y,mju,tau,sigu,s,lambda,nu)

    real(REAL64), dimension(:), intent(in) :: y,mju
    real(REAL64), intent(in) :: tau,sigu,s,lambda,nu

    real(REAL64), dimension(2) :: fd,fvec,eps,z
    real(REAL64), dimension(2,2) :: fjac, djac
    real(REAL64), dimension(-1:1,2) :: dg
    real(REAL64), dimension(-1:1) :: f
    integer :: i, j, l

    write(error_unit,*) 'yscaletest():'

    eps = epsilon(eps)**(1.0/3.0)
    fvec = ypsisum(mju,y,tau,sigu,s)
    fvec(2) = fvec(2) - lambda + nu*(s - scale)
    do i = 1, size(fvec)
       do l = -1,1,2
          z = [ sigu, s ]
          z(i) = z(i) + l*eps(i)
          f(l) = yeloglikes(mju,y,tau,z(1),z(2),lambda,nu)
       end do
       fd(i) = (f(1) - f(-1)) / (2*eps(i))
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
            real(fd(i)-fvec(i)),real(eps(i))
    end do

    fjac = ydpsisum(mju,y,tau,sigu,s)
    fjac(2,2) = fjac(2,2) + nu
    do i = 1,size(fjac,1)
       do j = i, size(fjac,2)
          do l = -1,1,2
             z = [ sigu, s ]
             z(i) = z(i) + l*eps(i)
             dg(l,:) = ypsisum(mju,y,tau,z(1),z(2))
             dg(l,2) = dg(l,2) - lambda + nu*(z(2) - scale)
          end do
          djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps(i))

          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do

  end subroutine yscaletest

  real(REAL64) function xeloglikes(mju,x,tau,sigv,s,lambda,nu)

    use rfun

    real(REAL64), dimension(:), intent(in) :: mju,x
    real(REAL64), intent(in) :: tau,sigv,s,lambda,nu

    real(REAL64) :: ss,rx,sigx
    integer :: n,i

    n = size(x)
    ss = 0
    do concurrent (i = 1:n)
       rx = x(i) - mju(i)
       sigx = sqrt(mju(i) + sigv**2/tau**2)
       if( robust ) then
          ss = ss + itukey(rx/sigx/s) + log(sigx)
       else
          ss = ss + ((rx/sigx)/s)**2/2 + log(sigx)
       end if
    end do

    xeloglikes = ss / n + log(s) - lambda*(s - scale) + (nu/2)*(s - scale)**2

  end function xeloglikes

  real(REAL64) function yeloglikes(mju,y,tau,sigu,s,lambda,nu)

    use rfun

    real(REAL64), dimension(:), intent(in) :: mju,y
    real(REAL64), intent(in) :: tau,sigu,s,lambda,nu

    real(REAL64) :: ss,ry,sigy
    integer :: n,i

    n = size(y)
    ss = 0
    do concurrent (i = 1:n)
       ry = y(i) - mju(i)*tau
       sigy = sqrt(mju(i)*tau + sigu**2)
       if( robust ) then
          ss = ss + itukey(ry/sigy/s) + log(sigy)
       else
          ss = ss + ((ry/sigy)/s)**2/2 + log(sigy)
       end if
    end do

    yeloglikes = ss / n + log(s) - lambda*(s - scale) + (nu/2)*(s - scale)**2

  end function yeloglikes

  subroutine graph_xscale(filename,sigv)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: sigv
    real(REAL64) :: v,s,mu
    integer :: i,j

    open(1,file=filename)
    do i = 0, 250, 5
!       v = 0.01*i*sigv
       v = 0.01*i*5e2
       do j = 10, 200, 5
          s = 0.01*j
!          write(1,*) v,s,xeloglikes(mju,pht,tauback,v,s,lambda(1),nu)
       end do
       write(1,*)
    end do
    close(1)

    s = 1
    open(1,file='/tmp/m')
    do i = 900, 1100, 1
       mu = 0.001*i*mju(1)
       do j = 50, 150, 1
          v = 0.01*j*sigv
!          write(1,*) mu,v,xeloglikes([mu],pht(1:1),tauback,v,s,lambda(1),nu)
       end do
       write(1,*)
    end do
    close(1)


  end subroutine graph_xscale


  subroutine graph_yscale(filename,sigu)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: sigu
    real(REAL64) :: u,s,mu,t
    integer :: i,j

    open(1,file=filename)
    do i = 0, 250, 5
       u = 0.01*i*sigu
       do j = 10, 200, 5
          s = 0.01*j
!          write(1,*) u,s,yeloglikes(mju,cts,tauback,u,s,lambda(2),nu)
       end do
       write(1,*)
    end do
    close(1)

    s = 1
    open(1,file='/tmp/n')
    do i = 900, 1100, 1
       mu = 0.001*i*mju(1)
       do j = 50, 150, 1
          u = 0.01*j*sigu
!          write(1,*) mu,u,yeloglikes([mu],cts(1:1),tauback,u,s,lambda(2),nu)
       end do
       write(1,*)
    end do

    s = 1
    open(1,file='/tmp/t')
    do i = 900, 1100, 1
       t = 0.001*i*tauback
       do j = 50, 150, 1
          u = 0.01*j*sigu
 !         write(1,*) t,u,yeloglikes(mju,cts,t,u,s,lambda(2),nu)
       end do
       write(1,*)
    end do


  end subroutine graph_yscale



  subroutine graph_ecdf(tau,rho,phi,cts,pht,delta)

    ! empirical CDFs for residuals in x,y

    use quantilmean

    real(REAL64), intent(in) :: tau,rho,phi
    real(REAL64), dimension(:), intent(in) :: cts, pht, delta

    real(REAL64), dimension(:), allocatable :: xpdf, ypdf, res, dlt
    real(REAL64) :: sigu, sigv
    character(len=80) :: errmsg
    integer :: n,i,stat

    sigu = rho*sin(phi)
    sigv = rho*cos(phi)/tau

    n = size(cts)
    allocate(xpdf(n),ypdf(n),res(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,*) 'Error: `',trim(errmsg),"'(graph_ecdf)"
       return
    end if

    ! pht
    res = 2.0 * delta / sqrt(pht + delta + sigv**2)
!    res = delta / sqrt((pht + delta)*(1+1/tau) + (rho/tau)**2)/2
    res = 2*delta/sqrt((pht+delta)*(1+1/tau) + (rho/tau)**2)
    call ecdf(res,xpdf,ypdf)

    open(1,file='/tmp/xcdf')
    do i = 1, size(xpdf)
       write(1,*) xpdf(i), ypdf(i)
    end do
    close(1)

    ! d,r check
!    write(*,*) n
    allocate(dlt(n))
    block
      real(REAL64) :: mag, mean

      do i = 1, n
         mag = 10 + (i-1) / 20.0
         mean = 10**((25 - mag)/2.5)
         dlt(i) = sqrt((pht(i) - mean)**2 + (cts(i)/tau - mean)**2)
         dlt(i) = (pht(i) - mean)**2 / (mean*(1 + 1/tau) + rho**2/tau**2) &
              + (cts(i)/tau - mean)**2 / (mean*(1 + 1/tau) + rho**2/tau**2)
         dlt(i) = sqrt(dlt(i)*2)
!         write(*,*) mag,mean,dlt(i)
      end do
    end block

!    res = (pht - cts/tau)/sqrt(mju*(1+1/tau) + (rho/tau)**2)
    res = (cts - tau*(pht + dlt))/sqrt(mju*(tau + tau**2) + rho**2)
    res = dlt / sqrt(mju*(tau + 0*tau**2) + rho**2/tau**2)
    res = (pht - cts/tau) / sqrt((pht+dlt)*(1 + 1/tau) + rho**2/tau**2)
    res = (pht - cts/tau) / sqrt(mju*(1 + 1/tau) + rho**2/tau**2)
    res = (cts - tau*(pht + dlt))/sqrt((pht+dlt)*(tau+tau**2) + rho**2)
!    res = (cts - tau*mju)/sqrt(mju*tau*(1+tau) + rho**2)
    !    res = (cts - tau*pht)/sqrt(pht*(tau+tau**2) + rho**2)
    res = (1+tau)*(cts - tau*(pht+delta))/sqrt((pht+delta)*(tau+tau**2) + rho**2)
    !         +(pht+delta)*(1+1/tau) + (rho/tau)**2)
    res = (cts - tau*(pht+delta)) / sqrt((pht+delta)*tau + sigu**2) * sqrt(2.0)/sqrt(tau)
    res = 2*(cts - tau*(pht+delta)) / sqrt((pht+delta)*tau*(1+tau) + rho**2)
    call ecdf(res,xpdf,ypdf)

    open(1,file='/tmp/ycdf')
    do i = 1, size(xpdf)
       write(1,*) xpdf(i), ypdf(i)
    end do
    close(1)


    ! distance residuals, Chi distribution: P(1/2,x**2/2)
    ! https://en.wikipedia.org/wiki/Chi_distribution
    res = 2**2*delta**2/sqrt((pht+delta)*(1+1/tau) + (rho/tau)**2)**2 &
         +2**2*(cts-tau*(pht+delta))**2/sqrt((pht+delta)*tau*(tau+1)+rho**2)**2
    res = 2**2*((tau*delta)**2 + (cts-tau*(pht+delta))**2) &
         /sqrt((pht+delta)*tau*(tau+1)+rho**2)**2
    call ecdf(sqrt(res/2),xpdf,ypdf)

    open(1,file='/tmp/dcdf')
    do i = 1, size(xpdf)
       write(1,*) xpdf(i), ypdf(i) ! igamma(0.5,x**2/2)
    end do
    close(1)


  end subroutine graph_ecdf

  subroutine graph_muline(tau,rho,delta,cts,pht)

    real(REAL64), intent(in) :: tau, rho
    real(REAL64), dimension(:), intent(in) :: cts, pht, delta

!    real(REAL64), dimension(:), allocatable :: xpdf, ypdf, res
    real(REAL64) :: t,s,slim,tmin,tmax,mu,d
    integer :: i

!    rho = sqrt(sigu**2 + tau**2*sigv**2)

    open(1,file='/tmp/ml')

    do i = 1, size(cts)

       d = delta(i)
       mu = pht(i) + d
       t = mu*tau
       s = sqrt(mu*tau*(tau+1) + rho**2)/2
       !       s = sqrt(pht(i)+delta + rho**2*sin(phi)**2)
       slim = 1*s
       tmin = t - slim
       tmax = t + slim
       write(1,*) pht(i),cts(i),t,tmin,tmax,cts(i)-t,slim,-slim,mu,tau*mu

    end do
    close(1)

    ! gnuplot> plot 'ml', 'ml' u 9:3 w l, 'ml' u 9:4 w l, 'ml' u 9:5 w l, 'ml' u 9:10
    ! gnuplot> set logscale x
    ! gnuplot> plot 'ml' u 1:6, 'ml' u 1:7 w l, 'ml' u 1:8 w l

  end subroutine graph_muline


  subroutine graph_rhos(tau,rho,cts,pht,delta,sx,sy)

    real(REAL64), intent(in) :: tau,rho,sx,sy
    real(REAL64), dimension(:), intent(in) :: cts,pht,delta

    real(REAL64), dimension(size(delta)+2) :: dd!,z
    real(REAL64) :: r,t,ll,s,d,f,d0
    integer :: i, j, l

!    z(1) = tau
!    z(2) = rho
!    z(3) = phi
!    z(3:) = delta


!!$    open(1,file='/tmp/rs')
!!$    do i = 20, 300
!!$       s = 0.01*i
!!$       do j = 0,5000, 100
!!$          r = j
!!$          z(2) = r
!!$          z(3) = s
!!$!          z(4) = s
!!$          lx = loglike(z)
!!$!          z(3) = s
!!$!          z(4) = s
!!$!          ly = loglike(z)
!!$          write(1,*) r,s,lx!,ly
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)


!    open(1,file='/tmp/tr')

!    z(3) = 1
!    z(4) = 1
!    z(:) = delta

!    do i = 1900, 2100,10
!       t = 0.001*i
!       z(1) = t
!       do j = 0,int(2*rho), max(int(2*rho/100),1)
!          r = j
!          z(2) = r
!          ll = loglike(z,rho,sx,sy)
!          write(1,*) t,r,ll
!       end do
!       write(1,*)
!    end do
!    close(1)
!    z(1) = tau
!    z(2) = rho

    d0 = sqrt(pht(1)*(1+1/tau) + (rho/tau)**2)/1.0
    open(1,file='/tmp/td')
    do i = 150, 250, 5
       t = 0.01*i
!       z(1) = t
       do j = -200,200,10
          d = j*d0
          do l = 1, size(delta)
!             if( pht(l) + d < 1 ) then
!                write(*,*) l,pht(l),d,pht(l) + d
             !             end if
!             z(3:) = 0
!             z(l+2) = d
!             z(2) = 1414
!             z(1) = 2
!             dd(l) = loglike(z(1:3),sx,sy)
             block
               use rfun
               real(REAL64) :: mu,r,sigx,sigy, rx,ry,rhox,rhoy

               mu = pht(l) + d
               r = cts(l) - t*mu
               sigy = sqrt(mu*t*(t+1) + rho**2)/2
               sigx = sigy / t

               rx = (d/sigx)/sx
               ry = (r/sigy)/sy
               if( robust ) then
                  rhox = itukey(rx)
                  rhoy = itukey(ry)
                  !          rhox = max(min(rx,rmax),rmin)**2/2
                  !          rhoy = max(min(ry,rmax),rmin)**2/2
               else
                  rhox = rx**2/2
                  rhoy = ry**2/2
               end if

               dd(l) = rhox + rhoy + log(sigx) + log(sigy) + log(sx) + log(sy)
             end block
          end do
          write(1,*) t,d,dd
       end do
       write(1,*)
    end do
    close(1)
!    z(1) = tau
!    z(4:) = delta


!    open(1,file='/tmp/rf')
!    do i = 100, 5000, 100
!       s = i
!       z(2) = s
!       do j = 0, 90
!          f = j / RAD
!          z(3) = f
!          ll = loglike(z,sx,sy)
!          write(1,*) s,j,ll
!       end do
!       write(1,*)
!    end do
!    close(1)


  end subroutine graph_rhos

  subroutine graph_rho(tau,rho,cts,pht,delta,lambda)

    real(REAL64), intent(in) :: tau,rho,lambda
    real(REAL64), dimension(:), intent(in) :: cts,pht,delta

    real(REAL64), dimension(size(delta)+2) :: fv
    real(REAL64), parameter :: sx = 1, sy = 1
    real(REAL64) :: r
    integer :: i

    open(1,file='/tmp/rmin')
    do i = 0,5000,100
       r = i
       fv = gradlike(tau,r,delta,sx,sy) + lambda*grad_disp(tau,r,delta)
       write(1,*) r,sum(fv**2),disp(tau,r,delta)
    end do
    close(1)

  end subroutine graph_rho

  subroutine graph_line(tau,rho,cts,pht,delta,lambda)

    real(REAL64), intent(in) :: tau,rho,lambda
    real(REAL64), dimension(:), intent(in) :: cts,pht,delta

    real(REAL64), parameter :: sx = 1, sy = 1
    real(REAL64), dimension(size(delta)) :: dd
    real(REAL64), dimension(size(delta)+2) :: fv
    real(REAL64) :: r,t,ll,s,f,d0,d
    integer :: i, j, l

    d0 = sqrt(pht(1)*(1+1/tau) + (rho/tau)**2)
    open(1,file='/tmp/lmin')
    do i = 199900, 200100, 1
       t = 0.00001*i
       do j = -5000,5000,100
          d = j*d0*0.01

          do l = 1, size(delta)
             block
               real(REAL64), dimension(:), allocatable :: xd

               xd = delta
               xd(l) = d
               fv = linegrad1(t,rho,xd) !+ lambda*grad_disp(t,rho,xd)
               dd(l) = sum(fv**2)/size(delta)
             end block
          end do
          write(1,*) t,d,dd
       end do
       write(1,*)
    end do
    close(1)


  end subroutine graph_line


!!$  subroutine graph_scale(tau,mju,cts,pht,rho)
!!$
!!$    real(REAL64), intent(in) :: tau,rho
!!$    real(REAL64), dimension(:), intent(in) :: mju,cts,pht
!!$
!!$    real(REAL64) :: r,s,ll
!!$    integer :: i, j, l
!!$
!!$    open(1,file='/tmp/rs')
!!$
!!$    do i = 20, 300
!!$       s = 0.01*i
!!$       do j = 0,5000, 100
!!$          r = j
!!$          ll = sloglike(mju,tau,r,s,lambda(1),nu)
!!$          write(1,*) r,s,ll
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graph_scale

  subroutine graph_disp(tau,rho,delta)

    real(REAL64), intent(in) :: tau,rho
    real(REAL64), dimension(:), intent(in) :: delta

    real(REAL64) :: r
    integer :: i

!    r = disp(tau,rho,delta)
!    stop

    open(1,file='/tmp/disp')

    do i = 500, 20000, 50
       r = i
       write(1,*) r,(disp(tau,r,delta)-scale)**2
    end do
    close(1)

  end subroutine graph_disp

  subroutine graph_like(filename,tau,dtau,rho,drho,mju)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: tau,dtau,rho,drho
    real(REAL64), dimension(:), intent(in) :: mju

    real(REAL64), dimension(2) :: fvec
    real(REAL64) :: x,y,c,llog
    integer :: i,j,iflag

    iflag = 1
    open(1,file=filename)
    do i = -500, 500, 5
       x = tau + 0.01*i*dtau
       do j = -500, 500, 5
          y = rho + 0.01*j*drho
          call theliked(2,2,[x,y],fvec,iflag)
          c = cvar(x,y) !/ ndata
          llog = loglike(x,y,mju)/ndata
          write(1,*) x,y,llog, sum(fvec**2), llog - lambda*c + (nu/2)*c**2,&
               c**2
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_like

  subroutine graph_rlike(filename,tau,dtau,rho,drho,mju)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: tau,dtau,rho,drho
    real(REAL64), dimension(:), intent(in) :: mju

    real(REAL64), dimension(3) :: fvec
    real(REAL64) :: x,y,c,llog
    integer :: i,j,iflag

    iflag = 1
    open(1,file=filename)
    do i = -500, 500, 5
       x = tau + 0.01*i*dtau
       do j = -500, 500, 5
          y = rho + 0.01*j*drho
          call rliked(3,[x,y,lambda],fvec,iflag)
          c = fvec(3)
          llog = rloglike(x,y,mju)/ndata
          write(1,*) x,y,llog, sum(fvec(1:2)**2), c**2, sum(fvec**2)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_rlike


  subroutine graph_fun(filename,tau,dtau,rho,drho,mju)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: tau,dtau,rho,drho
    real(REAL64), dimension(:), intent(in) :: mju

    character(len=666) :: name
    real(REAL64), dimension(2) :: fvec,cvec
    real(REAL64), dimension(2,2) :: fjac,cjac
    real(REAL64) :: x,y,lam,c,llog,f
    integer :: i,j,l

    do l = 1,1
       lam = l
       write(name,'(2a,i0)') filename,'_',l
       open(1,file=name)
       do i = -500, 500, 5
          x = tau + 0.01*i*dtau
          do j = -500, 500, 5
             y = rho + 0.01*j*drho
             call funky(x,y,lambda,f,fvec,fjac,c,cvec,cjac)
             write(1,*) x,y,f,c
!             write(1,*) x,y,f,sum(fvec**2),sum(matmul(transpose(fjac),fjac)**2),&
!                  c,sum(cvec**2),sum(matmul(transpose(cjac),cjac)**2)
          end do
          write(1,*)
       end do
       close(1)
    end do

  end subroutine graph_fun


end module photonrate_line
