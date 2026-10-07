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


module photonrate_mood

  use photonrate_common
  use iso_fortran_env

  implicit none
  private

  public :: mood, moodinit

  logical, parameter :: debug_graph = .true.
  logical, parameter :: debug_test = .false.

  ! Does to use Jacobian analytic, or numerical?
  logical, parameter :: analytic = .true.

  ! see, phopt
  integer, parameter :: nprint = 10

  ! use Huber's or Tukey's (descent) function
  logical, parameter :: descent = .false.

  ! limits
!  real(REAL64), parameter :: gmin = 0.001
!  real(REAL64), parameter :: gmax = 1000
!  real(REAL64), parameter :: phimin = 0
!  real(REAL64), parameter :: phimax = PI / 2

  integer, dimension(:), pointer :: iset
!  real(REAL64), dimension(:), pointer :: mju, cts, pht
  real(REAL64), dimension(:), pointer :: cts, pht, mju
  real(REAL64), dimension(:), allocatable :: ysel, xsel, mumin, mumax!, omju
  real(REAL64) :: tau, sx, sy, lambda,nu, scale, rhomax
  logical :: robustness

contains

    subroutine mood(y,x,isets,nset,sxs,sys,t,rho,drho,phi,dphi,mjus,kmean,robust)

    use photonrate_phiopt
    use medians
    use fmm
    use rfun

    real(REAL64), dimension(:), target, intent(in) :: x,y
    integer, dimension(:), target, intent(in) :: isets
    integer, intent(in) :: nset
    real(REAL64), intent(in) :: t
    real(REAL64), intent(in out) :: sxs,sys
    real(REAL64), intent(in out) :: rho, phi
    real(REAL64), intent(out) :: drho, dphi
    logical, intent(in) :: robust
    real(REAL64), dimension(:), target, intent(in out) :: mjus, kmean

    real(REAL64), dimension(size(kmean)) :: dkmean
    integer, dimension(:), allocatable :: idx
    character(len=80) :: errmsg
    real(REAL64) :: s, c, dmu, stol, sigu, sigv
!    real(REAL64), dimension(size(kmean)) :: dmus
    integer :: iter,i,n,k,m,stat
    logical :: reli, converge

    pht => x
    cts => y
    mju => mjus
    iset => isets
    tau = t
    sx = sxs
    sy = sys
    robustness = robust
    if( robustness ) then
       if( descent ) then
          scale = stukey
       else
          scale = shuber
       end if
    else
       scale = 1
    end if
    scale = 1

!    cts(1) = 1e10
!    pht(1) = 1

    !mumin = 0.9*min(minval(y)/tau,minval(x))
    !mumax = 1.1*min(maxval(y)/tau,maxval(x))

    allocate(mumin(nset),mumax(nset))
!    mumin = 0.9*minval(kmean)
    !    mumax = 1.1*maxval(kmean)
    mumin = max(kmean - 3*sqrt(kmean + rho**2/2)*6,sqrt(kmean))
    mumax = kmean + 3*sqrt(kmean + rho**2/2)*6
!    where( mumin < 0 ) mumin = 0
    rhomax = 10*max(rho,1.0_REAL64)

    if( debug ) then
       do i = 1, nset
          write(error_unit,*) 'mood: group=',i, &
               ' munim=',real(mumin(i)),' mumax=',real(mumax(i))
       end do
    end if


!    n = size(iset)
!    allocate(idx(n),stat=stat,errmsg=errmsg)
!    if( stat /= 0 ) then
!       if( debug ) write(error_unit,*) 'Error: `',trim(errmsg),"' (mood)"
!       return
!    end if

    ! setup mjus
    do k = 1, -nset
!       m = 0
!       do i = 1,n
!          if( iset(i) == k ) then
!             m = m + 1
!             idx(m) = i
!          end if
!       end do
!       write(*,*) idx(1:m)
       idx = picker(iset,k)
!       write(*,*) idx
       m = size(idx)
       if( m > 0 ) then
          kmean(k) = median((cts(idx)*tau + pht(idx)) / (1 + tau**2))
!          kmean(k) = median((cts(idx)/tau + pht(idx))) / 2
          mju(idx) =  kmean(k)
          if( debug ) write(error_unit,*) 'Init mu:',k,real(kmean(k))!, &
!               median((cts(idx)*tau + pht(idx)) / (1 + tau**2))
!          kmean(k) = median((cts(idx(1:m))*tau + pht(idx(1:m))) / (1 + tau**2))
!          mju(idx(1:m)) =  kmean(k)
       else
          stop 'An empty subset in mood.'
       end if
    end do

!    s = madscale(mju,cts,pht,tau,rho,sx,sy)
!    if( debug ) then
!       if( .not. (0.1 < s .and. s < 10) ) &
!         write(error_unit,*) 'Warning: the scale is out of range: ', &
!         real(s),' Data looks like the gain is uncorrected.'
!    end if
!    s = 1

    if( debug_test ) then
       k = 1
       idx = picker(iset,k)
       lambda = -1
       nu = 128
       s = 1.5
       phi = 30 /  RAD
       scale = 1
       kmean(k) = 9900
!       rho = 1414
       call rmoodtest(cts(idx),pht(idx),kmean(k),sx,sy,tau,rho,phi,s,lambda,nu)
       call xmoodtest(pht,mju,tau,rho,phi,s,lambda,nu)
       call ymoodtest(cts,mju,tau,rho,phi,s,lambda,nu)
       stop 'rmoodtest'
    end if

!!$    if( dog ) then
!!$       g = 1
!!$       phi = PI / 4
!!$       call gscale(cts,pht,mju,sx,sy,tau,rho,g,dg,phi,dphi,reli)
!!$!       if( .fail. ) then
!!$!          call ginit(s)
!!$    else
!!$       g = 1
!!$       ! optimal phi
!!$       phi = PI / 4
!!$       phi = phiopt(sx,sy,tau,rho,mju,cts,pht)
!!$    end if


!    call graph_mood(tau,rho*1.5,phi,sx,sy,cts,pht,mju)


    ! Penalty augmented method with the constrain: s == 1
    !phi = PI / 4
    stol = 0.01
!    scale = shuber
    s = scale
!    s = madscale(mju,cts,pht,tau,rho,phi,sx,sy)
    !    call phiopt(sx,sy,tau,rho,mju,cts,pht,s,phi,dphi)
!    stop

    block
      real(REAL64) :: dsigv, smin, smax, r
      logical :: w

      !      real(REAL64), dimension(size(kmean)) :: xkmean
      sigv = rho*cos(phi)/tau
      smin = 0
      smax = rhomax

      s = scale
      c = s - scale
      nu = 1
      lambda = nu*c
      converge = .false.

      do iter = 1, 17

         call xmood(mju,sigv,dsigv,s,reli)

         if( debug_iter ) then
103         format(a,i0,(2x,a,g0,a,es7.1),2x,2(1x,a,f0.4),1x,2(1x,a,f0.3),1x,l1)
            write(error_unit,103) 'xmood #',iter,'sigv=',real(sigv),'+-',dsigv,&
                 's=',s,' ds=',s-scale,'lambda=',lambda,'nu=',nu,reli
         end if

         if( debug_graph) call graph_xymood('/tmp/m','x',iter,sigu,sigv)

         c = s - scale
         converge = abs(c) < stol
         if( converge .or. .not. reli ) exit

         lambda = lambda - nu*c
         nu = 2*nu
      end do
      reli = .true.
      converge = .true.

      if( .not. (reli .and. converge) ) then
         write(error_unit,*) real(smin),real(smax)
         sigv = fmin(smin,smax,xmoodmin,1+0.01*rho)
         if( debug_iter ) then
            s = median(abs(pht - mju)/sqrt(mju + sigv**2)) / Qdev
            write(error_unit,'(a,g0,4x,a,f0.4)') &
                 'xmoodmad: sigv=',real(sigv),'s=',s
         end if
         reli = .true.
      end if

      do i = 1, size(pht)
         r = (pht(i) - mju(i))/sqrt(mju(i) + sigv**2)
         w = abs(r) < 1.5
         write(error_unit,'(i3,2(2x,g0),2x,f7.3,l3)') &
              i,real(pht(i)),real(mju(i)),real(r),w
      end do


    end block

    block
      real(REAL64) :: dsigu, smin, smax,r
      logical :: w

      !      real(REAL64), dimension(size(kmean)) :: xkmean
      sigu = rho*sin(phi)
      smin = 0
      smax = rhomax

      s = scale
      c = s - scale
      nu = 1
      lambda = nu*c
      converge = .false.

      do iter = 1, 17

         call ymood(mju,sigu,dsigu,s,reli)

         if( debug_iter ) then
104         format(a,i0,(2x,a,g0,a,es7.1),2x,2(1x,a,f0.4),1x,2(1x,a,f0.3),1x,l1)
            write(error_unit,104) 'ymood #',iter,'sigu=',real(sigu),'+-',dsigu,&
                 's=',s,' ds=',s-scale,'lambda=',lambda,'nu=',nu,reli
         end if

         if( debug_graph) call graph_xymood('/tmp/m','y',iter,sigu,sigv)

         c = s - scale
         converge = abs(c) < stol
         if( converge .or. .not. reli ) exit

         lambda = lambda - nu*c
         nu = 2*nu
      end do
!      reli = .false.

      if( .not. (reli .and. converge) ) then
         write(error_unit,*) real(smin),real(smax)
         sigu = fmin(smin,smax,ymoodmin,1+0.01*rho)
         if( debug_iter ) then
            s = median(abs(cts - mju*tau)/sqrt(mju*tau + sigu**2)) / Qdev
            write(error_unit,'(a,g0,4x,a,f0.4)') &
                 'ymoodmad: sigu=',real(sigu),'s=',s
         end if
         reli = .true.
      end if

      do i = 1, size(cts)
         r = (cts(i) - mju(i)*tau)/sqrt(mju(i)*tau + sigu**2)
         w = abs(r) < 1.5
         write(error_unit,'(i3,2(2x,g0),2x,f7.3,l2)') &
              i,real(cts(i)),real(mju(i)),real(r),w
      end do

      ! winsorize
!      do i = 1, size(cts)
!         r = (cts(i) - mju(i)*tau)/sqrt(mju(i)*tau + sigu**2)
!      end do

    end block

    rho = sqrt(sigu**2 + tau**2*sigv**2)
    phi = atan2(sigu,tau*sigv)

    if( debug ) then
       write(error_unit,*) 'Final rho=',real(rho),' phi=',real(RAD*phi)
    end if

    do k = 1, nset
       idx = picker(iset,k)
       m = size(idx)
       if( m > 0 ) then
          call rmood(sigu,sigv,tau,rho,phi,s,cts(idx),pht(idx),kmean(k),dkmean(k))
          !          kmean(k) = sum(cts(idx)/tau + pht(idx)) / (2*m)
          mju(idx) =  kmean(k)
       else
          stop 'Bad mood.'
       end if
    end do

    goto 70
!    stop


!    write(*,*) 'Initial phi=',phi

!    phi0 = phi
    c = s - scale
    nu = 1
    lambda = nu*c
    converge = .false.
    do iter = 1, 3*precision(phi)

       call allmood(kmean,dkmean,rho,drho,phi,dphi,s,reli)

       if( debug_iter ) then
102       format(a,i0,(2x,a,g0,a,es7.1),(2x,a,f4.1,a,f0.1),2(2x,a,f0.4), &
               3x,2(2x,a,f0.3),1x,l1)
          write(error_unit,102) 'mood #',iter,'rho=',real(rho),'+-',drho, &
               'phi=',RAD*phi,'+-',RAD*dphi,'s=',s,' ds=',s-scale, &
               'lambda=',lambda,'nu=',nu,reli
          write(error_unit,'(a,i0,*(2x,g0,a,es7.1))') 'mood #',iter,&
               (real(kmean(k)),'+-',dkmean(k),k=1,nset)
       end if

       c = s - scale
       converge = abs(c) < stol
       if( converge .or. .not. reli ) exit

       lambda = lambda - nu*c
       nu = 2*nu

    end do

    if( reli ) then
       mju = -1
       do k = 1, nset
          idx = picker(iset,k)
          mju(idx) = kmean(k)
       end do
    end if

    if( any(mju < 0) ) then
       write(error_unit,*) pack(mju,mju < 0)
       stop '-0'
    end if

    goto 70


    phi = PI / 4
    stol = 0.01
    scale = 1
    s = madscale(mju,cts,pht,tau,rho,phi,sx,sy)
!    call phiopt(sx,sy,tau,rho,mju,cts,pht,s,phi,dphi)
    c = s - scale
    nu = 1
    lambda = nu*c
    converge = .false.
    do iter = 1, precision(phi)

       do k = 1, nset
          idx = picker(iset,k)
          m = size(idx)
          if( m > 0 ) then
             call rmood(sx,sy,tau,rho,phi,s,cts(idx),pht(idx),kmean(k),dkmean(k))
             mju(idx) =  kmean(k)
          else
             stop 'Bad mood.'
          end if
       end do

       if( debug_iter ) then
101       format(a,i0,(2x,f4.1,a,f0.1),(2x,a,f0.4),3x,2(2x,f0.3),1x,l1)
          write(error_unit,101) 'mood #',iter, &
               RAD*phi,'+-',RAD*dphi,'s=',s,lambda,nu,reli
          write(error_unit,'(a,i0,*(2x,g0,a,es7.1))') 'mood #',iter,&
               (real(kmean(k)),'+-',dkmean(k),k=1,nset)
       end if

       !       s = madscale(mju,cts,pht,tau,rho,phi,sx,sy)

!       call rortho(mju,cts,pht,phi,dphi,s,reli)

       c = s - scale
!       converge = reli .and. abs(c) < stol
       converge = abs(c) < stol
       if( converge ) exit

       lambda = lambda - nu*c
       nu = 2*nu

    end do
    reli = converge

    !phi = PI / 4
    do iter = 1, 1

!       call phiopt(sx,sy,tau,rho,s,mju,cts,pht,phi,dphi)


       do k = 1, nset
!       m = 0
!       do i = 1, n
!          if( iset(i) == k ) then
!             m = m + 1
!             idx(m) = i
!          end if
       !       end do
          idx = picker(iset,k)
          m = size(idx)
          if( m > 2 ) then
             call robmood(sx,sy,s,cts(idx),pht(idx),kmean(k),dmu,phi,dphi)
!             kmean(k) = rmood(sx,sy,s,tau,rho,phi,cts(idx),pht(idx),kmean(k))
             mju(idx) =  kmean(k)
             !          kmean(k) = rmood(sx,sy,tau,rho,phi,g,cts(idx(1:m)),pht(idx(1:m)))
!          mju(idx(1:m)) =  kmean(k)
          else
             stop 'Bad mood.'
          end if
       end do

    end do

!    stop 666
70  continue

    if( debug_iter .and. size(mju) < 20 ) write(error_unit,*) 'mood:',real(mju)

    write(*,*) madscale(mju,cts,pht,tau,rho,phi,sx,sy)

    if( debug_graph) then
       call graph_mood(tau,rho,phi,sigu,sigv,cts,pht,mju)
       call graph_likemood('/tmp/likemood',mju(1),5*sqrt(mju(1)),rho)
       call graph_moodline(tau,rho*1.5,phi,sx,sy,cts,pht,mju,sigu,sigv)
    end if

    sxs = sx
    sys = sy

    deallocate(mumin,mumax)

  end subroutine mood

  function picker(iset,k) result(subset)

    integer, dimension(:), intent(in) :: iset
    integer, intent(in) :: k

    integer, dimension(:), allocatable :: idx, subset
    integer :: i,n,m

    n = size(iset)
    allocate(idx(n))

    m = 0
    do i = 1, n
       if( iset(i) == k ) then
          m = m + 1
          idx(m) = i
       end if
    end do

    subset = idx(1:m)

  end function picker

  function madscale(mju,cts,pht,tau,rho,phi,sx,sy) result(s)

    use medians

    real(REAL64) :: s
    real(REAL64), dimension(:), intent(in) :: mju,cts,pht
    real(REAL64), intent(in) :: tau,rho,phi
    real(REAL64), intent(out) :: sx,sy

    real(REAL64) :: xmad, ymad

    ! scale by MAD, phi is set to PI /4 (also the effective of cos()**2
    ! on the interval 0 .. PI/2)
!    xmad = median(abs(pht - mju)/sqrt(mju + (rho/tau)**2/2))
    xmad = median(abs(pht - mju)/sqrt(mju + (rho/tau*cos(phi))**2))
!    ymad = median(abs(cts - mju*tau)/sqrt(mju*(tau + 0*tau**2) + &
!         0*mju**2*dtau**2 + rho**2/2))
!    ymad = median(abs(cts - mju*tau)/sqrt(mju*(tau +0*tau**2/size(cts)) + 0*mju**2*dtau**2 + (rho)**2/2))
!    ymad = median(abs(cts - mju*tau)/sqrt(mju*tau + (rho)**2/2))
    ymad = median(abs(cts - mju*tau)/sqrt(mju*tau + (rho*sin(phi))**2))
    s = (xmad + ymad) / (2*Qdev)
    sx = xmad / Qdev
    sy = ymad / Qdev
    sx = 1
    sy = 1

    if( debug ) write(error_unit,*) 'mad scale: ', &
         real(s),'  xmad, ymad:',real(xmad/Qdev),real(ymad/Qdev)!,&
!         median((abs(pht - mju) + abs(cts - mju*tau)/tau)/sqrt(2*mju + (rho/tau)**2))/Qdev, &
!         median((abs(pht - mju)*tau + abs(cts - mju*tau)) &
!         /sqrt(mju*tau*(1 + tau) + rho**2))/Qdev

!    stop 0

  end function madscale

  subroutine rmood(sigu,sigv,tau,rho,phi,s,y,x,mu,stderr)

    use medians
    use rfun

    real(REAL64), intent(in) :: sigu,sigv,tau, rho, phi, s
    real(REAL64), dimension(:), intent(in) :: x,y
    real(REAL64), intent(in out) :: mu
    real(REAL64), intent(out) :: stderr

    real(REAL64), parameter :: mulim = 10
    real(REAL64), dimension(2) :: r, dr, sig, psi, dpsi, ddr, sc
    real(REAL64) :: mu0, f, df, d, tol, stdsig, mumin, mumax, sx, sy, fmin
    integer :: n, i, iter, nout
    character(len=1) :: method
    logical :: converge, point

    converge = .false.
    method = ' '
    n = size(y)

    if( n == 0 ) stop ' in rmood.'

!    sigu = rho*sin(phi)
!    sigv = rho*cos(phi)/tau
    !    sc = [ sx, sy ]
    sc = 1
!    sc = [ 7.0, 1.0 ]
    sx = 1
    sy  = 1

!    mu = median((y*tau + x) / (1 + tau**2))
    !    mu0 = mu
    mu0 = mu
    mumin = 0.9*min(minval(y)/tau,minval(x))
    mumax = 1.1*max(maxval(y)/tau,maxval(x))
    tol = max(10*epsilon(mu)*(1 + mu),0.1*sqrt(mu)/sqrt(real(n)))
    if( debug_iter ) write(error_unit,'(a,g0,a,es12.2,a,*(1x,g0))') &
         'mu start: ',real(mu),' mood tol: ',tol,' values:',real(x),real(y/tau)

!    mu = median(y/tau + x) / 2

!    if( abs(x(1) - 60.020e3) < 1e3 ) then
!!$    if( x(1) < 2 ) then
!!$    mu = mumin
!!$    mu0 = mu
!!$    fmin = huge(fmin)
!!$    d = (mumax - mumin) / 1e3
!!$    !    write(error_unit,*) mumin, mumax
!!$    open(1,file='/tmp/lp')
!!$    do i = 1, 1000
!!$       mu = mumin + d*i
!!$       f = (x(1) - mu)**2/(mu + sigv**2) + &
!!$            (y(1)/tau - mu)**2/(mu + sigu**2/tau**2)
!!$       f = (x(1) - mu)**2/(mu + sigv**2) + &
!!$            (y(1) - tau*mu)**2/(mu*tau + sigu**2)
!!$       sig = [ sqrt(mu/sx + sigv**2), sqrt(mu*tau/sy + sigu**2) ]
!!$       r = [ x(1) - mu, y(1) - mu*tau ] / sig
!!$       where( abs(r) > 1.5 )
!!$          r = sign(1.5_REAL64,r)
!!$       end where
!!$       f = sum(r**2)
!!$       write(1,*) mu,f
!!$!       f = (x(1) - mu)**2 + (y(1)/tau - mu)**2
!!$       if( f < fmin ) then
!!$          fmin = f
!!$          mu0 = mu
!!$ !         write(*,*) fmin,mu0
!!$       end if
!!$    end do
!!$    close(1)
!!$ end if
!    write(error_unit,*) fmin,mu0,(y(1)/tau + x(1))/2
!!$
!!$    mu = mu0
!!$    return

    point = .false.
    do iter = 1, precision(mu)

       sig = [ sqrt(mu/sx + sigv**2), sqrt(mu*tau/sy + sigu**2) ]

       f = 0
       df = 0
       do concurrent (i = 1:n)
          r = [ x(i) - mu, y(i) - mu*tau ] / sig

!          where( abs(r) > 1.345 )
!             r = sign(1.345_REAL64,r)
!          end where
          !          if( any(r > 1.34) ) write(*,*) real(r)
          nout = count(abs(r) > mulim)

          dr = -(1 + r/sig/sc) / sig
          ddr = (2 + 3*r/sig/sc)/sig**3/sc
          dr(2) = tau*dr(2)
          ddr(2) = tau**2*ddr(2)
!          call tukeys(r/s,psi,dpsi)
!          psi = huber(r/s)
!          dpsi = dhuber(r/s)
          psi = r/s
          dpsi = 1

          f = f + sum(psi*dr)
          df = df + sum(dpsi*dr**2) + sum(psi*ddr)
       end do

       if( n == 1 .and. nout == 2 ) then
          write(error_unit,*) 'r:',real(r)
          exit
       end if

!       if( n == 1 .and. all(abs(r) > 1.344) ) then
!          exit
!       end if

!       if( n == 1 .and. count(abs(r) > 1.344) == 1 ) then
       if( n == 1 .and. nout == 1 ) then
          write(error_unit,*) 'r:',real(r)
          if( abs(r(1)) < mulim ) then
             mu = x(1)
          else if( abs(r(2)) < mulim ) then
             mu = y(1) / tau
          end if
          point = .true.
          exit
       end if

       f = -f / n - (tau/sig(2)**2 + 1/sig(1)**2)
       df = -df / n + 2*((tau/sig(2)**2)**2 + 1/sig(1)**4)

!       block
!         real(REAL64), dimension(4) :: psi
!         real(REAL64), dimension(4,4) :: dpsi

!         psi = psisum(mu,y,x,tau,rho,phi,s,sx,sy,lambda,nu)
!         dpsi = dpsisum(mu,y,x,tau,rho,phi,s,sx,sy,nu)

!         f = -psi(1)
!         df = -dpsi(1,1)
!       end block

       if( .not. (df < 0) ) exit ! when a local minimum is encountered

       d = f / df
       mu = mu - d

       if( debug_iter ) write(error_unit,'(a,i2,2x,g0,2es12.2,es14.4)') &
            "mu, incr., f, f': ",iter,real(mu),d,f,df

       converge = abs(d) < tol .and. df < 0

       ! exit of iterations: the absolute error must be at least |d| < tol
       if( converge ) exit

       if( mu < mumin .or. mu > mumax ) exit

    end do

    if( converge ) then
       method = 'N'
       stdsig = 1 / sqrt(-df)
    else if( point ) then
       method = 'P'
       stdsig = 0
    else
       method = '-'
       mu = mu0
!       stdsig = 0
!!!       stdsig = median(abs((y*tau + x) / (1 + tau**2) - mu)) / Qdev
!!!       stdsig = median(abs((y/tau + x) / 2 - mu)) / Qdev
    end if

!    converge = .false.
    if( .not. converge .and. .false. ) then
       block
         real(REAL64) :: a,b,fa,fb,f
         integer :: nmax

         tol = max(10*epsilon(mu)*(1 + mu),0.5*sqrt(mu)/sqrt(real(n)))
         nmax = min(2*(nint(abs(log(tol/mu0))/log(2.0)) + 1),3*precision(mu))
         write(*,*) tol,nmax
         a = mumin
         b = mumax
         fa = psifun(a,sx,sy,s,tau,rho,phi,y,x)
         fb = psifun(b,sx,sy,s,tau,rho,phi,y,x)
         do iter = 1, nmax
            mu = (a + b) / 2
            f = psifun(mu,sx,sy,s,tau,rho,phi,y,x)

            if( debug_iter ) write(error_unit,'(a,i2,2x,3(g0:3x))') &
                 "bi-limits, a, b, tau: ",iter,real(a),real(b), real(mu)

            if( fa*f < 0 ) then
               b = mu
               fb = f
            else if( f*fb < 0 ) then
               a = mu
               fa = f
            else ! if( f == 0 .or. fa == 0 .or fb == 0 ) then
               if( debug ) write(error_unit,*) &
                    'Out of a root interval. ',real(f),real(fa),real(fb)
               exit
            end if

            converge = abs(b - a) < tol .and. abs(mumin - mu) > tol &
                 .and. abs(mumax - mu) > tol

            if( converge ) exit

         end do

         if( converge ) then
            method = 'B'
            stdsig = median(abs((y*tau + x) / (1 + tau**2) - mu)) / Qdev
         else
            mu = mu0
         end if

       end block
    end if


    if( debug ) then

       if( n > 1 ) then
          stderr = stdsig / sqrt(real(2*n-1))
       else
          stderr = stdsig
       end if
       write(error_unit,*) &
            ' mu,stderr,stdsig:',real(mu),real(stderr),real(stdsig),' ',method
    end if

  end subroutine rmood

  real(REAL64) function psifun(mu,sx,sy,s,tau,rho,phi,y,x)

    use rfun

    real(REAL64), intent(in) :: sx,sy,s,tau, rho, phi, mu
    real(REAL64), dimension(:), intent(in) :: x,y

    real(REAL64), dimension(2) :: r, dr, sig, psi, sc
    real(REAL64) :: sigu, sigv, f
    integer :: n, i

    n = size(y)
    sigu = rho*sin(phi)
    sigv = rho*cos(phi)/tau
    sc = [ sx, sy ]
    sig = [ sqrt(mu/sx + sigv**2), sqrt(mu*tau/sy + sigu**2) ]

    f = 0
    do concurrent (i = 1:n)
       r = [ x(i) - mu, y(i) - mu*tau ] / sig

!          where( abs(r) > 1.345 )
!             r = sign(1.345_REAL64,r)
!          end where
!          if( any(r > 1.34) ) write(*,*) real(r)


       dr = -(1 + r/sig/sc) / sig
       dr(2) = tau*dr(2)
!       psi = tukey(r/s)
!       psi = huber(r/s)
       psi = r/s
       f = f + sum(psi*dr)
    end do
    psifun = -f / n - (tau/sig(2)**2 + 1/sig(1)**2)

  end function psifun

  subroutine robmood(sx,sy,s,y,x,mu,dmu,phi,dphi)

    use medians
    use rfun

    real(REAL64), intent(in) :: sx,sy
    real(REAL64), dimension(:), intent(in) :: x,y
    real(REAL64), intent(in out) :: s, mu, phi, dmu, dphi

    real(REAL64) :: c, stol
    integer :: iter,n
    logical :: converge, reli

    n = size(y)
    allocate(ysel(n),xsel(n))
    ysel = y
    xsel = x
    mumin = 0.9*min(minval(y)/tau,minval(x))
    mumax = 1.1*min(maxval(y)/tau,maxval(x))
    stol = 0.01
    s = scale
    phi = PI / 4
    c = s - scale
    nu = 1
    lambda = nu*c

    converge = .false.
    do iter = 1, 3*precision(tau)

       call themood(mu,dmu,phi,dphi,s,reli)

       if( debug_iter ) then
100       format(a,i0,(2x,g0,a,es7.1),(2x,f4.1,a,f0.1),(2x,a,f0.4),3x, &
               2(2x,f0.3),1x,l1)
          write(error_unit,100) 'mood #',iter,real(mu),'+-',dmu, &
               RAD*phi,'+-',RAD*dphi,'s=',s,lambda,nu,reli
       end if

       c = s - scale
       converge = reli .and. abs(c) < stol
       if( converge ) exit

       lambda = lambda - nu*c
       nu = 2*nu

    end do

    deallocate(xsel, ysel)

  end subroutine robmood

  subroutine themood(mu,dmu,phi,dphi,s,reli)

    use minpacks

    real(REAL64), intent (in out) :: mu,phi,s
    real(REAL64), intent (out) :: dmu,dphi
    logical, intent(out) :: reli

    integer, parameter :: maxfev = 666
!    real(REAL64), parameter :: factor = 100
    integer, parameter :: npar = 3
    real(REAL64), dimension(npar) :: x, diag
    real(REAL64), dimension(npar,npar) :: cov
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: eps
    integer :: info, n, i

    n = size(xsel)
    x = [ mu, phi, s ]
    eps = 10*epsilon(x)*(1 + maxval(abs(x)))

    if( analytic ) then
       label = 'lmder'
       call lmder3(like,x,cov,eps,maxfev,nprint=nprint,info=info)
    else
       label = 'lmdif'
       call lmdif3(liked,x,cov,eps,maxfev,nprint=nprint,info=info)
    end if
    diag = [ (cov(i,i),i=1,npar) ]

    reli = mumin(1) < x(1) .and. x(1) < mumax(1) &
         .and. 0 <= x(2) .and. x(2) <= PI/2 &
         .and. 0.1 < x(3) .and. x(3) < 10 .and. all(diag > 0) &
         .and. (0 < info .and. info < 5)

    if( debug_iter ) then
       msg = ''
       if( 0 < info .and. info < 5 ) then
          msg = ' (success)'
       else if( info == -1 ) then
          msg = ' (out of range)'
       else if( info == 5 ) then
          msg = ' (no convergence)'
       else if( .not. all(diag > 0) ) then
          msg = ' (no minimum)'
       end if
       write(error_unit,'(3a,i0,a,a,l2)') ' themood() ',label,' info: ', &
            info,trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       mu = x(1)
       phi = x(2)
       s = x(3)

       dmu = sqrt(cov(1,1)/(2*n-3))
       dphi = sqrt(cov(2,2)/(2*n-3))
    else
!       if( debug_graph ) then
!          call graph_like_tr('/tmp/like_tr',tau,dtau,rho,phi)
!          call graph_like_tf('/tmp/like_tf',tau,dtau,rho)
!       end if
    end if
!call graph_maxlike(tau,rho,phi)

    if( debug_iter ) then
       block
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64) :: q
         integer :: i,j

         cor = 1
         do i = 1, npar
            do j = i+1, npar
               if( cov(i,i) > 0 .and. cov(j,j) > 0 ) then
                  q = sqrt(cov(i,i)*cov(j,j))
                  cor(i,j) = cov(i,j) / q
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         write(error_unit,'(a,3(es12.4,1x),t48,a,3(f7.4,1x))') &
              'Cov: ', real(cov(1,:)),'corr: ', cor(1,:)
         do i = 2,npar
            write(error_unit,'(5x,3(es12.4,1x),t48,6x,3(f7.4,1x))') &
                 real(cov(i,:)),cor(i,:)
         end do

       end block
    end if



  end subroutine themood


  subroutine allmood(kmean,dkmean,rho,drho,phi,dphi,s,reli)

    use minpacks

    real(REAL64), dimension(:), intent(in out) :: kmean
    real(REAL64), dimension(:), intent(out) :: dkmean
    real(REAL64), intent(in out) :: rho, phi, s
    real(REAL64), intent(out) :: drho, dphi
    logical, intent(out) :: reli

    integer, parameter :: maxfev = 666
    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(:), allocatable :: x, diag
    real(REAL64), dimension(:,:), allocatable :: cov
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: eps
    integer :: info, npar, n, i, K

    K = size(kmean)
    npar = 3 + K
    allocate(x(npar),diag(npar),cov(npar,npar))

    n = size(cts)
!    x = [ mu, phi, s ]
    x(1:K) = kmean
    x(K+1:npar) = [ rho, phi, s ]
    eps = 10*epsilon(x)*(1 + maxval(abs(x)))

    if( analytic ) then
       label = 'lmder'
       call lmder3(likeall,x,cov,eps,maxfev,nprint=nprint,info=info)
    else
       label = 'lmdif'
       call lmdif3(likedall,x,cov,eps,maxfev,nprint=nprint,info=info)
    end if
    diag = [ (cov(i,i),i=1,npar) ]

    reli = all(mumin < x(1:K) .and. x(1:K) < mumax) &
         .and. x(K+1) < rhomax &
         .and. -PI/4 < x(K+2) .and. x(K+2) < 3*PI/4 &
         .and. 0.1 < x(K+3) .and. x(K+3) < 10 &
         .and. all(diag > 0) &
         .and. (0 < info .and. info < 5)

    if( debug_iter ) then
       msg = ''
       if( 0 < info .and. info < 5 ) then
          msg = ' (success)'
       else if( info == -1 ) then
          msg = ' (out of range)'
       else if( info == 5 ) then
          msg = ' (no convergence)'
       else if( .not. all(diag > 0) ) then
          msg = ' (no minimum)'
       end if
       write(error_unit,'(3a,i0,a,a,l2)') ' allmood() ',label,' info: ', &
            info,trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       kmean = x(1:K)
       !       mu = x(1)
       rho = x(K+1)
       phi = x(K+2)
       s = x(K+3)
       if( phi < 0 ) then
          phi = -phi
       else if( phi > PI/2 ) then
          phi = phi
       end if

       if( n > npar ) then
          dkmean = [ (sqrt(cov(i,i)/(n - npar)),i=1,K) ]
          drho = sqrt(cov(K+1,K+1)/(n - npar))
          dphi = sqrt(cov(K+2,K+2)/(n - npar))
       else
          dkmean = -1
          drho = -1
          dphi = -1
       end if
    else
       dkmean = -1
       drho = -1
       dphi = -1
!       if( debug_graph ) then
!          call graph_like_tr('/tmp/like_tr',tau,dtau,rho,phi)
!          call graph_like_tf('/tmp/like_tf',tau,dtau,rho)
!       end if
    end if
!call graph_maxlike(tau,rho,phi)

    if( debug_iter ) then
       block
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64) :: q
         integer :: i,j

         cor = 1
         do i = 1, npar
            do j = i+1, npar
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
         do i = 1,npar
            write(error_unit,'(i2,a,*(es10.3,1x))') i,': ',cov(i,:)
         end do

         write(error_unit,*) 'Correlation:'
         do i = 1,npar
            write(error_unit,'(i2,a,*(f6.3,1x))') i,': ',cor(i,:)
         end do

!         write(error_unit,'(a,*(es12.4,1x),t48,a,*(f7.4,1x))') &
         !              'Cov: ', cov(1,:),'corr: ', cor(1,:)
!         write(error_unit,*) 'Cov, corr:'
!         do i = 1,npar
!            write(error_unit,'(5x,*(es12.4,1x),t48,6x,*(f7.4,1x))') &
!                 cov(i,:),cor(i,:)
!         end do

       end block
    end if



  end subroutine allmood


!!$  subroutine rortho(mju,cts,pht,phi,dphi,s,reli)
!!$
!!$    use minpacks
!!$
!!$    real(REAL64), dimension(:), intent (in) :: mju,cts,pht
!!$    real(REAL64), intent (in out) :: phi,s
!!$    real(REAL64), intent (out) :: dphi
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: maxfev = 666
!!$    real(REAL64), parameter :: factor = 100
!!$    integer, parameter :: npar = 2
!!$    real(REAL64), dimension(npar) :: x, diag
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    character(len=80) :: msg
!!$    character(len=5) :: label
!!$    real(REAL64) :: eps
!!$    integer :: info, n, i
!!$
!!$    n = size(cts) !size(xsel)
!!$    x = [ phi, s ]
!!$    eps = 1e-4
!!$
!!$    if( analytic ) then
!!$       label = 'lmder'
!!$       call lmder3(olike,x,cov,eps,maxfev,factor,nprint,info)
!!$    else
!!$       label = 'lmdif'
!!$       call lmdif3(oliked,x,cov,eps,maxfev,factor,nprint,info)
!!$    end if
!!$    diag = [ (cov(i,i),i=1,npar) ]
!!$
!!$    reli = 0 <= x(1) .and. x(1) <= PI/2 .and. 0.1 < x(2) .and. x(2) < 10 &
!!$         .and. all(diag > 0) .and. (0 < info .and. info < 5)
!!$
!!$    if( debug_iter ) then
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
!!$       write(error_unit,'(3a,i0,a,a,l2)') ' rortho() ',label,' info: ', &
!!$            info,trim(msg),',  reliable:',reli
!!$    end if
!!$
!!$    if( reli ) then
!!$       phi = x(1)
!!$       s = x(2)
!!$
!!$       ! we are incorrect for a single data twice, the estimate holds its order.
!!$!       m = max(n - 2,1)
!!$       dphi = sqrt(cov(1,1)/(n-2))
!!$    else
!!$!       if( debug_graph ) then
!!$!          call graph_like_tr('/tmp/like_tr',tau,dtau,rho,phi)
!!$!          call graph_like_tf('/tmp/like_tf',tau,dtau,rho)
!!$!       end if
!!$    end if
!!$!call graph_maxlike(tau,rho,phi)
!!$
!!$    if( debug_iter ) then
!!$       block
!!$         real(REAL64), dimension(npar,npar) :: cor
!!$         real(REAL64) :: q
!!$         integer :: i,j
!!$
!!$         cor = 1
!!$         do i = 1, npar
!!$            do j = i+1, npar
!!$               if( cov(i,i) > 0 .and. cov(j,j) > 0 ) then
!!$                  q = sqrt(cov(i,i)*cov(j,j))
!!$                  cor(i,j) = cov(i,j) / q
!!$               else
!!$                  cor(i,j) = huge(cor)
!!$               end if
!!$               cor(j,i) = cor(i,j)
!!$            end do
!!$         end do
!!$
!!$         write(error_unit,'(a,2(es12.4,1x),t48,a,2(f7.4,1x))') &
!!$              'Cov: ', real(cov(1,:)),'corr: ', cor(1,:)
!!$         do i = 2,npar
!!$            write(error_unit,'(5x,2(es12.4,1x),t48,6x,2(f7.4,1x))') &
!!$                 real(cov(i,:)),cor(i,:)
!!$         end do
!!$
!!$       end block
!!$    end if
!!$
!!$  end subroutine rortho



!!$  subroutine gscale(cts,pht,mju,sx,sy,tau,rho,g,dg,phi,dphi,reli)
!!$
!!$    use minpacks
!!$
!!$    real(REAL64), dimension(:), intent(in) :: cts, pht, mju
!!$    real(REAL64), intent(in) :: tau,rho,sx,sy
!!$    real(REAL64), intent(out) :: g, phi, dg, dphi
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: maxfev = 666
!!$    real(REAL64), parameter :: factor = 100
!!$    integer, parameter :: npar = 2
!!$    real(REAL64), dimension(npar) :: x
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    character(len=80) :: msg
!!$    character(len=5) :: label
!!$    real(REAL64) :: eps
!!$    integer :: info, n
!!$
!!$    reli = .false.
!!$    n = size(cts)
!!$
!!$    x = [ g, phi ]
!!$    eps = max(10*epsilon(x)*(1 + maxval(abs(x))),1e-3*RAD)
!!$    if( analytic ) then
!!$       label = 'lmder'
!!$       call lmder3(likescale,x,cov,eps,maxfev,factor,nprint,info)
!!$    else
!!$       label = 'lmdif'
!!$       call lmdif3(likescaled,x,cov,eps,maxfev,factor,nprint,info)
!!$    end if
!!$
!!$    reli = gmin < x(1) .and. x(1) < gmax .and. &
!!$        phimin <= x(2) .and. x(2) <= phimax .and. &
!!$        cov(1,1) > 0 .and. cov(2,2) > 0 .and. (0 < info .and. info < 5)
!!$
!!$    if( debug ) then
!!$       msg = ''
!!$       if( 0 < info .and. info < 5 ) then
!!$          msg = ' (success)'
!!$       else if( info == -1 ) then
!!$          msg = ' (out of range)'
!!$       else if( info == 5 ) then
!!$          msg = ' (no convergence)'
!!$       else if( .not. (cov(1,1) > 0 .and. cov(2,2) > 0) ) then
!!$          msg = ' (no minimum)'
!!$       end if
!!$       write(error_unit,'(3a,i0,a,a,l2)') ' gscale() ',label,' info: ', &
!!$            info,trim(msg),',  reliable:',reli
!!$    end if
!!$
!!$    if( reli ) then
!!$       g = x(1)
!!$       phi = abs(x(2))
!!$       ! the function is symetric in rho, we accept both +-rho,
!!$       ! the solutions can be tunneled to mirrored half-space
!!$       ! due to adjust of the ML lambda factor.
!!$
!!$       dg = sqrt(cov(1,1)/(n-2))
!!$       dphi = sqrt(cov(2,2)/(n-2))
!!$    else
!!$       dg = 100
!!$       if( debug_graph ) call graph_likescale('/tmp/likescale',mju(1),dg)
!!$    end if
!!$
!!$    if( debug ) then
!!$       block
!!$         real(REAL64), dimension(npar,npar) :: cor
!!$         real(REAL64) :: q
!!$         integer :: i,j
!!$
!!$         cor = 1
!!$         do i = 1, npar
!!$            do j = i+1, npar
!!$               if( cov(i,i) > 0 .and. cov(j,j) > 0 ) then
!!$                  q = sqrt(cov(i,i)*cov(j,j))
!!$                  cor(i,j) = cov(i,j) / q
!!$               else
!!$                  cor(i,j) = huge(cor)
!!$               end if
!!$               cor(j,i) = cor(i,j)
!!$            end do
!!$         end do
!!$
!!$         write(error_unit,'(a,2(es12.4,1x),t48,a,2(f7.4,1x))') &
!!$              'Cov: ', real(cov(1,:)),'corr: ', cor(1,:)
!!$         do i = 2,npar
!!$            write(error_unit,'(5x,2(es12.4,1x),t48,6x,2(f7.4,1x))') &
!!$                 real(cov(i,:)),cor(i,:)
!!$         end do
!!$
!!$       end block
!!$    end if
!!$
!!$
!!$  end subroutine gscale
!!$

  subroutine like(k,l,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: k, l, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(k), intent(in) :: x
    real(REAL64), dimension(l), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,k), intent(out) :: fjac

!    real(REAL64), dimension(3) :: psi
!    real(REAL64), dimension(3,3) :: dpsi
    real(REAL64) :: mu, rho, phi, s
!    integer :: n

!    n = size(cts)
    mu = x(1)
    rho = x(2)
    phi = x(3)
    s = x(4)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
          write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
          write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
          write(error_unit,*) 'fjac(3,:) = ',real(fjac(3,:))
          write(error_unit,*) 'fjac(4,:) = ',real(fjac(4,:))
       end if

       ! convergence check
       if( mu < mumin(1) .or. mu > mumax(1) .or. phi < 0 .or. phi > PI/2 &
            .or. s < 0 .or. s > 10 ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

       fvec = psisum(mu,ysel,xsel,tau,rho,phi,s,sx,sy,lambda,nu)

    else if( iflag == 2 ) then

       fjac = dpsisum(mu,ysel,xsel,tau,rho,phi,s,sx,sy,nu)

    end if

  end subroutine like

  subroutine liked(k,l,x,fvec,iflag)

    integer, intent(in) :: k, l
    integer, intent(in out) :: iflag
    real(REAL64), dimension(k), intent(in) :: x
    real(REAL64), dimension(l), intent(out) :: fvec
    real(REAL64) :: mu, rho, phi, s

    mu = x(1)
    rho = x(2)
    phi = x(3)
    s = x(4)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
       end if

       ! convergence check ?
       return
    end if

    fvec = psisum(mu,ysel,xsel,tau,rho,phi,s,sx,sy,lambda,nu)

  end subroutine liked

  subroutine likeall(npar,mpar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, mpar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

    real(REAL64), dimension(:), allocatable :: kmean
    integer, dimension(:), allocatable :: idx
    real(REAL64), dimension(4) :: psi
    real(REAL64), dimension(4,4) :: dpsi
    real(REAL64) :: rho, phi, s
    integer :: i,j,k,l,n,m

    n = size(cts)
!    nmju = npar - 2
    kmean = x(1:npar-3)
    rho = x(npar-2)
    phi = x(npar-1)
    s = x(npar)

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
       if( any(kmean < mumin) .or. any(kmean > mumax) &
            .or. abs(rho) > rhomax &
            .or. phi < -PI/4 .or. phi > 3*PI/4 &
            .or. s < 0.1 .or. s > 10 ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

       fvec = 0
       do k = 1, npar - 3
          idx = picker(iset,k)
          m = size(idx)
          psi = psisum(kmean(k),cts(idx),pht(idx),tau,rho,phi,s,sx,sy,lambda,nu)
          fvec(k) = fvec(k) + m*psi(1)
          do l = 0, 2
             fvec(npar-l) = fvec(npar-l) + m*psi(4-l)
          end do
       end do
       fvec = fvec / n

    else if( iflag == 2 ) then

       fjac = 0
       do k = 1, npar - 3
          idx = picker(iset,k)
          m = size(idx)
          dpsi = dpsisum(kmean(k),cts(idx),pht(idx),tau,rho,phi,s,sx,sy,nu)
          fjac(k,k) = m*dpsi(1,1)
          fjac(k,npar-2:npar) = fjac(k,npar-2:npar) + m*dpsi(1,2:4)
          do i = 0, 2
             fjac(npar-i,npar-i:npar) = fjac(npar-i,npar-i:npar) &
                  + m*dpsi(4-i,4-i:4)
          end do
       end do
       fjac = fjac / n

       do concurrent (i = 2:npar)
          do concurrent (j = 1:i-1)
             fjac(i,j) = fjac(j,i)
          end do
       end do

!       do i = 1, npar
!          write(error_unit,'(a,i0,a,*(1x,g0))') &
!               'fjac(',i,',:) = ',real(fjac(i,:)-dpsi(i,:))
!       end do
!       stop 0



    end if

  end subroutine likeall

  subroutine likedall(npar,mpar,x,fvec,iflag)

    integer, intent(in) :: npar, mpar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec

    real(REAL64), dimension(:), allocatable :: kmean
    real(REAL64), dimension(4) :: psi
    integer, dimension(:), allocatable :: idx
    real(REAL64) :: rho, phi, s
    integer :: n, m, k, l

    n = size(cts)
    kmean = x(1:npar-3)
    rho = x(npar-2)
    phi = x(npar-1)
    s = x(npar)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
       end if

       ! convergence check ?
       return
    end if

    fvec = 0
    do k = 1, npar - 3
       idx = picker(iset,k)
       m = size(idx)
       psi = psisum(kmean(k),cts(idx),pht(idx),tau,rho,phi,s,sx,sy,lambda,nu)
       fvec(k) = fvec(k) + m*psi(1)
       do l = 0, 2
          fvec(npar-l) = fvec(npar-l) + m*psi(4-l)
       end do
    end do
    fvec = fvec / n

  end subroutine likedall

!!$  subroutine olike(k,l,x,fvec,fjac,ldfjac,iflag)
!!$
!!$    integer, intent(in) :: k, l, ldfjac
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$    real(REAL64), dimension(ldfjac,k), intent(out) :: fjac
!!$
!!$    real(REAL64), dimension(3) :: psi
!!$    real(REAL64), dimension(3,3) :: dpsi
!!$    real(REAL64), dimension(2) :: sums
!!$    real(REAL64), dimension(2,2) :: dsums
!!$    real(REAL64) :: phi, s, sinfi, cosfi
!!$    integer :: i,n
!!$!    integer :: n
!!$
!!$!    n = size(cts)
!!$!    mu = x(1)
!!$    phi = x(1)
!!$    s = x(2)
!!$
!!$    if( iflag == 0 ) then
!!$       if( debug_iter ) then
!!$          write(error_unit,*) 'x = ',real(x)
!!$          write(error_unit,*) 'fvec = ',real(fvec)
!!$          write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
!!$          write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
!!$!          write(error_unit,*) 'fjac(3,:) = ',real(fjac(3,:))
!!$       end if
!!$
!!$       ! convergence check
!!$       if( phi < 0 .or. phi > PI/2 .or. s < 0.1 .or. s > 10 ) iflag = -1
!!$
!!$       return
!!$    end if
!!$
!!$    n = size(omju)
!!$    cosfi = cos(phi)
!!$    sinfi = sin(phi)
!!$
!!$    if( iflag == 1 ) then
!!$
!!$       sums = 0
!!$       do i = 1,n
!!$          psi = psilike(omju(i),cts(i),pht(i),tau,rho,s,cosfi,sinfi,sx,sy,lambda,nu)
!!$          sums = sums + psi(2:3)
!!$       end do
!!$       fvec = sums / n
!!$
!!$!       psi = psisum(mju,cts,pht,tau,rho,phi,s,sx,sy,lambda,nu)
!!$!       fvec = psi(2:3)
!!$
!!$    else if( iflag == 2 ) then
!!$
!!$!       dpsi = dpsisum(mju,cts,pht,tau,rho,phi,s,sx,sy,nu)
!!$!       fjac = dpsi(2:3,2:3)
!!$
!!$       dsums = 0
!!$       do i = 1,n
!!$          dpsi = dpsilike(omju(i),cts(i),pht(i),tau,rho,s,cosfi,sinfi,sx,sy,nu)
!!$          dsums = dsums + dpsi(2:3,2:3)
!!$       end do
!!$       fjac = dsums / n
!!$
!!$
!!$    end if
!!$
!!$  end subroutine olike
!!$
!!$  subroutine oliked(k,l,x,fvec,iflag)
!!$
!!$    integer, intent(in) :: k, l
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$    real(REAL64), dimension(3) :: psi
!!$    real(REAL64), dimension(2) :: sums
!!$    real(REAL64) :: phi, s, sinfi, cosfi
!!$    integer :: i,n
!!$
!!$!    mu = x(1)
!!$    phi = x(1)
!!$    s = x(2)
!!$
!!$    if( iflag == 0 ) then
!!$       if( debug_iter ) then
!!$          write(error_unit,*) 'x = ',real(x)
!!$          write(error_unit,*) 'fvec = ',real(fvec)
!!$       end if
!!$
!!$       ! convergence check ?
!!$       return
!!$    end if
!!$
!!$!    fvec = psisum(mu,ysel,xsel,tau,rho,phi,s,sx,sy,lambda,nu)
!!$!    psi = psisum(mju,cts,pht,tau,rho,phi,s,sx,sy,lambda,nu)
!!$!    fvec = psi(2:3)
!!$
!!$    n = size(omju)
!!$    cosfi = cos(phi)
!!$    sinfi = sin(phi)
!!$    sums = 0
!!$    do i = 1,n
!!$       psi = psilike(omju(i),cts(i),pht(i),tau,rho,s,cosfi,sinfi,sx,sy,lambda,nu)
!!$       sums = sums + psi(2:3)
!!$    end do
!!$    fvec = sums / n
!!$
!!$
!!$  end subroutine oliked
!!$




  function psisum(mu,y,x,tau,rho,phi,s,sx,sy,lambda,nu) result(sums)

    real(REAL64), dimension(:), intent(in) :: y,x
    real(REAL64), intent(in) :: mu, tau,rho,phi,sx,sy,s,lambda,nu
    real(REAL64), dimension(4) :: sums

    real(REAL64), dimension(4) :: psi
    real(REAL64) :: cosfi, sinfi
    integer :: i, n

    n = size(y)
    cosfi = cos(phi)
    sinfi = sin(phi)

    sums = 0
    !    do concurrent (i = 1:n)
    do i = 1,n
       psi = psilike(mu,y(i),x(i),tau,rho,s,cosfi,sinfi,sx,sy,lambda,nu)
       sums = sums + psi
    end do

    sums = sums / n

  end function psisum

  function dpsisum(mu,y,x,tau,rho,phi,s,sx,sy,nu) result(dsums)

    real(REAL64), dimension(:), intent(in) :: y,x
    real(REAL64), intent(in) :: mu,tau,rho,phi,sx,sy,s,nu
    real(REAL64), dimension(4,4) :: dsums

    real(REAL64), dimension(4,4) :: dpsi
    real(REAL64) :: cosfi, sinfi
    integer :: i, n

    n = size(y)
    cosfi = cos(phi)
    sinfi = sin(phi)

    dsums = 0
!    do concurrent (i = 1:n)
    do i = 1,n
       dpsi = dpsilike(mu,y(i),x(i),tau,rho,s,cosfi,sinfi,sx,sy,nu)
       dsums = dsums + dpsi
    end do

    dsums = dsums / n

  end function dpsisum



  function psilike(mu,y,x,tau,rho,s,cosfi,sinfi,sx,sy,lambda,nu) result(psi)

    use rfun

    real(REAL64), intent(in) :: mu,y,x,tau,rho,s,cosfi,sinfi,sx,sy,lambda,nu
    real(REAL64), dimension(4) :: psi
    real(REAL64) :: sigu, sigv, rx, ry, sigx, sigy, sigx_m, sigx_r, sigx_f, &
         sigy_m, sigy_r, sigy_f, psix, psiy

    sigv = (rho/tau)*cosfi
    sigu = rho*sinfi

    rx = x - mu
    ry = y - mu*tau
    sigx = sqrt(mu + sigv**2)
    sigy = sqrt(mu*tau + sigu**2)
    sigx_m = 1/(2*sigx)
    sigx_r = (rho/sigx)*(cosfi/tau)**2
    sigx_f = -(rho/tau**2)*(rho/sigx)*cosfi*sinfi
    sigy_m = tau/(2*sigy)
    sigy_r = (rho/sigy)*sinfi**2
    sigy_f = rho*(rho/sigy)*cosfi*sinfi

    if( robustness ) then
       if( descent ) then
          psix = tukey((rx / sigx) / s)
          psiy = tukey((ry / sigy) / s)
       else
          psix = huber((rx / sigx) / s)
          psiy = huber((ry / sigy) / s)
       end if
    else
       psix = (rx / sigx) / s
       psiy = (ry / sigy) / s
    end if

    psi(1) = -(psix/s)*(1/sigx + (rx/sigx)*(sigx_m/sigx)) + sigx_m / sigx &
         -(psiy/s)*(tau/sigy + (ry/sigy)*(sigy_m/sigy)) + sigy_m / sigy
    psi(2) = -(psix/s)*(rx/sigx)*(sigx_r/sigx) + sigx_r / sigx &
             -(psiy/s)*(ry/sigy)*(sigy_r/sigy) + sigy_r / sigy
    psi(3) = -(psix/s)*(rx/sigx)*(sigx_f/sigx) + sigx_f / sigx &
             -(psiy/s)*(ry/sigy)*(sigy_f/sigy) + sigy_f / sigy
    psi(4) = -(psix/s**2)*(rx/sigx) - (psiy/s**2)*(ry/sigy) &
         + 2/s - lambda + nu*(s - scale)

  end function psilike

  function dpsilike(mu,y,x,tau,rho,s,cosfi,sinfi,sx,sy,nu) result(dpsi)

    use rfun

    real(REAL64), intent(in) :: mu,y,x,tau,rho,cosfi,sinfi,s,sx,sy,nu
    real(REAL64), dimension(4,4) :: dpsi
    real(REAL64) :: sigu, sigv, rx, ry, sigx, sigy, sigx_m, sigx_r, sigx_f, &
         sigy_m, sigy_r, sigy_f, psix, psiy, sigx_mm, sigx_rr, sigx_ff, &
         sigx_mr, sigx_rf, sigx_mf, sigy_mm, sigy_rr, sigy_ff, sigy_rf, &
         sigy_mr, sigy_mf, dpsix, dpsiy
    integer :: i, j

    sigv = (rho/tau)*cosfi
    sigu = rho*sinfi

    rx = x - mu
    ry = y - mu*tau
    sigx = sqrt(mu + sigv**2)
    sigy = sqrt(mu*tau + sigu**2)

    sigx_m = 1/(2*sigx)
    sigx_r = (rho/sigx)*(cosfi/tau)**2
    sigx_f = -(rho/tau**2)*(rho/sigx)*cosfi*sinfi
    sigy_m = tau/(2*sigy)
    sigy_r = (rho/sigy)*sinfi**2
    sigy_f = rho*(rho/sigy)*cosfi*sinfi

    sigx_mm = -1/(4*sigx**3)
    sigx_rr = (cosfi/tau)**2*(1/sigx - (rho/sigx)*(sigx_r/sigx))
    sigx_ff = -(rho/tau**2)*(rho/sigx)* &
         (cosfi**2 - sinfi**2 - cosfi*sinfi*(sigx_f/sigx))
    sigx_mr = -(cosfi/tau)**2*(rho/sigx)*(sigx_m/sigx)
!    sigx_mr = -(sigx_r/sigx**2)/2
    sigx_rf = -(cosfi/tau**2)*(rho/sigx)*(2*sinfi + cosfi*(sigx_f/sigx))
    sigx_mf = -(sigx_f/sigx)*sigx_m

    sigy_mm = -tau**2/(4*sigy**3)
    sigy_rr = sinfi**2*(1/sigy - (rho/sigy)*(sigy_r/sigy))
    sigy_ff = rho*(rho/sigy)*(cosfi**2 - sinfi**2 - cosfi*sinfi*(sigy_f/sigy))
    sigy_mr = -sinfi**2*(rho/sigy)*(sigy_m/sigy)
!    sigy_mr = -tau*(sigy_r/sigy**2)/2
    sigy_rf = (rho/sigy)*(2*cosfi - sinfi*(sigy_f/sigy))*sinfi
    sigy_mf = -(sigy_f/sigy)*sigy_m

    if( robustness ) then
       if( descent ) then
          psix = tukey((rx / sigx) / s)
          psiy = tukey((ry / sigy) / s)
          dpsix = dtukey((rx / sigx) / s)
          dpsiy = dtukey((ry / sigy) / s)
       else
          psix = huber((rx / sigx) / s)
          psiy = huber((ry / sigy) / s)
          dpsix = dhuber((rx / sigx) / s)
          dpsiy = dhuber((ry / sigy) / s)
       end if
    else
       psix = (rx / sigx) / s
       psiy = (ry / sigy) / s
       dpsix = 1
       dpsiy = 1
    end if


    dpsi(1,1) = &
         (dpsix/s**2)*(1/sigx + (rx/sigx)*(sigx_m/sigx))**2 &
        +(dpsiy/s**2)*(tau/sigy + (ry/sigy)*(sigy_m/sigy))**2 &
         -(psix/s)*(-2*(sigx_m/sigx)/sigx + (rx/sigx)*(sigx_mm/sigx) &
                    -2*(rx/sigx)*(sigx_m/sigx)**2) &
         -(psiy/s)*(-2*tau*(sigy_m/sigy)/sigy + (ry/sigy)*(sigy_mm/sigy)&
                    -2*(ry/sigy)*(sigy_m/sigy)**2) &
         -(sigx_m/sigx)**2 + sigx_mm/sigx -(sigy_m/sigy)**2 + sigy_mm/sigy
    dpsi(2,2) = &
         (dpsix/s**2)*((rx/sigx)*(sigx_r/sigx))**2 &
        +(dpsiy/s**2)*((ry/sigy)*(sigy_r/sigy))**2 &
        +(psix/s)*(rx/sigx)*(2*(sigx_r/sigx)**2 - sigx_rr/sigx) &
        +(psiy/s)*(ry/sigy)*(2*(sigy_r/sigy)**2 - sigy_rr/sigy) &
        + sigx_rr/sigx - (sigx_r/sigx)**2 + sigy_rr/sigy - (sigy_r/sigy)**2
    dpsi(3,3) = &
         (dpsix/s**2)*((rx/sigx)*(sigx_f/sigx))**2 &
        +(dpsiy/s**2)*((ry/sigy)*(sigy_f/sigy))**2 &
        -(psix/s)*(rx/sigx)*((sigx_ff/sigx) - 2*(sigx_f/sigx)**2) &
        -(psiy/s)*(ry/sigy)*((sigy_ff/sigy) - 2*(sigy_f/sigy)**2) &
        -(sigx_f/sigx)**2 + sigx_ff/sigx -(sigy_f/sigy)**2 + sigy_ff/sigy
    dpsi(4,4) = &
        +dpsix*(rx/sigx)**2/s**4 + 2*(psix/s**3)*(rx/sigx) &
        +dpsiy*(ry/sigy)**2/s**4 + 2*(psiy/s**3)*(ry/sigy) &
        -2/s**2 + nu

    dpsi(1,2) = &
       (dpsix/s**2)*(rx/sigx)*(sigx_r/sigx)*(1/sigx + (rx/sigx)*(sigx_m/sigx))&
      +(dpsiy/s**2)*(ry/sigy)*(sigy_r/sigy)*(tau/sigy + (ry/sigy)*(sigy_m/sigy))&
      +(psix/s)*(sigx_r/sigx**2 - (rx/sigx)*(sigx_mr/sigx) &
      +2*(rx/sigx)*(sigx_m/sigx)*(sigx_r/sigx)) &
      +(psiy/s)*(tau*sigy_r/sigy**2 - (ry/sigy)*(sigy_mr/sigy) &
      +2*(ry/sigy)*(sigy_m/sigy)*(sigy_r/sigy)) &
      -(sigx_r/sigx)*(sigx_m/sigx) + sigx_mr/sigx &
      -(sigy_r/sigy)*(sigy_m/sigy) + sigy_mr/sigy


    dpsi(1,3) = &
       (dpsix/s**2)*(rx/sigx)*(sigx_f/sigx)*(1/sigx + (rx/sigx)*(sigx_m/sigx))&
      -(psix/s)*(-sigx_f/sigx**2 + (rx/sigx)*(sigx_mf/sigx) &
                 -2*(rx/sigx)*(sigx_m/sigx)*(sigx_f/sigx)) &
     +(dpsiy/s**2)*(tau/sigy + (ry/sigy)*(sigy_m/sigy))*(ry/sigy)*(sigy_f/sigy)&
      -(psiy/s)*(-tau*sigy_f/sigy**2 + (ry/sigy)*(sigy_mf/sigy) &
                 -2*(ry/sigy)*(sigy_m/sigy)*(sigy_f/sigy)) &
         -(sigx_f/sigx)*(sigx_m/sigx) + sigx_mf/sigx &
         -(sigy_f/sigy)*(sigy_m/sigy) + sigy_mf/sigy

    dpsi(1,4) = &
         (dpsix/s**3)*(rx/sigx)*(1/sigx + (rx/sigx)*(sigx_m/sigx)) &
         +(psix/s**2)*(1/sigx + (rx/sigx)*(sigx_m/sigx)) &
        +(dpsiy/s**3)*(ry/sigy)*(tau/sigy + (ry/sigy)*(sigy_m/sigy)) &
        +(psiy/s**2)*(tau/sigy + (ry/sigy)*(sigy_m/sigy))
    dpsi(2,3) = &
         (dpsix/s**2)*(rx/sigx**2)**2*(sigx_r*sigx_f) &
        +(dpsiy/s**2)*(ry/sigy**2)**2*(sigy_r*sigy_f) &
        +(psix/s)*(rx/sigx)*(2*(sigx_r/sigx)*(sigx_f/sigx) - sigx_rf/sigx) &
        +(psiy/s)*(ry/sigy)*(2*(sigy_r/sigy)*(sigy_f/sigy) - sigy_rf/sigy) &
        +sigx_rf/sigx - (sigx_r/sigx)*(sigx_f/sigx) &
        +sigy_rf/sigy - (sigy_r/sigy)*(sigy_f/sigy)

    dpsi(2,4) = &
         (dpsix/s**3)*(rx/sigx)**2*(sigx_r/sigx) &
        +(dpsiy/s**3)*(ry/sigy)**2*(sigy_r/sigy) &
        +(psix/s**2)*(rx/sigx)*(sigx_r/sigx) &
        +(psiy/s**2)*(ry/sigy)*(sigy_r/sigy)

    dpsi(3,4) = &
         (dpsix/s**3)*(rx/sigx)**2*(sigx_f/sigx) &
        + (psix/s**2)*(rx/sigx)*(sigx_f/sigx) &
        +(dpsiy/s**3)*(ry/sigy)**2*(sigy_f/sigy) &
        + (psiy/s**2)*(ry/sigy)*(sigy_f/sigy)

    do concurrent(i = 2:size(dpsi,1))
       do concurrent(j = 1:i-1)
          dpsi(i,j) = dpsi(j,i)
       end do
    end do
!    dpsi(2,1) = dpsi(1,2)
!    dpsi(3,1) = dpsi(1,3)
!    dpsi(3,2) = dpsi(2,3)

  end function dpsilike


  real(REAL64) function eloglike(mu,y,x,tau,rho,phi,s,lambda,nu)

    use rfun

    real(REAL64), dimension(:), intent(in) :: y,x
    real(REAL64), intent(in) :: mu,tau,rho,phi,s,lambda,nu

    real(REAL64) :: ss,rx,ry,sigx,sigy
    integer :: n,i

    n = size(y)
    sigx = sqrt(mu + ((rho/tau)*cos(phi))**2)
    sigy = sqrt(mu*tau + (rho*sin(phi))**2)

    ss = 0
    do concurrent (i = 1:n)
       rx = x(i) - mu
       ry = y(i) - mu*tau
       if( robustness ) then
          if( descent ) then
             ss = ss + itukey(rx/sigx/s) + itukey(ry/sigy/s)
          else
             ss = ss + ihuber((rx/sigx)/s) + ihuber((ry/sigy)/s)
          end if
       else
          ss = ss + ((rx/sigx)/s)**2/2 + ((ry/sigy)/s)**2/2
       end if
    end do

    eloglike = ss / n + log(sigx) + log(sigy) + 2*log(s) &
         - lambda*(s - scale) + (nu/2)*(s - scale)**2

  end function eloglike

  real(REAL64) function eloglikes(mju,y,x,tau,sigu,sigv,s,lambda,nu)

    use rfun

    real(REAL64), dimension(:), intent(in) :: mju,y,x
    real(REAL64), intent(in) :: tau,sigv,sigu,s,lambda,nu

    real(REAL64) :: ss,rx,ry,sigx,sigy
    integer :: n,i

    n = size(y)
    ss = 0
    do concurrent (i = 1:n)
       rx = x(i) - mju(i)
       ry = y(i) - mju(i)*tau
       sigx = sqrt(mju(i) + sigv**2)
       sigy = sqrt(mju(i)*tau + sigu**2)
       if( robustness ) then
          if( descent ) then
             ss = ss + itukey(rx/sigx/s) + itukey(ry/sigy/s)
          else
             ss = ss + ihuber((rx/sigx)/s) + ihuber((ry/sigy)/s)
          end if
       else
          ss = ss + ((rx/sigx)/s)**2/2 + ((ry/sigy)/s)**2/2
       end if
    end do

    eloglikes = ss / n + log(sigx) + log(sigy) + 2*log(s) &
         - lambda*(s - scale) + (nu/2)*(s - scale)**2

  end function eloglikes

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
       sigx = sqrt(mju(i) + sigv**2)
       if( robustness ) then
          if( descent ) then
             ss = ss + itukey(rx/sigx/s) + log(sigx)
          else
             ss = ss + ihuber((rx/sigx)/s) + log(sigx)
          end if
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
       if( robustness ) then
          if( descent ) then
             ss = ss + itukey(ry/sigy/s) + log(sigy)
          else
             ss = ss + ihuber((ry/sigy)/s) + log(sigy)
          end if
       else
          ss = ss + ((ry/sigy)/s)**2/2 + log(sigy)
       end if
    end do

    yeloglikes = ss / n + log(s) - lambda*(s - scale) + (nu/2)*(s - scale)**2

  end function yeloglikes


  subroutine rmoodtest(y,x,mu,sx,sy,tau,rho,phi,s,lambda,nu)

    real(REAL64), dimension(:), intent(in) :: y,x
    real(REAL64), intent(in) :: mu,sx,sy,tau,rho,phi,s,lambda,nu

    real(REAL64), dimension(4) :: fd,fvec,eps,z,z1,z2
    real(REAL64), dimension(4,4) :: fjac, djac
    real(REAL64), dimension(-1:1,4) :: dg
    real(REAL64), dimension(-1:1) :: f
    integer :: i, j, l

    eps = epsilon(eps)**(1.0/3.0)
    fvec = psisum(mu,y,x,tau,rho,phi,s,sx,sy,lambda,nu)
    do i = 1, size(fvec)
       do l = -1,1,2
          z = [ mu, rho, phi, s ]
          z(i) = z(i) + l*eps(i)
          f(l) = eloglike(z(1),y,x,tau,z(2),z(3),z(4),lambda,nu)
       end do
       fd(i) = (f(1) - f(-1)) / (2*eps(i))
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
            real(fd(i)-fvec(i)),real(eps(i))
    end do

    fjac = dpsisum(mu,y,x,tau,rho,phi,s,sx,sy,nu)
    do i = 1,size(fjac,1)
       do j = i, size(fjac,2)
          do l = -1,1,2
             z = [ mu, rho, phi, s ]
             z(i) = z(i) + l*eps(i)
             dg(l,:) = psisum(z(1),y,x,tau,z(2),z(3),z(4),sx,sy,lambda,nu)
          end do
          djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps(i))

          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do

  end subroutine rmoodtest


!!$  real(REAL64) function eloglikes(mju,cts,pht,tau,rho,phi,sx,sy)
!!$
!!$    real(REAL64), dimension(:), intent(in) :: mju,cts,pht
!!$    real(REAL64), intent(in) :: tau,rho,phi,sx,sy
!!$
!!$    real(REAL64) :: rx,ry,sigu,sigv,sigx,sigy,s
!!$    integer :: n,i
!!$
!!$    n = size(cts)
!!$    sigu = rho*sin(phi)
!!$    sigv = (rho/tau)*cos(phi)
!!$
!!$    s = 0
!!$    do concurrent (i = 1:n)
!!$       s = s + eloglike(mju(i),cts,pht,tau,rho,phi)
!!$    end do
!!$    eloglikes = s
!!$
!!$  end function eloglikes



  subroutine xmood(mju,sigv,dsigv,s,reli)

    use minpacks

    real(REAL64), dimension(:), intent(in) :: mju
!    real(REAL64), dimension(:), intent(out) :: dkmean
    real(REAL64), intent(in out) :: sigv, s
    real(REAL64), intent(out) :: dsigv
    logical, intent(out) :: reli

    integer, parameter :: npar = 2
    integer, parameter :: maxfev = 666
!    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(npar) :: x, diag
    real(REAL64), dimension(npar,npar) :: cov
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: eps
    integer :: info, n, i, K

!    K = size(kmean)
!    npar = 3 + K
!    allocate(x(npar),diag(npar),cov(npar,npar))

    n = size(cts)
    x = [ sigv, s ]
!    x(1:K) = kmean
!    x(K+1:npar) = [ rho, phi, s ]
    eps = 10*epsilon(x)*(1 + abs(sigv))

    if( analytic ) then
       label = 'lmder'
       call lmder3(xlikesig,x,cov,eps,maxfev,nprint=nprint,info=info)
    else
       label = 'lmdif'
       call lmdif3(xlikedsig,x,cov,eps,maxfev,nprint=nprint,info=info)
    end if
    diag = [ (cov(i,i),i=1,npar) ]

    reli = abs(x(1)) < rhomax .and. 0.1 < x(2) .and. x(2) < 10 &
         .and. all(diag > 0) .and. (0 < info .and. info < 5)

    if( debug_iter ) then
       msg = ''
       if( 0 < info .and. info < 5 ) then
          msg = ' (success)'
       else if( info == -1 ) then
          msg = ' (out of range)'
       else if( info == 5 ) then
          msg = ' (no convergence)'
       else if( .not. all(diag > 0) ) then
          msg = ' (no minimum)'
       end if
       write(error_unit,'(3a,i0,a,a,l2)') ' xmood() ',label,' info: ', &
            info,trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       ! Sometimes, the lambda parameter of LM-method is incidentally found
       ! in the mirrored part of the parameter space. Fortunately, the function
       ! is pretty symmetric, so I do to ignore the "tunnelling". The correct,
       ! in principle, solution should to keep the condition: sigv > 0.
       sigv = abs(x(1))
       s = x(2)

       dsigv = sqrt(cov(1,1)/(n - npar))
    else
       dsigv = -1
    end if

    if( debug_iter ) then
       block
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64) :: q
         integer :: i,j

         cor = 1
         do i = 1, npar
            do j = i+1, npar
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
         do i = 1,npar
            write(error_unit,'(i2,a,*(es10.3,1x))') i,': ',cov(i,:)
         end do

         write(error_unit,*) 'Correlation:'
         do i = 1,npar
            write(error_unit,'(i2,a,*(f6.3,1x))') i,': ',cor(i,:)
         end do

       end block
    end if

  end subroutine xmood

  function xpsilike(mu,x,tau,sigv,s,lambda,nu) result(psi)

    use rfun

    real(REAL64), intent(in) :: mu,x,tau,sigv,s,lambda,nu
    real(REAL64), dimension(2) :: psi
    real(REAL64) :: rx, sigx, psix, sigx_v

    rx = x - mu
    sigx = sqrt(mu + sigv**2)
    sigx_v = sigv/sigx

    if( robustness ) then
       if( descent ) then
          psix = tukey((rx / sigx) / s)
       else
          psix = huber((rx / sigx) / s)
       end if
    else
       psix = (rx / sigx) / s
    end if

    psi(1) = (1 - (psix/s)*(rx/sigx))*(sigx_v/sigx)
    psi(2) = -(psix/s**2)*(rx/sigx) + 1/s - lambda + nu*(s - scale)

  end function xpsilike


  function xdpsilike(mu,x,tau,sigv,s,nu) result(dpsi)

    use rfun

    real(REAL64), intent(in) :: mu,x,tau,sigv,s,nu
    real(REAL64), dimension(2,2) :: dpsi
    real(REAL64) :: rx, sigx, psix, dpsix, sigx_v, sigx_vv

    rx = x - mu
    sigx = sqrt(mu + sigv**2)
    sigx_v = sigv/sigx
    sigx_vv = -(sigv/sigx)*(sigx_v/sigx) + 1/sigx

    if( robustness ) then
       if( descent ) then
          psix = tukey((rx / sigx) / s)
          dpsix = dtukey((rx / sigx) / s)
       else
          psix = huber((rx / sigx) / s)
          dpsix = dhuber((rx / sigx) / s)
       end if
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
         - 1/s**2 + nu

    dpsi(2,1) = dpsi(1,2)

  end function xdpsilike


  function xpsisum(mju,x,tau,sigv,s,lambda,nu) result(sums)

    real(REAL64), dimension(:), intent(in) :: mju,x
    real(REAL64), intent(in) :: tau,sigv,s,lambda,nu
    real(REAL64), dimension(2) :: sums

    real(REAL64), dimension(2) :: psi
    integer :: i, n

    n = size(x)

    sums = 0
    !    do concurrent (i = 1:n)
    do i = 1,n
       psi = xpsilike(mju(i),x(i),tau,sigv,s,lambda,nu)
       sums = sums + psi
    end do
    sums = sums / n

  end function xpsisum

  function xdpsisum(mju,x,tau,sigv,s,nu) result(dsums)

    real(REAL64), dimension(:), intent(in) :: mju,x
    real(REAL64), intent(in) :: tau,sigv,s,nu
    real(REAL64), dimension(2,2) :: dsums

    real(REAL64), dimension(2,2) :: dpsi
    integer :: i, n

    n = size(x)

    dsums = 0
!    do concurrent (i = 1:n)
    do i = 1,n
       dpsi = xdpsilike(mju(i),x(i),tau,sigv,s,nu)
       dsums = dsums + dpsi
    end do
    dsums = dsums / n

  end function xdpsisum


  subroutine xmoodtest(x,mju,tau,rho,phi,s,lambda,nu)

    real(REAL64), dimension(:), intent(in) :: x,mju
    real(REAL64), intent(in) :: tau,rho,phi,s,lambda,nu

    real(REAL64), dimension(2) :: fd,fvec,eps,z,z1,z2
    real(REAL64), dimension(2,2) :: fjac, djac
    real(REAL64), dimension(-1:1,2) :: dg
    real(REAL64), dimension(-1:1) :: f
    real(REAL64) :: sigv, sigu
    integer :: i, j, l

    sigu = rho*sin(phi)
    sigv = rho*cos(phi)/tau

    eps = epsilon(eps)**(1.0/3.0)
    fvec = xpsisum(mju,x,tau,sigv,s,lambda,nu)
    do i = 1, size(fvec)
       do l = -1,1,2
          z = [ sigv, s ]
          z(i) = z(i) + l*eps(i)
          f(l) = xeloglikes(mju,x,tau,z(1),z(2),lambda,nu)
       end do
       fd(i) = (f(1) - f(-1)) / (2*eps(i))
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
            real(fd(i)-fvec(i)),real(eps(i))
    end do

    fjac = xdpsisum(mju,x,tau,sigv,s,nu)
    do i = 1,size(fjac,1)
       do j = i, size(fjac,2)
          do l = -1,1,2
             z = [ sigv, s ]
             z(i) = z(i) + l*eps(i)
             dg(l,:) = xpsisum(mju,x,tau,z(1),z(2),lambda,nu)
          end do
          djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps(i))

          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do

  end subroutine xmoodtest

  subroutine xlikesig(npar,mpar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, mpar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

    real(REAL64), dimension(2) :: psi
    real(REAL64), dimension(2,2) :: dpsi
    real(REAL64) :: sigv, s
    integer :: i,j,n

    n = size(cts)
    sigv = x(1)
    s = x(2)

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

       fvec = xpsisum(mju,pht,tau,sigv,s,lambda,nu)

    else if( iflag == 2 ) then

       fjac = xdpsisum(mju,pht,tau,sigv,s,nu)

       do concurrent (i = 2:npar)
          do concurrent (j = 1:i-1)
             fjac(i,j) = fjac(j,i)
          end do
       end do

    end if

  end subroutine xlikesig

  subroutine xlikedsig(npar,mpar,x,fvec,iflag)

    integer, intent(in) :: npar, mpar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec

    real(REAL64), dimension(2) :: psi
    real(REAL64) :: sigv, s

    sigv = x(1)
    s = x(2)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
       end if

       ! convergence check ?
       return
    end if

    fvec = xpsisum(mju,pht,tau,sigv,s,lambda,nu)

  end subroutine xlikedsig

  real(REAL64) function xmoodmin(sigv)

    use medians
    use quantilmean

    real(REAL64), intent(in) :: sigv

    real(REAL64), dimension(:), allocatable :: xcdf, ycdf
    real(REAL64), dimension(2), parameter :: q = [1, 3] / real(4.0,REAL64)
    real(REAL64), dimension(2) :: x
    character(len=80) :: errmsg
    real(REAL64) :: xmad
    integer :: i, n, stat

    xmad = median(abs(pht - mju)/sqrt(mju + sigv**2))
    xmoodmin = abs(xmad/Qdev - 1)

    n = size(pht)
    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,'(2a)') 'Error in xmoodmin: ',trim(errmsg)
       return
    end if

    ! empirical CDF
    call ecdf((pht - mju)/sqrt(mju + sigv**2),xcdf,ycdf)

    do i = 1, 2
       x(i) = quantile(q(i),xcdf,ycdf)
    end do
    xmad = (x(2) - x(1)) / 2
    xmoodmin = abs(xmad/Qdev - 1)

  end function xmoodmin


  subroutine ymood(mju,sigu,dsigu,s,reli)

    use minpacks

    real(REAL64), dimension(:), intent(in) :: mju
    real(REAL64), intent(in out) :: sigu, s
    real(REAL64), intent(out) :: dsigu
    logical, intent(out) :: reli

    integer, parameter :: npar = 2
    integer, parameter :: maxfev = 666
!    real(REAL64), parameter :: factor = 100
    real(REAL64), dimension(npar) :: par, diag
    real(REAL64), dimension(npar,npar) :: cov
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: eps
    integer :: info, n, i, K

    n = size(cts)
    par = [ sigu, s ]
    eps = 10*epsilon(par)*(1 + abs(sigu))

    if( analytic ) then
       label = 'lmder'
       call lmder3(ylikesig,par,cov,eps,maxfev,nprint=nprint,info=info)
    else
       label = 'lmdif'
       call lmdif3(ylikedsig,par,cov,eps,maxfev,nprint=nprint,info=info)
    end if
    diag = [ (cov(i,i),i=1,npar) ]

    reli = abs(par(1)) < rhomax .and. 0.1 < par(2) .and. par(2) < 10 &
         .and. all(diag > 0) .and. (0 < info .and. info < 5)

    if( debug_iter ) then
       msg = ''
       if( 0 < info .and. info < 5 ) then
          msg = ' (success)'
       else if( info == -1 ) then
          msg = ' (out of range)'
       else if( info == 5 ) then
          msg = ' (no convergence)'
       else if( .not. all(diag > 0) ) then
          msg = ' (no minimum)'
       end if
       write(error_unit,'(3a,i0,a,a,l2)') ' ymood() ',label,' info: ', &
            info,trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       sigu = abs(par(1))
       s = par(2)

       dsigu = sqrt(cov(1,1)/(n - npar))
    else
       dsigu = -1
    end if

    if( debug_iter ) then
       block
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64) :: q
         integer :: i,j

         cor = 1
         do i = 1, npar
            do j = i+1, npar
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
         do i = 1,npar
            write(error_unit,'(i2,a,*(es10.3,1x))') i,': ',cov(i,:)
         end do

         write(error_unit,*) 'Correlation:'
         do i = 1,npar
            write(error_unit,'(i2,a,*(f6.3,1x))') i,': ',cor(i,:)
         end do

       end block
    end if

  end subroutine ymood


  function ypsilike(mu,y,tau,sigu,s,lambda,nu) result(psi)

    use rfun

    real(REAL64), intent(in) :: mu,y,tau,sigu,s,lambda,nu
    real(REAL64), dimension(2) :: psi
    real(REAL64) :: ry, sigy, psiy, sigy_u

    ry = y - mu*tau
    sigy = sqrt(mu*tau + sigu**2)
    sigy_u = sigu/sigy

    if( robustness ) then
       if( descent ) then
          psiy = tukey((ry / sigy) / s)
       else
          psiy = huber((ry / sigy) / s)
       end if
    else
       psiy = (ry / sigy) / s
    end if

    psi(1) = (1 - (psiy/s)*(ry/sigy))*(sigy_u/sigy)
    psi(2) = -(psiy/s**2)*(ry/sigy) + 1/s - lambda + nu*(s - scale)

  end function ypsilike


  function ydpsilike(mu,y,tau,sigu,s,nu) result(dpsi)

    use rfun

    real(REAL64), intent(in) :: mu,y,tau,sigu,s,nu
    real(REAL64), dimension(2,2) :: dpsi
    real(REAL64) :: ry, sigy, psiy, dpsiy, sigy_u, sigy_uu

    ry = y - mu*tau
    sigy = sqrt(mu*tau + sigu**2)
    sigy_u = sigu/sigy
    sigy_uu = 1/sigy - (sigu/sigy)*(sigy_u/sigy)

    if( robustness ) then
       if( descent ) then
          psiy = tukey((ry / sigy) / s)
          dpsiy = dtukey((ry / sigy) / s)
       else
          psiy = huber((ry / sigy) / s)
          dpsiy = dhuber((ry / sigy) / s)
       end if
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
         - 1/s**2 + nu

    dpsi(2,1) = dpsi(1,2)

  end function ydpsilike


  function ypsisum(mju,y,tau,sigu,s,lambda,nu) result(sums)

    real(REAL64), dimension(:), intent(in) :: mju,y
    real(REAL64), intent(in) :: tau,sigu,s,lambda,nu
    real(REAL64), dimension(2) :: sums

    real(REAL64), dimension(2) :: psi
    integer :: i, n

    n = size(y)

    sums = 0
    !    do concurrent (i = 1:n)
    do i = 1,n
       psi = ypsilike(mju(i),y(i),tau,sigu,s,lambda,nu)
       sums = sums + psi
    end do
    sums = sums / n

  end function ypsisum

  function ydpsisum(mju,y,tau,sigu,s,nu) result(dsums)

    real(REAL64), dimension(:), intent(in) :: mju,y
    real(REAL64), intent(in) :: tau,sigu,s,nu
    real(REAL64), dimension(2,2) :: dsums

    real(REAL64), dimension(2,2) :: dpsi
    integer :: i, n

    n = size(y)

    dsums = 0
!    do concurrent (i = 1:n)
    do i = 1,n
       dpsi = ydpsilike(mju(i),y(i),tau,sigu,s,nu)
       dsums = dsums + dpsi
    end do
    dsums = dsums / n

  end function ydpsisum


  subroutine ylikesig(npar,mpar,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: npar, mpar, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec
    real(REAL64), dimension(ldfjac,npar), intent(out) :: fjac

    real(REAL64), dimension(2) :: psi
    real(REAL64), dimension(2,2) :: dpsi
    real(REAL64) :: sigu, s
    integer :: i,j

!    n = size(cts)
    sigu = x(1)
    s = x(2)

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
       if( abs(sigu) > rhomax .or. s < 0.1 .or. s > 10 ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

       fvec = ypsisum(mju,cts,tau,sigu,s,lambda,nu)

    else if( iflag == 2 ) then

       fjac = ydpsisum(mju,cts,tau,sigu,s,nu)

!       do concurrent (i = 2:npar)
!          do concurrent (j = 1:i-1)
!             fjac(i,j) = fjac(j,i)
!          end do
!       end do

    end if

  end subroutine ylikesig

  subroutine ylikedsig(npar,mpar,x,fvec,iflag)

    integer, intent(in) :: npar, mpar
    integer, intent(in out) :: iflag
    real(REAL64), dimension(npar), intent(in) :: x
    real(REAL64), dimension(mpar), intent(out) :: fvec

    real(REAL64), dimension(2) :: psi
    real(REAL64) :: sigu, s

    sigu = x(1)
    s = x(2)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
       end if

       ! convergence check ?
       return
    end if

    fvec = ypsisum(mju,cts,tau,sigu,s,lambda,nu)

  end subroutine ylikedsig


  real(REAL64) function ymoodmin(sigu)

    use medians
    use quantilmean

    real(REAL64), intent(in) :: sigu

    real(REAL64), dimension(:), allocatable :: xcdf, ycdf
    real(REAL64), dimension(2), parameter :: q = [1, 3] / real(4.0,REAL64)
    real(REAL64), dimension(2) :: y
    character(len=80) :: errmsg
    real(REAL64) :: ymad
    integer :: i, n, stat


    ymad = median(abs(cts - mju*tau)/sqrt(mju*tau + sigu**2))
    ymoodmin = abs(ymad/Qdev - 1)

    n = size(cts)
    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,'(2a)') 'Error in ymoodmin: ',trim(errmsg)
       return
    end if

    ! empirical CDF
    call ecdf((cts - mju*tau)/sqrt(mju*tau + sigu**2),xcdf,ycdf)

    do i = 1, 2
       y(i) = quantile(q(i),xcdf,ycdf)
    end do
    ymad = (y(2) - y(1)) / 2
    ymoodmin = abs(ymad/Qdev - 1)
    write(error_unit,*) sigu,ymoodmin

  end function ymoodmin




  subroutine ymoodtest(y,mju,tau,rho,phi,s,lambda,nu)

    real(REAL64), dimension(:), intent(in) :: y,mju
    real(REAL64), intent(in) :: tau,rho,phi,s,lambda,nu

    real(REAL64), dimension(2) :: fd,fvec,eps,z,z1,z2
    real(REAL64), dimension(2,2) :: fjac, djac
    real(REAL64), dimension(-1:1,2) :: dg
    real(REAL64), dimension(-1:1) :: f
    real(REAL64) :: sigu
    integer :: i, j, l

    sigu = rho*sin(phi)

    eps = epsilon(eps)**(1.0/3.0)
    fvec = ypsisum(mju,y,tau,sigu,s,lambda,nu)
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

    fjac = ydpsisum(mju,y,tau,sigu,s,nu)
    do i = 1,size(fjac,1)
       do j = i, size(fjac,2)
          do l = -1,1,2
             z = [ sigu, s ]
             z(i) = z(i) + l*eps(i)
             dg(l,:) = ypsisum(mju,y,tau,z(1),z(2),lambda,nu)
          end do
          djac(i,:) = (dg(1,:) - dg(-1,:)) / (2*eps(i))

          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))
       end do
    end do

  end subroutine ymoodtest


  subroutine graph_mood(tau,rho,phi,sigu,sigv,cts,pht,mju)

    use quantilmean

    real(REAL64), intent(in) :: tau,rho,phi,sigu,sigv
    real(REAL64), dimension(:), intent(in) :: cts, pht, mju

    real(REAL64), dimension(:), allocatable :: xpdf, ypdf, res
    character(len=80) :: errmsg
!    real(REAL64) :: sigu, sigv
    integer :: n,i,stat

    n = size(cts)
    allocate(xpdf(n),ypdf(n),res(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,*) 'Error: `',trim(errmsg),"'(graph_mood)"
       return
    end if

    ! empirical CDF
!    sigu = rho*sin(phi)
!    sigv = rho*cos(phi)/tau

    write(*,*) 'sigu, sigv:',sigu, sigv

    ! pht
    res = (pht - mju) / sqrt(mju + sigv**2)

    call ecdf(res,xpdf,ypdf)

    open(1,file='/tmp/xcdf')
    do i = 1, size(xpdf)
       write(1,*) xpdf(i), ypdf(i)
    end do
    close(1)

    res = (cts - tau*mju) / sqrt(mju*tau + sigu**2)
!    res = (cts - tau*mju) / sqrt(mju*(tau + 0*tau**2/n) + 0*mju**2*dtau**2 + sigu**2)

    call ecdf(res,xpdf,ypdf)

    open(1,file='/tmp/ycdf')
    do i = 1, size(xpdf)
       write(1,*) xpdf(i), ypdf(i)
    end do
    close(1)

!    res = ((cts - tau*mju)**2 + tau**2*(pht - mju)**2) &
!         / sqrt(mju*tau*(1 + tau) + rho**2)**2
!    call ecdf(res,xpdf,ypdf)
!    open(1,file='/tmp/rcdf')
!    do i = 1, size(xpdf)
!       write(1,*) xpdf(i), ypdf(i)
!    end do
!    close(1)



  end subroutine graph_mood


  subroutine graph_moodline(tau,rho,phi,sx,sy,cts,pht,mju,sigu,sigv)

    use medians
    use quantilmean

    real(REAL64), intent(in) :: tau,rho,phi,sx,sy, sigu, sigv
    real(REAL64), dimension(:), intent(in) :: cts, pht, mju

    real(REAL64), dimension(:), allocatable :: xpdf, ypdf, res
    character(len=80) :: errmsg
    real(REAL64) :: t,s,tmin, tmax
    integer :: i

    open(1,file='/tmp/ml')

    do i = 1, size(cts)

       t = pht(i)*tau
       s = sqrt(pht(i) + rho**2/2)
       tmin = t - s
       tmax = t + s
       write(1,*) pht(i),cts(i),t,tmin,tmax,cts(i)-t,s,-s,mju(i),tau*mju(i)

    end do
    close(1)

    ! gnuplot> plot 'ml', 'ml' u 1:3 w l, 'ml' u 1:4 w l, 'ml' u 1:5 w l, 'ml' u 9:10
    ! gnuplot> set logscale x
    ! gnuplot> plot 'ml' u 1:6, 'ml' u 1:7 w l, 'ml' u 1:8 w l

  end subroutine graph_moodline


  subroutine graph_likemood(filename,mu0,dmu,rho)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: rho,mu0,dmu
!    real(REAL64), dimension(2) :: psi
!    real(REAL64), dimension(2,2) :: dpsi
!    real(REAL64), dimension(:), allocatable :: mju
    real(REAL64) :: mu,phi,logf,s,lambda,nu
    integer :: i,j

    s = 1
    lambda = 0
    nu = 0
!    allocate(mju(size(cts)))
    open(1,file=filename)
    do i = -50,50
       mu = mu0 + 0.1*i*dmu
       do j = 0, 90
          phi = j / RAD
          logf = eloglike(mu,cts,pht,tau,rho,phi,s,lambda,nu)
          !   psi = psisum(mju,cts,pht,tau,rho,fi,g,sx,sy)
          !   dpsi = dpsisum(mju,cts,pht,tau,rho,fi,g,sx,sy)
          write(1,*) mu,RAD*phi,logf
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_likemood

  subroutine graph_xymood(filename,xy,iter,sigu,sigv)

    character(len=*), intent(in) :: filename, xy
    integer, intent(in) :: iter
    real(REAL64), intent(in) :: sigu,sigv
    real(REAL64) :: su,sv,s
    integer :: i,j
    character(len=len(filename)+10) :: xname, yname

    if( xy == 'x' ) then
       write(xname,'(a,a,i0)') trim(filename),'_x_',iter
       open(1,file=xname)
       do i = 50, 150
          sv = 0.01*i*sigv
          do j = 50, 150
             s = 0.01*j
             write(1,*) sv,s,xeloglikes(mju,pht,tau,sv,s,lambda,nu)
          end do
          write(1,*)
       end do
       close(1)
    end if

    if( xy == 'y' ) then
       write(yname,'(a,a,i0)') trim(filename),'_y_',iter
       open(1,file=yname)
       do i = 50, 150
          su = 0.01*i*sigu
          do j = 50, 150
             s = 0.01*j
             write(1,*) su,s,yeloglikes(mju,cts,tau,su,s,lambda,nu)
          end do
          write(1,*)
       end do
       close(1)
    end if

  end subroutine graph_xymood


  real(REAL64) function xeloglike(mu,y,x,tau,rho)

    real(REAL64), dimension(:), intent(in) :: y,x
    real(REAL64), intent(in) :: mu,tau,rho

    real(REAL64) :: s,z,rx,ry,r,r2,sig,dr,dz
    integer :: n,i

    n = size(y)
    sig = sqrt(2*mu + (rho/tau)**2)
    sig = sqrt(mu*tau*(1 + tau) + rho**2)

    s = 0
    do concurrent (i = 1:n)
       rx = x(i) - mu
       ry = y(i) - mu*tau
       r2 = rx**2 + (ry/tau)**2
       if( r2 > 0 ) then
          r = sqrt(r2)
       else
          r = 0
       end if
       r2 = rx**2*tau**2 + ry**2

       !       s = s + itukey(ry/sigy) + itukey(rx/sigx) + log(sigy) + log(sigx)
       dr = -2*(ry/tau + rx)
       dz = 2*((r/sig)*(dr/sig) - (r/sig**2)**2)
       z = (r / sig)**2
       s = s + z + log(z)
!       s = s + (1 - 1/z)*dz
       !s = s + r2/sig**2
    end do

    xeloglike = -s / (2*n) - log(sqrt(2*PI))
!    eloglike = -s / (2*n) - log(sig)

  end function xeloglike


!!$  subroutine graph_likescale(filename,mu0,dmu)
!!$
!!$    character(len=*), intent(in) :: filename
!!$    real(REAL64), intent(in) :: mu0,dmu
!!$    real(REAL64), dimension(2) :: psi
!!$    real(REAL64), dimension(2,2) :: dpsi
!!$    real(REAL64), dimension(:), allocatable :: mju
!!$    real(REAL64) :: f,fi,g
!!$    integer :: i,j,l
!!$
!!$    allocate(mju(size(cts)))
!!$    g = 1
!!$    open(1,file=filename)
!!$    do i = -50,50
!!$       mju = mu0 + 0.1*i*dmu
!!$       do j = 0, 90
!!$          fi = j / RAD
!!$          f = exp(eloglike(mju,cts,pht,tau,rho,fi,g,sx,sy))
!!$          psi = psisum(mju,cts,pht,tau,rho,fi,g,sx,sy)
!!$          dpsi = dpsisum(mju,cts,pht,tau,rho,fi,g,sx,sy)
!!$          write(1,*) mju(1),RAD*fi,f,psi,(dpsi(l,l),l=1,2)
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graph_likescale

  subroutine moodinit(iset,nset,t,rho,y,x,robust,kmean,phi,singular)

    use photonrate_phiopt
    use medians
    use rfun

!    real(REAL64), dimension(:), intent(in) :: cts, pht
    real(REAL64), dimension(:), target, intent(in) :: x,y
    integer, dimension(:), intent(in) :: iset
    integer, intent(in) :: nset
    real(REAL64), intent(in) :: t, rho
    logical, intent(in) :: robust
    real(REAL64), dimension(:), intent(out) :: kmean
    real(REAL64), intent(out) :: phi
    logical, intent(out) :: singular

    integer, dimension(:), allocatable :: idx
    real(REAL64), dimension(:), allocatable :: omju
    real(REAL64) :: dphi, c, stol, s
    integer :: k, m, iter
    logical :: reli, converge

    pht => x
    cts => y
    tau = t
!    rho = r
    sx = 1
    sy = 1
    robustness = robust
    if( robustness ) then
       if( descent ) then
          scale = stukey
       else
          scale = shuber
       end if
    else
       scale = 1
    end if

    allocate(omju(size(cts)))

    do k = 1, nset
       idx = picker(iset,k)
       m = size(idx)
       if( m > 0 ) then
!          kmean(k) = median((cts(idx)*tau + pht(idx)) / (1 + tau**2))
          kmean(k) = median((cts(idx)/tau + pht(idx))) / 2
!          kmean(k) = sum(pht(idx))/m
          omju(idx) =  kmean(k)
          if( debug ) write(error_unit,*) 'Init mu:',k,real(kmean(k))
       else
          stop 'An empty subset in mood.'
       end if
    end do

    s = madscale(omju,cts,pht,tau,rho,phi,sx,sy)
    call phiopt(sx,sy,tau,rho,omju,cts,pht,s,phi,dphi,singular)


!!$    stol = 0.01
!!$    s = scale
!!$    phi = PI / 4
!!$    c = s - scale
!!$    nu = 1
!!$    lambda = nu*c
!!$    converge = .false.

!    lambda = 0
!    nu = 0
!    s = 1
!    do m = 0, 90, 5
!       phi = m / RAD
!       write(*,*) m,eloglikes(omju,y,x,tau,rho,phi,s,lambda,nu)
!    end do


!    do iter = 1, precision(phi)

!       call rortho(omju,cts,pht,phi,dphi,s,reli)

!       if( debug_iter ) then
!102       format(a,i0,(2x,f4.1,a,f0.1),2(2x,a,f0.4),3x,2(2x,a,f0.3),1x,l1)
!!$          write(error_unit,102) 'moodinit #',iter, &
!!$               RAD*phi,'+-',RAD*dphi,'s=',s,' ds=',s-scale, &
!!$               'lambda=',lambda,'nu=',nu,reli
!!$       end if
!!$
!!$       c = s - scale
!!$       converge = abs(c) < stol
!!$       if( converge ) exit
!!$
!!$       lambda = lambda - nu*c
!!$       nu = 2*nu
!!$
!!$    end do

    deallocate(omju)

  end subroutine moodinit


end module photonrate_mood
