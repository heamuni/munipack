!
!  Robust estimation of factor in lambda of Poisson distribution
!     which is average number of events per time period.
!
!  Copyright © 2012 - 2019 F.Hroch (hroch@physics.muni.cz)
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
!  along with Resthmet.  If not, see <http://www.gnu.org/licenses/>.

module fluxmean

  ! This module implements a robust estimator of a ratio of events
  !
  !     t =  photon / counts
  !
  ! (including estimation of its scatter) appearing as the factor
  ! lambda parameter of Poisson distribution
  !
  !     Po(lambda) = lambda**k/k! * exp(-lambda),   lambda = t * C
  !
  !  which provides zero mean for difference of N1 - N2 events
  !  in Skellam's distribution
  !     https://en.wikipedia.org/wiki/Skellam_distribution
  !
  !  The estimation is valid only for large values of lambda >> 1
  !  in Gaussian regime when transformation
  !
  !      (photon - t*counts) / sqrt(photon_err**2 + t**2*count_err**2)
  !
  !  is valid. The approach usually requires proper estimates errors.
  !
  !  Both photons, counts should be non-integers > 1 as a product of
  !  previous computations.
  !
  !
  ! fmean should be called as:
  !
  !  call fmean(photons,photons_err,counts,counts_err,mean,stderr,stdsig,
  !             scale,reliable,poisson,verbose)
  !
  ! On input:
  !   photons - array of reference photon rates
  !   photons_err - array of statistical errors of photons
  !   counts - array of measured rates
  !   counts_err - array of statistical errors of counts
  !   verbose - print additional info
  !
  ! On output are estimated:
  !   mean - mean
  !   stderr (optional) - its standard error
  !   stdsig (optional) - estimate of standard deviation
  !   scale (optional) - estimate of scale
  !   reliable - indicates reliability of results (optional)
  !   poisson - indicates Poisson (.true.) or Normal (.false.)
  !             distribution used for determine of results (optional)
  !
  ! The given results  means that a true value T of the sample
  ! photons / counts can be, with 70% probability, found in interval
  !
  !         t - dt  <  T  <  t + dt
  !
  ! The estimate is performed as maximum of entropy of Poisson density
  ! for values of photons smaller than 'plimit' (666) or by robust estimate
  ! of parameters of Normal distribution. Both the estimates are robust.
  ! Poisson mean is estimated as the minimum of free energy:
  ! https://en.wikipedia.org/wiki/Helmholtz_free_energy

  use iso_fortran_env

  implicit none
  private

  integer, parameter, private :: kind = REAL64

  ! 50% quantil of N(0,1)
  real(kind), parameter, private :: Q50 = 0.6745

  ! debuging
  logical, private :: debug = .false.
  logical, private :: debug_graph = .false.

  ! The lower limit for Poisson statistics use (Normal otherwise)
  integer, parameter, private :: plimit = 666

  ! working arrays
  real(kind), dimension(:), pointer, private :: cts,dcts,pht,dpht
  real(kind), dimension(:), allocatable, private :: res,dcts2,dpht2,psi,sig
  real(kind), private :: zmean_sigma, dt_err

  public :: fmean

  interface fmean
     module procedure fmean64, fmean32
  end interface fmean


contains

  subroutine fmean64(photons,photons_err,counts,counts_err, &
       mean,stderr,stdsig,scale,reliable,poisson,flag,verbose)

    use rfun

    real(kind), dimension(:), target, intent(in) :: &
         photons,photons_err,counts,counts_err
    real(kind), intent(out) :: mean
    real(kind), intent(out), optional :: stderr,stdsig,scale
    logical, intent(out), optional :: reliable, poisson
    integer, intent(out), optional :: flag
    logical, intent(in), optional :: verbose

    integer :: n,i,nplus,nsig,nstep,stat,nflag,iflag
    real(kind) :: t,dt,d0,c,d,r,w,s,t0,s0,tmax,tmin,h11
    logical :: reli, poi
    character(len=80) :: errmsg

    debug = .false.
    if( present(verbose) ) debug = verbose

    ! data check
    if( .not. (all(photons > 0) .and. all(counts > 0) .and. &
         all(photons_err > 0) .and. all(counts > 0) ) ) then
       write(error_unit,*) 'fmean: some data has invalid values (needs all>0).'
       if( present(flag) ) flag = 4
       if( present(reliable) ) reliable = .false.
       return
    end if

    if(  size(counts) /= size(photons) .or. &
         size(counts) /= size(counts_err) .or. &
         size(counts_err) /= size(photons_err) ) then
       write(error_unit,*) 'fmean: input data of inconsistent length.'
       if( present(flag) ) flag = 4
       if( present(reliable) ) reliable = .false.
       return
    end if

    n = size(photons)
    pht => photons
    cts => counts
    dpht => photons_err
    dcts => counts_err

    if( n == 0 ) then ! no data
       mean = 1
       if( present(flag) ) flag = 2
       if( present(reliable) ) reliable = .false.
       if( present(stderr) ) stderr = 0
       if( present(stdsig) ) stdsig = 0
       if( present(scale) )  scale = 0
       return
    endif

    if( n == 1 ) then               ! a single point
       mean = photons(1) / counts(1)
       s = sqrt(photons_err(1)**2 + mean**2*counts_err(1)**2)
       if( present(stderr) ) stderr = s
       if( present(stdsig) ) stdsig = s
       if( present(scale) )  scale = 1
       if( present(poisson) ) poisson = .false.
       if( present(flag) ) flag = 2
       if( present(reliable) ) reliable = .true.
       return
    endif

    allocate(res(n),sig(n),psi(n),dpht2(n),dcts2(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'fmean: ',trim(errmsg)
       if( present(reliable) ) reliable = .false.
       if( present(flag) ) flag = 5
       return
    end if

    dpht2 = dpht**2
    dcts2 = dcts**2

    ! initial estimate of mean,scale.. by quantiles
    call rinit(t,s,d0,h11,tmin,tmax,iflag)
    dt_err = s*d0/sqrt(h11*n)

    ! listing
    if( debug ) then
       i = max(int(log10(real(n))) - 1,0)
       nstep = 10**i
       if( nstep > 1 ) &
         write(*,'(a,i0,a)') 'Info: printing only every ',nstep,'th star.'
       call residuals(t,res,sig)
       write(*,*) 'Data table: pht(i),dpht(i),cts(i),dcts(i), ratio, residuals'
       do i = 1,n,nstep
          write(*,'(i3,1p,2(g15.5,g10.3),0pf11.3,f10.3)') &
               i,pht(i),dpht(i),cts(i),dcts(i),pht(i)/cts(i),res(i)/sig(i)
       end do
       write(*,'(a50,g0.4,1x,es7.1,1x,f0.3)') &
            "Initial mean, stderr, scale by quantiles: ",t,s*d0/sqrt(h11*n),s
       write(*,'(a,3g15.5)') 'tmin, tmax:',tmin,tmax
       if( debug_graph ) call graph_ecdf(t,s)
    end if

    if( abs(s) < 10*epsilon(s) )then

       ! Over 75% of data contains identical values,
       ! results can be considered as reliable.

       mean = t
       if( present(stderr) )  stderr = d0 / sqrt(real(n - 1))
       if( present(stdsig) )  stdsig = s
       if( present(scale) )   scale = s
       if( present(reliable) ) reliable = .true.
       if( present(poisson) ) poisson = .false.
       if( present(flag) ) flag = 3
       deallocate(res,sig,psi,dpht2,dcts2)
       return

    endif

!    if( debug .and. t < plimit ) &
!         write(*,*) 'Low fluxes determination is imprecise.'
!    if( t < plimit ) then
    if( .false. .and. t < plimit ) then

       t0 = t
       t = pmean(s*d0/sqrt(h11*n),tmin,tmax,reli)
       if( reli .and. t > 0 ) then
          s = sqrt(t)
       else
          t = t0
       end if
       dt = s / sqrt(real(n - 1))


    else

       t0 = t
       s0 = s
       dt_err = s*d0/sqrt(h11*n)
       do i = 1, 2
          call rnorm(t,dt,s,d0,h11,tmin,tmax,reli,nflag)
          if( debug ) write(*,*) 'rnorm: ',t,dt,s,reli,nflag
          if( .not. reli ) then
             t = t0
             dt = s0 / sqrt(real(n))
             s = s0
             exit
          end if
          dt_err = dt
       end do

    end if

    ! residual sum
    if( debug ) then
       call residuals(t,res,sig)
       write(*,*) '# photons, computed counts, rel.err, expected rel.err, res., huber, dhuber:'
       do i = 1,n,nstep
          c = t*cts(i)
          d = (pht(i) - c)/((pht(i) + c)/2)
          w = sqrt((dpht(i)/t)**2 + dcts(i)**2)/((cts(i) + pht(i)/t)/2.0)
          r = res(i) / (sig(i)*s)
          write(*,'(i3,1p2g13.4,0p2f10.5,3f8.3)') &
               i,pht(i),c,d,w,r,huber(r),dhuber(r)
       enddo
       nsig = count(abs(res/sig) < s)
       nplus = count(res > 0 .and. abs(res/sig) < s)
       write(*,'(a,g0.5,a,1pg0.1,a,1pg0.2,5x,a,g0.3)') &
            'mean = ',t,' +- ',dt,'   stdsig = ',dt*sqrt(real(n)),'scale = ',s
       write(*,'(a,i0,a,i0,a,f0.1,a)') &
            'Points within 1-sigma error interval: ', &
            nsig,'/',n,' (',100.*real(nsig)/real(n),'%)'
       write(*,'(a,i0,"/",f0.1,", ",i0,"+-",f0.1)') &
            'Sign test in 1-sigma interval (total/expected, +): ',&
            nsig,nsig/2.0,nplus,sqrt(nplus*0.25 + 1e-3)
       ! test whatever the errors corresponds to values to prove Poisson origin
       ! it is same kind quantity as
       ! https://en.wikipedia.org/wiki/Index_of_dispersion
       write(*,'(a,2(1x,g0.3,1x,a,:","))') 'Mean Poissonity (best=1): ',&
            sqrt(sum(dpht**2/pht)/n),'(reference)', &
            sqrt(sum(dcts**2/cts)/n),'(data)'
       write(*,*) 'Poissonity:',poi,'   Reliable:',reli,'     Flag:',nflag
    end if

    mean = t
    if( present(stderr) ) stderr = dt
    if( present(stdsig) ) stdsig = dt*sqrt(real(n))
    if( present(scale) )  scale = s
    if( present(reliable) ) reliable = reli
    if( present(poisson) ) poisson = t < plimit
    if( present(flag) ) flag = nflag

    deallocate(res,sig,psi,dpht2,dcts2)

  end subroutine fmean64


  subroutine fmean32(photons,photons_err,counts,counts_err, &
       mean,stderr,stdsig,scale,reliable,poisson,flag,verbose)

    real(REAL32), dimension(:), target, intent(in) :: &
         photons,photons_err,counts,counts_err
    real(REAL32), intent(out) :: mean
    real(REAL32), intent(out), optional :: stderr,stdsig,scale
    logical, intent(in), optional :: verbose
    logical, intent(out), optional :: reliable, poisson
    integer, intent(out), optional :: flag

    real(kind) :: t,dt,s,c
    real(kind), dimension(:), allocatable :: x,dx,y,dy
    integer :: stat, iflag
    logical :: reli, poi
    character(len=80) :: errmsg

    allocate(x(size(photons)),dx(size(photons_err)), &
         y(size(counts)),dy(size(counts_err)),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'fmean32: ',trim(errmsg)
       if( present(reliable) ) reliable = .false.
       if( present(flag) ) flag = 5
       return
    end if

    x = photons
    dx = photons_err
    y = counts
    dy = counts_err

    if( present(verbose) ) then
       call fmean64(x,dx,y,dy,t,dt,s,c,reli,poi,iflag,verbose)
    else
       call fmean64(x,dx,y,dy,t,dt,s,c,reli,poi,iflag)
    end if

    mean = real(t)
    if( present(stderr) ) stderr = real(dt)
    if( present(stdsig) ) stdsig = real(s)
    if( present(scale) )   scale = real(c)
    if( present(reliable) ) reliable = reli
    if( present(poisson) ) poisson = poi
    if( present(flag) )       flag = iflag

  end subroutine fmean32


  ! --------------------------------------------------------

  subroutine rnorm(t,dt,s,d0,h11,tmin,tmax,reliable,flag)

    use rfun
    use medians
    use likescale
    use quantilmean

    real(kind), intent(in out) :: t,s
    real(kind), intent(out) :: dt
    real(kind), intent(in) :: d0,h11,tmin,tmax
    logical, intent(out) :: reliable
    integer, intent(out) :: flag

    real(kind) :: s0,t0,t1,dt1,d
    logical :: converge, treli, sreli
    integer :: n

    reliable = .false.
    t0 = t
    s0 = s
    n = size(pht)

    call rnewton(t,dt,s,treli)
    if( .not. treli ) then
       call zmean(s,s*d0/sqrt(h11*n),tmin,tmax,t,treli)
       if( .not. treli ) t = t0
       if( debug ) write(*,*) 'zmean:',t,treli,s*d0/sqrt(h11*n),dt
    end if

    if( treli ) then
       call residuals(t,res,sig)
       call iscale(res/sig,s,sreli,verbose=debug)
       if( .not. sreli ) s = s0
    else
       sreli = .false.
    end if

    if( debug ) write(*,'(a50,1pg15.5,1x,0pg8.2,1x,2l1)') &
         'Robust mean estimation, scale by information:',t,s,treli,sreli

    if( treli .and. sreli ) then

       t1 = t
       call rnewton(t1,dt1,s,converge)
       if( debug ) write(*,'(a20,1pg20.10,1p2g17.7,2l1)') &
            'Newton:',t1,dt1,converge,tmin < t1 .and. t1 < tmax

       converge = converge .and. (tmin < t1 .and. t1 < tmax)
       if( converge ) then
          t = t1
          dt = dt1
       end if
    else
       converge = .false.
       call residuals(t,res,sig)
       d = median(abs(res/(s*sig))) / Q50
       dt = d / sqrt(h11*n)
    end if

    reliable = treli
    if( treli .and. sreli .and. converge ) then
       flag = 0
    else
       flag = 1
    end if

    if( debug_graph ) then
       call rmean_graph(t - 3*dt,t + 3*dt,s)
       call poisson_graph
    end if

  end subroutine rnorm

  subroutine residuals(t,res,sig)

    real(kind), intent(in) :: t
    real(kind), dimension(:), intent(out) :: res, sig
    real(kind) :: t2

    t2 = t**2
    res = (pht - t*cts)
    sig = sqrt(t2*dcts2 + dpht2 + dt_err**2*cts**2)

  end subroutine residuals


  subroutine rnewton(t,dt,s,converge)

    use rfun

    real(kind), intent(in out) :: t
    real(kind), intent(out) :: dt
    real(kind), intent(in) :: s
    logical, intent(out) :: converge

    integer :: n,iter
    real(kind) :: tol,d,f,df,f2,dp,q
    real(kind), dimension(:), allocatable :: r,dpsi,ds,ds2,cs,dd,d1s,d2s

    converge = .false.

    if( .not. (s > 0) ) return

    n = size(pht)
    allocate(r(n),dpsi(n),ds(n),cs(n),ds2(n),dd(n),d1s(n),d2s(n))

    ! tolerance limit will depends on number of elements in sums
    ! because every element can include its rounding error...
    ! ... proportional to absolute value ...
    tol = 10*n*epsilon(t)*(abs(t) + 1)

    ! Number of iterations is limited by numerical precision.
    ! We belives, as optimists, that one order is reached by each iteration.
    ! But, the pesimist, inside us, says that the convergence can be slower.
    do iter = 1, precision(t)

       ! solution of f(t) = 0, where f(t) = sum(psi((n-t*c)/(s*ds))
       call residuals(t,r,ds)
       ds2 = ds**2
       call tukeys(r/(s*ds),psi,dpsi)
       dd = pht - t*cts
       d1s = t*dcts**2/ds
       d2s = dcts**2*(ds - t*d1s)/ds2
       cs = cts*ds + dd*d1s

       f = t*sum((psi*(r/ds)-1.0_kind)*(dcts/ds)**2) + sum(psi*cts/ds)
       df = sum(dpsi*(cts/ds + t*(r/ds)*dcts2/ds2)**2) &
           -sum(psi*(pht - 3*t*cts)*(dcts/ds)**2/ds) &
           +sum(psi*3*(r/ds)*t**2*(dcts/ds)**4) &
           +sum((dcts/ds)**2) - 2*t**2*sum((dcts2/ds2)**2)

       if( .not. (abs(df) > epsilon(df)) ) exit

       ! corrector for mean
       d = s * f / df

       ! update location
       t = t + d

       if( debug ) write(*,'(a,i2,1pg14.5,1p3g12.3)') &
            "mean, incr., f, f': ",iter,t,d,f,df

       f2 = sum(psi**2)  ! mean of psi**2
       dp = sum(dpsi)    ! mean of dpsi
       q = df*dp         ! df is an element of Jacobian
       if( f2 > 0 .and. q > 0 ) then
          dt = s*sqrt(f2/q*n/(n-1))
          converge = abs(d) < 0.001*dt
       else ! f2 == 0 or q == 0, t > |c|
          dt = 0
          if( abs(d) < tol ) exit
       end if

       ! exit of iterations: the absolute errot must be at least |d| < tol
       if( converge ) exit

    enddo

  end subroutine rnewton

  subroutine rinit(t,s,d0,h11,tmin,tmax,flag)

    ! This is the very first estimation of ratio parameter.
    use medians
    use quantilmean

    real(kind), intent(out) :: t,s,d0,h11,tmin,tmax
    integer, intent(out) :: flag

    integer, parameter :: medlim = 50 ! quick median limit
    real(kind) :: d,xmin,xmax

    if( size(pht) > medlim ) then
       t = qmedian(pht/cts)
       d = t*qmedian(sqrt(dcts2/cts**2 + dpht2/pht**2))
       h11 = qmedian(pht)
    else
       t = median(pht/cts)
       d = t*median(sqrt(dcts2/cts**2 + dpht2/pht**2))
       h11 = median(pht)
    end if
    h11 = h11 / (t**2*(1 + t)) ! /s**2
    dt_err = 0
    call residuals(t,res,sig)

    ! determine parameters by empirical CDF quantilies
    call qbracket(res/sig,xmin,xmax,flag)
    if( flag /= 0 ) return

    ! Medians of pht/cts and res slightly differs (ignored).
    tmin = t + d*xmin
    tmax = t + d*xmax
    s = (xmax - xmin) / (2*Q50)

    if( s < 10*epsilon(s) ) then
       ! identical data
       d0 = 0
       flag = 5
       return
    end if

    if( size(pht) > medlim ) then
       d0 = qmedian(abs(res/(s*sig))) / Q50 ! qbracket?
    else
       d0 = median(abs(res/(s*sig))) / Q50 ! qbracket?
    end if

  end subroutine rinit

  subroutine zmean(s,d,tmin,tmax,t,reli)

    use fmm

    real(kind), intent(in) :: s,d,tmin,tmax
    real(kind), intent(out) :: t
    logical, intent(out) :: reli
    real(kind) :: tol,h

    zmean_sigma = s
    tol = 0.001 * d
    t = zeroin(tmin,tmax,zfun,tol)

    h = 2*tol
    reli = tmin + h < t .and. t < tmax - h

    ! a realiable result is inside range tmin .. tmax,
    ! also not too close to the interval endpoints

  end subroutine zmean


  function zfun(t)

    use rfun

    real(kind), intent(in) :: t
    real(kind) :: zfun, sigma

    sigma = zmean_sigma
    call residuals(t,res,sig)
    psi = tukey(res/(sigma*sig))

    zfun = ( t*sum((psi*(res/sig) - 1.0_kind)*(dcts/sig)**2) + &
         sum(psi*cts/sig) ) / sigma

  end function zfun


  function pmean(d,tmin,tmax,reli) result(t)

    use fmm

    real(kind) :: t
    real(kind), intent(in) :: d,tmin, tmax
    logical, intent(out) :: reli
    real(kind) :: h,tol

    tol = 0.001 * d
    t = fmin(tmin,tmax,negfree,tol)

    h = 2*tol
    reli = tmin + h < t .and. t < tmax - h

  end function pmean


  function negfree(t)

    use rfun, A => etukey

    ! negative free energy

    ! the largest argument for which exp() should be computed without
    ! issuing FPE flag
    real(kind), parameter :: emax = log(huge(emax)) / A

    real(kind), intent(in) :: t
    real(kind) :: negfree
    real(kind), dimension(:), allocatable :: f

    allocate(f(size(pht)))

    ! The density for Poisson distribution is computed by numerically reliable
    ! way: https://en.wikipedia.org/wiki/Poisson_distribution#Definition

    f = pht*log(t*cts) - t*cts - log_gamma(pht + 1)
    where( f < emax )
       f = exp(f)
    elsewhere
       f = 0
    end where
    negfree = - sum(f)
    deallocate(f)

  end function negfree


  ! --------------------------------------------------------

  subroutine poisson_graph

    integer, parameter :: nmax = 1000
    integer :: n
    real(kind) :: t,dt,rmin,rmax
    real(kind), dimension(:), allocatable :: f

    allocate(f(size(pht)))

    rmax = maxval(pht/cts)
    rmin = minval(pht/cts)
    dt = (rmax - rmin) / nmax

    open(1,file='/tmp/f')
    do n = 1,nmax,2
       t = rmin + dt*n
       f = pht*log(t*cts) - t*cts - log_gamma(pht+1)
       write(1,*) t,sum(exp(f)),sum(f*exp(f))
    end do
    close(1)

    deallocate(f)

  end subroutine poisson_graph

  subroutine rmean_graph(t1,t2,s)

    use rfun

    real(kind), intent(in) :: t1,t2,s
    real(kind) :: t,dt
    real(kind), dimension(:), allocatable :: r,rho
    integer :: i,n

    if( .not. (s > 0) ) return

    n = size(pht)
    allocate(r(n),rho(n))

    dt = (t2 - t1) / 1000

    open(1,file='/tmp/r')
    do i = 1,1000,10
       t = t1 + dt*i

       call residuals(t,res,sig)
       rho = ihuber(res/(s*sig))

       write(1,*) t,sum(rho)
    end do
    close(1)

    deallocate(r,rho)

  end subroutine rmean_graph



  subroutine graph_ecdf(t,s)

    use quantilmean

    real(kind), intent(in) :: t,s
    real(kind), dimension(:), allocatable :: xcdf,ycdf
    integer :: n,l,u

    if( .not. s > 0 ) return

    n = size(pht)
    allocate(xcdf(n),ycdf(n))
    call residuals(t,res,sig)
    call ecdf(res/(s*sig),xcdf,ycdf)

    open(newunit=u,file='/tmp/cdf1')
    do l = 1, n
       write(u,*) xcdf(l),ycdf(l)
    end do
    close(u)

  end subroutine graph_ecdf

end module fluxmean
