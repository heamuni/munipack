!
!  Robust weighted mean
!
!  Copyright © 2016-9 F.Hroch (hroch@physics.muni.cz)
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
!  along with OakLeaf.  If not, see <http://www.gnu.org/licenses/>.
!
!

module weightedmean

  ! This module provides subroutines for estimation of both
  ! robust mean and its deviations.
  !
  ! rwmean should be called as:
  !
  !  call rwmean(data,errors,mean,stderr,stdsig,scale,reliable,flag,verbose)
  !
  ! on input:
  !   data - array of data values to be estimated
  !   errors - array of statistical errors of x, dx > 0
  !   verbose (optional) - print additional informations
  !
  ! on output are estimated by robust way:
  !   mean - robust mean
  !   stderr (optional) - standard error
  !   stdsig (optional) - standard deviation
  !   scale (optional) - scale of function normalisation
  !   reliable (optional) - indicates reliablity of results
  !   flag (optional) - indicates reliability of results:
  !          0 == success, results are both reliable and precise
  !          1 == consider rough estimate due convergence fail, try verbose
  !          2 == basic estimates, few datapoints available
  !          3 == data looks nerly identical, MAD is zero
  !          4 == improper input parameters (sizes of arrays errors, data
  !               are unequal, any element of error <= 0), try verbose
  !          5 == failed to allocate memory
  !
  ! The results gives the robust estimate of the true value X
  ! of the sample x with weights dx, with 68% probability, into the interval
  !
  !         mean - stderr  <  X  <  mean + stderr
  !
  ! Data distribution is approximated as Normal distribution
  !
  !          N(mean,stdsig).
  !
  !
  ! WARNING
  !
  ! This routine gives unexpected results when setup of input errors
  ! is incorrect. It is indicated by significant violation of scale /= 1.
  ! If you encountered it, run plain rmean() instead. If errors holds some
  ! additional property or condition, use it.
  !
  !
  ! Robust estimators has been prepared on base of
  !  * Hogg in Launer, Wilkinson: Robustness in Statistics
  !  * Hubber: Robust Statistics
  !  * my experiences

  use iso_fortran_env

  implicit none
  private

  ! print additional debuging results
  logical, private :: debug = .false.
  logical, private :: debug_graph = .false.

  ! numerical precision of real numbers
  integer, parameter, private :: kind = selected_real_kind(15)
  ! or use:
  ! integer, parameter, private :: kind = REAL64

  ! 50% quantil of N(0,1)
  real(kind), parameter, private :: Q50 = 0.6745

  ! private working variables
  real(kind), dimension(:), allocatable, private :: psi
  real(kind), dimension(:), pointer, private :: x, dx ! pointers on inputs
  real(kind) :: zmean_scale ! estimated scale for zmean

  ! Allocate of psi can be elimined by replacing of vector
  ! summations in rnewton and zfun by direct summations. It does
  ! slow-down computations (something is computed repeatelly and
  ! no masive parallelization is possible), so this way is prefered.

  private :: rnewton, zmean, zfun, zero_graph
  public :: rwmean_estimator

contains

  subroutine rwmean_estimator(data,errors,mean,stderr,stdsig,scale, &
       reliable,flag,verbose)

    use likescale

    real(kind), dimension(:), target, intent(in) :: data, errors
    real(kind), intent(out) :: mean
    real(kind), intent(out), optional :: stderr,stdsig,scale
    logical, intent(out), optional :: reliable
    integer, intent(out), optional :: flag
    logical, intent(in), optional :: verbose

    real(kind) :: t,dt,s,sig,s0,t0,t1,dt1,tmin,tmax,t1min,t1max
    integer :: n, stat, iflag
    logical :: converge, treli, sreli
    character(len=80) :: errmsg

    debug = .false.
    if( present(reliable) ) reliable = .false.
    if( present(verbose) ) debug = .true.

    x => data
    dx => errors
    n = size(x)

    if( size(x) /= size(dx) .or. .not. all(dx > 0) ) then
       if( debug ) then
          if( size(x) /= size(dx) ) &
               write(error_unit,*) 'rwmean: size(x) /= size(dx)'
          if( .not. all(dx > 0) ) write(error_unit,*) 'rwmean: any dx < 0'
       end if
       mean = 0
       if( present(reliable) ) reliable = .false.
       if( present(flag) ) flag = 4
       return
    end if

    if( n == 0 ) then

       mean = 0
       if( present(reliable) ) reliable = .false.
       if( present(stderr) ) stderr = 0
       if( present(stdsig) ) stdsig = 0
       if( present(scale) )  scale = 0
       if( present(flag) ) flag = 2
       return

    else if( n == 1 ) then

       mean = x(1)
       if( present(reliable) ) reliable = .true.
       if( present(stderr) )   stderr = dx(1)
       if( present(stdsig) )   stdsig = 1
       if( present(scale) )    scale = 1
       if( present(flag) )     flag = 2
       return

    end if

    !
    ! ... n >= 2 since this point ...
    !

    ! initial estimate of t,s
    call rwinit(x,dx,t0,s0,sig,tmin,tmax,iflag)
    if( iflag /= 0 ) then
       if( present(flag) ) flag = 5
       return
    end if

    ! Init
    t = t0
    s = s0
    dt = (s0*sig) / sqrt(real(n))

    if( debug ) then
       write(*,'(a,1p3g17.7)') 'rwinit:',t0,s0,sig
       write(*,'(a,3g15.5)') 'tmin, tmax:',tmin,tmax
    end if

    if( abs(s0) < 10*epsilon(s0)  )then

       ! If important part of the dataset are identical,
       ! results can be considered as reliable.

       mean = t0
       if( present(stderr) )  stderr = (s0*sig) / sqrt(real(n))
       if( present(stdsig) )  stdsig = s0
       if( present(scale) )   scale = s0
       if( present(flag) )     flag = 3
       if( present(reliable) ) reliable = .true.
       return

    endif

    ! Working arrays allocation: A memory can be allocated on the fly,
    ! but it may slow-down the computations.
    allocate(psi(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'rwmean(): ',trim(errmsg)
       if( present(flag) )     flag = 5
       if( present(reliable) ) reliable = .false.
       return
    end if

    if( debug_graph ) call zero_graph(x,dx,t0,s0)

    ! Newton's method (with derivatives) is used to estimate of t.
    ! If it does not converge, Brent's method (without derivatives)
    ! is its failback alternative. The approach is generally much faster
    ! than direct application of Brent's because an interval of
    ! non-convergence of Tukey's function is short.
    call rnewton(x,dx,s0,t,dt,treli)
    if( .not. treli ) then
       call zmean(s0,sig,tmin,tmax,t,treli)
       if( .not. treli ) t = t0
       if( debug ) write(*,*) 'zmean:',t,treli
    end if

    ! The scale parameter estimates average dispersion of residuals
    ! scaled by theirs, a priory known, statistical errors (dx).
    ! Normally, the value is expected near number one: scale ~= 1.
    if( treli ) then
       call iscale((x-t)/dx,s,sreli,verbose=debug)
       if( .not. sreli ) s = s0
    else
       sreli = .false.
    end if

    if( debug ) &
         write(*,'(a,1p2g17.7,1x,2l1)') 'mean, scale:',t,s,treli,sreli

    if( treli .and. sreli ) then

       ! Newton's iterations are used mainly to get more precise result
       ! and to estimate of std.err. and sig
       t1 = t
       t1min = t - 3*s*sig
       t1max = t + 3*s*sig
       call rnewton(x,dx,s,t1,dt1,converge)
       if( debug ) write(*,'(a,1pg20.10,1p2g17.7,2l1)') &
            'Newton:',t1,dt1,converge,t1min < t1 .and. t1 < t1max

       converge = converge .and. (t1min < t1 .and. t1 < t1max)
       if( converge ) then
          t = t1
          dt = dt1
       end if
    else
       converge = .false.
       dt = (s*sig) / sqrt(real(n))
    end if

    if( debug ) write(*,'(a,1pg20.15,3g10.3,1x,3l1)') &
         'rwmean: ',t,dt,s,sig,treli,sreli,converge

    mean = t
    if( present(stderr) ) stderr = dt
    if( present(stdsig) ) stdsig = dt*sqrt(real(n))
    if( present(scale) ) scale = s
    if( present(reliable) ) reliable = .true.
    if( present(flag) ) then
       if( treli .and. sreli .and. converge ) then
          flag = 0
       else !if( treli .and. (.not. sreli .or. .not. converge) ) then
          flag = 1
       end if
    end if

    deallocate(psi)

  end subroutine rwmean_estimator

  subroutine rnewton(x,dx,s,t,dt,converge)

    ! This routine implements Newton's method for solving of equation
    !
    !  f(t) = 0,  where f(t) = sum( psi((x-t)/(s*dx))/dx )
    !
    ! t   is estimation of mean,
    ! dt  is standard error of t
    ! s   is scale

    use rfun

    real(kind), dimension(:), intent(in) :: x, dx
    real(kind), intent(in out) :: t
    real(kind), intent(in) :: s
    real(kind), intent(out) :: dt
    logical, intent(out) :: converge

    real(kind), dimension(:), allocatable :: dpsi
    character(len=81) :: errmsg
    integer :: n,iter,stat
    real(kind) :: d,f,df,f2,dp,tol

    converge = .false.
    if( .not. (s > 0.0_kind) ) return

    n = size(x)
    allocate(dpsi(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'Newton rwmean: ',trim(errmsg)
       return
    end if

    ! Adjusting of tolerance
    tol = 10*n*epsilon(t)*(abs(t) + 1)

    ! Newton's iterations
    do iter = 1, precision(t)
       ! Number of iterations is limited by numerical precision.
       ! We belives in quadratic convergence of Newton's method;
       ! two orders are reached by each iteration, at least.

       ! derivation of derivatives
       call tukeys((x-t)/(s*dx),psi,dpsi)
       f = sum(psi/dx)
       df = sum(dpsi/dx**2)

       ! derivative (denominator) checkpoint:
       ! * If all values are out of range: |df| < epsilon
       ! * |df| < 0 indicates "re-descending M-estimate" function problem ...
       !   the problem usually appears when distribution strongly deviates
       !   from Normal distribution
       if( .not. (df > epsilon(df)) ) exit

       ! corrector
       d = s * f / df

       ! update location
       t = t + d

       if( debug ) write(*,'(a,i2,1pg20.10,1p3g10.2)') &
            "mean, incr., f, f': ",iter,t,d,-f/s,df/s**2

       ! exit of iterations, see rmean() for explanation
       f2 = sum(psi**2*dx**2)
       dp = sum(dpsi)
       if( f2 > 0 .and. dp > 0 ) then
          dt = s*sqrt(f2/dp**2*n/(n-1))
          converge = abs(d) < 0.001*dt
       else ! f2 == 0
          dt = 0
          if( abs(d) < tol ) exit
       end if

       if( converge ) exit
    enddo

  end subroutine rnewton

  subroutine zmean(s,sig,tmin,tmax,t,reli)

    use fmm

    real(kind), intent(in) :: s,sig,tmin,tmax
    real(kind), intent(out) :: t
    logical, intent(out) :: reli

    real(kind) :: d, tol
    integer :: n

    n = size(x)
    tol = 0.01 * s*sig / sqrt(real(n))
    zmean_scale = s
    t = zeroin(tmin,tmax,zfun,tol)

    ! tol gives precision of localisation by expected std.err.
    ! Number of calls of 0-fun will approximately ln((tmax-tmin)/tol).

    ! if result is too close to interval endpoints, it raises suspicion ...
    d = 2*tol
    reli = abs(t - tmin) > d .and. abs(tmax - t) > d

  end subroutine zmean

  function zfun(t)

    use rfun

    real(kind) :: zfun
    real(kind), intent(in) :: t
    real(kind) :: s

    s = zmean_scale
    psi = tukey((x-t)/(s*dx))
    zfun = sum(psi/dx)

  end function zfun


  subroutine rwinit(x,dx,t,s,sig,tmin,tmax,flag)

    ! This routine performs an initial estimate of both mean and scale
    ! using both median and absolute value.

    use medians
    use quantilmean

    real(kind), dimension(:),intent(in) :: x, dx
    real(kind), intent(out) :: t,s,sig,tmin,tmax
    integer, intent(out) :: flag
    real(kind) :: xmin,xmax

    if( size(x) > 50 ) then

       ! This choice of the threshold used for median computation method
       ! assumes that odd and even elements in sequence are the same (within
       ! requested precision). Note that qmed is the fastest known algorithm
       ! (~2*n) while median is simpler and slower (~n*log(n)).

       t = qmedian(x)
       sig = qmedian(dx)

    else

       t = median(x)
       sig = median(dx)

    end if

    ! determine parameters by empirical CDF quantilies
    call qbracket((x-t)/dx,xmin,xmax,flag)
    if( flag /= 0 ) return

    ! Both mean of x and mean of (x-t0)/dx will different in general (ignored).
    tmin = t + sig*xmin
    tmax = t + sig*xmax
    s = (xmax - xmin) / (2*Q50)

    if( debug_graph ) call graph_ecdf(t,sig)

  end subroutine rwinit



 ! diagnostics routines


  subroutine zero_graph(x,dx,t0,s)

    use rfun
    use medians

    real(kind), dimension(:),intent(in) :: x, dx
    real(kind), intent(in) :: t0,s

    integer :: i,n
    real(kind) :: t,d

    d = median(s*dx) / 0.6745

    n = size(x)
    zmean_scale = s

    open(1,file='/tmp/tw')
    do i = -500,500,10
       t = t0 + d * i / 100.0
       write(1,*) t,zfun(t)
    end do
    close(1)

  end subroutine zero_graph

  subroutine graph_ecdf(t,s)

    use quantilmean

    real(kind), intent(in) :: t,s
    real(kind), dimension(:), allocatable :: xcdf,ycdf
    integer :: n,l,u

    n = size(x)
    allocate(xcdf(n),ycdf(n))
    call ecdf((x-t)/s,xcdf,ycdf)
!    call ecdf(x-t,xcdf,ycdf)

    open(newunit=u,file='/tmp/cdf1')
    do l = 1, n
       write(u,*) xcdf(l),ycdf(l)
    end do
    close(u)

  end subroutine graph_ecdf


end module weightedmean
