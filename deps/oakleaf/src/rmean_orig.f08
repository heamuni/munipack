!
!  Robust mean
!
!  Copyright © 2001 - 2019 F.Hroch (hroch@physics.muni.cz)
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

module robustmean

  ! This module provides rmean subroutine for estimation of robust mean.
  !
  ! rmean should be called as:
  !
  !  call rmean(data,mean,stderr,stdsig,scale,reliable,flag,verbose)
  !
  ! on input:
  !   data - array of data values to be estimated
  !   verbose (optional) - print additional informations
  !
  ! on output are estimated:
  !   mean - robust mean
  !   stderr (optional) - standard error
  !   stdsig (optional) - standard deviation
  !   scale (optional) - scale of function normalisation
  !   reliable (optional) - realiable results?
  !   flag (optional) - indicates reliability of results:
  !          0 == success, results are both reliable and accurate
  !          1 == consider a rough estimate due convergence fail, try verbose
  !          2 == basic estimates, few datapoints available
  !          3 == data looks nearly identical, MAD is zero
  !          5 == failed to allocate memory, initial estimates are returned
  !
  ! The results gives the robust estimate of the true value X
  ! of the sample x, with 68% probability, into the interval
  !
  !         mean - stderr  <  X  <  mean + stderr
  !
  ! Data distribution is approximated as Normal distribution
  !
  !          N(mean,stdsig).
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

  ! 50% quantile of N(0,1)
  real(kind), parameter, private :: Q50 = 0.6745

  ! private data
  real(kind), dimension(:), pointer, private :: x
  real(kind), dimension(:), allocatable, private :: psi
  real(kind), private :: zmean_scale

  private :: rinit, zero_graph, zmean, zfun, rnewton
  public :: rmean_estimator

contains

  subroutine rmean_estimator(data,mean,stderr,stdsig,scale, &
       reliable,flag,verbose)

    ! This routine finds robust mean by looking for the maximum of likelihood.
    ! The scale is estimated from maximum of information.

    use likescale

    real(kind), dimension(:), target, intent(in) :: data
    real(kind), intent(out) :: mean
    real(kind), intent(out), optional :: stderr, stdsig, scale
    logical, intent(out), optional :: reliable
    integer, intent(out), optional :: flag
    logical, intent(in), optional :: verbose

    real(kind) :: t,dt,s,s0,t0,t1,dt1,tmin,tmax,t1min,t1max
    integer :: n, stat,iflag
    logical :: treli, sreli, converge
    character(len=80) :: errmsg

    debug = .false.
    if( present(reliable) ) reliable = .false.
    if( present(verbose) ) debug = verbose
    x => data
    n = size(x)

    if( n == 0 ) then

       ! results are undefined
       mean = 0
       if( present(reliable) ) reliable = .false.
       if( present(stderr) )   stderr = 0
       if( present(stdsig) )  stdsig = 0
       if( present(scale) )  scale = 0
       if( present(flag) ) flag = 2
       return

    else if( n == 1 ) then

       ! dispersion and std.err. are undefined
       mean = x(1)
       if( present(reliable) ) reliable = .false.
       if( present(stderr) )   stderr = 0
       if( present(stdsig) )   stdsig = 0
       if( present(scale) )    scale = 1
       if( present(flag) )     flag = 2
       return

    end if

    !
    ! ... n >= 2, since this point ...
    !

    ! values are initialised by median, MAD and quantiles
    call rinit(x,t0,s0,tmin,tmax,iflag)
    if( iflag /= 0 ) then
       if( present(flag) ) flag = 5
       return
    end if

    ! Initial estimation
    t = t0
    s = s0
    dt = s0 / sqrt(real(n))

    if( debug ) then
       write(*,*) "rinit: t, s= ",t0,s0
       write(*,'(a,3g15.5)') 'tmin, tmax:',tmin,tmax
    end if

    if( abs(s0) < 10*epsilon(s0) ) then
       ! If half of data, or more, has identical values,
       ! their MAD is zero. Any application of robust methods is
       ! uselles (and problematic while s == 0).
       !
       ! If all the data are supposed to be identical (but a few points
       ! can deviate), result can be considered as reliable.

       mean = t0
       if( present(stderr) ) stderr = 0
       if( present(stdsig) ) stdsig = 0
       if( present(scale) )  scale = s0
       if( present(flag) )   flag = 3
       if( present(reliable) ) reliable = .true.
       return

    end if

    ! allocate working arrays
    allocate(psi(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'rmean(): ',trim(errmsg)
       if( present(flag) ) flag = 5
       return
    end if

    if( debug_graph ) call zero_graph(t0,s0)

    ! The first step of estimate results by robust functions.
    ! Newton's method (with derivatives) is used to estimate of t.
    ! If it does not converge, Brent's method (without derivatives)
    ! is its failback alternative. The approach is generally much faster
    ! than preferential application of Brent's: non-convergence interval
    ! of Tukey's function domain is short.
    call rnewton(x,t,dt,s0,treli)
    if( .not. treli ) then
       call zmean(x,s0,tmin,tmax,t,treli)
       if( .not. treli ) t = t0
       if( debug ) write(*,*) 'zmean:',t,treli
    end if

    ! Scale estimation, if above determination of location is successfull
    if( treli ) then
       call iscale(x-t,s,sreli,verbose=debug)
       if( .not. sreli ) s = s0
    else
       sreli = .false.
    end if

    if( debug ) write(*,'(a,1pg20.10,1pg12.3,1x,2l1)') &
         'mean, scale:',t,s,treli,sreli

    if( treli .and. sreli ) then

       ! The second step.
       ! Newton's iterations are only used to obtain even precise result
       ! and to estimate of std.err.
       t1 = t
       t1min = t - 3*s
       t1max = t + 3*s
       call rnewton(x,t1,dt1,s,converge)
       if( debug ) write(*,'(a,1pg20.10,1p2g17.7,2l1)') &
            'Newton:',t1,dt1,converge,t1min < t1 .and. t1 < t1max

       converge = converge .and. (t1min < t1 .and. t1 < t1max)
       if( converge ) then
          t = t1
          dt = dt1
       end if
    else
       converge = .false.
       dt = s / sqrt(real(n))
    end if

    if( debug ) write(*,'(a,1g20.10,2g12.3,3l1)') &
         'rmean:',t,dt,s,treli,sreli,converge

    mean = t
    if( present(stderr) ) stderr = dt
    if( present(stdsig) ) stdsig = dt * sqrt(real(n))
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

  end subroutine rmean_estimator


  !-------------------------------------------------------------------
  subroutine rinit(x,t,s,tmin,tmax,flag)

    ! Initial estimate of parameters,
    ! it is extended version of qmean().

    use quantilmean

    real(kind), dimension(:), intent(in) :: x
    real(kind), intent(out) :: t,s,tmin,tmax
    integer, intent(out) :: flag
    real(kind), dimension(5), parameter :: q = [0.2, 0.25, 0.5, 0.75, 0.8]
    real(kind), dimension(size(q)) :: f

    call quantilefun(x,q,f,flag)
    if( flag /= 0 ) return

    t = (f(3) + (f(4) + f(2))/2)/2
    s = (f(4) - f(2)) / (2*Q50)
    tmin = f(1)
    tmax = f(5)

  end subroutine rinit


  !-----------------------------------------------------------------------

  subroutine rnewton(x,t,dt,s,converge)

    ! Newton's method has been the default robust estimator for a long time.
    ! Only the location parameter (t) is estimated (by Newton's iterations).
    ! The scale has its unchanged estimated value during all computations.
    !
    ! The reliability of result randomly depends on an initial estimate of t.
    ! There is no reliability if the estimate is out of convergence region
    ! of Newtons method or the scale s is far from its right value. The region
    ! of convergence depends on robust function in use. Two most frequent
    ! cases are Huber's function which is monotone and converges in full range
    ! and Tukey's (or Hampel's) function which converges only in domain
    ! determined by extremes of the functions.

    use rfun

    real(kind), dimension(:), intent(in) :: x
    real(kind), intent(in out) :: t
    real(kind), intent(out) :: dt
    real(kind), intent(in) :: s
    logical, intent(out) :: converge

    real(kind), dimension(:), allocatable :: dpsi
    character(len=81) :: errmsg
    integer :: n,it,stat
    real(kind) :: d,f,f2,df,tol

    converge = .false.
    if( .not. (s > epsilon(s)) ) return

    n = size(x)

    allocate(dpsi(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'Newton rmean: ',trim(errmsg)
       return
    end if

    ! set-up required precision
    tol = 10*n*epsilon(t)*(abs(t) + 1)

    ! Newton's iterations
    do it = 1, precision(t)
       ! Number of iterations is given by numerical precision.
       ! We belives, two orders are reached by each iteration, at least.

       call tukeys((x-t)/s,psi,dpsi)
       f = sum(psi)
       df = sum(dpsi)

       ! |df| will be vanishing when the important part of data has large
       ! deviations (above 3.4 for Hampel's piecewise or 1.6 for Tukey's
       ! bi-weight), this si "re-descending M-estimate" function problem.
       ! In this case, the computation is immediatelly interrupted and
       ! result should by considered as unreliable.
       ! The test is conservative, tiny() should be fine as well as.
       if( .not. (df > epsilon(df)) ) exit

       ! corrector
       d = s * f / df

       ! update location
       t = t + d

       if( debug ) write(*,'(a,i2,1pg20.10,1p3g10.2)') &
            "mean, inc., f, f': ",it,t,d,-f/s,df/s**2

       ! The exit of iterations is immediately done when updating
       ! takes less than 1/1000 of standard error indicating
       ! no further progress in iterations will give even precise result.
       ! As the fallback, the numerical tolerance is used.
       ! No convergence is indicated otherwise.

       f2 = sum(psi**2)
       if( f2 > 0 ) then  ! condition df > 0 is valid already
          ! estimation of standard error
          dt = s*sqrt(f2/df**2*n/(n-1)) ! the order does matter
          converge = abs(d) < 0.001*dt
       else ! f2 == 0, t is out of interval -c .. c of Tukey's function
          dt = 0
          if( abs(d) < tol ) exit
       end if

       if( converge ) exit
    enddo

  end subroutine rnewton

  subroutine zmean(x,s,tmin,tmax,t,reli)

    ! Zero mean estimates mean of robust function by Brent's method
    ! which is combination of bisection and inverse quadratic interpolation.
    ! It is important that it works without derivatives so there is
    ! no potential risk of divergence due to "re-descending" property
    ! of Tukey's function as in rnewton() with an improper initial.

    use fmm

    real(kind), dimension(:),intent(in) :: x
    real(kind), intent(in) :: s, tmin, tmax
    real(kind), intent(out) :: t
    logical, intent(out) :: reli

    real(kind) :: tol, d
    integer :: n

    n = size(x)
    tol = 0.01 * s / sqrt(real(n))
    zmean_scale = s
    t = zeroin(tmin,tmax,zfun,tol)

    ! tol gives precision of localisation by expected std.err.
    ! Number of calls of zfun will approximately log((tmax-tmin)/d).

    ! result too close to endpoints shows convergence fail
    d = 2*tol
    reli = abs(t - tmin) > d .and. abs(tmax - t) > d

  end subroutine zmean

  function zfun(t)

    use rfun

    real(kind) :: zfun
    real(kind), intent(in) :: t
    real(kind) :: s

    s = zmean_scale
    psi = tukey((x-t)/s)
    zfun = sum(psi) / size(psi)

  end function zfun



  !-------------------------------------
  !
  ! diagnostics routines
  !

  subroutine zero_graph(t0,s)

    real(kind), intent(in) :: t0,s

    integer :: i
    real(kind) :: t

    zmean_scale = s

    open(1,file='/tmp/t')
    do i = -500,500,10
       t = t0 + s * i / 100.0
       write(1,*) t,zfun(t)
    end do
    close(1)

  end subroutine zero_graph


end module robustmean
