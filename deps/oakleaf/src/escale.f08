!
!  Scale by entropy, or quantiles
!
!  Copyright © 2017-8 F.Hroch (hroch@physics.muni.cz)
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

module entropyscale

  ! This module provides escale subroutine for estimation of scale
  ! by maximalising of relative entropy.
  !
  ! escale should be called as:
  !
  !  call escale(r,s,reliable,flag,verbose)
  !
  ! on input:
  !   r - array of residuals: (data-mean), (data-mean)/errors, ...
  !   verbose (optional) - print additional informations
  !
  ! on output are estimated:
  !   s - scale
  !   reliable (optional) - indicates reliability of result
  !                         if .true., scale has been correctly
  !                         localised. if .false., s is undefined
  !   flag (optional) - indicates reliability of results:
  !          0 == success, results are both reliable and precise
  !          1 == consider rough estimate due convergence fail, try verbose
  !          2 == basic estimates, few datapoints available
  !          3 == data looks nearly identical
  !          5 == failed to allocate memory, initial estimates are returned
  !
  ! Caveats.
  !
  ! The implemented algorithm looks for minimum of negative entropy.
  ! Sometimes, commonly, for strongly non-Normal data, the function
  ! has more minimums. It can be impossible to distinguish between them
  ! (a false and the right minimus). The minimum  is selected with
  ! certain probability of false success.  As a consequence, more
  ! separated groups of results can be given which looks as a bifurcation.
  ! Result can be always verified by a K-S kind test.
  !
  ! Note.
  ! The scale parameter corrects oddly estimated statistical
  ! errors of passed data (dx). For most applications,
  ! the value is expected near number one: escale ~= 1.
  !
  ! https://en.wikipedia.org/wiki/Half-normal_distribution
  !

  use rfun
  use iso_fortran_env

  implicit none
  private

  ! print debug information ?
  logical, private :: debug = .false.
  logical, private :: debug_files = .false.

  ! numerical precision of real numbers
  integer, parameter, private :: kind = selected_real_kind(15)

  ! minimal precision of results
  real(kind), parameter, private :: etol = 0.1

  ! scale entropy parameter
  real(kind), parameter, private :: A = etukey

  ! the largest argument for which exp() can be computed without FPE flag issue
  real(kind), parameter, private :: rhomax = log(huge(rhomax)) / A

  ! data buffer
  real(kind), dimension(:), allocatable, private :: rho

  ! data pointer, just pointer
  real(kind), dimension(:), pointer, private :: res

  private :: negentropy, newton, grid, gmin

  public :: escale, escale_graph

contains

  subroutine escale(r,s,reliable,flag,verbose)

    ! An estimator of scale parameter as maximum of entropy
    ! (implemented as minimum of negative entropy).

    use quantilmean

    real(kind), dimension(:), target, intent(in) :: r
    real(kind), intent(out) :: s
    integer, intent(out), optional :: flag
    logical, intent(out), optional :: reliable
    logical, intent(in), optional :: verbose

    real(kind), dimension(:), allocatable :: xcdf, ycdf
    character(len=80) :: errmsg
    real(kind) :: s0,s1,smin,smax,qmin,qmax,qhalf
    integer :: n, stat
    logical :: found,converge

    debug = .false.
    if( present(verbose) ) debug = verbose

    n = size(r)
    res => r

    if( n == 0 ) then
       s = 0
       if( present(flag) ) flag = 2
       if( present(reliable) ) reliable = .false.
       return
    end if

    if( n == 1 ) then
       s = abs(res(1))
       if( present(flag) ) flag = 2
       if( present(reliable) ) reliable = .true.
       return
    end if

    if( debug ) write(*,'(a,2g15.3)') &
         'Maximum and minimum: ',minval(abs(res)),maxval(abs(res))

    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'escale: ',trim(errmsg)
       if( present(flag) ) flag = 5
       if( present(reliable) ) reliable = .false.
       return
    end if

    ! cummulative distribution function
    call ecdf(abs(res),xcdf,ycdf)

    ! keep eCDF
    block
      integer :: l
      if( debug_files ) then
         open(7,file='/tmp/c')
         do l = 1, n
            write(7,*) xcdf(l),ycdf(l)
         end do
         close(7)
      end if
    end block

    ! the very first, rought, estimate of scale at Q(0.6745), backup estimate
    qhalf = 0.6745
    s0 = quantile(qhalf,xcdf,ycdf)

    ! the intial interval is determined by qmin and qmax quantiles
    qmin = 0.4    ! = 0.5*s
    qmax = 0.87   ! = 1.5*s

    ! quantile for Q(0.87) for 1.5*x, or 1.5*s0, gives upper bound
    if( qmax < ycdf(n) ) then
       smax = quantile(qmax,xcdf,ycdf)
    else ! if( Q(0.87) is out of CDF values ) then
         ! the interval of q is uncovered for a few points
       smax = 1.5 * s0  ! = Q(4/5) * s0
    end if

    smin = max(quantile(qmin,xcdf,ycdf), epsilon(smin))

    deallocate(xcdf,ycdf)

    if( debug ) write(*,'(a,2g15.3,a,g15.3)') &
         'Quantiles bracketing:',smin,smax,'   Qmean:',s0

    ! check: smax > smin, s > 0, identical or zero value data
    if( (smax - smin) < 10*epsilon(s) .or. abs(smax) < 10*epsilon(s) ) then
       s = s0
       deallocate(rho)
       if( present(flag) ) flag = 3
       if( present(reliable) ) reliable = .true.
       return
    end if

    allocate(rho(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'escale: ',trim(errmsg)
       if( present(flag) ) flag = 5
       if( present(reliable) ) reliable = .false.
       return
    end if

    ! Optimum determination
    ! The extreme is searched by an optimisation method without derivatives.
    ! Determination is reliable but slower than follow-up Newton method.
    if( debug_files ) open(7,file='/tmp/s')
    if( n < 17 ) then

       ! If there are a few elements only, we will inspect full
       ! interval given by range of |res|. Any mixed distributions
       ! (Normal with some contamination) or plain random number
       ! generator will cause function with more minimums. We will
       ! select the most probable one, but is it the right one? Maybe?
       ! To distinguish the right minimum, any additional information
       ! should be included.

       call grid(smin,smax,s,found)

    else ! if( n is large ) then

       call gmin(smin,smax,s,found)

    end if
    if( debug_files ) close(7)

    if( found ) then

       ! We are attemping to do 10x more precise estimate of the scale
       ! by fast Newtons's method. The convergence sucess is given
       ! by starting point.
       s1 = s
       call newton(s1,0.1*etol,converge)
       ! Newton method is kind of root-find metod. We are looking for
       ! maximum of a function so we work with derivatives of the function.
       ! By second derivative, we can choice if it is truly maximum
       ! (and no only minimum).
       converge = converge .and. smin < s1 .and. s1 < smax
       if( converge  ) s = s1
    else

       ! if minimum is not found, recovers by backup solution
       ! this is common for a few datapoints with strongly non-Normal
       ! distribution
       converge = .false.
       s = s0
    end if

    if( present(flag) ) then
       if( converge ) then
          flag = 0
       else
          flag = 1
       end if
    end if
    if( present(reliable) ) reliable = .true.

    deallocate(rho)

  end subroutine escale

  function negentropy(s)

    ! negative of entropy per data point

    real(kind) :: negentropy
    real(kind), intent(in) :: s

    rho = itukey(res/s)
    negentropy = - sum(rho*exp(-A*rho), rho < rhomax) / size(rho)
    ! The sum includes only values within range acceptable for exp() function.
    ! Otherwise, FPE flag indicating underflow will be raised, although
    ! validity of results is untouched.

    if( debug_files ) write(7,*) s,negentropy

  end function negentropy

  subroutine newton(s,tol,converge)

    real(kind), intent(in) :: tol
    real(kind), intent(in out) :: s
    logical, intent(out) :: converge

    real(kind), dimension(:), allocatable :: psi, dpsi, erho
    integer :: n,it,stat
    real(kind) :: d,f,df,sum1,sum2
    character(len=80) :: errmsg

    converge = .false.

    n = size(res)
    allocate(psi(n),dpsi(n),erho(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,*) 'Netwon scale: ',trim(errmsg)
       return
    end if
    ! Newton needs to allocate additional memory.
    ! If the memory is unavailable, we exiting and fail-back
    ! solution must be used.

    ! Newton's iterations
    do it = 1, precision(s)

       if( .not. (s > epsilon(s)) ) return

       call tukeys(res/s,psi,dpsi,rho)

       where( rho < rhomax )
          erho = exp(-A*rho)
       elsewhere
          erho = 0.0_kind
       end where
       sum1 = sum(psi*(1.0_kind - A*rho)*res*erho)
       sum2 = sum((A*psi**2 + (1.0_kind - A*rho)*(A*psi**2 - dpsi))*res**2*erho)
       f = sum1
       df = (sum2 / s - 2.0_kind*sum1) / s

       if( .not. (abs(df) > epsilon(df) .and. -df < 0) ) exit

       ! corrector for mean
       d = f / df

       ! update location
       s = s - d

       if( debug ) write(*,'(a,i2,1pg20.6,1p3g10.2)') &
            "#, s, inc., f, f': ",it,s,d,-f/s**2,-df/s**2

       ! exit of iterations: the absolute errot must be at least |d| < tol
       converge = abs(d) < tol
       if( converge ) exit

    enddo

    ! second derivation of neg-entropy is positive for a stable local minimum
    converge = converge .and. -df < 0

  end subroutine newton

  subroutine grid(smin,smax,s,found)

    ! Grid search
    !
    ! The search looks for minimum on a grid within an interval
    ! of smin .. smax. Grid search examines all possible minimas
    ! of negentropy(). Every minimum reveals a group of data subset
    ! with similar values. The most probable solutions corresponds
    ! to the most numerous group.

    real(kind), intent(in) :: smin, smax
    real(kind), intent(out) :: s
    logical, intent(out) :: found

    ! grid elements
    real(kind), dimension(:), allocatable :: ord, tab
    real(kind) :: h,p,pmax
    integer :: n,ngrid,cmin,nmin

    ! develop grid
    ! common smin .. smax interval corresponds to 2-3 sigma,
    ! step is at least 1/10 of sigma
    ngrid = 23

    allocate(ord(ngrid+2),tab(ngrid+2))
    h = (smax - smin) / (ngrid+2)

    ! the first point is close to zero (but non-zero!)
    do n = 1, ngrid+2
       s = smin + (n-1)*h
       ord(n) = s
       tab(n) = negentropy(s)
       !          if( debug ) write(7,*) ord(n),tab(n) in negentropy....
    end do

    ! locate minimum
    ! minimal limit of probability, meaning that the maximum
    ! corresponds to group of over 50% of datapoints
    pmax = 0.5
    cmin = 0
    nmin = 0
    do n = 2, ngrid+1
       if( tab(n-1) > tab(n) .and. tab(n) < tab(n+1) ) then
          cmin = cmin + 1

          p = sum(dtukey(res/ord(n))) / size(res)
          if( p > pmax ) then
             nmin = n
             pmax = p
          end if

          if( debug ) write(*,'(a,i3,g15.3,f10.2)') &
               'minloc:',n,ord(n),p
       end if

    end do
    if( nmin > 0 ) then
       s = ord(nmin)
    else
       s = 0
    end if

    deallocate(tab,ord)

    ! found minimum with required data amount
    found = 2 <= nmin .and. nmin <= ngrid+1

    if( debug ) write(*,'(3(a,g10.3),2(a,i0,1x),l1)') &
         'Grid search s: ',smin," < ",s," < ",smax," n:",ngrid,' #',cmin,found

  end subroutine grid

  subroutine gmin(smin,smax,s,found)

    use fmm

    real(kind), intent(in) :: smin, smax
    real(kind), intent(out) :: s
    logical, intent(out) :: found

    real(kind), parameter :: e = 10*epsilon(s)

    ! Routine fmin localises a minimum by using of golden-ratio search,
    ! and parabola fit close to the extreme. Scale s, as a statistical
    ! parameter, can be known with low precision (2 significant digits
    ! are enought for our application), giving ln((smax - smin)/tol) < 7
    ! calls to negentropy().
    !
    ! The fmin has been choiced in precedence of zeroin() (root finding
    ! method), we are not able to distiguish between both maximum or
    ! minimum without derivatives

    s = fmin(smin,smax,negentropy,etol)

    ! Does the found minimum falls inside the interval smin < .. < smax ?
    found = smin+e < s .and. s < smax-e

    if( debug ) write(*,'(a,g20.6,a,f0.1,l2)') &
         'Minimum fmin location:',s," +- ",etol,found

  end subroutine gmin


  ! diagnostic, neg-entropy graph
  subroutine escale_graph(filename,r,smin,smax)

    character(len=*), intent(in) :: filename
    real(kind), dimension(:), target, intent(in) :: r
    real(kind), intent(in) :: smin, smax

    ! local rho overwrites global rho (!)
    real(kind), dimension(:), allocatable :: rho,psi, dpsi, erho
    integer, parameter :: m = 666
    integer :: i,n
    real(kind) :: s,ds,f,df,sum1,sum2

    n = size(res)

    if( n <= 0 ) return

    allocate(rho(n),psi(n),dpsi(n),erho(n))

    ds = (smax - smin) / m

    open(1,file=filename)
    res => r
    do i = 0,m
       s = smin + i * ds

       call tukeys(res/s,psi,dpsi,rho)
       where( rho < rhomax )
          erho = exp(-A*rho)
       elsewhere
          erho = 0.0_kind
       end where
       sum1 = sum(psi*(1.0_kind - A*rho)*res*erho)
       sum2 = sum((A*psi**2 + (1.0_kind - A*rho)*(A*psi**2 - dpsi))*res**2*erho)
       f = sum1
       df = (sum2 / s - 2.0_kind*sum1) / s

       write(1,*) s,negentropy(s),f,df,-sum((r/s)**2*exp(-(r/s)**2))/size(r)
    end do
    close(1)

  end subroutine escale_graph


end module entropyscale
