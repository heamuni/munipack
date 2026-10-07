!
!  Information scale:
!  scale by maximal likelihood, maximum of information, minimum of variance
!
!  Copyright © 2018 F.Hroch (hroch@physics.muni.cz)
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

module likescale

  ! This module provides the subroutine for estimation of scale
  ! by maximising of information, or minimising of variance (dispersion):
  !    https://en.wikipedia.org/wiki/Fisher_information
  !
  ! iscale should be called as:
  !
  !  call iscale(r,s,reliable,flag,verbose)
  !
  ! on input:
  !   r - array of residuals: (data-mean), (data-mean)/errors, ...
  !   verbose (optional) - verbose prints
  !
  ! on output are estimated:
  !   s - scale
  !   reliable (optional) - indicates reliability of result
  !                         if .true., scale has been correctly
  !                         localised. if .false., s is undefined
  !   flag (optional) - result code:
  !          0 == success, results are both reliable and precise
  !          1 == consider rough estimate due convergence fail, try verbose
  !          2 == basic estimates, few datapoints available
  !          3 == data looks nearly identical
  !          5 == failed to allocate memory, initial estimates are returned
  !
  ! Note.
  ! Scale parameter shouldn't confused with standard deviation.
  ! Their values are identical only in case of Normal distribution.
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

  ! data buffer
  real(kind), dimension(:), allocatable, private :: rho

  ! data pointer, just pointer
  real(kind), dimension(:), pointer, private :: res

  ! likelihood computations counter, debug purposes
  integer, private :: countlike

  private :: loglike, newton

  public :: iscale

contains

  subroutine iscale(r,s,reliable,flag,verbose)

    use quantilmean
    use fmm

    real(kind), dimension(:), target, intent(in) :: r
    real(kind), intent(out) :: s
    integer, intent(out), optional :: flag
    logical, intent(out), optional :: reliable
    logical, intent(in), optional :: verbose

    real(kind), parameter :: qhalf = 0.6745, qmin = 0.4, qmax = 0.87
    real(kind), dimension(:), allocatable :: xcdf, ycdf
    character(len=80) :: errmsg
    real(kind) :: s0,s1,tol,smin,smax
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

    if( debug ) write(*,'(a,2g13.3)') &
         'Maximum and minimum: ',minval(abs(res)),maxval(abs(res))

    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'scale: ',trim(errmsg)
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

    ! the very first, rought, estimate of scale at Q(0.6745)
    ! qhalf = 0.6745
    s0 = quantile(qhalf,xcdf,ycdf)

    ! the intial interval is determined by qmin and qmax quantiles
    ! qmin = 0.4    ! = 0.5*s
    ! qmax = 0.87   ! = 1.5*s

    ! quantile for Q(0.87) for 1.5*s0, gives upper bound
    if( qmax < ycdf(n) ) then
       smax = quantile(qmax,xcdf,ycdf)
    else ! if( Q(0.87) is out of CDF values ) then
       smax = 1.5 * s0  ! = Q(4/5) * s0
    end if

    ! quantile for Q(0.4), ycdf(1) = 1/n
    if( qmin > ycdf(1) ) then
       smin = quantile(qmin,xcdf,ycdf)
    else ! if( Q(0.4) < 1/n, interpolation undefined ) then
       smin = 0.5 * s0  ! = Q(2/5) * s0
    end if

    deallocate(xcdf,ycdf)

    if( debug ) write(*,'(a,2g13.3,a,g13.3)') &
         'Quantiles bracketing:',smin,smax,'   Initial scale:',s0

    ! check: smax > smin, s,s0 > 0, identical or zero value data
    if( (smax - smin) < 10*epsilon(s) .or. abs(smax) < 10*epsilon(s) .or. &
         abs(s0) < 10*epsilon(s) ) then
       s = s0
       if( present(flag) ) flag = 3
       if( present(reliable) ) reliable = s > epsilon(s)
       return
    end if

    allocate(rho(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'scale: ',trim(errmsg)
       if( present(flag) ) flag = 5
       if( present(reliable) ) reliable = .false.
       return
    end if

    ! Optimum determination
    ! The extreme is searched by an optimisation method without derivatives,
    ! it is reliable but slower than follow-up Newton method.
    if( debug_files ) then
       countlike = 0
       open(7,file='/tmp/s')
    end if
    tol = s0*etol
    s = fmin(smin,smax,loglike,tol)
    if( debug_files ) close(7)

    ! Does the minimum falls inside the interval smin < .. < smax ?
    found = smin + tol < s .and. s < smax - tol

    if( debug ) write(*,'(a,g20.6,a,g0.2,l2,i5)') &
         'Minimum fmin location:',s," +- ",tol,found,countlike

    if( found ) then
       ! We are attemping to do 10x more precise estimate of the scale
       ! by fast Newtons's method. The convergence success is mostly given
       ! by the proper starting point.
       s1 = s
       tol = 0.1*tol
       call newton(s1,tol,converge)
       ! Newton method is kind of root-finding metod. We are looking for
       ! maximum of a function so we work with derivatives of the function.
       ! By second derivative, we can choice if it is truly maximum
       ! (and no only minimum).
       converge = converge .and. smin + tol < s1 .and. s1 < smax - tol
       if( converge  ) s = s1
    else

       ! if minimum is not found, backup solution is used instead,
       ! it's common for few datapoints and strongly non-Normal distribution
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

  end subroutine iscale

  function loglike(s)

    ! negative of likelihood per data point

    real(kind) :: loglike
    real(kind), intent(in) :: s

    if( .not. s > epsilon(s) ) then
       loglike = min(1e3*maxval(abs(res))**2,0.1*huge(loglike))
       return
    end if

    rho = itukey(res/s)
    loglike = sum(rho)/size(rho) + log(s)

    if( debug_files ) then
       countlike = countlike + 1
       write(7,*) s,loglike
    end if

  end function loglike


  subroutine newton(s,tol,converge)

    real(kind), intent(in) :: tol
    real(kind), intent(in out) :: s
    logical, intent(out) :: converge

    real(kind), dimension(:), allocatable :: psi, dpsi, rs
    integer :: n,it,stat
    real(kind) :: d,f,df
    character(len=80) :: errmsg

    converge = .false.

    n = size(res)
    allocate(psi(n),dpsi(n),rs(n),stat=stat,errmsg=errmsg)
    ! Computations needs a lot of additional memory,
    ! a fail-back way should be used if it is unavailable.
    if( stat /= 0 ) then
       if( debug ) write(error_unit,*) 'Newton scale: ',trim(errmsg)
       return
    end if

    ! Newton's iterations
    do it = 1, precision(s)

       if( .not. (s > epsilon(s)) ) return

       rs = res / s
       call tukeys(rs,psi,dpsi)

       f = sum(psi*rs)/n - 1
       df = - f - sum((dpsi*rs + psi)*rs)/(n*s)

       if( .not. (df < -epsilon(df)) ) exit

       ! update scale
       d = f / df
       s = s - d

       if( debug ) write(*,'(a,i2,1pg20.6,1p3g10.2)') &
            "#, s, inc., f, f': ",it,s,d,-f/s,-df/s

       ! exit of iterations: the absolute error must be at least |d| < tol
       converge = abs(d) < tol
       if( converge ) exit

    enddo

    ! second derivation of likelog is negative for a stable local minimum
    converge = converge .and. df < 0

  end subroutine newton

end module likescale
