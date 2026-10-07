
module robnewton

  use sr1mean_function
  use iso_fortran_env
  implicit none

  logical, parameter, private :: debug = .false.
!  integer, dimension(precision(real(1.0,REAL64))) :: riter, triter

contains

  subroutine trnewton(fcn,t,dt,s,tol,converge,iters,trstep)

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

!    use rfun
!    use robust_data_function

!    class(robust_datafun), intent(in) :: fcn
    class(tsr1mean_function), intent(in) :: fcn
    !   real(kind), dimension(:), intent(in) :: x
!    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in out) :: t
    real(REAL64), intent(out) :: dt
    real(REAL64), intent(in) :: s, tol
    logical, intent(out) :: converge
    integer, intent(out) :: iters, trstep

!!$    interface
!!$       elemental pure function psi(x)
!!$         import :: REAL64
!!$         real(REAL64) :: psi
!!$         real(REAL64), intent(in) :: x
!!$       end function psi
!!$       elemental pure function dpsi(x)
!!$         import :: REAL64
!!$         real(REAL64) :: dpsi
!!$         real(REAL64), intent(in) :: x
!!$       end function dpsi
!!$    end interface


!    real(REAL64), dimension(:), allocatable :: psi,dpsi
!    character(len=81) :: errmsg
    integer :: n,iter!,stat
    real(REAL64) :: d,f,f2,df,g1,h11,dlt
    logical :: posdef
!    real(REAL64), dimension(2) :: g
!    real(REAL64), dimension(2,2) :: hess

    converge = .false.
    if( .not. (s > epsilon(s)) ) return
    dlt = s / 4
    trstep = 0
!    n = size(fcn%data)
!    n = size(data)

!    allocate(psi(n),dpsi(n),stat=stat,errmsg=errmsg)
!    allocate(res(n),stat=stat,errmsg=errmsg)
!    if( stat /= 0 ) then
!       write(error_unit,*) 'Newton rmean: ',trim(errmsg)
!       return
!    end if

    ! set-up required precision
!    tol = max(1e-3/sqrt(real(n,REAL64)),10*n*epsilon(t)*(abs(t) + 1))
!    write(*,*) 'tol=',tol

    ! Newton's iterations
    do iter = 1, 3*precision(t)
       ! Number of iterations is given by numerical precision.
       ! We belives, two orders are reached by each iteration, at least.
       iters = iter

       !       call tukeys((x-t)/s,psi,dpsi)
!       call fcn%grad(t,s,g)
       !       call fcn%hess(t,s,hess)
       !      call fcn%line(t,s,g1,h11)

!       call fcn%(t,s,g1,h11)
       !res = (data - t) / s
       !f = sum(psi(res))
       !df = sum(dpsi(res))
!       g1 = fcn%g1(t,s)
       !       h11 = fcn%h11(t,s)
       call fcn%rnewton(t,s,g1,h11)
       f = g1!g(1)
       df = h11
!       df = hess(1,1)

       ! |df| will be vanishing when the important part of data has large
       ! deviations (above 3.4 for Hampel's piecewise or 1.6 for Tukey's
       ! bi-weight), this si "re-descending M-estimate" function problem.
       ! In this case, the computation is immediatelly interrupted and
       ! result should by considered as unreliable.
       ! The test is conservative, tiny() should be fine as well as.
       !       if( .not. (df > epsilon(df)) ) exit
       !       if( hess(1,1) > 0 ) then
       posdef = h11 > 0
       if( posdef ) then
!          df = hess(1,1)
          df = h11
          d = -f / df
          d = -g1 / h11
          if( abs(d) > dlt ) then
             d = -sign(dlt,g1)
             trstep = trstep + 1
          end if
       else
          !          df = 1/s**2
          d = -sign(dlt,g1)
          trstep = trstep + 1
          !          d = dlt
       end if


       ! corrector
!       d = f / df

!       d = -f / df
!       if( abs(d) > s ) d = -sign(s,f)

       ! update location
       t = t + d

       if( debug ) &
            write(*,'(a,i2,1pg20.10,1p3g10.2,L2)') &
            "mean, inc., f, f': ",iter,t,d,-f/s,df/s**2,abs(abs(d)-dlt)<epsilon(dlt)

       ! The exit of iterations is immediately done when updating
       ! takes less than 1/1000 of standard error indicating
       ! no further progress in iterations will give even precise result.
       ! As the fallback, the numerical tolerance is used.
       ! No convergence is indicated otherwise.

       if( abs(d) < tol ) then
!       if( abs(d) < tol .and. posdef ) then
!          f2 = sum(psi**2)
          !          f2 = sum(psi(res)**2)
          if( posdef ) then
             dt = fcn%stderr(t,s,converge)
          end if
!       f2 = dt
!       if( f2 > 0 .and. df > 0) then  ! condition df > 0 is valid already
          ! estimation of standard error
!          dt = s*sqrt(f2/df**2*n/(n-1)) ! the order does matter
          !          converge = iflag == 0
          converge = .true.
          exit
 !         dt = fcn%stderr(t,s)
          !         converge = abs(d) < max(0.001*dt,10*epsilon(d)) .and. abs(dt) > 10*epsilon(dt)
          !converge = abs(dt) > 10*epsilon(dt)
!       else ! f2 == 0, t is out of global minimum
!          dt = 0
          ! flag = 2
          !          if( abs(d) < tol ) exit
!          exit
       end if
!    end if

!       if( converge ) exit
    enddo

  end subroutine trnewton


    !-----------------------------------------------------------------------

  subroutine rnewton(fcn,t,dt,s,tol,converge,iters)

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

    class(tsr1mean_function), intent(in) :: fcn
!   real(kind), dimension(:), intent(in) :: x
    real(REAL64), intent(in out) :: t
    real(REAL64), intent(out) :: dt
    real(REAL64), intent(in) :: s, tol
    logical, intent(out) :: converge
    integer, intent(out) :: iters

!    real(REAL64), dimension(:), allocatable :: psi,dpsi
    character(len=81) :: errmsg
    integer :: n,iter,stat
    real(REAL64) :: d,f,f2,df
    real(REAL64), dimension(2) :: g
    real(REAL64), dimension(2,2) :: hess

    converge = .false.
    if( .not. (s > epsilon(s)) ) return

    n = size(fcn%data)

!    allocate(psi(n),dpsi(n),stat=stat,errmsg=errmsg)
!    if( stat /= 0 ) then
!       write(error_unit,*) 'Newton rmean: ',trim(errmsg)
!       return
!    end if

    ! set-up required precision
!    tol = 10*n*epsilon(t)*(abs(t) + 1)
!    tol = max(1e-3/sqrt(real(n,REAL64)),10*n*epsilon(t)*(abs(t) + 1))

    ! Newton's iterations
    do iter = 1, precision(t)
       ! Number of iterations is given by numerical precision.
       ! We belives, two orders are reached by each iteration, at least.

       !       call tukeys((x-t)/s,psi,dpsi)
       call fcn%grad(t,s,g)
       call fcn%hess(t,s,hess)
!       f = sum(psi)
!       df = sum(dpsi)
       f = g(1)
       df = hess(1,1)

       ! |df| will be vanishing when the important part of data has large
       ! deviations (above 3.4 for Hampel's piecewise or 1.6 for Tukey's
       ! bi-weight), this si "re-descending M-estimate" function problem.
       ! In this case, the computation is immediatelly interrupted and
       ! result should by considered as unreliable.
       ! The test is conservative, tiny() should be fine as well as.
       if( .not. (df > epsilon(df)) ) exit

       ! corrector
       d = f / df

       ! update location
       t = t - d

       if( debug ) &
            write(*,'(a,i2,1pg20.10,1p3g10.2)') &
            "mean, inc., f, f': ",iter,t,d,-f/s,df/s**2

       ! The exit of iterations is immediately done when updating
       ! takes less than 1/1000 of standard error indicating
       ! no further progress in iterations will give even precise result.
       ! As the fallback, the numerical tolerance is used.
       ! No convergence is indicated otherwise.

       if( abs(d) < tol ) then
!       f2 = sum(psi**2)
!       f2 = 1
!       if( f2 > 0 ) then  ! condition df > 0 is valid already
          ! estimation of standard error
          !          dt = s*sqrt(f2/df**2*n/(n-1)) ! the order does matter
          dt = fcn%stderr(t,s,converge)
          !          riter(iter) = riter(iter) + 1
          iters = iter
          exit
          !          converge = abs(d) < max(0.001*dt,10*epsilon(d)) .and. abs(dt) > 10*epsilon(dt)

!       else ! f2 == 0, t is out of interval -c .. c of Tukey's function
!          exit
!          dt = 0
!          if( abs(d) < tol ) exit
       end if

!       if( converge ) exit
    enddo

  end subroutine rnewton


end module robnewton
