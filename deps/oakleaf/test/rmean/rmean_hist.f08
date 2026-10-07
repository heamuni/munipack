
program rmean_hist

  use toolbox
  use iso_fortran_env
  use robnewton
  implicit none

  integer, parameter :: ndata = 10
  integer, parameter :: nsamples = 999999
  integer, parameter :: step = 12
  real(REAL64), dimension(:), allocatable :: t,stderr,stdsig
  integer :: i, rfail, trfail, nfail

  allocate(t(nsamples),stderr(nsamples),stdsig(nsamples))

!  triter = 0
!  riter = 0
  rfail = 0
  trfail = 0
  nfail = 0
  do i = 1, nsamples
     call one_round(t(i),stderr(i),stdsig(i))
  end do

  call histogram(t,real(-3,REAL64),real(3,REAL64),150,'/tmp/rmean_t',.true.)
  call histogram(stderr,real(0,REAL64),real(5,REAL64),150,'/tmp/rmean_stderr',.true.)
  call histogram(stdsig,real(0,REAL64),real(5,REAL64),150,'/tmp/rmean_stdsig',.true.)

  open(1,file='/tmp/rmean_hist.par')
  write(1,*) ndata
  close(1)

  write(*,*) 'Convergence of rnewton failed in ',rfail,' case(s).'
  write(*,*) 'Convergence of trnewton failed in ',trfail,' case(s).'
  write(*,*) 'Results inconsistent in ',nfail,' case(s).'
!  do i = 1, size(riter)
!     write(*,*) i, triter(i), riter(i)
!     write(*,'(16I3)') triter
!     !     write(*,'(16I3)') riter
!  end do


contains

  subroutine one_round(t,dt,sig)

    use quantilmean
    use sr1mean_function
    use noise

    real(REAL64), intent(out) :: t,dt,sig
    !    type(sr1_mean) :: funmean
    type(tsr1mean_tukey), target :: tukeymean
    type(tsr1mean_huber), target :: hubermean
!    type(tsr1mean_square), target :: squaremean
    type(tsr1mean_experimental), target :: bathmean
    class(tsr1mean_function), pointer :: funmean

    real(REAL64), dimension(:), allocatable, target :: data
    integer :: i,iflag,iters,trstep
    real(REAL64) :: s, tmad, tol, t1, t2
    logical :: converge,reli

    allocate(data(ndata))
    call gnoise_fill(0.0_REAL64,1.0_REAL64,data)
    do i = 2, ndata, step
       data(i) = gnoise(-15.0_REAL64,1.0_REAL64)
!       write(*,*) data(i)
    end do

    call qmean(data,t,dt,s,.true.,reli,iflag)
    tmad = t
!    s = 1
    tol = 1e-2*s/sqrt(real(size(data)))

    call bathmean%init(data,iflag)
    call tukeymean%init(data,iflag)
    funmean => bathmean
!    funmean => tukeymean
!    funmean => hubermean
    funmean%data => data

    tmad = 0
    s = 1
    tol = 1e-2/sqrt(real(size(data)))
    t = tmad
    call rnewton(funmean,t,dt,s,tol,converge,iters)
    t1 = t
    if( .not. converge ) rfail = rfail + 1

    t = tmad
    call trnewton(funmean,t,dt,s,tol,converge,iters,trstep)
    t2 = t
    if( .not. converge ) trfail = trfail + 1

    if( converge .and. abs(t1 - t2) > 10*tol ) then
       nfail = nfail + 1
       write(*,*) t1,t2
!       trfail = trfail + 1
!       rfail = rfail + 1
    end if

    sig = sqrt(real(size(data)))*dt

    converge = .false.
    if( .not. converge ) then
!       nfail = nfail + 1

       ! save data
       open(3,file='/tmp/rmean_datafail',form='unformatted')
       write(3) size(data)
       write(3) data
       close(3)
!       stop

    end if

  end subroutine one_round

!!$  subroutine rinit(data,t,mad,s)
!!$
!!$    use quantilmean
!!$
!!$    real(REAL64), dimension(:), intent(in) :: data
!!$    real(REAL64), intent(out) :: t,mad,s
!!$
!!$    real(REAL64), dimension(3), parameter :: q = [0.25, 0.5, 0.75]
!!$    real(REAL64), dimension(size(q)) :: f
!!$    real(REAL64), dimension(:), allocatable :: x,y
!!$    real(REAL64) :: h
!!$    integer :: flag,n,i,m1,m2,r,m
!!$
!!$    n = size(data)
!!$    allocate(x(n),y(n))
!!$    call ecdf(data,x,y)
!!$
!!$    call quantilefun(data,q,f,flag)
!!$    t = f(2)
!!$!    t = (f(2) + (f(1) + f(3))/2)/2
!!$    s = (f(3) - f(1)) / (2*0.6745)
!!$    mad = s*0.6745
!!$
!!$    !    s = f(2) / 0.6745
!!$
!!$    do i = 1, n
!!$       write(*,*) x(i),y(i),abs(x(i)-t)
!!$    end do
!!$
!!$    if( mod(n,2) == 0 ) then
!!$       t = (x(n/2) + x(n/2+1))/2
!!$       h = 1 / real(n)
!!$!       m = real(n + 1)/2.0
!!$    else
!!$       h = 1 / real(n+1)
!!$!       m = real(n + 1) / 2
!!$       t = x(int((n+1)/2))
!!$    end if
!!$!    m = nint(0.5/h)
!!$    !    t = x(m)
!!$    write(*,*) h,m,n,mod(n,2)
!!$
!!$    h = 1 / real(n+1)
!!$    m = nint(0.5 / h)
!!$    r = int(0.25 / h)
!!$    if( h*r < 0.25 ) r = r + 1
!!$    m1 = n + 1 - int(0.75 / h)
!!$!    m1 = r !int(0.25 / h + 0.5)
!!$    write(*,*) r,0.25/h
!!$    r = int(0.75 / h)
!!$    m2 = r
!!$
!!$!    r = int(0.25 / h)
!!$!    m1 = int(m - r + 0.5)
!!$!    m2 = int(m + r - 0.5)
!!$
!!$!    m1 = nint(0.33/h)
!!$    !    m2 = nint(0.66/h)
!!$
!!$
!!$    s = (x(m2) - x(m1)) / (2*0.4125)
!!$    s = (x(m2) - x(m1)) / (2*0.6745)
!!$    s = (x(m2) - x(m1)) / (2*0.53)
!!$!    mad = (m2-m1)*h*(x(m2) - x(m1))
!!$
!!$    goto 33
!!$    call ecdf(abs(data-t),x,y)
!!$    if( mod(n,2) == 0 ) then
!!$       s = (x(n/2) + x(n/2+1))/(2*0.6745)
!!$       write(*,*) (x(n/2) + x(n/2+1))/2
!!$    else
!!$       s = x(int((n+1)/2)) / 0.6745
!!$       write(*,*) x(int((n+1)/2))
!!$    end if
!!$33 continue
!!$!    write(*,*) t,m,m1,m2,x(m1),x(m2),mad,(m2-m1)*h
!!$    if( mad < 10*epsilon(mad) ) stop
!!$!    stop
!!$
!!$  end subroutine rinit

end program rmean_hist
