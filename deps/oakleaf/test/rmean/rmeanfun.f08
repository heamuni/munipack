
program rmeanfun

  ! test an optimization

  use rfun
  use sr1mean_function
  use robnewton
  use quantilmean
  use noise
  use iso_fortran_env
  implicit none

  type(tsr1mean_tukey), pointer :: tukeymean
  type(tsr1mean_experimental), pointer :: bathmean
  class(tsr1mean_function), pointer :: funmean

  integer, parameter :: datamax = 5
  logical, parameter :: switch_on_generator = .true.
  character(len=666), parameter :: datafile = '/tmp/rmean_datafail'
  real(REAL64), dimension(:), allocatable, target :: data
  real(REAL64) :: frac, t, s, tmed, xmin, xmax
  integer :: i, step, ndat, iflag

  ! random generated data
  if( switch_on_generator ) then
     allocate(data(datamax))
     call gnoise_fill(0.0_REAL64,1.0_REAL64,data)
     frac = 0.3
     step = datamax / max(int(frac*datamax),1)
     step = 6
     write(*,*) 'step:',frac*datamax,step
     do i = 2, datamax, step
        data(i) = gnoise(-12.0_REAL64,1.3_REAL64)
        write(*,*) 'bad point:',real(data(i))
     end do
     write(*,*) 'data:',real(data)
  else
     ! saved data
     open(3,file=datafile,form='unformatted')
     read(3) ndat
     allocate(data(ndat))
     read(3) data
     close(3)
     write(*,*) 'data:',real(data)
  end if


!  data = [ 142, 141, 149, 149, 142, 148,149 ]
!  write(*,*) sum(data)/7

!  data = [2.25176118E-05, 4.98100162E-05, 2.31736030E-05,   9.98998730E-05]
  !  write(*,*) sum(data)/size(data)
!!$  data = [ 10109.000000000000 ,  10109.000000000000 ,   10109.000000000000 , &
!!$       10109.000000000000 ,   10109.000000000000  ,      10109.000000000000  , &
!!$       10109.000000000000 ,   10109.000000000000   ,     10109.000000000000  , &
!!$       729415.00000000000 ,   10109.000000000000  ,      10109.000000000000 , &
!!$       10109.000000000000 ,   10109.000000000000  ,      10109.000000000000 , &
!!$       10109.000000000000 ]
!!$  write(*,*) sum(data)/size(data)

  allocate(bathmean,tukeymean)
  call bathmean%init(data,iflag)
  call tukeymean%init(data,iflag)

!  tukeymean%data => data
!  bathmean%data => data

  block
    real(REAL64), dimension(5), parameter :: q = [0.2, 0.25, 0.5, 0.75, 0.8]
    real(REAL64), dimension(size(q)) :: f
    real(REAL64), dimension(:), allocatable :: x,y
    real(REAL64) :: tmin,tmax, s1, s2,h,q1,q2
    integer :: flag,m1,m2,n,i

    call quantilefun(data,q,f,flag)
    if( flag /= 0 ) stop

    !t = (f(3) + (f(4) + f(2))/2)/2
    t = f(3)
    tmin = f(2)
    tmax = f(4)
    s = (f(4) - f(2))/(2*0.6745)
!    s1 = f(4) / 0.6745
!    s2 = f(2) / 0.6745
    call quantilefun(abs(data-t),q,f,flag)
    s = f(3) / 0.6745
    write(*,*) s,s1,-s2
!    s = 1
    tmin = t - s
    tmax = t + s

    open(1,file='/tmp/rmeanfun_quant')
    write(1,*) t,-2,tmin,-2,tmax,-2
    write(1,*) t,2,tmin,2,tmax,2
    close(1)

    n = size(data)
    allocate(x(n),y(n))
    call ecdf(data,x,y)

    do i = 1, n
       write(*,*) x(i),y(i)
    end do

    if( mod(n,2) == 0 ) then
       t = (x(n/2) + x(n/2+1))/2
       h = 1 / real(n)
    else
       h = 1 / real(n+1)
       t = x(int((n+1)/2))
    end if

    m1 = n + 1 - int(0.75 / h)
    m2 = int(0.75 / h)
    m2 = int(0.75 * (n+1))
    m1 = n + 1 - m2
    q1 = m1*h
    q2 = m2*h
    write(*,*) q1,q2,m1,m2

!    s = (x(m2) - x(m1)) / (2*0.4125)
    s = (x(m2) - x(m1)) / (2*0.6745)

    tmed = t
    write(*,*) 'median:',t,s

    tmin = t - s
    tmax = t + s

    open(1,file='/tmp/rmeanfun_quant')
    write(1,*) t,-2,tmin,-2,tmax,-2
    write(1,*) t,2,tmin,2,tmax,2
    close(1)

  end block

  block
    real(REAL64) :: dt, tol
    logical :: converge
    integer :: iter,trstep

    tol = 1e-2/sqrt(real(datamax,REAL64))

    !    t = 5
    t = tmed
!    s = 1
    call rnewton(bathmean,t,dt,s,tol,converge,iter)
    write(*,*) 'rnewton hampel:',t,dt,converge

    t = tmed
!    t = 5
!    s = 1
    call trnewton(bathmean,t,dt,s,tol,converge,iter,trstep)
    write(*,*) 'trnewton hampel:',t,dt,converge

    t = tmed
    call rnewton(tukeymean,t,dt,s,tol,converge,iter)
    write(*,*) 'rnewton tukey:',t,dt,converge
  end block

  xmin = -20
  xmax = 20
  s = 0.5
  call graph(bathmean,s,xmin,xmax,'/tmp/rmeanfun')
  call graph(tukeymean,s,xmin,xmax,'/tmp/rmeanfun_tukey')
  call ecdf_graph(data,xmin,xmax,'/tmp/rmean_ecdf')
  call ecdf_graph(abs(data-tmed),xmin,xmax,'/tmp/rmean_ecdf_mad')
  !xmin = 130
  !xmax = 160
  !call graph(bathmean,s,xmin,xmax,'/tmp/trmean1d')
  !call graph(tukeymean,s,xmin,xmax,'/tmp/trmean1d_tukey')

  deallocate(bathmean,tukeymean)

contains

  subroutine graph(fcn,s,x0,x1,filename)

    class(tsr1mean_function), intent(in) :: fcn
    real(REAL64), intent(in) :: s,x0,x1
    character(len=*), intent(in) :: filename
    integer, parameter :: nx = 200
    real(REAL64) :: x,dx
    real(REAL64), dimension(2) :: g
    real(REAL64), dimension(2,2) :: hess
    integer :: i

    dx = (x1 - x0) / (nx - 1)

    open(1,file=filename)
    do i = 0, nx
       x = x0 + dx*i
       call fcn%grad(x,s,g)
       call fcn%hess(x,s,hess)
       write(1,*) x,fcn%val(x,s),g(1),hess(1,1)
    end do
    close(1)

  end subroutine graph

  subroutine ecdf_graph(data,x0,x1,filename)

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: x0,x1
    character(len=*), intent(in) :: filename
    real(REAL64), dimension(:), allocatable :: x,y
    integer :: i,n

    n = size(data)
    allocate(x(n),y(n))

    call ecdf(data,x,y)
    open(1,file=filename)
    do i = 1, n
       write(1,*) x(i),y(i)
    end do
    close(1)

  end subroutine ecdf_graph


end program rmeanfun
