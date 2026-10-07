
program rmean_region

  ! test an optimization

  use rfun
  use sr1mean_function
  use robnewton
  use quantilmean
  use noise
  use iso_fortran_env
  implicit none

!  type(robust_datafun) :: funmean
  type(tsr1mean_tukey), pointer :: tukeymean
  type(tsr1mean_experimental), pointer :: bathmean
  class(tsr1mean_function), pointer :: funmean

  character(len=666), parameter :: datafile = '/tmp/rmean_datafail'
  real(REAL64), dimension(:), allocatable, target :: data
  real(REAL64), dimension(:), allocatable :: t
  real(REAL64) :: tt, s, tmed, xmin, xmax, dt, ddt, tmin, tmax, tol
  logical, dimension(:), allocatable :: cnewton, ctrnewton
  integer :: i, n, step, ndat, iflag, trstep, rstep, iters
  logical :: reli

  ! saved data
  open(3,file=datafile,form='unformatted')
  read(3) ndat
  allocate(data(ndat))
  read(3) data
  close(3)
  write(*,*) 'data:',real(data)

  call qmean(data,tmed,dt,s,.true.,reli,iflag)

  open(1,file='/tmp/rmean_region_quant')
  write(1,*) tmed,0,s
  close(1)

  allocate(bathmean,tukeymean)
  call bathmean%init(data,iflag)
  call tukeymean%init(data,iflag)
  if( iflag /= 0 ) stop
  funmean => bathmean
  funmean => tukeymean

  !  bathmean%data => data
!  funmean%rho => ipshampel
!  funmean%psi => pshampel
!  funmean%dpsi => dpshampel
!  funmean%data => data

  tmin = -10
  tmax = 10
  dt = 0.1
  n = int((tmax - tmin) / dt)
  allocate(cnewton(n),ctrnewton(n),t(n))
  write(*,*) n

  tol = 1e-3/sqrt(real(ndat,REAL64))

  do i = 1, n
     t(i) = tmin + (i-1)*dt

     tt = t(i)
     call rnewton(funmean,tt,ddt,s,tol,cnewton(i),rstep)
     write(*,*) 'rnewton:',tt,ddt,cnewton(i)
     cnewton(i) = cnewton(i) .and. abs(tt - 0.33) <  0.2

     tt = t(i)
     call trnewton(funmean,tt,ddt,s,tol,ctrnewton(i),iters,trstep)
     write(*,*) 'trnewton:',tt,ddt,ctrnewton(i),iters,trstep
     ctrnewton(i) = ctrnewton(i) .and. abs(tt - 0.33) <  0.2

  end do

  open(1,file='/tmp/rmean_region_conv')
  do i = 1, n
!     write(1,*) t(i),cnewton(i),ctrnewton(i)
  end do
  do i = 2,n-1
     if( ctrnewton(i) .and. .not.ctrnewton(i-1) ) tmin = t(i)
     if( ctrnewton(i) .and. .not.ctrnewton(i+1) ) tmax = t(i)
  end do
  write(1,*) "set style rect fc lt -1 fs solid 0.1 noborder"
  write(1,*) "set obj rect from ",tmin,", graph 0 to ",tmax,", graph 1"

  do i = 2,n-1
     if( cnewton(i) .and. .not.cnewton(i-1) ) xmin = t(i)
     if( cnewton(i) .and. .not.cnewton(i+1) ) then
        xmax = t(i)
        write(1,*) "set style rect fc lt -1 fs solid 0.2 noborder"
        write(1,*) "set obj rect from ",xmin,", graph 0 to ",xmax,", graph 1"
     end if
  end do
  close(1)

  xmin = -10
  xmax = 10
  call graph(funmean,s,xmin,xmax,'/tmp/rmean_region')

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


end program rmean_region
