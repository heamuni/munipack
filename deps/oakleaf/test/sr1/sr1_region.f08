
program sr1_region

  use rfun
  use sr1mean_function
  use trustsr1mean
  use noise
  use iso_fortran_env
  implicit none

  real(REAL64), parameter :: tmin = -5
  real(REAL64), parameter :: tmax =  5
  real(REAL64), parameter :: smax =  5
  real(REAL64), parameter :: dt = 0.05_REAL64
  real(REAL64), parameter :: ds = 0.05_REAL64
  integer, parameter :: nmax = nint((tmax - tmin) / dt)
  integer, parameter :: mmax = nint((smax - ds) / ds)
  logical, dimension(0:nmax,mmax) :: reg

  type(tsr1mean_tukey), target :: tukeymean
  type(tsr1mean_huber), target :: hubermean
  type(tsr1mean_experimental), target :: bathmean
  class(tsr1mean_function), pointer :: funmean
  integer, parameter :: datamax = 7
  real(REAL64), dimension(:), allocatable, target :: data

  character(len=666) :: robfun = 'tukey'
  logical, parameter :: switch_on_generator = .false.
  character(len=666), parameter :: datafile = '/tmp/sr1region_data'
  real(REAL64) :: t,s,frac
  integer :: i,j,step,ndat

  call get_command_argument(1,robfun)

  if( switch_on_generator ) then
     ! random generated data
     allocate(data(datamax))
     call gnoise_fill(0.0_REAL64,1.0_REAL64,data)
     frac = 0.15
     step = datamax / max(int(frac*datamax),1)
     step = 99
     write(*,*) 'step:',frac*datamax,step
     do i = 2, datamax, step
        data(i) = gnoise(-7.0_REAL64,1.0_REAL64)
!        write(*,*) 'bad point:',real(data(i))
    end do
!    write(*,*) 'data:',real(data)
  else
     ! saved data
     open(3,file=datafile,form='unformatted')
     read(3) ndat
     allocate(data(ndat))
     read(3) data
     close(3)
!     write(*,*) 'data:',real(data)
  end if

  if( robfun == 'huber' ) then
     funmean => hubermean
  else if( robfun == 'hampel' ) then
     funmean => bathmean
  else
     funmean => tukeymean
  end if
  funmean%data => data

  do i = 0, nmax
     t = tmin + i*dt
     do j = 1, mmax
        s = j*ds

        block
          logical :: reli
          logical, parameter :: verbose = .true.
          logical, parameter :: quasi = .true.
          integer, parameter :: maxfev = 3*precision(t)
          real(REAL64), dimension(2,2) :: hess
          real(REAL64) :: tol,delta,stderr
          integer :: info
          character(len=666) :: msg

          tol = 0.01 * s / sqrt(real(size(data)))
          delta = 0.1*s

          call tsr1meansqp(funmean,t,s,hess,delta,tol,maxfev,quasi,verbose,info,msg)
          stderr = funmean%stderr(t,s,reli)

          reg(i,j) = info == 1 .and. abs(t - 0) < 0.5 .and. abs(s - 1) < 0.5
          reg(i,j) = info == 1 .and. reli
        end block

     end do
  end do


  call graph(funmean,tmin,tmax,smax,'/tmp/sr1region.dat')
  call graph_reg(reg,tmin,dt,ds,'/tmp/sr1region_region_reg.dat')
  call graph_reg2(reg,tmin,dt,ds,'/tmp/sr1region_region.dat')

  if( switch_on_generator ) then
     ! save data
     open(3,file=datafile,form='unformatted')
     write(3) datamax
     write(3) data
     close(3)
  end if

contains

    subroutine graph(fcn,x0,x1,y1,filename)

    class(tsr1mean_function), intent(in) :: fcn
    real(REAL64), intent(in) :: x0,x1,y1
    character(len=*), intent(in) :: filename
    integer, parameter :: nx = 200
    integer, parameter :: ny = 200
    real(REAL64) :: x,dx,y,dy
    integer :: i,j

    dx = (x1 - x0) / (nx - 1)
    dy = (y1 -  0) / (ny - 1)

    open(1,file=filename)
    do j = 1, ny
       y = 0 + dy*j
       do i = 0, nx
          x = x0 + dx*i
          write(1,*) x,y,fcn%val(x,y)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph


  subroutine graph_reg(reg,tmin,dt,ds,filename)

    logical, dimension(0:,:), intent(in) :: reg
    real(REAL64), intent(in) :: tmin,dt,ds
    character(len=*), intent(in) :: filename
    real(REAL64) :: t,s
    integer :: i,j,q

    open(1,file=filename)
    do i = 0, ubound(reg,1)
       t = tmin + dt*i
       do j = 1, ubound(reg,2)
          s = ds*j
          if( reg(i,j) ) then
             q = 1
          else
             q = 0
          end if
          write(1,*) t,s,q!,reg(i,j)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_reg

  subroutine graph_reg2(reg,tmin,dt,ds,filename)

    logical, dimension(0:,:), intent(in) :: reg
    real(REAL64), intent(in) :: tmin,dt,ds
    character(len=*), intent(in) :: filename
    real(REAL64) :: t,smin,smax
    integer :: i,j,q

    open(1,file=filename)
    do i = 0, ubound(reg,1)
       t = tmin + dt*i
       smin = -1
       smax = -1!ds*ubound(reg,2)
       do j = 1, ubound(reg,2)
          if( reg(i,j) ) then
             smin = ds*j
             exit
          end if
       end do
       do j = ubound(reg,2),1,-1
          if( reg(i,j) ) then
             smax = ds*j
             exit
          end if
       end do
       write(1,*) t,smin,smax
    end do
    close(1)

  end subroutine graph_reg2

end program sr1_region
