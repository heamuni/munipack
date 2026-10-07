
program trrho

  ! tr rho parameter study

  use rfun
  use sr1mean_function
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
  real(REAL64) :: tmin, tmax, dt, frac, s, p
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

  allocate(bathmean,tukeymean)
  call bathmean%init(data,iflag)
  call tukeymean%init(data,iflag)

  tmin =-10
  tmax = 10
  dt = 0.01
  s = 1
  p = 0.33

  call graph(bathmean,s,p,tmin,tmax,dt,'/tmp/trrho_pshampel')
  call graph(tukeymean,s,p,tmin,tmax,dt,'/tmp/trrho_tukey')
  call graphs(bathmean,s,tmin,tmax,dt,'/tmp/trrho_pshampels')

  deallocate(bathmean,tukeymean)

contains

  subroutine graph(fcn,s,p,tmin,tmax,dt,filename)

    class(tsr1mean_function), intent(in) :: fcn
    real(REAL64), intent(in) :: s,p,tmin,tmax,dt
    character(len=*), intent(in) :: filename
    real(REAL64) :: t,g1,h11,drho,dq,rho
    integer :: i

    open(1,file=filename)
    do i = 1, int((tmax - tmin) / dt)
       t = tmin + i*dt

       drho = fcn%val(t+p,s) - fcn%val(t,s)
       call fcn%rnewton(t,s,g1,h11)
       dq = g1*p + 0.5*h11*p**2
!       if( abs(dq) > 0.001*s .and. abs(drho) > 0.001*s ) then
       if( abs(dq) > 0.001*s ) then
          rho = drho/dq
          write(1,*) t,rho,drho,dq,g1,h11,fcn%val(t,s)
       end if
    end do
    close(1)

  end subroutine graph

  subroutine graphs(fcn,s,tmin,tmax,dt,filename)

    class(tsr1mean_function), intent(in) :: fcn
    real(REAL64), intent(in) :: s,tmin,tmax,dt
    character(len=*), intent(in) :: filename
    integer, parameter :: np = 5
    real(REAL64), dimension(np), parameter :: p = [0.1, 0.33, 0.5, 1.0, 2.0]
    real(REAL64), dimension(np) :: rho
    real(REAL64) :: t,g1,h11,drho,dq
    integer :: i,j

    open(1,file=filename)
    do i = 1, int((tmax - tmin) / dt)
       t = tmin + i*dt

       do j = 1, size(p)
          drho = fcn%val(t+p(j),s) - fcn%val(t,s)
          call fcn%rnewton(t,s,g1,h11)
          dq = g1*p(j) + 0.5*h11*p(j)**2
          if( abs(dq) > 0.001*s ) then
             rho(j) = drho / dq
          end if
       end do
       write(1,*) t,rho
    end do
    close(1)

  end subroutine graphs

end program trrho
