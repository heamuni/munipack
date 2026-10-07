
program rfun_multimin

  use rfun
  use sr1mean_function
  use noise
  use iso_fortran_env
  implicit none

  type(tsr1mean_tukey), pointer :: tukeymean
  type(tsr1mean_experimental), pointer :: bathmean
  class(tsr1mean_function), pointer :: funmean

  integer, parameter :: nsamples = 999
  integer, parameter :: datamax = 5
  logical, parameter :: switch_on_generator = .true.
  character(len=666), parameter :: datafile = '/tmp/rmean_datafail'
  real(REAL64), dimension(:), allocatable :: smin, topt, f
  real(REAL64), dimension(:), pointer :: data
  integer, dimension(:), allocatable :: hist
  real(REAL64) :: tmin, tmax, dt, s, t, tc, w
  integer :: i, j, l, n, m, iflag

  tmin = real(-20,REAL64)
  tmax = real( 20,REAL64)
  dt = 1d-2
  s = 1.0_REAL64
  w = 0.1_REAL64

!  write(*,*) tmin,tmax,dt,s
!  stop

  n = int((tmax - tmin) / dt)
  allocate(smin(n),topt(n),hist(n),f(n))
  hist = 0

  do i = 1, nsamples
!     write(*,*) i
     call one_sample(data)

     allocate(bathmean,tukeymean)
     call bathmean%init(data,iflag)
     call tukeymean%init(data,iflag)

     funmean => tukeymean
     funmean => bathmean

     do j = 1, n
        t = tmin + j*dt
        f(j) = funmean%val(t,s)
     end do

!     do j = 100, 100!n
     do j = 1, n
        tc = tmin + j*dt
        call one_round(funmean,w,tc,tmin,tmax,dt,s,f)
     end do

!     if( m /= 1 ) then
     if( m > 1 ) then
        do j = 1, m
           write(*,*) j,topt(j),smin(j)
        end do
        open(1,file='/tmp/rfun_multifun')
        do j = 1, n
           t = tmin + j*dt
           write(1,*) t,funmean%val(t,s),(ipshampel((t-data(l))/s),l=1,datamax),&
                ipshampel((t-tc)/s), &
                ((f(j)-log(s)) + w*ipshampel((t-tc)/s))/(real(1,REAL64) + w)+log(s)
!           write(1,*) t,funmean%val(t,s),(itukey((t-data(l))/s),l=1,datamax)
!           write(1,*) t,funmean%val(t,s)
        end do
        close(1)
        stop
     end if

     funmean => null()
     deallocate(data,bathmean,tukeymean)
  end do

  open(1,file='/tmp/rfun_multimin_10_hampel7_25')
  do i = 1, n
     t = tmin + i*dt
     write(1,*) t,hist(i)
  end do
  close(1)

contains

  subroutine one_sample(data)

    real(REAL64), dimension(:), pointer, intent(out) :: data
    real(REAL64) :: frac
    integer :: i, step, ndat

    ! random generated data
    if( switch_on_generator ) then
       allocate(data(datamax))
       call gnoise_fill(0.0_REAL64,1.0_REAL64,data)
       frac = 0.1
       step = datamax / max(int(frac*datamax),1)
!       step = 6
!       write(*,*) 'step:',frac*datamax,step
       do i = 2, datamax, step
          data(i) = gnoise(-12.0_REAL64,1.3_REAL64)
!          write(*,*) 'bad point:',real(data(i))
       end do
!       data(datamax/2) = gnoise(-5.0_REAL64,1.3_REAL64)
!       write(*,*) 'data:',real(data)
    else
       ! saved data
       open(3,file=datafile,form='unformatted')
       read(3) ndat
       allocate(data(ndat))
       read(3) data
       close(3)
       write(*,*) 'data:',real(data)
    end if

  end subroutine one_sample

  subroutine one_round(funmean,w,tc,tmin,tmax,dt,s,ff)

    class(tsr1mean_function), pointer, intent(in) :: funmean
    real(REAL64), intent(in) :: w,tc,tmin,tmax,dt,s
    real(REAL64), dimension(:), intent(in) :: ff
    real(REAL64), dimension(:), allocatable :: f,gt
    real(REAL64), dimension(2) :: g
    real(REAL64) :: tol
    integer :: i

    tol = 10*n*epsilon(f)
!    tol = 1e-7
    allocate(f(n),gt(n))

    do i = 1, n
       t = tmin + i*dt
!       f(i) = (ff(i) + w*itukey((t-tc)/s))/(real(1,REAL64) + w)
       f(i) = ((ff(i)-log(s)) + w*ipshampel((t-tc)/s))/(real(1,REAL64) + w) + log(s)
       !       f(i) = (funmean%val(t,s) + w*ipshampel(t-tc))/(real(1,REAL64) + w)
       call funmean%grad(t,s,g)
       gt(i) = g(2)
    end do

    m = 0
    do i = 2, n-1
       if( f(i-1) > f(i) .and. f(i) < f(i+1) .and. &
            abs(f(i-1) - f(i)) > tol .and. abs(f(i) - f(i+1)) > tol  ) then
          m = m + 1
          smin(m) = f(i)
          topt(m) = tmin + i*dt
          hist(i) = hist(i) + 1
!          write(*,*) '>',m,tmin + i*dt,f(i-1:i+1),f(i-1:i+1)-f(i),gt(i-1:i+1)
       end if
    end do

  end subroutine one_round

end program rfun_multimin
