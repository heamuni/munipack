
program sr1_anim

  use rfun
  use rlike
  use noise
  use iso_fortran_env
  implicit none

  type(robust_likelihood_hampel), target :: hampelfun
  type(robust_likelihood_tukey), target :: tukeyfun
  type(robust_likelihood_huber), target :: huberfun
  type(nonrobust_likelihood_square), target :: squarefun
  class(robust_likelihood), pointer :: funmean
  integer, parameter :: datamax = 13
  real(REAL64), dimension(:), allocatable, target :: data

  logical, parameter :: switch_on_generator = .true.
  character(len=666), parameter :: datafile = '/tmp/tsr1mean_data'
  character(len=666) :: dir
  integer :: i, n, step, ndat, iflag
  real(REAL64) :: frac, azimuth, zenit
  real(REAL64) :: tmin, tmax, smin, smax, ts

  call get_command_argument(1,dir)
  if( trim(dir) == '' ) stop 'empty dir'

  if( switch_on_generator ) then
     ! random generated data
     allocate(data(datamax))
     call gnoise_fill(0.0_REAL64,0.25_REAL64,data)
     frac = 0.25
     step = datamax / max(int(frac*datamax),1)
!     step = 3
     write(*,*) 'step:',frac*datamax,step
     do i = 2, datamax, step
        data(i) = gnoise(-10.0_REAL64,1.0_REAL64)
        write(*,*) 'bad point:',real(data(i))
!        if( abs(data(i)) < 1 ) ts = data(i)
!        if( i == 2 ) ts = data(i)
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

  funmean => hampelfun
  call funmean%init(data,iflag)

!  call old

  call new


contains

  subroutine new

    integer :: j

    open(2,file='/tmp/x')

    n = 0
    zenit = 60
    azimuth = 75

    tmin = -20
    tmax = 20
    smin = 0
    smax = 20
    do i = 0, 120
       azimuth = mod(70+i,360)
       zenit = 60 + 15*i/120.0
       write(2,*) tmin,tmax,smin,smax,azimuth,zenit
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

!!$    azimuth = 200
!!$    do i = 170,0,-1
!!$       tmin = -3 - 0.1*i
!!$       tmax =  3 + 0.1*i
!!$       smax =  3 + 0.1*i
!!$       zenit = 60 - 30*(i-170.0)/170.0
!!$       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
!!$       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
!!$       n = n + 1
!!$    end do
!!$
!!$    do i = 200,270
!!$       azimuth = mod(i,360)
!!$       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
!!$       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
!!$       n = n + 1
!!$    end do

    ts = data(1)
    do i = 1, size(data)
       if( -1 < data(i) .and. data(i) < 1 ) then
          ts = data(i)
       end if
    end do
    ts = data(9)

    !do i = 300, 30, -1
    !    do i = 90, 300
    do i = 300,0,-1
       j = 300 - i
!       tmin = -3 + (ts + 3 - 0.5)*(300-j)/270.0
!       tmax =  3 - (3 - ts - 0.5)*(300-j)/270.0
!       tmin = ts + (-20+-ts)/(20-0.3) *j/300
       !       tmin = ts +  (20-0.3-ts)/(20-0.3) *j/
       tmax =  0.3 + ts + i*(20 - (0.3 + ts))/300.0
       tmin = -0.3 + ts - i*(20 - (0.3 - ts))/300.0
       smax = 0.3 + i*(20-0.3)/300
       azimuth = 190 + (360-190)*j /300.0
!       zenit = 30 + (59.9)*(i - 30) /(300.0-30.0)
       zenit = 75 - 75*j/300.0
       write(2,*) tmin,tmax,smin,smax,azimuth,zenit
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    do i = 300, 1, -1
       tmin = ts - 0.001*i
       tmax = ts + 0.001*i
       smax =  0.001*i
       write(2,*) tmin,tmax,smin,smax,azimuth,zenit
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    close(2)

  end subroutine new


  subroutine old

    n = 0
    zenit = 60
    azimuth = 75

    tmin = -20
    tmax = 20
    smin = 0
    smax = 20
    do i = 75,360+75
       azimuth = mod(i,360)
       !     tmin = -5 - 0.1*i
       !     tmax =  5 + 0.1*i
       !     smax =  5 + 0.1*i
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    azimuth = 75
    do i = 150,0,-1
       tmin = -5 - 0.1*i
       tmax =  5 + 0.1*i
       smax =  5 + 0.1*i
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    do i = 75,360+75
       azimuth = mod(i,360)
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do
    azimuth = 75

    do i = 60,105
       zenit = i
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    tmin = -5
    tmax = 5
    smin = 0
    smax = 5
    do i = 75, 360+75
       azimuth = mod(i,360)
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       !     write(filename,'(2a,i5.5,a)') trim(dir),'/a',n,'.png'
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    azimuth = 75
    do i = 105,60,-1
       zenit = i
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    do i = 75, 190
       azimuth = i
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       !     write(filename,'(2a,i5.5,a)') trim(dir),'/a',n,'.png'
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    do i = 500, 100, -1
       tmin = -0.01*i
       tmax =  0.01*i
       smax =  0.01*i
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       !     write(filename,'(2a,i5.5,a)') trim(dir),'/a',n,'.png'
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    do i = 100, 500, 10
       tmin = -0.01*i
       tmax =  0.01*i
       smax =  0.01*i
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       !     write(filename,'(2a,i5.5,a)') trim(dir),'/a',n,'.png'
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    do i = 190, 350
       azimuth = i
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       !     write(filename,'(2a,i5.5,a)') trim(dir),'/a',n,'.png'
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

    !  ts = data(2)
    ts = data(1)
    do i = 1, size(data)
       if( -1 < data(i) .and. data(i) < 1 ) then
          ts = data(i)
       end if
    end do

    do i = 500, 50, -1
       tmin = -5 + (ts + 5 - 0.5)*(500-i)/450.0
       tmax =  5 - (5 - ts - 0.5)*(500-i)/450.0
       !     tmin = ts - ((5-ts)/450.0)*(500-i)
       !     tmax = ts + (ts/450.0)*(500-i)
       !     smax = 5 - (ts/450.0) *(500-i)
       smax = 0.01*i
       zenit = 60 - (59.0/450.0)*(500-i)
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       !     write(filename,'(2a,i5.5,a)') trim(dir),'/a',n,'.png'
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do


    !  zenit = 0
    do i = 500, 5, -1
       tmin = ts - 0.001*i
       tmax = ts + 0.001*i
       smax =  0.001*i
       call graph(funmean,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       call gplot(dir,n,tmin,tmax,smin,smax,azimuth,zenit)
       n = n + 1
    end do

  end subroutine old

  subroutine gplot(dirname,n,tmin,tmax,smin,smax,azimuth,zenit)

    character(len=*), intent(in) :: dirname
    integer, intent(in) :: n
    real(REAL64), intent(in) :: tmin,tmax,smin,smax,azimuth,zenit
    character(len=666) :: com,errmsg,filename
    integer :: estat

    errmsg = ''

    write(filename,'(2a,i5.5,a)') trim(dirname),'/a',n,'.png'

    write(com,'(2a,1x,6(g0,1x))') &
         'sh -x gnuplot_anim.sh ',trim(filename),tmin,tmax,smin,smax,azimuth,zenit
    call execute_command_line(com,exitstat=estat,cmdmsg=errmsg)
    if( estat /= 0 ) then
       write(*,*) 'Error for ',trim(com),": ",trim(errmsg)
       stop
    end if

  end subroutine gplot

  subroutine graph(fcn,x0,x1,y0,y1,filename)

    class(robust_likelihood), intent(in) :: fcn
    real(REAL64), intent(in) :: x0,x1,y0,y1
    character(len=*), intent(in) :: filename
    integer, parameter :: nx = 200
    integer, parameter :: ny = 200
    real(REAL64) :: x,dx,y,dy
    integer :: i,j

    dx = (x1 - x0) / (nx - 1)
    dy = (y1 - y0) / (ny - 1)

    open(1,file=filename)
    do j = 1, ny
       y = y0 + dy*j
       do i = 0, nx
          x = x0 + dx*i
          write(1,*) x,y,fcn%func(x,y)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph


end program sr1_anim
