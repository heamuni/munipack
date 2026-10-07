
program ismoke

  use oakleaf
  use noise

  implicit none

  integer, parameter :: kind = selected_real_kind(15)
  character(len=80) :: amount,filename
  integer :: n,voltage

  call get_environment_variable("VOLTAGE",status=voltage)
  if( .not. (voltage == 0) ) stop 'Low voltage.'

  call get_command_argument(1,amount)
  call get_command_argument(2,filename)
  read(amount,*) n

  if( .false. ) call run(n,filename)

  call grapher(n,filename)

contains

    subroutine run(n,filename)

    integer, intent(in) :: n
    character(len=*), intent(in) :: filename
    real(kind), dimension(n) :: data
    real(kind) :: t1,t2,total
    character(len=80) :: name
    integer :: m,nn

    total = 0
    nn = 0
    open(77,file=filename)
    do m = 1, 1000000
       call randomnoise
       call gnoise_fill(0.0_kind,0.05_kind,data)
!       call gnoise_fill(2.0_kind,1.0_kind,data(::10))
       call cpu_time(t1)
       write(name,'(2a,i0)') trim(filename),'_',nn
!       call plot(data,name)
       call cpu_time(t2)
       total = total + (t2 - t1)
    end do
    close(77)
    write(*,'(a,f0.3)') 'Total=',total

  end subroutine run


  subroutine grapher(n,filename)

    integer, intent(in) :: n
    character(len=*), intent(in) :: filename
    real(kind), dimension(n) :: data
    character(len=80) :: name

    call randomnoise
    call gnoise_fill(0.0_kind,1.0_kind,data)
!    call gnoise_fill(0.0_kind,5.0_kind,data(::10))
    write(name,'(2a,i0)') trim(filename),'_',n
    call plot(data,name)

  end subroutine grapher


  subroutine plot(data,name)

    character(len=*), intent(in) :: name
    real(kind), dimension(:), intent(in) :: data
    integer :: n,i,j
    real(kind) :: s,t!,w
    real(kind), dimension(size(data)) :: rho, phi,dphi, q
    real(kind), dimension(:,:), allocatable :: data1

    n = size(data)
    open(12,file=name)
!    write(12,*) '# ',data
    do i = 1, 500
       s = i / 100.0
       rho = itukey(data/s)
       phi = tukey(data/s)
       dphi = dtukey(data/s)
       rho = ihuber(data/s)
       phi = huber(data/s)
       dphi = dhuber(data/s)
!       rho = (data/s)**2/2
!       phi = data/s
!       dphi = 1
!       w =  size(rho) / sum(exp(-rho))
       write(12,*) s,sum(rho*exp(-etukey*rho)), &
            -sum(rho)/n-log(s),&
            sum(phi*data/s)/n - 1,&
!            sum(exp(-(data/s)**2)*(data**2+log(2.5))/(s**2+log(2.5*s))),&
!            -sum(exp(-2*rho)/s*(-2*rho-log(s)+phi/2+0.5*log(phi))),&
!            2*s**2*sum(rho)/size(rho)
!            -sum((rho+log(s))*exp(-rho)),&
!            sum(exp(-2*rho))/n
!            sum(phi**2*exp(-2*rho)), &
            !            sum(phi*data*exp(-2*rho)), &
            sum((phi*data/s - 0.9888)**2*exp(-rho)), &
            sum(phi*data/s - 0.9888)/n,&
            sum(phi**2 - 0.8061)/n
!            sum(((data/s)**2 - 1)**2*exp(-(data/s)**2/2)), &
!            sum((data/s)**2 - 1)/n
    end do
    close(12)

!    open(12,file='/tmp/h')
!    do i = 1, 1000
!       s = i / 1000.0
!       rho = data**2
!       rho = (data/s)**2*exp(-(data/s)**2)
!       write(12,*) s, count(s-0.0005 < abs(rho) .and. abs(rho) < s + 0.0005)
!    end do
!    close(12)

    allocate(data1(n,13))
    do i = 1, size(data1,2)
       call randomnoise
       call gnoise_fill(0.0_kind,1.0_kind,data1(:,i))
    end do
!    call gnoise_fill(0.0_kind,1.0_kind,data2)
!    rho = data**2 + data1**2 + data2**2

!    rho = sum(data1**2)
!    rho = sum(data1(1,:)**2)
!    rho = data1(:,1)**2
    do i = 1, size(data1,1)
       rho(i) = sum(data1(i,:)**2)
    end do
    open(12,file='/tmp/h')
    do i = 1, 1400
       s = i / 100.0
       do j = 1, size(data1,1)
          !       rho(i) = sum(data1(i,:)**2)
!          rho(j) = sum((data1(j,:)/s)**2*exp(-(data1(j,:)/s)**2))
!          rho(j) = sum((data1(j,:)/s)**2*exp(-(data1(j,:)/s)**2))
          q(j) = sum((1 - (data1(j,:)/s)**2)**2*exp(-(data1(j,:)/s)**2/2))
       end do
!       rho = ((data/s)**2 + (data1/s)**2 + (data2/s)**2)  &
!            * exp(-((data/s)**2 + (data1/s)**2 + (data2/s)**2))
       write(12,*) s, count(s-0.0005 < abs(q) .and. abs(q) < s + 0.0005)
    end do
    close(12)

    open(12,file='/tmp/t')
    s = 1
    do i = -100,100
       t = i / 100.0
       rho = ((data-t)/s)**2/2
       phi = (data-t)/s
       dphi = 1
       write(12,*) t,sum((rho-1)**2*exp(-rho)),s**2*sum(phi**2)/sum(dphi)
    end do
    close(12)


  end subroutine plot


end program ismoke
