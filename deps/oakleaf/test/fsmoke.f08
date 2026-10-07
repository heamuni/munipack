
program fsmoke

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

  call run(n,filename)

contains

    subroutine run(n,filename)

    integer, intent(in) :: n
    character(len=*), intent(in) :: filename
    real(kind), dimension(n) :: data1, errors1, data2, errors2
    real(kind) :: t,s,ds,c,t1,t2,total
    character(len=80) :: name
    integer :: m,nn,flag

    total = 0
    nn = 0
    open(77,file=filename)
    do m = 1, 1000000
       call randomnoise
       call pnoise_fill(10000.0_kind,data1)
       errors1 = sqrt(data1)
       where( errors1 < epsilon(errors1) )
          errors1 = 100
       end where
       call pnoise_fill(2500.0_kind,data2)
       errors2 = sqrt(data2)
       where( errors2 < epsilon(errors2) )
          errors2 = 50
       end where
       call cpu_time(t1)
       call fmean(data1,errors1,data2,errors2,t,ds,s,c,flag=flag)
       call cpu_time(t2)
       total = total + (t2 - t1)
       write(77,*) t,ds,s,c,flag
       if( flag /= 0 ) then
          open(22,file='/tmp/fail',position='append')
          write(22,*) data1
          close(22)
          nn = nn + 1
          write(name,'(2a,i0)') trim(filename),'_',nn
          call plot(data1,errors1,name,t)
          call surface(data1,errors1,trim(name)//'_s')
!          t = sum(data/errors**2)/sum(1/errors**2)
!          s = sqrt(sum((t-data)**2/errors**2)/(n-1))
!          c = sqrt(sum((t-data)**2)/(n-1))
!          ds = sqrt(s**2/sum(1/errors**2))
!          write(*,*) '# ',nn,' Arith:',t,s,ds,c
          call fmean(data1,errors1,data2,errors2,t,ds,s,c,verbose=.true.)
       end if
    end do
    close(77)
    write(*,'(a,f0.3)') 'Total=',total

  end subroutine run

  subroutine plot(data,errors,name,t)

    character(len=*), intent(in) :: name
    real(kind), dimension(:), intent(in) :: data, errors
    real(kind), intent(in) :: t
    integer :: n,i,j
    real(kind) :: x,s,e,w,q,e0,r,scale

    scale = 1

    n = size(data)
    open(12,file=name)
    write(12,*) '# ',data
    do i = -300, 300
       x = i / 100.0
       s = 0
       w = 0
       do j = 1, n
          s = s + tukey((x - data(j))/(errors(j)*scale))
          w = w + dtukey((x - data(j))/(errors(j)*scale))
       end do
       e = 0
       e0 = 0
       q = 0
       if( i > 0 ) then
          do j = 1, n
             r = itukey((data(j)-t)/(errors(j)*x))
             e = e + r*exp(-2.1227*r)
             r = itukey(data(j)/(errors(j)*x))
             e0 = e0 + r*exp(-2.1227*r)
             q = q + dtukey(data(j)/(errors(j)*scale))
          end do
       end if
       write(12,*) x, s/n,e/n,e0/n,w/n,q/n
    end do
    close(12)

  end subroutine plot


  subroutine surface(data,errors,name)

    real(kind), dimension(:), intent(in) :: data,errors
    character(len=*), intent(in) :: name

    real(kind), dimension(size(data)) :: r,x
    real(kind) :: t,s,d
    integer :: i,j,n

    n = size(data)

    open(48,file=name)
    do i = -30,30
       t = i / 10.0
       do j = 1,30
          s = j / 10.0
          x = (data-t)/(errors*s)
          r = itukey(x)
          d = sum(dtukey(x))
          if( abs(d) > 0*n/2.0 ) then
             write(48,*) t,s,sum(itukey(x))/n,sum(r*exp(-2.1227*r))/n, &
                  sum(tukey(x)**2)/sum(dtukey(x))**2
          end if
       end do
    end do
    close(48)

  end subroutine surface

end program fsmoke
