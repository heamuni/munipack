
program smoke

  use oakleaf
  use noise

  implicit none

  integer, parameter :: kind = selected_real_kind(15)
!  integer, dimension(1) :: n
  character(len=80) :: amount,filename
  integer :: n,voltage

  call get_environment_variable("VOLTAGE",status=voltage)
  if( .not. (voltage == 0) ) stop 'Low voltage.'

  call get_command_argument(1,amount)
  call get_command_argument(2,filename)
  read(amount,*) n

  call run(n,filename)

!  n = [ 2,3,4,5,6,7,8,9,10,20,50,10**2,10**3,10**4,10**6 ]
!  n = [ 2,3,4,5,6,7,8,9,10 ]
!  n = [ 7 ]
!  do l = 1, size(n)
!     call run(n(l))
!  end do

contains

  subroutine run(n,filename)

    integer, intent(in) :: n
    character(len=*), intent(in) :: filename
    real(kind), dimension(n) :: data
    real(kind) :: t,s,ds,c,t1,t2,total
    logical :: reliable
    character(len=80) :: name
    integer :: m,nn,flag

!    data = [5.2129640725329870E-002,  0.15293191139569814, &
!         -0.20754828909183076, -0.39700992178364081, 1.3411645446578198]
!    call rmean(data,t,ds,s,verbose=.true.)
!    stop

!    data = [ -0.32366372162119911,      -0.34013885339444594,  &
!         -0.28319988121743622 , -0.15217523163307081, &
!         -0.67745093584468041,      -0.29660261571936658, 1.2883265222850431 ]

!    data = [ -1.3666496394616230, 3.0716410066064586, 4.6510401721685517E-002, !&
!         1.2641804668388122E-002, 6.0265855187624381E-002, &
!         -0.77198678639721452, -0.38235178598268538 ]
!!$    data = [ -1.3666496394616230, 3.0716410066064586, 4.6510401721685517E-002, &
!!$         1.2641804668388122E-002, 6.0265855187624381E-002, &
!!$         -0.77198678639721452, -0.38235178598268538 ]
!    call randomnoise
!    call gnoise_fill(0.0_kind,1.0_kind,data)
!!$    call surface(data)
!    t = sum(data) / size(data)
!    s = sqrt(sum((data - t)**2)/(size(data) - 1))
!    write(*,*) t,s
!    call rmean(data,t,ds,s,verbose=.true.)
!!$    call plot(data,'/tmp/lll',-0.369961183d0)
!   stop
    !    write(name,'(a,i0)') '/tmp/smoke',n
    total = 0
    nn = 0
    open(77,file=filename)
    do m = 1, 1000000
       call randomnoise
       call gnoise_fill(0.0_kind,1.0_kind,data)
       call gnoise_fill(2.0_kind,10.0_kind,data(::11))
       call gnoise_fill(0.0_kind,0.01_kind,data(::13))
       call cpu_time(t1)
!       call rmean(data,t,ds,s,c,reliable=reliable,flag=flag)
       call rmean(data,t,ds,s,c,flag=flag)
       call cpu_time(t2)
       total = total + (t2 - t1)
       write(77,*) t,ds,s,c,flag
!       if( .not. reliable .or. (n == 4 .and. s < 0.07) ) then
!       if( .not. reliable ) then
       if( flag /= 0 ) then
          open(22,file='/tmp/fail',position='append')
          write(22,*) data
          close(22)
          nn = nn + 1
          write(name,'(2a,i0)') trim(filename),'_',nn
          call plot(data,name,t)
          call surface(data,trim(name)//'_s')
          t = sum(data) / size(data)
          s = sqrt(sum((data - t)**2)/(size(data) - 1))
          write(*,*) '# ',nn,' Arith:',t,s
          call rmean(data,t,ds,s,verbose=.true.)
!          stop
       end if
    end do
    close(77)
    write(*,'(a,f0.3)') 'Total=',total

  end subroutine run

  subroutine plot(data,name,t)

    character(len=*), intent(in) :: name
    real(kind), dimension(:), intent(in) :: data
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
          s = s + tukey((x - data(j))/scale)
          w = w + dtukey((x - data(j))/scale)
       end do
       e = 0
       e0 = 0
       q = 0
       if( i > 0 ) then
          do j = 1, n
             r = itukey((data(j)-t)/x)
             e = e + r*exp(-2.1227*r)
             r = itukey(data(j)/x)
             e0 = e0 + r*exp(-2.1227*r)
             q = q + dtukey(data(j)/scale)
          end do
       end if
       write(12,*) x, s/n,e/n,e0/n,w/n,q/n
    end do
    close(12)

  end subroutine plot


  subroutine surface(data,name)

    real(kind), dimension(:), intent(in) :: data
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
          x = (data-t)/s
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

end program smoke
