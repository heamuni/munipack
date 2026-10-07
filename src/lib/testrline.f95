

! gfortran -fcheck=all -g -p  -Wall testrline.f95 -L. -L../minpack -lnoise -lrstat -lsort -llmin -lminpacks -lminpack -lfmm -lm

program testrline

  use robustline
  use noise

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)
  integer, parameter :: nmax = 100
  real(dbl), dimension(:),allocatable :: x,y,dx,dy
  real(dbl) :: a,b,da,db,sig,t
  integer :: i,n
  logical :: reliable

  call rektoris
!  stop

  allocate(x(nmax),y(nmax),dx(nmax),dy(nmax))
  do i = 1,size(x)

     call random_number(t)
     x(i) = t
     if( mod(i,10) > 0 ) then
        y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,0.01_dbl)
     else
        y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,0.1_dbl)
     end if

  end do
  i = size(x)/3
!  y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,0.1_dbl)
  i = size(x)/2
!  y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,0.1_dbl)
  i = size(x)/4
!  y(i) = 1 + 1*x(i) + gnoise(0.0_dbl,0.1_dbl)

  dx = 0.01 !* 0.707
  dy = 0.01 !* 0.707
  dx = 0
  dy = 1

!  deallocate(x,y,dx,dy)
!  read(*,*) n
!  allocate(x(n),y(n),dx(n),dy(n))
!  do n = 1,size(x)
!     read(*,*) x(n),y(n),dx(n),dy(n)
!  enddo

  call rline(x,y,a,b,da,db,dx,dy,sig,reliable,verbose=.true.)

  write(*,*) 'a=',a,'+-',da
  write(*,*) 'b=',b,'+-',db
  write(*,*) 'sig=',sig,' reliable:',reliable

  open(1,file='/tmp/rline2')
  do i = 1, size(x)
     write(1,*) x(i),y(i),y(i) - (a + b*x(i))
  end do
  close(1)

  deallocate(x,y,dx,dy)


contains

  subroutine rektoris

    integer, parameter :: ndat = 14
    real(dbl), dimension(ndat) :: x,y
    real(dbl) :: a,b,da,db,sig
    logical :: reliable

    x = (/3,4,6,6,6,7,8,9,11,11,12,12,14,16 /)
    y = (/24.82, 23.26, 14.77, 19.06, 14.79, 17.66, 11.83, &
         14.82, 5.16, 11.12, 8.04, 2.72, 0.74, 1.21/)

    call rline(x,y,a,b,da,db,sigma=sig,reliable=reliable,verbose=.true.)

    write(*,*) "rline:", reliable
    write(*,*) "sig:",sig
    write(*,*) "a:",a,"+-",da
    write(*,*) "b:",b,"+-",db

  end subroutine rektoris

end program testrline
