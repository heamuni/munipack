
!
! gfortran -Wall testrfun.f95 -L. -lrstat -lm
!

program testrfun

  use rfun

  implicit none

  integer, parameter :: rp = selected_real_kind(15)

  integer :: i
  real(rp) :: x

  do i = -100,100
     x = i/10.0
!     write(*,*) x,huber(x),dhuber(x),hampel(x),dhampel(x),andrews(x),dandrews(x), &
!          tukey(x),dtukey(x)
     write(*,*) x,dhuber(x),(huber(x+0.1)-huber(x-0.1))/0.2, &
          dhampel(x),(hampel(x+0.1)-hampel(x-0.1))/0.2,&
          dandrews(x),(andrews(x+0.1)-andrews(x-0.1))/0.2,&
          dtukey(x),(tukey(x+0.1)-tukey(x-0.1))/0.2
  end do

end program testrfun
