
program plot_rfun

  use oakleaf
  use iso_fortran_env

  implicit none

  real(REAL64) :: x
  integer :: i

  open(1,file='hubers')
  do i = -1000,1000
     x = i / 100.0
     write(1,*) x,huber(x),dhuber(x),ihuber(x)
  end do
  close(1)

  open(1,file='tukeys')
  do i = -1000,1000
     x = i / 100.0
     write(1,*) x,tukey(x),dtukey(x),itukey(x)
  end do
  close(1)

end program plot_rfun
