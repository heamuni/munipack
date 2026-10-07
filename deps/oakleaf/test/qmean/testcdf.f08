
program cdf

  use quantilmean
  use medians
  use noise
  use toolbox

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)
  integer, parameter :: ndata = 9
  real(dbl), dimension(ndata) :: x,y,z,xcdf,pcdf
  integer :: i

  call gnoise_fill(0.0_dbl,1.0_dbl,x)
  call gnoise_fill(-15.0_dbl,1.0_dbl,y)
  call gnoise_fill(1.0_dbl,0.01_dbl,z)

  do i = 1, 0*ndata,10
     x(i) = y(i)
  end do

  do i = 5, 0*ndata,5
     x(i) = z(i)
  end do


  call ecdf(x,xcdf,pcdf)

  open(1,file='/tmp/testecdf')
  do i = 1, size(xcdf)
     write(1,*) xcdf(i),pcdf(i)
  end do
  close(1)

  write(*,*) xcdf((ndata+1)/2),qmedian(x)


end program cdf
