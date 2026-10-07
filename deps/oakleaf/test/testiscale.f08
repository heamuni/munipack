
program testlscale

  use oakleaf
  use noise
  use toolbox

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)
  integer, parameter :: ndata = 900000
  real(dbl), dimension(ndata) :: x,y,z
  real(dbl) :: s, t, stderr
  integer :: i
  logical :: reliable

  call gnoise_fill(0.0_dbl,10.0_dbl,x)
!  call gnoise_fill(0.0_dbl,100.1_dbl,y)
!  call gnoise_fill(0.0_dbl,0.01_dbl,z)

  do i = 1, ndata,5
     x(i) = y(i)
  end do

  do i = 5, ndata,5
     x(i) = z(i)
  end do


!  write(*,*) x

!  x = 1
!  x = [-0.936962512663888, 0.982829448493290, -0.045866935829402 ]
!  x = [ 177.0, 0.431, -0.104, 1.451, -1.165, 0.516, 0.0 ]

  t = sum(x) / ndata
  write(*,*) 'arith.mean:',t
  write(*,*) 'median:',median(x)
  t = median(x)
  t = 0

  write(*,*) 'std.dev. ',stddev(x-t)
  s = median(abs(x-t))/0.6745
  write(*,*) 'MAD/0.6745 ',s
!    s = 1
  !s = 2.18
  call iscale(x-t,s,reliable,verbose=.true.)
  stderr = s*sqrt(sum(tukey((x-t)/s)**2)/sum(dtukey((x-t)/s))**2)
  write(*,*) 'scale ',s,reliable, stderr, stderr*sqrt(real(size(x)))


  open(1,file='/tmp/x')
  write(1,*) '# x'
  write(1,'(f20.15)') x-t
  close(1)

  do i = 1, ndata
!     call plots(i,x(i)-t)
  end do

contains

  subroutine plots(i,r)

    use oakleaf

    integer, intent(in) :: i
    real(dbl), intent(in) :: r
    character(len=80) :: filename
    real(dbl) :: s,smin,ds,rho
    integer :: l

    write(filename,'(a,i0)') '/tmp/x',i
    open(1,file=filename)
    smin = tiny(s)
    ds = (5.0 - 0.0) / 66
    do l = 1, 66
       s = smin + l * ds
       rho = itukey(r/s)
       write(1,'(2f20.15)') s,- rho - log(s)
    end do
    close(1)

  end subroutine plots

end program testlscale
