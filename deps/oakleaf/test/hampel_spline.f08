
! gfortran -std=f2008 -Wall -g -p -no-pie -fcheck=all rfun_hampel.f08 spline.f seval.f hampel_spline.f08
!

! Spline interpolation of the function of Hampel

program hampel_spline

  use rfun_hampel_REAL64, ipshampel => ipshampel_REAL64, ihampel => ihampel_REAL64,&
       pshampel => pshampel_REAL64, dpshampel => dpshampel_REAL64, &
       ishampel => ishampel_REAL64, shampel => shampel_REAL64, &
       dshampel => dshampel_REAL64, hampelsqpf => hampelsqpf_REAL64, &
       hampel => hampel_REAL64, dhampel => dhampel_REAL64, ihampel => ihampel_REAL64
  use iso_fortran_env
  use, intrinsic :: ieee_arithmetic

  implicit none

  real(REAL64), dimension(:), allocatable :: x,y,b,c,d
  real(REAL64) :: u,t,f,g,seval,f1,g1,R,f2,g2,h,df,df2
  integer :: n, i, k
  external :: seval

  R = 0.275_REAL64
  N = nint(8.525_REAL64 / R) + 5 ! five is importnat for the boundary condition
  allocate(x(n),y(n),b(n),c(n),d(n))
  do i = 1, N
     x(i) = (i-1)*R
     y(i) = ipshampel(x(i))
  end do

  call spline(n,x,y,b,c,d)

  ! comment out the private in hampel's module
  do i = 0, -32
     write(*,*) i,d(i+1),spline_d(i),spline_d(i)-d(i+1)
!     write(*,*) i,c(i+1),spline_c(i),spline_c(i)-c(i+1)
!     write(*,*) i,b(i+1),spline_b(i),spline_b(i)-b(i+1)
!     write(*,*) i,y(i+1),spline_y(i),spline_y(i)-y(i+1)
  end do
  do i = 1, -N
     write(*,*) i,b(i),c(i),d(i),y(i)
  end do
!  stop

  do i = 1, -10*n
     call random_number(u)
!     t = 20*(u - 0.5)
     t = 0.1*R*(i-1) + 0.1*R*u
     f = ipshampel(t)
     f1 = pshampel(t)
     f2 = dpshampel(t)
     g = seval(n,t,x,y,b,c,d)
     k = int(t/R) + 1
     k = min(k,n)
     g = y(k) + (t - x(k))*(b(k) + (t - x(k))*(c(k) + (t - x(k))*d(k)))
     g1 = b(k) + (t - x(k))*(2*c(k) + 3*d(k)*((t - x(k))))
     g2 = 2*c(k) + 6*d(k)*((t - x(k)))
     write(*,'(10f25.16,i5)') t,f,g,f-g,f1,g1,f1-g1,f2,g2,f2-g2,k
  end do

  !  do k = 1, 1000
  h = 1e-3
  do i = -1000*n, 1000*n
     call random_number(u)
     t = h*R*(i-1) + h*R*u
     f = ipshampel(t)
     f1 = pshampel(t)
     f2 = dpshampel(t)
     df = (ipshampel(t+h) - ipshampel(t-h)) / (2*h)
     df2 = (pshampel(t+h) - pshampel(t-h)) / (2*h)
     g = ishampel(t)
     g1 = shampel(t)
     g2 = dshampel(t)
     call hampelsqpf(t,g,g1)
     write(*,'(14f25.16)') t,f,g,f-g,f1,g1,f1-g1,f2,g2,f2-g2,df,df2,f1-df,f2-df2
  end do
!enddo

  write(*,*) '#', min(abs(IEEE_VALUE(t, IEEE_QUIET_NAN)),10.0_REAL64)
  write(*,*) '#',ishampel(IEEE_VALUE(t, IEEE_POSITIVE_INF))
  write(*,*) '#',ishampel(IEEE_VALUE(t, IEEE_NEGATIVE_INF))
  write(*,*) '#',ishampel(IEEE_VALUE(t, IEEE_QUIET_NAN))
  write(*,*) '#',ishampel(10.0_REAL64)

end program hampel_spline
