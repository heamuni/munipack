

program hampel

  use quantilmean
  use oakleaf

  implicit none
  integer, parameter :: kind = selected_real_kind(15)
  real(kind), dimension(:), allocatable :: x, xcdf, ycdf, fcdf,tcdf
  real(kind), dimension(2,2) :: par
  real(kind), dimension(1) :: eps
  integer :: m, n
  real(kind) :: d,dh,t,dt,s

  eps(1) = 0.0
  par(1,:) = [0.0_kind, 1.0_kind]
  par(2,:) = [1.0_kind, 10.0_kind]

!  eps = 0.1
  do m = 1, 20
     n = 2**m
     x = gsample(n,eps,par)
     allocate(xcdf(n),ycdf(n),fcdf(n))
     call cdf(x,xcdf,ycdf)
     fcdf = gnormal(xcdf,par(1,1),par(1,2))
     d = maxval(abs(fcdf - ycdf))
     call rmean(x,t,dt,s,verbose=.true.)
     tcdf = gnormal(xcdf,t,s)
     dh = maxval(abs(fcdf - tcdf))
     write(*,*) n,d,dh,t,s
     deallocate(x,fcdf,xcdf,ycdf,tcdf)
  end do

  do m = 1, 20
     n = 2**m
     x = lsample(n,eps,par)
     allocate(xcdf(n),ycdf(n),fcdf(n))
     call cdf(x,xcdf,ycdf)
     fcdf = lnormal(xcdf,par(1,1),par(1,2))
     d = maxval(abs(fcdf - ycdf))
     call mmean(x,t,s)
     tcdf = lnormal(xcdf,t,s)
     dh = maxval(abs(fcdf - tcdf))
     write(*,*) n,d,dh,t,s
     deallocate(x,fcdf,xcdf,ycdf,tcdf)
  end do


contains

  function gsample(n,eps,par)

    use noise

    real(kind), dimension(:), allocatable :: gsample
    integer, intent(in) :: n
    real(kind), dimension(:), intent(in) :: eps
    real(kind), dimension(:,:), intent(in) :: par
    real(kind) :: p
    integer :: i,k

    allocate(gsample(n))
    do i = 1, n
       do k = 1, size(eps)
          call random_number(p)
          if( p > eps(k) ) then
             gsample(i) = gnoise(par(1,1),par(1,2))  ! primary
          else
             gsample(i) = gnoise(par(2,1),par(2,2))  ! contamination
          end if
       end do
    end do

  end function gsample


  function lsample(n,eps,par)

    use noise

    real(kind), dimension(:), allocatable :: lsample
    integer, intent(in) :: n
    real(kind), dimension(:), intent(in) :: eps
    real(kind), dimension(:,:), intent(in) :: par
    real(kind) :: p
    integer :: i,k

    allocate(lsample(n))
    do i = 1, n
       do k = 1, size(eps)
          call random_number(p)
          if( p > eps(k) ) then
             lsample(i) = lnoise(par(1,1),par(1,2))
          else
             lsample(i) = lnoise(par(2,1),par(2,2))
          end if
       end do
    end do

  end function lsample

  elemental function gnormal(x,t,s)

    real(kind), parameter :: sqrt2 = sqrt(2.0_kind)

    real(kind), intent(in) :: x,t,s
    real(kind) :: gnormal

    gnormal = (1.0_kind + erf((x-t)/(sqrt2*s)))/2.0_kind

  end function gnormal

  elemental function lnormal(x,t,s)

    real(kind), parameter :: sqrt2 = sqrt(2.0_kind)

    real(kind), intent(in) :: x,t,s
    real(kind) :: lnormal
    real(kind) :: u

    u = (x - t) / s
    if( x < t ) then
       lnormal = exp(u) / 2.0_kind
    else
       lnormal = 1.0_kind - exp(-u) / 2.0_kind
    end if

  end function lnormal

end program hampel
