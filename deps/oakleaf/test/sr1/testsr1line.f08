
module line_function

  use symmetric_rank_1
  use iso_fortran_env
  implicit none

  type, extends(sr1_function) :: sr1_line
     real(REAL64), dimension(:), pointer :: x,y
   contains
     procedure :: value => line_fval
     procedure :: gradient => line_fgrad
     procedure :: hessian => line_fhessian
     procedure :: print => line_print
  end type sr1_line

contains

  function line_fval(this,x) result(fval)

    use rfun

    real(REAL64) :: fval
    class(sr1_line) :: this
    real(REAL64), dimension(:), intent(in) :: x
!    real(REAL64), dimension(:), pointer :: x,y
    real(REAL64) :: a,b,s
    integer :: n

!    x => this%x
!    y => this%y
    n = size(this%x)
    a = x(1)
    b = x(2)
    s = x(3)

    associate( x => this%x, y => this%y )
!      fval = sum(itukey((y - (a + b*x))/s))/n + log(s)
      fval = sum(ihuber((y - (a + b*x))/s))/n + log(s)
!    fval = sum(((y - (a + b*x))/s)**2)/(2*n) + log(s)
    end associate

  end function line_fval

  subroutine line_fgrad(this,x,fgrad)

    use rfun

    class(sr1_line) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), intent(out) :: fgrad
!    real(REAL64), dimension(:), pointer :: x,y
    real(REAL64), dimension(:), allocatable :: res, psi
    real(REAL64) :: a, b, s
    integer :: n

!    x => this%x
!    y => this%y
    n = size(this%x)
    a = x(1)
    b = x(2)
    s = x(3)

    associate( x => this%x, y => this%y )

      allocate(res(n),psi(n))
      res = y - (a + b*x)
      psi = huber(res/s)
      !psi = tukey(res/s)

      fgrad(1) = -sum(psi)/(n*s)
      fgrad(2) = -sum(psi*x)/(n*s)
      fgrad(3) = -sum(psi*res)/(n*s**2) + 1/s
    end associate

!    fgrad(1) = -sum(res)/(n*s**2)
!    fgrad(2) = -sum(res*x)/(n*s**2)
!    fgrad(3) = -sum((res/s)**2)/(n*s) + 1/s

  end subroutine line_fgrad

  subroutine line_fhessian(this,x,fhess)

    use rfun

    class(sr1_line) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:,:), intent(out) :: fhess
    real(REAL64), dimension(:), allocatable :: res, psi, dpsi
    real(REAL64) :: a, b, s
    integer :: n, i, j

    n = size(this%x)
    a = x(1)
    b = x(2)
    s = x(3)

    associate( x => this%x, y => this%y )

      allocate(res(n),psi(n))
      res = y - (a + b*x)
      !psi = huber(res/s)
!      psi = tukey(res/s)
!      dpsi = dtukey(res/s)
      psi = huber(res/s)
      dpsi = dhuber(res/s)

      fhess(1,1) = sum(dpsi)/(n*s**2)
      fhess(1,2) = sum(dpsi*x)/(n*s**2)
      fhess(1,3) = (sum(dpsi*res)/s + sum(psi))/(n*s**2)
      fhess(2,2) = sum(dpsi*x**2)/(n*s**2)
      fhess(2,3) = (sum(dpsi*res*x)/s + sum(psi*x))/(n*s**2)
      fhess(3,3) = (sum(dpsi*res**2)/s + 2*sum(psi*res))/(n*s**3) - 1/s**2

      do i = 1, 3
         do j = i+1, 3
            fhess(j,i) = fhess(i,j)
         end do
      end do

    end associate

  end subroutine line_fhessian


  subroutine line_print(this,iter,x,fvalue,fgrad,fhess,step,delta,rho,sinfo)
    class(sr1_line) :: this
    integer, intent(in) :: iter, sinfo
    real(REAL64), intent(in) :: fvalue,delta,rho
    real(REAL64), dimension(:), intent(in) :: x,fgrad,step
    real(REAL64), dimension(:,:), intent(in) :: fhess

    write(2,*) iter,x,fvalue,delta,step

    ! rectangle
!!$    block
!!$      integer, parameter :: ngrid = 100
!!$      real(REAL64) :: u,v,h,f
!!$      real(REAL64), dimension(3) :: dr
!!$      integer :: i,j
!!$      character(len=666) :: filename
!!$
!!$      h = dlt / ngrid
!!$
!!$      write(filename,'(a,i0)') '/tmp/sr1line_approx_',iter
!!$      open(3,file=filename)
!!$      do i = -ngrid, ngrid
!!$         u = par(1) + i*h
!!$         do j = -ngrid, ngrid
!!$            v = par(2) + j*h
!!$            dr = [u - par(1), v - par(2), par(3)]
!!$            f = fval + dot_product(fgrad,dr) + 0.5*dot_product(dr,matmul(B,dr))
!!$            write(3,*) u,v,f
!!$         end do
!!$         write(3,*)
!!$      end do
!!$      close(3)
!!$    end block

  end subroutine line_print

end module line_function


program testsr1line

  ! test an optimization

  use line_function
  use symmetric_rank_1
  use noise
  use iso_fortran_env
  implicit none

  type(sr1_line) :: funline
  integer, parameter :: datamax = 6
  real(REAL64), dimension(:), allocatable, target :: x,y
  real(REAL64), dimension(3) :: p
  real(REAL64), dimension(3,3) :: cov
  real(REAL64) :: delta, maxdev
  integer :: maxfev, info, i
  character(len=80) :: msg
  logical :: verbose, quasi

  allocate(x(datamax),y(datamax))
  x = [ -4, -3, -2, -1, 0, 10 ]
  y = [ 2.48, 0.73, -0.04, -1.44, -1.32, 0.0 ]

  do i = 1, datamax*0
     x(i) = 1e-1*i
     y(i) = 2 + 1*x(i) + gnoise(0.0_REAL64,0.5_REAL64)
     if( mod(i,20) == 0 ) &
          y(i) = 2 + 1*x(i) + gnoise(0.0_REAL64,100.0_REAL64)
  end do
!  call gnoise_fill(0.0_REAL64,1.0_REAL64,data)
!  do i = 1, datamax,7
!     data(i) = gnoise(1.0_REAL64,10.0_REAL64)
!  end do

  funline%x => x
  funline%y => y

  p = [ -2.0, -1.0, 0.6 ]
  p = [ -1., -1., 1.0 ]
!  p = [ 1e10, 1e10, 1.0 ]
!  p = [0,0,1]
  delta = 1
  maxdev = 1e-7
  maxfev = 666
  verbose = .true.
  quasi = .false.
  open(2,file='/tmp/sr1line_step')

  call sr1rob(funline,p,cov,delta,maxdev,maxfev,quasi,verbose,info,msg)

  write(*,*) info,trim(msg)
  write(*,*) p

  close(2)

  call graph(funline,-5.0,5.0,-5.0,5.0,real(p(3)),'/tmp/sr1line')
  call graph_disp(funline,-5.0,5.0,-5.0,5.0,'/tmp/sr1line_disp')

contains

  subroutine graph(fcn,x0,x1,y0,y1,sig0,filename)

    class(sr1_function), intent(in) :: fcn
    real(REAL32), intent(in) :: x0,x1,y0,y1,sig0
    character(len=*), intent(in) :: filename
    integer, parameter :: nx = 500
    integer, parameter :: ny = 500
    real(REAL64) :: a,da,b,db,fval,sig
    integer :: i,j

    sig = sig0
!    sig = 0.6
!    sig = 0.1
    da = (x1 - x0) / (nx - 1)
    db = (y1 - y0) / (ny - 1)

    open(1,file=filename)
    do j = 1, ny
       b = y0 + db*j
       do i = 1, nx
          a = x0 + da*i
          fval = fcn%value([a,b,sig])
          write(1,*) a,b,fval,sum(abs(y - (a + b*x)),x > 3)/size(x)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph


  subroutine graph_disp(fcn,x0,x1,y0,y1,filename)

    use rfun

    class(sr1_function), intent(in) :: fcn
    real(REAL32), intent(in) :: x0,x1,y0,y1
    character(len=*), intent(in) :: filename
    integer, parameter :: nx = 500
    integer, parameter :: ny = 500
    real(REAL64), dimension(:), allocatable :: res, psi, dpsi
    real(REAL64) :: a,da,b,db,sig,disp,f2,df
    integer :: i,j,n

    sig = 0.6
    da = (x1 - x0) / (nx - 1)
    db = (y1 - y0) / (ny - 1)

    n = size(x)
    allocate(res(n),psi(n),dpsi(n))

    open(1,file=filename)
    do j = 1, ny
       b = y0 + db*j
       do i = 1, nx
          a = x0 + da*i
          res = (y - (a + b*x))
          psi = tukey(res/sig)
          dpsi = dtukey(res/sig)
!          psi = huber(res/sig)
!          dpsi = dhuber(res/sig)
          f2 = sum(psi**2)
          df = sum(dpsi)
          if( f2 > 0 ) then
             disp = df/sqrt(f2) / sig
!          if( df > 1e-5 ) then
!             disp = sig**2*f2/df**2
!             disp = 1/disp
          else
             disp = 0
          end if
          write(1,*) a,b,disp
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_disp


end program testsr1line
