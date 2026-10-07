!!$
!!$module robust_data_function
!!$
!!$  use rfun
!!$  use iso_fortran_env
!!$  implicit none
!!$
!!$  type, abstract, public :: robust_datafun
!!$     logical :: break = .false.
!!$     procedure(varrho), nopass, pointer :: rho
!!$     procedure(varpsi), nopass, pointer :: psi
!!$     procedure(vardpsi), nopass, pointer :: dpsi
!!$     real(REAL64), dimension(:), pointer :: data
!!$   contains
!!$     procedure(fvalue), deferred :: val
!!$     procedure(fgradient), deferred :: grad
!!$     procedure(fhessian), deferred :: hess
!!$     procedure(fstderr), deferred :: stderr
!!$     procedure(fprint), deferred :: print
!!$  end type robust_datafun
!!$
!!$
!!$contains
!!$
!!$  elemental pure function varrho(x)
!!$
!!$    real(REAL64) :: varrho
!!$    real(REAL64), intent(in) :: x
!!$
!!$  end function varrho
!!$
!!$  elemental pure function varpsi(x)
!!$
!!$    real(REAL64) :: varpsi
!!$    real(REAL64), intent(in) :: x
!!$
!!$  end function varpsi
!!$
!!$  elemental pure function vardpsi(x)
!!$
!!$    real(REAL64) :: vardpsi
!!$    real(REAL64), intent(in) :: x
!!$
!!$  end function vardpsi
!!$
!!$  function fvalue(this,t,s) result(fval)
!!$
!!$    real(REAL64) :: fval
!!$    class(robust_datafun) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$
!!$    procedure(varrho), pointer :: rho
!!$    real(REAL64), dimension(:), pointer :: data
!!$    integer :: n
!!$
!!$    rho => this%rho
!!$    data => this%data
!!$    n = size(data)
!!$
!!$    fval = sum(rho((data - t)/s))/n + log(s)
!!$
!!$  end function fvalue
!!$
!!$  function fstderr(this,t,s) result(stderr)
!!$
!!$    real(REAL64) :: stderr
!!$    class(robust_datafun) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$
!!$    procedure(varpsi), pointer :: psi
!!$    procedure(vardpsi), pointer :: dpsi
!!$    real(REAL64), dimension(:), pointer :: data
!!$    real(REAL64) :: f2,df,f
!!$    integer :: n
!!$
!!$    psi => this%psi
!!$    dpsi => this%dpsi
!!$    data => this%data
!!$    n = size(data)
!!$
!!$    f2 = sum(psi((data - t)/s)**2)
!!$    df = sum(dpsi((data - t)/s))
!!$
!!$    if( f2 > 0 .and. abs(df) > 0) then
!!$       f = sqrt(f2)
!!$!       stderr = s*sqrt(f2/df**2*n/(n-1))
!!$       stderr = s*(f/(n-1))*(n/df)
!!$    else
!!$       stderr = 0
!!$    end if
!!$
!!$  end function fstderr
!!$
!!$  subroutine fgradient(this,t,s,fgrad)
!!$
!!$    class(robust_datafun) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64), dimension(:), intent(out) :: fgrad
!!$
!!$    procedure(varpsi), pointer :: psi
!!$    real(REAL64), dimension(:), pointer :: data
!!$    real(REAL64), dimension(:), allocatable :: r
!!$    integer :: n
!!$
!!$    psi => this%psi
!!$    data => this%data
!!$    n = size(data)
!!$
!!$    allocate(r(n))
!!$    r = (data - t) / s
!!$
!!$    fgrad(1) = -sum(psi(r))/(n*s)
!!$    fgrad(2) = -sum(psi(r)*r)/(n*s) + 1/s
!!$
!!$  end subroutine fgradient
!!$
!!$  subroutine fhessian(this,t,s,fhess)
!!$
!!$    class(robust_datafun) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64), dimension(:,:), intent(out) :: fhess
!!$
!!$    procedure(varpsi), pointer :: psi
!!$    procedure(vardpsi), pointer :: dpsi
!!$    real(REAL64), dimension(:), pointer :: data
!!$    real(REAL64), dimension(:), allocatable :: r, p, dp
!!$    integer :: n
!!$
!!$    psi => this%psi
!!$    dpsi => this%dpsi
!!$    data => this%data
!!$    n = size(data)
!!$
!!$    allocate(r(n),p(n),dp(n))
!!$    r = (data - t) / s
!!$    p = psi(r)
!!$    dp = dpsi(r)
!!$
!!$    fhess(1,1) = sum(dp) / (n*s**2)
!!$    fhess(1,2) = (sum(dp*r) + sum(p)/s) / (n*s)
!!$    fhess(2,1) = fhess(1,2)
!!$    fhess(2,2) = (sum(dp*r**2) + 2*sum(p*r))/(n*s**2) - 1/s**2
!!$
!!$  end subroutine fhessian
!!$
!!$
!!$  subroutine fprint(this,iter,t,s,fvalue,fgrad,fhess,step,delta,rho,sinfo)
!!$    class(robust_datafun) :: this
!!$    integer, intent(in) :: iter, sinfo
!!$    real(REAL64), intent(in) :: t,s,fvalue,delta,rho
!!$    real(REAL64), dimension(:), intent(in) :: fgrad,step
!!$    real(REAL64), dimension(:,:), intent(in) :: fhess
!!$
!!$    write(2,*) iter,t,s,fvalue,delta,step
!!$
!!$  end subroutine fprint
!!$
!!$  subroutine fline(this,t,s,g1,h11)
!!$
!!$    class(robust_datafun) :: this
!!$
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64), intent(out) :: g1, h11
!!$
!!$    procedure(varpsi), pointer :: psi
!!$    procedure(vardpsi), pointer :: dpsi
!!$    real(REAL64), dimension(:), pointer :: data
!!$    real(REAL64), dimension(:), allocatable :: r
!!$    integer :: n
!!$
!!$    psi => this%psi
!!$    dpsi => this%dpsi
!!$    data => this%data
!!$    n = size(data)
!!$
!!$    allocate(r(n))
!!$    r = (data - t) / s
!!$
!!$    g1 = -sum(this%psi(r)) / (n*s)
!!$    h11 = sum(this%dpsi(r)) / (n*s**2)
!!$
!!$  end subroutine fline
!!$
!!$end module robust_data_function
!!$


module sr1mean_function

  use iso_fortran_env
  implicit none

  type, public :: tsr1mean_function
     logical :: break = .false.
     real(REAL64), dimension(:), pointer :: data
     real(REAL64), dimension(:), allocatable :: res
   contains
     procedure :: init => init_function
     procedure :: val => square_val
     procedure :: grad => square_fgrad
     procedure :: hess => square_fhessian
     procedure :: stderr => square_stderr
     procedure :: print => square_print
!     procedure :: g1 => square_g1
!     procedure :: h11 => square_h11
     procedure :: rnewton => square_rnewton
     procedure :: nfrac => square_nfrac
     final :: final_function
  end type tsr1mean_function

  abstract interface

     real(REAL64) function fvalue(this,t,s)
       import :: REAL64, tsr1mean_function
       class(tsr1mean_function) :: this
       real(REAL64), intent(in) :: t,s
     end function fvalue

     subroutine fgradient(this,t,s,fgrad)
       import :: REAL64, tsr1mean_function
       class(tsr1mean_function) :: this
       real(REAL64), intent(in) :: t,s
       real(REAL64), dimension(:), intent(out) :: fgrad
     end subroutine fgradient

     subroutine fhessian(this,t,s,fhess)
       import :: REAL64, tsr1mean_function
       class(tsr1mean_function) :: this
       real(REAL64), intent(in) :: t,s
       real(REAL64), dimension(:,:), intent(out) :: fhess
     end subroutine fhessian

     real(REAL64) function fstderr(this,t,s)
       import :: REAL64, tsr1mean_function
       class(tsr1mean_function) :: this
       real(REAL64), intent(in) :: t,s
     end function fstderr

     subroutine fprint(this,iter,t,s,fvalue,fgrad,fhess,step,delta,rho,sinfo)
       import :: REAL64, tsr1mean_function
       class(tsr1mean_function) :: this
       integer, intent(in) :: iter, sinfo
       real(REAL64), intent(in) :: t,s,fvalue,delta,rho
       real(REAL64), dimension(:), intent(in) :: fgrad,step
       real(REAL64), dimension(:,:), intent(in) :: fhess
     end subroutine fprint

  end interface

  type, extends(tsr1mean_function) :: tsr1mean_tukey
   contains
     procedure :: val => tukey_fval
     procedure :: grad => tukey_fgrad
     procedure :: hess => tukey_fhessian
     procedure :: stderr => tukey_stderr
     procedure :: print => tukey_print
!     procedure :: g1 => tukey_g1
!     procedure :: h11 => tukey_h11
     procedure :: rnewton => tukey_rnewton
     procedure :: nfrac => tukey_nfrac
  end type tsr1mean_tukey

  type, extends(tsr1mean_function) :: tsr1mean_huber
!     real(REAL64), dimension(:), pointer :: data
   contains
     procedure :: val => huber_fval
     procedure :: grad => huber_fgrad
     procedure :: hess => huber_fhessian
     procedure :: stderr => huber_stderr
     procedure :: print => huber_print
     procedure :: nfrac => huber_nfrac
  end type tsr1mean_huber

!!$  type, extends(tsr1mean_function) :: tsr1mean_square
!!$!     real(REAL64), dimension(:), pointer :: data
!!$   contains
!!$     procedure :: val => square_fval
!!$     procedure :: grad => square_fgrad
!!$     procedure :: hess => square_fhessian
!!$     procedure :: stderr => square_stderr
!!$     procedure :: print => square_print
!!$  end type tsr1mean_square

  type, extends(tsr1mean_function) :: tsr1mean_experimental
   contains
     procedure :: val => expr_fval
     procedure :: grad => expr_fgrad
     procedure :: hess => expr_fhessian
     procedure :: stderr => expr_stderr
     procedure :: print => expr_print
!     procedure :: g1 => expr_g1
!     procedure :: h11 => expr_h11
     procedure :: rnewton => expr_rnewton
     procedure :: nfrac => expr_nfrac
  end type tsr1mean_experimental

  type, extends(tsr1mean_function) :: tsr1mean_entropy
   contains
     procedure :: val => entropy_fval
     procedure :: grad => entropy_fgrad
!     procedure :: hess => entropy_fhessian
!     procedure :: stderr => eentropy_stderr
!     procedure :: print => entropy_print
!     procedure :: g1 => expr_g1
!     procedure :: h11 => expr_h11
!     procedure :: rnewton => expr_rnewton
!     procedure :: nfrac => expr_nfrac
  end type tsr1mean_entropy


  private :: final_function

contains

  subroutine init_function(this,data,iflag)

    class(tsr1mean_function) :: this
    real(REAL64), dimension(:), pointer, intent(in) :: data
    integer, intent(out) :: iflag

!    integer :: n
    character(len=81) :: errmsg

    this%data => data
!    n = size(data)

    allocate(this%res(size(data)),stat=iflag,errmsg=errmsg)
    if( iflag /= 0 ) then
       write(error_unit,*) 'init(tsr1mean): ',trim(errmsg)
       return
    end if

  end subroutine init_function

  subroutine final_function(this)

    type(tsr1mean_function) :: this

    if( allocated(this%res) ) deallocate(this%res)
!    write(*,*) '... Finished.'

  end subroutine final_function

  function tukey_fval(this,t,s) result(fval)

    use rfun

    real(REAL64) :: fval
    class(tsr1mean_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), pointer :: data
    integer :: n

    data => this%data
    n = size(data)

    fval = sum(itukey((data - t)/s))/n + log(s)

  end function tukey_fval

  function tukey_stderr(this,t,s,reli) result(stderr)

    use rfun

    real(REAL64) :: stderr
    class(tsr1mean_tukey) :: this
    real(REAL64), intent(in) :: t,s
    logical, intent(out) :: reli

    real(REAL64), dimension(:), pointer :: data
    real(REAL64) :: f2,df,f
    integer :: n

    data => this%data
    n = size(data)

    f2 = sum(tukey((data - t)/s)**2)
    df = sum(dtukey((data - t)/s))

!    if( f2 > 0 .and. df > 0.1*n ) then
    if( f2 > 0 .and. df > 0 ) then
       f = sqrt(f2)
!       stderr = s*sqrt(f2/df**2*n/(n-1))
       stderr = s*(f/(n-1))*(n/df)
       reli = .true.
    else
       stderr = 0!epsilon(stderr)
       reli = .false.
    end if

  end function tukey_stderr

  subroutine tukey_fgrad(this,t,s,fgrad)

    use rfun

    class(tsr1mean_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: fgrad

    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res
    integer :: n

    data => this%data
    n = size(data)

    allocate(res(n))
    res = (data - t) / s

    fgrad(1) = -sum(tukey(res))/(n*s)
    fgrad(2) = -sum(tukey(res)*res)/(n*s) + 1/s

  end subroutine tukey_fgrad

  subroutine tukey_fhessian(this,t,s,fhess)

    use rfun

    class(tsr1mean_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: fhess

    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res, psi, dpsi
    integer :: n

    data => this%data
    n = size(data)

    allocate(res(n),psi(n),dpsi(n))
    res = (data - t) / s
    psi = tukey(res)
    dpsi = dtukey(res)

    fhess(1,1) = sum(dpsi) / s**2 / n
    fhess(1,2) = (sum(dpsi*res)/s + sum(psi)/s**2) / n
    fhess(2,1) = fhess(1,2)
    fhess(2,2) = (sum(dpsi*res**2) + 2*sum(psi*res))/s**2/n - 1/s**2

  end subroutine tukey_fhessian


  subroutine tukey_print(this,iter,t,s,fvalue,fgrad,fhess,step,delta,rho,sinfo)
    class(tsr1mean_tukey) :: this
    integer, intent(in) :: iter, sinfo
    real(REAL64), intent(in) :: t,s,fvalue,delta,rho
    real(REAL64), dimension(:), intent(in) :: fgrad,step
    real(REAL64), dimension(:,:), intent(in) :: fhess

    write(2,*) iter,t,s,fvalue,delta,step

  end subroutine tukey_print


!!$  function tukey_g1(this,t,s) result(g1)
!!$
!!$    use rfun
!!$    class(tsr1mean_tukey) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64) :: g1
!!$    integer :: n
!!$
!!$    n = size(this%data)
!!$    this%res = (this%data - t) / s
!!$    g1 = -sum(tukey(this%res))/(n*s)
!!$
!!$  end function tukey_g1
!!$
!!$  function tukey_h11(this,t,s) result(h11)
!!$
!!$    use rfun
!!$    class(tsr1mean_tukey) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$
!!$    real(REAL64) :: h11
!!$    integer :: n
!!$
!!$    n = size(this%data)
!!$    this%res = (this%data - t) / s
!!$    h11 = sum(dtukey(this%res)) / (n*s**2)
!!$
!!$  end function tukey_h11

  subroutine tukey_rnewton(this,t,s,g1,h11)

    use rfun
    class(tsr1mean_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: g1,h11
    integer :: n

    n = size(this%data)
    this%res = (this%data - t) / s
    g1 = -sum(tukey(this%res)) / (n*s)
    h11 = sum(dtukey(this%res)) / (n*s**2)

  end subroutine tukey_rnewton

  function tukey_nfrac(this,t,s) result(nfrac)

    use rfun

    real(REAL64) :: nfrac
    class(tsr1mean_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), pointer :: data
    integer :: n

    data => this%data
    n = size(data)

    nfrac = sum(dtukey((data - t)/s))/n

  end function tukey_nfrac



  function huber_fval(this,t,s) result(fval)

    use rfun

    real(REAL64) :: fval
    class(tsr1mean_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), pointer :: data
    integer :: n

    data => this%data
    n = size(data)

    fval = sum(ipshuber((data - t)/s))/n + log(s)

  end function huber_fval

  subroutine huber_fgrad(this,t,s,fgrad)

    use rfun

    class(tsr1mean_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: fgrad

    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res
    integer :: n

    data => this%data
    n = size(data)

    allocate(res(n))
    res = (data - t) / s

    fgrad(1) = -sum(pshuber(res))/(n*s)
    fgrad(2) = -sum(pshuber(res)*res)/(n*s) + 1/s

  end subroutine huber_fgrad


  subroutine huber_fhessian(this,t,s,fhess)

    use rfun

    class(tsr1mean_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: fhess

    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res, psi, dpsi
    integer :: n

    data => this%data
    n = size(data)

    allocate(res(n),psi(n),dpsi(n))
    res = (data - t) / s
    psi = pshuber(res)
    dpsi = dpshuber(res)

    fhess(1,1) = sum(dpsi) / s**2 / n
    fhess(1,2) = (sum(dpsi*res)/s + sum(psi)/s**2) / n
    fhess(2,1) = fhess(1,2)
    fhess(2,2) = (sum(dpsi*res**2) + 2*sum(psi*res))/s**2/n - 1/s**2

  end subroutine huber_fhessian

  subroutine huber_print(this,iter,t,s,fvalue,fgrad,fhess,step,delta,rho,sinfo)
    class(tsr1mean_huber) :: this
    integer, intent(in) :: iter, sinfo
    real(REAL64), intent(in) :: t,s,fvalue,delta,rho
    real(REAL64), dimension(:), intent(in) :: fgrad,step
    real(REAL64), dimension(:,:), intent(in) :: fhess


    write(2,*) iter,t,s,fvalue,delta,step

  end subroutine huber_print

  function huber_stderr(this,t,s,reli) result(stderr)

    use rfun

    real(REAL64) :: stderr
    class(tsr1mean_huber) :: this
    real(REAL64), intent(in) :: t,s
    logical, intent(out) :: reli

    real(REAL64), dimension(:), pointer :: data
    real(REAL64) :: f2,df,f
    integer :: n

    data => this%data
    n = size(data)

    f2 = sum(pshuber((data - t)/s)**2)
    df = sum(dpshuber((data - t)/s)**2)

    if( f2 > 0 .and. df > 0 ) then
       f = sqrt(f2)
       stderr = s*sqrt(f2/df**2*n/(n-1))
       reli = .true.
    else
       stderr = 0
       reli = .false.
    end if

  end function huber_stderr

  function huber_nfrac(this,t,s) result(nfrac)

    use rfun

    real(REAL64) :: nfrac
    class(tsr1mean_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), pointer :: data
    integer :: n

    data => this%data
    n = size(data)

    nfrac = sum(dhuber((data - t)/s))/n

  end function huber_nfrac


  function square_val(this,t,s) result(fval)

    use rfun

    real(REAL64) :: fval
    class(tsr1mean_function) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), pointer :: data
    integer :: n

    data => this%data
    n = size(data)

    fval = sum(((data - t)/s)**2)/(2*n) + log(s)

  end function square_val

!!$  function square_g1(this,t,s) result(g1)
!!$
!!$    use rfun
!!$    class(tsr1mean_function) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64) :: g1
!!$    integer :: n
!!$
!!$    n = size(this%data)
!!$    this%res = (this%data - t) / s
!!$    g1 = -sum(this%res)/(n*s)
!!$
!!$  end function square_g1
!!$
!!$  function square_h11(this,t,s) result(h11)
!!$
!!$ !   use rfun
!!$    class(tsr1mean_function) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$
!!$    real(REAL64), dimension(:), allocatable :: dpsi
!!$    real(REAL64) :: h11
!!$    integer :: n
!!$
!!$    n = size(this%data)
!!$!    res = (this%data - t) / s
!!$    allocate(dpsi(n))
!!$    dpsi = 1
!!$    h11 = sum(dpsi) / s**2 / n
!!$
!!$  end function square_h11

  subroutine square_rnewton(this,t,s,g1,h11)

!    use rfun
    class(tsr1mean_function) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: g1,h11
    integer :: n

    n = size(this%data)
    this%res = (this%data - t) / s
    g1 = -sum(this%res) / (n*s)
    h11 = 1 / s**2

  end subroutine square_rnewton

  subroutine square_fgrad(this,t,s,fgrad)

    use rfun

    class(tsr1mean_function) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: fgrad

    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res
    integer :: n

    data => this%data
    n = size(data)

    allocate(res(n))
    res = (data - t) / s

    fgrad(1) = -sum(res)/(n*s)
    fgrad(2) = -sum(res**2)/(n*s) + 1/s

  end subroutine square_fgrad

  subroutine square_fhessian(this,t,s,fhess)

    use rfun

    class(tsr1mean_function) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: fhess

    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res, psi, dpsi
    integer :: n

    data => this%data
    n = size(data)

    allocate(res(n),psi(n),dpsi(n))
    res = (data - t) / s
    psi = res
    dpsi = 1

    fhess(1,1) = sum(dpsi) / s**2 / n
    fhess(1,2) = (sum(dpsi*res)/s + sum(psi)/s**2) / n
    fhess(2,1) = fhess(1,2)
    fhess(2,2) = (sum(dpsi*res**2) + 2*sum(psi*res))/s**2/n - 1/s**2

  end subroutine square_fhessian

  subroutine square_print(this,iter,t,s,fvalue,fgrad,fhess,step,delta,rho,sinfo)
    class(tsr1mean_function) :: this
    integer, intent(in) :: iter, sinfo
    real(REAL64), intent(in) :: t,s,fvalue,delta,rho
    real(REAL64), dimension(:), intent(in) :: fgrad,step
    real(REAL64), dimension(:,:), intent(in) :: fhess

    write(2,*) iter,t,s,fvalue,delta,step

  end subroutine square_print

  function square_stderr(this,t,s,reli) result(stderr)

    use rfun

    real(REAL64) :: stderr
    class(tsr1mean_function) :: this
    real(REAL64), intent(in) :: t,s
    logical, intent(out) :: reli
    real(REAL64), dimension(:), pointer :: data
    real(REAL64) :: f2,df,f
    integer :: n

    data => this%data
    n = size(data)

    f2 = sum(((data - t)/s)**2)
    df = n

    if( f2 > 0 .and. df > 0 ) then
       f = sqrt(f2)
       stderr = s*sqrt(f2/df**2*n/(n-1))
       reli = .true.
    else
       stderr = 0
       reli = .false.
    end if

  end function square_stderr

  function square_nfrac(this,t,s) result(nfrac)

    use rfun

    real(REAL64) :: nfrac
    class(tsr1mean_function) :: this
    real(REAL64), intent(in) :: t,s

    nfrac = 1.0_REAL64

  end function square_nfrac



  function expr_fval(this,t,s) result(fval)

    use rfun_hampel

    real(REAL64) :: fval
    class(tsr1mean_experimental) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), pointer :: data
    real(REAL64) :: r
    integer :: n,i

    data => this%data
    n = size(data)
    fval = sum(ipshampel((data - t)/s))/n + log(s)
!    fval = sum(ibathtub((data - t)/s))/n + log(s)
!    fval = sum(irat((data - t)/s))/n + log(s)

  end function expr_fval

  subroutine expr_fgrad(this,t,s,fgrad)

    use rfun_hampel

    class(tsr1mean_experimental) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: fgrad

    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res
    integer :: n

    data => this%data
    n = size(data)

    allocate(res(n))
    res = (data - t) / s

    fgrad(1) = -sum(pshampel(res))/(n*s)
    fgrad(2) = -sum(pshampel(res)*res)/(n*s) + 1/s

!    fgrad(1) = -sum(bathtub(res))/(n*s)
!    fgrad(2) = -sum(bathtub(res)*res)/(n*s) + 1/s
!    fgrad(1) = -sum(rat(res))/(n*s)
!    fgrad(2) = -sum(rat(res)*res)/(n*s) + 1/s

  end subroutine expr_fgrad

  subroutine expr_fhessian(this,t,s,fhess)

    use rfun_hampel

    class(tsr1mean_experimental) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: fhess

    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res, psi, dpsi
    integer :: n

    data => this%data
    n = size(data)

    allocate(res(n),psi(n),dpsi(n))
    res = (data - t) / s
    psi = pshampel(res)
    dpsi = dpshampel(res)
!    psi = bathtub(res)
!    dpsi = dbathtub(res)
!    psi = rat(res)
!    dpsi = drat(res)

    fhess(1,1) = sum(dpsi) / s**2 / n
    fhess(1,2) = (sum(dpsi*res)/s + sum(psi)/s**2) / n
    fhess(2,1) = fhess(1,2)
    fhess(2,2) = (sum(dpsi*res**2) + 2*sum(psi*res))/s**2/n - 1/s**2

  end subroutine expr_fhessian

  subroutine expr_print(this,iter,t,s,fvalue,fgrad,fhess,step,delta,rho,sinfo)
    class(tsr1mean_experimental) :: this
    integer, intent(in) :: iter, sinfo
    real(REAL64), intent(in) :: t,s,fvalue,delta,rho
    real(REAL64), dimension(:), intent(in) :: fgrad,step
    real(REAL64), dimension(:,:), intent(in) :: fhess

    write(2,*) iter,t,s,fvalue,delta,step

  end subroutine expr_print

  function expr_stderr(this,t,s,reli) result(stderr)

    use rfun_hampel

    real(REAL64) :: stderr
    class(tsr1mean_experimental) :: this
    real(REAL64), intent(in) :: t,s
    logical, intent(out) :: reli

    real(REAL64), dimension(:), pointer :: data
    real(REAL64) :: f2,df!,f
    integer :: n

    data => this%data
    n = size(data)

    f2 = sum(pshampel((data - t)/s)**2)
!    s0 = f2 / (n - 1)
    df = sum(dpshampel((data - t)/s))
!    nfrac = df / n
!    f2 = sum(bathtub((data - t)/s)**2)
!    df = sum(dbathtub((data - t)/s))
!    f2 = sum(rat((data - t)/s)**2)
!    df = sum(drat((data - t)/s))

    if( f2 > 0 .and. df > 0 ) then
!       f = sqrt(f2)
       stderr = s*sqrt(f2/df**2*n/(n-1))
!       stderr = s*sqrt(s0/nfrac/df)
       reli = .true.
    else
       stderr = 0
       reli = .false.
    end if

  end function expr_stderr

!!$  function expr_g1(this,t,s) result(g1)
!!$
!!$    use rfun
!!$    class(tsr1mean_experimental) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64) :: g1
!!$    integer :: n
!!$
!!$    n = size(this%data)
!!$    this%res = (this%data - t) / s
!!$    g1 = -sum(pshampel(this%res))/(n*s)
!!$
!!$  end function expr_g1
!!$
!!$  function expr_h11(this,t,s) result(h11)
!!$
!!$    use rfun
!!$    class(tsr1mean_experimental) :: this
!!$    real(REAL64), intent(in) :: t,s
!!$
!!$!    real(REAL64), dimension(:), allocatable :: dpsi
!!$    real(REAL64) :: h11
!!$    integer :: n
!!$
!!$    n = size(this%data)
!!$    this%res = (this%data - t) / s
!!$!    allocate(dpsi(n))
!!$!    dpsi = dpshampel(this%res)
!!$!    h11 = sum(dpsi) / s**2 / n
!!$    h11 = sum(dpshampel(this%res)) / (n*s**2)
!!$
!!$  end function expr_h11

  subroutine expr_rnewton(this,t,s,g1,h11)

    use rfun_hampel
    class(tsr1mean_experimental) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: g1,h11
    integer :: n

    n = size(this%data)
    this%res = (this%data - t) / s
    g1 = -sum(pshampel(this%res)) / (n*s)
    h11 = sum(dpshampel(this%res)) / (n*s**2)

  end subroutine expr_rnewton

  function expr_nfrac(this,t,s) result(nfrac)

    use rfun_hampel

    real(REAL64) :: nfrac
    class(tsr1mean_experimental) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), pointer :: data
    integer :: n

    data => this%data
    n = size(data)

!    write(*,*) 'nfrac:',t,s,n,sum(dpshampel((data - t)/s))
    if( s > 0 .and. n > 0 ) then
       nfrac = sum(dpshampel((data - t)/s))/n
    else
       nfrac = 0
    end if

  end function expr_nfrac


  function entropy_fval(this,t,s) result(fval)

    use rfun

    real(REAL64) :: fval
    class(tsr1mean_entropy) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res, f
    real(REAL64), parameter :: gamma = 2.650
    real(REAL64) :: r
    integer :: n,i

    data => this%data
    n = size(data)
    allocate(res(n),f(n))

    res = (data - t) / s
    f = exp(-itukey(res))/ (gamma*s)

    fval = sum(f*log(f),f>0) / n

  end function entropy_fval

  subroutine entropy_fgrad(this,t,s,fgrad)

    use rfun

    class(tsr1mean_entropy) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: fgrad

    real(REAL64), parameter :: gamma = 2.650
    real(REAL64), dimension(:), pointer :: data
    real(REAL64), dimension(:), allocatable :: res, q,f,phi
    integer :: n

    data => this%data
    n = size(data)

    allocate(res(n),q(n),phi(n),f(n))
    res = (data - t) / s
    phi = tukey(res)
    f = exp(-itukey(res))/ (gamma*s)
    where ( f > 0 )
       q = (log(f) + 1)*f
    elsewhere
       q = 0
    end where

!    fgrad(1) = -sum(q*phi) / (n*s)
!    fgrad(2) = -sum(q*(1 + phi*res)) / (n*s)

!    block
!      real(REAL64) :: h = 1e-7
!      fgrad(1) = (entropy_fval(this,t+h,s) - entropy_fval(this,t,s)) / h
!      fgrad(2) = (entropy_fval(this,t,s+h) - entropy_fval(this,t,s)) / h
!      write(*,*) 'numeric:',real(fgrad)
!    end block

    fgrad = [ sum(q*phi), sum(q*(phi*res - 1)) ] / (n*s)
!    write(*,*) 'analytic:',real(fgrad)


!    stop 'in entropy_fgrad'

  end subroutine entropy_fgrad


end module sr1mean_function
