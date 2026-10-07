!
!  A regression with Huber's robust function on base of BFGS minimiser
!
!  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Oakleaf.
!
!  Oakleaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Oakleaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!
! The subroutine rfit() implements robust estimations. As it utilises
! the bfgs method, it is suitable *** only for Huber function ***,
! and *** have not to be used for any non-increasing functions ***
! like Hampel or Tukey functions.
!
! The implementation is a wrapper for bfgs(). The scale is added separately
! to parameters of the fitting function; the methods of the derived type
! bfgs_robfun does it. The important characteristics of the algorithm
! is keeping scale s > 0, else log(s) or psi(r/s) in the likelihood function
! lead to a crash. It is ensured by maxalpha() implementation, as it,
! in every step, compute the step length, and if the step is infeasible,
! linear search is restricted on safe ranges --- it is possible only if
! the minimum never lie on the border; it is a trivial variant of
! interior-point methods. This condition s > 0 in linear search releases
! the most important drawback of Levendberg-Marquart method: the linear
! search in LM has no limit, making LM method unsuitable for the robust fitting.
!
! A handling of leverage points and winsorising are unimplemented.

module huber_curvefit

  use iso_fortran_env
  implicit none
  private

  type, abstract, public :: huberfit_function
   contains
     procedure(fvalue), deferred :: val
     procedure(fgradient), deferred :: grad
     procedure :: hess => fhessian
  end type huberfit_function

  abstract interface

     function fvalue(this,x,p) result(fval)
       import :: REAL64, huberfit_function
       class(huberfit_function) :: this
       real(REAL64), intent(in) :: x
       real(REAL64), dimension(:), intent(in) :: p
       real(REAL64) :: fval
     end function fvalue

     subroutine fgradient(this,x,p,fgrad)
       import :: REAL64, huberfit_function
       class(huberfit_function) :: this
       real(REAL64), intent(in) :: x
       real(REAL64), dimension(:), intent(in) :: p
       real(REAL64), dimension(:), intent(out) :: fgrad
     end subroutine fgradient

  end interface

  interface

     module subroutine huberfit(fcn,x,y,par,parerr,stdsig,scale,H,gtol, &
          maxiter,verbose,report_file,iflag)
       class(huberfit_function), intent(in), target :: fcn
       real(REAL64), dimension(:), intent(in), target :: x,y
       real(REAL64), dimension(:), intent(in out) :: par
       real(REAL64), dimension(:), intent(out) :: parerr
       real(REAL64), intent(in out) :: scale
       real(REAL64), intent(out) :: stdsig
       real(REAL64), dimension(:,:), intent(out) :: H
       real(REAL64), intent(in) :: gtol
       integer, intent(in) :: maxiter
       logical, intent(in) :: verbose
       character(len=*), intent(in) :: report_file
       integer, intent(out) :: iflag
     end subroutine huberfit

  end interface

  public :: huberfit

contains

  subroutine fhessian(this,x,p,fhess)
    class(huberfit_function) :: this
    real(REAL64), intent(in) :: x
    real(REAL64), dimension(:), intent(in) :: p
    real(REAL64), dimension(:,:), intent(out) :: fhess

    fhess = 0

  end subroutine fhessian

end module huber_curvefit


submodule(huber_curvefit) huberfun_functions

  use rfun
  use bfgs_method
  use iso_fortran_env
  implicit none

  type, abstract, extends(bfgs_function) :: bfgs_robfun
     class(huberfit_function), pointer :: fcn
     real(REAL64), dimension(:), pointer :: x,y
   contains
     procedure :: maxalf => maxalpha
  end type bfgs_robfun

  type, extends(bfgs_robfun) :: huber_bfgs
     real(REAL64) :: scale = 1
   contains
     procedure :: val => rvalue
     procedure :: grad => rgradient
     procedure :: hess => rhessian
     procedure :: fungrad => rfungrad
     procedure :: break => breaker
  end type huber_bfgs

  type, extends(bfgs_robfun) :: lsq_bfgs
   contains
     procedure :: val => lvalue
     procedure :: grad => lgradient
     procedure :: hess => lhessian
     procedure :: fungrad => lfungrad
  end type lsq_bfgs

  type, extends(bfgs_robfun) :: wsq_bfgs
     class(huber_bfgs), pointer :: fun
   contains
     procedure :: val => wvalue
     procedure :: grad => wgradient
     procedure :: hess => whessian
     procedure :: fungrad => wfungrad
  end type wsq_bfgs

contains

  function rvalue(this,x) result(fval)

    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64) :: fval
    real(REAL64), dimension(size(x)-1) :: p
    real(REAL64) :: s, r, ravg
    integer :: i,n

    s = x(0)
    p = x(1:)

    associate( fcn => this%fcn, x => this%x, y => this%y )
      n = size(x)

      ravg = 0
      do i = 1, n
         r = y(i) - fcn%val(x(i),p)
         ravg = ravg + ipshuber(r/s)
      end do

    end associate

    fval = ravg / n + log(s)

  end function rvalue

  subroutine rgradient(this,x,fgrad)

    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), dimension(0:), intent(out) :: fgrad
    real(REAL64), dimension(size(x)-1) :: p, fg
    real(REAL64) :: s, r, psi
    integer :: i,n

    p = x(1:)
    s = x(0)

    associate( fcn => this%fcn, x => this%x, y => this%y )
      n = size(x)
      fgrad = 0
      do i = 1, n
         r = y(i) - fcn%val(x(i),p)
         call fcn%grad(x(i),p,fg)
         psi = pshuber(r/s)
         fgrad(1:) = fgrad(1:) + psi*fg
         fgrad(0) = fgrad(0) + psi*r
      end do
    end associate

    fgrad(1:) = -fgrad(1:) / (s*n)
    fgrad(0) = -fgrad(0)/(s**2*n) + 1/s

  end subroutine rgradient

  subroutine rfungrad(this,x,fval,fgrad)

    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), intent(out) :: fval
    real(REAL64), dimension(0:), intent(out) :: fgrad
    real(REAL64), dimension(size(x)-1) :: p, fg
    real(REAL64) :: s, r, rs, ravg, psi
    integer :: i,n

    p = x(1:)
    s = x(0)
    this % scale = s

    associate( fcn => this%fcn, x => this%x, y => this%y )
      n = size(x)
      ravg = 0
      fgrad = 0
      do i = 1, n
         r = y(i) - fcn%val(x(i),p)
         call fcn%grad(x(i),p,fg)
         rs = r / s
         ravg = ravg + ipshuber(rs)
         psi = pshuber(rs)
         fgrad(1:) = fgrad(1:) + psi*fg
         fgrad(0) = fgrad(0) + psi*rs
      end do
    end associate

    fval = ravg / n + log(s)
    fgrad = -fgrad / (s*n)
    fgrad(0) = fgrad(0) + 1/s

  end subroutine rfungrad


  subroutine rhessian(this,x,fhess)
    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), dimension(0:,0:), intent(out) :: fhess
    real(REAL64), dimension(size(x)-1,size(x)-1) :: hess
    real(REAL64), dimension(size(x)-1) :: p, fg
    real(REAL64) :: s, r, psi, dpsi
    integer :: i,n,m,k,l

    m = size(x) - 1
    p = x(1:)
    s = x(0)

    associate( fcn => this%fcn, x => this%x, y => this%y )
      n = size(x)
      fhess = 0
      do i = 1, n
         r = y(i) - fcn%val(x(i),p)
         call fcn%grad(x(i),p,fg)
         call this%fcn%hess(x(i),p,hess)
         psi = pshuber(r/s)
         dpsi = dpshuber(r/s)
         do l = 1, m
            do k = 1, m
               fhess(k,l) = fhess(k,l) + dpsi/s * fg(k)*fg(l) - psi*hess(k,l)
            end do
         end do
         fhess(1:m,0) = fhess(1:m,0) + dpsi*(r/s)*fg(1:m) + psi*fg(1:m)
         fhess(0,0)   = fhess(0,0) + dpsi*(r/s)**2 + 2*psi*(r/s) - 1
      end do
    end associate

    fhess(1:m,1:m) = fhess(1:m,1:m) / (n*s)
    fhess(1:m,0) = fhess(1:m,0) / (n*s**2)
    fhess(0,1:m) = fhess(1:m,0)
    fhess(0,0) = fhess(0,0) / (n*s**2)

  end subroutine rhessian

  function wvalue(this,x) result(fval)

    class(wsq_bfgs) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64) :: fval
    real(REAL64), dimension(0:size(x)) :: p

    p(1:) = x
    p(0) = 1
    fval = this % fun % val(p)

  end function wvalue

  subroutine wgradient(this,x,fgrad)

    class(wsq_bfgs) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), intent(out) :: fgrad
    real(REAL64), dimension(0:size(x)) :: p

    p(1:) = x
    p(0) = 1
    call this % fun % grad(p,fgrad)

  end subroutine wgradient

  subroutine wfungrad(this,x,fval,fgrad)

    class(wsq_bfgs) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), intent(out) :: fval
    real(REAL64), dimension(:), intent(out) :: fgrad
    real(REAL64), dimension(0:size(x)) :: p, fg

    p(1:) = x
    p(0) = 1
    call this % fun % fungrad(p,fval,fg)
    fgrad = fg(1:)

  end subroutine wfungrad

  subroutine whessian(this,x,fhess)
    class(wsq_bfgs) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:,:), intent(out) :: fhess
    real(REAL64), dimension(0:size(x),0:size(x)) :: hess
    real(REAL64), dimension(0:size(x)) :: p

    p(1:) = x
    p(0) = 1
    call this % fun % hess(p,hess)
    fhess = hess(1:,1:)

  end subroutine whessian

  function lvalue(this,x) result(fval)

    class(lsq_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64) :: fval
    real(REAL64), dimension(size(x)-1) :: p
    real(REAL64) :: s, r, ravg, fval
    integer :: i,n

    p = x(1:)
    s = x(0)

    associate( fcn => this%fcn, x => this%x, y => this%y )
      n = size(x)

      ravg = 0
      do i = 1, n
         r = y(i) - fcn%val(x(i),p)
         ravg = ravg + (r/s)**2/2
      end do

    end associate

    fval = ravg / n + log(s)

  end function lvalue

  subroutine lgradient(this,x,fgrad)

    class(lsq_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), dimension(0:), intent(out) :: fgrad
    real(REAL64), dimension(size(x)-1) :: p, fg
    real(REAL64) :: s, r, psi
    integer :: i,n

    p = x(1:)
    s = x(0)

    associate( fcn => this%fcn, x => this%x, y => this%y )
      n = size(x)
      fgrad = 0
      do i = 1, n
         r = y(i) - fcn%val(x(i),p)
         call fcn%grad(x(i),p,fg)
         psi = r / s
         fgrad(1:) = fgrad(1:) + psi*fg
         fgrad(0) = fgrad(0) + psi*r
      end do
    end associate

    fgrad(1:) = -fgrad(1:) / (s*n)
    fgrad(0) = -fgrad(0)/(s**2*n) + 1/s

  end subroutine lgradient

  subroutine lfungrad(this,x,fval,fgrad)

    class(lsq_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), intent(out) :: fval
    real(REAL64), dimension(0:), intent(out) :: fgrad
    real(REAL64), dimension(size(x)-1) :: p, fg
    real(REAL64) :: s, r, rs, rs2, s2
    integer :: i,n

    p = x(1:)
    s = x(0)

    associate( fcn => this%fcn, x => this%x, y => this%y )
      n = size(x)
      s2 = 0
      fgrad = 0
      do i = 1, n
         r = y(i) - fcn%val(x(i),p)
         rs = r / s
         rs2 = rs**2
         call fcn%grad(x(i),p,fg)
         s2 = s2 + rs2/2
         fgrad(1:) = fgrad(1:) + rs*fg
         fgrad(0) = fgrad(0) + rs2
      end do
    end associate

    fval = s2 / n + log(s)
    fgrad = -fgrad / (s*n)
    fgrad(0) = fgrad(0) + 1 / s

  end subroutine lfungrad


  subroutine lhessian(this,x,fhess)
    class(lsq_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), dimension(0:,0:), intent(out) :: fhess
    real(REAL64), dimension(size(x)-1,size(x)-1) :: hess
    real(REAL64), dimension(size(x)-1) :: p, fg
    real(REAL64) :: s, r, psi, dpsi, fval
    integer :: i,n,m,k,l

    m = size(x) - 1
    p = x(1:)
    s = x(0)

    associate( fcn => this%fcn, x => this%x, y => this%y )
      n = size(x)
      fhess = 0
      do i = 1, n
         fval = fcn%val(x(i),p)
         call fcn%grad(x(i),p,fg)
         call this%fcn%hess(x(i),p,hess)
         r = y(i) - fval
         psi = r / s
         dpsi = 1
         do l = 1, m
            do k = 1, m
               fhess(k,l) = fhess(k,l) + dpsi/s * fg(k)*fg(l) - psi*hess(k,l)
            end do
         end do
         fhess(1:m,0) = fhess(1:m,0) + dpsi*(r/s)*fg(1:m) + psi*fg(1:m)
         fhess(0,0)   = fhess(0,0) + dpsi*(r/s)**2 + 2*psi*(r/s) - 1
      end do
    end associate

    fhess(1:m,1:m) = fhess(1:m,1:m) / (n*s)
    fhess(1:m,0) = fhess(1:m,0) / (n*s**2)
    fhess(0,1:m) = fhess(1:m,0)
    fhess(0,0) = fhess(0,0) / (n*s**2)

  end subroutine lhessian


  real(REAL64) function maxalpha(this,x,p)

    class(bfgs_robfun) :: this
    real(REAL64), dimension(0:), intent(in) :: x, p
    real(REAL64) :: s, ps

    s = x(0)
    ps = p(0)
    if( s + ps <= 0 ) then
       maxalpha = -0.9*(s / ps)
    else
       maxalpha = 1
    end if

  end function maxalpha

  logical function breaker(this)

    class(huber_bfgs) :: this
    real(REAL64) :: s

    s = this % scale
    breaker = s**2 < 1e3*epsilon(s)

  end function breaker

end submodule huberfun_functions


submodule(huber_curvefit : huberfun_functions) huberfit_huberfit

  logical :: debug = .true.

contains

  module subroutine huberfit(fcn,x,y,par,parerr,stdsig,scale,H,gtol,maxiter, &
       verbose,report_file,iflag)

    use bfgs_method
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN

    class(huberfit_function), intent(in), target :: fcn
    real(REAL64), dimension(:), intent(in), target :: x,y
    real(REAL64), dimension(:), intent(in out) :: par
    real(REAL64), dimension(:), intent(out) :: parerr
    real(REAL64), intent(in out) :: scale
    real(REAL64), intent(out) :: stdsig
    real(REAL64), dimension(:,:), intent(out) :: H
    real(REAL64), intent(in) :: gtol
    integer, intent(in) :: maxiter
    logical, intent(in) :: verbose
    character(len=*), intent(in) :: report_file
    integer, intent(out) :: iflag

    character(len=*), parameter :: robfun = 'huber'
    class(bfgs_robfun), allocatable :: fun, wfun
    real(REAL64), dimension(0:size(par)) :: p
    real(REAL64), dimension(0:size(par),0:size(par)) :: HH
    real(REAL64) :: s, psi2, dpsi, r,  rs
    integer :: i,m,n

    iflag = 6
    n = size(x)
    m = size(par)

    if( .not. (n > m) ) then
       write(error_unit,*) 'No data, no love in rfit().'
       return
    end if

    if( robfun == 'huber' ) then
       fun = huber_bfgs(fcn,x,y)
    else if( robfun == 'lsq' ) then
       fun = lsq_bfgs(fcn,x,y) ! for testing purposes
    else
       stop "Error: robfun parameter misspelled, try: `huber' or `lsq'."
    end if

    p(0) = scale
    p(1:) = par

    call bfgsmin(fun,p,HH,gtol,maxiter,verbose,report_file,iflag)
    par = p(1:)
    s = p(0)
    scale = s
    H = HH(1:,1:)

    if( verbose ) then
       write(error_unit,*) 'scale: ',s,'    pars: ', par,'  iflag: ',iflag
    end if

    ! the scale goes to zero, we are facing data without a scatter
    if( iflag == 6 .and. s**2 < 1e3*epsilon(s) ) then
       wfun = wsq_bfgs(fcn,x,y,fun)
       call bfgsmin(wfun,par,H,gtol,maxiter,verbose,report_file,iflag)
       s = 1
       if( verbose ) then
          write(error_unit,*) 'scale: 1','    pars: ', par,'  iflag: ',iflag
       end if
    end if

    ! stdsig
    psi2 = 0
    dpsi = 0
    do i = 1, n
       r = y(i) - fun%fcn%val(x(i),p(1:m))
       rs = r / s
       if( robfun == 'huber' ) then
          psi2 = psi2 + huber(rs)**2
          dpsi = dpsi + dhuber(rs)
       else if( robfun == 'lsq' ) then
          psi2 = psi2 + rs**2
          dpsi = dpsi + 1
       end if
    end do
    if( psi2 > 0 .and. dpsi > 0 ) then
       stdsig = s*sqrt(psi2*n)/dpsi
    else if( dpsi > 0 ) then
       stdsig = 0
    else
       stdsig = IEEE_VALUE(stdsig, IEEE_QUIET_NAN)
    end if

    do i = 1, m
       if( H(i,i) > 0 ) then
          parerr(i) = sqrt(H(i,i)/(n - m))
       else
          parerr(i) = IEEE_VALUE(parerr(i), IEEE_QUIET_NAN)
       end if
    end do

    if( debug ) then
       block

         use QRHouse

         real(REAL64), dimension(0:size(par),0:size(par)) :: eye, hess, Hinv
         real(REAL64) :: s0, s2, r, fval
         integer :: i

         call qrinv(HH,Hinv)

!         eye = matmul(HH,Hinv)
!         do i = 0, size(par)
!            write(error_unit,'(3f20.15)') eye(i,:)
!         end do

         s2 = 0
         do i = 1, n
            r = y(i) - fun%fcn%val(x(i),p(1:m))
            s2 = s2 + r**2
         end do

         fval = fun%val(p)
         s0 = 2*n*(fval - log(s))*s**2
         write(error_unit,*) 's2 = ',s2,fval,n*(fval-log(s))*s**2*2,s0

         call fun%hess(p,hess)
         do i = 0,2
            write(error_unit,'(a,3f25.15)') 'hess: ',hess(i,:)*(N*s**2)
         end do
         do i = 0,2
            write(error_unit,'(a,3f25.15)') 'Hinv: ',Hinv(i,:)*(N*s**2)
         end do

         s0 = 2*(fval - log(s))
         write(error_unit,*) s0
         write(error_unit,*) H(1:2,1:2)
         write(error_unit,*) 'S0 rob.: ', N*s**2
       end block

    end if

  end subroutine huberfit

end submodule huberfit_huberfit
