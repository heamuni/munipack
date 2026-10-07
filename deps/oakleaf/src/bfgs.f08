!
!  The BFGS method
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
!
!  This module implements the optimisation procedure based on quasi-Newton
!  approximation of the inverse Hessian by the BFGS method combined with the
!  linear search for the step length. The BSFG approximation gives strictly
!  positive definite matrixes representing the inverse Hessian H = B^-1.
!  The step is computed by matrix multiplication and the final inverse
!  can be used to estimate of the parameters errors without a matrix inversion.
!  The linear search for the Wolfe's conditions implements the enhanced
!  version of Alg. 3.5 and Alg. 3.6 preventing numerical uncertainities.
!
!  Based on Nocedal & Wright, Ch. 3. and 6.
!
!  This algorithm is considered as the robust counterpart of Levendberg-
!  Marquart algorithm implemented by Minpack. This implementation release
!  the restriction on the minimized function (LM is based on square roots),
!  improves procedures for the linear search parameter, and does not require
!  the labourly derived analytical Hessian. It is suitable for fully
!  determined systems without constrains. The minimised function
!  should be convex over the whole domain. Generally, the inverse Hessian
!  approximation is poor (~10%), but suitable for estimation of the parameters
!  statistical errors.
!
!  There is maxalpha() function which can be overwritten to control
!  the line search --- it prevents to calculate functions somewhere out of
!  a feasibility region; via that, constrains can be handled
!  if the minimum lie deep inside the feasibility region.
!
!  Computations can be interrupted by user via break() function.
!
!  There are two parameters to control convergence. The crutial is gtol,
!  a tolerance for the Euclidean norm of the gradient --- we are finding
!  the minimum. Its reccomended value is gtol ~= 1e3* npar * machine epsilon;
!  the value usually gives acceptable approximation of the inverse Hessian.
!  There is also the internal test for values |f_i - f_{i+1} | < mach. eps.
!  Both test produces valid results with iflag=10 and 11 respectively.
!  The maxiter is the last resort if no convergence is reached. Minpack
!  reccomends 200 iterations per parameter, thousands are usually enought.
!
!  Usage:
!   * specify function values and gradients in the abstract bfgs_function
!   * run bfgsmin()
!   * to improve efectivity, you can rewrite fungrad()
!
!  Examples are in test/bfgs/rosenbrock.f08 and test/bfgs/line.f08


module bfgs_method

  use iso_fortran_env
  implicit none
  private

  logical, parameter :: debug = .false.
  logical, parameter :: wolfe_search_activated = .false.
  real(REAL64), parameter :: macheps = epsilon(1.0_REAL64)

  ! parameters of Wolfe's conditions
  real, parameter :: c1 = 1e-4
  real, parameter :: c2 = 0.9

  type, abstract, public :: bfgs_function
   contains
     procedure(fvalue), deferred :: val
     procedure(fgradient), deferred :: grad
     procedure :: fungrad => fungradient
     procedure :: hess => fhessian
     procedure :: maxalf => maxalpha
     procedure :: break => breaker
  end type bfgs_function

  abstract interface

     real(REAL64) function fvalue(this,x)
       import :: REAL64, bfgs_function
       class(bfgs_function) :: this
       real(REAL64), dimension(:), intent(in) :: x
     end function fvalue

     subroutine fgradient(this,x,fgrad)
       import :: REAL64, bfgs_function
       class(bfgs_function) :: this
       real(REAL64), dimension(:), intent(in) :: x
       real(REAL64), dimension(:), intent(out) :: fgrad
     end subroutine fgradient

  end interface

  public :: bfgsmin

contains

  subroutine bfgsmin(fcn,x,H,gtol,maxiter,verbose,report_file,iflag)

    class(bfgs_function), intent(in) :: fcn
    real(REAL64), dimension(:), intent(in out) :: x
    real(REAL64), dimension(:,:), intent(out) :: H ! inverse quasi-Hessian
    real(REAL64), intent(in) :: gtol
    integer, intent(in) :: maxiter
    logical, intent(in) :: verbose
    character(len=*), intent(in) :: report_file
    integer, intent(out) :: iflag

    real(REAL64), parameter :: glim = sqrt(epsilon(gtol))
    real(REAL64) :: fval, fv, alf, maxalf, gnorm
    real(REAL64), dimension(size(x)) :: fgrad, p, fg, y, s
    integer :: n, info, iter, iostat, runit, nfev, ngrad
    logical :: curvature, report
    character(len=88) :: errmsg, xfmt,vfmt

    iflag = 6
    n = size(x)
    write(xfmt,'(a,i0,a)') '(a,',n**2,'(g0,1x))'
    write(vfmt,'(a,i0,a)') '(i5,a,',n,'(f9.3,1x),a,2es10.3,a,f7.3,es11.3,i4,l2)'
    call fcn%fungrad(x,fval,fgrad)
    H = eye(n)
    nfev = 1
    ngrad = 1
    p = 0

    report = report_file /= ''
    if( report ) then
       open(newunit=runit,file=report_file,iostat=iostat,iomsg=errmsg)
       if( iostat /= 0 ) then
          if( verbose ) write(error_unit,*) &
               'Error in bfgsmin(): ',trim(errmsg),': `',trim(report_file),"'."
          iflag = 6
          close(runit)
          return
       end if
       write(runit,'(a)') '# bfgsmin()'
       write(runit,'(a,i0)') '# npar = ',n
       write(runit,xfmt) '# x0 = ',x
       write(runit,'(a)') '# iter x fval alpha p fgrad H flag'
       write(runit,*) 0,x,fval,1,p,fgrad,H,0
    end if

    if( verbose ) then
       write(error_unit,*) 'iter: pars ... > fval, |grad| | alpha, fval-fv,  info'
       write(error_unit,vfmt) 0,': ',x," > ",fval,norm2(fgrad)
    end if

    curvature = .true.
    do iter = 1, maxiter

       ! the step
       p = -matmul(H,fgrad)
       if( debug ) write(error_unit,*) 'p=',real(p),', (matmul)'

       ! the linear search
       if( wolfe_search_activated ) then
          maxalf = fcn%maxalf(x,p)
          call wolfe_search(fcn,x,p,fval,fgrad,maxalf,alf,nfev,ngrad,info)
          if( debug ) write(error_unit,*) 'wolfe search =',alf,info
       else ! backtrace search prefered
          maxalf = min(1.0_REAL64,fcn%maxalf(x,p))
          call back_search(fcn,x,p,fval,fgrad,maxalf,alf,nfev,info)
          if( debug ) write(error_unit,*) 'back search =',alf,info
       end if
       if( info >= 60 ) then
          if( debug ) write(error_unit,*) &
               'Warning in bfgsmin: the linear search unsuccessful.'
       end if

       ! reports
       if( report ) write(runit,*) iter,x,fval,alf,alf*p,fgrad,H,info
       s = alf * p
       x = x + s
       call fcn%fungrad(x,fv,fg)
       nfev = nfev + 1
       ngrad = ngrad + 1

       if( verbose ) write(error_unit,vfmt) &
            iter,': ',x," > ",fv,norm2(fg)," |",alf,fval-fv,info,curvature

       ! interrupt on request
       if( fcn % break() ) then
          iflag = 9
          exit
       end if

       ! convergence check
       if( info < 60 ) then
          gnorm = norm2(fg)
          if( gnorm < gtol ) then
             iflag = 0
             exit
          else if( abs(fval - fv) <= macheps .and. gnorm < gtol ) then
             iflag = 1
             exit
          end if
       end if

       ! update
       y = fg - fgrad
       fgrad = fg
       fval = fv

       ! the first update of inverse Hessian by approximate value
       ! of the largest eigenvalue, NW Ch.6
       if( iter == 1 ) then
          H = dot_product(s,y) / dot_product(y,y) * H
          if( debug ) then
             block
               integer :: i
               character(len=88) :: fmt
               write(fmt,'(a,i0,a)') '(a,',n,'(f13.5,2x))'
               do i = 1, n
                  write(error_unit,fmt) 'init H:',H(i,:)
               end do
             end block
          end if
       end if

       ! BFGS update
       call bfgs_update(s,y,H,curvature)
       if( .not. curvature ) then
          ! The curvature check is ignored.
          ! For Wolfe search, it commonly fails if the line search fails,
          ! and it fails if the parameter range is too tight.
          ! For both Wolfe and back search fails close to the minimum
          ! due rounding as y and Bs are small, the replacing by the
          ! analytical Hessian usually does not improve accuracy.
          if( debug ) write(error_unit,*) 'Error in bfgsmin() on curvature.'
       end if

    end do

    if( report ) then
       write(runit,xfmt) '# solution = ',x
       write(runit,'(a,g0)') '# iflag = ',iflag
       write(runit,'(a,g0)') '# niter = ',iter
       write(runit,'(a,g0)') '# nfev = ',nfev
       write(runit,'(a,g0)') '# ngrad = ',ngrad
       write(runit,xfmt) '# inv-hess approx = ',H
       block
         real(REAL64), dimension(n,n) :: hess
         call fcn%hess(x,hess)
         write(runit,xfmt) '# hess = ',hess
       end block
       close(runit)
    end if

    if( debug ) then
       block
         use QRHouse
         integer :: i
         character(len=88) :: fmt
         real(REAL64), dimension(n,n) :: hess, B
         write(fmt,'(a,i0,a)') '(a,',n,'(f13.5,2x))'
         call qrinv(H,B)
         do i = lbound(H,1), ubound(H,1)
            write(error_unit,fmt) '   H: ',H(i,:)
         end do
         do i = lbound(B,1), ubound(B,1)
            write(error_unit,fmt) '   B: ',B(i,:)
         end do
         call fcn%hess(x,hess)
         do i = 1, n
            write(error_unit,fmt) 'hess: ',hess(i,:)
         end do
         write(error_unit,*) '|B - hess|/|hess|: ',norm2(B - hess)/norm2(hess)
       end block
    end if

  end subroutine bfgsmin


  pure subroutine bfgs_update(s,y,H,curvature)

    real(REAL64), dimension(:), intent(in) :: s,y
    real(REAL64), dimension(:,:), intent(in out) :: H
    logical, intent(out) :: curvature

    real(REAL64), dimension(size(y)) :: Hy
    real(REAL64) :: ys, rho, yHy, w
    integer :: i,j,n

    n = size(s)
    Hy = matmul(H,y)
    yHy = dot_product(y,Hy)
    ys = dot_product(y,s)
    curvature = ys > 0
    if( .not. curvature ) return
    rho = 1 / ys

    ! the formula by https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm

    w = rho*(ys + yHy)
    forall( i = 1:n, j = 1:n )
       H(i,j) = H(i,j) + rho*(w*s(i)*s(j) - (Hy(j)*s(i) + s(j)*Hy(i)))
    end forall

  end subroutine bfgs_update

  pure function eye(n) result(B)

    integer, intent(in) :: n
    real(REAL64), dimension(:,:), allocatable :: B
    integer :: i

    allocate(B(n,n))
    B = 0
    forall( i = 1:n ) B(i,i) = 1

  end function eye

  subroutine back_search(fcn,x,p,fval,fgrad,maxalf,alf,nfev,info)

    real(REAL64), parameter :: gold = (1 + sqrt(5.0)) / 2

    class(bfgs_function), intent(in), target :: fcn
    real(REAL64), dimension(:), intent(in) :: x,p,fgrad
    real(REAL64), intent(in) :: fval, maxalf
    real(REAL64), intent(out) :: alf
    integer, intent(in out) :: nfev
    integer, intent(out) :: info
    real(REAL64) :: dp, df, fv
    integer :: iter

    info = 60
    alf = min(1.0_REAL64, maxalf)
    df = dot_product(fgrad,p)

    do iter = 1, 5*precision(alf)
       fv = fcn % val(x  + alf*p)
       nfev = nfev + 1
       dp = c1*alf*df
!       write(*,*) iter, alf,fv-fval,fval + c1*alf*df,df,epsilon(df),c1*alf*df
       if( fv <= fval + dp ) then
          info = 11
          exit
       else if( abs(dp) < epsilon(df) ) then
          info = 12
          exit
       end if
       alf = alf / gold
    end do

  end subroutine back_search


  subroutine wolfe_search(fcn,x,p,fval,fgrad,maxalf,alf,nfev,ngrad,info)

    real(REAL64), parameter :: ftol = 10*macheps
    real(REAL64), parameter :: gold = (1 + sqrt(5.0)) / 2

    class(bfgs_function), intent(in), target :: fcn
    real(REAL64), dimension(:), intent(in) :: x,p,fgrad
    real(REAL64), intent(in) :: fval, maxalf
    real(REAL64), intent(out) :: alf
    integer, intent(in out) :: nfev, ngrad
    integer, intent(out) :: info
    real(REAL64) :: dp, df, fv, fv1, alf0, alf1
    real(REAL64), dimension(size(fgrad)) :: fg!, xp
    integer :: iter

    info = 60
    alf = min(1.0_REAL64, maxalf)
    fv1 = huge(fv1)

    dp = dot_product(fgrad,p)
    alf0 = 0
    do iter = 1,3*precision(alf)
       call fcn%fungrad(x + alf*p,fv,fg)
       nfev = nfev + 1
       ngrad = ngrad + 1

       if( debug ) write(error_unit,'(a,i3,2f8.3,4es13.3)') &
            'wolfe ',iter, alf,maxalf,fv,fval + c1*alf*dp,dp,dot_product(fg,p)

       if( fv > fval + c1*alf*dp .or. (fv >= fv1 .and. iter > 1) ) then
          alf1 = alf
          call zoom(fcn,x,p,fval,fgrad,dp,alf0,alf1,alf,nfev,ngrad,info)
          info = info + 1
          exit
       end if

       df = dot_product(fg,p)

       if( abs(df) <= -c2*dp .or. (abs(df) < ftol .and. abs(dp) < ftol) ) then
          info = 12
          exit
       end if

       if( df >= 0 ) then
          alf1 = alf
          call zoom(fcn,x,p,fval,fgrad,dp,alf1,alf0,alf,nfev,ngrad,info)
          info = info + 3
          exit
       end if

       alf0 = alf
       alf = alf * gold
       if( alf > maxalf ) then
          alf = min(1.0_REAL64, maxalf / 2)
          exit
       end if
       fv1 = fv
    end do

  end subroutine wolfe_search

  subroutine zoom(fcn,x,p,fval,fgrad,dp,alflo,alfhi,alf,nfev,ngrad,info)

    real(REAL64), parameter :: ftol = 10*macheps
    real(REAL64), parameter :: gtol = sqrt(macheps)

    class(bfgs_function), intent(in), target :: fcn
    real(REAL64), dimension(:), intent(in) :: x,p,fgrad
    real(REAL64), intent(in) :: fval,dp,alflo,alfhi
    real(REAL64), intent(out) :: alf
    integer, intent(in out) :: nfev, ngrad
    integer, intent(in out) :: info
    real(REAL64) :: fv, alfmin, alfmax, df, fmin
    real(REAL64), dimension(size(fgrad)) :: fg, xp
    integer :: iter

    info = 60
    alfmin = alflo
    alfmax = alfhi
    fmin = fcn%val(x + alfmin*p)
    nfev = nfev + 1

    if( debug ) write(error_unit,*) 'zoom',alfmin,alfmax

    do iter = 1, 3*precision(alf)

       alf = (alfmin + alfmax) / 2
       xp = x + alf*p
       fv = fcn%val(xp)
       nfev = nfev + 1

       if( debug ) then
          call fcn%grad(xp,fg)
          df = dot_product(fg,p)
          write(error_unit,'(a,i3,3f10.5,4l3,1x,4es12.2)') &
               ' zoom: ',iter,alf,alfmin,alfmax, &
               fv <= fval + c1*alf*dp + ftol, &
               fv <= fmin + ftol, &
               abs(df) <= -c2*dp,df*(alfmax - alfmin) >= 0, &
               fv-(fval + c1*alf*dp), fv - fmin, dp, df
       end if

       if( .not. (fv <= fval + c1*alf*dp + ftol .and. fv <= fmin + ftol) ) then
          alfmax = alf
       else
          call fcn%grad(xp,fg)
          ngrad = ngrad + 1
          df = dot_product(fg,p)
          if( abs(df) <= -c2*dp + gtol .or. (abs(df) < gtol .and. abs(dp) < gtol) ) then
             info = 10
             return
          end if
          if( df*(alfmax - alfmin) >= 0 ) then
             alfmax = alfmin
          end if
          alfmin = alf
          fmin = fv
       end if
    end do

  end subroutine zoom

  real(REAL64) function maxalpha(this,x,p)

    ! This gives maximal value of line search.
    ! The principial contribution is in case of
    ! robust fitting, see src/rfit.f08.

    class(bfgs_function) :: this
    real(REAL64), dimension(:), intent(in) :: x,p

    if( wolfe_search_activated ) then
       maxalpha = 999
    else
       maxalpha = 1
    end if

  end function maxalpha

  logical function breaker(this)

    ! This gives maximal value of line search.
    ! The principial contribution is in case of
    ! robust fitting, see src/rfit.f08.

    class(bfgs_function) :: this

    breaker = .false.

  end function breaker

  subroutine fhessian(this,x,fhess)

    ! The numerical approximation of Hessian.
    ! As bfgsmin() is based on quasi-Newton approximation,
    ! this routine is, in fact, only for testing, and similar, purposes.

    class(bfgs_function) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:,:), intent(out) :: fhess

    real(REAL64), dimension(size(x)) :: d,fv1,fv2
    real(REAL64) :: h
    integer :: i

    d = 0
    do i = lbound(x,1), ubound(x,1)
       h = sqrt(macheps)*(1 + abs(x(i)))
       d(i) = h
       call this%grad(x+d,fv1)
       call this%grad(x-d,fv2)
       fhess(i:,i) = (fv1(i:) - fv2(i:))/(2*h)
       fhess(i,i:) = fhess(i:,i)
       d(i) = 0
    end do

  end subroutine fhessian

  subroutine fungradient(this,x,fval,fgrad)
    class(bfgs_function) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), intent(out) :: fval
    real(REAL64), dimension(:), intent(out) :: fgrad

    fval = this % val(x)
    call this % grad(x,fgrad)

  end subroutine fungradient


! The commented-out routines are based on the Hessian,
! while the suroutine adapts the inverse Hessian and it quasi-approximation.
!!$  pure subroutine chstep(B,g,p,posdef)
!!$
!!$    use cholesky
!!$
!!$    real(REAL64), dimension(:,:), intent(in) :: B
!!$    real(REAL64), dimension(:), intent(in) :: g
!!$    real(REAL64), dimension(:), intent(out) :: p
!!$    logical, intent(out) :: posdef
!!$    real(REAL64), dimension(size(B,1),size(B,2)) :: R
!!$    real(REAL64), dimension(size(p)) :: q
!!$
!!$    call chol_rrfac(B,R,posdef)
!!$    if( posdef ) then
!!$       call chol_rrsolve(R,-g,p)
!!$       call chol_rrtsolve(R,p,q)
!!$    end if
!!$
!!$  end subroutine chstep
!!$
!!$  pure subroutine bfgs_update(s,y,B,curvature)
!!$
!!$    real(REAL64), dimension(:), intent(in) :: s,y
!!$    real(REAL64), dimension(:,:), intent(in out) :: B
!!$    logical, intent(out) :: curvature
!!$
!!$    real(REAL64), dimension(size(s)) :: Bs
!!$    real(REAL64) :: sBs, ys
!!$    integer :: i,j,n
!!$
!!$    n = size(s)
!!$    Bs = matmul(B,s)
!!$    sBs = dot_product(s,Bs)
!!$    ys = dot_product(y,s)
!!$    curvature = ys > 0
!!$    if( .not. curvature ) return
!!$    forall( i = 1:n, j = 1:n )
!!$       B(i,j) = B(i,j) - Bs(i)*Bs(j) / sBs + y(i)*y(j) / ys
!!$    end forall
!!$
!!$  end subroutine bfgs_update

end module bfgs_method
