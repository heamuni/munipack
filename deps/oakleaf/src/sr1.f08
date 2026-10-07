!
!  Trust-region method on base of SR1 and dgqt()
!
!  Copyright © 2021-2022 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of OakLeaf.
!
!  OakLeaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  OakLeaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with OakLeaf.  If not, see <http://www.gnu.org/licenses/>.
!
!

! This is the implementation of Symmetric-rank-1 method on base of
! the algorithm 6.2 (SR1 Trust-Region Method) published in the second
! edition of Numerical Optimization by Jorge Nocedal and Stephen J. Wright.
! Trust region step is computed by ideas of work More & Sorenson (1981).
! All details can be found in the descriptions of specfic routines.
!

! sr1 should be called as:
!
! call sr1(fcn,x,hess,delta,maxdev,maxfev,verbose,flag,msg)
!
! on input:
! fcn ... abstract class representing a function to minimise
! x ... initial estimate of minimum
! delta (optional) .. initial radius of the trust region
! maxdev (optional).. maximum deviation of results
! maxfev (optional).. maximum number of function evaluations
! verbose(optional).. causess to print an additional info about the solution path
!
! on output:
! x .. final estimate of the fcn minimum
! hess(optional).. final estimate of the Hessian matrix
! flag(optional).. info about the solution
! mgs(optional) .. description of the solution for humans
!
! fcn is of type sr1_function:
! the pure abstract class should be implemented in calling routine. All methods
! should be defined (but a dummy definition of print() is acceptable).
! fval() computes a value of the function; fgrad() computes the vector
! of the first partial derivatives. Hessian matrix (second derivatives)
! is computed by the SR1 quasi-Newton updates.
!
! The fcn contains the logical variable: break. The variable should be
! set on .true., only if the procedure should be imediatelly interupted.
!
! The easy way for implementing of the function is to develop a module
! containing the implementation, and possible data, or something else.
! test/sr1/*.f08 contains examples.
!
! sr1_function%fval(this,x)
!
! on input:
!  this .. the membership
!  x .. vector of the independer coordinates
! retuns derived function the point by x
!
! sr1_function%fgrad(this,x,fgrad)
!
! on input:
!  this .. the membership
!  x .. vector of the independer coordinates
!
! on output:
!  fgrad .. gradient of the function in point by x
!
!
! sr1_function%print(this,iter,x,fval,fgrad,Hess,step,delta,rho)
!
! on input, there're the current values of
!
!  this .. the membership
!  iter .. iteration counter
!  x .. point
!  fval .. function value
!  fgrad .. gradient
!  Hess .. Hessian
!  step .. step
!  delta .. trust radius
!




module symmetric_rank_1

  use iso_fortran_env
  implicit none
  private

  integer, parameter :: kind = REAL64

  type, abstract, public :: sr1_function
     logical :: break = .false.
   contains
     procedure(fvalue), deferred :: value
     procedure(fgradient), deferred :: gradient
     procedure(fhessian), deferred :: hessian
     procedure(fprint), deferred :: print
  end type sr1_function

  abstract interface

     real(kind) function fvalue(this,x)
       import :: kind, sr1_function
       class(sr1_function) :: this
       real(kind), dimension(:), intent(in) :: x
     end function fvalue

     subroutine fgradient(this,x,fgrad)
       import :: kind, sr1_function
       class(sr1_function) :: this
       real(kind), dimension(:), intent(in) :: x
       real(kind), dimension(:), intent(out) :: fgrad
     end subroutine fgradient

     subroutine fhessian(this,x,fhess)
       import :: kind, sr1_function
       class(sr1_function) :: this
       real(kind), dimension(:), intent(in) :: x
       real(kind), dimension(:,:), intent(out) :: fhess

     end subroutine fhessian

     subroutine fprint(this,iter,x,fvalue,fgrad,fhess,step,delta,rho,sinfo)
       import :: kind, sr1_function
       class(sr1_function) :: this
       integer, intent(in) :: iter, sinfo
       real(kind), intent(in) :: fvalue,delta,rho
       real(kind), dimension(:), intent(in) :: x,fgrad,step
       real(kind), dimension(:,:), intent(in) :: fhess

     end subroutine fprint


  end interface


  public :: sr1, sr1rob, sr1update

  logical :: debug = .true.
  real(kind), parameter :: eta = 0.5e-3  ! 0 < eta < 1e-3
  real(kind), parameter :: rpar = sqrt(epsilon(rpar))  ! 0 < r < 1
  real(kind), parameter :: rtol = 0.1
  real(kind), parameter :: atol = 0

contains


  subroutine sr1(fcn,x,hess,delta,maxdev,maxfev,quasi,verbose,info,msg)

    use trust_step

    class(sr1_function), intent(in) :: fcn
    real(kind), dimension(:), intent(in out) :: x
    real(kind), dimension(:,:), optional, intent(out) :: hess
    real(kind), optional, intent(in) :: delta, maxdev
    integer, optional, intent(in) :: maxfev
    logical, optional, intent(in) :: quasi, verbose
    integer, optional, intent(out) :: info
    character(len=*), optional, intent(out) :: msg

    real(kind), dimension(size(x),size(x)) :: B
    real(kind), dimension(size(x)) :: fgrad, fg, y, p, xp
    real(kind) :: tol, dlt, fval, fv, ared, pred, rho, pnorm
    integer :: n,i,iter,maxev,trinfo,iflag
    character(len=80) :: fmt
    logical :: qnewton, verb

    iflag = 666
    n = size(x)

    if( present(verbose) ) then
       verb = verbose
    else
       verb = .false.
    end if

    ! the convergence tolerance
    if( present(maxdev) ) then
       tol = maxdev
    else
       tol = 10*n*epsilon(x)
    end if

    ! Delta is both initial and maximal radius of the trust-region circle.
    if( present(delta) ) then
       if( .not. (delta > 0 ) ) stop 'Assumption (delta > 0) unsatisfied.'
       dlt = delta
    else
       dlt = 1 !norm2(x)
    end if

    ! the maximum number of iterations (the function evaluations)
    if( present(maxfev) ) then
       maxev = maxfev
    else
       maxev = 3*n*precision(x)
    end if

    ! quasi-Newton Hessian approximation
    if( present(quasi) ) then
       qnewton = quasi
    else
       qnewton = .true.
    end if


    ! the initial Hessian approximation
    if( qnewton ) then
       B = 0
       forall( i = 1:n ) B(i,i) = 1
    else
       call fcn%hessian(x,B)
    end if

    fval = fcn%value(x)
    call fcn%gradient(x,fgrad)

    if( verb ) then
       p = 0
       rho = 0
       call fcn%print(0,x,fval,fgrad,B,p,dlt,rho,0)
    end if

    if( debug ) then
       write(fmt,'("(i0,a,",i0,"f8.3,a,2es9.2,a,3f8.3,i3)")') size(x)
       write(*,fmt) 0,':',x," >",fval,norm2(fgrad)," |",dlt
    end if

    do iter = 1, maxev

       if( norm2(fgrad) < tol ) then
          iflag = 1
          exit
       end if

       call trstep(B,fgrad,dlt,rtol,atol,p,trinfo)
       ! what about convergnece?
       ! info = 0 ok
       ! info = 1 indefinite
       ! info = 3 detected no convergence

       if( verb ) call fcn%print(iter,x,fval,fgrad,B,p,dlt,rho,trinfo)
       if( debug ) write(*,fmt) &
            iter,': ',x," >",fval,norm2(fgrad)," |",dlt,norm2(p),rho,trinfo

       ! attempt the step
       xp = x + p
       fv = fcn%value(xp)
       if( fcn%break ) then
          iflag = 6
          exit
       end if
       call fcn%gradient(xp,fg)
       if( fcn%break ) then
          iflag = 6
          exit
       end if

       y = fg - fgrad
       ared = fval - fv
       pred = -(dot_product(fgrad,p) + 0.5*dot_product(p,matmul(B,p)))

!          write(*,*) 'ared,pred:',ared,pred

       if( abs(ared) < epsilon(ared) .and. abs(pred) < epsilon(pred) ) then
          rho = 1
       else if( abs(pred) < epsilon(pred) .and. abs(ared) > 0 ) then
          rho = 0
       else
          rho = ared / pred
       end if

       if( rho > eta ) then
          x = xp
          fval = fv
          fgrad = fg
       else
          x = x
       end if

       pnorm = norm2(p)
       if( rho > 0.75 ) then
          if( pnorm  <= 0.8*dlt ) then
             dlt = dlt
          else
             !                dlt = min(2*dlt,delta);
             dlt = 2*dlt;
          end if
       else if( 0.25 <= rho .and. rho <= 0.75 ) then
          dlt = dlt
       else ! if ( rho < 0.25 ) then
          dlt = dlt / 2;
          !             dlt = max(dlt/2, tol)
       end if

       if( qnewton ) then
          call sr1update(p,y,pnorm,B)
       else
          call fcn%hessian(x,B)
          if( fcn%break ) then
             iflag = 6
             exit
          end if
       end if

    end do

    if( present(hess) ) hess = B
    if( present(info) ) info = iflag
    if( present(msg) ) then
       if( iflag == 1 ) then
          msg = 'The optimum is localised with the required accuracy.'
       else if( iflag == 5 ) then
          msg = 'The initial point is placed out of trust region.'
       else if( iflag == 6 ) then
          msg = 'This optimalisation run interrupted on request.'
       else
          msg = ''
       end if
    end if

  end subroutine sr1


  subroutine sr1rob(fcn,x,hess,delta,maxdev,maxfev,quasi,verbose,info,msg)

    use trust_step

    class(sr1_function), intent(in) :: fcn
    real(kind), dimension(:), intent(in out) :: x
    real(kind), dimension(:,:), optional, intent(out) :: hess
    real(kind), optional, intent(in) :: delta, maxdev
    integer, optional, intent(in) :: maxfev
    logical, optional, intent(in) :: quasi, verbose
    integer, optional, intent(out) :: info
    character(len=*), optional, intent(out) :: msg

    real(kind), dimension(size(x),size(x)) :: B
    real(kind), dimension(size(x)) :: fgrad, fg, y, p, xp
    real(kind) :: tol, dlt, fval, fv, ared, pred, rho, pnorm, smin
    integer :: n,i,iter,maxev,trinfo,iflag
    character(len=80) :: fmt
    logical :: qnewton, verb

    iflag = 666
    n = size(x)

    if( present(verbose) ) then
       verb = verbose
    else
       verb = .false.
    end if

    ! the convergence tolerance
    if( present(maxdev) ) then
       tol = maxdev
    else
       tol = 10*n*epsilon(x)
    end if

    ! Delta is both initial and maximal radius of the trust-region circle.
    if( present(delta) ) then
       if( .not. (delta > 0 ) ) stop 'Assumption (delta > 0) unsatisfied.'
       dlt = delta
    else
       dlt = 1
    end if

    ! the maximum number of iterations (the function evaluations)
    if( present(maxfev) ) then
       maxev = maxfev
    else
       maxev = 3*n*precision(x)
    end if

    ! minimal scale parameter
    smin = 0.1*x(n)/sqrt(real(666.0))
!    smin = 0.5*x(n)/sqrt(real(n))
!    smin = 0

    dlt = min(dlt,0.1*(x(n) - smin))
    write(*,*) 'ssss',delta,0.5*(x(n) - smin),dlt,smin

    ! quasi-Newton Hessian approximation
    if( present(quasi) ) then
       qnewton = quasi
    else
       qnewton = .true.
    end if

    ! the initial Hessian approximation
    if( qnewton ) then
       B = 0
       forall( i = 1:n ) B(i,i) = 1
       B(1,1) = 1 / x(2)**2
       B(2,2) = 2 / x(2)**2
!       call fcn%hessian(x,B)
    else
       call fcn%hessian(x,B)
    end if

    fval = fcn%value(x)
    call fcn%gradient(x,fgrad)

    if( verb ) then
       p = 0
       rho = 0
       call fcn%print(0,x,fval,fgrad,B,p,dlt,rho,0)
    end if

    if( debug ) then
       write(fmt,'("(i0,a,",i0,"f8.3,a,2es9.2,a,2es8.1,f5.2i3)")') n
       write(*,fmt) 0,':',x," >",fval,norm2(fgrad)," |",dlt
    end if

    do iter = 1, maxev

       if( norm2(fgrad) < tol ) then
          iflag = 1
          exit
       end if

!       B(1,:) = [ 5.743E-01,      1.385E+00]
!       B(2,:) = [1.385E+00,      2.471E+00]
!       fgrad = [-8.668E-01  ,   -4.070E-01]
       call trstep(B,fgrad,dlt,rtol,atol,p,trinfo)
       ! what about convergnece?
       ! info = 0 ok
       ! info = 1 indefinite
       ! info = 3 detected no convergence

       if( verb ) call fcn%print(iter,x,fval,fgrad,B,p,dlt,rho,trinfo)
       if( debug ) write(*,fmt) &
            iter,': ',x," > ",fval,norm2(fgrad)," |",dlt,norm2(p),rho,trinfo

!       if( trinfo /= 1 ) stop 'trstep'

       ! attempt the step
       xp = x + p

       if( xp(n) > smin ) then

          fv = fcn%value(xp)
          if( fcn%break ) then
             iflag = 6
             exit
          end if
          call fcn%gradient(xp,fg)
          if( fcn%break ) then
             iflag = 6
             exit
          end if

          y = fg - fgrad
          ared = fval - fv
          pred = -(dot_product(fgrad,p) + 0.5*dot_product(p,matmul(B,p)))

!          write(*,*) 'ared,pred:',ared,pred

          if( abs(ared) < epsilon(ared) .and. abs(pred) < epsilon(pred) ) then
             rho = 1
          else if( abs(pred) < epsilon(pred) .and. abs(ared) > 0 ) then
             rho = 0
          else
             rho = ared / pred
          end if

       else ! s <= 0

          rho = 0

       end if

       if( rho > eta ) then
          x = xp
          fval = fv
          fgrad = fg
       else
          x = x
       end if

       pnorm = norm2(p)
       if( rho > 0.75 ) then
          if( pnorm  > 0.8*dlt ) then
             dlt = 2*dlt
          end if
       else if( 0.25 <= rho .and. rho <= 0.75 ) then
          dlt = dlt
       else ! if ( rho < 0.25 ) then
          dlt = dlt / 2;
       end if
       dlt = min(dlt,x(n)-smin);
       !       write(*,*) '****',dlt,x(n),pnorm,rho

       if( qnewton ) then
          call sr1update(p,y,pnorm,B)
       else
          call fcn%hessian(x,B)
          if( fcn%break ) then
             iflag = 6
             exit
          end if
       end if

       if( iter == maxev ) iflag = 5

    end do

    if( present(hess) ) hess = B
    if( present(info) ) info = iflag
    if( present(msg) ) then
       if( iflag == 1 ) then
          msg = 'The optimum is localised with the required accuracy.'
       else if( iflag == 6 ) then
          msg = 'This optimalisation interrupted on request.'
       else
          msg = ''
       end if
    end if

  end subroutine sr1rob


  subroutine sr1update(p,y,snorm,B)

    real(kind), dimension(:), intent(in) :: p,y
    real(kind), intent(in) :: snorm
    real(kind), dimension(:,:), intent(in out) :: B

    real(kind), dimension(size(p)) :: res
    real(kind) :: rnorm, q
    integer :: n,i,j

    n = size(p)
    res = y - matmul(B,p)
    rnorm = norm2(res)
    q = dot_product(p,res)
!    write(*,*) real(res)
!    write(*,*) real(y)
!    write(*,*) real(q)

    if( abs(q) >= rpar*snorm*rnorm .and. snorm*rnorm > 0 ) then
       forall( i = 1:n, j = 1:n ) B(i,j) = B(i,j) + res(i)*res(j)/q
    else
       B = B
    end if

  end subroutine sr1update


end module symmetric_rank_1
