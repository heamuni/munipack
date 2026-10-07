!
!  The trust region step by Newton-Raphson method, based on eigendecomposition
!
!  Copyright © 2023-2026 F.Hroch (hroch@physics.muni.cz)
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
!  This is an implementation of More & Sorenson's (1981) trust-region step
!  restricted into 2D space. It exploit the direct (non-iterative) factorisation
!  of the eigenproblem. For description, see the second edition Numerical
!  Optimization by Jorge Nocedal and Stephen J. Wright, Ch. 4.3.

module qtrstep_REAL64

  use iso_fortran_env
  implicit none
  private

  real(REAL64), parameter :: macheps = epsilon(1.0_REAL64)
  logical, parameter :: debug = .false.

  public :: trstep_REAL64

contains


  subroutine trstep_REAL64(B,g,delta,rtol,p,pnorm,posdef,info,verbose)

    ! trstep() returns the step (p), its norm (pnorm),
    ! the positive definitnes of B matrix (posdef)
    ! inside or on the edge of a trust-region with delta radius
    !
    ! Normal return is indicated by:
    ! info = 1 the step is inside TR, pnorm < delta, Newton-Raphson step
    ! info = 2 steepest descent for double eigenvalue, pnorm == delta
    ! info = 3 outside solution, pnorm ~= delta
    ! info = 4 indefinite solution, pnorm ~= delta
    ! info = 5 hard case solution, pnorm ~= delta
    !
    ! A.b.normal return for
    ! info = 6 detected no convergence

    use modQ2eig_REAL64

    real(REAL64), dimension(:,:), intent(in) :: B
    real(REAL64), dimension(:), intent(in) :: g
    real(REAL64), intent(in) :: delta, rtol
    real(REAL64), dimension(:), intent(out) :: p
    real(REAL64), intent(out) :: pnorm
    logical, intent(out) :: posdef
    integer, intent(out) :: info
    logical, intent(in) :: verbose

    real(REAL64), parameter :: omega = 1.0 / 360.0
    real(REAL64), dimension(2,2) :: Q
    real(REAL64), dimension(2) :: lam, Qg
    real(REAL64) :: gnorm, d, tau
    logical :: converge
    integer :: m, nroots, lmin
    character(len=*), parameter :: fmt = '(a,2f10.5,2x,a,f10.5)'

    info = 6

    if( verbose ) then
       write(error_unit,*) 'trstep init: delta=',real(delta),' g=',real(g)
       write(error_unit,*) 'rtol=',rtol
       write(error_unit,*) 'B:',real(B)
    end if

    if( debug ) then
       if( size(g) /= 2 .or. size(B,1) /= 2 .or. size(B,2) /= 2 ) &
            error stop 'trstep(): second order vectors/matrixes are required.'
       if( abs(B(1,2) - B(2,1)) > macheps ) then
          write(error_unit,*) 'B:',B, 'res:',B(1,2)-B(2,1)
          error stop 'trstep(): A symmetry of matrix B(1,2) == B(2,1) unsatisfied.'
       end if
    end if

    ! determine eigen values and vectors
    call q2eig_REAL64(B,lam,Q,nroots)

    ! Is B positive definite?
    posdef = all(lam > macheps)

    if( verbose ) then
       write(error_unit,*) '   lam: ',lam,' posdef: ',posdef
       write(error_unit,*) 'Q(1,:): ',real(Q(1,:)), ' nroots=',nroots
       write(error_unit,*) 'Q(2,:): ',real(Q(2,:))
    end if

    Qg = matmul(transpose(Q),g)
    gnorm = norm2(g)

    if( verbose ) write(error_unit,'(a,l2,3g13.5)') &
         'posdef, |Qg|/D-lam2, -lam1, |Qg|:',posdef, &
         real(gnorm/delta-maxval(abs(lam))),-real(minval(lam)),real(gnorm/delta)

    if( posdef .and. gnorm/delta < maxval(lam) ) then

       p = -matmul(Q,Qg/lam)
       pnorm = norm2(p)

       if( debug ) then
          if( maxval(abs(matmul(B,p)+g)) > 1e3*(maxval(abs(B))+macheps) ) then
             write(error_unit,*) matmul(B,p), -g, matmul(B,p)+g, &
                  maxval(abs(matmul(B,p)+g)),1e3*maxval(abs(B))+macheps,&
                  ',',lam,Q
             stop '.NOT. matmul(B,p) = -g'
          end if
       end if

       ! par = 0: |p| < delta
       if( pnorm < delta ) then
          if( verbose ) write(error_unit,fmt) &
               'trstep (*) newton sol.:',p,' norm:',pnorm
          info = 1
          return
       end if

    end if

    if( posdef .and. nroots == 1 ) then
       p = -(delta/gnorm)*g
       pnorm = delta ! == norm2(p)
       info = 2
       if( verbose ) write(error_unit,fmt) &
            'trstep (/) steepest descent sol.:',p,' norm:',pnorm
       if( debug ) then
          if( abs(pnorm - delta)/delta > rtol ) error stop 'pnorm /= delta'
          if( all(abs(g - matmul(Q,Qg)) > 10*macheps) ) then
             write(error_unit,*) p,-(delta/gnorm)*matmul(Q,Qg),g,matmul(Q,Qg),&
                  abs(g-matmul(Q,Qg)) < macheps,norm2(p)-delta
             error stop 'g /= matmul(Q,Qg)'
          end if
       end if
       return
    end if

    if( verbose ) write(error_unit,*) &
         'Qg=',real(Qg),real(gnorm),real(omega),macheps
    lmin = minloc(lam,1)
    if( .not. posdef .and. abs(Qg(lmin)) < max(omega*gnorm,42*macheps) ) then

       info = 5
       m = 3 - lmin
       d = Qg(m)/(lam(m) - lam(lmin))
       pnorm = abs(d)
       d = sign(min(pnorm,delta),d)
       p = -d * Q(:,m)
       if( debug )  write(error_unit,*) 'imin=',lmin,pnorm,Qg,p,d,gnorm
       if( pnorm < delta ) then
          tau = sqrt((delta - pnorm)*(delta + pnorm))
          p = p - sign(tau,Qg(lmin))*Q(:,lmin)
          pnorm = norm2(p)
          if( verbose ) write(error_unit,*) &
               'final:',p, pnorm,tau,Qg(m)/(lam(m) - lam(lmin)),lmin
       end if
       if( verbose ) write(error_unit,fmt) &
            'trstep ($) hard case sol.:',p,' norm:',pnorm
    else

       call bpar(lam,posdef,Q,Qg,gnorm,delta,rtol,p,pnorm,converge,verbose)

       if( converge ) then
          if( posdef ) then
             info = 3
             if( verbose ) write(error_unit,fmt) &
                  'trstep (|) outside sol.:',p,' norm:',pnorm
          else
             info = 4
             if( verbose ) write(error_unit,fmt) &
                  'trstep (-) indefinite sol.:',p,' norm:',pnorm
          end if
          if( debug ) then
             if( abs(pnorm - delta)/delta > 2*rtol ) error stop 'NOT pnorm == delta.'
          end if
       else
          if( verbose ) then
             if( posdef ) then
                write(error_unit,fmt) 'trstep (!) outside sol: NO CONVERGENCE'
             else
                write(error_unit,fmt) 'trstep (!) indefinite sol: NO CONVERGENCE'
             end if
          end if
       end if
    end if

  end subroutine trstep_REAL64

  subroutine bpar(lam,posdef,Q,Qg,gnorm,delta,rtol,p,pnorm,converge,verbose)

    real(REAL64), dimension(:,:), intent(in) :: Q
    real(REAL64), dimension(:), intent(in) :: lam,Qg
    real(REAL64), intent(in) :: gnorm, delta, rtol
    real(REAL64), dimension(:), intent(out) :: p
    real(REAL64), intent(out) :: pnorm
    logical, intent(out) :: converge
    logical, intent(in) :: verbose, posdef

    character(len=*), parameter :: fmt = '(a,i3,2x,g0,es15.3,l2)'

    real(REAL64), dimension(size(lam)) :: y,w,u
    real(REAL64) :: h,p2,q2,par
    integer :: iter, lmin

    converge = .false.

    if( posdef ) then
       par = max(gnorm / delta - maxval(lam), 0.0_REAL64)
    else ! negative definite or semi-definite
       lmin = minloc(lam,1)
       par = abs(Qg(lmin)) / delta - (1 + macheps)*lam(lmin)
       if( debug ) write(error_unit,*) &
            'indef init:',par, gnorm / delta - (1 + macheps)*minval(lam), &
            max(-(1 + macheps)*minval(lam),macheps)
    end if

    if( verbose ) write(error_unit,fmt) 'bpar:',0,par

    do iter = 1, precision(par)

       y = lam + par

       if( any(abs(y) == 0) .or. par < 0 ) then
          if( debug ) then
             if( par < 0 ) write(error_unit,*) 'bpar(): par < 0 reached.'
             if( any(abs(y) == 0) ) &
                  write(error_unit,*) 'y==0?',iter,y,lam,par,Qg,delta,macheps
          end if
          return
       end if

       u = Qg / y
       w = u**2
       p2 = sum(w)
       q2 = sum(w/y) ! p2 > 0 is allways satisfied
       pnorm = sqrt(p2)
       h = (pnorm - delta) / delta

       if( verbose ) write(error_unit,fmt) 'bpar:',iter,par,h, abs(h) < rtol

       if( abs(h) < rtol ) then
          p = -matmul(Q,u)
          converge = .true.
          return
       end if

       par = par + (p2 / q2) * h

    end do

  end subroutine bpar

end module qtrstep_REAL64
