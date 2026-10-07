!
!  Trust-region step  on base of dgqt()
!
!  Copyright © 2021 F.Hroch (hroch@physics.muni.cz)
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


module trust_step

  use iso_fortran_env
  implicit none
  private

  integer, parameter :: kind = REAL64

  public :: trstep

  logical :: debug = .true.


contains

  subroutine trstep(B,g,delta,rtol,atol,p,info)

    ! This subroutine re-implements subroutine dgqt.f published
    ! as a part of MINPACK-2 project, and downloaded from
    ! https://ftp.mcs.anl.gov/pub/MINPACK-2/.
    !
    ! Th subroutine (dgqt.f) is on base of ideas "Computing a Trust Region Step"
    ! by More & Sorenson (1981), SIAM J. Sci. Statist. Comput. 4 (1983), 553-572;
    ! refences of theorems & equations are included.
    !
    ! The procedure is also described in the second edition of Numerical
    ! Optimization by Jorge Nocedal and Stephen J. Wright: the algorithm 4.3.
    !
    ! Notes.
    ! The reimplementation keeps the lambda parameter safeguarding handling,
    ! some convergence criteria, and info meaning only. The structure
    ! and the steps are completely re-wroted:
    !  * The program is in Fortran 2008+, particulary: the vector operations,
    !    and functions, are frequently utilised.
    !  * I have implemented myself Cholesky factorisation (no pivoting) for
    !    the upper triangular matrix, which is uncommon (I derived the
    !    factorisation formulas myself). There is no dependence on LAPACK.
    !  * The estimate of z-step corresponds to destsv.f subroutine of MINPACK-2,
    !    but the implementation is new (see comments).
    !  * The intermediate computations of the procedure has been compared
    !    to instrumented original dgqt(). The results tightly corresponds;
    !    no full exact match is possible (I believes).
    !  * I've attempted to keep API, names of varibles, etc. the same.

    use cholesky

    real(kind), dimension(:,:), intent(in) :: B
    real(kind), dimension(:), intent(in) :: g
    real(kind), intent(in) :: delta, rtol, atol
    real(kind), dimension(:), intent(out) :: p
    integer, intent(out) :: info

    real(kind), dimension(size(p)) :: z, u, q, gbnorm, diag
    real(kind), dimension(size(b,1),size(b,2)) :: B1, R
    real(kind) :: gnorm, bnorm, pnorm, qnorm, rznorm, znorm, rpnorm, &
         par, parc, pars, parl, paru, tau, pz, t
    integer :: i,n,iter,maxiter,indef
    character(len=120) :: fmt
    logical :: rednc
!    logical :: hard = .true.

    ! initialisation
    n = size(B,1)
!    sig1 = 1e-3
    z = real(0,kind)
    p = real(0,kind)

!    if( debug ) &

!,size(p),'f10.5,",",es15.3)') '(a,2i3,1x,es15.3,",",',size(p),'f10.5,",",es15.3)'


    ! Norms (https://en.wikipedia.org/wiki/Norm_(mathematics))
    ! gnorm .. Euclidean norm L2 of g
    ! rbnorm(n) .. sum of absolute values of off-diagonal elements per row
    !              (https://en.wikipedia.org/wiki/Gershgorin_circle_theorem)
    ! bnorm .. modified L1,Linfty norm of the matrix B
    !          (dgqt() implementation computes maximum of L1 norms
    !           rather than their sums; or, maximum of rows rather than the maximum
    !           of all elements).

    gnorm = norm2(g)
    bnorm = real(0,kind)
    do i = 1, n
!       bnorm = max(sum(abs(B(:,i))),bnorm)
       gbnorm(i) = sum(abs(B(:,i)))
       bnorm = max(bnorm,gbnorm(i))
       gbnorm(i) = gbnorm(i) - abs(B(i,i))
    end do
!    bnorm = sum(abs(B))
!    write(*,*) gnorm,bnorm,real(gbnorm),sum(abs(B))

    !    bnorm = sum(abs(B))
    diag = [(B(i,i),i=1,n)]
    par = real(0,kind)
    pars = maxval(-diag)
!    pars = maxval(-[(B(i,i),i=1,n)])
    parl = max(real(0,kind),gnorm/delta - maxval(gbnorm+diag),pars)
    paru = max(real(0,kind),gnorm/delta + maxval(gbnorm-diag))
!    parl = max(real(0,kind),gnorm/delta - bnorm,pars)
!    paru = max(real(0,kind),gnorm/delta + bnorm)
    paru = max(paru,(1 + rtol)*parl)

    ! put parameter inside the bounds
    par = min(max(parl,par),paru)

    ! prevent a cataclysmatic break down
    paru = max(paru,(1+rtol)*parl)

    if( debug ) then
       write(fmt,'("(5x,a,2i3,1x,f10.2,2(a,",i0,"f8.2,a),f10.2)")') size(p)
!       write(*,*) fmt
       write(*,'(5x,a,4es13.3)') 'init par:',parl, paru, pars, par
       do i = 1, n
          write(*,'(5x,a,i0,a,100es18.6)') 'b(',i,',:):',b(i,:)
       end do
       write(*,'(5x,a,100es18.6)') 'g(:):',g(:)
    end if

    info = 0
    maxiter = precision(p)
    do iter = 1, maxiter

       ! Safeguard lambda
!       lambda = min(max(lambda,lambda_L),lambda_U)
!       if( lambda <= lambda_S ) &
       !            lambda = max(0.001*lambda_U,sqrt(lambda_L*lambda_U))

!       write(*,*) par,pars,parl,paru
       if( par <= pars .and. paru*parl > 0 ) then
          par = max(real(1e-3,kind),sqrt(parl*paru))
       end if

       ! update matrix
       B1 = B
       forall( i = 1:n ) B1(i,i) = B(i,i) + par

       if( debug ) then
          do i = 1, n
             write(*,'(5x,a,i0,a,100es18.6)') 'b1(',i,',:):',b1(i,:)
          end do
       end if

!       if( debug ) 'trstep hessian ',b1


       ! Cholesky decomposition, and the step estimation, if possible
       call chol_rrfac(B1,R,indef)
!       write(*,*) 'indef=',indef, lambda
       if( indef == 0 ) then
          call chol_rrsolve(R,-g,p)
          call chol_rrtsolve(R,p,q)
          qnorm = norm2(q)
          pnorm = norm2(p)

!       do i = 1, n
!          write(*,'(5x,a,i0,a,100es18.6)') 'R(',i,',:):',R(i,:)
!       end do

!       write(*,*) 'p:',p
!          write(*,*) 'q:',q,matmul(transpose(R),q)-p
!          write(*,*) 'pnorm',pnorm,p,qnorm,norm2(matmul(R,p)),pnorm/delta,&
!               abs(pnorm - delta)/delta

!          write(*,*) par,pnorm
          ! inside a hole having parabolic surface
          if( abs(par) < epsilon(par) .and. pnorm/delta < (1 + rtol) ) then
             info = 1
             exit
          end if

          ! on the edge
          if( par > 0  .and. abs(pnorm - delta)/delta < rtol ) then
             info = 1
             exit
          end if


!          if( abs(pnorm - delta) <= rtol*delta .or. &
!               (abs(lambda) < epsilon(lambda) .and. pnorm <= (1+rtol)*delta) ) &
!               info = 1


!          if( debug_qstep ) then
!             block
!               real(kind), dimension(size(z)) :: z1
!               call qstep(B1,delta,z1)
!               write(*,*) 'expected z:',real(z1)
!             end block
!          end if

!          rznorm = norm2(z)
!          tau = 0
!          z = 0
!          rznorm = 0
!          p = p + tau*z

          ! compute the z step (to update pars, and to solve the hard case)
          call zstep(R,z)

!          rznorm = norm2(z)

          ! DGQT comment
          ! Compute a direction of negative curvature and use this
          ! information to improve pars.
          rznorm = norm2(matmul(R,z))
          pars = max(pars,par - rznorm**2)

!          write(*,*) 'rznorm:',real(rznorm),real(pars)

          rednc = .false.
          if( pnorm < delta ) then ! .and. par >= 0

             pz = dot_product(p,z) / delta
             t = (delta - pnorm)*((delta + pnorm)/delta)
             tau = t/(abs(pz) + sqrt(pz**2 + t/delta))
             tau = sign(tau,pz)

             rznorm = abs(tau)*rznorm
!             rpnorm = norm2(matmul(R,p))
!             rednc = (rznorm/delta)**2 + lambda*(pnorm/delta)**2 <= lambda

             rednc = (rznorm/delta)**2 + par*(pnorm/delta)**2 < par

!             if( (rznorm/delta)**2 + par*(pnorm/delta)**2 < par ) then
             !if( (rznorm/delta)**2  < 0.5*((rpnorm/delta)**2 + par) )then
!                p = p + tau*z
!                exit
!                pnorm = norm2(p)
!             end if

             rpnorm = norm2(matmul(R,p))
             if( (rznorm/delta)**2 <= rtol*(2-rtol)*(par+(rpnorm/delta))**2 ) then
!             if( (rznorm/delta)**2 <= rtol*(2-rtol)*(par+(qnorm/delta))**2 ) then
                info = 1
                exit
             else if( 0.5*(par+(qnorm/delta)**2) <= (atol/delta)/delta ) then
                info = 2
                exit
             else if( pnorm < epsilon(pnorm) ) then
                info = 1
                exit
             end if
          end if

! c           Compute the Newton correction parc to par.
          if( pnorm > 0 ) then
             parc = (pnorm/qnorm)**2*((pnorm - delta)/delta)
          else
             parc = -par
          end if
!          write(*,'(a,4f15.5,a)') '>',pnorm,qnorm,pnorm - delta,parc,'<<<<<'

!          if( qnorm < epsilon(qnorm) ) then
!             dlam = -lambda
!          else

          ! update the parameter by Newton iterations
!          write(*,*) p
          !          write(*,*) q
!          do jiter = 1, precision(lambda)
!             call chol_rrsolve(transpose(R),p,q)
             !             call chol_rrsolve(R,p,q)

!             B1 = R
!             B1(2,1) = B(2,1)

!             write(*,*) B1(1,:),'...'
!             write(*,*) B1(2,:),'...'
!             write(*,*) p/pnorm

!             call chol_rsolve(R,p/pnorm,u)
             !call chol_solve(B1,p/pnorm,u)
!             qnorm = norm2(u)
!             write(*,*) u,qnorm,norm2(u)
!             dlam = ((pnorm - delta)/delta) / qnorm**2

!             dlam = (pnorm/qnorm)**2*((pnorm - delta)/delta)
!                lambda = lambda + (pnorm/qnorm)**2*(pnorm/delta - 1)
!                lambda_S = max(lambda_S,lambda-norm2(matmul(R,z))**2)
!          end if

!          if( pnorm > delta ) parl = max(parl,par)
!          if( pnorm < delta ) paru = min(paru,par)


!          write(*,*) 'xx',real(dlam),real(lambda_L),real(lambda_U)

!             write(*,*) jiter,lambda,qnorm
!          end do

       else ! the indefinite case

!          write(*,*) 'b1(1,:)',b1(1,:)
!          write(*,*) 'b1(2,:)',b1(2,:)

!          write(*,*) 'indefinite',indef

!!$          if( indef > 1 ) then
!!$
!!$             write(*,*) R(1,:),'before'
!!$             write(*,*) R(2,:)
!!$             R(1:indef-1,indef) = B(indef,1:indef-1)
!!$             R(indef,indef) = diag(indef) + par
!!$             write(*,*) R(1,:),'after'
!!$             write(*,*) R(2,:)
!!$
!!$             k = indef-1
!!$             write(*,*) 'b:',R(1:k,indef)
!!$             call chol_rsolve(transpose(R(1:k,1:k)),R(1:k,indef),u(1:k))
!!$             write(*,*) u(1:k)
!!$
!!$             R(indef,indef) = R(indef,indef) - norm2(u(1:k))**2
!!$             call chol_rsolve(R(1:k,1:k),u(1:k),gbnorm(1:k))
!!$             write(*,*) gbnorm(1:k)
!!$
!!$!             call chol_solve(transpose(B1(1:k,1:k)),B1(1:k,indef),u(1:k))
!!$!             call chol_rrsolve(R,-g,p)
!!$
!!$!             e = 0
!!$!             e(k) = 1
!!$!             call chol_rrsolve(R(1:k,1:k),e(1:k),u(1:k))
!!$!             mu = dot_product(e,matmul(B1,e))
!!$          end if
!!$
!!$
!!$!          k = indef-1
!!$!          R(indef,indef) = B(indef,indef) + lambda
!!$!          call chol_rrsolve(transpose(B1(1:k,1:k)),B1(1:k,indef),u(1:k))
!!$!          B1(indef,indef) = B(indef,indef) - norm2(u(1:k))**2
!!$!          call chol_rrtsolve(B1(1:k,1:k),u(1:k),u(1:k))
!!$
!!$          gbnorm(indef) = -1
!!$          unorm = norm2(gbnorm(1:indef))
!!$          parc = -(R(indef,indef)/unorm)/unorm


          do i = 1, indef - 1
             R(i,i:n) = R(i,i:n) / R(i,i)
          end do
          R(indef,indef) = real(1,kind)

          u = real(0,kind)
          u(indef) = real(1,kind)
          call chol_rsolve(R(1:indef,1:indef),u(1:indef),z(1:indef))
          znorm = norm2(z(1:indef))
          q(1:indef) = matmul(B1(1:indef,1:indef),z(1:indef))
          parc = -dot_product(z(1:indef),q(1:indef)) / znorm**2
          pars = max(pars,par,par+parc)


!          write(*,*) 'pars:',real(parc),real(pars)


!          write(*,*) 'indef...',indef,real(u),real(unorm),real(lambda_S)

!          write(*,*) 'R(1,:):',R(1,:)
!          write(*,*) 'R(2,:):',R(2,:)
!          write(*,*) 'gbnorm',gbnorm

          paru = max(paru,(1+rtol)*pars)
       end if

       if( iter == maxiter ) then
          info = 4

          block
            real(kind) :: dp
            integer :: j

            dp = paru / 100.0
            open(3,file='/tmp/p')
            do j = 0,100
               par = j*dp
               B1 = B
               forall( i = 1:n ) B1(i,i) = B(i,i) + par
               call chol_rrfac(B1,R,indef)
               if( indef == 0 ) then
                  call chol_rrsolve(R,-g,p)
                  write(3,*) par,norm2(p),1/delta - 1/norm2(p)
               end if
            end do

            close(3)
          end block

!          stop
          exit
       end if

       if( paru <= (1 + 0.5*rtol)*pars ) then
          info = 3
          exit
       end if

       if( abs(paru) < epsilon(paru) ) then
          info = 2
          exit
       end if

!       if( info == 0 ) then
!          if( iter == maxiter ) info = 4
!          if( lambda_U <= (1 + 0.5*rtol)*lambda_S ) info = 3
!          if( abs(lambda_U) < epsilon(lambda) ) info = 2
!       end if


!       if( info /= 0 ) then

!          if( .not. (dot_product(p,matmul(B1,p)) > 0) ) &
!               stop '.not. positive definite matrix'

          ! zvazit vypocet f
!          p = p + tau*z
!          exit
!       end if


       ! check convergence
!       d = norm2(p)
!       if( pnorm < delta .and. lambda < epsilon(lambda) ) return
!       if( abs(delta - pnorm) < sig1*delta ) return
!       if( hard ) then
!          p = p + tau*z
!       end if

       ! update safeguard parameters
       parl = max(parl,pars)
!       lambda_U = min(lambda_U,lambda)
       par = max(par + parc,parl)

       ! update lambda
!       if( iflag ) then
!          call 32
!       else
!          lambda =


       if( debug ) then

!          call trstep_debug(B1,delta)

!          write(fmt,'(a,i0,a)') '(a,2i3,1x,es15.3,",",',size(p),'f10.5,",",es15.3)'
!          write(*,*) trim(fmt)
          write(*,fmt) "trstep:",iter,indef,par,", ",p,", ","",z,", ",pnorm
!          write(*,'(4f15.5)') parl,paru,pars,par
!          write(*,'(4f15.5)') pnorm, qnorm,pnorm-delta,parc
       end if

    end do

    if( rednc .and. (info == 1 .or. info == 2) ) then
       p = p + tau*z
    end if

!    write(*,*) fmt
    if( debug ) write(*,fmt) "trstep:",iter,indef,par,", ",p,", ","",z,", ",pnorm


    ! Estimate lambda which solves ||p(lambda)|| = delta, where
    ! (B + lambda I) s = -g

!    call chol_rrfac(b1,r,indef)
!
!    if( indef == 0 ) then
!       call chol_rrsolve(r,-g,p)
!       if( enorm(p) < delta ) then
!          continue
!       end if
!    end if



!!$    call qrfac(b,q,r)
!!$
!!$    pdef = all(r > 0)
!!$
!!$    if( pdef ) then
!!$       call qrsolv(b,g,s)
!!$       if( hypot(s,s) <= delta ) return
!!$    end if
!!$
!!$
!!$    n = size(b,1)
!!$    b1 = b
!!$    lambda = 1
!!$    do iter = 0, 666
!!$       do i = 1, n
!!$          b1(i,i) = b(i,i) + lambda
!!$       end do
!!$
!!$       call qrsolve(b1,-g,p)
!!$
!!$       call llfac(b1,r,r1,iflag)
!!$       call qrsolve(r1,p,q)
!!$
!!$       phorm = hypot(p,p)
!!$       qhorm = hypot(q,q)
!!$       lambda = lambda + (pnorm/qnorm)**2*(pnorm - delta)/delta
!!$
!!$    end do

  end subroutine trstep


  subroutine zstep(R,z)

    use cholesky

    real(kind), dimension(:,:), intent(in) :: R
!    real(kind), intent(in) :: delta
    real(kind), dimension(:), intent(out) :: z

    real(kind), parameter :: zero = real(0,kind)
    real(kind), parameter :: one = real(1,kind)
    real(kind), dimension(size(z)) :: e
    real(kind), dimension(size(z)), target :: w1,w2
    real(kind), dimension(:), pointer :: w
    real(kind) :: r11,w1norm,w2norm
    integer :: n, i


!    write(*,*) 'R(1,:)',real(R(1,:))
!    write(*,*) 'R(2,:)',real(R(2,:))

    n = size(z)
    z = zero


    r11 = abs(R(1,1))
!    if( abs(r11) < epsilon(R) ) then
       ! this check is, perhaps, unnecessary: the z-step is only computed
       ! in case of any positive definite matrix, eg. for all R(i,i) > 0

       !       rznorm = one
!       z = zero
!       z(1) = one
!       return
!    end if

!    e = 1
!    u = zero
!    s1 = 0
!    s2 = 0
    do i = 1, n
       !       sw = dot_product(R(1:i-1,i),w(1:i-1))
       e(i) = r11
       call chol_rrtsolve(R(1:i,1:i),e(1:i),w1(1:i))
!       w1 = norm2(w(1:i))
!       w1 = sum(abs(w(1:i)))
       e(i) =-r11
       call chol_rrtsolve(R(1:i,1:i),e(1:i),w2(1:i))
!       w2 = norm2(w(1:i))
!       w2 = sum(abs(w(1:i)))

       w1norm = sum(abs(w1(1:i)))
       w2norm = sum(abs(w2(1:i)))
!       w1 = e - sw
!       w2 =-e - sw
       if( w1norm > w2norm ) then
          e(i) = r11
       else
          e(i) =-r11
!          w => w2
       end if
!       w(i) = w(i) / R(i,i)
    end do

    if( w1norm > w2norm ) then
       w => w1
    else
       w => w2
    end if

    !    call chol_rrtsolve(R,u,w)
!    if( w1norm < w2norm ) then
!       call chol_rsolve(R,w1,z)
!    else
!       call chol_rsolve(R,w2,z)
    !    end if

    call chol_rsolve(R,w,z)
    z = z / norm2(z)

!    write(*,*) 'zstep',real(z)!,real(rznorm)


  end subroutine zstep

end module trust_step
