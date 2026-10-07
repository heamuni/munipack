!
!  An implementation of SQP (sequence quadratic programming) minimisation method
!  with a constrain.
!
!  Copyright © 2020 F.Hroch (hroch@physics.muni.cz)
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
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!
!  This implementation is due to the algorithms of the second edition
!  of Numerical Optimization by Jorge Nocedal and Stephen J. Wright.
!
!  trsqp ... trust-region SQP
!  linesqp ... line search SQP
!  locsqp ... local SQP

module sqps

  use minpacks
  use iso_fortran_env

  implicit none
  private

  integer, parameter :: kind = REAL64

  type, abstract, public :: sqp_function
     logical :: break = .false.
   contains
     procedure(fvalue), deferred :: value
     procedure(fgradient), deferred :: gradient
     procedure(fhessian), deferred :: hessian
     procedure(fprint), deferred :: print
  end type sqp_function

  abstract interface

     real(kind) function fvalue(this,x)
       import :: kind, sqp_function
       class(sqp_function) :: this
       real(kind), dimension(:), intent(in) :: x
     end function fvalue

     subroutine fgradient(this,x,fgrad)
       import :: kind, sqp_function
       class(sqp_function) :: this
       real(kind), dimension(:), intent(in) :: x
       real(kind), dimension(:), intent(out) :: fgrad
     end subroutine fgradient

     subroutine fhessian(this,x,fhess)
       import :: kind, sqp_function
       class(sqp_function) :: this
       real(kind), dimension(:), intent(in) :: x
       real(kind), dimension(:,:), intent(out) :: fhess

     end subroutine fhessian

     subroutine fprint(this,iter,x,fvalue,fgrad,fhess,step,delta,rho,sinfo)
       import :: kind, sqp_function
       class(sqp_function) :: this
       integer, intent(in) :: iter, sinfo
       real(kind), intent(in) :: fvalue,delta,rho
       real(kind), dimension(:), intent(in) :: x,fgrad,step
       real(kind), dimension(:,:), intent(in) :: fhess

     end subroutine fprint

  end interface

!  public :: linesqp, trsqp, locsqp
  public :: trsqp, locsqp

  logical, parameter :: quasi = .true.
  logical, parameter :: debug = .true.
!  logical, parameter :: debug_graph = .true.
!  logical, parameter :: debug_test = .false.

contains

  subroutine trsqp(fcn,x,cov,maxdev,maxfev,nprint,info,msg)

    !  This implementation of trust-region SQP method thightly follows
    !  the algorithm 18.4 (Byrd-Omojokun Trust-Region SQL Method) published
    !  in the second edition of Numerical Optimization by Jorge Nocedal and
    !  Stephen J. Wright.

    interface
       subroutine fcn(x,f,fvec,fjac,c,cvec,cjac,iflag)
         use iso_fortran_env
         real(REAL64), dimension(:), intent(in) :: x
         real(REAL64), intent(out) :: f, c
         real(REAL64), dimension(:), intent(out) :: fvec, cvec
         real(REAL64), dimension(:,:), intent(out) :: fjac, cjac
         integer, intent(in out) :: iflag
       end subroutine fcn
    end interface

    real(REAL64), dimension(:), intent(in out) :: x
    real(REAL64), dimension(:,:), optional, intent(out) :: cov
    real(REAL64), optional, intent(in) :: maxdev
    integer, optional, intent(in) :: maxfev, nprint
    integer, optional, intent(out) :: info
    character(len=*), optional, intent(out) :: msg

    character(len=80) :: emsg
    real(REAL64), dimension(size(x)) :: xlam
    real(REAL64), dimension(size(x),size(x)) :: hes, fjac, cjac, mat, imat
    real(REAL64), dimension(size(x)) :: p, fvec, cvec, fvec_old, cvec_old,x_old
!    real(REAL64), dimension(:,:), allocatable :: b, fjac, cjac
!    real(REAL64), dimension(:), allocatable :: p, fvec, cvec, fvec_old, cvec_old, &
!         x_old,dx, r, fx
    real(REAL64) :: tol,f,c,cr,lambda,mu,eta, gamma, rho, tdelta
    integer :: n, i, iter, maxev, nfev, iflag

    eta = 0.5
    gamma = 0.5
    lambda = -1

    n = size(x)
    nfev = 0
    iflag = 666
    emsg = 'UNDEFINED'

!    if( debug_test ) then
!       call test_cgsolve
!       stop 'All tests.'
!    end if

    if( present(maxdev) ) then
       tol = maxdev
    else
       tol = sqrt(epsilon(x))
    end if

    if( present(maxfev) ) then
       maxev = maxfev
    else
       maxev = 666
    end if

    tdelta = 10
    mu = 10*tdelta

    do iter = 1, n*precision(x)
       write(*,*) 'iter=',iter,'   trust-region radius=', real(tdelta)

       ! compute functions and their gradients
       call fcn(x,f,fvec,fjac,c,cvec,cjac,iflag)
       nfev = nfev + 1
       if( iflag < 0  ) goto 666
       if( nfev > maxev ) then
          info = 2
          goto 666
       end if

       ! computer inverse of cvec gradients
       block
!         real(REAL64), dimension(size(x),size(x)) :: aa
         integer :: i,j

         do concurrent ( i = 1:n )
            do concurrent ( j = 1:n )
               mat(i,j) = cvec(i)*cvec(j)
            end do
!            r(i) = sum(a(i,:)*fvec)
         end do
         !         call qrsolve(a,r,xlam)
         call qrinv(mat,imat)

!         aa = matmul(mat,imat)
!         do i = 1, n
!            write(*,*) mat(i,:)
!         end do
!stop

       end block



       ! lambdas
       block
         real(REAL64), dimension(size(x)) :: r
!!$         real(REAL64), dimension(:,:), allocatable :: a
!!$         real(REAL64), dimension(:), allocatable :: r
!!$         integer :: i,j
!!$
!!$         allocate(a(n,n),r(n))
!!$         do concurrent ( i = 1:n )
!!$            do concurrent ( j = 1:n )
!!$               a(i,j) = cvec(i)*cvec(j)
!!$            end do
!!$            r(i) = sum(a(i,:)*fvec)
!!$         end do
!!$
!!$         call qrsolve(a,r,xlam)

!         do i = 1, n
!            write(*,*) a(i,:)
!         end do
!         write(*,*) r

!         write(*,*) 'xlam:',xlam

!         do concurrent ( i = 1:n )
!            r(i) = sum(a(i,:)*fvec)
         !         end do

         r = matmul(mat,-fvec)
         xlam = matmul(imat,r)
         write(*,*) 'xlam:',xlam

       end block

       ! convergence test
       write(*,*) fvec - xlam*cvec, c
       write(*,*) x
       if( maxval(abs(fvec - xlam*cvec)) < tol .and. abs(c) < tol ) then
          iflag = 1
          exit
       end if

       ! estimate relaxation factors by the dogleg method
       block
         real(REAL64), dimension(size(x)) :: g, v!, pb
!         g = fvec - lambda*cvec
         ! newton's step
         !         pb = -matmul(ainv,cvec*c)
!         pb = -matmul(transpose(cvec),matmul(imat,c))
!         call dogleg(mat,imat,cvec,c,tdelta,v)
!         call dogleg(cvec,c,tdelta,v)
         cr = sum(cvec*v) + c
         write(*,*) 'cr:',cr
         write(*,*) 'v:',v

!         if( iter == 1 ) p = v
         p = v
       end block
!       stop 2

       ! update Hessian
       if( quasi ) then
          ! by the quasi-Hessian BFGS method
          block
!            real(REAL64), dimension(size(x)) :: s, u
            integer :: i

            if( iter == 1 ) then
               hes = 0
               do concurrent (i = 1: n)
                  hes(i,i) = 1
               end do
            else
!               s = x - x_old
!               u = (fvec - fvec_old) - lambda*(cvec - cvec_old)
!               call bfgs(s,u,b)
               call bfgs(x,fvec,cvec,x_old,fvec_old,cvec_old,lambda,hes)
            end if
          end block
       else ! compute the analytic Hessian with the constrain of second derivatives
          block
            real(REAL64), dimension(size(x)) :: fvec,cvec
            real(REAL64) :: f,c
            integer :: i

            call fcn(x+p,f,fvec,fjac,c,cvec,cjac,iflag)
            nfev = nfev + 1
            if( iflag < 0  ) goto 666
            if( nfev > maxev ) then
               info = 2
               goto 666
            end if
            hes = fjac - lambda*cjac
!            do concurrent ( i = 1:n )
!               b(i,:) = fjac(i,:) - xlam(i)*cjac(i,:)
!            end do

          end block
       end if

       ! localise solution
!       goto 44
       block
         real(REAL64), dimension(size(x)+1,size(x)+1) :: mat
         real(REAL64), dimension(size(x)+1) :: r, q
         integer :: m, i

         m = n + 1

         mat(1:n,1:n) = hes
         mat(1:n,m) = -cvec
         mat(m,1:n) = cvec
         mat(m,m) = real(0,REAL64)
         r(1:n) = -fvec !+ xlam*cvec
         r(m) = -c + cr!r(1)
         q(1:n) = p
         q(m) = lambda
         if( iter == 1 ) then
            q(m) = 0
         else
            q(m) = lambda
         end if

!         write(*,*) '>',real(r)
         do i = 1, m
            write(*,*) real(mat(i,:)),' *x = ',real(r(i))
         end do
!         write(*,*) real(mat(2,:))
!         write(*,*) real(mat(3,:))
         write(*,*) '>',real(q)

         !         call cgsolve(v)
!         call cgsolve(mat,r,q,tol,tdelta)
         call qrsolve(mat,r,q)


         write(*,*) 'QR solve =>',real(q)

         p = q(1:n)
         lambda = q(m)
       end block
44     continue

       goto 45
       block
         real(REAL64), dimension(size(x)+1,size(x)+1) :: mat
         real(REAL64), dimension(size(x)+1) :: r, q
         integer :: m, i

         m = n + 1

         mat(1:n,1:n) = hes
         mat(1:n,m) = -cvec
         mat(m,1:n) = cvec
         mat(m,m) = real(0,REAL64)
         r(1:n) = -fvec !+ xlam*cvec
         r(m) = -c + cr!r(1)
         q(1:n) = p
         if( iter == 1 ) then
            q(m) = 0.01*tdelta
         else
            q(m) = lambda
         end if

!         write(*,*) '>',real(r)
         do i = 1, m
            write(*,*) real(mat(i,:)),' *x = ',real(r(i))
         end do
!         write(*,*) real(mat(2,:))
!         write(*,*) real(mat(3,:))
         write(*,*) '>',real(q)

         call cgsolve(mat,r,q,tol,1e6*tdelta)
         p = q(1:n)
         lambda = q(m)

!         call cgsolve(hes,-fvec,p,tol,tdelta)
!         call qrsolve(b,-fvec,p)


         write(*,*) 'CG solution =>',real(p),real(lambda)

       end block
45 continue

!       stop 3

!       mu = 1
       block
         real(REAL64), dimension(0:1) :: fx, cx
         real(REAL64), dimension(size(x)) :: fvec1, cvec1
         real(REAL64), dimension(size(x),size(x)) :: fjac1,cjac1
         real(REAL64) :: q0, q1


         do i = 0, 1
            call fcn(x+i*p,fx(i),fvec1,fjac1,cx(i),cvec1,cjac1,iflag)
            nfev = nfev + 1
            if( iflag < 0  ) goto 666
            if( nfev > maxev ) then
               info = 2
               goto 666
            end if
         end do

         q0 = fx(0) + mu*abs(cx(0))
         q1 = fx(1) + sum(fvec1*p) + 0.5*sum(p*matmul(hes,p)) + mu*norm2(cx(1) + cvec1*p)

        rho = (phi2(fx(0),cx(0),mu) - phi2(fx(1),cx(1),mu)) / (q0 - q1)
!         rho = -rho

         write(*,*) q0,q1
         write(*,*) fx
         write(*,*) rho
       end block



!       stop 4

       ! update solution
       if( rho > eta ) then
          x = x + p
          tdelta = 2*tdelta
       else
!          tdelta = gamma*sqrt(sum(p**2))
          tdelta = norm2(p)
       end if

    end do

666 continue

    if( present(info) ) info = iflag
    if( present(msg) ) then
       msg = '.. undefined ...'
       if( iflag == 1 ) &
            msg = 'An uncertainity of the solution is under the required precision.'
       if( iflag == -1 ) msg = 'An interrupt by user.'
       if( iflag == 2 ) msg = 'Number of calls to fcn reached maxfev.'
       if( iflag == 3 ) msg = 'Any functions or constrains diverges.'
       if( iflag == 6 ) msg = emsg
    end if

  end subroutine trsqp


  function phi2(f,c,mu)

    real(REAL64) :: phi2
    real(REAL64), intent(in) :: f, mu, c

    phi2 = f + mu*abs(c)
!    q = f + sum(fvec*p) + sum(p*matmul(b,P)) / 2 + mu * enorm(cvec + matmul(b,p))

  end function phi2


!  subroutine dogleg(b,grad,delta,p)
  subroutine dogleg(mat,imat,cvec,c,tdelta,p)

    ! see also
    ! https://en.wikipedia.org/wiki/Powell%27s_dog_leg_method

    real(REAL64), dimension(:,:), intent(in) :: mat, imat
    real(REAL64), dimension(:), intent(in) :: cvec
    real(REAL64), intent(in) :: c,tdelta
    real(REAL64), dimension(:), intent(out) :: p

    real(REAL64), dimension(size(p)) :: pb, pu, pt, d
    real(REAL64) :: sgnorm, t, fmin, f, s
    logical :: found
    integer :: i

    ! Steepest descent
    d = cvec*c
    s = sum(d*matmul(mat,d))
    if( s > 0 ) then
       pu = -sum(d**2)/s * d
    else
       pu = 0
!       stop 's>0'
    end if

    ! Gauss-Newton step
!    pb = matmul(imat,-d)
    call qrsolve(mat,-d,pb)
!    pb = -c/cvec


!    call qrsolve(b,-grad,pb)
!    call qrsolve(b,-grad,pb)


    write(*,*) 'newton step:',real(pb)
    write(*,*) 'cauchy step:',real(pu)

    !    tmin = 0
    p = 0
    fmin = c**2
!    fmin = huge(fmin)
    found = .true.
    do i = 1, 20
       t = 0.1 * i
       if( t <= 1 ) then
          pt = t*pu
       else
          pt = pu + (t - 1)*(pb - pu)
       end if
!       f = enorm(cvec*pt + c)**2
       f = sum((cvec*pt + c)**2)
       write(*,'(f4.1,1x,2(f8.1,1x),2x,g0)') real(t),pt,real(log10(f))
!       if( enorm(pt) > 0.8*delta .and. i > 1 ) then
!          pt = pt0
!          found = .true.
!          exit
!       end if
       if( f < fmin .and. norm2(pt) <= 0.8*tdelta ) then
          fmin = f
!          tmin = t
          p = pt
!          write(*,*) i,t
          found = .true.
       end if
!       pt0 = pt

    end do

    write(*,*) found,real(p),real(log10(fmin))

    if( .not. found ) then
       stop 'not found.'
    end if

!stop
!!$    ! if step is inside trust-region, we're done
!!$    if( enorm(pb) <= delta ) then
!!$       p = pb
!!$       return
!!$    end if
!!$
!!$    ! Steepest descent
!!$!    pu = -grad * (enorm(grad)/sum(grad*matmul(b,grad)))
!!$
!!$    sgnorm = enorm(pu)
!!$    if( sgnorm  > delta ) then
!!$       p = (delta/sgnorm) * pu
!!$       return
!!$    end if
!!$
!!$    ! intermediate region
!!$    block
!!$      real(REAL64) :: a,b,c,d,a2,s,t1,t2
!!$
!!$      a = sum((pb - pu)**2)
!!$      b = 2*sum(pu*(pb - pu))
!!$      c = sum(pu**2) - delta**2
!!$      d = b**2 - 4*a*c
!!$      if( d > 0 ) then
!!$         a2 = 2*a
!!$         s = sqrt(d)
!!$         t1 = -(b + s) / a2
!!$         t2 = -(b - s) / a2
!!$         if( t2 + 1 > 0 ) then
!!$            t = t2 + 1
!!$         else !if( t1 + 1 > 0 ) then
!!$            t = t1 + 1
!!$         end if
!!$      else
!!$         t = 0
!!$      end if
!!$
!!$      p = pu + (t - 1)*(pb - pu)
!!$    end block

  end subroutine dogleg





  subroutine bfgs(x,fgrad,cgrad,x_old,fgrad_old,cgrad_old,lambda,B)

    ! damped BFGS updating, Procedure 18.2

    !    s = x - x_old
    !    u = (fgrad - fgrad_old) - lambda*(cgrad - cgrad_old)

    real(REAL64), dimension(:), intent(in) :: x,fgrad,cgrad,x_old,fgrad_old,cgrad_old
    real(REAL64), intent(in) :: lambda
    real(REAL64), dimension(:,:), intent(in out) :: B

    real(REAL64), dimension(size(x)) :: s,u,v,y
    real(REAL64) :: th, alpha, beta, q
    integer :: i,j,n

    n = size(x)

    s = x - x_old
    y = cgrad - cgrad_old
    u = (fgrad - fgrad_old) - lambda*y

    v = matmul(B,s)
    beta = dot_product(s,y)
    q = dot_product(s,matmul(B,s))

    write(*,*) beta, q

    ! curvature condition (6.7)
!    if( .not. (beta > 0) ) stop 'The curvature condition abort encountered.'
    if( beta < 0.01*q .or. abs(q) < epsilon(q) ) return

    alpha = sum(s*u)

    if( alpha >= 0.2*beta ) then
       th = real(1,REAL64)
    else
       th = (0.8*beta)/(beta - alpha)
    end if
!    u = th*u + (real(1,REAL64) - th)*v

!    write(*,*) beta,alpha,v,u
    alpha = sum(s*u)  ! u updated
    do concurrent ( i = 1:n )
       do concurrent ( j = i:n )
          B(i,j) = B(i,j) + (u(i)*u(j)) / alpha - (v(i)*v(j)) / beta
       end do
    end do
    do concurrent ( i = 2:n )
       do concurrent ( j = 1:i-1 )
          B(i,j) = B(j,i)
       end do
    end do

  end subroutine bfgs

!!$  subroutine bfgs(s,y,b,iflag)
!!$
!!$    real(REAL64), dimension(:), intent(in) :: s,y
!!$    real(REAL64), dimension(:,:), intent(in out) :: b
!!$    integer, intent(out) :: iflag
!!$!    real(REAL64), dimension(size(x)) :: z
!!$
!!$!    z = matmul(b,s)
!!$!    b = b - matmul(z,transpose(z))/sum(s*z) + matmul(y,transpose(y))/sum(y*s)
!!$    !block
!!$    real(REAL64), dimension(:), allocatable :: z
!!$    !    real(REAL64), dimension(:,:), allocatable :: ymat,zmat
!!$    real(REAL64) :: u, v, alpha, beta
!!$    integer :: i,j,n
!!$
!!$    iflag = 0
!!$    n = size(s)
!!$    allocate(z(n))!,stat=stat,errmsg=emsg)
!!$!    allocate(z(n),zmat(n,n),ymat(n,n))!,stat=stat,errmsg=emsg)
!!$!    if( stat /= 0 ) then
!!$!       iflag = 6
!!$!       return
!!$!    end if
!!$
!!$!    s = x - x_old
!!$!         y = fvec - fvec_old
!!$    z = matmul(b,s)
!!$    do i = 1, n
!!$       do j = i,n
!!$          ymat(i,j) = y(i)*y(j)
!!$          ymat(j,i) = ymat(i,j)
!!$          zmat(i,j) = z(i)*z(j)
!!$          zmat(j,i) = zmat(i,j)
!!$       end do
!!$    end do
!!$
!!$!    write(*,*) 'zmat',real(zmat)
!!$!    write(*,*) 'ymat',real(ymat)
!!$!    write(*,*) sum(s*z),sum(s*y)
!!$
!!$    alpha = sum(s*z)
!!$    beta = sum(s*y)
!!$    do concurrent ( i = 1:n )
!!$       do concurrent ( j = 1:n )
!!$          u = y(i)*y(j)
!!$          v = z(i)*z(j)
!!$          b(i,j) = b(i,j) - v / alpha + u / beta
!!$!         b = b - matmul(z,transpose(z))/sum(s*z) + matmul(y,transpose(y))/sum(y*s)
!!$       end do
!!$    end do
!!$
!!$!    do concurrent (i = 2:n )
!!$!       do concurrent (j = i+1:n )
!!$!          b(j,i) = b(i,j)
!!$!       end do
!!$!    end do
!!$
!!$  end subroutine bfgs


  subroutine locsqp(fcn,cnstr,x,lambda,hess,maxdev,maxfev,quasihess,verbose,info,msg)

    ! local SQP (no region of a feasibility checking)
    ! Algorithm 18.1

    class(sqp_function), intent(in) :: fcn, cnstr


!!$    interface
!!$       subroutine fcn(x,f,fvec,fjac,c,cvec,cjac,iflag)
!!$         use iso_fortran_env
!!$         real(REAL64), dimension(:), intent(in) :: x
!!$         real(REAL64), intent(out) :: f, c
!!$         real(REAL64), dimension(:), intent(out) :: fvec, cvec
!!$         real(REAL64), dimension(:,:), intent(out) :: fjac, cjac
!!$         integer, intent(in out) :: iflag
!!$       end subroutine fcn
!!$    end interface

    real(REAL64), dimension(:), intent(in out) :: x
    real(REAL64), intent(in out) :: lambda
    real(REAL64), dimension(:,:), optional, intent(out) :: hess
    real(REAL64), optional, intent(in) :: maxdev
    integer, optional, intent(in) :: maxfev!, nprint
    logical, optional, intent(in) :: quasihess, verbose
    integer, optional, intent(out) :: info
    character(len=*), optional, intent(out) :: msg

    real(kind), parameter :: dlt = real(0,kind)
    real(kind), parameter :: rho = real(0,kind)
    integer, parameter :: trinfo = 0
    character(len=80) :: emsg
    real(REAL64), dimension(:,:), allocatable :: B, C
    real(REAL64), dimension(:), allocatable :: p,fgrad,cgrad,fgrad_old,x_old,cgrad_old
    real(REAL64) :: tol,f,a,fval,cval
    integer :: n,i,maxev,iter,iflag,stat
    logical :: quasi, verb

    iflag = 666
    n = size(x)

    if( present(maxdev) ) then
       tol = maxdev
    else
       tol = 10*n*epsilon(x)
    end if

    if( present(maxfev) ) then
       maxev = maxfev
    else
       maxev = 3*n*precision(x)
    end if

    if( present(quasihess) ) then
       quasi = quasihess
    else
       quasi = .true.
    end if

!    if( quasi ) then
!       qiflag = 1
!    else
!       qiflag = 2
!    end if

!    iprint = -1
!    if( present(nprint) ) then
!       if( nprint > 0 ) then
!          iprint = nprint
!       end if
!    end if

    allocate(B(n,n),C(n,n),fgrad(n),cgrad(n),p(n),fgrad_old(n), &
         x_old(n),cgrad_old(n),stat=stat,errmsg=emsg)
    if( stat /= 0 ) then
       iflag = 6
       goto 666
    end if

!    lambda = 1
!    nfev = 0
!    lagmin = 0.1*huge(f)

!    if( quasi ) then
!       B = 0
!       forall( i = 1:n ) B(i,i) = 1
!       do concurrent ( i = 1:n )
!          hes(i,i) = 1
       !       end do
!    else
!       call fcn%hessian(x,B)
!       call cnstr%hessian(x,C)
!    end if
!    hes(1,1) = 1e6
!    hes(2,2) = 1e-6

!    iter = 0
    do iter = 1, maxev
!       iter = iter + 1

       ! evaluate functions at x
       fval = fcn%value(x)
       if( fcn%break ) then
          iflag = 6
          exit
       end if
       call fcn%gradient(x,fgrad)
       if( fcn%break ) then
          iflag = 6
          exit
       end if

       ! evaluate constrains at x
       cval = cnstr%value(x)
       if( cnstr%break ) then
          iflag = 6
          exit
       end if
       call cnstr%gradient(x,cgrad)
       if( cnstr%break ) then
          iflag = 6
          exit
       end if

       ! debug prints
       if( debug ) write(*,*) 'iter=',iter,'convergence:',cval,' - ', &
            (fgrad - lambda*cgrad),' = ',f-lambda*cval

       if( verb ) then
          call fcn%print(iter,x,fval,fgrad,B,p,dlt,rho,trinfo)
          call cnstr%print(iter,x,cval,cgrad,C,p,dlt,rho,trinfo)
       end if

       ! check the convergence
       if( maxval(abs(fgrad - cgrad*lambda)) < tol .and. abs(cval) < tol ) then
          iflag = 0
          exit
       end if

!       if( iprint > 0 .and. mod(iter-1,iprint) == 0 ) then
!          iflag = 0
!          call fcn(x,f,fvec,fjac,c,cvec,cjac,iflag)
!          if( iflag < 0  ) goto 666
!       end if

       ! function
!       iflag = qiflag
!       call fcn(x,f,fvec,fjac,c,cvec,cjac,iflag)
!       nfev = nfev + 1
!       if( iflag < 0  ) goto 666
!       write(*,*) qiflag
!       if( nfev > maxev ) then
!          iflag = 2
!          goto 666
!       end if

       ! convergence check
!       lag = f - lambda*c
!       write(*,*) 'lag:',f,c,lag,lagmin
!       if( 10*lag < lagmin ) then
!          lagmin = lag
!       else
!          iflag = 3
!          goto 666
!       end if

!       if( iter == n*precision(x) ) then
!          iflag = 2
!          exit
!       end if

       if( quasi ) then ! quasi Hessian
          if( iter == 1 ) then
             B = 0
             forall( i = 1:n ) B(i,i) = 1
          else
             call bfgs(x,fgrad,cgrad,x_old,fgrad_old,cgrad_old,lambda,B)
          end if
       else ! compute the analytic Hessian with the constrain of second derivatives
          call fcn%hessian(x,B)
          call cnstr%hessian(x,C)
!          B = B !fjac !- lambda*cjac
       end if


       ! solution of the system of equations
       block

!         use KKTsolvers
         real(REAL64), dimension(:,:), allocatable :: mat, L, T
         real(REAL64), dimension(:), allocatable :: r, q
         integer :: m

         m = n + 1
         allocate(mat(m,m),L(m,m),T(m,m),q(m),r(m),stat=stat,errmsg=emsg)
         if( stat /= 0 ) then
            iflag = 6
            goto 666
         end if

         mat(1:n,1:n) = B - lambda*C
         mat(1:n,m) = -cgrad
         mat(m,1:n) = cgrad
         mat(m,m) = real(0,REAL64)
         r(1:n) = -fgrad
         r(m) = -cval
         !         call qrsolve(mat,r,q)
!         call aasen_fac(mat,L,T)
!         call aasen_solve(L,T,r,q)


         do i = 1, m
            write(*,*) real(mat(i,:)),' * x = ',real(r(i))
         end do
         write(*,*) '=',real(q)

!         write(*,*) r - matmul(mat,q)

         p = q(1:n)
         lambda = q(m)
       end block

       ! convergence test
!       diag = [ (1/sqrt(hes(i,i)),i = 1,n) ]
!       diag = 1
!       write(*,*) 'convergence:',c,' - ',(fvec - lambda*cvec),' = ',f-lambda*c
!       if( abs(c) < tol .and. norm2((fvec - lambda*cvec)) < tol ) then
!          iflag = 1
!          exit
!       end if

       ! keep results
       x_old = x
       fgrad_old = fgrad
       cgrad_old = cgrad

       ! update solution
       x = x + p

    end do

!    if( present(cov) ) call qrinv(hes,cov)

666 continue

    if( present(info) ) info = iflag
    if( present(msg) ) then
       msg = '.. undefined ...'
       if( iflag == 1 ) &
            msg = 'An uncertainity of the solution is under the required precision.'
       if( iflag == -1 ) msg = 'An interrupt by user.'
       if( iflag == 2 ) msg = 'Number of calls to fcn reached maxfev.'
       if( iflag == 3 ) msg = 'Some function or constrains diverges.'
       if( iflag == 6 ) msg = emsg
    end if
    if( present(hess) ) hess = B

  end subroutine locsqp


  subroutine linesqp(fcn,x,cov,maxdev,maxfev,nprint,info,msg)

! linear search SQP

!  This implementation thightly follows the algorithm 18.3 (Line search SQL)
!  in the second edition of Numerical Optimization by Jorge Nocedal and
!  Stephen J. Wright.
!

    interface
       subroutine fcn(x,f,fvec,fjac,c,cvec,cjac,iflag)
         use iso_fortran_env
         real(REAL64), dimension(:), intent(in) :: x
         real(REAL64), intent(out) :: f, c
         real(REAL64), dimension(:), intent(out) :: fvec, cvec
         real(REAL64), dimension(:,:), intent(out) :: fjac, cjac
         integer, intent(in out) :: iflag
       end subroutine fcn
!!$       subroutine con(x,c,cvec,iflag)
!!$         real(REAL64), dimension(:), intent(in) :: x
!!$         real(REAL64), intent(out) :: c
!!$         real(REAL64), dimension(:), intent(out) :: cvec
!!$         integer, intent(in out) :: iflag
!!$       end subroutine con
    end interface

    real(REAL64), parameter :: eta = 0.4
    real(REAL64), parameter :: tau = 0.5

    real(REAL64), dimension(:), intent(in out) :: x
    real(REAL64), dimension(:,:), optional, intent(out) :: cov
    real(REAL64), optional, intent(in) :: maxdev
    integer, optional, intent(in) :: maxfev, nprint
    integer, optional, intent(out) :: info
    character(len=*), optional, intent(out) :: msg

    character(len=80) :: emsg
    real(REAL64), dimension(:,:), allocatable :: b, fjac, cjac
    real(REAL64), dimension(:), allocatable :: p, fvec, cvec,fvec_old,x_old, &
         cvec_old, lagrange, lagrange_old
    real(REAL64) :: tol,f,c,lambda,a,mu,f_old,c_old,lag,lag_old
    integer :: n,i,maxev,nfev,iter,iflag,stat

    iflag = 666

    if( present(maxdev) ) then
       tol = maxdev
    else
       tol = sqrt(epsilon(x))
    end if

    if( present(maxfev) ) then
       maxev = maxfev
    else
       maxev = 666
    end if

    n = size(x)
    allocate(b(n,n),fjac(n,n),cjac(n,n),fvec(n),cvec(n),p(n),fvec_old(n),x_old(n),&
         cvec_old(n),lagrange(n),lagrange_old(n),stat=stat,errmsg=emsg)
    if( stat /= 0 ) then
       iflag = 6
       goto 666
    end if

    mu = -1
    lambda = 1
    call fcn(x,f,fvec,fjac,c,cvec,cjac,iflag)
    if( iflag < 0 ) goto 666
    nfev = 1
    if( abs(c) > 0 ) then
       lambda = f / c
    else
       lambda = 1
    end if
    x_old = x
    fvec_old = fvec
    cvec_old = cvec
    lagrange = fvec - lambda*cvec
    lagrange_old = lagrange
    f_old = f
    c_old = c
    lag = f - lambda*c
    lag_old = lag
!    call con(x,c,cvec,iflag)
!    if( iflag < 0 ) goto 666

    if( quasi ) then
       b = 0
       do concurrent ( i = 1:n )
          b(i,i) = 1
       end do
    else
       b = fjac - lambda*cjac
    end if

    do iter = 1, 3*precision(x)
       write(*,*) 'iter=',iter

       ! convergence test
       write(*,*) c,' - ',fvec - lambda*cvec,' = ',f
       if( abs(c) < tol .and. all(abs(fvec - lambda*cvec) < tol) ) then
          iflag = 1
          exit
       end if


       ! solution of the system of equations
!       call linsolve(b,fvec,c,cvec,p,lambda)
       block
         real(REAL64), dimension(:,:), allocatable :: mat
         real(REAL64), dimension(:), allocatable :: r, q
         integer :: m

         m = n + 1
         allocate(mat(m,m),q(m),r(m),stat=stat,errmsg=emsg)
         if( stat /= 0 ) then
            iflag = 6
            goto 666
         end if

         mat(1:n,1:n) = b
         mat(1:n,m) = -cvec
         mat(m,1:n) = cvec
         mat(m,m) = real(0,REAL64)
         r(1:n) = -fvec !+ lambda*cvec
         r(m) = -c
         call qrsolve(mat,r,q)

         write(*,*) '>',real(r)
         write(*,*) real(mat(1,:))
         write(*,*) real(mat(2,:))
         write(*,*) real(mat(3,:))
         write(*,*) '=',real(q)

         p = q(1:n)
         lambda = q(m)
       end block
!       lag = fvec - lambda*cvec

!       goto 33
       ! line search
       block
         real(REAL64) :: nu, q, q2, mrt, mrt0, dir, c1, f1
         real(REAL64), dimension(2) :: fvec1, cvec1
         real(REAL64), dimension(2,2) :: fjac1, cjac1
         integer :: miter

         if( mu < 0 ) mu = 2*abs(lambda)

         !       merit(x,mu) = f(x) +mu*abs(c)

         if( abs(c) > 0 ) then
            q = (1 - 0.5)*abs(c)
            q2 = max(sum(p*matmul(b,p)),real(0.0,REAL64))
            !         q2 = sum(p*matmul(b,p))
            !         q2 = 0
            nu = (sum(fvec*p) + 0.5*q2) / q
!         nu = sum(fvec*p) / q
            if( nu > mu ) mu = nu
         end if

!         minmer = f - lambda*c  + mu*abs(c)
         mrt0 = f + mu*abs(c)
         dir = sum(fvec*p) - mu*abs(c)
         if( .not. (dir < 0) ) exit
         write(*,*) 'mu, dir,',real(mu),real(nu),real(dir),real(q2),real(sum(fvec*p))
         a = 1
         do miter = 1, 8
            call fcn(x+a*p,f1,fvec1,fjac1,c1,cvec1,cjac1,iflag)
            nfev = nfev + 1
            if( iflag < 0  ) goto 666
            if( nfev > maxev ) then
               info = 2
               goto 666
            end if

!            mrit = f1 - lambda*c1  + mu*abs(c1)
            mrt = f1  + mu*abs(c1)
            write(*,*) 'merit:',real(a), real(mrt),real(mrt0 + eta*a*dir)
            if( mrt < mrt0 + eta*a*dir ) exit
            a = tau*a
         end do
       end block
!33     continue
!       a = 1

       ! update solution + function values, their gradients and jacobians
       x = x + a*p
       call fcn(x,f,fvec,fjac,c,cvec,cjac,iflag)
       nfev = nfev + 1
       if( iflag < 0  ) goto 666
       if( nfev > maxev ) then
          info = 2
          goto 666
       end if

       lagrange = fvec - lambda*cvec
       lag = f - lambda*c

       ! check for divergence
       write(*,*) lag,lag_old,sum(p*matmul(b,p))
!       if( lag > lag_old + max(sum(p*matmul(b,p)),0.01*abs(lag_old)) &
!            .or. abs(c) > 10*abs(c_old)  ) then
!       if( f > f_old + max(sum(p*matmul(b,p)),0.01*abs(f_old)) .or. abs(c) > 10*abs(c_old)  ) then
!          iflag = 3
!          exit
!       end if


!       write(*,*) fvec,fvec_old
!       write(*,*) cvec,cvec_old
!       call con(x,c,cvec,iflag)
!       if( iflag < 0 ) goto 666

       if( quasi ) then
          ! update by the quasi-Hessian BFGS method
!          call bfgs(x - x_old,(fvec - fvec_old) - lambda*(cvec - cvec_old),b,iflag)
          block
!            real(REAL64), dimension(size(x)) :: s, u

!            s = x - x_old
!            u = (fvec - fvec_old) - lambda*(cvec - cvec_old)
!            call bfgs(s,u,b)
            call bfgs(x,fvec,cvec,x_old,fvec_old,cvec_old,lambda,b)
          end block
       else ! compute the analytic Hessian with the constrain of second derivatives
          b = fjac - lambda*cjac
       end if

       ! keep values
       x_old = x
       fvec_old = fvec
       cvec_old = cvec
       lagrange_old = lagrange
       f_old = f
       c_old = c
       lag_old = lag

    end do

    if( present(cov) ) call qrinv(b,cov)

666 continue

    if( present(info) ) info = iflag
    if( present(msg) ) then
       msg = '.. undefined ...'
       if( iflag == 1 ) &
            msg = 'An uncertainity of the solution is under the required precision.'
       if( iflag == -1 ) msg = 'An interrupt by user.'
       if( iflag == 2 ) msg = 'Number of calls to fcn reached maxfev.'
       if( iflag == 3 ) msg = 'Any functions or constrains diverges.'
       if( iflag == 6 ) msg = emsg
    end if

  end subroutine linesqp


!!$  function que(p,f,fvec,c,cvec,hes,mu)
!!$
!!$    real(REAL64) :: que
!!$    real(REAL64), intent(in) :: f, mu, c
!!$    real(REAL64), dimension(:), intent(in) :: p, fvec, cvec
!!$    real(REAL64), dimension(:,:), intent(in) :: hes
!!$
!!$    que = f + sum(fvec*p) + sum(p*matmul(hes,p))/2 + mu*norm2(c + cvec*p)
!!$
!!$  end function que


!  function enorm(x)

!    real(REAL64) :: enorm
!    real(REAL64), dimension(:), intent(in) :: x
!    real(REAL64) :: s2

!    s2 = sum(x**2)
!    s2 = hypot(x)
!    if( s2 > 0 ) then
!       enorm = sqrt(s2)
!    else
!       enorm = real(0,REAL64)
!    end if

!  end function enorm

!  function norm1(x)

!    norm1 = sum(abs(x))

!  end function norm1

!  function inftynorm(x)

!    inftynorm = maxval(abs(x))

!  end function inftynorm


end module sqps
