!
!  QR factorisation
!
!  Copyright © 2021-2026 F.Hroch (hroch@physics.muni.cz)
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
! References: see doc/guide/oakleaf.pdf, Chapter: QR factorisation

module QRHouse

  use iso_fortran_env
  implicit none
  private

  integer, parameter, private :: kind = REAL64
  logical, parameter, private :: debug = .false.

  public :: qrfac, qform, rform, qrpivot,  qrsol, qrinv, qreig

contains

  subroutine qrfac(A,beta)

    real(kind), dimension(:,:), intent(in out) :: A
    real(kind), dimension(:), intent(out) :: beta
    real(kind), dimension(size(A,1)) :: v

    integer :: i,j,n,m

    m = size(A,1)
    n = size(A,2)

    do j = 1, n
       call housevec(A(j:m,j),v(j:m),beta(j))
       call housemul(v(j:m),beta(j),A(j:m,j:n))
       if( j < m ) &
            A(j+1:m,j) = v(j+1:m)

       if( debug ) then
          write(*,*) 'housevec, beta:',v(j:m),beta(j)
          do i = 1,n
             write(*,*) 'A:',A(i,:)
          end do
       end if

    end do

  end subroutine qrfac


  pure subroutine housevec(x,v,beta)

    real(kind), dimension(:), intent(in) :: x
    real(kind), dimension(:), intent(out) :: v
    real(kind), intent(out) :: beta

    integer :: m
    real(kind) :: rnorm2, xnorm

    m = size(x)

    v = [real(1,kind), x(2:m)]
    rnorm2 = dot_product(x(2:m),x(2:m))
    if( abs(rnorm2) > epsilon(rnorm2) ) then
       xnorm = sqrt(x(1)**2 + rnorm2)
       if( x(1) > 0 ) then
          v(1) = -rnorm2/(x(1) + xnorm)
       else
          v(1) = x(1) - xnorm
       end if
       v = v / v(1)
       beta = real(2,kind) / dot_product(v,v)
    else ! if( rnorm == 0 ) then
       if( x(1) >= 0 ) then
          beta = real(0,kind)
       else ! if( x(1) < 0 ) then
          beta = real(2,kind)
       end if
    end if

  end subroutine housevec


  pure subroutine housemul(v,beta,A)

    real(kind), dimension(:), intent(in) :: v
    real(kind), intent(in) :: beta
    real(kind), dimension(:,:), intent(in out) :: A

    real(kind), dimension(size(v)) :: u, w
    integer :: i,j,m,n

    m = size(A,1)
    n = size(A,2)
    u = beta*v
    forall( j = 1:n ) w(j) = dot_product(v,A(:,j))
    forall( i = 1:m, j = 1:n ) A(i,j) = A(i,j) - u(i)*w(j)

  end subroutine housemul

  pure subroutine houseapp(R,beta,b,u)

    real(kind), dimension(:,:), intent(in) :: R
    real(kind), dimension(:), intent(in) :: beta
    real(kind), dimension(:), intent(in) :: b
    real(kind), dimension(:), intent(out) :: u

    real(kind), dimension(size(R,1)) :: v
    integer :: j,m,n

    m = size(R,1)
    n = size(R,2)

    u = b
    do j = 1, n
       v(j:m) = [ real(1,kind), R(j+1:m,j) ]
       u(j:m) = u(j:m) - beta(j)*dot_product(v(j:m),u(j:m))*v(j:m)
    end do

  end subroutine houseapp


  ! Qform
  pure subroutine qform(A,beta,Q)

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), dimension(:), intent(in) :: beta
    real(kind), dimension(:,:), intent(out) :: Q

    real(kind), dimension(size(A,1)) :: v
    integer :: j,n,m

    m = size(A,1)
    n = size(A,2)

    Q = real(0,kind)
    forall( j = 1:n ) Q(j,j) = real(1,kind)
    do j = n,1,-1
       v(j:m) = [ real(1,kind), A(j+1:m,j) ]
       call housemul(v(j:m),beta(j),Q(j:m,j:n))
    end do

  end subroutine qform


  ! Rform
  pure subroutine rform(A,R)

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), dimension(:,:), intent(out) :: R

    integer :: i,j,m,n

    m = size(A,1)
    n = size(A,2)

    forall( i = 1:m, j = 1:n, i > j )  R(i,j) = real(0,kind)
    forall( i = 1:m, j = 1:n, i <= j ) R(i,j) = A(i,j)

  end subroutine rform


  pure subroutine rsol(R,qtb,x)

    ! R x = Q**T b, R is a upper triangular matrix

    real(kind), dimension(:,:), intent(in) :: R
    real(kind), dimension(:), intent(in) :: qtb
    real(kind), dimension(:), intent(out) :: x

    integer :: m,n,i,i1

    m = size(R,1)
    n = size(R,2)

    x(m) = qtb(m) / R(m,m)
    do i = m-1,1,-1
       i1 = i + 1
       x(i) = (qtb(i) - dot_product(R(i,i1:n),x(i1:n))) / R(i,i)
    end do

  end subroutine rsol

  pure subroutine lsol(L,qtb,x)

    ! L x = Q**T b, L is a lower triangular matrix

    real(kind), dimension(:,:), intent(in) :: L
    real(kind), dimension(:), intent(in) :: qtb
    real(kind), dimension(:), intent(out) :: x

    integer :: m,i,i1

    m = size(L,1)

    x(1) = qtb(1) / L(1,1)
    do i = 2,m
       i1 = i - 1
       x(i) = (qtb(i) - dot_product(L(1:i1,i),x(1:i1))) / L(i,i)
    end do

  end subroutine lsol


  subroutine qrpivot(A,beta,pivot)

    real(kind), dimension(:,:), intent(in out) :: A
    real(kind), dimension(:), intent(out) :: beta
    integer, dimension(:), intent(out) :: pivot

    real(kind), dimension(size(A,2)) :: v,c
    real(kind), dimension(size(A,1)) :: tmp
    integer :: k,l,n,m,r,r1
    real(kind) :: temp

    m = size(A,1)
    n = size(A,2)
    pivot = [ (r,r=1,n) ]
    forall( r = 1:n ) c(r) = dot_product(A(1:m,r),A(1:m,r))

    if( debug ) then
       do r = 1, m
          write(*,*) 'A(i,:):',A(r,:)
       end do
       write(*,*) pivot
       write(*,*) 'init c:',c
    end if

    r = 0
    do r = 1, n
       r1 = r + 1

       k = (r - 1) + maxloc(c(r:n),1)
       if( .not. (c(k) > 0) ) exit

       if( debug ) write(*,*) r,'<->',k,c(r:n),c(k)

       ! swap cols
       if( k /= r ) then
          l = pivot(r); pivot(r) = pivot(k); pivot(k) = l
          tmp(1:m) = A(1:m,r); A(1:m,r) = A(1:m,k); A(1:m,k) = tmp(1:m)
          temp = c(r); c(r) = c(k); c(k) = temp
       end if

       if( debug ) write(*,*) 'k,c, pivot:',k,c,pivot(1:r)

       call housevec(A(r:m,r),v(r:m),beta(r))
       call housemul(v(r:m),beta(r),A(r:m,r:n))
       A(r1:m,r) = v(r1:m)

       c(r1:n) = c(r1:n) - A(r,r1:n)**2

       if( debug ) then
          write(*,*) 'v:',v
          do l = 1, n
             write(*,*) 'A(i,:)',A(l,:)
          end do
          write(*,*) 'r,c:',r,c
       end if

    end do

    if( debug ) then
       do l = 1, m
          write(*,*) 'A(i,:)',A(l,:)
       end do
       write(*,*) 'pivot:',pivot
       write(*,*) 'beta:',beta
    end if

  end subroutine qrpivot

  pure function mrank(A,tol) result(r)

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), intent(in) :: tol
    integer :: i, r

    r = 0
    do i = 1, size(A,1)
       if( abs(A(i,i)) < tol ) exit
       r = i
    end do

  end function mrank

  subroutine qrsol(A,b,x,rtol)

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), dimension(:), intent(in) :: b
    real(kind), dimension(:), intent(out) :: x
    real(kind), intent(in), optional :: rtol

    real(kind), dimension(size(A,1),size(A,2)) :: R
    real(kind), dimension(size(A,2),size(A,2)) :: U
    integer, dimension(size(A,2)) :: pivot
    real(kind), dimension(size(b)) :: beta, qtb, z
    real(kind) :: tol
    integer :: i,m,n

    m = size(A,1)
    n = size(A,2)

    if( m > n ) then
       ! the least-square solution
       R = A
       call qrfac(R,beta)
       call houseapp(R,beta,b,qtb)
       call rsol(R(1:n,1:n),qtb,x)

    else

       R = A
       call qrpivot(R,beta,pivot)
       call houseapp(R,beta,b,qtb)

       if( present(rtol) ) then
          tol = rtol
       else
          tol = 1e2*epsilon(R)*abs(R(1,1))
       end if

       m = mrank(R,tol)

       if( debug ) write(*,*) 'rank=',m,'pivot:', pivot

       if( m == size(A,1) ) then
          ! a square matrix, full rank

          call rsol(R,qtb,z)
          x(pivot) = z

       else if( m < size(A,1) ) then
          ! a rank-deficient matrix

          if( debug ) then
             call rsol(R(1:m,1:m),qtb(1:m),z(1:m))
             z(m+1:) = 0
             z = z(pivot)
             write(*,*) 'the basic solution:',z
          end if

          ! R = R**T, nullify the upper trapezoid
          forall( i = 1:m ) R(i:n,i) = R(i,i:n)
          forall( i = 1:m ) R(1:i-1,i) = real(0,kind)

          if( debug ) then
             do i = 1, n
                write(*,*) 'R**T',R(i,1:m)
             end do
          end if

          ! the second factoristaion
          call qrfac(R(1:n,1:m),beta(1:m))
          call qform(R(1:n,1:m),beta(1:m),U(1:n,1:n))

          if( debug ) then
             write(*,*) m
             do i = 1, n
                write(*,*) 'R(i,:):',R(i,1:m)
             end do
             do i = 1, n
                write(*,*) 'U(i,:):',U(i,1:m)
             end do
          end if

          call lsol(R(1:m,1:m),qtb(1:m),z(1:m))

          if( debug ) then
             write(*,*) 'z:',z(1:m)
             write(*,*) matmul(U(:,1:m),z(1:m))
          end if

          z(1:m) = matmul(U(:,1:m),z(1:m))
          x(pivot) = z

       end if

    end if

  end subroutine qrsol

  subroutine qrinv(A,X,rtol)

    ! This routine computes of inverse of a squres matrix.
    ! All computations are done in least-square sense for
    ! possibly rank deficient matrixes.

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), dimension(:,:), intent(out) :: X
    real(kind), intent(in), optional :: rtol

    real(kind), dimension(size(A,1),size(A,2)) :: R
    real(kind), dimension(size(A,2),size(A,2)) :: Q,U
    real(kind), dimension(size(A,1)) :: z,beta
    integer, dimension(size(A,2)) :: pivot
    real(kind) :: tol
    integer :: i,n,m

    m = size(A,1)
    n = size(A,2)

    R = A
    call qrpivot(R,beta,pivot)
    call qform(R,beta,Q)

    if( present(rtol) ) then
       tol = rtol
    else
       tol = 1e2*epsilon(R)*abs(R(1,1))
    end if

    m = mrank(R,tol)

    if( debug ) write(*,*) 'rank=',m,'pivot:', pivot

    if( m == size(A,1) ) then
       ! a full rank matrix

       do i = 1, m
          call rsol(R,Q(i,:),z)
          X(pivot,i) = z
       end do

    else if( m < size(A,1) ) then

       ! a rank deficient matrix

       ! transpose R = R**T, also nullify the upper trapezoid
       forall( i = 1:m ) R(i:n,i) = R(i,i:n)
       forall( i = 1:m ) R(1:i-1,i) = real(0,kind)

       ! the second factorisation
       call qrfac(R(1:n,1:m),beta(1:m))
       call qform(R(1:n,1:m),beta(1:m),U(1:n,1:n))

       if( debug ) then
          do i = 1, n
             write(*,*) 'R(i,:)',R(i,:)
          end do
          do i = 1, n
             write(*,*) 'U(i,:)',U(i,:)
          end do
       end if

       do i = 1, n
          call lsol(R(1:m,1:m),Q(i,1:m),z(1:m))
          z = matmul(U(:,1:n),z)
          X(pivot,i) = z
       end do

    end if

  end subroutine qrinv


  subroutine qreig(D,a,Q)

    ! We're expecting a symmetric square matrix; that's garantied
    ! both real eigen values and vectors.

    real(kind), dimension(:,:), intent(in)  :: D
    real(kind), dimension(:), intent(out)   :: a
    real(kind), dimension(:,:), intent(out) :: Q

    real(kind), dimension(size(a)-1) :: b,tol
    integer :: i,j,k,l,m,n,iter

    n = size(D,1)

    if( size(D,2) /= n ) &
         stop 'Error: qreig solves only eigenvalue problems with a square matrix.'

    ! TEST for symmetry A(i,j) == A(j,i)
    do i = 1, n
       do j = i+1,n
          if( abs(D(i,j) - D(j,i)) > 10*epsilon(D) ) &
               stop 'Error: qreig solves only symmetric eigenvalue problems.'
       end do
    end do

    Q = D
    call housetridig(Q,a,b)

    if( debug ) then
       block
         real(kind), dimension(size(D,1),size(D,2)) :: D1, T

         do k = 1, n
            write(*,*) 'D:',D(k,:)
         end do
         do k = 1, n
            write(*,*) 'Q:',Q(k,:)
         end do
         write(*,*) 'a:',a
         write(*,*) 'b:',b
         T = 0
         forall( k = 1:n ) T(k,k) = a(k)
         forall( k = 1:n-1 )
            T(k,k+1) = b(k)
            T(k+1,k) = b(k)
         end forall

         D1 = matmul(Q,matmul(T,transpose(Q)))
         write(*,*) all(abs(D1 - D) < 1e3*epsilon(D)), maxval(abs(D1 - D))

         if( .not. all(abs(D1 - D) < 1e3*epsilon(D)) ) then
            do k = 1, n
               write(*,*) 'D1:',D1(k,:)
            end do
         end if
       end block
    end if

    l = 1
    m = n
    do iter = 1, 10*n*precision(D)

       forall( i = 1:n-1 ) tol(i) = 2*epsilon(D)*(abs(a(i)) + abs(a(i+1)))
       do i = 2, n-1
          if( abs(b(i)) > tol(i) .and. abs(b(i-1)) < tol(i-1) ) l = i
          if( abs(b(i)) < tol(i) .and. abs(b(i-1)) > tol(i-1) ) m = i
       end do

       if( debug ) then
          write(*,*) 'iter=',iter,'    interval:',l,m
          do i = 1, n-1
             write(*,*) tol(i),abs(b(i)),i
          end do
       end if

       if( m > l ) &
            call qrstep(a(l:m),b(l:m-1),Q(:,l:m))

       if( m-l == 1 ) exit

       if( debug ) then
          do k = 1, n
             write(*,*) 'Q:',Q(k,:)
          end do
       end if

    end do

    if( debug ) then
       do i = 1, n
          write(*,*) 'D:',D(i,:)
       end do
    end if

  end subroutine qreig

  subroutine housetridig(T,a,b)
    ! Householder tridiagonalisation

    real(kind), dimension(:,:), intent(in out) :: T
    real(kind), dimension(:), intent(out) :: a,b

    real(kind), dimension(size(T,1)) :: p,v,w,u,beta
    real(kind) :: r
    integer :: i,j,k,n,k1

    n = size(T,1)

    do k = 1, n-2
       k1 = k + 1
       call housevec(T(k1:n,k),v(k1:n),beta(k))
       p(k1:n) = beta(k) * matmul(T(k1:n,k1:n),v(k1:n))
       r = beta(k) * dot_product(p(k1:n),v(k1:n)) / real(2,kind)
       w(k1:n) = p(k1:n) - r*v(k1:n)
       T(k1,k) = norm2(T(k1:n,k))
       T(k,k1) = T(k1,k)
       forall( i = k1:n, j = k1:n )
          T(i,j) = T(i,j) - v(i)*w(j) - w(i)*v(j)
       end forall

       if( k > 1 ) T(k1:n,k-1) = u(k1:n)
       u(k1:n) = v(k1:n)

    end do
    if( n > 2 ) T(n,n-2) = u(n)

    if( debug ) then
       write(*,*) 'housetridig():'
       do i = 1, n
          write(*,*) 'T:',T(i,:)
       end do
    end if

    forall( i = 1:n ) a(i) = T(i,i)
    forall( i = 1:n-1 ) b(i) = T(i,i+1)

    T(n-1,n-1:n) = real([1, 0],kind)
    T(n  ,n-1:n) = real([0, 1],kind)

    do k = n-2,1,-1
       k1 = k + 1
       v(k1:n) = [ real(1,kind), T(k+2:n,k) ]

       T(k,k) = real(1,kind)
       T(k1:n,k) = real(0,kind)
       T(k,k1:n) = real(0,kind)

       forall( i = k1:n ) w(i) = dot_product(T(k1:n,i),v(k1:n))
       u(k1:n) = beta(k)*v(k1:n)
       forall( i = k1:n, j = k1:n ) T(i,j) = T(i,j) - w(j)*u(i)

       if( debug ) then
          do i = 1, n
             write(*,*) T(i,:)
          end do
       end if

    end do

  end subroutine housetridig


  subroutine qrstep(a,b,Q)

    real(kind), dimension(:), intent(in out) :: a,b
    real(kind), dimension(:,:), intent(in out) :: Q

    real(kind), dimension(2,2) :: G
    real(kind) :: d, mu, x, y, c, s, c2, s2, cs, u, v, w
    integer :: k,l,n,k1

    n = size(a,1)

    d = (a(n-1) - a(n)) / real(2,kind)
    mu = a(n) - b(n-1)**2 / (d + sign(norm2([d,b(n-1)]),d))
    x = a(1) - mu
    y = b(1)

    if( debug ) write(*,*) 'x,z,mu',x,y,mu

    do k = 1, n-1
       call givens(x,y,c,s)

       c2 = c**2
       s2 = s**2
       cs = c*s
       u = a(k)
       v = a(k+1)
       w = real(2,kind)*cs*b(k)
       a(k)   = c2*u + s2*v - w
       a(k+1) = c2*v + s2*u + w
       b(k)   = (c2 - s2)*b(k) + cs*(u - v)

       if( k > 1 ) b(k-1) = c*x - s*y
       if( k < n - 1 ) then
          x = b(k)
          y = -s*b(k+1)
          b(k+1) = c*b(k+1)
       end if

       k1 = k + 1
       G(1,1:2) = [ c, s ]
       G(2,1:2) = [-s, c ]
       Q(:,k:k1) = matmul(Q(:,k:k1),G)

       if( debug ) then
          write(*,*) 'x,y:',x,y
          write(*,*) 'a:',a
          write(*,*) 'b:',b
          do l = 1, 2
             write(*,*) 'G:',G(l,:)
          end do
       end if

    end do

  end subroutine qrstep


  pure subroutine givens(x,y,c,s)

    real(kind), intent(in) :: x,y
    real(kind), intent(out) :: c,s
    real(kind) :: t

    if( abs(y) < epsilon(y) ) then
       c = real(1,kind)
       s = real(0,kind)
    else
       if( abs(y) > abs(x) ) then
          t = -x / y
          s = real(1,kind) / hypot(real(1,kind),t)
          c = s * t
       else
          t = -y / x
          c = real(1,kind) / hypot(real(1,kind),t)
          s = c * t
       end if
    end if

  end subroutine givens

end module QRHouse
