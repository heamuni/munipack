!
!  QR factorisation, the debug version
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
!
!
! References: see doc/guide/oakleaf.pdf, Chapter: QR factorisation

module qrdebug

  use iso_fortran_env
  implicit none

  integer, parameter, private :: kind = REAL64

contains

  subroutine qrfac(A,Q,R)

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), dimension(:,:), intent(out) :: Q,R

    real(kind), dimension(size(A,1)) :: v
    real(kind), dimension(size(A,1),size(A,1)) :: P
    integer :: j,k,n,m
    real(kind) :: beta

    m = size(A,1)
    n = size(A,2)

    ! init the matrixes
    R = A
    Q = real(0,kind)
    forall( j = 1:m ) Q(j,j) = real(1,kind)

    do j = 1, n

       call housevec(R(:,j),j,v,beta)
       write(*,*) 'house',v(j:m),beta

       ! form Householder matrix
       call housemat(v,beta,P)

       do k = 1,m
          write(*,*) 'P:',P(k,:)
       end do

       ! compute next iteration of R and Q
!       R = matmul(P,R)
       R(j:m,j:n) = matmul(P(j:m,j:m),R(j:m,j:n))

       do k = 1,m
          write(*,*) 'R:',R(k,:)
       end do

       Q = matmul(Q,P)

    end do

  end subroutine qrfac


  subroutine housevec(x,j,v,beta)

    real(kind), dimension(:), intent(in) :: x
    integer, intent(in) :: j
    real(kind), dimension(:), intent(out) :: v
    real(kind), intent(out) :: beta

    integer :: m, j1
    real(kind) :: rnorm2, xnorm, xnorm2, vnorm2

    m = size(x)
    v(:j-1) = real(0,kind)

    j1 = j + 1
    rnorm2 = dot_product(x(j1:m),x(j1:m))
    xnorm2 = x(j)**2 + rnorm2
    if( xnorm2 > 0 ) then
       xnorm = sqrt(xnorm2)

       if( x(j) > 0 ) then
          v(j) = -rnorm2/(x(j) + xnorm)
       else
          v(j) = x(j) - xnorm
       end if
       v(j1:m) = x(j1:m)
       vnorm2 = v(j)**2 + rnorm2
       if( vnorm2 > epsilon(vnorm2) ) then
          beta = real(2,kind) *v(j)**2 / vnorm2
          v = v / v(j)
       else
          beta = real(0,kind)
       end if
    else
       beta = real(0,kind)
       v = real(0,kind)
    end if

  end subroutine housevec

  subroutine housemat(v,beta,H)

    real(kind), dimension(:), intent(in) :: v
    real(kind), intent(in) :: beta
    real(kind), dimension(:,:), intent(out) :: H
    integer :: n,m,k,l

    m = size(H,1)
    n = size(H,2)

    forall( k = 1:m, l = 1:n ) H(k,l) = - beta * v(k) * v(l)
    forall( k = 1:n ) H(k,k) = real(1,kind) + H(k,k)

  end subroutine housemat

  subroutine rsol(R,qtb,x)

    real(kind), dimension(:,:), intent(in) :: R
    real(kind), dimension(:), intent(in) :: qtb
    real(kind), dimension(:), intent(out) :: x

    integer :: m,i,i1

    m = size(R,1)

    do i = m,1,-1
       i1 = i + 1
       if( abs(R(i,i)) > epsilon(R) ) then
          x(i) = (qtb(i) - dot_product(R(i,i1:m),x(i1:m))) / R(i,i)
       else
          x(i) = real(1,kind)
       end if
    end do

  end subroutine rsol

  subroutine qrlssol(A,b,xls)

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), dimension(:), intent(in) :: b
    real(kind), dimension(:), intent(out) :: xls

    real(kind), dimension(size(A,1),size(A,1)) :: Q
    real(kind), dimension(size(A,1),size(A,2)) :: R
    real(kind), dimension(size(b)) :: qtb, res
    integer :: i,m,n

    m = size(A,1)
    n = size(A,2)

    call qrfac(A,Q,R)

    do i = 1, m
       write(*,*) 'Q(i,1:n)',Q(i,1:n)
    end do

    do i = 1, n
       write(*,*) 'R(i,:)',R(i,:)
    end do

    qtb = matmul(transpose(Q),b)

    call rsol(R(1:n,1:n),qtb(1:n),xls)

    res = matmul(A,xls) - b
    write(*,*) 'norms:',norm2(res), norm2(qtb(n+1:m))

  end subroutine qrlssol


  subroutine qrsol(A,b,x)

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), dimension(:), intent(in) :: b
    real(kind), dimension(:), intent(out) :: x

    real(kind), dimension(size(A,1),size(A,2)) :: R, Q
    real(kind), dimension(size(b)) :: y, qtb
    integer, dimension(size(b)) :: pivot

    call qrpivot(A,Q,R,pivot)

    qtb = matmul(transpose(Q),b)
    write(*,*) 'Q^T*B:', real(qtb)

    call rsol(R,qtb,y)
    write(*,*) 'y:',y
    x = y(pivot)

  end subroutine qrsol

  subroutine qrpivot(A,Q,R,pivot)

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), dimension(:,:), intent(out) :: Q,R
    integer, dimension(:), intent(out) :: pivot

    real(kind), dimension(size(A,2)) :: v,c
    real(kind), dimension(size(A,1)) :: tmp
    real(kind), dimension(size(A,1),size(A,2)) :: H
    integer :: i,j,k,l,n,m
    real(kind) :: tau, beta, temp

    m = size(a,1)
    n = size(a,2)

    R = A
    Q = 0
    forall( j = 1:m ) Q(j,j) = 1
    forall( j = 1:m ) pivot(j) = j

    do j = 1, n
       c(j) = dot_product(A(1:m,j),A(1:m,j))
    end do
    write(*,*) c
    tau = maxval(c)

    do j = 1, n

       if( .not. (tau > 0) ) exit
       k = (j - 1) + maxloc(c(j:n),1)


       ! switch the rows
       l = pivot(j); pivot(j) = pivot(k); pivot(k) = l
       tmp(1:m) = R(1:m,j); R(1:m,j) = R(1:m,k); R(1:m,k) = tmp(1:m)
       temp = c(j); c(j) = c(k); c(k) = temp
       write(*,*) k,c
       write(*,*) pivot

       call housevec(R(:,j),j,v,beta)
       call housemat(v,beta,H)

       write(*,*) v
       do i = 1, n
          write(*,*) 'H(i,:)',H(i,:)
       end do
       R = matmul(H,R)
       Q = matmul(Q,H)

       do i = 1, n
          write(*,*) 'R(i,:)',R(i,:)
       end do

       c(:) = c(:) - R(j,:)**2
       write(*,*) j,real(c)
       tau = maxval(c(j+1:n))
    end do

  end subroutine qrpivot

  subroutine qrinv(A,X)

    ! It computes the inverse matrix in least-square sense
    ! (possibly rank deficient)  by the use of the QR factorization
    ! Any efficiency does not matter. The input should be a square matrix.
    !
    ! http://en.wikipedia.org/wiki/QR_decomposition
    !     - see: "Solution of inverse problems"

    real(kind), dimension(:,:), intent(in) :: A
    real(kind), dimension(:,:), intent(out) :: X

    real(kind), dimension(size(A,1),size(A,2)) :: Q, R
    real(kind), dimension(size(A,2)) :: u
    integer, dimension(size(A,2)) :: pivot
    integer :: i, m

    m = size(A,1)

    call qrpivot(A,Q,R,pivot)

    do i = 1, m
       write(*,*) 'R:',R(i,:)
    end do
    do i = 1, m
       write(*,*) 'Q:',Q(i,:)
    end do

    do i = 1, m

       ! A shortcut for unit vectors: b = e(j):
       ! b = 0
       ! b(ipvt(i)) = 1
       ! qtb = matmul(transpose(q),b)

       call rsol(R,Q(pivot(i),:),u)
       X(:,pivot(i)) = u(pivot)
    end do

  end subroutine qrinv

  subroutine qreig(A,e,Z)

    ! We're expecting a symmetric square matrix; that's garantied
    ! both real eigen values and vectors.

    real(kind), dimension(:,:), intent(in)  :: A
    real(kind), dimension(:), intent(out)   :: e
    real(kind), dimension(:,:), intent(out) :: Z

    real(kind), dimension(size(A,1),size(A,2)) :: Q,R,L
    real(kind) :: d, tol, rnorm
    integer :: i,j,j1,n,iter

    n = size(A,1)
    tol = 1e-7

    ! TEST for symmetry A(i,j) == A(j,i)
    do i = 1, n
       do j = i+1,n
          if( A(i,j) /= A(j,i) ) &
               stop 'qreig() solves only symmetric eigenvalue problems.'
       end do
    end do

    ! Eigenvalues are computed interactively by A = Q * L * Q^T =>
    ! L * Q^T = R, L = R * Q
    L = A
    do iter = 1, 100*n*precision(A)

       call qrfac(L,Q,R)
       L = matmul(R,Q)

       do i = 1, n
          write(*,*) 'L:',L(i,:)
       end do

       ! off-diagonal elements should be negligible
       d = real(0,kind)
       do i = 1, n
          d = max(maxval(abs(L(i,1:i-1))),maxval(abs(L(i,i+1:n))),d)
       end do
       write(*,*) 'iter:',iter,d,tol
       if( d < tol ) exit
    end do
    forall( i = 1:n ) e(i) = L(i,i)

    ! Eigenvectors are computed by equation: (A - e_i*I) q_i = 0, i = 1, n.
    do i = 1, n

       L = A
       forall( j = 1:n ) L(j,j) = A(j,j) - e(i)
       call qrfac(L,Q,R)

       rnorm = norm2([(R(j,j),j=1,n-1)])
       Z(n,i) = rnorm
       do j = n-1,1,-1
          j1 = j + 1
          Z(j,i) = -dot_product(R(j,j1:n),Z(j1:n,i)) / R(j,j)
       end do
       d = norm2(Z(:,i))
       if( d > 0 ) then
          Z(:,i) = Z(:,i) / d
       end if

    end do

  end subroutine qreig


end module qrdebug
