!
!  Tests of QR related algorithms
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

! gfortran  -std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow ../src/qrhouse.f08 testqrhouse.f08 && ./a.out

program testqrhouse

  use QRHouse

  use iso_fortran_env
  implicit none
  integer, parameter :: kind = REAL64
  logical, dimension(13) :: test
  integer :: msuite, ntest
  real(kind) :: tol = 1e2*epsilon(real(1,kind))

  call get_environment_variable("AM_TEST_SUITE",status=msuite)

  if( msuite == 0 ) then
     ntest = 4

     ! selected tests for automake
     test(1) = wiki()
     test(2) = solve()
     test(3) = lssol()
     test(4) = ranksol()

  else
     ntest = size(test)

     ! all tests
     test(1) = wiki()
     test(2) = solve()
     test(3) = lssol()
     test(4) = ranksol()
     test(5) = pivo()
     test(6) = rankdefic()
     test(7) = invert()
     test(8) = invertrank()
     test(9) = eigen()
     test(10) = eigenb()
     test(11) = eigen2()
     test(12) = eigen4()
     test(13) = eigenX()

  end if
  write(*,*) 'All tests done: ',ntest,'/',ntest

  if( all(test(1:ntest)) ) then
     write(*,*) 'testqrhouse = ','PASS'
  else
     write(*,*) 'testqrhouse = ','FAIL:',test(1:ntest)
  end if


  ! additional tests
  ! call trieigen
  ! call zpivot()

contains

  logical function eigen()

    real(REAL64), dimension(3,3) :: A,D,B,Q,W, qres
    real(REAL64), dimension(3) :: lam
    integer :: i

    A(1,:) = [  4,  12, -16 ]
    A(2,:) = [ 12,  37, -43 ]
    A(3,:) = [-16, -43,  98 ]

    call qreig(A,lam,Q)

    do i = 1, 3
       write(*,*) i,real(lam(i)),':',real(Q(i,:))
    end do

    A(1,:) = [  3,  -1,   1 ]
    A(2,:) = [ -1,   5,  -1 ]
    A(3,:) = [  1,  -1,   3 ]

    W(1,:) = [ -1,   1,   1 ]
    W(2,:) = [  2,   0,   1 ]
    W(3,:) = [ -1,  -1,   1 ]

    D = A
    call qreig(D,lam,Q)

    B = matmul(transpose(Q),matmul(D,Q))
    do i = 1, 3
       write(*,*) 'B:',B(i,:)
    end do

    do i = 1, 3
       write(*,*) i,lam(i),':',Q(:,i)
    end do

    write(*,*) 'res lam:',lam-[6,2,3]

    do i = 1,3
       write(*,*) 'res e:',Q(:,i)/norm2(Q(:,i)) - W(:,i)/norm2(W(:,i))
    end do

    do i = 1,3
       qres(i,:) = matmul(A,Q(:,i))-lam(i)*Q(:,i)
       write(*,*) 'A*u-lam*u', qres(i,:)
    end do

    eigen = norm2(lam-[6,2,3]) < tol .and. norm2(qres) < tol

  end function eigen

  logical function eigenb()

    ! example by Bulrisch

    real(REAL64), dimension(5,5) :: T, Q, A0, qres
    real(REAL64), dimension(5) :: lam
    integer :: i

    T(1,:) = [ 12, 1, 0, 0, 0 ]
    T(2,:) = [  1, 9, 1, 0, 0 ]
    T(3,:) = [  0, 1, 6, 1, 0 ]
    T(4,:) = [  0, 0, 1, 3, 1 ]
    T(5,:) = [  0, 0, 0, 1, 0 ]

    A0 = T
    call qreig(T,lam,Q)

    do i = 1, 5
       write(*,*) i,real(lam(i)),':',real(Q(i,:))
    end do

    do i = 1,5
       qres(i,:) = matmul(A0,Q(:,i))-lam(i)*Q(:,i)
       write(*,*) i, qres(i,:)
    end do

    write(*,*) 'eigenb():',norm2(qres),tol
    eigenb = norm2(qres) < 1e2*tol

  end function eigenb

  logical function eigen2()

    real(kind), dimension(2,2) :: A,A0,Q,W,B,qres
    real(kind), dimension(2) :: lam
    real(kind), parameter :: two = real(2,kind)
    integer :: i

    write(*,*) 'eigen2():'

    A(1,:) = [ 2, -1]
    A(2,:) = [-1,  2]

    W(1,:) = [ sqrt(two)/two, sqrt(two)/two ]
    W(2,:) = [-sqrt(two)/two, sqrt(two)/two ]

    A0 = A
    call qreig(A,lam,Q)

    B = matmul(Q,matmul(A0,transpose(Q)))
    do i = 1, 2
       write(*,*) 'B:',B(i,:)
    end do
    do i = 1, 2
       write(*,*) 'res B:',B(i,:)-A0(i,:)
    end do

    do i = 1, 2
       write(*,*) i,lam(i),':',Q(:,i)
    end do

    write(*,*) 'res lam:',lam-[3,1]

    do i = 1,2
       qres(i,:) = Q(:,i)/norm2(Q(:,i)) - W(:,i)/norm2(W(:,i))
       write(*,*) 'res Q:',qres(i,:)
    end do

    write(*,*) norm2(lam-[3,1]),norm2(qres),tol
    eigen2 = norm2(lam-[3,1]) < tol .and. norm2(qres) < tol

  end function eigen2

  logical function eigen4()

    real(kind), dimension(4,4) :: A,A0,B,Q,W, qres, wres
    real(kind), dimension(4) :: lam, lam0
    integer :: i

    A(1,:) = [  4,  -1,   2,  1 ]
    A(2,:) = [ -1,   2,   0, -1 ]
    A(3,:) = [  2,   0,  -1,  3 ]
    A(4,:) = [  1,  -1,   3,  0 ]

    ! result of numpy.linalg.eigh(A) (LAPACK routines _syevd, _heevd)
    lam0(1) = 5.886433276120998315_kind
    lam0(2) =-3.748832947424356909_kind
    lam0(3) = 1.115376404806020183_kind
    lam0(4) = 1.747023266497336413_kind

    W(1,1) = -0.779681727939313785_kind
    W(1,2) = -0.135838867510508166_kind
    W(1,3) =  0.506391787260046744_kind
    W(1,4) =  0.342361742878208863_kind
    W(2,1) =  0.299414901547905177_kind
    W(2,2) = -0.130213554798958242_kind
    W(2,3) = -0.198678277107084189_kind
    W(2,4) =  0.924079049152959797_kind
    W(3,1) = -0.393714173008689550_kind
    W(3,2) =  0.767558121913135216_kind
    W(3,3) = -0.488640348004201808_kind
    W(3,4) =  0.130668625766089619_kind
    W(4,1) = -0.383974308802957409_kind
    W(4,2) = -0.612737106518989716_kind
    W(4,3) = -0.682147279041461330_kind
    W(4,4) = -0.108591243525245518_kind

    A0 = A
    call qreig(A,lam,Q)

    B = matmul(transpose(Q),matmul(A0,Q))
    do i = 1, 4
       write(*,*) 'B:',B(i,:)
    end do

    do i = 1, 4
       write(*,*) i,lam(i),':',Q(:,i)
    end do

    write(*,*) 'res lam:',lam-lam0

    do i = 1,4
       wres(i,:) = Q(:,i)/norm2(Q(:,i)) - W(:,i)/norm2(W(:,i))
       write(*,*) 'res Q:',wres(i,:)
    end do

    do i = 1,4
       qres(i,:) = matmul(A0,Q(:,i))-lam(i)*Q(:,i)
       write(*,*) i, qres(i,I)
    end do

    write(*,*) 'eigen4():',norm2(lam-lam0),norm2(qres),norm2(wres),tol
    eigen4 = norm2(lam-lam0) < tol .and. norm2(qres) < tol .and. norm2(wres) < tol

  end function eigen4

  logical function eigenX()

    real(kind), parameter :: pi = real(3.14159265358979323844,kind)
    integer, parameter :: n = 11
    real(kind), dimension(n,n) :: A,A0,B,Q, qres
    real(kind), dimension(n) :: lam
    integer :: i, j
    real(kind) :: x

    do i = 1, n
       do j = i, n
          do
             call random_number(x)
             if( x > 0.01 ) then
                A(i,j) = x
                A(j,i) = x
                exit
             end if
          end do
       end do
    end do

    do i = 1, n
       do j = 1, n
          A(i,j) = 2*min(i,j) - 1
       end do
    end do

    do i = 1, n
       write(*,*) A(i,:)
    end do

    A0 = A
    call qreig(A,lam,Q)

    B = matmul(transpose(Q),matmul(A0,Q))
    do i = 1, n
       write(*,*) 'B:',B(i,:)
    end do

    do i = 1, n
       write(*,*) i,lam(i),':',Q(:,i)
    end do

    do i = 1, n
       write(*,*) 'res lam:',lam(i) - 0.5 / cos ( (2*i-1)*pi / (4*n) )**2
    end do

    do i = 1,n
       qres(i,:) = matmul(A0,Q(:,i))-lam(i)*Q(:,i)
       write(*,*) i, qres(i,:)
    end do

    eigenX = norm2(qres) < 1e3*tol

  end function eigenX

!!$  subroutine trieigen
!!$
!!$    real(REAL64), dimension(4,4) :: A,T,W
!!$!    real(REAL64), dimension(4) :: lam
!!$    integer :: i
!!$
!!$    ! https://en.wikipedia.org/wiki/Householder_transformation
!!$    A(1,:) = real([  4,   2,  -2, 1 ],REAL64)
!!$    A(2,:) = real([  2,   3,   2, 1 ],REAL64)
!!$    A(3,:) = real([ -2,   2,   1, 0 ],REAL64)
!!$    A(4,:) = real([  1,   1,   0, 2  ],REAL64)
!!$
!!$    W(1,:) = real([  4.0,     3.0,     0.0,     0.0 ],REAL64)
!!$    W(2,:) = real([  3.0, 2.0/3.0, 5.0/3.0,     0.0 ],REAL64)
!!$    W(3,:) = real([  0.0, 5.0/3.0,     3.0,-4.0/3.0 ],REAL64)
!!$    W(4,:) = real([  0.0,     0.0,-4.0/3.0, 7.0/3.0 ],REAL64)
!!$    W(2,2:3) = [2,5] / real(3,REAL64)
!!$    W(3,2) = 5 / real(3,REAL64)
!!$    W(3,4) = -4 / real(3,REAL64)
!!$    W(4,3:4) = [-4, 7] / real(3,REAL64)
!!$
!!$!    call qreig(A,lam,Q)
!!$    call hausetridig(A,T)
!!$
!!$    do i = 1, 4
!!$       write(*,*) i,':',T(:,i)
!!$    end do
!!$
!!$    do i = 1, 4
!!$       write(*,*) i,':',T(:,i)-W(:,i)
!!$    end do
!!$
!!$
!!$  end subroutine trieigen
!!$


!!$  subroutine zpivot()
!!$
!!$    real(REAL64), dimension(3,3) :: A
!!$    real(REAL64), dimension(3) :: b,q,res
!!$
!!$    A(1,:) = [1.25, 0.0, 2.5]
!!$    A(2,:) = [0.0, 1.25, 2.5]
!!$    A(3,:) = [-2.5, -2.5, 0.0]
!!$    b = [-1.0, -1.0, -1.125]
!!$
!!$    call qrsol(A,b,q)
!!$
!!$    res = matmul(A,q) - b
!!$    write(*,*) 'zpivot(): norm(res):',norm2(res)
!!$
!!$  end subroutine zpivot


  logical function invert()

    real(REAL64), dimension(3,3) :: A,X,E,B
    integer :: i

    write(*,*) 'invert():'

    A(1,:) = [ 4, -1, 0 ]
    A(2,:) = [-1,  4,-1 ]
    A(3,:) = [ 0, -1, 4 ]

    A(1,:) = [ 1, 3, 1 ]
    A(2,:) = [ 0, 3, 1 ]
    A(3,:) = [ 4, 2, 0 ]

    B(1,:) = [ 2, -2, 0 ]
    B(2,:) = [-4,  4, 1 ]
    B(3,:) = [12,-10,-3 ]
    B = B / 2

    call qrinv(A,X)

    do i = 1, size(A,1)
       write(*,*) 'B(i,:)',B(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'X(i,:)',X(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'res X(i,:)',X(i,:) - B(i,:)
    end do
    write(*,*) 'norm:',norm2(X - B)

    E = matmul(A,X)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

    E = matmul(X,A)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

    invert = norm2(X - B) < 10*tol

  end function invert

  logical function invertrank()

    real(REAL64), dimension(3,3) :: A,X,B,E,AX,XA,T
    integer :: i

    A(1,:) = [  5,  7,   1 ]
    A(2,:) = [  5,  7,   1 ]
    A(3,:) = [ -1,  2,   1 ]
    A(2,:) = A(1,:) + A(3,:)

    call qrinv(A,X)

    do i = 1, size(A,1)
       write(*,*) 'X(i,:)',X(i,:)
    end do

    B(1,:) = [ 0.1952381 , -0.08095238, -0.27619048]
    B(2,:) = [-0.03428571,  0.09714286,  0.13142857]
    B(3,:) = [-0.06952381,  0.05809524,  0.12761905]

    write(*,*) 'Resudials against scipy.linalg.pinv(A):'
    do i = 1, size(A,1)
       write(*,*) 'res X(i,:)',X(i,:) - B(i,:)
    end do
    write(*,*) 'norm2(X - B):',norm2(X - B)

    write(*,*) "Unit matrixes:"

    E = matmul(A,X)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

    E = matmul(X,A)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

    write(*,*) "The conditions AXA=A, XAX=X are satisfied:"

    E = matmul(X,matmul(A,X))
    do i = 1, size(A,1)
       write(*,*) 'res XAX(i,:)',E(i,:) - X(i,:)
    end do
    write(*,*) 'norm2(XAX - X):',norm2(E - X)

    E = matmul(A,matmul(X,A))
    do i = 1, size(A,1)
       write(*,*) 'res AXA(i,:)',E(i,:) - A(i,:)
    end do
    write(*,*) 'norm2(AXA - A):',norm2(E - A)


    write(*,*) "The conditions (AX)**T=AX, (XA)**T=XA are also satisfied:"

    AX = matmul(A,X)
    T = transpose(AX)
    do i = 1, size(A,1)
       write(*,*) 'AX(i,:)**T - AX(i,:)',T(i,:) - AX(i,:)
    end do
    write(*,*) 'norm2(AX**T - AX):',norm2(T - AX)

    XA = matmul(X,A)
    T = transpose(XA)
    do i = 1, size(A,1)
       write(*,*) 'XA(i,:)**T - XA(i,:)',T(i,:) - XA(i,:)
    end do
    write(*,*) 'norm2(XA**T - XA):',norm2(T - XA)

    invertrank = norm2(X - B) < 1e-7

  end function invertrank


  logical function solve()

    real(kind), dimension(3,3) :: A
    real(kind), dimension(3) :: x,b,res

    A(1,:) = [ 4, -1, 0 ]
    A(2,:) = [-1,  4,-1 ]
    A(3,:) = [ 0, -1, 4 ]
    b = [ 7, 1, 3]

    call qrsol(A,b,x)

    write(*,*) 'x:',x

    res = x - [ 2,1,1]
    write(*,*) 'res, norm2:',res,norm2(res)

    solve = norm2(res) < tol

  end function solve

  logical function pivo()

    ! https://www.math.tamu.edu/~fnarc/m660/qr_pivot.html
    real(REAL64), dimension(4,4) :: A,Q,R,QR,E,P
    real(REAL64), dimension(4) :: beta
    integer, dimension(4) :: pivot
    integer :: i

    A(1,:) = [ 1, -1,  3, 1 ]
    A(2,:) = [-1,  3, -5, 2 ]
    A(3,:) = [ 2, -1,  5, 1 ]
    A(4,:) = [ 1, -1,  3, 2 ]

    do i = 1, size(A,1)
       write(*,*) 'A(i,:)',A(i,:)
    end do

    QR = A
    call qrpivot(QR,beta,pivot)
    call qform(QR,beta,Q)
    call rform(QR,R)

    write(*,*) pivot

    do i = 1, size(A,1)
       write(*,*) 'Q(i,:)',Q(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'R(i,:)',R(i,:)
    end do

    P = 0
    forall( i = 1:size(A,1) ) P(pivot(i),i) = 1
    do i = 1, size(A,1)
       write(*,*) 'P(i,:)',P(i,:)
    end do

    QR = matmul(Q,matmul(R,transpose(P)))
    do i = 1, size(A,1)
       write(*,*) 'res(i,:)',QR(i,:)-A(i,:)
    end do

    E = matmul(transpose(Q),Q)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

    pivo = norm2(QR - A) < 1e2*tol

  end function pivo

  logical function wiki()

    ! https://en.wikipedia.org/wiki/QR_decomposition

    real(kind), dimension(3,3) :: A,Q,R,QR,E,Qi,Ri,P, AP
    real(kind), dimension(3) :: beta
    integer, dimension(3) :: pivot
    integer :: i
    logical :: res, res1

    A(1,:) = [12,-51, 4 ]
    A(2,:) = [ 6,167,-68]
    A(3,:) = [-4, 24,-41]

    Qi(1,:) = [ 6._kind/7._kind,  -69._kind/175._kind, -58._kind/175._kind ]
    Qi(2,:) = [ 3._kind/7._kind,  158._kind/175._kind,   6._kind/175._kind ]
    Qi(3,:) = [-2._kind/7._kind,    6._kind/35._kind,  -33._kind/35._kind ]

    Ri(1,:) = [ 14,  21,-14 ]
    Ri(2,:) = [  0, 175,-70 ]
    Ri(3,:) = [  0,   0, 35 ]

    ! NO pivoting
    QR = A
    call qrfac(QR,beta)
    call rform(QR,R)
    call qform(QR,beta,Q)

    do i = 1, size(A,1)
       write(*,*) 'Q(i,:)',Q(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'R(i,:)',R(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'res Q(i,:)',Q(i,:)-Qi(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'res R(i,:)',R(i,:)-Ri(i,:)
    end do

    QR = matmul(Q,R)
    do i = 1, size(A,1)
       write(*,*) 'res(i,:)',QR(i,:)-A(i,:)
    end do

    E = matmul(transpose(Q),Q)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

    res = norm2(R-Ri) < tol .and. norm2(Q-Qi) < tol

    QR = A
    call qrpivot(QR,beta,pivot)
    call rform(QR,R)
    call qform(QR,beta,Q)

    write(*,*) pivot

    do i = 1, size(A,1)
       write(*,*) 'Q(i,:)',Q(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'R(i,:)',R(i,:)
    end do

    P = 0
    forall( i = 1:size(A,1)) P(pivot(i),i) = 1
    do i = 1, size(A,1)
       write(*,*) 'P(i,:)',P(i,:)
    end do
    AP = matmul(A,P)
    QR = matmul(Q,R)
    do i = 1, size(A,1)
       write(*,*) 'res(i,:)',QR(i,:)-AP(i,:)
    end do
    write(*,*) 'norm(res):',norm2(QR - AP),tol

    E = matmul(transpose(Q),Q)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

    write(*,*) 'wiki() -- qrpivot:',norm2(QR - AP)
    res1 = norm2(QR - AP) < tol*abs(R(1,1))

    wiki = res .and. res1

  end function wiki


  logical function rankdefic()

    ! https://www.math.usm.edu/lambers/mat610/sum10/lecture11.pdf

    real(kind), dimension(4,4) :: A,QR,Q,R,P,E,AP
    real(kind), dimension(3,3) :: R11
    integer, dimension(4) :: pivot
    real(kind), dimension(4) :: beta
    real(kind), dimension(3) :: rvec
    real(kind) :: rnorm
    integer :: i

    A(1,:) = [ 1, 3, 5, 1 ]
    A(2,:) = [ 2,-1, 2, 1 ]
    A(3,:) = [ 1, 4, 6, 1 ]
    A(4,:) = [ 4, 5,10, 1 ]

    do i = 1, size(A,1)
       write(*,*) 'A(i,:)',A(i,:)
    end do

    QR = A
    call qrpivot(QR,beta,pivot)
    call rform(QR,R)
    call qform(QR,beta,Q)

    write(*,*) pivot

    do i = 1, size(A,1)
       write(*,*) 'Q(i,:)',Q(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'R(i,:)',R(i,:)
    end do

    P = 0
    forall( i = 1:size(A,1)) P(pivot(i),i) = 1
    do i = 1, size(A,1)
       write(*,*) 'P(i,:)',P(i,:)
    end do
    AP = matmul(A,P)
    QR = matmul(Q,R)
    do i = 1, size(A,1)
       write(*,*) 'res(i,:)',QR(i,:)-AP(i,:)
    end do

    E = matmul(transpose(Q),Q)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

    call qrinv(R(1:3,1:3),R11)
    rvec = matmul(R11,R(1:3,4))
    rnorm = norm2(rvec)
    write(*,*) rnorm, sqrt(1+rnorm**2),rvec

    rankdefic = norm2(QR-AP) < 1e2*tol

  end function rankdefic

  logical function lssol()

    integer, parameter :: m = 14
    integer, parameter :: n = 2
    real(kind), dimension(m,n) :: A
    real(kind), dimension(m) :: x,y,b
    real(kind), dimension(n) :: xls

    ! Rektoris
    x = [ 3,4,6,6,6,7,8,9,11,11,12,12,14,16 ]
    y = [ 24.82, 23.26, 14.77, 19.06, 14.79, 17.66, 11.83, &
          14.82, 5.16, 11.12, 8.04, 2.72, 0.74, 1.21 ]

    A(:,1) = 1
    A(:,2) = x
!    A(:,3) = x**2
    b = y

    call qrsol(A,b,xls)
    write(*,*) 'xls:',xls
    write(*,*) 'res:',xls - [29.223, -1.913]

    lssol = norm2(xls - [29.223, -1.913]) < 1e-3

  end function lssol

  logical function ranksol()

    real(REAL64), dimension(3,3) :: A
    real(REAL64), dimension(3) :: b,x,res,d

    A(1,:) = [  5,  7,   1 ]
    A(3,:) = [ -1,  2,   1 ]
    A(2,:) = A(1,:) + A(3,:)
    b = [12, 12, 1]
    b(2) = b(1) + b(3)

    write(*,*) 'b:',b

    call qrsol(A,b,x)
    write(*,*) 'sol:',x

    res = x - [1,1,0]
    write(*,*) 'res, norm2:',res,norm2(res)

    res = matmul(A,x) - b
    write(*,*) 'norm(res):',norm2(res)

    d = [1.0142857142857147_kind , 0.98285714285714287_kind,  &
         4.8571428571428488e-002_kind]
    write(*,*) 'res norm (lapack):',norm2(d-x)

    ranksol = norm2(x - d) < max(10*epsilon(x),real(1e-15,kind))

  end function ranksol

end program testqrhouse
