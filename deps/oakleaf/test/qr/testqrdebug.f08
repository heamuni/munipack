!
!  Tests of QR factorisation, the debug version
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

program testqrdebug

  use qrdebug

  use iso_fortran_env
  implicit none
  integer, parameter :: kind = REAL64

!  call trieigen
!   call eigen
!  call eigenb
  call eigen2
!  call eigen4
!  call eigenX
!  call invert
!  call solve
!  call pivo
!  call wiki
!  call simple
!  call lssol


contains

  subroutine eigen

    real(REAL64), dimension(3,3) :: A,D,B,Q,W
    real(REAL64), dimension(3) :: e
    integer :: i

    A(1,:) = real([  3,  -1,   1 ],REAL64)
    A(2,:) = real([ -1,   5,  -1 ],REAL64)
    A(3,:) = real([  1,  -1,   3 ],REAL64)

    W(1,:) = real([  1,   1,  -1 ],REAL64)
    W(2,:) = real([ -2,   1,   0 ],REAL64)
    W(3,:) = real([  1,   1,   1 ],REAL64)

    call qreig(A,e,Q)

    B = matmul(transpose(Q),matmul(A,Q))
    do i = 1, 3
       write(*,*) 'B:',B(i,:)
    end do

    do i = 1, 3
       write(*,*) i,e(i),':',Q(:,i)
    end do

    write(*,*) 'res e:',e-real([6,3,2],REAL64)

    do i = 1,3
       write(*,*) 'res Q:',Q(:,i)/norm2(Q(:,i)) - W(:,i)/norm2(W(:,i))
    end do

    do i = 1,3
       write(*,*) 'A*u-e*u', matmul(A,Q(:,i))-e(i)*Q(:,i)
    end do

  end subroutine eigen

  subroutine eigenb

    ! example by Bulrisch

    real(REAL64), dimension(5,5) :: T, Q, A0
    real(REAL64), dimension(5) :: lam
    integer :: i

    T(1,:) = real([ 12, 1, 0, 0, 0 ],REAL64)
    T(2,:) = real([  1, 9, 1, 0, 0 ],REAL64)
    T(3,:) = real([  0, 1, 6, 1, 0 ],REAL64)
    T(4,:) = real([  0, 0, 1, 3, 1 ],REAL64)
    T(5,:) = real([  0, 0, 0, 1, 0 ],REAL64)

    A0 = T
    call qreig4(T,lam,Q)

    do i = 1, 5
       write(*,*) i,real(lam(i)),':',real(Q(i,:))
    end do

    do i = 1,5
       write(*,*) i, matmul(A0,Q(:,i))-lam(i)*Q(:,i)
    end do

  end subroutine eigenb

  subroutine eigen2

    real(kind), dimension(2,2) :: A,A0,Q,W,B
    real(kind), dimension(2) :: lam
    real(kind), parameter :: two = real(2,kind)
    integer :: i

    A(1,:) = real([ 2, -1],kind)
    A(2,:) = real([-1,  2],kind)

    A(1,:) = [1.0003937921998791  ,     -6.2761672849150837E-003]
    A(2,:) = [-6.2761672849150837E-003 ,  2.2264813751428720]

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

    write(*,*) 'res lam:',lam-real([3,1],kind)

    do i = 1,2
       write(*,*) 'res Q:',Q(:,i)/norm2(Q(:,i)) - W(:,i)/norm2(W(:,i))
    end do

  end subroutine eigen2

  subroutine eigen4

    real(kind), dimension(4,4) :: A,A0,B,Q,W
    real(kind), dimension(4) :: lam, lam0
    integer :: i

    A(1,:) = real([  4,  -1,   2,  1 ],kind)
    A(2,:) = real([ -1,   2,   0, -1 ],kind)
    A(3,:) = real([  2,   0,  -1,  3 ],kind)
    A(4,:) = real([  1,  -1,   3,  0 ],kind)

    ! result of numpy.linalg.eigh(A) (LAPACK routines _syevd, _heevd)
    ! .eig() gives result only on 7 digits!
    lam0(1) = 5.886433276120998315_kind
    lam0(2) =-3.748832947424356909_kind
    lam0(3) = 1.747023266497336413_kind
    lam0(4) = 1.115376404806020183_kind

    W(1,1) =  0.779681727939313785_kind
    W(1,2) =  0.135838867510508166_kind
    W(1,3) = -0.342361742878208863_kind
    W(1,4) = -0.506391787260046744_kind
    W(2,1) = -0.299414901547905177_kind
    W(2,2) =  0.130213554798958242_kind
    W(2,3) = -0.924079049152959797_kind
    W(2,4) =  0.198678277107084189_kind
    W(3,1) =  0.393714173008689550_kind
    W(3,2) = -0.767558121913135216_kind
    W(3,3) = -0.130668625766089619_kind
    W(3,4) =  0.488640348004201808_kind
    W(4,1) =  0.383974308802957409_kind
    W(4,2) =  0.612737106518989716_kind
    W(4,3) =  0.108591243525245518_kind
    W(4,4) =  0.682147279041461330_kind

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
       write(*,*) 'res Q:',Q(:,i)/norm2(Q(:,i)) - W(:,i)/norm2(W(:,i))
    end do

    do i = 1,4
       write(*,*) i, matmul(A0,Q(:,i))-lam(i)*Q(:,i)
    end do

  end subroutine eigen4

  subroutine eigenX

    real(kind), parameter :: pi = real(3.14159265358979323844,kind)
    integer, parameter :: n = 11
    real(kind), dimension(n,n) :: A,A0,B,Q
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
       write(*,*) i, matmul(A0,Q(:,i))-lam(i)*Q(:,i)
    end do

  end subroutine eigenX

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

  subroutine invert

    real(REAL64), dimension(3,3) :: A,X,E
    integer :: i

    A(1,:) = [ 4, -1, 0 ]
    A(2,:) = [-1,  4,-1 ]
    A(3,:) = [ 0, -1, 4 ]

    call qrinv(A,X)

    do i = 1, size(A,1)
       write(*,*) 'X(i,:)',X(i,:)
    end do

    E = matmul(A,X)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

    E = matmul(X,A)
    do i = 1, size(A,1)
       write(*,*) 'E(i,:)',E(i,:)
    end do

  end subroutine invert

  subroutine solve

    real(kind), dimension(3,3) :: A
    real(kind), dimension(3) :: x,b,res

    A(1,:) = [ 4, -1, 0 ]
    A(2,:) = [-1,  4,-1 ]
    A(3,:) = [ 0, -1, 4 ]
    b = [ 7, 1, 3]

    call qrsol(A,b,x)

    write(*,*) 'x:',x

    res = b - matmul(A,x)
    write(*,*) 'res:',res
    return

!    call qrpivot(A,Q,R,pivot)
!    call qrsolpivot(Q,R,pivot,b,x)
!    write(*,*) 'x:',x
!    res = b - matmul(A,x)
!    write(*,*) 'res:',res

  end subroutine solve
!!$
!!$  subroutine pivo
!!$
!!$    ! https://www.math.tamu.edu/~fnarc/m660/qr_pivot.html
!!$    real(REAL64), dimension(4,4) :: A,Q,R,QR,E,P
!!$    integer, dimension(4) :: pivot
!!$    integer :: i
!!$
!!$    A(1,:) = [ 1, -1,  3, 1 ]
!!$    A(2,:) = [-1,  3, -5, 2 ]
!!$    A(3,:) = [ 2, -1,  5, 1 ]
!!$    A(4,:) = [ 1, -1,  3, 2 ]
!!$
!!$    do i = 1, size(A,1)
!!$       write(*,*) 'A(i,:)',A(i,:)
!!$    end do
!!$
!!$    call qrpivot(A,Q,R,pivot)
!!$
!!$    write(*,*) pivot
!!$
!!$    do i = 1, size(A,1)
!!$       write(*,*) 'Q(i,:)',Q(i,:)
!!$    end do
!!$
!!$    do i = 1, size(A,1)
!!$       write(*,*) 'R(i,:)',R(i,:)
!!$    end do
!!$
!!$    P = 0
!!$    do i = 1, size(A,1)
!!$       P(pivot(i),i) = 1
!!$       write(*,*) 'P(i,:)',P(i,:)
!!$    end do
!!$
!!$    QR = matmul(Q,matmul(R,transpose(P)))
!!$    do i = 1, size(A,1)
!!$       write(*,*) 'res(i,:)',QR(i,:)-A(i,:)
!!$    end do
!!$
!!$    E = matmul(transpose(Q),Q)
!!$    do i = 1, size(A,1)
!!$       write(*,*) 'E(i,:)',E(i,:)
!!$    end do
!!$
!!$    call qrfac(A,Q,R)
!!$
!!$    do i = 1, size(A,1)
!!$       write(*,*) 'Q(i,:)',Q(i,:)
!!$    end do
!!$
!!$    do i = 1, size(A,1)
!!$       write(*,*) 'R(i,:)',R(i,:)
!!$    end do
!!$
!!$
!!$  end subroutine pivo

  subroutine wiki

    ! https://en.wikipedia.org/wiki/QR_decomposition

    real(kind), dimension(3,3) :: A,Q,R,QR,E,Qi,Ri,P,AP
    real(kind), dimension(3) :: beta
    integer, dimension(3) :: pivot
    integer :: i

    A(1,:) = [12,-51, 4 ]
    A(2,:) = [ 6,167,-68]
    A(3,:) = [-4, 24,-41]

    Qi(1,:) = [ 6._kind/7._kind,  -69._kind/175._kind, -58._kind/175._kind ]
    Qi(2,:) = [ 3._kind/7._kind,  158._kind/175._kind,   6._kind/175._kind ]
    Qi(3,:) = [-2._kind/7._kind,    6._kind/35._kind,  -33._kind/35._kind ]

    Ri(1,:) = [ 14,  21,-14 ]
    Ri(2,:) = [  0, 175,-70 ]
    Ri(3,:) = [  0,   0,-35 ]

    call qrfac(A,Q,R)

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

    call qrpivot(A,Q,R,pivot)

    write(*,*) pivot

    do i = 1, size(A,1)
       write(*,*) 'Q(i,:)',Q(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'R(i,:)',R(i,:)
    end do

    P = 0
    do i = 1, size(A,1)
       P(pivot(i),i) = 1
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

  end subroutine wiki

  subroutine simple

    real(kind), dimension(3,3) :: A,Q,R
    integer :: i

    A(1,:) = [ 5, 0, 0]
    A(2,:) = [ 0, 3*0, 0]
    A(3,:) = [ 0, 0,-1]

    call qrfac(A,Q,R)

    do i = 1, size(A,1)
       write(*,*) 'Q(i,:)',Q(i,:)
    end do

    do i = 1, size(A,1)
       write(*,*) 'R(i,:)',R(i,:)
    end do

  end subroutine simple

  subroutine lssol

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
    b = y

    call qrlssol(A,b,xls)

    write(*,*) 'xls:',xls
    write(*,*) 'res:',xls - [29.223, -1.913]

  end subroutine lssol


end program testqrdebug
