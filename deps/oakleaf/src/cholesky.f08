
module cholesky

  use iso_fortran_env
  implicit none

contains

  pure subroutine chol_rrfac(A,R,posdef)
    !!
    !! Cholesky factorisation of a square symetric real positive definite matrix
    !!
    !! $$ A = R^T R $$
    !!
    !! where $R$ is the upper triangular matrix.
    !!
    !! posdef = .true. for success
    !!
    ! Cholesky factorisation: A = R**T * R, where R is the upper triangular matrix

    real(REAL64), dimension(:,:), intent(in) :: A
    real(REAL64), dimension(:,:), intent(out) :: R
    logical, intent(out) :: posdef

    real(REAL64), parameter :: macheps = epsilon(A)
    real(REAL64) :: s
    integer :: i,k,n,l

    n = size(A,1)
    posdef = .false.

    ! Crout implementation
    ! https://en.wikipedia.org/wiki/Cholesky_decomposition
    do i = 1, n
       k = i - 1
       s = A(i,i) - dot_product(R(1:k,i),R(1:k,i))
       if( s < macheps ) return
       R(i,i) = sqrt(s)

       l = i + 1
       R(i,l:n) = (A(i,l:n) - matmul(R(1:k,i),R(1:k,l:n))) / R(i,i)
    end do
    posdef = .true.

  end subroutine chol_rrfac

  pure subroutine chol_rrsolve(R,b,x)

    ! Cholesky solution: R**T * R * x = b

    real(REAL64), dimension(:,:), intent(in) :: R
    real(REAL64), dimension(:), intent(in) :: b
    real(REAL64), dimension(:), intent(out) :: x

    integer :: i,k,n

    n = size(R,1)
    do i = 1,n
       k = i - 1
       x(i) = (b(i) - dot_product(R(1:k,i),x(1:k))) / R(i,i)
    end do
    do i = n, 1, -1
       k = i + 1
       x(i) = (x(i) - dot_product(R(i,k:n),x(k:n))) / R(i,i)
    end do

  end subroutine chol_rrsolve

  pure subroutine chol_rrtsolve(R,b,x)

    ! Cholesky solution: R**T * x = b

    real(REAL64), dimension(:,:), intent(in) :: R
    real(REAL64), dimension(:), intent(in) :: b
    real(REAL64), dimension(:), intent(out) :: x

    integer :: j,k,n

    n = size(R,1)
    do j = 1, n
       k = j - 1
       x(j) = (b(j) - dot_product(R(1:k,j),x(1:k))) / R(j,j)
    end do

  end subroutine chol_rrtsolve

  pure subroutine chol_rsolve(R,b,x)

    ! Cholesky solution: R * x = b

    real(REAL64), dimension(:,:), intent(in) :: R
    real(REAL64), dimension(:), intent(in) :: b
    real(REAL64), dimension(:), intent(out) :: x

    integer :: j,k,n

    n = size(R,1)
    do j = n, 1, -1
       k = j + 1
       x(j) = (b(j) - dot_product(R(j,k:n),x(k:n))) / R(j,j)
    end do

  end subroutine chol_rsolve

end module cholesky
