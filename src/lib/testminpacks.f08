
! gfortran -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace minpacks.f95 testminpacks.f08 -lminpack -lm

program testminpacks

  use minpacks
  use iso_fortran_env

  implicit none


  real(REAL64), dimension(3,3) :: a, a1, u
  real(REAL64), dimension(3) :: b,x
  integer :: i,j

  do i = 1, 3
     do j = 1,3
        a(i, j) = real(1.0,REAL64) / real(i + 3*j - 1,REAL64)
     enddo
  enddo

  write(*,*) 'Original matrix:'
  do i = 1,3
     write(*,*) a(i,:)
  enddo

  call qrinv(a,a1)

  write(*,*)
  write(*,*) 'Inverted matrix:'
  do i = 1,3
     write(*,*) a1(i,:)
  enddo

  write(*,*)
  write(*,*) 'Unit matrix:'
  u = matmul(a,a1)
  do i = 1,3
     write(*,*) u(i,:)
  enddo

  write(*,*)
  write(*,*) 'Solution of equations:'

  b = [ 1,0,0 ]
  call qrsolve(a,b,x)
  write(*,*) x

  b = [ 0,1,0 ]
  call qrsolve(a,b,x)
  write(*,*) x

  b = [ 0,0,1 ]
  call qrsolve(a,b,x)
  write(*,*) x

end program testminpacks
