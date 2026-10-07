
! gfortran -I/usr/include example_rmean.f08 -loakleaf

program example_rmean

   use oakleaf
   use iso_fortran_env

   implicit none

   real(REAL64), dimension(7) :: x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 ]
   real(REAL64) :: t,dt,sig
   logical :: reli

   call rmean(x,t,dt,sig,reli,verbose=.false.)
   write(*,*) 'Data: ',x
   write(*,*) 'Mean: ',t
   write(*,*) 'Std.err.: ',dt
   write(*,*) 'Std.dev: ',sig
   write(*,*) 'Reliable? ',reli

end program example_rmean
