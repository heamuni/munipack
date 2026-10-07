
! gfortran -I/usr/include example_fmean.f08 -loakleaf

program example_fmean

   use oakleaf
   use iso_fortran_env

   implicit none

   real(REAL64), dimension(7) :: x = [ 1.0, 2.0, 3.0, 4.0,  5.0,  6.0,  7.0 ]
   real(REAL64), dimension(7) :: y = [ 2.1, 4.0, 5.9, 7.9, 10.2, 12.0, 13.7 ]
   real(REAL64) :: t,dt,sig
   logical :: reli

   call fmean(x,sqrt(x),y,sqrt(y),t,dt,sig,reli)
   write(*,*) 'Data: ',x
   write(*,*) 'Mean: ',t
   write(*,*) 'Std.err.: ',dt
   write(*,*) 'Std.dev: ',sig
   write(*,*) 'Reliable? ',reli

 end program example_fmean
