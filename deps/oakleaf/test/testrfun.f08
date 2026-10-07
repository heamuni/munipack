

! gfortran -std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow ../src/rfun_real64.f08 ../src/rfun_hampel.f08  ../src/rfun.f08 testrfun.f08

program testrfun

  use rfun

  implicit none

  integer, parameter :: rp = selected_real_kind(15)

  integer :: i
  real(rp) :: x, dx = 1e-5_rp

  do i = -100,100
     x = i/10.0_rp
!     write(*,*) x,huber(x),dhuber(x),tukey(x),dtukey(x)
     write(*,*) x,dhuber(x),(huber(x+dx)-huber(x-dx))/(2*dx), &
          dpshuber(x),(pshuber(x+dx)-pshuber(x-dx))/(2*dx), &
          pshuber(x),(ipshuber(x+dx)-ipshuber(x-dx))/(2*dx), &
          dtukey(x),(tukey(x+dx)-tukey(x-dx))/(2*dx), &
          tukey(x),(itukey(x+dx)-itukey(x-dx))/(2*dx), &
          dhampel(x),(hampel(x+dx)-hampel(x-dx))/(2*dx), &
          hampel(x),(ihampel(x+dx)-ihampel(x-dx))/(2*dx), &
          dpshampel(x),(pshampel(x+dx)-pshampel(x-dx))/(2*dx), &
          pshampel(x),(ipshampel(x+dx)-ipshampel(x-dx))/(2*dx), &
          ihampel(x), ipshampel(x)
!          dbathtub(x),(bathtub(x+dx)-bathtub(x-dx))/(2*dx), &
!          bathtub(x),(ibathtub(x+dx)-ibathtub(x-dx))/(2*dx), &
!          drat(x),(rat(x+dx)-rat(x-dx))/(2*dx), &
!          rat(x),(irat(x+dx)-irat(x-dx))/(2*dx),&
!          irat(x), ibathtub(x)
  end do

end program testrfun
