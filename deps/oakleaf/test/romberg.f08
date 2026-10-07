module iromberg

  implicit none
  integer, parameter, private :: kind = selected_real_kind(15)

contains

  function romberg(f,a,b,tol,count,iflag)

    ! Rombergova metoda numerickeho vypoctu intergalu.
    !
    ! Znaceni:
    ! 2**MMAX = maximalni pocet vycisleni integrovane funkce F(X)
    ! Interval: <AINT,BINT>
    ! Absolutni presnost: EPS
    !
    ! iflag = 0 .. dosazena pozadovana presnost
    ! iflag = 1 .. nedosazena pozadovana presnost kvuli zaokrouhlovani
    ! iflag = 2 .. prekrocen maximalni pocet vycisleni funkce
    !

    interface
       function f(x)
         integer, parameter :: kind = selected_real_kind(15)
         real(kind) :: f
         real(kind), intent(in) :: x
       end function f
    end interface

    real(kind) :: romberg
    real(kind), intent(in) :: a,b
    real(kind), intent(in), optional :: tol
    integer, intent(in), optional :: count
    integer, intent(out), optional :: iflag

    real(kind),parameter :: macheps = epsilon(romberg)
    integer :: mmax,i,j,k
    real(kind) :: x,s,g,h,w,q,c,err,d,dmin,y,epsmin
    real(kind), dimension(:), allocatable :: t
!    double precision :: f

    if( present(tol) ) then
       err = tol
    else
       err = 42*macheps
    end if
    if( present(count) ) then
       mmax = count
    else
       mmax = 2 * nint(log((b - a) / err) / 0.6931) + 1
    end if
!    write(*,*) err,mmax

    allocate(t(0:mmax))

!    external f

! Integracni meze
!      AINT = 0.0
!      BINT = -1.0
! Zadana chyba
!      EPS = 1E-13
    !
    epsmin = 10*mmax*epsilon(1.0_kind)
    g = (b - a)/2
    h = g
    t(0) = g*(f(a) + f(b))
!    write(*,*) '0', t(0)
    k = 1
    dmin = 0
!    hmin = 0.0_kind
    do j = 1, mmax
       y = t(0)
       x = a + g
       s = 0
       c = 0
       do i = 0, k - 1
          w = f(x + h*i) - c
          q = s + w
          c = (q - s) - w
          s = q
       enddo
       t(j) = t(j - 1)/2 + s * g
       s = 1
       do i = j, 1, -1
          s = 4*s
          t(i - 1) = (s*t(i) - t(i - 1))/(s - 1)
       enddo
       !       write(*,*) k,real(t(0:k))
       d = abs(t(0) - y)
!       write(*,*) j,t(0),d,err
       if( j > 1 ) then
          if( d < j*err )then
             ! Pozadovana presnost dosazena
             if( present(iflag) ) iflag = 0
             romberg = t(0)
             return
          endif
!          if( d > epsmin )then
          if( d < macheps )then
             ! prenosti nedosahnuto kvuli zaokrouhlovani
             if( present(iflag) ) iflag = 1
             romberg = t(0)
             return
          endif
       endif
       k = 2*k
       h = g
       g = g / 2
       epsmin = d
    enddo
    if( present(iflag) ) iflag = 2
    romberg = t(0)

  end function romberg

end module iromberg
