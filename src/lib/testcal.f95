

! gfortran -Wall --check=all  -I../lib rratio.f95 testcal.f95 -L. -L../lib -L../minpack -lrstat -lsort -llmin -lminpacks -lminpack -lm


program testcal

  use robratio
  implicit none

  integer, parameter :: dbl = selected_real_kind(15)
  integer, parameter :: nmax = 15
  real(dbl), dimension(nmax) :: ph,dph,ct,dct
  real(dbl) :: r,dr,x,y,z,dx,dy
  integer :: i,k

!    write(*,*) 'sss'

!  x = inverf(-1.5_dbl)
!  x = inverf(-0.7_dbl)
!  x = inverf(0.0_dbl)
!  x = inverf(0.3_dbl)
!  x = inverf(1.5_dbl)
!stop 'sss'

  goto 33
  do i = 1,nmax
     ph(i) = pnoise(200.0_dbl)
     dph(i) = sqrt(ph(i))
     ct(i) = pnoise(100.0_dbl)
     dct(i) = sqrt(ct(i))
  end do

  call rcal(ph,dph,ct,dct,r,dr,.true.)


  do i = 1,nmax
     ph(i) = gdis(20000.0_dbl,141.42_dbl)
     dph(i) = sqrt(ph(i))
     ct(i) = gdis(10000.0_dbl,100.0_dbl)
     dct(i) = sqrt(ct(i))
  end do
  call rcal(ph,dph,ct,dct,r,dr,.true.)
33 continue

  do i = 1,nmax
     r = 1000 + 20*i
     !ph(i) = gdis(2*r,sqrt(2*r))
     ph(i) = pnoise(1*r)
     dph(i) = sqrt(ph(i))
!     ct(i) = gdis(r,sqrt(r))
     ct(i) = pnoise(r)
     dct(i) = sqrt(ct(i))
  end do
  call rcal(ph,dph,ct,dct,r,dr,.true.)

!  write(*,*) 'Non-robust:',sqrt(sum((ph - ct)**2) / nmax)

  stop

  open(1,file='/tmp/xdist')
  do k = 1,10000,10
     r = k
     do i = 1,nmax
        ph(i) = pnoise(r)
     end do
     x = sum(ph) / nmax
     dx = sqrt(sum((ph-x)**2)/nmax)
     do i = 1,nmax
        ph(i) = gdis(r,sqrt(r))
     end do
     y = sum(ph) / nmax
     dy = sqrt(sum((ph-y)**2)/nmax)
     write(1,*) k,x,dx,y,dy
  end do
  close(1)


  open(1,file='/tmp/err')
  do i = -5000,5000,10
     r = i / 1e3
     x = erf(r)
     write(1,*) r,x,ierf(x)-r,inverf(x)-r
  end do
  close(1)

  ! comparison rratios, Anscombe and magnitudes
  ! https://en.wikipedia.org/wiki/Anscombe_transform
  open(1,file='/tmp/ans')
  do k = 1,8
     r = 100 + 10**k
     do i = 1,nmax
        ph(i) = pnoise(r)
     end do
     x = sum(ph/sqrt(ph))/nmax
     y = sum(log10(ph))/nmax
     z = sum(2*sqrt(ph+3.0/8.0))/nmax
     write(*,*) r,x**2/r,10**y/r, &
          ((z/2)**2 - 1.0/8.0 + 0.306186218/z - 1.375/z**2 + 0.765465545/z**3 )/r, &
          (sum(ph)/nmax)/r
  end do
  close(1)

contains

  function pnoise(lam)

    real(dbl), intent(in) :: lam
    real(dbl) :: pnoise

    if( lam  < 500 ) then
       pnoise = pnoise_knuth(lam)
    else
       pnoise = gdis(lam,sqrt(lam))
    end if

  end function pnoise


  function pnoise_junhao(lam) result(k)

    ! Junhao, based on Knuth:
    !    http://en.wikipedia.org/wiki/Poisson_distribution

    real(dbl), intent(in) :: lam
    real(dbl), parameter :: e = exp(1.0_dbl)
    real(dbl), parameter :: lx = 500
    real(dbl), parameter :: elx = exp(lx)
    real(dbl) :: l,p,u
    integer :: k

    l = lam
    k = 0
    p = 1
    do
       k = k + 1
       call random_number(u)
       p = p*u
       if( p < e .and. l > 0 ) then
          if( l > lx ) then
             p = p * elx
             l = l - lx
          else
             p = p * exp(l)
             l = -1
          end if
       end if
       if( .not. (p > 1) ) exit
    enddo
    k = k - 1

  end function pnoise_junhao

  function pnoise_knuth(lam) result(k)

    ! Knuth's algorithm by http://en.wikipedia.org/wiki/Poisson_distribution

    real(dbl), intent(in) :: lam
    real(dbl) :: L,p,u
    integer :: k

    L = exp(-lam)
    k = 0
    p = 1
    do
       k = k + 1
       call random_number(u)
       p = p*u
       if( .not. (p > L) ) exit
    enddo
    k = k - 1

  end function pnoise_knuth

  function gdis(mean, sig)

    ! generate random data with gauss distribution

    real(dbl), intent(in) :: mean, sig
    real(dbl) :: gdis, x

    call random_number(x)
!    gdis = mean + 1.618*sqrt(2.0)*sig*ierf(2*x-1) ! Numerical Recipes (6.14.3)
    gdis = mean + sqrt(2.0)*sig*inverf(2*x-1)

  end function gdis

  function ierf(x)

    ! http://en.wikipedia.org/wiki/Error_function

    real(dbl), intent(in) :: x
    real(dbl) :: ierf, u, v, a, y, w

    a = 0.140012288686666
    u = log(1 - x**2)
    v = 2/3.14159265358979/a
    y = v + u/2

    w = sqrt(y**2 - u/a)
    ierf = sign(1.0d0,x)*sqrt(w - v)

  end function ierf


  ! by Numerical Recipes
  function inverf(p)

    real(dbl), intent(in) :: p
    real(dbl) :: inverf, pp,x,t,err
    integer :: j

    if( abs(p) < 1 ) then
       pp = 1 - abs(p)
       t = sqrt(-2.0_dbl*log(pp/2))
       x = - 0.70711*((2.30753+t*0.27061)/(1 + t*(0.99229+t*0.04481)) - t)
       do j = 1,2
          err = 1 - erf(x) - pp
          x = x + err/(1.128379167099551257 * exp(-x**2) - x*err)
       end do
       inverf = -sign(x,-p)

    else
       inverf = sign(huge(p),p)
    end if

  end function inverf


end program testcal
