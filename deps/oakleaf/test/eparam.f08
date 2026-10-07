
! gfortran -std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow -I$HOME/scratch/opt/include noise.f08 romberg.f08  ../src/fmm/fmm.f08  eparam.f08 -o eparam  -L$HOME/scratch/opt/lib  -loakleaf

program eparam

  use fmm
  use iromberg
  use oakleaf
  use iso_fortran_env

  implicit none
!  integer, parameter :: kind = selected_real_kind(15)
  integer, parameter :: kind = REAL64
  real(kind), parameter :: pi = 4.0_kind*atan(1.0_kind)
!  real(kind), parameter :: gamma_pshampel = 2.57602997677_REAL64
!  real(kind), parameter :: gamma_pshampel = 2.5729879_REAL64
  real(kind) :: gamma_pshampel = 2.57298793 !  2.57620926_REAL64
  real(kind), parameter :: beta_pshampel = 0.9484**2  !0.8712975324409_REAL64
  real(kind), parameter :: phr = 0.25_REAL64

  real(kind), dimension(:), allocatable :: likes
  real(kind) :: a,b,tol,s,f,psi2,dpsi, scale
  integer :: flag,n,m


  block
    integer :: i
    real(REAL64) :: a,b,tol

    a = 0.9
    b = 1.0
    tol = 1e-10
    scale = zeroin(a,b,fun_gt,tol)
    write(*,*) 'sqrt(nu) by fmin:',scale
  end block


!  a = 0
!  b = 2*atan(1.0_kind) ! pi/2
!  b = 1
!  tol = 1e-15
!  f = romberg(rtest,a,b,tol,iflag=flag)
!  write(*,*) "Romberg's integration = ",f,' res=',f-erf(b),' flag= ',flag
!  stop 0

!  write(*,*) 'ipshampel in infinity:',ipshampel(10.0_REAL64)

  a = 0.5
  b = 2
  tol = 1e-7

!!$  s = fmin(a,b,negentropy_gauss,tol)
!!$  write(*,'(a,2f7.4)') "Maximum of entropy for Normal distribution:",s,s**2
!!$
!!$  s = fmin(a,b,negentropy_huber,tol)
!!$  write(*,'(a,2f7.4)') "Maximum of entropy for Huber's minimax:",s,s**2
!!$
!!$  s = fmin(a,b,negentropy_tukey,tol)
!!$  write(*,'(a,2f7.4)') "Maximum of entropy for Tukey's biweight:",s,s**2

  a = 0
  b = 6
  tol = 1e-13
!  a = romberg(fun_tukey,a,b,tol,iflag=flag)
!  stop
  write(*,*) "Normalisation for Tukey:", norm_tukey(), &
       2*romberg(fun_tukey,a,b,tol,iflag=flag)
  write(*,*) "Normalisation for Huber:", norm_huber(), &
       2*romberg(fun_huber,a,b,tol,iflag=flag)
  write(*,*) "Normalisation for psHuber:", norm_pshuber(), &
       2*romberg(fun_pshuber,a,b,tol,iflag=flag)
  b = 8.5
  b = 8.5_REAL64 + phr
  gamma_pshampel = norm_pshampel()
  gamma_pshampel = 2*romberg(fun_pshampel,a,b,tol,iflag=flag)
  write(*,*) "Normalisation for psHampel:", gamma_pshampel, &
       2*romberg(fun_pshampel,a,b,tol,iflag=flag)
  write(*,*) "Normalisation for Hampel:", norm_hampel(), &
       2*romberg(fun_hampel,a,b,tol,iflag=flag)

  a = 0
  b = 10
  b = 8.5 + phr
  write(*,*) "Effective beta for Huber:",beta_huber(), &
       2*romberg(sig_huber,a,b,tol,iflag=flag)
  write(*,*) "Effective beta for pseudo Huber:",real(beta_pshuber()), &
       real(2*romberg(sig_pshuber,a,b,tol,iflag=flag))
  write(*,*) "Effective beta for Tukey:",beta_tukey(), &
       2*romberg(sig_tukey,a,b,tol,iflag=flag)
  b = 10
  b = 8.5 + phr
  tol = 1e-10
  write(*,*) "Effective beta for Hampel:",beta_hampel(), &
       2*romberg(sig_hampel,a,b,tol,iflag=flag)
  write(*,*) "Effective beta for psHampel:",fun_beta_pshampel(), &
       2*romberg(sig_pshampel,a,b,tol,iflag=flag)
  write(*,*) "Effective nu for Hampel:",nu_hampel(), &
       2*romberg(sig_nu_hampel,a,b,tol,iflag=flag)
  write(*,*) "Effective nu for psHampel:",fun_nu_pshampel(), &
       2*romberg(sig_nu_pshampel,a,b,tol,iflag=flag)

  a = 0
  b = 6
  tol = 1e-10
  write(*,*) "Effective beta for dHuber:",beta_dhuber(), &
       2*romberg(sig_dhuber,a,b,tol,iflag=flag)
  write(*,*) "Effective beta for pseudo dHuber:",real(beta_dpshuber()), &
       real(2*romberg(sig_dpshuber,a,b,tol,iflag=flag))
  write(*,*) "Effective beta for dTukey:",beta_dtukey(), &
       2*romberg(sig_dtukey,a,b,tol,iflag=flag)

  write(*,*) "Effective ϖ for psHampel:",fun_gamma_pshampel()
!  write(*,*) "Effective gamma1 for psHampel:",gamma1_pshampel()
!  write(*,*) "Effective nu for psHampel:",nu_pshampel()
  write(*,*) "Effective off for psHampel:",off_pshampel()

  b = 8.5+phr
  tol = 1e-10
  psi2 = romberg(fun_psi2,-b,b,tol,iflag=flag)
  dpsi = romberg(fun_dpsi,-b,b,tol,iflag=flag)
  write(*,*) "int(psi**2) vs. int(psi'):",psi2,dpsi,psi2-dpsi
  write(*,*) "int(psi**2)/int(psi') = 1:",psi2/dpsi,psi2/dpsi-1
  write(*,*) "int(psi') = sqrt(beta):",dpsi,dpsi-beta_pshampel

  f = romberg(fun_s,-b,b,tol,iflag=flag)
  write(*,*) "int(psi*x):",f
  write(*,*) "scale=",sqrt(psi2),sqrt(dpsi)

  write(*,*) "Test of the effective beta for classic, wins.,Huber,psHuber,Tukey:"
!  write(*,*)beta_huber_norm(1000000)
  do m = 3,3!8
     n = 10**m
     write(*,'(i10,5f10.5)') n, beta_sqr_norm(n), beta_wins_norm(n), &
          beta_huber_norm(n),beta_tukey_norm(n),beta_pshuber_norm(n)
  end do

!!$  write(*,'(a,2f7.4)') "Maximum of information for N(0,1):"
!!$  block
!!$    use noise_REAL64
!!$
!!$    a = 0.5
!!$    b = 1.5
!!$    tol = 1e-6
!!$
!!$    do m = 3,6
!!$       n = 10**m
!!$       allocate(likes(n))
!!$       call gnoise_fill(0.0_kind,1.0_kind,likes)
!!$       sh = fmin(a,b,likes_huber,tol)
!!$       st = fmin(a,b,likes_tukey,tol)
!!$       sp = fmin(a,b,likes_pshuber,tol)
!!$       write(*,*) n, real(sh), real(st), real(sp)
!!$       deallocate(likes)
!!$    end do
!!$  end block

  block
    integer :: i
    real(REAL64) :: a,b,tol
    write(*,*) 'Minimum of the scale for N(0,1) in Hampel:'
    do i = 940,950
       s = 1e-3_REAL64*i
       write(*,*) s,fun_gt(s)
    end do

    a = 0.9
    b = 1.0
    tol = 1e-7
    write(*,*) 'sqrt(nu) by fmin:',zeroin(a,b,fun_gt,tol)

  end block

contains

!!$  function negentropy_gauss(x)
!!$
!!$    real(kind) :: negentropy_gauss
!!$    real(kind), intent(in) :: x
!!$
!!$    negentropy_gauss = - x**2/2*exp(-x**2/2)
!!$
!!$  end function negentropy_gauss
!!$
!!$  function negentropy_huber(x)
!!$
!!$    real(kind) :: negentropy_huber
!!$    real(kind), intent(in) :: x
!!$
!!$    negentropy_huber = - ihuber(x)*exp(-ihuber(x))
!!$
!!$  end function negentropy_huber
!!$
!!$  function negentropy_tukey(x)
!!$
!!$    real(kind) :: negentropy_tukey
!!$    real(kind), intent(in) :: x
!!$
!!$    negentropy_tukey = - itukey(x)*exp(-itukey(x))
!!$
!!$  end function negentropy_tukey

  function likes_huber(s) result(f)

    real(kind), intent(in) :: s
    real(kind) :: f

    f = sum(ihuber(likes/s))/size(likes) + log(s)

  end function likes_huber

  function likes_pshuber(s) result(f)

    real(kind), intent(in) :: s
    real(kind) :: f

    f = sum(ipshuber(likes/s))/size(likes) + log(s)

  end function likes_pshuber

  function likes_tukey(s) result(f)

    real(kind), intent(in) :: s
    real(kind) :: f

    f = sum(itukey(likes/s))/size(likes) + log(s)

  end function likes_tukey

  function beta_huber() result(beta)

    real(kind) :: x, s, snorm, beta
    integer :: i

    s = 0d0
    snorm = 0!2.658945
    do i = 50000, 0, -1
       x = i / 1d4
!       s = s + huber(x)**2*exp(-ihuber(x))
       s = s + huber(x)**2*exp(-ihuber(x))
       snorm = snorm + exp(-x**2/2)  !exp(-ihuber(x))
    end do
!    snorm = 2 * snorm
!    s = 2 * s
    beta = (s  / snorm)
    snorm = 5 / 50000.0
    beta = 2*s / 2.5066282746310 * snorm
    beta = 2*s / 2.65957459 / 1d4

  end function beta_huber

  function beta_pshuber() result(beta)

    real(kind) :: x, s, snorm, beta
    integer :: i

    s = 0d0
    snorm = 0!2.658945
    do i = 100000, 0, -1
       x = i / 1d4
       s = s + pshuber(x)**2*exp(-ihuber(x))
       snorm = snorm + exp(-x**2/2)  !exp(-ihuber(x))
    end do
!    snorm = 2 * snorm
!    s = 2 * s
    beta = (s  / snorm)
    snorm = 5 / 50000.0
    beta = 2*s / 2.5066282746310 * snorm
    beta = 2*s / 2.72611755 / 1d4

  end function beta_pshuber


  function beta_tukey() result(beta)

    real(kind) :: x, s, snorm, beta
    integer :: i

    s = 0
    snorm = 0
    do i = 60000, 0, -1
       x = i / 1d4
       s = s + tukey(x)**2*exp(-itukey(x))
!       s = s + tukey(x)**2*exp(-tukey(x)**2/2)
!       snorm = snorm + exp(-tukey(x)**2/2)
!       s = s + itukey(x)*exp(-itukey(x))
!       snorm = snorm + exp(-itukey(x))
!       s = s + x**2*exp(-x**2/2)
!       snorm = snorm + exp(-x**2/2)
        snorm = snorm + 1d-3
    end do
!    snorm = sqrt(2*pi)
    beta = s  / snorm !/ 1e5
    snorm = 5 / 50000.0
!    beta = sqrt(2*s / 2.5066282746310 * snorm)
    beta = 2*s / 2.6496658716877421 * snorm
    beta = 2*s / 2.61992 * snorm

  end function beta_tukey

  function beta_hampel() result(beta)

    real(kind) :: x, s, beta
    integer :: i

    s = 0
    do i = 85000, 0, -1
       x = i / 1d4
       s = s + hampel(x)**2*exp(-ihampel(x))
    end do
    beta = 2*s / 2.562175 * 1d-4

  end function beta_hampel

  function fun_beta_pshampel() result(beta)

    real(kind) :: x, s, beta, c, y, t, dx, f0
    integer :: i
    integer, parameter :: M = 87500000

    s = 0
    c = 0
    f0 = exp(-ipshampel(9.0_REAL64))
    f0 = 0
    dx = 8.75_REAL64 / M
!    do i = 0, M
    do i = M, 0, -1
       x = i * dx
!       y = pshampel(x)**2*exp(-x**2/2) - c
!       y = pshampel(x)**2*exp(-x**2/2) - c
       y = pshampel(x)**2*(exp(-ipshampel(x)) -f0) - c
!       y = dpshampel(x)*(exp(-ipshampel(x)) -f0) - c
       t = s + y
       c = (t - s) - y
       s = t
!       s = s + pshampel(x)**2*exp(-ipshampel(x))
    end do
    beta = 2*s *dx / gamma_pshampel
!    beta = 2*s * dx / sqrt(2*pi)

  end function fun_beta_pshampel

  function fun_beta_pshampel1() result(beta)

    real(kind) :: x, s, beta, c, y, t, dx
    integer :: i
    integer, parameter :: M = 910000000

    s = 0
    c = 0
    dx = real(9.1,REAL64) / M
    do i = -M, M
       x = i * dx
!       y = dpshampel(x)*exp(-x**2/2) - c
       y = pshampel(x)**2*exp(-x**2/2) - c
       t = s + y
       c = (t - s) - y
       s = t
    end do
    beta = s / sqrt(2*pi) * dx

  end function fun_beta_pshampel1


  function nu_hampel() result(nu)

    real(kind) :: x, s, nu, h
    integer :: i

    s = 0
    h = 1e-6
    do i = 9100000, 0, -1
       x = i * h
       s = s + dhampel(x)*exp(-ihampel(x))
!       s = s + dhampel(x)*exp(-x**2/2)
    end do
    nu = 2*s / sqrt(2*pi) * h
!    nu = 2*s / 2.564 * h

  end function nu_hampel

  function fun_nu_pshampel() result(nu)

    real(kind) :: x, s, nu, c, y, t, dx, expc, hc
    integer :: i
    integer, parameter :: M = 87500000

    hc = 8.75_REAL64
    expc = exp(-ipshampel(hc))
    s = 0
    c = 0
    dx = hc / M
    do i = M, 0, -1
       x = i * dx
       y = dpshampel(x)*exp(-x**2/2) - c
       !       y = dpshampel(x)*exp(-ipshampel(x)) - c
!       y = dpshampel(x)*(exp(-ipshampel(x)) - expc) - c
       t = s + y
       c = (t - s) - y
       s = t
    end do
   nu = 2*s / sqrt(2*pi) * dx
!    nu = 2*s * dx / gamma_pshampel

  end function fun_nu_pshampel


  function beta_dhuber() result(beta)

    real(kind) :: x, s, snorm, beta
    integer :: i

    s = 0d0
    snorm = 0!2.658945
    do i = 50000, 0, -1
       x = i / 1d4
       s = s + dhuber(x)*exp(-ihuber(x))
       snorm = snorm + exp(-x**2/2)  !exp(-ihuber(x))
    end do
!    snorm = 2 * snorm
!    s = 2 * s
    beta = (s  / snorm)
    snorm = 5 / 50000.0
    beta = 2*s / 2.5066282746310 * snorm

  end function beta_dhuber

  function beta_dpshuber() result(beta)

    real(kind) :: x, s, snorm, beta
    integer :: i

    s = 0d0
    snorm = 0!2.658945
    do i = 50000, 0, -1
       x = i / 1d4
       s = s + dpshuber(x)*exp(-ihuber(x))
       snorm = snorm + exp(-x**2/2)  !exp(-ihuber(x))
    end do
!    snorm = 2 * snorm
!    s = 2 * s
    beta = (s  / snorm)
    snorm = 5 / 50000.0
    beta = 2*s / 2.5066282746310 * snorm

  end function beta_dpshuber


  function beta_dtukey() result(beta)

    real(kind) :: x, s, snorm, beta
    integer :: i

    s = 0
    snorm = 0
    do i = 60000, 0, -1
       x = i / 1d4
!       s = s + tukey(x)**2*exp(-itukey(x))
!       s = s + dtukey(x)*exp(-x**2/2)
       s = s + dtukey(x)*exp(-itukey(x))
!       snorm = snorm + exp(-tukey(x)**2/2)
!       s = s + itukey(x)*exp(-itukey(x))
!       snorm = snorm + exp(-itukey(x))
!       s = s + x**2*exp(-x**2/2)
!       snorm = snorm + exp(-x**2/2)
        snorm = snorm + 1d-3
    end do
!    snorm = sqrt(2*pi)
    beta = sqrt(s  / snorm) !/ 1e5
    snorm = 5 / 50000.0
    beta = 2*s / 2.6496658716877421 * snorm

  end function beta_dtukey


  function norm_huber() result(gama)

    real(kind) :: x, s, gama
    integer :: i

    s = 0
    do i = 0, 18000000
       x = i / 1d6
       s = s + exp(-ihuber(x))
!       s = s + exp(-x**2/2)
    end do
    gama = 2 * s / 1d6

  end function norm_huber

  function norm_pshuber() result(gama)

    real(kind) :: x, s, gama
    integer :: i

    s = 0
    do i = 0, 18000000
       x = i / 1d6
       s = s + exp(-ipshuber(x))
!       s = s + exp(-x**2/2)
    end do
    gama = 2 * s / 1d6

  end function norm_pshuber

  function norm_pshampel() result(gama)

    real(kind) :: x, s, gama, c, y, t, dx, f0
    integer :: i
    integer, parameter :: M = 8750000

!    f0 = exp(-ipshampel(100.0_REAL64))
    f0 = exp(-ipshampel(8.75_REAL64))
!    f0 = 0
    s = 0
    c = 0
!    dx = real(1e6,REAL64)
    dx = real(8.5,REAL64) / M
    dx = 8.75_REAL64 / M
!    do i = 0, M
    do i = M,0,-1
       x = i * dx
       y = (exp(-ipshampel(x)) - f0) - c
!       y = (exp(-ipshampel(x)) - 4*f0) - c
!       y = exp(-ipshampel(x)) - c
       t = s + y
       c = (t - s) - y
       s = t
!       s = s + exp(-ipshampel(x))
    end do
    gama = 2 * s * dx !- 2*8.75*exp(-ipshampel(8.75_REAL64))

  end function norm_pshampel


  function norm_hampel() result(gama)

    real(kind) :: x, s, gama
    integer :: i

    s = 0
    do i = 0, 850000
       x = i / 1e5
       s = s + exp(-ihampel(x))
    end do
    gama = 2 * s / 1e5

  end function norm_hampel



  function norm_tukey() result(gama)

    real(kind) :: x, s, gama
    integer :: i

    s = 0
    do i = 0, 600000
       x = i / 1e5
       s = s + (exp(-itukey(x)) - exp(-itukey(6.0001_REAL64))*0)
!       s = s + exp(-tukey(x)**2/2)
!       s = s + exp(-x**2/2)
    end do
    gama = 2 * s / 1e5

  end function norm_tukey


  function fun_tukey(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    if( abs(x) <= 6 ) then
       f = exp(-itukey(x)) - exp(-itukey(6.0001_REAL64))*0 !/ 3.73253078 !/ 2.6496858221126
!       f = exp(-x**2/2) !/ sqrt(2*pi)
!       f = exp(-tukey(x)**2/2)
    else
       f = 0
    end if

  end function fun_tukey

  function fun_huber(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

!    if( abs(x) <= 6 ) then
       f = exp(-ihuber(x)) !/ 3.73253078 !/ 2.6496858221126
!       f = exp(-x**2/2) !/ sqrt(2*pi)
!       f = exp(-tukey(x)**2/2)
!    else
!       f = 0
!    end if

  end function fun_huber

  function fun_pshuber(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    f = exp(-ipshuber(x))

  end function fun_pshuber

  function fun_pshampel(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    f = exp(-ipshampel(x)) !- exp(-ipshampel(100.0_REAL64))
    return

    if( abs(x) <= 8.5_REAL64 ) then
       f = exp(-ipshampel(x))
    else
       f = 0
    end if

  end function fun_pshampel

  function fun_hampel(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    f = exp(-ihampel(x))

  end function fun_hampel


  function sig_tukey(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    if( abs(x) <= 6 ) then
!       f = itukey(x)*exp(-itukey(x)) / 3.73253078
!       f = tukey(x)**2*exp(-itukey(x)) / 3.73253078
!       f = x**2*exp(-x**2/2) / sqrt(2*pi)
!       f = tukey(x)**2*exp(-x**2/2) / 2.50662827463100
       f = tukey(x)**2*exp(-itukey(x)) / 2.64968582
    else
       f = 0
    end if
!    f = x**2* exp(-x**2/2) / sqrt(2*pi)
    f = tukey(x)**2*exp(-itukey(x)) / 2.6496558604796374
    f = tukey(x)**2*exp(-itukey(x)) / 2.61992
!    f = tukey(x)**2*exp(-x**2/2) / 2.50662827463100

  end function sig_tukey

  function sig_hampel(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

!    if( abs(x) <= 8.5 ) then
!       f = hampel(x)**2*exp(-ihampel(x)) / 2.562175
!    else
!       f = 0
!    end if
    f = hampel(x)**2*exp(-ihampel(x)) / 2.562175

  end function sig_hampel

  function sig_pshampel(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x
    real(kind) :: w

    w = exp(-ipshampel(x)) !- exp(-ipshampel(100.0_REAL64))
    w = exp(-x**2/2)/sqrt(2*pi)
    f = pshampel(x)**2*w! / gamma_pshampel
    !    f = pshampel(x)**2*exp(-x**2/2) / sqrt(2*pi)
    f = pshampel(x/scale)**2*exp(-x**2/2)/sqrt(2*pi)
    return

    if( abs(x) <= 8.5_REAL64 ) then
       f = pshampel(x)**2*w / gamma_pshampel
    else
       f = 0
    end if
!    f = hampel(x)**2*exp(-ihampel(x)) / 2.562175

  end function sig_pshampel

  function sig_nu_hampel(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    f = dhampel(x)*exp(-x**2/2) / sqrt(2*pi)

  end function sig_nu_hampel

  function sig_nu_pshampel(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    f = dpshampel(x)*exp(-x**2/2) / sqrt(2*pi)
!    f = dpshampel(x)*exp(-ipshampel(x)) / gamma_pshampel
!    f = exp(-x**2/2) / sqrt(2*pi)

    f = dpshampel(x/scale)*exp(-x**2/2) / sqrt(2*pi)

  end function sig_nu_pshampel


  function sig_huber(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    f = huber(x)**2 * exp(-ihuber(x)) / 2.65957459
!    f = huber(x)**2 * exp(-x**2/2) / 2.50662827463100
!    f = x**2 *  exp(-x**2/2) / 2.50662827463100
!    f = exp(-ihuber(x)) / 2.658945

  end function sig_huber

  function sig_pshuber(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    f = pshuber(x)**2 * exp(-ipshuber(x)) / 2.7261175497405246

  end function sig_pshuber

  function sig_dtukey(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    if( abs(x) <= 6 ) then
!       f = itukey(x)*exp(-itukey(x)) / 3.73253078
!       f = tukey(x)**2*exp(-itukey(x)) / 3.73253078
!       f = x**2*exp(-x**2/2) / sqrt(2*pi)
!       f = tukey(x)**2*exp(-x**2/2) / 2.50662827463100
       f = dtukey(x)*exp(-itukey(x)) / 2.64968582
    else
       f = 0
    end if
!    f = x**2* exp(-x**2/2) / sqrt(2*pi)
    f = dtukey(x)*exp(-itukey(x)) / 2.6496558604796374
!    f = tukey(x)**2*exp(-x**2/2) / 2.50662827463100

  end function sig_dtukey

  function sig_dhuber(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    f = dhuber(x) * exp(-ihuber(x)) / 2.6578194259790333
!    f = huber(x)**2 * exp(-x**2/2) / 2.50662827463100
!    f = x**2 *  exp(-x**2/2) / 2.50662827463100
!    f = exp(-ihuber(x)) / 2.658945

  end function sig_dhuber

  function sig_dpshuber(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

    f = dpshuber(x) * exp(-ipshuber(x)) / 2.7261175497405246

  end function sig_dpshuber



  function rtest(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x

!    f = sin(x)
    f = exp(-x**2)*1.1283791670955126_REAL64   !*1.1283791

  end function rtest

  function beta_sqr_norm(n) result(s)

    use noise_REAL64

    integer, intent(in) :: n
    real(kind) :: s
    real(kind), dimension(:), allocatable :: x

    allocate(x(n))
    call gnoise_fill(0.0_kind,1.0_kind,x)
!    s = sqrt(sum(x**2) / n)
    s = sum(x**2) / n

  end function beta_sqr_norm

  function beta_wins_norm(n) result(s)

    use noise_REAL64

    integer, intent(in) :: n

    real(kind), parameter :: rmin = real(-1.618,kind)
    real(kind), parameter :: rmax = real( 1.618,kind)
    real(kind) :: s
    real(kind), dimension(:), allocatable :: x

    allocate(x(n))
    call gnoise_fill(0.0_kind,1.0_kind,x)
    x = max(rmin,min(x,rmax)) ! winsorizing
!    s = sqrt(sum(x**2) / n ) / 0.9068012 ! Phi(rmax)-Phi(rmin)
    s = sum(x**2) / n / 0.9068012**2 ! Phi(rmax)-Phi(rmin)

  end function beta_wins_norm

  function beta_huber_norm(n) result(s)

    use noise_REAL64

    integer, intent(in) :: n
    real(kind) :: s,s0
    real(kind), dimension(:), allocatable :: x
    integer :: iter


    allocate(x(n))
    call gnoise_fill(0.0_kind,1.0_kind,x)

    s0 = -1
    s = 1
    do iter = 1, 2*precision(s)
       s0 = s
       s = sqrt(sum(huber(x/s)**2)*s**2 / n / 0.843849785**2)
!       s = 2* sum(ihuber(x)) / n / 0.9326
!       s = sum(huber(x/s)**2) / n
       if( abs(s - s0) < 0.1/sqrt(real(n)) ) exit
!       write(*,*) iter,s,s*sqrt(sum(huber(x/s)**2)/sum(dhuber(x/s))**2*n)
    end do

  end function beta_huber_norm

  function beta_pshuber_norm(n) result(s)

    use noise_REAL64

    integer, intent(in) :: n
    real(kind) :: s,s0
    real(kind), dimension(:), allocatable :: x
    integer :: iter


    allocate(x(n))
    call gnoise_fill(0.0_kind,1.0_kind,x)

    s0 = -1
    s = 1
    do iter = 1, 2*precision(s)
       s0 = s
       s = sqrt(sum(pshuber(x/s)**2)*s**2 / n / 0.793123038**2)
!       s = 2* sum(ipshuber(x)) / n / 0.793123038**2
!       s = sum(huber(x/s)**2) / n
       if( abs(s - s0) < 0.1/sqrt(real(n)) ) exit
!       write(*,*) iter,s,s*sqrt(sum(huber(x/s)**2)/sum(dhuber(x/s))**2*n)
    end do

  end function beta_pshuber_norm

  function beta_tukey_norm(n) result(s)

    use noise_REAL64

    integer, intent(in) :: n
    real(kind) :: s,s0
    real(kind), dimension(:), allocatable :: x
    integer :: iter


    allocate(x(n))
    call gnoise_fill(0.0_kind,1.0_kind,x)

    s0 = -1
    s = 1
    do iter = 1, 2*precision(s)
       s0 = s
       s = sqrt(sum(tukey(x/s)**2)*s**2 / n / 0.8530367413**2)
       if( abs(s - s0) < 0.1/sqrt(real(n)) ) exit
!       write(*,*) iter,s,s*sqrt(sum(tukey(x/s)**2)/sum(dtukey(x/s))**2*n)
    end do

  end function beta_tukey_norm

  function fun_gamma_pshampel() result(gamma)

    real(kind) :: x, s, gamma, c, y, t, dx, w,u
    integer :: i
    integer, parameter :: M = 875000000

    s = 0
    c = 0
    dx = real(8.75,REAL64) / M
    do i = M, 0, -1
       x = i * dx
!       w = exp(-ipshampel(x)) - exp(-ipshampel(100.0_REAL64))
!       w = exp(-ipshampel(x))
       !       w = exp(-x**2/2)
       !y = dpshampel(x)*(x**2)*w - c
       u = x / scale
       y = dpshampel(u)*(u**2)*exp(-x**2/2) - c
       t = s + y
       c = (t - s) - y
       s = t
!       s = s + dpshampel(x)*(x**2)*exp(-ipshampel(x))
    end do
!    gamma = 2*s * dx / gamma_pshampel !/ beta_pshampel
    gamma = 2*s / sqrt(2*pi) * dx

  end function fun_gamma_pshampel

!!$  function gamma1_pshampel() result(gamma)
!!$
!!$    use rfun_hampel
!!$    real(kind) :: x, s, gamma, c, y, t, dx
!!$    integer :: i
!!$    integer, parameter :: M = 850000000
!!$
!!$    s = 0
!!$    c = 0
!!$    dx = real(8.5,REAL64) / M
!!$    do i = M, 0, -1
!!$       x = i * dx
!!$       y = exp(-ipshampel(x)) - c
!!$       t = s + y
!!$       c = (t - s) - y
!!$       s = t
!!$!       s = s + exp(-ipshampel(x))
!!$    end do
!!$    gamma = 2*s / gamma_pshampel * dx
!!$
!!$  end function gamma1_pshampel


!!$  function nu_pshampel() result(nu)
!!$
!!$    use rfun_hampel
!!$    real(kind) :: x, s, nu
!!$    integer :: i
!!$
!!$    s = 0
!!$    do i = 85000, 0, -1
!!$       x = i / 1d4
!!$       s = s + dpshampel(x)*exp(-ipshampel(x))
!!$    end do
!!$    nu = 2*s / gamma_pshampel * 1d-4
!!$
!!$  end function nu_pshampel

  function off_pshampel() result(off)

    real(kind) :: x, s, off
    integer :: i

    s = 0
    do i = -85000, 0, 1
       x = i / 1d4
!       s = s + dpshampel(x)*x*exp(-ipshampel(x))
       s = s + dpshampel(x)*x*exp(-x**2/2)
    end do
    do i = 85000, 0, -1
       x = i / 1d4
!       s = s + dpshampel(x)*x*exp(-ipshampel(x))
       s = s + dpshampel(x)*x*exp(-x**2/2)
    end do
!    off = s / gamma_pshampel * 1d-4
    off = s / sqrt(2*pi) * 1d-4

  end function off_pshampel


  function fun_psi2(x) result(psi2)

    real(kind) :: psi2
    real(kind), intent(in) :: x
    real(kind) :: w

    w = exp(-ipshampel(x))! - exp(-ipshampel(1000.0_REAL64))
    psi2 = pshampel(x)**2*w / gamma_pshampel
    !psi2 = pshampel(x)**2*exp(-x**2/2) / sqrt(2*pi)

  end function fun_psi2

  function fun_dpsi(x) result(dpsi)

    real(kind) :: dpsi
    real(kind), intent(in) :: x
    real(kind) :: w

    w = exp(-ipshampel(x))! - exp(-ipshampel(1000.0_REAL64))
    dpsi = dpshampel(x)*w / gamma_pshampel
!    dpsi = dpshampel(x)*exp(-x**2/2) / sqrt(2*pi)

  end function fun_dpsi

  function fun_s(x) result(f)

    real(kind) :: f
    real(kind), intent(in) :: x
    real(kind) :: w

!    w = exp(-ipshampel(x)) - exp(-ipshampel(1000.0_REAL64))
    w = exp(-ipshampel(x))
    f = pshampel(x)*x*w / gamma_pshampel

  end function fun_s

  real(REAL64) function fun_gt(scale)

    real(kind), intent(in) :: scale
    real(kind) :: x, s, c, y, t, dx
    integer :: i
    integer, parameter :: M = 87500000

    s = 0
    c = 0
    dx = real(8.75,REAL64) / M
    do i = M,0,-1
!    do i = M, 0, -1
       x = i * dx
       y = pshampel(x/scale)*(x/scale)*exp(-x**2/2) - c
!       y = pshampel(x/scale)*(x/scale)*exp(-ipshampel(x)) - c
!       y = dpshampel(x)*(exp(-ipshampel(x)) -f0) - c
       t = s + y
       c = (t - s) - y
       s = t
!       s = s + pshampel(x)**2*exp(-ipshampel(x))
    end do
    fun_gt = 2 * s * dx / sqrt(2*pi) - 1
!    fun_gt = 2 * s * dx / gamma_pshampel - 1

  end function fun_gt

end program eparam
