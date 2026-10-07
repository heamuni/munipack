
program newparam

  use fmm
  use iromberg
  use oakleaf
  use iso_fortran_env

  implicit none
  real(REAL64), parameter :: pi = 4.0_REAL64*atan(1.0_REAL64)
  real(REAL64), parameter :: c = 8.8_REAL64
  real(REAL64) :: scale

  block
    real(REAL64) :: a,b,tol

    a = 0.9
    b = 1.0
    tol = 1e-15
    scale = zeroin(a,b,eq_gt,tol)
    write(*,*) 'scale by (4b):',scale
  end block

  block
    use noise_REAL64

    real(REAL64), parameter :: t = 0
    real(REAL64), parameter :: sig = 1
    real(REAL64), dimension(:), allocatable :: data
    real(REAL64) :: f,s,fs,tol,a
    integer :: k,n

    s = scale
    do k = 3,9
       n = 10**k
       allocate(data(n))
       call gnoise_fill(t,sig,data)
       f = sum(dpshampel(data/s)) / n
       write(*,*) n,'  ν = ',f,sqrt(f)-s
       deallocate(data)
    end do

    a = 0
    tol = 1e-15
    fs = 2*romberg(dpsi,a,c,tol)
    write(*,*) " E[ψ'] = ",fs,fs-f
  end block

  block
    use noise_REAL64

    real(REAL64), parameter :: t = 0
    real(REAL64), parameter :: sig = 1
    real(REAL64), dimension(:), allocatable :: data
    real(REAL64) :: f,s,a,tol,fs
    integer :: k,n

    s = scale
    do k = 3,9
       n = 10**k
       allocate(data(n))
       call gnoise_fill(t,sig,data)
       f = sum(pshampel(data/s)**2) / n
       write(*,*) n,'  β = ',f
       deallocate(data)
    end do

    a = 0
    tol = 1e-15
    fs = 2*romberg(psi2,a,c,tol)
    write(*,*) " E[ψ²] = ",fs,fs-f
  end block

  block
    use noise_REAL64

    real(REAL64), parameter :: t = 0
    real(REAL64), parameter :: sig = 1
    real(REAL64), dimension(:), allocatable :: data
    real(REAL64) :: f,s,a,tol,fs
    integer :: k,n

    s = scale
    do k = 3,9
       n = 10**k
       allocate(data(n))
       call gnoise_fill(t,sig,data)
       f = sum(dpshampel(data/s)*(data/s)**2) / n
       write(*,*) n,"  ϖ = ",f
       deallocate(data)
    end do

    a = 0
    tol = 1e-15
    fs = 2*romberg(dpsir2,a,c,tol)
    write(*,*) " E[ψ'.r] = ",fs,fs-f
  end block

  write(*,*) 'rho_c=',ipshampel(c)

contains

  real(REAL64) function eq_gt(s)

    real(REAL64), intent(in) :: s
    real(REAL64) :: tol, f

    tol = 1e-15
    scale = s
    f = romberg(psir,-c,c,tol)

    eq_gt = 1 - f

  end function eq_gt

  real(REAL64) function psir(x)

    real(REAL64), intent(in) :: x
    real(REAL64), parameter :: norm = sqrt(2.0_REAL64*pi)

    psir = pshampel(x/scale)*(x/scale)*exp(-x**2/2) / norm

  end function psir

  real(REAL64) function dpsi(x)

    real(REAL64), intent(in) :: x
    real(REAL64), parameter :: norm = sqrt(2.0_REAL64*pi)

    dpsi = dpshampel(x/scale)*exp(-x**2/2) / norm

  end function dpsi

  real(REAL64) function psi2(x)

    real(REAL64), intent(in) :: x
    real(REAL64), parameter :: norm = sqrt(2.0_REAL64*pi)

    psi2 = pshampel(x/scale)**2*exp(-x**2/2) / norm

  end function psi2

  real(REAL64) function dpsir2(x)

    real(REAL64), intent(in) :: x
    real(REAL64), parameter :: norm = sqrt(2.0_REAL64*pi)

    dpsir2 = dpshampel(x/scale)*(x/scale)**2*exp(-x**2/2) / norm

  end function dpsir2


end program newparam
