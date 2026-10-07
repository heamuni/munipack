!
!  Light curves
!
!  Copyright © 2017-2026 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Munipack.
!
!  Munipack is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Munipack is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with Munipack.  If not, see <http://www.gnu.org/licenses/>.


module lcurves

  use iso_fortran_env

  implicit none
  private

  real(REAL64), parameter :: pi = 3.1415926535897931_REAL64

  type, public :: LightCurve
     real(REAL64) :: lcmag, lcamp, lcper, lcjd0
     real(REAL64), dimension(2) :: coo
   contains
     procedure :: flux => wave
     procedure :: magvar
  end type LightCurve

  type, extends(LightCurve), public :: LightCurveDeltaCep
   contains
     procedure :: flux => delta_cep
  end type LightCurveDeltaCep

  type, extends(LightCurve), public :: LightCurveAlgol
   contains
     procedure :: flux => beta_per
  end type LightCurveAlgol

  type, extends(LightCurve), public :: LightCurveFourier
     ! Fourier's coefficients
     complex(REAL64), dimension(:), allocatable, private :: fcoeff
   contains
     procedure :: SetFourCoefficients
     procedure :: flux => four_flux
  end type LightCurveFourier

  type SplineSerieType
     ! coeficients of the spline
     real(REAL64), dimension(:), allocatable, private :: time
     real(REAL64), dimension(:,:), allocatable, private :: taaux
   contains
     procedure :: Init, Spline
  end type SplineSerieType

  type, extends(LightCurve), public :: LightCurveTabular
     type(SplineSerieType) :: lcspline
   contains
     procedure :: SetLCTable
     procedure :: flux => spline_in_table
  end type LightCurveTabular

  public :: wave, delta_cep, beta_per, four_flux, spline_in_table, magvar, &
       SetFourCoefficients, SetLCTable

contains

  real(REAL64) function wave(this,phi)

    class(LightCurve), intent(in) :: this
    real(REAL64), intent(in) :: phi

    wave = cos(2*pi*phi)

  end function wave

  ! delta Cep

  real(REAL64) function delta_cep(this,phi)

    use lcparams

    class(LightCurveDeltaCep), intent(in) :: this
    real(REAL64), intent(in) :: phi

    delta_cep = fourfun(dcep,phi)

  end function delta_cep

  real(REAL64) function beta_per(this,phi)

    use lcparams

    class(LightCurveAlgol), intent(in) :: this
    real(REAL64), intent(in) :: phi

    beta_per = fourfun(algol,phi)

  end function beta_per

  pure function fourfun(c,phi)

    complex(REAL64), dimension(0:), intent(in) :: c
    real(REAL64), intent(in) :: phi
    real(REAL64) :: fourfun,t
    complex(REAL64) :: s,u,w
    integer :: j

    t = 2*pi*phi
    w = cmplx(cos(t),-sin(t),REAL64)
    u = cmplx(1.0_REAL64,0.0_REAL64,REAL64)
    s = 0.0_REAL64
    do j = 0,size(c)-1
       s = s + c(j)*u
       u = u * w
    end do
    fourfun = real(s)

  end function fourfun

  real(REAL64) pure function four_flux(this,phi)

    class(LightCurveFourier), intent(in) :: this
    real(REAL64), intent(in) :: phi

    four_flux = fourfun(this%fcoeff,phi)

  end function four_flux

  subroutine Init(this,x,y,dy)

    ! https://en.wikipedia.org/wiki/Smoothing_spline

    class(SplineSerieType) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), intent(in) :: y,dy

    real(REAL64), external :: smooth
    real(REAL64), dimension(:,:), allocatable :: aux, aaux

    integer :: n
    real(REAL64) :: w,s

    n = size(x)
    allocate(aux(n,7),aaux(n,4))
    allocate(this%time(n),this%taaux(4,n))

    ! interpolation
    s = epsilon(x)
    ! aproximation
    s = n + sqrt(2.0*n)

    w = smooth(x,y,dy,n,s,aux,aaux)
    this%taaux = transpose(aaux(:,1:4))
    this%time = x

    ! block
    ! integer :: i
    !    open(1,file='l')
    !    do i = -100,200
    !       write(1,*) i/100.0,real(lcspline(i/100.0_REAL64))
    !    end do
    !    close(1)
    ! end block

  end subroutine Init

  real(REAL64) function Spline(this,t)

    class(SplineSerieType), intent(in) :: this
    real(REAL64), intent(in) :: t
    integer :: n1
    real(REAL64), external :: ppvalu

    n1 = size(this%time) - 1
    Spline = ppvalu(this%time,this%taaux,n1,4,t,0)

  end function Spline

  real(REAL64) function magvar(this,jd)

    use astrosphere, only: phase

    class(LightCurve), intent(in) :: this
    real(REAL64), intent(in) :: jd
    real(REAL64) :: f

    f = 10**(-0.4 * this%lcamp * this%flux(phase(jd,this%lcjd0,this%lcper)))
    magvar = this%lcmag + 2.5_REAL64*log10(f)

  end function

  real(REAL64) function spline_in_table(this, phi)

    class(LightCurveTabular), intent(in) :: this
    real(REAL64), intent(in) :: phi

    spline_in_table = this%lcspline%Spline(phi)

  end function spline_in_table

  subroutine SetLCTable(this,time,flux,dflux)

    class(LightCurveTabular) :: this
    real(REAL64), dimension(:), intent(in) :: time, flux, dflux

    call this%lcspline%Init(time,flux,dflux)

  end subroutine SetLCTable

  subroutine SetFourCoefficients(this,fcoeff)

    class(LightCurveFourier) :: this
    complex(REAL64), dimension(:), intent(in) :: fcoeff

    allocate(this%fcoeff,source=fcoeff)

  end subroutine SetFourCoefficients

end module lcurves
