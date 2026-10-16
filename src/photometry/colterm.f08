!
!  Colour terms of the photometric calibration
!
!  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
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
!
!
!  A passband of an instrument differs from the passband of the catalogue
!  (or of the standard system), so the ratio of the expected photons and
!  the counts of a star depends on its colour. In magnitudes,
!
!      m_cat - m_inst = ZP + k1 (c - cref) + k2 (c - cref)**2,
!
!  c a colour index (a difference of two magnitudes), cref a reference
!  colour. Munipack calibrates the counts to photons by the ratio
!  t = ph/ct (CTPH); the same relation is
!
!      ln(ph/ct) = b0 + b1 (c - cref) + b2 (c - cref)**2,
!
!  t = exp(b0) the calibration of a star of the reference colour,
!  b1 = -k1/magln, b2 = -k2/magln (magln = 2.5/ln 10). The photons of
!  a star are ph = t exp(g(c)) ct, g(c) = b1 (c - cref) + b2 (c - cref)**2,
!  so every quantity derived from the photons (MAG, FLUX, ...) has
!  the colour term of the magnitudes.
!
!  The coefficients are fitted together with the excess scatter of
!  the stars over their errors, robustly to outliers, by mixture_fit
!  of exscatter.f08 (the likelihood of a mixture, the same as for
!  the zero point alone): the variance of a star includes the error of
!  its colour propagated by the slope g' = dg/dc (the effective variance)
!  and its covariance with the catalogue magnitude of the calibration
!  (the colour B-V has the error of V).
!

module colterm

  use iso_fortran_env

  implicit none
  private

  integer, parameter :: dbl = REAL64

  ! magnitudes per a unit of the natural logarithm of a flux, 2.5 / ln(10)
  real(dbl), parameter, public :: magln = 1.0857362047581294_dbl

  ! the reference colour is the mean colour of the calibration stars
  real(dbl), parameter, public :: cref_mean = huge(1.0_dbl)

  ! models of the colour dependence:
  ! COLT_NONE no colour term,
  ! COLT_CAT  the colour term of the calibration fitted on the catalogue
  !           colour (phcal --color-term), applied by the catalogue colour,
  ! COLT_TRA  the magnitude transformation of phfotran (--model mag),
  !           its coefficients fixed, applied by the colour of the
  !           calibrated magnitudes of the frames (phcal --tratab)
  integer, parameter, public :: COLT_NONE = 0, COLT_CAT = 1, COLT_TRA = 2

  type, public :: colour_term
     integer :: model = COLT_NONE
     ! the colour of every band: the names of the catalogue columns
     ! (COLT_CAT) or of the filters (COLT_TRA), c = A - B
     character(len=68), dimension(:), allocatable :: cola, colb
     ! COLT_TRA: the indices of the filters of the colour of every band
     integer, dimension(:), allocatable :: ia, ib
     ! the order of the polynomial, the number of the stars used
     integer, dimension(:), allocatable :: order, nstars
     ! ln(ph/ct) = beta(0) + beta(1) d + beta(2) d**2, d = c - cref
     ! (beta(0,:) is ln CTPH), the covariance of beta (in the order 0,1,2)
     real(dbl), dimension(:,:), allocatable :: beta
     real(dbl), dimension(:,:,:), allocatable :: cov
     ! the reference colour; the mean colour of the calibration stars,
     ! its scatter, the mean colour correction g of the calibration stars
     ! and its scatter (used for the stars without a colour)
     real(dbl), dimension(:), allocatable :: cref, cmean, cstd, gmean, gstd
     ! COLT_CAT: the catalogue (all its stars): positions [deg], colours,
     ! their errors (< 0 for an undefined colour), the search radius [deg]
     real(dbl), dimension(:), allocatable :: ra, dec, col, dcol
     real(dbl) :: tol = 0
  end type colour_term

  public :: colterm_init, colour_corr, colour_solve

contains

  subroutine colterm_init(ctm,model,nbands)

    type(colour_term), intent(out) :: ctm
    integer, intent(in) :: model, nbands

    ctm%model = model
    allocate(ctm%cola(nbands),ctm%colb(nbands),ctm%ia(nbands),ctm%ib(nbands),&
         ctm%order(nbands),ctm%nstars(nbands),ctm%beta(0:2,nbands), &
         ctm%cov(3,3,nbands),ctm%cref(nbands),ctm%cmean(nbands), &
         ctm%cstd(nbands),ctm%gmean(nbands),ctm%gstd(nbands))
    ctm%cola = ''
    ctm%colb = ''
    ctm%ia = 0
    ctm%ib = 0
    ctm%order = 1
    ctm%nstars = 0
    ctm%beta = 0
    ctm%cov = 0
    ctm%cref = 0
    ctm%cmean = 0
    ctm%cstd = 0
    ctm%gmean = 0
    ctm%gstd = 0
    allocate(ctm%ra(0),ctm%dec(0),ctm%col(0),ctm%dcol(0))

  end subroutine colterm_init


  subroutine colour_corr(ctm,n,c,dc,g,dg)

    ! The colour correction g = ln(ph/ct) - ln CTPH of a star of the colour
    ! c (with the error dc, dc < 0: an undefined colour) in the band n,
    ! and its error dg: the error of the calibration at the colour (by
    ! the covariance of beta, the error of CTPH included) and the error of
    ! the colour. A star without a colour gets the mean correction of
    ! the calibration stars, its error inflated by their scatter.

    type(colour_term), intent(in) :: ctm
    integer, intent(in) :: n
    real(dbl), intent(in) :: c, dc
    real(dbl), intent(out) :: g, dg

    real(dbl), dimension(3) :: x
    real(dbl) :: d, gp, v

    if( dc >= 0 ) then
       d = c - ctm%cref(n)
       g = ctm%beta(1,n)*d + ctm%beta(2,n)*d**2
       gp = ctm%beta(1,n) + 2*ctm%beta(2,n)*d
       x = [ 1.0_dbl, d, d**2 ]
       v = dot_product(x,matmul(ctm%cov(:,:,n),x)) + (gp*dc)**2
    else
       g = ctm%gmean(n)
       v = ctm%cov(1,1,n) + ctm%gstd(n)**2
    end if
    dg = sqrt(max(v,0.0_dbl))

  end subroutine colour_corr


  subroutine colour_solve(ctm,m0,dm0,m,dm,c,dc,status)

    ! The magnitudes of a star in the standard bands by the magnitude
    ! transformation (COLT_TRA): m0 are the magnitudes calibrated by CTPH
    ! (the colour terms not applied) with their errors dm0 (the counts,
    ! CTPH and the systematic error), m the solution of
    !
    !     m(n) = m0(n) + k1(n) d(n) + k2(n) d(n)**2,
    !     d(n) = m(ia(n)) - m(ib(n)) - cref(n),  k = -magln beta,
    !
    ! by Newton's method (the colours of the standard magnitudes), dm
    ! the errors of m by the linearised solution, (I - G) dm = dm0,
    ! the errors of the coefficients included; c, dc the colours d(n) +
    ! cref(n) and their errors. status /= 0 when it does not converge.

    type(colour_term), intent(in) :: ctm
    real(dbl), dimension(:), intent(in) :: m0, dm0
    real(dbl), dimension(:), intent(out) :: m, dm, c, dc
    integer, intent(out) :: status

    real(dbl), dimension(size(m0),size(m0)) :: jm, ji, vm
    real(dbl), dimension(size(m0)) :: f, d, k1, k2, gp, s2, delta
    integer :: nb, n, iter

    nb = size(m0)
    k1 = -magln*ctm%beta(1,:)
    k2 = -magln*ctm%beta(2,:)
    m = m0
    status = 1
    do iter = 1, 50
       d = m(ctm%ia) - m(ctm%ib) - ctm%cref
       f = m - m0 - k1*d - k2*d**2
       gp = k1 + 2*k2*d
       call jacobian
       call gsolve(jm,-f,delta,status)
       if( status /= 0 ) exit
       m = m + delta
       if( maxval(abs(delta)) < 1e-12_dbl ) then
          status = 0
          exit
       end if
       status = 1
    end do

    d = m(ctm%ia) - m(ctm%ib) - ctm%cref
    gp = k1 + 2*k2*d
    call jacobian
    c = d + ctm%cref

    ! the errors: dm = (I - G)**-1 (dm0 + errors of k1, k2)
    s2 = dm0**2 + magln**2*(d**2*ctm%cov(2,2,:) + d**4*ctm%cov(3,3,:))
    ji = 0
    do n = 1, nb
       delta = 0
       delta(n) = 1
       call gsolve(jm,delta,ji(:,n),iter)
       if( iter /= 0 ) status = 1
    end do
    vm = matmul(ji,matmul(diag(s2),transpose(ji)))
    do n = 1, nb
       dm(n) = sqrt(max(vm(n,n),0.0_dbl))
       dc(n) = sqrt(max(vm(ctm%ia(n),ctm%ia(n)) + vm(ctm%ib(n),ctm%ib(n)) &
            - 2*vm(ctm%ia(n),ctm%ib(n)),0.0_dbl))
    end do

  contains

    subroutine jacobian
      ! the Jacobian of f: I - G, G(n,ia) = gp, G(n,ib) = -gp
      integer :: i
      jm = 0
      do i = 1, nb
         jm(i,i) = 1
      end do
      do i = 1, nb
         jm(i,ctm%ia(i)) = jm(i,ctm%ia(i)) - gp(i)
         jm(i,ctm%ib(i)) = jm(i,ctm%ib(i)) + gp(i)
      end do
    end subroutine jacobian

    function diag(v) result(a)
      real(dbl), dimension(:), intent(in) :: v
      real(dbl), dimension(size(v),size(v)) :: a
      integer :: i
      a = 0
      do i = 1, size(v)
         a(i,i) = v(i)
      end do
    end function diag

  end subroutine colour_solve


  subroutine gsolve(a,b,x,status)

    ! the solution of a x = b by the Gauss elimination with partial pivoting

    real(dbl), dimension(:,:), intent(in) :: a
    real(dbl), dimension(:), intent(in) :: b
    real(dbl), dimension(:), intent(out) :: x
    integer, intent(out) :: status

    real(dbl), dimension(size(b),size(b)+1) :: w
    real(dbl), dimension(size(b)+1) :: t
    integer :: n, i, k, ip

    n = size(b)
    w(:,1:n) = a
    w(:,n+1) = b
    x = 0
    status = 1
    do k = 1, n
       ip = k - 1 + maxloc(abs(w(k:n,k)),1)
       if( .not. (abs(w(ip,k)) > tiny(1.0_dbl)) ) return
       if( ip /= k ) then
          t = w(k,:)
          w(k,:) = w(ip,:)
          w(ip,:) = t
       end if
       do i = k+1, n
          w(i,k:) = w(i,k:) - w(i,k)/w(k,k)*w(k,k:)
       end do
    end do
    do i = n, 1, -1
       x(i) = (w(i,n+1) - sum(w(i,i+1:n)*x(i+1:n))) / w(i,i)
    end do
    status = 0

  end subroutine gsolve

end module colterm
