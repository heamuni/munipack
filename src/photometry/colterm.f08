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
!  Several colour terms (COLT_CAT, at most ncolmax of exscatter.f08):
!
!      m_cat - m_inst = ZP + sum_j k1_j (c_j - cref_j) + k2_j (c_j - cref_j)**2,
!
!  a passband that depends on more than one colour (QUVIK 260-360 nm
!  against SDSS u: u - g and g - r); beta(2j-1), beta(2j) are b1, b2 of
!  the colour j, the covariance of beta in the same order. The errors of
!  the colours of a star are correlated when they share a band (u - g and
!  g - r share g: cov = -err(g)**2), the effective variance includes it.
!  A star without some of its colours gets them predicted by the colours
!  it has, by the colour-colour relation of the calibration stars (their
!  colours as a Gaussian, the conditional mean, colour_predict), with
!  the conditional scatter as their error (without any colour, the mean
!  colours and their scatter); the calibration uses the stars with all
!  the colours only.
!
!  The catalogue colours scatter about the true ones by their errors,
!  which biases the coefficients towards zero (the regression dilution:
!  by the reliability ratio lambda, the variance of the true colours over
!  the observed one) and, as the colour contains the catalogue magnitude
!  of the calibration, by the covariance of their errors. With the error
!  columns, the fit of one or several colours is corrected for them
!  (errcorr, the moments of Fuller 1987, mixfit of exscatter.f08), capped
!  when a combination of the terms is hardly wider than its errors.
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
     ! the number of the colours (terms), one for COLT_TRA
     integer :: nterm = 1
     ! the colours (j) of every band (n): the names of the catalogue
     ! columns (COLT_CAT) or of the filters (COLT_TRA), c = A - B, (j,n)
     character(len=68), dimension(:,:), allocatable :: cola, colb
     ! COLT_TRA: the indices of the filters of the colour of every band
     integer, dimension(:), allocatable :: ia, ib
     ! the order of the polynomial of every colour (j,n), the number of
     ! the stars used (with all the colours)
     integer, dimension(:,:), allocatable :: order
     integer, dimension(:), allocatable :: nstars
     ! ln(ph/ct) = beta(0) + sum_j beta(2j-1) d_j + beta(2j) d_j**2,
     ! d_j = c_j - cref_j (beta(0,:) is ln CTPH), the covariance of beta
     ! (in the order 0,1,2,...)
     real(dbl), dimension(:,:), allocatable :: beta
     real(dbl), dimension(:,:,:), allocatable :: cov
     ! the reference colours (j,n); the mean colours of the calibration
     ! stars, their scatter, the mean colour corrections g_j of
     ! the calibration stars and their scatter (used for the stars
     ! without a colour), the covariances of the colours of the calibration
     ! stars (j,k,n; several colours: the prediction of a missing colour)
     real(dbl), dimension(:,:), allocatable :: cref, cmean, cstd, gmean, gstd
     real(dbl), dimension(:,:,:), allocatable :: ccs
     ! COLT_CAT: the catalogue (all its stars): positions [deg], colours
     ! (i,j), their errors (< 0 for an undefined colour), the covariances
     ! of the colours of a star (i,j,k), the search radius [deg]
     real(dbl), dimension(:), allocatable :: ra, dec
     real(dbl), dimension(:,:), allocatable :: col, dcol
     real(dbl), dimension(:,:,:), allocatable :: ccov
     real(dbl) :: tol = 0
     ! COLT_CAT: the correction of the fit for the errors of the catalogue
     ! colours (the regression dilution, --color-err-correct), the colours
     ! with the error columns (j); the reliability ratios of the colours of
     ! the calibration stars, the variance of the true colours over
     ! the observed one (j,n; 1 without the errors), the status of
     ! the correction (n): 0 not applied, 1 applied, 2 applied with the cap
     ! (the corrected moments of a combination of the terms below lammin
     ! of exscatter.f08), 3 applied, the fit not converged; the fraction of
     ! the correction of the dilution applied (n; 1 full, < 1 capped)
     logical :: errcorr = .true.
     logical, dimension(:), allocatable :: colerr
     real(dbl), dimension(:,:), allocatable :: lam
     integer, dimension(:), allocatable :: eivst
     real(dbl), dimension(:), allocatable :: eivh
  end type colour_term

  public :: colterm_init, colour_corr, colour_predict, colour_solve

contains

  subroutine colterm_init(ctm,model,nbands,nterm)

    ! nterm (optional, 1 by default) the number of the colours

    type(colour_term), intent(out) :: ctm
    integer, intent(in) :: model, nbands
    integer, intent(in), optional :: nterm

    integer :: nt

    nt = 1
    if( present(nterm) ) nt = nterm
    ctm%model = model
    ctm%nterm = nt
    allocate(ctm%cola(nt,nbands),ctm%colb(nt,nbands),ctm%ia(nbands), &
         ctm%ib(nbands),ctm%order(nt,nbands),ctm%nstars(nbands), &
         ctm%beta(0:2*nt,nbands),ctm%cov(2*nt+1,2*nt+1,nbands), &
         ctm%cref(nt,nbands),ctm%cmean(nt,nbands),ctm%cstd(nt,nbands), &
         ctm%gmean(nt,nbands),ctm%gstd(nt,nbands),ctm%ccs(nt,nt,nbands), &
         ctm%colerr(nt),ctm%lam(nt,nbands),ctm%eivst(nbands),ctm%eivh(nbands))
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
    ctm%ccs = 0
    ctm%colerr = .false.
    ctm%lam = 1
    ctm%eivst = 0
    ctm%eivh = 0
    allocate(ctm%ra(0),ctm%dec(0),ctm%col(0,nt),ctm%dcol(0,nt),ctm%ccov(0,nt,nt))

  end subroutine colterm_init


  subroutine colour_corr(ctm,n,c,dc,g,dg,ccv)

    ! The colour correction g = ln(ph/ct) - ln CTPH of a star of
    ! the colours c (with the errors dc, dc < 0: an undefined colour) in
    ! the band n, and its error dg: the error of the calibration at
    ! the colours (by the covariance of beta, the error of CTPH included)
    ! and the errors of the colours, ccv (optional) their covariances
    ! (j /= k). A star without a colour gets the mean correction of
    ! the calibration stars for it, its error inflated by their scatter
    ! (one colour); with several colours, the missing ones are predicted
    ! by the known ones (colour_predict).

    type(colour_term), intent(in) :: ctm
    integer, intent(in) :: n
    real(dbl), dimension(:), intent(in) :: c, dc
    real(dbl), intent(out) :: g, dg
    real(dbl), dimension(:,:), intent(in), optional :: ccv

    real(dbl), dimension(3) :: x
    real(dbl), dimension(2*ctm%nterm+1) :: xm
    real(dbl), dimension(ctm%nterm) :: gpm, gk, cp
    real(dbl), dimension(ctm%nterm,ctm%nterm) :: cv, a
    logical, dimension(ctm%nterm) :: known
    real(dbl) :: d, gp, v
    integer :: j, k, nt

    nt = ctm%nterm
    if( nt == 1 ) then
       if( dc(1) >= 0 ) then
          d = c(1) - ctm%cref(1,n)
          g = ctm%beta(1,n)*d + ctm%beta(2,n)*d**2
          gp = ctm%beta(1,n) + 2*ctm%beta(2,n)*d
          x = [ 1.0_dbl, d, d**2 ]
          v = dot_product(x,matmul(ctm%cov(:,:,n),x)) + (gp*dc(1))**2
       else
          g = ctm%gmean(1,n)
          v = ctm%cov(1,1,n) + ctm%gstd(1,n)**2
       end if
       dg = sqrt(max(v,0.0_dbl))
       return
    end if

    ! several colours: the missing ones predicted by the known ones
    ! (colour_predict), the gradient of the correction by the known
    ! colours includes their part through the prediction
    call colour_predict(ctm,n,c,dc,cp,cv,a)
    known = dc(1:nt) >= 0
    g = 0
    xm = 0
    xm(1) = 1
    gpm = 0
    do j = 1, nt
       d = cp(j) - ctm%cref(j,n)
       g = g + ctm%beta(2*j-1,n)*d + ctm%beta(2*j,n)*d**2
       gpm(j) = ctm%beta(2*j-1,n) + 2*ctm%beta(2*j,n)*d
       xm(2*j) = d
       xm(2*j+1) = d**2
    end do
    v = dot_product(xm,matmul(ctm%cov(:,:,n),xm))
    ! the missing colours: the conditional covariance
    do j = 1, nt
       if( known(j) ) cycle
       do k = 1, nt
          if( .not. known(k) ) v = v + gpm(j)*gpm(k)*cv(j,k)
       end do
    end do
    ! the known colours: their errors, by the total gradient
    gk = 0
    do k = 1, nt
       if( .not. known(k) ) cycle
       gk(k) = gpm(k)
       do j = 1, nt
          if( .not. known(j) ) gk(k) = gk(k) + gpm(j)*a(j,k)
       end do
    end do
    do j = 1, nt
       if( .not. known(j) ) cycle
       v = v + (gk(j)*dc(j))**2
       if( present(ccv) ) then
          do k = j + 1, nt
             if( known(k) ) v = v + 2*gk(j)*gk(k)*ccv(j,k)
          end do
       end if
    end do
    dg = sqrt(max(v,0.0_dbl))

  end subroutine colour_corr


  subroutine colour_predict(ctm,n,c,dc,cp,cv,a)

    ! The colours of a star with the missing ones (dc < 0) predicted by
    ! the known ones: the colours of the calibration stars of the band n
    ! as a Gaussian (the mean cmean, the covariance ccs), the conditional
    ! mean of the missing colours M given the known K,
    !
    !    c_M = mu_M + A (c_K - mu_K),  A = S_MK S_KK**-1,
    !
    ! and their conditional covariance cv = S_MM - A S_KM (the scatter
    ! about the colour-colour relation of the calibration stars, small for
    ! u - g and g - r of the main sequence); without a known colour,
    ! the mean and the covariance. cp the colours (the known ones as
    ! given), cv (missing rows and columns only), a(j,k) the coefficients
    ! A of the missing j by the known k (zero elsewhere).

    type(colour_term), intent(in) :: ctm
    integer, intent(in) :: n
    real(dbl), dimension(:), intent(in) :: c, dc
    real(dbl), dimension(:), intent(out) :: cp
    real(dbl), dimension(:,:), intent(out) :: cv, a

    real(dbl), dimension(ctm%nterm,ctm%nterm) :: s, sk
    logical, dimension(ctm%nterm) :: known
    integer, dimension(ctm%nterm) :: ik
    real(dbl) :: det
    integer :: nt, nk, j, k, i, l

    nt = ctm%nterm
    s = ctm%ccs(:,:,n)
    known = dc(1:nt) >= 0
    cp = c(1:nt)
    cv = 0
    a = 0
    nk = 0
    do j = 1, nt
       if( known(j) ) then
          nk = nk + 1
          ik(nk) = j
       end if
    end do
    if( nk == nt ) return

    ! the inverse of S_KK (1 x 1 or 2 x 2)
    sk = 0
    if( nk == 1 ) then
       if( s(ik(1),ik(1)) > 0 ) sk(1,1) = 1 / s(ik(1),ik(1))
    else if( nk == 2 ) then
       det = s(ik(1),ik(1))*s(ik(2),ik(2)) - s(ik(1),ik(2))**2
       if( det > 0 ) then
          sk(1,1) = s(ik(2),ik(2)) / det
          sk(2,2) = s(ik(1),ik(1)) / det
          sk(1,2) = -s(ik(1),ik(2)) / det
          sk(2,1) = sk(1,2)
       end if
    end if

    do j = 1, nt
       if( known(j) ) cycle
       do k = 1, nk
          do l = 1, nk
             a(j,ik(k)) = a(j,ik(k)) + s(j,ik(l))*sk(l,k)
          end do
       end do
       cp(j) = ctm%cmean(j,n)
       do k = 1, nk
          cp(j) = cp(j) + a(j,ik(k))*(c(ik(k)) - ctm%cmean(ik(k),n))
       end do
    end do
    do j = 1, nt
       if( known(j) ) cycle
       do i = 1, nt
          if( known(i) ) cycle
          cv(j,i) = s(j,i)
          do k = 1, nk
             cv(j,i) = cv(j,i) - a(j,ik(k))*s(ik(k),i)
          end do
       end do
    end do

  end subroutine colour_predict


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
       d = m(ctm%ia) - m(ctm%ib) - ctm%cref(1,:)
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

    d = m(ctm%ia) - m(ctm%ib) - ctm%cref(1,:)
    gp = k1 + 2*k2*d
    call jacobian
    c = d + ctm%cref(1,:)

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
