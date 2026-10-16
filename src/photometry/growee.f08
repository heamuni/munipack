!
!  Growth-curve photometry by the enclosed energy of the apertures
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
!  The growth curve of the frame is the enclosed energy EE_k of the
!  apertures k, fitted on the bright isolated stars by aphot's
!  automatic aperture correction (apee_auto: the plateau of the growth
!  curves by its slope, EE varying over the field, the errors by the
!  covariance of the nested apertures). Every star has the counts C_k
!  of the apertures up to the plateau, the growth curve gives
!
!     C_k = F EE_k + (pixel noise) + A_k * (the error of the mean sky)
!                  + F * (the error of EE_k),
!
!  F the total flux (GCOUNT, the flux of the plateau apertures). The
!  covariance of the counts is that of nested apertures (they share
!  their pixels and the sky),
!
!     Cov(C_i,C_j) = Vp(min(i,j)) + A_i A_j SKYERR**2,
!
!  Vp the noise of the pixels of the aperture: the sky and read noise
!  of APCOUNTERR (APCOUNTERR**2 - A**2 SKYERR**2 - C/g) and the photon
!  noise of the star by the model F EE/g (g = APERRGN of APERPHOT; the
!  measured counts when unknown). The estimates of F are the weighted
!  sums of the apertures,
!
!     F = sum w_k C_k,   sum w_k EE_k = 1,
!
!  the best linear unbiased estimate of the apertures 1 .. m by the
!  covariance of the counts (it weights the pixels by the profile, as
!  PSF photometry: faint stars by the inner apertures), for every m up
!  to the plateau, and every single aperture (C_k/EE_k). The variance
!  of an estimate is w' Cov w + F**2 (sum |w_k| dEE_k)**2: the errors of
!  EE of the apertures are not independent (one growth curve, one
!  normalisation), the bound of fully correlated errors is used. The
!  estimate of the least variance is GCOUNT, its standard deviation
!  GCOUNTERR: faint stars get the profile-weighted sum (the noise of the
!  sky dominates), bright stars the larger apertures, where EE is known
!  better than the photon noise of the star.
!
!  The residual sky b (the error of the mean sky of the ring, common
!  to all apertures) is estimated together with F by all the apertures
!  to the plateau, with the prior b = 0 +- SKYERR: SKY + b is the sky of
!  the star by the ring and the pixels of the apertures, its error the
!  posterior error of b (smaller than SKYERR).
!

module grow_ee

  use iso_fortran_env

  implicit none
  private

  integer, parameter, private :: dbl = REAL64

  public :: grow_ee_phot, grow_ee_curve

contains

  subroutine grow_ee_phot(area,kend,gain,apcts,apcts_err,sky_err,ee,eeerr, &
       gcount,gcount_err,skycorr,skycorr_err,naper)

    ! The total fluxes gcount (and their errors) of all stars by the
    ! enclosed energy ee, eeerr of the apertures 1 .. kend; area the
    ! areas of the apertures, gain the gain of the photon noise of
    ! apcts_err (<= 0 unknown), sky_err the error of the mean sky;
    ! skycorr is the residual sky b [ADU/pix] and skycorr_err the error
    ! of the sky + b (-1 when unknown); naper the largest aperture of
    ! the estimate (0 for an invalid star)

    real(dbl), dimension(:), intent(in) :: area
    integer, intent(in) :: kend
    real(dbl), intent(in) :: gain
    real, dimension(:,:), intent(in) :: apcts, apcts_err
    real, dimension(:), intent(in) :: sky_err
    real(dbl), dimension(:,:), intent(in) :: ee, eeerr
    real(dbl), dimension(:), intent(out) :: gcount, gcount_err
    real, dimension(:), intent(out) :: skycorr, skycorr_err
    integer, dimension(:), intent(out) :: naper

    integer :: l

    do l = 1, size(apcts,1)
       call star_flux(l)
    end do

  contains

    subroutine star_flux(l)

      integer, intent(in) :: l

      integer, dimension(:), allocatable :: kk
      real(dbl), dimension(:), allocatable :: c, v, a, e, de, vb, vp, w, wbest
      real(dbl), dimension(:,:), allocatable :: cov
      real(dbl) :: s2, f, fnew, var, vbest
      integer :: n, m, k, it
      logical :: ok

      gcount(l) = -1
      gcount_err(l) = -1
      skycorr(l) = 0
      skycorr_err(l) = -1
      naper(l) = 0

      ! the valid apertures up to the plateau with EE
      kk = pack([(k, k = 1, min(kend,size(apcts,2)))], &
           [(apcts(l,k) > 0 .and. apcts_err(l,k) > 0 .and. &
           ee(l,k) > 0 .and. eeerr(l,k) >= 0, k = 1, min(kend,size(apcts,2)))])
      n = size(kk)
      if( n == 0 ) return

      allocate(c(n),v(n),a(n),e(n),de(n),vb(n),vp(n),w(n),wbest(n),cov(n,n))
      c = apcts(l,kk)
      v = real(apcts_err(l,kk),dbl)**2
      a = area(kk)
      e = ee(l,kk)
      de = eeerr(l,kk)
      s2 = max(real(sky_err(l),dbl),0.0_dbl)**2

      ! the noise of the pixels without the photon noise of the star
      ! (the sky, read noise), or with the measured one (gain unknown)
      if( gain > 0 ) then
         vb = max(v - a**2*s2 - c/gain,0.0_dbl)
      else
         vb = max(v - a**2*s2,0.0_dbl)
      end if

      ! the initial flux: the single aperture of the least variance
      ! by the measured counts
      k = minloc((v + (c*de/e)**2)/e**2,1)
      f = c(k) / e(k)

      do it = 1, 5

         ! the covariance of the nested apertures by the model flux
         if( gain > 0 ) then
            vp = vb + max(f*e,0.0_dbl)/gain
         else
            vp = vb
         end if
         do k = 2, n
            vp(k) = max(vp(k),vp(k-1))
         end do
         do m = 1, n
            do k = 1, n
               cov(k,m) = vp(min(k,m)) + a(k)*a(m)*s2
            end do
         end do

         ! the estimate of the least variance: single apertures and
         ! the best linear estimates of the apertures 1 .. m
         vbest = huge(vbest)
         wbest = 0
         do k = 1, n
            var = cov(k,k)/e(k)**2 + (f*de(k)/e(k))**2
            if( var < vbest ) then
               vbest = var
               wbest = 0
               wbest(k) = 1 / e(k)
            end if
         end do
         do m = 2, n
            call blue(cov(1:m,1:m),e(1:m),w(1:m),var,ok)
            if( .not. ok ) cycle
            var = var + f**2*sum(abs(w(1:m))*de(1:m))**2
            if( var < vbest ) then
               vbest = var
               wbest = 0
               wbest(1:m) = w(1:m)
            end if
         end do

         fnew = sum(wbest*c)
         if( abs(fnew - f) <= 1e-7_dbl*abs(f) ) then
            f = fnew
            exit
         end if
         f = fnew
      end do

      if( .not. (f > 0 .and. vbest > 0) ) return
      gcount(l) = f
      gcount_err(l) = sqrt(vbest)
      naper(l) = kk(findloc(wbest /= 0,.true.,dim=1,back=.true.))

      ! the residual sky b with the prior 0 +- SKYERR, together with F,
      ! by all the apertures (their pixel noise)
      if( s2 > 0 ) then
         block
           real(dbl), dimension(2,2) :: p, pinv
           real(dbl), dimension(2) :: q, theta
           real(dbl), dimension(n) :: ue, ua
           do m = 1, n
              do k = 1, n
                 cov(k,m) = vp(min(k,m))
              end do
           end do
           call cholsolve(cov,e,ue,ok)
           if( ok ) call cholsolve(cov,a,ua,ok)
           if( ok ) then
              p(1,1) = dot_product(e,ue)
              p(1,2) = dot_product(e,ua)
              p(2,1) = p(1,2)
              p(2,2) = dot_product(a,ua) + 1/s2
              q(1) = dot_product(ue,c)
              q(2) = dot_product(ua,c)
              call inv2(p,pinv,ok)
              if( ok ) then
                 theta = matmul(pinv,q)
                 if( pinv(2,2) > 0 ) then
                    skycorr(l) = real(theta(2))
                    skycorr_err(l) = real(sqrt(pinv(2,2)))
                 end if
              end if
           end if
         end block
      end if

    end subroutine star_flux

  end subroutine grow_ee_phot


  subroutine blue(cov,e,w,var,ok)

    ! The best linear unbiased estimate of F by C = F e + noise of the
    ! covariance cov: the weights w (F = w'C, w'e = 1) and the variance

    real(dbl), dimension(:,:), intent(in) :: cov
    real(dbl), dimension(:), intent(in) :: e
    real(dbl), dimension(:), intent(out) :: w
    real(dbl), intent(out) :: var
    logical, intent(out) :: ok

    real(dbl) :: s

    call cholsolve(cov,e,w,ok)
    if( .not. ok ) return
    s = dot_product(e,w)
    ok = s > 0
    if( .not. ok ) return
    w = w / s
    var = 1 / s

  end subroutine blue


  subroutine cholsolve(a,b,x,ok)

    ! the solution of a x = b, a symmetric positive definite (Cholesky)

    real(dbl), dimension(:,:), intent(in) :: a
    real(dbl), dimension(:), intent(in) :: b
    real(dbl), dimension(:), intent(out) :: x
    logical, intent(out) :: ok

    real(dbl), dimension(size(b),size(b)) :: g
    real(dbl) :: s
    integer :: n, i, j

    n = size(b)
    g = 0
    ok = .false.
    do j = 1, n
       s = a(j,j) - sum(g(j,1:j-1)**2)
       if( .not. (s > 1e-12_dbl*abs(a(j,j))) .or. .not. (s > 0) ) return
       g(j,j) = sqrt(s)
       do i = j + 1, n
          g(i,j) = (a(i,j) - sum(g(i,1:j-1)*g(j,1:j-1))) / g(j,j)
       end do
    end do
    do i = 1, n
       x(i) = (b(i) - sum(g(i,1:i-1)*x(1:i-1))) / g(i,i)
    end do
    do i = n, 1, -1
       x(i) = (x(i) - sum(g(i+1:n,i)*x(i+1:n))) / g(i,i)
    end do
    ok = .true.

  end subroutine cholsolve


  subroutine inv2(a,b,ok)

    real(dbl), dimension(2,2), intent(in) :: a
    real(dbl), dimension(2,2), intent(out) :: b
    logical, intent(out) :: ok
    real(dbl) :: d

    d = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    ok = abs(d) > 0
    if( .not. ok ) return
    b(1,1) = a(2,2) / d
    b(2,2) = a(1,1) / d
    b(1,2) = -a(1,2) / d
    b(2,1) = -a(2,1) / d

  end subroutine inv2


  subroutine grow_ee_curve(raper,area,ee,eeerr,curve,curve_err,prof)

    ! The growth curve of the frame (GROWCURVE): the median EE of the
    ! stars in every aperture and its error (-1 for the apertures
    ! without EE), the profile: the mean surface brightness in the
    ! central disc and in the annuli between the apertures [1/pix]

    real, dimension(:), intent(in) :: raper
    real(dbl), dimension(:), intent(in) :: area
    real(dbl), dimension(:,:), intent(in) :: ee, eeerr
    real(dbl), dimension(:), intent(out) :: curve, curve_err, prof

    real(dbl), dimension(:), allocatable :: x
    real(dbl) :: c0, a0
    integer :: k

    do k = 1, size(raper)
       x = pack(ee(:,k),ee(:,k) > 0 .and. eeerr(:,k) >= 0)
       if( size(x) > 0 ) then
          curve(k) = median(x)
          x = pack(eeerr(:,k),ee(:,k) > 0 .and. eeerr(:,k) >= 0)
          curve_err(k) = median(x)
       else
          curve(k) = -1
          curve_err(k) = -1
       end if
    end do

    c0 = 0
    a0 = 0
    do k = 1, size(raper)
       if( curve(k) >= 0 ) then
          prof(k) = max((curve(k) - c0) / (area(k) - a0),0.0_dbl)
          c0 = curve(k)
          a0 = area(k)
       else
          prof(k) = 0
       end if
    end do

  end subroutine grow_ee_curve


  real(dbl) function median(a)
    real(dbl), dimension(:), intent(in) :: a
    real(dbl), dimension(:), allocatable :: b
    real(dbl) :: t
    integer :: n, i, j
    n = size(a)
    b = a
    ! insertion sort (short arrays)
    do i = 2, n
       t = b(i)
       j = i - 1
       do while( j >= 1 )
          if( b(j) <= t ) exit
          b(j+1) = b(j)
          j = j - 1
       end do
       b(j+1) = t
    end do
    if( mod(n,2) == 1 ) then
       median = b(n/2+1)
    else
       median = (b(n/2) + b(n/2+1)) / 2
    end if
  end function median

end module grow_ee
