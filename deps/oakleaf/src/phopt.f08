!
!  Robust estimation of mean photon rates
!
!  Copyright © 2020 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of OakLeaf.
!
!  OakLeaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  OakLeaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!


module photonrate_opt

  use photonrate_common
  use iso_fortran_env

  implicit none
  private

  public :: phopt

  ! debugging
  logical, parameter :: debug_graph = .false.
  logical, parameter :: debug_graphpsi = .false.
  logical, parameter :: debug_test = .false.

  ! Does to use Jacobian analytic, or numerical?
  logical, parameter :: analytic = .true.

  ! I had changed meaning of nprint parameter of Minpack;
  ! the paramater calls a minimisation routine every nprint
  ! steps, by design.
  !   This routine uses the nprint calls for checking of variables
  ! on liyng in a desired region (tau is bracketed by quantiles,
  ! rho by its initial estimate). The checking requires nprint > 0.
  ! Regular debug prints should to set: debug_iter = .true.
  integer, parameter :: nprint = 10

  real(REAL64), dimension(:), pointer :: mju,tt, kmean,tmin,tmax
  integer, dimension(:), pointer :: iset
  real(REAL64) :: phifix, sx,sy,taumin, taumax, rmax

  character(len=*), parameter :: error_label = 'Error in phopt: '
  character(len=*), parameter :: warning_label = 'Warning in phopt: '


contains

  subroutine phopt(mjus,tts,kmeans,isets,tmins,tmaxs,sxs,sys, &
       tau,dtau,rho,drho,phi,dphi,reli)

    real(REAL64), dimension(:), target, intent(in) :: mjus,tts,kmeans,tmins,tmaxs
    integer , dimension(:), target, intent(in) :: isets
    real(REAL64), intent(in) :: sxs, sys
    real(REAL64), intent(in out) :: tau, dtau, rho, phi
    real(REAL64), intent(out) :: drho, dphi
    logical, intent(out) :: reli

    real(REAL64) :: tau0, rho0, phi0
    integer :: n

    reli = .false.

    mju => mjus
    tt => tts
    kmean => kmeans
    iset => isets
    tmin => tmins
    tmax => tmaxs
    sx = sxs
    sy = sys

    n = size(tt)
    tau0 = tau
    rho0 = rho
    phi0 = phi
    drho = -1
!    dphi = -1

    if( debug_test ) then
       call phoptest(tau,rho,phi)
       stop 'phoptest'
    end if

    ! hard limits of parameters
    taumin = max(0.9*minval(tmin),real(0.1,REAL64))
    taumax = 1.1*maxval(tmax)
    rmax = 10*(abs(rho) + 1)
    rmax = maxval(mju)
    if( debug ) write(error_unit,*) 'taumin: ',real(taumin), &
         ' taumax: ',real(taumax),'  rmax:',real(rmax)

    ! initial estimate of phi, grid search
!!$    block
!!$      real(REAL64) :: fi, f, fmin
!!$      integer :: i
!!$
!!$      fmin = huge(f)
!!$      do i = 15, 75, 20
!!$         fi = i / RAD
!!$         f = elogf(tau,rho,fi)
!!$         if( f < fmin ) then
!!$            fmin = f
!!$            phi = fi
!!$         end if
!!$         if( debug_iter ) write(error_unit,*) 'phi[deg] on grid:',i,real(f)
!!$      end do
!!$      phi0 = phi
!!$    end block


    if( n > 3 .and. .true. ) then

       call thelike(tau,dtau,rho,drho,phi,dphi,reli)
       if( debug ) write(error_unit, &
            '(a,(g0,a,es12.3,2x),(g0,a,es8.1,2x),(f0.1,a,es8.1,2x),l2)') &
            "The Like: ", &
            real(tau),'+-',dtau,real(rho),'+-',drho,RAD*phi,'+-',RAD*dphi,reli

       if( .not. reli ) then
          tau = tau0
          rho = rho0
          phi = phi0
       end if

    end if

    if( n > 2 .and. .not. reli .and. .true. ) then

       call onplane(tau,dtau,rho,drho,phi,reli)
       if( debug ) write(error_unit, &
            '(a,g0,a,es12.3,2x,1(f7.1,a,es8.1,2x),l2)') "OnPlane: ", &
            real(tau),'+-',dtau,rho,'+-',drho,reli

       if( .not. reli ) then
          tau = tau0
          rho = rho0
       end if

    end if

    if( .not. reli .and. .true. ) then

       call online(tau,dtau,rho,phi,reli)

       if( debug ) write(error_unit, &
            '(a,g0,a,es12.3,l2)') "Online: ",real(tau),'+-',dtau,reli

       if( .not. reli ) tau = tau0
    end if

    if( .not. reli ) then

       ! fail-back, last resort
       call bisection(tau,dtau,rho,phi,reli)

       if( debug ) write(error_unit,'(a,g0,l2)') "Bisection: ",real(tau),reli

       if( .not. reli ) tau = tau0
    end if


!!$    block
!!$      real(REAL64), dimension(3,size(kmeans)) :: psimin, psimax
!!$      real(REAL64), dimension(3,3,size(kmeans)) :: dpsimin, dpsimax
!!$      logical, dimension(3), parameter :: mask = [ .true., .false., .false. ]
!!$      integer :: i, nlo, nhi
!!$      real(REAL64) :: flo, fhi, df, delta
!!$      write(*,*) n, count(tt < tmin(iset)), count(tt > tmax(iset))
!!$      do i = 1,size(kmeans)
!!$         psimin(:,i) = psilike(tmin(i),kmean(i),tau,rho,cos(phi),sin(phi),mask)
!!$         psimax(:,i) = psilike(tmax(i),kmean(i),tau,rho,cos(phi),sin(phi),mask)
!!$         dpsimin(:,:,i) = dpsilike(tmin(i),kmean(i),tau,rho,cos(phi),sin(phi),mask)
!!$         dpsimax(:,:,i) = dpsilike(tmax(i),kmean(i),tau,rho,cos(phi),sin(phi),mask)
!!$         write(*,*) psimin(1,i),psimax(1,i),dpsimin(1,1,i),dpsimax(1,1,i)
!!$      end do
!!$      nlo = count(tt < tmin(iset))
!!$      nhi = count(tt > tmax(iset))
!!$      flo = psimin(1,1)
!!$      fhi = psimax(1,1)
!!$      df = (dpsimin(1,1,1) + dpsimax(1,1,1)) / 2
!!$      delta = (nlo*flo - nhi*fhi)/(df*n)
!!$      write(*,*) n, -38*count(tt < tmin(iset)), 30*count(tt > tmax(iset))
!!$      write(*,*) psilike(tau,kmean(1),tau,rho,cos(phi),sin(phi),mask), delta
!!$!      tau = tau - delta
!!$    end block

    if( debug_graphpsi ) then
       call graph_psi_ttau('/tmp/psi_ttau',mju(1),tau,rho,phi)
       call graph_psi_trho('/tmp/psi_trho',mju(1),tau,rho,phi)
       call graph_psi_tphi('/tmp/psi_tphi',mju(1),tau,rho,phi)
          call graph_like_tr('/tmp/like_tr',tau,dtau,rho,phi)
          call graph_like_tf('/tmp/like_tf',tau,dtau,rho)
    end if

  end subroutine phopt


  subroutine thelike(tau,dtau,rho,drho,phi,dphi,reli)

    use minpacks

    real(REAL64), intent(in out) :: tau, dtau, rho, phi
    real(REAL64), intent(out) :: drho,dphi
    logical, intent(out) :: reli

    integer, parameter :: maxfev = 666
!    real(REAL64), parameter :: factor = 100
    integer, parameter :: npar = 3
    real(REAL64), dimension(npar) :: x, diag
    real(REAL64), dimension(npar,npar) :: cov
    character(len=80) :: msg
    character(len=5) :: label
    real(REAL64) :: eps
    integer :: info, n, i

    reli = .false.
    n = size(tt)

    x = [ tau, rho, phi ]
    eps = max(10*epsilon(x)*(1 + maxval(abs(x))),1e-3*dtau)
    if( analytic ) then
       label = 'lmder'
       call lmder3(maxlike,x,cov,eps,maxfev,nprint=nprint,info=info)
    else
       label = 'lmdif'
       call lmdif3(maxliked,x,cov,eps,maxfev,nprint=nprint,info=info)
    end if
    diag = [ (cov(i,i),i=1,npar) ]

    reli = taumin < x(1) .and. x(1) < taumax .and. abs(x(2)) < rmax &
         .and. all(diag(1:2) > 0) .and. (0 < info .and. info < 5)

    if( debug ) then
       msg = ''
       if( 0 < info .and. info < 5 ) then
          msg = ' (success)'
       else if( info == -1 ) then
          msg = ' (out of range)'
       else if( info == 5 ) then
          msg = ' (no convergence)'
       else if( .not. all(diag > 0) ) then
          msg = ' (no minimum)'
       end if
       write(error_unit,'(3a,i0,a,a,l2)') ' thelike() ',label,' info: ', &
            info,trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       tau = x(1)
       rho = abs(x(2)) ! |r|

       ! 0 .. PI/2
       if( x(3) < 0 .or. x(3) > PI/2 ) then
          phi = mod(abs(x(3)),PI/2)
       else
          phi = x(3)
       end if

       dtau = sqrt(cov(1,1)/(n-3))
       drho = sqrt(cov(2,2)/(n-3))

       if( cov(3,3) > 0 ) then
          dphi = sqrt(cov(3,3)/(n-3))
       else
          dphi = 90 / RAD
       end if
    else
       if( debug_graph ) then
          call graph_like_tr('/tmp/like_tr',tau,dtau,rho,phi)
          call graph_like_tf('/tmp/like_tf',tau,dtau,rho)
       end if
    end if
!call graph_maxlike(tau,rho,phi)

    if( debug ) then
       block
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64) :: q
         integer :: i,j

         cor = 1
         do i = 1, npar
            do j = i+1, npar
               if( cov(i,i) > 0 .and. cov(j,j) > 0 &
                    .and. cov(i,i) < 1e10 .and. cov(j,j) < 1e10 ) then
                  q = sqrt(cov(i,i)*cov(j,j))
                  if( 1e-10 < q .and. q < 1e10 .and. abs(cov(i,j)) < 1e10 ) then
                     cor(i,j) = cov(i,j) / q
                  else
                     cor(i,j) = huge(cor)
                  end if
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         write(error_unit,'(a,3(es12.4,1x),t48,a,3(f7.4,1x))') &
              'Cov: ', cov(1,:),'corr: ', cor(1,:)
         do i = 2,npar
            write(error_unit,'(5x,3(es12.4,1x),t48,6x,3(f7.4,1x))') &
                 cov(i,:),cor(i,:)
         end do

       end block
    end if

  end subroutine thelike

  subroutine maxlike(k,l,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: k, l, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(k), intent(in) :: x
    real(REAL64), dimension(l), intent(in out) :: fvec
    real(REAL64), dimension(ldfjac,k), intent(in out) :: fjac

    logical, dimension(3), parameter :: mask = .true.

    real(REAL64) :: tau, rho, fi

    tau = x(1)
    rho = x(2)
    fi = x(3)

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
          write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
          write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
          write(error_unit,*) 'fjac(3,:) = ',real(fjac(3,:))
       end if

       if( tau < taumin .or. tau > taumax .or. abs(rho) > rmax ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

       fvec = psisum(tt,mju,tau,rho,fi,tmin,tmax,kmean,mask)

    else if( iflag == 2 ) then

       fjac = dpsisum(tt,mju,tau,rho,fi,tmin,tmax,kmean,mask)

    end if

  end subroutine maxlike


  subroutine maxliked(k,l,x,fvec,iflag)

    integer, intent(in) :: k, l
    integer, intent(in out) :: iflag
    real(REAL64), dimension(k), intent(in) :: x
    real(REAL64), dimension(l), intent(in out) :: fvec
    real(REAL64), dimension(3,3) :: fjac
    integer, parameter :: ldfjac = 3
    integer :: nflag

    if( iflag == 0 ) then
       nflag = 0
    else
       nflag = 1
    end if
    call maxlike(k,l,x,fvec,fjac,ldfjac,nflag)

  end subroutine maxliked



  subroutine onplane(tau,dtau,rho,drho,phi,reli)

    use minpacks

    real(REAL64), intent(in out) :: tau, dtau, rho
    real(REAL64), intent(in) :: phi
    real(REAL64), intent(out) :: drho
    logical, intent(out) :: reli

    integer, parameter :: maxfev = 666
!    real(REAL64), parameter :: factor = 100
    integer, parameter :: npar = 2
    real(REAL64), dimension(npar) :: x
    real(REAL64), dimension(npar,npar) :: cov
    character(len=80) :: msg
    real(REAL64) :: eps
    integer :: info, n

    reli = .false.
    n = size(tt)

    x = [ tau, rho ]
    phifix = phi
    eps = max(10*epsilon(x)*(1 + maxval(abs(x))),1e-3*dtau)
    call lmder3(planelike,x,cov,eps,maxfev,nprint=nprint,info=info)

    reli = taumin < x(1) .and. x(1) < taumax .and. abs(x(2)) < rmax &
         .and. cov(1,1) > 0 .and. cov(2,2) > 0 .and. (0 < info .and. info < 5)

    if( debug ) then
       msg = ''
       if( 0 < info .and. info < 5 ) then
          msg = ' (success)'
       else if( info == -1 ) then
          msg = ' (out of range)'
       else if( info == 5 ) then
          msg = ' (no convergence)'
       else if( .not. (cov(1,1) > 0 .and. cov(2,2) > 0) ) then
          msg = ' (no minimum)'
       end if
       write(error_unit,'(a,i0,a,a,l2)') ' onplane() lmder info: ', &
            info,trim(msg),',  reliable:',reli
    end if

    if( reli ) then
       tau = x(1)
       rho = abs(x(2))
       ! the function is symetric in rho, we accept both +-rho,
       ! the solutions can be tunneled to mirrored half-space
       ! due to adjust of the ML lambda factor.

       dtau = sqrt(cov(1,1)/(n-2))
       drho = sqrt(cov(2,2)/(n-2))
    else
       if( debug_graph ) call graph_like_tr('plane_tr',tau,dtau,rho,phi)
    end if
!call graph_like_tr('plane_tr',tau,dtau,rho,phi)

    if( debug ) then
       block
         real(REAL64), dimension(npar,npar) :: cor
         real(REAL64) :: q
         integer :: i,j

         cor = 1
         do i = 1, npar
            do j = i+1, npar
               if( cov(i,i) > 0 .and. cov(j,j) > 0 ) then
                  q = sqrt(cov(i,i)*cov(j,j))
                  cor(i,j) = cov(i,j) / q
               else
                  cor(i,j) = huge(cor)
               end if
               cor(j,i) = cor(i,j)
            end do
         end do

         write(error_unit,'(a,2(es12.4,1x),t48,a,2(f7.4,1x))') &
              'Cov: ', real(cov(1,:)),'corr: ', cor(1,:)
         do i = 2,npar
            write(error_unit,'(5x,2(es12.4,1x),t48,6x,2(f7.4,1x))') &
                 real(cov(i,:)),cor(i,:)
         end do

       end block
    end if

  end subroutine onplane

  subroutine planelike(k,l,x,fvec,fjac,ldfjac,iflag)

    integer, intent(in) :: k, l, ldfjac
    integer, intent(in out) :: iflag
    real(REAL64), dimension(k), intent(in) :: x
    real(REAL64), dimension(l), intent(in out) :: fvec
    real(REAL64), dimension(ldfjac,k), intent(in out) :: fjac

    logical, dimension(3), parameter :: mask = [ .true., .true., .false. ]
    real(REAL64), dimension(3) :: psi
    real(REAL64), dimension(3,3) :: dpsi
    real(REAL64) :: tau, rho, fi

    tau = x(1)
    rho = x(2)
    fi = phifix

    if( iflag == 0 ) then
       if( debug_iter ) then
          write(error_unit,*) 'x = ',real(x)
          write(error_unit,*) 'fvec = ',real(fvec)
          write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
          write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
       end if

       ! convergence check
       if( tau < taumin .or. tau > taumax .or. abs(rho) > rmax  ) iflag = -1

       return
    end if

    if( iflag == 1 ) then

       psi = psisum(tt,mju,tau,rho,fi,tmin,tmax,kmean,mask)
       fvec = psi(1:2)

    else if( iflag == 2 ) then

       dpsi = dpsisum(tt,mju,tau,rho,fi,tmin,tmax,kmean,mask)
       fjac = dpsi(1:2,1:2)

    end if

  end subroutine planelike


!!$  subroutine maxlike2(k,l,x,fvec,fjac,ldfjac,iflag)
!!$
!!$    integer, intent(in) :: k, l, ldfjac
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$    real(REAL64), dimension(ldfjac,k), intent(out) :: fjac
!!$
!!$    real(REAL64), dimension(3) :: sums, psi, psimin, psimax, bmax, bmin, u
!!$    real(REAL64), dimension(3,3) :: dsums, dpsi, dpsimin, dpsimax
!!$    real(REAL64) :: mu, t, tau, r, fi, sinfi, cosfi
!!$    integer :: n, i, j
!!$
!!$    mu = mju(1)
!!$
!!$    n = size(tt)
!!$    tau = x(1)
!!$    r = x(2)
!!$    fi = x(3)
!!$    cosfi = cos(fi)
!!$    sinfi = sin(fi)
!!$
!!$    if( .not. ( all(tmin <= tau) .and. all(tau <= tmax) .and. abs(r) < rmax) ) then
!!$       iflag = -1
!!$!       write(*,*) tmin,tau,tmax,r,rmax
!!$       return
!!$    end if
!!$
!!$
!!$
!!$!    if( .not. (tau > 0 .and. r >= 0) ) then
!!$!       if( debug_iter ) &
!!$!            write(error_unit,*) 'maxlike: assumptions violated.'
!!$!       write(*,*) tau,r
!!$!       iflag = -1
!!$!       return
!!$!    end if
!!$
!!$    if( iflag == 0 ) then
!!$       if( debug_iter ) then
!!$          write(error_unit,*) 'x = ',real(x)
!!$          write(error_unit,*) 'fvec = ',real(fvec)
!!$          write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
!!$          write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
!!$          write(error_unit,*) 'fjac(3,:) = ',real(fjac(3,:))
!!$!          write(error_unit,*) 'change:',x - maxlike_old
!!$       end if
!!$!       if( abs(x(1) - maxlike_old(1))/x(1) >
!!$
!!$       !       maxlike_old = x
!!$
!!$       return
!!$    end if
!!$
!!$    psimin = psilike(tmin,mu,tau,r,cosfi,sinfi) !/ tanh(bfac/2)
!!$    psimax = psilike(tmax,mu,tau,r,cosfi,sinfi) !/ tanh(bfac/2)
!!$!    bmin = (bfac/2)/(ttiny - tmin)
!!$!    bmax = (bfac/2)/(thuge - tmax)
!!$!    dpsimin = dpsilike(tmin,mu,tau,r,cosfi,sinfi)
!!$!    dpsimax = dpsilike(tmax,mu,tau,r,cosfi,sinfi)
!!$
!!$!    do i = 1, 3
!!$!       bmin(i) = abs(dpsimin(i,i) /psimin(i))
!!$!       bmax(i) = abs(dpsimax(i,i) /psimax(i))
!!$!    end do
!!$
!!$    if( iflag == 1 ) then
!!$
!!$       sums = real(0.0,REAL64)
!!$       do concurrent (i = 1:n)  !i = 1, n
!!$          t = tt(i)
!!$!          if( t <= ttiny ) then
!!$!             psi = real(0.0,REAL64)
!!$!             psi = psimin
!!$!          else if( t >= thuge ) then
!!$!             psi = real(0.0,REAL64)
!!$!             psi = psimax
!!$          !          else
!!$          if( t < tmin(iset(i)) ) then
!!$!             psi = psimin*tanh(bmin*(ttiny - t))
!!$             psi = psimin
!!$
!!$!             dpsi = 0
!!$!             u = bmin * (t - tmin)
!!$!             psi = min(psimin*tanh(u) + psimin,real(0.0,REAL64))
!!$!             write(*,*) real(u)
!!$!             write(*,*) real(psi)
!!$!             write(*,*) real(psimin)
!!$!             write(*,*) real(bmin)
!!$!             do j = 1, 3
!!$!                if( psi(j) <= 0 ) then
!!$!                   dpsi(j,j) = -bmin(j)*psimin(j)/cosh(u(j))**2
!!$!                else
!!$!                   psi(j) = 0
!!$!                end if
!!$!             end do
!!$!             where( psi <= 0 )
!!$!                dpsi = -bmin*psimin/cosh(u)**2
!!$!             elsewhere
!!$!                psi = 0
!!$!                dpsi = 0
!!$!             end where
!!$             psi = psimin
!!$
!!$          else if( t > tmax ) then
!!$!             psi = psimax*tanh(bmax*(thuge - t))
!!$
!!$!             dpsi = 0
!!$!             u = bmax * (tmax - t)
!!$!             psi = max(psimax*tanh(u) + psimax,real(0.0,REAL64))
!!$!             do j = 1, 3
!!$!                if( psi(j) >= 0 ) then
!!$!                   dpsi(j,j) = -bmax(j)*psimax(j)/cosh(u(j))**2
!!$!                else
!!$!                   psi(j) = 0
!!$!                end if
!!$!             end do
!!$!             where( psi >= 0 )
!!$!                dpsi = -bmax*psimax/cosh(u)**2
!!$!             elsewhere
!!$!                psi = 0
!!$!                dpsi = 0
!!$!             end where
!!$             psi = psimax
!!$
!!$          else
!!$             psi = psilike(t,mu,tau,r,cosfi,sinfi)
!!$          end if
!!$!          psi = psilike(t,mu,tau,r,cosfi,sinfi)
!!$
!!$          sums = sums + psi
!!$
!!$       end do
!!$
!!$       fvec = - sums / n
!!$
!!$    end if
!!$
!!$    if( iflag == 2 ) then
!!$
!!$!       dpsimin = dpsilike(tmin,mu,tau,r,cosfi,sinfi)
!!$!       dpsimax = dpsilike(tmax,mu,tau,r,cosfi,sinfi)
!!$
!!$       dsums = real(0.0,REAL64)
!!$       do concurrent (i = 1:n)  !i = 1, n
!!$          t = tt(i)
!!$!          if( t <= ttiny ) then
!!$!             dpsi = real(0.0,REAL64)
!!$!             psi = psimin
!!$!          else if( t >= thuge ) then
!!$!             dpsi = real(0.0,REAL64)
!!$!             psi = psimax
!!$          !          else
!!$          if( t < tmin(iset(i)) ) then
!!$             dpsi = real(0.0,REAL64)
!!$!            dpsi = psimin*tanh(bmin*(ttiny - t))
!!$             !            psi = psimin
!!$!             dpsi = dpsimin / cosh(bmin*(ttiny - t))**2 * (-bmin)
!!$!             do concurrent (j = 1:3)
!!$!                dpsi(:,j) = psimin(j) / cosh(bmin*(ttiny - t))**2 * (-bmin)
!!$!             end do
!!$
!!$!             dpsi = 0
!!$!             u = bmin * (t - tmin)
!!$!             psi = psimin*tanh(u) + psimin
!!$             do j = 1,-3
!!$                if( psi(j) <= 0 ) then
!!$                   dpsi(j,j) = -bmin(j)*psimin(j)/cosh(u(j))**2
!!$                else
!!$                   psi(j) = 0
!!$!                dpsi = 0
!!$                end if
!!$             end do
!!$
!!$          else if( t > tmax ) then
!!$             dpsi = real(0.0,REAL64)
!!$!            dpsi = psimax*tanh(bmax*(thuge - t))
!!$             !             psi = psimax
!!$             !             dpsi = psimax / cosh(bmax*(thuge - t))**2 * (-bmax)
!!$!             dpsi = dpsimax / cosh(bmax*(thuge - t))**2 * (-bmax)
!!$!             do concurrent (j = 1:3)
!!$!                dpsi(:,j) = psimax(j) / cosh(bmax*(thuge - t))**2 * (-bmax)
!!$             !             end do
!!$
!!$!             dpsi = 0
!!$             u = bmin * (t - tmin)
!!$!             psi = psimin*tanh(u) + psimin
!!$             do j = 1,-3
!!$                if( psi(j) <= 0 ) then
!!$                   dpsi(j,j) = -bmin(j)*psimin(j)/cosh(u(j))**2
!!$                else
!!$                   psi(j) = 0
!!$                end if
!!$                !                dpsi = 0
!!$             end do
!!$
!!$          else
!!$             dpsi = dpsilike(t,mu,tau,r,cosfi,sinfi)
!!$          end if
!!$!          dpsi = dpsilike(t,mu,tau,r,cosfi,sinfi)
!!$
!!$          dsums = dsums + dpsi
!!$       end do
!!$
!!$!       fjac(1,:) = dsums(1:3) / n
!!$!       fjac(2,2:3) = dsums(4:5) / n
!!$!       fjac(3,3) = dsums(6) / n
!!$
!!$       fjac = - dsums / n
!!$
!!$       do i = 2,3
!!$          do j = 1, i-1
!!$             fjac(i,j) = fjac(j,i)
!!$          end do
!!$       end do
!!$!       fjac = -fjac
!!$
!!$    end if
!!$
!!$  end subroutine maxlike2

  pure function psisum(t,mju,tau,rho,fi,tmin,tmax,kmean,mask) result(sums)

    real(REAL64), dimension(:), intent(in) :: t, mju, tmin,tmax, kmean
    real(REAL64), intent(in) :: tau,rho,fi
    logical, dimension(:), intent(in) :: mask
    real(REAL64), dimension(3) :: sums

    real(REAL64), dimension(:,:), allocatable :: psimin, psimax
    real(REAL64), dimension(3) :: psi
    real(REAL64) :: cosfi, sinfi
    integer :: i, n, k

    n = size(t)
    cosfi = cos(fi)
    sinfi = sin(fi)

    k = size(tmin)
    allocate(psimin(3,k),psimax(3,k))

    do concurrent (i = 1:k)
       psimin(:,i) = psilike(tmin(i),kmean(i),tau,rho,cosfi,sinfi,mask)
       psimax(:,i) = psilike(tmax(i),kmean(i),tau,rho,cosfi,sinfi,mask)
    end do

    sums = 0
    do concurrent (i = 1:n)
       if( t(i) < tmin(iset(i)) ) then
          psi = psimin(:,iset(i))
       else if( t(i) > tmax(iset(i)) ) then
          psi = psimax(:,iset(i))
       else
          psi = psilike(t(i),mju(i),tau,rho,cosfi,sinfi,mask)
       end if
!       psi = psilike(t(i),mju(i),tau,rho,cosfi,sinfi)

       sums = sums + psi
    end do

    sums = -sums / n

  end function psisum

  pure function dpsisum(t,mju,tau,rho,fi,tmin,tmax,kmean,mask) result(dsums)

    real(REAL64), dimension(:), intent(in) :: t, mju, tmin,tmax, kmean
    real(REAL64), intent(in) :: tau,rho,fi
    real(REAL64), dimension(3,3) :: dsums
    logical, dimension(:), intent(in) :: mask

    real(REAL64), dimension(:,:,:), allocatable :: dpsimin, dpsimax
    real(REAL64), dimension(3,3) :: dpsi
    real(REAL64) :: cosfi, sinfi
    integer :: i, j, k, n

    n = size(t)
    cosfi = cos(fi)
    sinfi = sin(fi)

    k = size(tmin)
    allocate(dpsimin(3,3,k),dpsimax(3,3,k))
    do concurrent (i = 1:k)
       dpsimin(:,:,i) = dpsilike(tmin(i),kmean(i),tau,rho,cosfi,sinfi,mask)
       dpsimax(:,:,i) = dpsilike(tmax(i),kmean(i),tau,rho,cosfi,sinfi,mask)
    end do

    dsums = 0
    do concurrent (i = 1:size(t))
       if( t(i) < tmin(iset(i)) ) then
          dpsi = dpsimin(:,:,iset(i))
       else if( t(i) > tmax(iset(i)) ) then
          dpsi = dpsimax(:,:,iset(i))
       else
          dpsi = dpsilike(t(i),mju(i),tau,rho,cosfi,sinfi,mask)
       end if
!       dpsi = dpsilike(t(i),mju(i),tau,rho,cosfi,sinfi)

       dsums = dsums + dpsi
    end do

    dsums = -dsums / n

    do concurrent (i = 2:3)
       do concurrent (j = 1:i-1)
          dsums(i,j) = dsums(j,i)
       end do
    end do

  end function dpsisum

  pure function psilike(t,mu,tau,rho,cosfi,sinfi,mask) result(psi)

    real(REAL64), intent(in) :: t,mu,tau,rho,cosfi,sinfi
    logical, dimension(:), intent(in) :: mask

    real(REAL64), dimension(3) :: psi
    real(REAL64) :: c2, cs, tt, sig, beta, res, sig_tau, beta_tau, sig_rho, &
         beta_rho, sig_fi, beta_fi

    tt = t / tau
    c2 = cosfi**2
    cs = cosfi*sinfi
    sig = sqrt(mu*(tau/sy + t**2/sx) + rho**2*(1 + (tt**2 - 1)*c2))
    beta = mu*tau*(t/sx + 1/sy) + rho**2*(1 + (tt - 1)*c2)
    res = (mu/sig)*(tau - t)

    if( mask(1) ) then
       sig_tau = mu/(2*sig*sy) - rho*((tt*cosfi)**2/tau)*(rho/sig)
       beta_tau = mu*(t/sx + 1/sy) - rho**2*tt/tau*c2
       psi(1) = -res*(mu/sig) + (res**2 - 3)*(sig_tau/sig) + beta_tau/beta
    end if

    if( mask(2) ) then
       sig_rho = (rho/sig)*(1 + (tt**2 - 1)*c2)
       beta_rho = 2*rho*(1 + (tt - 1)*c2)
       psi(2) = (res**2 - 3)*(sig_rho/sig) + beta_rho/beta
    end if

    if( mask(3) ) then
       sig_fi = -cs*(tt**2 - 1)*(rho/sig)*rho
       beta_fi = -2*rho**2*(tt - 1)*cs
       psi(3) = (res**2 - 3)*(sig_fi/sig) + beta_fi/beta
    end if

  end function psilike

  pure function dpsilike(t,mu,tau,rho,cosfi,sinfi,mask) result(dpsi)

    real(REAL64), intent(in) :: t,mu,tau,rho,cosfi,sinfi
    logical, dimension(:), intent(in) :: mask

    real(REAL64), dimension(3,3) :: dpsi
    real(REAL64) :: c2, cs, cf, sf, sig, beta, res, sig_t, beta_t, sig_r, &
         beta_r, sig_f, beta_f, sig_tt, sig_tr, sig_tf, sig_rr, &
         sig_rf, sig_ff, beta_tt, beta_tr, beta_tf, beta_rr, beta_rf, &
         beta_ff, dpsi_tt, dpsi_tr, dpsi_tf, dpsi_rr, dpsi_rf, dpsi_ff, &
         res_t, res_r, res_f, res_tt, res_tr, res_tf, res_rr, res_rf, &
         res_ff

    c2 = cosfi**2
    cs = cosfi*sinfi
    cf = cosfi
    sf = sinfi

    sig = sqrt(mu*(tau/sy + t**2/sx) + rho**2*(1 + ((t/tau)**2 - 1)*c2))
    beta = mu*tau*(t/sx + 1/sy) + rho**2*(1 + (t/tau - 1)*c2)
    res = (mu/sig)*(tau - t)

    if( mask(1) ) then
       sig_t = mu/(2*sig*sy) - rho*((t*cosfi/tau)**2/tau)*(rho/sig)
       beta_t = mu*(t/sx + 1/sy) - rho**2*t/tau**2*c2
       res_t = mu/sig - res*(sig_t/sig)
    end if

    if( mask(2) ) then
       sig_r = (rho/sig)*(1 + ((t/tau)**2 - 1)*c2)
       beta_r = 2*rho*(1 + (t/tau - 1)*c2)
       res_r = -res*(sig_r/sig)
    end if

    if( mask(3) ) then
       sig_f = -cs*((t/tau)**2 - 1)*(rho/sig)*rho
       beta_f = -2*rho**2*(t/tau - 1)*cs
       res_f = -res*(sig_f/sig)
    end if

    if( mask(1) ) then
       sig_tt = -sig_t**2/sig + 3*(rho**2/sig)*(t/tau**2)**2*c2
       beta_tt = 2*rho**2*t/tau**3*c2
       res_tt = -(mu/sig)*(sig_t/sig) &
            -res_t*(sig_t/sig) - res*(sig_tt/sig) + res*(sig_t/sig)**2
       dpsi_tt = -res*res_tt - res_t**2 &
            + beta_tt/beta - (beta_t/beta)**2 - 3*sig_tt/sig + 3*(sig_t/sig)**2
       dpsi(1,1) = dpsi_tt
    end if

    if( mask(2) ) then
       sig_tr = -sig_t*sig_r/sig - 2*(rho/sig)*(t/tau)**2/tau*c2
       sig_rr = (1 + ((t/tau)**2 - 1)*c2)*(1 - (rho/sig)*sig_r)/sig
       beta_rr = 2*(1 + (t/tau - 1)*c2)
       beta_tr = -2*rho*t/tau**2*c2
       res_tr = -(mu/sig)*(sig_r/sig) &
            -res_r*(sig_t/sig) - res*(sig_tr/sig) + res*(sig_t/sig)*(sig_r/sig)
       res_rr = -res_r*(sig_r/sig) - res*(sig_rr/sig) + res*(sig_r/sig)**2
       dpsi_tr = -res*res_tr - res_t*res_r &
            + beta_tr/beta - (beta_t/beta)*(beta_r/beta) &
            - 3*sig_tr/sig + 3*(sig_t/sig)*(sig_r/sig)
       dpsi_rr = -res*res_rr - res_r**2 &
            + beta_rr/beta - (beta_r/beta)**2 - 3*sig_rr/sig + 3*(sig_r/sig)**2
       dpsi(1,2) = dpsi_tr
       dpsi(2,2) = dpsi_rr
    end if

    if( mask(3) ) then
       sig_tf = -sig_t*sig_f/sig + 2*(rho**2/sig)*(t/tau)**2/tau*cs
       sig_rf = -sig_r*sig_f/sig - 2*(rho/sig)*((t/tau)**2 - 1)*cs
       sig_ff = -rho*(rho/sig)*((t/tau)**2-1)*(-(sig_f/sig)*cs + cf**2 - sf**2)
       beta_tf = 2*rho**2*t/tau**2*cs
       beta_rf = -4*rho*(t/tau - 1)*cs
       beta_ff = -2*rho**2*(t/tau - 1)*(cf**2 - sf**2)
       res_tf = -(mu/sig)*(sig_f/sig) &
            -res_f*(sig_t/sig) - res*(sig_tf/sig) + res*(sig_t/sig)*(sig_f/sig)
       res_rf = -res_f*(sig_r/sig) &
            - res*(sig_rf/sig) + res*(sig_r/sig)*(sig_f/sig)
       res_ff = -res_f*(sig_f/sig) - res*(sig_ff/sig) + res*(sig_f/sig)**2
       dpsi_tf = -res*res_tf - res_t*res_f &
            + beta_tf/beta - (beta_t/beta)*(beta_f/beta) &
            - 3*sig_tf/sig + 3*(sig_t/sig)*(sig_f/sig)
       dpsi_rf = -res*res_rf - res_r*res_f &
            + beta_rf/beta - (beta_r/beta)*(beta_f/beta) &
            - 3*sig_rf/sig + 3*(sig_r/sig)*(sig_f/sig)
       dpsi_ff = -res*res_ff - res_f**2 &
            + beta_ff/beta - (beta_f/beta)**2 - 3*sig_ff/sig + 3*(sig_f/sig)**2
       dpsi(1,3) = dpsi_tf
       dpsi(2,3) = dpsi_rf
       dpsi(3,3) = dpsi_ff
    end if

    ! only the upper triangle is filled

!    dpsi(1,1:3) = [ dpsi_tt, dpsi_tr, dpsi_tf ]
!    dpsi(2,2:3) = [ dpsi_rr, dpsi_rf ]
!    dpsi(3,3) = dpsi_ff

  end function dpsilike

!!$  pure function psilike1(t,mu,tau,r,cosfi,sinfi) result(psi)
!!$
!!$    real(REAL64), intent(in) :: t,mu,tau,r,cosfi,sinfi
!!$    real(REAL64), dimension(3) :: psi
!!$
!!$    real(REAL64) :: c2, cs, tt, sig, beta, res, sig_tau, beta_tau, sig_r, &
!!$         beta_r, sig_fi, beta_fi!, res_t, res_r, res_f
!!$
!!$    tt = t / tau
!!$    c2 = cosfi**2
!!$    cs = cosfi*sinfi
!!$    sig = sqrt(mu*(tau + t**2) + r**2*(1 + (tt**2 - 1)*c2))
!!$    beta = mu*tau*(t + 1) + r**2*(1 + (tt - 1)*c2)
!!$    res = (mu/sig)*(tau - t)
!!$
!!$    sig_tau = mu/(2*sig) - r*((tt*cosfi)**2/tau)*(r/sig)
!!$    beta_tau = mu*(t + 1) - r**2*tt/tau*c2
!!$    psi(1) = -res*(mu/sig) + (res**2 - 3)*(sig_tau/sig) + beta_tau/beta
!!$
!!$    sig_r = (r/sig)*(1 + (tt**2 - 1)*c2)
!!$    beta_r = 2*r*(1 + (tt - 1)*c2)
!!$    psi(2) = (res**2 - 3)*(sig_r/sig) + beta_r/beta
!!$
!!$    sig_fi = -cs*(tt**2 - 1)*(r/sig)*r
!!$    beta_fi = -2*r**2*(tt - 1)*cs
!!$    psi(3) = (res**2 - 3)*(sig_fi/sig) + beta_fi/beta
!!$
!!$!    res_t = mu/sig - res*(sig_tau/sig)
!!$!    res_r = -res*(sig_r/sig)
!!$!    res_f = -res*(sig_fi/sig)
!!$
!!$!    psi(1) = -res*res_t + beta_tau/beta - 3*(sig_tau/sig)
!!$!    psi(2) = -res*res_r + beta_r/beta - 3*(sig_r/sig)
!!$!    psi(3) = -res*res_f + beta_fi/beta - 3*(sig_fi/sig)
!!$
!!$  end function psilike1
!!$
!!$  pure function dpsilike1(t,mu,tau,r,cosfi,sinfi) result(dpsi)
!!$
!!$    real(REAL64), intent(in) :: t,mu,tau,r,cosfi,sinfi
!!$    real(REAL64), dimension(3,3) :: dpsi
!!$
!!$    real(REAL64) :: c2, cs, cf, sf, sig, beta, res, sig_t, beta_t, sig_r, &
!!$         beta_r, sig_f, beta_f, sig_tt, sig_tr, sig_tf, sig_rr, &
!!$         sig_rf, sig_ff, beta_tt, beta_tr, beta_tf, beta_rr, beta_rf, &
!!$         beta_ff, dpsi_tt, dpsi_tr, dpsi_tf, dpsi_rr, dpsi_rf, dpsi_ff, &
!!$         res_t, res_r, res_f, res_tt, res_tr, res_tf, res_rr, res_rf, &
!!$         res_ff
!!$
!!$    c2 = cosfi**2
!!$    cs = cosfi*sinfi
!!$    cf = cosfi
!!$    sf = sinfi
!!$
!!$    sig = sqrt(mu*(tau + t**2) + r**2*(1 + ((t/tau)**2 - 1)*c2))
!!$    beta = mu*tau*(t + 1) + r**2*(1 + (t/tau - 1)*c2)
!!$    res = (mu/sig)*(tau - t)
!!$
!!$    sig_t = mu/(2*sig) - r*((t*cosfi/tau)**2/tau)*(r/sig)
!!$    beta_t = mu*(t + 1) - r**2*t/tau**2*c2
!!$    res_t = mu/sig - res*(sig_t/sig)
!!$
!!$    sig_r = (r/sig)*(1 + ((t/tau)**2 - 1)*c2)
!!$    beta_r = 2*r*(1 + (t/tau - 1)*c2)
!!$    res_r = -res*(sig_r/sig)
!!$
!!$    sig_f = -cs*((t/tau)**2 - 1)*(r/sig)*r
!!$    beta_f = -2*r**2*(t/tau - 1)*cs
!!$    res_f = -res*(sig_f/sig)
!!$
!!$
!!$
!!$!    sig_tt = -(mu/sig)*(sig_t/sig)/2 + &
!!$!         (t/tau)**2/tau*c2*(3*sig/tau + sig_t)*(r/sig)**2
!!$    !    sig_tt = 3*r*(r/sig)*(t/tau**2)**2*c2 - (sig_t/sig)*sig_t*()
!!$    sig_tt = -(mu/sig)*(sig_t/sig)/2  + &
!!$         r*(r/sig)*(t**2/tau**3)*c2*(sig_t/sig + 3/tau)
!!$    sig_tr = -(mu/sig)*(sig_r/sig)/2  + &
!!$         (r/sig)*(t**2/tau**3)*c2*((r/sig)*sig_r - 2)
!!$    sig_tf = -(mu/sig)*(sig_f/sig)/2  + &
!!$         r*(r/sig)*(t**2/tau**3)*((sig_f/sig)*c2 + 2*cs)
!!$    sig_rr = (1 + ((t/tau)**2 - 1)*c2)*(1 - (r/sig)*sig_r)/sig
!!$    sig_rf = -(r/sig) &
!!$         * (2*((t/tau)**2 - 1)*cs + (1 + ((t/tau)**2 - 1)*c2)*(sig_f/sig))
!!$    sig_ff = -r*(r/sig)*((t/tau)**2-1)*(-(sig_f/sig)*cs + cf**2 - sf**2)
!!$
!!$!    block
!!$!      real(REAL64) :: d, tau1, sig_t1
!!$!      d = 1e-7
!!$!      tau1 = tau + d
!!$!      sig = sqrt(mu*(tau1 + t**2) + r**2*(1 + ((t/tau1)**2 - 1)*c2))
!!$!      beta = mu*tau1*(t + 1) + r**2*(1 + (t/tau1 - 1)*c2)
!!$!      res = (mu/sig)*(tau1 - t)
!!$!      sig_t1 = mu/(2*sig) - r*((t*cosfi/tau1)**2/tau1)*(r/sig)
!!$!      write(*,*) sig_tt, (sig_t1 - sig_t)/d
!!$!      stop 0
!!$!    end block
!!$
!!$
!!$!    sig_tr = -2*(t/tau)**2/tau*c2*(r/sig) - sig_t*(sig_r/sig)
!!$!    sig_tf = 2*cs*(t/tau)**2/tau*r*(r/sig) - sig_t*(sig_f/sig)
!!$!    sig_rr = (1+((t/tau)**2-1)*c2)*(1 - r*(sig_r/sig))/sig
!!$!    sig_rf = (1-(t/tau)**2)*cs*(2*(r/sig) - (r/sig)**2*sig_r)
!!$!    sig_ff = r*(r/sig)*(1-(t/tau)**2)*((cf**2-sf**2)-cf*sf*(sig_f/sig))
!!$
!!$    beta_tt = 2*r**2*t/tau**3*c2
!!$    beta_tr = -2*r*t/tau**2*c2
!!$    beta_tf = 2*r**2*t/tau**2*cs
!!$    beta_rr = 2*(1 + (t/tau - 1)*c2)
!!$    beta_rf = -4*r*(t/tau - 1)*cs
!!$    beta_ff = -2*r**2*(t/tau - 1)*(cf**2 - sf**2)
!!$
!!$    res_tt = -(mu/sig)*(sig_t/sig) &
!!$             -res_t*(sig_t/sig) - res*(sig_tt/sig) + res*(sig_t/sig)**2
!!$    res_tr = -(mu/sig)*(sig_r/sig) &
!!$             -res_r*(sig_t/sig) - res*(sig_tr/sig) + res*(sig_t/sig)*(sig_r/sig)
!!$    res_tf = -(mu/sig)*(sig_f/sig) &
!!$             -res_f*(sig_t/sig) - res*(sig_tf/sig) + res*(sig_t/sig)*(sig_f/sig)
!!$    res_rr = -res_r*(sig_r/sig) - res*(sig_rr/sig) + res*(sig_r/sig)**2
!!$    res_rf = -res_f*(sig_r/sig) - res*(sig_rf/sig) + res*(sig_r/sig)*(sig_f/sig)
!!$    res_ff = -res_f*(sig_f/sig) - res*(sig_ff/sig) + res*(sig_f/sig)**2

!!$    block
!!$      real(REAL64) :: d, tau1, sig_t1, res_t1, res_f1, fi1, sig_f1, fi, beta_f1
!!$      d = 1e-7
!!$      tau1 = tau! + d
!!$      fi = atan2(sinfi,cosfi)
!!$      fi1 = fi + d
!!$      c2 = cos(fi1)**2
!!$      cs = cos(fi1)*sin(fi1)
!!$      cf = cos(fi1)
!!$      sf = sin(fi1)
!!$
!!$      sig = sqrt(mu*(tau1 + t**2) + r**2*(1 + ((t/tau1)**2 - 1)*c2))
!!$      beta = mu*tau1*(t + 1) + r**2*(1 + (t/tau1 - 1)*c2)
!!$      res = (mu/sig)*(tau1 - t)
!!$!      sig_t1 = mu/(2*sig) - r*((t*cosfi/tau1)**2/tau1)*(r/sig)
!!$      !      res_t1 = mu/sig - res*(sig_t1/sig)!
!!$      sig_f1 = -cs*((t/tau)**2 - 1)*(r/sig)*r
!!$      res_f1 = -res*(sig_f1/sig)
!!$      beta_f1 = -2*r**2*(t/tau - 1)*cs
!!$      write(*,*) res_ff, (res_f1 - res_f)/d
!!$      write(*,*) beta_ff, (beta_f1 - beta_f)/d
!!$      stop 0
!!$    end block


!!$    dpsi_tt = -res*res_tt - res_t**2 &
!!$         + beta_tt/beta - (beta_t/beta)**2 - 3*sig_tt/sig + 3*(sig_t/sig)**2
!!$    dpsi_tr = -res*res_tr - res_t*res_r &
!!$         + beta_tr/beta - (beta_t/beta)*(beta_r/beta) &
!!$         - 3*sig_tr/sig + 3*(sig_t/sig)*(sig_r/sig)
!!$    dpsi_tf = -res*res_tf - res_t*res_f &
!!$         + beta_tf/beta - (beta_t/beta)*(beta_f/beta) &
!!$         - 3*sig_tf/sig + 3*(sig_t/sig)*(sig_f/sig)
!!$    dpsi_rr = -res*res_rr - res_r**2 &
!!$         + beta_rr/beta - (beta_r/beta)**2 - 3*sig_rr/sig + 3*(sig_r/sig)**2
!!$    dpsi_rf = -res*res_rf - res_r*res_f &
!!$         + beta_rf/beta - (beta_r/beta)*(beta_f/beta) &
!!$         - 3*sig_rf/sig + 3*(sig_r/sig)*(sig_f/sig)
!!$    dpsi_ff = -res*res_ff - res_f**2 &
!!$         + beta_ff/beta - (beta_f/beta)**2 - 3*sig_ff/sig + 3*(sig_f/sig)**2
!!$
!!$    ! Warning: only the upper triangle is filled !
!!$
!!$    dpsi(1,1:3) = [ dpsi_tt, dpsi_tr, dpsi_tf ]
!!$    dpsi(2,2:3) = [ dpsi_rr, dpsi_rf ]
!!$    dpsi(3,3) = dpsi_ff
!!$
!!$    dpsi_tt = -(mu/sig + res*(sig_t/sig))**2 + &
!!$              res*((mu/sig)*(sig_t/sig) + (mu/sig)*(sig_t/sig) &
!!$              +2*res*(sig_t/sig)**2 + res*(sig_tt/sig)) &
!!$              -(beta_t/beta)**2 + beta_tt/beta &
!!$              +3*((sig_t/sig)**2 - sig_tt/sig)
!!$
!!$    dpsi_tr = (mu/sig - res*(sig_t/sig))*res*(sig_r/sig) &
!!$            - res*(-(mu/sig)*(sig_r/sig) - &
!!$            res*(-2*(sig_t/sig)*(sig_r/sig) + sig_tr/sig)) &
!!$            -(beta_t/beta)*(beta_r/beta) + beta_tr/beta &
!!$            + 3*(sig_t/sig)*(sig_r/sig) - sig_tr/sig
!!$
!!$    dpsi_tf = (mu/sig - res*(sig_t/sig))*res*(sig_f/sig) &
!!$            +res*((mu/sig)*(sig_f/sig) + &
!!$            res*(-2*(sig_t/sig)*(sig_f/sig) + sig_tf/sig)) &
!!$            -(beta_t/beta)*(beta_f/beta) + beta_tf/beta &
!!$            + 3*((sig_t/sig)*(sig_f/sig) - sig_tf/sig)
!!$
!!$    dpsi_rr = -((res/sig)*(sig_r/sig))**2 &
!!$            +res*(-2*res*(sig_r/sig)**2 + res*(sig_rr/sig)) &
!!$            - (beta_r/beta)**2 + beta_rr/beta &
!!$            +3*((sig_r/sig)**2 - sig_rr/sig)
!!$
!!$    dpsi_rf = -res**2*(sig_r/sig)*(sig_f/sig) &
!!$            +res**2*(-2*(sig_r/sig)*(sig_f/sig) + sig_rf/sig) &
!!$            - (beta_r/beta)*(beta_f/beta) + beta_rf/beta &
!!$            + 3*((sig_r/sig)*(sig_f/sig) - sig_rf/sig)
!!$
!!$    dpsi_ff = -(res*(sig_f/sig))**2 &
!!$            +res*(-2*res*(sig_f/sig)**2 + res*(sig_ff/sig)) &
!!$            - (beta_f/beta)**2 + beta_ff/beta &
!!$            + 3*((sig_f/sig)**2 - sig_ff/sig)

!    dpsi = [ dpsi_tt, dpsi_tr, dpsi_rf, dpsi_rr, dpsi_rf, dpsi_ff ]

!!$  end function dpsilike1


!!$  subroutine maxlike0(k,l,x,fvec,fjac,ldfjac,iflag)
!!$
!!$    integer, intent(in) :: k, l, ldfjac
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$    real(REAL64), dimension(ldfjac,k), intent(out) :: fjac
!!$    real(REAL64), dimension(:), allocatable :: sig, psi, dpsi, beta, res, &
!!$         sig_t, sig_s, sig_f, beta_t, beta_s, beta_f, &
!!$         sig_tt,sig_ts,sig_tf,sig_ss,sig_sf,sig_ff, &
!!$         beta_tt,beta_ts,beta_tf,beta_ss,beta_sf,beta_ff
!!$    real(REAL64) :: t, s, f, c2, cs, sf, cf
!!$    integer :: n, i, j, stat
!!$    character(len=80) :: errmsg
!!$
!!$    n = size(tt)
!!$    t = x(1)
!!$    s = x(2)
!!$    f = x(3)
!!$    sf = sin(f)
!!$    cf = cos(f)
!!$    c2 = cf**2
!!$    cs = cf*sf
!!$    cs = sin(2*f)/2
!!$
!!$!    write(*,*) 'sss',t,s,iflag
!!$
!!$    if( .not. (t > 0 .and. s >= 0) ) then
!!$       if( debug_iter ) &
!!$            write(error_unit,*) 'maxlike: assumptions violated.'
!!$       write(*,*) t,s
!!$       iflag = -1
!!$       return
!!$    end if
!!$
!!$
!!$    allocate(sig(n),psi(n),dpsi(n),beta(n),res(n), &
!!$         sig_t(n), sig_s(n), sig_f(n), beta_t(n), beta_s(n), beta_f(n), &
!!$         stat=stat,errmsg=errmsg)
!!$    if( stat /= 0 ) then
!!$       write(error_unit,'(2a)') warning_label, trim(errmsg)
!!$       iflag = -1
!!$       return
!!$    end if
!!$
!!$    if( iflag == 0 .and. debug_iter ) then
!!$       write(error_unit,*) 'x = ',real(x)
!!$       write(error_unit,*) 'fvec = ',real(fvec)
!!$       write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
!!$       write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
!!$       write(error_unit,*) 'fjac(3,:) = ',real(fjac(3,:))
!!$       return
!!$    end if
!!$
!!$    sig = sqrt(mju*(t + tt**2) + s**2*(1 + ((tt/t)**2 - 1)*c2))
!!$    beta = (mju*t*(tt + 1) + s**2*(1 + (tt/t - 1)*c2))
!!$    sig_t = (mju - 2*s**2*(tt/t)**2/t*c2)/(2*sig)
!!$    sig_s = s*(1 + ((tt/t)**2 - 1)*c2)/sig
!!$    sig_f = -s**2*cs*((tt/t)**2 - 1)/sig
!!$    beta_t = mju*(tt + 1) - s**2*tt/t**2*c2
!!$    beta_s = 2*s*(1 + (tt/t - 1)*c2)
!!$    beta_f = -2*s**2*(tt/t - 1)*cs
!!$    res = mju*(t - tt)
!!$    psi = res/sig
!!$
!!$    if( iflag == 1 ) then
!!$
!!$       fvec(1) = -sum(psi*(mju/sig)) + sum(psi*(res/sig)*(sig_t/sig)) + &
!!$            sum(beta_t/beta) - real(3.0,REAL64)*sum(sig_t/sig)
!!$!       write(*,*) sum(psi*(mju/sig)),sum(psi*(res/sig)*(sig_t/sig)),&
!!$!            sum(beta_t/beta,mask),real(3.0,REAL64)*sum(sig_t/sig,mask)
!!$!       write(*,*) sum(psi*(mju/sig))-sum(psi*mju/sig), &
!!$!            sum(psi*res*sig_t/sig**2)-sum(psi*(res/sig)*sig_t/sig),&
!!$!            sum(beta_t/beta)- &
!!$!            sum((mju*(tt + 1) - s**2*tt/t**2*c2)/(mju*t*(tt + 1) + s**2*(1 + (tt/t - 1)*c2)))
!!$       fvec(2) = sum(psi*(res/sig)*sig_s/sig) + &
!!$            sum(beta_s/beta) - real(3.0,REAL64)*sum(sig_s/sig)
!!$       fvec(3) = sum(psi*(res/sig)*sig_f/sig) + &
!!$            sum(beta_f/beta) - real(3.0,REAL64)*sum(sig_f/sig)
!!$!       write(*,*) sum(psi*res*sig_f/sig**2), &
!!$!            sum(beta_f/beta) , 3*sum(sig_f/sig)
!!$       !       write(*,*) 'fvec:',real(fvec)
!!$
!!$       fvec = -fvec
!!$    end if
!!$
!!$    if( iflag == 2 ) then
!!$
!!$       dpsi = 1
!!$
!!$       allocate(sig_tt(n),sig_ts(n),sig_tf(n),sig_ss(n),sig_sf(n),sig_ff(n),&
!!$            beta_tt(n),beta_ts(n),beta_tf(n),beta_ss(n),beta_sf(n),beta_ff(n),&
!!$            stat=stat,errmsg=errmsg)
!!$       if( stat /= 0 ) then
!!$          write(error_unit,'(2a)') warning_label, trim(errmsg)
!!$          iflag = -1
!!$          return
!!$       end if
!!$
!!$       sig_tt = -mju*sig_t/(2*sig**2) + &
!!$            s**2*(tt/t)**2/t*c2*(3*sig/t+sig_t)/sig**2
!!$       sig_ts = -2*s*(tt/t)**2/t*c2/sig - sig_t*sig_s/sig
!!$       sig_tf = -mju*sig_t/(2*sig**2) + &
!!$            s**2*(tt/t)**2/t*(sig_t*cf - 2*sf)/sig**2
!!$       sig_tf = cf*sf*s**2/sig*((sig_t/sig)*((tt/t)**2-1) + 2*(tt/t)**2/t)
!!$       sig_tf = 2*s**2*cf*sf*(tt/t)**2/t/sig - sig_t*sig_f/sig
!!$       sig_ss = (1+((tt/t)**2-1)*c2)*(sig - s*sig_s)/sig**2
!!$       sig_sf = (1-(tt/t)**2)*cf*sf*(2*s/sig - (s/sig)**2*sig_s)
!!$       sig_ff = s**2*(1-(tt/t)**2)*((cf**2-sf**2)/sig-cf*sf*sig_f/sig**2)
!!$
!!$       beta_tt = 2*s**2*tt/t**3*c2
!!$       beta_ts = -2*s*tt/t**2*c2
!!$       beta_tf = 2*s**2*tt/t**2*cf*sf
!!$       beta_ss = 2*(1 + (tt/t - 1)*c2)
!!$       beta_sf = -4*s*(tt/t - 1)*cf*sf
!!$       beta_ff = 2*s**2*(tt/t - 1)*(cf**2 - sf**2)
!!$
!!$
!!$       fjac(1,1) = -sum(dpsi*(mju/sig-res*sig_t/sig**2)**2) + &
!!$            sum(psi*(mju*sig_t/sig**2 + mju*sig_t/sig**2 &
!!$            -2*res*(sig_t/sig)**2/sig + res*sig_tt/sig**2)) &
!!$            +sum(-(beta_t/beta)**2 + beta_tt/beta) &
!!$            -3*sum(-(sig_t/sig)**2 + sig_tt/sig)
!!$       fjac(1,2) = sum(dpsi*(mju/sig - res*sig_t/sig**2)*res*sig_s/sig**2) &
!!$            - sum(psi*(-mju*sig_s/sig**2 - &
!!$            res*(-2*sig_t*sig_s/sig**3 + sig_ts/sig**2))) &
!!$            + sum(-beta_t*beta_s/beta**2+beta_ts/beta) &
!!$            - 3*sum(-sig_t*sig_s/sig**2 + sig_ts/sig)
!!$
!!$       fjac(1,3) = sum(dpsi*(mju/sig - res*sig_t/sig**2)*res*sig_f/sig**2) &
!!$            -sum(psi*(-mju*sig_f/sig**2 - &
!!$            res*(-2*sig_t*sig_f/sig**3 + sig_tf/sig**2))) &
!!$            + sum(-beta_t*beta_f/beta**2+beta_tf/beta) &
!!$            - 3*sum(-sig_t*sig_f/sig**2 + sig_tf/sig)
!!$
!!$!       write(*,*) sum(dpsi*(mju/sig - res*sig_t/sig**2)*res*sig_f/sig**2), &
!!$!            -sum(psi*(-mju*sig_f/sig**2 - &
!!$!            res*(-2*sig_t*sig_f/sig**3 + sig_tf/sig**2))), &
!!$!             sum(-beta_t*beta_f/beta**2+beta_tf/beta), &
!!$!             3*sum(-sig_t*sig_f/sig**2 + sig_tf/sig)
!!$
!!$       fjac(2,2) = -sum(dpsi*(res*sig_s/sig**2)**2) &
!!$            +sum(psi*(-2*res*(sig_s/sig)**2/sig + res*sig_ss/sig**2)) &
!!$            +sum(-(beta_s/beta)**2+beta_ss/beta) &
!!$            -3*sum(-(sig_s/sig)**2 + sig_ss/sig)
!!$       fjac(2,3) = -sum(dpsi*(res/sig**2)**2*sig_s*sig_f) &
!!$            +sum(psi*(res/sig)*(-2*(sig_s/sig)*(sig_f/sig) + sig_sf/sig)) &
!!$            + sum(-(beta_s/beta)*(beta_f/beta)+beta_sf/beta) &
!!$            - 3*sum(-(sig_s/sig)*(sig_f/sig) + sig_sf/sig)
!!$       fjac(3,3) = -sum(dpsi*(res*sig_f/sig**2)**2) &
!!$            +sum(psi*(-2*res*(sig_f/sig)**2/sig + res*sig_ff/sig**2)) &
!!$            +sum(-(beta_f/beta)**2+beta_ff/beta) &
!!$            -3*sum(-(sig_f/sig)**2 + sig_ff/sig)
!!$
!!$       do i = 2,3
!!$          do j = 1, i-1
!!$             fjac(i,j) = fjac(j,i)
!!$          end do
!!$       end do
!!$
!!$       fjac = -fjac
!!$    end if
!!$
!!$  end subroutine maxlike0



  subroutine online(tau,dtau,rho,fi,converge)

    ! This routine does to search of minimum on tau line.
    ! It's implemented as the standard Newton-Raphson method.

    real(REAL64), intent(in out) :: tau, dtau
    real(REAL64), intent(in) :: rho,fi
    logical, intent(out) :: converge

    logical, dimension(3), parameter :: mask = [ .true., .false., .false. ]
    real(REAL64), dimension(3) :: psi
    real(REAL64), dimension(3,3) :: dpsi
    real(REAL64) :: f,df,d,tol
    integer :: n, iter

    converge = .false.
    n = size(tt)

    tol = max(10*epsilon(tau)*(abs(tau) + 1), 1e-3*dtau/sqrt(real(n)))
    if( debug_iter ) write(error_unit,'(a,es12.2)') 'online tol: ',tol

    do iter = 1, 3*precision(tau)

       ! the function, and the derivative
       psi = psisum(tt,mju,tau,rho,fi,tmin,tmax,kmean,mask)
       f = psi(1)
       dpsi = dpsisum(tt,mju,tau,rho,fi,tmin,tmax,kmean,mask)
       df = dpsi(1,1)

       if( abs(df) < epsilon(df) ) exit ! zero derivative

       ! update
       d = f / df
       tau = tau - d

       if( debug_iter ) write(error_unit,'(a,i2,2x,g0,2es12.2,es14.4)') &
            "rate, incr., f, f': ",iter,real(tau),d,f,df

       converge = abs(d) < tol .and. df > 0

       ! exit of iterations:
       if( converge ) exit   ! a local minimum is found, success !
       if( df < 0 ) exit     ! a local maximum encountered, bad luck
       if( tau < taumin .or. tau > taumax ) exit ! out of range, divergent

    end do

    if( debug ) then

       if( abs(df) < epsilon(df) ) write(error_unit,'(2a)') &
               warning_label," fail in online(): derivative df/dtau == 0"

       if( df < 0 ) write(error_unit,'(2a)') warning_label,&
               " fail in online(): a local maximum is encountered."

       if( tau < taumin .or. tau > taumax ) write(error_unit,'(2a)') &
            warning_label," fail in online(): out of range."
    end if

    if( converge ) then
       dtau = 1/sqrt(df*(n-1))
    end if

!!$    block
!!$      integer :: n1,n2,i
!!$      n1 = count(tt < tmin(iset))
!!$      n2 = count(tt > tmax(iset))
!!$!      write(*,*) n, n1,n2, (n2*30 - n1*40)/df/(n1 + n2)
!!$
!!$      do i = 1,size(tmin)
!!$         write(*,*) psilike(tmin(i),kmean(i),tau,rho,cos(fi),sin(fi),mask)
!!$         write(*,*) psilike(tmax(i),kmean(i),tau,rho,cos(fi),sin(fi),mask)
!!$      end do
!!$
!!$    end block

    if( .not. converge .and. debug_graph ) &
         call graph_like_t('/tmp/on',tau,dtau,rho,fi)
!    if( debug_graph ) then
!       block
!       end block
!    end if

  end subroutine online

  subroutine bisection(tau,dtau,rho,fi,converge)

    ! This routine implements bisection as the failback

    real(REAL64), intent(in out) :: tau
    real(REAL64), intent(in) :: dtau,rho,fi
    logical, intent(out) :: converge

    logical, dimension(3), parameter :: mask = [ .true., .false., .false. ]
    real(REAL64), dimension(3) :: psi
    real(REAL64) :: a,b,fa,fb,f,d,tol
    integer :: iter, nmax, n

    converge = .false.

    n = size(tt)
    tol = max(10*epsilon(tau)*(1+tau),0.1*dtau/sqrt(real(n)))
    nmax = nint(abs(log(0.1*tol/tau))/log(2.0)) + 1
    if( debug_iter ) &
         write(error_unit,'(a,es12.2,3x,i0)') 'bisec tol: ',tol,nmax

    a = taumin
    b = taumax
    psi = psisum(tt,mju,a,rho,fi,tmin,tmax,kmean,mask)
    fa = psi(1)
    psi = psisum(tt,mju,b,rho,fi,tmin,tmax,kmean,mask)
    fb = psi(1)
    do iter = 1, nmax

       tau = (a + b) / 2
       psi = psisum(tt,mju,tau,rho,fi,tmin,tmax,kmean,mask)
       f = psi(1)

       if( debug_iter ) write(error_unit,'(a,i2,2x,3(g0:3x))') &
            "bi-limits, a, b, tau: ",iter,real(a),real(b), real(tau)

       if( fa*f < 0 ) then
          b = tau
          fb = f
       else if( f*fb < 0 ) then
          a = tau
          fa = f
       else ! if( f == 0 .or. fa == 0 .or fb == 0 ) then
          !          continue
          if( debug ) write(error_unit,*) 'Out of a root interval. '
          exit
       end if

       converge = abs(b - a) < tol .and. abs(taumin - tau) > tol &
            .and. abs(taumax - tau) > tol

       if( converge ) exit
    end do

    if( .not. converge .and. debug_graph ) &
         call graph_like_t('/tmp/bi',tau,dtau,rho,fi)

  end subroutine bisection


  ! -------------------------------------------------------------------


  subroutine graph_like_tr(filename,tau,dtau,rho,phi)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: tau,dtau,rho,phi
    real(REAL64), dimension(3) :: psi
    real(REAL64), dimension(3,3) :: dpsi
    logical, dimension(3) :: mask = [ .true., .true., .true. ]
    real(REAL64) :: drho, t,r,f
    integer :: i,j,l

!    dtau = 0.02
    drho = 10

!    dtau = (taumax - taumin) / 150
!    drho = 2*rho / 200

    open(1,file=filename)
    do i = -60,60
       t = tau + 0.05*i*dtau
       do j = -50, 50
          r = rho + j*drho
          f = exp(elogf(t,r,phi))
          psi = psisum(tt,mju,t,r,phi,tmin,tmax,kmean,mask)
          dpsi = dpsisum(tt,mju,t,r,phi,tmin,tmax,kmean,mask)
          write(1,*) t,r,f,psi,(dpsi(l,l),l=1,3)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_like_tr


  subroutine graph_like_tf(filename,tau,dtau,rho)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: tau,dtau,rho
    real(REAL64), dimension(3) :: psi
    real(REAL64), dimension(3,3) :: dpsi
    logical, dimension(3) :: mask = [ .true., .true., .true. ]
    real(REAL64) :: t,fi,f
    integer :: i,j,l

!    dtau = 0.02

    open(1,file=filename)
    do i = -50,50
       t = tau + 0.1*i*dtau
       do j = -180, 180, 10
          fi = j / RAD
          f = exp(elogf(t,rho,fi))
          psi = psisum(tt,mju,t,rho,fi,tmin,tmax,kmean,mask)
          dpsi = dpsisum(tt,mju,t,rho,fi,tmin,tmax,kmean,mask)
          write(1,*) t,RAD*fi,f,psi,(dpsi(l,l),l=1,3)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_like_tf

  subroutine graph_like_t(filename,tau,dtau,rho,phi)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: tau, dtau,rho, phi

    real(REAL64), dimension(3) :: psi
    real(REAL64), dimension(3,3) :: dpsi
    logical, dimension(3) :: mask = [ .true., .true., .true. ]
    integer :: i
    real(REAL64) :: t

    open(1,file=filename)
    do i = -50, 50
       t = tau + 0.1*i*dtau
       psi = psisum(tt,mju,t,rho,phi,tmin,tmax,kmean,mask)
       dpsi = dpsisum(tt,mju,t,rho,phi,tmin,tmax,kmean,mask)
       write(1,*) t,psi(1),dpsi(1,1)
    end do
    close(1)

  end subroutine graph_like_t


!!$  subroutine graph_maxlike_tau(tau,rho,phi)
!!$
!!$    real(REAL64), intent(in) :: tau,rho,phi
!!$    real(REAL64), dimension(3) :: x,fvec,fvec1,fvec2
!!$    real(REAL64), dimension(3,3) :: fjac
!!$    real(REAL64) :: dtau, t,f1,h
!!$    integer :: i,iflag
!!$
!!$    dtau = 0.01
!!$    h = epsilon(tau)**(1.0/3.0)
!!$
!!$    open(1,file='/tmp/maxliketau')
!!$    do i = -100,100
!!$       t = tau + i*dtau
!!$       x = [ t, rho, phi ]
!!$       iflag = 1
!!$       call maxlike(3,3,x,fvec,fjac,3,iflag)
!!$       x = [ t + h, rho, phi ]
!!$       call maxlike(3,3,x,fvec1,fjac,3,iflag)
!!$       x = [ t - h, rho, phi ]
!!$       call maxlike(3,3,x,fvec2,fjac,3,iflag)
!!$       f1 = (fvec1(1) - fvec2(1))/(2*h)
!!$       iflag = 2
!!$       x = [ t, rho, phi ]
!!$       call maxlike(3,3,x,fvec,fjac,3,iflag)
!!$       write(1,*) t,fvec(1),fjac(1,1),f1
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graph_maxlike_tau

  subroutine graph_psi_ttau(filename,mu,tau0,rho,phi)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: mu,tau0,rho,phi
    real(REAL64), dimension(3) :: psi,psi1,psi2
    real(REAL64), dimension(3,3) :: dpsi
    logical, dimension(3), parameter :: mask = [ .true., .true., .true. ]
    real(REAL64) :: dtau, tau, t,psid,h
    integer :: i,j

    h = epsilon(tau)**(1.0/3.0)
    dtau = (taumax - taumin) / 150

    open(1,file=filename)
    do i = -50,150,5
       tau = tau0 + i*dtau
       do j = -50, 150,5
          t = tau0 + j*dtau

          psi = psisum([t],[mu],tau,rho,phi,tmin,tmax,[mu],mask)
          dpsi = dpsisum([t],[mu],tau,rho,phi,tmin,tmax,[mu],mask)
          psi1 = psisum([t],[mu],tau+h,rho,phi,tmin,tmax,[mu],mask)
          psi2 = psisum([t],[mu],tau-h,rho,phi,tmin,tmax,[mu],mask)
          psid = (psi1(1) - psi2(1))/(2*h)

          write(1,*) tau,t,psi(1),dpsi(1,1),psid
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_psi_ttau

  subroutine graph_psi_trho(filename,mu,tau,rho0,phi)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: mu,tau,rho0,phi
    real(REAL64), dimension(3) :: psi,psi1,psi2
    real(REAL64), dimension(3,3) :: dpsi
    logical, dimension(3), parameter :: mask = [ .true., .true., .true. ]
    real(REAL64) :: drho, rho, t,psid,h, dtau
    integer :: i,j

    dtau = (taumax - taumin) / 150
    drho = 3*rho0 / 55
    h = epsilon(tau)**(1.0/3.0)

    open(1,file=filename)
    do i = -10,50
       rho = rho0 + i*drho
       do j = -50, 150
          t = tau + j*dtau

          psi = psisum([t],[mu],tau,rho,phi,tmin,tmax,[mu],mask)
          dpsi = dpsisum([t],[mu],tau,rho,phi,tmin,tmax,[mu],mask)
          psi1 = psisum([t],[mu],tau,rho+h,phi,tmin,tmax,[mu],mask)
          psi2 = psisum([t],[mu],tau,rho-h,phi,tmin,tmax,[mu],mask)
          psid = (psi1(2) - psi2(2))/(2*h)

          write(1,*) rho,t,psi(2),dpsi(2,2),psid
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_psi_trho

  subroutine graph_psi_tphi(filename,mu,tau,rho,phi0)

    character(len=*), intent(in) :: filename
    real(REAL64), intent(in) :: mu,tau,rho,phi0
    real(REAL64), dimension(3) :: psi,psi1,psi2
    real(REAL64), dimension(3,3) :: dpsi
    logical, dimension(3), parameter :: mask = [ .true., .true., .true. ]
    real(REAL64) :: phi,t,psid,h,dtau
    integer :: i,j

    dtau = (taumax - taumin) / 150
    h = epsilon(tau)**(1.0/3.0)

    open(1,file=filename)
    do i = 0, 90
       phi = i / RAD
       do j = -50, 150
          t = tau + j*dtau

          psi = psisum([t],[mu],tau,rho,phi,tmin,tmax,[mu],mask)
          dpsi = dpsisum([t],[mu],tau,rho,phi,tmin,tmax,[mu],mask)
          psi1 = psisum([t],[mu],tau,rho,phi+h,tmin,tmax,[mu],mask)
          psi2 = psisum([t],[mu],tau,rho,phi-h,tmin,tmax,[mu],mask)
          psid = (psi1(3) - psi2(3))/(2*h)

          write(1,*) RAD*phi,t,psi(3),dpsi(3,3),psid
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_psi_tphi


  ! --------------------------------------------------------------------

!!$  subroutine thelike1(t,dt,sigu,dsigu,sigv,dsigv,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, sigu, sigv
!!$    real(REAL64), intent(out) :: dt, dsigu, dsigv
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 3
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$!    real(REAL64), parameter :: eps2 = sqrt(epsilon(t))
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps, s0, s
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$    diagonal = .false.
!!$    do i = 1, npar
!!$       diagonal(i,i) = .true.
!!$    end do
!!$
!!$    reli = .false.
!!$    nprint = 1
!!$    n = size(res)
!!$    allocate(dpsi(n),stat=stat,errmsg=errmsg)
!!$    ! Computations needs some additional memory,
!!$    ! a fail-back way should be used if it is unavailable.
!!$    if( stat /= 0 ) then
!!$       if( debug ) write(error_unit,*) warning_label,trim(errmsg),' (thelike)'
!!$!       reli = .false.
!!$       return
!!$    end if
!!$
!!$    eps = 10*n*epsilon(t)*(1 + abs(t))
!!$    eps = 0.01*sqrt(1/cts(1))
!!$!    eps = 10*n*sqrt(epsilon(t))*(1 + abs(t))
!!$!    rmax = maxval(abs(cts - t*pht))
!!$!    if( abs(r) > 10*n*sqrt(epsilon(r)) ) then
!!$!       rtol = 1e-5*r
!!$!    else
!!$!       rtol = 1e-5*rmax
!!$
!!$    !    end if
!!$!    ttol = max(1e-3*dt,10*n*sqrt(epsilon(t)))
!!$!    rtol = 0.01*dr
!!$!    stol = 0.5e-4
!!$!    ztol = ttol
!!$
!!$    dt = -1
!!$    dsigu = -1
!!$    dsigv = -1
!!$!    s0 = -1
!!$
!!$!    write(*,*) ttol, rtol, stol, ztol, (r / minval(cts))**2
!!$
!!$!    z = 0
!!$!    dz = (r / minval(cts))**2
!!$!    dz = ztol
!!$!    call reziduals(t,r,z,res,sig)
!!$!    c = s - stukey
!!$!    c = s - 1
!!$!    if( abs(c) > 0 ) then
!!$!       mu = sum(itukey(res/(s*sig))) / (10*c**2)
!!$!    else
!!$!       mu = sum(itukey(res/(s*sig))) / 10
!!$!    end if
!!$!    mu = sum(itukey(res/(s*sig))) / 10
!!$!    if( mu < 0.1 ) mu = 1
!!$!    lambda = mu*c
!!$
!!$!    write(*,*) mu,lambda,c
!!$!    stop 0
!!$
!!$    if( .true. ) then
!!$       block
!!$!         real(REAL64) :: sigv, sigu
!!$!         t = 2
!!$!         sigu = 2e3
!!$!         sigv = 2e3
!!$         call testmaxlike(t,sigu,sigv)
!!$!         stop 0
!!$       end block
!!$    end if
!!$
!!$    x = [ t, sigu, sigv ]
!!$    !call lmdif3(maxliked,x,cov,eps,666,100d0,nprint=nprint,info=info)
!!$    call lmder3(maxlike,x,cov,eps,666,1d0,nprint,info)
!!$
!!$!    write(*,*) x
!!$    t = x(1)
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$    reli = .true.
!!$
!!$!    write(*,*) 'diag:',pack(cov,diagonal)
!!$
!!$    if( any(x < 0) .or. any(pack(cov,diagonal) < 0) .or. &
!!$         info == 666 .or. info < 0 ) then
!!$       if( debug_iter ) write(error_unit,*) &
!!$            'lmder3 failed at ',real(x),'for',info,' or reached a maximum.'
!!$
!!$!       x = [ t, r, z, s ]
!!$!       dx = [ 0.1*dt, 0.1*dr, 0.01*dz, 1d-2 ]
!!$!       call nelmin1(minpeak,x,dx,f,info)
!!$!       kind = 'L'
!!$!       if( debug_iter ) write(error_unit,*) 'nelmin: ',real(x),info
!!$    end if
!!$
!!$    s = 1
!!$    if( all(pack(cov,diagonal) < 0) ) then
!!$       sig = sqrt(cts + sigu**2 + t**2*(pht + sigv**2))
!!$       res = (cts - t*pht)
!!$       s0 = s**2*sum(tukey(res/(s*sig))**2)/ &
!!$            sum(dtukey(res/(s*sig)))*real(n)**2/real(n-1)
!!$       s0 = sum((res/sig)**2)
!!$!          write(*,*) s0,s,sum(tukey(res/sig/s)**2),sum(dtukey(res/sig/s))
!!$       if( s0 > 0 ) then
!!$          dt = sqrt(s0*abs(cov(1,1))/(n-3))
!!$          dsigu = sqrt(s0*abs(cov(2,2))/(n-3))
!!$          dsigv = sqrt(s0*abs(cov(3,3))/(n-3))
!!$!          ttol = max(0.001*dt,ttol)
!!$!          rtol = max(0.01*dr,rtol)
!!$!          ztol = max(0.01*dz,ztol)
!!$!             stol = max(0.01*ds,stol)
!!$          s0 = sqrt(s0 / n)
!!$       end if
!!$    else
!!$       s0 = -1
!!$       dt = -1
!!$       dsigu = -1
!!$       dsigv = -1
!!$    end if




!!$    converge = .false.
!!$    do iter = 1, 3*precision(r)
!!$
!!$       x = [ t, r, z, s ]
!!$       !       call lmder3(peakder,x,cov,eps,666,1d0,nprint,info)
!!$       info = -1
!!$       kind = 'D'
!!$       call lmdif2(maxliked,x,eps,nprint,info)
!!$       cov = -1
!!$!       write(*,*) pack(cov,diagonal)
!!$       if( any(x < 0) .or. any( pack(cov,diagonal) < 0) .or. &
!!$            info == 666 .or. info < 0 ) then
!!$          if( debug_iter ) write(error_unit,*) &
!!$               'lmder3 failed at ',real(x),'for',info,' or reached a maximum.'
!!$
!!$          x = [ t, r, z, s ]
!!$          dx = [ 0.1*dt, 0.1*dr, 0.01*dz, 1d-2 ]
!!$          call nelmin1(minpeak,x,dx,f,info)
!!$          kind = 'L'
!!$!          if( debug_iter ) write(error_unit,*) 'nelmin: ',real(x),info
!!$       end if
!!$
!!$       if( debug_iter ) then
!!$          c = s - stukey
!!$          if( s > 0 ) then
!!$             call reziduals(x(1),x(2),x(3),res,sig)
!!$             f = sum(itukey(res/(s*sig)))+sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$          else
!!$             f = 0
!!$          end if
!!$          write(error_unit, &
!!$               '(3a,i0,2x,g0,2x,f0.1,2x,g0,2x,f6.4,1x,f10.5,2x,2(g0,2x),i0)') &
!!$               'peak(',kind,'): #',iter,real(x),lambda,real(mu),real(f),info
!!$       end if
!!$!       stop 0
!!$
!!$       if( (iter > 1 .and. (x(2) > rmax .or. x(4) < 0.1 .or. x(4) > 2 ))  &
!!$            .or. info < 0 ) then
!!$          converge = .false.
!!$          exit
!!$       end if
!!$
!!$!       write(*,*) x(1) - t, x(2) - r, x(3) - stukey
!!$       converge = abs(x(1) - t) < ttol  .and. abs(x(2) - r) < rtol .and. &
!!$            abs(x(3) - z) < ztol .and. abs(x(4) - stukey) < stol
!!$
!!$       t = x(1)
!!$       r = x(2)
!!$       z = x(3)
!!$       s = x(4)
!!$
!!$       call reziduals(t,r,z,res,sig)
!!$       if( s > 0 .and. all(pack(cov,diagonal) > 0) ) then
!!$          s0 = s**2*sum(tukey(res/(s*sig))**2)/ &
!!$               sum(dtukey(res/(s*sig)))*real(n)**2/real(n-1)
!!$!          write(*,*) s0,s,sum(tukey(res/sig/s)**2),sum(dtukey(res/sig/s))
!!$          if( s0 > 0 ) then
!!$             dt = sqrt(s0*cov(1,1)/(n-3))
!!$             dr = sqrt(s0*cov(2,2)/(n-3))
!!$             ds = sqrt(s0*cov(3,3)/(n-3))
!!$             ttol = max(0.001*dt,ttol)
!!$             rtol = max(0.01*dr,rtol)
!!$             ztol = max(0.01*dz,ztol)
!!$!             stol = max(0.01*ds,stol)
!!$             s0 = sqrt(s0 / n)
!!$          end if
!!$       else
!!$          s0 = -1
!!$          dt = -1
!!$          dr = -1
!!$          dz = -1
!!$          ds = -1
!!$       end if
!!$
!!$       if( converge ) exit
!!$
!!$       lambda = lambda - mu*(s - stukey)
!!$       mu = 2*mu
!!$
!!$    end do
!!$
!!$!    if( info == 0 ) then
!!$!       t = x(1)
!!$!       r = x(2)
!!$!       s = x(3)
!!$!    end if
!!$
!!$!    block
!!$!      call testders(t,s,r)
!!$!      stop 0
!!$!    end block
!!$
!!$
!!$
!!$    reli = converge
!!$!    r = abs(r)
!!$    deallocate(dpsi)
!!$

!!$    if( reli .and. debug_iter ) then
!!$       block
!!$         real(REAL64), dimension(npar,npar) :: cor
!!$         real(REAL64) :: q
!!$         integer :: i,j
!!$
!!$         if( all(pack(cov,diagonal) > 0) ) then
!!$            cor = 1
!!$            do i = 1, npar
!!$               do j = i+1, npar
!!$                  q = sqrt(cov(i,i)*cov(j,j))
!!$                  cor(i,j) = cov(i,j) / q
!!$                  cor(j,i) = cor(i,j)
!!$               end do
!!$            end do
!!$         else
!!$            cor = -1
!!$         end if
!!$         write(error_unit,'(a,3(es12.4,1x),t48,a,3(f7.4,1x))') &
!!$              'Cov: ', real(cov(1,:)),'corr: ', cor(1,:)
!!$         do i = 2,npar
!!$            write(error_unit,'(5x,3(es12.4,1x),t48,6x,3(f7.4,1x))') &
!!$                 real(cov(i,:)),cor(i,:)
!!$         end do
!!$       end block
!!$    end if
!!$
!!$  end subroutine thelike1
!!$
!!$  subroutine wmaxlike(k,l,x,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k, l, ldfjac
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$    real(REAL64), dimension(ldfjac,k), intent(out) :: fjac
!!$    real(REAL64), dimension(:), allocatable :: sig, psi, dpsi, oak, res, sigx2,&
!!$         wcts, wpht
!!$    logical, dimension(:), allocatable :: mask
!!$    real(REAL64) :: t, sigu, sigv, sumec, sig1
!!$    integer :: n, i, j
!!$
!!$    n = size(cts)
!!$    t = x(1)
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$
!!$    if( .not. (t > 0 .and. sigu < 1e10 .and. sigv < 1e10) ) then
!!$       iflag = -1
!!$       if( debug_iter ) &
!!$            write(error_unit,*) 'maxlike: assumptions violated.'
!!$       return
!!$    end if
!!$
!!$
!!$    allocate(sig(n),psi(n),dpsi(n),oak(n),res(n),sigx2(n),wcts(n),wpht(n))
!!$
!!$    if( iflag == 0 .and. debug_iter ) then
!!$       write(error_unit,*) 'x = ',real(x)
!!$       write(error_unit,*) 'fvec = ',real(fvec)
!!$       write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
!!$       write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
!!$       write(error_unit,*) 'fjac(3,:) = ',real(fjac(3,:))
!!$       return
!!$    end if
!!$
!!$    wcts = cts
!!$    wpht = pht
!!$    sigx2 = wpht + sigv**2
!!$    sig = sqrt(wcts + sigu**2 + t**2*sigx2)
!!$    res = wcts - t*wpht
!!$
!!$    mask = res/sig > -0.87 .and. res/sig < 1
!!$!    write(*,*) count(mask),n
!!$    mask = .true.
!!$
!!$    sigx2 = wpht + sigv**2
!!$    sig = sqrt(wcts + sigu**2 + t**2*sigx2)
!!$    sig1 = sum(1/sig**2,1,mask)
!!$    res = wcts - t*wpht
!!$    psi = huber(res/sig)
!!$    psi = res/sig
!!$    oak = wcts*(wpht*(1+t) + t*sigv**2) + wpht*sigu**2
!!$
!!$    if( iflag == 1 ) then
!!$       sumec = sum(psi*res/sig**3,1,mask)
!!$       fvec(1) = sum(psi*wpht/sig,1,mask) + t*sum(psi*res*sigx2/sig**3,1,mask) + &
!!$            sum(sigx2*wcts/oak,1,mask) - 3*t*sum(sigx2/sig**2,1,mask)
!!$       fvec(2) = (sumec   + 2*sum(wpht/oak,1,mask) - 3*sig1)*sigu
!!$       fvec(3) = (sumec*t + 2*sum(wcts/oak,1,mask) - 3*t*sig1)*t*sigv
!!$       write(*,*) real(fvec)
!!$    end if
!!$
!!$    if( iflag == 2 ) then
!!$
!!$       dpsi = dhuber(res/sig)
!!$       dpsi = 1
!!$       fjac(1,1) = -sum(dpsi*((wpht + t*res*sigx2/sig**2)/sig)**2,1,mask) - &
!!$            t*sum(psi*wpht*sigx2/sig**3,1,mask) + sum(psi*sigx2*(wcts-2*t*wpht)/sig**3,1,mask)+&
!!$            3*t*sum(psi*res*sigx2/sig**5,1,mask) - sum((wcts*sigx2/oak)**2,1,mask) + &
!!$            6*t*sum((sigx2/sig**2)**2,1,mask)
!!$       fjac(1,2) = -sigu*sum(dpsi*(wpht + t*res*(sigx2/sig**2))*res/sig**4,1,mask) - &
!!$            sigu*sum(psi*wpht/sig**3,1,mask) - 3*t*sigu*sum(psi*res*sigx2/sig**5,1,mask) - &
!!$            2*sigu*sum(wcts*wpht*sigx2/oak**2,1,mask) + 6*t*sigu*sum(sigx2/sig**4,1,mask)
!!$       fjac(1,3) = -sigv*t**2*sum(dpsi*(wpht+t*res*(sigx2/sig**2))*res/sig**4,1,mask)-&
!!$            t**2*sigv*sum(psi*wpht/sig**3,1,mask) + 2*t*sigv*sum(psi*res/sig**3,1,mask) - &
!!$            3*t**3*sigv*sum(psi*res*sigx2/sig**5,1,mask) + 2*sigv*sum(wcts/oak,1,mask) - &
!!$            2*t*sigv*sum(sigx2*(wcts/oak)**2,1,mask) - 6*t*sigv*sum(1/sig**2,1,mask) + &
!!$            6*t**3*sigv*sum(sigx2/sig**4,1,mask)
!!$       fjac(2,2) = -sigu**2*sum(dpsi*(res/sig**3)**2,1,mask) + &
!!$            sum(psi*res/sig**3,1,mask) - 3*sigu**2*sum(psi*res/sig**5,1,mask) + &
!!$            2*sum(wpht/oak,1,mask) - 4*sigu**2*sum((wpht/oak)**2,1,mask) - &
!!$            3*sum(1/sig**2,1,mask) + 6*sigu**2*sum(1/sig**4,1,mask)
!!$       fjac(2,3) = -t**2*sigu*sigv*sum(dpsi*(res/sig**3)**2,1,mask) - &
!!$            3*t**2*sigu*sigv*sum(psi*res/sig**5,1,mask) - &
!!$            4*t*sigu*sigv*sum(wpht*wcts/oak**2,1,mask) + 6*t**2*sigu*sigv*sum(1/sig**4,1,mask)
!!$       fjac(3,3) = -t**4*sigv**2*sum(dpsi*(res/sig**3)**2,1,mask) + &
!!$            t**2*sum(psi*res/sig**3,1,mask) - 3*t**4*sigv**2*sum(psi*res/sig**5,1,mask) + &
!!$            2*t*sum(wcts/oak,1,mask) - 4*t**2*sigv**2*sum((wcts/oak)**2,1,mask) - &
!!$            3*t**2*sum(1/sig**2,1,mask) + 6*t**4*sigv**2*sum(1/sig**4,1,mask)
!!$
!!$       do i = 2,3
!!$          do j = 1, i-1
!!$             fjac(i,j) = fjac(j,i)
!!$          end do
!!$       end do
!!$    end if
!!$
!!$  end subroutine wmaxlike
!!$
!!$
!!$  subroutine maxlike1(k,l,x,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k, l, ldfjac
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$    real(REAL64), dimension(ldfjac,k), intent(out) :: fjac
!!$    real(REAL64), dimension(:), allocatable :: sig, psi, dpsi, oak, res, sigx2
!!$    real(REAL64) :: t, sigu, sigv, sumec, sig1
!!$    integer :: n, i, j
!!$
!!$    n = size(cts)
!!$    t = x(1)
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$
!!$    if( .not. (t > 0 .and. sigu < 1e10 .and. sigv < 1e10) ) then
!!$       iflag = -1
!!$       if( debug_iter ) &
!!$            write(error_unit,*) 'maxlike: assumptions violated.'
!!$       return
!!$    end if
!!$
!!$
!!$    allocate(sig(n),psi(n),dpsi(n),oak(n),res(n),sigx2(n))
!!$
!!$    if( iflag == 0 .and. debug_iter ) then
!!$       write(error_unit,*) 'x = ',real(x)
!!$       write(error_unit,*) 'fvec = ',real(fvec)
!!$       write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
!!$       write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
!!$       write(error_unit,*) 'fjac(3,:) = ',real(fjac(3,:))
!!$       return
!!$    end if
!!$
!!$    sigx2 = pht + sigv**2
!!$    sig = sqrt(cts + sigu**2 + t**2*sigx2)
!!$    sig1 = sum(1/sig**2)
!!$    res = cts - t*pht
!!$    !psi = huber(res/sig)
!!$    psi = res/sig
!!$    oak = cts*(pht*(1+t) + t*sigv**2) + pht*sigu**2
!!$
!!$    if( iflag == 1 ) then
!!$       sumec = sum(psi*res/sig**3)
!!$       fvec(1) = sum(psi*pht/sig) + t*sum(psi*res*sigx2/sig**3) + &
!!$            sum(sigx2*cts/oak) - 3*t*sum(sigx2/sig**2)
!!$       fvec(2) = (sumec   + 2*sum(pht/oak) - 3*sig1)*sigu
!!$       fvec(3) = (sumec*t + 2*sum(cts/oak) - 3*t*sig1)*t*sigv
!!$       write(*,*) real(fvec)
!!$    end if
!!$
!!$    if( iflag == 2 ) then
!!$
!!$       !dpsi = dhuber(res/sig)
!!$       dpsi = 1
!!$       fjac(1,1) = -sum(dpsi*((pht + t*res*sigx2/sig**2)/sig)**2) - &
!!$            t*sum(psi*pht*sigx2/sig**3) + sum(psi*sigx2*(cts-2*t*pht)/sig**3)+&
!!$            3*t*sum(psi*res*sigx2/sig**5) - sum((cts*sigx2/oak)**2) + &
!!$            6*t*sum((sigx2/sig**2)**2)
!!$       fjac(1,2) = -sigu*sum(dpsi*(pht + t*res*(sigx2/sig**2))*res/sig**4) - &
!!$            sigu*sum(psi*pht/sig**3) - 3*t*sigu*sum(psi*res*sigx2/sig**5) - &
!!$            2*sigu*sum(cts*pht*sigx2/oak**2) + 6*t*sigu*sum(sigx2/sig**4)
!!$       fjac(1,3) = -sigv*t**2*sum(dpsi*(pht+t*res*(sigx2/sig**2))*res/sig**4)-&
!!$            t**2*sigv*sum(psi*pht/sig**3) + 2*t*sigv*sum(psi*res/sig**3) - &
!!$            3*t**3*sigv*sum(psi*res*sigx2/sig**5) + 2*sigv*sum(cts/oak) - &
!!$            2*t*sigv*sum(sigx2*(cts/oak)**2) - 6*t*sigv*sum(1/sig**2) + &
!!$            6*t**3*sigv*sum(sigx2/sig**4)
!!$       fjac(2,2) = -sigu**2*sum(dpsi*(res/sig**3)**2) + &
!!$            sum(psi*res/sig**3) - 3*sigu**2*sum(psi*res/sig**5) + &
!!$            2*sum(pht/oak) - 4*sigu**2*sum((pht/oak)**2) - &
!!$            3*sum(1/sig**2) + 6*sigu**2*sum(1/sig**4)
!!$       fjac(2,3) = -t**2*sigu*sigv*sum(dpsi*(res/sig**3)**2) - &
!!$            3*t**2*sigu*sigv*sum(psi*res/sig**5) - &
!!$            4*t*sigu*sigv*sum(pht*cts/oak**2) + 6*t**2*sigu*sigv*sum(1/sig**4)
!!$       fjac(3,3) = -t**4*sigv**2*sum(dpsi*(res/sig**3)**2) + &
!!$            t**2*sum(psi*res/sig**3) - 3*t**4*sigv**2*sum(psi*res/sig**5) + &
!!$            2*t*sum(cts/oak) - 4*t**2*sigv**2*sum((cts/oak)**2) - &
!!$            3*t**2*sum(1/sig**2) + 6*t**4*sigv**2*sum(1/sig**4)
!!$
!!$       do i = 2,3
!!$          do j = 1, i-1
!!$             fjac(i,j) = fjac(j,i)
!!$          end do
!!$       end do
!!$    end if
!!$
!!$  end subroutine maxlike1
!!$
!!$  subroutine testmaxlike1(t,sigu,sigv)
!!$
!!$    real(REAL64), intent(in) :: t,sigu,sigv
!!$    real(REAL64), dimension(3) :: fd,fvec,eps,x,x1,x2
!!$    real(REAL64), dimension(3,3) :: fjac, djac
!!$    real(REAL64) :: f0
!!$    integer :: n, i, j, iflag
!!$
!!$    n = size(cts)
!!$    eps = 1e-7
!!$
!!$    iflag = 1
!!$    x = [ t, sigu, sigv ]
!!$    call maxlike(3,3,x,fvec,fjac,3,iflag)
!!$    f0 = -n*ropt7min(x)
!!$    do i = 1, 3
!!$       x = [ t, sigu, sigv ]
!!$       x(i) = x(i) + eps(i)
!!$       fd(i) = (-n*ropt7min(x) - f0) / eps(i)
!!$       write(error_unit,*) i,real(fvec(i)),real(fd(i)),real((fd(i)-fvec(i)))
!!$    end do
!!$
!!$!    return
!!$
!!$    iflag = 2
!!$    x = [ t, sigu, sigv ]
!!$    call maxlike(3,3,x,fvec,fjac,3,iflag)
!!$    call fdjac1(testmaxlike_aux,3,x,fvec,djac,3,iflag,3,3,epsilon(x),x1,x2)
!!$    do i = 1,3
!!$       do j = i, 3
!!$          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
!!$               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))!/fjac(i,j))
!!$       end do
!!$    end do
!!$
!!$
!!$  end subroutine testmaxlike1
!!$
!!$  subroutine testmaxlike_aux1(k,x,fvec,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(k), intent(out) :: fvec
!!$
!!$    real(REAL64), dimension(3,3) :: fjac
!!$    integer :: nflag
!!$
!!$    nflag = 1
!!$    call maxlike(3,3,x,fvec,fjac,3,nflag)
!!$
!!$  end subroutine testmaxlike_aux1



!!$  subroutine maxliked1(k,l,x,fvec,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k, l
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$!    real(REAL64), dimension(:), allocatable :: tt, mu, sig, psi, oak
!!$    real(REAL64), dimension(:), allocatable :: sig, psi, oak, arc, res, sigx2
!!$    real(REAL64) :: t, sigu, sigv, sumec, sig1
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    t = x(1)
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$
!!$!    allocate(tt(n),sig(n),mu(n),psi(n),oak(n))
!!$    allocate(sig(n),psi(n),oak(n),arc(n),res(n),sigx2(n))
!!$!    tt = cts / pht
!!$!    mu = (pht + cts / t) / 2
!!$!    sig = sqrt(mu*(tt**2 + t) + tt**2*sigv**2 + sigu**2)
!!$!    ee = sum((mu*(tt - t)/sig)**2)/2
!!$!    psi = mu*(t - tt)/sig
!!$    !    oak = (t - tt)*tt*(mu + sigv**2)/sig**3 + 1/sig
!!$    sigx2 = pht + sigv**2
!!$    sig = sqrt(cts + sigu**2 + t**2*sigx2)
!!$    sig1 = sum(1/sig**2)
!!$!    arc = cts*(pht*(1+t) + t*sigv**2) + pht*sigu**2
!!$!    oak = arc / sig**3
!!$    res = cts - t*pht
!!$    psi = res/sig
!!$
!!$!    fvec(1) = sum(psi*(pht/sig + t*res*(pht + sigv**2)/sig**3)) + &
!!$!         sum((pht + sigv**2)*(cts - 3*t*arc/sig**2)/sig**3/oak)
!!$!    fvec(2) = sum(psi*res/sig**3)*sigu + &
!!$!         sum((2*pht-3*arc/sig**2)*sigu/sig**3/oak)
!!$!    fvec(3) = sum(psi*res/sig**3)*sigv*t**2 + &
!!$!         sum((2*cts-3*arc*t/sig**2)*t*sigv/sig**3/oak)
!!$
!!$    oak = cts*(pht*(1+t) + t*sigv**2) + pht*sigu**2
!!$    sumec = sum(psi*res/sig**3)
!!$
!!$    fvec(1) = sum(psi*(pht + t*res*(sigx2/sig**2))/sig) + &
!!$         sum(sigx2*cts/oak) - 3*t*sum(sigx2/sig**2)
!!$    fvec(2) = (sumec   + 2*sum(pht/oak) - 3*sig1)*sigu
!!$    fvec(3) = (sumec*t + 2*sum(cts/oak) - 3*t*sig1)*t*sigv
!!$
!!$
!!$!    fvec = fvec / n
!!$
!!$!    fvec(1) = -sum(psi*(mu/sig - (t-tt)*mu**2/sig**3/2)) + &
!!$!         sum((tt*(mu + sigv**2)/sig**3 - 1.5*(t-tt)*tt*(mu+sigv**2)/sig**5 - &
!!$!         mu/(2*sig))/oak)
!!$!    fvec(2) = sum(psi*mu*(t-tt)*sigu/sig**3) - &
!!$!         sum((3*(t - tt)*tt*(mu + sigv**2)/sig**2 + 1)*sigu/sig**3/oak)
!!$!    fvec(3) = sum(psi*mu*(t-tt)*tt**2*sigv/sig**3) - &
!!$!         sum((((t-tt)*tt*(2-3*tt**2*(mu+sigv**2)/sig**2))-tt**2)*sigv/sig**3/oak)
!!$
!!$    write(error_unit,*) 'maxliked:',real(x),real(fvec)
!!$
!!$  end subroutine maxliked1
!!$
!!$  subroutine testmaxliked(t,sigu,sigv)
!!$
!!$    real(REAL64), intent(in) :: t,sigu,sigv
!!$    real(REAL64), dimension(3) :: fd,fvec,eps,x,x1,x2
!!$    real(REAL64), dimension(3,3) :: fjac, djac
!!$!    real(REAL64), dimension(:), allocatable :: mu, tt
!!$    real(REAL64) :: f0
!!$    integer :: n, i, j, iflag
!!$
!!$!    n = size(cts)
!!$!    allocate(mu(n),tt(n))
!!$!    tt = cts / pht
!!$!    mu = (pht + cts / t) / 2
!!$
!!$    eps = 1e-7
!!$    iflag = 1
!!$    x = [ t, sigu, sigv ]
!!$    call maxliked(3,3,x,fvec,iflag)
!!$    fvec = fvec / size(cts)
!!$    f0 = -ropt7min(x)
!!$    do i = 1, 3
!!$       x = [ t, sigu, sigv ]
!!$       x(i) = x(i) + eps(i)
!!$       fd(i) = (-ropt7min(x) - f0) / eps(i)
!!$       write(error_unit,*) i,real(fvec(i)),real(fd(i)),real((fd(i)-fvec(i)))
!!$    end do
!!$
!!$  end subroutine testmaxliked

!!$  subroutine ropt5(t,dt,r,sigu,sigv,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, dt, r, sigu, sigv
!!$!    real(REAL64), intent(out) :: z, dr, dz, ds, s0
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 4
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$!    real(REAL64), parameter :: eps2 = sqrt(epsilon(t))
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$!    x = [ t, pht(1), r/2, r/2 ]
!!$!    dx = [ 0.1*dt, 0.1*sqrt(pht(1)), 0.01*r, 0.01*r ]
!!$    !    dx = [ 0.001*t, 0.1*sqrt(pht(1)), 0.01*r, 0.01*r ]
!!$    !    r = sum(sqrt((cts/pht)**2*cts + pht))/size(cts)
!!$!    r = 2e2
!!$    x = [ t, sum(pht)/size(pht), r/2, r/2 ]
!!$    dx = [ 0.001*t, 0.001*sqrt(pht(1)), 0.001*r, 1d-3*r ]
!!$    call nelmin1(ropt5min,x,dx,f,info)
!!$    reli = .true.
!!$
!!$    t = x(1)
!!$    mu = x(2)
!!$    sigu = x(3)
!!$    sigv = x(4)
!!$
!!$  end subroutine ropt5
!!$
!!$  function ropt5min(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64), dimension(:), allocatable :: sig, sigx2
!!$    real(REAL64) :: t,r,z,s,c,f,mu,nu,sigu, sigv, sign, sige
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    t = x(1)
!!$    mu = x(2)
!!$    r = x(3)
!!$!    nu = x(4)
!!$    sigu = x(3)
!!$    sigv = x(4)
!!$!    sigu = x(3)*x(4)/sqrt(t**2 + x(4)**2)
!!$!    sigv = x(3)/sqrt(t**2 + x(4)**2)
!!$!    sign = x(3)
!!$!    sige = x(4)
!!$    allocate(sig(n),sigx2(n))
!!$!    tt = cts / pht
!!$!    sigu = sqrt(sign**2 - sige**2)/2
!!$!    sigv = sqrt(sign**2 + sige**2)/2
!!$!    sigv = sqrt(r**2/(tt**2 + nu**2))
!!$!    sigu = nu*sigv
!!$!    sig = sqrt(tt**2*mu + t*mu + tt**2*sigv**2 + sigu**2)
!!$    sig = sqrt(tt**2*mu + t*mu + tt**2*sigv**2 + sigu**2)
!!$    !sig = sqrt(tt**2*mu + t*mu + r**2)
!!$!    sigx2 = mu + r**2/2
!!$    sigx2 = mu + sigv**2
!!$    !sig = r
!!$
!!$!    if( .not. (0.1 < s .and. s < 2) ) then
!!$!       f = mu*c**4
!!$!       return
!!$!    end if
!!$
!!$!    call reziduals(t,r,z,res,sig)
!!$!    f = sum(itukey(res/(s*sig))) + sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$    !    f = sum((res/(s*sig))**2/2) + sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$
!!$    ! accurate with **3/2
!!$!    f = mu**2/2*sum((tt - t)**2/(tt**2*sigx**2 + sigy**2))/n - &
!!$!         sum(log(mu*(t*tt*sigx**2 + sigy**2)/(tt**2*sigx**2 + sigy**2)**1.5))/n
!!$
!!$    ! low accuracy with **1/2
!!$    f = mu**2/2*sum((tt - t)**2/sig**2)/n - &
!!$         sum(log(mu/sig + mu*sigx2*(t - tt)*tt/sig**3))/n
!!$
!!$
!!$    write(error_unit,*) 'ropt5min:',real(x),real(f)
!!$
!!$  end function ropt5min


!!$  subroutine ropt6(t,dt,r,sigu,sigv,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, dt, r, sigu, sigv
!!$!    real(REAL64), intent(out) :: z, dr, dz, ds, s0
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 3
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$!    real(REAL64), parameter :: eps2 = sqrt(epsilon(t))
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$!    x = [ t, pht(1), r/2, r/2 ]
!!$!    dx = [ 0.1*dt, 0.1*sqrt(pht(1)), 0.01*r, 0.01*r ]
!!$    !    dx = [ 0.001*t, 0.1*sqrt(pht(1)), 0.01*r, 0.01*r ]
!!$    !    r = sum(sqrt((cts/pht)**2*cts + pht))/size(cts)
!!$!    r = 2e2
!!$    x = [ t, sum(pht)/size(pht), r ]
!!$    dx = [ 0.001*t, 0.001*sqrt(pht(1)), 0.001*r ]
!!$    call nelmin1(ropt6min,x,dx,f,info)
!!$    reli = .true.
!!$
!!$    t = x(1)
!!$    mu = x(2)
!!$    sigu = 0
!!$    sigv = x(3)
!!$
!!$  end subroutine ropt6
!!$
!!$  function ropt6min(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64), dimension(:), allocatable :: tt, sig, sigx2
!!$    real(REAL64) :: t,r,z,s,c,f,mu,sigu, sigv
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    t = x(1)
!!$    mu = x(2)
!!$    sigv = x(3)
!!$    sigu = 0
!!$
!!$    allocate(tt(n),sig(n),sigx2(n))
!!$    tt = cts / pht
!!$    sig = sqrt(tt**2*mu + t*mu + tt**2*sigv**2 + sigu**2)
!!$    sigx2 = mu + sigv**2
!!$
!!$    f = mu**2/2*sum((tt - t)**2/sig**2)/n - &
!!$         sum(log(mu/sig + mu*sigx2*(t - tt)*tt/sig**3))/n
!!$
!!$    write(error_unit,*) 'ropt6min:',real(x),real(f)
!!$
!!$  end function ropt6min

!!$  subroutine ropt7(t,dt,r,sigu,sigv,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, dt, r, sigu, sigv
!!$!    real(REAL64), intent(out) :: z, dr, dz, ds, s0
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 3
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$!    real(REAL64), parameter :: eps2 = sqrt(epsilon(t))
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps,s,w,sig,y1,x1
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$!    x = [ t, pht(1), r/2, r/2 ]
!!$!    dx = [ 0.1*dt, 0.1*sqrt(pht(1)), 0.01*r, 0.01*r ]
!!$    !    dx = [ 0.001*t, 0.1*sqrt(pht(1)), 0.01*r, 0.01*r ]
!!$    !    r = sum(sqrt((cts/pht)**2*cts + pht))/size(cts)
!!$    !    r = 2e2
!!$    r = 2e3
!!$    x = [ t,  r, r ]
!!$    dx = [ 0.0001*t, 0.001*r, 0.001*r ]
!!$    call nelmin1(ropt7alt,x,dx,f,info)
!!$    reli = .true.
!!$
!!$    t = x(1)
!!$!    mu = -1
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$
!!$    s = 0
!!$    w = 0
!!$    n = size(cts)
!!$    x1 = sum(pht)/n
!!$    y1 = sum(cts)/n
!!$    do i = 1, size(cts)
!!$       t = cts(i)/pht(i)
!!$       sig = sqrt(y1 + t**2*x1 + 0e3**2 + t**2*0e3**2)
!!$       f = (pht(i)*t*(y1+0e3**2)+cts(i)*(x1 + 0e3**2))/sig**3* &
!!$            exp(-((cts(i)-t*pht(i))/sig)**2/2)
!!$       s = s + f*t
!!$       w = w + f
!!$    end do
!!$    write(*,*) 'ratio:',s/w
!!$    write(*,*) x1, y1
!!$!stop 0
!!$!    t = 1.91
!!$!    sigu = 2e3
!!$!    sigv = 2e3
!!$
!!$    x = [ 2.0, 2e3, 2e3 ]
!!$    open(1,file='/tmp/x')
!!$    do i = 50, 5000,10
!!$       x(1) = i / 1e3
!!$       write(1,*) x(1), ropt7alt(x)
!!$    end do
!!$    close(1)
!!$!    stop 0
!!$
!!$  end subroutine ropt7
!!$
!!$  function ropt7alt(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64) :: t,r,z,s,c,f,sigu, sigv, sigx, sigy, sig, y1, x1
!!$    integer :: i,n
!!$
!!$    n = size(cts)
!!$    t = x(1)
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$
!!$!    x1 = sum(pht) / n
!!$    !    y1 = sum(cts) / n
!!$    x1 = 1e4
!!$    y1 = 2e4
!!$    sigx = sqrt(x1 + sigv**2)
!!$    sigy = sqrt(y1 + sigu**2)
!!$    sig = sqrt(sigx**2*t**2 + sigy**2)
!!$
!!$!    write(*,*) x1,sigx,sig,sigv
!!$!    stop 0
!!$
!!$
!!$    f = sum(((cts - t*pht)/sig)**2)/2/n - &
!!$         sum(log(pht/sig + t*(cts - t*pht)*sigx**2/sig**3))/n
!!$!         sum(log(max(pht/sig + t*(cts - t*pht)*sigx**2/sig**3,1d0)))/n
!!$
!!$    s = 0
!!$    n = 0
!!$    do i = 1, size(cts)
!!$       if( cts(i) > 0 .and. pht(i) > 0 ) then
!!$          s = s + ((cts(i) - t*pht(i))/sig)**2/2 - &
!!$               log((cts(i)*t*sigx**2 + pht(i)*sigy**2)/sig**3)
!!$!               log(pht(i)/sig + t*(cts(i) - t*pht(i))*sigx**2/sig**3)
!!$          n = n + 1
!!$       end if
!!$    end do
!!$    f = s / n
!!$
!!$    write(error_unit,*) 'ropt7alt:',real(x),real(f)
!!$
!!$  end function ropt7alt
!!$
!!$
!!$  function ropt7min(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64), dimension(:), allocatable :: tt, sig, sigx2, mu
!!$    real(REAL64) :: t,r,z,s,c,f,sigu, sigv,sigu2,sigv2
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    t = x(1)
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$
!!$!    allocate(tt(n),sig(n),sigx2(n),mu(n))
!!$!    tt = cts / pht
!!$!    mu = (pht + cts / t) / 2
!!$!    sig = sqrt(mu*(tt**2 + t) + tt**2*sigv**2 + sigu**2)
!!$!    sigx2 = mu + sigv**2
!!$!    f = sum(mu**2*(tt - t)**2/sig**2)/n/2 - &
!!$!         sum(log(mu/sig + mu*sigx2*(t - tt)*tt/sig**3))/n
!!$
!!$    allocate(sig(n))
!!$    sigu2 = sigu**2
!!$    sigv2 = sigv**2
!!$    sig = sqrt(cts + sigu2 + t**2*(pht + sigv2))
!!$    f = sum(((cts - t*pht)/sig)**2)/2/n - &
!!$         sum(log((cts*pht*(1+t)+cts*t*sigv2 + pht*sigu2)/sig**3))/n
!!$
!!$!    write(error_unit,*) 'ropt7min:',real(x),real(f)
!!$
!!$  end function ropt7min


!!$  subroutine ropt8(t,dt,sig,sigu,sigv,mu,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, dt, sig,sigu,sigv,mu
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 3
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps, &
!!$         t0
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$    n = size(cts)
!!$
!!$    t0 = t
!!$    sigu = sig
!!$    sigv = sig
!!$    do i = 1, 3
!!$       mju = (pht + cts/t)/2
!!$       mju = sum(pht)/n
!!$!       mju = 1e4
!!$!       mju = sum((pht + cts/t))/(2*n)
!!$       x = [ t, sigu, sigv ]
!!$       dx = x / 1000.0
!!$       t0 = t
!!$       call nelmin1(ropt8min,x,dx,f,info)
!!$       t = x(1)
!!$       sigu = x(2)
!!$       sigv = x(3)
!!$       write(*,*) i, real(t),real(t0),real(t-t0)
!!$       if( abs(t - t0)/t < 1e-5 ) exit
!!$    end do
!!$
!!$    reli = .true.
!!$    write(*,*) 'ropt8:',t,sigu,sigv
!!$    mu = 1e4
!!$
!!$  end subroutine ropt8
!!$
!!$  function ropt8min(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64), dimension(:), allocatable :: tt, sigx, sigy,r
!!$    real(REAL64) :: t,f,sigu,sigv
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    allocate(tt(n),sigx(n),sigy(n),r(n))
!!$!    allocate(tt(n),sig(n),sigx(n),sigy(n))
!!$    tt = cts / pht
!!$    t = x(1)
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$
!!$    where( mju > 0 )
!!$       sigx = sqrt(mju + sigv**2)
!!$       sigy = sqrt(t*mju + sigu**2)
!!$    elsewhere
!!$       sigx = sigv
!!$       sigy = sigu
!!$    end where
!!$    sig = sqrt(tt**2*sigx**2 + sigy**2)
!!$    r = sqrt(tt**2*sigv**2 + sigu**2)
!!$    f = sum((mju*(t - tt)/sig)**2,tt>0)/2/n - &
!!$         sum(log(mju*(t*tt*sigx**2 + sigy**2)/sig**3),tt>0)/n
!!$!    f = sum((mju*(t - tt)/sig)**2,tt>0)/2/n - &
!!$!         sum(log(mju*(t*(tt+1)+r**2 + tt*(t-tt)*sigv**2)/sig**3),tt>0)/n
!!$
!!$
!!$!    write(error_unit,*) 'ropt8min:',real(x),real(f)
!!$
!!$  end function ropt8min

!!$  subroutine ropt8a(t,dt,sig,sigu,sigv,mu,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, dt, sig,sigu,sigv,mu
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 4
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps, &
!!$         t0
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$    n = size(cts)
!!$    mju = pht
!!$
!!$    t0 = t
!!$    sigu = sig
!!$    sigv = sig
!!$    mu = sum((pht + cts/t))/(2*n)
!!$    x = [ t, sigu, sigv, mu ]
!!$    dx = x / 1000.0
!!$    call nelmin1(ropt8amin,x,dx,f,info)
!!$    t = x(1)
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$    mu = x(4)
!!$
!!$    reli = .true.
!!$    write(*,*) 'ropt8a:',t,sigu,sigv,mu
!!$    write(*,*) real(dx)
!!$!    mu = 1e4
!!$
!!$  end subroutine ropt8a
!!$
!!$  function ropt8amin(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64), dimension(:), allocatable :: tt, sigx, sigy
!!$    real(REAL64) :: t,f,sigu,sigv,mu
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    allocate(tt(n),sigx(n),sigy(n))
!!$!    allocate(tt(n),sig(n),sigx(n),sigy(n))
!!$    tt = cts / pht
!!$    t = x(1)
!!$    sigu = x(2)
!!$    sigv = x(3)
!!$    mu = x(4)
!!$
!!$    sigx = sqrt(mu + sigv**2)
!!$    sigy = sqrt(t*mu + sigu**2)
!!$    sig = sqrt(tt**2*sigx**2 + sigy**2)
!!$!    n = count(tt > 0)
!!$    f = sum((mu*(t - tt)/sig)**2,tt>0)/2/n - &
!!$!    f = sum(((cts - t*pht)/sig)**2,tt>0)/2/n - &
!!$         sum(log(mu*(t*tt*sigx**2 + sigy**2)/sig**3),tt>0)/n
!!$
!!$
!!$!    write(error_unit,*) 'ropt8amin:',real(x),real(f)
!!$
!!$  end function ropt8amin

!!$  subroutine ropt88(t,dt,sig,sigu,sigv,mu,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, dt, sig,sigu,sigv,mu
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 2
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps, &
!!$         t0,r
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$    n = size(cts)
!!$    mju = pht
!!$
!!$    t0 = t
!!$    r = sig
!!$!    sigu = sig
!!$!    sigv = sig
!!$    !    mu = sum((pht + cts/t))/(2*n)
!!$    do iter = 1, 3
!!$       mju = (pht + cts/t)/2
!!$       mju = 1e4
!!$       x = [ t, r ]
!!$       dx = x / 1000.0
!!$       t0 = t
!!$       call nelmin1(ropt88min,x,dx,f,info)
!!$       t = x(1)
!!$       r = x(2)
!!$       write(*,*) iter, real(t),real(t0),real(t-t0)
!!$       if( abs(t - t0)/t < 1e-5 ) exit
!!$    end do
!!$
!!$    reli = .true.
!!$    write(*,*) 'ropt88:',t,r
!!$!    write(*,*) real(dx)
!!$!    mu = 1e4
!!$
!!$  end subroutine ropt88
!!$
!!$  function ropt88min(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64), dimension(:), allocatable :: tt, sigx, sigy
!!$    real(REAL64) :: t,f,sigu,sigv,mu,r
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    allocate(tt(n),sigx(n),sigy(n))
!!$!    allocate(tt(n),sig(n),sigx(n),sigy(n))
!!$    tt = cts / pht
!!$    t = x(1)
!!$    r = x(2)
!!$!    sigu = 2e3
!!$!    sigv = 2e3
!!$
!!$!    sigx = sqrt(mu + sigv**2)
!!$!    sigy = sqrt(t*mu + sigu**2)
!!$    !    sig = sqrt(tt**2*sigx**2 + sigy**2)
!!$    sig = sqrt(mju*(t + tt**2) + r**2)
!!$!    n = count(tt > 0)
!!$    f = sum((mju*(t - tt)/sig)**2,tt>0)/2/n - &
!!$!    f = sum(((cts - t*pht)/sig)**2,tt>0)/2/n - &
!!$         sum(log(mju*(mju*t*(tt + 1) + r**2)/sig**3),tt>0)/n
!!$
!!$
!!$!    write(error_unit,*) 'ropt8amin:',real(x),real(f)
!!$
!!$  end function ropt88min


!!$  subroutine ropt7(t,dt,s,phi,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$    use fmm
!!$
!!$    real(REAL64), intent(in out) :: t, dt, s,phi
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 3
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps
!!$    integer :: iter, info, nprint, stat, n, i, j
!!$    logical :: converge
!!$    character :: kind
!!$
!!$    n = size(cts)
!!$
!!$    rmin = minval(pht)
!!$    rmax = maxval(pht)
!!$!    rmin = 5000
!!$!    rmax = 15000
!!$    ttol = 0.1
!!$    write(*,*) rmin, rmax, ttol
!!$
!!$    phi = 45.0/57.3
!!$    do iter = 1, 13
!!$       write(*,*) 'iter:',iter
!!$
!!$       global(1) = t
!!$       global(2) = s
!!$       global(3) = phi
!!$       do i = 1, -n
!!$          allocate(iset(1))
!!$          !          iset = [ (j, j=1,n) ]
!!$          iset(1) = i
!!$!          global(4) = i
!!$          mju = fmin(rmin,rmax,muopt7,ttol)
!!$          if( i < 10 ) write(*,*) mju(i),(pht(i) + t*cts(i))/(1 + t**2)
!!$          deallocate(iset)
!!$       end do
!!$       call mjuopt7(2*iter,t,s,phi)
!!$
!!$       block
!!$         real(REAL64) :: sigu, sigv
!!$         sigu = s*sin(phi)
!!$         sigv = (s/t)*cos(phi)
!!$         open(1,file='/tmp/h')
!!$         do i = 1, n
!!$            write(1,*) mju(i),(pht(i) + t*cts(i))/(1 + t**2), &
!!$                 (cts(i)-t*mju(i))/sqrt(mju(i)*t+sigu**2), &
!!$                 (pht(i)-mju(i))/sqrt(mju(i)+sigv**2),cts(i),pht(i)
!!$         end do
!!$         close(1)
!!$       end block
!!$
!!$       block
!!$         open(1,file='/tmp/fmu')
!!$         global(1) = t
!!$         global(2) = s
!!$         global(3) = phi
!!$         allocate(iset(1))
!!$         do i = 5000, 30000, 100
!!$            mju = i
!!$            iset = 1
!!$            write(1,*) i,ropt7min([t,s,phi]),muopt7(real(i,REAL64))
!!$         end do
!!$         close(1)
!!$         deallocate(iset)
!!$       end block

!!$       x = [ t, s, phi ]
!!$       dx = x / 1000.0
!!$       dx(3) = 1.0/57.3
!!$       call nelmin1(ropt7min,x,dx,f,info)
!!$       t = x(1)
!!$       s = x(2)
!!$       phi = x(3)
!!$       write(*,*) iter, real(t),real(s),real(57.3*phi)
!!$!       if( abs(t - t0)/t < 1e-5 ) exit
!!$
!!$    end do
!!$
!!$    reli = .true.
!!$    write(*,*) 'ropt7:',real(t),real(s),real(57.3*phi),real(mju(1))
!!$
!!$    if( .not. (0 <= phi .and. phi <= 1.57) ) phi = 45.0 / 57.3
!!$!    sigu = s*sin(phi)
!!$!    sigv = s*cos(phi)/t
!!$
!!$!    write(*,*) 'ropt888: sigu, sigv',real(sigu),real(sigv)
!!$
!!$  end subroutine ropt7
!!$
!!$  function ropt7min(x) result(f)
!!$
!!$!    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64), dimension(:), allocatable :: beta, res
!!$    real(REAL64) :: t,f,sigu,sigv,mu,r,q,phi,c2
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    t = x(1)
!!$    r = x(2)
!!$    phi = x(3)
!!$!    mu = mju(1)
!!$!    sigu = s*sin(phi)
!!$!    sigv = s*cos(phi)/t
!!$    c2 = cos(phi)**2
!!$
!!$    sig = sqrt(mju*(t + tt**2) + r**2*(1 + ((tt/t)**2 - 1)*c2))
!!$    beta = mju*t*(tt + 1) + r**2*(1 + (tt/t - 1)*c2)
!!$    res = (mju/sig)*(t - tt)
!!$    f = -sum(res**2)/(2*n) + sum(log(beta))/n - 3*sum(log(sig))/n
!!$!    f = sum(((mu/sig)*(t - tt))**2)/(2*n) - &
!!$!         sum(log((mu/sig)**2*t*(tt + 1)/sig + (mu/sig)*r**2*(1 + (tt/t - 1)*c2)/sig**2))/n
!!$!         sum(log((mju/sig)*(mju*t*(tt + 1) + s**2*(1 + (tt/t - 1)*c2))/sig**2),tt>0)/n
!!$
!!$    return
!!$
!!$    sig = sqrt(mju*(t + tt**2) + sigu**2 + tt**2*sigv**2)
!!$    f = sum((mju*(t - tt)/sig)**2,tt>0)/2/n - &
!!$         sum(log(mju*(mju*t*(tt + 1) + t*tt*sigv**2 + sigu**2)/sig**3),tt>0)/n
!!$
!!$    !    write(error_unit,*) 'ropt7amin:',real(x),real(f)
!!$
!!$  end function ropt7min
!!$
!!$  subroutine mjuopt7(ngroup,t,s,phi)
!!$
!!$    use fmm
!!$
!!$    integer, intent(in) :: ngroup
!!$    real(REAL64), intent(in) :: t,s,phi
!!$    integer, dimension(:), allocatable :: xset
!!$    real(REAL64) :: mumin, mumax, d, xmju, xmin, xmax, ttol, q
!!$    integer :: i,j,n,m
!!$
!!$    n = size(mju)
!!$    global(1) = t
!!$    global(2) = s
!!$    global(3) = phi
!!$
!!$    mumin = min(minval(pht),minval(cts/t))
!!$    mumax = max(maxval(pht),maxval(cts/t))
!!$    mumin = max(mumin,real(1.0,REAL64))
!!$    d = (mumax - mumin) / (ngroup + 1)
!!$    write(*,*) 'b:',mumin,mumax,d
!!$
!!$    if( d < 10*epsilon(d) ) then
!!$       allocate(iset(n))
!!$       iset = [(i,i=1,n)]
!!$       ttol = 0.1 * sqrt(mumin)
!!$       mju = fmin(mumin,mumax,muopt7,ttol)
!!$       deallocate(iset)
!!$       write(*,*) 0, n, mju(1)
!!$       return
!!$    end if
!!$
!!$    do i = 1, ngroup
!!$       m = 0
!!$       allocate(xset(n))
!!$       do j = 1, n
!!$          q = (pht(j) + t*cts(j))/(1+t**2)
!!$!          if( j < 5 ) write(*,*) int((q - mumin) / d),i
!!$          if( int((q - mumin) / d) == i - 1 ) then
!!$!             write(*,*) real(q),real(mumin),real(d),i-1,int((q - mumin) / d) == i - 1
!!$             m = m + 1
!!$             xset(m) = j
!!$          end if
!!$       end do
!!$
!!$       if( m > 0 ) then
!!$
!!$          allocate(iset(m))
!!$          iset = xset(1:m)
!!$
!!$          xmin = mumin + d*(i-1)
!!$          xmax = xmin + d
!!$          ttol = 0.01 * sqrt(xmin) + 10*epsilon(xmin)
!!$          xmju = fmin(xmin,xmax,muopt7,ttol)
!!$
!!$          mju(iset) = xmju
!!$!          do j = 1, size(iset)
!!$!             mju(iset(j)) = xmju
!!$!          end do
!!$          deallocate(iset)
!!$          write(*,*) i,m,xmju
!!$       end if
!!$
!!$       deallocate(xset)
!!$    end do
!!$!    write(*,*) 'mjuopt7:',real(mju)
!!$
!!$
!!$  end subroutine mjuopt7
!!$
!!$
!!$  function muopt7(mu)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: mu
!!$    real(REAL64) :: muopt7, s, t, phi, sigu, sigv, sig2x, sig2y
!!$    integer :: n
!!$
!!$    n = size(iset)
!!$    t = global(1)
!!$    s = global(2)
!!$    phi = global(3)
!!$    sigu = s*sin(phi)
!!$    sigv = s*cos(phi)/t
!!$    sig2y = mu*t + sigu**2
!!$    sig2x = mu + sigv**2
!!$
!!$    muopt7 = sum((cts(iset) - t*mu)**2/sig2y) + sum((pht(iset) - mu)**2/sig2x)
!!$    muopt7 = muopt7 / n
!!$
!!$!    muopt7 = sum((cts(iset) - t*mu)**2 + (pht(iset) - mu)**2)/n
!!$!    muopt7 = sum((cts - t*mu)**2 + (pht - mu)**2)/n
!!$
!!$!    write(*,*) 'muopt7:',real(mu),real(muopt7)
!!$!    if( n == -1 ) &
!!$!         write(*,*) real(cts(n)),real(pht(n)),real(mu), &
!!$!         real((pht(n) + t*cts(n))/(1+t**2)),real(muopt)
!!$
!!$  end function muopt7
!!$
!!$
!!$  function muopt7a(mu)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: mu
!!$    real(REAL64) :: muopt7a, s, t, phi, sigu, sigv, sig2x, sig2y
!!$    integer :: n
!!$
!!$    n = size(iset)
!!$    t = global(1)
!!$    s = global(2)
!!$    phi = global(3)
!!$    sigu = s*sin(phi)
!!$    sigv = s*cos(phi)/t
!!$    sig2y = mu*t + sigu**2
!!$    sig2x = mu + sigv**2
!!$
!!$    muopt7a = sum((cts(iset) - t*mu)**2/sig2y) + sum((pht(iset) - mu)**2/sig2x)
!!$
!!$!    muopt7a = muopt7 / n
!!$
!!$!    muopt7 = sum((cts(iset) - t*mu)**2 + (pht(iset) - mu)**2)/n
!!$!    muopt7 = sum((cts - t*mu)**2 + (pht - mu)**2)/n
!!$
!!$!    write(*,*) 'muopt7:',real(mu),real(muopt7)
!!$!    if( n == -1 ) &
!!$!         write(*,*) real(cts(n)),real(pht(n)),real(mu), &
!!$!         real((pht(n) + t*cts(n))/(1+t**2)),real(muopt)
!!$
!!$  end function muopt7a
!!$
!!$
!!$  subroutine ropt888(t,dt,sig,sigu,sigv,mu,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$    use fmm
!!$
!!$    real(REAL64), intent(in out) :: t, dt, sig,sigu,sigv,mu
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 3
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps, &
!!$         t0,r,phi,teff,s,w
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$    n = size(cts)
!!$    mju = pht
!!$
!!$    t0 = t
!!$    r = sig
!!$    teff = sum(cts/pht)/n
!!$    phi = atan(1.0)
!!$    rmin = minval(pht)
!!$    rmax = maxval(pht)
!!$!    rmin = 5000
!!$!    rmax = 15000
!!$    ttol = 0.1
!!$    write(*,*) rmin, rmax, ttol
!!$
!!$
!!$    ! update mu
!!$    block
!!$      real(REAL64) :: tmax,fmax,lnf,f,mu,tx
!!$      integer :: m
!!$
!!$      tmax = -1
!!$      fmax = -huge(fmax)
!!$      phi = 45.0/57.3
!!$!      t = 2
!!$      sigu = 2e3 !r*sin(phi)
!!$      sigv = 2e3 !r*cos(phi)/t
!!$
!!$      open(3,file='/tmp/mmax')
!!$      ! function tmax = 2-2*(2.6e3/x)**2
!!$      do m = 20,-60
!!$         mu = 10**(m/10.0)
!!$         do i = 1000, 2500
!!$            tx = i/1000.0
!!$            !               mu = mju(1)
!!$            sig = sqrt(mu*(t + tx**2) + sigu**2 + tx**2*sigv**2)
!!$            f = mu*(mu*t*(tx+1)+tx*t*sigv**2+sigu**2)/sig**3 * &
!!$                 exp(-(mu*(tx-t)/sig)**2/2)
!!$            !            f = exp(lnf)
!!$            !            f = lnf
!!$            if( f > fmax ) then
!!$               fmax = f
!!$               tmax = tx
!!$            end if
!!$         end do
!!$         s = 0
!!$         w = 0
!!$         do i = 10,90000
!!$            tx = i / 100.0
!!$            sig = sqrt(mu*(t + tx**2) + sigu**2 + tx**2*sigv**2)
!!$            f = mu*(mu*t*(tx+1)+tx*t*sigv**2+sigu**2)/sig**3 * &
!!$                 exp(-(mu*(tx-t)/sig)**2/2)
!!$            s = s + tx*f
!!$            w = w + f
!!$         end do
!!$         write(*,*) 'max:',fmax,tmax
!!$         if( tmax > 0 ) then
!!$            write(3,*) mu, tmax, s/w
!!$            mju = (sum(pht) + sum(cts)/tmax)/(2*n)
!!$         else
!!$            mju = 1e4
!!$         end if
!!$      end do
!!$      close(3)
!!$
!!$    end block
!!$
!!$
!!$!    sigu = sig
!!$!    sigv = sig
!!$    !    mu = sum((pht + cts/t))/(2*n)
!!$    !    mju = sum(pht)/n
!!$    mju = sum(pht)/n
!!$!    mju = 1e4
!!$    phi = 45.0/57.3
!!$    do iter = 1, 5
!!$!       mju = (pht + cts/t)/2
!!$!       mju = sum(pht)/n
!!$       !       mju = 1e
!!$
!!$       global(1) = t
!!$       global(2) = r
!!$       global(3) = phi
!!$       rmin = minval(pht)
!!$       rmax = maxval(pht)
!!$       sigu = r*sin(phi)
!!$       sigv = (r/t)*cos(phi)
!!$       do i = 1, n
!!$          global(4) = i
!!$          mju(i) = fmin(rmin,rmax,muopt,ttol)
!!$          if( i < 10 ) write(*,*) mju(i),(pht(i) + t*cts(i))/(1 + t**2)
!!$       end do
!!$
!!$       open(1,file='/tmp/h')
!!$       do i = 1, n
!!$          write(1,*) mju(i),(pht(i) + t*cts(i))/(1 + t**2), &
!!$               (cts(i)-t*mju(i))/sqrt(mju(i)*t+sigu**2), &
!!$               (pht(i)-mju(i))/sqrt(mju(i)+sigv**2),cts(i),pht(i)
!!$       end do
!!$       close(1)
!!$
!!$
!!$!       mju = (pht + t*cts)/(1 + t**2)
!!$       mju = sum(mju)/n
!!$!       write(*,*) 'E(mju):',mju(1)
!!$
!!$
!!$
!!$
!!$ !      mju = 1e4
!!$!       write(*,*) 'mu:',sum(pht)/n, sum(cts)/n,(sum(pht)+sum(cts)/t)/(2*n)
!!$       x = [ t, r, phi ]
!!$       dx = x / 1000.0
!!$       dx(3) = 1.0/57.3
!!$       t0 = t
!!$       call nelmin1(ropt888min,x,dx,f,info)
!!$       t = x(1)
!!$       r = x(2)
!!$       phi = x(3)
!!$       write(*,*) iter, real(t),real(t0),real(t-t0),real(mju(1))
!       if( abs(t - t0)/t < 1e-5 ) exit

       ! update mu
!!$       global(1) = t
!!$       global(2) = r
!!$       global(3) = phi
!!$       do i = 5, -5
!!$          global(4) = cts(i)/pht(i)
!!$          block
!!$            integer :: j
!!$            real(REAL64) :: tmut,tmin,tmax
!!$            write(*,*) '****',cts(i),pht(i),tt(i)
!!$            tmin = 1.0
!!$            tmax = 2.5
!!$            open(4,file='/tmp/k')
!!$            do j = 5000,15000,100
!!$               global(4) = j
!!$               tmut = fmin(tmin,tmax,tmutmin,ttol)
!!$               write(4,*) j,cts(i) - tmut*(2*j-pht(i)),tmut
!!$            end do
!!$            close(4)
!!$          end block
!!$!          mju(i) = fmin(rmin,rmax,mumin,ttol)
!!$!          if( i < 5 ) write(*,*) i,real(global(4)),real(mju(i)),real(pht(i))
!!$          stop 0
!!$       end do
!       mju = sum(mju)/n

!!$    end do
!!$
!!$    reli = .true.
!!$    write(*,*) 'ropt888:',real(t),real(r),real(57.3*phi),real(mju(1))
!!$
!!$    if( .not. (0 <= phi .and. phi < 1.57) ) phi = 45.0 / 57.3
!!$    sigu = r*sin(phi)
!!$    sigv = r*cos(phi)/t
!!$
!!$    write(*,*) 'ropt888: sigu, sigv',real(sigu),real(sigv)
!!$!    write(*,*) real(dx)
!!$!    mu = 1e4
!!$
!!$    open(1,file='/tmp/o')
!!$    do i = 1, n
!!$       write(1,*) pht(i),tt(i),(pht(i) + t*cts(i))/(1+t**2)
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine ropt888
!!$
!!$  function ropt888min(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64) :: t,f,sigu,sigv,mu,r,q,phi
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    t = x(1)
!!$    r = x(2)
!!$    phi = x(3)
!!$    sigu = r*sin(phi)
!!$    sigv = r*cos(phi)/t
!!$
!!$    sig = sqrt(mju*(t + tt**2) + sigu**2 + tt**2*sigv**2)
!!$    f = sum((mju*(t - tt)/sig)**2,tt>0)/2/n - &
!!$         sum(log(mju*(mju*t*(tt + 1) + t*tt*sigv**2 + sigu**2)/sig**3),tt>0)/n
!!$
!!$!    write(error_unit,*) 'ropt888min:',real(x),real(f)
!!$
!!$  end function ropt888min
!!$
!!$  function muopt(mu)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: mu
!!$    real(REAL64) :: muopt, r, t, phi, sigu, sigv
!!$    integer :: n
!!$
!!$    t = global(1)
!!$    r = global(2)
!!$    phi = global(3)
!!$    n = nint(global(4))
!!$    sigu = r*sin(phi)
!!$    sigv = r*cos(phi)/t
!!$
!!$    muopt = (cts(n) - t*mu)**2/(mu*t + sigu**2) + &
!!$         (pht(n) - mu)**2/(mu + sigv**2)
!!$
!!$!    muopt = (cts(n) - t*mu)**2 + (pht(n) - mu)**2
!!$
!!$    if( n == -1 ) &
!!$         write(*,*) real(cts(n)),real(pht(n)),real(mu), &
!!$         real((pht(n) + t*cts(n))/(1+t**2)),real(muopt)
!!$
!!$  end function muopt

!!$  function mumin(mu)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: mu
!!$    real(REAL64) :: mumin, mumax, r, t, phi, tt, sigu, sigv, sig
!!$!    integer :: n
!!$
!!$    t = global(1)
!!$    r = global(2)
!!$    phi = global(3)
!!$    tt = global(4)
!!$    sigv = r*sin(phi)
!!$    sigu = (r/t)*cos(phi)
!!$    sig = sqrt(mu*(t+tt**2) + sigu**2 + tt**2*sigv**2)
!!$
!!$    mumax = mu*(mu*t*(tt+1)+t*tt*sigv**2 +sigu**2)/sig**3 * &
!!$         exp(-(mu*(t-tt)/sig)**2/2)
!!$    mumin = -mumax
!!$!    write(*,*) real(t),real(r),real(tt),real(mu),real(mumax)
!!$
!!$  end function mumin

!!$  function tmutmin(tt)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: tt
!!$    real(REAL64) :: tmutmin, f, r, phi, t, sigu, sigv, sig
!!$!    integer :: n
!!$
!!$    t = global(1)
!!$    r = global(2)
!!$    phi = global(3)
!!$    mu = global(4)
!!$    sigv = r*sin(phi)
!!$    sigu = (r/t)*cos(phi)
!!$    sig = sqrt(mu*(t+tt**2) + sigu**2 + tt**2*sigv**2)
!!$
!!$    f = mu*(mu*t*(tt+1)+t*tt*sigv**2 +sigu**2)/sig**3 * &
!!$         exp(-(mu*(t-tt)/sig)**2/2)
!!$    tmutmin = -f
!!$!    write(*,*) real(tt),real(r),real(t),real(mu),real(tmutmin)
!!$
!!$  end function tmutmin



!!$  subroutine ropt9(t,dt,sig,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, dt, sig
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 4
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$    real(REAL64) :: f,t0,mu,phi
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$!    kappa = 0.1
!!$!    r = 2e2
!!$    !    mux = 1e4
!!$
!!$    n = size(cts)
!!$
!!$    phi = 1
!!$    mu = sum(pht)/n
!!$!    mu = 1e4
!!$    x = [ t, sig, phi, mu ]
!!$    dx = x / 1e3
!!$    dx(3) = 0.01
!!$    call nelmin1(ropt9min,x,dx,f,info)
!!$    t = x(1)
!!$    sig = x(2)
!!$    phi = x(3)
!!$    mu = x(4)
!!$!    write(*,*) i, real(t),real(t0),real(t-t0)
!!$
!!$    reli = .true.
!!$
!!$!    t = x(1)
!!$!    sig = x(2)
!!$
!!$    write(*,*) 'ropt9:',real(t),real(sig),real(57.3*phi),real(mu)
!!$
!!$  end subroutine ropt9
!!$
!!$  function ropt9min(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64), dimension(:), allocatable :: tt
!!$    real(REAL64) :: t,f,r,sigu,sigv,phi
!!$    integer :: n
!!$
!!$    n = size(cts)
!!$    allocate(tt(n))
!!$    tt = cts / pht
!!$    t = x(1)
!!$    r = x(2)
!!$    phi = x(3)
!!$    mu = x(4)
!!$
!!$    sigu = r*sin(phi)
!!$    sigv = r*cos(phi)/t
!!$    sig = sqrt(mu*(t + tt**2) + sigu**2 + tt**2*sigv**2)
!!$    f = sum((mu*(t - tt)/sig)**2,tt>0)/2/n - &
!!$         sum(log(mu*(mu*t*(tt + 1) + t*tt*sigv**2 + sigu**2)/sig**3),tt>0)/n
!!$!    write(error_unit,*) 'ropt9min:',real(x),real(f)
!!$
!!$  end function ropt9min


  ! *** TEST ***


  subroutine phoptest(tau,rho,phi)

    real(REAL64), intent(in) :: tau,rho,phi

    call testmaxlike(tau,rho,phi)

  end subroutine phoptest


  subroutine testmaxlike(tau,rho,phi)

    real(REAL64), intent(in) :: tau,rho,phi
    real(REAL64), dimension(3) :: fd,fvec,eps,x,x1,x2
    real(REAL64), dimension(3,3) :: fjac, djac
    real(REAL64) :: f1,f2
    integer :: n, i, j, iflag

    n = size(tt)

    iflag = 1
    x = [ tau, rho, phi ]
!    eps = (epsilon(eps)*(1 + abs(x)))**(1.0/3.0)
    eps = epsilon(eps)**(1.0/3.0)
    write(*,*) 'sx, sy:',real(sx),real(sy)
    write(*,*) '>',real(x)
    call maxlike(3,3,x,fvec,fjac,3,iflag)
    do i = 1, 3
       x = [ tau, rho, phi ]
       x(i) = x(i) + eps(i)
       f1 = elogf(x(1),x(2),x(3))
       x = [ tau, rho, phi ]
       x(i) = x(i) - eps(i)
       f2 = elogf(x(1),x(2),x(3))
       fd(i) = (f1 - f2) / (2*eps(i))
       write(error_unit,*) i,real(fvec(i)),real(fd(i)),&
            real(fd(i)-fvec(i)),real(eps(i))
!            real((fd(i)-fvec(i))/abs(1+fvec(i))),real(eps(i))
    end do

!    return

    iflag = 2
    x = [ tau, rho, phi ]
    eps = epsilon(eps)
    call maxlike(3,3,x,fvec,fjac,3,iflag)
    call fdjac1(testmaxlike_aux,3,x,fvec,djac,3,iflag,3,3,eps,x1,x2)
    do i = 1,3
       do j = i, 3
          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
               real(djac(i,j)), real((fjac(i,j)-djac(i,j))/abs(1+fjac(i,j)))
       end do
    end do

!!$    if( debug_der ) then
!!$       dt = 1e-2
!!$       mu = mju(1)
!!$       open(1,file='/tmp/tder')
!!$
!!$       do j = 50, 600
!!$          tau = j * dt
!!$
!!$          s = 0
!!$          ds = 0
!!$          psi1 = 0
!!$          do i = 1, n
!!$             s = s + psilike(t,mu,tau,r,cos(phi),sin(phi))
!!$             ds = ds + dpsilike(t,mu,tau,r,cos(phi),sin(phi))
!!$          end do
!!$          write(1,*) tau, s / n, ds / n, (s/n - psi1)/dt
!!$          psi1 = s / n
!!$       end do
!!$       close(1)
!!$    end if

  end subroutine testmaxlike

  real(REAL64) function elogf(tau,rho,phi)

    real(REAL64), intent(in) :: tau, rho, phi
    real(REAL64), dimension(:), allocatable :: beta, res, sig
    real(REAL64) :: c2
    integer :: n

    n = size(tt)
    c2 = cos(phi)**2
    sig = sqrt(mju*(tau/sy + tt**2/sx) + rho**2*(1 + ((tt/tau)**2 - 1)*c2))
    beta = mju*tau*(tt/sx + 1/sy) + rho**2*(1 + (tt/tau - 1)*c2)
    res = (mju/sig)*(tau - tt)
    elogf = -sum(res**2)/(2*n) + sum(log(beta/sig**2))/n &
         + sum(log(mju/sig))/n + log(2*PI)/2
    elogf = -elogf

  end function elogf


  subroutine testmaxlike_aux(k,x,fvec,iflag)

    integer, intent(in) :: k
    integer, intent(in out) :: iflag
    real(REAL64), dimension(k), intent(in) :: x
    real(REAL64), dimension(k), intent(out) :: fvec

    real(REAL64), dimension(3,3) :: fjac
    integer :: nflag

    nflag = 1
    call maxlike(3,3,x,fvec,fjac,3,nflag)

  end subroutine testmaxlike_aux




end module photonrate_opt
