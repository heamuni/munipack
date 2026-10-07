!
!  Fortran 95+ Minpack module
!
!  Copyright © 2013-8 F.Hroch (hroch@physics.muni.cz)
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
!  This module implements:
!
!    * Fotran 90+ interfaces for original Minpack routines
!    * Simplified wrappers to original functions
!    * Add-ons: a linear equations solver, the inverse matrix routine
!
!  Adaptations:
!
!  * Power of modern Fortran is utilised: Dimensions of arrays are no more
!    required as arguments, working arrays are created transparently.
!
!  * Tolerance limits are set for maximal accuracy (should be slower).
!
!  * Only principal parameters of routines are visible.
!
!  * Added the convenience functions:
!    qrsolve - a solver of linear systems equations
!    qrinv - an inverse matrix routine
!
!    These are implemented via QR factorised matrix, so they will
!    no fail for a singular matrix; the approach is equivalent to SVD.
!
!  Notes:
!
!  * While the original Minpack's function covar.f is not available in all
!    distribution packages (perhaps, due to some mess), I cancel its support.
!

module minpacks

  implicit none

  ! precision for double real
  integer, parameter, private :: dbl = selected_real_kind(15)

contains


  subroutine lmdif2(fcn,x,tol,nprint,info)

    integer, intent(out) :: nprint
    integer, intent(out) :: info
    real(dbl), intent(in) :: tol
    real(dbl), dimension(:), intent(in out) :: x

    interface
       subroutine fcn(m,n,x,fvec,iflag)
         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n, m
         integer, intent(in out) :: iflag
         real(dbl), dimension(m), intent(in) :: x
         real(dbl), dimension(n), intent(out) :: fvec
       end subroutine fcn
    end interface

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,diag
    real(dbl), dimension(size(x),size(x)) :: fjac
    integer, dimension(size(x)) :: ipvt
    real(dbl) :: xtol,ftol, gtol
    integer :: npar,nfev,maxfev

    npar = size(x)
    ftol = tol
    xtol = tol
    gtol = epsilon(gtol)
    maxfev = 200*(npar+1)

    call lmdif(fcn,npar,npar,x,fvec,ftol,xtol,gtol,maxfev,epsilon(0.0_dbl), &
         diag,1,100.0_dbl,nprint,info,nfev,fjac,npar,ipvt,qtf,wa1,wa2,wa3,wa4)

  end subroutine lmdif2


  subroutine lmdif3(fcn,x,cov,tol,nprint,info)

    real(dbl), dimension(:), intent(in out) :: x
    real(dbl), dimension(:,:), intent(out), optional :: cov
    real(dbl), intent(in), optional :: tol
    integer, intent(in), optional :: nprint
    integer, intent(out), optional :: info

    interface
       subroutine fcn(m,n,x,fvec,iflag)
         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n, m
         integer, intent(in out) :: iflag
         real(dbl), dimension(m), intent(in) :: x
         real(dbl), dimension(n), intent(out) :: fvec
       end subroutine fcn
    end interface

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,diag
    real(dbl), dimension(size(x),size(x)) :: fjac
    integer, dimension(size(x)) :: ipvt
    real(dbl) :: xtol,ftol, gtol,factor,epsfcn
    integer :: npar,nfev,maxfev,nprints, infos,mode, iflag

    npar = size(x)
    if( present(tol) ) then
       ftol = tol
       xtol = tol
    else
       ftol = epsilon(ftol)
       xtol = ftol
    end if
    gtol = 0.0_dbl
    epsfcn = 0.0_dbl
    maxfev = 200*(npar+1)
    factor = 100
    mode = 1
    if( present(nprint) ) then
       nprints = nprint
    else
       nprints = 0
    end if

    call lmdif(fcn,npar,npar,x,fvec,ftol,xtol,gtol,maxfev,epsfcn, &
         diag,mode,factor,nprints,infos,nfev,fjac,npar,ipvt,qtf,wa1,wa2,wa3,wa4)

    if( present(cov) )then
       iflag = 2
       call fdjac2(fcn,npar,npar,x,fvec,fjac,npar,iflag,epsfcn,wa4)
       call qrinv(fjac,cov)
    end if

    if( present(info) ) info = infos

  end subroutine lmdif3


  subroutine lmder2(fcn,x,tol,nprint,info)

    integer, intent(out) :: nprint
    integer, intent(out) :: info
    real(dbl), intent(in) :: tol
    real(dbl), dimension(:), intent(in out) :: x

    interface
       subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n, m, ldfjac
         integer, intent(in out) :: iflag
         real(dbl), dimension(n), intent(in) :: x
         real(dbl), dimension(m), intent(out) :: fvec
         real(dbl), dimension(ldfjac,n), intent(out) :: fjac
       end subroutine fcn
    end interface

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,diag
    real(dbl), dimension(size(x),size(x)) :: fjac
    integer, dimension(size(x)) :: ipvt
    real(dbl) :: xtol,ftol, gtol
    integer :: npar,nfev,njev,maxfev

    npar = size(x)
    ftol = tol
    xtol = tol
    gtol = epsilon(gtol)
    maxfev = 200*(npar+1)

    call lmder(fcn,npar,npar,x,fvec,fjac,npar,ftol,xtol,gtol,maxfev, &
         diag,1,100.0_dbl,nprint,info,nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)

  end subroutine lmder2

  subroutine lmder3(fcn,x,cov,tol,nprint,info)

    real(dbl), dimension(:), intent(in out) :: x
    real(dbl), dimension(:,:), optional, intent(out) :: cov
    real(dbl), optional, intent(in) :: tol
    integer, optional, intent(in) :: nprint
    integer, optional, intent(out) :: info

    interface
       subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n, m, ldfjac
         integer, intent(in out) :: iflag
         real(dbl), dimension(n), intent(in) :: x
         real(dbl), dimension(m), intent(out) :: fvec
         real(dbl), dimension(ldfjac,n), intent(out) :: fjac
       end subroutine fcn
    end interface

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,diag
    real(dbl), dimension(size(x),size(x)) :: fjac
    integer, dimension(size(x)) :: ipvt
    real(dbl) :: xtol,ftol, gtol, factor
    integer :: npar,nfev,njev,maxfev,iflag,infos,nprints,mode

    npar = size(x)
    if( present(tol) ) then
       ftol = tol
       xtol = tol
    else
       ftol = epsilon(ftol)
       xtol = ftol
    end if
    gtol = 0.0_dbl
    factor = 100
    maxfev = 100*(npar+1)
    infos = 0
    mode = 1

    if( present(nprint) ) then
       nprints = nprint
    else
       nprints = 0
    end if

    call lmder(fcn,npar,npar,x,fvec,fjac,npar,ftol,xtol,gtol,maxfev, &
         diag,mode,factor,nprints,infos,nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)

    if( present(cov) ) then
       iflag = 2
       call fcn(npar,npar,x,fvec,fjac,npar,iflag)
       call qrinv(fjac,cov)
    end if

    if( present(info) ) info = infos

  end subroutine lmder3


  ! Functions hybr[d,j]2 are currently unused,
  ! because Levendberg-Marquart (lm*2) provides the regularised way).
  subroutine hybrd2(fcn,x,info)

    integer, intent(out) :: info
    real(dbl), dimension(:), intent(in out) :: x

    interface
       subroutine fcn(n,x,fvec,iflag)
         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n
         integer, intent(in out) :: iflag
         real(dbl), dimension(n), intent(in) :: x
         real(dbl), dimension(n), intent(out) :: fvec
       end subroutine fcn
    end interface

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,diag
    real(dbl), dimension(size(x),size(x)) :: fjac
    real(dbl), dimension(size(x)**2) :: xr
    real(dbl) :: xtol
    integer :: npar,nfev,maxfev

    npar = size(x)
    xtol = epsilon(x)
    maxfev = 200*(npar+1)

    call hybrd(fcn,npar,x,fvec,xtol,maxfev,npar-1,npar-1,0.0_dbl, &
         diag,1,100.0_dbl,1,info,nfev,fjac,npar,xr,size(xr),qtf,wa1,wa2,wa3,wa4)

  end subroutine hybrd2


  subroutine hybrj2(fcn,x,info)

    integer, intent(out) :: info
    real(dbl), dimension(:), intent(in out) :: x

    interface
       subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n, ldfjac
         integer, intent(in out) :: iflag
         real(dbl), dimension(n), intent(in) :: x
         real(dbl), dimension(n), intent(out) :: fvec
         real(dbl), dimension(ldfjac,n), intent(out) :: fjac
       end subroutine fcn
    end interface

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,diag
    real(dbl), dimension(size(x),size(x)) :: fjac
    real(dbl), dimension((size(x)*(size(x)+1))/2) :: xr
    real(dbl) :: xtol
    integer :: npar,nfev,njev,maxfev

    npar = size(x)
    xtol = epsilon(x)
    maxfev = 200*(npar+1)

    call hybrj(fcn,npar,x,fvec,fjac,npar,xtol,maxfev, &
         diag,1,100.0_dbl,1,info,nfev,njev,xr,size(xr),qtf,wa1,wa2,wa3,wa4)

  end subroutine hybrj2


! ----------------------------------------------------------------------
!
!


  subroutine qrsolve(a,b,x)

    real(dbl), dimension(:,:), intent(in) :: a
    real(dbl), dimension(:), intent(in) :: b
    real(dbl), dimension(:), intent(out) :: x

    integer :: m,n,i
    real(dbl), dimension(size(a,1),size(a,2)) :: r, q
    integer, dimension(size(a,1)) :: ipvt
    real(dbl), dimension(size(a,2)) :: rdiag,wa,qtb,diag

    m = size(a,1)
    n = size(a,2)

    ! form the r matrix, r is upper trinagle (without diagonal)
    ! of the factorized a, diagonal is presented in rdiag
    q = a
    call qrfac(m,n,q,m,.true.,ipvt,n,rdiag,diag,wa)

!    write(*,*) 'qrfac:',q,rdiag,diag,ipvt

    ! form R, upper triangular
    r = q
    forall( i = 1:n ) r(i,i) = rdiag(i)

    ! form Q orthogonal matrix
    call qform(m,n,q,m,wa)
    qtb = matmul(transpose(q),b)

    ! accurate up to machine epsilon
    diag = diag*epsilon(diag)    ! lmpar.f:206,

    call qrsolv(n,r,m,ipvt,diag,qtb,x,rdiag,wa)

  end subroutine qrsolve


  subroutine qrinv(a,ainv)

    ! Compute inverse matrix in least-square sense by the use
    ! of the QR factorization. An efficiency does not matter.
    ! A singular matrix check is implemented.
    !
    ! http://en.wikipedia.org/wiki/QR_decomposition
    !     - see: "Solution of inverse problems"

    implicit none

    real(dbl), dimension(:,:), intent(in) :: a
    real(dbl), dimension(:,:), intent(out) :: ainv

    integer :: m,n,i,j
    real(dbl), dimension(size(a,1),size(a,2)) :: r, q
    integer, dimension(size(a,1)) :: ipvt
    real(dbl), dimension(size(a,1)) :: rdiag,acnorm,wa,x,b

    m = size(a,1)
    n = size(a,2)

    ! form the r matrix, r is upper triangle (without diagonal)
    ! of the factorized a, its diagonal is stored in rdiag
    q = transpose(a)
    call qrfac(m,n,q,m,.false.,ipvt,n,rdiag,acnorm,wa)

!    write(*,*) 'qrfac:',q,rdiag,ipvt

    ! form R, upper triangular
    r = q
    forall( i = 1:n ) r(i,i) = rdiag(i)

    ! form Q orthogonal matrix
    call qform(m,n,q,n,wa)

!    do i = 1,n
!       write(*,'(a,4f15.7)') 'q:',q(i,:)
!    end do
!    do i = 1,n
!       write(*,'(a,4f15.7)') 'r:',r(i,:)
!    end do

    ! determine the inverse matrix by substitution
    do i = 1, n
       b = 0.0_dbl
       b(i) = 1.0_dbl

       ! forward subtitution with checking on singular values
       ! (for overdetermined problem, the inverse matrix has
       ! set elements to zero when no inverse is possible).
       ! We assumes: no element index like x(1:0) is acceptable.
       do j = 1,n
          if( abs(r(j,j)) > epsilon(r) ) then
             x(j) = (b(j) - sum(x(1:j-1)*r(1:j-1,j))) / r(j,j)
          else
             x(j) = 0.0_dbl
          end if
       end do
       ainv(:,i) = matmul(q,x)
    end do

!    write(*,*) 'ainv:',ainv

  end subroutine qrinv


end module minpacks
