!
!  Fortran 95+ Minpack module
!
!  Copyright © 2013-8, 2020 F.Hroch (hroch@physics.muni.cz)
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
    gtol = real(0.0,dbl)
    maxfev = 200*(npar+1)

    call lmdif(fcn,npar,npar,x,fvec,ftol,xtol,gtol,maxfev,epsilon(0.0_dbl), &
         diag,1,100.0_dbl,nprint,info,nfev,fjac,npar,ipvt,qtf,wa1,wa2,wa3,wa4)

  end subroutine lmdif2


  subroutine lmdif3(fcn,x,cov,tol,maxev,diag,factor,nprint,info,msg)

    real(dbl), dimension(:), intent(in out) :: x
    real(dbl), dimension(:,:), intent(out), optional :: cov
    real(dbl), dimension(:), intent(in), optional :: diag
    real(dbl), intent(in), optional :: tol, factor
    integer, intent(in), optional :: nprint, maxev
    integer, intent(out), optional :: info
    character(len=*), intent(out), optional :: msg

    interface
       subroutine fcn(m,n,x,fvec,iflag)
         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n, m
         integer, intent(in out) :: iflag
         real(dbl), dimension(m), intent(in) :: x
         real(dbl), dimension(n), intent(in out) :: fvec
       end subroutine fcn
    end interface

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,xdiag
    real(dbl), dimension(size(x),size(fvec)) :: fjac
    integer, dimension(size(x)) :: ipvt
    real(dbl) :: xtol,ftol, gtol,xfactor
    integer :: npar,nfev,maxfev,nprints, infos,mode, iflag

    npar = size(x)
    if( present(tol) ) then
       ftol = tol
       xtol = tol
    else
       ftol = sqrt(epsilon(ftol))
       xtol = ftol
    end if
    gtol = real(0.0,dbl)
    if( present(factor) ) then
       xfactor = factor
    else
       xfactor = 100
    end if
    if( present(maxev) ) then
       maxfev = maxev
    else
       maxfev = 200*(npar+1)
    end if
    if( present(diag) ) then
       xdiag = diag
       mode = 2
    else
       mode = 1
    end if
    if( present(nprint) ) then
       nprints = nprint
    else
       nprints = 0
    end if

    call lmdif(fcn,npar,npar,x,fvec,ftol,xtol,gtol,maxfev,epsilon(x), &
         xdiag,mode,xfactor,nprints,infos,nfev,fjac,npar,ipvt,qtf, &
         wa1,wa2,wa3,wa4)

    if( present(cov) ) then
       if( infos > 0 )then

!!$       block
!!$         integer :: i
!!$         real(dbl),dimension(size(x)) :: acnorm
!!$         real(dbl) :: tols
!!$
!!$!         fjac1 = fjac
!!$
!!$         iflag = 1
!!$!         epsfcn = epsilon(epsfcn)
!!$         call fdjac2(fcn,npar,npar,x,fvec,fjac,npar,iflag,epsilon(x),wa4)
!!$!         fjac = transpose(fjac)
!!$         call qrfac(npar,npar,fjac,npar,.true.,ipvt,npar,xdiag,acnorm,wa2)
!!$         do i = 1, npar
!!$            fjac(i,i) = xdiag(i)
!!$         end do
!!$         write(*,*) ipvt
!!$         tols = epsilon(tol)
!!$         call covars(npar,fjac,npar,ipvt,tols)
!!$!         call covar(npar,fjac,npar,ipvt,tols,wa2)
!!$!         cov = fjac
!!$        write(*,*) 'covar:'
!!$         do i = 1, size(fjac,1)
!!$           write(*,*) real(fjac(i,:))
!!$        end do
!!$     end block

       iflag = 1
!       epsfcn =
       call fdjac2(fcn,npar,npar,x,fvec,fjac,npar,iflag,epsilon(x),wa4)
       !       fjac = fjac1
!       fjac(2,2) = 0

!       fjac1(1,:) = [ 1, 2, 0]
!       fjac1(2,:) = [ 3, 4, 0]
       call qrinv(fjac,cov)

!!$       block
!!$         integer :: i
!!$         real(dbl),dimension(3,3) :: fjac1
!!$         real(dbl),dimension(3,3) :: cov1
!!$!         real(dbl), dimension(size(x),size(x)) :: fjac1
!!$
!!$!       iflag = 1
!!$!       epsfcn = epsilon(epsfcn)
!!$!       call fdjac2(fcn,npar,npar,x,fvec,fjac1,npar,iflag,epsfcn,wa4)
!!$
!!$       fjac1(1,:) = [ 1, 2, 0]
!!$       fjac1(2,:) = [ 3, 4, 0]
!!$       fjac1(2,:) = [ 0, 0, 0]
!!$       call qrinv(fjac1,cov1)
!!$
!!$         write(*,*) 'qrinv:'
!!$         do i = 1, 2
!!$            write(*,*) real(cov1(i,1:2))
!!$         end do
!!$         fjac1 = matmul(cov1,fjac1)
!!$         do i = 1, 2
!!$            write(*,*) real(fjac1(i,1:2))
!!$         end do
!!$!         stop
!!$       end block

    else
       cov = 0
    end if

    end if

    if( present(info) ) info = infos

    if( present(msg) ) then
       if( info == 0 ) then
          msg = 'improper input parameters.'
       else if( info == 1 ) then
          msg = 'both actual and predicted relative reductions' // &
               ' in the sum of squares are at most ftol.'
       else if( info == 2 ) then
          msg = 'relative error between two consecutive iterates' // &
               ' is at most xtol.'
       else if( info == 3 ) then
          msg = 'conditions for info = 1 and info = 2 both hold.'
       else if( info == 4 ) then
          msg = 'the cosine of the angle between fvec and any' // &
               ' column of the jacobian is at most gtol in absolute value. '
       else if( info == 5 ) then
          msg = 'number of calls to fcn has reached or exceeded maxfev.'
       else if( info == 6 ) then
          msg = 'ftol is too small. no further reduction in' // &
               ' the sum of squares is possible.'
       else if( info == 7 ) then
          msg = 'xtol is too small. no further improvement in' // &
               ' the approximate solution x is possible.'
       else if( info == 8 ) then
          msg = 'gtol is too small. fvec is orthogonal to the' // &
               ' columns of the jacobian to machine precision.'
       else if( info < 0 ) then
          msg = 'terminated on call.'
       else
          msg = ''
       end if
    end if

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
    gtol = real(0.0,dbl)
    maxfev = 200*(npar+1)

    call lmder(fcn,npar,npar,x,fvec,fjac,npar,ftol,xtol,gtol,maxfev, &
         diag,1,100.0_dbl,nprint,info,nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)

  end subroutine lmder2

  subroutine lmder3(fcn,x,cov,tol,maxev,diag,factor,nprint,info,msg)

    real(dbl), dimension(:), intent(in out) :: x
    real(dbl), dimension(:,:), optional, intent(out) :: cov
    real(dbl), dimension(:), optional, intent(in) :: diag
    real(dbl), optional, intent(in) :: tol, factor
    integer, optional, intent(in) :: nprint, maxev
    integer, optional, intent(out) :: info
    character(len=*), intent(out), optional :: msg

    interface
       subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n, m, ldfjac
         integer, intent(in out) :: iflag
         real(dbl), dimension(n), intent(in) :: x
         real(dbl), dimension(m), intent(in out) :: fvec
         real(dbl), dimension(ldfjac,n), intent(in out) :: fjac
       end subroutine fcn
    end interface

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,xdiag
    real(dbl), dimension(size(x),size(fvec)) :: fjac
    integer, dimension(size(x)) :: ipvt
    real(dbl) :: xtol,ftol, gtol, xfactor
    integer :: npar,nfev,njev,maxfev,iflag,infos,nprints,mode

    npar = size(x)
    if( present(tol) ) then
       ftol = tol
       xtol = tol
    else
       ftol = sqrt(epsilon(ftol))
       xtol = ftol
    end if
    gtol = real(0.0,dbl)
    if( present(factor) ) then
       xfactor = factor
    else
       xfactor = 100
    end if
    if( present(maxev) ) then
       maxfev = maxev
    else
       maxfev = 200*(npar+1)
    end if
    if( present(diag) ) then
       xdiag = diag
       mode = 2
    else
       mode = 1
    end if
    if( present(nprint) ) then
       nprints = nprint
    else
       nprints = 0
    end if
    infos = 0

    call lmder(fcn,npar,npar,x,fvec,fjac,npar,ftol,xtol,gtol,maxfev, &
         xdiag,mode,xfactor,nprints,infos,nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)

    if( present(cov) ) then
       if( infos > 0  ) then
!       call qrfac(npar,npar,fjac,npar,.true.,ipvt,npar,xdiag,acnorm,wa2)
!       call covars(npar,fjac,npar,ipvt,tol)
!       cov = fjac
          iflag = 2
          call fcn(npar,npar,x,fvec,fjac,npar,iflag)
          call qrinv(fjac,cov)
       else
          cov = 0
       end if
    end if

    if( present(info) ) info = infos

    if( present(msg) ) then
       if( info == 0 ) then
          msg = 'improper input parameters.'
       else if( info == 1 ) then
          msg = 'both actual and predicted relative reductions' // &
               ' in the sum of squares are at most ftol.'
       else if( info == 2 ) then
          msg = 'relative error between two consecutive iterates' // &
               ' is at most xtol.'
       else if( info == 3 ) then
          msg = 'conditions for info = 1 and info = 2 both hold.'
       else if( info == 4 ) then
          msg = 'the cosine of the angle between fvec and any' // &
               ' column of the jacobian is at most gtol in absolute value. '
       else if( info == 5 ) then
          msg = 'number of calls to fcn has reached or exceeded maxfev.'
       else if( info == 6 ) then
          msg = 'ftol is too small. no further reduction in' // &
               ' the sum of squares is possible.'
       else if( info == 7 ) then
          msg = 'xtol is too small. no further improvement in' // &
               ' the approximate solution x is possible.'
       else if( info == 8 ) then
          msg = 'gtol is too small. fvec is orthogonal to the' // &
               ' columns of the jacobian to machine precision.'
       else if( info < 0 ) then
          msg = 'terminated on call.'
       else
          msg = ''
       end if
    end if

  end subroutine lmder3


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


  subroutine hybrj3(fcn,x,tol,maxev,factor,diag,nprint,info,msg)

    real(dbl), dimension(:), intent(in out) :: x
    real(dbl), intent(in), optional :: tol, factor
    real(dbl), dimension(:), intent(in), optional :: diag
    integer, intent(in), optional :: nprint, maxev
    integer, intent(out), optional :: info
    character(len=*), intent(out), optional :: msg

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

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,xdiag
    real(dbl), dimension(size(x),size(x)) :: fjac
    real(dbl), dimension((size(x)*(size(x)+1))/2) :: xr
    real(dbl) :: xtol,ftol, gtol,xfactor
    integer :: npar,nfev,njev,maxfev,nprints, infos, mode

    npar = size(x)

    if( present(tol) ) then
       ftol = tol
       xtol = tol
    else
       ftol = epsilon(ftol)
       xtol = ftol
    end if
    if( present(factor) ) then
       xfactor = factor
    else
       xfactor = 100
    end if
    if( present(diag) ) then
       xdiag = diag
    else
       xdiag = 1
    end if
    if( present(maxev) ) then
       maxfev = maxev
    else
       maxfev = 200*(npar+1)
    end if
    if( present(nprint) ) then
       nprints = nprint
    else
       nprints = 0
    end if
    gtol = 10*epsilon(gtol)
    mode = 2

    call hybrj(fcn,npar,x,fvec,fjac,npar,xtol,maxfev, &
         xdiag,mode,xfactor,nprints,infos,nfev,njev,xr,size(xr), &
         qtf,wa1,wa2,wa3,wa4)

    if( present(info) ) info = infos

    if( present(msg) ) then
       if( info == 0 ) then
          msg = 'improper input parameters.'
       else if( info == 1 ) then
          msg = 'relative error between two consecutive iterates is' // &
               ' at most xtol.'
       else if( info == 2 ) then
          msg = 'number of calls to fcn with iflag = 1 has reached maxfev.'
       else if( info == 3 ) then
          msg = 'xtol is too small. no further improvement in ' // &
               ' the approximate solution x is possible.'
       else if( info == 4 ) then
          msg = 'iteration is not making good progress, as' // &
               ' measured by the improvement from the last' // &
               ' five jacobian evaluations.'
       else if( info == 5 ) then
          msg = ' iteration is not making good progress, as' // &
               ' measured by the improvement from the last' // &
               ' ten iterations.'
       else if( info < 0 ) then
          msg = 'terminated on call.'
       else
          msg = ''
       end if
    end if

  end subroutine hybrj3

  subroutine hybrd3(fcn,x,tol,maxev,factor,diag,nprint,info,msg)

    real(dbl), dimension(:), intent(in out) :: x
    real(dbl), optional, intent(in) :: tol, factor
    integer, optional, intent(in) :: nprint, maxev
    real(dbl), dimension(:), intent(in), optional :: diag
    integer, optional, intent(out) :: info
    character(len=*), intent(out), optional :: msg

    interface
       subroutine fcn(n,x,fvec,iflag)
         integer, parameter :: dbl = selected_real_kind(15)
         integer, intent(in) :: n
         integer, intent(in out) :: iflag
         real(dbl), dimension(n), intent(in) :: x
         real(dbl), dimension(n), intent(out) :: fvec
       end subroutine fcn
    end interface

    real(dbl), dimension(size(x)) :: fvec,qtf,wa1,wa2,wa3,wa4,xdiag
    real(dbl), dimension(size(x),size(x)) :: fjac
    real(dbl), dimension(size(x)**2) :: xr
    real(dbl) :: xtol,xfactor
    integer :: npar,nfev,maxfev,infos,nprints,mode

    npar = size(x)

    if( present(tol) ) then
       xtol = tol
    else
       xtol = epsilon(xtol)
    end if
    if( present(factor) ) then
       xfactor = factor
    else
       xfactor = 100
    end if
    if( present(maxev) ) then
       maxfev = maxev
    else
       maxfev = 200*(npar+1)
    end if
    if( present(diag) ) then
       xdiag = diag
    else
       xdiag = 1
    end if
    infos = 0
    mode = 2


    if( present(nprint) ) then
       nprints = nprint
    else
       nprints = 0
    end if

    call hybrd(fcn,npar,x,fvec,xtol,maxfev,npar-1,npar-1,epsilon(x), &
         xdiag,mode,xfactor,nprints,infos,nfev,fjac,npar,xr,size(xr), &
         qtf,wa1,wa2,wa3,wa4)

    if( present(info) ) info = infos

    if( present(msg) ) then
       if( info == 0 ) then
          msg = 'improper input parameters.'
       else if( info == 1 ) then
          msg = 'relative error between two consecutive iterates is at most xtol.'
       else if( info == 2 ) then
          msg = 'number of calls to fcn with iflag = 1 has reached maxfev.'
       else if( info == 3 ) then
          msg = 'xtol is too small. no further improvement in ' // &
               'the approximate solution x is possible.'
       else if( info == 4 ) then
          msg = 'iteration is not making good progress, as ' // &
               'measured by the improvement from the last ' // &
               'five jacobian evaluations.'
       else if( info == 5 ) then
          msg = 'iteration is not making good progress, as ' // &
               'measured by the improvement from the last ' // &
               'ten iterations.'
       else if( info < 0 ) then
          msg = 'terminated on call.'
       else
          msg = ''
       end if
    end if

  end subroutine hybrd3


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

!   write(*,*) 'qrfac:',q,rdiag,diag,ipvt

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

    ! It computes the inverse matrix in least-square sense
    ! (possibly rank deficient)  by the use of the QR factorization
    ! Any efficiency does not matter. The input should be a square matrix.
    !
    ! http://en.wikipedia.org/wiki/QR_decomposition
    !     - see: "Solution of inverse problems"

    real(dbl), dimension(:,:), intent(in) :: a
    real(dbl), dimension(:,:), intent(out) :: ainv
!    real(dbl), optional, intent(in) :: tol

    integer :: m,n,i!,j
    real(dbl), dimension(size(a,1),size(a,2)) :: r, q
    integer, dimension(size(a,2)) :: ipvt
    real(dbl), dimension(size(a,2)) :: rdiag,acnorm,wa,x,qtb
!    real(dbl) :: rtol

!    if( present(tol) ) then
!       rtol = tol
!    else
!       rtol = epsilon(rtol)
!    end if

    m = size(a,1)
    n = size(a,2)

    ! form the r matrix, r is upper triangle (without diagonal)
    ! of the factorized a, its diagonal is stored in rdiag
 !   q = transpose(a)
    q = a
    call qrfac(m,n,q,m,.true.,ipvt,n,rdiag,acnorm,wa)

!    write(*,*) 'qrfac:',q,rdiag,ipvt

    ! form R, upper triangular only
    r = q
    forall( i = 1:n ) r(i,i) = rdiag(i)

    ! form Q orthogonal matrix
    call qform(m,n,q,n,wa)

!    do i = 1,n
!       write(*,'(a,*(f15.7))') 'q:',q(i,:)
!    end do
!    do i = 1,n
!       write(*,'(a,*(f15.7))') 'r:',r(i,:)
!    end do

    ! accurate up to machine epsilon
    acnorm = acnorm*epsilon(acnorm)    ! lmpar.f:206,

    ! Q is replaced by transposed Q since this point
!    q = transpose(q)
    do i = 1, m
!       b = 0.0_dbl
!       b(ipvt(i)) = 1.0_dbl

       !       qtb = matmul(qt,b)
       qtb = q(ipvt(i),:)
       ! this is a shortcut for:
       ! b = 0
       ! b(ipvt(i)) = 1
       ! qtb = matmul(transpose(q),b)

       call qrsolv(n,r,m,ipvt,acnorm,qtb,x,rdiag,wa)
       ainv(:,ipvt(i)) = x !matmul(q,x)
    end do

!    block
!      integer :: i,j
!      do i = 1, n
!         do j = 1, m
!            write(*,*) sum(a(i,:)*ainv(:,j))
!         end do
!      end do
!    end block

!    return


!!$    ! determine the inverse matrix by substitution
!!$    do i = 1, n
!!$       b = 0.0_dbl
!!$       b(i) = 1.0_dbl
!!$
!!$       ! forward subtitution with checking on singular values
!!$       ! (for overdetermined problem, the inverse matrix has
!!$       ! set elements to zero when no inverse is possible).
!!$       ! We assumes: no element index like x(1:0) is acceptable.
!!$       do j = 1,n
!!$          if( abs(r(j,j)) > rtol ) then
!!$             x(j) = (b(j) - sum(x(1:j-1)*r(1:j-1,j))) / r(j,j)
!!$          else
!!$             x(j) = 0.0_dbl
!!$          end if
!!$       end do
!!$       ainv(:,i) = matmul(q,x)
!!$    end do

!    write(*,*) 'ainv:',ainv

  end subroutine qrinv




  subroutine covars(n,r,ldr,ipvt,tol)

    ! covars is a F90+ implementation of covar() subroutine
    ! by Minpack manual
    !

    integer, intent(in) :: n,ldr
    integer, dimension(:), intent(in) :: ipvt
    real(dbl), intent(in) :: tol
    real(dbl), dimension(ldr,n), intent(in out) :: r
    real(dbl), dimension(n) :: wa

    real(dbl), parameter :: zero = real(0.0,dbl)
    real(dbl), parameter :: one = real(1.0,dbl)
    real(dbl) :: temp, tolr
    integer :: i, ii, j, jj, k, l
    logical :: sing

    tolr = tol*abs(r(1,1))
    l = 0
    do k = 1, n
       if( abs(r(k,k)) < tolr ) exit
       r(k,k) = one / r(k,k)
       do j = 1, k - 1
          temp = r(k,k)*r(j,k)
          r(j,k) = zero
          do i = 1, j
             r(i,k) = r(i,k) - temp*r(i,j)
          end do
       end do
       l = k
    end do

    do k = 1, l
       do j = 1, k - 1
          temp = r(j,k)
          do i = 1, j
             r(i,j) = r(i,j) + temp*r(i,k)
          end do
       end do
!       temp = r(k,k)
!       do i = 1, k
!          r(i,k) = temp*r(i,k)
!       end do
       r(1:k,k) = r(k,k)*r(1:k,k)
    end do

    do j = 1, n
       jj = ipvt(j)
       sing = j > l
       do i = 1, j
          if( sing ) r(i,j) = zero
          ii = ipvt(i)
          if( ii > jj ) r(ii,jj) = r(i,j)
          if( ii < jj ) r(jj,ii) = r(i,j)
       end do
       wa(jj) = r(j,j)
    end do

    do j = 1, n
!       do i = 1, j
!          r(i,j) = r(j,i)
!       end do
       r(1:j,j) = r(j,1:j)
       r(j,j) = wa(j)
    end do

  end subroutine covars


end module minpacks
