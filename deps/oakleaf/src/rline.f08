!
!  Robust estimate of parameters of line
!
!  Copyright © 2014-5, 2017-8, 2020 F.Hroch (hroch@physics.muni.cz)
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
!  along with OakLeaf.  If not, see <http://www.gnu.org/licenses/>.
!
!


module robustline

  ! This module provides subroutines for estimation
  ! of parameters (a,b) for a linear function in the form
  !
  !     y = a + b*x                        (no errors of x,y)
  !
  ! or
  !                a + b*x
  !     y = ------------------------       (errors included)
  !         sqrt(dy**2 + b**2*dx**2)
  !
  !
  ! with optional estimate of theirs statistical errors (da, db).
  !
  ! rline should be called as:
  !
  !  call rline(x,y,a,b,da,db,dx,dy,sigma,reliable,robustness,verbose)
  !
  ! on input:
  !   x,y - array of data values to be estimated
  !   dx,dy - array of statistical errors of x,y (optional), dx can be omitted,
  !           if dy is presented.
  !   robustness - use or don't use robust methods (optional, .true.)
  !   verbose - verbose output (optional, .false.)
  !
  ! on output are estimated:
  !   a,b - parameters of the line
  !   da,db - standard error of the parameters (optional)
  !   sigma - standard deviation (optional)
  !   reliable - indicates reliability of results (optional)
  !
  ! The weight-free alternative uses "leverage points" smoothing:
  ! all points has added a weight proportional to the deviation
  ! of the measured and expected data.

  use iso_fortran_env

  implicit none
  private

  ! numerical precision of real numbers
  integer, parameter :: dbl = selected_real_kind(15)

  ! print debug messages
  logical :: debug = .false.

  ! a graph of the log-like debug messages
  logical :: debug_graph = .false.

  ! Does to estimate the jacobian by differences or derivations ?
  logical, parameter :: analytic = .false.

  ! Does to use the robust estimators?
  logical :: robust = .true.

  ! weights?
  logical :: weights = .false.

  integer, private :: ndat
  real(dbl), dimension(:), pointer :: xdata, ydata
  real(dbl), dimension(:), allocatable :: xsig, ysig, res, psi, dpsi, w
  real(dbl), private :: rsig

  public :: rline

contains

  subroutine rline(x,y,a,b,da,db,dx,dy,sigma,reliable,robustness,verbose)

    real(dbl), dimension(:), target, intent(in) :: x,y
    real(dbl), intent(in out) :: a,b
    real(dbl), dimension(:), intent(in), optional :: dx,dy
    real(dbl), intent(out), optional :: da,db
    real(dbl), intent(out), optional :: sigma
    logical, intent(out), optional :: reliable
    logical, intent(in), optional :: robustness, verbose

    real(dbl), dimension(:), allocatable :: h
    real(dbl) :: s,d,d1,d2,sig

    logical :: abreli, sreli, reli
    integer :: n

    rsig = 1

    weights = present(dx) .or. present(dy)

    if( present(verbose) ) debug = verbose
    if( present(robustness) ) robust = robustness

    if( weights .and. .not. present(dy) ) stop 'rline(): dy must be present.'

    if( size(x) /= size(y) ) stop 'rline(): Size of x,y data sets are different.'

    if( present(dx) ) then
       if( size(dx) /= size(x) ) stop 'rline(): Size of x and dx are different.'
       if( .not. all(dx >= 0) ) stop 'rline(): all dx >= 0, it is required.'
    end if

    if( present(dy) ) then
       if( size(dy) /= size(x) ) stop 'rline(): Size of y and dy are different.'
       if( .not. all(dy > 0) ) stop 'rline(): all dy > 0, it is required.'
    end if

    n = size(x)
    ndat = n
    xdata => x
    ydata => y

    if( n < 2 ) then
       a = 0
       b = 0
       if( present(da) ) da = 0
       if( present(db) ) db = 0
       if( present(sigma) ) sigma = 1
       if( present(reliable) ) reliable = .false.
       return

    else if( n == 2 ) then
       d = x(2) - x(1)
       if( abs(d) > 0 ) then
          b = (y(2) - y(1)) / (x(2) - x(1))
          a = (y(1) + y(1)) / 2 + b*(x(1) + x(2)) / 2
       else
          a = 0
          b = 0
       end if
       if( present(da) ) da = 0
       if( present(db) ) db = 0
       if( present(sigma) ) sigma = 1
       if( present(reliable) ) reliable = .false.
       return

    end if

    allocate(xsig(n),ysig(n),res(n),psi(n),dpsi(n),w(n))

    if( present(dx) ) then
       xsig = dx
    else
       xsig = 0.0_dbl
    end if

    if( present(dy) ) then
       ysig = dy
    else
       ysig = 1.0_dbl
    end if

    ! the inital estimate by absolute deviations
    call absline(a,b,s,abreli)
    if( robust ) rsig = s
    if( debug ) write(error_unit,'(a,3(1x,g0),1x,2l1)') &
         'Robust line by abs.dev., sig by MAD:', &
         real(a),real(b),real(s),abreli,sreli

    reli = abreli
    if( reli ) then

       ! determine hat matrix, diagonal elements revelates
       ! leverage points, if h > 0.5
       if( robust ) then
          allocate(h(n))
          call leverage(x,h)
          w = sqrt(1 - h)
          deallocate(h)
       else
          w = 1
       end if

       ! update estimates by a robust method
       call logline(a,b,abreli)
       call nscale(a,b,s,sreli)
       if( debug ) write(error_unit,'(a,3(1x,g0),1x,2l1)') &
            'Robust line by likelihood, sig by information:', &
            real(a),real(b),real(s),abreli,sreli

       if( debug_graph ) call graph_like(1.0_dbl,1.0_dbl)

       if( sreli .and. robust ) rsig = s
       reli = abreli .and. sreli
       if( reli ) then

          ! the high-accurate solution, including statistical errors
          call newline(a,b,d1,d2,sig,reli)
          if( reli ) s = sig

       end if

    end if

    if( .not. reli ) then
       d2 = rsig/sqrt(real(n))
       if( b > 0 ) then
          d1 = d2 / b
       else
          d1 = huge(d1)
       end if
    end if

    if( present(da) ) da = d1
    if( present(db) ) db = d2
    if( present(reliable) ) reliable = reli
    if( present(sigma) ) sigma = s

    deallocate(xsig,ysig,res,psi,dpsi,w)

  end subroutine rline

  subroutine leverage(x,h)

    use minpacks

    real(dbl), dimension(:), intent(in) :: x
    real(dbl), dimension(:), intent(out) :: h

    real(dbl), parameter :: maxhat = 1 - epsilon(x)
    real(dbl), dimension(ndat,ndat) :: hat
    real(dbl), dimension(2,2) :: xjac, xjac1
    real(dbl), dimension(ndat,2) :: xdat
    integer :: i

    xdat(:,1) = 1
    xdat(:,2) = x
    xjac = matmul(transpose(xdat),xdat)
    call qrinv(xjac,xjac1)
    hat = matmul(matmul(xdat,xjac1),transpose(xdat))
    h = [ (min(hat(i,i),maxhat),i=1,size(hat,1)) ]

    if( debug ) then
       do i = 1, size(xjac,1)
          write(error_unit,'(a,*(1x,g0))') 'xjac:',real(xjac(i,:))
       end do
       do i = 1, size(hat,1)
          write(error_unit,'(a,*(1x,f8.4))') 'hat:',hat(i,:)
       end do
    end if

  end subroutine leverage



  subroutine absline(a,b,s,reli)

    use NelderMead
    use medians

    real(dbl), intent(in out) :: a,b
    real(dbl), intent(out) :: s
    real(dbl), dimension(2) :: t, dt
    logical, intent(out) :: reli
    real(dbl) :: ams, mad
    integer :: ifault

    t = [ a, b ]
    dt = 0.1*(1 + t)
    call nelmin1(absfun,t,dt,ams,ifault)
    reli = ifault == 0

    if( debug .and. ifault > 0 ) &
         write(error_unit,*) "No convergence in rline: absline()."

    res = residuals(t(1),t(2))

    if( ndat < 50 ) then
       mad = median(abs(res))
    else
       mad = qmedian(abs(res))
    end if

    a = t(1)
    b = t(2)
    s = mad / 0.6745

  end subroutine absline

  function absfun(t) result(s)

    real(dbl), dimension(:), intent(in) :: t
    real(dbl) :: s

    res = residuals(t(1),t(2))
    s = sum(abs(res))

  end function absfun


  subroutine logline(a,b,reli)

    use NelderMead

    real(dbl), intent(in out) :: a,b
    real(dbl), dimension(2) :: t, dt
    logical, intent(out) :: reli
    real(dbl) :: ams
    integer :: ifault

    t = [ a, b ]
    dt = rsig
    call nelmin1(loglikely,t,dt,ams,ifault)
    reli = ifault == 0

    if( debug .and. ifault > 0 ) &
         write(error_unit,*) "No convergence in rline: logline."

    a = t(1)
    b = t(2)

  end subroutine logline


  subroutine newline(a,b,da,db,sig,reli)

    use rfun
    use minpacks

    real(dbl), intent(in out) :: a,b
    real(dbl), intent(out) :: da,db,sig
    real(dbl), dimension(2) :: t
    logical, intent(out) :: reli
    real(dbl), dimension(2,2) :: cov
    real(dbl) :: df, f2, tol
    integer :: nprint,info

    if( debug ) then
       nprint = 1
    else
       nprint = 0
    end if

    t = [ a, b ]
    tol = max(10*ndat*epsilon(t)*(abs(maxval(t)) + 1), 1e-3*rsig)

    ! Try to find the best solution. Our initial estimate is very nearly
    ! to the right one.
    if( analytic ) then
       call lmder3(funder,t,cov,tol,nprint=nprint,info=info)
    else
       call lmdif3(fundif,t,cov,tol,nprint=nprint,info=info)
    end if

    reli = (0 < info .and. info < 5) .and. cov(1,1) > 0 .and. cov(2,2) > 0

    if( reli ) then
       ! the minimum localised successfully, update estimates
       a = t(1)
       b = t(2)

       ! estimate uncerainities
       res = residuals(a,b)
       if( robust ) then
          call tukeys(res/rsig,psi,dpsi)
       else
          psi = res/rsig
          dpsi = 1
       end if
       df = sum(dpsi)
       f2 = sum(psi**2)

       sig = rsig * sqrt(f2/df*ndat/(ndat-2))
       da = sig * sqrt(cov(1,1))
       db = sig * sqrt(cov(2,2))

       if( debug ) then
          write(error_unit,*) 'cov:',real(cov(1,:))
          write(error_unit,*) 'cov:',real(cov(2,:))
       end if

    end if

  end subroutine newline


  function residuals(a,b) result(r)

    real(dbl), intent(in) :: a,b
    real(dbl), dimension(:), allocatable :: r

    if( weights ) then
       r = (ydata - (a + b*xdata)) / sqrt(ysig**2 + b**2*xsig**2)
    else
       r = ydata - (a + b*xdata)
    end if

  end function residuals

  real(dbl) function loglikely(t)

    use rfun

    real(dbl), dimension(:), intent(in) :: t
    real(dbl), dimension(:), allocatable :: rho

    res = residuals(t(1),t(2))

    if( robust ) then
       rho = itukey(res/(rsig*w))
    else
       rho = (res/(rsig*w))**2/2
    end if

    loglikely = sum(rho*(w*rsig)**2) + sum(log(w*rsig))

  end function loglikely


  subroutine funder(m,np,p,fvec,fjac,ldfjac,iflag)

    use rfun

    integer, intent(in) :: m,np,ldfjac
    integer, intent(inout) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(in out) :: fvec
    real(dbl), dimension(ldfjac,np), intent(in out) :: fjac
    real(dbl), dimension(:), allocatable :: ds,dra,drb
    real(dbl) :: a,b,s

    if( iflag == 0 ) then

       if( debug ) then
          write(error_unit,'(4g15.5)') p,fvec
          write(error_unit,*) 'fjac:',fjac(1,:)
          write(error_unit,*) 'fjac:',fjac(2,:)

          block
            real(dbl), dimension(2,2) :: dfjac

            call difjac(p(1),p(2),dfjac)
            write(error_unit,*) 'djac:',dfjac(1,:)
            write(error_unit,*) 'djac:',dfjac(2,:)
          end block
       end if

       return
    end if

    a = p(1)
    b = p(2)
    s = rsig

    if( iflag == 1 ) then

       if( weights ) then

          ds = sqrt(ysig**2 + b**2*xsig**2)
          res = residuals(a,b) / ds
          if( robust ) then
             psi = tukey(res/s)
          else
             psi = res/s
          end if

          fvec(1) = -sum(psi/ds)
          fvec(2) = -sum(psi*(xdata/ds + res*b*xsig**2/ds**2))

       else

          res = residuals(a,b)
          if( robust ) then
             psi = tukey(res/(s*w))
          else
             psi = res / (s*w)
          end if

          fvec(1) = -sum(psi*w)*s
          fvec(2) = -sum(psi*xdata*w)*s

       end if

    else if( iflag == 2 ) then

       res = residuals(a,b)
       if( robust ) then
          dpsi = dtukey(res/(s*w))
       else
          dpsi = 1
       end if

       if( weights ) then

          ds = sqrt(ysig**2 + b**2*xsig**2)
          dra = - 1 / ds
          drb = - (xdata*ds + res*b*xsig**2) / ds**2

          fjac(1,1) = sum(dpsi/ds**2)
          fjac(1,2) = sum(dpsi*dra*drb) - s*sum(psi*b*xsig**2/ds**3)
          fjac(2,1) = fjac(1,2)
          fjac(2,2) = sum(dpsi*drb**2) - &
             s*sum(psi*(res*(1-3*b**2**xsig**2/ds**2)-2*xdata/ds)*xsig**2/ds**2)

       else

          fjac(1,1) = sum(dpsi)
          fjac(1,2) = sum(dpsi*xdata)
          fjac(2,1) = fjac(1,2)
          fjac(2,2) = sum(dpsi*xdata**2)

       end if

    end if

  end subroutine funder

  subroutine difjac(a,b,jac)

    real(dbl), intent(in) :: a,b
    real(dbl), dimension(:,:), intent(out) :: jac
    real(dbl), parameter :: d = epsilon(d)**(1.0/3.0)
    real(dbl), dimension(2) :: fv1,fv2
    integer :: iflag

    iflag = 1

    call fundif(2,2,[a+d,b],fv1,iflag)
    call fundif(2,2,[a-d,b],fv2,iflag)
    jac(1,:) = (fv1 - fv2)/(2*d)

    call fundif(2,2,[a,b+d],fv1,iflag)
    call fundif(2,2,[a,b-d],fv2,iflag)
    jac(2,:) = (fv1 - fv2)/(2*d)

  end subroutine difjac

  subroutine fundif(m,np,p,fvec,iflag)

    use rfun

    integer, intent(in) :: m,np
    integer, intent(in out) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(in out) :: fvec
    real(dbl), dimension(:), allocatable :: ds
    real(dbl) :: a,b,s

    if( iflag == 0 ) then

       if( debug ) then
          write(error_unit,'(a,6g13.5)') ' p,fvec:',p,fvec

          block
            real(dbl), dimension(2,2) :: jac
            real(dbl), dimension(2) :: fv
            integer :: n

            n = 2
            call funder(2,2,p,fv,jac,2,n)
            write(error_unit,*) ' jac:',jac(1,:)
            write(error_unit,*) ' jac:',jac(2,:)
          end block

       end if

       return
    end if

    a = p(1)
    b = p(2)
    s = rsig

    if( weights ) then

       ds = sqrt(ysig**2 + b**2*xsig**2)
       res = residuals(a,b) /  ds
       if( robust ) then
          psi = tukey(res/s)
       else
          psi = res/s
       end if

       fvec(1) = -sum(psi/ds)
       fvec(2) = -sum(psi*(xdata*ds + res*b*xsig**2)/ds**2)

    else

       res = residuals(a,b)

       if( robust ) then
          psi = tukey(res/(w*s))
       else
          psi = res/(w*s)
       end if

       fvec(1) = -sum(psi*w)*s
       fvec(2) = -sum(psi*xdata*w)*s

    end if

  end subroutine fundif

  subroutine nscale(a,b,s,reli)

    use likescale

    real(dbl), intent(in) :: a,b
    real(dbl), intent(in out) :: s
    logical, intent(out) :: reli

    res = residuals(a,b)
    call iscale(res,s,reli)

  end subroutine nscale

  subroutine graph_like(a,b)

    real(dbl), intent(in) :: a,b

    real(dbl) :: x,y
    integer :: i,j

    open(1,file='/tmp/rline')
    do i = -200, 200, 5
       x = a + 0.01*i
       do j = -200, 200, 5
          y = b + 0.01*j
          write(1,*) x,y,loglikely([x,y])
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_like


end module robustline
