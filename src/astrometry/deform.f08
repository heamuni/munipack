!
! deform - the plate distortions
!
! Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
!
! This file is part of Munipack.
!
! Munipack is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Munipack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
!
! The reference on TPV projection:
! https://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html
!
! There is implemented the TPV deformation by polynomial
! of 3rd degree, the last term is PVi_11.
!

module astrometry_deformfun

  use oakleaf
  use iso_fortran_env
  implicit none

  logical, parameter, private :: debug = .false.
  real(REAL64), parameter, private :: rad = 57.295779513082322865_REAL64
  integer, parameter, private :: npar = 12

  type, extends(bfgs_function) :: deform_bfgs
     real(REAL64), dimension(:), pointer :: x,y,u,v
     real(REAL64) :: sc,rot,refl
     real(REAL64), dimension(:), pointer :: u1,v1,rdist
   contains
     procedure :: val => pvalue
     procedure :: grad => pgradient
     procedure :: hess => phessian
     procedure :: stdsig => pstdsig
!     procedure :: fungrad => rfungradient
     procedure :: maxalf => palpha
  end type deform_bfgs

contains

  subroutine dist_fun(pv,x,y,u1,v1,rdist)

    real(REAL64), dimension(:,0:), intent(in) :: pv
    real(REAL64), dimension(:), intent(in) :: x,y
    real(REAL64), dimension(:), intent(out) :: u1,v1,rdist
    integer :: i
    real(REAL64) :: r,ux,uy,vy,vx,xy

    do i = 1, size(x)

       r = hypot(x(i),y(i))
       rdist(i) = r

       ux = x(i)*(pv(1,1) + x(i)*(pv(1,4) + x(i)*pv(1,7)))
       uy = y(i)*(pv(1,2) + y(i)*(pv(1,6) + y(i)*pv(1,10)))
       xy = x(i)*y(i)*(pv(1,5) + pv(1,8)*x(i) + pv(1,9)*y(i))
       u1(i) = pv(1,0) + ux + uy + xy + r*(pv(1,3) + pv(1,11)*r**2)

       vy = y(i)*(pv(2,1) + y(i)*(pv(2,4) + y(i)*pv(2,7)))
       vx = x(i)*(pv(2,2) + x(i)*(pv(2,6) + x(i)*pv(2,10)))
       xy = x(i)*y(i)*(pv(2,5) + pv(2,8)*y(i) + pv(2,9)*x(i))
       v1(i) = pv(2,0) + vy + vx + xy + r*(pv(2,3) + pv(2,11)*r**2)

    end do

  end subroutine dist_fun

  real(REAL64) function pvalue(this,x)

    class(deform_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    integer :: i,n
    real(REAL64) :: du,dv,ravg,s
    real(REAL64), dimension(2,0:npar-1) :: pv

    s = x(0)
    pv(1,:) = x(1:npar)
    pv(2,:) = x(npar+1:)

    n = size(this % x)
    ravg = 0
    associate( u1 => this % u1, v1 => this % v1, u => this % u, v => this % v)
      call dist_fun(pv,this % x, this % y, u1, v1, this % rdist)
      do i = 1, n
         du = u(i) - u1(i)
         dv = v(i) - v1(i)
         ravg = ravg + ipshuber(du/s) + ipshuber(dv/s)
       !         write(error_unit,*) x(i),u(i),du,y(i),v(i),dv
      end do
    end associate
    pvalue = ravg / n + log(s)

  end function pvalue

  subroutine pgradient(this,x,fgrad)

    class(deform_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), dimension(0:), intent(out) :: fgrad
    integer :: i, n
    real(REAL64) :: du,dv,fu,fv,s
    real(REAL64), dimension(2,0:npar-1) :: pv

    s = x(0)
    pv(1,:) = x(1:npar)
    pv(2,:) = x(npar+1:)

    associate( u => this % u, v => this % v, u1 => this % u1, v1 => this % v1, &
         r => this % rdist, x => this % x, y => this % y )
      call dist_fun(pv,x,y,u1,v1,r)
      n = size(x)
      fgrad = 0
      do i = 1, n
         du = u(i) - u1(i)
         dv = v(i) - v1(i)
         fu = pshuber(du / s)
         fv = pshuber(dv / s)
         fgrad(0) = fgrad(0) + fu*du + fv*dv
         fgrad(1) = fgrad(1) + fu
         fgrad(2) = fgrad(2) + fu*x(i)
         fgrad(3) = fgrad(3) + fu*y(i)
         fgrad(4) = fgrad(4) + fu*r(i)
         fgrad(5) = fgrad(5) + fu*x(i)**2
         fgrad(6) = fgrad(6) + fu*x(i)*y(i)
         fgrad(7) = fgrad(7) + fu*y(i)**2
         fgrad(8) = fgrad(8) + fu*x(i)**3
         fgrad(9) = fgrad(9) + fu*x(i)**2*y(i)
         fgrad(10) = fgrad(10)+fu*x(i)*y(i)**2
         fgrad(11) = fgrad(11)+fu*y(i)**3
         fgrad(12) = fgrad(12)+fu*r(i)**3
         fgrad(npar+1) = fgrad(npar+1) + fv
         fgrad(npar+2) = fgrad(npar+2) + fv*y(i)
         fgrad(npar+3) = fgrad(npar+3) + fv*x(i)
         fgrad(npar+4) = fgrad(npar+4) + fv*r(i)
         fgrad(npar+5) = fgrad(npar+5) + fv*y(i)**2
         fgrad(npar+6) = fgrad(npar+6) + fv*x(i)*y(i)
         fgrad(npar+7) = fgrad(npar+7) + fv*x(i)**2
         fgrad(npar+8) = fgrad(npar+8) + fv*y(i)**3
         fgrad(npar+9) = fgrad(npar+9) + fv*y(i)**2*x(i)
         fgrad(npar+10) = fgrad(npar+10)+fv*y(i)*x(i)**2
         fgrad(npar+11) = fgrad(npar+11)+fv*x(i)**3
         fgrad(npar+12) = fgrad(npar+12)+fv*r(i)**3
      end do
    end associate

    fgrad(1:)= -fgrad(1:) / (n*s)
    fgrad(0) = -fgrad(0)  / (n*s**2) + 1 / s

    if( debug ) then
       block
         real(REAL64) :: h = 1e-7
         real(REAL64), dimension(0:ubound(x,1)) :: dp
         dp = 0
         do i = 0,ubound(x,1)
            h = sqrt(epsilon(x))*(1 + abs(x(i)))
            dp(i) = h
            write(error_unit,*) 'fgrad:',(this%val(x+dp) - this%val(x-dp)) / (2*h),fgrad(i)
            dp(i) = 0
         end do
       end block
    end if

  end subroutine pgradient


  subroutine phessian(this,x,fhess)

    ! the hessian in the minimum, the gradients are supposed nullified

    class(deform_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), dimension(0:,0:), intent(out) :: fhess
!    real(REAL64), dimension(0:size(x)) :: fc
!    type(AstrometryGeometry) :: tg
!    integer :: i,n
!    real(REAL64) :: xx,yy,c,s,xoff,yoff,rx,ry,fx,fy,dfx,dfy,dx,dy

!!$    goto 33
!!$
!!$    s = x(0)
!!$    xoff = x(1)
!!$    yoff = x(2)
!!$    c = x(3)
!!$!    call tg % init(c,x(4),this % refl)
!!$!    call this % grad(x,fc)
!!$
!!$    associate( x => this % x, y => this % y, u => this % u, v => this % v)
!!$
!!$      fhess = 0
!!$      n = size(x)
!!$      do i = 1, n
!!$!         call tg % uv2pix(u(i),v(i),xx,yy)
!!$         dx = x(i) - (xx + xoff)
!!$         dy = y(i) - (yy + yoff)
!!$         rx = dx / s
!!$         ry = dy / s
!!$         fx = huber(rx)
!!$         fy = huber(ry)
!!$         dfx = dhuber(rx)
!!$         dfy = dhuber(ry)
!!$         fhess(0,0) = fhess(0,0) + dfx*rx**2 + dfy*ry**2
!!$         fhess(0,1) = fhess(0,1) + dfx*dx
!!$         fhess(0,2) = fhess(0,2) + dfy*dy
!!$         fhess(0,3) = fhess(0,3) + dfx*dx*xx + dfy*dy*yy
!!$         fhess(0,4) = fhess(0,4) + dfx*dx*yy - dfy*dy*xx
!!$         fhess(1,1) = fhess(1,1) + dfx
!!$         fhess(1,3) = fhess(1,3) + dfx*xx
!!$         fhess(1,4) = fhess(1,4) + dfx*yy
!!$         fhess(2,2) = fhess(2,2) + dfy
!!$         fhess(2,3) = fhess(2,3) + dfy*yy
!!$         fhess(2,4) = fhess(2,4) + dfy*xx
!!$         fhess(3,3) = fhess(3,3) + dfx*yy**2 + dfy*xx**2
!!$         fhess(3,4) = fhess(3,4) + (dfx -  dfy)*xx*yy
!!$         fhess(4,4) = fhess(4,4) + dfx*xx**2 + dfy*yy**2
!!$      end do
!!$      fhess = fhess / (n*s**2)
!!$      fhess(0,0) = fhess(0,0) + 1 / s**2
!!$      fhess(0,3) = fhess(0,3) / c
!!$      fhess(0,4) = fhess(0,4) / rad
!!$      fhess(0,1:) = fhess(0,1:) / s
!!$      fhess(1,3) = fhess(1,3) / c
!!$      fhess(1,4) = fhess(1,4) / rad
!!$      fhess(2,3) = fhess(2,3) / c
!!$      fhess(2,4) =-fhess(2,4) / rad
!!$      fhess(3,3) = fhess(3,3) / c**2
!!$      fhess(3,4) = fhess(3,4) / c / rad
!!$      fhess(4,4) = fhess(4,4) / rad**2
!!$
!!$      do i = lbound(fhess,1), ubound(fhess,1)
!!$         fhess(i+1:,i) = fhess(i,i+1:)
!!$      end do
!!$
!!$    end associate
!!$
!!$    fhess = 0
!!$    33 continue
!    if( debug ) then
       block
         real(REAL64), dimension(lbound(x,1):ubound(x,1)) :: d,fv1,fv2
         character(len=80) :: fmt
         real(REAL64) :: h
         integer :: i,npar

         npar = ubound(x,1) + 1
         write(fmt,'(a,i0,a)') '(a,',npar,'es15.5)'

!         do i = lbound(x,1), ubound(x,1)
!            write(error_unit,fmt) 'Analytic Hessian: ', fhess(i,:)
!         end do

         ! this is a numerical approximation of the Hessian,
         ! the inverse of this matrix can fail or be undetermined

         d = 0
         do i = lbound(x,1), ubound(x,1)
            h = sqrt(epsilon(h))*(1 + abs(x(i)))
            d(i) = h
            call this%grad(x+d,fv1)
            call this%grad(x-d,fv2)
            fhess(i:,i) = (fv1(i:) - fv2(i:))/(2*h)
            fhess(i,i:) = fhess(i:,i)
            d(i) = 0
         end do

         if( debug ) then
            do i = lbound(x,1), ubound(x,1)
               write(error_unit,fmt) 'Numeric Hessian: ', fhess(i,:)
            end do
         end if

       end block
!      end if

  end subroutine phessian


  function pstdsig(this,x)  result(sig)

    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN

    class(deform_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64) :: s,psi2,dpsi,su,sv,sig
    integer :: n,i
    real(REAL64), dimension(2,0:npar-1) :: pv

    s = x(0)
    pv(1,:) = x(1:npar)
    pv(2,:) = x(npar+1:)

    associate( u => this % u, v => this % v, u1 => this % u1, v1 => this % v1)

      call dist_fun(pv,this % x,this % y,u1,v1,this % rdist)
      n = size(u)
      psi2 = 0
      dpsi = 0
      do i = 1, n
         su = (u(i) - u1(i)) / s
         sv = (v(i) - v1(i)) / s
         psi2 = psi2 + huber(su)**2 + huber(sv)**2
         dpsi = dpsi + dhuber(su) + dhuber(sv)
      end do

    end associate

    if( psi2 > 0 .and. dpsi > 0 ) then
       sig = s * (sqrt(2*psi2*n) / dpsi)
    else if( dpsi > 0 ) then
       sig = 0
    else
       sig = IEEE_VALUE(sig, IEEE_QUIET_NAN)
    end if

  end function pstdsig


  real(REAL64) function palpha(this,x,p)

    class(deform_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x, p
    real(REAL64) :: s, ps

    palpha = 1
    s = x(0)
    ps = p(0)

    ! the statistical scale: s > 0
    if( s + ps <= 0 ) then
       palpha = -0.9*(s / ps)
    end if

  end function palpha

end module astrometry_deformfun

module deform

  use iso_fortran_env
  implicit none
  private

  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: rad = 57.295779513082322865_dbl

  public :: defmin

contains


  subroutine defmin(type,rec,dec,refl,xpix,ypix,s,acen,dcen, &
       sc,rot,stdsig,pv,verbose,reliable)

    ! fitting of the transformation by the robust method

    use astrometry_deformfun
    use oakleaf
    use astrotrafo
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN

    character(len=*),intent(in) :: type
    real(dbl),dimension(:), intent(in) :: rec, dec, xpix,ypix
    real(dbl),intent(in) :: refl
    real(dbl),intent(in) :: acen,dcen,sc,rot
    real(dbl),intent(out) :: s,stdsig
    real(dbl), dimension(:,0:), intent(out) :: pv
    logical, intent(in) :: verbose
    logical, intent(out) :: reliable

    real(dbl), parameter :: macheps = epsilon(1.0_dbl)
    integer, parameter :: npar = 12
    integer, parameter :: maxiter = 1000 + npar*200

    real(dbl),dimension(0:size(pv)) :: p,dp
    real(dbl),dimension(0:size(pv),0:size(pv)) :: B
    real(dbl),dimension(:), allocatable, target :: u,v,u1,v1,rdist
    real(dbl), dimension(:), allocatable, target :: uref,vref, x, y
    integer :: i,j,iflag,ndat
    real(dbl) :: gtol
    class(deform_bfgs), allocatable :: fun
    type(AstrometryProjection) :: tp
    type(AstrometryGeometry) :: tg

    reliable = .false.

    if( verbose ) &
       write(error_unit,'(a)') '#  === Robust fitting with PV deformation ==='

    ndat = size(rec)

    if( .not. (npar == size(pv,2) ) ) &
         stop 'In deform, npar == size(pv,2) must hold.'
    if( .not. (ndat > 2*npar+1) ) &
         stop 'Determination of deformations requires more stars.'
    if( .not. s > 0  ) &
         stop 'Robust fitting needs statistical scale > 0.'

    allocate(uref(ndat),vref(ndat),u(ndat),v(ndat),x(ndat),y(ndat), &
         u1(ndat),v1(ndat),rdist(ndat))

    ! setup
    gtol = sqrt(macheps)
    p = 0
    p(0) = s
    p(2) = 1
    p(npar+2) = 1

    call tp % init(type,acen,dcen)
    call tp % proj(rec,dec,uref,vref)
    call tg % init(sc,rot,refl)
    call tg % pix2uv(xpix,ypix,x,y)

    fun = deform_bfgs(x,y,uref,vref,sc,rot,refl,u1,v1,rdist)
    call bfgsmin(fun,p,B,gtol,maxiter,verbose,'',iflag)
    reliable = iflag == 0 .or. iflag == 1

    if( verbose ) then
       write(error_unit,'(a,i2,f10.3,l2)') &
            '# deformfit ',iflag,s*sc,reliable
    end if

    if( .not. reliable ) return

    s = p(0)
    pv(1,:) = p(1:npar)
    pv(2,:) = p(npar+1:)
    dp = IEEE_VALUE(dp, IEEE_QUIET_NAN)
    do i = 0,2*npar
       if( B(i,i) > macheps ) dp(i) = sqrt(B(i,i)/(2*ndat - 2*npar-1))
    end do
    stdsig = fun%stdsig(p)

    if( verbose ) then

       block
         character(len=20) :: fmt = '(a,5g15.5)'
         real(dbl),dimension(0:size(pv),0:size(pv)) :: cov,hess

         call fun % hess(p,hess)
         hess = hess * (s / stdsig)
         call qrinv(hess,cov)

         write(fmt,'(a,i0,a)') '(a,',2*npar+1,'g15.5)'
         write(error_unit,'(a)') '# quasi-inverse-Hessian at minimum (estimate of covariance):'
         write(error_unit,fmt) ('#',B(i,:),i=0,2*npar)
         write(error_unit,'(a)') '# Covariance matrix:'
         write(error_unit,fmt) ('#',cov(i,:),i=0,2*npar)
         write(error_unit,'(a)') '# Correlation matrix:'
         do i = 0,2*npar
            write(error_unit,'(a)',advance="no")  '#'
            do j = 0,2*npar
               if( cov(i,i)*cov(j,j) > 0 ) then
                  write(error_unit,'(f9.3)',advance="no") &
                       cov(i,j)/sqrt(cov(i,i)*cov(j,j))
               else
                  write(error_unit,'(a9)',advance="no") ' *?* '
               end if
            end do
            write(error_unit,*)
         end do
         write(fmt,'(a,i0,a)') '(a,',2*npar+1,'es13.3)'
         write(error_unit,fmt) '# solution:  ',p
         write(error_unit,fmt) '# deviations:',dp
         write(error_unit,'(a,es10.2,a,f0.3,a)') &
              '# stdsig: ',stdsig,' [deg]    ',3600*stdsig,' [arcsec]'
         write(error_unit,'(a,es10.2,a,f0.3,a)') &
            '# stat. scale: ',s,' [deg]    ',3600*s,' [arcsec]'
       end block
    end if

  end subroutine defmin

end module deform
