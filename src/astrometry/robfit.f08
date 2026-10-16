!
! robfit - robust fitting
!
! Copyright © 2011-2026 F.Hroch (hroch@physics.muni.cz)
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

module astrometry_affinefun

  use oakleaf
  use iso_fortran_env
  implicit none

  logical, parameter, private :: debug = .false.
  real(REAL64), parameter, private :: rad = 57.295779513082322865_REAL64

  type, abstract, extends(bfgs_function) :: astro_bfgs
     real(REAL64), dimension(:), pointer :: x,y,u,v
     real(REAL64) :: refl
     ! the weights of the stars (the inverse errors relative to their
     ! mean square), none: equal weights
     real(REAL64), dimension(:), pointer :: w => null()
   contains
     procedure(fstdsig), deferred :: stdsig
     procedure :: wt
  end type astro_bfgs

  type, extends(astro_bfgs) :: huber_bfgs
   contains
     procedure :: val => rvalue
     procedure :: grad => rgradient
     procedure :: hess => rhessian
     procedure :: stdsig => rstdsig
!     procedure :: fungrad => rfungradient
     procedure :: maxalf => maxalpha
  end type huber_bfgs

  type, extends(astro_bfgs) :: lsq_bfgs
   contains
     procedure :: val => lvalue
     procedure :: grad => lgradient
     procedure :: hess => lhessian
     procedure :: stdsig => lstdsig
  end type lsq_bfgs

  abstract interface

     real(REAL64) function fstdsig(this,x)
       import :: REAL64, astro_bfgs
       class(astro_bfgs) :: this
       real(REAL64), dimension(:), intent(in) :: x
     end function fstdsig

  end interface

contains

  pure real(REAL64) function wt(this,i)

    ! the weight of the star i

    class(astro_bfgs), intent(in) :: this
    integer, intent(in) :: i

    if( associated(this % w) ) then
       wt = this % w(i)
    else
       wt = 1
    end if

  end function wt

  real(REAL64) function lvalue(this,x)

    use astrotrafo

    class(lsq_bfgs) :: this
    real(REAL64), dimension(:), intent(in) :: x
    type(AstrometryGeometry) :: tg
    real(REAL64) :: xx,yy,dx,dy,xoff,yoff,s2
    integer :: i

    xoff = x(1)
    yoff = x(2)
    call tg % init(x(3),x(4),this % refl)

    associate( x => this % x, y => this % y, u => this % u, v => this % v )
      s2 = 0
      do i = 1, size(x)
         call tg % uv2pix(u(i),v(i),xx,yy)
         dx = x(i) - (xx + xoff)
         dy = y(i) - (yy + yoff)
         s2 = s2 + this%wt(i)**2*(dx**2 + dy**2)
      end do
    end associate

    lvalue = s2

  end function lvalue

  subroutine lgradient(this,x,fgrad)

    use astrotrafo

    class(lsq_bfgs) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), intent(out) :: fgrad
    type(AstrometryGeometry) :: tg
    integer :: i
    real(REAL64) :: xx,yy,dx,dy,c,xoff,yoff

    xoff = x(1)
    yoff = x(2)
    c = x(3)
    call tg % init(c,x(4),this % refl)

    associate( x => this % x, y => this % y, u => this % u, v => this % v)

      fgrad = 0
      do i = 1, size(x)
         call tg % uv2pix(u(i),v(i),xx,yy)
         dx = this%wt(i)**2*(x(i) - (xx + xoff))
         dy = this%wt(i)**2*(y(i) - (yy + yoff))
         fgrad(1) = fgrad(1) + dx
         fgrad(2) = fgrad(2) + dy
         fgrad(3) = fgrad(3) + dx*xx + dy*yy
         fgrad(4) = fgrad(4) + dx*yy - dy*xx
      end do
    end associate

    ! d xx/d rot = refl*yy, d yy/d rot = -refl*xx
    fgrad(3) = fgrad(3) / c
    fgrad(4) = this % refl * fgrad(4) / rad
    fgrad = -2*fgrad

  end subroutine lgradient


  subroutine lhessian(this,x,fhess)
    ! the hessian in the minimum, the gradients are supposed nullified

    use astrotrafo

    class(lsq_bfgs) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:,:), intent(out) :: fhess
    real(REAL64), dimension(size(x)) :: fc
    type(AstrometryGeometry) :: tg
    integer :: i
    real(REAL64) :: xx,yy,c,xoff,yoff,dfx,dfy

    xoff = x(1)
    yoff = x(2)
    c = x(3)
    call tg % init(c,x(4),this % refl)
    call this % grad(x,fc)

    associate( x => this % x, y => this % y, u => this % u, v => this % v)

      fhess = 0
      do i = 1, size(x)
         call tg % uv2pix(u(i),v(i),xx,yy)
         dfx = this%wt(i)**2
         dfy = dfx
         fhess(1,1) = fhess(1,1) + dfx
         fhess(1,3) = fhess(1,3) + dfx*xx
         fhess(1,4) = fhess(1,4) + dfx*yy
         fhess(2,2) = fhess(2,2) + dfy
         fhess(2,3) = fhess(2,3) + dfy*yy
         fhess(2,4) = fhess(2,4) + dfy*xx
         fhess(3,3) = fhess(3,3) + dfx*xx**2 + dfy*yy**2
         fhess(4,4) = fhess(4,4) + dfx*yy**2 + dfy*xx**2
      end do

    end associate

    fhess = 2*fhess
    fhess(1,3) = fhess(1,3) / c
    fhess(1,4) = this % refl * fhess(1,4) / rad
    fhess(2,3) = fhess(2,3) / c
    fhess(2,4) =-this % refl * fhess(2,4) / rad
    fhess(3,3) = fhess(3,3) / c**2
    fhess(3,4) = fhess(3,4) / c
    fhess(4,4) = fhess(4,4) / rad**2

    do i = 1, size(fhess,1)
       fhess(i+1:,i) = fhess(i,i+1:)
    end do

  end subroutine lhessian

  function lstdsig(this,x)  result(sig)

    class(lsq_bfgs) :: this
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64) :: sig, s2
    integer :: n

    n = size(this % x)
    s2 = this % val(x)
    if( s2 > 0 ) then
       sig = sqrt( s2 / (2*n - 4))
    else
       sig = 0
    end if

  end function lstdsig

  real(REAL64) function rvalue(this,x)

    use astrotrafo

    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    type(AstrometryGeometry) :: tg
    integer :: i,n
    real(REAL64) :: xx,yy,dx,dy,c,xoff,yoff,ravg,s

    s = x(0)
    xoff = x(1)
    yoff = x(2)
    c = x(3)
    call tg % init(c,x(4),this % refl)

    associate( x => this % x, y => this % y, u => this % u, v => this % v)
      n = size(x)
      ravg = 0
      do i = 1, n
         call tg % uv2pix(u(i),v(i),xx,yy)
         dx = this%wt(i)*(x(i) - (xx + xoff))
         dy = this%wt(i)*(y(i) - (yy + yoff))
         ravg = ravg + ipshuber(dx/s) + ipshuber(dy/s)
      end do
    end associate

    rvalue = ravg / n + log(s)

  end function rvalue

  subroutine rgradient(this,x,fgrad)

    use astrotrafo

    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), dimension(0:), intent(out) :: fgrad
    integer :: i, n
    type(AstrometryGeometry) :: tg
    real(REAL64) :: xx,yy,dx,dy,c,xoff,yoff,fx,fy,s

    s = x(0)
    xoff = x(1)
    yoff = x(2)
    c = x(3)
    call tg % init(c,x(4),this % refl)

    associate( x => this % x, y => this % y, u => this % u, v => this % v)
      n = size(x)
      fgrad = 0
      do i = 1, n
         call tg % uv2pix(u(i),v(i),xx,yy)
         dx = x(i) - (xx + xoff)
         dy = y(i) - (yy + yoff)
         ! the weighted residual w*dx: psi(w*dx/s)*w by the parameters
         fx = this%wt(i)*pshuber(this%wt(i)*dx / s)
         fy = this%wt(i)*pshuber(this%wt(i)*dy / s)
         fgrad(0) = fgrad(0) + fx*dx + fy*dy
         fgrad(1) = fgrad(1) + fx
         fgrad(2) = fgrad(2) + fy
         fgrad(3) = fgrad(3) + fx*xx + fy*yy
         fgrad(4) = fgrad(4) + fx*yy - fy*xx
      end do
    end associate

    fgrad(3:4) = fgrad(3:4) / c
    fgrad(4) = this % refl * c*fgrad(4) / rad
    fgrad(1:)= -fgrad(1:) / (n*s)
    fgrad(0) = -fgrad(0)  / (n*s**2) + 1 / s

    if( debug ) then
       block
         real(REAL64) :: h = 1e-7
         real(REAL64), dimension(0:4) :: dp
         dp = 0
         do i = 0,4
            h = sqrt(epsilon(x))*(1 + abs(x(i)))
            dp(i) = h
            write(error_unit,*) 'fgrad:',(this%val(x+dp) - this%val(x-dp)) / (2*h),fgrad(i)
            dp(i) = 0
         end do
       end block
    end if

  end subroutine rgradient


  subroutine rhessian(this,x,fhess)

    ! the hessian in the minimum, the gradients are supposed nullified

    use astrotrafo

    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), dimension(0:,0:), intent(out) :: fhess
    real(REAL64), dimension(0:size(x)) :: fc
    type(AstrometryGeometry) :: tg
    integer :: i,n
    real(REAL64) :: xx,yy,c,s,xoff,yoff,rx,ry,fx,fy,dfx,dfy,dx,dy

    s = x(0)
    xoff = x(1)
    yoff = x(2)
    c = x(3)
    call tg % init(c,x(4),this % refl)
    call this % grad(x,fc)

    associate( x => this % x, y => this % y, u => this % u, v => this % v)

      fhess = 0
      n = size(x)
      do i = 1, n
         call tg % uv2pix(u(i),v(i),xx,yy)
         dx = x(i) - (xx + xoff)
         dy = y(i) - (yy + yoff)
         rx = dx / s
         ry = dy / s
         fx = this%wt(i)*pshuber(this%wt(i)*rx)
         fy = this%wt(i)*pshuber(this%wt(i)*ry)
         dfx = this%wt(i)**2*dpshuber(this%wt(i)*rx)
         dfy = this%wt(i)**2*dpshuber(this%wt(i)*ry)
         fhess(0,0) = fhess(0,0) + dfx*rx**2 + dfy*ry**2
         fhess(0,1) = fhess(0,1) + dfx*dx
         fhess(0,2) = fhess(0,2) + dfy*dy
         fhess(0,3) = fhess(0,3) + dfx*dx*xx + dfy*dy*yy
         fhess(0,4) = fhess(0,4) + dfx*dx*yy - dfy*dy*xx
         fhess(1,1) = fhess(1,1) + dfx
         fhess(1,3) = fhess(1,3) + dfx*xx
         fhess(1,4) = fhess(1,4) + dfx*yy
         fhess(2,2) = fhess(2,2) + dfy
         fhess(2,3) = fhess(2,3) + dfy*yy
         fhess(2,4) = fhess(2,4) + dfy*xx
         fhess(3,3) = fhess(3,3) + dfx*xx**2 + dfy*yy**2
         fhess(3,4) = fhess(3,4) + (dfx -  dfy)*xx*yy
         fhess(4,4) = fhess(4,4) + dfx*yy**2 + dfy*xx**2
      end do
      ! the derivatives by rotation: d xx/d rot = refl*yy, d yy/d rot = -refl*xx
      fhess = fhess / (n*s**2)
      fhess(0,0) = fhess(0,0) + 1 / s**2
      fhess(0,3) = fhess(0,3) / c
      fhess(0,4) = this % refl * fhess(0,4) / rad
      fhess(0,1:) = fhess(0,1:) / s
      fhess(1,3) = fhess(1,3) / c
      fhess(1,4) = this % refl * fhess(1,4) / rad
      fhess(2,3) = fhess(2,3) / c
      fhess(2,4) =-this % refl * fhess(2,4) / rad
      fhess(3,3) = fhess(3,3) / c**2
      fhess(3,4) = this % refl * fhess(3,4) / c / rad
      fhess(4,4) = fhess(4,4) / rad**2

      do i = lbound(fhess,1), ubound(fhess,1)
         fhess(i+1:,i) = fhess(i,i+1:)
      end do

    end associate

    if( debug ) then
       block
         real(REAL64), dimension(lbound(x,1):ubound(x,1)) :: d,fv1,fv2
         real(REAL64) :: h
         integer :: i

         do i = lbound(x,1), ubound(x,1)
            write(error_unit,'(a,5es15.5)') 'Analytic Hessian: ', fhess(i,:)
         end do

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

         do i = lbound(x,1), ubound(x,1)
            write(error_unit,'(a,5es15.5)') 'Numeric Hessian: ', fhess(i,:)
         end do

       end block
      end if

  end subroutine rhessian

  function rstdsig(this,x)  result(sig)

    use astrotrafo
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN

    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    type(AstrometryGeometry) :: tg
    real(REAL64) :: s,psi2,dpsi,sx,sy,xx,yy,sig
    integer :: n,i

    s = x(0)
    call tg % init(x(3),x(4),this % refl)

    associate( x => this % x, y => this % y, u => this % u, v => this % v)

      n = size(x)
      psi2 = 0
      dpsi = 0
      do i = 1, n
         call tg % uv2pix(u(i),v(i),xx,yy)
         ! both xoff, yoff should be nullified in the minimum
         sx = this%wt(i)*(x(i) - xx) / s
         sy = this%wt(i)*(y(i) - yy) / s
         psi2 = psi2 + huber(sx)**2 + huber(sy)**2
         dpsi = dpsi + dhuber(sx) + dhuber(sy)
      end do

    end associate

    if( psi2 > 0 .and. dpsi > 0 ) then
       sig = s * (sqrt(2*psi2*n) / dpsi)
    else if( dpsi > 0 ) then
       sig = 0
    else
       sig = IEEE_VALUE(sig, IEEE_QUIET_NAN)
    end if

  end function rstdsig

  real(REAL64) function maxalpha(this,x,p)

    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x, p
    real(REAL64), parameter :: maxalf = 1
    real(REAL64) :: s, ps, sc, psc, alf

    maxalpha = maxalf
    s = x(0)
    ps = p(0)
    sc = x(3)
    psc = p(3)

    ! the statistical scale: s > 0
    if( s + ps <= 0 ) then
       alf = -0.9*(s / ps)
       maxalpha = min(alf,maxalpha)
    end if

    ! the scale of affine projection: sc > 0
    if( sc + psc <= 0 ) then
       alf = -0.9*(sc / psc)
       maxalpha = min(alf,maxalpha)
    end if

  end function maxalpha

  subroutine rfungradient(this,x,fval,fgrad)
    ! for now, one is left unused,

    use astrotrafo

    class(huber_bfgs) :: this
    real(REAL64), dimension(0:), intent(in) :: x
    real(REAL64), intent(out) :: fval
    real(REAL64), dimension(0:), intent(out) :: fgrad
    integer :: i, n
    type(AstrometryGeometry) :: tg
    real(REAL64) :: xx,yy,dx,dy,c,xoff,yoff,fx,fy,s,sx,sy,ravg

    s = x(0)
    xoff = x(1)
    yoff = x(2)
    c = x(3)
    call tg % init(c,x(4),this % refl)

    associate( x => this % x, y => this % y, u => this % u, v => this % v)
      n = size(x)
      ravg = 0
      fgrad = 0
      do i = 1, n
         call tg % uv2pix(u(i),v(i),xx,yy)
         dx = x(i) - (xx + xoff)
         dy = y(i) - (yy + yoff)
         sx = dx / s
         sy = dy / s
         fx = pshuber(sx)
         fy = pshuber(sy)
         ravg = ravg + ipshuber(sx) + ipshuber(sy)
         fgrad(0) = fgrad(0) + fx*dx + fy*dy
         fgrad(1) = fgrad(1) + fx
         fgrad(2) = fgrad(2) + fy
         fgrad(3) = fgrad(3) + fx*xx + fy*yy
         fgrad(4) = fgrad(4) + fx*yy - fy*xx
      end do
    end associate

    fval = ravg / n + log(s)
    fgrad(3:4) = fgrad(3:4) / c
    fgrad(4) = this % refl * c*fgrad(4) / rad
    fgrad(1:)= -fgrad(1:) / (n*s)
    fgrad(0) = -fgrad(0)  / (n*s**2) + 1 / s

  end subroutine rfungradient

end module astrometry_affinefun



module robfit

  use iso_fortran_env
  implicit none
  private

  real(REAL64), parameter, private :: rad = 57.295779513082322865_REAL64

  public :: robmin

contains

  subroutine robmin(type,fit,a,d,refl,x,y,s,acen,dacen,dcen,ddcen, &
       sc,dsc,rot,drot,stdsig,verbose,reliable,w)

    ! fitting of the transformation; w (optional) are the weights
    ! of the stars (the inverse errors, relative to their mean square)

    use astrometry_affinefun
    use oakleaf
    use astrotrafo
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN

    character(len=*),intent(in) :: type, fit
    real(REAL64),dimension(:), intent(in) :: a, d
    real(REAL64),dimension(:), target, intent(in) :: x,y
    real(REAL64),intent(in) :: refl
    real(REAL64),intent(in out) :: acen,dcen,sc,rot,s
    real(REAL64),intent(out) :: dacen,ddcen,dsc,drot,stdsig
    logical, intent(in) :: verbose
    logical, intent(out) :: reliable
    real(REAL64), dimension(:), target, intent(in), optional :: w

    real(REAL64), parameter :: macheps = epsilon(1.0_REAL64)
    integer, parameter :: npar = 4

    real(REAL64),dimension(0:npar) :: p,dp
    real(REAL64),dimension(0:npar,0:npar) :: B
    real(REAL64),dimension(:), allocatable :: u,v
    real(REAL64), dimension(:), allocatable, target :: uref,vref
    integer :: iter,i,j,iflag,ndat,maxiter
    real(REAL64) :: xoff,yoff,uoff,voff,gtol
    class(astro_bfgs), allocatable :: fun
    type(AstrometryProjection) :: tp
    type(AstrometryGeometry) :: tg
    logical :: reli

    reliable = .false.

    if( verbose ) then
       if( fit == 'SQUARES' ) then
          write(error_unit,'(a)') '#  === Least-squares fitting ==='
       else
          write(error_unit,'(a)') '#  === Robust fitting ==='
       end if
       write(error_unit,'(a)') &
            '# iter iflag  s[pix]  acen         dcen          offset[pix]        s[pix/deg] rot[deg] reliable'
    end if

    ndat = size(a)

    ! too few stars, or no scale for the robust fit: no solution
    ! (the frame fails, the next frames are processed)
    if( .not. (ndat > npar) ) then
       write(error_unit,'(a,i0,a)') 'Warning: the astrometry fit of ',ndat, &
            ' stars: at least five stars are required.'
       return
    end if
    if( .not. s > 0 .and. fit == 'ROBUST' ) then
       write(error_unit,'(a)') &
            'Warning: the robust astrometry fit needs a statistical scale > 0.'
       return
    end if

    allocate(uref(ndat),vref(ndat),u(ndat),v(ndat))

    ! setup
    s = s * sc
    gtol = sqrt(macheps)

    maxiter = precision(p)
    if( type == '' ) maxiter = 1

    do iter = 1, maxiter

       call tp % init(type,acen,dcen)
       call tp % proj(a,d,uref,vref)

       p(0) = s
       p(1:2) = 0
       p(3:4) = [ sc, rot ]

       if( fit == 'SQUARES' ) then

          fun = lsq_bfgs(x,y,uref,vref,refl)
          if( present(w) ) fun % w => w
          call bfgsmin(fun,p(1:),B(1:,1:),gtol,200*npar,verbose,'',iflag)
          B(0,:) = 0
          B(:,0) = 0
          B(0,0) = 1
          ! s is the root mean square
          s = fun % stdsig(p(1:))
          p(0) = s
          xoff = p(1)
          yoff = p(2)
          sc = p(3)
          rot = p(4)

       else ! fit == 'ROBUST' ) then

          fun = huber_bfgs(x,y,uref,vref,refl)
          if( present(w) ) fun % w => w
          call bfgsmin(fun,p,B,gtol,200*npar,verbose,'',iflag)
          ! s is the scale parameter of a robust function
          s = p(0)
          xoff = p(1)
          yoff = p(2)
          sc = p(3)
          rot = p(4)

          if( iflag == 6 .and. s**2 < 1e3*epsilon(s) ) write(error_unit,*) &
               'Warning in robfit: it looks like data have no statistical errors, try --fit squares.'
       end if
       reli = iflag == 0 .or. iflag == 1
       if( .not. reli ) return

       ! update the centre of projection
       call tg % init(sc,rot,refl)
       call tg % pix2uv(-xoff,-yoff,uoff,voff)
       call tp % invproj(uoff,voff,acen,dcen)

       if( verbose ) &
            write(error_unit,'(a,2i2,f10.3,2f11.5,2f10.3,f10.1,f7.1,1x,l1,a1)') &
            '# rfit ',iter,iflag,s,acen,dcen,xoff,yoff,sc,rot,reli

       if( hypot(xoff,yoff) < 0.01*s .and. reli ) then
          reliable = .true.
          exit
       end if

    end do
    if( .not. reliable ) return

    dp = IEEE_VALUE(dp, IEEE_QUIET_NAN)
    do i = 0, npar
       if( B(i,i) > macheps ) dp(i) = sqrt(B(i,i)/(2*ndat - npar))
    end do

    ! deviations of output parameters, NaNs are quietly accepted
    dacen = dp(1) / sc !/ cos(dcen / rad) ! cos() sqeezes metric (distance) on sphere
    ddcen = dp(2) / sc
    dsc = sc*(dp(3)/p(3))
    drot = dp(4)

    s = s / sc ! [pix] -> [deg]
    stdsig = fun % stdsig(p) / sc

    if( verbose ) then

       block
         character(len=20), parameter :: fmt = '(a,5g15.5)'
         real(REAL64),dimension(0:npar,0:npar) :: cov,hess

         if( fit == 'SQUARES' ) then
            call fun % hess(p(1:),hess(1:,1:))
            hess(0,:) = 0
            hess(:,0) = 0
            hess(0,0) = 1
         else
            call fun % hess(p,hess)
            ! the correction for the precise estimate of stdsig for Huber's function
            hess = hess * (s / stdsig)
         end if
         call qrinv(hess,cov)

         write(error_unit,'(a)') '# quasi-inverse-Hessian at minimum (estimate of covariance):'
         write(error_unit,fmt) ('#',B(i,:),i=0,npar)
         write(error_unit,'(a)') '# Covariance matrix:'
         write(error_unit,fmt) ('#',cov(i,:),i=0,npar)
         write(error_unit,'(a)') '# Correlation matrix:'
         do i = 0,npar
            write(error_unit,'(a)',advance="no")  '#'
            do j = 0,npar
               if( cov(i,i)*cov(j,j) > 0 ) then
                  write(error_unit,'(f9.3)',advance="no") &
                       cov(i,j)/sqrt(cov(i,i)*cov(j,j))
               else
                  write(error_unit,'(a9)',advance="no") ' *?* '
               end if
            end do
            write(error_unit,*)
         end do
         write(error_unit,fmt) '# solution:  ',p
         write(error_unit,fmt) '# deviations:',dp
         write(error_unit,'(a,es10.2,a,f0.3,a)') &
            '# stdsig: ',stdsig,' [deg]    ',3600*stdsig,' [arcsec]'
         write(error_unit,'(a,es10.2,a,f0.3,a)') &
              '# stat. scale: ',s,' [deg]    ',3600*s,' [arcsec]'
       end block
    end if

  end subroutine robmin



end module robfit
