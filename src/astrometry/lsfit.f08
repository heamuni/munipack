!
! lsfit - least square fit
!
!
! Copyright © 2011-3, 2015-6, 2018 F.Hroch (hroch@physics.muni.cz)
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

module lsfit

  implicit none

  logical, private :: debug = .false.
  logical, private :: analytic = .true.
  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: rad = 57.295779513082322865_dbl
  real(dbl), dimension(:), allocatable, private :: x,y,u,v
  integer, private :: ndat
  real(dbl), private :: s0, refl, xcen, ycen

  private :: lsminfun,funder,hessian_min

contains

  subroutine lsmin(type,a,d,xc,yc,rf,xx,yy,acen,dacen,dcen,ddcen,sc,dsc,pa,dpa,&
       rms,verbose)

    use astrotrafo
    use minpack
    use minpacks

    ! parameters
    character(len=*),intent(in) :: type
    real(dbl),dimension(:),intent(in) :: a,d,xx,yy
    real(dbl),intent(in) :: xc,yc,rf
    real(dbl),intent(in out) :: acen,dcen,sc,pa
    real(dbl),intent(out) :: dacen,ddcen,dsc,dpa,rms
    logical, intent(in) :: verbose

    real(dbl), parameter :: eps = 1e2*epsilon(eps)

    ! fitting parameters
    integer, parameter :: npar = 4
    real(dbl),dimension(npar) :: p,dp
    real(dbl),dimension(npar,npar) :: fjac,cov

    integer :: i,iter,info,nprint
    type(AstroTrafoProj) :: t,ti

    debug = verbose
    refl = rf
    xcen = xc
    ycen = yc

    if( size(a) < npar ) stop 'N=4 is absolute minimum of provided data.'
    if( size(a) /= size(d) .or. size(a) /= size(xx) .or. size(xx) /= size(yy) )&
         stop 'Lsmin bad dimensions.'

    ! normalization
    ndat = size(a)
    allocate(x(ndat),y(ndat),u(ndat),v(ndat))
    x = xx
    y = yy

    ! setup projection
    call trafo_init(t,type,acen,dcen)
    call proj(t,a,d,u,v)

    ! fitting parameters
    p = [ 0.0_dbl, 0.0_dbl, 1.0_dbl, -pa/rad ]
    p(3) = 1.0_dbl / sc
    p(4) = pa / rad

    if( debug) write(*,*) "#   info   rms  x0    y0    scale   phi   acen   dcen"
    do iter = 1, 10
       nprint = 0 ! or 1

       ! transformation
       if( analytic ) then
          call lmder2(funder,p,eps,nprint,info)
       else
          call lmdif2(lsminfun,p,eps,nprint,info)
       end if
       if( info == 0 ) stop 'lsfit: Improper fit parameters.'

       rms = sqrt(s0 / (ndat - npar))

       ! correction of center of projection
       call trafo_init(ti,type,acen,dcen,xcen=p(1),ycen=p(2),scale=p(3), &
            rot=rad*p(4),refl=refl)
       call invtrafo(ti,0.0_dbl,0.0_dbl,acen,dcen)

       if( debug ) write(*,'(i2,1x,i1,es9.2,2es10.2,es11.3,3f11.5)') &
            iter,info,rms,p(1:2),p(3),pa,acen,dcen

       ! finish when parameters are appropriate small
       ! and iterations couldn't make further progress
       if( iter > 1 .and. all(abs(p(1:2)) < eps)  ) exit

       ! new projection by using of fitted parameters
       call trafo_init(t,type,acen,dcen)
       call proj(t,a,d,u,v)
       p(1:2) = 0.0_dbl

    end do

    sc = 1.0_dbl/p(3)
    pa = rad*p(4)

    ! statistical uncertainities of parameters
    call hessian_min(p,fjac)
    call qrinv(fjac,cov)
    if( debug ) then
       write(*,*) 'Hessian matrix:'
       write(*,'(4g15.3)') (fjac(i,:),i=1,npar)
       write(*,*) 'Covariance matrix:'
       write(*,'(4g15.3)') (cov(i,:),i=1,npar)
    end if
    forall( i = 1:npar ) dp(i) = sqrt(s0*cov(i,i)/(ndat - npar))

    ! output parameters with deviations
    dsc = sc*(dp(3)/p(3))
    dacen = dp(1) / cos(dcen / rad)
    ddcen = dp(2)
    dpa = rad*dp(4)

    if( debug ) then
       write(*,*) "Final solution:"
       write(*,'(a,4g15.5)') '#lsfit solution:',p
       write(*,'(a,4en15.2)') '#lsfit uncerts.: ',dp
       write(*,*) '#lsfit scale [pix/deg]: ',sc,'+-',dsc
       write(*,*) '#lsfit s0: ',s0
    end if

    deallocate(x,y,u,v)

  end subroutine lsmin



  subroutine funder(m,np,p,fvec,fjac,ldfjac,iflag)

    use astrotrafo

    integer, intent(in) :: m,np,ldfjac
    integer, intent(inout) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(out) :: fvec
    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
    real(dbl), dimension(:), allocatable :: du,dv,uu,vv
    real(dbl), dimension(m) :: fc
    real(dbl) :: u1,v1,r,su,sv,s2
    type(AstroTrafoProj) :: t
    integer :: i,j

    if( iflag == 0 ) then
       write(*,'(4g15.5)') p
       return
    end if

    call trafo_init(t,xcen=xcen,ycen=ycen, &
         scale=p(3),rot=rad*p(4),refl=refl)

    allocate(du(ndat),dv(ndat),uu(ndat),vv(ndat))

    do i = 1, ndat
       call invaffine(t,x(i),y(i),u1,v1,uu(i),vv(i))

       ! u,v - projected coordinates of reference stars
       ! x,y - coordinates of stars on CCD
       ! uu,vv - rotation and scale of x,y coordinates
       ! du,dv - differences

       du(i) = u(i) - u1
       dv(i) = v(i) - v1

!       write(*,'(9f12.5)') u(i),v(i),x(i),y(i),3600*du,3600*dv
    end do

    r = p(3)
    fc(1) = sum(du)
    fc(2) = sum(dv)
    fc(3) = sum(du*uu + dv*vv)
    fc(4) = sum(dv*uu - du*vv)

    if( iflag == 1 ) then

       fvec(1:3) = fc(1:3)
       fvec(4) = fc(4)*r
       fvec = 2*fvec

    else if( iflag == 2 ) then

       su = sum(uu)
       sv = sum(vv)
       s2 = sum(uu**2 + vv**2)

       fjac = 0.0_dbl
       fjac(1,1) = ndat
       fjac(1,3) = su
       fjac(1,4) = -sv*r
       fjac(2,2) = ndat
       fjac(2,3) = sv
       fjac(2,4) = su*r
       fjac(3,3) = s2
       fjac(3,4) = fc(4)
       fjac(4,4) = r**2*s2 - fc(3)*r

       do i = 1,size(p)
          do j = 1,i-1
             fjac(i,j) = fjac(j,i)
          end do
       end do
       fjac = 2*fjac

    end if

    s0 = sum(du**2 + dv**2)

    deallocate(du,dv,uu,vv)

  end subroutine funder

  subroutine hessian_min(p,hess)

    real(dbl), dimension(:), intent(in) :: p
    real(dbl), dimension(:,:), intent(out) :: hess
    real(dbl), dimension(size(p)) :: fvec
    integer :: iflag,m

    iflag = 2
    m = size(p)
    call funder(m,m,p,fvec,hess,m,iflag)

  end subroutine hessian_min

  subroutine lsminfun(m,np,p,fvec,iflag)

    use astrotrafo

    integer, intent(in) :: m,np
    integer, intent(inout) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(out) :: fvec
    real(dbl) :: sx,sy,sc,sf,rx,ry,xx,yy,uu,vv
    type(AstroTrafoProj) :: t
    integer :: i

    if( iflag == 0 .and. debug ) write(*,'(a,4en15.3)') '#lsfit t:',p(1:np)

    call trafo_init(t,xcen=xcen,ycen=ycen, &
         scale=p(3),rot=rad*p(4),refl=refl)

    sx = 0.0_dbl; sy = 0.0_dbl; sc = 0.0_dbl; sf = 0.0_dbl; s0 = 0.0_dbl;
    do i = 1, ndat
       call invaffine(t,x(i),y(i),uu,vv,xx,yy)
       rx = u(i) - uu
       ry = v(i) - vv

       sx = sx + rx
       sy = sy + ry
       sc = sc + (rx*xx + ry*yy)
       sf = sf + (rx*yy - ry*xx)

       s0 = s0 + rx**2 + ry**2
!       write(*,'(9f10.5)') u(i),v(i),x(i),y(i),rx,ry
    end do
    fvec = [ -sx,-sy,-sc,sf*p(3) ]

!    if( iflag == 0 .and. debug ) write(*,'(a,6es10.2)') '#lsfit: fvec:',fvec(1:np),sum(fvec**2),s0

  end subroutine lsminfun


end module lsfit
