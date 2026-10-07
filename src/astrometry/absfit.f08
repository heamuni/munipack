!
! absfit - fitting absolute deviations
!
!
! Copyright © 2011-2, 2015-8 F.Hroch (hroch@physics.muni.cz)
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

module absfit

  use astrotrafo
  use iso_fortran_env

  implicit none

  logical, private :: debug = .true.
  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: rad = 57.295779513082322865_dbl
  real(dbl), dimension(:), allocatable, private :: x,y,u,v
  real(dbl), private :: refl, xcen, ycen
  integer, private :: ndat

  private :: medfun

contains


  function medfun(p)

    ! This function computes mean Manhattan metrics (sum of absolute values)
    ! per point: http://en.wikipedia.org/wiki/Taxicab_geometry

    real(dbl) :: medfun
    real(dbl), dimension(:), intent(in) :: p
    type(AstroTrafoProj) :: t
    real(dbl) :: uu,vv,xx,yy
    integer :: i

!    call trafo_init(t,xcen=xcen,ycen=ycen, &
    call trafo_init(t, &
         scale=p(3),rot=p(4),refl=refl)
!         ucen=p(1),vcen=p(2),

    medfun = 0.0_dbl
    do i = 1, ndat
!       call invaffine(t,x(i),y(i),uu,vv)
       call affine(t,u(i),v(i),xx,yy)
!       medfun = medfun + abs(u(i) - uu) + abs(v(i) - vv)
       medfun = medfun + abs(x(i) - xx - p(1)) + abs(y(i) - yy - p(2))
!       write(error_unit,*) x(i),xx,y(i),yy
    end do
    medfun = medfun / (2*ndat)

!    stop
  end function medfun


!  subroutine  absmin(type,a,d,xc,yc,xx,yy,acen,dcen,sc,dsc,rot,drot,xrefl, &
!       xoff,dxoff,yoff,dyoff,mad,verbose)
  subroutine  absmin(type,a,d,xx,yy,acen,dcen,sc,dsc,rot,drot,xrefl, &
       mad,verbose)

    use NelderMead

    character(len=*),intent(in) :: type
!    real(dbl),intent(in) :: xc,yc
    real(dbl),dimension(:),intent(in) :: a,d,xx,yy
    real(dbl),intent(inout) :: acen,dcen,sc,dsc,rot,drot
    real(dbl),intent(in out) :: mad
    real(dbl),intent(in) :: xrefl!,xoff,dxoff,yoff,dyoff
    logical, intent(in) :: verbose

    real(dbl), parameter :: sqreps = sqrt(epsilon(mad))
    real(dbl), parameter :: reqmin = epsilon(mad)
    integer, parameter :: npar = 4

    real(dbl),dimension(npar) :: p,p0,dp ! = [x0, y0, sc, rot]
    real(dbl) :: uu,vv,dd,pa,pd
    integer :: iter,icount, numres, ifault, i, nra, ndec
    type(AstroTrafoProj) :: t,ti

    if( size(a) < 1 ) stop 'Absfit: no star available.'

    if(verbose) write(error_unit,*) "=== Absolute deviations fitting ==="

    debug = verbose
!    xcen = xc
!    ycen = yc

    ! x,y,u,v arrays are a common working arrays for both absmin and medfun !!
    ndat = size(a)
    allocate(x(ndat),y(ndat),u(ndat),v(ndat))

    ! origin of rectangular coordinates is moved to the center of frame
    ! because that's point the frame is rotated around
    x = xx
    y = yy

!!$    block
!!$!      real(REAL64) :: a0, d0
!!$      write(error_unit,*) '====',xoff,yoff,acen,dcen,sc,rot
!!$      call trafo_init(ti,type,acen=acen,dcen=dcen, &
!!$           xcen=xcen,ycen=ycen,scale=sc,rot=rot,refl=xrefl)
!!$      call invproj(ti,sc*xoff,sc*yoff,acen,dcen)
!!$      write(error_unit,*) '====',xoff,yoff,acen,dcen
!!$!      stop
!!$    end block
!    sc = 1 / sc
!    mad = 1.0 / 3600.0
!    return

    call trafo_init(t,type,acen,dcen)
    call proj(t,a,d,u,v)

    !    p = [ 0*xoff, 0*yoff, sc, rot ]
    p = 0
    p(3:4) = [ sc, rot ]
    dd = 0.1*sc ! star centroids are usually located with the precision 0.1 pix
!    dp = [ max(dxoff,dd), max(dyoff,dd), dsc, drot ]
    dp = [ max(mad,dd), max(mad,dd), dsc, drot ]
    refl = xrefl

    if( debug ) write(error_unit,'(a)') &
         '# ifault  mad[pix]  acen         dcen       offset[pix]    s[pix/deg] rot[deg]'

    do iter = 1, precision(p)

       ! estimate transformation by minimizing of absolute deviations
       p0 = p
       dp = max(dp,sqreps)

       call nelmin(medfun,npar,p0,p,mad,reqmin,dp,1,9999,icount,numres,ifault)
       ! a typical run needs about one thousand of medfun calls
       ! ifault=0 means convergence, 2 means slow convergence

       if( debug ) write(error_unit,'(2i2,f10.2,2f11.5,2es11.2,f10.1,f7.1,2i12)') &
            iter,ifault,mad,acen,dcen,p(1:2),1/p(3),p(4),icount,numres

       ! iterations couldn't make further progress, with precision 1/10 [pix]
!       if( iter > 1 .and. all(abs(p(1:2)/p(3)) < 0.1) ) exit
       if( iter > 1 .and. all(abs(p(1:2)) < 0.1) ) exit

       ! update center of projection and projected coordinates
!       call trafo_init(ti,type,acen=acen,dcen=dcen, &
!            xcen=xcen,ycen=ycen,scale=p(3),rot=p(4),refl=refl)
       call trafo_init(ti,type,acen=acen,dcen=dcen, &
            scale=1/p(3),rot=p(4),refl=refl)
       call invproj(ti,-p(1)*p(3),-p(2)*p(3),acen,dcen)
!       acen = acen + pa
!       dcen = dcen + pd
       call trafo_init(t,type,acen,dcen)
       call proj(t,a,d,u,v)
       p(1:2) = 0.0_dbl

    end do

    call trafo_init(ti,type,acen=acen,dcen=dcen,refl=refl)
    call invproj(ti,-p(1)*p(3),-p(2)*p(3),acen,dcen)

    rot = p(4)
    sc = 1.0_dbl / p(3)

    if( debug ) then
       nra = 0
       ndec = 0
       call trafo_init(t,xcen=xcen,ycen=ycen,scale=p(3),rot=p(4),refl=refl)
       do i = 1, ndat
          call invaffine(t,x(i),y(i),uu,vv)
          if( uu - u(i) > 0 ) nra  = nra + 1
          if( vv - v(i) > 0 ) ndec = ndec + 1
       end do
       write(error_unit,'(a,i7,"/",f0.1,2(i5,"+-",f0.1))') &
            '# absfit sign test (total/expected, RA+, Dec+): ',&
            ndat,ndat/2.0,nra,sqrt(nra*0.25),ndec,sqrt(ndec*0.25)
       write(error_unit,'(a,en12.2,en13.3,f8.1,2f11.5,f5.1)') &
            "# absfit final:",mad,sc,rot,acen,dcen,refl
    end if

    deallocate(x,y,u,v)

  end subroutine absmin

end module absfit
