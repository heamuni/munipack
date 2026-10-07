!
! Astrometry related transformations using WCSLIB
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

module astrotrafo

  use wcslib
  use iso_fortran_env
  implicit none
  private

  ! Transformation and projection:
  !
  ! type ... type of projection (GNOMONIC,..)
  !      ... mostly affine transformation ...
  ! scale... scale [pixel/deg]
  ! rot  ... rotation angle [deg]
  ! acen ... center projection in alpha [deg]
  ! dcen ... center projection in delta [deg]
  ! xcen ... horizontal center projection [pixel]
  ! ycen ... vertical center projection [pixel]
  ! refl ... -1 for reflected frame
  ! err  ... random statistical error of acen, dcen [deg]

  type, public :: AstrometryGeometry
     real(REAL64), private :: scale,rot,refl
     real(REAL64), dimension(2,2) :: cd, icd
   contains
     procedure :: init => initgeom
     procedure :: pix2uv => touv
     procedure :: uv2pix => topix
     procedure :: get_fcd
  end type AstrometryGeometry

  type, public ::  AstrometryProjection
     type(wcspar), private :: wcs
   contains
     procedure :: init => initbyhand
     procedure, private :: proj1 => cel2uv1
     procedure, private :: projX => cel2uv
     generic :: proj => proj1, projX
     procedure, private :: invprojX => uv2cel
     procedure, private :: invproj1 => uv2cel1
     generic :: invproj => invproj1, invprojX
     final :: deproject
  end type AstrometryProjection

  type, public :: AstroWCS
     type(wcspar), private :: wcs
   contains
     procedure, private :: bypar => initbypar
     procedure, private :: byhead => initwcs
     generic :: init => bypar, byhead
     procedure :: pix2img
     generic :: pix2cel => pix2celX, pix2cel1
     procedure :: pix2celX
     procedure :: pix2cel1
     generic :: cel2pix => cel2pixX, cel2pix1
     procedure :: cel2pixX
     procedure :: cel2pix1
     procedure :: get_code
     procedure :: get_affines
     procedure :: get_crval
     procedure :: get_crpix
     procedure :: get_crder
     procedure :: get_radesys
     procedure :: get_epoch
     procedure :: print_header
     final :: dealloc
  end type AstroWCS

  logical, parameter, private :: debug = .false.
  real(REAL64), parameter, private :: rad = 57.295779513082322865_REAL64
  real(REAL64), parameter, private :: UNDEFINED = 987654321.0e99_REAL64

contains

  ! AstrometryGeometry

  subroutine initgeom(t,scale,rot,refl)

    class(AstrometryGeometry), intent(out) :: t
    real(REAL64), intent(in) :: scale,rot,refl

    real(REAL64), dimension(2,2) :: a,b
    real(REAL64) :: f,c,s

    t%scale= scale
    t%rot  = rot
    t%refl = refl

    f = rot / rad
    c = cos(f)
    s = sin(f)

    ! the affine projection matrix
    a(1,:) = [    refl, 0.0_REAL64 ]
    a(2,:) = [ 0.0_REAL64, 1.0_REAL64 ]
    b(1,:) = [ c, s ]
    b(2,:) = [-s, c ]
    t%cd = scale * matmul(a,b)

    ! .. and the inverse: (A*B)^-1 = B^T * A^-1 = B^T * A
    b(1,:) = [ c,-s]
    b(2,:) = [ s, c]
    t%icd = matmul(b,a) / scale

  end subroutine initgeom

  pure subroutine get_fcd(this, cd)

    ! get x-axis fliped cd

    class(AstrometryGeometry), intent(in) :: this
    real(REAL64), dimension(:,:), intent(out) :: cd

    real(REAL64), dimension(2,2) :: m

    m(1,:) = real([-1, 0 ],REAL64)
    m(2,:) = real([ 0, 1 ],REAL64)
    cd = matmul(m,this % icd)

  end subroutine get_fcd

  pure elemental subroutine touv(t,x,y,u,v)

    class(AstrometryGeometry), intent(in) :: t
    real(REAL64), intent(in) :: x,y
    real(REAL64), intent(out) :: u,v

    real(REAL64), dimension(2) :: r

    r = matmul(t%icd,[x,y])
    u = r(1)
    v = r(2)

  end subroutine touv

  pure elemental subroutine topix(t,u,v,x,y)

    class(AstrometryGeometry), intent(in) :: t
    real(REAL64), intent(in) :: u,v
    real(REAL64), intent(out) :: x,y

    real(REAL64), dimension(2) :: r

    r = matmul(t%cd,[u,v])
    x = r(1)
    y = r(2)

  end subroutine topix


! AstrometryProjection

  subroutine initbyhand(this,type,acen,dcen)

    class(AstrometryProjection), intent(out) :: this
    character(len=*), intent(in) :: type
    real(REAL64), intent(in) :: acen,dcen
    character(len=3) :: code
    integer :: status1, status2, status

    code = type(1:3)
    status = wcserr_enable(1)

    associate( wcs => this % wcs)

      wcs = wcs_alloc()
      status = wcs_ini(2,wcs)
      if( status /= 0 ) goto 666

      if( type /= '' ) then

         status1 = wcs_put(wcs,WCS_CRVAL,acen,1,0)
         status2 = wcs_put(wcs,WCS_CRVAL,dcen,2,0)
         if( status1 + status2 > 0 ) goto 666

         status1 = wcs_put(wcs,WCS_CTYPE,'RA---'//code,1,0)
         status2 = wcs_put(wcs,WCS_CTYPE,'DEC--'//code,2,0)
         if( status1 + status2 > 0 ) goto 666

         status = wcs_put(wcs,WCS_PC,-1.0_REAL64,1,1) ! flip of x-axis
      else
         ! no projection is applied
         continue
      end if

      status = wcs_set(wcs)
      if( status /= 0 ) goto 666

      if( debug ) then
         block
           integer, parameter :: WCSHDO_all = 255
           integer, parameter :: LEN = 80
           character(len=:), allocatable :: header
           integer :: status, nkeyrec, i,j

           status = wcs_header_output(WCSHDO_all,this % wcs,nkeyrec,header)
           if( status == 0 ) then
              do i = 1, nkeyrec
                 j = LEN*(i - 1)
                 write(error_unit,*) header(j+1:j+LEN)
              end do
           end if
           deallocate(header)
         end block
      end if

    end associate
    return

666 continue

    call wcs_perr(this % wcs)
    stop 'Attempt to initialise WCSLIB by hand failed.'

  end subroutine initbyhand

  subroutine cel2uv(this,a,d,u,v)

    class(AstrometryProjection), intent(in) :: this
    real(REAL64), dimension(:), intent(in) :: a,d
    real(REAL64), dimension(:), intent(out) :: u,v
    real(REAL64), dimension(size(a),2) :: world,img,pix
    real(REAL64), dimension(size(a)) :: phi, theta
    integer, dimension(size(a)) :: stat
    integer :: status

    world(:,1) = a
    world(:,2) = d
    status = wcs_world2pixel(this % wcs,world,phi,theta,img,pix,stat)
    u = pix(:,1)
    v = pix(:,2)

  end subroutine cel2uv

  subroutine cel2uv1(this,a,d,u,v)

    class(AstrometryProjection), intent(in) :: this
    real(REAL64), intent(in) :: a,d
    real(REAL64), intent(out) :: u,v
    real(REAL64), dimension(1,2) :: world,img,pix
    real(REAL64), dimension(1) :: phi, theta
    integer, dimension(1) :: stat
    integer :: status

    world(1,:) = [ a, d ]
    status = wcs_world2pixel(this % wcs,world,phi,theta,img,pix,stat)
    u = pix(1,1)
    v = pix(1,2)

  end subroutine cel2uv1

  subroutine uv2cel(this,u,v,a,d)

    class(AstrometryProjection), intent(in) :: this
    real(REAL64), dimension(:), intent(out) :: a,d
    real(REAL64), dimension(:), intent(in) :: u,v
    real(REAL64), dimension(size(a),2) :: world,img,pix
    real(REAL64), dimension(size(a)) :: phi, theta
    integer, dimension(size(a)) :: stat
    integer :: status

    pix(:,1) = u
    pix(:,2) = v
    status = wcs_pixel2world(this % wcs,pix,img,phi,theta,world,stat)
    a = world(:,1)
    d = world(:,2)

  end subroutine uv2cel

  subroutine uv2cel1(this,u,v,a,d)

    class(AstrometryProjection), intent(in) :: this
    real(REAL64), intent(out) :: a,d
    real(REAL64), intent(in) :: u,v
    real(REAL64), dimension(1,2) :: world,img,pix
    real(REAL64), dimension(1) :: phi, theta
    integer, dimension(1) :: stat
    integer :: status

    pix(1,:) = [ u, v ]
    status = wcs_pixel2world(this % wcs,pix,img,phi,theta,world,stat)
    a = world(1,1)
    d = world(1,2)

  end subroutine uv2cel1

  subroutine deproject(this)
    type(AstrometryProjection), intent(in out) :: this
    integer :: status

    call wcs_free(this % wcs)
    status = wcserr_enable(0)

  end subroutine deproject

  ! -- AstroWCS

  subroutine initbypar(this,code,crval,crpix,scale,rot,refl,pv)

    class(AstroWCS), intent(out) :: this
    character(len=*), intent(in) :: code
    real(REAL64), dimension(:), intent(in) :: crval,crpix
    real(REAL64), intent(in) :: scale,rot,refl
    real(REAL64), dimension(:,:), intent(in), optional :: pv
    real(REAL64), dimension(2,2) :: a,b,pc
    real(REAL64) :: f,s,c
    integer :: i,j,status, status1, status2

    if( .not. (size(crpix) == 2 .and. size(crval) == size(crpix)) ) &
         stop 'AstroWCS: dimensions in initbypar() does not match.'
    if( .not. (scale > 0) ) stop 'AstroWCS: scale > 0 unsatisfied'

!    if( code == '' ) stop 'AstroWCS: the projection code is empty.'

    status = wcserr_enable(1)

    associate( wcs => this % wcs )
      wcs = wcs_alloc()
      status = wcs_ini(2,wcs)

      if( status /= 0 ) goto 666

      status1 = wcs_put(wcs,WCS_CTYPE,'RA---'//code(1:3),1,0)
      status2 = wcs_put(wcs,WCS_CTYPE,'DEC--'//code(1:3),2,0)
      if( status1 + status2 > 0 ) goto 666

      status = 0
      do i = 1, size(crpix)
         status = status + wcs_put(wcs,WCS_CRPIX,crpix(i),i,0)
      end do
      if( status > 0 ) goto 666

      status = 0
      do i = 1, size(crval)
         status = status + wcs_put(wcs,WCS_CRVAL,crval(i),i,0)
      end do
      if( status > 0 ) goto 666

      f = rot / rad
      c = cos(f)
      s = sin(f)

      ! the affine projection matrix, the inverse of AstrometryGeometry%cd,
      ! pc = diag(-1,1) * B^T * diag(refl,1) / scale;
      ! the x-axis is flipped (taken negative) due the left-hand convention in RA
      a(1,:) = [       refl, 0.0_REAL64 ]
      a(2,:) = [ 0.0_REAL64, 1.0_REAL64 ]
      b(1,:) = [ c,-s ]
      b(2,:) = [ s, c ]
      pc = matmul(b,a) / scale
      pc(1,:) = - pc(1,:)
      do i = 1, 2
         do j = 1,2
            status = wcs_put(wcs,WCS_PC,pc(i,j),i,j)
         end do
      end do

      if( present(pv) ) then
         do i = 1, size(pv,1)
            do j = 1, size(pv,2)
               status = wcs_put(wcs,WCS_PV, pv(i,j), i, j-1)
            end do
         end do
         if( .not. (code == 'TPV') ) &
              write(error_unit,*) 'Warning: PV must have set the code=TPV.'
      end if

      status = wcs_set(wcs)
      if( status /= 0 ) goto 666

    end associate

    return

666 continue

    call wcs_perr(this % wcs)
    stop 'Attempt to initialise WCSLIB by parameters failed.'

  end subroutine initbypar


  subroutine initwcs(this,header,nkeys,status)

    class(AstroWCS), intent(out) :: this
    character(len=*), intent(in) :: header
    integer, intent(in) :: nkeys
    integer, intent(out) :: status
    integer, parameter :: ctrl = 0
    integer :: nreject, nwcs

    status = wcserr_enable(1)

    status = wcs_image_header(header,nkeys,WCSHDR_all,ctrl, nreject,nwcs,this % wcs)
    if( status /= 0 ) goto 666

    status = wcs_set(this % wcs)
    if( status /= 0 ) goto 666

    return

666 continue

    call wcs_perr(this % wcs)
    stop 'Attempt to initialise WCSLIB by header failed.'

  end subroutine initwcs

  subroutine pix2img(this,x,y,ximg,yimg)

    class(AstroWCS), intent(in) :: this
    real(REAL64), dimension(:), intent(out) :: ximg,yimg
    real(REAL64), dimension(:), intent(in) :: x,y
    real(REAL64), dimension(size(x),2) :: world,img,pix
    real(REAL64), dimension(size(x)) :: phi, theta
    integer, dimension(size(x)) :: stat
    integer :: status

    pix(:,1) = x
    pix(:,2) = y
    status = wcs_pixel2world(this % wcs,pix,img,phi,theta,world,stat)
    ximg = img(:,1)
    yimg = img(:,2)

  end subroutine pix2img

  subroutine pix2celX(this,x,y,a,d)

    class(AstroWCS), intent(in) :: this
    real(REAL64), dimension(:), intent(out) :: a,d
    real(REAL64), dimension(:), intent(in) :: x,y
    real(REAL64), dimension(size(a),2) :: world,img,pix
    real(REAL64), dimension(size(a)) :: phi, theta
    integer, dimension(size(a)) :: stat
    integer :: status

    pix(:,1) = x
    pix(:,2) = y
    status = wcs_pixel2world(this % wcs,pix,img,phi,theta,world,stat)
    a = world(:,1)
    d = world(:,2)

  end subroutine pix2celX

  subroutine pix2cel1(this,x,y,a,d)

    class(AstroWCS), intent(in) :: this
    real(REAL64), intent(out) :: a,d
    real(REAL64), intent(in) :: x,y
    real(REAL64), dimension(1,2) :: world,img,pix
    real(REAL64), dimension(1) :: phi, theta
    integer, dimension(1) :: stat
    integer :: status

    pix(1,:) = [ x, y ]
    status = wcs_pixel2world(this % wcs,pix,img,phi,theta,world,stat)
    a = world(1,1)
    d = world(1,2)

  end subroutine pix2cel1

  subroutine cel2pixX(this,a,d,x,y)

    class(AstroWCS), intent(in) :: this
    real(REAL64), dimension(:), intent(in) :: a,d
    real(REAL64), dimension(:), intent(out) :: x,y
    real(REAL64), dimension(size(a),2) :: world,img,pix
    real(REAL64), dimension(size(a)) :: phi, theta
    integer, dimension(size(a)) :: stat
    integer :: status

    world(:,1) = a
    world(:,2) = d
    status = wcs_world2pixel(this % wcs,world,phi,theta,img,pix,stat)
    x = pix(:,1)
    y = pix(:,2)

  end subroutine cel2pixX

  subroutine cel2pix1(this,a,d,x,y)

    class(AstroWCS), intent(in) :: this
    real(REAL64), intent(in) :: a,d
    real(REAL64), intent(out) :: x,y
    real(REAL64), dimension(1,2) :: world,img,pix
    real(REAL64), dimension(1) :: phi, theta
    integer, dimension(1) :: stat
    integer :: status

    world(1,:) = [ a, d]
    status = wcs_world2pixel(this % wcs,world,phi,theta,img,pix,stat)
    x = pix(1,1)
    y = pix(1,2)

  end subroutine cel2pix1

  subroutine get_code(this,code)

    class(AstroWCS), intent(in out) :: this
    character(len=*), intent(out) :: code
    character(len=88), dimension(2) :: ctype
    character(len=88) :: a
    integer :: status

    status = wcs_get(this % wcs,WCS_CTYPE,ctype)
    if( status == 0 ) then
       a = ctype(1)
       code = a(6:8)  ! the 4-3 convention
    else
       code = ''
    end if

  end subroutine get_code

  subroutine get_affines(this,scale,rot,refl)
    class(AstroWCS), intent(in out) :: this
    real(REAL64), intent(out) :: scale,rot,refl
    real(REAL64), dimension(2,2) :: pc, mi, m
    real(REAL64), dimension(2) :: cdelt
    real(REAL64) :: c,s
    integer :: status

    status = wcs_get(this % wcs,WCS_PC,pc)
    if( status /= 0 ) then
       pc(1,:) = [ 1, 0 ]
       pc(2,:) = [ 0, 1 ]
       write(error_unit,*) 'Warning: PC undefined in get_affines().'
    end if

    status = wcs_get(this % wcs,WCS_CDELT,cdelt)
    if( status /= 0 ) then
       cdelt = 1
       write(error_unit,*) 'Warning: CDELT undefined in get_affines().'
    end if
    if( any(cdelt == 0) ) stop 'Error: CDELT zero in get_affines().'

    ! un-flip spherical to rectangular coordinates
    mi(1,:) = real([ -1, 0 ],REAL64)
    mi(2,:) = real([  0, 1 ],REAL64)
    m = matmul(mi,pc)

    ! the reflection (flip, mirror) is determined by diagonal elements
    ! having non-corresponding signs
    refl = sign(1.0_REAL64,m(1,1)*m(2,2))

    ! transform-out reflex
    mi(1,:) = [    refl, 0.0_REAL64 ]
    mi(2,:) = [ 0.0_REAL64, 1.0_REAL64 ]
    m = matmul(m,mi)

    c = (m(1,1) + m(2,2)) / 2
    s = (m(2,1) - m(1,2)) / 2

    scale = hypot(s,c) * (sum(abs(cdelt)) / size(cdelt))
    if( scale > 0 ) then
       scale = 1 / scale ! in pix per deg
    else
       stop 'Error in get_affines(): the scale is undefined.'
    end if

    if( abs(c) == 0 ) then
       rot = 0
    else
       rot = rad*(atan2(-m(1,2),m(1,1)) + atan2(m(2,1),m(2,2))) / 2
    end if

  end subroutine get_affines

  subroutine get_crval(this,crval)

    class(AstroWCS), intent(in out) :: this
    real(REAL64), dimension(:), intent(out) :: crval
    integer :: status

    status = wcs_get(this % wcs,WCS_CRVAL,crval)
    if( status /= 0 ) crval = 0

  end subroutine get_crval

  subroutine get_crpix(this,crpix)

    class(AstroWCS), intent(in out) :: this
    real(REAL64), dimension(:), intent(out) :: crpix
    integer :: status

    status = wcs_get(this % wcs,WCS_CRPIX,crpix)
    if( status /= 0 ) crpix = 0

  end subroutine get_crpix

  subroutine get_crder(this,crder)

    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN

    class(AstroWCS), intent(in out) :: this
    real(REAL64), dimension(:), intent(out) :: crder
    integer :: status

    status = wcs_get(this % wcs,WCS_CRDER,crder)
    if( status /= 0 ) crder = IEEE_VALUE(crder, IEEE_QUIET_NAN)

  end subroutine get_crder

  subroutine get_epoch(this,epoch)

    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN

    class(AstroWCS), intent(in out) :: this
    real(REAL64), intent(out) :: epoch
    integer :: status

    status = wcs_get(this % wcs,WCS_JEPOCH,epoch)
    if( status /= 0 .or. abs(epoch) >= UNDEFINED ) &
         epoch = IEEE_VALUE(epoch, IEEE_QUIET_NAN)

  end subroutine get_epoch

  subroutine get_radesys(this,radesys)

    class(AstroWCS), intent(in out) :: this
    character(len=*), intent(out) :: radesys
    character(len=88) :: a
    integer :: status

    status = wcs_get(this % wcs,WCS_RADESYS,a)
    if( status == 0 ) then
       radesys = a
    else
       radesys = ''
    end if

  end subroutine get_radesys

  subroutine print_header(this)

    class(AstroWCS), intent(in) :: this
    integer, parameter :: WCSHDO_all = 255
    integer, parameter :: LEN = 80
    character(len=:), allocatable :: header
    integer :: status, nkeyrec, i,j

    status = wcs_header_output(WCSHDO_all,this % wcs,nkeyrec,header)
    if( status == 0 ) then
       do i = 1, nkeyrec
          j = LEN*(i - 1)
          write(error_unit,*) header(j+1:j+LEN)
       end do
    end if
    deallocate(header)

  end subroutine print_header

  subroutine dealloc(this)
    type(AstroWCS), intent(in out) :: this
    integer :: status

    call wcs_free(this % wcs)
    status = wcserr_enable(0)

  end subroutine dealloc


end module astrotrafo
