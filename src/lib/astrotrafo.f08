!
! Astrometry related transformations
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

module astrotrafo

  use iso_fortran_env
  implicit none

  logical, parameter, private :: debug = .false.
  real(REAL64), parameter, private :: rad = 57.295779513082322865_REAL64

  ! Transformation and projection:
  !
  ! type ... type of projection (GNOMONIC,..)
  !      ... mostly affine transformation ...
!!  ! scale... scale [deg/pixel]
  ! scale... scale [pixel/deg]
  ! rot  ... rotation angle [deg]
  ! acen ... center projection in alpha [deg]
  ! dcen ... center projection in delta [deg]
  ! xcen ... horizontal center projection [pixel]
  ! ycen ... vertical center projection [pixel]
  ! refl ... -1 for reflected frame
!!  ! err  ... typical statistical error of acen, dcen [deg]

  type AstrometryGeometry
     real(REAL64), private :: scale,rot,refl
     real(REAL64), dimension(2,2) :: cd, icd
   contains
     procedure :: init => initgeom
     procedure :: pix2uv => touv
     procedure :: uv2pix => topix
  end type AstrometryGeometry

  type AstrometryProjection
     character(len=80), private :: type = 'GNOMONIC'
     real(REAL64), private :: acen, dcen
   contains
     procedure :: init => initproj
     procedure :: proj => cell2uv
     procedure :: invproj => uv2cell
  end type AstrometryProjection

  type AstroTrafoProj
     character(len=80) :: type
!     real(REAL64) :: acen,dcen,xcen,ycen,ucen,vcen,scale,rot,err,refl
     real(REAL64) :: acen,dcen,xcen,ycen,scale,rot,err,refl
     real(REAL64), dimension(2,2) :: mrot,smat,mat,mat1
!     real(REAL64), dimension(2) :: xy0,uv0
  end type AstroTrafoProj

  interface trafo
     module procedure trafo_point, trafo_array
  end interface trafo

  interface invtrafo
     module procedure invtrafo_point, invtrafo_array
  end interface invtrafo

  interface affine
     module procedure affine_point, affine_array
  end interface affine

  interface invaffine
     module procedure invaffine_point, invaffine_array
  end interface invaffine

  interface proj
     module procedure proj_point, proj_array
  end interface proj

  interface invproj
     module procedure invproj_point, invproj_array
  end interface invproj

  private :: trafo_point, trafo_array,invtrafo_point, invtrafo_array, &
       affine_point, affine_array, invaffine_point, invaffine_array, &
       proj_point, proj_array, invproj_point, invproj_array, &
       gproj, invgproj

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

    ! .. and the inverse
    b(1,:) = [ c,-s]
    b(2,:) = [ s, c]
    t%icd = matmul(a,b) / scale

  end subroutine initgeom

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


  subroutine initproj(this,type,acen,dcen)

    use projections

    class(AstrometryProjection), intent(out) :: this
    character(len=*), intent(in) :: type
    real(REAL64), intent(in) :: acen,dcen

    this % type = type
    this % acen = acen
    this % dcen = dcen

  end subroutine initproj

  pure elemental subroutine cell2uv(t,a,d,u,v)

    use projections

    class(AstrometryProjection), intent(in) :: t
    real(REAL64), intent(in) :: a,d
    real(REAL64), intent(out) :: u,v

    if( t%type == "GNOMONIC" ) then
       call gnomon(a,d,t%acen,t%dcen,u,v)
    else !if( t%type == "IDENTITY" ) then
       call carree(a,d,t%acen,t%dcen,u,v)
    end if

  end subroutine cell2uv

  pure elemental subroutine uv2cell(t,u,v,a,d)

    use projections

    class(AstrometryProjection), intent(in) :: t
    real(REAL64), intent(out) :: a,d
    real(REAL64), intent(in) :: u,v

    if( t%type == "GNOMONIC" ) then
       call invgnomon(u,v,t%acen,t%dcen,a,d)
    else !if( t%type == "IDENTITY" ) then
       call invcarree(u,v,t%acen,t%dcen,a,d)
    end if

  end subroutine uv2cell



! --- init ---

!  subroutine trafo_init(t,type,acen,dcen,xcen,ycen,ucen,vcen,scale,rot,refl,err)
  subroutine trafo_init(t,type,acen,dcen,xcen,ycen,scale,rot,refl,err)

    type(AstroTrafoProj), intent(out) :: t
    character(len=*), intent(in), optional :: type
!    real(REAL64), intent(in), optional :: acen,dcen,xcen,ycen,ucen,vcen,scale,rot,err,refl
    real(REAL64), intent(in), optional :: acen,dcen,xcen,ycen,scale,rot,err,refl

    t%type = 'GNOMONIC'
    t%acen = 0.0_REAL64
    t%dcen = 0.0_REAL64
    t%xcen = 0.0_REAL64
    t%ycen = 0.0_REAL64
!    t%ucen = 0.0_REAL64
!    t%vcen = 0.0_REAL64
    t%scale= 1.0_REAL64
    t%rot  = 0.0_REAL64
    t%refl = 1.0_REAL64
    t%err  = epsilon(t%err)

    if( present(type) ) t%type = type
    if( present(acen) ) t%acen = acen
    if( present(dcen) ) t%dcen = dcen
    if( present(xcen) ) t%xcen = xcen
    if( present(ycen) ) t%ycen = ycen
!    if( present(ucen) ) t%ucen = ucen
!    if( present(vcen) ) t%vcen = vcen
    if( present(scale)) t%scale= scale
    if( present(rot)  ) t%rot  = rot
    if( present(refl) ) t%refl = refl
    if( present(err)  ) t%err  = err

    call trafo_refresh(t)

  end subroutine trafo_init


  subroutine trafo_refresh(t)

    type(AstroTrafoProj), intent(inout) :: t
    real(REAL64), dimension(2,2) :: m,sm
    real(REAL64) :: f,c,s!,det

    ! the affine projection matrix
    f = t%rot / rad
    c = cos(f)
    s = sin(f)

    t%smat(1,:) = [ t%scale * t%refl, 0.0_REAL64 ]
    t%smat(2,:) = [ 0.0_REAL64,          t%scale ]
    t%mrot(1,:) = [ c,-s ]
    t%mrot(2,:) = [ s, c ]
    t%mat = matmul(t%mrot,t%smat)

    ! .. and the inverse
    sm(1,:) = [ t%refl / t%scale,     0.0_REAL64 ]
    sm(2,:) = [          0.0_REAL64, 1 / t%scale ]
    m(1,:) = [ c, s]
    m(2,:) = [-s, c]
    t%mat1 = matmul(m,sm)

    if( debug ) then
       block
         real(REAL64), dimension(2,2) :: eye
         eye = matmul(t%mat,t%mat1)
         write(error_unit,*) eye(1,:)
         write(error_unit,*) eye(2,:)
         eye = matmul(t%mat1,t%mat)
         write(error_unit,*) eye(1,:)
         write(error_unit,*) eye(2,:)

         write(error_unit,*) t%mat(1,:)
         write(error_unit,*) t%mat(2,:)
         write(error_unit,*) t%mat1(1,:)
         write(error_unit,*) t%mat1(2,:)
       end block
    end if

!    t%xy0 = (/t%xcen, t%ycen/)
!    t%uv0 = (/t%ucen, t%vcen/)

  end subroutine trafo_refresh

! --- trafo ---

  subroutine trafo_point(t, a, d, x,y)

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), intent(in) :: a,d
    real(REAL64), intent(out) :: x,y
    real(REAL64), dimension(1) :: xx,yy

    call trafo(t, (/a/), (/d/), xx,yy)
    x = xx(1)
    y = yy(1)

  end subroutine trafo_point

  subroutine trafo_array(t, a, d, x,y)


    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), dimension(:), intent(in) :: a,d
    real(REAL64), dimension(:), intent(out) :: x,y
    real(REAL64), dimension(:), allocatable :: u,v
    integer :: n

    n = size(a)
    allocate(u(n),v(n))

    call proj(t,a,d,u,v)
    call affine(t,u,v,x,y)

    deallocate(u,v)

  end subroutine trafo_array

! --- invtrafo ---

  subroutine invtrafo_point(t,x,y,a,d)

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), intent(in) :: x,y
    real(REAL64), intent(out) :: a,d
    real(REAL64), dimension(1) :: aa,dd

    call invtrafo(t,(/x/),(/y/),aa,dd)
    a = aa(1)
    d = dd(1)

  end subroutine invtrafo_point

  subroutine invtrafo_array(t,x,y,a,d)

    use projections

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), dimension(:), intent(in) :: x,y
    real(REAL64), dimension(:), intent(out) :: a,d
    real(REAL64), dimension(:), allocatable :: u,v
    integer :: n

    n = size(a)
    allocate(u(n),v(n))

    call invaffine(t,x,y,u,v)
    call invproj(t,u,v,a,d)

    deallocate(u,v)

  end subroutine invtrafo_array


! --- proj ---

  subroutine proj_array(t,a,d,u,v)

    use projections

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), dimension(:), intent(in) :: a,d
    real(REAL64), dimension(:), intent(out) :: u,v

    if( t%type == "GNOMONIC" ) then
       call gproj(gnomond,t%acen,t%dcen,a,d,u,v)
    else if( t%type == " " ) then
       call gproj(identity,t%acen,t%dcen,a,d,u,v)
    end if

  end subroutine proj_array

  subroutine proj_point(t,a,d,u,v)

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), intent(in) :: a,d
    real(REAL64), intent(out) :: u,v
    real(REAL64), dimension(1) :: uu,vv

    call proj(t,(/a/),(/d/),uu,vv)
    u = uu(1)
    v = vv(1)

  end subroutine proj_point

! --- invproj ---

  subroutine invproj_array(t,u,v,a,d)

    use projections

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), dimension(:), intent(in) :: u,v
    real(REAL64), dimension(:), intent(out) :: a,d

    if( t%type == "GNOMONIC" ) then
       call invgproj(invgnomond,t%acen,t%dcen,u,v,a,d)
    else if( t%type == " " ) then
       call invgproj(invidentity,t%acen,t%dcen,u,v,a,d)
    end if

  end subroutine invproj_array

  subroutine invproj_point(t,u,v,a,d)

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), intent(in) :: u,v
    real(REAL64), intent(out) :: a,d
    real(REAL64), dimension(1) :: aa,dd

    call invproj(t,(/u/),(/v/),aa,dd)
    a = aa(1)
    d = dd(1)

  end subroutine invproj_point

! --- invaffine ---


  subroutine invaffine_point(t,x,y,u,v,g,h)

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), intent(in) :: x,y
    real(REAL64), intent(out) :: u,v
    real(REAL64), intent(out), optional :: g,h

    real(REAL64), dimension(1) :: uu,vv,gg,hh

    if( present(g) .and. present(h) ) then
       call invaffine(t,(/x/),(/y/),uu,vv,gg,hh)
       u = uu(1)
       v = vv(1)
       g = gg(1)
       h = hh(1)
    else
       call invaffine(t,(/x/),(/y/),uu,vv)
       u = uu(1)
       v = vv(1)
    end if

  end subroutine invaffine_point

  subroutine invaffine_array(t,x,y,u,v,g,h)

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), dimension(:), intent(in) :: x,y
    real(REAL64), dimension(:), intent(out) :: u,v
    real(REAL64), dimension(:), intent(out), optional :: g,h

    real(REAL64), dimension(2) :: r,s
    integer :: i

    do i = 1,size(x)
       r = [ x(i) - t%xcen, y(i) - t%ycen ] !- [t%xcen, t%ycen]  !t%xy0
       s = matmul(t%mat,r) !+ t%uv0
       u(i) = s(1)
       v(i) = s(2)

       if( present(g) .and. present(h) ) then
          s = matmul(t%mrot,r)
          g(i) = s(1)
          h(i) = s(2)
       end if

    end do

  end subroutine invaffine_array

! --- affine ---

  subroutine affine_point(t,u,v,x,y)

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), intent(in) :: u,v
    real(REAL64), intent(out) :: x,y

    real(REAL64), dimension(1) :: xx,yy

    call affine(t,(/u/),(/v/),xx,yy)
    x = xx(1)
    y = yy(1)

  end subroutine affine_point

  subroutine affine_array(t,u,v,x,y)

    type(AstroTrafoProj), intent(in) :: t
    real(REAL64), dimension(:), intent(in) :: u,v
    real(REAL64), dimension(:), intent(out) :: x,y

    real(REAL64), dimension(2) :: r!,s
    integer :: i

    do i = 1,size(u)
!       s = [ u(i),v(i) ] - [ t%ucen, t%vcen ]
!       r = matmul(t%mat1,s) + [ t%xcen, t%ycen ]
       r = matmul(t%mat1,[u(i),v(i)]) + [ t%xcen, t%ycen ]
       x(i) = r(1)
       y(i) = r(2)
    end do

  end subroutine affine_array


! --- projections ---

  subroutine gproj(proj, acen,dcen, a, d, u,v)

    interface
       subroutine proj(a,d,a0,d0,u,v)
         implicit none
         integer, parameter :: dp = selected_real_kind(15)
         real(dp),intent(in) :: a,d,a0,d0
         real(dp),intent(out) :: u,v
       end subroutine proj
    end interface

    real(REAL64), intent(in) :: acen,dcen
    real(REAL64), dimension(:),intent(in) :: a,d
    real(REAL64), dimension(:),intent(out) :: u,v
    integer :: i

    do i = 1, size(a)
       call proj(a(i),d(i),acen,dcen,u(i),v(i))
    end do

  end subroutine gproj

  subroutine invgproj(invproj, acen,dcen, u, v, a,d)

    interface
       subroutine invproj(a,d,a0,d0,u,v)
         implicit none
         integer, parameter :: dp = selected_real_kind(15)
         real(dp),intent(in) :: a,d,a0,d0
         real(dp),intent(out) :: u,v
       end subroutine invproj
    end interface

    real(REAL64), intent(in) :: acen,dcen
    real(REAL64), dimension(:),intent(in) :: u,v
    real(REAL64), dimension(:),intent(out) :: a,d
    integer :: i

    do i = 1, size(a)
       call invproj(u(i),v(i),acen,dcen,a(i),d(i))
    end do

  end subroutine invgproj


end module astrotrafo
