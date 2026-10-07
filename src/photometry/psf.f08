!
!  fitspphot
!
!  Copyright © 2013 F.Hroch (hroch@physics.muni.cz)
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


module psfengine

  use titsio
  use iso_fortran_env

  implicit none

  integer, parameter, private :: rp = selected_real_kind(15)

contains

  subroutine psfmodel(xcens,ycens,asky,acts,frame,ring,lobad,hibad,phpadu, &
       verbose,plog,cts,cts_err,sky,sky_err,rframe,status)

    real(rp), dimension(:), intent(in) :: xcens,ycens,asky,acts
    real(rp), dimension(:,:), intent(in) :: frame
    real(rp), intent(in) :: ring,lobad,hibad,phpadu
    logical, intent(in) :: verbose,plog
    real(rp), dimension(:), intent(out) :: cts,cts_err,sky,sky_err
    real(rp), dimension(:,:), intent(out) :: rframe
    integer, intent(in out) :: status

    logical, dimension(:,:), allocatable :: psfmask
    real(rp), dimension(:,:), allocatable :: psf,tpsf,convol
    real(rp), dimension(:), allocatable :: r
    integer, dimension(:), allocatable :: icen,jcen
    integer :: npsf,nring,nstars,nx,ny,i,j,n,k,l
    real(rp) :: s,x,y
    type(fitsfiles) :: fits


    nring = nint(ring)
    npsf = 2*nring + 1
    nstars = size(xcens)
    allocate(psf(npsf,npsf),tpsf(npsf,npsf),convol(npsf,npsf), &
         r(nstars),icen(nstars),jcen(nstars))

    ! central pixels
    icen = nint(xcens)
    jcen = nint(ycens)
!    do i = 1,size(xcens)
!       icen(i) = nint(xcens(i))
!       jcen(i) = nint(ycens(i))
!    end do

    do n = 1,nstars
       if( abs(icen(n) - 840) < ring .and. abs(jcen(n) - 866) < ring ) then

          tpsf = frame(icen(n)-nring:icen(n)+nring,jcen(n)-nring:jcen(n)+nring)-asky(n)
          do i = 1,npsf
             do j = 1,npsf

                x = i + (xcens(n) - icen(n))
                y = j + (ycens(n) - jcen(n))
!                psf(i,j) =  bicubic(x,y,tpsf)
                psf(i,j) =  bilinear(x,y,tpsf)

             end do
          end do

!          psf = frame(icen(n)-nring:icen(n)+nring,jcen(n)-nring:jcen(n)+nring) - asky(n)
          psf = psf / acts(n)
       end if
    end do

!    s = 0
!    do i = 1,npsf
!       do j = 1,npsf
!          s = s + max(psf(i,j),0.0_rp)
!       end do
!    end do
!    psf = psf / s
!    psf /



    ! save PSF
    call fits_create_file(fits,'!psf.fits',status)
    call fits_insert_img(fits,-32,2,[npsf,npsf],status)
    call fits_write_image(fits,0,psf,status)
    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)



    ! synthetic stars
    nx = size(frame,1)
    ny = size(frame,2)
    allocate(psfmask(nx,ny))
    psfmask = .false.
    rframe = 1
    do n = 1,nstars

       do i = 1, npsf
          do j = 1,npsf

             l = icen(n)-nring -1 + i
             k = jcen(n)-nring -1 + j
             if( 1 <= l .and. l <= nx .and. 1 <= k .and. k <= ny ) then
                x = i + (xcens(n) - icen(n))
                y = j + (ycens(n) - jcen(n))
!                tpsf(,j) = acts(n)*bicubic(x,y,psf)
                tpsf(i,j) = acts(n)*bilinear(x,y,psf)

!                tpsf(i,j) = bicubic(x,y,psf)
!                rframe(l,k) = acts(n)*s + rframe(l,k)
!                rframe(l,k) = (acts(n)*psf(i,j) + 0*asky(n)) + rframe(l,k)
!                rframe(l,k) = acts(n)*psf(i,j) * rframe(l,k)
!                psfmask(l,k) = .true.
             end if
          end do
       end do

       do i = 1, npsf
          do j = 1,npsf

             !rframe(i,j) = sum(tpsf(i,:)*rframe(:,j))

             l = icen(n)-nring -1 + i
             k = jcen(n)-nring -1 + j
             if( npsf < l .and. l < nx-npsf .and. npsf < k .and. k < ny-npsf ) then
                !write(*,*) n,l,k,i,j
                convol(i,j) = tpsf(i,j) + rframe(l,k)
!                convol(i,j) = sum(tpsf(i,:)*rframe(l:l+npsf-1,k))
!                call convolution(tpsf,rframe(l:l+npsf-1,k),convol)
!                rframe(l,k) = sum(tpsf(i,:)*rframe(l:l+npsf-1,k))
             end if
          end do
       end do


       l = icen(n)-nring -1
       k = jcen(n)-nring -1
       if( npsf < l .and. l < nx-npsf .and. npsf < k .and. k < ny-npsf ) then
!          call convolution(tpsf,rframe(l:l+npsf-1,k:k+npsf-1),convol)
          rframe(l:l+npsf-1,k:k+npsf-1) = convol
       end if
    end do

!    return

!    where( psfmask )
       rframe = frame - rframe
!    elsewhere
!       rframe = frame
!    end where

    deallocate(psf)

  end subroutine psfmodel

  subroutine convolution(a,b,c)

    real(rp), dimension(:,:), intent(in) :: a,b
    real(rp), dimension(:,:), intent(out) :: c

    integer :: i,j,k,l,m,n,sn
    real(rp) :: s

    n = size(a,1)
    m = size(a,2)

    do i = 1,n
      do j = 1,m
!    do i = -n/2,n/2
!       do j = -m/2,m/2

          s = 0
          sn = 0

          do l = -n/2+1,n/2
             do k = -m/2+1,m/2
!               write(*,*)  n/2+i-l-1,m/2+j-k-1,l,k,i,j,m,n
 !               s = s + a(n/2+i-l,m+j-k)*b(l,k)
               if( 0 < i-l .and. i-l < n .and. 0 < j-k .and. j-k < m .and. 0 < l &
                    .and. l < n .and. 0 < k .and. k < m ) then
                s = s + a(i-l,j-k)*b(l,k)
                sn = sn + 1
             end if
             end do
          end do


!!$          do l = 1,0*i
!!$             do k = 1,0*j
!!$!                write(*,*) (i-l)+1,(j-k)+1,l,k,i,j,m,n
!!$                s = s + a((i-l)+1,(j-k)+1)*b(l,k)
!!$                sn = sn + 1
!!$             end do
!!$          end do
!!$          do l = i+1,n
!!$
!!$             do k = j+1,n
!!$!                write(*,*) n+(i-l),m+(j-k),l,k,i,j,m,n
!!$                s = s + a(n+(i-l),m+(j-k))*b(l,k)
!!$                sn = sn + 1
!!$             end do
!!$
!!$          end do
          c(i,j) = s/sn

       end do
    end do

  end subroutine convolution


  function bilinear(x,y,array) result(fun)

    real(rp), intent(in) :: x,y
    real(rp), dimension(:,:), intent(in) :: array
    real(rp) :: fun
    integer :: i,j
    real(rp) :: u,v

    i = int(x)
    j = int(y)

    if( (lbound(array,1) <= i .and. i < ubound(array,1)) .and. &
        (lbound(array,2) <= j .and. j < ubound(array,2)) ) then

       u = (x - i)/(i+1 - i)
       v = (y - j)/(j+1 - j)

       fun = (1.0 - u)*(1.0 - v)*array(i,j) + u*v*array(i+1,j+1) &
            + u*(1.0 - v)*array(i+1,j) + (1.0 - u)*v*array(i,j+1)


    else
       fun = 0.0
    end if

  end function bilinear

  function bicubic(x,y,array)

    implicit none

    real(rp), intent(in) :: x,y
    real(rp), dimension(:,:), intent(in) :: array
    real(rp) :: bicubic

    integer :: i,j,nx,ny
    real(rp) :: d, z, z1, z2, x1l, x1u, x2l, x2u
    real(rp), dimension(4) :: f, f1, f2, f12

    nx = size(array,1)
    ny = size(array,2)
    i = int(x)
    j = int(y)

    d = 1.0  ! coordinate step

    if( 1 < i .and. i < nx - 1 .and. 1 < j .and. j < ny - 1 ) then

       f(1) = array(i,j)
       f(2) = array(i+1,j)
       f(3) = array(i+1,j+1)
       f(4) = array(i,j+1)

       f1(1) = (array(i+1,j)   - array(i-1,j))  /(2.0*d)
       f1(2) = (array(i+2,j)   - array(i,j))    /(2.0*d)
       f1(3) = (array(i+2,j+1) - array(i,j+1))  /(2.0*d)
       f1(4) = (array(i+1,j+1) - array(i-1,j+1))/(2.0*d)

       f2(1) = (array(i,j+1)   - array(i,j-1))  /(2.0*d)
       f2(2) = (array(i+1,j+1) - array(i+1,j-1))/(2.0*d)
       f2(3) = (array(i+1,j+2) - array(i+1,j))  /(2.0*d)
       f2(4) = (array(i,j+2)   - array(i,j))    /(2.0*d)

       f12(1) = (array(i+1,j+1)-array(i+1,j-1)-array(i-1,j+1)+array(i-1,j-1))/(2.0*d)**2
       f12(2) = (array(i+2,j+1)-array(i+2,j-1)-array(i,j+1)+array(i,j-1))/(2.0*d)**2
       f12(3) = (array(i+2,j+2)-array(i+2,j)-array(i,j+2)+array(i,j))/(2.0*d)**2
       f12(4) = (array(i+1,j+2)-array(i+1,j)-array(i-1,j+2)+array(i-1,j))/(2.0*d)**2

       x1l = i
       x1u = x1l + d
       x2l = j
       x2u = x2l + d

       call bcuint(f,f1,f2,f12,x1l,x1u,x2l,x2u,x,y,z,z1,z2)
       bicubic = z

    else

       bicubic = 0

    end if

  end function bicubic

  subroutine bcucof(y,y1,y2,y12,d1,d2,c)

    implicit none

    real(rp), intent(in) :: d1,d2
    real(rp), dimension(4), intent(in) :: y,y1,y2,y12
    real(rp), dimension(4,4), intent(out) :: c

    real(rp), dimension(16) :: x

    real(rp), parameter, dimension(16,16) :: wt = reshape( (/ &
     1,0,-3,2,0,0,0,0,-3,0,9,-6,2,0,-6,4,0,0,0,0,0,0,0,0,3,0,-9,6,-2,0,6,-4, &
     0,0,0,0,0,0,0,0,0,0,9,-6,0,0,-6,4,0,0,3,-2,0,0,0,0,0,0,-9,6,0,0,6,-4, &
     0,0,0,0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,0,0,0,0,0,0,0,0,-1,0,3,-2,1,0,-3,2,&
     0,0,0,0,0,0,0,0,0,0,-3,2,0,0,3,-2,0,0,0,0,0,0,3,-2,0,0,-6,4,0,0,3,-2,&
     0,1,-2,1,0,0,0,0,0,-3,6,-3,0,2,-4,2,0,0,0,0,0,0,0,0,0,3,-6,3,0,-2,4,-2,&
     0,0,0,0,0,0,0,0,0,0,-3,3,0,0,2,-2,0,0,-1,1,0,0,0,0,0,0,3,-3,0,0,-2,2,&
     0,0,0,0,0,1,-2,1,0,-2,4,-2,0,1,-2,1,0,0,0,0,0,0,0,0,0,-1,2,-1,0,1,-2,1, &
     0,0,0,0,0,0,0,0,0,0,1,-1,0,0,-1,1,0,0,0,0,0,0,-1,1,0,0,2,-2,0,0,-1,1 /),&
     (/16,16/) )

    x(1:4) = y
    x(5:8) = y1*d1
    x(9:12) = y2*d2
    x(13:16) = y12*d1*d2
    x = matmul(wt,x)
    c = reshape(x,(/4,4/),order=(/2,1/))

  end subroutine bcucof

  subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)

    implicit none

    real(rp), dimension(4), intent(in) :: y,y1,y2,y12
    real(rp), intent(in) :: x1l,x1u,x2l,x2u,x1,x2
    real(rp), intent(out) :: ansy,ansy1,ansy2

    integer :: i
    real(rp) :: t,u
    real(rp), dimension(4,4) :: c

    call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)

    if( x1u == x1l .or. x2u == x2l ) &
         stop 'bcuint: problem with input values - boundary pair equal?'

    t = (x1 - x1l)/(x1u - x1l)
    u = (x2 - x2l)/(x2u - x2l)
    ansy = 0.0
    ansy1 = 0.0
    ansy2 = 0.0
    do i = 4,1,-1
       ansy = t*ansy + ((c(i,4)*u + c(i,3))*u + c(i,2))*u + c(i,1)
       ansy1 = u*ansy1 + (3.0*c(4,i)*t + 2.0*c(3,i))*t + c(2,i)
       ansy2 = t*ansy2 + (3.0*c(i,4)*u + 2.0*c(i,3))*u + c(i,2)
    end do
    ansy1 = ansy1/(x1u - x1l)
    ansy2 = ansy2/(x2u - x2l)

  end subroutine bcuint

end module psfengine
