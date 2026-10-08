!
!  Generate artifical sky frame
!
!
!  Copyright © 2016-2026 F.Hroch (hroch@physics.muni.cz)
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


module generator

  use iso_fortran_env
  implicit none

  private :: stardist

contains

  subroutine spreadstars(spsf,xcoo,ycoo,flux,crpix,sky)

    ! Adds the noiseless (mean) star images to the frame,
    ! the noise is added later on the whole frame.

    use spray

    class(SprayType), intent(in) :: spsf
    real(REAL64), dimension(:), intent(in) :: xcoo,ycoo,flux,crpix
    real(REAL64), dimension(:,:), intent(in out) :: sky
    integer, dimension(2) :: naxes
    real(REAL64), dimension(:,:), allocatable :: psf
    real(REAL64) :: f,dx,dy
    integer :: nw,i,j,n,x,y,k,l

    naxes = shape(sky)

    nw = spsf%nbox
    allocate(psf(-nw:nw,-nw:nw))

    do n = 1, size(xcoo)
       x = nint(xcoo(n))
       y = nint(ycoo(n))
       dx = x - xcoo(n)
       dy = y - ycoo(n)

       call spsf%Pixelize(xcoo(n)-crpix(1),ycoo(n)-crpix(2),dx,dy,psf,nw)
       psf = flux(n)*psf

       do i = x - nw, x + nw
          k = i - x
          do j = y - nw, y + nw
             l = j - y
             if(  1 <= i .and. i <= naxes(1) .and. &
                  1 <= j .and. j <= naxes(2) ) then
                f = psf(k,l)
                if( f > 0 ) sky(i,j) = sky(i,j) + f
             end if
          end do
       end do

    end do

    deallocate(psf)

  end subroutine spreadstars

  subroutine rain(prjcode,alpha,delta,mag,scale,angle,crval,crpix,phsys, &
       filter,xcoo,ycoo,flux)

    use astrotrafo
    use phsysfits
    use photoconv

    real(REAL64), dimension(:), intent(in) :: alpha,delta,mag
    character(len=*), intent(in) :: prjcode,filter
    real(REAL64), dimension(:), intent(in) :: crpix, crval
    real(REAL64), intent(in) :: scale,angle
    real(REAL64), dimension(:), allocatable, intent(out) :: xcoo,ycoo,flux
    type(type_phsys), intent(in) :: phsys
    type(AstroWCS) :: wcs
    real(REAL64), dimension(:), allocatable :: dmag, dflux
    real(REAL64), parameter :: refl = 1
    integer :: n

    call wcs % init(prjcode,crval,crpix,scale,angle,refl)

    n = size(alpha)
    allocate(xcoo(n),ycoo(n),flux(n),dmag(n),dflux(n))
    call wcs % cel2pix(alpha,delta,xcoo,ycoo)

    dmag = 0
    call phsysmagph1(phsys,filter,mag,dmag,flux,dflux)

    deallocate(dmag, dflux)

  end subroutine rain


  subroutine catstars(cat,labels,prjcode,scale,angle,crval,crpix,phsys, &
       filter,xcoo,ycoo,flux)

    use catio
    use phsysfits

    character(len=*), intent(in) :: cat, filter, prjcode
    type(type_phsys), intent(in) :: phsys
    character(len=*), dimension(:), intent(in) :: labels
    real(REAL64), dimension(:), intent(in) :: crpix, crval
    real(REAL64), intent(in) :: scale,angle
    real(REAL64), dimension(:), allocatable, intent(out) :: xcoo,ycoo,flux
    real(REAL64), dimension(:), allocatable :: alpha,delta,mag
    character(len=80) :: catid

    call catalogue(cat,labels,alpha,delta,mag,catid)

    call rain(prjcode,alpha,delta,mag,scale,angle,crval,crpix,phsys, &
         filter,xcoo,ycoo,flux)

    deallocate(alpha,delta,mag)

  end subroutine catstars


  subroutine genstars(maglim,naxes,xcoo,ycoo,flux)

    integer, intent(in) :: maglim
    integer, dimension(:), intent(in) :: naxes
    real(REAL64), dimension(:), allocatable, intent(out) :: xcoo,ycoo,flux
    integer, dimension(0:maglim) :: nsq
    real :: x,y
    integer :: mag,n,i

    call stardist(maglim,nsq)

    n = sum(nsq)
    allocate(xcoo(n),ycoo(n),flux(n))

    n = 0

    do mag = 0, maglim

       do i = 1,nsq(mag)

          n = n + 1

          ! coordinates
          call random_number(x)
          call random_number(y)
          xcoo(n) = naxes(1) * x
          ycoo(n) = naxes(2) * y

          ! flux
          call random_number(x)
          flux(n) = 5e4*10**(-0.4*(mag + x - maglim + 1))

       end do

    end do

  end subroutine genstars


  pure subroutine stardist(maglim,nsq)

    integer, intent(in) :: maglim
    integer, dimension(0:), intent(out) :: nsq
    integer :: mag

    ! mean magnitude distribution by Gaia DR1:
    ! https://www.cosmos.esa.int/web/gaia/dr1
    ! approximated as sigma(m) = 1e-4*exp(1.1*mag)

    do mag = 0, maglim
       nsq(mag) = int(1e-4*exp(1.1*mag)) ! per square degree
    end do

  end subroutine stardist


  subroutine background(b,grad,sky)

    ! The noiseless (mean) sky background in electrons per pixel:
    !
    !   B(x,y) = b * 10**(-0.4*(grad(1)*(x - x0) + grad(2)*(y - y0)))
    !
    ! b is the level at the centre (x0,y0), grad is the gradient
    ! of the sky brightness in mag/arcsec2/pix (positive values
    ! makes the sky fainter in the direction).

    real(REAL64), intent(in) :: b
    real(REAL64), dimension(:), intent(in) :: grad
    real(REAL64), dimension(:,:), intent(out) :: sky
    integer :: i,j,ic,jc
    real(REAL64) :: dx,dy

    ic = size(sky,1) / 2
    jc = size(sky,2) / 2
    do i = 1,size(sky,1)
       dx = (i - ic)*grad(1)
       do j = 1,size(sky,2)
          dy = (j - jc)*grad(2)
          sky(i,j) = b * 10**(-0.4_REAL64*(dx + dy))
       end do
    end do

  end subroutine background


  subroutine detector_noise(sky,readnoise)

    ! The detector noise: the photon (Poisson) noise of the mean
    ! number of electrons (all sources: sky, stars, dark current)
    ! in every pixel, and the Gaussian read-out noise N(0,readnoise)
    ! added afterwards, both in electrons.

    use noise

    real(REAL64), dimension(:,:), intent(in out) :: sky
    real(REAL64), intent(in) :: readnoise
    integer :: i,j

    do j = 1,size(sky,2)
       do i = 1,size(sky,1)
          sky(i,j) = pnoise(max(sky(i,j),0.0_REAL64))
       end do
    end do

    if( readnoise > 0 ) then
       do j = 1,size(sky,2)
          do i = 1,size(sky,1)
             sky(i,j) = sky(i,j) + gnoise(0.0_REAL64,readnoise)
          end do
       end do
    end if

  end subroutine detector_noise

end module generator
