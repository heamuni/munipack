!
!  Spreaded profiles
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

module spray

  use iso_fortran_env

  implicit none
  private

  integer, parameter, private :: defzoom = 2**4
  real(REAL64), parameter, private :: rad = 45 / atan(1.0_REAL64)
  real(REAL64), parameter, private :: sqrt2 = sqrt(2.0_REAL64)

  type, public ::SprayType
     integer :: zoom = defzoom
     integer :: nbox = 0
     real(REAL64), dimension(:,:), allocatable :: psf
   contains
     procedure :: SprayMe
     procedure :: Pixelize => PixelizePsf
     procedure :: gauss, moffat, seeing, cseeing, fseeing
  end type SprayType

  type, extends(SprayType), public :: FileSprayType
     real(REAL64) :: scale = 0
     real(REAL64), dimension(:,:,:), allocatable :: mpsf
     real(REAL64), dimension(:), allocatable :: rcoo, acoo
   contains
     procedure :: SprayCan
     procedure :: Pixelize => PixelizeFile
  end type FileSprayType

  public :: SprayMe, SprayCan, PixelizePsf, PixelizeFile

contains

  subroutine SprayMe(this,profile,spread,hwhm,eccentricity,inclination, &
       airy,beta)

    class(SprayType) :: this
    character(len=*), intent(in) :: profile, spread
    real(REAL64), intent(in) :: hwhm, airy, beta, eccentricity,inclination
    integer :: n,m,n2

    if( profile == 'GAUSS' ) then

       this%nbox = nint(7*max(hwhm,1.0))
       n = this%zoom * this%nbox
       allocate(this%psf(-n:n,-n:n))
       call this%gauss(hwhm,eccentricity,inclination,n,this%psf)

    else if( profile == 'MOFFAT' ) then

       this%nbox = nint(30*max(hwhm,1.0))
       n = this%zoom * this%nbox
       allocate(this%psf(-n:n,-n:n))
       call this%moffat(hwhm,beta,eccentricity,inclination,n,this%psf)

    else if( profile == 'SEEING' ) then

       if( spread == 'RANDOM' ) then

          this%nbox = nint(13*max(hwhm,1.0)*max(airy,5.0))
          n = this%nbox*this%zoom
          allocate(this%psf(-n:n,-n:n))
          call this%seeing(hwhm,airy,eccentricity,inclination,n,this%psf)

       else if( spread == 'TEST' ) then

          this%nbox = nint(13*max(hwhm,1.0)*max(airy,5.0))
          n = this%nbox*this%zoom
          allocate(this%psf(-n:n,-n:n))
          call this%cseeing(hwhm,airy,eccentricity,inclination,n,this%psf)

       else ! if (spread == 'FFT' ) then

          n = nint(13*max(hwhm,1.0)*max(airy,5.0))
          m = int(log(real(n*this%zoom))/log(2.0) + 1)  ! power o 2 for FFT
          n = 2**m
          this%nbox = n / this%zoom
          n = this%nbox * this%zoom
          n2 = n / 2
          allocate(this%psf(-n2+1:n2-1,-n2+1:n2-1))
          call this%fseeing(hwhm,airy,eccentricity,inclination,n,m,n2,this%psf)

       end if

    else
       write(error_unit,*) "Error: unknown PSF `",trim(profile),"'."
       stop 'Error: unknown PSF.'
    end if

    ! normalisation
    this%psf = this%psf * this%zoom**2 / sum(this%psf)

  end subroutine SprayMe


  subroutine seeing(this,hwhm,airy,e,incl,n,psf)

    ! The seeing is simulated by summing of Airy profiles with random
    ! offsets. The distribution of offsets is Gaussian N(0,hwhm).
    ! This method is very illustrative; it simulates the real process
    ! of an image generation. On the other side, it is slow for many
    ! iterations. It gives an asymmetric and randomly noised images.

    use noise

    integer, parameter :: maxiter = 5000

    class(SprayType), intent(in) :: this
    real(REAL64), intent(in) :: hwhm, airy, e,incl
    integer, intent(in) :: n
    real(REAL64), dimension(-n:,-n:), intent(out) :: psf

    real(REAL64), dimension(:,:), allocatable :: bessel
    integer :: m,i,j,iter,l1,l2,k1,k2
    real(REAL64) :: x,y,u,v,a,b,c,s,sig
    real(REAL64), dimension(2) :: coo

    m = n / 2
    allocate(bessel(-m:m,-m:m))

    ! fill-up by Bessel
    call fill_bessel(m,this%zoom*airy,bessel)

    ! spread the profile by random shifts with N(0,sig),
    ! gnoise2(0,1) gives N(0,1), the same sigma as fill_gaussian
    sig = hwhm*this%zoom
    a = sig
    b = a*sqrt(1 - e**2)
    c = cos(incl/rad)
    s = sin(incl/rad)

    psf = 0
    do iter = 1, maxiter
       coo = gnoise2(0.0_REAL64,1.0_REAL64)
       u = coo(1)*a
       v = coo(2)*b
       x = u*c + v*s
       y =-u*s + v*c
       i = nint(x)
       j = nint(y)

       k1 = i - m
       k2 = i + m
       l1 = j - m
       l2 = j + m
       psf(k1:k2,l1:l2) = psf(k1:k2,l1:l2) + bessel
    end do

    deallocate(bessel)

  end subroutine seeing


  pure subroutine cseeing(this,hwhm,airy,e,incl,n,psf)

    ! The seeing is simulated by convolution of Airy disk function and
    ! gaussian giving spread. The convolution is computed by direct
    ! sum, so by very slow way.

    class(SprayType), intent(in) :: this
    real(REAL64), intent(in) :: hwhm, airy, e,incl
    integer, intent(in) :: n
    real(REAL64), dimension(-n:,-n:), intent(out) :: psf

    real(REAL64), dimension(:,:), allocatable :: bessel, kernel
    integer :: m,i,j,k,l
    real(REAL64) :: s

    m = n / 2
    allocate(bessel(-m:m,-m:m),kernel(-n:n,-n:n))

    ! fill-up the Bessel intensity
    call fill_bessel(m,this%zoom*airy,bessel)

    ! fill kernel
    call fill_gaussian(m,this%zoom*hwhm,e,incl,kernel)

    ! spread profile by direct convolution
    psf = 0
    do i = -m,m
       do j = -m,m
          s = 0
          do l = -m,m
             do k = -m,m
                s = s + bessel(k,l)*kernel(i-l,j-k)
             end do
          end do
          psf(i,j) = s
       end do
    end do

  end subroutine cseeing

  pure subroutine fseeing(this,hwhm,airy,e,incl,n,m,n2,psf)

    ! The seeing is simulated by convolution of Airy disk function and
    ! gaussian giving spread. The convolution is computed by using
    ! os Fourier transformation. It is most fast and precise method.

    use ftransform

    class(SprayType), intent(in) :: this
    real(REAL64), intent(in) :: hwhm, airy,e,incl
    integer, intent(in) :: n,m,n2
    real(REAL64), dimension(-n2+1:,-n2+1:), intent(out) :: psf

    complex(REAL64), dimension(:,:), allocatable :: bessel, kernel, qpsf, &
         zbessel, zkernel, zpsf
    real(REAL64), dimension(:,:), allocatable :: rkernel, rbessel
    integer, dimension(:), allocatable :: icoo
    integer :: i,j

    ! The FFT convolution is circular: both the inputs are stored
    ! with origin at the first element (the coordinates 0..n2-1
    ! followed by -n2..-1, as ifftshift does), so the origin
    ! of the result is at the first element too.
    allocate(icoo(n))
    do i = 1, n
       if( i - 1 < n2 ) then
          icoo(i) = i - 1
       else
          icoo(i) = i - 1 - n
       end if
    end do

    ! fill-up Bessel
    allocate(bessel(n,n),rbessel(-n2:n2,-n2:n2))
    call fill_bessel(n2,this%zoom*airy,rbessel)
    do j = 1, n
       do i = 1, n
          bessel(i,j) = rbessel(icoo(i),icoo(j))
       end do
    end do

    ! fill-up Gaussian spread
    allocate(kernel(n,n),rkernel(-n2:n2,-n2:n2))
    call fill_gaussian(n2,this%zoom*hwhm,e,incl,rkernel)
    ! conversion Complex to Real numbers
    do j = 1, n
       do i = 1, n
          kernel(i,j) = rkernel(icoo(i),icoo(j))
       end do
    end do

    ! forward FFT
    allocate(zbessel(n,n),zkernel(n,n))
    call ftra(bessel,zbessel,m,1)
    call ftra(kernel,zkernel,m,1)

    ! the convolution
    allocate(zpsf(n,n),qpsf(n,n))
    zpsf = zbessel * zkernel

    ! back FFT
    call ftra(zpsf,qpsf,m,-1)

    ! fill PSF by re-arranged quarters (the origin at the first element)
    psf = 0
    psf(    0:n2-1,    0:n2-1) = real(qpsf(1:n2,1:n2),REAL64)
    psf(-n2+1:-1,      0:n2-1) = real(qpsf(n2+2:n,1:n2),REAL64)
    psf(    0:n2-1,-n2+1:-1)   = real(qpsf(1:n2,n2+2:n),REAL64)
    psf(-n2+1:-1,  -n2+1:-1)   = real(qpsf(n2+2:n,n2+2:n),REAL64)

  end subroutine fseeing

  pure subroutine fill_bessel(m,sig,bessel)

    integer, intent(in) :: m
    real(REAL64), intent(in) :: sig
    real(REAL64), dimension(-m:,-m:), intent(out) :: bessel

    real(REAL64) :: x,y,r
    integer :: i,j

    ! fill-up the top-right quarter of Airy disk
    do i = 0,m
       x = i / sig
       do j = 0,m
          y = j / sig
          ! https://en.wikipedia.org/wiki/Airy_disk
          if( i == 0 .and. j == 0 ) then
             bessel(i,j) = 1
          else
             !r = sqrt(x**2 + y**2)
             r = hypot(x,y)
             bessel(i,j) = (bessel_j1(r)/r)**2
          end if
       end do
    end do

    ! fill-up rest of the kernel
    forall( i = -m:0, j = 0:m ) bessel(i,j) = bessel(-i,j)
    forall( i = -m:m, j = -m:-1 ) bessel(i,j) = bessel(i,-j)

  end subroutine fill_bessel


  pure subroutine fill_gaussian(m,sig,e,incl,kernel)

    integer, intent(in) :: m
    real(REAL64), intent(in) :: sig, e, incl
    real(REAL64), dimension(-m:,-m:), intent(out) :: kernel

    real(REAL64) :: x,y,a,b,c,s
    integer :: i,j

    a = sqrt2*sig
    b = a*sqrt(1 - e**2)
    c = cos(incl/rad)
    s = sin(incl/rad)

    do i = -m, m
       do j = -m, m
          x = ( i*c + j*s) / a
          y = (-i*s + j*c) / b
          kernel(i,j) = exp(-(x**2 + y**2))
       end do
    end do

  end subroutine fill_gaussian

  pure subroutine gauss(this,hwhm,e,incl,n,psf)

    class(SprayType), intent(in) :: this
    real(REAL64), intent(in) :: hwhm, e, incl
    integer, intent(in) :: n
    real(REAL64), dimension(-n:,-n:), intent(out) :: psf

    call fill_gaussian(n,this%zoom*hwhm,e,incl,psf)

  end subroutine gauss


  pure subroutine moffat(this, hwhm, beta, e,incl,n,psf)

    class(SprayType), intent(in) :: this
    real(REAL64), intent(in) :: hwhm, beta, e, incl
    integer, intent(in) :: n
    real(REAL64), dimension(-n:,-n:), intent(out) :: psf
    real(REAL64) :: a, b, c, s, x, y
    integer :: i,j

    a = this%zoom * hwhm
    b = a*sqrt(1 - e**2)
    c = cos(incl/rad)
    s = sin(incl/rad)

    do i = -n, n
       do j = -n, n
          x = ( i*c + j*s) / a
          y = (-i*s + j*c) / b
          psf(i,j) = (1 + (x**2 + y**2))**(-beta)
       end do
    end do

  end subroutine moffat


  subroutine PixelizePsf(this,x,y,dx,dy,xpsf,n)

    ! The PSF integrated over pixels; the pixel centres are at dx,dy
    ! (any real numbers) from the star, see psf_pixelize.

    use psfmodel

    class(SprayType) :: this
    real(REAL64), intent(in) :: x,y,dx, dy
    real(REAL64), dimension(-n:n,-n:n), intent(out) :: xpsf
    integer, intent(in) :: n

    call psf_pixelize(this%psf,ubound(this%psf,1),this%zoom,dx,dy,n,xpsf)

  end subroutine PixelizePsf

  subroutine PixelizeFile(this,x,y,dx,dy,xpsf,n)

    ! x,y are offsets from the field centre in pixels

    use psfmodel

    class(FileSprayType) :: this
    real(REAL64), intent(in) :: x,y, dx, dy
    real(REAL64), dimension(-n:n,-n:n), intent(out) :: xpsf
    integer, intent(in) :: n

    call psf_field_interpolate(this%mpsf,this%rcoo,this%acoo,x,y,this%scale, &
         this%psf)
    call PixelizePsf(this,x,y,dx,dy,xpsf,n)

  end subroutine PixelizeFile

  subroutine SprayCan(this,mpsf,rcoo,acoo,scale,zoom,nbox)

    ! mpsf(-nbox:nbox,-nbox:nbox,:) are PSFs loaded by psf_file_load

    class(FileSprayType) :: this
    real(REAL64), dimension(:,:,:), intent(in) :: mpsf
    real(REAL64), dimension(:), intent(in) :: rcoo,acoo
    real(REAL64), intent(in) :: scale
    integer, intent(in) :: zoom, nbox
    integer :: n,m

    this%zoom = zoom
    this%scale = scale
    ! the half-size of the pixelised PSF in pixels
    this%nbox = nbox / zoom + 2

    m = size(mpsf,3)
    n = nbox
    allocate(this%psf(-n:n,-n:n),this%mpsf(-n:n,-n:n,m),this%rcoo(m),this%acoo(m))
    this%psf = 0
    this%mpsf = mpsf
    this%rcoo = rcoo
    this%acoo = acoo

  end subroutine SprayCan


end module spray
