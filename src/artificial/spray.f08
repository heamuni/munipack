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
     real(REAL64) :: rgrid = 0, agrid = 0, scale = 0
     real(REAL64), dimension(:,:,:), allocatable :: mpsf
     real(REAL64), dimension(:), allocatable :: rcoo, acoo
   contains
     procedure :: SprayCan
     procedure :: Pixelize => PixelizeFile
  end type FileSprayType

  public :: SprayMe, SprayCan, LoadSprayFile, PixelizePsf, PixelizeFile

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

    class(SprayType) :: this
    real(REAL64), intent(in) :: x,y,dx, dy
    real(REAL64), dimension(-n:n,-n:n), intent(out) :: xpsf
    integer, intent(in) :: n
    ! The default lower dimension in subroutine is 1 (allocatable array?),
    ! we pass the dimensions which we needs.

    integer :: x0,y0,k,l,i,j,k1,k2,l1,l2,m,m2,nmax
    real(REAL64), dimension(:), allocatable :: wx, wy
    real(REAL64) :: we

    nmax = ubound(this%psf,1)
    x0 = nint(this%zoom*dx)
    y0 = nint(this%zoom*dy)
    m = this%zoom
    m2 = m / 2

    ! A pixel is integrated over the window k-m2 .. k+m2 of sub-pixels.
    ! For even zoom, the window has zoom+1 elements, the edge sub-pixels
    ! are shared by neighbouring pixels, so the trapezoid rule (weights 1/2
    ! on the edges) is used to conserve the total flux.
    if( mod(m,2) == 0 ) then
       we = 0.5_REAL64
    else
       we = 1
    end if
    allocate(wx(-m2:m2),wy(-m2:m2))

    xpsf = 0
    do i = -n,n
       k = i*m + x0
       k1 = max(k - m2,-nmax)
       k2 = min(k + m2, nmax)
       if( k1 > k2 ) cycle
       wx = 1
       wx(-m2) = we
       wx(m2) = we
       do j = -n,n
          l = j*m + y0
          l1 = max(l - m2,-nmax)
          l2 = min(l + m2, nmax)
          if( l1 > l2 ) cycle
          wy = 1
          wy(-m2) = we
          wy(m2) = we
          xpsf(i,j) = dot_product(wx(k1-k:k2-k), &
               matmul(this%psf(k1:k2,l1:l2),wy(l1-l:l2-l))) / m**2
       end do
    end do

    deallocate(wx,wy)

  end subroutine PixelizePsf

  subroutine PixelizeFile(this,x,y,dx,dy,xpsf,n)

    class(FileSprayType) :: this
    real(REAL64), intent(in) :: x,y, dx, dy
    real(REAL64), dimension(-n:n,-n:n), intent(out) :: xpsf
    integer, intent(in) :: n

    call psfinterpol(this%mpsf,x,y,this%rgrid,this%agrid,this%scale,this%psf)
    call PixelizePsf(this,x,y,dx,dy,xpsf,n)

  end subroutine PixelizeFile

  pure subroutine psfinterpol(mpsf,x,y,rgrid,agrid,scale,obr)

    real(REAL64), dimension(:,:,:), intent(in) :: mpsf
    real(REAL64), intent(in) :: x,y,rgrid,agrid,scale
    real(REAL64), dimension(:,:), intent(out) :: obr
    real(REAL64) :: r,a,df,dr
    integer :: m,l,k,i0,j1,j2

    r = hypot(x,y) / scale
    m = size(mpsf,1)
    l = min(int(r / rgrid),0)
    i0 = 1 + l
    a = rad*atan2(y,x)
    if( a < 0 ) a = a + 360
    k = int(a / agrid)
!    write(*,*) a,k,k*agrid
    j1 = k+2
    if( k > 9 ) j1 = 1 + mod(k,9)
    j2 = j1 + 1
    if( j2 > 9 ) j2 = 2
!    if( idx(3) > 9 ) idx(3) = 2 + mod(k,9)
!    write(*,*) 'idx=',idx,mpsf(idx)%acoo,k,x,y

    dr = (r - l*rgrid) / rgrid
    df = (a - k*agrid) / agrid

    obr = (1 - df)*mpsf(j1,:,:) + df*mpsf(j2,:,:)
!    write(*,*) x,y,i0,j1,j2,k,dr,df,r,a

    if( dr < 1 ) then
       obr = (1 - dr)*mpsf(i0,:,:) + dr*obr
    end if

  end subroutine psfinterpol

  subroutine SprayCan(this,mpsf,rcoo,acoo,rgrid,agrid,scale,zoom,nbox)

    class(FileSprayType) :: this
    real(REAL64), dimension(:,:,:), intent(in) :: mpsf
    real(REAL64), dimension(:), intent(in) :: rcoo,acoo
    real(REAL64), intent(in) :: rgrid, agrid, scale
    integer, intent(in) :: zoom, nbox
    integer :: n,m
    integer, dimension(3) :: dims

    this%rgrid = rgrid
    this%agrid = agrid
    this%zoom = zoom
    this%nbox = nbox
    this%scale = scale

    dims = shape(mpsf)
    m = dims(1)
    n = nbox
    allocate(this%psf(-n:n,-n:n),this%mpsf(m,-n:n,-n:n),this%rcoo(m),this%acoo(m))
    this%psf = 0
    this%mpsf = mpsf
    this%rcoo = rcoo
    this%acoo = acoo

  end subroutine SprayCan

  subroutine LoadSprayFile(filename,mpsf,rcoo,acoo,rgrid,agrid,zoom,nbox)

    use titsio

    character(len=*), intent(in) :: filename
    real(REAL64), dimension(:,:,:), allocatable, intent(out) :: mpsf
    real(REAL64), dimension(:), allocatable, intent(out) :: rcoo, acoo
    real(REAL64), intent(out) :: rgrid, agrid
    integer, intent(out) :: zoom,nbox
    type(fitsfiles) :: fits
    integer :: naxis,status,k,n,n1,nhdu, hdutype
    integer, dimension(2) :: naxes
    logical :: anynull

    status = 0
    call fits_open_file(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) goto 666

    call fits_read_key(fits,'RGRID',rgrid,status)
    call fits_read_key(fits,'AGRID',agrid,status)
    call fits_read_key(fits,'ZOOM',zoom,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) goto 666

    nhdu = nint(360 / agrid) + 1

    allocate(rcoo(nhdu),acoo(nhdu))

    do k = 1, nhdu

       call fits_movrel_hdu(fits,1,hdutype,status)
       if( status == FITS_BAD_HDU_NUM ) goto 666
       call fits_read_key(fits,'RADCOO',rcoo(k),status)
       call fits_read_key(fits,'ANGCOO',acoo(k),status)
       if( status == FITS_KEYWORD_NOT_FOUND ) goto 666

       call fits_get_img_dim(fits,naxis,status)
       if( naxis /= 2 ) stop 'spray.Load(): PSF should by 2D image array.'
       call fits_get_image_size(fits,naxes,status)
       if( naxes(1) /= naxes(2) ) stop 'spray.Load(): PSF should be a square.'

       n = naxes(1) / 2
       n1 = naxes(1) - (n + 1)
       if( k == 1 ) then
          nbox = n
          allocate(mpsf(nhdu,-n:n1,-n:n1))
       else
          if( nbox /= n ) then
             write(error_unit,*) 'spray.Load() is speaking:', &
                  ' All PSF pieces should have equal dimensions.'
             goto 666
          end if
       end if

       call fits_read_image(fits,1,0.0_REAL64,mpsf(k,:,:),anynull,status)
    end do

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    return

666 continue
    call fits_report_error(error_unit,status)
    if( status /= 0 ) stop 'Error: spray.Load() failed to read a FITS file with PSF.'

  end subroutine LoadSprayFile

end module spray
