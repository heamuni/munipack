!
!  PSF models shared by the simulation (artificial) and the PSF photometry
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
!
!  A PSF is represented by a sampled (zoomed) image: ZOOM samples per
!  pixel along each axis, the sample (0,0) is the centre of the star.
!  The samples are point values of the intensity normalised so that
!  their sum is ZOOM**2 (the whole flux); the image of a star in
!  detector pixels is the integral of the (piece-wise linear) interpolant
!  of the samples over pixels, see psf_pixelize.
!
!  The same code renders stars in artificial and fits them in psfphot,
!  so the simulation and the analysis use identical PSFs.
!
!  The field-dependent PSF file (see dataform_artpsf):
!   * primary header: ZOOM (integer, samples per pixel),
!   * each image extension with keywords RADCOO [deg] (radial distance
!     from the field centre) and ANGCOO [deg] (position angle,
!     atan2(y,x) of the pixel offsets) is one PSF of the grid,
!     RADCOO = 0 is the PSF at the centre; extensions without both
!     keywords are ignored,
!   * the optional CRPIX1, CRPIX2 of an extension give the (integer,
!     1-based) pixel of the PSF centre, NAXISn/2 by default,
!   * a file without such extensions and a 2D primary image is the
!     older single-PSF format (the PSF is the same over the field).
!

module psfmodel

  use iso_fortran_env

  implicit none
  private

  real(REAL64), parameter, private :: rad = 45 / atan(1.0_REAL64)
  real(REAL64), parameter, private :: pi = 4 * atan(1.0_REAL64)

  ! the analytic PSF types
  integer, parameter, public :: PSF_GAUSS = 1, PSF_MOFFAT = 2, PSF_FILE = 3

  ! the sampling of the analytic PSFs, the same as in artificial
  integer, parameter, public :: PSF_DEFZOOM = 16

  public :: psf_file_load, psf_file_save, psf_field_interpolate, &
       psf_pixelize, psf_pixel_weights, psf_analytic_sample, &
       psf_analytic_fwhm, psf_sampled_fwhm, psf_moffat_alpha, &
       psf_moffat_usample

contains

  subroutine psf_file_load(filename,stamps,rcoo,acoo,zoom,nbox,verbose,status)

    ! Loads the field-dependent PSF. The stamps are ordered by radius
    ! and then by the angle (in 0..360), their centres are at the element
    ! (0,0); the half-size of stamps is nbox samples.

    use titsio

    character(len=*), intent(in) :: filename
    real(REAL64), dimension(:,:,:), allocatable, intent(out) :: stamps
    real(REAL64), dimension(:), allocatable, intent(out) :: rcoo, acoo
    integer, intent(out) :: zoom, nbox
    logical, intent(in) :: verbose
    integer, intent(in out) :: status

    type(fitsfiles) :: fits
    integer, dimension(:), allocatable :: ihdu, idx
    integer, dimension(:,:), allocatable :: naxes, refs
    real(REAL64), dimension(:), allocatable :: rs, as
    real(REAL64), dimension(:,:), allocatable :: tmp
    real(REAL64) :: s, cx, cy
    integer :: nhdu, hdutype, naxis, m, k, l, i, j, n1, n2
    integer, dimension(2) :: nx
    logical :: anyf, primary

    if( status /= 0 ) return

    call fits_open_file(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) "Error: failed to open the PSF file `", &
            trim(filename),"'."
       return
    end if

    call fits_read_key(fits,'ZOOM',zoom,status)
    if( status /= 0 ) then
       write(error_unit,*) "Error: the keyword ZOOM not found in `", &
            trim(filename),"'."
       goto 666
    end if
    if( zoom < 1 ) then
       write(error_unit,*) "Error: ZOOM < 1 in `",trim(filename),"'."
       status = 1
       goto 666
    end if

    call fits_get_num_hdus(fits,nhdu,status)
    if( status /= 0 ) goto 666
    allocate(ihdu(nhdu),naxes(2,nhdu),refs(2,nhdu),rs(nhdu),as(nhdu))

    ! the grid of PSFs in extensions
    m = 0
    do k = 2, nhdu
       call fits_movabs_hdu(fits,k,hdutype,status)
       if( status /= 0 ) goto 666
       if( hdutype /= FITS_IMAGE_HDU ) cycle
       call fits_get_img_dim(fits,naxis,status)
       if( naxis /= 2 ) cycle

       call fits_write_errmark
       call fits_read_key(fits,'RADCOO',rs(m+1),status)
       call fits_read_key(fits,'ANGCOO',as(m+1),status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
          if( verbose ) write(error_unit,'(a,i0,a)') &
               ' PSF file: extension ',k-1, &
               ' without RADCOO, ANGCOO keywords ignored.'
          cycle
       end if
       if( status /= 0 ) goto 666
       m = m + 1
       ihdu(m) = k
       call fits_get_image_size(fits,nx,status)
       naxes(:,m) = nx
       call psf_ref_pixel(fits,naxes(:,m),refs(:,m),status)
       if( status /= 0 ) goto 666
    end do

    ! the single PSF in the primary HDU
    primary = m == 0
    if( primary ) then
       call fits_movabs_hdu(fits,1,hdutype,status)
       call fits_get_img_dim(fits,naxis,status)
       if( status /= 0 .or. naxis /= 2 ) then
          write(error_unit,*) "Error: no PSF found in `",trim(filename),"'."
          if( status == 0 ) status = 1
          goto 666
       end if
       m = 1
       ihdu(1) = 1
       rs(1) = 0
       as(1) = 0
       call fits_get_image_size(fits,nx,status)
       naxes(:,1) = nx
       call psf_ref_pixel(fits,naxes(:,1),refs(:,1),status)
       if( status /= 0 ) goto 666
    end if

    if( any(rs(1:m) < 0) ) then
       write(error_unit,*) "Error: negative RADCOO in `",trim(filename),"'."
       status = 1
       goto 666
    end if
    as(1:m) = modulo(as(1:m),360.0_REAL64)
    where( rs(1:m) == 0 ) as(1:m) = 0

    ! order by radius, by angle
    allocate(idx(m))
    idx = [(k, k = 1, m)]
    do k = 2, m
       l = idx(k)
       j = k - 1
       do while( j >= 1 )
          if( rs(idx(j)) < rs(l) .or. &
               (rs(idx(j)) == rs(l) .and. as(idx(j)) <= as(l)) ) exit
          idx(j+1) = idx(j)
          j = j - 1
       end do
       idx(j+1) = l
    end do
    do k = 2, m
       if( rs(idx(k)) == rs(idx(k-1)) .and. as(idx(k)) == as(idx(k-1)) ) then
          write(error_unit,*) "Error: two PSFs at the same RADCOO, ANGCOO in `", &
               trim(filename),"'."
          status = 1
          goto 666
       end if
    end do

    ! the common box (centre at the element 0)
    nbox = 0
    do k = 1, m
       nbox = max(nbox,refs(1,k)-1,naxes(1,k)-refs(1,k), &
            refs(2,k)-1,naxes(2,k)-refs(2,k))
    end do

    allocate(stamps(-nbox:nbox,-nbox:nbox,m),rcoo(m),acoo(m))
    stamps = 0
    do k = 1, m
       l = idx(k)
       call fits_movabs_hdu(fits,ihdu(l),hdutype,status)
       allocate(tmp(naxes(1,l),naxes(2,l)))
       call fits_read_image(fits,0,0.0_REAL64,tmp,anyf,status)
       if( status /= 0 ) goto 666
       n1 = refs(1,l)
       n2 = refs(2,l)
       stamps(1-n1:naxes(1,l)-n1,1-n2:naxes(2,l)-n2,k) = tmp
       deallocate(tmp)
       rcoo(k) = rs(l)
       acoo(k) = as(l)
    end do

    ! checks of normalisation and centring
    do k = 1, m
       s = sum(stamps(:,:,k))
       if( abs(s / zoom**2 - 1) > 1e-3 ) &
            write(error_unit,'(a,2(f0.3,a),es10.3,a,i0,a,f0.5,a)') &
            ' Warning: PSF at RADCOO=',rcoo(k),', ANGCOO=',acoo(k), &
            ' sums to ',s,' instead of ZOOM**2=',zoom**2, &
            ': a star gets ',s/zoom**2,' of its flux.'
       if( rcoo(k) == 0 .and. s > 0 ) then
          cx = 0
          cy = 0
          do j = -nbox, nbox
             do i = -nbox, nbox
                cx = cx + i*stamps(i,j,k)
                cy = cy + j*stamps(i,j,k)
             end do
          end do
          cx = cx / s
          cy = cy / s
          if( verbose ) write(error_unit,'(a,2f8.3,a)') &
               ' PSF file: centroid of the central PSF at',cx/zoom,cy/zoom, &
               ' pix from its reference pixel.'
          if( hypot(cx,cy) > 0.25 ) write(error_unit,'(a,2f8.3,a)') &
               ' Warning: the centroid of the central PSF is offset by',cx,cy, &
               ' samples from its reference pixel (CRPIX1, CRPIX2 or NAXIS/2).'
       end if
    end do

    if( verbose ) then
       write(error_unit,'(a,i0,a,i0,a,i0,a)') ' PSF file: ',m,' PSFs, ZOOM=', &
            zoom,', box ',2*nbox+1,' samples.'
       do k = 1, m
          write(error_unit,'(a,f10.5,a,f8.2,a,f8.3,a)') '   RADCOO=',rcoo(k), &
               ' deg, ANGCOO=',acoo(k),' deg, FWHM=', &
               psf_sampled_fwhm(stamps(:,:,k),zoom),' pix'
       end do
    end if

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( status /= 0 ) then
       if( allocated(stamps) ) deallocate(stamps)
       if( allocated(rcoo) ) deallocate(rcoo,acoo)
    end if

  end subroutine psf_file_load


  subroutine psf_file_save(filename,stamps,rcoo,acoo,zoom,keys,vals,coms, &
       comments,status)

    ! Writes the PSFs stamps(-nbox:nbox,-nbox:nbox,m) (centres at the
    ! element (0,0)) at the field positions rcoo, acoo [deg] in the format
    ! read by psf_file_load: the primary header with ZOOM (and the optional
    ! real keywords keys = vals / coms and comments), one image extension
    ! with RADCOO, ANGCOO, CRPIX1, CRPIX2 for each PSF. An existing file
    ! is overwritten.

    use titsio

    character(len=*), intent(in) :: filename
    real(REAL64), dimension(:,:,:), intent(in) :: stamps
    real(REAL64), dimension(:), intent(in) :: rcoo, acoo
    integer, intent(in) :: zoom
    character(len=*), dimension(:), intent(in) :: keys, coms, comments
    real(REAL64), dimension(:), intent(in) :: vals
    integer, intent(in out) :: status

    type(fitsfiles) :: fits
    integer :: k, n

    if( status /= 0 ) return

    call fits_create_file(fits,'!'//trim(filename),status)
    if( status /= 0 ) then
       write(error_unit,*) "Error: failed to create the PSF file `", &
            trim(filename),"'."
       return
    end if

    call fits_insert_img(fits,-32,0,[0],status)
    call fits_write_key(fits,'ZOOM',zoom,'PSF samples per pixel',status)
    do k = 1, size(keys)
       call fits_write_key(fits,keys(k),vals(k),-10,coms(k),status)
    end do
    do k = 1, size(comments)
       ! (in pieces which fit to the cards)
       do n = 1, len_trim(comments(k)), 70
          call fits_write_comment(fits,comments(k)(n:min(n+69, &
               len_trim(comments(k)))),status)
       end do
    end do
    call fits_write_comment(fits,'Every extension is the PSF at the field '// &
         'position RADCOO, ANGCOO,',status)
    call fits_write_comment(fits,'the sum of samples is ZOOM**2 for the '// &
         'unit flux (see dataform_artpsf).',status)

    n = size(stamps,1)
    do k = 1, size(stamps,3)
       call fits_insert_img(fits,-32,2,[n,n],status)
       call fits_write_key(fits,'RADCOO',rcoo(k),-10, &
            '[deg] radial distance from the field centre',status)
       call fits_write_key(fits,'ANGCOO',acoo(k),-10, &
            '[deg] position angle, atan2(y,x)',status)
       call fits_write_key(fits,'CRPIX1',(n + 1)/2,'centre of the PSF',status)
       call fits_write_key(fits,'CRPIX2',(n + 1)/2,'centre of the PSF',status)
       call fits_write_image(fits,0,real(stamps(:,:,k)),status)
    end do

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine psf_file_save


  subroutine psf_ref_pixel(fits,naxes,refs,status)

    ! the reference pixel (the PSF centre) of the current HDU

    use titsio

    type(fitsfiles), intent(in) :: fits
    integer, dimension(:), intent(in) :: naxes
    integer, dimension(:), intent(out) :: refs
    integer, intent(in out) :: status

    character(len=FLEN_KEYWORD) :: key
    real(REAL64) :: crpix
    integer :: n

    ! The default is NAXISn/2 (integer division, the 1-based pixel),
    ! the convention of the PSF generators artpsf.f08 and deformpsf
    ! as well as of CRPIX of artificial frames.
    refs = naxes / 2
    do n = 1, 2
       call fits_make_keyn('CRPIX',n,key,status)
       call fits_write_errmark
       call fits_read_key(fits,key,crpix,status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
       else if( status == 0 ) then
          refs(n) = nint(crpix)
          if( abs(crpix - refs(n)) > 1e-6 ) write(error_unit,*) &
               'Warning: PSF file: non-integer ',trim(key),' rounded.'
       end if
    end do
    if( any(refs < 1 .or. refs > naxes) ) then
       write(error_unit,*) 'Error: PSF file: the reference pixel out of image.'
       status = 1
    end if

  end subroutine psf_ref_pixel


  subroutine psf_field_interpolate(stamps,rcoo,acoo,x,y,scale,obr)

    ! The PSF at the field position given by the pixel offsets x,y from
    ! the field centre; scale is in [pix/deg]. The PSFs (ordered by
    ! radius and angle) on rings of equal RADCOO are interpolated linearly
    ! in the angle (periodic), then linearly in the radius between rings;
    ! the centre PSF (RADCOO = 0) is used inside of the first ring and
    ! the outermost ring is used beyond it.

    real(REAL64), dimension(:,:,:), intent(in) :: stamps
    real(REAL64), dimension(:), intent(in) :: rcoo, acoo
    real(REAL64), intent(in) :: x,y,scale
    real(REAL64), dimension(:,:), intent(out) :: obr

    real(REAL64), dimension(:,:), allocatable :: obr2
    real(REAL64) :: r, a, dr
    integer :: m, ic, k1, k2, l1, l2

    m = size(rcoo)
    r = hypot(x,y) / scale
    a = rad*atan2(y,x)
    if( a < 0 ) a = a + 360

    ! the centre PSF
    if( rcoo(1) == 0 ) then
       ic = 1
    else
       ic = 0
    end if

    ! a PSF only at the centre
    if( ic == m ) then
       obr = stamps(:,:,1)
       return
    end if

    ! the ring k1:k2 just above the radius, or the outermost one
    k1 = ic + 1
    do
       k2 = ring_end(k1)
       if( r < rcoo(k1) .or. k2 == m ) exit
       k1 = k2 + 1
    end do

    if( r < rcoo(k1) ) then
       if( k1 == ic + 1 ) then
          ! inside of the first ring
          call ring_psf(k1,k2,obr)
          if( ic > 0 ) then
             dr = r / rcoo(k1)
             obr = (1 - dr)*stamps(:,:,ic) + dr*obr
          end if
       else
          ! between the rings
          l2 = k1 - 1
          l1 = ring_begin(l2)
          allocate(obr2(size(obr,1),size(obr,2)))
          call ring_psf(l1,l2,obr2)
          call ring_psf(k1,k2,obr)
          dr = (r - rcoo(l1)) / (rcoo(k1) - rcoo(l1))
          obr = (1 - dr)*obr2 + dr*obr
          deallocate(obr2)
       end if
    else
       ! beyond the outermost ring
       call ring_psf(k1,k2,obr)
    end if

  contains

    integer function ring_end(k)
      integer, intent(in) :: k
      ring_end = k
      do while( ring_end < m )
         if( rcoo(ring_end+1) /= rcoo(k) ) exit
         ring_end = ring_end + 1
      end do
    end function ring_end

    integer function ring_begin(k)
      integer, intent(in) :: k
      ring_begin = k
      do while( ring_begin > 1 )
         if( rcoo(ring_begin-1) /= rcoo(k) ) exit
         ring_begin = ring_begin - 1
      end do
    end function ring_begin

    subroutine ring_psf(k1,k2,obr)

      ! the angular interpolation on the ring k1:k2

      integer, intent(in) :: k1,k2
      real(REAL64), dimension(:,:), intent(out) :: obr
      real(REAL64) :: df, a1, a2
      integer :: j1, j2

      if( k1 == k2 ) then
         obr = stamps(:,:,k1)
         return
      end if

      ! the angles are ordered in 0 .. 360
      j1 = k2
      do while( j1 >= k1 )
         if( acoo(j1) <= a ) exit
         j1 = j1 - 1
      end do
      if( j1 < k1 ) then
         ! below the first angle: the sector from the last one
         j1 = k2
         j2 = k1
         a1 = acoo(j1) - 360
         a2 = acoo(j2)
      else if( j1 == k2 ) then
         j2 = k1
         a1 = acoo(j1)
         a2 = acoo(j2) + 360
      else
         j2 = j1 + 1
         a1 = acoo(j1)
         a2 = acoo(j2)
      end if
      df = (a - a1) / (a2 - a1)
      obr = (1 - df)*stamps(:,:,j1) + df*stamps(:,:,j2)

    end subroutine ring_psf

  end subroutine psf_field_interpolate


  pure subroutine psf_pixelize(psf,nmax,zoom,u,v,n,xpsf,dxpsf,dypsf)

    ! Integrates the sampled PSF over detector pixels. The pixel (i,j)
    ! of xpsf(-n:n,-n:n) has its centre at (i+u,j+v) pixels from the star
    ! centre. The samples psf(-nmax:nmax,-nmax:nmax) are interpolated
    ! linearly (even zoom) or as constants over the sample cells (odd
    ! zoom) and integrated over the pixels exactly, so the total flux
    ! sum(psf)/zoom**2 is conserved for any sub-pixel position.
    ! For positions on the sample grid (zoom*u, zoom*v integer), the
    ! integration is the trapezoid (even zoom) or the midpoint (odd zoom)
    ! rule over zoom+1 (zoom) samples.
    ! The optional dxpsf, dypsf are the derivatives by u and v
    ! (by the star position, they are with the opposite sign).

    integer, intent(in) :: nmax, zoom, n
    real(REAL64), dimension(-nmax:nmax,-nmax:nmax), intent(in) :: psf
    real(REAL64), intent(in) :: u,v
    real(REAL64), dimension(-n:n,-n:n), intent(out) :: xpsf
    real(REAL64), dimension(-n:n,-n:n), intent(out), optional :: dxpsf, dypsf

    real(REAL64), dimension(:,:), allocatable :: wx, wy, dwx, dwy, tmp, dtmp
    integer, dimension(-n:n) :: sx, lx, sy, ly
    integer :: i, j, l, s1, s2
    logical :: deriv

    deriv = present(dxpsf) .and. present(dypsf)

    allocate(wx(0:zoom+2,-n:n),wy(0:zoom+2,-n:n), &
         dwx(0:zoom+2,-n:n),dwy(0:zoom+2,-n:n))
    call psf_pixel_weights(nmax,zoom,n,u,sx,lx,wx,dwx)
    call psf_pixel_weights(nmax,zoom,n,v,sy,ly,wy,dwy)

    ! the range of samples needed along x
    s1 = nmax + 1
    s2 = -nmax - 1
    do i = -n, n
       if( lx(i) > 0 ) then
          s1 = min(s1,sx(i))
          s2 = max(s2,sx(i)+lx(i)-1)
       end if
    end do

    xpsf = 0
    if( deriv ) then
       dxpsf = 0
       dypsf = 0
    end if
    if( s1 > s2 ) return

    allocate(tmp(s1:s2,-n:n))
    if( deriv ) allocate(dtmp(s1:s2,-n:n))

    ! integration along y
    do j = -n, n
       tmp(:,j) = 0
       do l = 0, ly(j) - 1
          tmp(:,j) = tmp(:,j) + psf(s1:s2,sy(j)+l)*wy(l,j)
       end do
       if( deriv ) then
          dtmp(:,j) = 0
          do l = 0, ly(j) - 1
             dtmp(:,j) = dtmp(:,j) + psf(s1:s2,sy(j)+l)*dwy(l,j)
          end do
       end if
    end do

    ! integration along x
    do j = -n, n
       if( ly(j) == 0 ) cycle
       do i = -n, n
          if( lx(i) == 0 ) cycle
          xpsf(i,j) = dot_product(wx(0:lx(i)-1,i), &
               tmp(sx(i):sx(i)+lx(i)-1,j)) / zoom**2
          if( deriv ) then
             dxpsf(i,j) = dot_product(dwx(0:lx(i)-1,i), &
                  tmp(sx(i):sx(i)+lx(i)-1,j)) / zoom**2
             dypsf(i,j) = dot_product(wx(0:lx(i)-1,i), &
                  dtmp(sx(i):sx(i)+lx(i)-1,j)) / zoom**2
          end if
       end do
    end do

  end subroutine psf_pixelize


  pure subroutine psf_pixel_weights(nmax,zoom,n,u,sw,lw,w,dw)

    ! The weights of the samples -nmax..nmax along one axis for the pixels
    ! -n..n, whose centres are at i+u pixels from the star centre: the
    ! pixel i is the sum of samples sw(i) .. sw(i)+lw(i)-1 times w(0:,i)
    ! (divided by zoom along each axis, see psf_pixelize); dw are the
    ! derivatives of the weights by u. The interpolant of samples is
    ! linear (even zoom) or constant over the sample cells (odd zoom).
    ! The pixelisation is separable: the image of a star is
    !   p(i,j) = sum_k sum_l wx(k-sx(i),i) wy(l-sy(j),j) psf(k,l) / zoom**2.

    integer, intent(in) :: nmax, zoom, n
    real(REAL64), intent(in) :: u
    integer, dimension(-n:), intent(out) :: sw, lw
    real(REAL64), dimension(0:,-n:), intent(out) :: w, dw

    real(REAL64) :: a, b, h
    integer :: i, s, k1, k2

    ! the half-width of the interpolation kernel in samples
    if( mod(zoom,2) == 0 ) then
       h = 1
    else
       h = 0.5_REAL64
    end if

    do i = -n, n
       ! the pixel borders in samples
       a = zoom*(i + u - 0.5_REAL64)
       b = zoom*(i + u + 0.5_REAL64)
       k1 = max(ceiling(a - h),-nmax)
       k2 = min(floor(b + h),nmax)
       ! remove samples of zero weights on both ends
       do while( k1 <= k2 )
          if( b - k1 > -h .and. a - k1 < h ) exit
          k1 = k1 + 1
       end do
       do while( k2 >= k1 )
          if( b - k2 > -h .and. a - k2 < h ) exit
          k2 = k2 - 1
       end do
       sw(i) = k1
       lw(i) = max(k2 - k1 + 1,0)
       do s = k1, k2
          w(s-k1,i) = cumkernel(b - s) - cumkernel(a - s)
          dw(s-k1,i) = zoom*(kernel(b - s) - kernel(a - s))
       end do
    end do

  contains

    pure real(REAL64) function cumkernel(x)

      ! the integral of the interpolation kernel from -infinity to x

      real(REAL64), intent(in) :: x

      if( mod(zoom,2) == 0 ) then
         ! the triangle (linear interpolation)
         if( x <= -1 ) then
            cumkernel = 0
         else if( x <= 0 ) then
            cumkernel = (x + 1)**2 / 2
         else if( x < 1 ) then
            cumkernel = 1 - (1 - x)**2 / 2
         else
            cumkernel = 1
         end if
      else
         ! the box (constant over the sample cell)
         cumkernel = min(max(x + 0.5_REAL64,0.0_REAL64),1.0_REAL64)
      end if

    end function cumkernel

    pure real(REAL64) function kernel(x)

      real(REAL64), intent(in) :: x

      if( mod(zoom,2) == 0 ) then
         kernel = max(1 - abs(x),0.0_REAL64)
      else
         if( abs(x) < 0.5_REAL64 ) then
            kernel = 1
         else
            kernel = 0
         end if
      end if

    end function kernel

  end subroutine psf_pixel_weights


  pure subroutine psf_analytic_sample(ptype,width,beta,q,pa,zoom,nmax,psf)

    ! The analytic PSF sampled with zoom samples per pixel, normalised
    ! so that its integral over the plane is zoom**2 (the unit flux):
    !  GAUSS:  exp(-(x/a)**2 - (y/b)**2), a = sqrt(2)*width (width is
    !          sigma along the major axis),
    !  MOFFAT: (1 + (x/a)**2 + (y/b)**2)**(-beta), a = width (alpha),
    ! where b = q*a (q the axis ratio) and x,y are along the major
    ! and the minor axis, the major axis is inclined by pa [deg] from
    ! the x-axis (artificial: q = sqrt(1 - eccentricity**2), pa = inclination).

    integer, intent(in) :: ptype, zoom, nmax
    real(REAL64), intent(in) :: width, beta, q, pa
    real(REAL64), dimension(-nmax:nmax,-nmax:nmax), intent(out) :: psf

    real(REAL64) :: a, b, c, s, x, y, norm
    integer :: i, j

    c = cos(pa / rad)
    s = sin(pa / rad)

    if( ptype == PSF_GAUSS ) then
       a = sqrt(2.0_REAL64)*width*zoom
       b = a*q
       norm = zoom**2 / (pi*a*b)
       do j = -nmax, nmax
          do i = -nmax, nmax
             x = ( i*c + j*s) / a
             y = (-i*s + j*c) / b
             psf(i,j) = norm*exp(-(x**2 + y**2))
          end do
       end do
    else
       a = width*zoom
       b = a*q
       norm = zoom**2 * (beta - 1) / (pi*a*b)
       do j = -nmax, nmax
          do i = -nmax, nmax
             x = ( i*c + j*s) / a
             y = (-i*s + j*c) / b
             psf(i,j) = norm*(1 + (x**2 + y**2))**(-beta)
          end do
       end do
    end if

  end subroutine psf_analytic_sample


  elemental real(REAL64) function psf_analytic_fwhm(ptype,width,beta,q)

    ! FWHM [pix] of an analytic PSF, the geometric mean of both axes

    integer, intent(in) :: ptype
    real(REAL64), intent(in) :: width, beta, q

    if( ptype == PSF_GAUSS ) then
       psf_analytic_fwhm = 2*sqrt(2*log(2.0_REAL64))*width*sqrt(q)
    else
       psf_analytic_fwhm = 2*width*sqrt(2**(1/beta) - 1)*sqrt(q)
    end if

  end function psf_analytic_fwhm


  elemental real(REAL64) function psf_moffat_alpha(fwhm,beta)

    ! the Moffat alpha for the given FWHM

    real(REAL64), intent(in) :: fwhm, beta

    psf_moffat_alpha = fwhm / (2*sqrt(2**(1/beta) - 1))

  end function psf_moffat_alpha


  pure subroutine psf_moffat_usample(fwhm,u,q,pa,zoom,nmax,psf)

    ! The Moffat PSF by its FWHM along the major axis and u = 1/beta,
    ! sampled and normalised as psf_analytic_sample (MOFFAT):
    !
    !   (1 + 4 (2**u - 1) r**2 / fwhm**2)**(-1/u),  r**2 = x**2 + y**2/q**2,
    !
    ! written as exp(-s L(u s)), s = 4 E r**2 / fwhm**2, E = (2**u - 1)/u,
    ! L(t) = log(1 + t)/t. For u -> 0 (beta -> infinity), E -> ln 2 and
    ! L -> 1: the Gaussian of the same FWHM, exp(-4 ln 2 r**2 / fwhm**2),
    ! reached smoothly (u = 0 is the Gaussian itself). The shape fit of
    ! psfphot uses it: the Moffat of the parameters (alpha, beta) runs
    ! towards beta, alpha -> infinity on Gaussian stars, while u = 0 is
    ! a finite boundary of (fwhm, u) where the fit converges.

    real(REAL64), intent(in) :: fwhm, u, q, pa
    integer, intent(in) :: zoom, nmax
    real(REAL64), dimension(-nmax:nmax,-nmax:nmax), intent(out) :: psf

    real(REAL64) :: e, f2, c, sn, x, y, r2, sr, norm
    integer :: i, j

    c = cos(pa / rad)
    sn = sin(pa / rad)
    e = log(2.0_REAL64) * expm1x(u*log(2.0_REAL64))
    f2 = (fwhm*zoom)**2
    norm = zoom**2 * 4 * (1 - u) * e / (pi*f2*q)
    do j = -nmax, nmax
       do i = -nmax, nmax
          x = ( i*c + j*sn)
          y = (-i*sn + j*c) / q
          r2 = x**2 + y**2
          sr = 4*e*r2 / f2
          psf(i,j) = norm*exp(-sr*log1px(u*sr))
       end do
    end do

  contains

    pure real(REAL64) function expm1x(x)

      ! (exp(x) - 1)/x, accurate for small x (Kahan)

      real(REAL64), intent(in) :: x
      real(REAL64) :: w

      w = exp(x)
      if( w == 1 ) then
         expm1x = 1
      else
         expm1x = (w - 1) / log(w)
      end if

    end function expm1x

    pure real(REAL64) function log1px(t)

      ! log(1 + t)/t, accurate for small t (Kahan)

      real(REAL64), intent(in) :: t
      real(REAL64) :: w

      w = 1 + t
      if( w == 1 ) then
         log1px = 1
      else
         log1px = log(w) / (w - 1)
      end if

    end function log1px

  end subroutine psf_moffat_usample


  pure real(REAL64) function psf_sampled_fwhm(psf,zoom)

    ! FWHM [pix] of a sampled PSF: the diameter of the circle with the area
    ! of samples above half of the maximum

    real(REAL64), dimension(:,:), intent(in) :: psf
    integer, intent(in) :: zoom

    real(REAL64) :: h

    h = maxval(psf) / 2
    psf_sampled_fwhm = 2*sqrt(count(psf >= h) / pi) / zoom

  end function psf_sampled_fwhm

end module psfmodel
