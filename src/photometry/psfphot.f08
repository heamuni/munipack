!
!  PSF photometry
!
!  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
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
!  Stars detected by find are fitted by a PSF: an analytic one (GAUSS,
!  MOFFAT; their width can be fitted on bright isolated stars), or the
!  field-dependent PSF by a file (the same file and the same code as
!  in artificial). The stars are fitted one by one from the brightest,
!  the fitted neighbours are subtracted; the next passes refit every
!  star with all others subtracted.
!

program psfphot

  use titsio
  use psfmodel
  use psffit
  use fitspsfphot
  use iso_fortran_env

  implicit none

  real(REAL64), parameter :: rad = 45 / atan(1.0_REAL64)
  real(REAL64), parameter :: fwhm2sig = 1 / (2*sqrt(2*log(2.0_REAL64)))

  character(len=4*FLEN_FILENAME) :: record,key,val
  character(len=FLEN_FILENAME) :: filename, output, psffile = '', residual = ''
  character(len=FLEN_KEYWORD), dimension(3) :: fkeys
  character(len=FLEN_VALUE) :: psfkind = 'GAUSS', skymode = 'FIXED'
  character(len=80) :: msg
  real(REAL64) :: fwhm_par = -1, beta_par = -1, fitrad_par = -1, &
       saturate_par = -1, gain_par = -1, readns_par = -1, scale_par = -1
  real(REAL64), dimension(2) :: psfcen_par = 0
  logical :: psfcen_set = .false., elliptic = .false., verbose = .false.
  logical :: plog = .false., exitus = .true., ex
  integer :: maxiter = 50, npass = 3, eq, stat, status, nfiles

  ! the field-dependent PSF
  real(REAL64), dimension(:,:,:), allocatable :: stamps
  real(REAL64), dimension(:), allocatable :: rcoo, acoo
  integer :: zoom, nbox

  ! the current frame, PSF and stars
  type(PsfFrame) :: frame
  type(FrameType) :: fr
  type(PsfType) :: psf
  type(StarFit), dimension(:), allocatable :: stars
  real(REAL64), dimension(:), allocatable :: fwhms, snr
  integer, dimension(:), allocatable :: order, cell_head, cell_next
  real(REAL64), dimension(:), allocatable :: satrad
  real(REAL64) :: cell
  integer :: nstars, ncx, ncy, nfind

  fkeys = [ character(len=FLEN_KEYWORD) :: &
       FITS_KEY_SATURATE, FITS_KEY_READNS, FITS_KEY_GAIN ]
  nfiles = 0

  do
     read(*,'(a)',iostat=stat,iomsg=msg) record
     if( stat == IOSTAT_END ) exit
     if( stat > 0 ) then
        write(error_unit,*) trim(msg)
        error stop 'Some input error.'
     end if

     eq = index(record,'=')
     if( eq == 0 ) stop 'Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'PIPELOG' ) then

        read(val,*) plog

     else if( key == 'FITS_KEY_SATURATE' ) then

        read(val,*) fkeys(1)

     else if( key == 'FITS_KEY_READNOISE' ) then

        read(val,*) fkeys(2)

     else if( key == 'FITS_KEY_GAIN' ) then

        read(val,*) fkeys(3)

     else if( key == 'PSF' ) then

        read(val,*) psfkind
        if( .not. (psfkind == 'GAUSS' .or. psfkind == 'MOFFAT' .or. &
             psfkind == 'FILE') ) then
           write(error_unit,*) "Error: unknown PSF `",trim(psfkind), &
                "' (GAUSS, MOFFAT, FILE)."
           stop 'Error: unknown PSF.'
        end if

     else if( key == 'PSFFILE' ) then

        read(val,*) psffile

     else if( key == 'FWHM' ) then

        read(val,*) fwhm_par
        if( .not. (fwhm_par > 0) ) stop 'Error: FWHM > 0 unsatisfied.'

     else if( key == 'BETA' ) then

        read(val,*) beta_par
        if( .not. (beta_par > 1) ) stop 'Error: beta > 1 unsatisfied.'

     else if( key == 'ELLIPTIC' ) then

        read(val,*) elliptic

     else if( key == 'FITRAD' ) then

        read(val,*) fitrad_par
        if( .not. (fitrad_par > 0) ) stop 'Error: fitrad > 0 unsatisfied.'

     else if( key == 'SKY' ) then

        read(val,*) skymode
        if( .not. (skymode == 'FIXED' .or. skymode == 'LOCAL') ) &
             stop 'Error: sky is FIXED or LOCAL.'

     else if( key == 'MAXITER' ) then

        read(val,*) maxiter

     else if( key == 'NPASS' ) then

        read(val,*) npass
        if( npass < 1 ) stop 'Error: npass >= 1 unsatisfied.'

     else if( key == 'SATURATE' ) then

        read(val,*) saturate_par

     else if( key == 'GAIN' ) then

        read(val,*) gain_par
        if( .not. (gain_par > 0) ) stop 'Error: gain > 0 unsatisfied.'

     else if( key == 'READNOISE' ) then

        read(val,*) readns_par

     else if( key == 'SCALE' ) then

        read(val,*) scale_par

     else if( key == 'PSFCENTRE' ) then

        read(val,*) psfcen_par
        psfcen_set = .true.

     else if( key == 'RESIDUAL' ) then

        read(val,*) residual

     else if( key == 'FILE' ) then

        read(val,*) filename, output
        nfiles = nfiles + 1

        if( psfkind == 'FILE' .and. .not. allocated(stamps) ) then
           if( psffile == '' ) stop 'Error: PSF=FILE requires a PSF file.'
           status = 0
           call psf_file_load(psffile,stamps,rcoo,acoo,zoom,nbox,verbose,status)
           if( status /= 0 ) stop 'Error: failed to read the PSF file.'
        end if

        if( verbose ) write(error_unit,*) trim(filename)//": "
        call the_psfphot(ex)
        exitus = exitus .and. ex

     end if

  end do

  if( exitus ) then
     stop 0
  else
     stop 'Some error(s) occurred during this PSFPHOT run.'
  end if

contains

  subroutine the_psfphot(exitus)

    logical, intent(out) :: exitus

    integer, parameter :: ncols = 13
    type(PsfHeader) :: hdr
    real(REAL64), dimension(:,:), allocatable :: table
    integer, dimension(:), allocatable :: flags, npix
    character(len=FLEN_VALUE), dimension(ncols) :: ttype, tunit
    character(len=FLEN_VALUE) :: cunit
    real(REAL64) :: fwhm, t0, t1, rmod
    integer :: n, k, l, pass, status, nmod, nshape
    logical :: shape_fit, localsky
    type(FrameType) :: fr0
    type(PsfType) :: psf0

    ! reset of the frame data
    fr = fr0
    psf = psf0
    if( allocated(stars) ) deallocate(stars,order,fwhms,snr,cell_head,cell_next)
    if( allocated(satrad) ) deallocate(satrad)

    exitus = .false.
    status = 0
    call cpu_time(t0)

    call fits_psfphot_read(filename,fkeys,frame,status)
    if( status /= 0 ) return

    nstars = size(frame%x)
    localsky = skymode == 'LOCAL'

    ! the detector
    if( saturate_par > 0 ) frame%saturate = saturate_par
    if( gain_par > 0 ) then
       frame%gain = gain_par
       frame%gaininfo = 'option'
    end if
    if( readns_par >= 0 ) then
       frame%readns = readns_par
       frame%rninfo = 'option'
    end if
    fr%gain = frame%gain
    fr%rn2 = frame%readns**2
    fr%bias = frame%bias
    if( frame%bitpix > 0 ) then
       fr%qnoise = 1.0_REAL64 / 12
    else
       fr%qnoise = 0
    end if
    allocate(fr%mstar(frame%naxes(1),frame%naxes(2)), &
         fr%valid(frame%naxes(1),frame%naxes(2)), &
         fr%satur(frame%naxes(1),frame%naxes(2)))
    call move_alloc(frame%data,fr%data)
    fr%mstar = 0
    ! like find and aphot: valid pixels are lobad < data < saturation
    fr%satur = .not. (fr%data < frame%saturate)
    fr%valid = fr%data < frame%saturate .and. fr%data > frame%lobad &
         .and. fr%data == fr%data

    ! saturated stars are not detected by find: they are added, so they
    ! are modelled and subtracted from their neighbours
    nfind = nstars
    call add_saturated

    if( verbose ) then
       write(error_unit,'(a,i0,a,i0,a,i0,a,i0,a)') ' Frame ',frame%naxes(1),'x', &
            frame%naxes(2),', stars: ',nfind,' (+',nstars-nfind, &
            ' saturated not found by find)'
       write(error_unit,'(a,g0.6,3a)') ' Gain: ',frame%gain,' e-/ADU (', &
            trim(frame%gaininfo),')'
       write(error_unit,'(a,g0.6,3a)') ' Read noise: ',frame%readns,' ADU (', &
            trim(frame%rninfo),')'
       write(error_unit,'(a,g0.6,a)') ' Bias: ',frame%bias,' ADU'
       write(error_unit,'(a,g0.10,a)') ' Saturation: ',frame%saturate,' ADU'
       write(error_unit,'(2a)') ' Sky: ',trim(skymode)
    end if

    allocate(stars(nstars),order(nstars),fwhms(nstars),snr(nstars))

    ! the PSF
    if( psfkind == 'FILE' ) then
       psf%ptype = PSF_FILE
       psf%zoom = zoom
       psf%nbox = nbox
       allocate(psf%stamps,source=stamps)
       allocate(psf%rcoo,source=rcoo)
       allocate(psf%acoo,source=acoo)
       if( psfcen_set ) then
          psf%xcen = psfcen_par(1)
          psf%ycen = psfcen_par(2)
       else
          ! the same reference pixel as artificial
          psf%xcen = frame%naxes(1) / 2
          psf%ycen = frame%naxes(2) / 2
       end if
       if( scale_par > 0 ) then
          psf%scale = scale_par
       else
          psf%scale = frame%scale
       end if
       if( .not. (psf%scale > 0) .and. size(rcoo) > 1 ) then
          write(error_unit,*) 'Error: the scale [pix/deg] is undefined ', &
               '(no WCS in the header): use --scale.'
          goto 666
       end if
       fwhm = psf%fwhm(psf%xcen,psf%ycen)
       if( verbose ) write(error_unit,'(a,2f10.2,a,g0.10,a,f0.3,a)') &
            ' PSF by file, field centre:',psf%xcen,psf%ycen,' pix, scale: ', &
            psf%scale,' pix/deg, FWHM(centre): ',fwhm,' pix'
    else
       if( psfkind == 'GAUSS' ) then
          psf%ptype = PSF_GAUSS
       else
          psf%ptype = PSF_MOFFAT
       end if
       psf%zoom = PSF_DEFZOOM
       if( beta_par > 0 ) psf%beta = beta_par
       if( fwhm_par > 0 ) then
          fwhm = fwhm_par
       else if( frame%fwhm > 0 ) then
          fwhm = frame%fwhm
       else if( frame%aphwhm > 0 ) then
          fwhm = 2*frame%aphwhm
       else
          fwhm = 3
       end if
       call set_width(fwhm)
    end if

    ! fitting radii
    do k = 1, nstars
       stars(k)%x = frame%x(k)
       stars(k)%y = frame%y(k)
       stars(k)%flux = 0
       stars(k)%flag = 0
       stars(k)%satrad = satrad(k)
    end do
    call set_radii

    ! the initial sky
    do k = 1, nstars
       if( frame%have_aper .and. frame%sky(k) > 0 .and. &
            frame%skyerr(k) > 0 ) then
          stars(k)%sky = frame%sky(k)
          stars(k)%skyerr = frame%skyerr(k)
       else
          call local_sky(stars(k))
       end if
    end do

    ! the spatial grid for neighbours
    cell = max(maxval(stars%fitrad,1,nstars > 0),4.0_REAL64)
    ncx = int(frame%naxes(1) / cell) + 1
    ncy = int(frame%naxes(2) / cell) + 1
    allocate(cell_head(ncx*ncy),cell_next(nstars))
    cell_head = 0
    do k = 1, nstars
       l = cellof(stars(k)%x,stars(k)%y)
       cell_next(k) = cell_head(l)
       cell_head(l) = k
    end do

    ! the shape of the analytic PSF by bright isolated stars
    shape_fit = psf%ptype /= PSF_FILE .and. fwhm_par < 0 .and. nstars > 0
    hdr%fwhm_fitted = shape_fit
    nshape = 0
    if( shape_fit ) then
       call psf%prepare(ceiling(max(3*fwhm,6.0_REAL64)) + 4)
       call estimate_shape(nshape)
       call set_radii
    end if

    ! the final PSF
    if( psf%ptype == PSF_FILE ) then
       nmod = nbox / zoom + 1
    else
       fwhm = psf_analytic_fwhm(psf%ptype,psf%width,psf%beta,psf%q)
       rmod = max(2*maxval(stars%fitrad,1,nstars > 0),3*fwhm)
       ! Moffat wings: the radius with 0.5% of the flux outside
       if( psf%ptype == PSF_MOFFAT ) &
            rmod = max(rmod,psf%width*sqrt(0.005_REAL64**(1/(1-psf%beta)) - 1))
       nmod = min(nint(rmod) + 1,64)
       call psf%prepare(nmod)
    end if
    do k = 1, nstars
       stars(k)%modrad = nmod
    end do

    if( verbose .and. psf%ptype /= PSF_FILE ) then
       write(error_unit,'(3a,f0.3,a)',advance='no') ' PSF: ',trim(psfkind), &
            ', FWHM=',fwhm,' pix'
       if( psf%ptype == PSF_GAUSS ) then
          write(error_unit,'(a,f0.4)',advance='no') ', sigma=',psf%width
       else
          write(error_unit,'(2(a,f0.4))',advance='no') ', alpha=',psf%width, &
               ', beta=',psf%beta
       end if
       write(error_unit,'(a,f0.4,a,f0.1,a)',advance='no') ', axis ratio=',psf%q, &
            ', angle=',psf%pa,' deg'
       if( shape_fit ) then
          write(error_unit,'(a,i0,a)') ' (fitted on ',nshape,' stars)'
       else
          write(error_unit,*)
       end if
    end if

    ! the initial fluxes at detected positions, the order by brightness
    do k = 1, nstars
       call linear_flux(fr,psf,stars(k),localsky)
    end do
    call sort_order(stars%flux)

    ! passes over all stars: subtract the fitted neighbours
    do pass = 1, npass
       ! the fixed sky in the annulus of the residual frame (all fitted
       ! stars subtracted, their light does not bias the sky)
       if( .not. localsky .and. pass > 1 ) then
          do k = 1, nstars
             call annulus_sky(stars(k))
          end do
       end if
       do n = 1, nstars
          k = order(n)
          if( stars(k)%modelled ) call add_model(fr,psf,stars(k),-1)
          stars(k)%modelled = .false.
          call fit_star(fr,psf,stars(k),localsky,maxiter,3.0_REAL64,rcoreof(k))
          if( mod(stars(k)%flag/FLAG_FAILED,2) == 0 ) then
             call add_model(fr,psf,stars(k),1)
             stars(k)%modelled = .true.
          end if
       end do
       if( verbose ) then
          l = count(mod(stars%flag/FLAG_FAILED,2) == 0)
          write(error_unit,'(a,i0,a,i0,a,i0,a,f0.3)') ' Pass ',pass,': ', &
               l,' of ',nstars,' stars fitted, median reduced chi2: ', &
               median(pack(stars%chi2,mod(stars%flag/FLAG_FAILED,2) == 0))
       end if
       if( pass < npass ) call sort_order(stars%flux)
    end do

    ! a neighbour within the fitting radius
    do k = 1, nstars
       if( neighbour(k,stars(k)%fitrad) > 0 ) &
            stars(k)%flag = stars(k)%flag + FLAG_CROWDED
       if( k > nfind ) stars(k)%flag = stars(k)%flag + FLAG_ADDED
    end do

    ! output
    allocate(table(nstars,ncols-2),flags(nstars),npix(nstars))
    do k = 1, nstars
       table(k,:) = [stars(k)%x, stars(k)%y, stars(k)%xerr, stars(k)%yerr, &
            stars(k)%flux, stars(k)%fluxerr, stars(k)%sky, stars(k)%skyerr, &
            stars(k)%chi2, stars(k)%sharp, stars(k)%round]
       flags(k) = stars(k)%flag
       npix(k) = stars(k)%npix
    end do

    if( frame%bunit /= '' ) then
       cunit = frame%bunit
    else
       cunit = 'adu'
    end if
    ttype = [ character(len=FLEN_VALUE) :: FITS_COL_X, FITS_COL_Y, &
         FITS_COL_XERR, FITS_COL_YERR, FITS_COL_PSFCOUNT, FITS_COL_PSFCOUNTERR, &
         FITS_COL_SKY, FITS_COL_SKYERR, FITS_COL_CHI2, FITS_COL_SHARP, &
         FITS_COL_ROUND, FITS_COL_NPIX, FITS_COL_PSFFLAG ]
    tunit = [ character(len=FLEN_VALUE) :: 'pix', 'pix', 'pix', 'pix', &
         cunit, cunit, cunit, cunit, '', '', '', '', '' ]

    hdr%psftype = psfkind
    if( psfkind == 'FILE' ) hdr%psffile = psffile
    hdr%skymode = skymode
    hdr%fwhm = fwhm
    hdr%width = psf%width
    hdr%beta = psf%beta
    hdr%q = psf%q
    hdr%pa = psf%pa
    hdr%fitrad = fitrad_par
    hdr%gain = fr%gain
    hdr%readns = frame%readns
    hdr%bias = fr%bias
    hdr%saturate = frame%saturate
    hdr%psfcen = [psf%xcen, psf%ycen]
    hdr%scale = psf%scale
    hdr%zoom = psf%zoom
    hdr%nshape = nshape
    hdr%npass = npass

    call fits_psfphot_save(filename,output,hdr,ttype,tunit,table,npix,flags,status)

    if( residual /= '' .and. status == 0 ) then
       if( nfiles > 1 ) then
          write(error_unit,*) 'Warning: the residual image is saved ', &
               'for the first file only.'
       else
          call fits_psfphot_residual(filename,residual,fr%data - fr%mstar,status)
       end if
    end if

    call cpu_time(t1)
    if( verbose ) write(error_unit,'(a,f0.2,a)') ' Time: ',t1 - t0,' s'

666 continue

    exitus = status == 0

  end subroutine the_psfphot

    subroutine set_width(fwhm)

      real(REAL64), intent(in) :: fwhm

      if( psf%ptype == PSF_GAUSS ) then
         psf%width = fwhm * fwhm2sig / sqrt(psf%q)
      else
         psf%width = psf_moffat_alpha(fwhm,psf%beta) / sqrt(psf%q)
      end if

    end subroutine set_width

    subroutine set_radii

      integer :: k
      real(REAL64) :: f

      do k = 1, nstars
         if( psf%ptype == PSF_FILE ) then
            f = psf%fwhm(stars(k)%x,stars(k)%y)
         else
            f = psf_analytic_fwhm(psf%ptype,psf%width,psf%beta,psf%q)
         end if
         fwhms(k) = f
         if( fitrad_par > 0 ) then
            stars(k)%fitrad = fitrad_par
         else
            stars(k)%fitrad = max(1.5_REAL64*f,3.0_REAL64)
         end if
         ! a saturated star is fitted by the wings around the saturated area
         if( stars(k)%satrad > 0 ) &
              stars(k)%fitrad = max(stars(k)%fitrad,stars(k)%satrad + 1.5_REAL64*f)
      end do

    end subroutine set_radii

    subroutine add_saturated

      ! Connected areas of saturated pixels; those without a detected star
      ! are added to the list (at their centroid). The radius of the area
      ! is kept for all stars on it.

      real(REAL64), parameter :: cs = 16
      integer, dimension(:), allocatable :: stack, chead, cnext, px, py
      logical, dimension(:,:), allocatable :: seen
      real(REAL64), dimension(:), allocatable :: tx, ty, ts, te, tr
      real(REAL64) :: xc, yc, r, d, dmin
      integer :: nx, ny, i, j, k, l, m, n, ns, top, i0, j0, gx, gy, ci, cj, &
           nadd, kmin

      nx = frame%naxes(1)
      ny = frame%naxes(2)
      allocate(satrad(nstars))
      satrad = 0
      if( .not. any(fr%satur) ) return

      ! a grid of detected stars
      gx = int(nx / cs) + 1
      gy = int(ny / cs) + 1
      allocate(chead(gx*gy),cnext(nstars))
      chead = 0
      do k = 1, nstars
         l = 1 + min(max(int(frame%x(k)/cs),0),gx-1) + gx*min(max(int(frame%y(k)/cs),0),gy-1)
         cnext(k) = chead(l)
         chead(l) = k
      end do

      allocate(seen(nx,ny),stack(nx*ny/4+16),px(nx*ny/4+16),py(nx*ny/4+16))
      allocate(tx(0),ty(0),tr(0))
      seen = .false.
      nadd = 0
      do j0 = 1, ny
         do i0 = 1, nx
            if( .not. fr%satur(i0,j0) .or. seen(i0,j0) ) cycle
            ! flood fill (8-connected)
            ns = 0
            top = 1
            stack(1) = i0 + nx*(j0 - 1)
            seen(i0,j0) = .true.
            do while( top > 0 )
               m = stack(top)
               top = top - 1
               i = 1 + mod(m - 1,nx)
               j = 1 + (m - 1) / nx
               if( ns < size(px) ) then
                  ns = ns + 1
                  px(ns) = i
                  py(ns) = j
               end if
               do l = max(j-1,1), min(j+1,ny)
                  do k = max(i-1,1), min(i+1,nx)
                     if( fr%satur(k,l) .and. .not. seen(k,l) .and. &
                          top < size(stack) ) then
                        seen(k,l) = .true.
                        top = top + 1
                        stack(top) = k + nx*(l - 1)
                     end if
                  end do
               end do
            end do
            xc = sum(px(1:ns)) / real(ns,REAL64)
            yc = sum(py(1:ns)) / real(ns,REAL64)
            r = 0.5_REAL64
            do n = 1, ns
               r = max(r,hypot(px(n) - xc,py(n) - yc) + 0.5_REAL64)
            end do

            ! the detected stars on the area
            dmin = huge(dmin)
            kmin = 0
            ci = min(max(int(xc/cs),0),gx-1)
            cj = min(max(int(yc/cs),0),gy-1)
            m = ceiling((r + 2) / cs)
            do l = max(cj - m,0), min(cj + m,gy-1)
               do k = max(ci - m,0), min(ci + m,gx-1)
                  n = chead(1 + k + gx*l)
                  do while( n > 0 )
                     d = hypot(frame%x(n) - xc,frame%y(n) - yc)
                     if( d < r + 2 ) satrad(n) = max(satrad(n),r)
                     if( d < dmin ) then
                        dmin = d
                        kmin = n
                     end if
                     n = cnext(n)
                  end do
               end do
            end do
            if( dmin < r + 2 ) cycle

            nadd = nadd + 1
            tx = [tx, xc]
            ty = [ty, yc]
            tr = [tr, r]
         end do
      end do

      if( nadd == 0 ) return

      ! append the added stars
      n = nstars + nadd
      allocate(ts(n),te(n))
      ts(1:nstars) = frame%sky
      te(1:nstars) = frame%skyerr
      ts(nstars+1:) = -1
      te(nstars+1:) = -1
      call move_alloc(ts,frame%sky)
      call move_alloc(te,frame%skyerr)
      frame%x = [frame%x, tx]
      frame%y = [frame%y, ty]
      satrad = [satrad, tr]
      nstars = n

    end subroutine add_saturated

    real(REAL64) function rcoreof(k)

      ! the core of a star: half of FWHM, at least one pixel

      integer, intent(in) :: k

      rcoreof = max(fwhms(k) / 2,1.0_REAL64)

    end function rcoreof

    subroutine annulus_sky(star)

      ! The sky by the clipped mean in the annulus (of aphot) of the frame
      ! with all fitted stars subtracted.

      type(StarFit), intent(in out) :: star

      real(REAL64), dimension(:), allocatable :: s
      logical, dimension(:), allocatable :: keep
      real(REAL64) :: r1, r2, r, mean, sig
      integer :: i, j, n, m, iter

      if( frame%ring(1) > 0 .and. frame%ring(2) > frame%ring(1) ) then
         r1 = frame%ring(1)
         r2 = frame%ring(2)
      else
         r1 = star%fitrad + 2
         r2 = r1 + max(5.0_REAL64,star%fitrad)
      end if
      m = ceiling(r2)
      allocate(s((2*m+1)**2))
      n = 0
      do j = nint(star%y) - m, nint(star%y) + m
         do i = nint(star%x) - m, nint(star%x) + m
            if( i < 1 .or. j < 1 .or. i > size(fr%data,1) .or. &
                 j > size(fr%data,2) ) cycle
            if( .not. fr%valid(i,j) ) cycle
            r = hypot(i - star%x,j - star%y)
            if( r1 <= r .and. r <= r2 ) then
               n = n + 1
               s(n) = fr%data(i,j) - fr%mstar(i,j)
            end if
         end do
      end do
      if( n < 10 ) return

      ! the mean clipped at 3 sigma (the residual sky is nearly symmetric,
      ! the mean is not biased like the median of Poisson counts)
      mean = median(s(1:n))
      sig = 1.4826_REAL64*median(abs(s(1:n) - mean))
      allocate(keep(n))
      keep = .true.
      do iter = 1, 5
         if( .not. (sig > 0) ) exit
         keep = abs(s(1:n) - mean) < 3*sig
         if( count(keep) < 10 ) exit
         mean = sum(s(1:n),keep) / count(keep)
         sig = sqrt(sum((s(1:n) - mean)**2,keep) / (count(keep) - 1))
      end do
      star%sky = mean
      star%skyerr = sig / sqrt(real(count(keep),REAL64))

    end subroutine annulus_sky

    subroutine local_sky(star)

      ! a robust sky in the annulus outside the fitting radius

      type(StarFit), intent(in out) :: star

      real(REAL64), dimension(:), allocatable :: s
      real(REAL64) :: r1, r2, r, med, mad
      integer :: i, j, n, m

      r1 = star%fitrad + 2
      r2 = r1 + max(5.0_REAL64,star%fitrad)
      m = ceiling(r2)
      allocate(s((2*m+1)**2))
      n = 0
      do j = nint(star%y) - m, nint(star%y) + m
         do i = nint(star%x) - m, nint(star%x) + m
            if( i < 1 .or. j < 1 .or. i > size(fr%data,1) .or. &
                 j > size(fr%data,2) ) cycle
            if( .not. fr%valid(i,j) ) cycle
            r = hypot(i - star%x,j - star%y)
            if( r1 <= r .and. r <= r2 ) then
               n = n + 1
               s(n) = fr%data(i,j)
            end if
         end do
      end do
      if( n > 5 ) then
         med = median(s(1:n))
         mad = median(abs(s(1:n) - med))
         ! clip stars in the annulus
         med = median(pack(s(1:n),abs(s(1:n) - med) < 3*1.4826_REAL64*mad &
              .or. mad == 0))
         star%sky = med
         star%skyerr = 1.2533_REAL64*1.4826_REAL64*max(mad,1e-10_REAL64)/sqrt(real(n))
      else
         star%sky = fr%bias
         star%skyerr = -1
      end if

    end subroutine local_sky

    subroutine sort_order(key)

      ! the order by the key (flux), descending, by a merge sort

      real(REAL64), dimension(:), intent(in) :: key
      integer, dimension(:), allocatable :: tmp
      integer :: k

      order = [(k, k = 1, nstars)]
      call msort(order,key,tmp)

    end subroutine sort_order

    subroutine msort(a,key,tmp)

      integer, dimension(:), intent(in out) :: a
      real(REAL64), dimension(:), intent(in) :: key
      integer, dimension(:), allocatable, intent(in out) :: tmp

      if( .not. allocated(tmp) ) allocate(tmp(size(a)))
      call msort_rec(a,key,tmp,1,size(a))

    end subroutine msort

    recursive subroutine msort_rec(a,key,tmp,l,r)

      integer, dimension(:), intent(in out) :: a, tmp
      real(REAL64), dimension(:), intent(in) :: key
      integer, intent(in) :: l, r
      integer :: m, i, j, k

      if( r <= l ) return
      m = (l + r) / 2
      call msort_rec(a,key,tmp,l,m)
      call msort_rec(a,key,tmp,m+1,r)
      i = l
      j = m + 1
      k = l
      do while( i <= m .and. j <= r )
         if( key(a(i)) >= key(a(j)) ) then
            tmp(k) = a(i)
            i = i + 1
         else
            tmp(k) = a(j)
            j = j + 1
         end if
         k = k + 1
      end do
      do while( i <= m )
         tmp(k) = a(i)
         i = i + 1
         k = k + 1
      end do
      do while( j <= r )
         tmp(k) = a(j)
         j = j + 1
         k = k + 1
      end do
      a(l:r) = tmp(l:r)

    end subroutine msort_rec

    integer function cellof(x,y)

      real(REAL64), intent(in) :: x,y
      integer :: i, j

      i = min(max(int(x / cell),0),ncx-1)
      j = min(max(int(y / cell),0),ncy-1)
      cellof = 1 + i + ncx*j

    end function cellof

    integer function neighbour(k,r)

      ! a star within the distance r from the star k (by the initial
      ! positions), zero when none

      integer, intent(in) :: k
      real(REAL64), intent(in) :: r

      integer :: i, j, i0, j0, m, l

      neighbour = 0
      m = ceiling(r / cell)
      i0 = min(max(int(frame%x(k) / cell),0),ncx-1)
      j0 = min(max(int(frame%y(k) / cell),0),ncy-1)
      do j = max(j0 - m,0), min(j0 + m,ncy-1)
         do i = max(i0 - m,0), min(i0 + m,ncx-1)
            l = cell_head(1 + i + ncx*j)
            do while( l > 0 )
               if( l /= k ) then
                  if( hypot(frame%x(l) - frame%x(k),frame%y(l) - frame%y(k)) &
                       < r ) then
                     neighbour = l
                     return
                  end if
               end if
               l = cell_next(l)
            end do
         end do
      end do

    end function neighbour

    subroutine estimate_shape(nshape)

      ! fits of bright, isolated and unsaturated stars with free PSF shape

      integer, intent(out) :: nshape

      integer, parameter :: maxshape = 30
      real(REAL64), dimension(:,:), allocatable :: shapes
      real(REAL64), dimension(4) :: pshape
      type(StarFit) :: star
      real(REAL64) :: snrmin, c1, c2, rs
      integer :: n, k, i, j, m
      logical :: ok, sat

      do k = 1, nstars
         star = stars(k)
         call linear_flux(fr,psf,star,.true.)
         if( star%fluxerr > 0 ) then
            snr(k) = star%flux / star%fluxerr
         else
            snr(k) = 0
         end if
      end do
      call sort_order(snr)

      allocate(shapes(4,maxshape))
      nshape = 0
      snrmin = 50
      do n = 1, nstars
         k = order(n)
         if( nshape >= maxshape ) exit
         if( snr(k) < snrmin ) cycle
         ! the shape is fitted in the radius of 3 FWHM (the wings)
         rs = max(3*fwhms(k),6.0_REAL64)
         ! inside of the frame
         if( stars(k)%x - rs - 2 < 1 .or. stars(k)%y - rs - 2 < 1 .or. &
              stars(k)%x + rs + 2 > frame%naxes(1) .or. &
              stars(k)%y + rs + 2 > frame%naxes(2) ) cycle
         ! isolated
         if( neighbour(k,rs + 2*fwhms(k)) > 0 ) cycle
         ! no saturated or bad pixels
         m = ceiling(rs) + 1
         sat = .false.
         do j = nint(stars(k)%y) - m, nint(stars(k)%y) + m
            do i = nint(stars(k)%x) - m, nint(stars(k)%x) + m
               if( i < 1 .or. j < 1 .or. i > frame%naxes(1) .or. &
                    j > frame%naxes(2) ) then
                  sat = .true.
               else if( .not. fr%valid(i,j) ) then
                  sat = .true.
               end if
            end do
         end do
         if( sat ) cycle

         star = stars(k)
         star%fitrad = rs
         call linear_flux(fr,psf,star,.true.)
         pshape = [psf%width, psf%beta, 0.0_REAL64, 0.0_REAL64]
         call fit_shape(fr,psf,star,psf%ptype == PSF_MOFFAT .and. beta_par < 0, &
              elliptic,maxiter,pshape,ok)
         if( .not. ok ) cycle
         nshape = nshape + 1
         shapes(:,nshape) = pshape
         if( verbose ) write(error_unit,'(a,2f9.2,a,es10.3,a,f7.4,a,f7.3,a,2f8.4)') &
              ' Shape star at',star%x,star%y,' flux',star%flux,' width',pshape(1), &
              ' beta',pshape(2),' c1,c2',pshape(3:4)
      end do

      if( nshape == 0 ) then
         write(error_unit,*) 'Warning: no star for the PSF shape fit, ', &
              'the FWHM by find used.'
         return
      end if

      psf%width = median(shapes(1,1:nshape))
      if( psf%ptype == PSF_MOFFAT .and. beta_par < 0 ) &
           psf%beta = median(shapes(2,1:nshape))
      if( elliptic ) then
         c1 = median(shapes(3,1:nshape))
         c2 = median(shapes(4,1:nshape))
         psf%q = 1 - hypot(c1,c2)
         psf%pa = rad*atan2(c2,c1) / 2
      end if

    end subroutine estimate_shape

  function median(x)

    real(REAL64) :: median
    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), dimension(:), allocatable :: y
    real(REAL64) :: t
    integer :: n, i, j

    n = size(x)
    if( n == 0 ) then
       median = 0
       return
    end if
    allocate(y,source=x)
    ! insertion sort is enough for short arrays, a shell sort otherwise
    j = 1
    do while( j < n / 3 )
       j = 3*j + 1
    end do
    do while( j >= 1 )
       do i = j + 1, n
          t = y(i)
          block
            integer :: k
            k = i
            do while( k > j )
               if( y(k-j) <= t ) exit
               y(k) = y(k-j)
               k = k - j
            end do
            y(k) = t
          end block
       end do
       j = j / 3
    end do
    if( mod(n,2) == 1 ) then
       median = y(n/2 + 1)
    else
       median = (y(n/2) + y(n/2 + 1)) / 2
    end if

  end function median

end program psfphot
