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
!  MOFFAT; their width can be fitted on bright isolated stars), the
!  field-dependent PSF by a file (the same file and the same code as
!  in artificial), or the empirical PSF built on the stars of the frame
!  (an analytic core, or the PSF of a file, plus residual tables varying
!  over the field, see psfempir). The stars are fitted from the
!  brightest, the fitted neighbours are subtracted; stars closer than
!  the link distance are fitted simultaneously in groups. The next
!  passes refit every star with all others subtracted. Stars missed by
!  find can be detected in the residuals and added.
!

program psfphot

  use titsio
  use psfmodel
  use psffit
  use psfempir
  use fitspsfphot
  use skymean, only: sky_mean
  use iso_fortran_env

  implicit none

  real(REAL64), parameter :: rad = 45 / atan(1.0_REAL64)
  real(REAL64), parameter :: fwhm2sig = 1 / (2*sqrt(2*log(2.0_REAL64)))

  ! the largest group fitted simultaneously
  integer, parameter :: maxgroup = 15

  character(len=4*FLEN_FILENAME) :: record,key,val
  character(len=FLEN_FILENAME) :: filename, output, psffile = '', residual = '', &
       psflist = '', psfsave = ''
  character(len=FLEN_KEYWORD), dimension(3) :: fkeys
  character(len=FLEN_VALUE) :: psfkind = 'GAUSS', skymode = 'FIXED', psfcore = 'MOFFAT'
  character(len=80) :: msg
  real(REAL64) :: fwhm_par = -1, beta_par = -1, fitrad_par = -1, &
       saturate_par = -1, gain_par = -1, readns_par = -1, scale_par = -1, &
       psfrad_par = -1, group_par = 1.5_REAL64, psfsmooth = 0.01_REAL64, &
       resthresh_par = -1, resthresh = -1
  real(REAL64), dimension(2) :: psfcen_par = 0
  logical :: psfcen_set = .false., elliptic = .false., verbose = .false.
  logical :: plog = .false., exitus = .true., ex
  integer :: maxiter = 50, npass = 3, eq, stat, status, nfiles, &
       psfzoom_par = -1, psfvar = 0, psfstars = 100, resdetect = 0, nbuild = 3

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
  integer :: nstars, ncx, ncy, nfind, nsat, npassdone

  ! the Moffat shape fit reached beta -> infinity: the Gaussian limit used
  logical :: gauss_limit

  ! the empirical PSF: its stars, the radius of the tables [pix]
  logical :: empirical
  integer, dimension(:), allocatable :: psfsel
  integer :: npsfsel
  real(REAL64) :: psfrad

  ! groups: the members of the group g are glist(gfirst(g):gfirst(g+1)-1)
  integer, dimension(:), allocatable :: gfirst, glist
  integer :: ngroups

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
             psfkind == 'FILE' .or. psfkind == 'EMPIRICAL') ) then
           write(error_unit,*) "Error: unknown PSF `",trim(psfkind), &
                "' (GAUSS, MOFFAT, FILE, EMPIRICAL)."
           stop 'Error: unknown PSF.'
        end if

     else if( key == 'PSFFILE' ) then

        read(val,*) psffile

     else if( key == 'PSFCORE' ) then

        read(val,*) psfcore
        if( .not. (psfcore == 'GAUSS' .or. psfcore == 'MOFFAT') ) &
             stop 'Error: the core of the empirical PSF is GAUSS or MOFFAT.'

     else if( key == 'PSFZOOM' ) then

        read(val,*) psfzoom_par
        if( psfzoom_par < 1 ) stop 'Error: PSF zoom >= 1 unsatisfied.'

     else if( key == 'PSFRADIUS' ) then

        read(val,*) psfrad_par
        if( .not. (psfrad_par > 0) ) stop 'Error: PSF radius > 0 unsatisfied.'

     else if( key == 'PSFVAR' ) then

        read(val,*) psfvar
        if( psfvar < 0 .or. psfvar > 2 ) stop 'Error: PSF variation is 0, 1 or 2.'

     else if( key == 'PSFSTARS' ) then

        read(val,*) psfstars
        if( psfstars < 1 ) stop 'Error: PSF stars >= 1 unsatisfied.'

     else if( key == 'PSFLIST' ) then

        read(val,*) psflist

     else if( key == 'PSFSAVE' ) then

        read(val,*) psfsave

     else if( key == 'PSFSMOOTH' ) then

        read(val,*) psfsmooth
        if( psfsmooth < 0 ) stop 'Error: PSF smoothing >= 0 unsatisfied.'

     else if( key == 'PSFITER' ) then

        read(val,*) nbuild
        if( nbuild < 1 ) stop 'Error: PSF iterations >= 1 unsatisfied.'

     else if( key == 'GROUP' ) then

        read(val,*) group_par
        if( group_par < 0 ) stop 'Error: group distance >= 0 unsatisfied.'

     else if( key == 'RESDETECT' ) then

        read(val,*) resdetect
        if( resdetect < 0 ) stop 'Error: residual detection >= 0 unsatisfied.'

     else if( key == 'RESTHRESH' ) then

        read(val,*) resthresh_par
        if( .not. (resthresh_par > 0) ) &
             stop 'Error: residual detection threshold > 0 unsatisfied.'

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

        if( (psfkind == 'FILE' .or. (psfkind == 'EMPIRICAL' .and. &
             psffile /= '')) .and. .not. allocated(stamps) ) then
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

    integer, parameter :: ncols = 15
    type(PsfHeader) :: hdr
    real(REAL64), dimension(:,:), allocatable :: table
    real(REAL64), dimension(:), allocatable :: chi2
    integer, dimension(:), allocatable :: flags, npix, gid, gsz
    character(len=FLEN_VALUE), dimension(ncols) :: ttype, tunit
    character(len=FLEN_VALUE) :: cunit
    real(REAL64) :: fwhm, t0, t1, rmod, resrms, psfchi2
    integer :: n, k, l, pass, status, nmod, nshape, it, niter, nvar, nres
    logical :: shape_fit, localsky
    type(FrameType) :: fr0
    type(PsfType) :: psf0

    ! reset of the frame data
    fr = fr0
    psf = psf0
    if( allocated(stars) ) deallocate(stars,order,fwhms,snr,cell_head,cell_next)
    if( allocated(satrad) ) deallocate(satrad)
    if( allocated(psfsel) ) deallocate(psfsel)
    if( allocated(gfirst) ) deallocate(gfirst,glist)
    ngroups = 0
    npsfsel = 0
    npassdone = 0
    resrms = -1
    psfchi2 = -1
    niter = 0
    nres = 0
    empirical = psfkind == 'EMPIRICAL'
    gauss_limit = .false.

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
    nsat = nstars - nfind

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
    if( psfkind == 'FILE' .or. (empirical .and. psffile /= '') ) then
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
          status = 1
          goto 666
       end if
       fwhm = psf%fwhm(psf%xcen,psf%ycen)
       if( verbose ) write(error_unit,'(a,2f10.2,a,g0.10,a,f0.3,a)') &
            ' PSF by file, field centre:',psf%xcen,psf%ycen,' pix, scale: ', &
            psf%scale,' pix/deg, FWHM(centre): ',fwhm,' pix'
    else
       if( psfkind == 'GAUSS' .or. (empirical .and. psfcore == 'GAUSS') ) then
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
       ! the field geometry, for --psf-save
       if( psfcen_set ) then
          psf%xcen = psfcen_par(1)
          psf%ycen = psfcen_par(2)
       else
          psf%xcen = frame%naxes(1) / 2
          psf%ycen = frame%naxes(2) / 2
       end if
       if( scale_par > 0 ) then
          psf%scale = scale_par
       else
          psf%scale = frame%scale
       end if
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
    call make_grid

    ! the shape of the analytic PSF by bright isolated stars
    shape_fit = psf%ptype /= PSF_FILE .and. fwhm_par < 0 .and. nstars > 0
    hdr%fwhm_fitted = shape_fit
    nshape = 0
    ! the empirical PSF: the analytic core is sampled as the tables
    ! (the zoom by the initial FWHM), the core is fitted with it
    if( empirical .and. psf%ptype /= PSF_FILE ) then
       if( psfzoom_par > 0 ) then
          psf%zoom = psfzoom_par
       else if( fwhm < 5 ) then
          psf%zoom = 4
       else
          psf%zoom = 2
       end if
    end if
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
       if( psf%ptype == PSF_MOFFAT .and. .not. empirical ) &
            rmod = max(rmod,psf%width*sqrt(0.005_REAL64**(1/(1-psf%beta)) - 1))
       nmod = min(nint(rmod) + 1,64)
       if( .not. empirical ) call psf%prepare(nmod)
    end if

    ! the empirical PSF: the base with zero residual tables
    if( empirical ) then
       if( psf%ptype == PSF_FILE .and. psfzoom_par > 0 .and. &
            psfzoom_par /= zoom ) then
          write(error_unit,'(a,i0,a)') ' Warning: --psf-zoom ignored, the zoom ', &
               zoom,' of the PSF file is used.'
       end if
       if( psfrad_par > 0 ) then
          psfrad = psfrad_par
       else
          psfrad = 4*fwhm
          ! a file: the whole stamps (the PSF of the file may have a halo
          ! far beyond 4 FWHM; the sky of PSF stars is fitted outside)
          if( psf%ptype == PSF_FILE ) psfrad = max(psfrad,real(nbox / zoom,REAL64))
       end if
       psf%ntab = psf%zoom*ceiling(psfrad)
       nvar = psfvar
       psf%npoly = npoly_of(nvar)
       psf%pxc = (frame%naxes(1) + 1) / 2.0_REAL64
       psf%pyc = (frame%naxes(2) + 1) / 2.0_REAL64
       psf%pxs = max(frame%naxes(1) / 2.0_REAL64,1.0_REAL64)
       psf%pys = max(frame%naxes(2) / 2.0_REAL64,1.0_REAL64)
       allocate(psf%res(-psf%ntab:psf%ntab,-psf%ntab:psf%ntab,psf%npoly))
       psf%res = 0
       if( psf%ptype == PSF_FILE ) then
          nmod = max(nbox,psf%ntab) / zoom + 1
       else
          ! the radius with 0.1% of the flux of the core outside
          rmod = max(rmod,psfrad + 1)
          if( psf%ptype == PSF_MOFFAT ) &
               rmod = max(rmod,psf%width*sqrt(0.001_REAL64**(1/(1-psf%beta)) - 1))
          if( psf%ptype == PSF_GAUSS ) &
               rmod = max(rmod,psf%width*sqrt(2*log(1000.0_REAL64)))
          nmod = min(nint(rmod) + 1,64)
          nmod = max(nmod,ceiling(psfrad) + 1)
          call psf%prepare(nmod)
       end if
       psf%empirical = .true.
       call set_radii
    end if
    do k = 1, nstars
       stars(k)%modrad = nmod
    end do

    if( verbose .and. psf%ptype /= PSF_FILE ) then
       if( empirical ) then
          write(error_unit,'(3a,f0.3,a)',advance='no') ' PSF core: ',trim(psfcore), &
               ', FWHM=',fwhm,' pix'
       else
          write(error_unit,'(3a,f0.3,a)',advance='no') ' PSF: ',trim(psfkind), &
               ', FWHM=',fwhm,' pix'
       end if
       if( psf%ptype == PSF_GAUSS ) then
          write(error_unit,'(a,f0.4)',advance='no') ', sigma=',psf%width
          if( gauss_limit ) write(error_unit,'(a)',advance='no') &
               ' (beta -> infinity, the Gaussian limit)'
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

    ! the empirical PSF: a pass with the base, the PSF stars, the tables
    ! by the stars fitted (and their neighbours subtracted) by the
    ! current PSF, iterated
    if( empirical ) then
       call run_pass(localsky)
       call select_psf_stars
       if( npsfsel == 0 ) then
          write(error_unit,*) 'Error: no star for the empirical PSF.'
          status = 1
          goto 666
       end if
       ! enough stars for the polynomial (three per term)
       do while( nvar > 0 .and. npsfsel < 3*npoly_of(nvar) )
          nvar = nvar - 1
       end do
       if( nvar /= psfvar ) then
          write(error_unit,'(a,i0,a,i0,a,i0,a)') ' Warning: ',npsfsel, &
               ' PSF stars are too few for the variation of the order ',psfvar, &
               ', the order ',nvar,' used.'
          psf%npoly = npoly_of(nvar)
          deallocate(psf%res)
          allocate(psf%res(-psf%ntab:psf%ntab,-psf%ntab:psf%ntab,psf%npoly))
          psf%res = 0
       end if
       if( verbose ) write(error_unit,'(a,i0,a,f0.2,a,i0,a,i0)') &
            ' Empirical PSF: ',npsfsel,' stars, table radius ',psfrad, &
            ' pix, zoom ',psf%zoom,', variation order ',nvar

       do it = 1, nbuild
          allocate(chi2(npsfsel))
          call empirical_build(fr,psf,stars(psfsel(1:npsfsel)),psfrad, &
               it > 1,psfsmooth,chi2,resrms,niter)
          psfchi2 = median(chi2)
          if( verbose ) write(error_unit,'(a,i0,a,i0,a,f0.5,a,f0.3,a,i0,a)') &
               ' PSF iteration ',it,': ',npsfsel,' stars, table rms ',resrms, &
               ', median chi2 ',psfchi2,' (',niter,' CG iterations)'
          ! (after the second build: the first one is not converged yet)
          if( it >= 2 .and. it < nbuild ) &
               call reject_psf_stars(chi2,3*npoly_of(nvar))
          deallocate(chi2)
          ! the models by the new PSF
          call set_radii
          call rebuild_models
          if( it < nbuild ) call run_pass(localsky)
       end do
       fwhm = psf%fwhm(psf%xcen,psf%ycen)
       if( verbose ) write(error_unit,'(a,f0.3,a)') &
            ' Empirical PSF: FWHM(centre)=',fwhm,' pix'
    end if

    ! passes over all stars: subtract the fitted neighbours
    do pass = 1, npass
       call run_pass(localsky)
    end do

    ! stars in the residuals
    do it = 1, resdetect
       call residual_detect(n)
       nres = nres + n
       if( verbose ) write(error_unit,'(a,i0,a,i0,a)') ' Residuals: ',n, &
            ' stars found (detection ',it,')'
       if( n == 0 ) exit
       do pass = 1, 2
          call run_pass(localsky)
       end do
    end do

    ! a neighbour within the fitting radius
    do k = 1, nstars
       if( neighbour(k,stars(k)%fitrad) > 0 ) &
            stars(k)%flag = stars(k)%flag + FLAG_CROWDED
       if( k > nfind + nsat ) then
          stars(k)%flag = stars(k)%flag + FLAG_RESID
       else if( k > nfind ) then
          stars(k)%flag = stars(k)%flag + FLAG_ADDED
       end if
    end do

    ! output
    allocate(table(nstars,ncols-4),flags(nstars),npix(nstars),gid(nstars), &
         gsz(nstars))
    do k = 1, nstars
       table(k,:) = [stars(k)%x, stars(k)%y, stars(k)%xerr, stars(k)%yerr, &
            stars(k)%flux, stars(k)%fluxerr, stars(k)%sky, stars(k)%skyerr, &
            stars(k)%chi2, stars(k)%sharp, stars(k)%round]
       flags(k) = stars(k)%flag
       npix(k) = stars(k)%npix
       gid(k) = stars(k)%group
       gsz(k) = stars(k)%gsize
    end do

    if( frame%bunit /= '' ) then
       cunit = frame%bunit
    else
       cunit = 'adu'
    end if
    ttype = [ character(len=FLEN_VALUE) :: FITS_COL_X, FITS_COL_Y, &
         FITS_COL_XERR, FITS_COL_YERR, FITS_COL_PSFCOUNT, FITS_COL_PSFCOUNTERR, &
         FITS_COL_SKY, FITS_COL_SKYERR, FITS_COL_CHI2, FITS_COL_SHARP, &
         FITS_COL_ROUND, FITS_COL_NPIX, FITS_COL_PSFFLAG, FITS_COL_GROUP, &
         FITS_COL_GSIZE ]
    tunit = [ character(len=FLEN_VALUE) :: 'pix', 'pix', 'pix', 'pix', &
         cunit, cunit, cunit, cunit, '', '', '', '', '', '', '' ]

    hdr%psftype = psfkind
    if( psffile /= '' .and. (psfkind == 'FILE' .or. empirical) ) &
         hdr%psffile = psffile
    hdr%skymode = skymode
    hdr%fwhm = fwhm
    hdr%width = psf%width
    hdr%beta = psf%beta
    hdr%gauss_limit = gauss_limit
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
    hdr%group = group_par
    hdr%ngroups = ngroups
    hdr%maxgroup = maxval(stars%gsize,1,nstars > 0)
    if( nstars == 0 ) hdr%maxgroup = 0
    hdr%resdetect = resdetect
    hdr%nresid = nres
    hdr%resthresh = resthresh
    if( empirical ) then
       if( psf%ptype == PSF_FILE ) then
          hdr%core = 'FILE'
       else if( psf%ptype == PSF_GAUSS ) then
          hdr%core = 'GAUSS'
       else
          hdr%core = 'MOFFAT'
       end if
       hdr%nvar = nvar
       hdr%npsfstar = npsfsel
       hdr%psfrad = psfrad
       hdr%resrms = resrms
       hdr%psfchi2 = psfchi2
       hdr%nbuild = nbuild
       hdr%psfsave = psfsave
    end if

    call fits_psfphot_save(filename,output,hdr,ttype,tunit,table,npix,flags, &
         gid,gsz,status)

    if( residual /= '' .and. status == 0 ) then
       if( nfiles > 1 ) then
          write(error_unit,*) 'Warning: the residual image is saved ', &
               'for the first file only.'
       else
          call fits_psfphot_residual(filename,residual,fr%data - fr%mstar,status)
       end if
    end if

    if( empirical .and. psfsave /= '' .and. status == 0 ) then
       if( nfiles > 1 ) then
          write(error_unit,*) 'Warning: the empirical PSF is saved ', &
               'for the first file only.'
       else
          call save_psf(nvar,nshape,status)
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
         if( psf%ptype == PSF_FILE .or. psf%empirical ) then
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
      ! with all fitted stars subtracted. The pixels with the light of the
      ! fitted stars are left out: the photon noise of the subtracted
      ! light is symmetric, while the clip of the noise model of the sky
      ! is not (a Poisson tail above, 5 sigma below), so it raised the sky
      ! (+4.6% of a sky of 0.42 e-/pixel among the stars of a 32 pixel
      ! grid, faint stars -0.2 sigma low). The light is significant when
      ! its photon variance exceeds a quarter of the variance of the sky
      ! (all pixels are used when fewer than 10 are left).

      type(StarFit), intent(in out) :: star

      real(REAL64), dimension(:), allocatable :: s
      logical, dimension(:), allocatable :: lit
      real(REAL64) :: r1, r2, r, mean, sig, sdev, v4, g
      integer :: i, j, n, m

      if( frame%ring(1) > 0 .and. frame%ring(2) > frame%ring(1) ) then
         r1 = frame%ring(1)
         r2 = frame%ring(2)
      else
         r1 = star%fitrad + 2
         r2 = r1 + max(5.0_REAL64,star%fitrad)
      end if
      m = ceiling(r2)
      allocate(s((2*m+1)**2),lit((2*m+1)**2))
      v4 = fr%variance(star%sky) / 4
      g = max(fr%gain,epsilon(g))
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
               lit(n) = fr%mstar(i,j) / g > v4
            end if
         end do
      end do
      if( n < 10 ) return
      if( count(.not. lit(1:n)) >= 10 ) then
         s(1:count(.not. lit(1:n))) = pack(s(1:n),.not. lit(1:n))
         n = count(.not. lit(1:n))
      end if

      ! the mean clipped by the noise model (skymean): a symmetric clip
      ! by the scatter (at 3 sigma by the median absolute deviation, as
      ! before) rejects the photons of a sky of a few photons per pixel
      call sky_mean(s(1:n),fr%gain,fr%rn2,fr%bias,mean,sig,sdev,m)
      star%sky = mean
      star%skyerr = sig

    end subroutine annulus_sky

    subroutine local_sky(star)

      ! a robust sky in the annulus outside the fitting radius

      type(StarFit), intent(in out) :: star

      real(REAL64), dimension(:), allocatable :: s
      real(REAL64) :: r1, r2, r, med, mad, sdev
      integer :: i, j, n, m, nu

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
         ! the mean clipped by the noise model (the median of a sky of
         ! a few photons per pixel is biased, BUGS P23)
         call sky_mean(s(1:n),fr%gain,fr%rn2,fr%bias,med,mad,sdev,nu)
         star%sky = med
         star%skyerr = max(mad,1e-10_REAL64)
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
      real(REAL64) :: snrmin, c1, c2, rs, fw, u
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
         ! Moffat: the FWHM along the major axis and u = 1/beta
         if( psf%ptype == PSF_MOFFAT ) then
            pshape = [psf_analytic_fwhm(PSF_MOFFAT,psf%width,psf%beta,1.0_REAL64), &
                 1 / psf%beta, 0.0_REAL64, 0.0_REAL64]
         else
            pshape = [psf%width, psf%beta, 0.0_REAL64, 0.0_REAL64]
         end if
         if( empirical ) then
            ! the core: beta up to 10 (the tables correct the shape)
            call fit_shape(fr,psf,star,psf%ptype == PSF_MOFFAT .and. beta_par < 0, &
                 .true.,maxiter,pshape,ok,10.0_REAL64)
         else
            call fit_shape(fr,psf,star,psf%ptype == PSF_MOFFAT .and. beta_par < 0, &
                 elliptic,maxiter,pshape,ok)
         end if
         if( .not. ok ) cycle
         nshape = nshape + 1
         shapes(:,nshape) = pshape
         if( verbose ) then
            if( psf%ptype == PSF_MOFFAT ) then
               write(error_unit,'(a,2f9.2,a,es10.3,a,f7.4,a,f7.4,a,2f8.4)') &
                    ' Shape star at',star%x,star%y,' flux',star%flux,' fwhm',pshape(1), &
                    ' 1/beta',pshape(2),' c1,c2',pshape(3:4)
            else
               write(error_unit,'(a,2f9.2,a,es10.3,a,f7.4,a,f7.3,a,2f8.4)') &
                    ' Shape star at',star%x,star%y,' flux',star%flux,' width',pshape(1), &
                    ' beta',pshape(2),' c1,c2',pshape(3:4)
            end if
         end if
      end do

      if( nshape == 0 ) then
         write(error_unit,*) 'Warning: no star for the PSF shape fit, ', &
              'the FWHM by find used.'
         return
      end if

      if( psf%ptype == PSF_MOFFAT ) then
         ! the median FWHM and u = 1/beta; u < 0.001 (beta > 1000, Gaussian
         ! stars give 0 to 1e-3) is the Gaussian of the FWHM: the Moffat
         ! profile differs from it by u/10 of the flux (fitted on Gaussian
         ! stars: +0.09% for beta = 100, +0.47% for beta = 20)
         fw = median(shapes(1,1:nshape))
         if( beta_par < 0 ) then
            u = median(shapes(2,1:nshape))
         else
            u = 1 / psf%beta
         end if
         if( u >= 1e-3_REAL64 .or. beta_par > 0 ) then
            psf%beta = 1 / u
            psf%width = psf_moffat_alpha(fw,psf%beta)
         else
            psf%ptype = PSF_GAUSS
            psf%width = fw * fwhm2sig
            gauss_limit = .true.
         end if
      else
         psf%width = median(shapes(1,1:nshape))
      end if
      if( elliptic .or. empirical ) then
         c1 = median(shapes(3,1:nshape))
         c2 = median(shapes(4,1:nshape))
         psf%q = 1 - hypot(c1,c2)
         psf%pa = rad*atan2(c2,c1) / 2
      end if

    end subroutine estimate_shape

    pure integer function npoly_of(nvar)

      ! the terms of the polynomial of the order nvar

      integer, intent(in) :: nvar

      npoly_of = (nvar + 1)*(nvar + 2) / 2

    end function npoly_of

    subroutine make_grid

      ! the spatial grid of the detected positions

      integer :: k, l

      if( allocated(cell_head) ) deallocate(cell_head,cell_next)
      cell = max(maxval(stars%fitrad,1,nstars > 0),4.0_REAL64)
      ncx = int(frame%naxes(1) / cell) + 1
      ncy = int(frame%naxes(2) / cell) + 1
      allocate(cell_head(ncx*ncy),cell_next(nstars))
      cell_head = 0
      do k = 1, nstars
         l = cellof(frame%x(k),frame%y(k))
         cell_next(k) = cell_head(l)
         cell_head(l) = k
      end do

    end subroutine make_grid

    subroutine fit_single(k)

      ! the fit of the star k alone, the others subtracted

      integer, intent(in) :: k

      if( stars(k)%modelled ) call add_model(fr,psf,stars(k),-1)
      stars(k)%modelled = .false.
      call fit_star(fr,psf,stars(k),skymode == 'LOCAL',maxiter,3.0_REAL64, &
           rcoreof(k))
      if( mod(stars(k)%flag/FLAG_FAILED,2) == 0 ) then
         call add_model(fr,psf,stars(k),1)
         stars(k)%modelled = .true.
      end if

    end subroutine fit_single

    subroutine run_pass(localsky)

      ! A pass over all stars (groups) from the brightest, the models of
      ! the fitted stars are subtracted from the others.

      logical, intent(in) :: localsky

      type(StarFit), dimension(:), allocatable :: gs
      real(REAL64), dimension(:), allocatable :: rc
      logical, dimension(:), allocatable :: done
      integer :: n, k, l, g, m, i
      logical :: ok

      ! the fixed sky in the annulus of the residual frame (all fitted
      ! stars subtracted, their light does not bias the sky)
      if( .not. localsky .and. npassdone > 0 ) then
         do k = 1, nstars
            call annulus_sky(stars(k))
         end do
      end if
      if( npassdone > 0 ) call sort_order(stars%flux)
      call make_groups

      allocate(done(nstars))
      done = .false.
      do n = 1, nstars
         k = order(n)
         if( done(k) ) cycle
         g = stars(k)%group
         if( g == 0 ) then
            call fit_single(k)
            done(k) = .true.
            cycle
         end if

         ! the group, the members by brightness
         m = gfirst(g+1) - gfirst(g)
         allocate(gs(m),rc(m))
         do i = 1, m
            l = glist(gfirst(g)+i-1)
            if( stars(l)%modelled ) call add_model(fr,psf,stars(l),-1)
            stars(l)%modelled = .false.
            gs(i) = stars(l)
            rc(i) = rcoreof(l)
         end do
         call fit_group(fr,psf,gs,localsky,maxiter,3.0_REAL64,rc,ok)
         do i = 1, m
            l = glist(gfirst(g)+i-1)
            if( ok ) then
               stars(l) = gs(i)
               if( mod(stars(l)%flag/FLAG_FAILED,2) == 0 ) then
                  call add_model(fr,psf,stars(l),1)
                  stars(l)%modelled = .true.
               end if
            else
               ! singular: the members alone
               stars(l)%group = 0
               stars(l)%gsize = 1
               call fit_single(l)
            end if
            done(l) = .true.
         end do
         deallocate(gs,rc)
      end do
      npassdone = npassdone + 1

      if( verbose ) then
         l = count(mod(stars%flag/FLAG_FAILED,2) == 0)
         write(error_unit,'(a,i0,a,i0,a,i0,a,f0.3)',advance='no') ' Pass ', &
              npassdone,': ',l,' of ',nstars,' stars fitted, median reduced chi2: ', &
              median(pack(stars%chi2,mod(stars%flag/FLAG_FAILED,2) == 0))
         if( ngroups > 0 ) then
            write(error_unit,'(a,i0,a,i0,a)') ', ',ngroups,' groups of ', &
                 count(stars%group > 0),' stars'
         else
            write(error_unit,*)
         end if
      end if

    end subroutine run_pass

    subroutine make_groups

      ! Groups of stars closer than the link distance (group_par times
      ! the mean FWHM of the pair) by the current positions; the groups
      ! larger than maxgroup are split by shorter links (by 0.8 times)
      ! until they fit, the stars left alone are fitted alone.

      integer, dimension(:), allocatable :: parent, fin, head, nxt, cnt, pos
      real(REAL64) :: cs
      integer :: gx, gy, k, l, i, j, i0, j0, n, g

      stars%group = 0
      stars%gsize = 1
      ngroups = 0
      if( allocated(gfirst) ) deallocate(gfirst,glist)
      allocate(gfirst(1),glist(0))
      gfirst = 1
      if( .not. (group_par > 0) .or. nstars < 2 ) return

      ! the grid of the current positions, the cell of the longest link
      cs = max(group_par*maxval(fwhms(1:nstars)),1.0_REAL64)
      gx = int(frame%naxes(1) / cs) + 1
      gy = int(frame%naxes(2) / cs) + 1
      allocate(head(gx*gy),nxt(nstars),parent(nstars),fin(nstars))
      head = 0
      do k = 1, nstars
         l = 1 + min(max(int(stars(k)%x / cs),0),gx-1) + &
              gx*min(max(int(stars(k)%y / cs),0),gy-1)
         nxt(k) = head(l)
         head(l) = k
         parent(k) = k
      end do
      do k = 1, nstars
         i0 = min(max(int(stars(k)%x / cs),0),gx-1)
         j0 = min(max(int(stars(k)%y / cs),0),gy-1)
         do j = max(j0 - 1,0), min(j0 + 1,gy-1)
            do i = max(i0 - 1,0), min(i0 + 1,gx-1)
               l = head(1 + i + gx*j)
               do while( l > 0 )
                  if( l > k ) then
                     if( linked(k,l,group_par) ) call unite(parent,k,l)
                  end if
                  l = nxt(l)
               end do
            end do
         end do
      end do

      ! the components, split when too large
      fin = 0
      allocate(cnt(nstars),pos(nstars))
      cnt = 0
      do k = 1, nstars
         l = root(parent,k)
         cnt(l) = cnt(l) + 1
      end do
      ! the members of all components, by offsets
      pos = 0
      n = 0
      do k = 1, nstars
         l = root(parent,k)
         if( cnt(l) < 2 ) cycle
         if( pos(l) == 0 ) then
            n = n + 1
            pos(l) = n
         end if
      end do
      block
        integer, dimension(:), allocatable :: first, list, fill
        allocate(first(n+1),list(nstars),fill(n))
        first = 0
        do k = 1, nstars
           l = root(parent,k)
           if( cnt(l) < 2 ) cycle
           first(pos(l)) = first(pos(l)) + 1
        end do
        ! counts to offsets
        j = 1
        do g = 1, n
           i = first(g)
           first(g) = j
           j = j + i
        end do
        first(n+1) = j
        fill = first(1:n)
        do k = 1, nstars
           l = root(parent,k)
           if( cnt(l) < 2 ) cycle
           list(fill(pos(l))) = k
           fill(pos(l)) = fill(pos(l)) + 1
        end do
        do g = 1, n
           call assign_group(list(first(g):first(g+1)-1),group_par,fin)
        end do
      end block

      ! the lists of members ordered by brightness
      deallocate(cnt)
      allocate(cnt(ngroups+1))
      cnt = 0
      do k = 1, nstars
         if( fin(k) > 0 ) cnt(fin(k)) = cnt(fin(k)) + 1
      end do
      deallocate(gfirst,glist)
      allocate(gfirst(ngroups+1),glist(sum(cnt(1:ngroups))))
      gfirst(1) = 1
      do g = 1, ngroups
         gfirst(g+1) = gfirst(g) + cnt(g)
      end do
      cnt = 0
      do n = 1, nstars
         k = order(n)
         g = fin(k)
         if( g == 0 ) cycle
         glist(gfirst(g) + cnt(g)) = k
         cnt(g) = cnt(g) + 1
      end do
      do k = 1, nstars
         stars(k)%group = fin(k)
         if( fin(k) > 0 ) stars(k)%gsize = cnt(fin(k))
      end do

    end subroutine make_groups

    logical function linked(k,l,f)

      ! the stars k, l closer than f times their mean FWHM

      integer, intent(in) :: k,l
      real(REAL64), intent(in) :: f

      linked = hypot(stars(k)%x - stars(l)%x,stars(k)%y - stars(l)%y) < &
           f*(fwhms(k) + fwhms(l)) / 2

    end function linked

    recursive subroutine assign_group(list,f,fin)

      ! the group of the stars of the list (fin), or its parts by shorter
      ! links when it is larger than maxgroup

      integer, dimension(:), intent(in) :: list
      real(REAL64), intent(in) :: f
      integer, dimension(:), intent(in out) :: fin

      integer, dimension(:), allocatable :: par, sub
      integer :: m, a, b, r

      m = size(list)
      if( m < 2 ) return
      if( m <= maxgroup ) then
         ngroups = ngroups + 1
         fin(list) = ngroups
         return
      end if
      allocate(par(m))
      par = [(a, a = 1, m)]
      do a = 1, m
         do b = a + 1, m
            if( linked(list(a),list(b),0.8_REAL64*f) ) call unite(par,a,b)
         end do
      end do
      do a = 1, m
         r = root(par,a)
         if( r /= a ) cycle
         sub = pack(list,[(root(par,b) == r, b = 1, m)])
         if( size(sub) == m .and. f < 1e-3_REAL64 ) then
            ! no split possible: by chunks
            ngroups = ngroups + 1
            fin(list(1:maxgroup)) = ngroups
            call assign_group(list(maxgroup+1:),f,fin)
            return
         end if
         call assign_group(sub,0.8_REAL64*f,fin)
      end do

    end subroutine assign_group

    integer function root(parent,k)

      integer, dimension(:), intent(in out) :: parent
      integer, intent(in) :: k
      integer :: l, m

      root = k
      do while( parent(root) /= root )
         root = parent(root)
      end do
      ! path compression
      l = k
      do while( parent(l) /= root )
         m = parent(l)
         parent(l) = root
         l = m
      end do

    end function root

    subroutine unite(parent,k,l)

      integer, dimension(:), intent(in out) :: parent
      integer, intent(in) :: k,l
      integer :: a, b

      a = root(parent,k)
      b = root(parent,l)
      if( a /= b ) parent(max(a,b)) = min(a,b)

    end subroutine unite

    subroutine rebuild_models

      ! the sum of models of all stars by the current PSF

      integer :: k

      fr%mstar = 0
      do k = 1, nstars
         if( stars(k)%modelled ) call add_model(fr,psf,stars(k),1)
      end do

    end subroutine rebuild_models

    subroutine select_psf_stars

      ! The stars of the empirical PSF: by the list, or the brightest
      ! (by S/N, at least 30) unsaturated stars with no bad pixels in
      ! the radius of the table, inside of the frame, without a brighter
      ! neighbour in the table and a neighbour brighter than 10% within
      ! 3 FWHM, of normal sharpness and roundness.

      real(REAL64), dimension(:), allocatable :: key, v
      integer, dimension(:), allocatable :: idx, tmp
      logical, dimension(:), allocatable :: cand
      real(REAL64) :: med, mad
      integer :: k, n

      allocate(cand(nstars),key(nstars))
      cand = .false.
      do k = 1, nstars
         if( stars(k)%fluxerr > 0 ) then
            snr(k) = stars(k)%flux / stars(k)%fluxerr
         else
            snr(k) = 0
         end if
      end do

      if( psflist /= '' ) then
         call read_psf_list(cand)
      else
         do k = 1, nfind
            if( stars(k)%satrad > 0 ) cycle
            if( iand(stars(k)%flag,FLAG_SATURATED + FLAG_EXCLUDED + &
                 FLAG_POSFIXED + FLAG_FAILED) /= 0 ) cycle
            if( snr(k) < 30 ) cycle
            if( .not. clean_table(k) ) cycle
            if( .not. isolated(k) ) cycle
            cand(k) = .true.
         end do
         ! normal sharpness and roundness (4 sigma)
         if( count(cand) >= 5 ) then
            v = pack(stars%sharp,cand)
            med = median(v)
            mad = 1.4826_REAL64*median(abs(v - med))
            cand = cand .and. abs(stars%sharp - med) <= max(4*mad,0.02_REAL64)
            v = pack(stars%round,cand)
            med = median(v)
            mad = 1.4826_REAL64*median(abs(v - med))
            cand = cand .and. abs(stars%round - med) <= max(4*mad,0.02_REAL64)
         end if
      end if

      ! the brightest
      where( cand )
         key = snr
      elsewhere
         key = -1
      end where
      allocate(idx(nstars))
      idx = [(k, k = 1, nstars)]
      call msort(idx,key,tmp)
      n = min(count(cand),psfstars)
      npsfsel = n
      psfsel = idx(1:n)

      if( verbose ) then
         write(error_unit,'(a,i0,a,i0,a)') ' PSF stars: ',n,' (of ',count(cand), &
              ' candidates)'
         do k = 1, n
            associate(s => stars(psfsel(k)))
              write(error_unit,'(a,2f9.2,a,es10.3,a,f8.1,a,f7.3,a,f7.3)') &
                   '   PSF star at',s%x,s%y,' flux',s%flux,' S/N',snr(psfsel(k)), &
                   ' sharp',s%sharp,' round',s%round
            end associate
         end do
      end if

    end subroutine select_psf_stars

    logical function clean_table(k)

      ! the table of the star in the frame, no bad (saturated) pixels

      integer, intent(in) :: k
      integer :: i, j, m

      clean_table = .false.
      if( stars(k)%x - psfrad - 1 < 1 .or. stars(k)%y - psfrad - 1 < 1 .or. &
           stars(k)%x + psfrad + 1 > frame%naxes(1) .or. &
           stars(k)%y + psfrad + 1 > frame%naxes(2) ) return
      m = ceiling(psfrad)
      do j = nint(stars(k)%y) - m, nint(stars(k)%y) + m
         do i = nint(stars(k)%x) - m, nint(stars(k)%x) + m
            if( hypot(i - stars(k)%x,j - stars(k)%y) > psfrad ) cycle
            if( .not. fr%valid(i,j) ) return
         end do
      end do
      clean_table = .true.

    end function clean_table

    logical function isolated(k)

      ! no brighter neighbour in the table, no neighbour brighter than
      ! 10% within 3 FWHM, no saturated star reaching the table

      integer, intent(in) :: k
      real(REAL64) :: r, d
      integer :: i, j, i0, j0, m, l

      isolated = .false.
      r = max(psfrad,3*fwhms(k)) + maxval(satrad,1,nstars > 0) + 2*fwhms(k)
      m = ceiling(r / cell)
      i0 = min(max(int(frame%x(k) / cell),0),ncx-1)
      j0 = min(max(int(frame%y(k) / cell),0),ncy-1)
      do j = max(j0 - m,0), min(j0 + m,ncy-1)
         do i = max(i0 - m,0), min(i0 + m,ncx-1)
            l = cell_head(1 + i + ncx*j)
            do while( l > 0 )
               if( l /= k ) then
                  d = hypot(stars(l)%x - stars(k)%x,stars(l)%y - stars(k)%y)
                  if( d < 3*fwhms(k) .and. stars(l)%flux > 0.1_REAL64*stars(k)%flux ) &
                       return
                  if( d < psfrad .and. stars(l)%flux > stars(k)%flux ) return
                  if( satrad(l) > 0 .and. d < psfrad + satrad(l) + fwhms(k) ) return
               end if
               l = cell_next(l)
            end do
         end do
      end do
      isolated = .true.

    end function isolated

    subroutine read_psf_list(cand)

      ! The PSF stars by the list: x y [pix] per line (# comments), or
      ! the columns X, Y of the first table of a FITS file. The nearest
      ! fitted star within 2 pixels is used.

      logical, dimension(:), intent(in out) :: cand

      type(fitsfiles) :: fits
      real(REAL64), dimension(:), allocatable :: x, y
      real(REAL64) :: u, v, d, dmin
      character(len=FLEN_CARD) :: line
      integer :: iu, ios, n, k, l, kmin, nrows, xcol, ycol, st
      logical :: anyf, isfits

      n = len_trim(psflist)
      isfits = .false.
      do k = max(n - 6,1), n
         if( psflist(k:k) == '.' ) then
            line = psflist(k+1:n)
            isfits = line == 'fits' .or. line == 'fit' .or. line == 'fts' .or. &
                 line == 'FITS' .or. line == 'FIT' .or. line == 'FTS'
         end if
      end do

      allocate(x(0),y(0))
      if( isfits ) then
         st = 0
         call fits_open_table(fits,psflist,FITS_READONLY,st)
         call fits_get_num_rows(fits,nrows,st)
         call fits_get_colnum(fits,.false.,'X',xcol,st)
         call fits_get_colnum(fits,.false.,'Y',ycol,st)
         if( st == 0 ) then
            deallocate(x,y)
            allocate(x(nrows),y(nrows))
            call fits_read_col(fits,xcol,1,0.0_REAL64,x,anyf,st)
            call fits_read_col(fits,ycol,1,0.0_REAL64,y,anyf,st)
         end if
         call fits_close_file(fits,st)
         call fits_report_error(error_unit,st)
         if( st /= 0 ) then
            write(error_unit,*) "Warning: failed to read X, Y of `", &
                 trim(psflist),"'."
            return
         end if
      else
         open(newunit=iu,file=psflist,status='old',action='read',iostat=ios)
         if( ios /= 0 ) then
            write(error_unit,*) "Warning: failed to open `",trim(psflist),"'."
            return
         end if
         do
            read(iu,'(a)',iostat=ios) line
            if( ios /= 0 ) exit
            line = adjustl(line)
            if( line == '' .or. line(1:1) == '#' ) cycle
            read(line,*,iostat=ios) u, v
            if( ios /= 0 ) cycle
            x = [x, u]
            y = [y, v]
         end do
         close(iu)
      end if

      do l = 1, size(x)
         dmin = huge(dmin)
         kmin = 0
         do k = 1, nstars
            d = hypot(stars(k)%x - x(l),stars(k)%y - y(l))
            if( d < dmin ) then
               dmin = d
               kmin = k
            end if
         end do
         if( dmin > 2 .or. kmin == 0 ) then
            write(error_unit,'(a,2f9.2,a)') ' Warning: no star at',x(l),y(l), &
                 ' of the PSF list.'
            cycle
         end if
         if( stars(kmin)%satrad > 0 .or. iand(stars(kmin)%flag,FLAG_SATURATED + &
              FLAG_FAILED) /= 0 .or. .not. (snr(kmin) > 0) ) then
            write(error_unit,'(a,2f9.2,a)') ' Warning: the star at',x(l),y(l), &
                 ' of the PSF list is saturated or failed, skipped.'
            cycle
         end if
         cand(kmin) = .true.
      end do

    end subroutine read_psf_list

    subroutine reject_psf_stars(chi2,nmin)

      ! PSF stars with chi2 above the median + 5 sigma (and 1.5 times
      ! the median) are not used, when at least nmin stars remain

      real(REAL64), dimension(:), intent(in) :: chi2
      integer, intent(in) :: nmin

      logical, dimension(:), allocatable :: keep
      real(REAL64) :: med, mad
      integer :: k

      med = median(chi2)
      mad = 1.4826_REAL64*median(abs(chi2 - med))
      allocate(keep(npsfsel))
      keep = chi2 <= max(med + 5*mad,1.5_REAL64*med)
      if( count(keep) < max(nmin,1) .or. all(keep) ) return
      if( verbose ) then
         do k = 1, npsfsel
            if( .not. keep(k) ) write(error_unit,'(a,2f9.2,a,f0.2)') &
                 '   PSF star at',stars(psfsel(k))%x,stars(psfsel(k))%y, &
                 ' rejected, chi2=',chi2(k)
         end do
      end if
      psfsel = pack(psfsel(1:npsfsel),keep)
      npsfsel = size(psfsel)

    end subroutine reject_psf_stars

    subroutine residual_detect(nnew)

      ! Stars in the residuals (the frame with all fitted stars
      ! subtracted): the amplitude h of the lowered Gaussian (as find)
      ! of the median FWHM, significant above the threshold
      ! against the noise of the model, which includes the PSF mismatch
      ! perr*model (perr measured on the bright stars), local maxima
      ! of normal sharpness (the limits of find), farther than 0.75 FWHM
      ! from all stars and not on saturated stars.
      !
      ! The threshold (resthresh_par, by default) is the one of find, at
      ! least the Normal deviate exceeded by one pixel of the frame by
      ! chance (P = 1/npix: 5.0 for 2048**2 pixels, 4.75 for 1024**2,
      ! 5.3 for 4096**2): the residuals are searched in every pixel,
      ! and the skewed noise of a sky of a few photons per pixel gives
      ! more noise peaks than the Normal tail (about 110 per 2048**2
      ! frame at 4 sigma, 4 at 5 sigma on a sky of 0.4 e-/pixel).

      integer, intent(out) :: nnew

      integer, parameter :: blk = 64
      real(REAL64), dimension(:,:), allocatable :: res, sky, h, sg, ker
      logical, dimension(:,:), allocatable :: mask
      real(REAL64), dimension(:), allocatable :: e2, xs, ys, ss, buf
      integer, dimension(:), allocatable :: idx, tmp
      real(REAL64) :: fw, sig, rk, sk2, perr, s1, s2, v, dmin, a, c, dx, dy, sharp
      integer :: nx, ny, nk, i, j, k, l, i0, j0, ib, jb, nbx, nby, nc, m, nsb
      logical :: ok

      nnew = 0
      nx = frame%naxes(1)
      ny = frame%naxes(2)
      if( resthresh_par > 0 ) then
         resthresh = resthresh_par
      else
         resthresh = max(frame%threshold,normal_deviate(1 / real(nx*ny,REAL64)))
      end if
      if( nstars == 0 ) return
      fw = median(fwhms(1:nstars))
      sig = fw*fwhm2sig
      rk = max(2.0_REAL64,1.5_REAL64*sig)
      nk = floor(rk)
      dmin = 0.75_REAL64*fw

      ! the lowered Gaussian
      allocate(ker(-nk:nk,-nk:nk),mask(-nk:nk,-nk:nk))
      do j = -nk, nk
         do i = -nk, nk
            mask(i,j) = i**2 + j**2 <= rk**2
            ker(i,j) = exp(-(i**2 + j**2) / (2*sig**2))
         end do
      end do
      a = sum(ker,mask) / count(mask)
      where( mask )
         ker = ker - a
      elsewhere
         ker = 0
      end where
      sk2 = sum(ker**2)

      ! the residuals, the sky by clipped means in blocks
      allocate(res(nx,ny),sky(nx,ny))
      res = fr%data - fr%mstar
      nbx = (nx - 1) / blk + 1
      nby = (ny - 1) / blk + 1
      allocate(buf(blk*blk))
      s1 = median(pack(res,fr%valid))
      do jb = 1, nby
         do ib = 1, nbx
            m = 0
            do j = (jb - 1)*blk + 1, min(jb*blk,ny)
               do i = (ib - 1)*blk + 1, min(ib*blk,nx)
                  if( fr%valid(i,j) ) then
                     m = m + 1
                     buf(m) = res(i,j)
                  end if
               end do
            end do
            if( m > blk**2 / 4 ) then
               ! the mean (for the noise model), clipped
               call sky_mean(buf(1:m),fr%gain,fr%rn2,fr%bias,a,s2,v,nsb)
            else
               a = s1
            end if
            sky((ib - 1)*blk + 1:min(ib*blk,nx),(jb - 1)*blk + 1:min(jb*blk,ny)) = a
         end do
      end do

      ! the PSF mismatch by bright isolated stars
      allocate(e2(0))
      do k = 1, nstars
         if( .not. (snr(k) > 100) .or. stars(k)%satrad > 0 ) cycle
         if( iand(stars(k)%flag,FLAG_SATURATED + FLAG_EXCLUDED + FLAG_FAILED) /= 0 ) cycle
         if( neighbour(k,2*stars(k)%fitrad) > 0 ) cycle
         s1 = 0
         s2 = 0
         m = ceiling(fwhms(k))
         do j = nint(stars(k)%y) - m, nint(stars(k)%y) + m
            do i = nint(stars(k)%x) - m, nint(stars(k)%x) + m
               if( i < 1 .or. j < 1 .or. i > nx .or. j > ny ) cycle
               if( .not. fr%valid(i,j) ) cycle
               if( hypot(i - stars(k)%x,j - stars(k)%y) > fwhms(k) ) cycle
               s1 = s1 + (res(i,j) - sky(i,j))**2 - fr%variance(sky(i,j) + fr%mstar(i,j))
               s2 = s2 + fr%mstar(i,j)**2
            end do
         end do
         if( s2 > 0 ) e2 = [e2, s1 / s2]
      end do
      if( size(e2) > 0 ) then
         perr = sqrt(max(median(e2),0.0_REAL64))
      else
         perr = 0.02_REAL64
      end if
      perr = max(perr,0.005_REAL64)
      if( verbose ) write(error_unit,'(a,f0.4,a,i0,a)') ' Residuals: PSF mismatch ', &
           perr,' (by ',size(e2),' stars)'

      ! the amplitude and its significance
      allocate(h(nx,ny),sg(nx,ny))
      h = 0
      sg = 0
      do j = 1 + nk, ny - nk
         do i = 1 + nk, nx - nk
            if( .not. all(fr%valid(i-nk:i+nk,j-nk:j+nk) .or. .not. mask) ) cycle
            s1 = 0
            v = 0
            do l = -nk, nk
               do k = -nk, nk
                  if( .not. mask(k,l) ) cycle
                  s1 = s1 + ker(k,l)*(res(i+k,j+l) - sky(i+k,j+l))
                  v = v + ker(k,l)**2*(fr%variance(sky(i+k,j+l) + fr%mstar(i+k,j+l)) + &
                       (perr*fr%mstar(i+k,j+l))**2)
               end do
            end do
            h(i,j) = s1 / sk2
            sg(i,j) = h(i,j) / (sqrt(v) / sk2)
         end do
      end do

      ! the local maxima
      allocate(xs(0),ys(0),ss(0))
      do j = 1 + nk, ny - nk
         do i = 1 + nk, nx - nk
            if( sg(i,j) < resthresh .or. .not. (h(i,j) > 0) ) cycle
            ok = .true.
            do l = -nk, nk
               do k = -nk, nk
                  if( .not. mask(k,l) .or. (k == 0 .and. l == 0) ) cycle
                  ! ties: the first one in the scan order
                  if( h(i+k,j+l) > h(i,j) .or. (h(i+k,j+l) == h(i,j) .and. &
                       (l < 0 .or. (l == 0 .and. k < 0))) ) ok = .false.
               end do
            end do
            if( .not. ok ) cycle
            ! sharpness as find: the central pixel above the mean of others
            a = 0
            m = 0
            do l = -nk, nk
               do k = -nk, nk
                  if( .not. mask(k,l) .or. (k == 0 .and. l == 0) ) cycle
                  a = a + res(i+k,j+l) - sky(i+k,j+l)
                  m = m + 1
               end do
            end do
            sharp = (res(i,j) - sky(i,j) - a/m) / h(i,j)
            if( sharp < frame%shrplo .or. sharp > frame%shrphi ) cycle
            ! no saturated pixels near
            if( any(fr%satur(max(i-nk-1,1):min(i+nk+1,nx), &
                 max(j-nk-1,1):min(j+nk+1,ny))) ) cycle
            ! the local scatter of the amplitude (2-4 FWHM), observed
            ! relative to predicted: structured residuals near bright
            ! stars (the PSF mismatch) raise the threshold locally
            block
              real(REAL64), dimension(:), allocatable :: obs, pred
              real(REAL64) :: r, q, s3, s4
              integer :: ii, jj, mm, nn
              mm = ceiling(4*fw)
              allocate(obs((2*mm+1)**2),pred((2*mm+1)**2))
              nn = 0
              do jj = max(j - mm,1 + nk), min(j + mm,ny - nk)
                 do ii = max(i - mm,1 + nk), min(i + mm,nx - nk)
                    r = hypot(real(ii - i,REAL64),real(jj - j,REAL64))
                    if( r < 2*fw .or. r > 4*fw .or. .not. (sg(ii,jj) /= 0) ) cycle
                    nn = nn + 1
                    obs(nn) = h(ii,jj)
                    pred(nn) = h(ii,jj) / sg(ii,jj)
                 end do
              end do
              q = 1
              if( nn >= 20 ) then
                 s3 = median(obs(1:nn))
                 s3 = 1.4826_REAL64*median(abs(obs(1:nn) - s3))
                 s4 = median(pred(1:nn))
                 if( s4 > 0 ) q = max(s3 / s4,1.0_REAL64)
              end if
              if( sg(i,j) < resthresh*q ) cycle
            end block
            ! the sub-pixel position
            c = h(i-1,j) - 2*h(i,j) + h(i+1,j)
            dx = 0
            if( c < 0 ) dx = max(min(0.5_REAL64*(h(i-1,j) - h(i+1,j)) / c, &
                 0.5_REAL64),-0.5_REAL64)
            c = h(i,j-1) - 2*h(i,j) + h(i,j+1)
            dy = 0
            if( c < 0 ) dy = max(min(0.5_REAL64*(h(i,j-1) - h(i,j+1)) / c, &
                 0.5_REAL64),-0.5_REAL64)
            xs = [xs, i + dx]
            ys = [ys, j + dy]
            ss = [ss, sg(i,j)]
         end do
      end do
      nc = size(xs)
      if( nc == 0 ) return

      ! far from stars and from each other (the more significant first)
      allocate(idx(nc))
      idx = [(k, k = 1, nc)]
      call msort(idx,ss,tmp)
      block
        real(REAL64), dimension(:), allocatable :: ax, ay
        allocate(ax(0),ay(0))
        do l = 1, nc
           k = idx(l)
           if( near_star(xs(k),ys(k),dmin) ) cycle
           if( size(ax) > 0 ) then
              if( minval(hypot(ax - xs(k),ay - ys(k))) < dmin ) cycle
           end if
           ax = [ax, xs(k)]
           ay = [ay, ys(k)]
        end do
        nnew = size(ax)
        if( nnew > 0 ) call append_stars(ax,ay)
      end block

    end subroutine residual_detect

    real(REAL64) function normal_deviate(p)

      ! the deviate t of the Normal distribution exceeded with the
      ! probability p: erfc(t/sqrt(2))/2 = p (bisection)

      real(REAL64), intent(in) :: p
      real(REAL64) :: a, b, t
      integer :: i

      a = 0
      b = 38
      do i = 1, 100
         t = (a + b) / 2
         if( erfc(t / sqrt(2.0_REAL64)) / 2 > p ) then
            a = t
         else
            b = t
         end if
      end do
      normal_deviate = (a + b) / 2

    end function normal_deviate

    logical function near_star(x,y,d)

      ! a star (current position) closer than d to x,y

      real(REAL64), intent(in) :: x, y, d
      integer :: i, j, i0, j0, m, l

      near_star = .true.
      m = ceiling((d + 2*cell) / cell)
      i0 = min(max(int(x / cell),0),ncx-1)
      j0 = min(max(int(y / cell),0),ncy-1)
      do j = max(j0 - m,0), min(j0 + m,ncy-1)
         do i = max(i0 - m,0), min(i0 + m,ncx-1)
            l = cell_head(1 + i + ncx*j)
            do while( l > 0 )
               if( hypot(stars(l)%x - x,stars(l)%y - y) < d ) return
               if( satrad(l) > 0 .and. hypot(stars(l)%x - x,stars(l)%y - y) < &
                    satrad(l) + d ) return
               l = cell_next(l)
            end do
         end do
      end do
      near_star = .false.

    end function near_star

    subroutine append_stars(x,y)

      ! new stars at x,y: their sky, initial flux, the grid

      real(REAL64), dimension(:), intent(in) :: x, y

      type(StarFit), dimension(:), allocatable :: tmp
      real(REAL64), dimension(:), allocatable :: rtmp
      integer :: n0, n, k

      n0 = nstars
      n = n0 + size(x)
      allocate(tmp(n))
      tmp(1:n0) = stars
      call move_alloc(tmp,stars)
      frame%x = [frame%x, x]
      frame%y = [frame%y, y]
      allocate(rtmp(size(x)))
      rtmp = -1
      frame%sky = [frame%sky, rtmp]
      frame%skyerr = [frame%skyerr, rtmp]
      rtmp = 0
      satrad = [satrad, rtmp]
      fwhms = [fwhms, rtmp]
      snr = [snr, rtmp]
      deallocate(order)
      allocate(order(n))
      nstars = n
      do k = n0 + 1, n
         stars(k)%x = x(k-n0)
         stars(k)%y = y(k-n0)
         stars(k)%flux = 0
         stars(k)%flag = 0
         stars(k)%satrad = 0
         stars(k)%modrad = stars(1)%modrad
         stars(k)%modelled = .false.
      end do
      call set_radii
      call make_grid
      do k = n0 + 1, n
         if( skymode == 'LOCAL' ) then
            call local_sky(stars(k))
         else
            call local_sky(stars(k))
            call annulus_sky(stars(k))
         end if
         call linear_flux(fr,psf,stars(k),skymode == 'LOCAL')
      end do
      call sort_order(stars%flux)

    end subroutine append_stars

    subroutine save_psf(nvar,nshape,status)

      ! The empirical PSF in the format of PSF files: the stamp at the
      ! centre (a constant PSF with an analytic core), or the PSFs
      ! on rings covering the frame (2 rings by 16 angles for the order 0
      ! or 1, 12 rings by 32 angles for 2, plus the rings of the PSF file
      ! of the base inside of the frame); the radius of the field is the distance
      ! from the field centre (NAXIS/2 or --psf-centre) divided by the
      ! scale [pix/deg] (by WCS, --scale, or 1 when undefined: RADCOO is
      ! in pixels then).

      integer, intent(in) :: nvar, nshape
      integer, intent(in out) :: status

      real(REAL64), dimension(:,:,:), allocatable :: st
      real(REAL64), dimension(:,:), allocatable :: sub
      real(REAL64), dimension(:), allocatable :: rr, rc, ac
      real(REAL64) :: scale, rmax, x, y, a
      character(len=200), dimension(3) :: comments
      integer :: nsub, nr, na, k, l, m, ns

      scale = psf%scale
      na = 16
      if( .not. (scale > 0) ) then
         scale = 1
         if( nvar > 0 .or. psf%ptype == PSF_FILE ) write(error_unit,*) &
              'Warning: no scale (WCS) of the frame: RADCOO of the saved ', &
              'PSF in pixels (use it with --scale 1).'
      end if

      if( nvar == 0 .and. psf%ptype /= PSF_FILE ) then
         allocate(rr(0))
      else
         rmax = 0
         do k = 0, 1
            do l = 0, 1
               x = 1 + k*(frame%naxes(1) - 1)
               y = 1 + l*(frame%naxes(2) - 1)
               rmax = max(rmax,hypot(x - psf%xcen,y - psf%ycen))
            end do
         end do
         if( nvar <= 1 ) then
            nr = 2
            na = 16
         else
            nr = 12
            na = 32
         end if
         rr = [(rmax*k/nr / scale, k = 1, nr)]
         ! the rings of the base inside of the frame
         if( psf%ptype == PSF_FILE ) then
            do k = 1, size(psf%rcoo)
               if( psf%rcoo(k) > 0 .and. psf%rcoo(k) < rr(nr) ) then
                  if( all(abs(rr - psf%rcoo(k)) > 1e-9_REAL64) ) rr = [rr, psf%rcoo(k)]
               end if
            end do
            call sort_real(rr)
         end if
      end if

      nr = size(rr)
      m = 1 + nr*na
      allocate(rc(m),ac(m))
      call psf%sampled(psf%xcen,psf%ycen,sub,nsub)
      ns = nsub
      allocate(st(-ns:ns,-ns:ns,m))
      st(:,:,1) = sub
      rc(1) = 0
      ac(1) = 0
      l = 1
      do k = 1, nr
         do m = 1, na
            l = l + 1
            a = 360.0_REAL64*(m - 1) / na
            x = psf%xcen + rr(k)*scale*cos(a / rad)
            y = psf%ycen + rr(k)*scale*sin(a / rad)
            deallocate(sub)
            call psf%sampled(x,y,sub,nsub)
            st(:,:,l) = sub
            rc(l) = rr(k)
            ac(l) = a
         end do
      end do

      write(comments(1),'(2a)') 'Empirical PSF by psfphot from ', &
           trim(shortname(filename))
      if( psf%ptype == PSF_FILE ) then
         write(comments(2),'(3a,i0,a,i0,a)') 'Base: FILE ',trim(shortname(psffile)), &
              ', variation order ',nvar,', ',npsfsel,' PSF stars.'
      else
         write(comments(2),'(3a,i0,a,i0,a)') 'Base: ',trim(psfcore), &
              ', variation order ',nvar,', ',npsfsel,' PSF stars.'
      end if
      write(comments(3),'(a)') 'RADCOO by the centre PSFCEN1, PSFCEN2 [pix], ' // &
           'the scale PSFSCALE [pix/deg].'
      call psf_file_save(psfsave,st,rc,ac,psf%zoom, &
           [character(len=FLEN_KEYWORD) :: 'PSFCEN1','PSFCEN2','PSFSCALE', &
           'FWHM','PSFRAD'], &
           [psf%xcen, psf%ycen, scale, psf%fwhm(psf%xcen,psf%ycen), psfrad], &
           [character(len=FLEN_COMMENT) :: '[pix] field centre, x', &
           '[pix] field centre, y','[pix/deg] scale of RADCOO', &
           '[pix] FWHM at the field centre','[pix] radius of the residual table'], &
           comments,status)
      if( verbose .and. status == 0 ) write(error_unit,'(3a,i0,a)') &
           ' Empirical PSF saved to `',trim(psfsave),"' (",size(rc),' PSFs).'

    end subroutine save_psf

    subroutine sort_real(x)

      real(REAL64), dimension(:), intent(in out) :: x
      real(REAL64) :: t
      integer :: i, j

      do i = 2, size(x)
         t = x(i)
         j = i - 1
         do while( j >= 1 )
            if( x(j) <= t ) exit
            x(j+1) = x(j)
            j = j - 1
         end do
         x(j+1) = t
      end do

    end subroutine sort_real

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
