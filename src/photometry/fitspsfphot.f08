!
!  fitspsfphot - FITS input and output of the PSF photometry
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

module fitspsfphot

  use titsio
  use iso_fortran_env

  implicit none

  ! keywords of the PSF photometry extension
  character(len=*), parameter :: FITS_KEY_PSFTYPE = 'PSFTYPE'
  character(len=*), parameter :: FITS_KEY_PSFFILE = 'PSFFILE'
  character(len=*), parameter :: FITS_KEY_FITRAD = 'FITRAD'

  ! the columns
  character(len=*), parameter :: FITS_COL_XERR = 'XERR'
  character(len=*), parameter :: FITS_COL_YERR = 'YERR'
  character(len=*), parameter :: FITS_COL_CHI2 = 'CHI2'
  character(len=*), parameter :: FITS_COL_NPIX = 'NPIX'
  character(len=*), parameter :: FITS_COL_PSFFLAG = 'FLAG'

  type :: PsfFrame
     real(REAL64), dimension(:,:), allocatable :: data
     integer, dimension(2) :: naxes
     integer :: bitpix
     ! saturation, gain [e-/ADU], read noise [ADU], bias [ADU]
     real(REAL64) :: saturate = -1, gain = 1, readns = 0, bias = 0
     character(len=FLEN_VALUE) :: bunit = ''
     character(len=FLEN_CARD) :: gaininfo = '', rninfo = ''
     ! the reference pixel and the scale [pix/deg]
     real(REAL64), dimension(2) :: crpix = 0
     real(REAL64) :: scale = -1
     ! detected stars (FIND), their sky (APERPHOT)
     real(REAL64), dimension(:), allocatable :: x, y, sky, skyerr
     real(REAL64) :: fwhm = -1, lobad = -huge(1.0_REAL64), hibad = huge(1.0_REAL64)
     real(REAL64) :: aphwhm = -1
     real(REAL64), dimension(2) :: ring = -1
     logical :: have_aper = .false.
  end type PsfFrame

  type :: PsfHeader
     character(len=FLEN_VALUE) :: psftype = '', skymode = ''
     character(len=FLEN_FILENAME) :: psffile = ''
     real(REAL64) :: fwhm = -1, width = -1, beta = -1, q = 1, pa = 0
     real(REAL64) :: fitrad = -1, gain = 1, readns = 0, bias = 0, saturate = -1
     real(REAL64), dimension(2) :: psfcen = 0
     real(REAL64) :: scale = -1
     integer :: zoom = 0, nshape = 0, npass = 0
     logical :: fwhm_fitted = .false.
  end type PsfHeader

contains

  subroutine fits_psfphot_read(filename,fkeys,frame,status)

    ! fkeys: keywords of saturation, read noise (ADU), gain

    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: fkeys
    type(PsfFrame), intent(out) :: frame
    integer, intent(in out) :: status

    type(fitsfiles) :: fits
    integer, parameter :: extver = 0
    real(REAL64), dimension(2,2) :: cd
    real(REAL64), dimension(2) :: cdelt
    real(REAL64) :: rdnoise
    integer :: naxis, nrows, i, l, srows, frow, xcol, ycol, scol, ecol
    logical :: anyf, have_rdnoise

    if( status /= 0 ) return

    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read an image in the file `',trim(filename),"'."
       return
    end if

    call fits_get_img_dim(fits,naxis,status)
    if( naxis /= 2 .and. status == 0 ) then
       write(error_unit,*) 'Error in psfphot: Only 2D frames are supported.'
       status = 1
       goto 666
    end if
    call fits_get_img_size(fits,2,frame%naxes,status)
    call fits_get_img_type(fits,frame%bitpix,status)
    if( status /= 0 ) goto 666

    allocate(frame%data(frame%naxes(1),frame%naxes(2)))
    call fits_read_image(fits,0,0.0_REAL64,frame%data,anyf,status)
    if( status /= 0 ) then
       write(error_unit,*) trim(filename),": Failed to read data."
       goto 666
    end if

    ! detector parameters
    call fits_write_errmark
    call fits_read_key(fits,fkeys(1),frame%saturate,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       status = 0
       if( frame%bitpix > 0 ) then
          frame%saturate = 2.0_REAL64**frame%bitpix - 1
       else
          frame%saturate = huge(frame%saturate)
       end if
    end if

    call fits_write_errmark
    call fits_read_key(fits,FITS_KEY_BUNIT,frame%bunit,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       status = 0
       frame%bunit = ''
    end if

    call fits_write_errmark
    call fits_read_key(fits,fkeys(3),frame%gain,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       status = 0
       frame%gain = 1
       frame%gaininfo = 'no '//trim(fkeys(3))//' keyword, 1 assumed'
    else
       frame%gaininfo = trim(fkeys(3))//' keyword'
    end if
    if( frame%bunit == 'COUNTS' ) then
       ! converted to electrons by phcorr
       frame%gain = 1
       frame%gaininfo = 'BUNIT = COUNTS (electrons)'
    end if
    if( .not. (frame%gain > 0) ) then
       write(error_unit,*) 'Warning: gain <= 0 in the header, 1 used.'
       frame%gain = 1
    end if

    ! read noise: READNS [ADU] (Munipack), or RDNOISE [e-] (artificial)
    call fits_write_errmark
    call fits_read_key(fits,fkeys(2),frame%readns,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       status = 0
       call fits_write_errmark
       call fits_read_key(fits,'RDNOISE',rdnoise,status)
       have_rdnoise = status == 0
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
       end if
       if( have_rdnoise ) then
          frame%readns = rdnoise / frame%gain
          frame%rninfo = 'RDNOISE keyword [e-] / gain'
       else
          frame%readns = 0
          frame%rninfo = 'no read noise keyword, 0 assumed'
       end if
    else
       frame%rninfo = trim(fkeys(2))//' keyword [ADU]'
    end if

    ! the bias level (a pedestal), only for the noise model
    call fits_write_errmark
    call fits_read_key(fits,'BIASLEV',frame%bias,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       status = 0
       frame%bias = 0
    end if

    ! the reference pixel and the scale
    call fits_write_errmark
    call fits_read_key(fits,'CRPIX1',frame%crpix(1),status)
    call fits_read_key(fits,'CRPIX2',frame%crpix(2),status)
    if( status /= 0 ) then
       call fits_clear_errmark
       status = 0
       frame%crpix = frame%naxes / 2
    end if
    call fits_write_errmark
    call fits_read_key(fits,'CD1_1',cd(1,1),status)
    call fits_read_key(fits,'CD1_2',cd(1,2),status)
    call fits_read_key(fits,'CD2_1',cd(2,1),status)
    call fits_read_key(fits,'CD2_2',cd(2,2),status)
    if( status == 0 ) then
       frame%scale = 1 / sqrt(abs(cd(1,1)*cd(2,2) - cd(1,2)*cd(2,1)))
    else
       call fits_clear_errmark
       status = 0
       call fits_write_errmark
       call fits_read_key(fits,'CDELT1',cdelt(1),status)
       call fits_read_key(fits,'CDELT2',cdelt(2),status)
       if( status == 0 ) then
          frame%scale = 1 / sqrt(abs(cdelt(1)*cdelt(2)))
       else
          call fits_clear_errmark
          status = 0
          frame%scale = -1
       end if
    end if

    ! detected stars
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,FINDEXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       write(error_unit,*) "Error: ",trim(FINDEXTNAME)//" extension not found."
       write(error_unit,*) "       Has been stars detected by `munipack find " &
            //trim(filename)//"' ?"
       goto 666
    end if

    call fits_get_num_rows(fits,nrows,status)
    call fits_read_key(fits,FITS_KEY_LOWBAD,frame%lobad,status)
    call fits_read_key(fits,FITS_KEY_HIGHBAD,frame%hibad,status)
    call fits_read_key(fits,FITS_KEY_FWHM,frame%fwhm,status)
    call fits_get_colnum(fits,.true.,FITS_COL_X,xcol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_Y,ycol,status)
    if( status /= 0 ) goto 666

    allocate(frame%x(nrows),frame%y(nrows),frame%sky(nrows),frame%skyerr(nrows))
    frame%sky = -1
    frame%skyerr = -1
    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows-1,nrows)
       call fits_read_col(fits,xcol,frow,0.0_REAL64,frame%x(i:l),anyf,status)
       call fits_read_col(fits,ycol,frow,0.0_REAL64,frame%y(i:l),anyf,status)
    end do
    if( status /= 0 ) goto 666

    ! the sky by aperture photometry, if available
    call fits_write_errmark
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,APEREXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       call fits_clear_errmark
       status = 0
       frame%have_aper = .false.
    else if( status == 0 ) then
       call fits_get_num_rows(fits,l,status)
       call fits_read_key(fits,FITS_KEY_HWHM,frame%aphwhm,status)
       call fits_read_key(fits,trim(FITS_KEY_ANNULUS)//'1',frame%ring(1),status)
       call fits_read_key(fits,trim(FITS_KEY_ANNULUS)//'2',frame%ring(2),status)
       call fits_get_colnum(fits,.true.,FITS_COL_SKY,scol,status)
       call fits_get_colnum(fits,.true.,FITS_COL_SKYERR,ecol,status)
       if( status /= 0 ) goto 666
       if( l == nrows ) then
          frame%have_aper = .true.
          do i = 1, nrows, srows
             frow = i
             l = min(i+srows-1,nrows)
             call fits_read_col(fits,scol,frow,-1.0_REAL64,frame%sky(i:l),anyf,status)
             call fits_read_col(fits,ecol,frow,-1.0_REAL64,frame%skyerr(i:l),anyf,status)
          end do
       else
          write(error_unit,*) 'Warning: ',trim(APEREXTNAME), &
               ' and ',trim(FINDEXTNAME),' differ in rows, ', &
               trim(APEREXTNAME),' ignored.'
       end if
    end if

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine fits_psfphot_read


  subroutine fits_psfphot_save(filename,output,hdr,ttype,tunit,table,npix, &
       flags,status)

    ! The table columns are real: X, Y, XERR, YERR, PSFCOUNT, PSFCOUNTERR,
    ! SKY, SKYERR, CHI2, SHARP, ROUND, and integer: NPIX, FLAG

    character(len=*), intent(in) :: filename, output
    type(PsfHeader), intent(in) :: hdr
    character(len=*), dimension(:), intent(in) :: ttype, tunit
    real(REAL64), dimension(:,:), intent(in) :: table
    integer, dimension(:), intent(in) :: flags, npix
    integer, intent(in out) :: status

    type(fitsfiles) :: fits
    integer, parameter :: extver = 0
    character(len=FLEN_VALUE), dimension(:), allocatable :: tform
    character(len=FLEN_COMMENT) :: com
    integer :: hdutype, i, nrows, ncols

    if( status /= 0 ) return

    if( output == '' ) then
       call fits_open_file(fits,filename,FITS_READWRITE,status)
       if( status /= 0 ) then
          write(error_unit,*) &
               'Error: failed to open the file `',trim(filename),"' for a table update."
          return
       end if
    else
       call fits_precopy_file(fits,filename,output,FITS_READWRITE,.true.,status)
       if( status /= 0 ) then
          write(error_unit,*) &
               'Error: failed to create the file `',trim(output),"'."
          return
       end if
    end if

    ! already presented ? remove it !
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,PSFEXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       status = 0
    else
       call fits_delete_hdu(fits,hdutype,status)
    end if
    if( status /= 0 ) goto 666

    ncols = size(ttype)
    allocate(tform(ncols))
    tform = '1D'
    tform(ncols-1:ncols) = '1J'
    nrows = size(table,1)
    call fits_insert_btbl(fits,nrows,ttype,tform,tunit,PSFEXTNAME,status)

    call fits_write_key(fits,FITS_KEY_PSFTYPE,hdr%psftype, &
         'PSF: GAUSS, MOFFAT or FILE',status)
    if( hdr%psffile /= '' ) &
         call fits_write_key(fits,FITS_KEY_PSFFILE,hdr%psffile,'PSF file',status)
    if( hdr%fwhm_fitted ) then
       com = '[pix] FWHM of PSF (fitted)'
    else
       com = '[pix] FWHM of PSF'
    end if
    if( hdr%psftype == 'FILE' ) com = '[pix] FWHM of PSF at the field centre'
    call fits_write_key(fits,FITS_KEY_FWHM,hdr%fwhm,-5,com,status)
    call fits_write_key(fits,FITS_KEY_HWHM,hdr%fwhm/2,-5, &
         '[pix] half width at half of maximum',status)
    if( hdr%psftype == 'GAUSS' ) then
       call fits_write_key(fits,'PSFSIGMA',hdr%width,-6, &
            '[pix] Gaussian sigma (major axis)',status)
    else if( hdr%psftype == 'MOFFAT' ) then
       call fits_write_key(fits,'PSFALPHA',hdr%width,-6, &
            '[pix] Moffat alpha (major axis)',status)
       call fits_write_key(fits,'PSFBETA',hdr%beta,-6,'Moffat beta',status)
    end if
    if( hdr%psftype /= 'FILE' ) then
       call fits_write_key(fits,'PSFAXRAT',hdr%q,-5, &
            'axis ratio of PSF (minor/major)',status)
       call fits_write_key(fits,'PSFPA',hdr%pa,-3, &
            '[deg] position angle of the major axis from x',status)
       if( hdr%fwhm_fitted ) call fits_write_key(fits,'PSFNSTAR',hdr%nshape, &
            'stars used to fit the PSF shape',status)
    else
       call fits_write_key(fits,'PSFCEN1',hdr%psfcen(1),-5, &
            '[pix] field centre of the PSF grid, x',status)
       call fits_write_key(fits,'PSFCEN2',hdr%psfcen(2),-5, &
            '[pix] field centre of the PSF grid, y',status)
       call fits_write_key(fits,'PSFSCALE',hdr%scale,-10, &
            '[pix/deg] scale for the PSF grid',status)
    end if
    call fits_write_key(fits,'PSFZOOM',hdr%zoom,'PSF samples per pixel',status)
    if( hdr%fitrad > 0 ) then
       call fits_write_key(fits,FITS_KEY_FITRAD,hdr%fitrad,-4, &
            '[pix] fitting radius',status)
    else
       call fits_write_key(fits,FITS_KEY_FITRAD,-1, &
            'fitting radius 1.5*FWHM (min. 3) of local PSF',status)
    end if
    call fits_write_key(fits,'SKYMODE',hdr%skymode,'sky: FIXED or LOCAL',status)
    call fits_write_key(fits,'NPASS',hdr%npass,'passes over stars',status)
    call fits_write_key(fits,FITS_KEY_GAIN,hdr%gain,-6, &
         '[e-/ADU] gain used for weights',status)
    call fits_write_key(fits,FITS_KEY_READNS,hdr%readns,-6, &
         '[ADU] read noise used for weights',status)
    call fits_write_key(fits,'BIASLEV',hdr%bias,-6, &
         '[ADU] bias level used for weights',status)
    call fits_write_key(fits,FITS_KEY_SATURATE,hdr%saturate,-8, &
         '[ADU] saturation (pixels excluded)',status)

    call fits_write_comment(fits,'PSFCOUNT is the total flux of the PSF ' // &
         '(the same units as APCOUNT, GCOUNT).',status)
    call fits_write_comment(fits,'FLAG (sum): 1 saturated core, ' // &
         '2 pixels excluded (bad, saturated, edge),',status)
    call fits_write_comment(fits,'  4 position fixed (faint or run away), ' // &
         '8 not converged, 16 failed (no valid flux),',status)
    call fits_write_comment(fits,'  32 a neighbour within the fitting radius, ' // &
         '64 a saturated star added (not found by find).',status)
    call fits_write_comment(fits,'SHARP: relative residual in the core ' // &
         '(>0 sharper than PSF), ROUND: ellipticity data - model.',status)

    if( nrows > 0 ) then
       do i = 1, ncols - 2
          call fits_write_col(fits,i,1,table(:,i),status)
       end do
       call fits_write_col(fits,ncols-1,1,npix,status)
       call fits_write_col(fits,ncols,1,flags,status)
    end if

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine fits_psfphot_save


  subroutine fits_psfphot_residual(filename,resname,residual,status)

    ! the residual image (data - stars), sky included, with the header
    ! of the original frame

    character(len=*), intent(in) :: filename, resname
    real(REAL64), dimension(:,:), intent(in) :: residual
    integer, intent(in out) :: status

    type(fitsfiles) :: fits, rfits
    integer :: i

    if( status /= 0 ) return

    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) return
    call fits_create_file(rfits,resname,status)
    if( status /= 0 ) then
       write(error_unit,*) "Error: failed to create the residual file `", &
            trim(resname),"' (exists?)."
       call fits_close_file(fits,status)
       return
    end if
    call fits_copy_header(fits,rfits,status)
    call fits_close_file(fits,status)

    call fits_update_key(rfits,'BITPIX',-32,'',status)
    do i = 1, 2
       call fits_write_errmark
       if( i == 1 ) call fits_delete_key(rfits,'BSCALE',status)
       if( i == 2 ) call fits_delete_key(rfits,'BZERO',status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
       end if
    end do
    call fits_write_comment(rfits,'Residuals of PSF photometry: the frame '// &
         'with the fitted stars subtracted.',status)
    call fits_write_image(rfits,0,real(residual),status)
    call fits_close_file(rfits,status)
    call fits_report_error(error_unit,status)

  end subroutine fits_psfphot_residual

end module fitspsfphot
