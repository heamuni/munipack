!
!  Phcorr - photometric corrections
!  Copyright (C) 2012-4, 2016-20, 2023  Filip Hroch, Masaryk University, Brno, CZ
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

program phcorr

  use fitscorr
  use oakleaf
  use titsio
  use iso_fortran_env

  implicit none

  character(len=FLEN_FILENAME) :: sciname, output
  character(len=FLEN_FILENAME) :: flatname = ''
  character(len=FLEN_FILENAME) :: biasname=''
  character(len=FLEN_FILENAME) :: darkname=''
  character(len=FLEN_FILENAME) :: maskname=''

  character(len=4*FLEN_FILENAME) :: key, val, record
  character(len=80) :: msg

  integer :: eq, stat
  logical :: verbose = .false.       ! vebosity
  logical :: normalise = .false.     ! normalise flat prior to use
  real    :: gain_user = 0           ! gain
  logical :: gain_set = .false.
  real    :: tol_exptime = 1e-6      ! tolerances
  real    :: tol_temperat = 1.0
  real    :: xdark = -1.0            ! dark time factor
  real    :: flat_mean = 1.0, flat_err = 0.0 ! flat parameters
  integer :: bitpix = -32            ! bitpix
  integer :: box = 5                 ! box size for median mask
  logical :: ex, exitus = .true.
  logical :: overwrite = .false.     ! do not overwrite files
  logical :: gain_warning = .false.  ! reports unset gain on finish

  character(len=FLEN_KEYWORD) :: tmask = 'MEAN'
  character(len=FLEN_KEYWORD), dimension(6) :: keys = [ &
       character(len=FLEN_KEYWORD) :: &
       FITS_KEY_DATEOBS, &
       FITS_KEY_EXPTIME, &
       FITS_KEY_FILTER, &
       FITS_KEY_SATURATE, &
       FITS_KEY_TEMPERATURE, &
       FITS_KEY_GAIN ]

  ! FITS wrappers
  type(CorrFits), allocatable :: flatfits, darkfits, biasfits, maskfits

  do
     read(*,'(a)',iostat=stat,iomsg=msg) record
     if( stat == IOSTAT_END ) exit
     if( stat > 0 ) then
        write(error_unit,*) trim(msg)
        error stop 'Error: an general input error.'
     end if

     eq = index(record,'=')
     if( eq == 0 ) error stop 'Error: Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'VERBOSE' ) then
        read(val,*) verbose
     endif

     if( key == 'OVERWRITE' ) then
        read(val,*) overwrite
     endif

     if( key == 'BITPIX' ) then
        read(val,*) bitpix
     endif

     if( key == 'FITS_KEY_DATEOBS' ) then
        read(val,*) keys(1)
     endif

     if( key == 'FITS_KEY_EXPTIME' ) then
        read(val,*) keys(2)
     endif

     if( key == 'FITS_KEY_FILTER' ) then
        read(val,*) keys(3)
     endif

     if( key == 'FITS_KEY_SATURATE' ) then
        read(val,*) keys(4)
     endif

     if( key == 'FITS_KEY_TEMPERATURE' ) then
        read(val,*) keys(5)
     endif

     if( key == 'FITS_KEY_GAIN' ) then
        read(val,*) keys(6)
     endif

     if( key == 'GAIN' ) then
        read(val,*) gain_user
        if( .not. (gain_user > 0) ) &
             error stop 'Error: Assertion: gain > 0 unsatisfied.'
        gain_set = .true.
     endif

     if( key == 'MUNIPACK_TEMPERATURE_TOLERANCE' ) then
        read(val,*) tol_temperat
     endif

     if( key == 'MUNIPACK_EXPTIME_TOLERANCE' ) then
        read(val,*) tol_exptime
     endif

     if( key == 'NORMALISE' ) then
        read(val,*) normalise
     endif

     if( key == 'FLAT' ) then

        read(val,*) flatname

        if( verbose ) write(error_unit,'(a)',advance="no") &
             "Flat-field frame: "//trim(flatname)//", "

        allocate(flatfits)
        call flatfits%Load(flatname,keys)
        if( .not. flatfits%status ) &
             error stop 'Error: Failed to load flat-field frame.'

        if( verbose ) &
             write(error_unit,'(2a)') ' filter = ',trim(flatfits%filter)

        if( flatfits%filter == '' ) &
             write(error_unit,*) &
             "Warning: a photometric filter of flat-field unknown."

        if( .not. all(flatfits%image > 0) ) &
           error stop "Error: any value in the flat-field is non-positive."

        if( normalise .or. flatfits%bitpix > 0 ) then

           ! average over the flat
           call rmean(pack(flatfits%image, flatfits%image > 0),flat_mean,flat_err)

           if( verbose ) &
                write(error_unit,*) 'Flat-frame statistics: mean=', &
                flat_mean,'+-',flat_err

           if( flat_mean < epsilon(flat_mean) ) &
                error stop 'Error: Suspicious value(s) of mean of flat-field.'
        end if

     end if

     if( key == 'BIAS' ) then

        read(val,*) biasname
        if( verbose ) write(error_unit,'(a)',advance="no") &
             "Bias frame: "//trim(biasname)//", "

        allocate(biasfits)
        call biasfits%Load(biasname,keys)
        if( .not. biasfits%status ) &
             error stop 'Error: Failed to load bias frame.'

        if( verbose ) write(error_unit,'(a,1pg0.3,a,0pf0.1)') &
             ' exp.time[s] = ',biasfits%exptime,', T[degC] = ',biasfits%temper
     end if

     if( key == 'DARK' ) then

        read(val,*) darkname

        if( verbose ) write(error_unit,'(a)',advance="no") &
             "Dark frame: "//trim(darkname)//", "

        allocate(darkfits)
        call darkfits%Load(darkname,keys)
        if( .not. darkfits%status ) error stop 'Error: Failed to load dark frame.'

        if( verbose ) write(error_unit,'(a,1pg0.3,a,0pf0.1)') &
             ' exp.time[s] = ',darkfits%exptime,', T[degC] = ',darkfits%temper

        if( .not. darkfits%exptime_set .and. verbose ) &
             write(error_unit,*) "Warning: Exposure time for dark is undefined."

        if( darkfits%temper < -273 .and. verbose ) &
             write(error_unit,*) "Warning: Temperature of dark is undefined."

     end if

     if( key == 'XDARK' ) then
        read(val,*) xdark
     end if

     if( key == 'MASK' ) then

        read(val,*) maskname
        if( verbose ) write(error_unit,*) "Mask frame: ",trim(maskname)

        allocate(maskfits)
        call maskfits%Load(maskname)
        if( .not. maskfits%status ) &
             error stop 'Error: Failed to load the mask frame.'

     end if

     if( key == 'XMASK' ) then
        read(val,*) tmask
     end if

     if( key == 'BOX' ) then
        read(val,*) box
     end if

     if( key == 'FILE' ) then

        read(val,*) sciname,output
        ex = corrector(sciname,output)
        exitus = exitus .and. ex

     end if

  end do

  if( gain_warning ) write(error_unit,*) &
       "Warning: Gain is unset for one or more frames (use -gain)."

  if( exitus ) then
     stop 0
  else
     error stop 'Error: an error occurred during phcorr session.'
  end if

contains

  logical function corrector(sciname,output)

    character(len=*), parameter :: phid = 'PHCORR'

    character(len=*), intent(in) :: sciname, output
    character(len=FLEN_CARD) :: buf
    character(len=4*FLEN_FILENAME) :: flattext, darktext, biastext
    real, dimension(:,:), allocatable :: FLAT, BIAS, DARK, MASK, OUT, &
         DFLAT, DBIAS, DDARK, DOUT, VAR
    type(CorrFits) :: scifits
    type(fitsfiles) :: fits, cfits
    character(len=80) :: msg
    integer :: stat,status,hdutype,hdunum,chdu
    real :: t, xerr, xmask, saturate, gain, gain2, maxout, rnadu, qnoise, tb
    logical :: rnset, electrons
    integer :: n, m

    corrector = .false.
    call scifits%Load(sciname,keys)
    if( .not. scifits%status ) then
       write(error_unit,*) &
            "Warning: Failed to read `",trim(sciname),"'. Skippy."
       return
    end if

    if( gain_set ) then
       gain = gain_user
    else
       if( scifits%gain_set ) then
          gain = scifits%gain
       else
          gain_warning = .true.
          gain = 1
       end if
    end if
    gain2 = gain**2

    if( verbose ) &
         write(error_unit,'(1x,f0.2,a)',advance="no") gain,' * ('//trim(sciname)

    if( .not. (gain > 0) ) then
       write(error_unit,*) "Warning: assertion gain > 0 unsatisfied. Skip."
       return
    end if

    n = scifits%naxes(1)
    m = scifits%naxes(2)
    allocate(OUT(n,m),DOUT(n,m),FLAT(n,m),DFLAT(n,m),BIAS(n,m),DBIAS(n,m), &
         DARK(n,m),DDARK(n,m),MASK(n,m),VAR(n,m),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) 'allocate() in phcorr: ',trim(msg)
       error stop 'Error: an insufficient memory.'
    end if

    t = 0.0
    FLAT = 1.0
    DFLAT = 0.0
    BIAS = 0.0
    DBIAS = 0.0
    DARK = 0.0
    DDARK = 0.0

    if( allocated(flatfits) ) then
       if( .not. all(flatfits%naxes == scifits%naxes) )then
          write(error_unit,*) &
               "Warning: Flat-field size does not match corrected frame `",&
               trim(sciname),"', Skip."
          return
       end if

       if( verbose .and. flatfits%filter /= '' .and. scifits%filter /= '' )then
          if( flatfits%filter /= scifits%filter ) &
               write(error_unit,*) "Warning: filters does not match."
       end if

       FLAT =  flatfits%image / flat_mean
       DFLAT = flatfits%imgerr / flat_mean

    endif

    if( allocated(biasfits) ) then
       if( .not. all(biasfits%naxes == scifits%naxes) )then
          write(error_unit,*) &
               "Warning: Bias size does not match corrected frame `",&
               trim(sciname),"', Skip."
          return
       end if
       BIAS = biasfits%image
       DBIAS = biasfits%imgerr

    endif

    if( allocated(maskfits) ) then
       if( .not. all(maskfits%naxes == scifits%naxes) )then
          write(error_unit,*) &
               "Warning: Mask size does not match corrected frame `",&
               trim(sciname),"', Skip."
          return
       end if
    endif

    if( allocated(darkfits) ) then
       if( .not. all(darkfits%naxes == scifits%naxes) )then
          write(error_unit,*) &
               "Warning: Dark sizes does not match corrected frame `",&
               trim(sciname),"', Skip."
          return
       end if

       if( verbose .and. darkfits%temper > -273.0 .and. &
            scifits%temper > -273.0 )then
          if( abs(darkfits%temper - scifits%temper) > tol_temperat ) &
               write(error_unit,*) &
               "Warning: Temperature difference limit reached."
       end if

       if( xdark > 0.0 ) then
          t = xdark
       else if( darkfits%exptime_set .and. scifits%exptime_set .and. &
            darkfits%exptime > 0 .and. scifits%exptime > 0 ) then
          t = scifits%exptime / darkfits%exptime
       else
          t = 1
       end if

       DARK = darkfits%image
       DDARK = darkfits%imgerr

    endif

    if( verbose ) then

       flattext = ''
       biastext = ''
       darktext = ''

       if( flatname /= '' ) write(flattext,'(3a,1pg0.6,a)') &
            '/ (',trim(flatname),'/',flat_mean,')'

       if( biasname /= '' ) &
            biastext = " - "//trim(biasname)

       if( darkname /= '' ) &
            write(darktext,'(a,f0.3,2a)') ' - ',t,'*',trim(darkname)

       write(error_unit,'(a)',advance="no") &
            trim(darktext)//trim(biastext)//') '//trim(flattext)//' ='
    end if

    ! the corrections, we are sure that FLAT > 0
    OUT = gain*(scifits%image - (t*DARK + BIAS)) / FLAT

    ! The read noise of the frame [ADU]: READNS [ADU], or RDNOISE [e-]
    ! (artificial) divided by the gain
    rnset = scifits%readns_set .or. scifits%rdnoise_set
    if( scifits%readns_set ) then
       rnadu = scifits%readns
    else if( scifits%rdnoise_set ) then
       rnadu = scifits%rdnoise / gain
    else
       rnadu = 0
    end if
    ! the quantisation noise of integer frames [ADU**2]
    if( scifits%bitpix > 0 ) then
       qnoise = 1.0 / 12
    else
       qnoise = 0
    end if
    ! The error of the bias enters by (1 - t) when the dark was
    ! bias-subtracted by the same master bias (its STDERR has the scatter
    ! of the darks only): OUT = g*(S - t*D_raw - (1 - t)*B).
    tb = 1
    if( allocated(darkfits) .and. allocated(biasfits) ) then
       if( darkfits%darkbias /= '' .and. &
            basename(darkfits%darkbias) == basename(biasname) ) tb = 1 - t
    end if

    ! The variance of electrons before flat-fielding (first order, Merline
    ! & Howell 1995, Newberry 1991): the Poisson noise of all detected
    ! electrons (the sky, the stars and the dark current, S - B) by a model
    ! of the pixel, the mean of its eight neighbours (not of the pixel
    ! itself, so the weights 1/DOUT**2 do not depend on its noise and an
    ! inverse-variance weighted mean is not biased; not floored at one
    ! electron), the read noise, the quantisation noise of integer frames,
    ! the errors of the bias and the dark frames. The flat-field divides
    ! the noise and adds its own relative error. DOUT > 0 for every pixel.
    call neighbour_mean(scifits%image - BIAS,VAR)
    VAR = gain*max(VAR,0.0) + gain2*(rnadu**2 + qnoise + (tb*DBIAS)**2 + &
         (t*DDARK)**2)
    DOUT = sqrt(max(VAR / FLAT**2 + (OUT*DFLAT/FLAT)**2,1e-4*gain2))
    ! YES, this small piece of code is the full correction!!!
    ! All other staff is just a practise in reading and writing of FITS files.

    ! masking
    if( allocated(maskfits) ) then

       MASK = maskfits%image

       if( tmask == "MEAN" .or. tmask == "MEDIAN" ) then
          ! median of the image
          call rmean(pack(OUT,OUT>0),xmask,xerr)
          if( verbose ) write(error_unit,*) 'Mask frame statistics: mean=',&
               xmask,'+-',xerr
       else ! if( tmask == "ZERO" ) then
          xmask = 0
          xerr = 0
       end if

       where( MASK < epsilon(MASK) )
          OUT = xmask
          DOUT = xerr
       end where

       if( tmask == 'MEDIAN' ) then

          block
            integer :: i,j,h,i1,i2,j1,j2
            h = max(box / 2, 1)

            do i = 1, size(MASK,1)
               do j = 1, size(MASK,2)
                  if( MASK(i,j) < epsilon(MASK) ) then
                     i1 = max(1,i-h)
                     i2 = min(i+h,size(MASK,1))
                     j1 = max(1,j-h)
                     j2 = min(j+h,size(MASK,2))
                     if( count(MASK(i1:i2,j1:j2) > 0) > 0 ) then
                        call rmean(pack(OUT(i1:i2,j1:j2), &
                             MASK(i1:i2,j1:j2) > 0),OUT(i,j),DOUT(i,j))
                     end if
                  end if
               end do
            end do
          end block
       end if

    end if

    ! write out corrected image
    status = 0
    call fits_open_image(fits,sciname,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) 'Error: failed to read the file `',trim(sciname),"'."
       corrector = .false.
       return
    end if

    call fits_create_scratch(cfits,status)
    call fits_get_num_hdus(fits,hdunum,status)
    if( status /= 0 ) goto 666

    block
      character(len=FLEN_VALUE) :: extname
      logical :: found
      found = .false.

      do chdu = 1, hdunum
         call fits_movabs_hdu(fits,chdu,hdutype,status)

         call fits_write_errmark
         call fits_read_key(fits,'EXTNAME',extname,status)
         if( status == FITS_KEYWORD_NOT_FOUND ) then
            call fits_clear_errmark
            status = 0
            extname = ''
         end if

         ! we are looking for the very first IMAGE HDU to apply the corrections
         if( .not. found .and. hdutype == FITS_IMAGE_HDU ) then
            found = .true.

            call fits_copy_header(fits,cfits,status)

            call fits_update_key(cfits,'BITPIX',bitpix,'',status)
            if( bitpix > 0 ) then
               maxout = 2.0**bitpix - 1
               OUT = max(0.0,min(OUT,maxout))
               DOUT = max(1.0,min(DOUT,maxout))
            else
               call fits_write_errmark
               call fits_delete_key(cfits,'BSCALE',status)
               call fits_delete_key(cfits,'BZERO',status)
               if( status == FITS_KEYWORD_NOT_FOUND ) then
                  call fits_clear_errmark
                  status = 0
               end if
            end if

            if( scifits%saturate_set ) then
               saturate = scifits%get_saturate(gain)
               call fits_update_key(cfits,keys(4),saturate,-5, &
                    '[counts] saturation level',status)
            end if

            ! set units of images; the noise model of the output (GAIN,
            ! READNS, RDNOISE, BIASLEV) in its units, so find, aphot and
            ! psfphot read a correct one: the values in electrons have the
            ! gain 1, the read noise in electrons, no bias pedestal
            electrons = gain_set .or. scifits%gain_set
            if( electrons ) then
               call fits_update_key(cfits,FITS_KEY_BUNIT,'COUNTS', &
                    'captured electrons',status)
               call fits_update_key(cfits,keys(6),1.0,5, &
                    '[e-/ADU] the values are electrons (by PHCGAIN)',status)
               call fits_update_key(cfits,'PHCGAIN',gain,-7, &
                    '[e-/ADU] gain of the conversion to electrons',status)
               write(buf,'(g0.3)') gain
               call fits_write_history(cfits,phid//" gain: "//trim(buf),status)
            else
               call fits_update_key(cfits,FITS_KEY_BUNIT,'ADU', &
                    'values directly provided by camera',status)
            end if
            if( rnset ) then
               if( electrons ) then
                  call fits_update_key(cfits,FITS_KEY_READNS,gain*rnadu,-7, &
                       '[e-] read noise (the values are electrons)',status)
                  call fits_update_key(cfits,'RDNOISE',gain*rnadu,-7, &
                       '[e-] read-out noise',status)
               else
                  call fits_update_key(cfits,FITS_KEY_READNS,rnadu,-7, &
                       '[ADU] read noise',status)
               end if
            end if
            ! the pedestal is subtracted with the bias (or with a dark
            ! including it), else it is in the units of the output
            call fits_write_errmark
            if( biasname /= '' .or. darkname /= '' ) then
               call fits_delete_key(cfits,'BIASLEV',status)
            else
               block
                 real :: biaslev
                 call fits_read_key(cfits,'BIASLEV',biaslev,status)
                 if( status == 0 ) call fits_update_key(cfits,'BIASLEV', &
                      gain*biaslev,-7,'bias level (pedestal) in the units',status)
               end block
            end if
            if( status == FITS_KEYWORD_NOT_FOUND ) then
               call fits_clear_errmark
               status = 0
            end if

            if( flatname /= '' ) then
               write(buf,'(1pg0.6)') flat_mean
               call fits_write_history(cfits,phid//" flat-field: '" //  &
                    trim(flatname)//"' / "//trim(buf),status)
            end if

            if( darkname /= '' ) then
               write(buf,'(f0.5)') t
               call fits_write_history(cfits,phid//" dark: '"//trim(darkname) &
                    //"' * "//trim(buf),status)
            end if

            if( biasname /= '' ) then
               call fits_write_history(cfits,phid//" bias: '"//trim(biasname)//"'",&
                    status)
            end if

            if( maskname /= '' ) then
               call fits_write_history(cfits,phid//" mask: '"//trim(maskname)//"'",&
                    status)
            end if

            call fits_update_key(cfits,FITS_KEY_CREATOR,FITS_VALUE_CREATOR, &
                 FITS_COM_CREATOR,status)
            call fits_write_comment(cfits,MUNIPACK_VERSION,status)

            ! corrected data
            call fits_write_image(cfits,0,OUT,status)

         else if( extname == EXT_STDERR ) then
            ! skip all previous std. error images
            continue

         else
            ! other HDUs are just copied
            call fits_copy_hdu(fits,cfits,0,status)
         end if

         if( status /= 0 ) goto 666
      end do
    end block

    ! stderr
    call fits_insert_img(cfits,bitpix,scifits%naxis,scifits%naxes,status)
    call fits_update_key(cfits,'EXTNAME',EXT_STDERR,'',status)
    call fits_write_comment(cfits,'The estimation of standard error.',status)
    call fits_write_image(cfits,0,DOUT,status)

    if( verbose ) write(error_unit,*) trim(output)

666 continue

    if( status == 0 ) then
       if( fits_file_exist(output) .and. overwrite ) &
            call fits_file_delete(output)
       call fits_file_duplicate(cfits,output,status)
    end if

    call fits_close_file(fits,status)
    call fits_delete_file(cfits,status)
    call fits_report_error(error_unit,status)

    corrector = status == 0

  end function corrector

  subroutine neighbour_mean(x,m)

    ! the mean of the eight neighbours of every pixel (fewer on the edges)

    real, dimension(:,:), intent(in) :: x
    real, dimension(:,:), intent(out) :: m
    real, dimension(:,:), allocatable :: s, c
    integer :: nx, ny, i, j

    nx = size(x,1)
    ny = size(x,2)
    allocate(s(0:nx+1,0:ny+1),c(0:nx+1,0:ny+1))
    s = 0
    c = 0
    s(1:nx,1:ny) = x
    c(1:nx,1:ny) = 1
    m = 0
    do j = -1, 1
       do i = -1, 1
          if( i == 0 .and. j == 0 ) cycle
          m = m + s(1+i:nx+i,1+j:ny+j)
       end do
    end do
    s = c
    c = 0
    do j = -1, 1
       do i = -1, 1
          if( i == 0 .and. j == 0 ) cycle
          c(1:nx,1:ny) = c(1:nx,1:ny) + s(1+i:nx+i,1+j:ny+j)
       end do
    end do
    m = m / max(c(1:nx,1:ny),1.0)

  end subroutine neighbour_mean

  function basename(path) result(b)

    ! the file name without the directory

    character(len=*), intent(in) :: path
    character(len=len(path)) :: b
    integer :: l

    l = index(path,'/',back=.true.)
    b = path(l+1:)

  end function basename

end program phcorr
