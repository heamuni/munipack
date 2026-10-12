!
!  Dark     Average of a set of dark or bias frames.
!
!  Copyright © 1998-2024 F.Hroch (hroch@physics.muni.cz)
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

program dark

  use fitscorr
  use titsio
  use iso_fortran_env

  implicit none

  character(len=FLEN_FILENAME) :: darkname = 'dark.fits'   ! output
  character(len=FLEN_FILENAME) :: biasname = '' ! bias
  character(len=FLEN_FILENAME) :: maskname = '' ! mask

  character(len=80) :: msg
  character(len=4*FLEN_FILENAME) :: key, val, record
  character(len=FLEN_VALUE) :: dateobs, imagetyp
  character(len=FLEN_FILENAME), dimension(:), allocatable :: darknames
  character(len=FLEN_KEYWORD), dimension(6) :: keys = [ &
       character(len=FLEN_KEYWORD) :: &
       FITS_KEY_DATEOBS, &
       FITS_KEY_EXPTIME, &
       FITS_KEY_FILTER, &
       FITS_KEY_SATURATE, &
       FITS_KEY_TEMPERATURE, &
       FITS_KEY_GAIN ]

  integer :: ndark = 0           ! total count of darks
  integer :: bitpix = -32        ! default output dark bitpix
  integer :: status = 0          ! FITS status
  logical :: verbose = .false.
  integer :: eq, stat, naxis
  integer, dimension(2) :: naxes ! code is restricted on 2D frames

  ! check limits
  real :: tol_exptime = 1e-6
  real :: tol_temper = 1.0

  ! lower and upper limits
  real :: saturate = -1
  real :: threshold = 1
  logical :: threshold_set = .false.

  ! dark image and its statistical error
  real, dimension(:,:), allocatable :: image, imgerr
  real :: exptime_mean = 0,  exptime_err = 0  ! exposure time
  real :: temper_mean = -666, temper_err  = 0  ! temperature
  real :: dark_mean = 0,      dark_err = 0     ! mean level

  ! FITS wrappers
  type(CorrFits), allocatable :: biasfits, maskfits
  type(CorrFits), dimension(:), allocatable :: darks


  ! Section: Input  ---------------------------------------
  do
     read(*,'(a)',iostat=stat,iomsg=msg) record
     if( stat == IOSTAT_END ) exit
     if( stat > 0 ) then
        write(error_unit,*) trim(msg)
        error stop 'An input error.'
     end if

     eq = index(record,'=')
     if( eq == 0 ) error stop 'Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'OUTPUT' ) then
        read(val,*) darkname
     endif

     if( key == 'VERBOSE' ) then
        read(val,*) verbose
     endif

     if( key == 'BITPIX' ) then
        read(val,*) bitpix
     endif

     if( key == 'SATURATE' ) then
        read(val,*) saturate
        if( .not. (saturate > 0) ) stop 'Saturation > 0 is required.'
     endif

     if( key == 'THRESHOLD' ) then
        read(val,*) threshold
        if( .not. (threshold > 0) ) stop 'Threshold > 0 is required.'
        threshold_set = .true.
     endif

     if( key == 'FITS_KEY_DATEOBS' ) then
        read(val,*) keys(1)
     endif

     if( key == 'FITS_KEY_EXPTIME' ) then
        read(val,*) keys(2)
     endif

     if( key == 'FITS_KEY_SATURATE' ) then
        read(val,*) keys(4)
     endif

     if( key == 'FITS_KEY_TEMPERATURE' ) then
        read(val,*) keys(5)
     endif

     if( key == 'MUNIPACK_TEMPERATURE_TOLERANCE' ) then
        read(val,*) tol_temper
     endif

     if( key == 'MUNIPACK_EXPTIME_TOLERANCE' ) then
        read(val,*) tol_exptime
     endif

     if( key == 'BIAS' ) then
        read(val,*) biasname
     end if

     if( key == 'MASK' ) then
        read(val,*) maskname
     end if

     if( key == 'NFILES' ) then

        read(val,*) ndark
        allocate(darknames(ndark),stat=stat,errmsg=msg)
        if( stat /= 0 ) then
           write(error_unit,*) trim(msg)
           error stop 'Insufficient memory.'
        end if
        ndark = 0

     end if

     if( key == 'FILE' ) then

        ndark = ndark + 1
        if( ndark > size(darknames) ) stop 'NFILES unspecified?'
        read(val,*) darknames(ndark)

     end if

  enddo

  if( ndark == 0 ) stop 'Missing data'
  if( ndark /= size(darknames) ) error stop 'n /= size(darknames)'

  ! Section: FITS files read ------------------------------------------
  block

    integer :: n

    ! bias
    if( biasname /= '' ) then
       if( verbose ) write(error_unit,'(a)',advance="no") &
            "Bias frame: "//trim(biasname)//", "

       allocate(biasfits)
       call biasfits%Load(biasname,keys)
       if( .not. biasfits%status ) stop 'Failed to load bias frame.'

       if( verbose ) write(error_unit,'(a,1pg0.3,a,0pf0.1)') &
            ' exp.time[s] = ',biasfits%exptime, &
            ', T[degC] = ',biasfits%temper
    end if

    ! mask
    if( maskname /= '' ) then
       if( verbose ) write(error_unit,*) "Mask frame: ",trim(maskname)

       allocate(maskfits)
       call maskfits%Load(maskname)
       if( .not. maskfits%status ) stop 'Failed to load the mask frame.'
    end if

    if( verbose ) &
         write(error_unit,*) "Filename, exposure time[s], temperature [degC]:"

    ! darks
    allocate(darks(ndark),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Insufficient memory.'
    end if

    do n = 1, size(darknames)

       if( verbose ) write(error_unit,'(a)',advance="no") &
            trim(darknames(n))//":"

       call darks(n)%Load(darknames(n),keys)
       if( .not. darks(n)%status ) stop 'Frame read failed.'

       if( verbose ) write(error_unit,'(2x,1pg12.3,1x,0pf7.1)') &
            darks(n)%exptime,darks(n)%temper

       if( n == 1 .and. allocated(biasfits) ) then
          if( .not. all(biasfits%naxes == darks(n)%naxes) ) then
             stop "Dimensions of bias and the frame does not corresponds."
          end if
       end if

       if( n == 1 .and. allocated(maskfits) ) then
          if( .not. all(maskfits%naxes == darks(n)%naxes) ) then
             stop "Dimensions of bitmask and the frame does not corresponds."
          end if
       end if

       if( n > 1 ) then
          if( .not. all(darks(n-1)%naxes == darks(n)%naxes) ) then
             stop "Dimensions of images mutually does not corresponds."
          end if
       endif

       ! setup saturation
       if( saturate > 0 ) darks(n)%saturate = saturate
    end do

  end block

  ! setup common parameters
  naxis = darks(1)%naxis
  naxes = darks(1)%naxes
  dateobs = darks(1)%dateobs
  imagetyp = darks(1)%imagetyp

  if( verbose ) then
     write(error_unit,*)
     write(error_unit,*) 'Total count of frame(s):',ndark
     write(error_unit,*) 'Frame dimensions:',naxes(1),'x',naxes(2)
     write(error_unit,*) 'A bias or a dark frame is being prepared ...'
  end if

  ! Section: mean dark determination --------------------------------
  block

    use oakleaf

    real, dimension(:,:), allocatable :: bias
    logical, dimension(:,:), allocatable :: mask, bitmask
    real, dimension(:), allocatable :: x
    integer :: i,j,n,m
    real :: f

    allocate(bitmask(naxes(1),naxes(2)),mask(naxes(1),naxes(2)),x(ndark), &
         image(naxes(1),naxes(2)),imgerr(naxes(1),naxes(2)), &
         bias(naxes(1),naxes(2)),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Insufficient memory.'
    end if

    image = 0.0
    imgerr = 0.0

    if( allocated(biasfits) ) then
       bias = biasfits%image
       deallocate(biasfits)
    else
       bias = 0.0
    end if

    if( allocated(maskfits) ) then
       bitmask = maskfits%image > 0.5
       deallocate(maskfits)
    else
       bitmask = .true.
    end if

    if( ndark > 1 ) then

       do j = 1, naxes(2)
          do i = 1, naxes(1)

             if( bitmask(i,j) ) then

                n = 0
                do m = 1, ndark
                   f = darks(m)%image(i,j)
                   if( threshold < f .and. f < darks(m)%saturate ) then
                      n = n + 1
                      x(n) = f - bias(i,j)
                   end if
                end do

                if( n > 1 ) then
                   call rmean(x(1:n),image(i,j),imgerr(i,j))
                else
                   bitmask(i,j) = .false.
                end if
             end if
          enddo
       enddo

       ! mean exposure time
       n = 0
       do m = 1, ndark
          if( darks(m)%exptime_set ) then
             n = n + 1
             x(n) = darks(m)%exptime
          end if
       end do
       if( n > 0 ) then
          call rmean(x(1:n),exptime_mean,exptime_err)
       else
          exptime_mean = 0
          exptime_err = 0
       end if

       ! mean temperature
       n = 0
       do m = 1, ndark
          if( darks(m)%temper > -273 ) then
             n = n + 1
             x(n) = darks(m)%temper
          end if
       end do
       if( n > 0 ) then
          call rmean(x(1:n),temper_mean,temper_err)
       else
          temper_mean = -666
          temper_err = 0
       end if

       ! mean dark value
       mask = imgerr > 0
       call rmean(pack(image,mask),dark_mean,dark_err)

    else if( ndark == 1 ) then

       image = darks(1)%image
       temper_mean = darks(1)%temper
       temper_err = 0.0
       exptime_mean = darks(1)%exptime
       exptime_err = 0.0
       mask = threshold < image .and. image < darks(1)%saturate .and. bitmask
       call rmean(pack(image,mask),dark_mean,dark_err)

    end if

    ! replace undefined pixels by the mean
    where( .not. bitmask )
       image = dark_mean
       imgerr = dark_err
    end where

    deallocate(bias,mask,bitmask,x)
  end block

  ! checks
  block

    real :: t, e, maximg
    integer :: n

    ! check exptimes
    do n = 1, ndark
       e = darks(n)%exptime
       if( darks(n)%exptime_set .and. abs(e - exptime_mean) > tol_exptime ) &
            write(error_unit,*) "Warning: Exposure time out of limit for `",&
            trim(darks(n)%filename),"': ",e
    end do

    ! check temperatures
    do n = 1, ndark
       t = darks(n)%temper
       if( t > -273 .and. abs(t-temper_mean) > tol_temper ) &
            write(error_unit,*) "Warning: Temperature out of limit for `", &
            trim(darks(n)%filename),"': ",t
    enddo

    ! cut-off
    if( bitpix > 0 ) then
       maximg = 2.0**bitpix - 1
       image = max(0.0,min(image,maximg))
       imgerr = max(0.0,min(imgerr,maximg))
    end if

  end block

  if( verbose ) then
     write(error_unit,*) 'Average level: ',dark_mean,'+-',dark_err
     write(error_unit,*) 'Average exposure time: ',exptime_mean,'+-',exptime_err
     write(error_unit,*) 'Average temperature: ',temper_mean,'+-',temper_err
     write(error_unit,*) 'Saving to output file: ',trim(darkname)
  end if

  ! Section: FITS save -----------------------------------------------
  block

    integer, parameter :: group = 1
    character(len=FLEN_CARD) :: buf
    integer :: n
    type(fitsfiles) :: fits

    call fits_create_scratch(fits,status)
    call fits_insert_img(fits,bitpix,naxis,naxes,status)

    if( dateobs /= '' ) &
       call fits_write_key(fits,keys(1),dateobs,'UTC of the first on input',status)
    call fits_write_key(fits,keys(2),exptime_mean,6, &
         'The average of exposure times',status)
    call fits_write_key(fits,keys(5),temper_mean,4, &
         'The average of camera temperatures',status)
    if( imagetyp /= '' ) &
         call fits_write_key(fits,FITS_KEY_IMAGETYP,imagetyp,'The image type',status)
    write(buf,'(a,1pg0.6,a,1pg0.2)') &
         ' Average exposure time = ',exptime_mean,' +- ',exptime_err
    call fits_write_comment(fits,buf(2:),status)
    write(buf,'(a,f0.2,a,f0.2,a)') &
         ' Average temperature = ',temper_mean,' +- ',temper_err
    call fits_write_comment(fits,buf(2:),status)
    write(buf,'(a,1pg0.5,a,1pg0.1)') ' Average level = ',dark_mean,' +- ',&
         dark_err
    call fits_write_comment(fits,buf(2:),status)

    if( allocated(biasfits) ) &
         call fits_write_comment(fits,'BIAS: '//trim(biasfits%filename),status)

    if( allocated(maskfits) ) &
         call fits_write_comment(fits,'MASK: '//trim(maskfits%filename),status)

    if( ndark > 0 ) then
       write(buf,'(a,i0,a)') 'Result of average of ',ndark,'  exposure(s):'
       call fits_write_comment(fits,buf,status)
       do n = 1, ndark
          call fits_write_comment(fits,"'"//trim(darknames(n))//"'",status)
       end do
    endif

    if( biasname /= '' ) then
       call fits_write_history(fits,"DARK bias: '"//trim(biasname)//"'",status)
       ! the bias subtracted (phcorr: the error of the bias then enters
       ! the corrected frame by (1 - t), t the dark time factor; the
       ! STDERR of the dark has the scatter of the darks only)
       block
         character(len=FLEN_VALUE) :: val
         integer :: l
         l = len_trim(biasname)
         if( l > 64 ) then
            val = '...'//biasname(l-60:l)
         else
            val = biasname
         end if
         call fits_write_key(fits,'DARKBIAS',val, &
              'master bias subtracted from the dark',status)
       end block
    end if

    if( maskname /= '' ) &
         call fits_write_history(fits,"DARK bitmask: '"//trim(maskname)//"'", &
         status)

    if( threshold_set ) then
       write(buf,*) threshold
       call fits_write_history(fits,"DARK threshold: "//trim(buf),status)
    end if

    if( saturate > 0 ) then
       write(buf,*) saturate
       call fits_write_history(fits,"DARK saturation: "//trim(buf),status)
    end if

    call fits_update_key(fits,FITS_KEY_CREATOR,FITS_VALUE_CREATOR, &
         FITS_COM_CREATOR,status)
    call fits_write_comment(fits,MUNIPACK_VERSION,status)

    ! dark frame data
    call fits_write_image(fits,group,image,status)

    ! stderr
    call fits_insert_img(fits,bitpix,naxis,naxes,status)
    call fits_update_key(fits,'EXTNAME',EXT_STDERR,'',status)
    call fits_write_comment(fits, &
         'An estimation of standard error of mean of pixels of dark frame.',&
         status)
    call fits_write_image(fits,group,imgerr,status)

    if( status == 0 ) then
       if( fits_file_exist(darkname) ) call fits_file_delete(darkname)
       call fits_file_duplicate(fits,darkname,status)
    end if
    call fits_delete_file(fits,status)
    call fits_report_error(error_unit,status)

  end block

  deallocate(darks,darknames,image,imgerr)

  if( status == 0 ) then
     stop 0
  else
     stop 'An error occurred during some dark or bias frames determination.'
  end if

end program dark
