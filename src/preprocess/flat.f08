!
!  Flat     Average of a set of flat-fields.
!  Copyright (C) 1997 - 2024  Filip Hroch, Masaryk University, Brno, CZ
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

program flatmaker

  use fitscorr
  use oakleaf
  use titsio
  use iso_fortran_env

  implicit none

  ! debuging
  logical, parameter :: debug = .false.

  ! verbosity
  logical :: verbose = .false.

  ! Default Output image name:
  character(len=FLEN_FILENAME) :: flatname = 'flat.fits'
  character(len=FLEN_FILENAME) :: biasname = ''
  character(len=FLEN_FILENAME) :: darkname = ''
  character(len=FLEN_FILENAME) :: maskname = ''

  ! lower and upper limits
  real :: saturate =  -1
  real :: threshold = epsilon(threshold)
  logical :: saturate_set = .false.
  logical :: threshold_set = .false.

  integer :: eq, stat, naxis
  integer, dimension(2) :: naxes
  integer :: nflat = 0     ! counter of flats
  integer :: bitpix = -32
  integer :: maxiter = 7
  integer :: status = 0
  character(len=8) :: approximation = 'STANDARD'

  real :: flat_mean = 1, flat_err = 0
  real :: gain_mean = 1, gain_err = 0

  real :: time = 1
  real :: xdark = -1.0
  real :: gain = 1
  logical :: gain_set = .false.

  ! mean flat with statistical error
  real, dimension(:,:), allocatable :: flat, flaterr
  logical, dimension(:,:), allocatable :: bitmask

  character(len=80) :: msg
  character(len=4*FLEN_FILENAME) :: record, key, val
  character(len=FLEN_VALUE) :: dateobs, filter, imagetyp

  character(len=FLEN_FILENAME), dimension(:), allocatable :: flatnames
  character(len=FLEN_KEYWORD), dimension(6) :: keys = [ &
       character(len=FLEN_KEYWORD) :: &
       FITS_KEY_DATEOBS, &
       FITS_KEY_EXPTIME, &
       FITS_KEY_FILTER, &
       FITS_KEY_SATURATE, &
       FITS_KEY_TEMPERATURE, &
       FITS_KEY_GAIN ]

  ! FITS wrappers
  type(CorrFits), allocatable :: darkfits, biasfits, maskfits
  type(CorrFits), dimension(:), allocatable :: flats

  ! Section: Input  ---------------------------------------
  do
     read(*,'(a)',iostat=stat,iomsg=msg) record
     if( stat == IOSTAT_END ) exit
     if( stat > 0 ) then
        write(error_unit,*) trim(msg)
        error stop 'Input error.'
     end if

     eq = index(record,'=')
     if( eq == 0 ) error stop 'Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'VERBOSE' ) then
        read(val,*) verbose
     endif

     if( key == 'OUTPUT' ) then
        read(val,*) flatname
     endif

     if( key == 'BITPIX' ) then
        read(val,*) bitpix
     endif

     if( key == 'SATURATE' ) then
        read(val,*) saturate
        if( .not. (saturate > 0) ) stop 'Saturation > 0 is required.'
        saturate_set = .true.
     endif

     if( key == 'THRESHOLD' ) then
        read(val,*) threshold
        if( .not. (threshold > 0) ) stop 'Threshold > 0 is required.'
        threshold_set = .true.
     endif

     if( key == 'GAIN' ) then
        read(val,*) gain
        gain_set = .true.
     endif

     if( key == 'APPROXIMATION' ) then
        read(val,*) approximation
        if( approximation == 'BASIC' ) then
           maxiter = 0
        else !if( approximation == 'STANDARD' ) then
           maxiter = 3
        end if
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

     if( key == 'BIAS' ) then
        read(val,*) biasname
     end if

     if( key == 'XDARK' ) then
        read(val,*) xdark
     end if

     if( key == 'DARK' ) then
        read(val,*) darkname
     end if

     if( key == 'MASK' ) then
        read(val,*) maskname
     end if

     if( key == 'NFILES' ) then

        read(val,*) nflat
        allocate(flatnames(nflat),stat=stat,errmsg=msg)
        if( stat /= 0 ) then
           write(error_unit,*) trim(msg)
           error stop 'Insufficient memory.'
        end if
        nflat = 0

     end if

     if( key == 'FILE' ) then

        nflat = nflat + 1
        if( nflat > size(flatnames) ) stop 'NFILES unspecified?'
        read(val,*) flatnames(nflat)

     end if

  enddo

  if( nflat == 0 ) stop 'No input image(s).'
  if( nflat /= size(flatnames) ) error stop 'n /= size(flatnames)'

  ! Section: FITS files input ---------------------------------------
  block

    integer :: n
    character :: gflag
    logical :: filter_match = .true.
    logical :: filter_empty = .false.
    logical :: gain_default = .false.

    ! bias
    if( biasname /= '' ) then
       if( verbose ) write(error_unit,'(a)',advance="no") &
            "Bias frame: "//trim(biasname)//","

       allocate(biasfits)
       call biasfits%Load(biasname,keys)
       if( .not. biasfits%status ) stop 'Failed to load bias frame.'

       if( verbose ) write(error_unit,'(a,1pg0.3,a,0pf0.1)') &
            ' exp.time[s] = ',biasfits%exptime, &
            ', T[degC] = ',biasfits%temper
    end if

    ! dark
    if( darkname /= '' ) then
       if( verbose ) write(error_unit,'(a)',advance="no") &
            "Dark frame: "//trim(darkname)//", "

       allocate(darkfits)
       call darkfits%Load(darkname,keys)
       if( .not. darkfits%status ) stop 'Failed to load dark frame.'

       if( verbose ) write(error_unit,'(a,1pg0.3,a,0pf0.1)') &
            ' exp.time[s] = ',darkfits%exptime, &
            ', T[degC] = ',darkfits%temper
    end if

    ! mask
    if( maskname /= '' ) then
       if( verbose ) write(error_unit,*) "Mask frame: ",trim(maskname)

       allocate(maskfits)
       call maskfits%Load(maskname)
       if( .not. maskfits%status ) stop 'Failed to load the mask frame.'
    end if

    if( verbose ) write(error_unit,*) &
         "Filename, filter, exptime[s], gain[ct/adu], saturation[ct], T[degC]:"

    ! flats
    allocate(flats(nflat),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Insufficient memory.'
    end if

    do n = 1, size(flatnames)

       if( verbose ) write(error_unit,'(a)',advance="no") &
            trim(flatnames(n))//":"

       call flats(n)%Load(flatnames(n),keys)
       if( .not. flats(n)%status ) stop 'Flat-field frame read failed.'

       if( n == 1 .and. allocated(biasfits) ) then
          if( .not. all(biasfits%naxes == flats(n)%naxes) ) then
             stop "Dimensions of bias and the frame does not corresponds."
          end if
       end if

       if( n == 1 .and. allocated(darkfits) ) then
          if( .not. all(darkfits%naxes == flats(n)%naxes) ) then
             stop "Dimensions of dark and the frame does not corresponds."
          end if
       end if

       if( n == 1 .and. allocated(maskfits) ) then
          if( .not. all(maskfits%naxes == flats(n)%naxes) ) then
             stop "Dimensions of bitmask and the frame does not corresponds."
          end if
       end if

       if( n > 1 ) then
          if( .not. all(flats(n-1)%naxes == flats(n)%naxes) ) then
             stop "Dimensions of images mutually does not corresponds."
          end if
       endif

       ! setup saturation
       if( saturate_set ) flats(n)%saturate = saturate

       ! set gain by the provided value
       if( gain_set ) then
          flats(n)%gain = gain
          flats(n)%gain_set = .true.
       end if

       if( verbose ) then
          if( flats(n)%gain_set .or. gain_set ) then
             gflag = ' '
          else
             gflag = '!'
             gain_default = .true.
          end if

          write(error_unit,'(2x,a,2x,1pg0.3,2x,0pf0.3,a,2x,1p,g0.2,2x,0pf0.1)')&
               trim(flats(n)%filter),flats(n)%exptime,flats(n)%gain,gflag, &
               flats(n)%saturate,flats(n)%temper

          ! detect filter mismatch
          if( flats(n)%filter /= flats(1)%filter ) filter_match = .false.
          if( flats(n)%filter == '' ) filter_empty = .true.

       end if
    end do

    if( .not. filter_match ) write(error_unit,*) &
         "Warning: incompatible filters detected (try --verbose)."

    if( filter_empty ) write(error_unit,*) &
         "Warning: empty filter value is encountered (try --verbose). "

    if( verbose .and. gain_default ) &
       write(error_unit,*) "Gain default setup is indicated by '!' flag."

  end block

  ! setup common parameters
  naxis = flats(1)%naxis
  naxes = flats(1)%naxes
  dateobs = flats(1)%dateobs
  imagetyp = flats(1)%imagetyp
  filter = flats(1)%filter

  allocate(flat(naxes(1),naxes(2)),flaterr(naxes(1),naxes(2)), &
       bitmask(naxes(1),naxes(2)),stat=stat,errmsg=msg)
  if( stat /= 0 ) then
     write(error_unit,*) trim(msg)
     error stop 'Insufficient memory.'
  end if

  if( allocated(maskfits) ) then
     bitmask = maskfits%image > 0.5
     deallocate(maskfits)
  else
     bitmask = .true.
  end if

  if( verbose ) then
     write(error_unit,*)
     write(error_unit,*) 'Flat-field frames:',nflat
     write(error_unit,*) 'Flat frame dimensions:',naxes(1),'x',naxes(2)
     write(error_unit,*) 'Filter: ',trim(filter)
     write(error_unit,*) 'Accuracy of approximation: ',trim(approximation)
  end if

  ! Section: flat-field preparation --------------------------------
  ! All flats are corrected for gain, bias, darks.
  block

    real, dimension(:,:), allocatable :: bias, ebias, dark, edark
    logical, dimension(:,:), allocatable :: mask
    real :: exptime
    logical :: exptime_set
    integer :: n,m

    n = naxes(1)
    m = naxes(2)
    allocate(bias(n,m),ebias(n,m),dark(n,m),edark(n,m),mask(n,m), &
         stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Insufficient memory.'
    end if

    if( allocated(biasfits) ) then
       bias = biasfits%image
       ebias = biasfits%imgerr
       deallocate(biasfits)
    else
       bias = 0
       ebias = 0
    end if
    if( allocated(darkfits) ) then
       dark = darkfits%image
       edark = darkfits%imgerr
       exptime = darkfits%exptime
       exptime_set = darkfits%exptime_set .and. darkfits%exptime > 0
       deallocate(darkfits)
    else
       dark = 0
       edark = 0
       exptime_set = .false.
    end if

    if( verbose ) write(error_unit,*) &
         'Pre-corrections by gain, bias, or dark of raw flats ...'

    do n = 1, size(flats)

       ! dark frame multiplicator
       if( xdark > 0 )then
          time = xdark
       else if( flats(n)%exptime_set .and. exptime_set ) then
          time = flats(n)%exptime / exptime
       else
          time = 1
       end if

       ! Preparatory correct input flats (updates images in FitsCorr).
       ! Standard deviation is set with assumption of Poisson distribution
       ! of data. It requires large light fluxes, around half of a full range.
       associate (image => flats(n)%image, imgerr => flats(n)%imgerr, &
            saturate => flats(n)%saturate, gain => flats(n)%gain )

         mask = threshold < image .and. image < saturate .and. bitmask
         where( mask )
            image = gain*(image - (bias + time*dark))
            mask = image > epsilon(image)
         elsewhere
            ! Important. Pixels out of the mask has set negative values
            ! during all computations below, use of individual masks
            ! for every frame takes a lot of memory.
            image = -1
         end where
         where( mask )
            imgerr = sqrt(image + gain**2*(ebias**2 + time**2*edark**2))
         elsewhere
            imgerr = -1
         end where
         ! IMPORTANT
         ! imgerr means standard deviation: stdsig
         ! similarity stderr and imgerr, or meaning in dark.f08,
         ! is pure coincidental (and memory saving)

         ! determine of mean level of our current frame
         call rmean(pack(image,mask),flats(n)%mean,flats(n)%stderr,flats(n)%sig)

       end associate
    end do

    deallocate(mask,bias,ebias,dark,edark)
  end block


  ! Section: Initial flat-field estimate as mean of scaled flats  ------
  !          by averadge levels.
  block

    real, dimension(:), allocatable :: x
    integer :: i,j,n,m
    real :: thresh

    n = naxes(1)
    m = naxes(2)
    allocate(x(nflat))

    ! Initial mean flat
    if( verbose ) write(error_unit,'(a)') &
         'Calculating the initial flat-field frame (iter. #0) ...'
    do j = 1,naxes(2)
       do i = 1,naxes(1)

          ! this is initial estimate only
          ! lighter frames has higher influence due Poisson statistics
          n = 0
          do m = 1, nflat
             associate( image => flats(m)%image, imgerr => flats(m)%imgerr, &
                  mean => flats(m)%mean, &
                  saturate => flats(m)%saturate, gain => flats(m)%gain )

               thresh = gain*threshold
               if( thresh < image(i,j) .and. image(i,j) < saturate .and. &
                    bitmask(i,j) .and. mean > 0 ) then
                  n = n + 1
                  x(n) = image(i,j) / mean
               end if
             end associate
          end do

          if( n > 1 ) then
             call rmean(x(1:n),flat(i,j),flaterr(i,j))
          else
             flat(i,j) = 1
             flaterr(i,j) = 0
          end if

       enddo
    enddo

    deallocate(x)
  end block

  ! Section: debug prints, write out the first estimate -----------
  if( debug ) then

     block
       real, allocatable, dimension(:,:) :: res
       integer :: i,j

       allocate(res(naxes(1),naxes(2)),stat=stat,errmsg=msg)
       if( stat /= 0 ) then
          write(error_unit,*) trim(msg)
          error stop 'Insufficient memory.'
       end if

       where( bitmask .and. flaterr > 0 )
          res = (flat - 1) / (flaterr * sqrt(real(nflat)))
       end where

       open(1,file='/tmp/flatdebug_zero.dat')
       do i = 1,size(flat,1),2
          do j = 1,size(flat,2),2
             if( bitmask(i,j) .and. flaterr(i,j) > 0 .and. &
                  abs(res(i,j)) < 25 ) then
                write(1,*) flat(i,j)-1,res(i,j), flat(i,j)
             end if
          end do
       end do
       close(1)
       deallocate(res)
     end block
  end if

  ! Section: estimate the flat-field with scaling by individual -----
  !          ratios against previously estimated flat-field.
  block

    logical, dimension(:,:), allocatable :: mask
    real, dimension(:,:), allocatable  :: res
    real, dimension(:), allocatable :: u,ue,v,ve,gains
    integer :: i,j,m,n,iter
    logical :: terminate, reliable
    real :: t, dt, d, sig, avg, thresh

    allocate(res(naxes(1),naxes(2)),u(nflat),ue(nflat),v(nflat),ve(nflat), &
         gains(nflat),mask(naxes(1),naxes(2)),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Insufficient memory.'
    end if

    terminate = .false.
    do iter = 1, maxiter

       ! Now, we're improving accuracy of approximation. The number of
       ! iterations is controled by `terminate' variable which tests
       ! convergence of subsequent estimates of the created flat.

       if( terminate ) exit

       if( verbose ) then
          write(error_unit,'(a)') &
               'Scaling individual frames by the new flat ...'
          write(error_unit,'(a)') &
               "Filename, mean level[ct], std.err., std.dev., gain, reliable:"
       end if

       ! Update mean levels for individual frames
       do n = 1, nflat

          associate( image => flats(n)%image, imgerr => flats(n)%imgerr, &
               mean => flats(n)%mean, stderr => flats(n)%stderr, &
               stdsig => flats(n)%sig, &
               saturate => flats(n)%saturate, gain => flats(n)%gain )

            thresh = flats(n)%gain * threshold
            mask = thresh < image .and. image < saturate .and. imgerr > 0 &
                 .and. flat > 0 .and. flaterr > 0 .and. bitmask
            call fmean(pack(image,mask),pack(imgerr,mask),pack(flat,mask), &
                 pack(flaterr,mask),t,dt,sig,reliable=reliable)
!            call fmean(pack(flat,mask),pack(flaterr,mask), &
!                 pack(image,mask),pack(imgerr,mask),t,dt,sig,reliable=reliable)
            ! We assumes Poisson distribution of flat pixels, their dispersion
            ! is bound to the mean value. Data with strong fluence of
            ! non-Poisson noise component should be avoided already.

            ! update only when our estimate is realiable
            if( reliable ) then
               mean = t
               stderr = dt
               stdsig = sig
            end if

            if( verbose ) then

               gains(n) = t / sig**2
               write(error_unit, &
                    '(2a,2x,1pg0.5,2x,1pg0.3,1x,1pg0.5,2x,0pf0.3,2x,l1)') &
                    trim(flatnames(n)),": ",t,dt,sig,gains(n),reliable

            end if

          end associate

       enddo ! over all frames

       ! update flat
       if( verbose ) write(error_unit,'(a,i0,a)') &
            'Calculating accurate flat-field frame (iter. #',iter,') ...'
       res = -1
       do j = 1,naxes(2)
          do i = 1,naxes(1)

             n = 0
             do m = 1, nflat
                associate( image => flats(m)%image, imgerr => flats(m)%imgerr, &
                     mean => flats(m)%mean, sig => flats(m)%sig, &
                     saturate => flats(m)%saturate, gain => flats(m)%gain )

                  thresh = gain*threshold
                  if( thresh < image(i,j) .and. image(i,j) < saturate .and. &
                       imgerr(i,j) > 0 .and. bitmask(i,j) .and. &
                       mean > 0 .and. sig > 0 ) then
                     n = n + 1
                     u(n) = image(i,j)
                     ue(n) = imgerr(i,j)
                     v(n) = mean
                     ve(n) = sig
                  end if
                end associate
             enddo

             if( n > 0 ) then

                avg = flat(i,j)

                if( debug .and. i == naxes(1)/2 .and. j == naxes(2)/2 ) &
!                     call fmean(v(1:n),ve(1:n),u(1:n),ue(1:n), &
!                     flat(i,j),flaterr(i,j),reliable=reliable,verbose=.true.)
                     call fmean(u(1:n),ue(1:n),v(1:n),ve(1:n), &
                     flat(i,j),flaterr(i,j),reliable=reliable,verbose=.true.)

                call fmean(u(1:n),ue(1:n),v(1:n),ve(1:n), &
                     flat(i,j),flaterr(i,j))
!                call fmean(v(1:n),ve(1:n),u(1:n),ue(1:n), &
!                     flat(i,j),flaterr(i,j))

                ! Non reliable pixels are silently ignored. Some data
                ! has strongly non-gaussian distribution: the case of
                ! bad columns and pixels, an overscan data and etc.
                ! if( .not. reliable .and. verbose ) write(*,*) i,j

                ! the absolute difference between the result of previous
                ! computation and the current one controls termination
                if( reliable ) res(i,j) = abs(avg - flat(i,j))

             else
                flat(i,j) = 1
                flaterr(i,j) = 0
             end if

          enddo
       enddo

       ! condition for terminate: the mean difference of running flats
       ! between two latest subsequent iterations is smaller
       ! the mean std.err. limit
       d = qmedian(pack(res,res>0))
       flat_err = qmedian(pack(flaterr,flaterr>0))
       terminate = d < flat_err
       if( verbose ) &
            write(error_unit,'(a,3x,1p,2(g0.1,2x),l1)') &
            'Mean residual and std.dev., terminate:',d,flat_err,terminate

    end do ! iter

    ! final mean over the whole area
    call rmean(pack(flat,bitmask),flat_mean)

    if( verbose ) &
       call rmean(gains,gain_mean,gain_err)

    if( maxiter == 0 ) flat_err = qmedian(pack(flaterr,flaterr>0))

    deallocate(mask,res,u,ue,v,ve,gains)
  end block

  ! Section: diagnostics  ----------------------------------------
  block

    character(len=FLEN_FILENAME) :: buf
    logical, dimension(:,:), allocatable :: mask
    real, dimension(:,:), allocatable  :: res,des
    real :: thresh
    integer :: i,j,n

    allocate(mask(naxes(1),naxes(2)),res(naxes(1),naxes(2)), &
         des(naxes(1),naxes(2)),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Insufficient memory.'
    end if

    if( debug ) then
       ! (**)
       ! Diagnostics. The second column of the files
       ! are residuals intended for Normality testing.
       !
       ! https://stackoverflow.com/questions/2471884/
       ! or:
       ! gnuplot> binwidth=0.05
       ! gnuplot> bin(x,width)=width*floor(x/width)
       ! gnuplot> plot '/tmp/flatdebug_666.fits.dat' \
       !          using (bin($2,binwidth)):(1.0) smooth freq with boxes

       do n = 1, nflat

          associate( image => flats(n)%image, imgerr => flats(n)%imgerr, &
               mean => flats(n)%mean, saturate => flats(n)%saturate, &
               gain => flats(n)%gain  )

            thresh = gain*threshold
            mask = thresh < image .and. image < saturate .and. &
                 flat > 0 .and. flaterr > 0 .and. bitmask
            where( mask )
               des = sqrt(imgerr**2 + mean**2*flaterr**2)
               res = (image - mean*flat) / des
            end where

            write(buf,'(a,i0,a)') '/tmp/flatdebug_',n,'.dat'
            open(1,file=buf)
            write(1,'(2a)') '# ',trim(flatnames(n))
            do i = 1,size(image,1),2
               do j = 1,size(image,2),2
                  if( mask(i,j) .and. abs(res(i,j)) < 5 ) then
                     write(1,*) image(i,j)-mean, res(i,j), flat(i,j)
                  end if
               end do
            end do
            close(1)
          end associate

        end do
     end if ! debug
     deallocate(mask,res,des)
   end block

   ! Section: save in integers implies to scale the flat ------------------
   block

     integer :: waterline
     real :: maxflat

     ! Representation of flat by integers can't be recommended in any case,
     ! but it can be useful for a compatibility
     if( bitpix > 0 ) then

        ! scaled to waterline
        waterline = nint(10.0**(int(log10(2.0**(bitpix-1)))))
        ! waterline updates mean levels on 1e2,1e4 and 1e9
        flat_mean = waterline * flat_mean
        flat = waterline * flat
        flaterr = waterline * flaterr
        maxflat = 2.0**bitpix - 1

        ! range cut-off
        flat = max(0.0,min(flat,maxflat))
        flaterr = max(0.0,min(flaterr,maxflat))
        if( verbose ) write(error_unit,*) &
             'Warning: Numerical accurate degraded by conversion to integers.'
     end if
   end block

  if( verbose ) then
     write(error_unit,'(2a)') ' Output image: ',trim(flatname)
     write(error_unit,'(a,3x,1pg0.7)') ' Final mean:',flat_mean
     write(error_unit,'(a,3x,1pg0.1)') &
          ' Expected photometry standard error per pixel:',flat_err
     write(error_unit,'(a,2(2x,f0.3),a,f0.3,a)') &
          ' Estimated relative gain, std.err:',gain_mean,gain_err, &
          ' (original gain was ',gain,').'
  end if

  ! Section: FITS save ----------------------------------------------
  block

    integer, parameter :: group = 1
    character(len=*), parameter :: afid = 'FLAT'
    character(len=FLEN_CARD) :: buf
    integer :: n
    type(fitsfiles) :: fits

    call fits_create_scratch(fits,status)
    call fits_insert_img(fits,bitpix,naxis,naxes,status)
    if( dateobs /= '' ) &
         call fits_write_key(fits,keys(1),dateobs, &
         'UTC of the first on input',status)
    if( filter /= '' ) &
         call fits_write_key(fits,keys(3),filter,'filter',status)
    if( imagetyp /= '' ) &
         call fits_write_key(fits,FITS_KEY_IMAGETYP,imagetyp, &
         'Image type',status)
    if( verbose ) then
       call fits_write_key(fits,'GAIN_AVG',gain_mean,6, &
            '[ct/ADU] estimated gain',status)
       call fits_write_key(fits,'GAIN_STD',gain_err,2, &
            '[ct/ADU] std.dev of estimated gain',status)
    end if

    if( nflat > 0 ) then
       write(buf,'(a,i0,a)') 'Result of flat-fielding of ',nflat,' frames(s):'
       call fits_write_comment(fits,buf,status)
       do n = 1, nflat
          call fits_write_comment(fits,"'"//trim(flatnames(n))//"'",status)
       enddo
    endif

    if( gain_set ) then
       write(buf,'(f0.3)') gain
       call fits_write_history(fits,afid//" gain: "//trim(buf),status)
    end if

    if( darkname /= '' ) then
       write(buf,'(a,g0.5)') afid//" dark: '"//trim(darkname)//"' *",time
       call fits_write_history(fits,buf,status)
    end if

    if( biasname /= '' ) then
       buf = afid//" bias: '"//trim(biasname)//"'"
       call fits_write_history(fits,buf,status)
    end if

    if( maskname /= '' ) then
       buf = afid//" bitmask: '"//trim(maskname)//"'"
       call fits_write_history(fits,buf,status)
    end if

    if( threshold_set ) then
       write(buf,*) afid//" threshold: ", threshold, " (no gain applied)"
       call fits_write_history(fits,buf(2:),status)
    end if

    if( saturate_set ) then
       write(buf,*) afid//" saturation: ",saturate," (no gain applied)"
       call fits_write_history(fits,buf(2:),status)
    end if

    write(buf,*) afid//" mean level: ",flat_mean
    call fits_write_history(fits,buf(2:),status)

    call fits_update_key(fits,FITS_KEY_CREATOR,FITS_VALUE_CREATOR, &
         FITS_COM_CREATOR,status)
    call fits_write_comment(fits,MUNIPACK_VERSION,status)

    ! flat-field data
    call fits_write_image(fits,group,flat,status)

    ! standard error of mean
    call fits_insert_img(fits,bitpix,naxis,naxes,status)
    call fits_update_key(fits,'EXTNAME',EXT_STDERR,'',status)
    call fits_write_comment(fits,&
         'The estimation of standard error of mean of pixels of flat-field.',&
         status)

    call fits_write_image(fits,group,flaterr,status)

    if( status == 0 ) then
       if( fits_file_exist(flatname) ) call fits_file_delete(flatname)
       call fits_file_duplicate(fits,flatname,status)
    end if
    call fits_delete_file(fits,status)

    call fits_report_error(error_unit,status)

  end block

  deallocate(flats,flatnames,flat,flaterr,bitmask)

  if( status == 0 ) then
     stop 0
  else
     stop 'An error occurred during flat-field determination.'
  end if

end program flatmaker
