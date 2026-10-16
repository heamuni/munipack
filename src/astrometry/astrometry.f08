!
! astrometry
!
! Copyright © 2011-2026 F.Hroch (hroch@physics.muni.cz)
!
! This file is part of Munipack.
!
! Munipack is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Munipack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
!

! TODO:
!
!   * add all WCS related projections (done)
!   * add catalogue/reference matching identifier (done)
!   * add correlation matching by FFT with help of artificial
!
!  testing:
!   * improve high-precision algorithm to fit all parameters simulatenously (done)
!   * rewrite to use Penalty function to non-fitting parameters (centre) (cancelled)

! TODO
!   * add mirror projection (done)
!   * add FITS extension with astrometry results (done)
!   * filtering input for double stars (on photography or digital camera) (done)
!   * try to select better triangles firstly (abandoned)
!   * crossmatching by backtracking
!   * provide manual sequences as a reference parameter (done)
!   * estimate dispersion simultaneously
!   * use Rayleigh statistics (done)
!   * selection stars by magnitude limits (upper and lower, done)

program astrometry

  use titsio
  use astrofits
  use astrofitting
  use astromatcher
  use iso_fortran_env

  implicit none

  character(len=4*FLEN_FILENAME) :: record, key, val, seq1, seq2
  character(len=FLEN_FILENAME) :: ref = '', cat = '', rel = '', &
       output = '', filename = ''
  character(len=FLEN_FILENAME), dimension(:), allocatable :: filenames, outputs
  character(len=FLEN_VALUE) :: col_ra = 'RAJ2000', col_dec = 'DEJ2000', &
       col_pmra = '', col_pmdec = '', col_mag = '', col_epra = '', col_epdec = ''
  character(len=FLEN_VALUE) :: ptype = 'TAN', fit = 'ROBUST', aunits = '', &
       mode = 'MATCH', radesys = 'ICRS'
  character(len=FLEN_VALUE) :: key_dateobs = FITS_KEY_DATEOBS
  real(REAL64), dimension(2) :: crpix = 0.0_REAL64, crval = 0.0_REAL64
  real(REAL64) :: scale = 1.0_REAL64, angle = 0.0_REAL64, refl=1.0_REAL64
  ! the scale prior: the relative tolerance of the scale given by --scale,
  ! or by the astrometry calibration of a frame, for the match (0 = free)
  real(REAL64) :: scaletol = 0.05_REAL64
  ! copies of --scale, --angle, --reflex for the frames
  real(REAL64) :: user_scale, user_angle, user_refl
  real(REAL64) :: sig2 = 1.0, sig1 = 1.0/3600.0, xsig = 5
  real(REAL64) :: flxmin = 10.0**(-0.5/2.5), flxmax = 10.0**(0.5/2.5)
  real(REAL64) :: manepoch
  integer :: i,n,eq,status,istat,nseq=0,nfiles
  character(len=80) :: errmsg
  integer :: minmatch = 5, maxmatch = 42, maxstars = 250
  integer :: luckymatch = 3
  logical :: fcrpix = .false., fcrval = .false., reflex = .false. ,&
       fscale = .false., fangle = .false., &
       initpar = .true., rmscheck = .true., &
       fluxtest = .true., verbose = .false., wcssave = .true., &
       plog = .false., init_match = .false., report = .false., success
  logical, codimension[*] :: finished = .false.
  type(event_type), codimension[*] :: events
  integer, codimension[*] :: ifile = 0
  integer, codimension[*] :: nerr = 0
  integer, dimension(:), codimension[:], allocatable :: exitus

  ! a file to accumulate filenames
  open(22,status='scratch',form='unformatted',iostat=istat,iomsg=errmsg)
  if( istat /= 0 ) then
     write(error_unit,*) 'Error: ',trim(errmsg)
     stop 'Astrometry is unable to create a scratch file.'
  end if

  nfiles = 0

  if( this_image() == 1 ) then

     ! read input parameters
     do
        read(*,'(a)',iostat=istat,iomsg=errmsg) record
        if( is_iostat_end(istat) ) exit
        if( istat /= 0 ) then
           write(error_unit,*) 'Error: ',trim(errmsg)
           stop 'An std. input error in astrometry.'
        end if

        eq = index(record,'=')
        if( eq == 0 ) stop 'Malformed input record.'
        key = record(:eq-1)
        val = record(eq+1:)

        if( key == 'MODE' ) then

           read(val,*) mode

        elseif( key == 'PROJECTION' ) then

           read(val,*) ptype

        else if( key == 'FIT' ) then

           read(val,*) fit

        else if( key == 'AUNITS' ) then

           read(val,*) aunits

        else if( key == 'WCSSAVE' ) then

           read(val,*) wcssave

        else if( key == 'VERBOSE' ) then

           read(val,*) verbose

        else if( key == 'REPORT' ) then

           read(val,*) report

        else if( key == 'PIPELOG' ) then

           read(val,*) plog

        else if( key == 'CRPIX' ) then

           read(val,*) crpix
           fcrpix = .true.

        else if( key == 'CRVAL' ) then

           read(val,*) crval
           fcrval = .true.

        else if( key == 'SCALE' ) then

           read(val,*) scale
           fscale = .true.

        else if( key == 'ANGLE' ) then

           read(val,*) angle
           fangle = .true.

        else if( key == 'REFLEX' ) then

           read(val,*) reflex
           if( reflex ) refl = -1.0_REAL64

        else if( key == 'RADESYS' ) then

           read(val,*) radesys

        else if( key == 'EPOCH' ) then

           read(val,*) manepoch

        else if( key == 'SCALETOL' ) then

           read(val,*) scaletol

        else if( key == 'SIG' ) then

           read(val,*) sig2

        else if( key == 'SIGCAT' ) then

           read(val,*) sig1

        else if( key == 'FLXMIN' ) then

           read(val,*) flxmin

        else if( key == 'FLXMAX' ) then

           read(val,*) flxmax

        else if( key == 'FLUXCHECK' ) then

           read(val,*) fluxtest

        else if( key == 'INITPAR' ) then

           read(val,*) initpar

        else if( key == 'MINMATCH' ) then

           read(val,*) minmatch
           init_match = .true.

        else if( key == 'MAXMATCH' ) then

           read(val,*) maxmatch
           init_match = .true.

        else if( key == 'MAXSTARS' ) then

           read(val,*) maxstars

        else if( key == 'LUCKYMATCH' ) then

           read(val,*) luckymatch
           if( luckymatch < 0 ) stop 'Lucky match requires one or more stars.'

        else if( key == 'RMSCHECK' ) then

           read(val,*) rmscheck

        else if( key == 'COL_RA' ) then

           read(val,*) col_ra

        else if( key == 'COL_DEC' ) then

           read(val,*) col_dec

        else if( key == 'COL_PMRA' ) then

           read(val,*) col_pmra

        else if( key == 'COL_PMDEC' ) then

           read(val,*) col_pmdec

        else if( key == 'COL_MAG' ) then

           read(val,*) col_mag

        else if( key == 'COL_EPRA' ) then

           read(val,*) col_epra

        else if( key == 'COL_EPDEC' ) then

           read(val,*) col_epdec

        else if( key == 'FITS_KEY_DATEOBS' ) then

           read(val,*) key_dateobs

        else if( key == 'NSEQ' ) then

           read(val,*) nseq

        else if( key == 'SEQ1' ) then

           read(val,*) seq1

        else if( key == 'SEQ2' ) then

           read(val,*) seq2

        else if( key == 'REL' ) then

           read(val,*) rel

        else if( key == 'REF' ) then

           read(val,*) ref

        else if( key == 'CAT' ) then

           read(val,*) cat

        else if( key == 'FILE' ) then

           read(val,*) filename, output
           nfiles = nfiles + 1
           write(22) filename,output

        end if

     end do

     rewind(22)
     allocate(filenames(nfiles),outputs(nfiles))
     do n = 1, nfiles
        read(22) filenames(n),outputs(n)
     end do
     close(22)

  end if ! this_image() == 1

  call co_broadcast(verbose,1)
  call co_broadcast(mode,1)
  call co_broadcast(ptype,1)
  call co_broadcast(fit,1)
  call co_broadcast(aunits,1)
  call co_broadcast(wcssave,1)
  call co_broadcast(report,1)
  call co_broadcast(plog,1)
  call co_broadcast(crpix,1)
  call co_broadcast(fcrpix,1)
  call co_broadcast(crval,1)
  call co_broadcast(fcrval,1)
  call co_broadcast(scale,1)
  call co_broadcast(angle,1)
  call co_broadcast(fscale,1)
  call co_broadcast(fangle,1)
  call co_broadcast(reflex,1)
  if( reflex ) refl = -1.0_REAL64
  user_scale = scale
  user_angle = angle
  user_refl = refl
  call co_broadcast(scaletol,1)
  ! the solution needs five stars, the sequence as well (the front end
  ! rejects smaller values)
  if( minmatch < astrominstars ) stop 'Error: --minmatch must be five at least.'
  if( maxmatch <= minmatch ) stop 'Error: --maxmatch must be larger than --minmatch.'
  call co_broadcast(sig1,1)
  call co_broadcast(sig2,1)
  call co_broadcast(flxmin,1)
  call co_broadcast(flxmax,1)
  call co_broadcast(fluxtest,1)
  call co_broadcast(initpar,1)
  call co_broadcast(minmatch,1)
  call co_broadcast(maxmatch,1)
  call co_broadcast(maxstars,1)
  call co_broadcast(luckymatch,1)
  call co_broadcast(rmscheck,1)
  call co_broadcast(col_ra,1)
  call co_broadcast(col_dec,1)
  call co_broadcast(col_mag,1)
  call co_broadcast(col_pmra,1)
  call co_broadcast(col_pmdec,1)
  call co_broadcast(col_epra,1)
  call co_broadcast(col_epdec,1)
  call co_broadcast(key_dateobs,1)
  call co_broadcast(nseq,1)
  call co_broadcast(seq1,1)
  call co_broadcast(seq2,1)
  call co_broadcast(rel,1)
  call co_broadcast(ref,1)
  call co_broadcast(cat,1)
  call co_broadcast(nfiles,1)

  sync all

  if( this_image() > 1 ) then
     allocate(filenames(nfiles),outputs(nfiles))
  end if

  call co_broadcast(filenames,1)
  call co_broadcast(outputs,1)
  allocate(exitus(nfiles)[*])
  nerr = 0

  sync all

  if( mode == 'MATCH' .or. mode == 'SEQUENCE' .or. mode == 'CROSS' ) then

     block
       real(REAL64), allocatable, dimension(:) :: alpha,delta,pmalpha,pmdelta, &
            epalpha,epdelta,rflux
       real(REAL64), dimension(2) :: xcrval, xcrpix
       real(REAL64) :: repoch
       character(len=FLEN_VALUE) :: rframe = '', idcat = ''
       integer :: i,m,ic,n,iproc
       logical :: launch, fin, success

       ! lauching processes
       if( verbose .and. (this_image() == 1 .and. num_images() > 1) ) then
          write(error_unit,'(a,i0,a)') ' Launching the astrometry in ',num_images(), &
               ' parallel processes ...'
       end if

       if( cat /= '' ) then
          call read_cat(cat,alpha,delta,pmalpha,pmdelta,epalpha,epdelta,rflux,&
               rframe,idcat,repoch,xcrval,xcrpix)
       else if( ref /= '' ) then
          call read_ref(ref,alpha,delta,pmalpha,pmdelta,epalpha,epdelta,rflux,&
               rframe,idcat,repoch,xcrval,xcrpix)
       else if( rel /= '' ) then
          call read_rel(rel,alpha,delta,pmalpha,pmdelta,epalpha,epdelta,rflux,&
               rframe,idcat,repoch,xcrval,xcrpix)
       end if

       nerr = 0
       iproc = 0
       do n = 1, nfiles
          critical
            ic = 0
            do m = 1, num_images()
               ic = max(ic,ifile[m])
            end do
            launch = n > ic
            if( launch ) ifile = n
          end critical

          if( launch ) then
             iproc = iproc + 1
             call calibrate(filenames(n),outputs(n),alpha,delta,pmalpha,pmdelta,epalpha,epdelta,rflux,&
               rframe,idcat,repoch,xcrval,xcrpix,success)
             if( .not. success ) then
                nerr = nerr + 1
                exitus(nerr) = n
             end if
          end if
       end do
       finished = .true.

       if( verbose ) write(error_unit,'(2(a,i0))') 'Finished process #',this_image(),&
            ', total images processed: ', iproc

       ! the last finished process send an event to all images
       fin = .true.
       do i = 1, num_images()
          if( .not. finished[i] ) fin = .false.
       end do
       if ( fin ) then
          do i = 1, num_images()
             event post(events[i])
          end do
       end if

       event wait(events)

     end block

  elseif( mode == 'MANUAL' ) then

     if( this_image() == 1 ) then
        do n = 1, nfiles
           call mancal(filenames(n),outputs(n),success)
           if( .not. success ) then
              nerr = nerr + 1
              exitus(nerr) = n
           end if
        end do
     end if

     do i = 1, num_images()
        event post(events[i])
     end do
     event wait(events)

  elseif( mode == 'REMOVE' ) then

     if( this_image() == 1 ) then
        do n = 1, nfiles
           call remove(filenames(n),outputs(n),success)
           if( .not. success ) then
              nerr = nerr + 1
              exitus(nerr) = n
           end if
        end do
     end if

     do i = 1, num_images()
        event post(events[i])
     end do
     event wait(events)

  elseif( mode == 'DEFAULTS' ) then

     if( this_image() == 1 ) &
          call show_defaults

  end if

  sync all

  if( this_image() == 1 ) then

     block
       logical :: success = .true.
       integer :: k

       do i = 1, num_images()
          if( nerr[i] > 0 ) success = .false.
       end do

       if( success ) then
          !     stop 0
          write(error_unit,*) 'STOP 0'
       else

          if( verbose ) then
             write(error_unit,'(a)') 'The astrometry calibration of the file(s) was unsuccessful:'
             do i = 1, num_images()
                do k = 1, nerr[i]
                   write(error_unit,*) trim(filenames(exitus(k)[i]))
                end do
             end do
          end if
          stop 'Some error(s) occurred during the astrometry.'

       end if

     end block
  end if

contains

  subroutine calibrate(filename, output, alpha,delta,pmalpha,pmdelta,epalpha,epdelta,rflux,&
               rframe,idcat,repoch,xcrval,xcrpix,success)

    character(len=*), intent(in) :: filename, output
    real(REAL64), dimension(:), intent(in) :: alpha,delta,pmalpha,pmdelta,&
         epalpha,epdelta,rflux
    character(len=*), intent(in) :: rframe,idcat
    real(REAL64), intent(in) :: repoch
    real(REAL64), dimension(:), intent(in) :: xcrval,xcrpix
    logical, intent(out) :: success

    real(REAL64), parameter :: rad = 57.295779513082322865_REAL64

    type(fitsfiles) :: fits
    real(REAL64), allocatable, dimension(:) :: rec,dec,ares,dres,x,y,flux
    integer, allocatable, dimension(:) :: id1,id2
    integer, allocatable, dimension(:) :: i1,i2
    real(REAL64), dimension(2,0:11) :: pv
    character(len=2*FLEN_COMMENT), allocatable, dimension(:) :: com
    character(len=FLEN_VALUE) :: pcode, xframe
    real(REAL64), dimension(2) :: zcrval, dcrval, zcrpix
    real(REAL64) :: dscale = epsilon(1.0_REAL64), dangle=epsilon(1.0_REAL64), &
         a,da,b,db,s,ds, sm, epoch, rms=0.0_REAL64, zscale, zangle, zrefl
    ! the solution (the parameters given by the user, --scale, --angle,
    ! --reflex, are user_scale, user_angle, user_refl, unchanged
    ! for the next frames)
    real(REAL64) :: scale, angle, refl, scale0, ntrial
    integer, dimension(2) :: naxes
    integer :: n,nmin,nskip,nattempt
    logical :: matched = .false., reliable, zwcs
    real(REAL64), dimension(2) :: crval0
    ! the maximum of attempts to find a reliable match
    integer, parameter :: maxattempt = 20

    success = .false.
    scale = user_scale
    angle = user_angle
    refl = user_refl
    if( ptype == 'TPV' ) then
       pcode = 'TAN' ! for the base projection
    else
       pcode = ptype
    end if

    if( verbose ) &
       write(error_unit,*) "Astrometry calibration of `"//trim(filename)//"'."

    if( plog ) write(*,'(a)') "=ASTROMETRY> Start `"//trim(filename)//"'."

    status = 0
    call readfile(filename,key_dateobs,maxstars,x,y,flux,zcrpix,epoch,xframe,&
         zcrval,status,scale=zscale,rot=zangle,refl=zrefl,wcsok=zwcs,frame=naxes)
    if( status /= 0 ) then
       write(error_unit,*) "Error: Failed to open  `" //trim(filename)//"'"
       goto 77
    end if

    if( .not. fcrval ) crval = xcrval
    if( .not. fcrpix ) crpix = zcrpix

    ! a solution needs five stars at least, on the frame and in the catalogue
    if( size(x) < astrominstars .or. size(alpha) < astrominstars ) then
       write(error_unit,'(5a,i0,a,i0,a,i0,a)') "Error: astrometry of `",trim(filename), &
            "' with `",trim(cat)//trim(ref)//trim(rel),"': too few stars (",size(x), &
            " on the frame, ",size(alpha)," in the reference), at least ", &
            astrominstars," are required."
       goto 77
    end if

    ! the scale prior for the match: --scale, or the astrometry
    ! calibration of the frame, known within --scale-tol
    scale0 = -1
    if( scaletol > 0 .and. mode == 'MATCH' .and. rel == '' ) then
       if( fscale ) then
          scale0 = user_scale
       else if( zwcs ) then
          scale0 = zscale
       end if
       if( verbose .and. scale0 > 0 ) write(error_unit,'(a,f0.2,a,f0.1,a)') &
            ' The scale prior: ',scale0,' [pix/deg] +- ',100*scaletol,'%.'
    end if

    ! equatorial coordinates at current epoch
    if( IEEE_IS_NORMAL(epoch) ) then
       rec = alpha + (epoch - epalpha)*pmalpha !*  cos(delta / rad)
       dec = delta + (epoch - epdelta)*pmdelta
    else
       ! missing proper motions or epochs
       rec = alpha
       dec = delta
    end if

    ! When the solution by a match is rejected (it is unreliable, or not
    ! significant; a false match), the next match is tried.
    nskip = 0
    nattempt = 0
    ntrial = 0
    crval0 = crval
55  continue
    nattempt = nattempt + 1

    if( mode == 'MATCH' ) then

       if( plog ) then
          write(*,'(a)') "=ASTROMATCH> Start"
          call uvhisto(x,y)
       end if

       call astromatch(minmatch,maxmatch,pcode,rec,dec,rflux,x,y,flux, &
            crval(1),crval(2),crpix(1),crpix(2), &
            id1,id2,sig1,sig2,flxmin,flxmax,fluxtest,&
            luckymatch,verbose,matched,nskip,scale0,scaletol,ntrial)
!           write(*,'(a,f10.7)') "=MATCH> End "

    else if( mode == 'SEQUENCE' ) then

       block
         integer :: iostat1, iostat2
         allocate(id1(nseq),id2(nseq))
         read(seq1,*,iostat=iostat1) id1
         read(seq2,*,iostat=iostat2) id2
         if( .not. (nseq > 0 .and. iostat1 + iostat2 == 0) ) then
            write(error_unit,*) 'astrometry: sequences defined improperly.'
            return
         end if
         matched = nseq > 0
       end block

    else if( mode == 'CROSS' ) then

       block
         use nearmatch
         real(REAL64) :: xscale, xangle, xrefl, tol

         ! The approximate calibration is given by --scale, --angle, --reflex
         ! (and --rcen, --dcen, else the reference centre is used), or,
         ! when neither scale nor angle is given, it is taken from
         ! the current astrometry calibration of the frame.
         xscale = scale
         xangle = angle
         xrefl = refl
         if( .not. (fscale .or. fangle) ) then
            if( zwcs ) then
               xscale = zscale
               xangle = zangle
               xrefl = zrefl
               if( .not. fcrval ) crval = zcrval
               if( .not. fcrpix ) crpix = zcrpix
            else
               write(error_unit,*) "Error: cross-match of `",trim(filename), &
                    "' requires --scale and --angle, or an astrometry calibration."
               goto 77
            end if
         end if

         ! the cross-match tolerance [deg], xsig*sig2 [pix]
         tol = xsig*sig2 / xscale
         if( verbose ) write(error_unit,'(a,2f13.7,a,2f9.2,a,f0.3,a,f0.3,a,f3.0,a,f0.2,a)') &
              ' Cross-match by crval=',crval,' crpix=',crpix,' scale=',xscale, &
              ' [pix/deg] angle=',xangle,' refl=',xrefl,' tolerance=',tol*xscale,' [pix]'
         call listmatch(pcode,rec,dec,x,y,crval(1),crval(2),crpix(1),crpix(2), &
              xscale,xangle,xrefl,tol,id1,id2,matched)
       end block
    else
       write(error_unit,*) 'astrometry: unspecified match method.'
       return
    end if

    if( .not. matched ) then
       if(plog) write(*,'(a)') "=ASTROMATCH> Finish Fail"

       write(error_unit,'(7a)') "Warning: Match for files `", &
            trim(filename),"' and `",trim(cat),trim(ref),trim(rel), &
            "' failed. No common stars found."
    else
       if( plog ) write(*,'(a)') "=ASTROMATCH> Finish Success"

       if( plog ) write(*,'(a)') "=ASTROFIT> Start"
 !          if( verbose ) write(*,'(a)') "=> Fitting ..."
       call astrofit(fit,ptype,pcode,rec,dec,x,y,sig2,xsig, &
            crpix(1),crpix(2),crval(1),dcrval(1),crval(2),dcrval(2),&
            scale,dscale,refl,rms,angle,dangle,pv,ares,dres,id1,id2,i1,i2, &
            rmscheck,verbose,plog,reliable,scale0,scaletol,naxes,ntrial)

       if( .not. reliable .and. mode == 'MATCH' .and. nattempt < maxattempt ) then
          ! try a next match, nskip is the chunk pair of the rejected one
          crval = crval0
          if( verbose ) write(error_unit,'(a,i0,a)') &
               ' The solution is rejected, trying the next match (attempt #',nattempt+1,').'
          goto 55
       end if

       if( .not. reliable ) write(error_unit,*) &
            'Warning: The astrometry solution is unreliable for `',trim(filename),"'"

       if( reliable ) then

          if( plog ) then
             write(*,'(a,e15.7)') "=ASTROFIT> Final ",rms
             write(*,'(a)') "=ASTROFIT> Finish Success"
          end if

          allocate(com(17+size(i1)))

          a = cos(angle/rad)
          da = abs((dangle/rad)*sin(angle/rad))
          b = sin(angle/rad)
          db = abs((dangle/rad)*cos(angle/rad))
          s = 3600.0_REAL64*scale
          ds = s*dscale/scale

          if( aunits == '' ) then
             if( rel /= '' ) then
                aunits = 'pix'
             else if( log10(3600.0e6*rms) < 0.5 ) then
                aunits = 'uas'
             else if( log10(3600.0e3*rms) < 0.5 ) then
                aunits = 'mas'
             else if( log10(3600.0*rms) < 0.5 ) then
                aunits = 'arcsec'
             else if( log10(60.0*rms) < 0.5 ) then
                aunits = 'arcmin'
             else
                aunits = 'deg'
             end if
          end if

          if( rel /= '' ) then
             write(com(1),'(a)') "Type: relative"
             com(2) = "Reference frame: "//trim(rel)
          else
             write(com(1),'(a)') "Type: absolute"
             if( ref /= '' ) then
                com(2) = "Reference frame: "//trim(ref)
             else
                com(2) = "Reference catalogue: "//trim(idcat)//" (VizieR)"
             end if
          end if
          write(com(3),'(a)') "Projection: "//trim(ptype)
          write(com(4),'(a,i0)') "Match sequence length = ",size(id1)
          write(com(5),'(a,i0)') "Objects used = ",size(i1)
          write(com(6),'(a)') "Fitting method: "//trim(fit)
          if( aunits == 'arcsec' .or. aunits == 'arcmin' .or. &
               aunits == 'deg' .or. aunits == 'mas' .or. aunits == 'uas') then
             if( aunits == 'uas' ) then
                sm = 3600.0e6_REAL64
             else if( aunits == 'mas' ) then
                sm = 3600.0e3_REAL64
             else if( aunits == 'arcsec' ) then
                sm = 3600.0_REAL64
             else if( aunits == 'arcmin' ) then
                sm = 60.0_REAL64
             else if( aunits == 'deg' ) then
                sm = 1.0_REAL64
             end if
             s = sm/scale
             ds = s*dscale/scale
             write(com(7),'(a,en11.1,3a)') "RMS = ",sm*rms,' [',trim(aunits),']'
             write(com(8),'(a,g15.7,a,es7.1,3a)') "Scale = ",s," +- ",ds, &
                  ' [',trim(aunits),'/pix]'
          else if( aunits == "pix" ) then
             write(com(7),'(a,g15.7,a)') "RMS = ",rms,' [pix]'
             write(com(8),'(a,g15.7,a,es7.1)') "Scale = ",scale," +- ",dscale
          end if
          write(com(9),'(a,f12.7,a,es7.1)') "cos(angle) = ",a," +- ",da
          write(com(10),'(a,f12.7,a,es7.1)') "sin(angle) = ",b," +- ",db
          write(com(11),'(a,l5)') "Reflexion = ", refl < 0
          write(com(12),'(a,f11.5,a,es7.1,a)') "Angle of rotation(angle) = ", &
               angle," +- ",dangle," [deg]"
          if( aunits == 'pix' ) then
             write(com(13), '(a,f12.7,a,es7.1)') &
                 "Alpha center projection (CRVAL1) = ",crval(1)," +- ",dcrval(1)
             write(com(14),'(a,f12.7,a,es7.1)') &
                 "Delta center projection (CRVAL2) = ",crval(2)," +- ",dcrval(2)
          else
             write(com(13),'(a,f12.7,a,es7.1,a)') &
                "Alpha center projection (CRVAL1) = ",crval(1)," +- ",dcrval(1)," [deg]"
             write(com(14),'(a,f12.7,a,es7.1,a)') "Delta center projection (CRVAL2) = ",&
                  crval(2)," +- ",dcrval(2)," [deg]"
          end if
          write(com(15),'(a,f9.3,a)') "Horizontal center (CRPIX1) = ",crpix(1),' [pix]'
          write(com(16),'(a,f9.3,a)') "Vertical   center (CRPIX2) = ",crpix(2),' [pix]'

          if( aunits == 'uas' .or. aunits == 'mas' .or. aunits == 'arcsec' .or. &
               aunits == 'arcmin' .or. aunits == 'deg' ) then
             write(com(17),'(3a)') '    Equatorial RA,DEC [deg]     Data X,Y [pix]     Residuals [',trim(aunits),']'
          else if( aunits == 'pix' ) then
             write(com(17),'(a)') '   Reference X,Y [pix]       Data X,Y [pix]       Residuals [pix]'
          end if

          if( rel /= '' ) then
             do n = 1, size(ares)
                write(com(17+n),'(2f9.3,6x,2f10.3,1x,2f9.3)') &
                     alpha(i1(n))+xcrpix(1),delta(i1(n))+xcrpix(2),x(i2(n)),y(i2(n)),ares(n),dres(n)
             end do
          else
             do n = 1, size(ares)
                write(com(17+n),'(2f13.7,2f10.3,2en12.1)') &
                     alpha(i1(n)),delta(i1(n)),x(i2(n)),y(i2(n)),sm*ares(n),sm*dres(n)
             end do
          end if

!              if( plog ) write(*,'("=C> ",a)') (trim(com(n)),n=1,size(com))
          if( verbose ) then
             nmin = min(25,size(com))
             write(error_unit,'(a)') (trim(com(n)),n=1,nmin)
             if( nmin < size(com) ) write(error_unit,*) '...'
          end if

          if( wcssave ) then
             if( plog ) write(*,'(a)') "=WCSSAVE> Start"

             status = 0

             if( output == '' ) then
                call fits_open_file(fits,filename,FITS_READWRITE,status)
             else
                call fits_precopy_file(fits,filename,output,FITS_READWRITE, &
                     .true.,status)
             end if

             call wcsupdate(fits,ptype,crval,scale,angle,crpix,refl, &
                  pv,rms,rframe,epoch,com,status)

             call wcs_delete_report(fits,status)
             if( report ) then
                n = size(ares)
                call wcsreport(fits,epoch,com,alpha(i1(1:n)),delta(i1(1:n)), &
                     x(i2(1:n)),y(i2(1:n)),ares(1:n),dres(1:n),i1(1:n),id1,status)
             end if

             if( status == 0 ) then
                if( plog ) write(*,'(a)') "=WCSSAVE> Finish Success"
             else
                if( plog ) write(*,'(a)') "=WCSSAVE> Finish Fail"
                write(error_unit,*) &
                     "Warning: WCS save failed for `"//trim(filename)//"'."
             end if
             call fits_close_file(fits,status)
             call fits_report_error(error_unit,status)
             success = status == 0

          else
             success = .true.
          end if

          deallocate(com)

       else
          write(error_unit,'(a)') "Warning: Astrometry of `" // &
               trim(filename)//"' failed."
          if( plog ) write(*,'(a)') "=ASTROFIT> Finish Fail"
       end if

       if( allocated(ares)  ) deallocate(ares,dres,i1,i2)

    end if

77  continue
    if( allocated(x) ) deallocate(x,y,flux)
    if( allocated(id1) ) deallocate(id1,id2)

    if( plog ) write(*,'(a)') "=ASTROMETRY> Finish"

  end subroutine calibrate


  subroutine mancal(filename, output, success)

    use astrofits

    character(len=*), intent(in) :: filename, output
    logical, intent(out) :: success

    type(fitsfiles) :: fits
    character(len=FLEN_COMMENT), dimension(0) :: com
    real(REAL64), dimension(2,0:11) :: pv
    real(REAL64) :: rms = 0!, epoch = 2000.0 !!!!! need adjust
!    character(len=FLEN_VALUE) :: rframe = ''
    integer :: status

    if( verbose ) write(error_unit,*) &
         "Manual astrometry calibration of `"//trim(filename)//"'."

    pv = 0
    pv(1,1) = 1
    pv(2,1) = 1

    status = 0
    if( output == '' ) then
       call fits_open_file(fits,filename,FITS_READWRITE,status)
    else
       call fits_precopy_file(fits,filename,output,FITS_READWRITE,.true.,status)
    end if
    call wcsupdate(fits,ptype,crval,scale,angle,crpix,refl,pv,rms,&
         radesys,manepoch,com,status)
    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    success = status == 0

  end subroutine mancal


  subroutine remove(filename, output, success)

    use astrofits

    character(len=*), intent(in) :: filename, output
    logical, intent(out) :: success

    type(fitsfiles) :: fits
    integer :: status

    if( verbose ) write(error_unit,*) &
         "Removing astrometry calibration of `"//trim(filename)//"'."

    status = 0
    if( output == '' ) then
       call fits_open_file(fits,filename,FITS_READWRITE,status)
    else
       call fits_precopy_file(fits,filename,output,FITS_READWRITE,.true.,status)
    end if
    call wcsremove(fits,status)
    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)
    success = status == 0

  end subroutine remove


  subroutine read_rel(rel,alpha,delta,pmalpha,pmdelta,epalpha,epdelta,rflux, &
       rframe,idcat,repoch,xcrval,xcrpix)

    character(len=*), intent(in) :: rel
    real(REAL64), dimension(:), allocatable, intent(out) :: alpha,delta,pmalpha,pmdelta,&
         epalpha,epdelta,rflux
    character(len=*), intent(out) :: rframe,idcat
    real(REAL64), intent(out) :: repoch
    real(REAL64), dimension(:), intent(out) :: xcrval,xcrpix

    character(len=88) :: xframe
    real(REAL64), dimension(:), allocatable :: x,y
    integer :: n

    call readfile(rel,key_dateobs,maxstars,x,y,rflux,xcrpix,repoch,xframe, &
         xcrval,status)
    if( status /= 0 ) stop 'Error: a relative reference frame open failed.'

    ptype = ' '
    xcrval = 0 ! undefined
    scale = 1
    angle = 0
    n = size(x)
    allocate(alpha(n),delta(n),pmalpha(n),pmdelta(n),epalpha(n),epdelta(n))
    alpha = x - xcrpix(1)
    delta = y - xcrpix(2)
    pmalpha = 0.0_REAL64
    pmdelta = 0.0_REAL64
    epalpha = 0
    epdelta = 0
    idcat = ''
    rframe = ''

  end subroutine read_rel

  subroutine read_ref(ref,alpha,delta,pmalpha,pmdelta,epalpha,epdelta,rflux, &
       rframe,idcat,repoch,xcrval,xcrpix)

    character(len=*), intent(in) :: ref
    real(REAL64), dimension(:), allocatable, intent(out) :: alpha,delta,pmalpha,pmdelta,&
         epalpha,epdelta,rflux
    character(len=*), intent(out) :: rframe,idcat
    real(REAL64), intent(out) :: repoch
    real(REAL64), dimension(:), intent(out) :: xcrval,xcrpix
    integer :: n

    real(REAL64), dimension(:), allocatable :: x,y

    call readfile(ref,key_dateobs,maxstars,x,y,rflux,xcrpix,repoch,rframe, &
         xcrval,status,alpha,delta)
    if( status /= 0 .or. .not. allocated(alpha) ) &
       stop 'Error: a reference frame read failed or missing astrometry calibration.'

    n = size(alpha)
    allocate(pmalpha(n),pmdelta(n),epalpha(n),epdelta(n))
    pmalpha = 0.0_REAL64
    pmdelta = 0.0_REAL64
    epalpha = 0
    epdelta = 0
    idcat = ''

  end subroutine read_ref

  subroutine read_cat(cat,alpha,delta,pmalpha,pmdelta,epalpha,epdelta,rflux, &
       rframe,idcat,repoch,xcrval,xcrpix)

    use oakleaf

    character(len=*), intent(in) :: cat
    real(REAL64), dimension(:), allocatable, intent(out) :: alpha,delta,pmalpha,pmdelta,&
         epalpha,epdelta,rflux
    character(len=*), intent(out) :: rframe,idcat
    real(REAL64), intent(out) :: repoch
    real(REAL64), dimension(:), intent(out) :: xcrval,xcrpix

    call readtable(cat,col_ra,col_dec,col_pmra,col_pmdec,col_mag,col_epra, &
         col_epdec,maxstars,alpha,delta,pmalpha,pmdelta,epalpha,epdelta,rflux, &
         rframe,idcat,repoch,verbose,status)
    if( status /= 0 ) stop "Failed to open the catalogue file."
    if( size(alpha) < 1 ) stop "Error: no valid stars in the catalogue."

    ! a default centre of projection, --rcen, --dcen takes precedence
    xcrval(1) = median(alpha)
    xcrval(2) = median(delta)
    xcrpix = 0

  end subroutine read_cat


  subroutine show_defaults

    write(*,*) 'Default values of parameters :'
    write(*,*) '  Celestial projection code: ',trim(ptype)
    write(*,*) '  Right Ascension column label: ',trim(col_ra)
    write(*,*) '  Declination column label: ',trim(col_dec)
    write(*,*) '  Method for fitting: ',trim(fit)
    write(*,*) '  Maximum stars: ',maxstars
    write(*,*) '  Minimum stars for matching: ',minmatch
    write(*,*) '  Maximum stars for matching: ',maxmatch
    write(*,*) '  Number stars for lucky matching: ',luckymatch
    write(*,'(a,f0.3,a)') '   Relative tolerance of the scale prior: ',scaletol, &
         ' (--scale or the WCS of the frame, 0: free)'

    write(*,*) '  Check fluxes: ',fluxtest
    write(*,*) '  Check RMS:    ',rmscheck
    write(*,*) '  Save results: ',wcssave
    write(*,'(a,f0.1,a)') &
         '   Statistical error of centers stars on frame: ',sig2,' pix'
    write(*,'(a,f0.1,a)') &
         '   Statistical error of stars in catalogue: ',sig1*3600.0,' arcsec'
    write(*,'(a,4(f4.1,a))') '   Fluxes ratio range: ',flxmin,' .. ',flxmax,&
         '   in magnitudes:  ',2.5*log10(flxmin),' .. ',2.5*log10(flxmax)

  end subroutine show_defaults

end program astrometry
