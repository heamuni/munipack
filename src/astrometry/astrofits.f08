!
! astrofits - write WCS to FITS image & read star table
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

module astrofits

  use titsio
  use astrotrafo
  use trajd
  use iso_fortran_env
  use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN, IEEE_IS_NORMAL

  implicit none

contains

  subroutine readtable(filename, label_alpha, label_delta, label_pmalpha, &
       label_pmdelta, label_mag, label_epoch_alpha, label_epoch_delta, &
       maxstar, alpha, delta, pmalpha, pmdelta, epalpha, epdelta, &
       refflux, rframe, catid, epoch, verbose,status)

    character(len=*), intent(in) :: filename, label_alpha, label_delta, &
         label_pmalpha, label_pmdelta, label_mag, &
         label_epoch_alpha, label_epoch_delta
    integer, intent(in) :: maxstar
    real(REAL64), allocatable, dimension(:), intent(out) :: alpha, delta, &
         pmalpha, pmdelta, refflux, epalpha, epdelta
    character(len=*), intent(out) :: rframe, catid
    real(REAL64), intent(out) :: epoch
    logical, intent(in) :: verbose
    integer, intent(out) :: status

    real(REAL64), parameter :: macheps = epsilon(1.0_REAL64)
    real(REAL64), parameter :: nullcoo = 1e33_REAL64
    real(REAL64), parameter :: rlim = 360 + macheps
    real(REAL64), parameter :: dlim = 90 + macheps
    real(REAL64), parameter :: pmlim = 3600e3 ! [degree in mas]
    real(REAL64), parameter :: nullmag = 99.999
    integer, parameter :: ncols = 7
    integer, dimension(ncols) :: col, statuses
    real(REAL64), allocatable, dimension(:) :: ra,dec,pmra,pmdec,mag,epra,epdec
    integer :: nrows,i,l,n,frow,srows,stat
    real(REAL64) :: nullep
    logical :: anyf
    character(len=FLEN_VALUE), dimension(ncols) :: label,colname
    character(len=80) :: errmsg
    type(fitsfiles) :: fits

    label(1) = label_alpha
    label(2) = label_delta
    label(3) = label_pmalpha
    label(4) = label_pmdelta
    label(5) = label_mag
    label(6) = label_epoch_alpha
    label(7) = label_epoch_delta

    status = 0

    ! open, and move, to a table extension
    call fits_open_table(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read a table in the file `',trim(filename),"'."
       return
    end if

    call fits_get_num_rows(fits,nrows,status)
    if( status /= 0 ) goto 666
    if( .not. (nrows > 0) )then
       write(error_unit,*) 'Error: empty table in the file `',trim(filename),"'."
       goto 666
    end if

    ! define reference frame and identification of the catalogue
    rframe = ''
    catid = ''
    call fits_read_key(fits,'EXTNAME',catid,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       catid = 'unknown cat' ! somethig is terribly bad
       status = 0
    end if

    call fits_read_key(fits,FITS_KEY_RADESYS,rframe,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       rframe = ''
       status = 0
    end if

    call fits_read_key(fits,FITS_KEY_EPOCH,epoch,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       epoch = IEEE_VALUE(epoch, IEEE_QUIET_NAN)
       if( verbose ) &
            write(error_unit,*) "Warning: the keyword `",trim(FITS_KEY_EPOCH), &
            "' not found in `",trim(filename),"'. Proper motions are ignored."
    end if

    ! find columns by labels
    statuses = 0
    do i = 1, ncols
       colname(i) = ''
       if( label(i) /= '' ) then
          call fits_get_colname(fits,.false.,label(i),colname(i),col(i),statuses(i))
       else
          statuses(i) = -1
       end if
    end do
    if( any(statuses(1:2) /= 0) ) then
       write(error_unit,*) &
            'Error: Right Ascension and Declination columns not found.'
       status = maxval(statuses(1:2))
       goto 666
    end if

    allocate(ra(nrows),dec(nrows),pmra(nrows),pmdec(nrows),mag(nrows),&
         epra(nrows),epdec(nrows),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,'(2a)') 'Error:', trim(errmsg)
       goto 666
    end if

    epra = epoch
    epdec = epoch
    pmra = 0.0_REAL64
    pmdec = 0.0_REAL64

    call fits_get_rowsize(fits,srows,status)

    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       ! coordinates
       call fits_read_col(fits,col(1),frow,nullcoo,ra(i:l),anyf,status)
       call fits_read_col(fits,col(2),frow,nullcoo,dec(i:l),anyf,status)

       ! proper motion
       if( all(statuses(3:4) == 0) ) then
          call fits_read_col(fits,col(3),frow,nullcoo,pmra(i:l),anyf,status)
          call fits_read_col(fits,col(4),frow,nullcoo,pmdec(i:l),anyf,status)
       end if

       ! magnitudes
       if( statuses(5) == 0 ) then
          call fits_read_col(fits,col(5),frow,nullmag,mag(i:l),anyf,status)
       end if

       ! epochs for individual stars
       if( all(statuses(6:7) == 0) ) then
          nullep = IEEE_VALUE(nullep,IEEE_QUIET_NAN)
          call fits_read_col(fits,col(6),frow,nullep,epra(i:l),anyf,status)
          call fits_read_col(fits,col(7),frow,nullep,epdec(i:l),anyf,status)
       end if

       if( status /= 0 ) goto 666
    end do

    ! select valid stars
    n = 0
    do i = 1, nrows
       if( abs(ra(i)) < rlim .and. abs(dec(i)) < rlim .and. mag(i) < 99 &
            .and. n < maxstar ) then
          n = n + 1
          ra(n) = ra(i)
          dec(n) = dec(i)
          mag(n) = mag(i)
          if( abs(pmra(i)) < pmlim .and. abs(pmdec(i)) < pmlim ) then
             pmra(n) = pmra(i) !/ 3.6e6_REAL64  ! [mas/year] -> [deg/year]
             pmdec(n) = pmdec(i) !/ 3.6e6_REAL64
          else
             pmra(n) = 0
             pmdec(n) = 0
          end if
          if( IEEE_IS_NORMAL(epra(i)) .and. IEEE_IS_NORMAL(epdec(i)) ) then
             epra(n) = epra(i)
             epdec(n) = epdec(i)
          else if( IEEE_IS_NORMAL(epoch) ) then
             epra(n) = epoch
             epdec(n) = epoch
          else
             ! no time info, the proper motion is ignored
             epra(n) = 0
             epdec(n) = 0
             pmra(n) = 0
             pmdec(n) = 0
          end if
!          write(*,*) ra(n),dec(n),pmra(n),pmdec(n),mag(n),epra(n),epdec(n)
       end if
    end do

    allocate(alpha(n),delta(n),pmalpha(n),pmdelta(n),refflux(n), &
         epalpha(n), epdelta(n), stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,'(2a)') 'Error:', trim(errmsg)
       goto 666
    end if

    alpha = ra(1:n)
    delta = dec(1:n)
    pmalpha = pmra(1:n) / 3.6e6_REAL64  ! [mas/year] -> [deg/year]
    pmdelta = pmdec(1:n) / 3.6e6_REAL64
    epalpha = epra(1:n)
    epdelta = epdec(1:n)
    refflux = 10.0**(0.4*(25 - mag(1:n)))

    call fits_close_file(fits,status)
    return

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine readtable


  subroutine readfile(filename, key_dateobs, maxstar, &
       xcoo, ycoo, flux, crpix, epoch, radesys, crval,status,alpha,delta, &
       scale,rot,refl,wcsok)

    ! On output, crval are spherical coordinates of the reference point crpix
    ! (the centre of the frame) by its astrometry calibration (if any);
    ! scale, rot and refl describe the calibration, wcsok indicates its presence.

    character(len=*), intent(in) :: filename, key_dateobs
    integer, intent(in) :: maxstar
    real(REAL64), allocatable, dimension(:), intent(out) :: xcoo, ycoo, flux
    real(REAL64), dimension(:), intent(out) :: crpix, crval
    character(len=*), intent(out) :: radesys
    real(REAL64), intent(out) :: epoch
    integer,intent(out) :: status
    real(REAL64), allocatable, dimension(:), intent(out), optional :: alpha, delta
    real(REAL64), intent(out), optional :: scale, rot, refl
    logical, intent(out), optional :: wcsok

    integer :: naxis, bitpix, nrows, n, ns, i, l, srows, frow, wstatus
    integer, dimension(2) :: naxes
    logical :: anyf
    integer, parameter :: extver = 0
    real(REAL64), allocatable, dimension(:) :: x,y,cts
    real(REAL64), parameter :: nullval = 0.0_REAL64
    character(len=FLEN_CARD) :: colname1, colname2, dateobs, aplabel
    integer :: col1, col2, nkeys
    character(len=FLEN_CARD), dimension(0) :: exclist
    character(len=:), allocatable :: header
    character(len=FLEN_VALUE) :: code
    type(AstroWCS) :: wcs
    type(fitsfiles) :: fits

    status = 0
    wstatus = -1
    if( present(wcsok) ) wcsok = .false.
    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read an image in the file `',trim(filename),"'."
       return
    end if

    call fits_get_img_param(fits,bitpix,naxis,naxes,status)
    if( status /= 0 .or. naxis /= 2 ) then
       if( naxis /= 2 ) then
            write(error_unit,*) trim(filename), &
                 ": Sorry, just two dimensional images are implemented yet."
            status = -1
         end if
       goto 666
    end if
    crpix = naxes / 2.0_REAL64

    call fits_read_key(fits,key_dateobs,dateobs,status)
!    call fits_get_dateobs(fits,key_dateobs,dateobs,status)
    if( status == 0 ) then
       epoch = jdepoch(fits_jd(dateobs,status))
    else if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       epoch = IEEE_VALUE(epoch, IEEE_QUIET_NAN)
       write(error_unit,*) "Warning: the keyword `",trim(key_dateobs), &
            "' not found in `",trim(filename),"'. Proper motions are ignored."
    end if
    if( status /= 0 ) goto 666

    ! get WCS calibration
    call fits_hdr2str(fits,1,exclist,nkeys,header,status)
    if( status == 0 ) then
       call wcs % init(header,nkeys,wstatus)
       if( wstatus == 0 ) then
          ! a celestial calibration is required
          call wcs % get_code(code)
          if( code == '' ) wstatus = -1
       end if
       if( wstatus == 0 ) then
          call wcs % pix2cel(crpix(1),crpix(2),crval(1),crval(2))
          call wcs % get_radesys(radesys)
          if( present(scale) .and. present(rot) .and. present(refl) ) &
               call wcs % get_affines(scale,rot,refl)
          if( present(wcsok) ) wcsok = .true.
       else
          crval = -999.999
       end if
    else
       status = 0
       crval = -999.999
    end if
    deallocate(header)

    ! FIND extension
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,FINDEXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       write(error_unit,*) trim(FINDEXTNAME)// &
            ': extension not found. Coordinates of objects are not available.'
       goto 666
    end if
    call fits_get_num_rows(fits,nrows,status)
    if( status /= 0 ) goto 666

    allocate(x(nrows),y(nrows))

    ! the column labels are specific for MUNIPACK extension only
    colname1 = ''
    colname2 = ''
    call fits_get_colname(fits,.false.,FITS_COL_X,colname1,col1,status)
    call fits_get_colname(fits,.false.,FITS_COL_Y,colname2,col2,status)

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,col1,frow,nullval,x(i:l),anyf,status)
       call fits_read_col(fits,col2,frow,nullval,y(i:l),anyf,status)

       if( status /= 0 ) goto 666
    end do

    ! photometry extension
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,APEREXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       status = 0
       ! aperture photometry is not available, no problem,
       ! we are continuing still and all ...?
       ns = size(x)
       allocate(xcoo(ns),ycoo(ns),flux(ns))
       xcoo = x
       ycoo = y
       flux = 1
    else if( status == 0 ) then
       ! aperture photometry is available, we can use photometry informations
       ! to improve reliability
       aplabel = trim(FITS_COL_APCOUNT)//"1"
       call fits_get_colname(fits,.false.,aplabel,colname1,col1,status)
       allocate(cts(nrows))
       frow = 1
       call fits_read_col(fits,col1,frow,nullval,cts,anyf,status)

       ! we're selecting stars having valid photometry only
       ns = min(count(cts > 0.0),maxstar)
       allocate(xcoo(ns),ycoo(ns),flux(ns))
       ns = 0
       do n = 1,nrows
          if( cts(n) > 0.0 .and. ns < size(xcoo) .and. ns <= maxstar ) then
             ns = ns + 1
             xcoo(ns) = x(n)
             ycoo(ns) = y(n)
             flux(ns) = cts(n)
          end if
       end do
       deallocate(cts)

    end if
    if( status /= 0 ) goto 666

    if( present(alpha) .and. present(delta) .and. wstatus == 0 ) then
       n = size(xcoo)
       allocate(alpha(n),delta(n))
       call wcs % pix2cel(xcoo,ycoo,alpha,delta)
    end if

    call fits_close_file(fits,status)

    deallocate(x,y)
    return

666 continue

    if( allocated(x) ) deallocate(x,y)
    if( allocated(cts) ) deallocate(cts)

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine readfile

!-----------------------------------------------------------------------------

  subroutine wcsupdate(fits,type,crval,sc,pa,crpix,refl,pv,rms,radesys,&
       epoch,com,status)

    use astrotrafo

    integer, parameter :: ndim = 2
    real(REAL64), parameter :: rad = 57.295779513082322865_REAL64

    type(fitsfiles) :: fits
    character(len=*), intent(in) :: type, radesys
    real(REAL64), dimension(:), intent(in) :: crpix,crval
    real(REAL64), intent(in) :: sc,pa,rms,refl,epoch
    real(REAL64), dimension(:,0:), intent(in) :: pv
!    logical, intent(in) :: fcrpix
    character(len=*), dimension(:), intent(in) :: com
    integer, intent(in out) :: status

    character(len=3), dimension(2,0:11) :: pvc
!    integer :: naxis, bitpix
!    integer, dimension(ndim) :: naxes
!    real(REAL64), dimension(ndim) :: crpix
    character(len=6), dimension(4), parameter :: keys = &
         [ 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2' ]
    real(REAL64), dimension(2,2) :: cd
    type(AstrometryGeometry) :: tg
    character(len=5) :: un
    character(len=FLEN_KEYWORD) :: pvkey
    integer :: i,j

    pvc(1,:) = ['1  ', 'x  ', 'y  ', 'r  ', 'xx ', 'xy ', 'yy ', &
         'xxx', 'xxy', 'xyy', 'yyy', 'rrr' ]
    pvc(2,:) = ['1  ', 'y  ', 'x  ', 'r  ', 'yy ', 'yx ', 'xx ', &
         'yyy', 'yyx', 'xxy', 'xxx', 'rrr' ]

    if( type == " " ) then
       un = '[pix]'
    else
       un = '[deg]'
    end if

    ! clear a previous calibration
    call wcsremove(fits,status)

!    if( fcrpix ) then
!       crpix = crpixels
!    else
!       call fits_get_img_param(fits,bitpix,naxis,naxes,status)
!       if( status /= 0 ) return
!       crpix = naxes / 2
!    end if

    call fits_update_key(fits,FITS_KEY_CREATOR,FITS_VALUE_CREATOR,FITS_COM_CREATOR,status)

    if( type /= " " )then
       call fits_update_key(fits,'CTYPE1', &
            'RA---'//type(1:3),'the coordinate type for the first axis',status)
       call fits_update_key(fits,'CTYPE2', &
            'DEC--'//type(1:3),'the coordinate type for the second axis',status)
       call fits_update_key(fits,'CRVAL1',crval(1),16, &
            un//' first axis value at reference pixel',status)
       call fits_update_key(fits,'CRVAL2',crval(2),16, &
            un//' second axis value at reference pixel',status)
    else !if( type == " " ) then
       call fits_write_errmark
       do i = 1,size(keys)
          call fits_delete_key(fits,keys(i),status)
          if( status == FITS_KEYWORD_NOT_FOUND ) status = 0
       end do
       call fits_clear_errmark
    end if

    call fits_update_key(fits,'CRDER1',rms,1, &
         un//' random statistical error in first axis',status)
    call fits_update_key(fits,'CRDER2',rms,1, &
         un//' random statistical error in second axis',status)

    if( type == " " ) then
       call fits_update_key(fits,'CUNIT1','pix','units of first axis',status)
       call fits_update_key(fits,'CUNIT2','pix','units of second axis',status)
    else
       call fits_update_key(fits,'CUNIT1','deg','units of first axis',status)
       call fits_update_key(fits,'CUNIT2','deg','units of second axis',status)
    end if

    ! flip x-axis due to conversion rectangular -> spherical coordinates
    call tg % init(sc,pa,refl)
    call tg % get_fcd(cd)

    if( type /= '' ) then
       call fits_update_key(fits,'CRPIX1',crpix(1),16, &
            'x-coordinate of reference pixel',status)
       call fits_update_key(fits,'CRPIX2',crpix(2),16, &
            'y-coordinate of reference pixel',status)
    else
       call fits_update_key(fits,'CRPIX1',crpix(1)+crval(1),16, &
            'x-coordinate of reference pixel',status)
       call fits_update_key(fits,'CRPIX2',crpix(2)+crval(2),16, &
            'y-coordinate of reference pixel',status)
    end if

    call fits_update_key(fits,'CD1_1',cd(1,1),16, &
         'partial of first axis coordinate w.r.t. x',status)
    call fits_update_key(fits,'CD1_2',cd(1,2),16, &
         'partial of first axis coordinate w.r.t. y',status)
    call fits_update_key(fits,'CD2_1',cd(2,1),16, &
         'partial of second axis coordinate w.r.t. x',status)
    call fits_update_key(fits,'CD2_2',cd(2,2),16, &
         'partial of second axis coordinate w.r.t. y',status)

    if( type == 'TPV' ) then
       do j = 1, 2
          do i = lbound(pv,2), ubound(pv,2)
             write(pvkey,'(a,i1,a,i0)') 'PV',j,'_',i
             call fits_update_key(fits,pvkey,pv(j,i),15, &
                  'TPV coefficient: '//pvc(j,i),status)
          end do
       end do
    end if

    call fits_update_key(fits,'RADESYS',radesys,'reference coordinate system',status)

    if( IEEE_IS_NORMAL(epoch) ) call fits_update_key(fits,FITS_KEY_EPOCH,epoch,-6, &
         'epoch of equatorial coordinates',status)

    call fits_write_comment(fits,BEGIN_ASTROMETRY,status)
    do i = 1, size(com)
       if( index(com(i),'Residuals') > 0 ) exit
       call fits_write_comment(fits,com(i),status)
    end do

    ! no fitting was done
    if( size(com) < 1 ) &
         call fits_write_comment(fits,"Astrometry calibration provided by user.",status)

    call fits_write_comment(fits,MUNIPACK_VERSION,status)
    call fits_write_comment(fits,END_ASTROMETRY,status)

  end subroutine wcsupdate


!-----------------------------------------------------------------------------


  subroutine wcsremove(fits,status)

    type(fitsfiles) :: fits
    integer, intent(in out) :: status
    integer :: i, j, nbegin, nend
    character(len=FLEN_CARD) :: card
    character(len=7), dimension(16) :: keys = [ &
         'CTYPE1 ', 'CTYPE2 ','CRVAL1 ','CRVAL2 ','CRDER1 ','CRDER2 ', &
         'CRPIX1 ', 'CRPIX2 ','CUNIT1 ','CUNIT2 ','CD1_1  ','CD1_2  ', &
         'CD2_1  ', 'CD2_2  ','RADESYS','JEPOCH ' ]
    character(len=FLEN_KEYWORD) :: pvkey

    call fits_write_errmark
    do i = 1, size(keys)
       call fits_delete_key(fits,keys(i),status)
       if( status == FITS_KEYWORD_NOT_FOUND ) status = 0
    end do
    do j = 1, 2
       do i = 0, 39
          write(pvkey,'(a,i1,a,i0)') 'PV',j,'_',i
          call fits_delete_key(fits,pvkey,status)
          if( status == FITS_KEYWORD_NOT_FOUND ) status = 0
       end do
    end do
    call fits_clear_errmark

    i = 0
    nbegin = 0
    nend = 0
    card = ''
    do
       i = i + 1
       call fits_read_record(fits,i,card,status)
       if( status /= 0) exit
       if( card(9:) == BEGIN_ASTROMETRY ) nbegin = i
       if( card(9:) == END_ASTROMETRY ) nend = i
    end do
    if( status == 203 ) status = 0

    if( nbegin > 0 ) then
       do i = nbegin, nend
          ! looks strangle, but it's correct --
          ! records in header are moved up imediately following delete
          call fits_delete_record(fits,nbegin,status)
       end do
    end if

  end subroutine wcsremove


  subroutine wcsreport(fits,epoch,com,alpha,delta,x,y,ares,dres,i1,id1,status)

    type(fitsfiles), intent(in) :: fits
    real(REAL64), intent(in) :: epoch
    character(len=*), dimension(:), intent(in) :: com
    real(REAL64), dimension(:), intent(in) :: alpha,delta,x,y,ares,dres
    integer, dimension(:), intent(in) :: i1,id1
    integer, intent(in out) :: status

    integer, parameter :: ncols = 7
    integer, parameter :: frow = 1
    character(len=FLEN_VALUE), dimension(ncols) :: ttype, tform, tunit
    integer, dimension(:), allocatable :: idx
    integer :: i, j, n

    n = size(alpha)
    allocate(idx(n))
    ttype = [ character(len=FLEN_VALUE) :: &
         FITS_COL_RA, FITS_COL_DEC, FITS_COL_X, FITS_COL_Y, 'RA_OFF', 'DEC_OFF', &
         'MATCH' ]
    tform = '1D'
    tform(7) = '1J'
    tunit = 'deg'
    tunit(3:4) = ''
    tunit(7) = ''

    ! MATCH of a star (row j, catalogue star i1(j)) is its position
    ! in the match sequence id1, zero for the stars added by the final
    ! cross-match; the sequence can be longer than the rows (stars
    ! of the sequence rejected by the fit), it was used as the row index
    idx = 0
    do j = 1, n
       do i = 1, size(id1)
          if( id1(i) == i1(j) ) then
             idx(j) = i
             exit
          end if
       end do
    end do

    call fits_insert_btbl(fits,ncols,ttype,tform,tunit,ASTROREPORT,status)

    if( IEEE_IS_NORMAL(epoch) ) call fits_write_key(fits,FITS_KEY_EPOCH,epoch,-6, &
         'epoch of equatorial coordinates',status)

    call fits_write_comment(fits,BEGIN_ASTROMETRY,status)
    do i = 1, size(com)
       call fits_write_comment(fits,com(i),status)
    end do
    call fits_write_comment(fits,MUNIPACK_VERSION,status)
    call fits_write_comment(fits,END_ASTROMETRY,status)

    call fits_write_col(fits,1,frow,alpha,status)
    call fits_write_col(fits,2,frow,delta,status)
    call fits_write_col(fits,3,frow,x,status)
    call fits_write_col(fits,4,frow,y,status)
    call fits_write_col(fits,5,frow,ares,status)
    call fits_write_col(fits,6,frow,dres,status)
    call fits_write_col(fits,7,frow,idx,status)

  end subroutine wcsreport

  subroutine wcs_delete_report(fits,status)

    type(fitsfiles), intent(in) :: fits
    integer, intent(in out) :: status
    integer :: hdutype

    call fits_movnam_hdu(fits,FITS_BINARY_TBL,ASTROREPORT,0,status)
    if( status == 0 ) then
       call fits_delete_hdu(fits,hdutype,status)
    else if( status == FITS_BAD_HDU_NUM ) then
       ! no problem, we add the new one
       status = 0
    else
       ! this is a serious problem
       return
    end if

  end subroutine wcs_delete_report

end module astrofits
