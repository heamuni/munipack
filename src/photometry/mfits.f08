!
!  FITS I/O for photometric calibration
!
!  Copyright © 2012-2026 F.Hroch (hroch@physics.muni.cz)
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

module mfits

  use iso_fortran_env
  use titsio

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: rad = 57.295779513082322865_dbl
  real(dbl), parameter, private :: pi = 3.14159265359_dbl

  ! the photometry method selected by phcal: APER, GROW, PSF, or empty
  ! (GROW when its extension is present, otherwise APER)
  character(len=FLEN_VALUE) :: phmethod = ''

  ! the aperture photometry is divided by the enclosed energy of
  ! the aperture (APEEi of APERPHOT, when present): the total flux
  logical :: apee_use = .true.


contains


  subroutine readcat(cat,labels,label_mag,label_magerr, &
       alpha,delta,mag,magerr,catid,status)

    character(len=*), intent(in) :: cat
    character(len=*), dimension(:), intent(in) :: labels,label_mag,label_magerr
    real(dbl), dimension(:), allocatable, intent(out) :: alpha,delta
    real(dbl), dimension(:,:), allocatable, intent(out) :: mag,magerr
    character(len=*), intent(out) :: catid
    integer, intent(out) :: status

    real(dbl), parameter :: nullcoo = -999
    real(dbl), parameter :: nullmag = 99.99999
    real(dbl), parameter :: dullmag = 9.99999

    integer :: nrows, srows, ncols, i, k, l, frow
    integer, dimension(size(labels)) :: cols
    integer, dimension(size(label_mag)) :: col_mag
    integer, dimension(size(label_magerr)) :: col_magerr
    type(fitsfiles) :: fits
    logical :: anyf


    status = 0

    ! open and move to a table extension
    call fits_open_table(fits,cat,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read a table in the file `',trim(cat),"'."
       return
    end if

    call fits_get_num_rows(fits,nrows,status)
    if( status /= 0 ) goto 666
    if( .not. (nrows > 0) ) then
       write(error_unit,*) 'Error: an empty table in the file `',trim(cat),"'."
       goto 666
    end if

    ! define reference frame and identification of catalogue
    call fits_read_key(fits,'EXTNAME',catid,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       catid = ''
       status = 0
    end if

    ! find columns by labels
    do i = 1, size(labels)
       call fits_get_colnum(fits,.true.,labels(i),cols(i),status)
    end do
    do i = 1, size(label_mag)
       call fits_get_colnum(fits,.true.,label_mag(i),col_mag(i),status)
    end do
    do i = 1, size(label_magerr)
       call fits_get_colnum(fits,.true.,label_magerr(i),col_magerr(i),status)
    end do
    if( status /= 0 ) goto 666

    ncols = size(label_mag)
    allocate(alpha(nrows),delta(nrows),mag(nrows,ncols),magerr(nrows,ncols))

    magerr = 9.99999

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,cols(1),frow,nullcoo,alpha(i:l),anyf,status)
       call fits_read_col(fits,cols(2),frow,nullcoo,delta(i:l),anyf,status)

       do k = 1, size(col_mag)
          call fits_read_col(fits,col_mag(k),frow,nullmag,mag(i:l,k),anyf,status)
       end do
       do k = 1, size(col_magerr)
          call fits_read_col(fits,col_magerr(k),frow,dullmag,magerr(i:l,k),anyf,status)
       end do
       if( status /= 0 ) goto 666
    end do


    call fits_close_file(fits,status)
    return

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( allocated(mag) ) deallocate(mag,magerr,alpha,delta)

  end subroutine readcat


  subroutine readframe(filename,keys,threshold,maxerr,alpha,delta,cts,dcts, &
       ftol,exptime,area,photosys,filter,init_area,saper,aper,verbose,status,&
       amass)

    use astrotrafo
    use astrosphere
    use trajd
!    use phio

    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: keys
    real(dbl), intent(in) :: threshold, maxerr
    real(dbl), dimension(:), allocatable, intent(out) :: alpha,delta,cts,dcts
    real(dbl), intent(out) :: ftol,exptime, area
    character(len=*), intent(out) :: photosys,filter
    logical, intent(in) :: init_area
    integer, intent(in) :: saper
    real(dbl), intent(out) :: aper
    logical, intent(in) :: verbose
    integer, intent(out) :: status
    real(dbl), optional, intent(out) :: amass

    integer, parameter :: DIM = 2
    integer, dimension(DIM) :: naxes
    integer :: naxis, nrows, srows, naper, maper, s1, s2, s4, s5, i, j, l, nid, &
         frow, nkeys
    real(dbl) :: jd, longitude, latitude, spot, scale, rot, refl
    real(dbl), dimension(DIM) :: crder, crval
    character(len=FLEN_CARD) :: key, keycts, keyerr
    character(len=FLEN_CARD), dimension(0) :: exclist
    real(dbl), dimension(:), allocatable :: x,y,sky,ct,dct
    integer, dimension(:), allocatable :: id
    integer, dimension(5) :: colnum
    character(len=:), allocatable :: header
    character(len=FLEN_VALUE) :: dateobs
    real(dbl) :: nullval
    real :: hwhm
    logical :: anyf
    logical :: undef_airmass, have_grow, have_ap, have_psf, aperture
!    type(AstroTrafoProj) :: tproj
    type(fitsfiles) :: fits
    type(AstroWCS) :: wcs


    ! input FITS file
    status = 0
    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read an image in the file `',trim(filename),"'."
       return
    end if

    call fits_get_img_dim(fits,naxis,status)
    call fits_get_img_size(fits,DIM,naxes,status)
    if( status /= 0 ) goto 666

    ! read astrometric calibration
    call fits_hdr2str(fits,1,exclist,nkeys,header,status)
    if( status == 0 ) then
       call wcs % init(header,nkeys,status)
!       block
!         real(REAL64) :: rot, refl
!       end block
!       if( status /= 0 ) then
    end if
    if( allocated(header) ) deallocate(header)
!    else
!    call wcsget(fits,tproj,status)
!    if( status == FITS_KEYWORD_NOT_FOUND ) then
    if( status /= 0 ) then
       write(error_unit,*) 'Astrometry keywords not found in header.'
       goto 666
    end if
!    ftol = 5.0*tproj%err
    call wcs % get_affines(scale,rot,refl)
    call wcs % get_crder(crder)
    ftol = 5 * maxval(crder)
!    stop

    ! Citizen date
    call fits_get_dateobs(fits,[keys(3),keys(8)],dateobs,status)
    if( status == 0 ) then
       jd = fits_jd(dateobs,status)
    else
       write(error_unit,*) "Warning: Failed to get date or time."
       goto 666
    end if

    ! Geographic coordinates
    call fits_read_key(fits,keys(6),longitude,status)
    call fits_read_key(fits,keys(7),latitude,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       if( verbose ) write(error_unit,*) &
            "Warning: Geographical longitude or latitude identified by FITS keywords `",&
            trim(keys(6)),"',`",trim(keys(7)),"' not found. Airmass is undefined."
       status = 0
       undef_airmass = .true.
    else
       undef_airmass = .false.
    end if

    s1 = 0
    call fits_read_key(fits,keys(1),exptime,s1)
    if( s1 == FITS_KEYWORD_NOT_FOUND ) then
       write(error_unit,*) "Warning: Exposure time identified by FITS keyword `",&
            trim(keys(1)),"' not found. Set to default 1 sec."
       exptime = 1
       s1 = 0
    end if

    s2 = 0
    if( .not. init_area ) then
       call fits_read_key(fits,keys(2),area,s2)
       if( s2 == FITS_KEYWORD_NOT_FOUND ) then
          write(error_unit,*) "Warning: Area identified by FITS keyword `", &
               trim(keys(2)),"' not found. Set to default 1 sq. meter."
          area = 1
          s2 = 0
       end if
    end if

!    s3 = 0
!    if( present(airmass) ) then
!       call ftgkyd(20,keys(3),airmass,com,s3)
!       if( s3 == KEYWORD_NOT_FOUND ) then
!          write(error_unit,*) "Warning: Airmass identified by FITS keyword `", &
!               trim(keys(3)),"' not found. Default is 0 (extra-atmospheric)."
!          airmass = 0
!          s3 = 0
!       end if
!    end if

    s4 = 0
    call fits_read_key(fits,keys(4),photosys,s4)
    if( s4 == FITS_KEYWORD_NOT_FOUND ) then
!       write(error_unit,*) "Photometry system by FITS keyword `", &
!            trim(keys(4)),"' not found, leaving it empty."
       photosys = ''
       s4 = 0
    end if

    s5 = 0
    call fits_read_key(fits,keys(5),filter,s5)
    if( s5 == FITS_KEYWORD_NOT_FOUND ) then
       write(error_unit,*) "Filter identified by FITS keyword `", &
            trim(keys(5)),"' not found. Set to default `' (empty)."
       filter = ''
       s5 = 0
    end if

    status = s1 + s2 + s4 + s5
    if( status /= 0 ) goto 666

    ! Julian date at midpoint of exposure
    jd = jd + exptime / 2.0_dbl / 86400.0_dbl

    if( present(amass) ) then
       if( .not. undef_airmass ) then
          ! airmass at the time
          call wcs % get_crval(crval)
          amass = xairmass(jd,longitude,latitude,crval(1),crval(2))
!          amass = xairmass(jd,longitude,latitude,tproj%acen,tproj%dcen)
       else
          amass = -1
       end if
    end if

    ! detect all possible photometry tables
!    phexts(1) = APEREXTNAME
!    phexts(2) = GROWEXTNAME
!    have_phext = .false.
!    do i = 1,size(phexts)
!       call ftmnhd(20,BINARY_TBL,phexts(i),0,status)
!       if( status == 0 ) have_phext(i) = .true.
!       if( status == BAD_HDU_NUM ) status = 0
!       if( status /= 0 ) goto 666
!    end do



    have_psf = phmethod == 'PSF'
    if( have_psf ) then

       ! the PSF photometry extension
       call fits_movnam_hdu(fits,FITS_BINARY_TBL,PSFEXTNAME,0,status)
       if( status == FITS_BAD_HDU_NUM ) then
          write(error_unit,*) &
               "Error: Failed to find a PSF photometry extension: ", &
               trim(filename)
          goto 666
       end if
       have_grow = .false.
       have_ap = .false.

    else

    ! try to select grow-curve photometry table
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,GROWEXTNAME,0,status)
    have_grow = .not. (status == FITS_BAD_HDU_NUM)
    if( phmethod == 'GROW' .and. .not. have_grow ) then
       write(error_unit,*) &
            "Error: Failed to find a growth-curve photometry extension: ", &
            trim(filename)
       goto 666
    end if
    if( phmethod == 'APER' ) have_grow = .false.

    ! try to select the aperture photometry extension
    if( .not. have_grow ) then
       status = 0
       call fits_movnam_hdu(fits,FITS_BINARY_TBL,APEREXTNAME,0,status)
       have_ap = .not. (status == FITS_BAD_HDU_NUM)
    else
       have_ap = .true.
    end if

    if( .not. (have_grow .or. have_ap) )then
       write(error_unit,*) &
            "Error: Failed to find an aperture or a growth photometry extension: ", &
            trim(filename)
       goto 666
    end if

    end if

    ! key for HWHM is presented in all photometry extensions
    ! (methods for their estimates are different)
    call fits_read_key(fits,FITS_KEY_HWHM,hwhm,status)
    !    ftol = hwhm * tproj%scale
    ftol = hwhm /  scale

    keycts = ''
    keyerr = ''
    aperture = .false.

    if( have_psf ) then
       keycts = FITS_COL_PSFCOUNT
       keyerr = FITS_COL_PSFCOUNTERR
       ! the area of the selection threshold: the circle of radius FWHM
       aper = 2*hwhm
    end if

    if( have_grow .and. (saper == 0 .or. phmethod == 'GROW') ) then
       keycts = FITS_COL_GCOUNT
       keyerr = FITS_COL_GCOUNTERR
       call fits_read_key(fits,FITS_KEY_RF90,aper,status)
    end if

    if( (.not. have_grow .and. have_ap) .or. &
         (have_ap .and. saper > 0 .and. phmethod /= 'GROW') ) then

       call fits_movnam_hdu(fits,FITS_BINARY_TBL,APEREXTNAME,0,status)

       if( saper > 0 ) then
          maper = saper
          call fits_make_keyn(FITS_KEY_APER,maper,key,status)
          call fits_read_key(fits,key,aper,status)
       else
          call fits_read_key(fits,FITS_KEY_NAPER,naper,status)
          do i = 1, naper
             call fits_make_keyn(FITS_KEY_APER,i,key,status)
             call fits_read_key(fits,key,aper,status)
             ! the first aperture larger than 2*sigma of Gaussian profile
             if( aper > 2*hwhm / sqrt(2*log(2.0)) ) then
                maper = i
                exit
             end if
          end do
       end if

       write(keycts,'(a,i0)') FITS_COL_APCOUNT,maper
       write(keyerr,'(a,i0)') FITS_COL_APCOUNTERR,maper
       aperture = .true.

    end if

    call fits_get_num_rows(fits,nrows,status)
    if( status /= 0 ) goto 666

    allocate(x(nrows),y(nrows),sky(nrows),ct(nrows),dct(nrows))

    ! table
    call fits_get_colnum(fits,.true.,FITS_COL_X,colnum(1),status)
    call fits_get_colnum(fits,.true.,FITS_COL_Y,colnum(2),status)
    call fits_get_colnum(fits,.true.,FITS_COL_SKY,colnum(3),status)
    call fits_get_colnum(fits,.true.,keycts,colnum(4),status)
    call fits_get_colnum(fits,.true.,keyerr,colnum(5),status)

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,colnum(1),frow,nullval,x(i:l),anyf,status)
       call fits_read_col(fits,colnum(2),frow,nullval,y(i:l),anyf,status)
       call fits_read_col(fits,colnum(3),frow,nullval,sky(i:l),anyf,status)
       call fits_read_col(fits,colnum(4),frow,nullval,ct(i:l),anyf,status)
       call fits_read_col(fits,colnum(5),frow,nullval,dct(i:l),anyf,status)
       if( status /= 0 ) goto 666
    end do

    ! the total flux by the enclosed energy of the aperture
    if( aperture .and. apee_use ) then
       call apee_apply(fits,maper,ct,dct,status)
       if( status /= 0 ) goto 666
    end if

    ! PSF photometry: stars with a saturated core or a failed fit
    ! (FLAG 1, 16) are not used
    if( have_psf ) then
       block
         integer, dimension(:), allocatable :: flag
         integer :: fcol
         allocate(flag(nrows))
         call fits_get_colnum(fits,.true.,'FLAG',fcol,status)
         do i = 1, nrows, srows
            frow = i
            l = min(i+srows,nrows)
            call fits_read_col(fits,fcol,frow,0,flag(i:l),anyf,status)
         end do
         if( status /= 0 ) goto 666
         where( iand(flag,17) /= 0 ) ct = -1
       end block
    end if

    ! select stars acceptable for photometry
    spot = pi*aper**2
    allocate(id(nrows))
    nid = 0
    do i = 1, nrows
       if(  ct(i) > 0 .and. dct(i) > 0 .and. sky(i) > 0 .and. &
            ct(i) > threshold * (sky(i)*spot) .and. &
            dct(i) < maxerr * ct(i) ) then
          nid = nid + 1
          id(nid) = i
       end if
    end do

    allocate(alpha(nid),delta(nid),cts(nid),dcts(nid))
    do i = 1, nid
       j = id(i)
       x(i) = x(j)
       y(i) = y(j)
       cts(i) = ct(j)
       dcts(i) = dct(j)
    end do
!    call invtrafo(tproj,x(1:nid),y(1:nid),alpha,delta)
    call wcs % pix2cel(x(1:nid),y(1:nid),alpha,delta)
    deallocate(x,y,ct,dct,sky)

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine readframe

  subroutine apee_apply(fits,maper,ct,dct,status,found,source)

    ! The counts in the aperture maper (of the current APERPHOT HDU)
    ! are divided by its enclosed energy APEE (when the column is present):
    ! the total flux; the relative error of APEE (APEEERR) is added
    ! in quadrature. Stars without a valid APEE are not used.

    type(fitsfiles), intent(in out) :: fits
    integer, intent(in) :: maper
    real(dbl), dimension(:), intent(in out) :: ct, dct
    integer, intent(in out) :: status
    logical, intent(out), optional :: found
    character(len=*), intent(out), optional :: source

    real(dbl), dimension(:), allocatable :: ee, eeerr
    character(len=FLEN_CARD) :: key
    integer :: c1, c2, srows, frow, i, l, nrows
    real(dbl) :: nullval
    logical :: anyf

    if( present(found) ) found = .false.
    if( present(source) ) source = ''
    if( status /= 0 ) return

    write(key,'(a,i0)') FITS_COL_APEE,maper
    call fits_write_errmark
    call fits_get_colnum(fits,.true.,key,c1,status)
    if( status == FITS_COLUMN_NOT_FOUND ) then
       call fits_clear_errmark
       status = 0
       return
    end if
    if( present(source) ) then
       call fits_read_key(fits,'EESOURCE',source,status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
          source = ''
       end if
    end if
    write(key,'(a,i0)') FITS_COL_APEEERR,maper
    call fits_get_colnum(fits,.true.,key,c2,status)
    if( status /= 0 ) return

    nrows = size(ct)
    allocate(ee(nrows),eeerr(nrows))
    nullval = -1
    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)
       call fits_read_col(fits,c1,frow,nullval,ee(i:l),anyf,status)
       call fits_read_col(fits,c2,frow,nullval,eeerr(i:l),anyf,status)
    end do
    if( status /= 0 ) return

    where( ct > 0 .and. dct > 0 .and. ee > 0 .and. eeerr >= 0 )
       dct = sqrt((dct / ee)**2 + (ct * eeerr / ee**2)**2)
       ct = ct / ee
    elsewhere
       ct = -1
       dct = -1
    end where
    if( present(found) ) found = .true.

  end subroutine apee_apply


  subroutine checkframe(filename,keys,filter,nrows,status)

    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: keys
    character(len=*), intent(out) :: filter
    integer, intent(out) :: nrows
    integer, intent(in out) :: status

    type(fitsfiles) :: fits

    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) return

    call fits_read_key(fits,keys(5),filter,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       write(error_unit,*) "Filter identified by FITS keyword `", &
            trim(keys(5)),"' not found. Default is `' (empty)."
       filter = ''
       status = 0
    end if

    call fits_movnam_hdu(fits,FITS_BINARY_TBL,APEREXTNAME,0,status)
    if( status == FITS_BAD_HDU_NUM ) then
       write(error_unit,*) "Failed to find a photometry extension: ",trim(filename)
       goto 666
    end if

    call fits_get_num_rows(fits,nrows,status)

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine checkframe


  subroutine readframes(filenames,keys,filters_ref,photosys,threshold, &
       maxerr,ra,dec,cts,dcts,utol,exptime,area,filter,tol,init_area, &
       saper,aper,verbose,status)

    character(len=*), dimension(:), intent(in) :: filenames,keys,filters_ref
    character(len=*), intent(out) :: photosys
    real(dbl), intent(in) :: threshold, maxerr
    real(dbl), dimension(:,:), allocatable, intent(out) :: ra,dec,cts,dcts
    real(dbl), dimension(:), allocatable, intent(out) :: exptime,area,tol
    character(len=FLEN_VALUE), dimension(:), allocatable, intent(out) :: filter
    real(dbl), intent(in) :: utol
    integer, intent(in) :: saper
    logical, intent(in) :: init_area
    real(dbl), dimension(:), allocatable, intent(out) :: aper
    logical, intent(in) :: verbose
    integer, intent(out) :: status

    real(dbl) :: ftol
    integer :: n, nd, nfiles!, i, j
    real(dbl), dimension(:), allocatable :: qra,qdec,dn,ddn
!    integer, allocatable, dimension(:) :: ndat, order
    integer, allocatable, dimension(:) :: order
    character(len=FLEN_VALUE) :: phsys

    nfiles = size(filenames)

    call file_order(filenames,keys,filters_ref,order,nd)

!!$    allocate(filter(nfiles),ndat(nfiles),order(nfiles))
!!$
!!$    ! before loading of the (big) data, we are discovering their sizes and filters
!!$    status = 0
!!$    do n = 1, nfiles
!!$       call checkframe(filenames(n),keys,filter(n),ndat(n),status)
!!$       if( status /= 0 ) then
!!$          write(error_unit,*) "File: `",trim(filenames(n)),"'"
!!$          stop 'Failed to read a frame.'
!!$       end if
!!$    end do
!!$
!!$    ! array allocation dimension
!!$    nd = maxval(ndat)
!!$
!!$    ! and using the info to arrange order by filters
!!$    n = 0
!!$    do i = 1,size(filters_ref)
!!$       do j = 1,size(filter)
!!$          if( filters_ref(i) == filter(j) ) then
!!$             n = n + 1
!!$             order(n) = j
!!$          end if
!!$       end do
!!$    end do
!!$    if( n /= nfiles ) &
!!$         stop 'Filters of frames does not unique corresponds to ones of standard set.'

    ! now, we are prepared to read the data
    allocate(exptime(nfiles),area(nfiles),tol(nfiles),filter(nfiles), &
         cts(nd,nfiles),dcts(nd,nfiles),ra(nd,nfiles),dec(nd,nfiles), &
         aper(nfiles))

    cts = -1.0_dbl
    dcts = -1.0_dbl
    ra = 0.0_dbl
    dec = 0.0_dbl
    aper = -1.0_dbl

    status = 0
    do n = 1, nfiles

       call readframe(filenames(order(n)),keys,threshold,maxerr, &
            qra,qdec,dn,ddn,ftol,exptime(n),area(n),phsys,filter(n),&
            init_area,saper,aper(n),verbose,status)

       if( status /= 0 ) then
          write(error_unit,*) "File: `",trim(filenames(n)),"'"
          stop 'Failed to read a frame.'
       end if

       if( utol > epsilon(utol) ) then
          tol(n) = utol
       else
          tol(n) = ftol
       end if

       if( n == 1 ) then
          photosys = phsys
       else
          if( phsys /= photosys ) write(error_unit,*) &
               "Uncompatible photometry systems: `",trim(phsys),"' and `", &
               trim(photosys),"'."
       end if

       nd = size(dn)

       cts(1:nd,n) = dn
       dcts(1:nd,n) = ddn
       ra(1:nd,n) = qra
       dec(1:nd,n) = qdec

       deallocate(qra,qdec,dn,ddn)

    end do
    deallocate(order)

  end subroutine readframes

  subroutine readref(ref,keys,maxerr,threshold, &
       refra,refdec,refph,drefph,filter,exptime, area, aper, status)

    !
    ! WARNING: Only one band is supported !
    !
    use astrotrafo
!    use phio

    real(dbl), parameter :: rad = 57.295779513082322865_dbl
    real(dbl), parameter :: pi = 3.14159

    character(len=*), intent(in) :: ref
    character(len=*), dimension(:), intent(in) :: keys
    real(dbl), intent(in) :: maxerr,threshold
    real(dbl), dimension(:), allocatable, intent(out) :: refra,refdec
    real(dbl), dimension(:,:), allocatable, intent(out) :: refph,drefph
    character(len=*), intent(out) :: filter
    real(dbl), intent(out) :: exptime, area, aper
    integer, intent(in out) :: status

    real(dbl), dimension(:), allocatable :: ra,dec,ph,dph,sky
    integer, dimension(:), allocatable :: id
    character(len=FLEN_CARD), dimension(0) :: exclist
    character(len=:), allocatable :: header
    character(len=FLEN_CARD) :: key
    character(len=FLEN_VALUE) :: phkind
    real(dbl), parameter :: nullval = 0.0_dbl
    real(dbl) :: sep, r, cosd2, scale, spot, rot, refl
    logical :: anyf, found
    integer :: nrows,srows,nid,rcol,dcol,scol,pcol,ecol,frow,i,j,l,nkeys
    !    type(AstroTrafoProj) :: tastr
    type(fitsfiles) :: fits
    type(AstroWCS) :: wcs

    if( status /= 0 ) return

    ! open first image extension
    call fits_open_image(fits,ref,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read an image in the file `',trim(ref),"'."
       return
    end if

    call fits_read_key(fits,keys(5),filter,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) &
         stop 'A filter keyword not found in reference frame.'

    ! WCS
!    call wcsget(fits,tastr,status)
    call fits_hdr2str(fits,1,exclist,nkeys,header,status)
    if( status == 0 ) then
       call wcs % init(header,nkeys,status)
    end if
    if( allocated(header) ) deallocate(header)
!    if( status == FITS_KEYWORD_NOT_FOUND ) &
    if( status /= 0 ) stop 'Error: Astrometry keywords not found in header.'
    call wcs % get_affines(scale,rot,refl)
    !scale = tastr%scale

    call fits_read_key(fits,keys(1),exptime,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       write(error_unit,*) &
            "Warning: Exposure time identified by FITS keyword `", &
            trim(keys(1)),"' not found (default to 1 sec)."
       exptime = 1
       status = 0
    end if

    call fits_read_key(fits,keys(2),area,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       !       write(error_unit,*)
       !"Warning: Area identified by FITS keyword `",trim(keys(2)), &
!            "' not found (default to 1 [m2])."
       area = 1
       status = 0
    end if

    ! move to photometry table
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,PHOTOEXTNAME,0,status)
    if( status == FITS_BAD_HDU_NUM ) &
         stop 'Failed to find a photometry extension in reference frame.'

    call fits_get_num_rows(fits,nrows,status)

    call fits_make_keyn(FITS_KEY_ANNULUS,1,key,status)
    call fits_read_key(fits,key,sep,status)

    call fits_read_key(fits,FITS_KEY_PHOTOTYP,phkind,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) &
         stop 'PHOTOTYP not found in reference frame.'

    if( phkind == 'APHOT' ) then
       !    sep = 2*sep
       call fits_read_key(fits,FITS_KEY_APER,aper,status)
       aper = aper / scale
    else if( phkind == 'GPHOT' ) then
       call fits_read_key(fits,FITS_KEY_RF90,aper,status)
       aper = aper / scale
    end if

    call fits_get_colnum(fits,.true.,FITS_COL_RA,rcol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_DEC,dcol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_SKY,scol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_PHOTON,pcol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_PHOTONERR,ecol,status)
    if( status /= 0 ) goto 666

    allocate(ra(nrows),dec(nrows),ph(nrows),dph(nrows),sky(nrows),id(nrows))

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,rcol,frow,nullval,ra(i:l),anyf,status)
       call fits_read_col(fits,dcol,frow,nullval,dec(i:l),anyf,status)
       call fits_read_col(fits,scol,frow,nullval,sky(i:l),anyf,status)
       call fits_read_col(fits,pcol,frow,nullval,ph(i:l),anyf,status)
       call fits_read_col(fits,ecol,frow,nullval,dph(i:l),anyf,status)
       if( status /= 0 ) goto 666
    end do

    call fits_close_file(fits,status)

    ! selection of suitable calibration stars
    spot = pi*aper**2    ! sky is per square pixel

    ! selection of isolated stars
    nid = 0
    do i = 1, nrows
       found = .false.
       cosd2 = cos(dec(i) / rad)**2
       do j = 1, nrows
          if( i /= j ) then
             r = sqrt((ra(i) - ra(j))**2*cosd2 + (dec(i) - dec(j))**2)
             if( r < sep ) then
                found = .true.
                goto 90
             end if
          end if
       end do
90     continue

       ! accept only those isolated stars with valid measurements  ...
!       if( .not. found .and. ph(i) > 0 .and. dph(i) > 0 .and. &
!          ! ... and only bright ones above the limit given by sky
!            ph(i) / (spot*sky(i)) > threshold .and. &
!          ! ... with the maximum error
!            dph(i) / ph(i) < maxerr ) then

       if( .not. found .and. ph(i) > 0 .and. dph(i) > 0 )then

          nid = nid + 1
          id(nid) = i

       end if

    end do

    if( nid == 0 ) then
       write(error_unit,*) &
            'Error: There are no suitable and accurate stars',&
            ' on the frame. Try decrease threshold by -th or increase',&
            ' the magnitude error by -e.'
       stop 'No calibration stars available.'
    end if

    allocate(refra(nid),refdec(nid),refph(nid,1),drefph(nid,1))
    do i = 1, nid
       j = id(i)
       refra(i) = ra(j)
       refdec(i) = dec(j)
       refph(i,1) = ph(j)
       drefph(i,1) = dph(j)
    end do

    return

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( allocated(ra) ) deallocate(ra,dec,ph,dph,sky,id)

  end subroutine readref

  subroutine file_order(filenames,keys,filters_ref,order,maxdat)

    character(len=*), dimension(:), intent(in) :: filenames,keys,filters_ref
    integer, dimension(:), allocatable, intent(out) :: order
    integer, intent(out), optional :: maxdat
    character(len=FLEN_VALUE), dimension(:), allocatable :: filter
    integer, dimension(:), allocatable :: ndat
    integer :: i,j,n,status,nfiles,nfilters

    nfilters = size(filters_ref)
    nfiles = size(filenames)
    allocate(filter(nfiles),ndat(nfiles),order(nfiles))

    ! prior to loading of the (big) data,
    ! we are discovering their sizes and filters
    status = 0
    do n = 1, nfiles
       call checkframe(filenames(n),keys,filter(n),ndat(n),status)
       if( status /= 0 ) then
          write(error_unit,*) "File: `",trim(filenames(n)),"'"
          stop 'Failed to read a frame.'
       end if
    end do

    ! and using the info to arrange order by the filters
    n = 0
    do i = 1,nfilters
       do j = 1,size(filter)
          if( filters_ref(i) == filter(j) ) then
             n = n + 1
             order(n) = j
          end if
       end do
    end do

    if( n /= nfiles ) then
       write(error_unit,*)
       write(error_unit,*) 'Summary of filter survey '
       write(error_unit,*) 'Reference filters: ', &
            (trim(filters_ref(i)),i=1,nfilters)
       write(error_unit,*) 'List of files and filters:'
       do i = 1, nfiles
          write(error_unit,*) trim(filenames(i)),': `',trim(filter(i)),"'"
       end do
       stop 'Filters of frames does not corresponds to the standard filter set.'
    end if

    if( present(maxdat) ) then
       ! array allocation dimension
       maxdat = maxval(ndat)
    end if

    deallocate(filter,ndat)

  end subroutine file_order

end module mfits
