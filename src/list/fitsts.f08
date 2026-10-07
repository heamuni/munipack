!
!  FITS tool for timeseries
!
!  Copyright © 2018 - 2026 F.Hroch (hroch@physics.muni.cz)
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

module fitsts

  use titsio
  use iso_fortran_env

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: rad = 57.295779513082322865_dbl

  type :: TSFits
     character(len=FLEN_FILENAME) :: filename
     character(len=FLEN_VALUE) :: dateobs = ''
     character(len=FLEN_VALUE) :: radesys = 'ICRS'
     character(len=FLEN_VALUE), dimension(:), allocatable :: &
          headvals, cols, units, keyunits
     character(len=FLEN_KEYWORD), dimension(:), allocatable :: &
          headkeys
     character, dimension(:), allocatable :: keytype
     character(len=FLEN_VALUE), dimension(2) :: ctype
     real(dbl), dimension(2) :: crval, crpix, crerr, geo
     real(dbl), dimension(2,2) :: cd
     real(dbl), dimension(:,:), allocatable :: data
     logical, dimension(:), allocatable :: valid
     real :: exptime = -1
     real :: scale = 0, hwhm = 1
     integer :: year
     real(dbl) :: epoch = 2000.0
     real(dbl) :: jd = 2451545.0_dbl   ! J2000.0
     real(dbl) :: azimuth = 0, zenitd = 0, airmass = 1
     logical :: status = .false.
     logical :: cmatch = .true.
     logical :: key_found = .true.
   contains
     procedure :: Load
     procedure, nopass, private :: FitsOpen
     procedure, private :: updatejd
     procedure :: gettime
     procedure :: gettol
     procedure, private :: updatehorizon
     procedure, private :: updateairmass
     procedure, private :: crossmatch
     procedure :: report
  end type TSFits

  private :: tits_read_key

contains

  subroutine Load(this,filename,extname,fitskeys,headkeys,coonames,cootype, &
       colnames,tolerance,cfits,verbose)

    use fitscat

    class(TSFits) :: this
    character(len=*), intent(in) :: filename, extname, cootype
    character(len=*), dimension(:), intent(in) :: fitskeys, headkeys, &
         coonames,colnames
    real(dbl) :: tolerance
    type(CatFits), intent(in) :: cfits
    logical, intent(in) :: verbose

    real(dbl), dimension(:,:), allocatable :: coo, data

    call FitsOpen(this,filename,extname,fitskeys,headkeys,coonames,cootype, &
         colnames,coo,data,verbose)
    if( .not. this%status ) return

    call this%updatejd
    call this%updatehorizon
    call this%updateairmass

    ! select stars
    block
      character(len=80) :: msg
      character(len=3) :: units
      integer, dimension(:), allocatable :: idx
      real(dbl), dimension(:), allocatable :: xref,yref
      integer :: i,j,n,stat,ncols
      real(dbl) :: tol

      n = size(cfits%coo,1)
      ncols = size(colnames)
      allocate(this%data(n,ncols),this%valid(n),xref(n),yref(n),idx(n), &
           stat=stat,errmsg=msg)
      if( stat /= 0 ) then
         write(error_unit,*) trim(msg)
         error stop 'Failed to allocate images.'
      end if

      tol = this%gettol(cootype,tolerance)
      call cfits%getcurrent(this%jd,xref,yref) ! coordinates at current time
      call this%crossmatch(cootype,tol,xref,yref,coo(:,1),coo(:,2),idx)

      if( cootype == 'RECT' ) then
         units = 'pix'
      else
         units = 'deg'
      end if

      do i = 1, size(idx)
         j = idx(i)
         if( j > 0 ) then
            this%data(i,:) = data(j,:)
            this%valid(i) = .true.
         else
            this%cmatch = .false.
            this%data(i,:) = 0 !data_undefined
            this%valid(i) = .false.
            if( verbose ) then
               write(error_unit,'(2a,2(f0.3,1x),a,es8.1,2a)') &
                    trim(filename),': * @ ',xref(i),yref(i), &
                    'not found (>',tol,units,').'
            end if
         end if
      end do
    end block

  end subroutine Load

  subroutine FitsOpen(fits,filename,extname,fitskeys,headkeys,coonames, &
       cootype,colnames,coo,data,verbose)

    type(TSFits),intent(out) :: fits
    character(len=*), intent(in) :: filename, extname, cootype
    character(len=*), dimension(:), intent(in) :: fitskeys, headkeys, &
         coonames,colnames
    real(dbl), dimension(:,:), allocatable, intent(out) :: coo, data
    logical, intent(in) :: verbose

    real(REAL64), parameter :: null = real(0.0,REAL64)
    logical, dimension(size(headkeys)) :: headkey_found
    integer, parameter :: extver = 0
    character(len=FLEN_KEYWORD) :: key
    character(len=80) :: msg
    real(dbl) :: sqscale
    integer :: i, ncols, nrows, ndata, colnum, status, stat, nkeys
    type(fitsfiles) :: fitsfile
    logical :: anyf

    status = 0
    fits%filename = filename
    allocate(fits%cols, source=colnames)
    allocate(fits%headkeys, source=headkeys)
    nkeys = size(headkeys)
    ncols = size(colnames)
    allocate(fits%headvals(nkeys),fits%keytype(nkeys),fits%keyunits(nkeys), &
         fits%units(ncols))
    fits%units = ''
    fits%headvals = ''
    fits%keytype = 'C'
    fits%keyunits = ''
    headkey_found = .true.

    call fits_open_image(fitsfile,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) 'Error: failed to read the file `',trim(filename),"'."
       fits%status = .false.
       return
    end if

    ! read date and time, dateobs
    call fits_write_errmark
    call fits_get_dateobs(fitsfile,fitskeys(1:2),fits%dateobs,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       if( verbose ) write(error_unit,*) 'Warning: Date of observation by ', &
            trim(fitskeys(1)),',',trim(fitskeys(2)),' keywords not found.'
       status = 0
       fits%dateobs = ''
    end if

    ! exposure time
    call fits_write_errmark
    call fits_read_key(fitsfile,fitskeys(3),fits%exptime,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       if( verbose ) write(error_unit,*) &
            'Warning: Exposure time by ',trim(fitskeys(3)),' keyword not found.'
       status = 0
       fits%exptime = -1
    end if

    if( cootype == 'EQT' ) then

       call fits_read_wcs(fitsfile,fits%radesys,fits%epoch,fits%ctype, &
            fits%crval,fits%crpix,fits%cd,fits%crerr,status)
       if( status == 0 ) then
          sqscale = fits%cd(1,1)**2 + fits%cd(1,2)**2
          if( sqscale > 0 )then
             fits%scale = real(1.0 / sqrt(sqscale))
          else
             fits%scale = 0
          end if
       end if

    else ! if( cootype == 'RECT' ) then
       fits%ctype = ''
       fits%crval = 0
       fits%crpix = 0
       fits%cd = 0
       fits%crerr = 0
    end if

    ! geographical coordinates
    call fits_write_errmark
    call fits_read_key(fitsfile,fitskeys(4),fits%geo(1),status)
    call fits_read_key(fitsfile,fitskeys(5),fits%geo(2),status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       if( verbose ) write(error_unit,*) 'Warning: Geographic coordinates by ',&
            trim(fitskeys(4)),",",trim(fitskeys(5)),' keywords not found.'
       status = 0
       fits%geo = 0
    end if

    ! user required keywords
    call fits_write_errmark
    do i = 1, size(headkeys)

       call tits_read_key(fitsfile,headkeys(i),fits%headvals(i),fits%keyunits(i),&
            fits%keytype(i),status)
       headkey_found(i) = status == 0
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          fits%headvals(i) = ''
          fits%keyunits(i) = ''
          fits%keytype(i) = 'C'
          status = 0
       end if

    end do
    call fits_clear_errmark

    ! select data extension
    call fits_write_errmark
    call fits_movnam_hdu(fitsfile,FITS_BINARY_TBL,extname,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       call fits_clear_errmark
       if( verbose .and. size(colnames) > 0 ) write(error_unit,*)  &
            "Warning: A photometry extension is unavailable: ",trim(filename)
       goto 666
    end if

    call fits_write_errmark
    call fits_read_key(fitsfile,FITS_KEY_HWHM,fits%hwhm,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       status = 0
       fits%hwhm = 1
    end if

    ! user required keywords
    call fits_write_errmark
    do i = 1, size(headkeys)

       ! if the same keyword has been already found in any previous HDU,
       ! it is omitted
       if( .not. headkey_found(i) ) then
          call tits_read_key(fitsfile,headkeys(i),fits%headvals(i),fits%keyunits(i),&
               fits%keytype(i),status)
          headkey_found(i) = status == 0
          if( status == FITS_KEYWORD_NOT_FOUND ) then
             fits%headvals(i) = ''
             fits%keyunits(i) = ''
             fits%keytype(i) = 'C'
             status = 0
          end if
       end if
    end do
    call fits_clear_errmark

    ! read table data
    ncols = size(coonames)
    ndata = size(colnames)

    call fits_get_num_rows(fitsfile,nrows,status)
    if( status /= 0 ) goto 666

    ! units
    call fits_write_errmark
    do i = 1, ndata
       call fits_get_colnum(fitsfile,.true.,colnames(i),colnum,status)
       call fits_make_keyn('TUNIT',colnum,key,status)
       call fits_read_key(fitsfile,key,fits%units(i),status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          fits%units(i) = ''
          status = 0
       end if
    end do
    call fits_clear_errmark

    allocate(coo(nrows,ncols),data(nrows,ndata),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Failed to allocate volume for data.'
    end if

    do i = 1, ncols
       call fits_get_colnum(fitsfile,.true.,coonames(i),colnum,status)
       call fits_read_col(fitsfile,colnum,1,null,coo(:,i),anyf,status)
    end do

    do i = 1, ndata
       call fits_get_colnum(fitsfile,.true.,colnames(i),colnum,status)
       call fits_read_col(fitsfile,colnum,1,null,data(:,i),anyf,status)
    end do

666 continue

    call fits_close_file(fitsfile,status)
    call fits_report_error(error_unit,status)

    fits%status = status == 0
    fits%key_found = all(headkey_found .eqv. .true.)

    ! inform user when any keyword is empty
    if( verbose .and. .not. fits%key_found ) then
       write(error_unit,'(a)',advance='no') 'Warning: keys not found: '
       do i = 1, size(headkeys)
          if( .not. headkey_found(i) ) &
               write(error_unit,"(a,' ')",advance='no') trim(headkeys(i))
       end do
       write(error_unit,*)
    end if

  end subroutine FitsOpen

  subroutine updatejd(fits)

    use trajd

    class(TSFits) :: fits
    integer ::year, month, day, hour, minute
    real(dbl) :: second, d
    integer :: status

    if( fits%dateobs /= '' ) then

       ! decode date-time string
       status = 0
       call fits_str2date(fits%dateobs,year,month,day, &
            hour,minute,second,status)
       if( status /= 0 ) return

       ! Evaluate Julian date
       d = day + (hour + (minute + second/60.0_dbl) / 60.0_dbl) / 24.0_dbl
       fits%jd = datjd(real(year,dbl),real(month,dbl),d)
       fits%year = year

    end if

  end subroutine updatejd


  function gettime(fits,timetype,timestamp,min0,per) result(t)

    use trajd
    use astrosphere

    class(TSFits), intent(in) :: fits
    character(len=*), intent(in) :: timetype, timestamp
    real(dbl), intent(in) :: min0,per

    real(dbl) :: t,jd,etime,ls,d

    if( fits%dateobs /= '' ) then

       jd = fits%jd

       ! time-stamps
       if( fits%exptime >= 0) then

          etime = fits%exptime / 86400.0_dbl
          if( timestamp == 'MID' ) then
             jd = jd + etime / 2.0_dbl
          else if ( timestamp == 'END' ) then
             jd = jd + etime
          end if

       end if

       ! time (-like) determination

       if( timetype == 'JD' ) then
          t = jd
       else if( timetype == 'MJD' ) then
          t = mjd(jd)
       else if( timetype == 'HJD' ) then
          d = jd - datjd(real(fits%year,dbl),1.0_dbl,1.0_dbl)
          ls = longsun(d)
          t = hjd(jd,helcor(fits%crval(1),fits%crval(2),ls))
       else if( timetype == 'PHASE' ) then
          t = phase(jd,min0,per)
       else
          t = 0
       end if

    else
       t = 0
    end if

  end function gettime

  function gettol(fits,cootype,tolerance) result(tol)

    class(TSFits), intent(in) :: fits
    character(len=*), intent(in) :: cootype
    real(dbl), intent(in) :: tolerance

    real(dbl) :: tol

    if( cootype == 'EQT' ) then
       if( tolerance > 0 ) then
          tol = tolerance
       else if( maxval(fits%crerr) > 0 .and. fits%scale > 0 ) then
          ! If tol is not given, maximum of both astrometry error
          ! and pixel size is used. We are suppose that user
          ! specifies coordinates of objects with corresponding
          ! precision.
          tol = 5*max(maxval(fits%crerr),1/fits%scale)
       else if( fits%hwhm > 0 .and. fits%scale > 0 ) then
          tol = 5*fits%hwhm / fits%scale
       else if( fits%scale > 0 ) then !
          tol = 5/fits%scale
       else ! fail-back, 1 arcsec
          tol = 1.0 / 3600.0
       end if
    else if( cootype == 'RECT' ) then
       if( tolerance > 0 ) then
          tol = tolerance
       else
          tol = 5 ! pixels
       end if
    else
       write(error_unit,*) 'Error: tolerance is unset.'
    end if

  end function gettol

  subroutine updatehorizon(fits)

    use astrosphere

    class(TSFits) :: fits
    real(dbl) :: ra, dec, t, ha, a, h, longitude, latitude

    ra = fits%crval(1)
    dec = fits%crval(2)
    longitude = fits%geo(1)
    latitude = fits%geo(2)
    t = lmst(fits%jd,longitude)
    ha = hangle(15*t, ra)
    call eq2hor(ha, dec, latitude, a, h)
    fits%azimuth = a
    fits%zenitd = real(90.0,dbl) - h

  end subroutine updatehorizon

  subroutine updateairmass(fits)

    use astrosphere

    class(TSFits) :: fits

    fits%airmass = airmass(fits%zenitd)

  end subroutine updateairmass


  subroutine crossmatch(fits,cootype,tol,u,v,x,y,idx)

    ! Returns idx with valid elements > 0 if the given
    ! reference object (u,v) has been found at coordinates (x,y).
    ! The search is made inside a circle with tol radius.
    ! A closer object is selected when multiple objects falls
    ! into the circle.

    class(TSFits), intent(in) :: fits
    character(len=*), intent(in) :: cootype
    real(dbl), intent(in) :: tol
    real(dbl), dimension(:), intent(in) :: u,v,x,y
    integer, dimension(:), intent(out) :: idx

    real(dbl) :: g, rmin, r, dx, dy
    integer :: i, j

    idx = -1

    ! crossmatch is based on distance
    ! metrics for rectangular coordinates is
    !  dr**2 = dx**2 + dy**2
    ! while for spherical
    !  dr**2 = dalpha**2*cos(delta)**2 + ddelta**2
    ! the metric works satisfactory for small angular differences
    if( cootype == 'EQT' ) then
       g = cos(fits%crval(2)/rad)
    else
       g = 1
    end if
    do i = 1, size(u)
       rmin = tol
       do j = 1,size(x)
          dx = u(i) - x(j)
          dy = v(i) - y(j)
          r = hypot(dx*g, dy)
          if( r < rmin ) then
             rmin = r
             idx(i) = j
          end if
       end do
    end do

  end subroutine crossmatch


  subroutine report(fits,timetype,timestamp,lc_epoch,lc_period, &
       printname,printtime,printhorizon,printairmass)

    class(TSFits), intent(in) :: fits
    character(len=*), intent(in) :: timetype, timestamp
    real(dbl), intent(in) :: lc_epoch,lc_period
    logical, intent(in) :: printname, printtime, printhorizon, printairmass

    real(dbl), dimension(:,:), allocatable :: data
    character(len=666) :: fmt = "", form = ""
    integer :: i,n,m
    real(dbl) :: tm


    ! filename
    if( printname ) write(*,'(a,1x)',advance="no") trim(fits%filename)

    ! time
    if( printtime ) then
       tm = fits%gettime(timetype,timestamp,lc_epoch,lc_period)
       write(*,'(f0.10,1x)',advance="no") tm
    end if

    ! header values
    if( allocated(fits%headvals) ) then
       associate( vals => fits%headvals )
         if( size(vals) > 0 ) then
            write(fmt,'(a,i0,a)') '(',size(vals),'(a,1x))'
            write(*,fmt,advance="no") (trim(vals(i)),i=1,size(vals))
         end if
       end associate
    end if

    ! horizontal coordinates and airmass
    if( printhorizon ) &
       write(*,'(2(f0.10,1x))',advance="no") fits%azimuth,fits%zenitd
    if( printairmass ) &
       write(*,'(f0.10,1x)',advance="no") fits%airmass

    if( allocated(fits%data) ) then

       n = size(fits%data,1) ! objects
       m = size(fits%data,2) ! columns

       if( m > 0 .and. n > 0 ) then

          ! format for one object derived by columns
          fmt = ''
          do i = 1, m
             if( index(fits%cols(i),'MAG') > 0 ) then
                fmt = trim(fmt) // 'f0.5'
             else
                fmt = trim(fmt) // 'g0'
             end if
             fmt = trim(fmt) // ',1x'
             if( i < m ) fmt = trim(fmt) // ','
          end do
          write(form,'(a,i0,3a)') '(1x,',n,'(',trim(fmt),'))'

          ! temporary memoty is allocated due undefined values
          allocate(data(n,m))
          data = fits%data
          do i = 1, n
             if( .not. fits%valid(i) ) then
                where( fits%cols == 'MAG' )
                   data(i,:) = 99.999
                elsewhere( fits%cols == 'MAGERR' )
                   data(i,:) = 9.999
                elsewhere
                   data(i,:) = huge(data)*0
                end where
             end if
          end do
          write(*,form,advance="no") (data(i,:),i=1,n)

       end if
    end if
    write(*,*)

  end subroutine report

  subroutine tits_read_key(fitsfile,key,value,units,type,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: key
    character(len=*), intent(out) :: value,units,type
    integer, intent(in out) :: status

    character(len=FLEN_COMMENT) :: com
    character(len=FLEN_VALUE) :: val

    call fits_read_key(fitsfile,key,value,status)
    if( status == 0 ) then
       call fits_read_key_unit(fitsfile,key,units,status)
       call fits_read_keyword(fitsfile,key,val,com,status)
       call fits_get_keytype(val,type,status)
    else
       value = ''
       units = ''
       type = 'C'
    end if

  end subroutine tits_read_key

end module fitsts
