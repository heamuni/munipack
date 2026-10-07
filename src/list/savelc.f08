!
!  Save timeserie to FITS table
!
!  Copyright © 2018-20 F.Hroch (hroch@physics.muni.cz)
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

module savelcfits

  use fitsts
  use fitscat
  use titsio
  use iso_fortran_env

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)

  private :: maxlen, tits_write_col_logical, tits_write_col_integer, &
       tits_write_col_double

contains

  subroutine savelc(output,fitses,timetype,timestamp,cat,cols, &
       lc_epoch,lc_period,printname,printtime,printhorizon,printairmass,status)

    character(len=*), intent(in) :: output, timetype, timestamp
    type(TSFits), dimension(:), intent(in) :: fitses
    type(CatFits), intent(in) :: cat
    character(len=*), dimension(:), intent(in) :: cols
    real(dbl), intent(in) :: lc_epoch,lc_period
    logical, intent(in) :: printname, printtime, printhorizon, printairmass
    logical, intent(out) :: status

    character(len=FLEN_VALUE), dimension(:), allocatable :: labels, units
    character(len=FLEN_VALUE), dimension(:,:), allocatable :: headvals
    character(len=FLEN_FILENAME), dimension(size(fitses)) :: filenames
    real(dbl), dimension(size(fitses)) :: timelike,azimuth,zenitd,airmass
    real(dbl), dimension(:,:), allocatable :: data
    character(len=80) :: msg, fmt

    integer :: i,j,n,m,ncols,nrows,nkeys,nvalid,nstars,ndata,stat

    ! investigating valid fitses
    nrows = 0
    nvalid = 0 ! the first valid fits
    do i = 1, size(fitses)
       if( fitses(i)%status ) nrows = nrows + 1
       if( nvalid == 0 ) nvalid = i
    end do

    nstars = size(cat%coo,1) ! also size(fitses(nvalid)%data,1)
    ncols = size(cols)
    nkeys = size(fitses(nvalid)%headkeys)
    ndata = size(fitses(nvalid)%data)
    n = nstars*ncols
    allocate(labels(n),units(n),headvals(nrows,nkeys),data(nrows,ndata), &
         stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Failed to allocate temporary arrays.'
    end if

    units = [ (fitses(nvalid)%units,i=1,nstars) ]
    do n = 1,nstars
       m = ncols*(n-1)
       write(fmt,'(i0)' ) n
       do i = 1,ncols
          labels(m+i) = trim(cols(i)) // "_" // fmt
       end do
    end do

    n = 0
    do i = 1, size(fitses)
       if( fitses(i)%status ) then
          n = n + 1
          filenames(n) = fitses(i)%filename
          timelike(n) = fitses(i)%gettime(timetype,timestamp,lc_epoch,lc_period)
          headvals(n,:) = fitses(i)%headvals
          azimuth(n) = fitses(i)%azimuth
          zenitd(n) = fitses(i)%zenitd
          airmass(n) = fitses(i)%airmass
          data(n,:) = [ (fitses(i)%data(j,:),j=1,nstars) ]
          ! pack() stores data in unprefered column-order (Fortran way)
       end if
    end do

    call savefits(output,filenames(:n),printtime,printname,printhorizon, &
         printairmass,timetype,timelike(:n),labels,units, &
         fitses(nvalid)%keyunits,fitses(nvalid)%headkeys, &
         fitses(nvalid)%keytype,headvals(:n,:),cat%jd, &
         azimuth(:n),zenitd(:n),airmass(:n),cat%coo(:,1),cat%coo(:,2), &
         data(:n,:),stat)

    status = stat == 0

  end subroutine savelc


  subroutine savefits(output,filenames,printtime,printname,printhorizon, &
       printairmass,timetype,timelike,labels,units,hunits,headkeys,keytype,&
       headvals,epoch,azimuth,zenitd,airmass,xcoo,ycoo,data,status)

    character(len=*), intent(in) :: output,timetype
    character(len=*), dimension(:), intent(in) :: filenames
    logical, intent(in) :: printtime, printname, printhorizon, printairmass
    character(len=*), dimension(:), intent(in) :: headkeys, labels,units, &
         hunits, keytype
    character(len=*), dimension(:,:), intent(in) :: headvals
    real(dbl), dimension(:), intent(in) :: timelike,xcoo,ycoo,azimuth,zenitd,&
         airmass
    real(dbl), dimension(:,:), intent(in) :: data
    real(dbl), intent(in) :: epoch
    integer, intent(out) :: status

    character(len=FLEN_VALUE), dimension(:), allocatable :: ttype, tform, tunit
    character(len=FLEN_VALUE), dimension(size(keytype)) :: hform
    integer :: i,l,n,m,num,nadd,nrows,srows,ntypes,nkeys,frow
    type(fitsfiles) :: fitsfile

    nadd = 0
    if( printtime ) nadd = nadd + 1
    if( printname ) nadd = nadd + 1
    if( printhorizon ) nadd = nadd + 2
    if( printairmass ) nadd = nadd + 1

    nkeys = size(headkeys)
    ntypes = nadd + nkeys + size(labels)
    allocate(ttype(ntypes),tform(ntypes),tunit(ntypes))
    tform = '1D'
    tunit = ''

    num = 0
    if( printname ) then
       num = num + 1
       ttype(1) = FITS_COL_FILENAME
       write(tform(1),"(i0,'A')") max(maxlen(filenames),1)
       tunit(1) = ''
    end if

    if( printtime ) then
       num = num + 1
       ttype(num) = FITS_COL_TIME
       tform(num) = '1D'
       tunit(num) = 'day'   ! jd
       if( timetype == 'PHASE' ) tunit(num) = ''
    end if

    if( nkeys > 0 ) then
       do i = 1, nkeys
          select case(keytype(i))
          case('L')
             hform(i) = '1L'
          case('I')
             hform(i) = '1J'
          case('F')
             hform(i) = '1D'
          case default
             write(hform(i),"(i0,'A')") max(maxlen(headvals(:,i)),1)
          end select
       end do
       n = num + 1
       m = num + nkeys
       num = num + nkeys
       ttype(n:m) = headkeys
       tform(n:m) = hform
       tunit(n:m) = hunits
    end if

    if( printhorizon ) then
       num = num + 2
       ttype(num-1:num) = [ FITS_COL_AZIMUTH, FITS_COL_ZENITD ]
       tform(num-1:num) = '1D'
       tunit(num-1:num) = 'deg'
    end if

    if( printairmass ) then
       num = num + 1
       ttype(num) = FITS_COL_AIRMASS
       tform(num) = '1D'
       tunit(num) = ''
    end if

    ! all other columns
    ttype(num+1:) = labels
    tform(num+1:) = '1D'
    tunit(num+1:) = units

    status = 0
    nrows = size(filenames)
    if( fits_file_exist(output) ) call fits_file_delete(output)
    call fits_create_file(fitsfile,output,status)
    call fits_insert_img(fitsfile,8,0,[0],status)
    call fits_write_key(fitsfile,'HDUNAME',MEXTNAMETS,'time series table',status)
    call fits_write_key(fitsfile,FITS_KEY_CREATOR,FITS_VALUE_CREATOR, &
         FITS_COM_CREATOR,status)
    call fits_write_comment(fitsfile,MUNIPACK_VERSION,status)
    call fits_write_comment(fitsfile,&
         'Description: http://munipack.physics.muni.cz/dataform_tmserie.html',&
         status)

    if( status /= 0  ) goto 666

    ! LC extension
    call fits_insert_btbl(fitsfile,nrows,ttype,tform,tunit,MEXTNAMETS,status)

    ! time info: http://dotastro.org/simpletimeseries/
    call fits_write_key(fitsfile,'TIMETYPE',timetype,'JD,MJD,HJD,phase',status)
    call fits_write_key(fitsfile,'TIMESYS','UTC','',status)
    call fits_write_key(fitsfile,'TIMEREF','GEOCENTER','',status)
    call fits_write_key(fitsfile,'TIMESTMP','MIDPOINT', &
         'time stamp = MIDPOINT,BEGIN,END',status)

    call fits_write_comment(fitsfile,'File names (by input order):',status)
    do i = 1, size(filenames)
       call fits_write_comment(fitsfile,filenames(i),status)
    end do

    ! data vomit
    call fits_get_rowsize(fitsfile,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       num = 0
       if( printname ) then
          num = num + 1
          call fits_write_col(fitsfile,num,frow,filenames(i:l),status)
       end if

       if( printtime ) then
          num = num + 1
          call fits_write_col(fitsfile,num,frow,timelike(i:l),status)
       end if

       ! keys
       do n = 1, nkeys
          select case(keytype(n))
          case('L')
             call tits_write_col_logical(fitsfile,num+n,frow,headvals(i:l,n),status)
          case('I')
             call tits_write_col_integer(fitsfile,num+n,frow,headvals(i:l,n),status)
          case('F')
             call tits_write_col_double(fitsfile,num+n,frow,headvals(i:l,n),status)
          case default
             call fits_write_col(fitsfile,num+n,frow,headvals(i:l,n),status)
          end select
       end do
       num = num + nkeys

       ! horizon
       if( printhorizon ) then
          call fits_write_col(fitsfile,num+1,frow,azimuth(i:l),status)
          call fits_write_col(fitsfile,num+2,frow,zenitd(i:l),status)
          num = num + 2
       end if

       ! airmass
       if( printairmass ) then
          num = num + 1
          call fits_write_col(fitsfile,num,frow,airmass(i:l),status)
       end if

       ! data
       do n = 1,size(labels)
          call fits_write_col(fitsfile,num+n,frow,data(i:l,n),status)
       end do

       if( status /= 0 ) goto 666
    end do
    deallocate(ttype,tform,tunit)

    ! OBJECT extension
    nrows = size(xcoo)
    ntypes = 2
    allocate(ttype(ntypes),tform(ntypes),tunit(ntypes))
    tform = '1D'
    tunit = 'deg'
    ttype(1) = FITS_COL_RA
    ttype(2) = FITS_COL_DEC

    call fits_insert_btbl(fitsfile,nrows,ttype,tform,tunit,MEXTNAMETSC,status)
    call fits_write_key(fitsfile,"REFRAME",'ICRS', &
         'reference celestial coordinate system',status)
    call fits_write_key(fitsfile,"EPOCH",epoch,15,'catalogue reference epoch',status)

    call fits_get_rowsize(fitsfile,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_write_col(fitsfile,1,frow,xcoo(i:l),status)
       call fits_write_col(fitsfile,2,frow,ycoo(i:l),status)
       if( status /= 0 ) goto 666
    end do

    deallocate(ttype,tform,tunit)

666 continue

    call fits_close_file(fitsfile,status)
    call fits_report_error(error_unit,status)

  end subroutine savefits


  integer function maxlen(a)

    ! returns maximum lenght of an element of the passed character array
    ! minimal possible lenght is 1 due FITS table format

    character(len=*), dimension(:) :: a
    integer :: i,n,m

    m = 0
    do i = 1, size(a)
       n = len_trim(a(i))
       if( n > m ) m = n
    end do
    maxlen = m

  end function maxlen

  subroutine tits_write_col_logical(fitsfile,colnum,frow,values,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, frow
    character(len=*), dimension(:), intent(in) :: values
    integer, intent(in out) :: status

    logical, dimension(size(values)) :: vals
    integer :: i

    do i = 1, size(values)
       read(values(i),*) vals(i)
    end do
    call fits_write_col(fitsfile,colnum,frow,vals,status)

  end subroutine tits_write_col_logical

  subroutine tits_write_col_integer(fitsfile,colnum,frow,values,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, frow
    character(len=*), dimension(:), intent(in) :: values
    integer, intent(in out) :: status

    integer, dimension(size(values)) :: vals
    integer :: i

    do i = 1, size(values)
       read(values(i),*) vals(i)
    end do
    call fits_write_col(fitsfile,colnum,frow,vals,status)

  end subroutine tits_write_col_integer

  subroutine tits_write_col_double(fitsfile,colnum,frow,values,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, frow
    character(len=*), dimension(:), intent(in) :: values
    integer, intent(in out) :: status

    real(dbl), dimension(size(values)) :: vals
    integer :: i

    do i = 1, size(values)
       read(values(i),*) vals(i)
    end do
    call fits_write_col(fitsfile,colnum,frow,vals,status)

  end subroutine tits_write_col_double

end module savelcfits
