!
!  fitspphot
!
!  Copyright © 2013 F.Hroch (hroch@physics.muni.cz)
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


module fitspphot

  use titsio
  use iso_fortran_env

  implicit none

  integer, parameter, private :: rp = selected_real_kind(15)

contains

subroutine fpphot(filename,verbose,plog,fkeys)

  use psfengine

  integer, parameter :: DIM = 2

  character(len=*),intent(in) :: filename
  logical, intent(in) :: verbose, plog
  character(len=*),dimension(:),intent(in) :: fkeys

  integer :: status,naxis,bitpix,nrows,xcol,ycol,scol,ccol,hdutype,i,j,n
  integer, parameter :: group = 1, extver = 0, frow = 1, nbegin = 4
  real(rp), parameter :: nullval = 0.0_rp
  integer, dimension(DIM) :: naxes
  logical :: anyf
  real(rp), dimension(:,:), allocatable :: data, rframe
  real(rp), dimension(:), allocatable :: xcens,ycens,sky,sky_err,asky,acts,cts,cts_err,q
  logical, dimension(:), allocatable :: id
  character(len=FLEN_VALUE), dimension(:), allocatable :: ttype, tform, tunit
  real(rp) :: lobad, hibad, phpadu, fwhm, ring
  type(fitsfiles) :: fits

  status = 0
  call fits_open_image(fits,filename,FITS_READWRITE,status)
  if( status /= 0 ) then
     write(error_unit,*) &
          'Error: failed to open an image in the file `',trim(filename),"'."
     return
  end if

  call fits_get_img_param(fits,bitpix,naxis,naxes,status)
  if( status /= 0 ) goto 666

  call fits_read_key(fits,fkeys(1),phpadu,status)
  if( status == FITS_KEYWORD_NOT_FOUND ) then
     write(error_unit,*) trim(filename),": Gain identified by `",trim(fkeys(1)),"' keyword not found."
     goto 666
  end if

  allocate(data(naxes(1),naxes(2)))
  call fits_read_image(fits,group,nullval,data,anyf,status)

  if( status /= 0 ) then
     write(error_unit,*) trim(filename),": Failed to read data."
     goto 666
  end if

  call fits_movnam_hdu(fits,FITS_BINARY_TBL,FINDEXTNAME,extver,status)
  if( status == FITS_BAD_HDU_NUM ) then
     write(error_unit,*) trim(filename),": ", &
          trim(FINDEXTNAME)//' extension not found. Has been stars detected already?'
     goto 666
  end if

  call fits_read_key(fits,FITS_KEY_LOWBAD,lobad,status)
  call fits_read_key(fits,FITS_KEY_HIGHBAD,hibad,status)
  call fits_read_key(fits,FITS_KEY_FWHM,fwhm,status)

  if( status /= 0 ) then
     write(error_unit,*) trim(filename),": Required keywords ",trim(FITS_KEY_LOWBAD),&
          ",",trim(FITS_KEY_HIGHBAD)," or ",trim(FITS_KEY_FWHM)," not found."
     goto 666
  end if

  call fits_movnam_hdu(fits,FITS_BINARY_TBL,APEREXTNAME,extver,status)
  if( status == FITS_BAD_HDU_NUM ) then
     write(error_unit,*) trim(filename),": ",trim(APEREXTNAME) // &
          ' extension not found. An aperture photometry is missing.'
     goto 666
  end if

  call fits_get_num_rows(fits,nrows,status)

  call fits_read_key(fits,trim(FITS_KEY_ANNULUS)//'1',ring,status)
  if( status /= 0 ) then
     write(error_unit,*) trim(filename),": Required keywords ",trim(FITS_KEY_ANNULUS),&
          " not found."
     goto 666
  end if
  ring = 10*fwhm

  allocate(xcens(nrows),ycens(nrows),asky(nrows),acts(nrows))

  call fits_get_colnum(fits,.true.,FITS_COL_X,xcol,status)
  call fits_get_colnum(fits,.true.,FITS_COL_Y,ycol,status)
  call fits_get_colnum(fits,.true.,FITS_COL_SKY,scol,status)
  call fits_get_colnum(fits,.true.,trim(FITS_COL_APCOUNT)//'1',ccol,status)

  call fits_read_col(fits,xcol,frow,nullval,xcens,anyf,status)
  call fits_read_col(fits,ycol,frow,nullval,ycens,anyf,status)
  call fits_read_col(fits,scol,frow,nullval,asky,anyf,status)
  call fits_read_col(fits,ccol,frow,nullval,acts,anyf,status)
  if( status /= 0 ) goto 666

  allocate(cts(nrows),cts_err(nrows),sky(nrows),sky_err(nrows), &
       rframe(size(data,1),size(data,2)))

  call psfmodel(xcens,ycens,asky,acts,data,ring,lobad,hibad,phpadu, &
       verbose,plog,cts,cts_err,sky,sky_err,rframe,status)

  if( status /= 0 ) goto 666

  ! store results to next extension
  call fits_movnam_hdu(fits,FITS_BINARY_TBL,PSFEXTNAME,extver,status)
  if( status == FITS_BAD_HDU_NUM ) then
     status = 0
  else
     ! already presented ? remove it !
     call fits_delete_hdu(fits,hdutype,status)
     if( status /= 0 ) goto 666
  end if

  n = 6
  allocate(ttype(n), tform(n), tunit(n))

  tform = '1D'
  tunit = ''
  ttype(1) = FITS_COL_X
  ttype(2) = FITS_COL_Y
  ttype(3) = FITS_COL_SKY
  ttype(4) = FITS_COL_SKYERR
  ttype(5) = FITS_COL_PSFCOUNT
  ttype(6) = FITS_COL_PSFCOUNTERR

  ! PSF photometry
  call fits_insert_btbl(fits,0,ttype,tform,tunit,PSFEXTNAME,status)
  call fits_update_key(fits,FITS_KEY_FWHM,fwhm,-5,'[pix] standard FWHM of objects',status)

  call fits_write_col(fits,1,frow,xcens,status)
  call fits_write_col(fits,2,frow,ycens,status)
  call fits_write_col(fits,3,frow,sky,status)
  call fits_write_col(fits,4,frow,sky_err,status)
  call fits_write_col(fits,5,frow,cts,status)
  call fits_write_col(fits,6,frow,cts_err,status)
  deallocate(ttype,tform,tunit)

  ! photometry table
  call fits_movnam_hdu(fits,FITS_BINARY_TBL,PHOTOEXTNAME,extver,status)
  if( status == FITS_BAD_HDU_NUM ) then
     status = 0
  else
     ! already presented ? remove it !
     call fits_delete_hdu(fits,hdutype,status)
     if( status /= 0 ) goto 666
  end if

  allocate(ttype(6), tform(6), tunit(6))

  tform = '1D'
  tunit = ''
  ttype(1) = FITS_COL_X
  ttype(2) = FITS_COL_Y
  ttype(3) = FITS_COL_SKY
  ttype(4) = FITS_COL_SKYERR
  ttype(5) = FITS_COL_COUNT
  ttype(6) = FITS_COL_COUNTERR
  call fits_insert_btbl(fits,0,ttype,tform,tunit,PHOTOEXTNAME,status)
  call fits_update_key(fits,FITS_KEY_ORIGHDU,APEREXTNAME,'copied from',status)
  call fits_update_key(fits,FITS_KEY_FWHM,fwhm,-5,'[pix] standard FWHM of objects',status)

  ! select valid records only
  allocate(id(size(xcens)))
  id = .false.
  n = 0
  do i = 1, size(xcens)
     if( cts(i) >= 0 ) then
        n = n + 1
        id(i) = .true.
     end if
  end do

  allocate(q(n))

  q = pack(xcens,id)
  call fits_write_col(fits,1,frow,q,status)

  q = pack(ycens,id)
  call fits_write_col(fits,2,frow,q,status)

  q = pack(sky,id)
  call fits_write_col(fits,3,frow,q,status)

  q = pack(sky_err,id)
  call fits_write_col(fits,4,frow,q,status)

  q = pack(cts,id)
  call fits_write_col(fits,5,frow,q,status)

  q = pack(cts_err,id)
  call fits_write_col(fits,6,frow,q,status)

  call fits_close_file(fits,status)

  ! save residuals
  call fits_create_file(fits,'!r.fits',status)
  call fits_insert_img(fits,-32,naxis,naxes,status)
  call fits_write_image(fits,group,rframe,status)
  call fits_close_file(fits,status)
  call fits_report_error(error_unit,status)

  deallocate(id,q)

  deallocate(ttype,tform,tunit)

666 continue

  if( allocated(data) ) deallocate(data)
  if( allocated(xcens) ) deallocate(xcens,ycens,asky,acts)
  if( allocated(cts) ) deallocate(cts,cts_err,sky,sky_err)

  call fits_close_file(fits,status)
  call fits_report_error(error_unit,status)

end subroutine fpphot

end module fitspphot
