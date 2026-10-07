!
!  Fortran 2008+ interface for (c)FITSIO library
!
!  Purpose of this module is to provide the interface to cFITSIO
!  in native Fortran data types. The long-names provides generic
!  names (selected by actual passed arguments).
!
!  Copyright © 2020-2026 F.Hroch (hroch@physics.muni.cz)
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

module fitsio

  use cfitsio
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  implicit none

  ! basic constants (fitsio.h)
  integer, parameter :: FLEN_FILENAME = 1025 ! max length of a filename
  integer, parameter :: FLEN_KEYWORD = 75    ! max length of a keyword
  integer, parameter :: FLEN_CARD = 81       ! length of a FITS header card
  integer, parameter :: FLEN_VALUE = 71      ! max length of a keyword value string
  integer, parameter :: FLEN_COMMENT = 73    ! max length of a keyword comment string
  integer, parameter :: FLEN_ERRMSG = 81     ! max length of a FITSIO error message
  integer, parameter :: FLEN_STATUS = 31     ! max length of a FITSIO status text


  ! Codes for FITS extension types
  integer, parameter :: &
       FITS_IMAGE_HDU = 0, &
       FITS_ASCII_TBL = 1, &
       FITS_BINARY_TBL = 2, &
       FITS_ANY_HDU = 2

  integer, parameter :: &
       FITS_READONLY = 0, &
       FITS_READWRITE = 1

  ! table labels sensitivity
  integer, parameter :: &
       FITS_CASEINSEN = 0, &
       FITS_CASESEN = 1

  ! error codes
  integer, parameter :: FITS_FILE_NOT_OPENED = 104
  integer, parameter :: FITS_FILE_NOT_CREATED = 105
  integer, parameter :: FITS_READ_ERROR = 108
  integer, parameter :: FITS_BAD_HDU_NUM = 301
  integer, parameter :: FITS_KEYWORD_NOT_FOUND = 202
  integer, parameter :: FITS_COLUMN_NOT_FOUND = 219
  integer, parameter :: FITS_NUMERICAL_OVERFLOW = 412
  integer, parameter :: FITS_MULTIPLE_MATCH = 237


  type :: fitsfiles
     type(c_ptr) :: fptr
     character(len=FLEN_FILENAME) :: filename
  end type fitsfiles


  interface fits_read_key
     module procedure fits_read_key_str, fits_read_key_str_null, &
          fits_read_key_lng, fits_read_key_lng_null, &
          fits_read_key_flt, fits_read_key_flt_null, &
          fits_read_key_dbl, fits_read_key_dbl_null
  end interface fits_read_key

  interface fits_write_key
     module procedure fits_write_key_str, fits_write_key_lng, &
          fits_write_key_flt, fits_write_key_dbl
  end interface fits_write_key

  interface fits_update_key
     module procedure fits_update_key_str, fits_update_key_lng, &
          fits_update_key_flt, fits_update_key_dbl
  end interface fits_update_key

  interface fits_read_col
     module procedure fits_read_col_str, fits_read_col_log, &
          fits_read_col_lng, fits_read_col_flt, fits_read_col_dbl
  end interface fits_read_col

  interface fits_write_col
     module procedure fits_write_col_str, fits_write_col_log, &
          fits_write_col_lng, fits_write_col_flt, fits_write_col_dbl
  end interface fits_write_col

  interface fits_read_img
     module procedure fits_read_img_chr, fits_read_img_int, &
          fits_read_img_flt, fits_read_img_dbl
  end interface fits_read_img

  interface fits_write_img
     module procedure fits_write_img_chr, fits_write_img_int, &
          fits_write_img_flt, fits_write_img_dbl
  end interface fits_write_img

  interface fits_read_2d
     module procedure fits_read_2d_chr, fits_read_2d_int, &
          fits_read_2d_flt, fits_read_2d_dbl
  end interface fits_read_2d

  interface fits_read_3d
     module procedure fits_read_3d_chr, fits_read_3d_int, &
          fits_read_3d_flt, fits_read_3d_dbl
  end interface fits_read_3d

  interface fits_write_2d
     module procedure fits_write_2d_chr, fits_write_2d_int, &
          fits_write_2d_flt, fits_write_2d_dbl
  end interface fits_write_2d

  interface fits_write_3d
     module procedure fits_write_3d_chr, fits_write_3d_int, &
          fits_write_3d_flt, fits_write_3d_dbl
  end interface fits_write_3d

  private :: fchar, cchar

contains


  ! File Access Routines

  subroutine fits_open_file(fitsfile,filename,iomode,status)

    type(fitsfiles), intent(out) :: fitsfile
    character(len=*), intent(in) :: filename
    integer, intent(in) :: iomode
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(len_trim(filename)+1) :: name
    integer(kind=C_INT) :: rwmode, stat

    fitsfile%filename = filename
    call cchar(filename,name)
    stat = status
    rwmode = iomode
    status = ffopen(fitsfile%fptr,name,rwmode,stat)

  end subroutine fits_open_file

  subroutine fits_open_table(fitsfile,filename,iomode,status)

    type(fitsfiles), intent(out) :: fitsfile
    character(len=*), intent(in) :: filename
    integer, intent(in) :: iomode
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(len_trim(filename)+1) :: name
    integer(kind=C_INT) :: rwmode, stat

    fitsfile%filename = filename
    call cchar(filename,name)
    stat = status
    rwmode = iomode
    status = fftopn(fitsfile%fptr,name,rwmode,stat)

  end subroutine fits_open_table

  subroutine fits_open_image(fitsfile,filename,iomode,status)

    type(fitsfiles), intent(out) :: fitsfile
    character(len=*), intent(in) :: filename
    integer, intent(in) :: iomode
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(len_trim(filename)+1) :: name
    integer(kind=C_INT) :: rwmode, stat

    fitsfile%filename = filename
    call cchar(filename,name)
    stat = status
    rwmode = iomode
    status = ffiopn(fitsfile%fptr,name,rwmode,stat)

  end subroutine fits_open_image


  subroutine fits_create_file(fitsfile,filename,status)

    type(fitsfiles), intent(out) :: fitsfile
    character(len=*), intent(in) :: filename
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(len_trim(filename)+1) :: name
    integer(kind=C_INT) :: stat

    fitsfile%filename = filename
    call cchar(filename,name)
    stat = status
    status = ffinit(fitsfile%fptr,name,stat)

  end subroutine fits_create_file


  subroutine fits_close_file(fitsfile,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat

    stat = status
    status = ffclos(fitsfile%fptr,stat)

  end subroutine fits_close_file

  subroutine fits_delete_file(fitsfile,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat

    stat = status
    status = ffdelt(fitsfile%fptr,stat)

  end subroutine fits_delete_file


  ! HDU Access Routines

  subroutine fits_movabs_hdu(fitsfile,hdunum,hdutype,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: hdunum
    integer, intent(out) :: hdutype
    integer, intent(in out) :: status

    integer(kind=C_INT) :: n,m,stat

    stat = status
    n = hdunum
    status = ffmahd(fitsfile%fptr,n,m,stat)
    hdutype = m

  end subroutine fits_movabs_hdu

  subroutine fits_movrel_hdu(fitsfile,nmove,hdutype,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: nmove
    integer, intent(out) :: hdutype
    integer, intent(in out) :: status

    integer(kind=C_INT) :: n,m,stat

    stat = status
    n = nmove
    status = ffmrhd(fitsfile%fptr,n,m,stat)
    hdutype = m

  end subroutine fits_movrel_hdu

  subroutine fits_movnam_hdu(fitsfile,hdutype,extname,extver,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: hdutype, extver
    character(len=*), intent(in) :: extname
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(len_trim(extname)+1) :: name
    integer(kind=C_INT) :: n,m,stat

    call cchar(extname,name)
    stat = status
    n = hdutype
    m = extver
    status = ffmnhd(fitsfile%fptr,n,name,m,stat)

  end subroutine fits_movnam_hdu

  subroutine fits_get_num_hdus(fitsfile,hdunum,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: hdunum
    integer, intent(in out) :: status

    integer(kind=C_INT) :: n,stat

    stat = status
    status = ffthdu(fitsfile%fptr,n,stat)
    hdunum = n

  end subroutine fits_get_num_hdus

  subroutine fits_get_hdu_num(fitsfile,hdunum)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: hdunum

    integer(kind=C_INT) :: n

    hdunum = ffghdn(fitsfile%fptr,n)

  end subroutine fits_get_hdu_num

  subroutine fits_copy_file(infits,outfits,previous,current,following,status)

    type(fitsfiles), intent(in) :: infits, outfits
    integer, intent(in) :: previous,current,following
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat,n,m,k

    stat = status
    n = previous
    m = current
    k = following
    status = ffcpfl(infits%fptr,outfits%fptr,n,m,k,stat)

  end subroutine fits_copy_file

  subroutine fits_copy_hdu(infits,outfits,morekeys,status)

    type(fitsfiles), intent(in) :: infits, outfits
    integer, intent(in) :: morekeys
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat,m

    stat = status
    m = morekeys
    status = ffcopy(infits%fptr,outfits%fptr,m,stat)

  end subroutine fits_copy_hdu

  subroutine fits_copy_header(infits,outfits,status)

    type(fitsfiles), intent(in) :: infits, outfits
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat

    stat = status
    status = ffcphd(infits%fptr,outfits%fptr,stat)

  end subroutine fits_copy_header

  subroutine fits_delete_hdu(fitsfile,hdutype,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: hdutype
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat, type

    stat = status
    status = ffdhdu(fitsfile%fptr,type,stat)
    hdutype = type

  end subroutine fits_delete_hdu

  subroutine fits_insert_img(fitsfile,bitpix,naxis,naxes,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: bitpix, naxis
    integer, dimension(:), intent(in) :: naxes
    integer, intent(in out) :: status

    integer(kind=C_INT) :: bpix, maxis, stat
    integer(kind=C_LONG), dimension(size(naxes)) :: maxes

    stat = status
    bpix = bitpix
    maxis = naxis
    maxes = naxes
    status = ffiimg(fitsfile%fptr,bpix,maxis,maxes,stat)

  end subroutine fits_insert_img

  subroutine fits_insert_btbl(fitsfile,nrows,ttype,tform,tunit,extname,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: nrows
    character(len=*), dimension(:), intent(in) :: ttype,tform,tunit
    character(len=*) :: extname
    integer, intent(in out) :: status

    integer :: i
    character(kind=C_CHAR), dimension(len_trim(extname)+1) :: name
    integer(kind=C_LONG_LONG) :: mrows, pcount
    integer(kind=C_INT) :: tfields, stat
    character(kind=C_CHAR), dimension(len(ttype)+1,size(ttype)), target :: mtype
    character(kind=C_CHAR), dimension(len(tform)+1,size(ttype)), target :: mform
    character(kind=C_CHAR), dimension(len(tunit)+1,size(ttype)), target :: munit
    type(c_ptr), dimension(size(ttype)) :: type_ptr, form_ptr, unit_ptr

    stat = status
    mrows = nrows
    tfields = size(ttype)
    pcount = 0

    do i = 1, tfields
       call cchar(ttype(i),mtype(:,i))
       type_ptr(i) = C_LOC(mtype(:,i))
       call cchar(tform(i),mform(:,i))
       form_ptr(i) = C_LOC(mform(:,i))
       call cchar(tunit(i),munit(:,i))
       unit_ptr(i) = C_LOC(munit(:,i))
    end do

    call cchar(extname,name)

    status = ffibin(fitsfile%fptr,mrows,tfields,type_ptr,form_ptr,unit_ptr, &
         name,pcount,stat)

  end subroutine fits_insert_btbl


  ! Header Keyword Read/Write Routines

  subroutine fits_get_hdrspace(fitsfile,keysexist,morekeys,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: keysexist,morekeys
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat, nk, mk

    stat = status
    status = ffghsp(fitsfile%fptr,nk,mk,stat)
    keysexist = nk
    morekeys = mk

  end subroutine fits_get_hdrspace

  subroutine fits_read_keyword(fitsfile,keyname,val,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    character(len=*), intent(out) :: val, com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k,v
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgkey(fitsfile%fptr,k,v,c,stat)
    call fchar(v,val)
    call fchar(c,com)

  end subroutine fits_read_keyword

  subroutine fits_read_key_str(fitsfile,keyname,val,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    character(len=*), intent(out) :: val, com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k,v
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgkys(fitsfile%fptr,k,v,c,stat)
    call fchar(v,val)
    call fchar(c,com)

  end subroutine fits_read_key_str

  subroutine fits_read_key_str_null(fitsfile,keyname,val,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    character(len=*), intent(out) :: val
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k,v
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgkys(fitsfile%fptr,k,v,c,stat)
    call fchar(v,val)

  end subroutine fits_read_key_str_null


  subroutine fits_read_key_lng(fitsfile,keyname,val,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    integer, intent(out) :: val
    character(len=*), intent(out) :: com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    integer(kind=C_LONG) :: n
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgkyj(fitsfile%fptr,k,n,c,stat)
    call fchar(c,com)
    val = int(n)

  end subroutine fits_read_key_lng

  subroutine fits_read_key_lng_null(fitsfile,keyname,val,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    integer, intent(out) :: val
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    integer(kind=C_LONG) :: n
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgkyj(fitsfile%fptr,k,n,c,stat)
    val = int(n)

  end subroutine fits_read_key_lng_null


  subroutine fits_read_key_flt(fitsfile,keyname,val,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    real(REAL32), intent(out) :: val
    character(len=*), intent(out) :: com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    real(kind=C_FLOAT) :: x
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgkye(fitsfile%fptr,k,x,c,stat)
    call fchar(c,com)
    val = x

  end subroutine fits_read_key_flt

  subroutine fits_read_key_flt_null(fitsfile,keyname,val,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    real(REAL32), intent(out) :: val
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    real(kind=C_FLOAT) :: x
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgkye(fitsfile%fptr,k,x,c,stat)
    val = x

  end subroutine fits_read_key_flt_null

  subroutine fits_read_key_dbl(fitsfile,keyname,val,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    real(REAL64), intent(out) :: val
    character(len=*), intent(out) :: com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    real(kind=C_DOUBLE) :: x
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgkyd(fitsfile%fptr,k,x,c,stat)
    call fchar(c,com)
    val = x

  end subroutine fits_read_key_dbl

  subroutine fits_read_key_dbl_null(fitsfile,keyname,val,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    real(REAL64), intent(out) :: val
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    real(kind=C_DOUBLE) :: x
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgkyd(fitsfile%fptr,k,x,c,stat)
    val = x

  end subroutine fits_read_key_dbl_null

  subroutine fits_write_key_str(fitsfile,keyname,val,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname, val, com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k,v
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    call cchar(val,v)
    call cchar(com,c)
    status = ffpkys(fitsfile%fptr,k,v,c,stat)

  end subroutine fits_write_key_str

  subroutine fits_write_key_lng(fitsfile,keyname,val,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    integer, intent(in) :: val
    character(len=*), intent(in) :: com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    integer(kind=C_LONG_LONG) :: n
    integer(kind=C_INT) :: stat

    stat = status
    n = val
    call cchar(keyname,k)
    call cchar(com,c)
    status = ffpkyj(fitsfile%fptr,k,n,c,stat)

  end subroutine fits_write_key_lng

  subroutine fits_write_key_flt(fitsfile,keyname,val,decimals,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    real(REAL32), intent(in) :: val
    integer, intent(in) :: decimals
    character(len=*), intent(in) :: com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    real(kind=C_FLOAT) :: x
    integer(kind=C_INT) :: stat, d

    x = val
    d = decimals
    stat = status
    call cchar(keyname,k)
    call cchar(com,c)
    status = ffpkye(fitsfile%fptr,k,x,d,c,stat)

  end subroutine fits_write_key_flt

  subroutine fits_write_key_dbl(fitsfile,keyname,val,decimals,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    real(REAL64), intent(in) :: val
    integer, intent(in) :: decimals
    character(len=*), intent(in) :: com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    real(kind=C_DOUBLE) :: x
    integer(kind=C_INT) :: stat, d

    x = val
    d = decimals
    stat = status
    call cchar(keyname,k)
    call cchar(com,c)
    status = ffpkyd(fitsfile%fptr,k,x,d,c,stat)

  end subroutine fits_write_key_dbl

  subroutine fits_update_key_str(fitsfile,keyname,val,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname, val, com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k,v
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    call cchar(val,v)
    call cchar(com,c)
    status = ffukys(fitsfile%fptr,k,v,c,stat)

  end subroutine fits_update_key_str

  subroutine fits_update_key_lng(fitsfile,keyname,val,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    integer, intent(in) :: val
    character(len=*), intent(in) :: com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    integer(kind=C_LONG_LONG) :: n
    integer(kind=C_INT) :: stat

    stat = status
    n = val
    call cchar(keyname,k)
    call cchar(com,c)
    status = ffukyj(fitsfile%fptr,k,n,c,stat)

  end subroutine fits_update_key_lng

  subroutine fits_update_key_flt(fitsfile,keyname,val,decimals,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    real(REAL32), intent(in) :: val
    integer, intent(in) :: decimals
    character(len=*), intent(in) :: com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    real(kind=C_FLOAT) :: x
    integer(kind=C_INT) :: stat, d

    x = val
    d = decimals
    stat = status
    call cchar(keyname,k)
    call cchar(com,c)
    status = ffukye(fitsfile%fptr,k,x,d,c,stat)

  end subroutine fits_update_key_flt

  subroutine fits_update_card(fitsfile,keyname,card,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname, card
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    call cchar(card,c)
    status = ffucrd(fitsfile%fptr,k,c,stat)

  end subroutine fits_update_card

  subroutine fits_update_key_dbl(fitsfile,keyname,val,decimals,com,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    real(REAL64), intent(in) :: val
    integer, intent(in) :: decimals
    character(len=*), intent(in) :: com
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: c,k
    real(kind=C_DOUBLE) :: x
    integer(kind=C_INT) :: stat, d

    x = val
    d = decimals
    stat = status
    call cchar(keyname,k)
    call cchar(com,c)
    status = ffukyd(fitsfile%fptr,k,x,d,c,stat)

  end subroutine fits_update_key_dbl

  subroutine fits_read_key_unit(fitsfile,keyname,unit,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    character(len=*), intent(out) :: unit
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: k,u
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffgunt(fitsfile%fptr,k,u,stat)
    call fchar(u,unit)

  end subroutine fits_read_key_unit

  subroutine fits_delete_key(fitsfile,keyname,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: keyname
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: k
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(keyname,k)
    status = ffdkey(fitsfile%fptr,k,stat)

  end subroutine fits_delete_key

  subroutine fits_write_comment(fitsfile,comment,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: comment
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: com
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(comment,com)
    status = ffpcom(fitsfile%fptr,com,stat)

  end subroutine fits_write_comment

  subroutine fits_write_history(fitsfile,history,status)

    type(fitsfiles), intent(in) :: fitsfile
    character(len=*), intent(in) :: history
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: hist
    integer(kind=C_INT) :: stat

    stat = status
    call cchar(history,hist)
    status = ffpcom(fitsfile%fptr,hist,stat)

  end subroutine fits_write_history

  subroutine fits_read_record(fitsfile,keynum,card,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: keynum
    character(len=*), intent(out) :: card
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(FLEN_CARD) :: line
    integer(kind=C_INT) :: stat, n

    stat = status
    n = keynum
    status = ffgrec(fitsfile%fptr,n,line,stat)
    call fchar(line,card)

  end subroutine fits_read_record

  subroutine fits_delete_record(fitsfile,keynum,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: keynum
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat, n

    stat = status
    n = keynum
    status = ffdrec(fitsfile%fptr,n,stat)

  end subroutine fits_delete_record

  subroutine fits_hdr2str(fitsfile,nocomments,exclist,nkeys,header,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: nocomments
    character(len=*), dimension(:), intent(in) :: exclist
    integer, intent(out) :: nkeys
    character(len=:), allocatable, intent(out) :: header
    integer, intent(in out) :: status

    integer, parameter :: LEN_CARD = 80
    character(kind=C_CHAR), dimension(len(exclist)+1,size(exclist)), target :: excl

    type(c_ptr) :: phdr = C_NULL_PTR
    character, dimension(:), pointer :: hdr
    type(c_ptr), dimension(size(exclist)) :: ptrs_excl
    integer(kind=C_INT) :: stat, nc, nexc, nk
    integer :: i

    stat = status
    nc = nocomments
    nexc = size(exclist)
    do i = 1, size(exclist)
       call cchar(exclist(i),excl(:,i))
       ptrs_excl(i) = C_LOC(excl(:,i))
    end do
    status = ffhdr2str(fitsfile%fptr,nc,ptrs_excl,nexc,phdr,nk,stat)
    if( status == 0 ) then
       nkeys = nk
       call c_f_pointer(phdr,hdr,[nkeys*LEN_CARD+1])
       ! this allocation is for a non-local variable, the calling
       ! routine should deallocate `header' by hand!
       allocate(character(len=nkeys*LEN_CARD+1) :: header)
       call fchar(hdr,header)
    else
       nkeys = 0
       allocate(character(len=1) :: header)
       header = ''
    end if

    stat = status
    status = fffree(phdr,stat)

  end subroutine fits_hdr2str


  subroutine fits_parse_template(templt,card,keytype,status)

    character(len=*), intent(in) :: templt
    character(len=*), intent(out) :: card
    integer, intent(out) :: keytype
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat, n
    character(kind=C_CHAR), dimension(FLEN_CARD) :: line, tmp

    stat = status
    call cchar(templt,tmp)
    status = ffgthd(tmp,line,n,stat)
    call fchar(line,card)
    keytype = n

  end subroutine fits_parse_template

  subroutine fits_get_keyname(card,keyname,keylenght,status)

    character(len=*), intent(in) :: card
    character(len=*), intent(out) :: keyname
    integer, intent(out) :: keylenght
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat, n
    character(kind=C_CHAR), dimension(FLEN_CARD) :: line, tmp

    stat = status
    call cchar(card,tmp)
    status = ffgknm(tmp,line,n,stat)
    call fchar(line,keyname)
    keylenght = n

  end subroutine fits_get_keyname


  ! Primary Array or IMAGE Extension I/O Routines

  subroutine fits_get_img_type(fitsfile,bitpix,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: bitpix
    integer, intent(in out) :: status

    integer(kind=C_INT) :: bpix, stat

    stat = status
    status = ffgidt(fitsfile%fptr,bpix,stat)
    status = stat
    bitpix = bpix

  end subroutine fits_get_img_type

  subroutine fits_get_img_dim(fitsfile,naxis,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: naxis
    integer, intent(in out) :: status

    integer(kind=C_INT) :: maxis, stat

    stat = status
    status = ffgidm(fitsfile%fptr,maxis,stat)
    naxis = maxis

  end subroutine fits_get_img_dim

  subroutine fits_get_img_size(fitsfile,naxis,naxes,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: naxis
    integer, dimension(:), intent(out) :: naxes
    integer, intent(in out) :: status

    integer(kind=C_INT) :: maxdim, stat
    integer(kind=C_LONG), dimension(naxis) :: maxes

    stat = status
    maxdim = naxis
    status = ffgisz(fitsfile%fptr,maxdim,maxes,stat)
    naxes = int(maxes)

  end subroutine fits_get_img_size

  subroutine fits_get_img_param(fitsfile,bitpix,naxis,naxes,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: bitpix, naxis
    integer, dimension(:), intent(out) :: naxes
    integer, intent(in out) :: status

    integer(kind=C_INT) :: maxdim, bpix, maxis, stat
    integer(kind=C_LONG), dimension(size(naxes)) :: maxes

    stat = status
    maxdim = size(maxes)
    status = ffgipr(fitsfile%fptr,maxdim,bpix,maxis,maxes,stat)
    bitpix = bpix
    naxis = maxis
    naxes = int(maxes)

  end subroutine fits_get_img_param

  subroutine fits_read_img_chr(fitsfile,group,firstelem,nelements,nullval, &
       array,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, firstelem, nelements
    character, intent(in) :: nullval
    character, dimension(:), intent(out) :: array
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: fpixel, npixels
    integer(kind=C_INT) :: stat, anynull
    character(kind=C_CHAR) :: nval
    character(kind=C_CHAR), dimension(size(array)) :: col

    fpixel = firstelem
    npixels = nelements
    nval = nullval
    grp = group
    stat = status
    status = ffgpvb(fitsfile%fptr,grp,fpixel,npixels,nval,col,anynull,stat)
    array = col
    anyf = anynull == 1

  end subroutine fits_read_img_chr

  subroutine fits_read_img_int(fitsfile,group,firstelem,nelements,nullval, &
       array,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, firstelem, nelements
    integer, intent(in) :: nullval
    integer, dimension(:), intent(out) :: array
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: fpixel, npixels
    integer(kind=C_INT) :: stat, anynull
    integer(kind=C_INT) :: nval
    integer(kind=C_INT), dimension(size(array)) :: col

    fpixel = firstelem
    npixels = nelements
    nval = nullval
    grp = group
    stat = status
    status = ffgpvk(fitsfile%fptr,grp,fpixel,npixels,nval,col,anynull,stat)
    array = col
    anyf = anynull == 1

  end subroutine fits_read_img_int

  subroutine fits_read_img_flt(fitsfile,group,firstelem,nelements,nullval, &
       array,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, firstelem, nelements
    real(REAL32), intent(in) :: nullval
    real(REAL32), dimension(:), intent(out) :: array
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: fpixel, npixels
    integer(kind=C_INT) :: stat, anynull
    real(kind=C_FLOAT) :: nval
    real(kind=C_FLOAT), dimension(size(array)) :: col

    fpixel = firstelem
    npixels = nelements
    nval = nullval
    grp = group
    stat = status
    status = ffgpve(fitsfile%fptr,grp,fpixel,npixels,nval,col,anynull,stat)
    array = col
    anyf = anynull == 1

  end subroutine fits_read_img_flt

  subroutine fits_read_img_dbl(fitsfile,group,firstelem,nelements,nullval, &
       array,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, firstelem, nelements
    real(REAL64), intent(in) :: nullval
    real(REAL64), dimension(:), intent(out) :: array
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: fpixel, npixels
    integer(kind=C_INT) :: stat, anynull
    real(kind=C_DOUBLE) :: nval
    real(kind=C_DOUBLE), dimension(size(array)) :: col

    fpixel = firstelem
    npixels = nelements
    nval = nullval
    grp = group
    stat = status
    status = ffgpvd(fitsfile%fptr,grp,fpixel,npixels,nval,col,anynull,stat)
    array = col
    anyf = anynull == 1

  end subroutine fits_read_img_dbl


  subroutine fits_write_img_chr(fitsfile,group,firstelem,nelements,array,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, firstelem, nelements
    character, dimension(:), intent(in) :: array
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: fpixel, npixels
    integer(kind=C_INT) :: stat
    character(kind=C_CHAR), dimension(size(array)) :: col

    fpixel = firstelem
    npixels = nelements
    grp = group
    stat = status
    col = array
    status = ffpprb(fitsfile%fptr,grp,fpixel,npixels,col,stat)

  end subroutine fits_write_img_chr

  subroutine fits_write_img_int(fitsfile,group,firstelem,nelements,array,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, firstelem, nelements
    integer, dimension(:), intent(in) :: array
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: fpixel, npixels
    integer(kind=C_INT) :: stat
    integer(kind=C_INT), dimension(size(array)) :: col

    fpixel = firstelem
    npixels = nelements
    grp = group
    stat = status
    col = array
    status = ffpprk(fitsfile%fptr,grp,fpixel,npixels,col,stat)

  end subroutine fits_write_img_int

  subroutine fits_write_img_flt(fitsfile,group,firstelem,nelements,array,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, firstelem, nelements
    real(REAL32), dimension(:), intent(in) :: array
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: fpixel, npixels
    integer(kind=C_INT) :: stat
    real(kind=C_FLOAT), dimension(size(array)) :: col

    fpixel = firstelem
    npixels = nelements
    grp = group
    stat = status
    col = array
    status = ffppre(fitsfile%fptr,grp,fpixel,npixels,col,stat)

  end subroutine fits_write_img_flt

  subroutine fits_write_img_dbl(fitsfile,group,firstelem,nelements,array,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, firstelem, nelements
    real(REAL64), dimension(:), intent(in) :: array
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: fpixel, npixels
    integer(kind=C_INT) :: stat
    real(kind=C_DOUBLE), dimension(size(array)) :: col

    fpixel = firstelem
    npixels = nelements
    grp = group
    stat = status
    col = array
    status = ffpprd(fitsfile%fptr,grp,fpixel,npixels,col,stat)

  end subroutine fits_write_img_dbl


  subroutine fits_read_2d_chr(fitsfile,group,nullval,dim1,naxis1,naxis2, &
       image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, dim1,naxis1,naxis2
    character, intent(in) :: nullval
    character, dimension(:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,n1,n2
    integer(kind=C_INT) :: stat, anynull
    character(kind=C_CHAR) :: nval
    character(kind=C_CHAR), dimension(size(image)) :: array

    d1 = dim1
    n1 = naxis1
    n2 = naxis2
    nval = nullval
    grp = group
    stat = status
    status = ffg2db(fitsfile%fptr,grp,nval,d1,n1,n2,array,anynull,stat)
    image = reshape(array,[size(image,1),size(image,2)])
    anyf = anynull == 1

  end subroutine fits_read_2d_chr

  subroutine fits_read_2d_int(fitsfile,group,nullval,dim1,naxis1,naxis2, &
       image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, dim1,naxis1,naxis2
    integer, intent(in) :: nullval
    integer, dimension(:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,n1,n2
    integer(kind=C_INT) :: stat, anynull
    integer(kind=C_INT) :: nval
    integer(kind=C_INT), dimension(size(image)) :: array

    d1 = dim1
    n1 = naxis1
    n2 = naxis2
    nval = nullval
    grp = group
    stat = status
    status = ffg2dk(fitsfile%fptr,grp,nval,d1,n1,n2,array,anynull,stat)
    image = reshape(array,[size(image,1),size(image,2)])
    anyf = anynull == 1

  end subroutine fits_read_2d_int

  subroutine fits_read_2d_flt(fitsfile,group,nullval,dim1,naxis1,naxis2, &
       image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, dim1,naxis1,naxis2
    real(REAL32), intent(in) :: nullval
    real(REAL32), dimension(:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,n1,n2
    integer(kind=C_INT) :: stat, anynull
    real(kind=C_FLOAT) :: nval
    real(kind=C_FLOAT), dimension(size(image)) :: array

    d1 = dim1
    n1 = naxis1
    n2 = naxis2
    nval = nullval
    grp = group
    stat = status
    status = ffg2de(fitsfile%fptr,grp,nval,d1,n1,n2,array,anynull,stat)
    image = reshape(array,[size(image,1),size(image,2)])
    anyf = anynull == 1

  end subroutine fits_read_2d_flt

  subroutine fits_read_2d_dbl(fitsfile,group,nullval,dim1,naxis1,naxis2, &
       image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, dim1,naxis1,naxis2
    real(REAL64), intent(in) :: nullval
    real(REAL64), dimension(:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,n1,n2
    integer(kind=C_INT) :: stat, anynull
    real(kind=C_DOUBLE) :: nval
    real(kind=C_DOUBLE), dimension(size(image)) :: array

    d1 = dim1
    n1 = naxis1
    n2 = naxis2
    nval = nullval
    grp = group
    stat = status
    status = ffg2dd(fitsfile%fptr,grp,nval,d1,n1,n2,array,anynull,stat)
    image = reshape(array,[size(image,1),size(image,2)])
    anyf = anynull == 1

  end subroutine fits_read_2d_dbl

  subroutine fits_write_2d_chr(fitsfile,group,dim1,naxis1,naxis2,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group,dim1,naxis1,naxis2
    character, dimension(:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,n1,n2
    integer(kind=C_INT) :: stat
    character(kind=C_CHAR), dimension(size(image)) :: array

    d1 = dim1
    n1 = naxis1
    n2 = naxis2
    grp = group
    stat = status
    array = pack(image,.true.)
    status = ffp2db(fitsfile%fptr,grp,d1,n1,n2,array,stat)

  end subroutine fits_write_2d_chr

  subroutine fits_write_2d_int(fitsfile,group,dim1,naxis1,naxis2,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group,dim1,naxis1,naxis2
    integer, dimension(:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,n1,n2
    integer(kind=C_INT) :: stat
    integer(kind=C_INT), dimension(size(image)) :: array

    d1 = dim1
    n1 = naxis1
    n2 = naxis2
    grp = group
    stat = status
    array = pack(image,.true.)
    status = ffp2dk(fitsfile%fptr,grp,d1,n1,n2,array,stat)

  end subroutine fits_write_2d_int

  subroutine fits_write_2d_flt(fitsfile,group,dim1,naxis1,naxis2,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group,dim1,naxis1,naxis2
    real(REAL32), dimension(:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,n1,n2
    integer(kind=C_INT) :: stat
    real(kind=C_FLOAT), dimension(size(image)) :: array

    d1 = dim1
    n1 = naxis1
    n2 = naxis2
    grp = group
    stat = status
    array = pack(image,.true.)
    status = ffp2de(fitsfile%fptr,grp,d1,n1,n2,array,stat)

  end subroutine fits_write_2d_flt

  subroutine fits_write_2d_dbl(fitsfile,group,dim1,naxis1,naxis2,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group,dim1,naxis1,naxis2
    real(REAL64), dimension(:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,n1,n2
    integer(kind=C_INT) :: stat
    real(kind=C_DOUBLE), dimension(size(image)) :: array

    d1 = dim1
    n1 = naxis1
    n2 = naxis2
    grp = group
    stat = status
    array = pack(image,.true.)
    status = ffp2dd(fitsfile%fptr,grp,d1,n1,n2,array,stat)

  end subroutine fits_write_2d_dbl

  subroutine fits_read_3d_chr(fitsfile,group,nullval,dim1,dim2,naxis1,naxis2,naxis3,&
       image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, dim1,dim2,naxis1,naxis2,naxis3
    character, intent(in) :: nullval
    character, dimension(:,:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,d2,n1,n2,n3
    integer(kind=C_INT) :: stat, anynull
    character(kind=C_CHAR) :: nval
    character(kind=C_CHAR), dimension(size(image)) :: array

    d1 = dim1
    d2 = dim2
    n1 = naxis1
    n2 = naxis2
    n3 = naxis3
    nval = nullval
    grp = group
    stat = status
    status = ffg3db(fitsfile%fptr,grp,nval,d1,d2,n1,n2,n3,array,anynull,stat)
    image = reshape(array,[size(image,1),size(image,2),size(image,3)])
    anyf = anynull == 1

  end subroutine fits_read_3d_chr

  subroutine fits_read_3d_int(fitsfile,group,nullval,dim1,dim2,naxis1,naxis2,naxis3,&
       image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, dim1,dim2,naxis1,naxis2,naxis3
    integer, intent(in) :: nullval
    integer, dimension(:,:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,d2,n1,n2,n3
    integer(kind=C_INT) :: stat, anynull
    integer(kind=C_INT) :: nval
    integer(kind=C_INT), dimension(size(image)) :: array

    d1 = dim1
    d2 = dim2
    n1 = naxis1
    n2 = naxis2
    n3 = naxis3
    nval = nullval
    grp = group
    stat = status
    status = ffg3dk(fitsfile%fptr,grp,nval,d1,d2,n1,n2,n3,array,anynull,stat)
    image = reshape(array,[size(image,1),size(image,2),size(image,3)])
    anyf = anynull == 1

  end subroutine fits_read_3d_int

  subroutine fits_read_3d_flt(fitsfile,group,nullval,dim1,dim2,naxis1,naxis2,naxis3,&
       image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, dim1,dim2,naxis1,naxis2,naxis3
    real(REAL32), intent(in) :: nullval
    real(REAL32), dimension(:,:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,d2,n1,n2,n3
    integer(kind=C_INT) :: stat, anynull
    real(kind=C_FLOAT) :: nval
    real(kind=C_FLOAT), dimension(size(image)) :: array

    d1 = dim1
    d2 = dim2
    n1 = naxis1
    n2 = naxis2
    n3 = naxis3
    nval = nullval
    grp = group
    stat = status
    status = ffg3de(fitsfile%fptr,grp,nval,d1,d2,n1,n2,n3,array,anynull,stat)
    image = reshape(array,[size(image,1),size(image,2),size(image,3)])
    anyf = anynull == 1

  end subroutine fits_read_3d_flt

  subroutine fits_read_3d_dbl(fitsfile,group,nullval,dim1,dim2,naxis1,naxis2,naxis3,&
       image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group, dim1,dim2,naxis1,naxis2,naxis3
    real(REAL64), intent(in) :: nullval
    real(REAL64), dimension(:,:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,d2,n1,n2,n3
    integer(kind=C_INT) :: stat, anynull
    real(kind=C_DOUBLE) :: nval
    real(kind=C_DOUBLE), dimension(size(image)) :: array

    d1 = dim1
    d2 = dim2
    n1 = naxis1
    n2 = naxis2
    n3 = naxis3
    nval = nullval
    grp = group
    stat = status
    status = ffg3dd(fitsfile%fptr,grp,nval,d1,d2,n1,n2,n3,array,anynull,stat)
    image = reshape(array,[size(image,1),size(image,2),size(image,3)])
    anyf = anynull == 1

  end subroutine fits_read_3d_dbl

  subroutine fits_write_3d_chr(fitsfile,group,dim1,dim2,naxis1,naxis2,naxis3, &
       image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group,dim1,dim2,naxis1,naxis2,naxis3
    character, dimension(:,:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,d2,n1,n2,n3
    integer(kind=C_INT) :: stat
    character(kind=C_CHAR), dimension(size(image)) :: array

    d1 = dim1
    d2 = dim2
    n1 = naxis1
    n2 = naxis2
    n3 = naxis3
    grp = group
    stat = status
    array = pack(image,.true.)
    status = ffp3db(fitsfile%fptr,grp,d1,d2,n1,n2,n3,array,stat)

  end subroutine fits_write_3d_chr

  subroutine fits_write_3d_int(fitsfile,group,dim1,dim2,naxis1,naxis2,naxis3, &
       image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group,dim1,dim2,naxis1,naxis2,naxis3
    integer, dimension(:,:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,d2,n1,n2,n3
    integer(kind=C_INT) :: stat
    integer(kind=C_INT), dimension(size(image)) :: array

    d1 = dim1
    d2 = dim2
    n1 = naxis1
    n2 = naxis2
    n3 = naxis3
    grp = group
    stat = status
    array = pack(image,.true.)
    status = ffp3dk(fitsfile%fptr,grp,d1,d2,n1,n2,n3,array,stat)

  end subroutine fits_write_3d_int

  subroutine fits_write_3d_flt(fitsfile,group,dim1,dim2,naxis1,naxis2,naxis3, &
       image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group,dim1,dim2,naxis1,naxis2,naxis3
    real(REAL32), dimension(:,:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,d2,n1,n2,n3
    integer(kind=C_INT) :: stat
    real(kind=C_FLOAT), dimension(size(image)) :: array

    d1 = dim1
    d2 = dim2
    n1 = naxis1
    n2 = naxis2
    n3 = naxis3
    grp = group
    stat = status
    array = pack(image,.true.)
    status = ffp3de(fitsfile%fptr,grp,d1,d2,n1,n2,n3,array,stat)

  end subroutine fits_write_3d_flt

  subroutine fits_write_3d_dbl(fitsfile,group,dim1,dim2,naxis1,naxis2,naxis3, &
       image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group,dim1,dim2,naxis1,naxis2,naxis3
    real(REAL64), dimension(:,:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: grp
    integer(kind=C_LONG_LONG) :: d1,d2,n1,n2,n3
    integer(kind=C_INT) :: stat
    real(kind=C_DOUBLE), dimension(size(image)) :: array

    d1 = dim1
    d2 = dim2
    n1 = naxis1
    n2 = naxis2
    n3 = naxis3
    grp = group
    stat = status
    array = pack(image,.true.)
    status = ffp3dd(fitsfile%fptr,grp,d1,d2,n1,n2,n3,array,stat)

  end subroutine fits_write_3d_dbl



  ! ASCII and Binary Table Routines

  subroutine fits_get_num_rows(fitsfile,nrows,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: nrows
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: mrows
    integer(kind=C_INT) :: stat

    stat = status
    status = ffgnrw(fitsfile%fptr,mrows,stat)
    nrows = int(mrows)

  end subroutine fits_get_num_rows

  subroutine fits_get_num_cols(fitsfile,ncols,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: ncols
    integer, intent(in out) :: status

    integer(kind=C_INT) :: mcols, stat

    stat = status
    status = ffgncl(fitsfile%fptr,mcols,stat)
    ncols = mcols

  end subroutine fits_get_num_cols

  subroutine fits_get_rowsize(fitsfile,nrows,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(out) :: nrows
    integer, intent(in out) :: status

    integer(kind=C_LONG) :: mrows
    integer(kind=C_INT) :: stat

    stat = status
    status = ffgrsz(fitsfile%fptr,mrows,stat)
    nrows = int(mrows)

  end subroutine fits_get_rowsize

  subroutine fits_get_colnum(fitsfile,casesens,templt,colnum,status)

    type(fitsfiles), intent(in) :: fitsfile
    logical, intent(in) :: casesens
    character(len=*), intent(in) :: templt
    integer, intent(out) :: colnum
    integer, intent(in out) :: status

    integer(kind=C_INT) :: mcase, mcolnum, stat
    character(kind=C_CHAR), dimension(len(templt)+1) :: templ

    stat = status
    if( casesens ) then
       mcase = FITS_CASESEN
    else
       mcase = FITS_CASEINSEN
    end if
    call cchar(templt,templ)

    status = ffgcno(fitsfile%fptr,mcase,templ,mcolnum,stat)
    colnum = mcolnum

  end subroutine fits_get_colnum


  subroutine fits_get_colname(fitsfile,casesens,templt,colname,colnum,status)

    type(fitsfiles), intent(in) :: fitsfile
    logical, intent(in) :: casesens
    character(len=*), intent(in) :: templt
    character(len=*), intent(out) :: colname
    integer, intent(out) :: colnum
    integer, intent(in out) :: status

    integer(kind=C_INT) :: mcase, mcolnum, stat
    character(kind=C_CHAR), dimension(len(templt)+1) :: templ
    character(kind=C_CHAR), dimension(len(colname)+1) :: name

    stat = status
    if( casesens ) then
       mcase = FITS_CASESEN
    else
       mcase = FITS_CASEINSEN
    end if
    call cchar(templt,templ)

    status = ffgcnn(fitsfile%fptr,mcase,templ,name,mcolnum,stat)
    call fchar(name,colname)
    colnum = mcolnum

  end subroutine fits_get_colname

  subroutine fits_insert_col(fitsfile,colnum,ttype,tform,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum
    character(len=*), intent(in) :: ttype,tform
    integer, intent(in out) :: status

    integer(kind=C_INT) :: icol, stat
    character(kind=C_CHAR), dimension(len(ttype)+1) :: mtype
    character(kind=C_CHAR), dimension(len(tform)+1) :: mform

    stat = status
    icol = colnum

    call cchar(ttype,mtype)
    call cchar(tform,mform)

    status = fficol(fitsfile%fptr,icol,mtype,mform,stat)

  end subroutine fits_insert_col

  subroutine fits_insert_cols(fitsfile,colnum,ncols,ttype,tform,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, ncols
    character(len=*), dimension(:), intent(in) :: ttype,tform
    integer, intent(in out) :: status

    integer :: i
    integer(kind=C_INT) :: icol, tfields, stat
    character(kind=C_CHAR), dimension(len(ttype)+1,size(ttype)), target :: mtype
    character(kind=C_CHAR), dimension(len(tform)+1,size(tform)), target :: mform
    type(c_ptr), dimension(size(ttype)) :: type_ptr
    type(c_ptr), dimension(size(tform)) :: form_ptr

    stat = status
    icol = colnum
    tfields = ncols

    do i = 1, ncols
       call cchar(ttype(i),mtype(:,i))
       type_ptr(i) = C_LOC(mtype(:,i))
       call cchar(tform(i),mform(:,i))
       form_ptr(i) = C_LOC(mform(:,i))
    end do

    status = fficls(fitsfile%fptr,icol,tfields,type_ptr,form_ptr,stat)

  end subroutine fits_insert_cols

  subroutine fits_delete_col(fitsfile,colnum,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum
    integer, intent(in out) :: status

    integer(kind=C_INT) :: icol, stat

    stat = status
    icol = colnum

    status = ffdcol(fitsfile%fptr,icol,stat)

  end subroutine fits_delete_col

  subroutine fits_read_col_str(fitsfile,colnum,firstrow,nullstr,array,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    character(len=*), intent(in) :: nullstr
    character(len=*), dimension(:), intent(out) :: array
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer :: i
    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer(kind=C_INT) :: icol, anynull,  stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    character(kind=C_CHAR), dimension(len(nullstr)+1) :: nulls
    character(kind=C_CHAR), dimension(len(array)+1,size(array)), target :: col
    type(c_ptr), dimension(size(array)) :: ptrs

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)

    call cchar(nullstr,nulls)
    do i = 1, size(array)
       ptrs(i) = C_LOC(col(:,i))
    end do
    status = ffgcvs(fitsfile%fptr,icol,frow,felement,nelements,nulls, &
         ptrs,anynull,stat)

    if( status == 0 ) then
       do i = 1, size(array)
          call fchar(col(:,i),array(i))
       end do
    end if
    anyf = anynull == 1

  end subroutine fits_read_col_str

  subroutine fits_read_col_log(fitsfile,colnum,firstrow,nullval,array,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    logical, intent(in) :: nullval
    logical, dimension(:), intent(out) :: array
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer(kind=C_INT) :: icol, anynull,  stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    character(kind=C_CHAR) :: nulls
    character(kind=C_CHAR), dimension(size(array)) :: bools

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)
    if( nullval ) then
       nulls = char(1)
    else
       nulls = char(0)
    end if
    status = ffgcvl(fitsfile%fptr,icol,frow,felement,nelements,nulls, &
         bools,anynull,stat)
    array = ichar(bools) == 1
    anyf = anynull == 1

  end subroutine fits_read_col_log

  subroutine fits_read_col_lng(fitsfile,colnum,firstrow,nullval,array,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    integer, intent(in) :: nullval
    integer, dimension(:), intent(out) :: array
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer(kind=C_INT) :: icol, anynull,  stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    integer(kind=C_LONG) :: nulls
    integer(kind=C_LONG), dimension(size(array)) :: col

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)
    nulls = nullval
    status = ffgcvj(fitsfile%fptr,icol,frow,felement,nelements,nulls, &
         col,anynull,stat)
    array = int(col)
    anyf = anynull == 1

  end subroutine fits_read_col_lng

  subroutine fits_read_col_flt(fitsfile,colnum,firstrow,nullval,array,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    real(REAL32), intent(in) :: nullval
    real(REAL32), dimension(:), intent(out) :: array
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer(kind=C_INT) :: icol, anynull,  stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    real(kind=C_FLOAT) :: nulls
    real(kind=C_FLOAT), dimension(size(array)) :: col

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)
    nulls = nullval
    status = ffgcve(fitsfile%fptr,icol,frow,felement,nelements,nulls, &
         col,anynull,stat)
    array = col
    anyf = anynull == 1

  end subroutine fits_read_col_flt

  subroutine fits_read_col_dbl(fitsfile,colnum,firstrow,nullval,array,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    real(kind=REAL64), intent(in) :: nullval
    real(kind=REAL64), dimension(:), intent(out) :: array
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer(kind=C_INT) :: icol, anynull,  stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    real(kind=C_DOUBLE) :: nulls
    real(kind=C_DOUBLE), dimension(size(array)) :: col

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)
    nulls = nullval
    status = ffgcvd(fitsfile%fptr,icol,frow,felement,nelements,nulls, &
         col,anynull,stat)
    array = col
    anyf = anynull == 1

  end subroutine fits_read_col_dbl


  subroutine fits_write_col_str(fitsfile,colnum,firstrow,array,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    character(len=*), dimension(:), intent(in) :: array
    integer, intent(in out) :: status

    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer :: i
    integer(kind=C_INT) :: icol, stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    character(kind=C_CHAR), dimension(len(array)+1,size(array)), target :: strings
    type(c_ptr), dimension(size(array)) :: ptrs

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)

    do i = 1, size(array)
       call cchar(array(i),strings(:,i))
       ptrs(i) = C_LOC(strings(:,i))
    end do

    status = ffpcls(fitsfile%fptr,icol,frow,felement,nelements,ptrs,stat)

  end subroutine fits_write_col_str

  subroutine fits_write_col_log(fitsfile,colnum,firstrow,array,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    logical, dimension(:), intent(in) :: array
    integer, intent(in out) :: status

    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer(kind=C_INT) :: icol, stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    character(kind=C_CHAR), dimension(size(array)) :: col

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)
    where( array )
       col = char(1)
    elsewhere
       col = char(0)
    end where

    status = ffpcll(fitsfile%fptr,icol,frow,felement,nelements,col,stat)

  end subroutine fits_write_col_log

  subroutine fits_write_col_lng(fitsfile,colnum,firstrow,array,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    integer, dimension(:), intent(in) :: array
    integer, intent(in out) :: status

    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer(kind=C_INT) :: icol, stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    integer(kind=C_LONG), dimension(size(array)) :: col

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)
    col = array

    status = ffpclj(fitsfile%fptr,icol,frow,felement,nelements,col,stat)

  end subroutine fits_write_col_lng

  subroutine fits_write_col_flt(fitsfile,colnum,firstrow,array,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    real(kind=REAL32), dimension(:), intent(in) :: array
    integer, intent(in out) :: status

    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer(kind=C_INT) :: icol, stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    real(kind=C_FLOAT), dimension(size(array)) :: col

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)
    col = array

    status = ffpcle(fitsfile%fptr,icol,frow,felement,nelements,col,stat)

  end subroutine fits_write_col_flt

  subroutine fits_write_col_dbl(fitsfile,colnum,firstrow,array,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: colnum, firstrow
    real(kind=REAL64), dimension(:), intent(in) :: array
    integer, intent(in out) :: status

    integer(kind=C_LONG_LONG), parameter :: felement = 1
    integer(kind=C_INT) :: icol, stat
    integer(kind=C_LONG_LONG) :: frow, nelements
    real(kind=C_DOUBLE), dimension(size(array)) :: col

    stat = status
    icol = colnum
    frow = firstrow
    nelements = size(array)
    col = array

    status = ffpcld(fitsfile%fptr,icol,frow,felement,nelements,col,stat)

  end subroutine fits_write_col_dbl

  subroutine fits_delete_rows(fitsfile,firstrow,nrows,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: firstrow, nrows
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat
    integer(kind=C_LONG_LONG) :: frow, nelements

    stat = status
    frow = firstrow
    nelements = nrows

    status = ffdrow(fitsfile%fptr,frow,nelements,stat)

  end subroutine fits_delete_rows


  !  Error Status Routines

  subroutine fits_report_error(unit,status)

    integer, intent(in) :: unit, status

    integer(kind=C_INT) :: stat
    character(len=FLEN_ERRMSG) :: errmsg
    character(len=FLEN_STATUS) :: status_str
    character(kind=C_CHAR), dimension(0:FLEN_STATUS) :: err_text
    character(kind=C_CHAR), dimension(0:FLEN_ERRMSG) :: err_msg

    if( status == 0 ) return

    stat = status
    call ffgerr(stat,err_text)
    call fchar(err_text,status_str)
    write(unit,'(3a,i0,a)') 'FITSIO error: "',trim(status_str), &
         '" (status = ',status,')'

    do
       if( ffgmsg(err_msg) == 0 ) exit
       call fchar(err_msg,errmsg)
       write(unit,'(2a)') 'FITSIO: ',trim(errmsg)
    end do

  end subroutine fits_report_error

  subroutine fits_write_errmark
    call ffpmrk
  end subroutine fits_write_errmark

  subroutine fits_clear_errmark
    call ffcmrk
  end subroutine fits_clear_errmark

  subroutine fits_clear_errmsg
    call ffcmsg
  end subroutine fits_clear_errmsg


  ! Utility Routines

  subroutine fits_get_keytype(val,dtype,status)

    character(len=*), intent(in) :: val
    character, intent(out) :: dtype
    integer, intent(in out) :: status

    integer(kind=C_INT) :: stat
    character(kind=C_CHAR), dimension(len(val)+1) :: string
    character(kind=C_CHAR) :: type

    stat = status
    call cchar(val,string)
    status = ffdtyp(string,type,stat)
    dtype = type

  end subroutine fits_get_keytype

  subroutine fits_make_keyn(keyroot,val,keyname,status)

    character(len=*), intent(in) :: keyroot
    integer, intent(in) :: val
    character(len=*), intent(out) :: keyname
    integer, intent(in out) :: status

    integer(kind=C_INT) :: n, stat
    character(kind=C_CHAR), dimension(len(keyroot)+1) :: root
    character(kind=C_CHAR), dimension(len(keyname)+1) :: name

    stat = status
    n = val
    call cchar(keyroot,root)
    status = ffkeyn(root,n,name,stat)
    call fchar(name,keyname)

  end subroutine fits_make_keyn

  subroutine fits_str2date(datestr,year,month,day,hour,minute,second,status)

    character(len=*), intent(in) :: datestr
    integer, intent(out) :: year,month,day,hour,minute
    real(REAL64), intent(out) :: second
    integer, intent(in out) :: status

    character(kind=C_CHAR), dimension(len_trim(datestr)+1) :: str
    integer(kind=C_INT) :: y,m,d,h,mnt,stat
    real(kind=C_DOUBLE) :: s

    stat = status
    call cchar(datestr,str)
    status = ffs2tm(str,y,m,d,h,mnt,s,stat)
    year = y
    month = m
    day = d
    hour = h
    minute = mnt
    second = s

  end subroutine fits_str2date


  ! Auxiliary utilities for string conversion: Fortran <-> C .

  pure subroutine fchar(cstring,fstring)

    character(kind=C_CHAR), dimension(*), intent(in) :: cstring
    character(len=*), intent(out) :: fstring

    integer :: i,n

    n = 1
    do i = 1, len(fstring)
       n = i
       if( cstring(n) == C_NULL_CHAR ) exit
       fstring(n:n) = cstring(n)
    end do
    fstring(n:) = ' '

  end subroutine fchar

  pure subroutine cchar(fstring,cstring)

    character(len=*), intent(in) :: fstring
    character(kind=C_CHAR), dimension(*), intent(out) :: cstring
    integer :: i, n

    n = len_trim(fstring)
    do i = 1, n
       cstring(i) = fstring(i:i)
    end do
    cstring(n+1) = C_NULL_CHAR

  end subroutine cchar


end module fitsio
