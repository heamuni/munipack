!
!  Purpose of this module is to provide API interfaces to C routines.
!  The short name versions of C routines are included only.
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

module cfitsio

  use, intrinsic :: iso_c_binding
  implicit none

  interface

     ! File Access Routines

     integer(kind=C_INT) function ffopen(fptr,filename,iomode,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), intent(out) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: filename
       integer(kind=C_INT), value, intent(in) :: iomode
       integer(kind=C_INT), intent(in out) :: status
     end function ffopen

     integer(kind=C_INT) function fftopn(fptr,filename,iomode,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), intent(out) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: filename
       integer(kind=C_INT), value, intent(in) :: iomode
       integer(kind=C_INT), intent(in out) :: status
     end function fftopn

     integer(kind=C_INT) function ffiopn(fptr,filename,iomode,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), intent(out) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: filename
       integer(kind=C_INT), value, intent(in) :: iomode
       integer(kind=C_INT), intent(in out) :: status
     end function ffiopn

     integer(kind=C_INT) function ffinit(fptr,filename,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), intent(out) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: filename
       integer(kind=C_INT), intent(in out) :: status
     end function ffinit

     integer(kind=C_INT) function ffclos(fptr,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), intent(in out) :: status
     end function ffclos

     integer(kind=C_INT) function ffdelt(fptr,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), intent(in out) :: status
     end function ffdelt

     ! HDU Access Routines

     integer(kind=C_INT) function ffmahd(fptr,hdunum,hdutype,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: hdunum
       integer(kind=C_INT), intent(out) :: hdutype
       integer(kind=C_INT), intent(in out) :: status
     end function ffmahd

     integer(kind=C_INT) function ffmrhd(fptr,nmove,hdutype,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: nmove
       integer(kind=C_INT), intent(out) :: hdutype
       integer(kind=C_INT), intent(in out) :: status
     end function ffmrhd

     integer(kind=C_INT) function ffmnhd(fptr,hdutype,extname,extver,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: hdutype, extver
       character(kind=C_CHAR), dimension(*), intent(in) :: extname
       integer(kind=C_INT), intent(in out) :: status
     end function ffmnhd

     integer(kind=C_INT) function ffthdu(fptr,hdunum,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), intent(out) :: hdunum
       integer(kind=C_INT), intent(in out) :: status
     end function ffthdu

     integer(kind=C_INT) function ffghdn(fptr,hdunum) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), intent(out) :: hdunum
     end function ffghdn

     integer(kind=C_INT) function ffcpfl(infptr,outfptr,previous,current,following,&
          status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: infptr,outfptr
       integer(kind=C_INT), value, intent(in) :: previous,current,following
       integer(kind=C_INT), intent(in out) :: status
     end function ffcpfl

     integer(kind=C_INT) function ffcopy(infptr,outfptr,morekeys,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: infptr,outfptr
       integer(kind=C_INT), value, intent(in) :: morekeys
       integer(kind=C_INT), intent(in out) :: status
     end function ffcopy

     integer(kind=C_INT) function ffcphd(infptr,outfptr,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: infptr,outfptr
       integer(kind=C_INT), intent(in out) :: status
     end function ffcphd

     integer(kind=C_INT) function ffdhdu(fptr,hdutype,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), intent(out) :: hdutype
       integer(kind=C_INT), intent(in out) :: status
     end function ffdhdu

     integer(kind=C_INT) function ffiimg(fptr,bitpix,naxis,naxes,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: bitpix, naxis
       integer(kind=C_LONG), dimension(*), intent(in) :: naxes
       integer(kind=C_INT), intent(in out) :: status
     end function ffiimg

     integer(kind=C_INT) function ffibin(fptr,nrows,tfields,ttype,tform,tunit,&
          extname,pcount,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG_LONG), value, intent(in) :: nrows
       integer(kind=C_INT), value, intent(in) :: tfields
       type(c_ptr), dimension(*), intent(in) :: ttype, tform, tunit
       character(kind=C_CHAR), dimension(*), intent(in) :: extname
       integer(kind=C_LONG_LONG), value, intent(in) :: pcount
       integer(kind=C_INT), intent(in out) :: status
     end function ffibin


     ! Header Keyword Read/Write Routines

     integer(kind=C_INT) function ffghsp(fptr,keysexist,morekeys,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), intent(out) :: keysexist,morekeys
       integer(kind=C_INT), intent(in out) :: status
     end function ffghsp

     integer(kind=C_INT) function ffgkey(fptr,keyname,val,com,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname
       character(kind=C_CHAR), dimension(*), intent(out) :: val,com
       integer(kind=C_INT), intent(in out) :: status
     end function ffgkey

     integer(kind=C_INT) function ffgkys(fptr,keyname,val,com,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname
       character(kind=C_CHAR), dimension(*), intent(out) :: val,com
       integer(kind=C_INT), intent(in out) :: status
     end function ffgkys

     integer(kind=C_INT) function ffgkyj(fptr,keyname,numval,com,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname
       integer(kind=C_LONG), intent(out) :: numval
       character(kind=C_CHAR), dimension(*), intent(out) :: com
       integer(kind=C_INT), intent(in out) :: status
     end function ffgkyj

     integer(kind=C_INT) function ffgkye(fptr,keyname,numval,com,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname
       real(kind=C_FLOAT), intent(out) :: numval
       character(kind=C_CHAR), dimension(*), intent(out) :: com
       integer(kind=C_INT), intent(in out) :: status
     end function ffgkye

     integer(kind=C_INT) function ffgkyd(fptr,keyname,numval,com,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname
       real(kind=C_DOUBLE), intent(out) :: numval
       character(kind=C_CHAR), dimension(*), intent(out) :: com
       integer(kind=C_INT), intent(in out) :: status
     end function ffgkyd

     integer(kind=C_INT) function ffpkys(fptr,keyname,val,com,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname,val,com
       integer(kind=C_INT), intent(in out) :: status
     end function ffpkys

     integer(kind=C_INT) function ffpkyj(fptr,keyname,numval,com,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname, com
       integer(kind=C_LONG_LONG), value, intent(in) :: numval
       integer(kind=C_INT), intent(in out) :: status
     end function ffpkyj

     integer(kind=C_INT) function ffpkye(fptr,keyname,numval,decimals,com,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname, com
       real(kind=C_FLOAT), value, intent(in) :: numval
       integer(kind=C_INT), value, intent(in) :: decimals
       integer(kind=C_INT), intent(in out) :: status
     end function ffpkye

     integer(kind=C_INT) function ffpkyd(fptr,keyname,numval,decimals,com,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname, com
       real(kind=C_DOUBLE), value, intent(in) :: numval
       integer(kind=C_INT), value, intent(in) :: decimals
       integer(kind=C_INT), intent(in out) :: status
     end function ffpkyd

     integer(kind=C_INT) function ffukys(fptr,keyname,val,com,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname,val,com
       integer(kind=C_INT), intent(in out) :: status
     end function ffukys

     integer(kind=C_INT) function ffukyj(fptr,keyname,numval,com,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname, com
       integer(kind=C_LONG_LONG), value, intent(in) :: numval
       integer(kind=C_INT), intent(in out) :: status
     end function ffukyj

     integer(kind=C_INT) function ffukye(fptr,keyname,numval,decimals,com,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname, com
       real(kind=C_FLOAT), value, intent(in) :: numval
       integer(kind=C_INT), value, intent(in) :: decimals
       integer(kind=C_INT), intent(in out) :: status
     end function ffukye

     integer(kind=C_INT) function ffukyd(fptr,keyname,numval,decimals,com,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname, com
       real(kind=C_DOUBLE), value, intent(in) :: numval
       integer(kind=C_INT), value, intent(in) :: decimals
       integer(kind=C_INT), intent(in out) :: status
     end function ffukyd

     integer(kind=C_INT) function ffucrd(fptr,keyname,card,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname, card
       integer(kind=C_INT), intent(in out) :: status
     end function ffucrd

     integer(kind=C_INT) function ffgunt(fptr,keyname,unit,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname
       character(kind=C_CHAR), dimension(*), intent(out) :: unit
       integer(kind=C_INT), intent(in out) :: status
     end function ffgunt

     integer(kind=C_INT) function ffdkey(fptr,keyname,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: keyname
       integer(kind=C_INT), intent(in out) :: status
     end function ffdkey

     integer(kind=C_INT) function ffpcom(fptr,comment,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: comment
       integer(kind=C_INT), intent(in out) :: status
     end function ffpcom

     integer(kind=C_INT) function ffphis(fptr,history,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       character(kind=C_CHAR), dimension(*), intent(in) :: history
       integer(kind=C_INT), intent(in out) :: status
     end function ffphis

     integer(kind=C_INT) function ffgrec(fptr,keynum,card,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(C_INT), value, intent(in) :: keynum
       character(kind=C_CHAR), dimension(*), intent(out) :: card
       integer(kind=C_INT), intent(in out) :: status
     end function ffgrec

     integer(kind=C_INT) function ffdrec(fptr,keynum,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(C_INT), value, intent(in) :: keynum
       integer(kind=C_INT), intent(in out) :: status
     end function ffdrec

     integer(kind=C_INT) function ffgthd(templt,card,keytype,status) bind(c)
       use, intrinsic :: iso_c_binding
       character(kind=C_CHAR), dimension(*), intent(in) :: templt
       character(kind=C_CHAR), dimension(*), intent(out) :: card
       integer(C_INT), intent(out) :: keytype
       integer(kind=C_INT), intent(in out) :: status
     end function ffgthd

     integer(kind=C_INT) function ffgknm(card,keyname,keylength,status) bind(c)
       use, intrinsic :: iso_c_binding
       character(kind=C_CHAR), dimension(*), intent(in) :: card
       character(kind=C_CHAR), dimension(*), intent(out) :: keyname
       integer(C_INT), intent(out) :: keylength
       integer(kind=C_INT), intent(in out) :: status
     end function ffgknm

     integer(kind=C_INT) function ffhdr2str(fptr,nocomments,exclist,nexc, &
          header,nkeys,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: nocomments
       type(c_ptr), value, intent(in) :: exclist
       integer(kind=C_INT), value, intent(in) :: nexc
       type(c_ptr), intent(out) :: header
       integer(kind=C_INT), intent(out) :: nkeys
       integer(kind=C_INT), intent(in out) :: status
     end function ffhdr2str

     integer(kind=C_INT) function fffree(phdr,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: phdr
       integer(kind=C_INT), intent(in out) :: status
     end function fffree

     ! Primary Array or IMAGE Extension I/O Routines

     integer(kind=C_INT) function ffgidt(fptr,bitpix,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), intent(out) :: bitpix
       integer(kind=C_INT), intent(in out) :: status
     end function ffgidt

     integer(kind=C_INT) function ffgidm(fptr,naxis,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), intent(out) :: naxis
       integer(kind=C_INT), intent(in out) :: status
     end function ffgidm

     integer(kind=C_INT) function ffgisz(fptr,maxdim,naxes,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: maxdim
       integer(kind=C_LONG), dimension(*), intent(out) :: naxes
       integer(kind=C_INT), intent(in out) :: status
     end function ffgisz

     integer(kind=C_INT) function ffgipr(fptr,maxdim,bitpix,naxis,naxes,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: maxdim
       integer(kind=C_INT), intent(out) :: bitpix, naxis
       integer(kind=C_LONG), dimension(*), intent(out) :: naxes
       integer(kind=C_INT), intent(in out) :: status
     end function ffgipr

     integer(kind=C_INT) function ffgpvb(fptr,group,fpixel,nelements,nulval,&
          array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: fpixel, nelements
       character(kind=C_CHAR), value, intent(in) :: nulval
       character(kind=C_CHAR), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffgpvb

     integer(kind=C_INT) function ffgpvk(fptr,group,fpixel,nelements,nulval,&
          array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: fpixel, nelements
       integer(kind=C_INT), value, intent(in) :: nulval
       integer(kind=C_INT), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffgpvk

     integer(kind=C_INT) function ffgpve(fptr,group,fpixel,nelements,nulval,&
          array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: fpixel, nelements
       real(kind=C_FLOAT), value, intent(in) :: nulval
       real(kind=C_FLOAT), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffgpve

     integer(kind=C_INT) function ffgpvd(fptr,group,fpixel,nelements,nulval,&
          array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: fpixel, nelements
       real(kind=C_DOUBLE), value, intent(in) :: nulval
       real(kind=C_DOUBLE), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffgpvd

     integer(kind=C_INT) function ffpprb(fptr,group,fpixel,nelements,array,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: fpixel, nelements
       character(kind=C_CHAR), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffpprb

     integer(kind=C_INT) function ffpprk(fptr,group,fpixel,nelements,array,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: fpixel, nelements
       integer(kind=C_INT), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffpprk

     integer(kind=C_INT) function ffppre(fptr,group,fpixel,nelements,array,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: fpixel, nelements
       real(kind=C_FLOAT), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffppre

     integer(kind=C_INT) function ffpprd(fptr,group,fpixel,nelements,array,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: fpixel, nelements
       real(kind=C_DOUBLE), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffpprd

     integer(kind=C_INT) function ffg2db(fptr,group,nullval,dim1, &
          naxis1,naxis2,array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,naxis1,naxis2
       character(kind=C_CHAR), value, intent(in) :: nullval
       character(kind=C_CHAR), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffg2db

     integer(kind=C_INT) function ffg2dk(fptr,group,nullval, &
          dim1,naxis1,naxis2,array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,naxis1,naxis2
       integer(kind=C_INT), value, intent(in) :: nullval
       integer(kind=C_INT), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffg2dk

     integer(kind=C_INT) function ffg2de(fptr,group,nullval, &
          dim1,naxis1,naxis2,array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,naxis1,naxis2
       real(kind=C_FLOAT), value, intent(in) :: nullval
       real(kind=C_FLOAT), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffg2de

     integer(kind=C_INT) function ffg2dd(fptr,group,nullval, &
          dim1,naxis1,naxis2,array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,naxis1,naxis2
       real(kind=C_DOUBLE), value, intent(in) :: nullval
       real(kind=C_DOUBLE), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffg2dd

     integer(kind=C_INT) function ffp2db(fptr,group,dim1,naxis1,naxis2,array,status)&
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,naxis1,naxis2
       character(kind=C_CHAR), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffp2db

     integer(kind=C_INT) function ffp2dk(fptr,group,dim1,naxis1,naxis2,array,status)&
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,naxis1,naxis2
       integer(kind=C_INT), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffp2dk

     integer(kind=C_INT) function ffp2de(fptr,group,dim1,naxis1,naxis2,array,status)&
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,naxis1,naxis2
       real(kind=C_FLOAT), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffp2de

     integer(kind=C_INT) function ffp2dd(fptr,group,dim1,naxis1,naxis2,array,status)&
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,naxis1,naxis2
       real(kind=C_DOUBLE), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffp2dd

     integer(kind=C_INT) function ffg3db(fptr,group,nullval,dim1,dim2, &
          naxis1,naxis2,naxis3,array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,dim2,naxis1,naxis2,naxis3
       character(kind=C_CHAR), value, intent(in) :: nullval
       character(kind=C_CHAR), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffg3db

     integer(kind=C_INT) function ffg3dk(fptr,group,nullval, &
          dim1,dim2,naxis1,naxis2,naxis3,array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,dim2,naxis1,naxis2,naxis3
       integer(kind=C_INT), value, intent(in) :: nullval
       integer(kind=C_INT), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffg3dk

     integer(kind=C_INT) function ffg3de(fptr,group,nullval, &
          dim1,dim2,naxis1,naxis2,naxis3,array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,dim2,naxis1,naxis2,naxis3
       real(kind=C_FLOAT), value, intent(in) :: nullval
       real(kind=C_FLOAT), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffg3de

     integer(kind=C_INT) function ffg3dd(fptr,group,nullval, &
          dim1,dim2,naxis1,naxis2,naxis3,array,anynul,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,dim2,naxis1,naxis2,naxis3
       real(kind=C_DOUBLE), value, intent(in) :: nullval
       real(kind=C_DOUBLE), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynul
       integer(kind=C_INT), intent(in out) :: status
     end function ffg3dd

     integer(kind=C_INT) function ffp3db(fptr,group,dim1,dim2,naxis1,naxis2,naxis3,&
          array,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,dim2,naxis1,naxis2,naxis3
       character(kind=C_CHAR), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffp3db

     integer(kind=C_INT) function ffp3dk(fptr,group,dim1,dim2,naxis1,naxis2,naxis3,&
          array,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,dim2,naxis1,naxis2,naxis3
       integer(kind=C_INT), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffp3dk

     integer(kind=C_INT) function ffp3de(fptr,group,dim1,dim2,naxis1,naxis2,naxis3,&
          array,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,dim2,naxis1,naxis2,naxis3
       real(kind=C_FLOAT), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffp3de

     integer(kind=C_INT) function ffp3dd(fptr,group,dim1,dim2,naxis1,naxis2,naxis3,&
          array,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), value, intent(in) :: group
       integer(kind=C_LONG_LONG), value, intent(in) :: dim1,dim2,naxis1,naxis2,naxis3
       real(kind=C_DOUBLE), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffp3dd

     ! ASCII and Binary Table Routines

     integer(kind=C_INT) function ffgnrw(fptr,nrows,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), intent(out) :: nrows
       integer(kind=C_INT), intent(in out) :: status
     end function ffgnrw

     integer(kind=C_INT) function ffgncl(fptr,ncols,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), intent(out) :: ncols
       integer(kind=C_INT), intent(in out) :: status
     end function ffgncl

     integer(kind=C_INT) function ffgrsz(fptr,nrows,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG), intent(out) :: nrows
       integer(kind=C_INT), intent(in out) :: status
     end function ffgrsz

     integer(kind=C_INT) function ffgcno(fptr,casesen,templt,colnum,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: casesen
       character(kind=C_CHAR), dimension(*), intent(out) :: templt
       integer(kind=C_INT), intent(out) :: colnum
       integer(kind=C_INT), intent(in out) :: status
     end function ffgcno

     integer(kind=C_INT) function ffgcnn(fptr,casesen,templt,colname,colnum,status)&
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: casesen
       character(kind=C_CHAR), dimension(*), intent(in) :: templt
       character(kind=C_CHAR), dimension(*), intent(out) :: colname
       integer(kind=C_INT), intent(out) :: colnum
       integer(kind=C_INT), intent(in out) :: status
     end function ffgcnn

     integer(kind=C_INT) function fficol(fptr,colnum,ttype,tform,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       character(kind=C_CHAR), dimension(*), intent(in) :: ttype, tform
       integer(kind=C_INT), intent(in out) :: status
     end function fficol

     integer(kind=C_INT) function fficls(fptr,colnum,ncols,ttype,tform,status) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum, ncols
       type(c_ptr), dimension(*), intent(in) :: ttype, tform
       integer(kind=C_INT), intent(in out) :: status
     end function fficls

     integer(kind=C_INT) function ffdcol(fptr,colnum,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_INT), intent(in out) :: status
     end function ffdcol

     integer(kind=C_INT) function ffgcvs(fptr,colnum,firstrow,firstelem,nelements,&
          nullstr,array,anynull,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       character(kind=C_CHAR), dimension(*), intent(in) :: nullstr
       type(c_ptr), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynull
       integer(kind=C_INT), intent(in out) :: status
     end function ffgcvs

     integer(kind=C_INT) function ffgcvl(fptr,colnum,firstrow,firstelem,nelements,&
          nullval,array,anynull,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       character(kind=C_CHAR), value, intent(in) :: nullval
       character(kind=C_CHAR), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynull
       integer(kind=C_INT), intent(in out) :: status
     end function ffgcvl

     integer(kind=C_INT) function ffgcvj(fptr,colnum,firstrow,firstelem,nelements,&
          nullval,array,anynull,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       integer(kind=C_LONG), value, intent(in) :: nullval
       integer(kind=C_LONG), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynull
       integer(kind=C_INT), intent(in out) :: status
     end function ffgcvj

     integer(kind=C_INT) function ffgcve(fptr,colnum,firstrow,firstelem,nelements,&
          nullval,array,anynull,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       real(kind=C_FLOAT), value, intent(in) :: nullval
       real(kind=C_FLOAT), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynull
       integer(kind=C_INT), intent(in out) :: status
     end function ffgcve

     integer(kind=C_INT) function ffgcvd(fptr,colnum,firstrow,firstelem,nelements,&
          nullval,array,anynull,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       real(kind=C_DOUBLE), value, intent(in) :: nullval
       real(kind=C_DOUBLE), dimension(*), intent(out) :: array
       integer(kind=C_INT), intent(out) :: anynull
       integer(kind=C_INT), intent(in out) :: status
     end function ffgcvd

     integer(kind=C_INT) function ffpcls(fptr,colnum,firstrow,firstelem,nelements,&
          array,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       type(c_ptr), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffpcls

     integer(kind=C_INT) function ffpcll(fptr,colnum,firstrow,firstelem,nelements,&
          array,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       character(kind=C_CHAR), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffpcll

     integer(kind=C_INT) function ffpclj(fptr,colnum,firstrow,firstelem,nelements,&
          array,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       integer(kind=C_LONG), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffpclj

     integer(kind=C_INT) function ffpcle(fptr,colnum,firstrow,firstelem,nelements,&
          array,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       real(kind=C_FLOAT), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffpcle

     integer(kind=C_INT) function ffpcld(fptr,colnum,firstrow,firstelem,nelements,&
          array,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_INT), value, intent(in) :: colnum
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, firstelem, nelements
       real(kind=C_DOUBLE), dimension(*), intent(in) :: array
       integer(kind=C_INT), intent(in out) :: status
     end function ffpcld

     integer(kind=C_INT) function ffdrow(fptr,firstrow,nelements,status) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: fptr
       integer(kind=C_LONG_LONG), value, intent(in) :: firstrow, nelements
       integer(kind=C_INT), intent(in out) :: status
     end function ffdrow

     !  Error Status Routines

     subroutine ffgerr(status,err_text) bind(c)
       use, intrinsic :: iso_c_binding
       integer(kind=C_INT), value, intent(in) :: status
       character(kind=C_CHAR), dimension(*), intent(out) :: err_text
     end subroutine ffgerr

     integer(kind=C_INT) function ffgmsg(err_msg) bind(c)
       use, intrinsic :: iso_c_binding
       character(kind=C_CHAR), dimension(*), intent(out) :: err_msg
     end function ffgmsg

     subroutine ffpmrk() bind(c)
     end subroutine ffpmrk

     subroutine ffcmrk() bind(c)
     end subroutine ffcmrk

     subroutine ffcmsg() bind(c)
     end subroutine ffcmsg

     ! Utility Routines

     integer(kind=C_INT) function ffdtyp(val,dtype,status) bind(c)
       use, intrinsic :: iso_c_binding
       character(kind=C_CHAR), dimension(*), intent(in) :: val
       character(kind=C_CHAR), intent(out) :: dtype
       integer(kind=C_INT), intent(in out) :: status
     end function ffdtyp

     integer(kind=C_INT) function ffkeyn(keyroot,val,keyname,status) bind(c)
       use, intrinsic :: iso_c_binding
       character(kind=C_CHAR), dimension(*), intent(in) :: keyroot
       integer(kind=C_INT), value, intent(in) :: val
       character(kind=C_CHAR), dimension(*), intent(out) :: keyname
       integer(kind=C_INT), intent(in out) :: status
     end function ffkeyn

     integer(kind=C_INT) function ffs2tm(datestr,year,month,day,hour,minute,second,&
          status) bind(c)
       use, intrinsic :: iso_c_binding
       character(kind=C_CHAR), dimension(*), intent(in) :: datestr
       integer(kind=C_INT), intent(out) :: year,month,day,hour,minute
       real(kind=C_DOUBLE), intent(out) :: second
       integer(kind=C_INT), intent(in out) :: status
     end function ffs2tm

  end interface

end module cfitsio
