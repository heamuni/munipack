!
!  Fortran 2008+ interface for (c)FITSIO library
!  This module provides some high-level convenience functions.
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


module fitsio_mmviii

  use fitsio
  use iso_fortran_env
  implicit none


  interface fits_read_image
     module procedure fits_read_image_chr, fits_read_image_int, &
          fits_read_image_flt, fits_read_image_dbl
  end interface fits_read_image

  interface fits_write_image
     module procedure fits_write_image_chr, fits_write_image_int, &
          fits_write_image_flt, fits_write_image_dbl
  end interface fits_write_image

  interface fits_read_cube
     module procedure fits_read_cube_chr, fits_read_cube_int, &
          fits_read_cube_flt, fits_read_cube_dbl
  end interface fits_read_cube

  interface fits_write_cube
     module procedure fits_write_cube_chr, fits_write_cube_int, &
          fits_write_cube_flt, fits_write_cube_dbl
  end interface fits_write_cube


contains


  ! Convenience functions

  subroutine fits_insert_image(fitsfile,bitpix,naxis,naxes,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: bitpix, naxis
    integer, dimension(:), intent(in) :: naxes
    integer, intent(in out) :: status

    character(len=*), parameter :: com = 'value = BSCALE*image + BZERO'
    real(REAL64) :: bzero

    call fits_insert_img(fitsfile,bitpix,naxis,naxes,status)

    ! The integer image is scaled to the unsigned numerical range by
    ! the FITS convention (Pence et al. 2010, Table 11): BZERO = 2**(bitpix-1)
    ! for 16 and 32 bit integers. BITPIX = 8 is unsigned by definition,
    ! so no offset is applied (BZERO = 128 would turn it to signed bytes).
    if( bitpix > 8 ) then
       bzero = 2.0_REAL64**(bitpix-1)   ! 32768, 2147483648 exactly
       call fits_write_key(fitsfile,'BSCALE',1,com,status)
       call fits_write_key(fitsfile,'BZERO',bzero,-20,com,status)
    endif

  end subroutine fits_insert_image

  subroutine fits_get_image_size(fitsfile,naxes,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, dimension(:), intent(out) :: naxes
    integer, intent(in out) :: status

    call fits_get_img_size(fitsfile,size(naxes),naxes,status)

  end subroutine fits_get_image_size

  subroutine fits_read_image_int(fitsfile,group,nullval,image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    integer, intent(in) :: nullval
    integer, dimension(:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer :: dim1, naxis1, naxis2

    dim1 = size(image,1)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    call fits_read_2d_int(fitsfile,group,nullval,dim1,naxis1,naxis2, &
         image,anyf,status)

  end subroutine fits_read_image_int

  subroutine fits_read_image_chr(fitsfile,group,nullval,image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    character, intent(in) :: nullval
    character, dimension(:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer :: dim1, naxis1, naxis2

    dim1 = size(image,1)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    call fits_read_2d_chr(fitsfile,group,nullval,dim1,naxis1,naxis2, &
         image,anyf,status)

  end subroutine fits_read_image_chr

  subroutine fits_read_image_flt(fitsfile,group,nullval,image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    real(REAL32), intent(in) :: nullval
    real(REAL32), dimension(:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer :: dim1, naxis1, naxis2

    dim1 = size(image,1)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    call fits_read_2d_flt(fitsfile,group,nullval,dim1,naxis1,naxis2, &
         image,anyf,status)

  end subroutine fits_read_image_flt

  subroutine fits_read_image_dbl(fitsfile,group,nullval,image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    real(REAL64), intent(in) :: nullval
    real(REAL64), dimension(:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer :: dim1, naxis1, naxis2

    dim1 = size(image,1)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    call fits_read_2d_dbl(fitsfile,group,nullval,dim1,naxis1,naxis2, &
         image,anyf,status)

  end subroutine fits_read_image_dbl

  subroutine fits_write_image_chr(fitsfile,group,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    character, dimension(:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer :: dim1, naxis1, naxis2

    dim1 = size(image,1)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    call fits_write_2d_chr(fitsfile,group,dim1,naxis1,naxis2,image,status)

  end subroutine fits_write_image_chr

  subroutine fits_write_image_int(fitsfile,group,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    integer, dimension(:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer :: dim1, naxis1, naxis2

    dim1 = size(image,1)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    call fits_write_2d_int(fitsfile,group,dim1,naxis1,naxis2,image,status)

  end subroutine fits_write_image_int

  subroutine fits_write_image_flt(fitsfile,group,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    real(REAL32), dimension(:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer :: dim1, naxis1, naxis2

    dim1 = size(image,1)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    call fits_write_2d_flt(fitsfile,group,dim1,naxis1,naxis2,image,status)

  end subroutine fits_write_image_flt

  subroutine fits_write_image_dbl(fitsfile,group,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    real(REAL64), dimension(:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer :: dim1, naxis1, naxis2

    dim1 = size(image,1)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    call fits_write_2d_dbl(fitsfile,group,dim1,naxis1,naxis2,image,status)

  end subroutine fits_write_image_dbl


  subroutine fits_read_cube_chr(fitsfile,group,nullval,image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    character, intent(in) :: nullval
    character, dimension(:,:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer :: dim1, dim2, naxis1, naxis2, naxis3

    dim1 = size(image,1)
    dim2 = size(image,2)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    naxis3 = size(image,3)
    call fits_read_3d_chr(fitsfile,group,nullval,dim1,dim2,naxis1,naxis2,naxis3,&
         image,anyf,status)

  end subroutine fits_read_cube_chr

  subroutine fits_read_cube_int(fitsfile,group,nullval,image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    integer, intent(in) :: nullval
    integer, dimension(:,:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer :: dim1, dim2, naxis1, naxis2, naxis3

    dim1 = size(image,1)
    dim2 = size(image,2)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    naxis3 = size(image,3)
    call fits_read_3d_int(fitsfile,group,nullval,dim1,dim2,naxis1,naxis2,naxis3,&
         image,anyf,status)

  end subroutine fits_read_cube_int

  subroutine fits_read_cube_flt(fitsfile,group,nullval,image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    real(REAL32), intent(in) :: nullval
    real(REAL32), dimension(:,:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer :: dim1, dim2, naxis1, naxis2, naxis3

    dim1 = size(image,1)
    dim2 = size(image,2)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    naxis3 = size(image,3)
    call fits_read_3d_flt(fitsfile,group,nullval,dim1,dim2,naxis1,naxis2,naxis3,&
         image,anyf,status)

  end subroutine fits_read_cube_flt

  subroutine fits_read_cube_dbl(fitsfile,group,nullval,image,anyf,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    real(REAL64), intent(in) :: nullval
    real(REAL64), dimension(:,:,:), intent(out) :: image
    logical, intent(out) :: anyf
    integer, intent(in out) :: status

    integer :: dim1, dim2, naxis1, naxis2, naxis3

    dim1 = size(image,1)
    dim2 = size(image,2)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    naxis3 = size(image,3)
    call fits_read_3d_dbl(fitsfile,group,nullval,dim1,dim2,naxis1,naxis2,naxis3,&
         image,anyf,status)

  end subroutine fits_read_cube_dbl

  subroutine fits_write_cube_chr(fitsfile,group,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    character, dimension(:,:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer :: dim1, dim2, naxis1, naxis2, naxis3

    dim1 = size(image,1)
    dim2 = size(image,2)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    naxis3 = size(image,3)
    call fits_write_3d_chr(fitsfile,group,dim1,dim2,naxis1,naxis2,naxis3, &
         image,status)

  end subroutine fits_write_cube_chr

  subroutine fits_write_cube_int(fitsfile,group,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    integer, dimension(:,:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer :: dim1, dim2, naxis1, naxis2, naxis3

    dim1 = size(image,1)
    dim2 = size(image,2)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    naxis3 = size(image,3)
    call fits_write_3d_int(fitsfile,group,dim1,dim2,naxis1,naxis2,naxis3, &
         image,status)

  end subroutine fits_write_cube_int

  subroutine fits_write_cube_flt(fitsfile,group,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    real(REAL32), dimension(:,:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer :: dim1, dim2, naxis1, naxis2, naxis3

    dim1 = size(image,1)
    dim2 = size(image,2)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    naxis3 = size(image,3)
    call fits_write_3d_flt(fitsfile,group,dim1,dim2,naxis1,naxis2,naxis3, &
         image,status)

  end subroutine fits_write_cube_flt

  subroutine fits_write_cube_dbl(fitsfile,group,image,status)

    type(fitsfiles), intent(in) :: fitsfile
    integer, intent(in) :: group
    real(REAL64), dimension(:,:,:), intent(in) :: image
    integer, intent(in out) :: status

    integer :: dim1, dim2, naxis1, naxis2, naxis3

    dim1 = size(image,1)
    dim2 = size(image,2)
    naxis1 = size(image,1)
    naxis2 = size(image,2)
    naxis3 = size(image,3)
    call fits_write_3d_dbl(fitsfile,group,dim1,dim2,naxis1,naxis2,naxis3, &
         image,status)

  end subroutine fits_write_cube_dbl


  ! FITS file handling routines

  function fits_file_exist(filename) result(exist)

    character(len=*), intent(in) :: filename
    logical :: exist

    integer :: status
    type(fitsfiles) :: fitsfile

    status = 0
    call ftpmrk
    call fits_open_file(fitsfile,filename,FITS_READONLY,status)
    exist = status == 0
    call fits_close_file(fitsfile,status)
    call ftcmrk

  end function fits_file_exist

  subroutine fits_file_copy(source,destination,status)

    character(len=*), intent(in) :: source,destination
    integer, intent(in out) :: status

    type(fitsfiles) :: src,dst

    call fits_open_file(src,source,FITS_READONLY,status)
    call fits_create_file(dst,destination,status)
    call fits_copy_file(src,dst,1,1,1,status)
    call fits_close_file(src,status)
    if( status == 0 ) then
       call fits_close_file(dst,status)
    else
       call fits_delete_file(dst,status)
    end if
!    call fits_report_error(error_unit,status)

  end subroutine fits_file_copy

  subroutine fits_file_duplicate(fitsfile,destination,status)

    type(fitsfiles), intent(out) :: fitsfile
    character(len=*), intent(in) :: destination
    integer, intent(in out) :: status

    type(fitsfiles) :: dst

    call fits_create_file(dst,destination,status)
    call fits_copy_file(fitsfile,dst,1,1,1,status)
    if( status == 0 ) then
       call fits_close_file(dst,status)
    else
       call fits_delete_file(dst,status)
    end if

  end subroutine fits_file_duplicate

  subroutine fits_precopy_file(fitsfile,source,destination,rwmode,overwrite,status)

    type(fitsfiles), intent(out) :: fitsfile
    character(len=*), intent(in) :: source,destination
    integer, intent(in) :: rwmode
    logical, intent(in) :: overwrite
    integer, intent(in out) :: status

    type(fitsfiles) :: fits

    call fits_open_file(fits,source,FITS_READONLY,status)
    if( fits_file_exist(destination) .and. overwrite ) &
         call fits_file_delete(destination)
    call fits_file_duplicate(fits,destination,status)
    call fits_close_file(fits,status)
    call fits_open_file(fitsfile,destination,rwmode,status)

  end subroutine fits_precopy_file

  subroutine fits_file_delete(filename)

    character(len=*), intent(in) :: filename
    integer :: unit, iostat
    character(len=80) :: msg

    open(newunit=unit,file=filename,status='old',iostat=iostat,iomsg=msg)
    if( iostat == 0 ) then
       close(unit,status='DELETE')
    else
       write(error_unit,*) 'fits_file_delete: ',trim(msg)
    end if

  end subroutine fits_file_delete

  subroutine fits_create_scratch(fitsfile,status)

    type(fitsfiles), intent(out) :: fitsfile
    integer, intent(in out) :: status

    real(REAL64) :: x
    character(len=42) :: scratch

    call random_seed()
    status = 0
    do
       call random_number(x)
       write(scratch,'(a,f22.20,a)') 'fitsio_scratch_',x,'.fits'

       call fits_create_file(fitsfile,scratch,status)
       if( status == 0 ) exit
       if( status == FITS_FILE_NOT_CREATED ) status = 0
       if( status /= 0 ) exit
    end do

  end subroutine fits_create_scratch


end module fitsio_mmviii
