!
!  FITS subroutines for pre-corrections
!
!  Copyright © 2018-2022 F.Hroch (hroch@physics.muni.cz)
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

module fitscorr

  use titsio
  use iso_fortran_env

  implicit none

  type :: CorrFits
     character(len=FLEN_FILENAME) :: filename
     character(len=FLEN_VALUE) :: dateobs='', imagetyp='', filter=''
     real, dimension(:,:), allocatable :: image, imgerr
     integer :: naxis, bitpix
     integer, dimension(2) :: naxes
     real :: exptime = 0
     real :: saturate = huge(1.0)*(1.0-10*epsilon(1.0))
     real :: temper = -666
     real :: gain = 1
     real :: readns = 0              ! readout noise [ADU]
     logical :: saturate_set = .false.
     logical :: gain_set = .false.
     logical :: exptime_set = .false.
     logical :: status = .false.
     real :: mean, stderr, sig ! flat only
   contains
     procedure :: Load, get_saturate
  end type CorrFits

contains

  subroutine Load(fits,filename,keys)

    class(CorrFits) :: fits
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in), optional :: keys

    integer, parameter :: extver = 0, group = 1
    real, parameter :: null = 0.0
    integer :: status, stat
    logical :: anyf
    character(len=80) :: msg
    type(fitsfiles) :: fitsfile

    status = 0
    fits%filename = filename
    call fits_open_image(fitsfile,filename,FITS_READONLY,status)
    if( status /= 0  ) then
       write(error_unit,*) 'Error: failed to read the file `',trim(filename),"'."
       fits%status = .false.
       return
    end if

    call fits_get_img_type(fitsfile,fits%bitpix,status)
    call fits_get_img_dim(fitsfile,fits%naxis,status)
    if( status /= 0 ) goto 666

    if( fits%naxis /= 2 ) then
       write(error_unit,*) trim(filename),": Assertion failed: naxis /= 2"
       goto 666
    end if

    call fits_get_img_size(fitsfile,2,fits%naxes,status)

    if( present(keys) ) then

       ! dateobs
       call fits_write_errmark
       call fits_read_key(fitsfile,keys(1),fits%dateobs,status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
          fits%dateobs = ''
       end if

       ! exposure time
       call fits_write_errmark
       call fits_read_key(fitsfile,keys(2),fits%exptime,status)
       if( status == 0 ) then
          fits%exptime_set = .true.
       else if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
          fits%exptime = -1
       end if

       ! filter
       call fits_write_errmark
       call fits_read_key(fitsfile,keys(3),fits%filter,status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
          fits%filter = ''
       end if

       call fits_write_errmark
       call fits_read_key(fitsfile,keys(4),fits%saturate,status)
       if( status == 0 ) then
          fits%saturate_set = .true.
       else if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
       end if

       ! device temperature
       call fits_write_errmark
       call fits_read_key(fitsfile,keys(5),fits%temper,status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
          fits%temper = -666
       end if

       ! gain
       call fits_write_errmark
       call fits_read_key(fitsfile,keys(6),fits%gain,status)
       if( status == 0 ) then
          fits%gain_set = .true.
       else if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
          fits%gain = 1
       end if

       ! readout noise
       call fits_write_errmark
       call fits_read_key(fitsfile,FITS_KEY_READNS,fits%readns,status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
          fits%readns = 0
       end if

       ! imagetyp
       call fits_write_errmark
       call fits_read_key(fitsfile,FITS_KEY_IMAGETYP,fits%imagetyp,status)
       if( status == FITS_KEYWORD_NOT_FOUND ) then
          call fits_clear_errmark
          status = 0
          fits%imagetyp = ''
       end if

    end if ! keys
    if( status /= 0 ) goto 666

    ! data
    allocate(fits%image(fits%naxes(1),fits%naxes(2)), &
         fits%imgerr(fits%naxes(1),fits%naxes(2)),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Failed to allocate images.'
    end if

    call fits_read_image(fitsfile,group,null,fits%image,anyf,status)
    if( status /= 0 ) goto 666

    ! errors
    call fits_write_errmark
    call fits_movnam_hdu(fitsfile,FITS_IMAGE_HDU,EXT_STDERR,extver,status)
    if( status == 0 ) then
       call fits_read_image(fitsfile,group,null,fits%imgerr,anyf,status)
       if( status /= 0 ) goto 666
    else if( status == FITS_BAD_HDU_NUM ) then
       ! if the information about standard deviation is not available,
       ! we are continuing with zeros
       fits%imgerr = 0
       status = 0
       call fits_clear_errmark
    end if

    ! initial saturation setup
    if( .not. fits%saturate_set ) then
       if( status == 0 .and. fits%bitpix > 0 ) then
          fits%saturate = 2.0**fits%bitpix - 1
       else
          fits%saturate = huge(fits%image) * (1.0 - 10*epsilon(fits%image))
       end if
    end if

666 continue

    if( status /= 0 ) then
       if( allocated(fits%image) ) deallocate(fits%image)
       if( allocated(fits%imgerr) ) deallocate(fits%imgerr)
    end if

    call fits_close_file(fitsfile,status)
    call fits_report_error(error_unit,status)

    fits%status = status == 0

  end subroutine Load


  real function get_saturate(this,gain)

    class(CorrFits) :: this
    real, intent(in) :: gain

    associate( saturate => this%saturate )

      get_saturate = saturate

      if( gain > 1 ) then

         if( saturate < huge(saturate) / (1.001*gain) ) &
              get_saturate = gain * saturate

      else if( gain < 1 ) then

         if( saturate < huge(saturate) * (1.001*gain) ) &
              get_saturate = gain * saturate

      end if

    end associate

  end function get_saturate

end module fitscorr
