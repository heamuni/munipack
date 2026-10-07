!
!  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
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

module flatmodule

  implicit none

contains

  subroutine flatmodel_load(filename,image,status)

    use titsio
    use iso_fortran_env

    character(len=*), intent(in) :: filename
    real, dimension(:,:), intent(out) :: image
    integer, intent(out) :: status

    type(fitsfiles) :: fits
    integer, parameter :: group = 1
    real, parameter :: null = 0.0
    integer :: naxis
    integer, dimension(2) :: naxes
    logical :: anyf

    status = 0
    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) goto 661

    call fits_get_img_dim(fits,naxis,status)
    if( status /= 0 ) goto 661
    if( .not. (naxis == 2) ) then
       write(error_unit,*) 'Error: The dimensions of the flat-model should be 2'
       goto 661
    end if

    call fits_get_img_size(fits,2,naxes,status)
    if( status /= 0 ) goto 661
    if( .not. all(naxes == shape(image)) ) then
       write(error_unit,*) 'Error: Dimensions of the flat-model does not match.'
       goto 661
    end if

    call fits_read_image(fits,group,null,image,anyf,status)

661 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( status  /= 0 )  write(error_unit,*) &
         "Warning: flat-model `",trim(filename),"' load failed."

  end subroutine flatmodel_load

end module flatmodule
