!
!  Read catalogue
!
!
!  Copyright © 2016-2026 F.Hroch (hroch@physics.muni.cz)
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

!  *  proper motion is omited


module catio

  use iso_fortran_env
  implicit none

contains

  subroutine catalogue(cat,labels,alpha,delta,mag,catid)

    use titsio

    character(len=*), intent(in) :: cat
    character(len=*), dimension(:), intent(in) :: labels
    real(REAL64), dimension(:), allocatable, intent(out) :: alpha,delta,mag
    character(len=*), intent(out) :: catid

    real(REAL64), parameter :: nullval = 99.99999_REAL64

    integer :: nrows, srows, status, i, l, frow
    integer, dimension(size(labels)) :: cols
    logical :: anyf
    type(fitsfiles) :: fits

    status = 0

    ! open and move to a table extension
    call fits_open_table(fits,cat,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read a table in the file `',trim(cat),"'."
       stop 'CATSIO'
    end if

    call fits_get_num_rows(fits,nrows,status)
    if( status /= 0 ) goto 666
    if( nrows == 0 ) stop 'catio: no stars available.'

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
    if( status /= 0 ) goto 666

    allocate(alpha(nrows),delta(nrows),mag(nrows))

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,cols(1),frow,nullval,alpha(i:l),anyf,status)
       call fits_read_col(fits,cols(3),frow,nullval,delta(i:l),anyf,status)
       call fits_read_col(fits,cols(5),frow,nullval,mag(i:l),anyf,status)
       if( status /= 0 ) goto 666
    end do

    call fits_close_file(fits,status)
    return

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( allocated(mag) ) deallocate(mag,alpha,delta)

    stop 'CATSIO'

  end subroutine catalogue

end module catio
