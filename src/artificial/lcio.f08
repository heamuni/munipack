!
!  Light curve parameters by FITS files
!
!
!  Copyright © 2017-26 F.Hroch (hroch@physics.muni.cz)
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


module lcio

  use titsio
  use iso_fortran_env
  implicit none

contains

  subroutine LoadLCTable(tablename,time,flux,dflux)

    character(len=*), intent(in) :: tablename
    real(REAL64), dimension(:), allocatable, intent(out) :: time,flux,dflux

    real(REAL64), parameter :: nullval = 0.0_REAL64

    integer :: nrows, srows, status, i, l, frow
    character(len=FLEN_CARD) :: qname
    logical :: anyf, check
    type(fitsfiles) :: fits

    status = 0

    ! open and move to a table extension
    call fits_open_table(fits,tablename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read a table in the file `',trim(tablename),"'."
       goto 666
    end if

    call fits_get_num_rows(fits,nrows,status)
    if( status /= 0 ) goto 666
    if( nrows == 0 ) stop 'LoadLCTable: There is no a timeserie.'

    allocate(time(nrows),flux(nrows),dflux(nrows))

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,1,frow,nullval,time(i:l),anyf,status)
       if( status /= 0 ) then
          write(error_unit,*) "Error: the read of the first column (time) of `", &
               trim(tablename),"' failed."
          goto 666
       end if

       call fits_read_col(fits,2,frow,nullval,flux(i:l),anyf,status)
       call fits_read_col(fits,3,frow,nullval,dflux(i:l),anyf,status)
       if( status /= 0 ) then
          write(error_unit,*) "Error: the read of the second, or third,", &
               "column of '",trim(tablename)," failed."
          goto 666
       end if
    end do

    if( .not. all(dflux > 0) ) then
       write(error_unit,*) "Error: the flux error (third) column of '", &
            trim(tablename)," should be positive numbers."
       stop 'LoadLCTable: .not. all(dflux > 0)'
    end if

    check = .false.
    do i = 2, size(time)
       if( time(i) < time(i-1) ) check = .true.
    end do
    if( check ) write(error_unit,*) "Warning: non-increasing time sequence."

    call fits_close_file(fits,status)
    return

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    stop 'LoadLCTable'

  end subroutine LoadLCTable


  subroutine LoadFourier(tablename,c)

    character(len=*), intent(in) :: tablename
    complex(REAL64), dimension(:), allocatable, intent(out) :: c

    real(REAL64), parameter :: nullval = 0.0_REAL64

    integer :: nrows, srows, status, frow, i, l
    real(REAL64), dimension(:), allocatable :: a,b
    logical :: anyf
    type(fitsfiles) :: fits

    status = 0

    ! open and move to a table extension
    call fits_open_table(fits,tablename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) 'Error: failed to read a table in `',trim(tablename),"'."
       goto 666
    end if

    call fits_get_num_rows(fits,nrows,status)
    if( status /= 0 ) goto 666

    if( nrows == 0 ) stop 'LoadFourier: There are no Fourier coefficients.'

    allocate(a(nrows),b(nrows),c(0:nrows-1))

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,1,frow,nullval,a(i:l),anyf,status)
       call fits_read_col(fits,2,frow,nullval,b(i:l),anyf,status)
       if( status /= 0 ) goto 666
    end do

    c = cmplx(a,b,REAL64)

    call fits_close_file(fits,status)
    return

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    stop 'LoadFourier'

  end subroutine LoadFourier

end module lcio
