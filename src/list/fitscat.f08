!
!  FITS catalogue for timeseries
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

module fitscat

  use titsio
  use iso_fortran_env

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: J2000 = 2451545.0_dbl   ! J2000.0

  type :: CatFits
     character(len=FLEN_FILENAME) :: filename = ''
     real(dbl) :: jd = 0
     real(dbl), dimension(:,:), allocatable :: coo, pm
     logical :: status = .false.
   contains
     procedure :: Catalogue,Coords,ZeroInit,getcurrent
  end type CatFits

contains

  subroutine ZeroInit(this)

    class(CatFits) :: this

    this%filename = ''
    this%jd = J2000
    allocate(this%coo(0,0),this%pm(0,0))
    this%status = .false.

  end subroutine ZeroInit

  subroutine Coords(this,x,y)

    class(CatFits) :: this
    real(dbl), dimension(:), intent(in) :: x,y

    integer :: n, stat
    character(len=80) :: msg

    n = size(x)

    allocate(this%coo(n,2),this%pm(n,2),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) trim(msg)
       error stop 'Failed to allocate arrays.'
    end if

    this%coo(:,1) = x
    this%coo(:,2) = y
    this%pm = 0
    this%jd = J2000
    this%status = .true.

  end subroutine Coords


  subroutine Catalogue(fits,cat)

    class(CatFits) :: fits
    character(len=*), intent(in) :: cat

    real(dbl), parameter :: xpm = 365.25_dbl * 3600.0_dbl
    integer, parameter :: nlabels = 4
    character(len=FLEN_VALUE), dimension(nlabels), parameter :: &
         labels = [ FITS_COL_RA, FITS_COL_DEC, FITS_COL_PMRA, FITS_COL_PMDEC ]

    real(REAL64), parameter :: null = real(0.0,REAL64)
    integer :: nrows, ncols, i, l, status, stat, srows, frow
    integer, dimension(nlabels) :: col
    character(len=FLEN_CARD) :: runits,dunits,buf
    character(len=80) :: msg
    type(fitsfiles) :: fitsfile
    logical :: anyf

    fits%filename = cat

    status = 0

    call fits_open_table(fitsfile,cat,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) 'Error: failed to read the file `',trim(cat),"'."
       fits%status = .false.
       return
    end if

    call fits_write_errmark
    call fits_read_key(fitsfile,FITS_KEY_EPOCH,fits%jd,buf,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       call fits_clear_errmark
       status = 0
       fits%jd = J2000
    end if

    call fits_get_num_cols(fitsfile,ncols,status)
    call fits_get_num_rows(fitsfile,nrows,status)
    if( ncols /= size(labels) ) stop 'ncols /= size(labels)'

    ! find columns by labels
    do i = 1, size(labels)
       call fits_get_colnum(fitsfile,.true.,labels(i),col(i),status)
    end do
    if( status /= 0 ) goto 666

    ! get units of PM_*
    call fits_read_key(fitsfile,'TUNIT3',runits,status)
    call fits_read_key(fitsfile,'TUNIT4',dunits,status)

    allocate(fits%coo(nrows,2),fits%pm(nrows,2),stat=stat,errmsg=msg)
    if( stat /= 0 ) then
       write(error_unit,*) 'Error: ',trim(msg)
       error stop 'Failed to allocate data memory'
    end if

    call fits_get_rowsize(fitsfile,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fitsfile,col(1),frow,null,fits%coo(i:l,1),anyf,status)
       call fits_read_col(fitsfile,col(2),frow,null,fits%coo(i:l,2),anyf,status)
       call fits_read_col(fitsfile,col(3),frow,null,fits%pm(i:l,1),anyf,status)
       call fits_read_col(fitsfile,col(4),frow,null,fits%pm(i:l,2),anyf,status)
       if( status /= 0 ) goto 666
    end do

    if( runits == 'arcsec/year' .and. dunits == 'arcsec/year')  then
       fits%pm = fits%pm  / xpm
    else if( runits == 'deg/day' .and. dunits == 'deg/day')  then
       fits%pm = fits%pm
    else
       stop 'Error: unrecognised proper motion units (allowed: arcsec/year, deg/day).'
    end if

666 continue

    call fits_close_file(fitsfile,status)
    call fits_report_error(error_unit,status)

    fits%status = status == 0

  end subroutine Catalogue

  subroutine getcurrent(this,jd,x,y)

    class(CatFits) :: this
    real(dbl), intent(in) :: jd
    real(dbl), dimension(:), intent(out) :: x,y

    real(dbl) :: dt

    if( this%status ) then

       dt = (jd - this%jd)!/365.25_dbl
       x = this%coo(:,1) + dt * this%pm(:,1)
       y = this%coo(:,2) + dt * this%pm(:,2)

    end if

  end subroutine getcurrent

end module fitscat
