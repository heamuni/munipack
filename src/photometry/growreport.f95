!
!  Make a report of growth-curves processing
!
!  Copyright © 2016 F.Hroch (hroch@physics.muni.cz)
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
!

module grow_report

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)

  type grow_reporter

     integer :: naper, nstars

     real(dbl), dimension(:,:), allocatable :: grow, grow_err, resgrow
     real(dbl), dimension(:), allocatable :: radius, curve, prof, rapers

  end type grow_reporter


contains


  subroutine grow_report_dump(fits,r,status)

    use titsio

    type(fitsfiles) :: fits
    type(grow_reporter), intent(in) :: r
    integer, intent(in out) :: status

    integer, parameter :: group = 1, extver = 0, frow = 1
    integer :: hdutype,n,i
    character(len=FLEN_VALUE), dimension(:), allocatable :: ttype, tform, tunit
    character(len=FLEN_VALUE) :: key

    if( status /= 0 ) return

    ! HDU for growth-curve (function)
    call fits_movnam_hdu(fits,FITS_ANY_HDU,GROWFUNCEXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       status = 0
    else
       ! already presented ? remove it !
       call fits_delete_hdu(fits,hdutype,status)
       if( status /= 0 ) goto 666
    end if

    n = 3
    allocate(ttype(n), tform(n), tunit(n))

    tform = '1E'
    tunit = ''
    ttype(1) = FITS_COL_R
    ttype(2) = FITS_COL_GROW
    ttype(3) = FITS_COL_RPROF

    call fits_insert_btbl(fits,0,ttype,tform,tunit,GROWFUNCEXTNAME,status)

    call fits_write_col(fits,1,frow,r%radius,status)
    call fits_write_col(fits,2,frow,r%curve,status)
    call fits_write_col(fits,3,frow,r%prof,status)
    deallocate(ttype,tform,tunit)
    if( status /= 0 ) goto 666

    ! HDU for growth-data (empirical grow + residuals) of individual stars
    call fits_movnam_hdu(fits,FITS_ANY_HDU,GROWDATEXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       status = 0
    else
       ! already presented ? remove it !
       call fits_delete_hdu(fits,hdutype,status)
       if( status /= 0 ) goto 666
    end if

    n = 3 * r%naper
    allocate(ttype(n), tform(n), tunit(n))
    tform = '1E'
    tunit = ''

    do i = 1, r%naper
       write(ttype(i          ),'(a,i0)') trim(FITS_COL_GROW),i
       write(ttype(i+  r%naper),'(a,i0)') trim(FITS_COL_GROWERR),i
       write(ttype(i+2*r%naper),'(a,i0)') trim(FITS_COL_RESGROW),i
    end do

    call fits_insert_btbl(fits,0,ttype,tform,tunit,GROWDATEXTNAME,status)
    call fits_write_key(fits,FITS_KEY_NAPER,r%naper,'Count of apertures',status)
    do i = 1, r%naper
       call fits_make_keyn(FITS_KEY_APER,i,key,status)
       call fits_write_key(fits,key,r%rapers(i),-5,'[pix] aperture radius',status)
    end do
    do i = 1,r%naper
       call fits_write_col(fits,i,frow,r%grow(:,i),status)
    end do
    do i = 1,r%naper
       call fits_write_col(fits,i+r%naper,frow,r%grow_err(:,i),status)
    end do
    do i = 1,r%naper
       call fits_write_col(fits,i+2*r%naper,frow,r%resgrow(:,i),status)
    end do
    deallocate(ttype,tform,tunit)

666 continue

  end subroutine grow_report_dump

  subroutine grow_report_init(r,naper,nstars,rapers)

    integer, parameter :: nr = 10
    type(grow_reporter), intent(out) :: r
    integer, intent(in) :: naper,nstars
    real, dimension(:), intent(in) :: rapers

    integer :: nc,n
    real :: rmax

    rmax = rapers(naper)
    n = nstars*naper
    nc = nint(nr*rmax) + 1

    allocate(r%radius(nc),r%curve(nc),r%prof(nc),r%grow(nstars,naper),&
         r%grow_err(nstars,naper),r%resgrow(nstars,naper),r%rapers(naper))

    do n = 1,nc
       r%radius(n) = real(n-1) / real(nr)
    end do
    r%naper = naper
    r%nstars = nstars
    r%rapers = rapers

  end subroutine grow_report_init

  subroutine grow_report_terminate(r)

    type(grow_reporter), intent(in out) :: r

    deallocate(r%radius,r%curve,r%prof,r%grow,r%grow_err,r%resgrow,r%rapers)

  end subroutine grow_report_terminate

end module grow_report
