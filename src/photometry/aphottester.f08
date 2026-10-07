!
!  aperture photometry tester
!
!  Copyright © 2019-20 F.Hroch (hroch@physics.muni.cz)
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

program testaphot

  use fitsaphot
  use daofotometr

  implicit none

  integer :: naper, i
  real :: r0, e, f, g, g2
  real, dimension(:), allocatable :: raper, apsum
  character(len=666) :: arg

  call get_command_argument(1,arg)
  read(arg,*) r0
  call get_command_argument(2,arg)
  read(arg,*) e

  call rfits('art.fits',naper,raper,apsum)

  do i = 1, naper
     f = 1 - exp(-raper(i)**2/2/r0**2)
     g = apsum(i) / apsum(naper)
     g2 = g2d(raper(i),r0,e)
     write(*,'(6f10.3)') raper(i),g,f,g-f,g2,g2-f
!          grow(raper(i)/r0)*grow(raper(i)/r0/sqrt(1-e**2))
  end do

  deallocate(raper,apsum)

contains

  real function grow(r)

    real, intent(in) :: r

    grow = (1 + erf(r/sqrt(2.0))) / 2

  end function grow

  real function g2d(raper,r0,e)

    real, intent(in) :: raper,r0,e

    integer, parameter :: ndim = 2500
    !    real, dimension(-ndim:ndim,-ndim,ndim) :: g
    integer :: n,i,j
    real :: s,c,g,r2

    n = 0
    s = 0
    c = 100
    do i = -ndim, ndim
       do j = -ndim, ndim
          r2 = (i**2 + j**2) / c**2
!          if( .true. ) then
          if( r2 <= raper**2 ) then
             g = exp(-r2/r0**2/2)
             s = s + g
             n = n + 1
          end if
       end do
    end do
!    write(*,*) s,n,s / n / c,s/250591.047
    g2d = s / n / c
    g2d = s/250591.047
    g2d = s / (c*ndim)

  end function g2d



  subroutine rfits(filename,naper,raper,apsum)

    use titsio

    character(len=*), intent(in) :: filename
    integer, intent(out) :: naper
    real, dimension(:), allocatable, intent(out) :: raper, apsum

    character(len=80) :: key
    real, dimension(:), allocatable :: col
    integer :: status,nrows,irow,i
    type(fitsfiles) :: fits
    logical :: anyf

    irow = 1

    status = 0
    call fits_open_table(fits,filename,FITS_READONLY,status)
    call fits_get_num_rows(fits,nrows,status)
    allocate(col(nrows))

    call fits_read_key(fits,'NAPER',naper,status)
    allocate(raper(naper),apsum(naper))
    do i = 1, naper
       write(key,'(a,i0)') 'APER',i
       call fits_read_key(fits,key,raper(i),status)
    end do

    do i = 1, naper
       call fits_read_col(fits,3+2*i,1,0.0,col,anyf,status)
       apsum(i) = col(irow)
    end do
    deallocate(col)
    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine rfits

end program testaphot
