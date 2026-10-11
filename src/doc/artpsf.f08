!
! gfortran -Wall -g -p -no-pie -fcheck=all cfitsio.f08 fitsio.f08 \
!	 -o artpsf artpsf.f08 -lcfitsio -lm

!  Copyright © 2022 F.Hroch (hroch@physics.muni.cz)
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

program artpsf

  use fitsio
  use iso_fortran_env

  implicit none
  integer, parameter :: width = 16
  integer, parameter :: height = width
  integer, parameter :: i0 = width / 2
  integer, parameter :: j0 = height / 2
  integer, parameter :: zoom = 4
  integer :: status, nelements, i, j
  real, dimension(width,height) :: psf
  type(fitsfiles) :: fits
  real(REAL64) :: s
  real, parameter :: sigma = 1

  ! generate gaussian
  do i = 1, width
     do j = 1, height
        psf(i,j) = exp(-((i - i0)**2 + (j - j0)**2)/sigma**2)
     end do
  end do

  ! normalisation: the sum of samples is zoom**2 (the unit flux),
  ! the centre of the star is the pixel (width/2, height/2)
  s = sum(psf)
  psf = real(psf / s * zoom**2)

  status = 0
  nelements = width*height
  call fits_create_file(fits,'!artpsf.fits',status)
  call fits_insert_img(fits,-32,2,[width,height],status)
  call fits_write_key(fits,'ZOOM',zoom,'',status)
  call fits_write_img(fits,1,1,nelements,pack(psf,.true.),status)
  call fits_close_file(fits,status)
  call fits_report_error(error_unit,status)

end program artpsf
