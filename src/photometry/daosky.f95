!
!  Sky - average sky value
!
!  Copyright © 1999, 2010, 2013, 2016, 2018-9 F.Hroch (hroch@physics.muni.cz)
!  Copyright (C) 1991 P.B. Stetson, Dominon Astrophysical Observatory
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
!====================================================================
!
!  This source is on base of DAOPHOT II package by P.B.Stetson
!
!=======================================================================
!
! This subroutine estimates an average sky value for a picture by taking
! individual pixels scattered over the picture.  The brightness values
! are sorted, and the modal value is estimated using the MMM subroutine.
!
!               OFFICIAL DAO VERSION:  1991 April 18
!
!=======================================================================
!

module mdaosky

contains

subroutine  daosky (d, nstep, verbose, hibad, sky, skyerr, skysig)

  use oakleaf

  implicit none

  integer, intent(in) :: nstep
  logical, intent(in) :: verbose
  real, intent(in) :: hibad
  real, intent(in), dimension(:,:) :: d
  real, intent(out) :: sky, skyerr, skysig

!
! NMAX    is the maximum number of sky pixels we can deal with,
!         given the limited amount of working space.
!

  integer :: ncol, nrow, n
  real, dimension(:), allocatable :: s

  ncol = size(d,1)
  nrow = size(d,2)
!
!-----------------------------------------------------------------------
!
! The spacing between pixels that will be included in the sample is
! estimated by the ratio of the total number of pixels in the picture to
! the maximum number of pixels that can be accomodated in the vector S.
!
!  istep = max((ncol*nrow)/nmax, 1)
!
! Go through the disk file reading a row at a time and extracting every
! ISTEP-th pixel.  If ISTEP is not equal to 1, make sure that the
! starting pixel for each row is staggered.
!

  allocate(s(size(d)))
  s = pack(d(1:ncol:nstep,1:nrow:nstep), d < hibad)
  n = count(d(1:ncol:nstep,1:nrow:nstep) < hibad)

!
! Sort these values, then estimate the mode.
!

  call rmean(s(1:n),sky,skyerr,skysig)

  if( verbose ) &
       write (*,"(A,F9.1,A,F10.2)") &
       ' Approximate sky value for this frame =', sky,' +- ',skyerr

  deallocate(s)

end subroutine daosky

end module mdaosky
