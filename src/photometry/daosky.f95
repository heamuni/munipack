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

subroutine  daosky (d, nstep, verbose, hibad, sky, skyerr, skysig, gain, rn2, bias)

  ! The mean sky level of the frame by every nstep-th pixel (in both
  ! directions) below hibad: the clipped mean of skymean (the robust mean
  ! used before is biased on a sky of a few photons per pixel, BUGS P23).
  ! The noise model (gain [e-/ADU], read noise squared [ADU**2], bias
  ! [ADU]) is optional.

  use iso_fortran_env
  use skymean

  implicit none

  integer, intent(in) :: nstep
  logical, intent(in) :: verbose
  real, intent(in) :: hibad
  real, intent(in), dimension(:,:) :: d
  real, intent(out) :: sky, skyerr, skysig
  real(REAL64), intent(in), optional :: gain, rn2, bias

  real(REAL64), dimension(:), allocatable :: s
  real(REAL64) :: g, rn, b, mean, err, sig
  integer :: ncol, nrow, n

  ncol = size(d,1)
  nrow = size(d,2)

  g = -1
  rn = 0
  b = 0
  if( present(gain) ) g = gain
  if( present(rn2) ) rn = rn2
  if( present(bias) ) b = bias

  ! the mask of the same (sub-sampled) shape as the pixels
  s = pack(real(d(1:ncol:nstep,1:nrow:nstep),REAL64), &
       d(1:ncol:nstep,1:nrow:nstep) < hibad)

  call sky_mean(s,g,rn,b,mean,err,sig,n)
  sky = real(mean)
  skyerr = real(err)
  skysig = real(sig)

  if( verbose ) &
       write (*,"(A,F9.1,A,F10.2)") &
       ' Approximate sky value for this frame =', sky,' +- ',skyerr

end subroutine daosky

end module mdaosky
