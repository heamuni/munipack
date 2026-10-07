!
! Common for read/write photometry-related FITSes
!
!
! Copyright © 2013,2015-6,2018-20 F.Hroch (hroch@physics.muni.cz)
!
! This file is part of Munipack.
!
! Munipack is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Munipack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
!

module phio

  use titsio
  use astrotrafo

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: rad = 57.295779513082322865_dbl

contains


  subroutine wcsget(fits,wcs,status)

    type(fitsfiles), intent(in) :: fits
    integer, intent(in out) :: status
!    type(AstroTrafoProj), intent(out) :: t
    type(AstroWCS), intent(out) :: wcs

    integer, parameter :: DIM = 2
    integer, dimension(DIM) :: naxes
    character(len=FLEN_VALUE),dimension(2) :: ctype
    real(dbl), dimension(2,2) :: cd
    real(dbl), dimension(2) :: crval,crpix,crder

    ! read image dimensions
    call fits_get_img_size(fits,DIM,naxes,status)

    ! read astrometric calibration
    call fits_read_key(fits,'CTYPE1',ctype(1),status)
    call fits_read_key(fits,'CTYPE2',ctype(2),status)
    call fits_read_key(fits,'CRVAL1',crval(1),status)
    call fits_read_key(fits,'CRVAL2',crval(2),status)
    call fits_read_key(fits,'CRPIX1',crpix(1),status)
    call fits_read_key(fits,'CRPIX2',crpix(2),status)
    call fits_read_key(fits,'CD1_1',cd(1,1),status)
    call fits_read_key(fits,'CD1_2',cd(1,2),status)
    call fits_read_key(fits,'CD2_1',cd(2,1),status)
    call fits_read_key(fits,'CD2_2',cd(2,2),status)
    if( status /= 0 ) return

    ! optional keywords
    call fits_read_key(fits,'CRDER1',crder(1),status)
    call fits_read_key(fits,'CRDER2',crder(2),status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       crder = 0
    end if

    if( status /= 0 ) return

    block
      character(len=FLEN_CARD), dimension(0) :: exclist
      character(len=:), allocatable :: header
      integer :: nkeys, status

      call fits_hdr2str(fits,1,exclist,nkeys,header,status)
      call wcs % init(header,nkeys,status)
      deallocate(header)

!      call wcs % get_affines(ctype,crval,crpix
    end block

    ! interpret the calibration
!    call trafo_fromwcs(t,ctype,crval,crpix,cd,crder)

  end subroutine wcsget

!!$  subroutine trafo_fromwcs(t,ctype,crval,crpix,cd,crder)
!!$
!!$    type(AstroTrafoProj), intent(out) :: t
!!$    character(len=*), dimension(:), intent(in) :: ctype
!!$    real(dbl), dimension(:),intent(in) :: crval,crpix,crder
!!$    real(dbl), dimension(:,:),intent(in) :: cd
!!$    character(len=FLEN_VALUE) :: type
!!$    real(dbl) :: c,s,sc,rot,refl,err
!!$    real(dbl), dimension(2,2) :: mi,m
!!$
!!$    if( size(ctype) /= 2 ) stop 'Bad dimensions of CTYPE (no 2).'
!!$
!!$    if( index(ctype(1),"-TAN") > 0 ) then
!!$       type = 'GNOMONIC'
!!$    else
!!$       type = ' '
!!$    end if
!!$
!!$    ! un-flip spherical to rectangular coordinates
!!$    mi(1,:) = [ -1.0_dbl, 0.0_dbl ]
!!$    mi(2,:) = [  0.0_dbl, 1.0_dbl ]
!!$    m = matmul(mi,cd)
!!$
!!$    ! the reflection (flip, mirror) is determined from
!!$    ! diagonal elements with (non-)corresponding signs
!!$    refl = sign(1.0_dbl,m(1,1)*m(2,2))
!!$
!!$    ! transform-out reflex
!!$    mi(1,:) = [    refl, 0.0_dbl ]
!!$    mi(2,:) = [ 0.0_dbl, 1.0_dbl ]
!!$    m = matmul(m,mi)
!!$
!!$    c = (m(1,1) + m(2,2)) / 2.0_dbl
!!$    s = (m(2,1) - m(1,2)) / 2.0_dbl
!!$    sc = 1.0_dbl / sqrt(c**2 + s**2)   ! in pix per deg
!!$    rot = rad*(atan2(-m(1,2),m(1,1)) + atan2(m(2,1),m(2,2)))/2.0_dbl
!!$
!!$    if( sum(abs(crder)) > size(crder)*epsilon(1.0_dbl) ) then
!!$       err = sqrt(sum(crder**2)/size(crder))
!!$    else
!!$       err = 1.0/sc
!!$    end if
!!$
!!$    call trafo_init(t,type,crval(1),crval(2),crpix(1),crpix(2), &
!!$         scale=1/sc,rot=rot,refl=refl,err=err)
!!$
!!$  end subroutine trafo_fromwcs


  function fits_jd(dateobs,status) result(jd)

    use trajd

    real(dbl) :: jd
    character(len=*), intent(in) :: dateobs
    integer ::year, month, day, hour, minute
    real(dbl) :: second, d
    integer :: status

    jd = 0
    if( status /= 0 ) return

    if( dateobs /= '' ) then

       ! decode date-time string
       call fits_str2date(dateobs,year,month,day,hour,minute,second,status)
       if( status == 0 ) then

          ! Evaluate Julian date
          d = day + (hour + (minute + second/60.0_dbl) / 60.0_dbl) / 24.0_dbl
          jd = datjd(real(year,dbl),real(month,dbl),d)
       end if

    end if

  end function fits_jd


end module phio
