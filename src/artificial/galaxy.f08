!
!  Galaxies
!
!  Copyright © 2021-2026 F.Hroch (hroch@physics.muni.cz)
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


module artgalaxy

  use iso_fortran_env

  implicit none
  private

  type :: galtype
     character(len=80) :: type
     real(REAL64) :: alpha,delta,Reff,mueff,PA,e,sersic
     ! mueff .. effective surface magnitude at Reff in [mag/arcsec2],
     !          determine effective ieff in counts
     ! Reff .. effective radius [deg]
     ! PA ... position angle [deg, o'clock orientation]
     ! e .. eccentricity
     ! sersic .. Sersic exponent in R**(1/n)
  end type galtype

  type, public :: ArtGalaxies
    type(galtype), dimension(:), allocatable :: gs
  contains
    procedure :: Render => galaxy
    procedure :: Load => LoadGalaxy
  end type ArtGalaxies

  public :: Galaxy, LoadGalaxy

contains

  subroutine galaxy(this,prjcode,crval,crpix,scale,angle,exptime,area,qe, &
       phsys,filter,sky)

    use phsysfits
    use photoconv
    use astrotrafo

    class(ArtGalaxies) :: this
    character(len=*), intent(in) :: prjcode,filter
    real(REAL64), dimension(:), intent(in) :: crval,crpix
    real(REAL64), intent(in) :: scale,angle,exptime,area,qe
    real(REAL64), dimension(:,:), intent(in out) :: sky

    type(AstroWCS) :: wcs
    type(type_phsys), intent(in) :: phsys
    real(REAL64), parameter :: refl = 1
    integer :: i

    call wcs % init(prjcode,crval,crpix,scale,angle,refl)

    do i = 1, size(this%gs)

       if( this%gs(i)%type == 'ELLIPTICAL' ) then
          call elliptical(this%gs(i),wcs,phsys,scale,exptime,area,qe,filter,sky)
       end if
    end do

  end subroutine galaxy

  subroutine elliptical(g,wcs,phsys,scale,exptime,area,qe,filter,sky)

    use phsysfits
    use photoconv
    use astrotrafo

    type(galtype), intent(in) :: g
    type(AstroWCS), intent(in) :: wcs
    type(type_phsys), intent(in) :: phsys
    character(len=*), intent(in) :: filter
    real(REAL64), intent(in) :: scale,exptime,area,qe
    real(REAL64), dimension(:,:), intent(in out) :: sky

    real(REAL64), parameter :: rad = 45 / atan(1.0_REAL64)
    real(REAL64) :: x,y,xc,yc,c,s,u,v,a,b,reff,n,beff,r,q,ceff, dmag
    real(REAL64), dimension(1) :: ct, dct
    integer :: i,j

    call wcs % cel2pix(g%alpha,g%delta,xc,yc)
    reff = g%reff * scale
    c = cos(g%PA / rad)
    s = sin(g%PA / rad)
    a = 1
    b = a*sqrt(1 - g%E**2)

    dmag = 9.999
    call phsysmagph1(phsys,filter,[g%mueff],[dmag],ct,dct)
    ceff = ct(1)
    ceff = ceff * exptime * area * qe * (3600/scale)**2

    n = g%sersic
    beff = 2*n - 1.0/3.0 + (9.8765e-3 + (1.8e-3 + (1.14e-4 - 7e-5/n)/n)/n)/n

    do i = 1, size(sky,1)
       do j = 1, size(sky,2)
          x = i - xc
          y = j - yc
          u = a*(c*x - s*y)
          v = b*(s*x + c*y)
          r = hypot(u,v)
          q = ceff * exp(-beff*((r/reff)**(1.0/n) - 1))
          sky(i,j) = sky(i,j) + q
       end do
    end do

  end subroutine elliptical

  subroutine LoadGalaxy(this,filename)

    use titsio

    class(ArtGalaxies) :: this
    character(len=*), intent(in) :: filename

    real(REAL64), parameter :: nullval = 99.99999

    character(len=FLEN_VALUE) :: extname
    character(len=FLEN_KEYWORD), dimension(7) :: labels
    real(REAL64), dimension(:), allocatable :: alpha,delta,mueff,reff,pa,e,sersic
    integer :: nrows, srows, status, i, l, frow
    integer, dimension(size(labels)) :: cols
    logical :: anyf
    type(fitsfiles) :: fits

    labels = [ character(len=FLEN_KEYWORD) :: &
         'RAJ2000', 'DEJ2000', 'MUEFF', 'REFF', 'PA', 'E', 'SERSIC' ]
    status = 0

    ! open and move to a table extension
    call fits_open_table(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: unable to read a table in the file `',trim(filename),"'."
       stop 'load_galaxy'
    end if

    call fits_get_num_rows(fits,nrows,status)
    if( status /= 0 ) goto 666
    if( nrows == 0 ) stop 'load_galaxy: empty table.'

    ! define reference frame and identification of catalogue
    call fits_read_key(fits,'EXTNAME',extname,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       extname = ''
       status = 0
    end if

    if( extname /= 'ELLIPTICAL' ) stop 'load_galaxy: unknown type galaxy requested.'

    ! find columns by labels
    do i = 1, size(labels)
       call fits_get_colnum(fits,.true.,labels(i),cols(i),status)
    end do
    if( status /= 0 ) goto 666

    allocate(alpha(nrows),delta(nrows),mueff(nrows),reff(nrows),pa(nrows),&
         e(nrows),sersic(nrows))

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,cols(1),frow,nullval,alpha(i:l),anyf,status)
       call fits_read_col(fits,cols(2),frow,nullval,delta(i:l),anyf,status)
       call fits_read_col(fits,cols(3),frow,nullval,mueff(i:l),anyf,status)
       call fits_read_col(fits,cols(4),frow,nullval,reff(i:l),anyf,status)
       call fits_read_col(fits,cols(5),frow,nullval,pa(i:l),anyf,status)
       call fits_read_col(fits,cols(6),frow,nullval,e(i:l),anyf,status)
       call fits_read_col(fits,cols(7),frow,nullval,sersic(i:l),anyf,status)
       if( status /= 0 ) goto 666
    end do

    call fits_close_file(fits,status)

    allocate(this%gs(nrows))
    do i = 1, nrows
       this%gs(i) = galtype(extname,alpha(i),delta(i),reff(i),mueff(i),&
            pa(i),e(i),sersic(i))
    end do
    return

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    stop 'load_galaxy'

  end subroutine LoadGalaxy

end module artgalaxy
