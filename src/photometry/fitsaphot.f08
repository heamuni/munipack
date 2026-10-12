!
!  fitsaphot
!
!  Copyright © 2013-6, 2018-2024 F.Hroch (hroch@physics.muni.cz)
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

! Parameters by:
! http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?psfmeasure.hlp


module fitsaphot

  use titsio
  use iso_fortran_env

  implicit none

contains

  subroutine fits_aphot_image(filename,fkeys,data,stderr,saturate,status)

    integer, parameter :: naxis = 2

    character(len=*),intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: fkeys
    real, dimension(:,:), allocatable, intent(out) :: data,stderr
    real, intent(out) :: saturate
    integer, intent(in out) :: status

    type(fitsfiles) :: fits
    integer, parameter :: extver = 0
    integer, dimension(naxis) :: naxes
    integer :: bitpix
    logical :: anyf

    if( status /= 0 ) return

    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read an image in the file `',trim(filename),"'."
       return
    end if

    call fits_get_img_size(fits,naxis,naxes,status)
    if( status /= 0 ) goto 666

    call fits_read_key(fits,fkeys(1),saturate,status)
    if( status == FITS_KEYWORD_NOT_FOUND ) then
       status = 0
       call fits_get_img_type(fits,bitpix,status)
       if( status /= 0 ) goto  666
       if( bitpix > 0 ) then
          saturate = 2.0**bitpix - 1
       else
          saturate = huge(saturate)
       end if
    end if

    allocate(data(naxes(1),naxes(2)),stderr(naxes(1),naxes(2)))
    call fits_read_image(fits,0,0.0,data,anyf,status)

    call fits_movnam_hdu(fits,FITS_IMAGE_HDU,EXT_STDERR,extver,status)
    if( status == 0 ) then
       call fits_read_image(fits,0,0.0,stderr,anyf,status)
    else if ( status == FITS_BAD_HDU_NUM ) then
       ! if the information about standard errors is not available,
       ! we are continuing with the Poisson component only
       where( data > 0 )
          stderr = sqrt(data)
       elsewhere
          stderr = -1
       end where
       status = 0
    end if

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( status /= 0 ) then
       if( allocated(data) ) deallocate(data,stderr)
    end if

  end subroutine fits_aphot_image


  subroutine fits_aphot_read(filename,data,stderr,xcens,ycens, &
       lobad,hibad,fwhm,ecc,incl,status)

    integer, parameter :: DIM = 2

    character(len=*),intent(in) :: filename
    real, dimension(:,:), allocatable, intent(out) :: data,stderr
    real, dimension(:), allocatable, intent(out) :: xcens,ycens
    real, intent(out) :: lobad,hibad,fwhm,ecc,incl
    integer, intent(in out) :: status

    type(fitsfiles) :: fits
    integer, parameter :: extver = 0
    integer, dimension(DIM) :: naxes
    integer :: naxis,nrows,srows,frow,i,l,xcol,ycol
    logical :: anyf

    if( status /= 0 ) return

    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) &
            'Error: failed to read an image in the file `',trim(filename),"'."
       return
    end if

    call fits_get_img_dim(fits,naxis,status)
    if( naxis /= 2 .and. status == 0 ) then
       write(error_unit,*) 'Error in aphot: Only 2D frames are supported.'
       goto 666
    end if
    call fits_get_img_size(fits,DIM,naxes,status)
    if( status /= 0 ) goto 666

    allocate(data(naxes(1),naxes(2)),stderr(naxes(1),naxes(2)))
    call fits_read_image(fits,0,0.0,data,anyf,status)

    call fits_movnam_hdu(fits,FITS_IMAGE_HDU,EXT_STDERR,extver,status)
    if( status == 0 ) then
       call fits_read_image(fits,0,0.0,stderr,anyf,status)
    else if ( status == FITS_BAD_HDU_NUM ) then
       ! if the information about standard errors is not available,
       ! we continues with Poisson component only
       where( data > 0 )
          stderr = sqrt(data)
       elsewhere
          stderr = -1
       end where
       status = 0
    end if

    if( status /= 0 ) then
       write(error_unit,*) trim(filename),": Failed to read data."
       goto 666
    end if

    call fits_movnam_hdu(fits,FITS_BINARY_TBL,FINDEXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       write(error_unit,*) "Error: ",trim(FINDEXTNAME)//" extension not found."
       write(error_unit,*) "       Has been stars detected by `munipack find " &
            //trim(filename)//"' ?"
       goto 666
    end if

    call fits_get_num_rows(fits,nrows,status)

    call fits_read_key(fits,FITS_KEY_LOWBAD,lobad,status)
    call fits_read_key(fits,FITS_KEY_HIGHBAD,hibad,status)
    call fits_read_key(fits,FITS_KEY_FWHM,fwhm,status)
    call fits_read_key(fits,FITS_KEY_ECCENTRICITY,ecc,status)
    call fits_read_key(fits,FITS_KEY_INCLINATION,incl,status)
    if( status /= 0 ) then
       write(error_unit,*) trim(filename),": Required keywords ", &
            trim(FITS_KEY_LOWBAD),", ",trim(FITS_KEY_HIGHBAD),", ", &
            trim(FITS_KEY_FWHM),", ",trim(FITS_KEY_ECCENTRICITY),", or ",&
            trim(FITS_KEY_INCLINATION)," not found."
       goto 666
    end if

    call fits_get_colnum(fits,.true.,FITS_COL_X,xcol,status)
    call fits_get_colnum(fits,.true.,FITS_COL_Y,ycol,status)

    allocate(xcens(nrows),ycens(nrows))

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_read_col(fits,xcol,frow,0.0,xcens(i:l),anyf,status)
       call fits_read_col(fits,ycol,frow,0.0,ycens(i:l),anyf,status)
    end do

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

    if( status /= 0 ) then
       if( allocated(data) ) deallocate(data,stderr)
       if( allocated(xcens) ) deallocate(xcens,ycens)
    end if

  end subroutine fits_aphot_read

  subroutine fits_aphot_info(filename,naxes,scale,gain,status)

    ! the frame size, the scale [pix/deg] by WCS (-1 when undefined)
    ! and the gain [e-/ADU] (1 when undefined)

    character(len=*),intent(in) :: filename
    integer, dimension(:), intent(out) :: naxes
    real(REAL64), intent(out) :: scale, gain
    integer, intent(in out) :: status

    type(fitsfiles) :: fits
    real(REAL64), dimension(2,2) :: cd
    real(REAL64), dimension(2) :: cdelt

    if( status /= 0 ) return

    call fits_open_image(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) return
    call fits_get_img_size(fits,2,naxes,status)
    if( status /= 0 ) goto 666

    call fits_write_errmark
    call fits_read_key(fits,'CD1_1',cd(1,1),status)
    call fits_read_key(fits,'CD1_2',cd(1,2),status)
    call fits_read_key(fits,'CD2_1',cd(2,1),status)
    call fits_read_key(fits,'CD2_2',cd(2,2),status)
    if( status == 0 ) then
       scale = 1 / sqrt(abs(cd(1,1)*cd(2,2) - cd(1,2)*cd(2,1)))
    else
       call fits_clear_errmark
       status = 0
       call fits_write_errmark
       call fits_read_key(fits,'CDELT1',cdelt(1),status)
       call fits_read_key(fits,'CDELT2',cdelt(2),status)
       if( status == 0 ) then
          scale = 1 / sqrt(abs(cdelt(1)*cdelt(2)))
       else
          call fits_clear_errmark
          status = 0
          scale = -1
       end if
    end if

    call fits_write_errmark
    call fits_read_key(fits,FITS_KEY_GAIN,gain,status)
    if( status /= 0 .or. .not. (gain > 0) ) then
       call fits_clear_errmark
       status = 0
       gain = 1
    end if

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine fits_aphot_info

  subroutine fits_aphot_save(filename, output, hwhm, ecc, incl, raper, ring,  &
       xcens, ycens, apcts,apcts_err,sky,sky_err, errgain, status, ee, eeerr)

    ! The optional ee, eeerr (the enclosed energy of apertures and its
    ! error) are saved as the columns APEEi, APEEERRi following
    ! the original ones. errgain is the gain [e-/ADU] of the Poisson
    ! noise of the stars in APCOUNTERR (the keyword APERRGN).

    character(len=*), intent(in) :: filename, output
    real, intent(in) :: hwhm, ecc, incl
    real, dimension(:), intent(in) :: raper, ring, xcens, ycens, sky,sky_err
    real, dimension(:,:), intent(in) :: apcts,apcts_err
    real(REAL64), intent(in) :: errgain
    integer, intent(in out) :: status
    real(REAL64), dimension(:,:), intent(in), optional :: ee, eeerr

    type(fitsfiles) :: fits
    integer, parameter :: extver = 0, nbegin = 4
    character(len=FLEN_VALUE), dimension(:), allocatable :: ttype, tform, tunit
    character(len=FLEN_VALUE) :: key
    integer :: hdutype,i,j,k,l,n,nrows,srows,frow,nee
    logical :: have_ee


    if( status /= 0 ) return

    if( output == '' ) then
       call fits_open_file(fits,filename,FITS_READWRITE,status)
       if( status /= 0 ) then
          write(error_unit,*) &
            'Error: failed to open the file `',trim(filename),"' for a table update."
          return
       end if
    else
       call fits_precopy_file(fits,filename,output,FITS_READWRITE,.true.,status)
       if( status /= 0 ) then
          write(error_unit,*) &
            'Error: failed to create the file `',trim(output),"'."
          return
       end if
    end if

    ! store results to the aperture photometry extension
    call fits_movnam_hdu(fits,FITS_BINARY_TBL,APEREXTNAME,extver,status)
    if( status == FITS_BAD_HDU_NUM ) then
       status = 0
    else
       ! already presented ? remove it !
       call fits_delete_hdu(fits,hdutype,status)
    end if
    if( status /= 0 ) goto 666

    have_ee = present(ee) .and. present(eeerr)
    nee = nbegin + 2*size(raper)
    if( have_ee ) then
       n = nbegin + 4*size(raper)
    else
       n = nbegin + 2*size(raper)
    end if
    allocate(ttype(n), tform(n), tunit(n))

    tform = '1D'
    tunit = ''
    ttype(1) = FITS_COL_X
    ttype(2) = FITS_COL_Y
    ttype(3) = FITS_COL_SKY
    ttype(4) = FITS_COL_SKYERR

    do i = 1, size(raper)
       j = nbegin - 1 + 2*i
       write(ttype(j),'(a,i0)') trim(FITS_COL_APCOUNT),i
       write(ttype(j+1),'(a,i0)') trim(FITS_COL_APCOUNTERR),i
    end do
    if( have_ee ) then
       do i = 1, size(raper)
          j = nee - 1 + 2*i
          write(ttype(j),'(a,i0)') trim(FITS_COL_APEE),i
          write(ttype(j+1),'(a,i0)') trim(FITS_COL_APEEERR),i
       end do
    end if

    ! aperture photometry table
    nrows = size(xcens)
    call fits_insert_btbl(fits,nrows,ttype,tform,tunit,APEREXTNAME,status)

    call fits_write_key(fits,FITS_KEY_HWHM,hwhm,-4, &
         '[pix] half width at half of maximum',status)
    call fits_update_key(fits,FITS_KEY_ECCENTRICITY,ecc,-2, &
         ' eccentricity',status)
    call fits_update_key(fits,FITS_KEY_INCLINATION,nint(incl), &
         ' inclination',status)
    call fits_write_key(fits,FITS_KEY_NAPER,size(raper), &
         'Count of apertures',status)
    do i = 1, size(raper)
       call fits_make_keyn(FITS_KEY_APER,i,key,status)
       call fits_write_key(fits,key,raper(i),-5,'[pix] aperture radius',status)
    end do

    call fits_write_key(fits,trim(FITS_KEY_ANNULUS)//'1',ring(1),-5, &
         '[pix] inner sky annulus radius',status)
    call fits_write_key(fits,trim(FITS_KEY_ANNULUS)//'2',ring(2),-5, &
         '[pix] outer sky annulus radius',status)
    ! (a free-format value, so the whole comment fits the card)
    block
      character(len=FLEN_VALUE) :: val
      integer :: ie, iz
      write(val,'(g0.12)') errgain
      ie = scan(val,'Ee')
      if( ie == 0 ) ie = len_trim(val) + 1
      iz = ie - 1
      do while( iz > 1 .and. val(iz:iz) == '0' )
         iz = iz - 1
      end do
      if( val(iz:iz) == '.' ) iz = iz + 1
      val = val(1:iz)//val(ie:)
      call fits_update_card(fits,'APERRGN','APERRGN = '//trim(val)// &
           ' / [e-/ADU] gain applied to the star term of APCOUNTERR',status)
    end block

    call fits_get_rowsize(fits,srows,status)
    do i = 1, nrows, srows
       frow = i
       l = min(i+srows,nrows)

       call fits_write_col(fits,1,frow,xcens(i:l),status)
       call fits_write_col(fits,2,frow,ycens(i:l),status)
       call fits_write_col(fits,3,frow,sky(i:l),status)
       call fits_write_col(fits,4,frow,sky_err(i:l),status)

       do k = 1,size(apcts,2)
          j = nbegin-1+2*k
          call fits_write_col(fits,j,frow,apcts(i:l,k),status)
          call fits_write_col(fits,j+1,frow,apcts_err(i:l,k),status)
       end do
       if( have_ee ) then
          do k = 1,size(ee,2)
             j = nee-1+2*k
             call fits_write_col(fits,j,frow,ee(i:l,k),status)
             call fits_write_col(fits,j+1,frow,eeerr(i:l,k),status)
          end do
       end if
    end do
    deallocate(ttype,tform,tunit)

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine fits_aphot_save

  subroutine estim_hwhm(data,xcens,ycens,sky,fwhm,lobad,hibad,hwhm)

    use oakleaf

    real, dimension(:,:), intent(in) :: data
    real, dimension(:), intent(in) :: xcens,ycens,sky
    real, intent(in) :: fwhm,lobad,hibad
    real, intent(out) :: hwhm

    ! estimate of HWHM (half width at half of maximum) of stars:
    ! the second moments of the brightest stars give the Gaussian sigma,
    ! HWHM = sqrt(2*ln(2)) * sigma

    real, dimension(:), allocatable :: xhwhm, w0
    real :: sx,sy,w,sw,sigma
    integer :: nx,ny,i,j,l,m,n,nmin,i0,j0

    allocate(xhwhm(size(xcens)),w0(13))

    nx = size(data,1)
    ny = size(data,2)

    n = 0
    m = nint(1.5*fwhm)
    m = nint(fwhm / 2) * 3
    ! neighborhood is 3*hwhm of expected which prefers important
    ! parts of profile

    do l = 1, size(xcens)

       sx = 0
       sy = 0
       sw = 0
       i0 = nint(xcens(l))
       j0 = nint(ycens(l))
       if( data(i0,j0) > sky(l) .and. sky(l) > 0 ) then
          do i = i0-m,i0+m
             do j = j0-m,j0+m
                if( 0 < i .and. i <= nx .and. 0 < j .and. j <= ny ) then
                   w = data(i,j) - sky(l)
                   if( lobad < data(i,j) .and. data(i,j) < hibad .and. w > 0) then
                      sx = sx + w*(i - xcens(l))**2
                      sy = sy + w*(j - ycens(l))**2
                      sw = sw + w
                   end if
                end if
             end do
          end do
          if( sw > 0 ) then
             ! estimation of hwhm is sensitive on noise in data (w),
             ! we're selecting only the brighest stars
             if( n < size(w0) ) then
                n = n + 1
                xhwhm(n) = (sqrt(sx/sw) + sqrt(sy/sw)) / 2
                w0(n) = data(i0,j0) - sky(l)
             else
                nmin = minloc(w0,1)
                if( data(i0,j0) - sky(l) > w0(nmin) ) then
                   w0(nmin) = data(i0,j0) - sky(l)
                   xhwhm(nmin) = (sqrt(sx/sw) + sqrt(sy/sw)) / 2
                end if
             end if
          end if
       end if
    end do

    if( n > 1 ) then
       call rmean(xhwhm(1:n),sigma,w)
       hwhm = sqrt(2*log(2.0)) * sigma
    else
       hwhm = -1
    end if

    deallocate(xhwhm)

  end subroutine estim_hwhm

end module fitsaphot
