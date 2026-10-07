!
!  transformation of an instrumental to standard photometry system
!
!  Copyright © 2013-6 F.Hroch (hroch@physics.muni.cz)
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


program phfotran

  use titsio
  use mfits
  use phsysfits

  implicit none

  integer, parameter :: dbl = selected_real_kind(15)
  character(len=FLEN_KEYWORD), dimension(8) :: keys
  character(len=4*FLEN_FILENAME) :: record, key, val, output
  logical :: verbose = .false., plog = .false., list = .false.
  character(len=FLEN_FILENAME), dimension(:), allocatable :: filenames
  character(len=FLEN_VALUE), dimension(:), allocatable :: col_mag, col_magerr, filters
  character(len=FLEN_FILENAME) :: phsystable = 'photosystems.fits', cat = ''
  character(len=FLEN_VALUE) :: col_ra = FITS_COL_RA, col_dec = FITS_COL_DEC, &
       photsys_instr = '', photsys_ref = ''
  real(dbl) :: utol = -1.0/3600.0
  real(dbl), dimension(:), allocatable :: ctph, extin
  logical :: init_area = .false.
  real(dbl) :: area
  integer :: nfile,n,eq,nmag,nmagerr,nfilters,nctph,nextin

  output = '!phfotran.fits'
  area = 1

  keys(1) = FITS_KEY_EXPTIME
  keys(2) = FITS_KEY_AREA
  keys(3) = FITS_KEY_DATEOBS
  keys(4) = FITS_KEY_PHOTSYS
  keys(5) = FITS_KEY_FILTER
  keys(6) = FITS_KEY_LONGITUDE
  keys(7) = FITS_KEY_LATITUDE
  keys(8) = FITS_KEY_TIMEOBS

  allocate(filenames(0))
  nfile = 0
  allocate(col_mag(0), col_magerr(0), filters(0), ctph(0),extin(0))
  nfilters = 0
  nmag = 0
  nmagerr = 0
  nctph = 0
  nextin = 0

  do
     read(*,'(a)',end=20) record

     eq = index(record,'=')
     if( eq == 0 ) stop 'Improper control data on input.'

     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'PIPELOG' ) then

        read(val,*) plog

     else if( key == 'FITS_KEY_EXPTIME' ) then

        read(val,*) keys(1)

     else if( key == 'FITS_KEY_AREA' ) then

        read(val,*) keys(2)

     else if( key == 'FITS_KEY_DATEOBS' ) then

        read(val,*) keys(3)

     else if( key == 'FITS_KEY_PHOTOSYS' ) then ! remove

        read(val,*) keys(4)

     else if( key == 'FITS_KEY_FILTER' ) then

        read(val,*) keys(5)

     else if( key == 'FITS_KEY_LONGITUDE' ) then

        read(val,*) keys(6)

     else if( key == 'FITS_KEY_LATITUDE' ) then

        read(val,*) keys(7)

     else if( key == 'FITS_KEY_TIMEOBS' ) then

        read(val,*) keys(8)

     else if( key == 'COL_RA' ) then

        read(val,*) col_ra

     else if( key == 'COL_DEC' ) then

        read(val,*) col_dec

     else if( key == 'COL_NMAG' ) then

        read(val,*) nmag

        deallocate(col_mag)
        allocate(col_mag(nmag))

     else if( key == 'COL_MAG' ) then

        ! a set of columns of filters
        read(val,*) col_mag

     else if( key == 'COL_NMAGERR' ) then

        read(val,*) nmagerr

        deallocate(col_magerr)
        allocate(col_magerr(nmagerr))

     else if( key == 'COL_MAGERR' ) then

        read(val,*) col_magerr

     else if( key == 'NFILTERS' ) then

        read(val,*) nfilters

        deallocate(filters)
        allocate(filters(nfilters))

     else if( key == 'FILTERS' ) then

        ! a set of filters
        read(val,*) filters

     else if( key == 'NCTPH' ) then

        read(val,*) nctph

        deallocate(ctph)
        allocate(ctph(nctph))

     else if( key == 'CTPH' ) then

        ! efficiency
        read(val,*) ctph

     else if( key == 'NEXTIN' ) then

        read(val,*) nextin

        deallocate(extin)
        allocate(extin(nextin))

     else if( key == 'EXTIN' ) then

        read(val,*) extin

     else if( key == 'PHOTSYS_INSTR' ) then

        read(val,*) photsys_instr

     else if( key == 'PHOTSYS_REF' ) then

        read(val,*) photsys_ref

     else if( key == 'PHSYSTABLE' ) then

        read(val,*) phsystable

     else if( key == 'LIST' ) then

        read(val,*) list

     else if( key == 'TOL' ) then

        read(val,*) utol

     else if( key == 'AREA' ) then

        read(val,*) area
        init_area = .true.

     else if( key == 'CAT' ) then

        read(val,*) cat

     else if( key == 'OUTPUT' ) then

        read(val,*) output

     else if( key == 'NFILE' ) then

        read(val,*) nfile
        deallocate(filenames)
        allocate(filenames(nfile))
        n = 0

     else if( key == 'FILE' ) then

        n = n + 1
        if( n > size(filenames) ) stop 'Too many files.'

        read(val,*) filenames(n)

     end if

  end do
20 continue

  if( photsys_instr == '' .or. photsys_ref == '' ) &
       stop 'Reference or instrumental photometry system undefined.'

  if( list ) then
     call listphsys(phsystable)
  else if ( nfile == 1 .and. cat == '' ) then
     call table
  else
     call frames
  end if

  deallocate(filenames,col_mag,col_magerr,filters,ctph,extin)

  stop 0

contains


  ! transformation
  subroutine frames

    use fits_fotran
    use fotran
    use jamming

    character(len=FLEN_VALUE) :: catid, photsys_frames
    character(len=FLEN_VALUE), dimension(:), allocatable :: filter
    real(dbl), allocatable, dimension(:) :: ra, dec, airmass
    real(dbl), allocatable, dimension(:,:) :: tr,tr1,trerr,tr1err, ph,dph,cts,dcts
    integer, dimension(:,:), allocatable :: pairs
    real(dbl) :: e
    integer :: k,nph

    if( size(filters) == 0 ) stop 'Filters undefined. Use -f option.'
    if( size(filenames) < 2 ) stop 'One frame only? No more frames - no more love.'
    if( nextin > 0 .and. nextin /= size(filenames) ) &
         stop 'Number of extinction elements does not match number of files.'
    if( nctph > 0 .and. nctph /= size(filenames) ) &
         stop 'Number of ctph elements does not match number of files.'
    if( size(filters) /= size(filenames) ) &
         stop 'Number of filters elements does not match number of files.'

    if( size(filters) > 0 .and. size(col_mag) == 0 ) then
       deallocate(col_mag)
       allocate(col_mag(size(filters)))
       col_mag = filters
    end if

    call jamcatx(cat,phsystable,keys,col_ra,col_dec,col_mag,col_magerr,&
       utol, area, init_area, filters, filenames, catid, filter, ra, dec, airmass, &
       ph,dph, cts,dcts, photsys_ref, photsys_frames, pairs)

    if( size(ra) < 2 ) &
         stop 'Determination of parameters of line requires more than two points.'

    if( photsys_frames /= '' ) photsys_instr = photsys_frames

    nph = nfilters
    allocate(tr(nph,nph),trerr(nph,nph),tr1(nph,nph),tr1err(nph,nph))

    do k = 1,nctph
       cts(:,k) = cts(:,k)*ctph(k)
       dcts(:,k) = dcts(:,k)*ctph(k)
    end do

    if( nextin > 0 .and. verbose ) write(*,*) 'Extinction: filter, airmass, attenuation:'
    do k = 1,nextin
       e = exp(extin(k)*airmass(k))
       cts(:,k) = e*cts(:,k)
       dcts(:,k) = e*dcts(:,k)
       if( verbose ) write(*,'(1x,a,2f7.3)') trim(filters(k)),airmass(k),e
    end do

    ! testing
!    tr(1,:) = (/0.25,-0.03 /)
!    tr(2,:) = (/-0.02,0.5 /)

!    tr(1,:) = (/11.0,0.1,0.0,0.0,0.0 /)
!    tr(2,:) = (/-0.1,3.0,0.1,0.0,0.0 /)
!    tr(3,:) = (/ 0.0,0.2,1.0,-0.1,0.0 /)
!    tr(4,:) = (/ 0.0,0.0,-0.2,2.0,0.2 /)
!    tr(5,:) = (/ 0.0,0.0,0.0,-0.1,3.0 /)

!    do k = 1,size(ph,1)
!       cts(k,:) = matmul(tr,ph(k,:))
!       dcts(k,:) = matmul(tr,dph(k,:))
!    end do

    call trafo(pairs,ph,dph,cts,dcts,tr,trerr,tr1,tr1err,verbose)

    call trawrite(output,photsys_ref,photsys_instr,filters,filters, &
         ra,dec,ci0,ph,dph,cts,dcts,tr,trerr,tr1,tr1err,filenames,airmass,ctph,extin)

    deallocate(ra,dec,airmass,cts,dcts,ph,dph,filter,tr,trerr,tr1,tr1err)

  end subroutine frames


  ! transformation
  subroutine table

    use fits_fotran
    use fotran
    use phsysfits

    real(dbl), allocatable, dimension(:) :: ra, dec
    real(dbl), allocatable, dimension(:,:) :: tr,tr1,trerr,tr1err, ph,dph,cts,dcts
    integer, dimension(:,:), allocatable :: pairs
    type(type_phsys) :: phsys
    integer :: n

    call phselect(phsystable,photsys_ref,phsys)
    call phsyspairs(phsys,filters,pairs)

    call tradata(filenames(1),ra,dec,ph,dph,cts,dcts)

    n = size(ph,2)
    allocate(tr(n,n),trerr(n,n),tr1(n,n),tr1err(n,n))

    call trafo(pairs,ph,dph,cts,dcts,tr,trerr,tr1,tr1err,verbose)

    call trawrite(output,photsys_ref,photsys_instr,filters,filters, &
         ra,dec,ci0,ph,dph,cts,dcts,tr,trerr,tr1,tr1err)

    deallocate(ra,dec,cts,dcts,ph,dph,tr,trerr,tr1,tr1err)

  end subroutine table


  subroutine trafo(pairs,ph,dph,cts,dcts,tr,trerr,tr1,tr1err,verbose)

    use fotran

    integer, dimension(:,:), intent(in) :: pairs
    real(dbl), dimension(:,:), intent(in) :: ph,dph,cts,dcts
    real(dbl), dimension(:,:), intent(out) :: tr,trerr,tr1,tr1err
    logical, intent(in) :: verbose

    real(dbl), dimension(:), allocatable :: ct,dct
    integer :: k,l,nph

    nph = size(ph,2)
    allocate(ct(nph),dct(nph))

    call stdtra(pairs,ph,dph,cts,dcts,tr,trerr,verbose)
    if( verbose ) then
       write(*,*) 'Forward matrix (ph-rate to cts-rate):'
       do k = 1,size(tr,1)
          write(*,'(5f12.4)') tr(k,:)
       end do
    end if


    if( verbose ) then
       write(*,*) 'Relative residuals in filters:'
       write(*,'(tr5,5a10)') filters
       do l = 1,size(ph,1)
          call fotra(tr,pairs,ph(l,:),dph(l,:),ct,dct)
          write(*,'(5f10.5)') (ct - cts(l,:))/ct
       end do
    end if

    call stdtra(pairs,cts,dcts,ph,dph,tr1,tr1err,verbose)
    if( verbose ) then
       write(*,*) 'Backward matrix (cts-rate to ph-rate):'
       do k = 1,size(tr1,1)
          write(*,'(5f12.4)') tr1(k,:)
       end do
    end if

    if( verbose ) then
       write(*,*) 'Relative residuals in filters:'
       write(*,'(tr5,5a10)') filters
       do l = 1,size(ph,1)
          call fotra(tr1,pairs,cts(l,:),dcts(l,:),ct,dct)
          write(*,'(5f10.5)') (ph(l,:) - ct)/ph(l,:)
       end do
    end if

    deallocate(ct,dct)

  end subroutine trafo

end program phfotran
