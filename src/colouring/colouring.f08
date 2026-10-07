!
! colouring
!
! Copyright © 2010-2, 2018-9 F.Hroch (hroch@physics.muni.cz)
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

program colouring

  use colouring_book
  use fitscolour
  use ctables
  use titsio
  use iso_fortran_env

  implicit none

  character(len=4*FLEN_FILENAME) :: record, key, val
  character(len=FLEN_FILENAME) :: filename
  type(ColourFits), dimension(:), allocatable :: bands
  character(len=FLEN_KEYWORD), dimension(4) :: fitskeys
  real, dimension(:,:), allocatable :: cmatrix
  real, dimension(:), allocatable :: weights, ctphs, backs
  character(len=80) :: msg
  integer :: eq,stat,status
  real :: x, y, r = 7  ! for white spot, star
  integer :: bitpix = -32
  integer :: nbands = 0
  integer :: nweights = 0
  integer :: nctphs = 0
  integer :: nbacks = 0
  logical :: verbose = .false.   ! verbose output
  logical :: estimbacks = .true. ! estimate backgrounds
  integer :: white = 0           ! zero means undefined
  logical :: list
  character(len=FLEN_FILENAME) :: output = 'colour.fits'
  character(len=FLEN_FILENAME) :: ctable = 'ctable.dat'
  character(len=FLEN_FILENAME) :: phsystable = 'photosystems.fits'
  character(len=FLEN_VALUE) :: cspace = ''
  character(len=FLEN_VALUE) :: outspace = 'CIE 1931 XYZ'

  ! mandatory keys
  fitskeys(1) = FITS_KEY_FILTER
  fitskeys(2) = FITS_KEY_EXPTIME
  fitskeys(3) = FITS_KEY_AREA
  fitskeys(4) = FITS_KEY_OBJECT

  do
     read(*,'(a)',iostat=stat,iomsg=msg) record
     if( stat == IOSTAT_END ) exit
     if( stat > 0 ) then
        write(error_unit,*) trim(msg)
        error stop 'An input error.'
     end if

     eq = index(record,'=')
     if( eq == 0 ) stop 'Malformed input record.'

     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'OUTPUT' ) then

        read(val,*) output

     else if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'FITS_KEY_FILTER' ) then

        read(val,*) fitskeys(1)

     else if( key == 'FITS_KEY_EXPTIME' ) then

        read(val,*) fitskeys(2)

     else if( key == 'FITS_KEY_AREA' ) then

        read(val,*) fitskeys(3)

     else if( key == 'FITS_KEY_OBJECT' ) then

        read(val,*) fitskeys(4)

     else if( key == 'BITPIX' ) then

        read(val,*) bitpix

    else if( key == 'PHSYSTABLE' ) then

        read(val,*) phsystable

     else if( key == 'CTABLE' ) then

        read(val,*) ctable

     else if( key == 'COLOURSPACE' ) then

        read(val,*) cspace

     else if( key == 'OUTERSPACE' ) then

        read(val,*) outspace

     else if( key == 'ESTIMBACKS' ) then

        read(val,*) estimbacks

     else if( key == 'WRADIUS' ) then

        read(val,*) r

     else if( key == 'WSPOT' ) then

        read(val,*) x,y
        white = WHITE_SPOT

     else if( key == 'WSTAR' ) then

        read(val,*) x,y
        white = WHITE_STAR

     else if( key == 'NWEIGHTS' ) then

        read(val,*) nweights
        allocate(weights(nweights))
        white = WHITE_WEIGHTS

     else if( key == 'WEIGHTS' ) then

        read(val,*) weights

     else if( key == 'NCTPHS' ) then

        read(val,*) nctphs
        allocate(ctphs(nctphs))
        white = WHITE_CTPH

     else if( key == 'CTPHS' ) then

        read(val,*) ctphs

     else if( key == 'NBACKS' ) then

        read(val,*) nbacks
        allocate(backs(nbacks))

     else if( key == 'BACKS' ) then

        read(val,*) backs

     else if( key == 'NBANDS' ) then

        read(val,*) nbands
        allocate(bands(nbands))
        nbands = 0

     else if( key == 'BAND' ) then

        read(val,*) filename

        nbands = nbands + 1
        call bands(nbands)%Load(filename,fitskeys)
        if( .not. bands(nbands)%status ) then
           write(error_unit,*) "Error: file `",trim(filename),"' read failed."
           stop 'FITS read fail.'
        end if

        if( nbands > 1 ) then
           if( .not. all(bands(nbands-1)%naxes == bands(nbands)%naxes) ) then
              if( verbose ) write(error_unit,*) "Current file: `",filename,"'."
              stop "Dimensions of images does not corresponds each other."
           end if
        endif

     end if

     if( key == 'LIST' ) then
        read(val,*) list
        if( list ) then
           call ctable_list(ctable)
           stop 0
        end if
     end if

  end do

  if( size(bands) == 0 ) stop 'No bands, no love.'
  if( .not. (size(bands) == nbands) ) stop 'A full set of FITS files required.'

  ! input colour space guess
  if( cspace == '' ) then
     call colour_oracle(bands,cspace)
     if( cspace == '' ) &
          stop 'Error: the built in oracle is unable to determine a colour space.'
  end if

  ! check calibration
  if( white == 0 .and. all(bands(:)%calibrated) ) white = WHITE_CAL

  ! load transformation table
  call ctable_load(ctable,cspace,outspace,cmatrix,verbose)

  ! load wavelengths
  call colour_waves(bands,cspace,phsystable)

  ! backs
  if( nbacks > 0 ) estimbacks = .false.

  call colour(cspace,outspace,bands,bitpix,cmatrix,x,y,r,white,estimbacks,&
       weights,ctphs,backs,output,verbose,status)

  deallocate(bands,cmatrix)

  if( status == 0 )then
     stop 0
  else
     stop 'An error occurred.'
  end if

contains

  subroutine colour_oracle(bands,cspace)

    type(ColourFits), dimension(:), intent(in) :: bands
    character(len=*), intent(in out) :: cspace

    integer :: k

    if( size(bands) == 3 ) then

       if( all(bands(:)%filter == ['B', 'V', 'R']) ) then
          cspace = 'Johnson BVR'
       else if ( all(bands(:)%filter == ['Z', 'Y', 'X']) ) then
          cspace = 'CIE 1931 XYZ'
       end if

    end if

    if( verbose ) then
       if( cspace /= '' ) then
          write(error_unit,*) "Info: Colour space of inputs: '", &
               trim(cspace),"'"
       else
          write(error_unit,*) 'Warning: Failed to guess a colour space',&
               ' for filters: ',("'", &
               trim(bands(k)%filter),"' ",k=1,size(bands)),&
               '(in this order).'
       end if
    end if

  end subroutine colour_oracle

  subroutine colour_waves(bands,cspace,phsystable)

    use phsysfits

    type(ColourFits), dimension(:), intent(in out) :: bands
    character(len=*), intent(in) :: cspace, phsystable

    type(type_phsys) :: phsys
    character(len=*), parameter :: Johnson = 'Johnson'
    integer :: l,k
    logical :: found

    if( index(cspace,Johnson) > 0 ) then

       call phselect(phsystable,Johnson,phsys)
       do l = 1, size(bands)
          found = .false.
          do k = 1, size(phsys%filter)
             if( bands(l)%filter == phsys%filter(k) ) then
                found = .true.
                bands(l)%leff = real(phsys%lam_eff(k))
                bands(l)%fref = &
                 real(2.50663 * phsys%flam_ref(k) * 1e9*phsys%lam_fwhm(k) / 2.0)
             end if
          end do
          if( .not. found ) then
             write(error_unit,*) 'Error: For filter `',trim(bands(l)%filter),"'"
             stop 'Error: filter not found.'
          end if
       end do

!!$       if( verbose ) then
!!$          write(error_unit,*) "Info: Effective wavelengts [nm] of filters: '"
!!$          do k = 1, size(bands)
!!$             write(error_unit,'(3a,f7.1)') &
!!$                  " Info: ",trim(bands(k)%filter),'   ',1e9*bands(k)%leff
!!$          end do
!!$       end if

    end if

  end subroutine colour_waves

end program colouring
