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
!
!  Two models (MODEL):
!   MAG  (default) colour terms of the magnitudes of every filter by
!        the colour of its neighbour, m_std - m_inst = ZP + k1 (c - cref)
!        [+ k2 (c - cref)**2] (the standard transformation; magtrafo),
!   FLUX the linear mixing of the photon rates of the standard filters,
!        counts = T photons (a nearly tridiagonal matrix; trafo).


program phfotran

  use titsio
  use mfits
  use phsysfits
  use iso_fortran_env

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
  ! the model of the transformation: MAG (colour terms of the magnitudes,
  ! of the order color_order) or FLUX (the linear mixing of photon rates)
  character(len=FLEN_VALUE) :: model = 'MAG'
  integer :: color_order = 1

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

     else if( key == 'MODEL' ) then

        read(val,*) model
        if( .not. (model == 'MAG' .or. model == 'FLUX') ) &
             stop 'Error: the model is MAG or FLUX.'

     else if( key == 'COLOR_ORDER' ) then

        read(val,*) color_order
        if( color_order < 1 .or. color_order > 2 ) &
             stop 'Error: the order of the colour terms is 1 or 2.'

     else if( key == 'TOL' ) then

        read(val,*) utol

     else if( key == 'AREA' ) then

        read(val,*) area
        init_area = .true.

     else if( key == 'CAT' ) then

        read(val,*) cat

     else if( key == 'OUTPUT' ) then

        read(val,*) output

     else if( key == 'NFILE' .or. key == 'NFILES' ) then

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

  if( list ) then
     call listphsys(phsystable)
  else if ( nfile == 1 .and. cat == '' ) then
     if( photsys_instr == '' .or. photsys_ref == '' ) &
          stop 'Reference or instrumental photometry system undefined.'
     call table
  else
     if( photsys_ref == '' ) stop 'Reference photometry system undefined.'
     if( cat == '' ) stop 'Reference catalogue undefined.'
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
    character(len=FLEN_VALUE), dimension(:), allocatable :: filter, filter_frames
    real(dbl), allocatable, dimension(:) :: ra, dec, airmass
    real(dbl), allocatable, dimension(:,:) :: tr,tr1,trerr,tr1err, ph,dph,cts,dcts
    integer, dimension(:,:), allocatable :: pairs
    integer, dimension(:), allocatable :: order, perm
    type(type_phsys) :: phsys
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
    if( size(col_mag) /= size(filters) ) &
         stop 'Magnitude columns and filters dimensions differs.'
    if( size(col_magerr) > 0 .and. size(col_magerr) /= size(filters) ) &
         stop 'Magnitude error columns and filters dimensions differs.'

    ! The transformation supposes the filters in the order of the photometry
    ! system (by wavelength): the frames are arranged so, and the filters,
    ! the catalogue columns and the -C, -E values (in the order of -f)
    ! are paired with the frames by the filters of the frames (as phcal).
    call phselect(phsystable,photsys_ref,phsys)
    call file_order(filenames,keys,phsys%filter,order,filters=filter_frames)
    allocate(perm(size(filenames)))
    call filters_arrange(filter_frames,filters,perm)
    filenames = filenames(order)
    filters = filters(perm)
    col_mag = col_mag(perm)
    if( size(col_magerr) > 0 ) col_magerr = col_magerr(perm)
    if( nctph > 0 ) ctph = ctph(perm)
    if( nextin > 0 ) extin = extin(perm)
    deallocate(order,perm,filter_frames)

    call jamcatx(cat,phsystable,keys,col_ra,col_dec,col_mag,col_magerr,&
       utol, area, init_area, filters, filenames, catid, filter, ra, dec, airmass, &
       ph,dph, cts,dcts, photsys_ref, photsys_frames, pairs)

    if( size(ra) < 2 ) &
         stop 'Determination of parameters of line requires more than two points.'

    ! the instrumental system: --photsys-instr, or PHOTSYS of the frames
    if( photsys_instr == '' ) photsys_instr = photsys_frames

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

    if( model == 'MAG' ) then
       call magtrafo(phsys,ph,dph,cts,dcts,verbose,ra,dec,filenames,airmass)
    else
       call trafo(pairs,ph,dph,cts,dcts,tr,trerr,tr1,tr1err,verbose)

       call trawrite(output,photsys_ref,photsys_instr,filters,filters, &
            ra,dec,ci0,ph,dph,cts,dcts,tr,trerr,tr1,tr1err,filenames,airmass,ctph,extin)
    end if

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
    if( size(filters) /= n ) &
         stop 'Filters (-f) do not match the data of the table.'

    if( model == 'MAG' ) then
       call magtrafo(phsys,ph,dph,cts,dcts,verbose,ra,dec)
    else
       allocate(tr(n,n),trerr(n,n),tr1(n,n),tr1err(n,n))

       call trafo(pairs,ph,dph,cts,dcts,tr,trerr,tr1,tr1err,verbose)

       call trawrite(output,photsys_ref,photsys_instr,filters,filters, &
            ra,dec,ci0,ph,dph,cts,dcts,tr,trerr,tr1,tr1err)
       deallocate(tr,trerr,tr1,tr1err)
    end if

    deallocate(ra,dec,cts,dcts,ph,dph)

  end subroutine table


  subroutine magtrafo(phsys,ph,dph,cts,dcts,verbose,ra,dec,filenames,airmass)

    ! The magnitude transformation: for every filter n (in the order of
    ! the photometry system),
    !
    !    m_std(n) - m_inst(n) = ZP(n) + k1(n) c(n) + k2(n) c(n)**2,
    !
    ! m_std the catalogue magnitudes, m_inst = -2.5 log10(count rate),
    ! c(n) the colour of the standard magnitudes of the filter and its
    ! neighbour: m(n-1) - m(n), for the first filter m(1) - m(2) (B-V for
    ! B and V, V-R for R of Johnson BVR). The fit is the one of the colour
    ! term of phcal (colterm): the mixture likelihood with the systematic
    ! error, the errors of the colours propagated, outliers without weight.

    use colterm
    use fits_fotran

    type(type_phsys), intent(in) :: phsys
    real(dbl), dimension(:,:), intent(in) :: ph,dph,cts,dcts
    logical, intent(in) :: verbose
    real(dbl), dimension(:), intent(in) :: ra,dec
    character(len=*), dimension(:), intent(in), optional :: filenames
    real(dbl), dimension(:), intent(in), optional :: airmass

    character(len=FLEN_VALUE), dimension(:), allocatable :: cola, colb
    real(dbl), dimension(:,:), allocatable :: m, dm
    real(dbl), dimension(:), allocatable :: y, e, x, ex, cxy, q, r, cref, zp, dzp, &
         k1, dk1, k2, dk2, sys, rms
    integer, dimension(:), allocatable :: nstars
    real(dbl), dimension(0:2) :: b
    real(dbl), dimension(3,3) :: cov
    real(dbl), dimension(1) :: m1, dm1, p1, dp1
    real(dbl) :: sig, fout, d, s0, s1
    integer :: nb, nd, n, ia, ib

    nd = size(ph,1)
    nb = size(ph,2)
    if( nb < 2 ) stop 'The magnitude transformation needs two filters at least.'
    if( nd < 3 ) stop 'Too few stars for the magnitude transformation.'

    allocate(m(nd,nb),dm(nd,nb),y(nd),e(nd),x(nd),ex(nd),cxy(nd),q(nd),r(nd), &
         cref(nb),zp(nb),dzp(nb),k1(nb),dk1(nb),k2(nb),dk2(nb),sys(nb),rms(nb), &
         nstars(nb),cola(nb),colb(nb))

    ! the catalogue magnitudes (the inverse of the photons of jamcatx)
    do n = 1, nb
       call phsysphmag1(phsys,filters(n),ph(:,n),dph(:,n),m(:,n),dm(:,n))
    end do

    if( verbose ) write(*,*) 'Magnitude transformation: m_std - m_inst = ', &
         'ZP + k1*(c - cref) + k2*(c - cref)**2:'

    do n = 1, nb

       if( n == 1 ) then
          ia = 1
          ib = 2
       else
          ia = n - 1
          ib = n
       end if
       cola(n) = filters(ia)
       colb(n) = filters(ib)

       where( ph(:,n) > 0 .and. cts(:,n) > 0 .and. dcts(:,n) > 0 )
          y = log(ph(:,n) / cts(:,n))
          e = sqrt(max(dph(:,n),0.0_dbl)**2/ph(:,n)**2 + (dcts(:,n)/cts(:,n))**2)
       elsewhere
          y = 0
          e = -1
       end where
       where( m(:,ia) > 98 .or. m(:,ib) > 98 ) e = -1
       if( count(e > 0) < 3 ) stop 'Too few stars for the magnitude transformation.'
       ! the reference colour: the mean colour of the stars, so the zero
       ! point and the coefficients are nearly uncorrelated (phcal adds
       ! the errors of the coefficients of the table to the errors
       ! of the stars as independent)
       cref(n) = sum(m(:,ia) - m(:,ib),mask=e > 0) / count(e > 0)
       x = m(:,ia) - m(:,ib) - cref(n)
       where( dm(:,ia) < 9 .and. dm(:,ib) < 9 )
          ex = sqrt(dm(:,ia)**2 + dm(:,ib)**2)
       elsewhere
          ex = 0
       end where
       d = 0
       if( ia == n ) d = d - 1
       if( ib == n ) d = d + 1
       where( dm(:,ia) < 9 .and. dm(:,ib) < 9 )
          cxy = d * dm(:,n)**2 / magln
       elsewhere
          cxy = 0
       end where

       call mixture_fit(y,e,x,ex,cxy,color_order,-1.0_dbl,b,cov,sig,fout,q)

       ! the zero point: m_std - m_inst = magln (ln ph0 - b0)
       m1 = 0
       dm1 = 0
       call phsysmagph1(phsys,filters(n),m1,dm1,p1,dp1)
       zp(n) = magln * (log(p1(1)) - b(0))
       dzp(n) = magln * sqrt(max(cov(1,1),0.0_dbl))
       k1(n) = -magln * b(1)
       dk1(n) = magln * sqrt(max(cov(2,2),0.0_dbl))
       k2(n) = -magln * b(2)
       dk2(n) = magln * sqrt(max(cov(3,3),0.0_dbl))
       sys(n) = magln * sig
       nstars(n) = count(e > 0)

       ! the scatter of the calibration stars [mag]: with the colour
       ! terms and with the zero point alone (the colour terms ignored)
       r = magln * (y - b(0) - b(1)*x - b(2)*x**2)
       s1 = sqrt(sum(q*r**2,mask=e > 0) / max(sum(q,mask=e > 0),1.0_dbl))
       rms(n) = s1
       if( verbose ) then
          r = magln * y
          s0 = sum(q*r,mask=e > 0) / max(sum(q,mask=e > 0),1.0_dbl)
          s0 = sqrt(sum(q*(r - s0)**2,mask=e > 0) / max(sum(q,mask=e > 0),1.0_dbl))
          write(*,'(1x,6a,f0.4,a,f0.4,a,f0.4,a,f0.5,a,f0.5)') trim(filters(n)), &
               ' (c = ',trim(cola(n)),'-',trim(colb(n)),'): cref = ',cref(n), &
               ', ZP = ',zp(n),' +- ',dzp(n),', k1 = ',k1(n),' +- ',dk1(n)
          if( color_order >= 2 ) write(*,'(3x,a,f0.5,a,f0.5)') 'k2 = ',k2(n), &
               ' +- ',dk2(n)
          write(*,'(3x,a,f0.4,a,f0.3,a,i0,a,f0.4,a,f0.4,a)') 'systematic error: ', &
               sys(n),' mag, outliers: ',fout,', stars: ',nstars(n), &
               ', rms ',s0,' -> ',s1,' mag'
       end if
    end do

    if( present(filenames) ) then
       call trawrite_mag(output,photsys_ref,photsys_instr,filters,filters, &
            cola,colb,cref,zp,dzp,k1,dk1,k2,dk2,sys,nstars,rms, &
            ra,dec,ph,dph,cts,dcts,filenames,airmass,ctph,extin)
    else
       call trawrite_mag(output,photsys_ref,photsys_instr,filters,filters, &
            cola,colb,cref,zp,dzp,k1,dk1,k2,dk2,sys,nstars,rms, &
            ra,dec,ph,dph,cts,dcts)
    end if

  end subroutine magtrafo


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
