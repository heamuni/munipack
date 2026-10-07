!
!  Photometric systems table
!
!  Copyright © 2013 - 2015, 2018-26 F.Hroch (hroch@physics.muni.cz)
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


module phsysfits

  use titsio
  use iso_fortran_env

  implicit none
  private

  integer, parameter :: maxbands = 13
  ! The photometry systems are implemented as static tables.
  ! They offer much more simpler handling, particulary for coarray
  ! programs. It is not significant limitation as common
  ! photometric systems has few filters only.

  type, public :: type_phsys
     character(len=FLEN_VALUE) :: phsystem = '', filter_ref = ''
     character(len=FLEN_VALUE), dimension(maxbands) :: filter
     real(REAL64), dimension(maxbands) :: lam_eff,lam_fwhm,nu_eff,nu_fwhm, &
          fnu_ref,flam_ref
     integer :: nbands = 0
  end type type_phsys

  public :: listphsys, phselect, selphsystem, phsysmagph, &
       phsysphmag1, phsysmagph1, phsysphmag, phsysmagflux, phsysconv, &
       phstkeys, phsyspairs

contains


  subroutine phsysread(filename,phsys,status)

    character(len=*), intent(in) :: filename
    type(type_phsys), dimension(:), allocatable, intent(in out) :: phsys
    integer, intent(out) :: status

    real(REAL64), parameter :: nullval = 0.0_REAL64

    integer :: nrows, ncols, hdutype, nhdu, i, n, m
    integer, parameter :: maxcols = 7
    integer, dimension(maxcols) :: col
    character(len=FLEN_VALUE), dimension(maxcols) :: cols
    character(len=FLEN_VALUE) :: photosys
    logical :: anyf
    type(fitsfiles) :: fits

    cols(1) = FITS_COL_FILTER
    cols(2) = FITS_COL_LAMEFF
    cols(3) = FITS_COL_LAMFWHM
    cols(4) = FITS_COL_NUEFF
    cols(5) = FITS_COL_NUFWHM
    cols(6) = FITS_COL_FNUREF
    cols(7) = FITS_COL_FLAMREF

    status = 0

    ! open, and get number of HDU
    call fits_open_file(fits,filename,FITS_READONLY,status)
    if( status /= 0 ) then
       write(error_unit,*) 'Error: failed to read the file `',trim(filename),"'."
       return
    end if

    call fits_get_num_hdus(fits,nhdu,status)
    if( status /= 0 ) goto 666

    ! checkpoint
    call fits_read_key(fits,'HDUNAME',photosys,status)
    if( status == FITS_KEYWORD_NOT_FOUND .or. photosys /= FHDUNAME ) then
       status = 0
       write(error_unit,*) "Warning: Photometry system file `",trim(filename), &
            "' has wrong identificator."
    end if
    if( status /= 0 ) goto 666

    allocate(phsys(nhdu-1))

    do n = 1, nhdu-1

       call fits_movabs_hdu(fits,n+1,hdutype,status)
       if( hdutype /= FITS_BINARY_TBL .or. status /= 0 ) goto 666

       call fits_read_key(fits,'EXTNAME',phsys(n)%phsystem,status)
       call fits_read_key(fits,FITS_KEY_FILTREF,phsys(n)%filter_ref,status)
       if( status == FITS_KEYWORD_NOT_FOUND .or. status /= 0 ) goto 666

       call fits_get_num_cols(fits,ncols,status)
       call fits_get_num_rows(fits,nrows,status)

       ! find columns by cols
       do i = 1, size(cols)
          call fits_get_colnum(fits,.true.,cols(i),col(i),status)
       end do
       if( status /= 0 ) goto 666
       if( ncols /= size(cols) ) stop 'ncols /= size(cols)'

       if( nrows > maxbands ) &
            write(error_unit,*) 'Warning in phsys: Maximum of filters reached.'
       m = min(nrows,maxbands)

       phsys(n)%filter = ''
       call fits_read_col(fits,col(1),1,'',phsys(n)%filter(1:m),anyf,status)
       call fits_read_col(fits,col(2),1,nullval,phsys(n)%lam_eff(1:m),anyf,status)
       call fits_read_col(fits,col(3),1,nullval,phsys(n)%lam_fwhm(1:m),anyf,status)
       call fits_read_col(fits,col(4),1,nullval,phsys(n)%nu_eff(1:m),anyf,status)
       call fits_read_col(fits,col(5),1,nullval,phsys(n)%nu_fwhm(1:m),anyf,status)
       call fits_read_col(fits,col(6),1,nullval,phsys(n)%fnu_ref(1:m),anyf,status)
       call fits_read_col(fits,col(7),1,nullval,phsys(n)%flam_ref(1:m),anyf,status)
       if( status /= 0 ) goto 666

       phsys(n)%nbands = m
    end do

    call fits_close_file(fits,status)
    return

666 continue

    call fits_close_file(fits,status)
    call fits_report_error(error_unit,status)

  end subroutine phsysread

  subroutine listphsys(phsystable)

    character(len=*), intent(in) :: phsystable
    type(type_phsys), dimension(:), allocatable :: phsyscal
    integer :: status,i,j

    call phsysread(phsystable,phsyscal,status)
    if( status /= 0 ) stop 'Failed to read photometry system table.'

    write(*,*)
    do i = 1,size(phsyscal)
       write(*,'(2a)',advance="no") trim(phsyscal(i)%phsystem),": "
       do j = 1, phsyscal(i)%nbands
          write(*,'(a)',advance="no") trim(phsyscal(i)%filter(j))
          if( j == phsyscal(i)%nbands ) then
             write(*,*)
          else
             write(*,'(a)',advance="no") ","
          end if
       end do
    end do

    deallocate(phsyscal)

  end subroutine listphsys

  subroutine phselect(phsystable,phsystem,phsys)

    character(len=*), intent(in) :: phsystable,phsystem
    type(type_phsys), intent(out) :: phsys

    type(type_phsys), dimension(:), allocatable :: phsyscal
    integer :: i,n,status

    if( phsystem == '' ) stop 'Photometry system undefined.'

    call phsysread(phsystable,phsyscal,status)
    if( status /= 0 ) stop 'Failed to read the photometry system table.'

    n = -1
    do i = 1, size(phsyscal)
       if( phsyscal(i)%phsystem == phsystem ) then
          n = i
          exit
       end if
    end do

    if( n > 0 ) then
       phsys = phsyscal(n)
    else
       stop 'Specified photometry system not found.'
    end if

  end subroutine phselect


  subroutine selphsystem(phsyscal,phsystem,filter,leff,lfwhm,feff,ffwhm,flam_ref, &
       fnu_ref,status)

    type(type_phsys), dimension(:), intent(in) :: phsyscal
    character(len=*), intent(in) :: phsystem,filter
    real(REAL64), intent(out) :: leff,lfwhm,feff,ffwhm,flam_ref,fnu_ref
    integer, intent(out) :: status

    integer :: i,j

    status = -1

    do i = 1,size(phsyscal)
       if( phsyscal(i)%phsystem == phsystem ) then
          do j = 1, phsyscal(i)%nbands
             if( phsyscal(i)%filter(j) == filter ) then
                status = 0
                feff = phsyscal(i)%nu_eff(j)
                ffwhm = phsyscal(i)%nu_fwhm(j)
                leff = phsyscal(i)%lam_eff(j)
                lfwhm = phsyscal(i)%lam_fwhm(j)
                flam_ref = phsyscal(i)%flam_ref(j)
                fnu_ref = phsyscal(i)%fnu_ref(j)
                return
             end if
          end do
       end if
    end do

  end subroutine selphsystem


  subroutine phsysmagph(phsys,filter,pairs,mag,dmag,ph,dph)

    use photoconv

    type(type_phsys), intent(in) :: phsys
    character(len=*), dimension(:), intent(in) :: filter
    integer, dimension(:,:), intent(in) :: pairs
    real(REAL64), dimension(:,:), intent(in) :: mag,dmag
    real(REAL64), dimension(:,:), intent(out) :: ph,dph

    real(REAL64), parameter :: eps = epsilon(ph)
    real(REAL64), dimension(size(filter)) :: feff,ffwhm,leff,lfwhm,flamref,fnuref
    integer :: n,i,j,k,status
!    real(REAL64), dimension(:,:), allocatable :: ph0,dph0,ph3,dph3,mag0,dmag0, &
!         magx,dmagx,mag3,dmag3,flux,dflux
    real(REAL64), dimension(:,:), allocatable :: flux,dflux
    real(REAL64) :: w

    if( size(mag,2) /= size(filter) ) stop 'Counts of filters and magnitudes different.'

    j = 1
    do n = 1, size(mag,2)
       call selphsystem([phsys],phsys%phsystem,filter(n), &
            leff(n),lfwhm(n),feff(n),ffwhm(n),flamref(n),fnuref(n),status)

       if( status /= 0 ) then
          write(error_unit,*) 'Filter `',trim(filter(n)),"' not found in system `",&
               trim(phsys%phsystem),"'."
          stop 'Filter not found.'
       end if
    end do

!    allocate(ph0(size(ph,1),size(ph,2)),dph0(size(ph,1),size(ph,2)))

!    write(*,*) leff(2),lfwhm(2),flamref(2)
!    do j = 1,size(mag,2)
!       call mag2ph0(leff(j),lfwhm(j),flamref(j),mag(:,j),dmag(:,j),ph(:,j),dph(:,j))
!    end do
!    write(*,'(5e10.3)') ph(1,:)
!    call mag2ph(leff(j),lfwhm(j),flamref(j),mag(:,j),dmag(:,j),ph(:,j),dph(:,j))


!    do n = 1,size(mag,1)
!       call mag2ph3(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph(n,:),dph(n,:))
!    end do

!    do j = 1,size(mag,1)
!       write(*,*) j,mag(j,3),dmag(j,3),ph(j,3),dph(j,3)
!    end do

!    allocate(mag0(size(ph,1),size(ph,2)),dmag0(size(ph,1),size(ph,2)))

!    call ph2mag0(leff(1),lfwhm(1),flamref(1),ph(:,1),dph(:,1),mag0(:,1),dmag0(:,1))
!    do j = 1,size(mag,1)
!       write(*,'(5g15.5)') real(mag(j,3)),real(ph(j,3)),real(ph0(j,3)),real(ph0(j,2)/ph0(j,3)),&
!            real(ph0(j,3)/ph(j,3))
!    end do

    allocate(flux(size(ph,1),size(ph,2)),dflux(size(ph,1),size(ph,2)))
    call phsysmagflux(phsys,filter,pairs,mag,dmag,flux,dflux)

    k = pairs(1,1)
    call mag2ph0(leff(k),lfwhm(k),flamref(k),mag(:,k),dmag(:,k),ph(:,k),dph(:,k))

    do n = 2, size(pairs,1)
       i = pairs(n,1)
       j = pairs(n,2)
       if( pairs(n,1) < k ) then
          i = pairs(n,1)
          j = pairs(n,2)
          w = leff(i) / leff(j)
          where( flux(:,i) > 0 .and. flux(:,j) > 0 .and. ph(:,j) > 0 )
             ph(:,i) = ph(:,j)*(flux(:,i)/flux(:,j))*w
          elsewhere
             ph(:,i) = -1
          end where
          where( dflux(:,i) > 0 .and. ph(:,i) > 0 .and. flux(:,i) > 0)
             dph(:,i) = dflux(:,i) * ph(:,i) / flux(:,i)
          elsewhere
             dph(:,i) = 0
          end where
       else if( pairs(n,2) > k ) then
          i = pairs(n,1)
          j = pairs(n,2)
          w = leff(j) / leff(i)
          where( flux(:,i) > 0 .and. flux(:,j) > 0 .and. ph(:,i) > 0 )
             ph(:,j) = ph(:,i)*(flux(:,j)/flux(:,i))*w
          elsewhere
             ph(:,j) = -1
          end where
          where( dflux(:,j) > 0 .and. ph(:,j) > 0 .and. flux(:,j) > 0)
             dph(:,j) = dflux(:,j) * ph(:,j) / flux(:,j)
          elsewhere
             dph(:,j) = 0
          end where
       end if
    end do
    deallocate(flux,dflux)



!    allocate(ph0(size(ph,1),size(ph,2)),dph0(size(ph,1),size(ph,2)), &
!         ph3(size(ph,1),size(ph,2)),dph3(size(ph,1),size(ph,2)), &
!         mag3(size(ph,1),size(ph,2)),dmag3(size(ph,1),size(ph,2)),&
!         mag0(size(ph,1),size(ph,2)),dmag0(size(ph,1),size(ph,2)),&
!         magx(size(ph,1),size(ph,2)),dmagx(size(ph,1),size(ph,2)),&
!         )



!!$    goto 666
!!$
!!$
!!$!    do n = 1,size(mag,1)
!!$!       call mag2ph(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph(n,:),dph(n,:))
!!$!       call mag2ph0(leff,flamref,mag(n,:),dmag(n,:),ph0(n,:),dph0(n,:))
!!$!       write(*,*) real(mag(n,:)),real(ph(n,:))
!!$!    end do
!!$
!!$    do n = 1,size(mag,1)
!!$!       call mag2ph(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph(n,:),dph(n,:))
!!$!       call mag2ph3(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph3(n,:),dph3(n,:))
!!$!       write(*,*) real(mag(n,:)),real(ph(n,:))
!!$    end do
!!$
!!$!    do n = 1,size(mag,1)
!!$!       call mag2ph(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph(n,:),dph(n,:))
!!$!       call mag2ph(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph(n,:),dph(n,:))
!!$!       write(*,*) real(mag(n,:)),real(ph(n,:))
!!$!    end do
!!$
!!$!    do n = 1,size(mag,1)
!!$!       write(*,'(2f10.3,3g15.5)') mag(n,3),mag(n,2)-mag(n,3),ph0(n,3),ph(n,3),ph3(n,3)
!!$!    end do
!!$
!!$    do j = 1,size(mag,2)
!!$!       call mag2ph0(leff(j),flamref(j),mag(:,j),dmag(:,j),ph0(:,j),dph0(:,j))
!!$!       call ph2mag0(leff(j),flamref(j),ph0(:,j),dph0(:,j),mag0(:,j),dmag0(:,j))
!!$    end do
!!$    do n = 1,size(mag,1)
!!$       call mag2ph3(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph3(n,:),dph3(n,:))
!!$!       call ph2mag3(leff,lfwhm,flamref,ph3(n,:),dph3(n,:),mag3(n,:),dmag3(n,:))
!!$    end do
!!$
!!$    do j = 1,size(mag,2)-2
!!$!       write(*,*) trim(filter(j))
!!$       do n = 1,size(mag,1)
!!$          write(*,'(2f10.3,3g15.5)') mag(n,j),mag(n,j)-mag(n,j+1),ph0(n,j),ph3(n,j)
!!$       end do
!!$       do n = 1,size(mag,1)
!!$          write(*,'(5f10.3)') mag(n,j),mag(n,j)-mag(n,j+1),mag0(n,j),mag3(n,j)
!!$       end do
!!$    end do
!!$    do j = size(mag,2)-1,size(mag,2)
!!$       write(*,*) trim(filter(j))
!!$!       do n = 1,size(mag,1)
!!$!          call mag2ph3(leff,lfwhm,flamref,mag(j,:),dmag(j,:),ph3(j,:),dph3(j,:))
!!$!          call ph2mag3(leff,lfwhm,flamref,ph3(j,:),dph3(j,:),mag3(j,:),dmag3(j,:))
!!$!       end do
!!$       do n = 1,size(mag,1)
!!$          write(*,'(2f10.3,3g15.5)') mag(n,j),mag(n,j-1)-mag(n,j),ph0(n,j),ph3(n,j)
!!$       end do
!!$       do n = 1,size(mag,1)
!!$          write(*,'(5f10.3)') mag(n,j),mag(n,j-1)-mag(n,j),mag0(n,j),mag3(n,j)
!!$       end do
!!$    end do
!!$
!!$    do n = 1,size(mag,1)
!!$!       call ph2mag3(leff,lfwhm,flamref,ph(n,:),dph(n,:),magx(n,:),dmagx(n,:))
!!$!       call ph2mag0(lfwhm(3),flamref(3),ph0(:,3),dph0(:,3),mag0(:,3),dmag0(:,3))
!!$!       call ph2mag0(leff(3),lfwhm(3),flamref(3),ph(:,3),dph(:,3),mag0(:,3),dmag0(:,3))
!!$!       call ph2mag3(leff,lfwhm,flamref,ph3(n,:),dph3(n,:),mag3(n,:),dmag3(n,:))
!!$    end do
!!$
!!$
!!$    do j = 1,size(mag,2)-2
!!$       write(*,*) trim(filter(j))
!!$       do n = 1,size(mag,1)
!!$          write(*,'(5f10.3)') mag(n,j),mag(n,j)-mag(n,j+1),mag0(n,j),mag3(n,j)
!!$       end do
!!$    end do
!!$    do j = size(mag,2)-1,size(mag,2)
!!$       write(*,*) trim(filter(j))
!!$       do n = 1,size(mag,1)
!!$          write(*,'(5f10.3)') mag(n,j),mag(n,j-1)-mag(n,j),mag0(n,j),mag3(n,j)
!!$       end do
!!$    end do
!!$
!!$!    stop
!!$
!!$!    do n = 1,size(mag,1)
!!$!       call mag2ph(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph(n,:),dph(n,:))
!!$!       call mag2ph3(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph(n,:),dph(n,:))
!!$!       call mag2ph3(leff,lfwhm,flamref,mag(n,:),dmag(n,:),ph(n,:),dph(n,:))
!!$!       write(*,*) real(mag(n,:)),real(ph(n,:))
!!$!    end do
!!$
!!$    ph = ph0
!!$    dph = dph0
!!$
!!$
!!$    666 continue

    ! when magnitude error is unavailable, we are use the most optimistics
    ! variant with non-photon noise neglecting
!    where( dmag > 9 .or. dmag < epsilon(dmag) )
!       dph = sqrt(ph)
!    end where

!    deallocate(ph0,dph0,ph3,dph3,mag0,dmag0,magx,dmagx,mag3,dmag3)

  end subroutine phsysmagph


  subroutine phsysphmag1(phsys,filter,ph,dph,mag,dmag)

    use photoconv

    type(type_phsys), intent(in) :: phsys
    character(len=*), intent(in) :: filter
    real(REAL64), dimension(:), intent(in) :: ph,dph
    real(REAL64), dimension(:), intent(out) :: mag,dmag

    real(REAL64) :: feff,ffwhm,leff,lfwhm,flamref,fnuref
    integer :: status

    call selphsystem([phsys],phsys%phsystem,filter, &
         leff,lfwhm,feff,ffwhm,flamref,fnuref,status)

    call ph2mag0(leff,lfwhm,flamref,ph,dph,mag,dmag)

  end subroutine phsysphmag1

  subroutine phsysmagph1(phsys,filter,mag,dmag,ph,dph)

    use photoconv

    type(type_phsys), intent(in) :: phsys
    character(len=*), intent(in) :: filter
    real(REAL64), dimension(:), intent(in) :: mag,dmag
    real(REAL64), dimension(:), intent(out) :: ph,dph

    real(REAL64) :: feff,ffwhm,leff,lfwhm,flamref,fnuref
    integer :: status

    call selphsystem([phsys],phsys%phsystem,filter, &
         leff,lfwhm,feff,ffwhm,flamref,fnuref,status)

    call mag2ph0(leff,lfwhm,flamref,mag,dmag,ph,dph)

  end subroutine phsysmagph1

  subroutine phsysphmag(phsys,filters,pairs,ph,dph,mag,dmag)

    use photoconv

    type(type_phsys), intent(in) :: phsys
    character(len=*), dimension(:), intent(in) :: filters
    integer, dimension(:,:), intent(in) :: pairs
    real(REAL64), dimension(:,:), intent(in) :: ph,dph
    real(REAL64), dimension(:,:), intent(out) :: mag,dmag

    real(REAL64), dimension(:,:), allocatable :: flux,dflux
    integer, dimension(size(filters)) :: idx
    integer :: n,i,j,k
    real(REAL64) :: w,f

    idx = 0
    do i = 1, size(filters)
       do j = 1, phsys%nbands
          if( phsys%filter(j) == filters(i) ) then
             n = n + 1
             idx(i) = j
          end if
       end do
    end do
    if( any(idx == 0 ) ) stop 'Any filter not found (in phsysphmag).'

    allocate(flux(size(ph,1),size(ph,2)),dflux(size(ph,1),size(ph,2)))

    k = pairs(1,1)
    n = idx(k)
    call ph2flux0(phsys%lam_eff(n),ph(:,k),dph(:,k),flux(:,k),dflux(:,k))
    call flux2mag0(phsys%lam_fwhm(n),phsys%flam_ref(n),flux(:,k),dflux(:,k), &
         mag(:,k),dmag(:,k))

    do n = 2, size(pairs,1)

       if( pairs(n,1) < k ) then
          i = pairs(n,1)
          j = pairs(n,2)
          w = phsys%lam_eff(idx(j)) / phsys%lam_eff(idx(i))
          f = (phsys%lam_fwhm(idx(i))*phsys%flam_ref(idx(i)))/ &
               (phsys%lam_fwhm(idx(j))*phsys%flam_ref(idx(j)))
          where( ph(:,i) > 0 .and. ph(:,j) > 0 .and. flux(:,j) > 0 )
             flux(:,i) = flux(:,j)*(ph(:,i)/ph(:,j))*w/f
             mag(:,i) = mag(:,j) - 2.5*log10(flux(:,i) / flux(:,j))
          elsewhere
             flux(:,i) = -1
             mag(:,i) = 99.99999
          end where
          where( dph(:,i) > 0 .and. ph(:,i) > 0 .and. flux(:,i) > 0)
             dflux(:,i) = flux(:,i) * dph(:,i) / ph(:,i)
             dmag(:,i) = 1.086 * dflux(:,i) / flux(:,i)
          elsewhere
             dflux(:,i) = 0
             dmag(:,i) = 9.99999
          end where
       else if( pairs(n,2) > k ) then
          i = pairs(n,1)
          j = pairs(n,2)
          w = phsys%lam_eff(idx(i)) / phsys%lam_eff(idx(j))
          f = (phsys%lam_fwhm(idx(j))*phsys%flam_ref(idx(j)))/ &
               (phsys%lam_fwhm(idx(i))*phsys%flam_ref(idx(i)))
          where( ph(:,i) > 0 .and. ph(:,j) > 0 .and. flux(:,i) > 0 )
             flux(:,j) = flux(:,i)*(ph(:,j)/ph(:,i))*w / f
             mag(:,j) = mag(:,i) - 2.5*log10(flux(:,j) / flux(:,i))
          elsewhere
             flux(:,j) = -1
             mag(:,j) = 99.99999
          end where
          where( dph(:,j) > 0 .and. flux(:,j) > 0 .and. ph(:,j) > 0)
             dflux(:,j) = dph(:,j) * flux(:,j) / ph(:,j)
             dmag(:,j) = 1.086 * dflux(:,j) / flux(:,j)
          elsewhere
             dflux(:,j) = 0
             dmag(:,j) = 9.99999
          end where
       end if
    end do

    deallocate(flux,dflux)

  end subroutine phsysphmag

  subroutine phsysmagflux(phsys,filter,pairs,mag,dmag,flux,dflux)

    use photoconv

    type(type_phsys), intent(in) :: phsys
    character(len=*), dimension(:), intent(in) :: filter
    integer, dimension(:,:), intent(in) :: pairs
    real(REAL64), dimension(:,:), intent(in) :: mag,dmag
    real(REAL64), dimension(:,:), intent(out) :: flux,dflux

    real(REAL64), dimension(size(filter)) :: feff,ffwhm,leff,lfwhm,flamref,fnuref
    integer :: n,i,j,k,status
    real(REAL64) :: w

    if( size(mag,2) /= size(filter) ) stop 'Counts of filters and magnitudes different.'

    j = 1
    do n = 1, size(mag,2)
       call selphsystem([phsys],phsys%phsystem,filter(n), &
            leff(n),lfwhm(n),feff(n),ffwhm(n),flamref(n),fnuref(n),status)

       if( status /= 0 ) then
          write(error_unit,*) 'Filter `',trim(filter(n)),"' not found in system `",&
               trim(phsys%phsystem),"'."
          stop 'Filter not found.'
       end if
    end do

    do n = 1,size(mag,1)
       call mag2flux0(lfwhm,flamref,mag(n,:),dmag(n,:),flux(n,:),dflux(n,:))
    end do

    k = pairs(1,1)
    do n = 2,size(pairs,1)
       i = pairs(n,1)
       j = pairs(n,2)
       if( pairs(n,1) < k ) then
          w = (lfwhm(i)*flamref(i)) / (lfwhm(j)*flamref(j))
          where( mag(:,i) < 99 .and. mag(:,j) < 99 )
             flux(:,i) = w*flux(:,j)*10.0_REAL64**(-0.4_REAL64*(mag(:,i) - mag(:,j)))
          elsewhere
             flux(:,i) = -1
          end where
          where( dmag(:,i) < 9 )
             dflux(:,i) = flux(:,i) * dmag(:,i) / 1.0857
          elsewhere
             dflux(:,i) = -1
          end where
       else if( pairs(n,2) > k ) then
          w = (lfwhm(j)*flamref(j)) / (lfwhm(i)*flamref(i))
          where( mag(:,i) < 99 .and. mag(:,j) < 99 )
             flux(:,j) = w*flux(:,i)*10.0_REAL64**( 0.4_REAL64*(mag(:,i) - mag(:,j)))
          elsewhere
             flux(:,j) = -1
          end where
          where( dmag(:,j) < 9 )
             dflux(:,j) = flux(:,j) * dmag(:,i) / 1.0857
          elsewhere
             dflux(:,j) = -1
          end where
       end if
    end do

  end subroutine phsysmagflux


  subroutine phsysconv(key,filter,phsyscal,area,exptime,arcscale,ph,dph,q,dq)

    use photoconv

    character(len=*), intent(in) :: key, filter
    type(type_phsys), intent(in) :: phsyscal
    real(REAL64), intent(in) :: exptime,area,arcscale
    real(REAL64), dimension(:), intent(in) :: ph,dph
    real(REAL64), dimension(:), intent(out) :: q,dq
    real(REAL64), dimension(:), allocatable :: r,dr
    real(REAL64) :: feff,ffwhm,leff,lfwhm,fnuref,flamref
    integer :: i
    logical :: defined

    defined = .false.
    do i = 1, phsyscal%nbands
       if( filter == phsyscal%filter(i) ) then
          feff = phsyscal%nu_eff(i)
          ffwhm = phsyscal%nu_fwhm(i)
          leff = phsyscal%lam_eff(i)
          lfwhm = phsyscal%lam_fwhm(i)
          fnuref = phsyscal%fnu_ref(i)
          flamref = phsyscal%flam_ref(i)
          defined = .true.
       end if
    end do

    if( .not. defined ) then
       write(error_unit,*) "Filter `",trim(filter),"' not defined."
       stop 'Undefined filter of photometry system.'
    end if

    allocate(r(size(q)),dr(size(dq)))

    call phrate(area,exptime,arcscale,ph,dph,r,dr)

    if( key == 'FNU' ) then
       call ph2fnu(feff,ffwhm,r,dr,q,dq)
    else if( key == 'FLAM' ) then
       call ph2flam0(leff,lfwhm,r,dr,q,dq)
    else if( key == 'PHOTNU' ) then
       call ph2photnu(ffwhm,r,dr,q,dq)
    else if( key == 'PHOTLAM' ) then
       call ph2photlam(lfwhm,r,dr,q,dq)
    else if( key == 'ABMAG' ) then
       call ph2abmag(feff,ffwhm,r,dr,q,dq)
    else if( key == 'STMAG' ) then
       call ph2stmag(leff,lfwhm,r,dr,q,dq)
    else if( key == 'MAG' ) then
       call ph2mag0(leff,lfwhm,flamref,r,dr,q,dq)
    else if( key == 'FLUX' ) then
       call ph2flux(feff,r,dr,q,dq)
    else if( key == 'PHRATE' ) then
       q = r
       dq = dr
    else if( key(1:6) == 'PHOTON' ) then
       q = ph
       dq = dph
    else
       stop 'An unknown quantity for calibrated data requested.'
    end if

    deallocate(r,dr)

  end subroutine phsysconv

  subroutine phstkeys(phsyscal,filter,photflam,photzpt,photplam,photbw)

    use photoconv

    character(len=*), intent(in) :: filter
    type(type_phsys), intent(in) :: phsyscal
    real(REAL64), intent(out) :: photflam,photzpt,photplam,photbw
    integer :: i

    do i = 1, phsyscal%nbands
       if( filter == phsyscal%filter(i) ) then
          call hstparam(phsyscal%lam_eff(i),phsyscal%lam_fwhm(i), &
               photflam,photzpt,photplam,photbw)
       end if
    end do

  end subroutine phstkeys

  subroutine phsyspairs(phsys,filters,pairs)

    ! arrange filters for photometry transformation

    type(type_phsys), intent(in) :: phsys
    character(len=*), dimension(:), intent(in) :: filters
    integer, dimension(:,:), allocatable :: pairs

    integer :: i,j,l,m,n,npairs

    if( size(filters) == 1 ) then
       allocate(pairs(1,2))
       pairs = 1
       return
    end if

    n = 0
    do j = 1,size(filters)
       if( filters(j) == phsys%filter_ref ) n = j
    end do
    if( n == 0 ) stop 'Reference filter not included in filter set.'

    m = 0
    do j = 1, phsys%nbands
       if( phsys%filter(j) == phsys%filter_ref ) m = j
    end do
    if( m == 0 ) stop 'Reference filter not found in standard set.'

    allocate(pairs(size(filters),2))
    npairs = 1
    pairs(1,:) = n

    ! we are supposes that the filters are sorted by wavelength order
    n = pairs(1,1)
    do j = m-1,1,-1
       l = 0
       do i = 1, size(filters)
          if( phsys%filter(j) == filters(i) ) l = i
       end do
       if( l /= 0 ) then
          npairs = npairs + 1
          pairs(npairs,:) = (/l,n/)
          n = l
       end if
    end do
    n = pairs(1,1)
    do j = m+1, phsys%nbands
       l = 0
       do i = 1,size(filters)
          if( phsys%filter(j) == filters(i) ) l = i
       end do
       if( l /= 0 ) then
          npairs = npairs + 1
          pairs(npairs,:) = (/n,l/)
          n = l
       end if
    end do

  end subroutine phsyspairs

end module phsysfits
