!
!  photometric calibration
!
!  Copyright © 2012-15, 2017-8 F.Hroch (hroch@physics.muni.cz)
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


module calibre

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)

  type photores

     integer :: ndat
     real(dbl) :: ctph,dctph
     real(dbl), dimension(:), allocatable :: ra,dec,pht,cts,res

  end type photores


contains

  subroutine photores_init(phres,ndat)

    type(photores), intent(out) :: phres
    integer, intent(in) :: ndat

    allocate(phres%ra(ndat),phres%dec(ndat),phres%pht(ndat),phres%cts(ndat), &
         phres%res(ndat))
    phres%ndat = ndat
    phres%ctph = 0
    phres%dctph = 0

  end subroutine photores_init

  subroutine photores_destroy(phres)

    type(photores), intent(in out) :: phres

    deallocate(phres%ra,phres%dec,phres%pht,phres%cts,phres%res)
    phres%ndat = 0

  end subroutine photores_destroy


  subroutine syscal(ph,dph,ct,dct,mode,t,dt,sys,verbose)

    ! The calibration t = ph / ct (fmean) with the systematic error:
    ! mode < 0 (auto): the excess scatter of the relative residuals
    ! (ph - t*ct) / (t*ct) over their errors (exscatter), mode >= 0:
    ! the given value [mag]. When the systematic error is positive, it is
    ! added in quadrature to the relative errors of ct and t is fitted
    ! again; sys is the systematic error in magnitudes.

    use oakleaf
    use exscatter

    real(dbl), dimension(:), intent(in) :: ph,dph,ct,dct
    real(dbl), intent(in) :: mode
    real(dbl), intent(out) :: t,dt,sys
    logical, intent(in) :: verbose

    real(dbl), parameter :: mag = 1.0857362047581294_dbl  ! 2.5 / ln(10)
    real(dbl), dimension(:), allocatable :: r,e,dc
    real(dbl) :: srel
    integer :: iter

    call fmean(ph,dph,ct,dct,t,dt,verbose=verbose)
    sys = 0
    if( size(ph) < 3 .or. mode == 0 ) return

    allocate(r(size(ph)),e(size(ph)),dc(size(ph)))
    dc = dct
    do iter = 1, 3
       if( mode < 0 ) then
          r = (ph - t*ct) / (t*ct)
          e = sqrt(dph**2 + t**2*dct**2) / (t*ct)
          call excess_scatter(r,e,srel)
       else
          srel = mode / mag
       end if
       if( .not. (srel > 0) ) exit
       dc = sqrt(dct**2 + (srel*ct)**2)
       call fmean(ph,dph,ct,dc,t,dt,verbose=verbose)
       if( mode >= 0 ) exit
    end do
    sys = srel * mag
    if( verbose ) write(*,'(a,f0.5,a)') ' Systematic error: ',sys,' mag'

  end subroutine syscal


  subroutine caliber(pairs,filters,tr,tr1,pht,dpht,cts,dcts,ctph,dctph,phres,&
       verbose,syserr,sys)

    use fotran
    use oakleaf

    integer, dimension(:,:), intent(in) :: pairs
    character(len=*), dimension(:), intent(in) :: filters
    real(dbl), dimension(:,:), intent(in) :: tr,tr1
    real(dbl), dimension(:,:), intent(in) :: pht,dpht,cts,dcts
    real(dbl), dimension(:), intent(in out) :: ctph,dctph
    type(photores), dimension(:), intent(in out) :: phres
    logical, intent(in) :: verbose
    real(dbl), intent(in), optional :: syserr
    real(dbl), dimension(:), intent(out), optional :: sys

    real(dbl),dimension(:,:),allocatable :: ct, dct
    real(dbl),dimension(:,:),allocatable :: relerr
    real(dbl),dimension(:),allocatable :: c, dc
    real(dbl) :: mode, s
    integer :: i,n,ndat,npht

    ndat = size(pht,1)
    npht = size(pht,2)
    allocate(ct(ndat,npht),dct(ndat,npht),relerr(ndat,npht),c(ndat),dc(ndat))

    ! the catalogue photons transformed to the instrumental filters
    do i = 1, ndat
       call fotra(tr,pairs,pht(i,:),dpht(i,:),ct(i,:),dct(i,:))
    end do

    ! every filter as in calibr: the missing errors by Poisson
    mode = 0
    if( present(syserr) ) mode = syserr
    do n = 1, npht
       if( verbose ) write(*,*) 'Filter: ',trim(filters(n))
       c = cts(:,n)
       dc = dcts(:,n)
       where( ct(:,n) > 0 .and. dct(:,n) < epsilon(dct) )
          dct(:,n) = sqrt(ct(:,n))
       end where
       where( c > 0 .and. dc < epsilon(dc) )
          dc = sqrt(c)
       end where
       call syscal(ct(:,n),dct(:,n),c,dc,mode,ctph(n),dctph(n),s,verbose)
       if( present(sys) ) sys(n) = s
    end do

    ! the counts converted to the photons of the catalogue filters
    do i = 1, ndat
       call fotra(tr1,pairs,cts(i,:)*ctph,dcts(i,:),ct(i,:),dct(i,:))
    end do

    do n = 1, npht
       relerr(:,n) = (pht(:,n) - ct(:,n)) / ((pht(:,n) + ct(:,n)) / 2)
    end do

    phres(:)%ctph = ctph
    phres(:)%dctph = dctph

    do i = 1,npht
       phres(i)%pht = pht(:,i)
       phres(i)%cts = ct(:,i)
       phres(i)%res = relerr(:,i)
    end do

    deallocate(ct,dct,relerr,c,dc)

  end subroutine caliber


  subroutine calibr(pht,dpht,cts,dcts,ctph,dctph,phres,verbose,syserr,sys)

    ! syserr (optional) is the systematic error [mag], < 0 for its
    ! estimate by the residuals (see syscal), 0 by default;
    ! sys are the systematic errors applied [mag]

    use oakleaf

    real(dbl), dimension(:,:), intent(in) :: pht,dpht,cts,dcts
    real(dbl), dimension(:), intent(out) :: ctph,dctph
    type(photores), dimension(:), intent(in out) :: phres
    logical, intent(in) :: verbose
    real(dbl), intent(in), optional :: syserr
    real(dbl), dimension(:), intent(out), optional :: sys

    real(dbl),dimension(:),allocatable :: ct, dct,ph,dph
    real(dbl),dimension(:,:),allocatable :: relerr
    real(dbl) :: mode, s
    integer :: i,n,ndat,ncol

    if( size(pht,1) /= size(cts,1) ) stop 'Non-conformal data.'

    ndat = size(pht,1)
    ncol = size(pht,2)
    allocate(ct(ndat),dct(ndat),ph(ndat),dph(ndat),relerr(ndat,ncol))

    do n = 1, size(cts,2)

       ct(:) = cts(:,n)
       dct(:) = dcts(:,n)

       ph(:) = pht(:,n)
       dph(:) = dpht(:,n)

       where( ph > 0 .and. dph < epsilon(dph) )
          dph = sqrt(ph)
       end where

       where( ct > 0 .and. dct < epsilon(dct) )
          dct = sqrt(ct)
       end where

!       call phrate(ct,ph,ctph(n),dctph(n),verbose=verbose)
!       call phrate(ph,ct,ctph(n),dctph(n),verbose=verbose)
!       write(*,*) 'phrate result=',ctph(n),dctph(n)
       mode = 0
       if( present(syserr) ) mode = syserr
       call syscal(ph,dph,ct,dct,mode,ctph(n),dctph(n),s,verbose)
       if( present(sys) ) sys(n) = s
!       call fmean(ct,dct,ph,dph,ctph(n),dctph(n),verbose=verbose)
       relerr(:,n) = (ph - ctph(n)*ct) / ((ph + ctph(n)*ct) / 2)
!       relerr(:,n) = (ct - ctph(n)*ph) / ct !((ct + ctph(n)*ph) / 2)

    end do

    phres(:)%ctph = ctph
    phres(:)%dctph = dctph

    do i = 1, size(ctph)
       phres(i)%pht = pht(:,i)
       phres(i)%cts = cts(:,i)
       phres(i)%res = relerr(:,i)
    end do

    deallocate(ct,dct,ph,dph)

  end subroutine calibr


end module calibre
