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

    ! The calibration t = ph / ct with the systematic error:
    ! mode = 0: t by fmean, no systematic error;
    ! mode > 0: the given value [mag], added in quadrature to the relative
    ! errors of ct, t by fmean with these errors;
    ! mode < 0 (auto): t and the systematic error (the excess scatter of
    ! the stars over their errors) together by the likelihood of a mixture
    ! of the calibration stars and outliers (mixture_mean of exscatter):
    ! log t is the mean of log(ph/ct) weighted by the probabilities of
    ! the stars to be inliers and by their errors with the excess scatter,
    ! the excess scatter by the restricted likelihood of the inliers
    ! (unbiased for few stars), so outliers (misidentified stars, blends,
    ! variables) have no weight in either.
    ! sys is the systematic error in magnitudes.

    use oakleaf
    use exscatter

    real(dbl), dimension(:), intent(in) :: ph,dph,ct,dct
    real(dbl), intent(in) :: mode
    real(dbl), intent(out) :: t,dt,sys
    logical, intent(in) :: verbose

    real(dbl), parameter :: mag = 1.0857362047581294_dbl  ! 2.5 / ln(10)
    real(dbl), dimension(:), allocatable :: y,e,dc
    logical, dimension(:), allocatable :: valid
    real(dbl) :: srel, mu, dmu, fout

    call fmean(ph,dph,ct,dct,t,dt,verbose=verbose)
    sys = 0
    if( size(ph) < 3 .or. mode == 0 ) return

    if( mode < 0 ) then
       allocate(y(size(ph)),e(size(ph)),valid(size(ph)))
       valid = ph > 0 .and. ct > 0 .and. dph >= 0 .and. dct >= 0
       where( valid )
          y = log(ph / ct)
          e = sqrt((dph/ph)**2 + (dct/ct)**2)
       elsewhere
          y = 0
          e = -1
       end where
       if( count(valid .and. e > 0) < 3 ) return
       call mixture_mean(y,e,mu,dmu,srel,fout)
       t = exp(mu)
       dt = t * dmu
       sys = srel * mag
       if( verbose ) write(*,'(a,f0.5,a,f0.3,a,1pg0.6,a,1pg0.2)') &
            ' Systematic error: ',sys,' mag, outlier fraction: ',fout, &
            ', calibration: ',t,' +- ',dt
       return
    end if

    allocate(dc(size(ph)))
    srel = mode / mag
    dc = sqrt(dct**2 + (srel*ct)**2)
    call fmean(ph,dph,ct,dc,t,dt,verbose=verbose)
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


  subroutine calibr_colour(pht,dpht,cts,dcts,col,dcol,cyc,ctm,ctph,dctph,phres,&
       verbose,syserr,sys)

    ! The calibration of every band (column) with a colour term (colterm):
    ! ln(ph/ct) = ln ctph + b1 (c - cref) + b2 (c - cref)**2. col are
    ! the colours of the stars, dcol their errors (< 0: no colour, the star
    ! is not used), cyc the covariance of ln(ph) and of the colour
    ! (the colour contains the catalogue magnitude of the band).
    ! ctm%model = COLT_CAT: the coefficients are fitted with ctph (of
    ! the order ctm%order), cref is the mean colour of the calibration
    ! stars when ctm%cref = cref_mean; COLT_TRA: the coefficients
    ! ctm%beta(1:2,:) and ctm%cref are fixed, ctph is fitted.
    ! syserr as in calibr (< 0: the systematic error fitted by the mixture
    ! likelihood, >= 0: the value [mag], fixed); sys the values applied.
    ! On output, ctm holds the coefficients, their covariance (ln ctph
    ! included), the mean colour, the mean colour correction of
    ! the calibration stars and their scatter.

    use colterm
    use exscatter

    real(dbl), dimension(:,:), intent(in) :: pht,dpht,cts,dcts,col,dcol,cyc
    type(colour_term), intent(in out) :: ctm
    real(dbl), dimension(:), intent(out) :: ctph,dctph
    type(photores), dimension(:), intent(in out) :: phres
    logical, intent(in) :: verbose
    real(dbl), intent(in) :: syserr
    real(dbl), dimension(:), intent(out) :: sys

    real(dbl), dimension(:), allocatable :: ct,dct,ph,dph,y,e,x,ex,cx,q,g,gp,w
    real(dbl), dimension(0:2) :: b
    real(dbl), dimension(3,3) :: cv, jm
    real(dbl) :: fixsig, sig, fout, cbar, d, sw
    integer :: n, ndat, i
    logical, dimension(:), allocatable :: valid

    ndat = size(pht,1)
    allocate(ct(ndat),dct(ndat),ph(ndat),dph(ndat),y(ndat),e(ndat),x(ndat), &
         ex(ndat),cx(ndat),q(ndat),g(ndat),gp(ndat),w(ndat),valid(ndat))

    fixsig = -1
    if( syserr >= 0 ) fixsig = syserr / magln

    do n = 1, size(cts,2)

       ct = cts(:,n)
       dct = dcts(:,n)
       ph = pht(:,n)
       dph = dpht(:,n)
       where( ph > 0 .and. dph < epsilon(dph) )
          dph = sqrt(ph)
       end where
       where( ct > 0 .and. dct < epsilon(dct) )
          dct = sqrt(ct)
       end where

       valid = ph > 0 .and. ct > 0 .and. dcol(:,n) >= 0
       where( valid )
          y = log(ph / ct)
          e = sqrt((dph/ph)**2 + (dct/ct)**2)
          ex = dcol(:,n)
          cx = cyc(:,n)
       elsewhere
          y = 0
          e = -1
          ex = -1
          cx = 0
       end where
       ctm%nstars(n) = count(valid)
       if( ctm%nstars(n) == 0 ) stop 'No calibration star with a colour.'

       if( ctm%model == COLT_CAT ) then

          ! the colours centred for the numerics, then shifted to cref
          cbar = sum(col(:,n),mask=valid) / count(valid)
          x = col(:,n) - cbar
          call mixture_fit(y,e,x,ex,cx,ctm%order(n),fixsig,b,cv,sig,fout,q)
          if( ctm%cref(n) == cref_mean ) then
             w = q / (e**2 + sig**2)
             sw = sum(w,mask=valid)
             if( sw > 0 ) then
                ctm%cref(n) = sum(w*col(:,n),mask=valid) / sw
             else
                ctm%cref(n) = cbar
             end if
          end if
          ! b(c - cbar) = b'(c - cref): d = cref - cbar
          d = ctm%cref(n) - cbar
          jm = 0
          jm(1,:) = [ 1.0_dbl, d, d**2 ]
          jm(2,:) = [ 0.0_dbl, 1.0_dbl, 2*d ]
          jm(3,:) = [ 0.0_dbl, 0.0_dbl, 1.0_dbl ]
          ctm%beta(:,n) = matmul(jm,b)
          ctm%cov(:,:,n) = matmul(jm,matmul(cv,transpose(jm)))

       else

          ! the coefficients fixed: the zero point of y - g(c)
          x = col(:,n) - ctm%cref(n)
          g = ctm%beta(1,n)*x + ctm%beta(2,n)*x**2
          gp = ctm%beta(1,n) + 2*ctm%beta(2,n)*x
          where( valid )
             y = y - g
             e = sqrt(max(e**2 + gp**2*ex**2 - 2*gp*cx, 0.01_dbl*e**2))
          end where
          ex = 0
          cx = 0
          call mixture_fit(y,e,x,ex,cx,0,fixsig,b,cv,sig,fout,q)
          ctm%beta(0,n) = b(0)
          ! the coefficients of the table: their errors only
          ctm%cov(1,1,n) = cv(1,1)

       end if

       ctph(n) = exp(ctm%beta(0,n))
       dctph(n) = ctph(n) * sqrt(ctm%cov(1,1,n))
       sys(n) = sig * magln

       ! the mean colour and the mean correction of the calibration stars
       x = col(:,n) - ctm%cref(n)
       g = ctm%beta(1,n)*x + ctm%beta(2,n)*x**2
       sw = sum(q,mask=valid)
       if( sw > 0 ) then
          ctm%cmean(n) = sum(q*col(:,n),mask=valid) / sw
          ctm%cstd(n) = sqrt(sum(q*(col(:,n) - ctm%cmean(n))**2,mask=valid)/sw)
          ctm%gmean(n) = sum(q*g,mask=valid) / sw
          ctm%gstd(n) = sqrt(sum(q*(g - ctm%gmean(n))**2,mask=valid) / sw)
       end if

       if( verbose ) then
          write(*,'(a,i0,a,i0,a)') ' Colour term, band ',n,' (',ctm%nstars(n), &
               ' stars with a colour):'
          write(*,'(a,f0.4,a,f0.4,a,f0.4)') '  reference colour: ',ctm%cref(n), &
               ', mean colour: ',ctm%cmean(n),' +- ',ctm%cstd(n)
          write(*,'(a,f0.5,a,f0.5,a)') '  k1 = ',-magln*ctm%beta(1,n),' +- ', &
               magln*sqrt(max(ctm%cov(2,2,n),0.0_dbl)),' mag/mag'
          if( ctm%order(n) >= 2 .or. ctm%beta(2,n) /= 0 ) &
               write(*,'(a,f0.5,a,f0.5,a)') '  k2 = ',-magln*ctm%beta(2,n), &
               ' +- ',magln*sqrt(max(ctm%cov(3,3,n),0.0_dbl)),' mag/mag**2'
          write(*,'(a,f0.5,a,f0.3,a,1pg0.6,a,1pg0.2)') '  systematic error: ', &
               sys(n),' mag, outlier fraction: ',fout,', calibration: ', &
               ctph(n),' +- ',dctph(n)
       end if

       ! the residuals by the colour of the stars (the mean correction
       ! for the stars without a colour); the counts are corrected
       ! to the reference colour
       do i = 1, ndat
          if( .not. valid(i) ) g(i) = ctm%gmean(n)
       end do
       ct = cts(:,n) * exp(g)
       phres(n)%pht = pht(:,n)
       phres(n)%cts = ct
       phres(n)%res = (pht(:,n) - ctph(n)*ct) / ((pht(:,n) + ctph(n)*ct) / 2)

    end do

    phres(:)%ctph = ctph
    phres(:)%dctph = dctph

  end subroutine calibr_colour


end module calibre
