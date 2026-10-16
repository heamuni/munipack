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

  ! a systematic error estimated from fewer stars is uncertain (a warning)
  integer, parameter :: nsyswarn = 10

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


  subroutine syscal(ph,dph,ct,dct,mode,t,dt,sys,verbose,sysint,nsys)

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
    ! sys is the systematic error in magnitudes; sysint (optional) its
    ! interval [mag] by the profile of the restricted likelihood (the lower
    ! bound zero for an upper limit; sys for a fixed value, -1 when it
    ! is not estimated), nsys (optional) the number of the stars used.

    use oakleaf
    use exscatter

    real(dbl), dimension(:), intent(in) :: ph,dph,ct,dct
    real(dbl), intent(in) :: mode
    real(dbl), intent(out) :: t,dt,sys
    logical, intent(in) :: verbose
    real(dbl), dimension(2), intent(out), optional :: sysint
    integer, intent(out), optional :: nsys

    real(dbl), parameter :: mag = 1.0857362047581294_dbl  ! 2.5 / ln(10)
    real(dbl), dimension(:), allocatable :: y,e,dc
    logical, dimension(:), allocatable :: valid
    real(dbl), dimension(2) :: sint
    real(dbl) :: srel, mu, dmu, fout
    integer :: nf

    call fmean(ph,dph,ct,dct,t,dt,verbose=verbose)
    sys = 0
    if( present(sysint) ) sysint = -1
    if( present(nsys) ) nsys = count(ph > 0 .and. ct > 0)
    if( mode == 0 ) then
       if( present(sysint) ) sysint = 0
       return
    end if
    if( size(ph) < 3 ) return

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
       if( present(nsys) ) nsys = count(valid .and. e > 0)
       if( count(valid .and. e > 0) < 3 ) return
       call mixture_mean(y,e,mu,dmu,srel,fout,sigint=sint,nfit=nf)
       t = exp(mu)
       dt = t * dmu
       sys = srel * mag
       if( present(sysint) ) then
          sysint = sint * mag
          if( sint(2) < 0 ) sysint = -1
       end if
       if( present(nsys) ) nsys = nf
       if( verbose ) write(*,'(a,f0.5,a,f0.5,a,f0.5,a,i0,a,f0.3,a,1pg0.6,a,1pg0.2)') &
            ' Systematic error: ',sys,' mag (',mag*sint(1),' .. ',mag*sint(2), &
            ', ',nf,' stars), outlier fraction: ',fout, &
            ', calibration: ',t,' +- ',dt
       return
    end if

    allocate(dc(size(ph)))
    srel = mode / mag
    dc = sqrt(dct**2 + (srel*ct)**2)
    call fmean(ph,dph,ct,dc,t,dt,verbose=verbose)
    sys = srel * mag
    if( present(sysint) ) sysint = sys
    if( verbose ) write(*,'(a,f0.5,a)') ' Systematic error: ',sys,' mag'

  end subroutine syscal


  subroutine caliber(pairs,filters,tr,tr1,pht,dpht,cts,dcts,ctph,dctph,phres,&
       verbose,syserr,sys,sysint,nsys)

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
    real(dbl), dimension(:,:), intent(out), optional :: sysint
    integer, dimension(:), intent(out), optional :: nsys

    real(dbl),dimension(:,:),allocatable :: ct, dct
    real(dbl),dimension(:,:),allocatable :: relerr
    real(dbl),dimension(:),allocatable :: c, dc
    real(dbl), dimension(2) :: si
    real(dbl) :: mode, s
    integer :: i,n,ndat,npht,ns

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
       call syscal(ct(:,n),dct(:,n),c,dc,mode,ctph(n),dctph(n),s,verbose,si,ns)
       if( present(sys) ) sys(n) = s
       if( present(sysint) ) sysint(:,n) = si
       if( present(nsys) ) nsys(n) = ns
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


  subroutine calibr(pht,dpht,cts,dcts,ctph,dctph,phres,verbose,syserr,sys, &
       sysint,nsys)

    ! syserr (optional) is the systematic error [mag], < 0 for its
    ! estimate by the residuals (see syscal), 0 by default;
    ! sys are the systematic errors applied [mag], sysint their intervals
    ! and nsys the numbers of the stars (see syscal)

    use oakleaf

    real(dbl), dimension(:,:), intent(in) :: pht,dpht,cts,dcts
    real(dbl), dimension(:), intent(out) :: ctph,dctph
    type(photores), dimension(:), intent(in out) :: phres
    logical, intent(in) :: verbose
    real(dbl), intent(in), optional :: syserr
    real(dbl), dimension(:), intent(out), optional :: sys
    real(dbl), dimension(:,:), intent(out), optional :: sysint
    integer, dimension(:), intent(out), optional :: nsys

    real(dbl),dimension(:),allocatable :: ct, dct,ph,dph
    real(dbl),dimension(:,:),allocatable :: relerr
    real(dbl), dimension(2) :: si
    real(dbl) :: mode, s
    integer :: i,n,ndat,ncol,ns

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
       call syscal(ph,dph,ct,dct,mode,ctph(n),dctph(n),s,verbose,si,ns)
       if( present(sys) ) sys(n) = s
       if( present(sysint) ) sysint(:,n) = si
       if( present(nsys) ) nsys(n) = ns
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
       verbose,syserr,sys,sysint,nsys,ccv)

    ! The calibration of every band (column) with colour terms (colterm):
    ! ln(ph/ct) = ln ctph + sum_j b1_j (c_j - cref_j) + b2_j (c_j - cref_j)**2.
    ! col(:,n,j) are the colours j of the stars for the band n, dcol their
    ! errors (< 0: no colour, the star is not used), cyc the covariance of
    ! ln(ph) and of the colour (the colour contains the catalogue magnitude
    ! of the band), ccv (optional, several colours) the covariances of
    ! the colours j /= k of a star, (:,j,k).
    ! ctm%model = COLT_CAT: the coefficients are fitted with ctph (of
    ! the orders ctm%order), cref is the mean colour of the calibration
    ! stars when ctm%cref = cref_mean; COLT_TRA (one colour):
    ! the coefficients ctm%beta(1:2,:) and ctm%cref are fixed, ctph is
    ! fitted.
    ! syserr as in calibr (< 0: the systematic error fitted by the mixture
    ! likelihood, >= 0: the value [mag], fixed); sys the values applied,
    ! sysint (optional) their intervals, nsys (optional) the numbers of
    ! the stars (see syscal).
    ! On output, ctm holds the coefficients, their covariance (ln ctph
    ! included), the mean colours, the mean colour corrections of
    ! the calibration stars and their scatter.

    use colterm
    use exscatter

    real(dbl), dimension(:,:), intent(in) :: pht,dpht,cts,dcts
    real(dbl), dimension(:,:,:), intent(in) :: col,dcol,cyc
    type(colour_term), intent(in out) :: ctm
    real(dbl), dimension(:), intent(out) :: ctph,dctph
    type(photores), dimension(:), intent(in out) :: phres
    logical, intent(in) :: verbose
    real(dbl), intent(in) :: syserr
    real(dbl), dimension(:), intent(out) :: sys
    real(dbl), dimension(:,:), intent(out), optional :: sysint
    integer, dimension(:), intent(out), optional :: nsys
    real(dbl), dimension(:,:,:), intent(in), optional :: ccv

    real(dbl), dimension(:), allocatable :: ct,dct,ph,dph,y,e,x,ex,cx,q,g,gp,w,bm,cbm
    real(dbl), dimension(:,:), allocatable :: xm,exm,cxm,gj,jmm,cvm
    real(dbl), dimension(:,:,:), allocatable :: exx
    real(dbl), dimension(0:2) :: b
    real(dbl), dimension(3,3) :: cv, jm
    real(dbl), dimension(2) :: si
    real(dbl) :: fixsig, sig, fout, cbar, d, sw
    integer :: n, ndat, i, j, k, nf, nt
    logical, dimension(:), allocatable :: valid
    logical, dimension(:,:), allocatable :: known

    ndat = size(pht,1)
    nt = ctm%nterm
    allocate(ct(ndat),dct(ndat),ph(ndat),dph(ndat),y(ndat),e(ndat),x(ndat), &
         ex(ndat),cx(ndat),q(ndat),g(ndat),gp(ndat),w(ndat),valid(ndat))
    if( nt > 1 ) allocate(xm(ndat,nt),exm(ndat,nt),cxm(ndat,nt), &
         exx(ndat,nt,nt),gj(ndat,nt),known(ndat,nt),bm(0:2*nt),cbm(nt), &
         jmm(2*nt+1,2*nt+1),cvm(2*nt+1,2*nt+1))

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

       valid = ph > 0 .and. ct > 0
       do j = 1, nt
          valid = valid .and. dcol(:,n,j) >= 0
       end do
       where( valid )
          y = log(ph / ct)
          e = sqrt((dph/ph)**2 + (dct/ct)**2)
          ex = dcol(:,n,1)
          cx = cyc(:,n,1)
       elsewhere
          y = 0
          e = -1
          ex = -1
          cx = 0
       end where
       ctm%nstars(n) = count(valid)
       if( ctm%nstars(n) == 0 ) then
          if( nt > 1 ) stop 'No calibration star with all the colours.'
          stop 'No calibration star with a colour.'
       end if

       if( ctm%model == COLT_CAT .and. nt == 1 ) then

          ! the colours centred for the numerics, then shifted to cref
          cbar = sum(col(:,n,1),mask=valid) / count(valid)
          x = col(:,n,1) - cbar
          call mixture_fit(y,e,x,ex,cx,ctm%order(1,n),fixsig,b,cv,sig,fout,q, &
               sigint=si,nfit=nf)
          if( ctm%cref(1,n) == cref_mean ) then
             w = q / (e**2 + sig**2)
             sw = sum(w,mask=valid)
             if( sw > 0 ) then
                ctm%cref(1,n) = sum(w*col(:,n,1),mask=valid) / sw
             else
                ctm%cref(1,n) = cbar
             end if
          end if
          ! b(c - cbar) = b'(c - cref): d = cref - cbar
          d = ctm%cref(1,n) - cbar
          jm = 0
          jm(1,:) = [ 1.0_dbl, d, d**2 ]
          jm(2,:) = [ 0.0_dbl, 1.0_dbl, 2*d ]
          jm(3,:) = [ 0.0_dbl, 0.0_dbl, 1.0_dbl ]
          ctm%beta(:,n) = matmul(jm,b)
          ctm%cov(:,:,n) = matmul(jm,matmul(cv,transpose(jm)))

       else if( ctm%model == COLT_CAT ) then

          ! several colours: centred for the numerics, then shifted to
          ! the references; the Jacobian of the shift of every colour
          ! as above
          jmm = 0
          jmm(1,1) = 1
          do j = 1, nt
             cbar = sum(col(:,n,j),mask=valid) / count(valid)
             where( valid )
                xm(:,j) = col(:,n,j) - cbar
                exm(:,j) = dcol(:,n,j)
                cxm(:,j) = cyc(:,n,j)
             elsewhere
                xm(:,j) = 0
                exm(:,j) = -1
                cxm(:,j) = 0
             end where
             do k = 1, nt
                exx(:,j,k) = 0
                if( present(ccv) .and. j /= k ) then
                   where( valid ) exx(:,j,k) = ccv(:,j,k)
                end if
             end do
             cbm(j) = cbar
          end do
          call mixture_fit_multi(y,e,xm,exm,exx,cxm,ctm%order(:,n),fixsig, &
               bm,cvm,sig,fout,q,sigint=si,nfit=nf)
          w = q / (e**2 + sig**2)
          sw = sum(w,mask=valid)
          do j = 1, nt
             cbar = cbm(j)
             if( ctm%cref(j,n) == cref_mean ) then
                if( sw > 0 ) then
                   ctm%cref(j,n) = sum(w*col(:,n,j),mask=valid) / sw
                else
                   ctm%cref(j,n) = cbar
                end if
             end if
             d = ctm%cref(j,n) - cbar
             jmm(1,2*j) = d
             jmm(1,2*j+1) = d**2
             jmm(2*j,2*j) = 1
             jmm(2*j,2*j+1) = 2*d
             jmm(2*j+1,2*j+1) = 1
          end do
          ctm%beta(:,n) = matmul(jmm,bm)
          ctm%cov(:,:,n) = matmul(jmm,matmul(cvm,transpose(jmm)))

       else

          ! the coefficients fixed: the zero point of y - g(c)
          x = col(:,n,1) - ctm%cref(1,n)
          g = ctm%beta(1,n)*x + ctm%beta(2,n)*x**2
          gp = ctm%beta(1,n) + 2*ctm%beta(2,n)*x
          where( valid )
             y = y - g
             e = sqrt(max(e**2 + gp**2*ex**2 - 2*gp*cx, 0.01_dbl*e**2))
          end where
          ex = 0
          cx = 0
          call mixture_fit(y,e,x,ex,cx,0,fixsig,b,cv,sig,fout,q,sigint=si,nfit=nf)
          ctm%beta(0,n) = b(0)
          ! the coefficients of the table: their errors only
          ctm%cov(1,1,n) = cv(1,1)

       end if

       ctph(n) = exp(ctm%beta(0,n))
       dctph(n) = ctph(n) * sqrt(ctm%cov(1,1,n))
       sys(n) = sig * magln
       if( present(sysint) ) then
          sysint(:,n) = si * magln
          if( si(2) < 0 ) sysint(:,n) = -1
          if( syserr == 0 ) sysint(:,n) = 0
       end if
       if( present(nsys) ) nsys(n) = nf

       ! the mean colour and the mean correction of the calibration stars
       if( nt == 1 ) then
          x = col(:,n,1) - ctm%cref(1,n)
          g = ctm%beta(1,n)*x + ctm%beta(2,n)*x**2
          sw = sum(q,mask=valid)
          if( sw > 0 ) then
             ctm%cmean(1,n) = sum(q*col(:,n,1),mask=valid) / sw
             ctm%cstd(1,n) = sqrt(sum(q*(col(:,n,1) - ctm%cmean(1,n))**2, &
                  mask=valid)/sw)
             ctm%gmean(1,n) = sum(q*g,mask=valid) / sw
             ctm%gstd(1,n) = sqrt(sum(q*(g - ctm%gmean(1,n))**2,mask=valid) / sw)
          end if
       else
          ! of every colour, the covariances of the corrections
          do j = 1, nt
             x = col(:,n,j) - ctm%cref(j,n)
             gj(:,j) = ctm%beta(2*j-1,n)*x + ctm%beta(2*j,n)*x**2
             known(:,j) = dcol(:,n,j) >= 0
          end do
          sw = sum(q,mask=valid)
          if( sw > 0 ) then
             do j = 1, nt
                ctm%cmean(j,n) = sum(q*col(:,n,j),mask=valid) / sw
                ctm%cstd(j,n) = sqrt(sum(q*(col(:,n,j) - ctm%cmean(j,n))**2, &
                     mask=valid)/sw)
                ctm%gmean(j,n) = sum(q*gj(:,j),mask=valid) / sw
                ctm%gstd(j,n) = sqrt(sum(q*(gj(:,j) - ctm%gmean(j,n))**2, &
                     mask=valid) / sw)
             end do
             do j = 1, nt
                do k = 1, nt
                   ctm%ccs(j,k,n) = sum(q*(col(:,n,j) - ctm%cmean(j,n))* &
                        (col(:,n,k) - ctm%cmean(k,n)),mask=valid) / sw
                end do
             end do
          end if
          ! the corrections of the stars: the missing colours predicted
          ! by the known ones (colour_corr)
          do i = 1, ndat
             if( valid(i) ) then
                g(i) = sum(gj(i,:))
             else
                call colour_corr(ctm,n,col(i,n,:),dcol(i,n,:),g(i),d)
             end if
          end do
       end if

       if( verbose ) then
          if( nt == 1 ) then
             write(*,'(a,i0,a,i0,a)') ' Colour term, band ',n,' (',ctm%nstars(n), &
                  ' stars with a colour):'
          else
             write(*,'(a,i0,a,i0,a)') ' Colour terms, band ',n,' (',ctm%nstars(n), &
                  ' stars with all the colours):'
          end if
          do j = 1, nt
             if( nt > 1 ) write(*,'(a,i0,4a)') '  colour ',j,': ', &
                  trim(ctm%cola(j,n)),' - ',trim(ctm%colb(j,n))
             write(*,'(a,f0.4,a,f0.4,a,f0.4)') '  reference colour: ',ctm%cref(j,n), &
                  ', mean colour: ',ctm%cmean(j,n),' +- ',ctm%cstd(j,n)
             write(*,'(a,f0.5,a,f0.5,a)') '  k1 = ',-magln*ctm%beta(2*j-1,n),' +- ', &
                  magln*sqrt(max(ctm%cov(2*j,2*j,n),0.0_dbl)),' mag/mag'
             if( ctm%order(j,n) >= 2 .or. ctm%beta(2*j,n) /= 0 ) &
                  write(*,'(a,f0.5,a,f0.5,a)') '  k2 = ',-magln*ctm%beta(2*j,n), &
                  ' +- ',magln*sqrt(max(ctm%cov(2*j+1,2*j+1,n),0.0_dbl)),' mag/mag**2'
          end do
          write(*,'(a,f0.5,a,f0.5,a,f0.5,a,i0,a,f0.3,a,1pg0.6,a,1pg0.2)') &
               '  systematic error: ',sys(n),' mag (',magln*si(1),' .. ', &
               magln*si(2),', ',nf,' stars), outlier fraction: ',fout, &
               ', calibration: ',ctph(n),' +- ',dctph(n)
       end if

       ! the residuals by the colour of the stars (the mean correction
       ! for the stars without a colour); the counts are corrected
       ! to the reference colour
       if( nt == 1 ) then
          do i = 1, ndat
             if( .not. valid(i) ) g(i) = ctm%gmean(1,n)
          end do
       end if
       ct = cts(:,n) * exp(g)
       phres(n)%pht = pht(:,n)
       phres(n)%cts = ct
       phres(n)%res = (pht(:,n) - ctph(n)*ct) / ((pht(:,n) + ctph(n)*ct) / 2)

    end do

    phres(:)%ctph = ctph
    phres(:)%dctph = dctph

  end subroutine calibr_colour


end module calibre
