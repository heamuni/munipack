!
!  An implementation of growth-curves in aperture photometry
!  of point-like objects.
!
!     This piece of software is dedicated to my colleague
!     Vladimír Štefl. because the growth-curve method is his
!     favourite method for spectral line analysis. :)
!
!  Copyright © 2015-7 F.Hroch (hroch@physics.muni.cz)
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

! grow curve construction:
!   * localy (en every point), the variartion principle
!     must be satisfied, base funtions be linear?


module grow_curve

  use grow_report

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: pi = 3.14159265358979312_dbl

  type(grow_reporter), pointer, private :: reporter
  logical, private :: report

!  private :: growcurve, select_valid_isolated, curve_init, grow_deepest, &
  private :: growcurve, select_valid_isolated, curve_init, model_curve, &
       curve_spline, grow_phot, estim_rflux90, estim_hwhm

contains


  subroutine growphot(xcens,ycens,apcts,apcts_err,sky,sky_sig,raper,hwhm,sep,threshold, &
       gcount,gcount_err,skycorr,skyerrcorr,curve,curve_err,growflag,prof, &
       rflux90,ghwhm,make_report,xreport)

    real, dimension(:), intent(in) :: xcens,ycens
    real, dimension(:,:), intent(in) :: apcts,apcts_err
    real, dimension(:),intent(in) :: sky,sky_sig,raper
    real, intent(in) :: hwhm, sep, threshold
    real(dbl), dimension(:), intent(out) :: gcount,gcount_err,curve,curve_err,prof
    real, dimension(:), intent(out) :: skycorr,skyerrcorr
    integer, dimension(:), intent(out) :: growflag
    real, intent(out) :: rflux90, ghwhm
    logical, intent(in) :: make_report
    type(grow_reporter), target, intent(in out) :: xreport

    integer :: i,j,n

    report = make_report
    reporter => xreport

    n = size(gcount)

    ! construct growth-curve
    call growcurve(xcens,ycens,apcts,apcts_err,sky,sky_sig,raper,hwhm,sep,threshold, &
         curve,curve_err,growflag,prof,rflux90,ghwhm)

    ! photometry by using of the curve
    call grow_phot(raper,apcts,apcts_err,curve,curve_err,sky_sig,hwhm, &
         gcount,gcount_err,skycorr,skyerrcorr)

    if( report ) then

       call curve_spline(raper,curve,curve_err)

       n = 0
       do i = 1, size(apcts,1)
          do j = 1, size(apcts,2)
             n = n + 1
             if( apcts(i,j) > 0 .and. curve(j) > 0 .and. gcount(i) > 0 ) then
                reporter%grow(i,j) = (apcts(i,j) - skycorr(i)*pi*raper(j)**2)/ gcount(i)
                reporter%resgrow(i,j) =  reporter%grow(i,j) - curve(j)
                reporter%grow_err(i,j) = reporter%grow(i,j)* &
                     sqrt((apcts_err(i,j)/apcts(i,j))**2 + (curve_err(j)/curve(j))**2)
             else
                reporter%resgrow(i,j) = -1
                reporter%grow(i,j) = -1
                reporter%grow_err(i,j) = -1
             end if
          end do
       end do
    end if

  end subroutine growphot


  subroutine growcurve(xcens,ycens,apcts,apcts_err,sky,sky_err,raper,hwhm,sep, &
       threshold,curve,curve_err,growflag,prof,rflux90,ghwhm)

    use grow_model
    use grow_fit

    real, dimension(:), intent(in) :: xcens,ycens
    real, dimension(:,:), intent(in) :: apcts,apcts_err
    real, dimension(:),intent(in) :: sky,sky_err,raper
    real, intent(in) :: hwhm, sep, threshold
    real(dbl), dimension(:), intent(out) :: curve,curve_err, prof
    integer, dimension(:), intent(out) :: growflag
    real, intent(out) :: rflux90, ghwhm

    real(dbl), dimension(:,:), allocatable :: grow,grow_err!,grows,grows_err
    integer, dimension(:), allocatable :: id, idn
    integer :: nid,nrows,naper,i,nhwhm,maper
!    real(dbl) :: y1,y2,x2,x1,a,b,gx

    nrows = size(apcts,1)
    naper = size(apcts,2)
    growflag = 0
    nhwhm = optimal_index(raper,hwhm)

    allocate(id(nrows),idn(nrows),grow(nrows,naper),grow_err(nrows,naper))

    ! generate initial growth-curve estimate on base of Gaussian
    ! the estimate of HWHM is important (especialy when
    ! there are no valid stars to determine the accurate one).
    call curve_init(raper,hwhm,curve,curve_err,prof)

!    do i = 1,12
!       write(*,*) raper(i),curve(i)
!    end do
!    stop 0

    ! select aperture optimal for star selection procedure
    maper = 0
    do i = 1,naper
       if( raper(i) >= 2*hwhm ) then
          maper = i
          exit
       end if
    end do
    if( maper == 0 ) goto 666

    ! select isolated stars with all aperture measurements valid
    call select_valid_isolated(xcens,ycens,sep,threshold,apcts,apcts_err,sky, &
         raper(maper),maper,nid,id)
!    write(*,*) nid,maper
    if( nid == 0 ) goto 666

    ! ** WARNING **
    ! The order in apcts,apcts_err field is unchanged.
    ! All valid grow(_err) are spreaded over all the field keeping their ranks.
    ! It makes possibility to identify any object at any time of processing.

    ! flags for stars with an accepted profile and objects with non-stellar ones
    forall( i = 1:nid  ) growflag(id(i))  = 1

    ! construct the principial growth curve
    call model_curve(nid,id,hwhm,nhwhm,raper,apcts,apcts_err,sky,curve,curve_err)

!!$    ! construct growth-curves of individual (selected) stars
!!$    call grow_deepest(nid,id,hwhm,nhwhm,raper,apcts,apcts_err,sky_err,grow,grow_err)
!!$!    write(*,*) nid
!!$    if( nid == 0 ) goto 666
!!$
!!$    ! ** WARNING **
!!$    ! The order in apcts,apcts_err field is unchanged.
!!$    ! All valid grow(_err) are spreaded over all the field keeping their ranks.
!!$    ! It makes possibility to identify any object at any time of processing.
!!$
!!$    ! flags for stars with an accepted profile and objects with non-stellar ones
!!$    forall( i = 1:nid  ) growflag(id(i))  = 1
!!$
!!$    allocate(grows(nid,naper),grows_err(nid,naper))
!!$    do i = 1,nid
!!$       k = id(i)
!!$       grows(i,:) = grow(k,:)
!!$       grows_err(i,:) = grow_err(k,:)
!!$    end do
!!$
!!$    call growmodel1(grows,grows_err,curve,curve_err)
!!$    goto 33
!!$
!!$    ! correct for finite radius
!!$    y1 = 1/(1 - curve(naper))
!!$    y2 = 1/(1 - curve((nhwhm+naper)/2))
!!$    x1 = raper(naper)
!!$    x2 = raper((nhwhm+naper)/2)
!!$    a = (y2 - y1)/(x2 - x1)
!!$    b = (y1 + y2 - a*(x1 + x2)) / 2
!!$    gx = 1 - 1/(a*raper(naper) + b)
!!$    write(*,*) a,b,gx
!!$
!!$    call grow_deepest(nid,id,hwhm,nhwhm,raper,apcts,apcts_err,sky_err,grow,grow_err,a,b)
!!$    call growmodel1(grows,grows_err,curve,curve_err)
!!$
!!$    do i = 1,size(curve)
!!$       write(*,*) raper(i),curve(i) - 1,1/(1-curve(i))
!!$    end do
!!$!    curve = curve / gx
!!$33  continue
!!$
!!$    call grow_deepest1(nid,id,hwhm,nhwhm,raper,apcts,apcts_err,curve,curve_err)

!    deallocate(grows,grows_err)

666 continue

    ! estimate HWHM and energy within a radius as an indicator of actual spread
    rflux90 = estim_rflux90(raper,curve)

    ! profile is derived from rectangle rule: the mean surface brightness
    ! in the central disc (area pi*r**2) and in the annuli
    prof(1) = curve(1)/raper(1)**2/pi
    do i = 2,naper
       prof(i) = (curve(i) - curve(i-1)) / (pi * (raper(i)**2 - raper(i-1)**2))
       prof(i) = max(prof(i),0.0)
    end do
    ghwhm = estim_hwhm(raper,prof)

!    do i = 1,naper
!       write(*,'(f5.1,3f10.4)') raper(i),curve(i),curve_err(i),prof(i)
!    end do

    deallocate(id,idn,grow,grow_err)

  end subroutine growcurve


  subroutine select_valid_isolated(xcens,ycens,sep,threshold,apcts,apcts_err,&
       sky,raper,maper,nid,id)

    ! select isolated stars suitable for estimate of growth-curve
    !
    ! the selection is on base of rules:
    !   * isolated star by 'sep' parameters
    !   * all measurements are valid
    !   * bright stars

    real, dimension(:), intent(in) :: xcens,ycens,sky
    real, intent(in) :: sep, threshold, raper
    real, dimension(:,:), intent(in) :: apcts, apcts_err
    integer, intent(in) :: maper
    integer, intent(out) :: nid
    integer, dimension(:), intent(out) :: id

    integer :: n,i,j,naper
    logical :: found
    real :: r,aper

    n = size(xcens)
    naper = size(apcts,2)
    aper = 3.14*raper**2

    ! selection of isolated stars
    nid = 0
    do i = 1, n
       found = .false.
       do j = 1, n
          if( i /= j ) then
             r = sqrt((xcens(i) - xcens(j))**2 + (ycens(i) - ycens(j))**2)
             if( r < sep ) then
                found = .true.
                goto 90
             end if
          end if
       end do
90     continue

       ! accept only valid bright stars highly above sky limit
!       write(*,*) apcts(i,maper)/(aper*sky(i)),threshold,found
       if( .not. found .and. &
            all(apcts(i,:) > 0) .and. all(apcts_err(i,:) > 0) .and. &
            apcts(i,maper) / (aper*sky(i))  > threshold ) then

          nid = nid + 1
          id(nid) = i

       end if
    end do

  end subroutine select_valid_isolated

  subroutine curve_init(r,hwhm,curve,curve_err,prof)

    real, dimension(:), intent(in) :: r
    real, intent(in) :: hwhm
    real(dbl), dimension(:), intent(out) :: curve, curve_err, prof

    prof = exp(-(r/hwhm)**2/2.0_dbl)
    curve = 1.0_dbl - prof
    curve_err = 1e-5

    if( report ) then
       reporter%prof = exp(-(reporter%radius/hwhm)**2/2.0_dbl)
       reporter%curve = 1.0_dbl - reporter%prof
    end if

  end subroutine curve_init


  subroutine model_curve(nid,id,hwhm,nhwhm,raper,apcts,apcts_err,sky,grow,grow_err)

    use grow_model
    use grow_fit

    integer, intent(in out) :: nid
    integer, dimension(:), allocatable, intent(in out) :: id
    real, intent(in) :: hwhm
    integer, intent(in) :: nhwhm
    real, dimension(:), intent(in) :: raper, sky
    real, dimension(:,:), intent(in) :: apcts,apcts_err
    real(dbl), dimension(:), intent(out) :: grow,grow_err

    real(dbl), dimension(:,:), allocatable :: cts, dcts
    real, dimension(:), allocatable :: sky1
    real(dbl) :: a,aerr,b,berr,sig
    integer :: i,k,n,info

    n = size(raper)
    allocate(cts(nid,n),dcts(nid,n),sky1(nid))

    n = 0
    do i = 1,nid
       k = id(i)
       n = n + 1
       cts(n,:) = apcts(k,:)
       dcts(n,:) = apcts_err(k,:)
       sky1(n) = sky(k)
    end do

    call growmodel(raper,cts,dcts,hwhm,nhwhm,grow,grow_err)

!stop 0
    goto 55
    grow(1) = 0.378608555
    grow(2) = 0.556426227
    grow(3) = 0.760242879
    grow(4) = 0.897968650
    grow(5) = 0.972795010
    grow(6) = 0.988545179
    grow(7) = 0.992449760
    grow(8) = 0.994584739
    grow(9) = 0.996115625
    grow(10) = 0.996945143
    grow(11) = 0.997353375
    grow(12) = 0.997653544

    55 continue

    do i = 1,nid*0
       k = id(i)
       a = apcts(k,6)
       b = 0
       call growfit(raper,1d0*apcts(k,:),1d0*apcts_err(k,:),grow,grow_err,real(30), &
            a,aerr,b,berr,sig,info)
!       write(*,*) k,real(a),real(b),sky(k),'ccc'
!       stop
    end do

!stop
    deallocate(cts,dcts)

  end subroutine model_curve



!!$  subroutine grow_deepest(nid,id,hwhm,nhwhm,raper,apcts,apcts_err,sky_err,grow,grow_err,a,b)
!!$
!!$    ! The First Cut Is the Deepest
!!$
!!$    use grow_init
!!$
!!$    integer, intent(in out) :: nid
!!$    integer, dimension(:), allocatable, intent(in out) :: id
!!$    real, intent(in) :: hwhm
!!$    real(dbl), intent(in), optional :: a,b
!!$    integer, intent(in) :: nhwhm
!!$    real, dimension(:), intent(in) :: raper,sky_err
!!$    real, dimension(:,:), intent(in) :: apcts,apcts_err
!!$    real(dbl), dimension(:,:), intent(out) :: grow,grow_err
!!$
!!$    integer, dimension(:), allocatable :: mid
!!$    real(dbl), dimension(size(apcts,2)) :: cts, dcts
!!$    real(dbl) :: skyerr
!!$    integer :: i,k,j,n
!!$    logical :: valid
!!$
!!$    n = 0
!!$    allocate(mid(nid))
!!$
!!$    do i = 1,nid
!!$       k = id(i)
!!$       cts = apcts(k,:)
!!$       dcts = apcts_err(k,:)
!!$       skyerr = sky_err(k)
!!$       if( present(a) .and. present(b) ) then
!!$          call growinit(raper,cts,dcts,skyerr,hwhm,nhwhm,grow(k,:),grow_err(k,:),valid,a,b)
!!$       else
!!$          call growinit(raper,cts,dcts,skyerr,hwhm,nhwhm,grow(k,:),grow_err(k,:),valid)
!!$       end if
!!$       if( valid ) then
!!$          do j = 1,12
!!$!             write(*,'(i4,f5.1,2g20.5)') k,raper(j),grow(k,j)
!!$          end do
!!$
!!$          n = n + 1
!!$          mid(n) = k
!!$!          stop 0
!!$
!!$       end if
!!$    end do
!!$
!!$    deallocate(id)
!!$    allocate(id(n))
!!$    id = mid(1:n)
!!$    nid = n
!!$
!!$    deallocate(mid)
!!$
!!$!stop 0
!!$  end subroutine grow_deepest



  subroutine curve_spline(raper,curve,curve_err)

    real, dimension(:), intent(in) :: raper
    real(dbl), dimension(:), intent(in) :: curve,curve_err

    real(dbl), external :: smooth, ppvalu
    real(dbl), dimension(size(raper)+1) :: x,y,d
    real(dbl), dimension(size(raper)+1,7) :: aux
    real(dbl), dimension(size(raper)+1,4) :: aaux
    real(dbl), dimension(4,size(raper)+1) :: taaux

    integer :: i,n1,naper
    real(dbl) :: t,s,dnoise

    naper = size(raper)
    n1 = size(x) - 1

    ! approximate it
    dnoise = 1e-5       ! minimal error for spline fit (and display)
    x(1) = 0.0_dbl
    x(2:naper+1) = raper
    y(1) = 0.0_dbl
    y(2:naper+1) = curve
    d(1) = dnoise
    d(2:naper+1) = max(curve_err,dnoise)

    ! interpolation
    s = epsilon(x)
    ! aproximation
    s = naper
!    do i = 1,naper+1
!       write(*,*) i,real(x(i)),real(y(i)),real(d(i))
!    end do
    t = smooth(x,y,d,size(x),s,aux,aaux)
    taaux = transpose(aaux(:,1:4))

    if( report ) then
       do i = 1,size(reporter%radius)
          t = reporter%radius(i)
          reporter%curve(i) = ppvalu(x,taaux,n1,4,t,0)
          reporter%prof(i) = ppvalu(x,taaux,n1,4,t,1)
       end do
    end if

  end subroutine curve_spline

  subroutine grow_phot(raper,apcts,apcts_err,curve,curve_err,sky_sig, hwhm, &
       apcount,apcount_err, back, back_err)

    use grow_fit
    use oakleaf

    real, dimension(:), intent(in) :: raper,sky_sig
    real, dimension(:,:), intent(in) :: apcts,apcts_err
    real(dbl), dimension(:), intent(in) :: curve,curve_err
    real, intent(in) :: hwhm
    real(dbl), dimension(:), intent(out) :: apcount,apcount_err
    real, dimension(:), intent(out) :: back, back_err

    real(dbl), dimension(size(apcts,2)) :: cts, dcts, c,dc
    real, dimension(size(apcts,2)) :: r
    integer :: i,j,n,info,nhwhm,m
    real(dbl) :: sig,b,berr,a,aerr, rbest

    nhwhm = optimal_index(raper,hwhm)

    do i = 1, size(apcts,1)
       n = 0
       m = 0
       rbest = raper(size(raper)) - raper(1)
       do j = 1, size(apcts,2)
          if( apcts(i,j) > 0 ) then
             n = n + 1
             cts(n) = apcts(i,j)
             dcts(n) = apcts_err(i,j)
             c(n) = curve(j)
             dc(n) = curve_err(j)
             r(n) = raper(j)
             ! the valid aperture nearest to hwhm
             if( abs(r(n) - hwhm) < rbest ) then
                m = n
                rbest = abs(r(n) - hwhm)
             end if
          end if
       end do
       back(i) = 0
       back_err(i) = -1
       apcount(i) = -1
       apcount_err(i) = -1
       if( n > 0 .and. m > 0) then

          ! The values are initialized on values near HWHM (having optimal
          ! noise) and corrected for aperture correction. When fitting
          ! procedure will unsuccessfull, the values are considered as results.

          apcount(i) = cts(m)/c(m)
          apcount_err(i) = dcts(m)/c(m)

          if( n > nhwhm + 1 ) then
!             a = apcount(i)
             !             aerr = apcount_err(i)
             call fmean(cts(:n),dcts(:n),c(:n),dc(:n),a,aerr)
             ! only 1:nhwhm ?
             call growfit(r(1:n),cts(1:n),dcts(1:n),c(1:n),dc(1:n),sky_sig(i), &
                  a,aerr,b,berr,sig,info)
             if( info == 0 ) then
                apcount(i) = a
                apcount_err(i) = aerr
                back(i) = real(b)
                back_err(i) = real(berr)
             else
                back(i) = 0
                back_err(i) = 0
             end if
          end if

       end if
    end do

  end subroutine grow_phot


  function estim_rflux90(raper,curve)

    ! estimate a radius with 90% energy within

    real :: estim_rflux90
    real, parameter :: q90 = 0.9
    real, dimension(:), intent(in) :: raper
    real(dbl), dimension(:), intent(in) :: curve
    real(dbl) :: d
    integer :: i

    if( curve(1) > q90 ) then
       estim_rflux90 = real(q90*raper(1)/curve(1))
       return
    end if

    do i = 2,size(raper)
       if( curve(i-1) <= q90 .and. q90 < curve(i) ) then
          d = (raper(i) - raper(i-1)) / (curve(i) - curve(i-1))
          estim_rflux90 = real(raper(i-1) + d*(q90 - curve(i-1)))
          return
       end if
    end do

    estim_rflux90 = -1

  end function estim_rflux90


  function estim_hwhm(raper,prof)

    ! estimate a half of width in hight maximum

    real :: estim_hwhm
    real, parameter :: q = 0.5
    real, dimension(:), intent(in) :: raper
    real(dbl), dimension(:), intent(in) :: prof
    real, dimension(size(prof)) :: p, rm
    real :: d
    real(dbl) :: prof0
    integer :: i

    ! The profile values are mean surface brightnesses in the central
    ! disc and in the annuli between apertures; they are assigned to
    ! radii splitting the area of the disc and the annuli in halves
    ! (assigning them to the outer radius biases HWHM).
    rm(1) = raper(1) / sqrt(2.0)
    do i = 2,size(raper)
       rm(i) = sqrt((raper(i)**2 + raper(i-1)**2) / 2)
    end do

    ! extrapolate to raper == 0 with first and second point
    ! and the profile  prof(r) = a*r**2 + c
    prof0 = (prof(1)*rm(2)**2 - prof(2)*rm(1)**2) / (rm(2)**2 - rm(1)**2)

    p = real(prof / prof0)

    if( p(1) < q ) then
       estim_hwhm = real(q*rm(1)/p(1))
       return
    end if

    do i = 2,size(raper)
       if( p(i-1) >= q .and. q > p(i) ) then
          d = (rm(i) - rm(i-1)) / (p(i) - p(i-1))
          estim_hwhm = rm(i-1) + d*(q - p(i-1))
          return
       end if
    end do

    estim_hwhm = -1

  end function estim_hwhm


  function optimal_index(raper,hwhm) result(nhwhm)

    real, dimension(:), intent(in) :: raper
    real, intent(in) :: hwhm

    integer :: nhwhm,naper,i

    naper = size(raper)

    ! the optimal aperture for starting point, see section
    ! Optimal Aperture Choice in On Calibration ...
    nhwhm = naper / 2
    do i = 1, naper
       ! the factor 3 is for bright stars with gaussian profile
       if( raper(i) > 3*hwhm ) then
          nhwhm = i
          exit
       end if
    end do

  end function optimal_index


end module grow_curve
