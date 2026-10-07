!
!  Robust construction of a growth-curve.
!
!  Copyright © 2016-7 F.Hroch (hroch@physics.muni.cz)
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
!

module grow_model

  implicit none

  ! numerical precision of real numbers
  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: pi = 3.14159265358979312_dbl

  ! assumed error of approximation by asymptote
  real(dbl), parameter, private :: approxerr = 1e-5

  ! print debug informations ?
  logical, parameter, private :: debug = .false.

!  real(dbl), dimension(:), allocatable, private :: aper,raper,cts_, dcts_, back, &
!       grow_,dgrow_, tx, fvec, dvec, radius, phix
  real(dbl), dimension(:), allocatable, private :: aper,raper,radius,grow,dgrow, sigs
  real(dbl), dimension(:,:), allocatable, private :: cts,dcts,grows,dgrows
!  real(dbl), private :: hwhm, par_a, par_b, par_e, par_c, par_t, sig, lambda, rback
!  real(dbl), private :: hwhm, sig, lambda, rback
  real(dbl), private :: hwhm, sig, lambda
  integer, private :: nhwhm, naper, nstar, nref
!  integer, private :: type

  !private :: minall, the_constructor

contains


  subroutine growmodel(xraper,xcts,xdcts,rhwhm,mhwhm,curve,dcurve)

    use oakleaf

    real, dimension(:), intent(in) :: xraper
    real(dbl), dimension(:,:), intent(in) :: xcts,xdcts
    real, intent(in) :: rhwhm
!    real(dbl), intent(in) :: skyerr
    integer, intent(in) :: mhwhm
    real(dbl), dimension(:), intent(out) :: curve,dcurve

!    real(dbl), dimension(3) :: p
!    real(dbl), dimension(size(xraper)) :: ff,aa
!    real(dbl), dimension(:), allocatable :: back1,tx1,flux0,aas,daas,betas,qs,phi,&
!         bs
    real(dbl), dimension(:), allocatable :: flux,betas,phi
!    real(dbl) :: a,b,c,e,a1,b1,e1,t,d,s,t1,dt,db,tol,q,c1,ctsback,backmed,beta
    real(dbl) :: a,b
!    real(dbl), dimension(2,2) :: mat
!    real(dbl), dimension(2) :: vec
!    integer :: i, iter, j, n, n1, n2
    integer :: i, iter, n
!    logical :: valid

    nstar = size(xcts,1)
    naper = size(xraper)
!    allocate(raper(naper),aper(naper),back(nstar),back1(nstar),cts_(naper),dcts_(naper),&
!         cts(nstar,naper),dcts(nstar,naper),grows(nstar,naper),dgrows(nstar,naper),&
!         tx(nstar),tx1(nstar),grow_(naper),dgrow_(naper),radius(naper), &
!         flux0(nstar),aas(nstar),daas(nstar),betas(nstar),qs(nstar),phi(nstar),&
!         bs(nstar),phix(nstar))
    allocate(raper(naper),aper(naper),radius(naper),flux(nstar),betas(nstar),phi(nstar),&
         cts(nstar,naper),dcts(nstar,naper),grows(nstar,naper),dgrows(nstar,naper),&
         grow(naper),dgrow(naper),sigs(nstar))

    ! apertures radii in pixels
    raper = xraper

    ! aperture areas
    aper = pi*raper**2

    ! radius is scaled to interval 0..1, to improve numerical stability
    radius = raper / raper(naper)

    cts = xcts
    dcts = xdcts
    hwhm = rhwhm
    nhwhm = mhwhm
!    tol = sqrt(epsilon(tol))
!    tol = 10*epsilon(tol)
!    tol = 1e-8
!    nref = (nhwhm + naper)/2 ! min 5 points required?
    nref = nhwhm

    ! osetript pripad, ze nstar < 3 ?



!    if( debug ) then ! write initial grows

!!$    flux0 = cts(:,nhwhm)
!!$    do i = 1,nstar
!!$       call arc_init(cts(i,:),dcts(i,:),cts(i,nhwhm),aas(i),bs(i),betas(i))
!!$    end do
!!$
!!$! nutne:   call rwmean(aas,daas,a,s)
!!$    call rmean(aas,a,s)
!!$    write(*,*) 'aas:',real(a),real(s)
!!$    call rmean(bs,b,s)
!!$    write(*,*) 'bs:',real(b),real(s)

    a = 0
    b = 0
    betas = 0
    phi = 0
    flux = cts(:,nref)

!    call arc_inits1(flux0,b,e,betas)
!    write(*,*) b
    call arc_init(flux,b)
!    call arc_init1(flux,a,b)
    if( debug ) write(*,'(a,3f15.5)') 'Very first a,b:',a,b

    !    write(*,*) flux0(1),flux0(1)*(1 + e)
!    flux0 = flux0*(1 - e)


!    do n = 1,nstar
!       do i = nref, naper
!          write(*,*) raper(i),cts(n,i)/flux(n),-a/raper(i)-b/raper(i)**2
!       end do
!    end do
!    stop
!    call rmean(tx(1:5),t,s) ! nutne rwmean
!    write(*,*) 't:',real(t),real(s)
!    stop

    flux = flux * (1 - b/raper(naper)**2) / (1 - b/raper(nref)**2)

!    phi = 0
!    a = 0.0045
!    b = 0.13
!    call anchor_inits1(flux0,a,b,e)
!    write(*,*) 1+t

!stop
!    write(*,*) flux0(1),flux0(1)*(1 - t)
!    flux0 = flux0*(1 - t)

    ! checkpoint: first estimate of arc parameters and anchor

!    t = 1 + t
    ! t anchors grow curve at nref point

    ! update flux0
    !    flux0 = flux0 * (1 + phi*0) / t

    do iter = 1,50

       ! determine beta,phi for every star
       do n = 1,nstar
          call back_update(cts(n,:),dcts(n,:),flux(n),a,b,betas(n),phi(n))
          if( debug ) &
               write(*,'(i3,2es12.2,2f15.1)') n,betas(n),phi(n),flux(n)*(1+phi(n)),flux(n)
!                 if( n == 2 ) stop
       end do

       call arc_update(flux,phi,betas,a,b)
       !call arc_update1(flux,phi,betas,b)
       if( debug ) write(*,'(a,i5,3f15.5)') 'iter,a,b:',iter,a,b
       !stop

       flux = flux*(1 + phi)

       if( all(abs(phi) < 0.1*sqrt(flux)/flux ) ) exit

    end do

    do n = 1,nstar*0
       do i = 1,naper
          write(*,*) raper(i),cts(n,i)/flux(n),cts(n,i)/flux(n)-betas(n)*aper(i), &
               i, n, 'FF'
       end do
    end do

    ! grow curve
    do i = naper,nref,-1
       curve(i) = 1 - (a/raper(i) + b/raper(i)**2)
       do n = 1,nstar
          grows(n,i) = cts(n,i)/flux(n) - betas(n)*aper(i)
          dgrows(n,i) = grows(n,i) * (dcts(n,i) / cts(n,i))
       end do
!       call rwmean(grows(:,i),dgrows(:,i),grow(i),dgrow(i))
       call rmean(grows(:,i),grow(i),dgrow(i))
       dcurve(i) = dgrow(i)
       if( debug ) write(*,'(i3,3e15.5)') i,grow(i),dgrow(i),grow(i)-curve(i)
    end do
!    t = grow(nref)
!    grows(:,nref+1) = grow(nref+1)
    betas = betas*flux
    do i = nref,1,-1
       do n = 1,nstar
          grows(n,i) = grows(n,i+1) * &
               (cts(n,i) - betas(n)*aper(i)) / (cts(n,i+1) - betas(n)*aper(i+1))
!          grows(n,i) = (cts(n,i) - betas(n)*aper(i)) /flux(n)
          dgrows(n,i) = grows(n,i) * (dcts(n,i) / cts(n,i))
!          write(*,*) i,n,real(grows(n,i))
       end do
!       call rwmean(grows(:,i),dgrows(:,i),curve(i),dcurve(i))
       call rmean(grows(:,i),curve(i),dcurve(i))
       dcurve(i) = dcurve(i) + approxerr
    end do


!!$goto 999
!!$    stop
!!$
!!$
!!$
!!$
!!$
!!$
!!$    ! Initial estimate of individual backgrounds of stars.
!!$    ! Its estimate is on base of slope of every individual
!!$    ! growth curve close to largest aperture.
!!$    !
!!$    ! rback is radius where the slope is determined
!!$    ! ctsback is mean of counts at rback (change by slope is negligible)
!!$    rback = sum(radius(nref:naper))/(naper - nref + 1)
!!$!    write(*,*) rback,real(radius(nref:naper)),naper - nref
!!$    rback = 1
!!$
!!$    do i = 1,nstar*0
!!$!       call back_init(cts(i,:),dcts(i,:),back(i))
!!$       ctsback = sum(cts(i,nref:naper)) / (naper - nref + 1)
!!$!       write(*,*) real(cts(i,nref:naper))
!!$       back(i) = ctsback * back(i) / (2*pi*rback*raper(naper)**2)
!!$       write(*,*) real(sky(i)-1e3),real(back(i)), real(back(i)-(sky(i)-1e3)),real(ctsback)
!!$    end do
!!$    back = sky-1e3
!!$!stop
!!$!    if( debug ) then ! write grows with corrected sky
!!$
!!$
!!$    ! Initial estimate of asympotic parameters
!!$    ! There is a lot of curves with potentially false points
!!$    ! and the parameters are determined from a few points.
!!$    ! Thats is why the data are filtered by median.
!!$    !
!!$    ! The curve is anchored on number one at naper, so the
!!$    ! curve needs additional normalisation.
!!$    !
!!$    t = 1
!!$!    back = 0
!!$    do i = 1,nstar
!!$       ! starting construction from hwhm (last element > 1)
!!$       call a_constructor(cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$       if( .not. valid ) stop 'a_constructor: not valid'
!!$    end do
!!$    do i = 1,naper
!!$       call rinit(grows(:,i),grow_(i),dgrow_(i))
!!$!       write(*,'(3f15.5)') raper(i),grow_(i),dgrow_(i)
!!$    end do
!!$    dgrow_(nhwhm) = (dgrow_(nhwhm-1) + dgrow_(nhwhm+1))/2 ! reference has zero noise
!!$    grow_ = grow_ / grow_(naper)    ! normalize
!!$    t = grow_(nref) - grow_(nhwhm)
!!$    do i = 1,naper
!!$       write(*,'(3f15.5)') raper(i),grow_(i),dgrow_(i)
!!$    end do
!!$!    call asymptote_init(a,b,c,q)
!!$!    write(*,'(a,4f12.5)') 'asymptote_init (q=1?):',a,b,c,q
!!$
!!$    do i = 1,nstar
!!$       call rcal(cts(i,1:nhwhm),dcts(i,1:nhwhm),grow_(1:nhwhm),dgrow_(1:nhwhm), &
!!$            flux0(i),s,.false.)
!!$       s = sum(abs(cts(i,1:nhwhm)/grow_(1:nhwhm)-flux0(i)))/(nhwhm-1)
!!$       write(*,*) i,real(flux0(i)),real(s),real(cts(i,nhwhm))
!!$    end do
!!$
!!$    a = 0
!!$    ! inidividual grows
!!$    do i = 1,nstar
!!$       grow_ = (cts(i,:) - back(i)*pi*raper(:)**2)/ flux0(i)
!!$       dgrow_ = dcts(i,:) / flux0(i)
!!$       call asymptote_init5(a,q,beta,daas(i),b,c)
!!$       write(*,*) i,a,beta*flux0(i),q,flux0(i),daas(i),1e3-sky(i),'init5'
!!$!       flux0(i) = q*flux0(i)
!!$       qs(i) = q
!!$       aas(i) = a/q
!!$       daas(i) = daas(i)/q
!!$       betas(i) = beta/q
!!$!       stop
!!$    end do
!!$
!!$
!!$    call rwmean(aas,daas,a,s)
!!$    write(*,*) 'aas:',real(a),real(s)
!!$
!!$    call rmean(qs,q,s)
!!$    write(*,*) 'qs:',real(q),real(s)
!!$
!!$    ! update individual fluxes and backgrounds
!!$    a = 0
!!$    do i = 1,nstar
!!$       call update_backflux(a,cts(i,:),dcts(i,:),flux0(i),back(i))
!!$       write(*,*) i,real(flux0(i)),real(back(i)),real(1e3-sky(i))
!!$    end do
!!$    betas = back
!!$
!!$stop
!!$
!!$    ! update grow curve
!!$    do i = nref,naper
!!$       grow_(i) = 1 - a/raper(i)! - t
!!$       write(*,*) real(raper(i)),real(grow_(i))
!!$    end do
!!$    grows(:,nref) = grow_(nref)
!!$    do j = nref-1,1,-1
!!$       do i = 1,nstar
!!$          grows(i,j) = & !grows(i,j+1)* &
!!$               (cts(i,j) - betas(i)*raper(j)**2)/(cts(i,j+1) - betas(i)*raper(j+1)**2)
!!$       end do
!!$       call rwmean(grows(:,j),dgrows(:,j),t,dgrow_(j))
!!$       grow_(j) = t*grow_(j+1)
!!$       write(*,*) real(raper(j)),real(grow_(j)),real(dgrow_(j))
!!$    end do
!!$
!!$    grow = grow_
!!$    dgrow = dgrow_
!!$!    stop
!!$    return
!!$
!!$
!!$
!!$
!!$
!!$    do i = 1,nstar
!!$       grow_ = cts(i,:) / flux0(i)
!!$       dgrow_ = dcts(i,:) / flux0(i)
!!$       call asymptote_init5(a,q,beta,daas(i),b,c)
!!$!       write(*,*) i,a,beta*flux0(i),q,flux0(i),daas(i),1e3-sky(i),'init5'
!!$       qs(i) = q
!!$       aas(i) = a/q
!!$       daas(i) = daas(i)/q
!!$       betas(i) = beta/q
!!$    end do
!!$    call rwmean(aas,daas,a,s)
!!$    write(*,*) 'aas1:',real(a),real(s)
!!$
!!$
!!$    do i = 1,nstar
!!$       call update_backflux(a,cts(i,:),dcts(i,:),flux0(i),back(i))
!!$       write(*,*) i,real(flux0(i)),real(back(i))
!!$    end do
!!$
!!$
!!$
!!$stop
!!$
!!$
!!$
!!$    ! update fluxes
!!$    do i = 1,nstar
!!$       t = flux0(i)
!!$       call rcal(cts(i,1:nhwhm),dcts(i,1:nhwhm),grow_(1:nhwhm),dgrow_(1:nhwhm), &
!!$            flux0(i),s,.false.)
!!$       write(*,*) i,real(flux0(i)),real(t),real(s),'fff'
!!$    end do
!!$stop
!!$    ! update inidividual grows
!!$    do i = 1,nstar
!!$       grow_ = cts(i,:) / flux0(i)
!!$       dgrow_ = dcts(i,:) / flux0(i)
!!$       call asymptote_init5(a,q,beta,daas(i),b,c)
!!$       write(*,*) i,a,beta*flux0(i),q,flux0(i),daas(i),1e3-sky(i),'init7'
!!$!       flux0(i) = q*flux0(i)
!!$       qs(i) = q
!!$       aas(i) = a/q
!!$       daas(i) = daas(i)/q
!!$       betas(i) = beta/q
!!$    end do
!!$
!!$    call rwmean(aas,daas,a,s)
!!$    write(*,*) 'aas:',real(a),real(s)
!!$
!!$    ! update grow curve
!!$    do i = nref,naper
!!$       grow_(i) = 1 - a/raper(i)! - t
!!$       write(*,*) real(raper(i)),real(grow_(i))
!!$    end do
!!$    grows(:,nref) = grow_(nref)
!!$    do j = nref-1,1,-1
!!$       do i = 1,nstar
!!$          grows(i,j) = grows(i,j+1)* &
!!$               (cts(i,j) - betas(i)*raper(j)**2)/(cts(i,j+1) - betas(i)*raper(j+1)**2)
!!$       end do
!!$       call rwmean(grows(:,j),dgrows(:,j),grow_(j),dgrow_(j))
!!$       write(*,*) real(raper(j)),real(grow_(j)),real(dgrow_(j))
!!$    end do
!!$stop
!!$    grow = grow_
!!$    dgrow = dgrow_
!!$
!!$    return
!!$
!!$!    stop
!!$
!!$    do i = 1,nstar
!!$       call rmean((cts(i,nref:)-flux0(i)*grow_(nref:))/aper(nref:),back(i),s)
!!$!       call rmean((cts(i,nref:)/tx(i)-grow_(nref:))/aper(nref:),back(i),s)
!!$!       write(*,*) real((cts(i,nref:)-0*tx(i)*grow_(nref:))/aper(nref:))
!!$!       call rmean(cts(i,:)/grow_,tx(i),s)
!!$       write(*,*) i,real(back(i)/3.14),real(s),real(sky(i)-1e3)
!!$       if( i == nstar ) stop
!!$    end do
!!$    ! v teto fazi je jeste rozptyl prilis velky a pozadi urcovat nejde
!!$    ! a urceni f0 je taky zbytecne...
!!$    back = 0
!!$
!!$
!!$!    do i = 1,nstar
!!$!       grows(i,:) = cts(i,:) / tx(i)
!!$    !    end do
!!$    t = 1
!!$    back(1) = 5
!!$    do i = 1,nstar
!!$       call the_constructor(t,back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$    end do
!!$    do i = 1,naper
!!$!       call rinit(grows(:,i),grow_(i),dgrow_(i))
!!$       call rwmean(grows(:,i),dgrows(:,i),grow_(i),dgrow_(i))
!!$!       write(*,'(3f15.5)') raper(i),grow_(i),dgrow_(i)
!!$       write(*,'(3f15.5)') raper(i),cts(1,i),dgrows(1,i)
!!$    end do
!!$!        stop
!!$
!!$    grow_ = cts(1,:) / 1.603e7
!!$    dgrow_ = dcts(i,:) / 1.603e7
!!$    grow_ = cts(1,:) / 1.7e7
!!$    dgrow_ = dcts(i,:) / 1.7e7
!!$    p(2) = 1
!!$    do i = -1000,1000,100
!!$       p(1) = i/5e3
!!$       do j = -500,500,50
!!$          p(3) = j/1e3
!!$          write(*,*) p(1),p(3),fainit5(p),'qqq'
!!$       end do
!!$    end do
!!$    stop
!!$
!!$!    call asymptote_init3(a,q)
!!$
!!$
!!$    do i = 1,nstar
!!$       grow_ = grows(i,:)
!!$       dgrow_ = dgrows(i,:)
!!$       call asymptote_init3(a,q)
!!$       write(*,*) i,a,q,'aaa'
!!$       stop
!!$    end do
!!$
!!$
!!$    c = 0
!!$    stop
!!$
!!$    ! The curve's anchor point is given by q as 1-q at maximum aperture.
!!$    ! Both a,b parameters are only sligthly different against to final (true) ones.
!!$
!!$!stop
!!$
!!$    do iter = 1,13
!!$
!!$       ! Next iteration
!!$       call asymptote_init2(a,b,q)
!!$       write(*,'(a,4f12.5)') 'asymptote_init2 (q->1):',a,b,q
!!$       stop
!!$
!!$
!!$       !       t = 1 - q
!!$!       t = 1 + (1-q) !- (a + b/raper(naper))
!!$       t = 1 - (a + b/raper(naper))/raper(naper)
!!$       write(*,*) t
!!$!       stop
!!$       do i = nref,naper
!!$!          grow_(i) = asymptote(raper(i),a,b,c)
!!$          grow_(i) = t - (a + b/raper(i))/raper(i)
!!$          ! dgrow_?
!!$          write(*,'(3f15.5)') raper(i),grow_(i),dgrow_(i)
!!$       end do
!!$       grows(:,nref) = grow_(nref)
!!$       do j = nref-1,1,-1
!!$          do i = 1,nstar
!!$             grows(i,j) = grows(i,j+1)*(cts(i,j) - back(i)*aper(j))/(cts(i,j+1) - back(i)*aper(j+1))
!!$          end do
!!$       end do
!!$
!!$       do i = 1,nref-1
!!$          call rwmean(grows(:,i),dgrows(:,i),grow_(i),dgrow_(i))
!!$          write(*,'(3f15.5)') raper(i),grow_(i),dgrow_(i)
!!$       end do
!!$
!!$       do i = 1,nstar
!!$!          call rcal(cts(i,1:nref),dcts(i,1:nref),grow_(1:nref),dgrow_(1:nref),tx(i),s,.false.)
!!$          call rcal(cts(i,:),dcts(i,:),grow_,dgrow_,tx(i),s,.false.)
!!$          write(*,*) i,real(tx(i)),real(s),real(sky(i)-1e3)
!!$       end do
!!$
!!$       do i = 1,nstar
!!$          call rmean((cts(i,nref:)-tx(i)*grow_(nref:))/aper(nref:),back(i),s)
!!$!          write(*,*) real((cts(i,nref:)-tx(i)*grow_(nref:))/aper(nref:))
!!$          write(*,*) i,real(back(i)),real(s),real(sky(i)-1e3)
!!$!          if( i == 2 ) stop
!!$       end do
!!$!       stop
!!$
!!$!       t = 1 - q
!!$       do i = 1,nstar
!!$          call the_constructor(t,back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$       end do
!!$       do i = 1,naper
!!$          call rwmean(grows(:,i),dgrows(:,i),grow_(i),dgrow_(i))
!!$          write(*,'(3f15.5)') raper(i),grow_(i),dgrow_(i)
!!$       end do
!!$
!!$    !       call asymptote_update(a,b,c)
!!$!       write(*,'(a,4f12.5)') 'asymptote_update:',a,b,c
!!$
!!$    end do
!!$
!!$stop
!!$
!!$    ! Update of the initials of asymptote with anchor value as mean at naper.
!!$    ! update grow curve by the parameters (correct solution is
!!$    ! update of t as t = 1/q and repeat the procedure while q will
!!$    ! approach 1, but we are belives, that the all deviations are small...)
!!$!    t = 1/q
!!$!    t = 0.997
!!$    t = asymptote(raper(naper),a,b,c)
!!$!    write(*,*) t
!!$    do i = 1,nstar
!!$       call the_constructor(t,back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$    end do
!!$    do i = 1,naper
!!$!       call rinit(grows(:,i),grow_(i),dgrow_(i))
!!$       call rwmean(grows(:,i),dgrows(:,i),grow_(i),dgrow_(i))
!!$    end do
!!$    dgrow_(naper) = dgrow_(naper-1)
!!$!    call asymptote_init(a,b,c)
!!$    call asymptote_init(a,b,c,q)
!!$    write(*,'(a,4f12.5)') 'asymptote_init:',a,b,c,q
!!$
!!$    stop
!!$    ! choose determination of method to minimise
!!$    ! setup norm (absolute deviations)
!!$    type = 0
!!$
!!$    ! Now, its time to do more precise estimation of background..
!!$    call rinit(abs(back)/cts(:,naper),backmed,s)
!!$    do i = 1,nstar
!!$       call anchor(cts(i,:),dcts(i,:),a,b,c,t,backmed,back(i))
!!$       write(*,*) real(sky(i)-1e3),real(back(i)), real(back(i)-(sky(i)-1e3)), sqrt(cts(i,naper))
!!$    end do
!!$
!!$    ! .. anchor t
!!$!    call anchor_update(a,b,c,back,t)
!!$!    write(*,'(a,4f12.5)') 'anchor_update:',t
!!$
!!$    ! .. and asymptote
!!$    call asymptote_update(a,b,c,t)
!!$    t = asymptote(raper(naper),a,b,c)
!!$    write(*,'(a,4f12.5)') 'asymptote_update:',a,b,c,t
!!$
!!$
!!$    ! Recompute with least square meaning statistics
!!$    ! setup norm (Huber's function)
!!$    type = 1
!!$
!!$    ! determine mean sig in tail
!!$    call asymptote_sig(a,b,c,t,back,sig)
!!$
!!$    ! Now, its time to do more precise estimation of background..
!!$    call rinit(abs(back)/cts(:,naper),backmed,s)
!!$    do i = 1,nstar
!!$       call anchor(cts(i,:),dcts(i,:),a,b,c,t,backmed,back(i))
!!$       write(*,*) real(sky(i)-1e3),real(back(i)), real(back(i)-(sky(i)-1e3)), i
!!$    end do
!!$
!!$    ! .. anchor t
!!$!    call anchor_update(a,b,c,back,t)
!!$!    write(*,'(a,4f12.5)') 'anchor_update:',t
!!$
!!$    ! .. and asymptote
!!$    call asymptote_update(a,b,c,t)
!!$    t = asymptote(raper(naper),a,b,c)
!!$    write(*,'(a,4f12.5)') 'asymptote_update:',a,b,c,t
!!$
!!$
!!$
!!$
!!$    stop 0
!!$
!!$!    write(*,*) 'nref=',nref
!!$
!!$
!!$!    lambda = 1
!!$!    do j = 1,0
!!$!       call par_init(a,b,c,back)
!!$!       lambda = 2*lambda
!!$!    end do
!!$!    back =
!!$
!!$
!!$    open(1,file='/tmp/1')
!!$!    back = 0
!!$    tx = 1
!!$    t = 1 !0.9984
!!$    !    stop
!!$    t = 1
!!$    do i = 1,nstar
!!$       call the_constructor(t,back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$    end do
!!$    do i = 1,naper
!!$!       write(*,*) 'i=',i
!!$!       call rmean(grows(:,i),grow_(i),dgrow_(i))
!!$       call rinit(grows(:,i),grow_(i),dgrow_(i))
!!$!       call rmean(grows(:,i),a,b)
!!$!       call rwmean(grows(:,i),dgrows(:,i),grow_(i),dgrow_(i))
!!$       write(*,*) 'i=',i,real(raper(i)),real(grow_(i)),real(dgrow_(i))
!!$       write(1,*) raper(i),grow_(i),dgrow_(i)
!!$!       write(*,'(4g15.5,a)') a,b,grow_(i),dgrow_(i),' www'
!!$!       if( i == 11 ) then
!!$!          do j = 1,nstar
!!$!             write(*,*) grows(j,i),dgrows(j,i)**2
!!$!          end do
!!$       !       end if
!!$!       if( i == 9  ) stop 0
!!$    end do
!!$    close(1)
!!$
!!$    dgrow_(naper) = dgrow_(naper-1)
!!$    write(*,*) dgrow_(naper),dgrow_(naper-1),dgrow_(naper-2),naper
!!$
!!$
!!$    ! initial estimate of asympotic parameters
!!$    lambda = 0
!!$    do i = 1,1
!!$       call asymptote_init(a,b,c,q)
!!$       lambda = 2*lambda
!!$    end do
!!$    write(*,'(a,4f12.5)') 'asymptote_init:',a,b,c,q
!!$
!!$
!!$    ! update grow curve by the parameters (correct solution is
!!$    ! update of t as t = 1/q and repeat the procedure while q will
!!$    ! approach 1, but we are belives, that the all deviations are small...)
!!$    if( q > 1 ) grow_ = grow_ / q
!!$
!!$    ! compute asymptote
!!$    lambda = 0
!!$    do i = 1,1
!!$       call asymptote_init(a,b,c)
!!$       lambda = 2*lambda
!!$    end do
!!$
!!$    call asymptote_init(a,b,c,q)
!!$    write(*,'(a,4f12.5)') 'asymptote_init1:',a,b,c,q
!!$
!!$    stop
!!$    ! setup norm (absolute deviations)
!!$    type = 0
!!$
!!$
!!$    do j = 1,1
!!$
!!$    ! update individual t,b-parameters
!!$    tx = asymptote(raper(naper),a,b,c)
!!$    write(*,*) 'tx(1)',tx(1)
!!$    do i = nref, naper
!!$       grow_(i) = asymptote(raper(i),a,b,c)
!!$    end do
!!$    do i = 1,nstar
!!$       cts_ = cts(i,:)
!!$       dcts_ = dcts(i,:)
!!$!       call anchor(a,b,c,tx(i),back(i))
!!$!       write(*,*) i,tx(i),tx(i)-asymptote(raper(naper),a,b,c)
!!$       write(*,*) i,back(i)
!!$    end do
!!$
!!$
!!$    ! update grow curve
!!$    do i = 1,nstar
!!$       call the_constructor(tx(i),back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$    end do
!!$    do i = 1,naper
!!$       call rinit(grows(:,i),grow_(i),dgrow_(i))
!!$!       write(*,*) 'i=',i,real(raper(i)),real(grow_(i)),real(dgrow_(i))
!!$    end do
!!$    dgrow_(naper) = dgrow_(naper-1)
!!$    write(*,*) dgrow_(naper),dgrow_(naper-1),dgrow_(naper-2),naper
!!$    lambda = 0
!!$    call asymptote_init(a,b,c,q)
!!$    write(*,'(a,4f12.5)') 'asymptote_init:',a,b,c,q
!!$
!!$ enddo
!!$
!!$stop
!!$
!!$    grow_(9) = 0.99611562
!!$    grow_(10) = 0.996945143
!!$    grow_(11) = 0.997353375
!!$    grow_(12) = 0.997653544
!!$
!!$    open(1,file='/tmp/3')
!!$    back = 0
!!$    do i = 1,nstar
!!$       call the_constructor(tx(i),back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$       do j = nref,naper
!!$          write(1,*) raper(j),grows(i,j),(grows(i,j)-grow_(j))/dgrows(i,j)
!!$       end do
!!$    end do
!!$    close(1)
!!$
!!$    ! update asymptote
!!$    do i = 1,nstar
!!$       call the_constructor(tx(i),back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$    end do
!!$
!!$    lambda = 1
!!$    a1 = 100
!!$    b1 = 100
!!$    c1 = 100
!!$    do i = 1,10
!!$!       call asymptote_update(a,b,c)
!!$       write(*,'(i3,4f15.5)') i,a,b,c
!!$          if( abs(c - c1) < tol .and. abs(a - a1) < tol .and. abs(b - b1) < tol ) exit
!!$       lambda = 2*lambda
!!$       a1 = a
!!$       b1 = b
!!$       c1 = c
!!$    end do
!!$
!!$
!!$    ! setup norm (Huber's function)
!!$    type = 1
!!$
!!$    ! determine mean sig in tail
!!$!    call asymptote_sig(a,b,c,tx,sig)
!!$
!!$    ! update t-par
!!$    do i = nref, naper
!!$       grow_(i) = asymptote(raper(i),a,b,c)
!!$    end do
!!$    do i = 1,nstar
!!$       cts_ = cts(i,:)
!!$       dcts_ = dcts(i,:)
!!$!       call anchor(a,b,c,tx(i),back(i))
!!$       write(*,*) i,tx(i),tx(i)-asymptote(raper(naper),a,b,c)
!!$    end do
!!$
!!$    ! update asymptote, with Huber's help
!!$    do i = 1,nstar
!!$       call the_constructor(tx(i),back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$    end do
!!$
!!$    lambda = 1
!!$!    a1 = 100
!!$!    b1 = 100
!!$!    c1 = 100
!!$    do i = 1,10
!!$!       call asymptote_update(a,b,c)
!!$       write(*,'(i3,4f15.5)') i,a,b,c
!!$          if( abs(c - c1) < tol .and. abs(a - a1) < tol .and. abs(b - b1) < tol ) exit
!!$       lambda = 2*lambda
!!$       a1 = a
!!$       b1 = b
!!$       c1 = c
!!$    end do
!!$
!!$    goto 99
!!$    stop 0
!!$
!!$!    call rmean(dgrows(nhwhm,:),t,dt)
!!$    !    dgrow(nhwhm) = t
!!$
!!$!    call background3(back)
!!$!    write(*,*) real(back)
!!$
!!$
!!$!    stop 0
!!$    t = grow_(nref) / grow_(naper)
!!$
!!$    ! Initial value must be 0 < t < 1.
!!$    if( 0 < t .and. t < 1 ) then
!!$       continue
!!$    else
!!$       t = 0.95
!!$    end if
!!$
!!$!    t = t - 3*dgrow_(naper)
!!$    t = t - 0.002
!!$!    t = 0.9*t
!!$
!!$!    t = 0.5
!!$!    e = 1 - t
!!$!    e = e / 2
!!$!    a = ((1 - e) - t) * raper(nref)
!!$!    a = (1 - t) * raper(nref)
!!$    a = (1 - t) / 2
!!$!    c = a * raper(nref)**3
!!$    b = a * raper(nref)**2
!!$    a = a * raper(nref)
!!$    e = 0
!!$    write(*,*) t,asymptote(raper(naper),a,b,c)
!!$!    t = t * asymptote(raper(naper),e,a,b)
!!$
!!$    write(*,*) grow_(naper) - grow_(naper-1), grow_(naper-2) - grow_(naper-1)
!!$    write(*,*) 1/raper(naper) - 1/raper(naper-1)
!!$    write(*,*) 1/raper(naper-2) - 1/raper(naper-1)
!!$    write(*,*) 1/raper(naper)**2 - 1/raper(naper-1)**2
!!$    write(*,*) 1/raper(naper-2)**2 - 1/raper(naper-1)**2
!!$
!!$    mat(1,1) = 1/raper(naper) - 1/raper(naper-1)
!!$    mat(1,2) = 1/raper(naper)**2 - 1/raper(naper-1)**2
!!$    mat(2,1) = 1/raper(naper-2) - 1/raper(naper-1)
!!$    mat(2,2) = 1/raper(naper-2)**2 - 1/raper(naper-1)**2
!!$    mat = -mat
!!$    vec(1) = grow_(naper) - grow_(naper-1)
!!$    vec(2) = grow_(naper-2) - grow_(naper-1)
!!$
!!$    d = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
!!$
!!$    a = (vec(1)*mat(2,2) - vec(2)*mat(1,2)) / d
!!$!    a = -a / d
!!$    b = (vec(1)*mat(1,1) - vec(2)*mat(2,1)) / d
!!$!    b = b / d
!!$    e = 0
!!$    !    t =  grow_(nref) / grow_(naper) * asymptote(raper(naper),e,a,b)
!!$!    t = asymptote(raper(nref),e,a,b)
!!$    a = vec(1) / mat(1,1)
!!$    b = 0
!!$    call part(e,a,b,t)
!!$
!!$!    t = asymptote(raper(nref),a,b)
!!$!    write(*,*) a,b,t
!!$!    c = b /  !1 - t
!!$!    e = 1 - t - a/raper(nref) - b/raper(nref)**2
!!$!    a = a / 2
!!$!    e = 1 - a
!!$!    write(*,*) real(t),real(a),real(e)
!!$
!!$!    write(*,*) 1 - e - a / raper(nref) - b / raper(nref)**2
!!$    write(*,*) 'Initial estimate:',real(t),real(a),real(b),real(e)
!!$
!!$!    call tnoise(t)
!!$!    write(*,*) 'tnoise=',t
!!$
!!$!    open(11,file='/tmp/t')
!!$!    do i = 5000,10000
!!$!       t = i/1e4
!!$!       write(11,*) t,noisefun(t)
!!$!    end do
!!$!    close(11)
!!$
!!$!    stop 0
!!$
!!$!    stop 0
!!$!    back1 = 0! maxval(cts)
!!$    back = 0
!!$    back1 = back
!!$    tx1 = t
!!$    tx = t
!!$    t1 = t
!!$    dt = 1 - t
!!$
!!$    par_t = t
!!$
!!$    goto 33
!!$    open(11,file='/tmp/s')
!!$    do i = 1,1000,10
!!$       do j = 1,1000,10
!!$          b = i*1e-3
!!$          a = j*1e-3
!!$          e = j*1e-4
!!$!          s = mincom((/e,0.0_dbl,b/))
!!$!          s = mincom((/e,a,0.0_dbl/))
!!$          s = mincom((/0.0_dbl,a,b/))
!!$          write(11,*) a,b,s
!!$       end do
!!$    end do
!!$    close(11)
!!$33 continue
!!$
!!$    open(11,file='/tmp/t')
!!$    par_a = a
!!$    par_b = b
!!$    do i = 9600,10000
!!$       t = i/1e4
!!$!       write(11,*) t,tfun(t)
!!$    end do
!!$    close(11)
!!$
!!$
!!$    !    call towel(a,b,t)
!!$    call part(e,a,b,t)
!!$    write(*,*) t
!!$
!!$!    stop 0
!!$
!!$    lambda = 1
!!$    do iter = 1, precision(t)
!!$
!!$
!!$       ! estimate parameters common to all stars: a,b,e
!!$
!!$       back = 0
!!$       type = 0
!!$       lambda = 1
!!$       e1 = 100
!!$       a1 = 100
!!$       b1 = 100
!!$       do i = 1,20
!!$          call commons1(e,a,b)
!!$          !          call part(e,a,b,t)
!!$!          t = asymptote(raper(nref),e,a,b)
!!$          write(*,'(i3,4f15.5)') i,a,b,e,t
!!$!          if( abs(e - e1) < tol .and. abs(a - a1) < tol .and. abs(b - b1) < tol ) exit
!!$          lambda = 2*lambda
!!$          e1 = e
!!$          a1 = a
!!$          b1 = b
!!$       end do
!!$
!!$!       call part(e,a,b,t)
!!$       par_t = t
!!$
!!$
!!$       write(*,'(a,4f15.5)') 'abs:',a,b,c,t
!!$
!!$       exit
!!$       tx = t
!!$
!!$!       stop 0
!!$
!!$goto 45
!!$       lambda = 0.1
!!$       do i = 1,10
!!$          call commons(e,a,b)
!!$          write(*,*) real(e),real(a),real(b)
!!$          lambda = 2*lambda
!!$       end do
!!$
!!$       open(1,file='/tmp/5')
!!$       do i = nhwhm,naper
!!$          write(1,*) raper(i),(1-e-a/raper(i)-b/raper(i)**2)
!!$       end do
!!$       close(1)
!!$!stop 0
!!$45 continue
!!$
!!$       ! estimate background of individual stars
!!$       do i = 1,nstar
!!$          cts_ = cts(i,:)
!!$          dcts_ = dcts(i,:)
!!$!          call background(e,a,b,back(i),tx(i))
!!$!          call background3(e,a,b,t,back(i))
!!$!          back(i) = xbackground(e,a,b,t)
!!$!          write(*,*) i,real(back(i))
!!$!stop 0
!!$       end do
!!$!       stop 0
!!$!       call rmean(tx,t,dt)
!!$!       write(*,*) t,dt
!!$!       t = t1
!!$
!!$goto 44
!!$
!!$       call background2(e,a,b,t,back)
!!$       write(*,*) real(t),real(back)
!!$       tx = t
!!$
!!$       open(1,file='/tmp/4')
!!$       do i = 1,nstar
!!$          call the_constructor(tx(i),back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$          do j = 1,naper
!!$             write(1,*) raper(j),grows(i,j)
!!$          end do
!!$       end do
!!$       close(1)
!!$       open(1,file='/tmp/5')
!!$       do i = nhwhm,naper
!!$          call rinit(grows(:,i),grow_(i),dgrow_(i))
!!$          write(1,*) raper(i),grow_(i),dgrow_(i),grow_(i)-(1+e-a/raper(i)-b/raper(i)**2)
!!$       end do
!!$       close(1)
!!$
!!$!stop 0
!!$       tx = t
!!$
!!$       !       call part(e,a,b,t,dt)
!!$
!!$       ! estimate parameters common to all stars: a,b,e
!!$       lambda = 0.1
!!$       do i = 1,10
!!$          call commons(e,a,b)
!!$          write(*,*) real(e),real(a),real(b)
!!$          lambda = 2*lambda
!!$       end do
!!$
!!$44 continue
!!$
!!$       db = sum(abs(back - back1)/cts(:,nhwhm))/nstar
!!$!       dt = sum(abs(tx - tx1))/nstar
!!$       dt = abs(t - t1)
!!$
!!$       write(*,'(i2,3f10.5,5x,1p,2g10.2)') iter,a,b,t,db,dt
!!$
!!$       if( db < tol .and. dt < tol ) exit
!!$
!!$!       do i = 1,nstar
!!$!          do j = 8,naper
!!$!             if( abs(cts(i,j)/cts(i,naper) - (1-a/raper(j)-b/raper(j)**2)) > 0.01 ) then
!!$!                cts(i,j) = cts(i,naper)*(1-a/raper(j)-b/raper(j)**2)
!!$!             end if
!!$!          end do
!!$!       end do
!!$
!!$       lambda = 2*lambda
!!$       back1 = back
!!$       tx1 = tx
!!$       dt = abs(t1-t)
!!$       t1 = t
!!$!       stop
!!$    end do
!!$
!!$
!!$    do i = 1,nstar
!!$       call the_constructor(t,back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$    end do
!!$    allocate(fvec(naper*nstar),dvec(naper*nstar))
!!$    n = 0
!!$!    do i = (nhwhm+naper)/2,naper
!!$    do i = nref,naper
!!$       n1 = n + 1
!!$       n = n + nstar
!!$       n2 = n
!!$       d = asymptote(raper(i),a,b,c)
!!$       fvec(n1:n2) = abs(grows(:,i) - d) !/ dgrows(:,i)
!!$       dvec(n1:n2) = dgrows(:,i)
!!$!       d = (1-e-a/raper(i)-b/raper(i)**2)
!!$    end do
!!$!    write(*,*) n,nstar,real(fvec(1:n))
!!$    call rinit(fvec(1:n),d,s)
!!$    call rwmean(fvec(1:n),dvec(1:n),d,s,sig)
!!$!    sig = d / 0.6745
!!$    write(*,*) 'sig:',d, sig
!!$    deallocate(fvec,dvec)
!!$
!!$    back = 0
!!$    type = 1
!!$    lambda = 1
!!$    e1 = 100
!!$    a1 = 100
!!$    b1 = 100
!!$    do i = 1,20
!!$       call commons1(e,a,b)
!!$       !       call part(e,a,b,t)
!!$!       t = asymptote(raper(nref),e,a,b)
!!$       write(*,'(i3,4f15.5)') i,a,b,e,t
!!$!       if( abs(e - e1) < tol .and. abs(a - a1) < tol .and. abs(b - b1) < tol ) exit
!!$       lambda = 2*lambda
!!$       e1 = e
!!$       a1 = a
!!$       b1 = b
!!$    end do
!!$!       t = 1 - a/raper(nhwhm) - b/raper(nhwhm)**2
!!$
!!$!    call part(e,a,b,t)
!!$    write(*,*) 'commons1:',real(a),real(b),real(t),real(e)
!!$
!!$    ! estimate background of individual stars
!!$    do i = 1,nstar
!!$       cts_ = cts(i,:)
!!$       dcts_ = dcts(i,:)
!!$!       back(i) = ybackground(e,a,b,t)
!!$!       write(*,*) i,real(back(i))
!!$    end do
!!$
!!$99  continue
!!$
!!$    ! final curve
!!$    do i = 1,nstar
!!$!       call the_constructor(e,a,cts(i,:),dcts(i,:),back(i),grows(i,:),dgrows(i,:),valid)
!!$!       call the_constructor(tx(i),back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$       call the_constructor(tx(i),back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$!       write(*,*) i,real(grows(i,:)),real(dgrows(i,:))
!!$    end do
!!$    do i = 1,naper
!!$!       call rmean(grows(:,i),grow(i),dgrow(i))
!!$       call rwmean(grows(:,i),dgrows(:,i),grow(i),dgrow(i))
!!$!       write(*,*) i,real(grows(:,i))
!!$    end do
!!$!    do i = (nhwhm+naper)/2,naper
!!$!       grow(i) = (1-a/raper(i)-b/raper(i)**2)
!!$    !    end do
!!$
!!$999 continue

    ! check for identical curves
!    where( abs(dgrow) < epsilon(dgrow) )
!       dgrow = 1e-6
!    end where

!    open(1,file='/tmp/8')
!    do i = 1,naper
!       write(1,*) raper(i),grow(i)-(1-e-a/raper(i)-b/raper(i)**2),grow(i),dgrow(i),&
!            (1-e-a/raper(i)-b/raper(i)**2)
!    end do
!    close(1)

    if( debug ) then
       do i = 1,naper
          write(*,*) i,real(raper(i)),real(curve(i)),real(dcurve(i))
       end do
    end if

!    do i = nref+1,naper
!    t1 = (grow(naper) - grow(nref))*0.998
!    b1 = a*(1/raper(nref) -1/raper(naper)) +  2*b*(1/raper(nref)**2 - 1/raper(naper)**2)
!    write(*,'(4g15.5)') t1-b1,t1,b1
!   stop 0

!    deallocate(raper,aper,cts,dcts,back,back1,grows,dgrows,tx,tx1,grow_,dgrow_, &
!         cts_,dcts_)
    deallocate(raper,aper,radius,flux,betas,phi,cts,dcts,grows,dgrows,grow,dgrow,sigs)

  end subroutine growmodel



  subroutine back_update(cts,dcts,flux,a,b,beta,phi)

    use NelderMead
    use oakleaf
    use minpacks

    real(dbl), dimension(:), intent(in) :: cts, dcts
    real(dbl), intent(in) :: flux,a,b
    real(dbl), intent(out) :: beta, phi
    real(dbl), dimension(2) :: p,dp
    real(dbl), dimension(size(cts)) :: res
    real(dbl) :: pmin, s, g
    integer :: i,ifault,info, nprint

    do i = nref,naper
       grow(i) = cts(i)/flux - (1 - (a + b/raper(i))/raper(i))
!       write(*,'(3e15.5,a)') raper(i),grow_(i),dgrow_(i),' EE'
    end do
    s = median(dcts(nref:naper))
    sig = s / flux + approxerr

    p = 0
    dp = 1e-6

    call nelmin1(backfun,p,dp,pmin,ifault)

    if( ifault == 0 ) then
       phi = p(1)
       beta = p(2)
    else
       if( debug ) write(*,*) 'back_init:',ifault,real(p)
       phi = 0
       beta = 0
       return
    end if

    do i = nref, naper
       g = phi + beta*aper(i)
       res(i) = abs(grow(i) - g) / s
    end do
    sig = median(res(nref:naper)) / 0.6745 + approxerr

!    write(*,*) real(p),real(sig)

    nprint = 0
    call lmder2(backder,p,epsilon(p),nprint,info)

!    write(*,*) real(p),info
    if( info == 2 ) then
       phi = p(1)
       beta = p(2)
    end if

  end subroutine back_update

  function backfun(p) result(s)

    real(dbl), dimension(:), intent(in) :: p

    real(dbl) :: s,g,beta,phi
    integer :: i

    phi = p(1)
    beta = p(2)

    s = 0
    do i = nref, naper
       g = phi + beta*aper(i)
       s = s + abs(grow(i) - g) / sig
!       write(*,*) radius(i),f,x,grow_(i),dgrow_(i)
    end do

  end function backfun

  subroutine backder(m,np,p,fvec,fjac,ldfjac,iflag)

    use oakleaf

    integer, intent(in) :: m,np,ldfjac
    integer, intent(inout) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(out) :: fvec
    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
    real(dbl), dimension(:), allocatable :: rs,f,df,r
    integer :: n,i,l
    real(dbl) :: phi, beta

    if( iflag == 0 ) then
       write(*,'(6g12.3)') p,fvec

       if( .true. ) then
          write(*,*) ' jac:',real(fjac(1,:))
          write(*,*) ' jac:',real(fjac(2,:))
          write(*,*) ' jac:',real(fjac(3,:))
       end if

       return
    end if

    n = naper - nref + 1
    allocate(rs(n),f(n),df(n),r(n))

    phi = p(1)
    beta = p(2)

    do i = 1,n
       l = i + nref - 1
       r(i) = raper(l)
       rs(i) = (phi + beta*aper(l) - grow(l)) / sig
       f(i) = huber(rs(i))
    end do
!    f = r

    if( iflag == 1 ) then

       fvec(1) = sum(f)
       fvec(2) = sum(f*r**2)*pi
       fvec = fvec / sig

    else if( iflag == 2 ) then

       df = dhuber(rs)
!       df = 1

       fjac(1,1) = sum(df)
       fjac(1,2) = sum(df*r**2)*pi
       fjac(2,2) = sum(df*r**4)*pi**2

       fjac(2,1) = fjac(1,2)
       fjac = fjac / sig**2

    end if

    deallocate(rs,f,df,r)

  end subroutine backder


  subroutine arc_update(flux,phi,beta,a,b)

    use oakleaf
    use neldermead
    use minpacks

    real(dbl), dimension(:), intent(in) :: flux, phi, beta
    real(dbl), intent(in out) :: a,b
    real(dbl), dimension(3) :: p,dp,p1
    real(dbl), dimension(:), allocatable :: res
    integer :: ifault,info,nprint,i,l,n
    real(dbl) :: s,g,e

    do n = 1,nstar
       do i = nref,naper
          grows(n,i) = cts(n,i) / flux(n) - phi(n) - beta(n)*aper(i)
       end do
       s = median(dcts(n,nref:naper))
       sigs(n) = s / flux(n) + approxerr
!       write(*,*) n,s,sigs(n)
    end do

    p(1) = a / raper(naper)
    p(2) = b / raper(naper)**2
    p(3) = 0

    dp = 1e-4
    lambda = 1e-4
    do i = 1,50
       call nelmin1(arcupdate,p,dp,s,ifault)
       if( debug ) &
            write(*,'(2f12.5,es12.2,i3)') p(1)*raper(naper),p(2)*raper(naper)**2,&
            p(3),ifault
       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) ) exit
       p1 = p
       lambda = 2*lambda
    end do
!    a = p(1)*raper(naper)
!    b = p(2)*raper(naper)**2
!    return

    if( ifault == 0 ) then
       a = p(1)
       b = p(2)
       e = p(3)
    else
       a = p(1)*raper(naper)
       b = p(2)*raper(naper)**2
       e = p(3)
       return
    end if

    l = 0
    allocate(res(naper*nstar))
    do n = 1,nstar
       do i = nref, naper
          g = 1 - e - (a + b/radius(i))/radius(i)
          l = l + 1
          res(l) = abs(grows(n,i) - g) / sigs(n)
       end do
    end do
    s = median(res(1:l))
    sig = s / 0.6745 + approxerr
    deallocate(res)
    if( debug ) write(*,*) 'sig,s:',real(sig),real(s)

    lambda = 1e-4
    do i = 1,50
       call nelmin1(archuber,p,dp,s,ifault)
       if( debug ) &
            write(*,'(2f12.5,es12.2,i3)') p(1)*raper(naper),p(2)*raper(naper)**2,&
            p(3),ifault
       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) ) exit
       p1 = p
       lambda = 2*lambda
    end do
    if( ifault == 0 ) then
       a = p(1)*raper(naper)
       b = p(2)*raper(naper)**2
       return
    end if

    nprint = 0
    lambda = 1e-4
    do i = 1,50*0
       call lmder2(arcder,p,epsilon(p),nprint,info)
       !    call lmdif2(backdif,p,epsilon(p),nprint,info)
       write(*,'(2f12.5,es12.2,i3)') p(1)*raper(naper),p(2)*raper(naper)**2,&
            p(3),info
       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) ) exit
       p1 = p
       lambda = 2*lambda
    end do

    if( info == 2 ) then
       a = p(1)*raper(naper)
       b = p(2)*raper(naper)**2
    end if

  end subroutine arc_update


 function arcupdate(p) result(s)

    real(dbl), dimension(:), intent(in) :: p

    real(dbl) :: s,g,a,b,e
    integer :: i,n

    a = p(1)
    b = p(2)
    e = p(3)

    s = 0
    do n = 1,nstar
       do i = nref, naper
          g = 1 - e - (a + b/radius(i))/radius(i)
          s = s + abs(g - grows(n,i)) / sigs(n)
       end do
    end do
    s = s / (nstar*(naper - nref + 1)) + lambda*(max(0.0,-a) + max(0.0,-b) + max(0.0,-e) + abs(e))

  end function arcupdate

  function archuber(p) result(s)

    use oakleaf

    real(dbl), dimension(:), intent(in) :: p

    real(dbl) :: s,g,a,b,e
    integer :: i,n

    a = p(1)
    b = p(2)
    e = p(3)

    s = 0
    do n = 1,nstar
       do i = nref, naper
          g = 1 - e - (a + b/radius(i))/radius(i)
          !          s = s + abs(grows(n,i) - g) / sigs(n)
          s = s + ihuber((g - grows(n,i)) / (sigs(n)*sig))
!          write(*,*) (g - grows(n,i)) / sigs(n),n,i,'666'
       end do
    end do
!    s = s / (nstar*(naper - nref + 1)) + lambda*(max(0.0,-a)**2 + max(0.0,-b)**2 + max(0.0,-e)**2 + e**2)
    s = s / (nstar*(naper - nref + 1)) &
         + lambda*(smooth_max(a)**2 + smooth_max(b)**2 + smooth_max(e)**2)
!stop
  end function archuber


  subroutine arcder(mp,np,p,fvec,fjac,ldfjac,iflag)

    use oakleaf

    integer, intent(in) :: mp,np,ldfjac
    integer, intent(inout) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(mp), intent(out) :: fvec
    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
    real(dbl), dimension(:,:), allocatable :: rs,f,df
    real(dbl), dimension(:), allocatable :: r,s
!    real(dbl), dimension(3,3) :: dfjac
    integer :: n,m,i,j,l
    real(dbl) :: a,b,e

    if( iflag == 0 ) then
       write(*,'(6g12.3)') p,fvec

       if( .true. ) then
          write(*,*) ' jac:',real(fjac(1,:))
          write(*,*) ' jac:',real(fjac(2,:))
          write(*,*) ' jac:',real(fjac(3,:))

!          call difjac(p(1),p(2),dfjac)
!          write(*,*) 'djac:',dfjac(1,:)
!          write(*,*) 'djac:',dfjac(2,:)

       end if

       return
    end if

    m = naper - nref + 1
    n = nstar
    allocate(rs(n,m),f(n,m),df(n,m),r(m),s(n))
    a = p(1)
    b = p(2)
    e = p(3)

    r = radius(nref:naper)

    do n = 1,nstar
       do l = 1, m
          i = l + nref - 1
          s(n) = sigs(n)*sig
          rs(n,l) = (1 - e - (a + b/r(l))/r(l) - grows(n,i)) / s(n)
          f(n,l) = huber(rs(n,l))
       end do
    end do
    f = rs

    if( iflag == 1 ) then

       fvec = 0
       do n = 1,nstar
          fvec(1) = fvec(1) + sum(f(n,:)/r) / s(n)
          fvec(2) = fvec(2) + sum(f(n,:)/r**2) / s(n)
          fvec(3) = fvec(3) + sum(f(n,:)) / s(n)
       end do
       fvec = - fvec - lambda*(/ smooth_max(a), smooth_max(b), smooth_max(e) /)

!       write(*,*) real(fvec)

    else if( iflag == 2 ) then

       do n = 1, nstar
          df(n,:) = dhuber(rs(n,:))
       end do
       df = 1

       fjac = 0
       do n = 1,nstar
          fjac(1,1) = fjac(1,1) + sum(df(n,:)/r**2) / s(n)**2
          fjac(1,2) = fjac(1,2) + sum(df(n,:)/r**3) / s(n)**2
          fjac(1,3) = fjac(1,3) + sum(df(n,:)/r) / s(n)**2
          fjac(2,2) = fjac(2,2) + sum(df(n,:)/r**4) / s(n)**2
          fjac(2,3) = fjac(2,3) + sum(df(n,:)/r**2) / s(n)**2
          fjac(3,3) = fjac(3,3) + sum(df(n,:)) / s(n)**2
       end do
       do i = 1,3
          do j = 1,i-1
             fjac(i,j) = fjac(j,i)
          end do
       end do

       fjac(1,1) = fjac(1,1) + lambda*smooth_step(a)
       fjac(2,2) = fjac(2,2) + lambda*smooth_step(b)
       fjac(3,3) = fjac(3,3) + lambda*smooth_step(e)

    end if

    deallocate(rs,f,df,r,s)

  end subroutine



    subroutine arc_update1(flux,phi,beta,b)

    use neldermead
    use minpacks
    use oakleaf

    real(dbl), dimension(:), intent(in) :: flux, phi, beta
    real(dbl), intent(in out) :: b
    real(dbl), dimension(2) :: p,dp,p1
    real(dbl), dimension(:), allocatable :: res
    integer :: ifault,info,nprint,i,l,n
    real(dbl) :: s,g,e

    do n = 1,nstar
       do i = nref,naper
          grows(n,i) = cts(n,i) / flux(n) - phi(n) - beta(n)*aper(i)
       end do
       s = median(dcts(n,nref:naper))
       sigs(n) = s / flux(n) + approxerr
!       write(*,*) n,s,sigs(n)
    end do

    p(1) = b / raper(naper)**2
    p(2) = 0

    dp = 1e-4
    lambda = 1e-4
    do i = 1,50
       call nelmin1(arcupdate1,p,dp,s,ifault)
       if( debug ) &
            write(*,'(f12.5,es12.2,i3)') p(1)*raper(naper)**2,p(2),ifault
       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) ) exit
       p1 = p
       lambda = 2*lambda
    end do

    if( ifault == 0 ) then
       b = p(1)
       e = p(2)
    else
       b = p(1)*raper(naper)**2
       e = p(2)
       return
    end if

    l = 0
    allocate(res(naper*nstar))
    do n = 1,nstar
       do i = nref, naper
          g = 1 - e - b/radius(i)**2
          l = l + 1
          res(l) = abs(grows(n,i) - g) / sigs(n)
       end do
    end do
    s = median(res(1:l))
    sig = s / 0.6745 + approxerr
    deallocate(res)
    if( debug ) write(*,*) 'sig,s:',real(sig),real(s)

    lambda = 1e-4
    do i = 1,50
       call nelmin1(archuber1,p,dp,s,ifault)
       if( debug ) &
            write(*,'(f12.5,es12.2,i3)') p(1)*raper(naper)**2,p(2),ifault
       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) ) exit
       p1 = p
       lambda = 2*lambda
    end do
    b = p(1)*raper(naper)**2
    return

!!$    nprint = 0
!!$    lambda = 1e-4
!!$    do i = 1,50*0
!!$       call lmder2(arcder,p,epsilon(p),nprint,info)
!!$       !    call lmdif2(backdif,p,epsilon(p),nprint,info)
!!$       write(*,'(2f12.5,es12.2,i3)') p(1)*raper(naper),p(2)*raper(naper)**2,&
!!$            p(3),info
!!$       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) ) exit
!!$       p1 = p
!!$       lambda = 2*lambda
!!$    end do
!!$
!!$    if( info == 2 ) then
!!$       a = p(1)*raper(naper)
!!$       b = p(2)*raper(naper)**2
!!$    end if

  end subroutine arc_update1


 function arcupdate1(p) result(s)

    real(dbl), dimension(:), intent(in) :: p

    real(dbl) :: s,g,b,e
    integer :: i,n

    b = p(1)
    e = p(2)

    s = 0
    do n = 1,nstar
       do i = nref, naper
          g = 1 - e - b/radius(i)**2
          s = s + abs(g - grows(n,i)) / sigs(n)
       end do
    end do
    s = s / (nstar*(naper - nref + 1)) + lambda*(max(0.0,-b) + max(0.0,-e))

  end function arcupdate1

  function archuber1(p) result(s)

    use oakleaf

    real(dbl), dimension(:), intent(in) :: p

    real(dbl) :: s,g,b,e
    integer :: i,n

    b = p(1)
    e = p(2)

    s = 0
    do n = 1,nstar
       do i = nref, naper
          g = 1 - e - b/radius(i)**2
          !          s = s + abs(grows(n,i) - g) / sigs(n)
          s = s + ihuber((g - grows(n,i)) / (sigs(n)*sig) )
!          write(*,*) (g - grows(n,i)) / sigs(n),n,i,'666'
       end do
    end do
!    s = s / (nstar*(naper - nref + 1)) + lambda*(max(0.0,-b)**2 + max(0.0,-e)**2)
    s = s / (nstar*(naper - nref + 1)) &
         + lambda*(smooth_max(b)**2 + smooth_max(e)**2)

  end function archuber1





  subroutine arc_init(flux,b)

    ! initialise arc of grow curve

    use fmm
    use oakleaf

    real(dbl), dimension(:), intent(in) :: flux
    real(dbl), intent(out) :: b

    real(dbl) :: tol,bmax,bmin,s,x
    integer :: n,i

    bmax = 0
    do n = 1,nstar
       do i = naper-1,nref,-1
          x = 1/radius(i+1)**1 - 1/radius(i)**1
          grows(n,i) = (cts(n,i+1) - cts(n,i)) / flux(n)
          s = grows(n,i) / x
          if( s > bmax ) bmax = s
       end do
       s = median(dcts(n,nref:naper)) / 0.6745
       sigs(n) = s / flux(n) + approxerr
       ! added function approximation error
    end do

    tol = 1e-7
    bmin = 0 !epsilon(b)
!    bmax = 1.0/20.**2!00
    !    bmax = maxval(grows(n,nref:naper-1)/
    bmax = max(bmax,bmin+tol)
    b = fmin(bmin,bmax,arcfun,tol)
!    write(*,*) bmin,b,bmax
    b = b*raper(naper)**2

  end subroutine arc_init

  function arcfun(b) result(s)

    real(dbl), intent(in) :: b

    real(dbl) :: s,g,x,y
    integer :: i,n

    s = 0
    do n = 1,nstar
       do i = nref, naper - 1
          x = 1/radius(i+1)**1 - 1/radius(i)**1
          y = 1/radius(i+1)**2 - 1/radius(i)**2
          g = -b*x
          s = s + abs(grows(n,i) - g) / sigs(n)
          !          write(*,*) grows(n,i) - g,(grows(n,i) - g) / sigs(n),i,n
!          write(*,*) -x,y,grows(n,i),i,n,'RR'
       end do
    end do
    s = s / (nstar*(naper - nref))

  end function arcfun


  subroutine arc_init1(flux,a,b)

    ! initialise arc of grow curve

    use neldermead
    use oakleaf

    real(dbl), dimension(:), intent(in) :: flux
    real(dbl), intent(out) :: a,b

    real(dbl) :: s,x
    real(dbl), dimension(3) :: p,dp,p1
    integer :: n,i,ifault

    do n = 1,nstar
       do i = naper-1,nref,-1
          x = 1/radius(i+1)**2 - 1/radius(i)**2
          grows(n,i) = (cts(n,i+1) - cts(n,i)) / flux(n)
       end do
       s = median(dcts(n,nref:naper)) / 0.6745
       sigs(n) = s / flux(n) + approxerr
       ! added function approximation error
    end do

    p(1) = 0
    p(2) = 0
    p(3) = 0

    dp = 1e-4
    lambda = 1e-4
    p1 = huge(p)
    do i = 1,50
       call nelmin1(arcinit1,p,dp,s,ifault)
!       write(*,*) real(p)
       if( all(p > 0) .and. all(abs(p - p1) < 1e-7) ) exit
       lambda = 2*lambda
       p1 = p
    end do

    a = p(1)*raper(naper)
    b = p(2)*raper(naper)**2

  end subroutine arc_init1

  function arcinit1(p) result(s)

    real(dbl), dimension(:), intent(in) :: p

    real(dbl) :: s,g,x,y,b,a,e
    integer :: i,n

    a = p(1)
    b = p(2)
    e = p(3)

    s = 0
    do n = 1,nstar
       do i = nref, naper - 1
          x = 1/radius(i+1)    - 1/radius(i)
          y = 1/radius(i+1)**2 - 1/radius(i)**2
          g = -b*y - a*x + e
          s = s + abs(grows(n,i) - g) / sigs(n)
          !          write(*,*) grows(n,i) - g,(grows(n,i) - g) / sigs(n),i,n
!          write(*,*) -x,grows(n,i),i,n,'RR'
       end do
    end do
    s = s / (nstar*(naper - nref)) + lambda*(max(0.0,-a) + max(0.0,-b)+max(0.0,-e))
!    stop

  end function arcinit1


  function smooth_step(x) result(f)

    ! negative logistics function: 1 - 1/(1 + exp(-x)):  _
    !         it is smooth equivalent of step function:   |_

    real(dbl), intent(in) :: x
    real(dbl), parameter :: s = epsilon(1.0) ! the same as in isigmoid
    real(dbl), parameter :: ehuge = log(huge(1.0)) / 10
    real(dbl) :: f

    if( abs(x)/s < ehuge ) then
       f = 1 / (1 + exp((x-s)/s))
    else
       if( x < 0 ) then
          f = 1
       else
          f = 0
       end if
    end if

  end function smooth_step

  function smooth_max(x) result(f)

    ! integral of negative logistics function: int(1 - 1/(1 + exp(-x))) = log(1+exp(-x))
    ! it is smooth equivalent of max(0,-x): \_

    real(dbl), intent(in) :: x
    real(dbl), parameter :: s = epsilon(1.0) ! the same as in nsigmoid
    real(dbl), parameter :: ehuge = log(huge(1.0)) / 10
    real(dbl) :: f

    if( abs(x)/s < ehuge ) then
       f = log(1 + exp(-(x-s)/s))
    else
       if( x < 0 ) then
          f = -x
       else
          f = 0
       end if
    end if

  end function smooth_max



!---------------------------------------------------------------------------


!!$  subroutine arc_inits1(flux0,b,e,beta)
!!$
!!$    ! initialise arc of grow curve
!!$
!!$    use minpacks
!!$    use robustmean
!!$    use neldermead
!!$
!!$    real(dbl), dimension(:), intent(in) :: flux0
!!$    real(dbl), intent(in out) :: b,e
!!$    real(dbl), dimension(:), intent(in out) :: beta
!!$    real(dbl), dimension(3) :: p,dp,p1
!!$    real(dbl), dimension(:), allocatable :: dd,x,y,z
!!$    real(dbl) :: c,d,s,g,fmax,fmin,a
!!$    integer :: info,nprint,i,n,ifault,m,l
!!$
!!$!    nstar = 10
!!$
!!$    fmin = minval(cts(:,nref)/flux0)
!!$    fmax = maxval(cts(:,naper)/flux0)
!!$    c = fmax - fmin
!!$!    c = 1
!!$    ! remove c
!!$
!!$    do n = 1,nstar
!!$       do i = naper-1,nref,-1
!!$          grows(n,i) = (cts(n,i+1) - cts(n,i)) / flux0(n) / c
!!$!          grows(n,i) = sqrt(cts(n,i+1)/flux0(n)) - sqrt(cts(n,i) / flux0(n))
!!$!          write(*,*) i,i+1,grows(n,i),'QQ'
!!$       end do
!!$       call rinit(dcts(n,nref:naper)/c,d,s)
!!$!       d = d / 0.6745
!!$!       write(*,*) dcts(n,nhwhm) / flux0(n),sqrt(cts(1,nhwhm))/flux0(n),d/flux0(n)
!!$       !       dgrows(n,:) = sqrt(dcts(n,nref)**2 + cts(n,nref)) / flux0(n)
!!$       dgrows(n,:) = d / flux0(n) + approxerr
!!$       ! added function approximation error
!!$!       write(*,*) real(d/flux0(n)),real(dgrows(n,1))
!!$    end do
!!$!    call rinit(pack(abs(grows(:,nref:naper-1)),.true.),sig,s)
!!$!    dgrows(:,nhwhm) = sqrt(sig**2 + dgrows(:,nhwhm)**2)
!!$
!!$!    fmin = minval(grows(:,nref:naper-1))
!!$!    fmax = maxval(grows(:,nref:naper-1))
!!$!    grows(:,nref:naper-1) = (grows(:,nref:naper-1) - fmin) / d
!!$!    dgrows(:,nref:naper-1) = dgrows(:,nref:naper-1) / d
!!$
!!$    p(1) = a / raper(naper)
!!$    p(2) = b / raper(naper)**2
!!$!    p(3) = 0
!!$!    p(4) = 0
!!$!    p(3:) = beta*(pi*raper(naper)**2)/flux0
!!$
!!$    p = 0
!!$
!!$    dp = 1e-4
!!$!    dp(3) = 1e-7
!!$!    dp(1:2) = 1e-4
!!$    lambda = 1e-4
!!$!    back = 0
!!$    do i = 1,50
!!$       call nelmin1(arcfuns2,p,dp,s,ifault)
!!$       write(*,'(2f15.5,es10.2,f15.5,f10.3)') p(1)*raper(naper)**2*c,p(2)*c!raper(naper)*c, &
!!$!            p(2)*c!,sum(abs(p(1:)))/nstar
!!$       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) ) exit
!!$!       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) .and. &
!!$!            sum(abs(back-p(3:)))/nstar < 0.001 ) exit
!!$!       back = p(3:)
!!$       p1 = p
!!$       lambda = 2*lambda
!!$    end do
!!$!    e = p(2)
!!$!    b = p(1)
!!$
!!$    m = naper - nref
!!$    allocate(dd(m),x(m),y(m),z(m))
!!$
!!$    do l = 1, m
!!$       i = l + nref - 1
!!$       x(l) = 1/radius(i) -    1/radius(i+1)
!!$       y(l) = 1/radius(i)**2 - 1/radius(i+1)**2
!!$!       z(l) = radius(i+1)**2 - radius(i)**2
!!$    end do
!!$
!!$    do n = 1,nstar
!!$       do i = 1,m
!!$          g = a*x(i) + b*y(i)
!!$          dd(i) = abs(grows(n,i+nref-1) - g)
!!$!          write(*,*) g,grows(n,i+nref-1)
!!$       end do
!!$       call rinit(dd,d,s)
!!$!       write(*,*) real(dd)
!!$!       stop
!!$       dgrows(n,:) = d / 0.6745 + approxerr
!!$    end do
!!$    deallocate(dd,x,y,z)
!!$    !stop
!!$
!!$    ! je to rozumne vzheledem k systematickym chybam?
!!$    nprint = 0
!!$    p1 = p
!!$    lambda = lambda / 10
!!$    do i = 1,50*0
!!$       !call lmder2(initfun,p,epsilon(p),nprint,info)
!!$       call lmdif2(inidfun1,p,epsilon(p),nprint,info)
!!$       write(*,'(e15.5,3f19.7)') p(1)*raper(naper)*c,p(2)*raper(naper)**2*c, &
!!$            p(2)/(pi*raper(naper)**2)*c
!!$       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-5) ) exit
!!$       lambda = 2*lambda
!!$       p1 = p
!!$!       back = p(3:)
!!$    end do
!!$
!!$    a = p(1)*raper(naper)*c
!!$    b = p(2)*raper(naper)**2*c
!!$    a = 0
!!$    b = p(1)*raper(naper)**2*c
!!$    a = p(2)*c
!!$    b = p(1)*raper(naper)**2*c
!!$
!!$
!!$  end subroutine arc_inits1
!!$
!!$  function arcfuns2(p) result(s)
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(:), allocatable :: x,y,z,sig
!!$    real(dbl) :: s,a,b,g,t,f2,r,sbeta,q
!!$    integer :: l,i,n,m
!!$
!!$    m = naper - nref
!!$    n = nstar
!!$    allocate(x(m),y(m),z(m),sig(n))
!!$
!!$    b = p(1)
!!$    a = p(2)
!!$    q = p(3)
!!$!    sbeta = p(2)
!!$    !    sbeta = 0
!!$!    q = 0
!!$    sig = dgrows(1:nstar,nhwhm)
!!$
!!$    do l = 1, m
!!$       i = l + nref - 1
!!$       x(l) = 1/radius(i) -    1/radius(i+1)
!!$       y(l) = 1/radius(i)**2 - 1/radius(i+1)**2
!!$       z(l) = radius(i+1)**2 - radius(i)**2
!!$    end do
!!$
!!$    s = 0
!!$    do n = 1,nstar
!!$       do i = 1, m
!!$!          g = b*y(i) + sbeta*z(i) + q!*radius(i+nref-1)
!!$          g = a*x(i) + b*y(i) + q
!!$!          write(*,'(6g12.2)') x(i),y(i),grows(n,i+nref-1) - g,sig(n)
!!$          s = s + abs(grows(n,i+nref-1) - g) / sig(n)
!!$       end do
!!$    end do
!!$
!!$    s = s / (m*n) + lambda*(max(0.0,-b) + max(0.0,-a) + abs(sbeta)*0)
!!$!    s = s / (m*n) + lambda*max(0.0,-b)
!!$
!!$    deallocate(x,y,z,sig)
!!$
!!$ end function arcfuns2
!!$
!!$
!!$ function arcfuns1(p) result(s)
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(:), allocatable :: x,y,z,sig
!!$    real(dbl) :: s,a,b,g,t,f2,r,sbeta,q
!!$    integer :: l,i,n,m
!!$
!!$    m = naper - nref
!!$    n = nstar
!!$    allocate(x(m),y(m),z(m),sig(n))
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    sbeta = p(3)
!!$    q = p(4)
!!$!    beta = p(3:)
!!$    sig = dgrows(1:nstar,nhwhm)
!!$
!!$    do l = 1, m
!!$       i = l + nref - 1
!!$       x(l) = 1/radius(i) -    1/radius(i+1)
!!$       y(l) = 1/radius(i)**2 - 1/radius(i+1)**2
!!$       z(l) = radius(i+1)**2 - radius(i)**2
!!$    end do
!!$
!!$    s = 0
!!$    do n = 1,nstar
!!$       do i = 1, m
!!$          g = a*x(i) + b*y(i) + sbeta*z(i) + q*radius(i+nref-1)
!!$!          write(*,'(6g12.2)') x(i),y(i),grows(n,i+nref-1) - g,sig(n)
!!$          s = s + abs(grows(n,i+nref-1) - g) / sig(n)
!!$       end do
!!$    end do
!!$
!!$    f2 = 0
!!$    do i = nref+1,naper
!!$       r = (radius(i) + radius(i-1))/2
!!$       f2 = f2 + 2*(a + 3*b/r)/r**3 * (radius(i) - radius(i-1))
!!$    end do
!!$
!!$    s = s / (m*n) + lambda*(max(0.0,-a) + max(0.0,-b) + max(0.0,-q))
!!$
!!$    deallocate(x,y,z,sig)
!!$
!!$ end function arcfuns1
!!$
!!$  subroutine inidfun1(mm,np,p,fvec,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: mm,np
!!$    integer, intent(in out) :: iflag
!!$    real(dbl), dimension(np), intent(in) :: p
!!$    real(dbl), dimension(mm), intent(out) :: fvec
!!$    real(dbl), dimension(:), allocatable :: beta,x,y,z,s
!!$    real(dbl), dimension(:,:), allocatable :: rs,f
!!$    real(dbl) :: a,b,g
!!$    integer :: n,m,i,j,l
!!$
!!$    if( iflag == 0 ) then
!!$       write(*,'(2g13.5)') p
!!$       write(*,'(2g13.5)') fvec
!!$       return
!!$    end if
!!$
!!$    m = naper - nref
!!$    n = nstar
!!$    allocate(rs(n,m),f(n,m),x(m),y(m),z(m),beta(n),s(n))
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    s = dgrows(1:nstar,nhwhm)
!!$
!!$    do l = 1, m
!!$       i = l + nref - 1
!!$       x(l) = 1/radius(i) -    1/radius(i+1)
!!$       y(l) = 1/radius(i)**2 - 1/radius(i+1)**2
!!$!       z(l) = radius(i+1)**2 - radius(i)**2
!!$    end do
!!$
!!$    do n = 1,nstar
!!$       do i = 1, m
!!$          g = a*x(i) + b*y(i)
!!$          rs(n,i) = (grows(n,i+nref-1) - g) / s(n)
!!$          f(n,i) = huber(rs(n,i))
!!$!          if( n == 1 ) &
!!$!          write(*,'(i2,8f15.5)') i,rs(n,i),f(n,i),g,grows(n,i+nref-1),s(n)
!!$       end do
!!$    end do
!!$    f = rs
!!$
!!$    fvec(1:2) = 0
!!$    do n = 1,nstar
!!$       fvec(1) = fvec(1) + sum(f(n,:)*x) / s(n)
!!$       fvec(2) = fvec(2) + sum(f(n,:)*y) / s(n)
!!$    end do
!!$!    forall( n = 1:nstar )
!!$!       fvec(2+n) = sum(f(n,:)*z) / s(n) - lambda*(beta(n) - back(n))
!!$!    end forall
!!$
!!$    fvec(1) = fvec(1) - lambda*isigmoid(a)
!!$    fvec(2) = fvec(2) - lambda*isigmoid(b)
!!$
!!$    fvec = - fvec
!!$
!!$    deallocate(rs,f,x,y,z,beta,s)
!!$
!!$  end subroutine inidfun1
!!$
!!$  subroutine arc_inits(flux0,beta,a,b)
!!$
!!$    ! initialise arc of grow curve
!!$
!!$    use minpacks
!!$    use robustmean
!!$    use neldermead
!!$
!!$    real(dbl), dimension(:), intent(in) :: flux0
!!$    real(dbl), dimension(:), intent(in) :: beta
!!$    real(dbl), intent(in out) :: a,b
!!$    real(dbl), dimension(2) :: p,dp,p1
!!$    real(dbl), dimension(:), allocatable :: dd,x,y,z
!!$    real(dbl) :: c,d,s,g,fmax,fmin
!!$    integer :: info,nprint,i,n,ifault,m,l
!!$
!!$!    nstar = 10
!!$
!!$    fmin = minval(cts(:,nref)/flux0)
!!$    fmax = maxval(cts(:,naper)/flux0)
!!$    c = fmax - fmin
!!$!    c = 1
!!$    ! remove c
!!$
!!$    do n = 1,nstar
!!$       do i = nref,naper
!!$          grows(n,i) = (cts(n,i) / flux0(n) - beta(n)*aper(i))  / c
!!$!          grows(n,i) = ((cts(n,i+1) - cts(n,i)) / flux0(n) - beta(n)*(aper(i+1) - aper(i)))  / c
!!$!          grows(n,i) = sqrt(cts(n,i+1)/flux0(n)) - sqrt(cts(n,i) / flux0(n))
!!$!                    write(*,*) raper(i),grows(n,i),n,'ee'
!!$       end do
!!$       call rinit(dcts(n,nref:naper)/c,d,s)
!!$!       d = d / 0.6745
!!$!       write(*,*) dcts(n,nhwhm) / flux0(n),sqrt(cts(1,nhwhm))/flux0(n),d/flux0(n)
!!$       !       dgrows(n,:) = sqrt(dcts(n,nref)**2 + cts(n,nref)) / flux0(n)
!!$       dgrows(n,:) = d / flux0(n) + approxerr
!!$       ! added function approximation error
!!$!       write(*,*) real(d/flux0(n)),real(dgrows(n,1))
!!$    end do
!!$!    call rinit(pack(abs(grows(:,nref:naper-1)),.true.),sig,s)
!!$!    dgrows(:,nhwhm) = sqrt(sig**2 + dgrows(:,nhwhm)**2)
!!$
!!$!    fmin = minval(grows(:,nref:naper-1))
!!$!    fmax = maxval(grows(:,nref:naper-1))
!!$!    grows(:,nref:naper-1) = (grows(:,nref:naper-1) - fmin) / d
!!$!    dgrows(:,nref:naper-1) = dgrows(:,nref:naper-1) / d
!!$
!!$!    p(1) = a / raper(naper) / c
!!$!    p(2) = b / raper(naper)**2 / c
!!$!    p(3) = 0
!!$!    p(4) = 0
!!$!    p(3:) = beta*(pi*raper(naper)**2)/flux0
!!$
!!$    p(1) = 0
!!$    p(2) = b / raper(naper)**2 / c
!!$!    p(3) = 0
!!$
!!$    par_t = c
!!$    dp = 1e-4
!!$!    dp(1:2) = 1e-4
!!$    lambda = 1e-4
!!$!    back = 0
!!$    do i = 1,50
!!$       call nelmin1(arcfuns,p,dp,s,ifault)
!!$       write(*,'(2f15.5,es10.2,f15.5,f10.3)') p(1)*raper(naper)*c,p(2)*raper(naper)**2*c!, &
!!$!            p(3)*c,sum(abs(back-p(3:)))/nstar,lambda
!!$!       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) .and. &
!!$!            sum(abs(back-p(3:)))/nstar < 0.001 ) exit
!!$       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-7) ) exit
!!$!       back = p(3:)
!!$       p1 = p
!!$       lambda = 2*lambda
!!$    end do
!!$    write(*,'(2f15.5)') p(1)*c,p(2)*raper(naper)**2*c
!!$    b = p(2)*raper(naper)**2
!!$    return
!!$    stop
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$!    beta = p(3:)
!!$
!!$    m = naper - nref
!!$    allocate(dd(m),x(m),y(m),z(m))
!!$
!!$    do l = 1, m
!!$       i = l + nref - 1
!!$       x(l) = 1/radius(i) -    1/radius(i+1)
!!$       y(l) = 1/radius(i)**2 - 1/radius(i+1)**2
!!$       z(l) = radius(i+1)**2 - radius(i)**2
!!$    end do
!!$
!!$    do n = 1,nstar
!!$       do l = 1,m
!!$          i = l + nref - 1
!!$          g = beta(n)*z(l) + a*x(l) + b*y(l)
!!$          dd(l) = abs(grows(n,i) - g)
!!$       end do
!!$       call rinit(dd,d,s)
!!$       dgrows(n,:) = d / 0.6745 + approxerr
!!$    end do
!!$    deallocate(dd,x,y,z)
!!$
!!$    nprint = 0
!!$    p1 = p
!!$    lambda = lambda / 10
!!$    do i = 1,50
!!$       !call lmder2(initfun,p,epsilon(p),nprint,info)
!!$       call lmdif2(inidfun,p,epsilon(p),nprint,info)
!!$!       write(*,'(e15.5,3f19.7)') p(1)*raper(naper)*c,p(2)*raper(naper)**2*c, &
!!$!            p(3)/(pi*raper(naper)**2)*c,sum(abs(back-p(3:)))/nstar
!!$!       if( all(p(1:2) > 0) .and. all(abs(p(1:2)-p1(1:2)) < 1e-5) .and. &
!!$!            sum(abs(back-p(3:)))/nstar < 0.001 ) exit
!!$       lambda = 2*lambda
!!$       p1 = p
!!$       back = p(3:)
!!$    end do
!!$
!!$    a = p(1)*raper(naper)*c
!!$    b = p(2)*raper(naper)**2*c
!!$!    beta = p(3:)*flux0/(pi*raper(naper)**2)*c
!!$!    beta = p(3:)/(pi*raper(naper)**2)*c
!!$
!!$!    write(*,'(2f12.5,i3,a)') a,b,info,' rrr'
!!$!    write(*,*) real(p(3)*flux0/pi/raper(naper)**2)
!!$!    write(*,'(5000f10.2)') beta(1:nstar)
!!$!    write(*,*) ifault
!!$!    stop
!!$
!!$  end subroutine arc_inits
!!$
!!$
!!$ function arcfuns(p) result(s)
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(:), allocatable :: beta,x,y,z,sig
!!$    real(dbl) :: s,a,b,g,f2,r,sbeta,q
!!$    integer :: l,i,n,m
!!$
!!$    m = naper - nref
!!$    n = nstar
!!$    allocate(x(m),y(m),z(m),beta(n),sig(n))
!!$
!!$    a = 0
!!$    q = p(1)
!!$    b = p(2)
!!$    !    sbeta = p(3)
!!$    sbeta = 0
!!$!    q = p(4)
!!$!    beta = p(3:)
!!$    sig = dgrows(1:nstar,nhwhm)
!!$
!!$    do l = 1, m
!!$       i = l + nref - 1
!!$       x(l) = 1/radius(i) -    1/radius(i+1)
!!$       y(l) = 1/radius(i)**2 - 1/radius(i+1)**2
!!$!       z(l) = radius(i+1)**2 - radius(i)**2
!!$    end do
!!$
!!$    s = 0
!!$    do n = 1,nstar
!!$       !       do i = 1, m
!!$       do i = nref, naper
!!$          !          g = a*x(i) + b*y(i)
!!$          g = 1 - q - (a/radius(i) + b/radius(i)**2) + sbeta*radius(i)**1
!!$!          write(*,'(6g12.2)') z(i),x(i),y(i),grows(n,i+nref-1) - g,sig(n)
!!$!          s = s + abs(grows(n,i+nref-1) - g) / sig(n)
!!$          s = s + abs(grows(n,i) - g /par_t) / sig(n)
!!$       end do
!!$    end do
!!$
!!$    f2 = 0
!!$    do i = nref+1,naper
!!$       r = (radius(i) + radius(i-1))/2
!!$       f2 = f2 + 2*(a + 3*b/r)/r**3 * (radius(i) - radius(i-1))
!!$    end do
!!$
!!$    s = s / (m*n) + lambda*(max(0.0,-a)*0 + max(0.0,-b) + max(0.0,-q))
!!$
!!$!    s = s / (m*n) + lambda*(max(0.0,-a) + max(0.0,-b)) !+ sum(abs(beta-back)))
!!$
!!$    deallocate(x,y,z,beta,sig)
!!$
!!$ end function arcfuns
!!$
!!$
!!$
!!$  subroutine initfun(mm,np,p,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: mm,np,ldfjac
!!$    integer, intent(inout) :: iflag
!!$    real(dbl), dimension(np), intent(in) :: p
!!$    real(dbl), dimension(mm), intent(out) :: fvec
!!$    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
!!$    real(dbl), dimension(:), allocatable :: beta,x,y,z,s
!!$    real(dbl), dimension(:,:), allocatable :: rs,f,df
!!$    real(dbl) :: a,b,g
!!$    integer :: n,m,i,l
!!$
!!$    if( iflag == 0 ) then
!!$       write(*,'(6g12.3)') p,fvec
!!$
!!$       if( .true. ) then
!!$          write(*,*) ' jac:',real(fjac(1,:))
!!$          write(*,*) ' jac:',real(fjac(2,:))
!!$          write(*,*) ' jac:',real(fjac(3,:))
!!$
!!$       end if
!!$
!!$       return
!!$    end if
!!$
!!$    m = naper - nref
!!$    n = nstar
!!$    allocate(rs(n,m),f(n,m),df(n,m),x(m),y(m),z(m),beta(n),s(n))
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    beta = p(3:)
!!$    s = dgrows(1:nstar,nhwhm)
!!$
!!$    do l = 1, m
!!$       i = l + nref - 1
!!$       x(l) = 1/radius(i) -    1/radius(i+1)
!!$       y(l) = 1/radius(i)**2 - 1/radius(i+1)**2
!!$       z(l) = radius(i+1)**2 - radius(i)**2
!!$    end do
!!$
!!$    do n = 1,nstar
!!$       do l = 1, m
!!$          i = l + nref - 1
!!$          g = beta(n)*z(l) + a*x(l) + b*y(l)
!!$          rs(n,l) = (grows(n,i) - g) / s(n)
!!$          f(n,l) = huber(rs(n,l))
!!$!          write(*,'(i2,4f15.5)') i,rs(n,i),f(n,i),grows(n,i+nref-1) - beta(n)*z(i) - a*x(i) - b*y(i),s(n)
!!$       end do
!!$    end do
!!$!    stop
!!$!    f = rs
!!$
!!$    if( iflag == 1 ) then
!!$
!!$       fvec(1:2) = 0
!!$       do n = 1,nstar
!!$          fvec(1) = fvec(1) + sum(f(n,:)*x) / s(n)
!!$          fvec(2) = fvec(2) + sum(f(n,:)*y) / s(n)
!!$       end do
!!$       forall( n = 1:nstar )
!!$          fvec(2+n) = sum(f(n,:)*z) / s(n) - lambda*(beta(n) - back(n))
!!$       end forall
!!$       fvec(1) = fvec(1) - lambda*isigmoid(a)
!!$       fvec(2) = fvec(2) - lambda*isigmoid(b)
!!$       fvec = - fvec
!!$
!!$    else if( iflag == 2 ) then
!!$
!!$       do n = 1, nstar
!!$          call dhubers(rs(n,:),df(n,:))
!!$       end do
!!$!       df = 1
!!$
!!$       fjac = 0
!!$       do n = 1,nstar
!!$          fjac(1,1) = fjac(1,1) + sum(df(n,:)*x**2) / s(n)**2
!!$          fjac(1,2) = fjac(1,2) + sum(df(n,:)*x*y) / s(n)**2
!!$          fjac(2,2) = fjac(2,2) + sum(df(n,:)*y**2) / s(n)**2
!!$       end do
!!$       forall( n = 1:nstar )
!!$          fjac(1,2+n) = sum(df(n,:)*z*x) / s(n)**2
!!$          fjac(2,2+n) = sum(df(n,:)*z*y) / s(n)**2
!!$          fjac(2+n,2+n) = sum(df(n,:)*z**2) / s(n)**2
!!$       end forall
!!$
!!$       fjac(2,1) = fjac(1,2)
!!$       forall( n = 1:nstar )
!!$          fjac(2+n,1) = fjac(1,2+n)
!!$          fjac(2+n,2) = fjac(2,2+n)
!!$       end forall
!!$
!!$       fjac(1,1) = fjac(1,1) - lambda*nsigmoid(a)
!!$       fjac(2,2) = fjac(2,2) - lambda*nsigmoid(b)
!!$       forall( n = 1:nstar ) fjac(2+n,2+n) = fjac(2+n,2+n) + lambda
!!$
!!$    end if
!!$
!!$    deallocate(rs,f,df,x,y,z,beta,s)
!!$
!!$  end subroutine initfun


!!$  subroutine inidfun(mm,np,p,fvec,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: mm,np
!!$    integer, intent(in out) :: iflag
!!$    real(dbl), dimension(np), intent(in) :: p
!!$    real(dbl), dimension(mm), intent(out) :: fvec
!!$    real(dbl), dimension(:), allocatable :: beta,x,y,z,s
!!$    real(dbl), dimension(:,:), allocatable :: rs,f
!!$    real(dbl) :: a,b,g
!!$    integer :: n,m,i,j,l
!!$
!!$    if( iflag == 0 ) then
!!$       write(*,'(2g13.5)') p
!!$       write(*,'(2g13.5)') fvec
!!$       return
!!$    end if
!!$
!!$    m = naper - nref
!!$    n = nstar
!!$    allocate(rs(n,m),f(n,m),x(m),y(m),z(m),beta(n),s(n))
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    beta = p(3:)
!!$    s = dgrows(1:nstar,nhwhm)
!!$
!!$    do l = 1, m
!!$       i = l + nref - 1
!!$       x(l) = 1/radius(i) -    1/radius(i+1)
!!$       y(l) = 1/radius(i)**2 - 1/radius(i+1)**2
!!$       z(l) = radius(i+1)**2 - radius(i)**2
!!$    end do
!!$
!!$    do n = 1,nstar
!!$       do i = 1, m
!!$          g = beta(n)*z(i) + a*x(i) + b*y(i)
!!$          rs(n,i) = (grows(n,i+nref-1) - g) / s(n)
!!$          f(n,i) = huber(rs(n,i))
!!$!          if( n == 1 ) &
!!$!          write(*,'(i2,8f15.5)') i,rs(n,i),f(n,i),g,grows(n,i+nref-1),s(n)
!!$       end do
!!$    end do
!!$!    f = rs
!!$
!!$    fvec(1:2) = 0
!!$    do n = 1,nstar
!!$       fvec(1) = fvec(1) + sum(f(n,:)*x) / s(n)
!!$       fvec(2) = fvec(2) + sum(f(n,:)*y) / s(n)
!!$    end do
!!$    forall( n = 1:nstar )
!!$       fvec(2+n) = sum(f(n,:)*z) / s(n) - lambda*(beta(n) - back(n))
!!$    end forall
!!$
!!$    fvec(1) = fvec(1) - lambda*isigmoid(a)
!!$    fvec(2) = fvec(2) - lambda*isigmoid(b)
!!$
!!$    fvec = - fvec
!!$
!!$    deallocate(rs,f,x,y,z,beta,s)
!!$
!!$  end subroutine inidfun
!!$
!!$
!!$  subroutine anchor_inits1(flux0,beta,a,b,e)
!!$
!!$    use robustmean
!!$    use neldermead
!!$
!!$    real(dbl), dimension(:), intent(in) :: flux0, beta
!!$    real(dbl), intent(in) :: a,b
!!$    real(dbl), intent(in out) :: e
!!$    real(dbl), dimension(1) :: p,dp
!!$    integer :: ifault,i,n,m
!!$    real(dbl) :: d,s,g
!!$
!!$!    nstar = 2
!!$
!!$    do n = 1,nstar
!!$       do i = nref,naper
!!$          grows(n,i) = 1 - cts(n,i) / flux0(n) - (a/raper(i) + b/raper(i)**2) - &
!!$               beta(n)*aper(i)
!!$!          write(*,*) raper(i),grows(n,i),'LL'
!!$       end do
!!$       call rinit(dcts(n,nref:naper),d,s)
!!$       dgrows(n,:) = d / flux0(n) + approxerr
!!$!       write(*,*) n,dgrows(n,nhwhm)
!!$    end do
!!$!stop
!!$    p(1) = 0
!!$
!!$    dp = 1e-4
!!$    call nelmin1(anchorfuns1,p,dp,s,ifault)
!!$!    write(*,'(4f15.5)') p(1),p(2)*flux0(1),sum(abs(phix-p(2:)))/nstar
!!$
!!$    e = p(1) !*raper(naper)
!!$    write(*,'(f12.5,i3,a)') e,ifault,' anchor_inits1'
!!$!    write(*,*) real(p(3)*flux0/pi/raper(naper)**2)
!!$    !    write(*,'(5000f10.5)') phi(1:nstar)
!!$!    write(*,*) real(phi(1:nstar))
!!$!    write(*,*) real(p(2:))
!!$!    write(*,*) ifault
!!$!    stop
!!$
!!$  end subroutine anchor_inits1
!!$
!!$
!!$ function anchorfuns1(p) result(s)
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(:), allocatable :: sig
!!$    real(dbl) :: s,t,g
!!$    integer :: i,n,m
!!$
!!$    m = naper - nref + 1
!!$    n = nstar
!!$    allocate(sig(n))
!!$
!!$    t = p(1)
!!$    sig = dgrows(1:nstar,nhwhm)
!!$
!!$    s = 0
!!$    do n = 1,nstar
!!$       do i = 1, m
!!$          s = s + abs(grows(n,i+nref-1) - t) / sig(n)
!!$       end do
!!$    end do
!!$    s = s / (m*n)
!!$
!!$    deallocate(sig)
!!$
!!$  end function anchorfuns1
!!$
!!$
!!$  subroutine anchor_inits(flux0,a,b,beta,t,phi)
!!$
!!$    use minpacks
!!$    use robustmean
!!$    use neldermead
!!$
!!$    real(dbl), dimension(:), intent(in) :: flux0,beta
!!$    real(dbl), intent(in) :: a,b
!!$    real(dbl), dimension(:), intent(in out) :: phi
!!$    real(dbl), intent(in out) :: t
!!$    real(dbl), dimension(size(phi)+1) :: p,dp,p1
!!$    real(dbl), dimension(:), allocatable :: dabs
!!$    integer :: ifault,info,nprint,i,n,m
!!$    real(dbl) :: d,s,g
!!$
!!$!    nstar = 2
!!$
!!$    do n = 1,nstar
!!$       do i = nref,naper
!!$          grows(n,i) = 1 - cts(n,i) / flux0(n) + 0*beta(n)*aper(i) & !/flux0(n) &
!!$               - (a/raper(i) + b/raper(i)**2)
!!$!          write(*,*) raper(i),1 - cts(n,i) / flux0(n),(a/raper(i) + b/raper(i)**2)
!!$!          grows(n,i) = sqrt(cts(n,i+1)/flux0(n)) - sqrt(cts(n,i) / flux0(n))
!!$!          write(*,*) i
!!$       end do
!!$       call rinit(dcts(n,nref:naper),d,s)
!!$       dgrows(n,:) = d / flux0(n) + approxerr
!!$!       write(*,*) n,dgrows(n,nhwhm)
!!$    end do
!!$!stop
!!$    p(1) = t
!!$    p(2:) = phi !/ flux0(1:nstar)
!!$!    p = 0*phi / flux0(1:nstar)
!!$
!!$
!!$    dp = 1e-6
!!$    dp(1) = 1e-4
!!$    lambda = 1e-4
!!$    phix = 0
!!$    do i = 1,50
!!$       call nelmin1(anchorfuns,p,dp,s,ifault)
!!$       write(*,'(4f15.5)') p(1),p(2)*flux0(1),sum(abs(phix-p(2:)))/nstar
!!$       if( abs(p(1)-p1(1)) < 1e-5 .and. sum(abs(phix-p(2:)))/nstar < 1e-7 ) exit
!!$       phix = p(2:)
!!$       p1 = p
!!$       lambda = 2*lambda
!!$    end do
!!$!stop
!!$
!!$    t = p(1)
!!$    phi = p(2:)
!!$
!!$    m = naper - nref + 1
!!$    allocate(dabs(m))
!!$    do n = 1,nstar
!!$       do i = 1,m
!!$          g = t + phi(n)
!!$          dabs(i) = abs(grows(n,i+nref-1) - g)
!!$       end do
!!$       call rinit(dabs,d,s)
!!$       dgrows(n,:) = d / 0.6745 + approxerr
!!$    end do
!!$    deallocate(dabs)
!!$
!!$!    stop
!!$    nprint = 0
!!$
!!$    p1 = p
!!$!    lambda = 1e-6
!!$    do i = 1,20
!!$       call lmder2(anchorder,p,epsilon(p),nprint,info)
!!$!       call lmdif2(anchordif,p,epsilon(p),nprint,info)
!!$       write(*,*) real(p(1:2))
!!$!       if( all(abs(p - p1) < 1e-5) ) exit
!!$       if( abs(p(1)-p1(1)) < 1e-7 .and. sum(abs(phix-p(2:)))/nstar < 1e-7 ) exit
!!$       lambda = 2*lambda
!!$       p1 = p
!!$       phix = p(2:)
!!$    end do
!!$
!!$    t = p(1) !*raper(naper)
!!$    phi = p(2:)!*flux0(1:nstar)
!!$
!!$    write(*,'(f12.5,i3,a)') t,info,' anchor_inits'
!!$!    write(*,*) real(p(3)*flux0/pi/raper(naper)**2)
!!$    !    write(*,'(5000f10.5)') phi(1:nstar)
!!$!    write(*,*) real(phi(1:nstar))
!!$!    write(*,*) real(p(2:))
!!$!    write(*,*) ifault
!!$!    stop
!!$
!!$  end subroutine anchor_inits
!!$
!!$
!!$ function anchorfuns(p) result(s)
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(:), allocatable :: phi,sig
!!$    real(dbl) :: s,t,g
!!$    integer :: i,n,m
!!$
!!$    m = naper - nref + 1
!!$    n = nstar
!!$    allocate(phi(n),sig(n))
!!$
!!$    t = p(1)
!!$    phi = p(2:)
!!$    sig = dgrows(1:nstar,nhwhm)
!!$
!!$    s = 0
!!$    do n = 1,nstar
!!$       do i = 1, m
!!$          g = t + phi(n)
!!$          s = s + abs(grows(n,i+nref-1) - g) / sig(n)
!!$       end do
!!$    end do
!!$    s = s / (m*n) + lambda*sum(abs(phi-phix)) + lambda*sum(abs(phi))
!!$
!!$    deallocate(phi,sig)
!!$
!!$  end function anchorfuns
!!$
!!$
!!$  subroutine anchordif(mm,np,p,fvec,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: mm,np
!!$    integer, intent(in out) :: iflag
!!$    real(dbl), dimension(np), intent(in) :: p
!!$    real(dbl), dimension(mm), intent(out) :: fvec
!!$    real(dbl), dimension(:), allocatable :: s,phi
!!$    real(dbl), dimension(:,:), allocatable :: rs,f
!!$!    real(dbl), dimension(2,2) :: dfjac
!!$    real(dbl) :: t, g
!!$    integer :: n,m,i
!!$
!!$    if( iflag == 0 ) then
!!$       write(*,'(2g13.5)') p
!!$       write(*,'(2g13.5)') fvec
!!$       return
!!$    end if
!!$
!!$    m = naper - nref + 1
!!$    n = nstar
!!$    allocate(rs(n,m),f(n,m),phi(n),s(n))
!!$
!!$    t = p(1)
!!$    phi = p(2:)
!!$    s = dgrows(1:nstar,nhwhm)
!!$
!!$    do n = 1,nstar
!!$       do i = 1, m
!!$          g = t + phi(n)
!!$          rs(n,i) = (grows(n,i+nref-1) - g) / s(n)
!!$          f(n,i) = huber(rs(n,i))
!!$!          write(*,'(i2,8f15.5)') i,rs(n,i),f(n,i),g,grows(n,i+nref-1),s(n)
!!$       end do
!!$    end do
!!$!    stop
!!$!    f = rs !*s(1)
!!$
!!$    fvec(1) = 0
!!$    do n = 1,nstar
!!$       fvec(1) = fvec(1) + sum(f(n,:)) / s(n)
!!$!       write(*,*) real(sum(f(n,:)*x))
!!$    end do
!!$!    write(*,*) real(t)
!!$    forall( n = 1:nstar )
!!$!       fvec(1+n) = sum(f(n,:)) / s(n) + lambda*(phi(n)*(phi(n) - 1))
!!$       fvec(1+n) = sum(f(n,:)) / s(n) + lambda*(phi(n) - phix(n))
!!$    end forall
!!$    !    write(*,*) real(sum(f(1,:)*z)) / s(1)
!!$
!!$!    fvec(1) = fvec(1) + lambda*(log(t) + 1)
!!$    fvec = - fvec
!!$
!!$!    write(*,'(6g12.3)') fvec,rs!sqrt(sum(rs**2)/(m-4))
!!$
!!$    deallocate(rs,f,phi,s)
!!$
!!$  end subroutine anchordif
!!$
!!$  subroutine anchorder(mm,np,p,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: mm,np,ldfjac
!!$    integer, intent(inout) :: iflag
!!$    real(dbl), dimension(np), intent(in) :: p
!!$    real(dbl), dimension(mm), intent(out) :: fvec
!!$    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
!!$    real(dbl), dimension(:), allocatable :: phi,s
!!$    real(dbl), dimension(:,:), allocatable :: rs,f,df
!!$!    real(dbl), dimension(2,2) :: dfjac
!!$    real(dbl) :: t,g
!!$    integer :: n,m,i,j
!!$
!!$    if( iflag == 0 ) then
!!$       write(*,'(6g12.3)') p,fvec
!!$
!!$       if( .true. ) then
!!$          write(*,*) ' jac:',real(fjac(1,:))
!!$          write(*,*) ' jac:',real(fjac(2,:))
!!$          write(*,*) ' jac:',real(fjac(3,:))
!!$
!!$!          call difjac(p(1),p(2),dfjac)
!!$!          write(*,*) 'djac:',dfjac(1,:)
!!$!          write(*,*) 'djac:',dfjac(2,:)
!!$
!!$       end if
!!$
!!$       return
!!$    end if
!!$
!!$    m = naper - nref + 1
!!$    n = nstar
!!$    allocate(rs(n,m),f(n,m),df(n,m),phi(n),s(n))
!!$
!!$    t = p(1)
!!$    phi = p(2:)
!!$    s = dgrows(1:nstar,nhwhm)
!!$
!!$    do n = 1,nstar
!!$       do i = 1, m
!!$          g = t + phi(n)
!!$          rs(n,i) = (grows(n,i+nref-1) - g) / s(n)
!!$          f(n,i) = huber(rs(n,i))
!!$!          write(*,'(i2,4f15.5)') i,rs(n,i),f(n,i),grows(n,i+nref-1) - beta(n)*z(i) - a*x(i) - b*y(i),s(n)
!!$       end do
!!$    end do
!!$!    stop
!!$!    f = rs
!!$
!!$    if( iflag == 1 ) then
!!$
!!$       fvec(1) = 0
!!$       do n = 1,nstar
!!$          fvec(1) = fvec(1) + sum(f(n,:)) / s(n)
!!$       end do
!!$       forall( n = 1:nstar )
!!$!          fvec(1+n) = sum(f(n,:)) / s(n) + lambda*(phi(n)*(phi(n) - 1))
!!$          fvec(1+n) = sum(f(n,:)) / s(n) + lambda*(phi(n) - phix(n))
!!$       end forall
!!$!       forall( n = 1:nstar ) fvec(1+n) = sum(f(n,:)) / s(n)
!!$!       fvec(1) = fvec(1) + lambda*max(0.0,-a)
!!$!       fvec(2) = fvec(2) + lambda*max(0.0,-b)
!!$       fvec = - fvec
!!$
!!$    else if( iflag == 2 ) then
!!$
!!$       do n = 1, nstar
!!$          call dhubers(rs(n,:),df(n,:))
!!$       end do
!!$!       df = 1
!!$
!!$       fjac = 0
!!$       do n = 1,nstar
!!$          fjac(1,1) = fjac(1,1) + sum(df(n,:)) / s(n)**2
!!$       end do
!!$       forall( n = 1:nstar )
!!$          fjac(1,1+n) = sum(df(n,:)) / s(n)**2 + lambda
!!$          fjac(1+n,1+n) = sum(df(n,:)) / s(n)**2
!!$       end forall
!!$
!!$       forall( n = 1:nstar )
!!$          fjac(1+n,1) = fjac(1,1+n)
!!$       end forall
!!$
!!$    end if
!!$
!!$!    write(*,'(6g12.3)') p,fvec
!!$
!!$
!!$    deallocate(rs,f,df,phi,s)
!!$
!!$  end subroutine anchorder
!!$
!!$
!!$
!!$
!!$  !-----------------------------------------------------------------------
!!$
!!$
!!$
!!$
!!$  subroutine update_backflux(a,cts,dcts,flux,back)
!!$
!!$    use NelderMead
!!$    use minpacks
!!$
!!$    real(dbl), intent(in) :: a
!!$    real(dbl), dimension(:), intent(in) :: cts, dcts
!!$    real(dbl), intent(in out) :: flux, back
!!$
!!$    real(dbl), dimension(2) :: p,dp
!!$    real(dbl) :: pmin
!!$    integer :: ifault, nprint, info
!!$
!!$    p = 0
!!$    dp = 1e-6
!!$    cts_(nref:naper) = cts(nref:naper) / flux - (1 - a/raper(nref:naper))
!!$!    cts_(nref:naper) = sqrt(cts(nref:naper) / flux) - (1 - a/raper(nref:naper)/2)
!!$    dcts_(nref:naper) = dcts(nref:naper) / flux
!!$
!!$    call nelmin1(backflux,p,dp,pmin,ifault)
!!$
!!$
!!$    nprint = 0
!!$    sig = 1
!!$    call lmder2(fackblux,p,epsilon(p),nprint,info)
!!$
!!$    back = flux*p(2)/pi
!!$    flux = flux*(1 + p(1))
!!$!    back = 2*flux*p(2)/pi
!!$!    flux = flux*(1 + p(1)/2)
!!$
!!$!    write(*,*) real(flux),real(back),'ppp'
!!$!    write(*,*) real(p)
!!$
!!$  end subroutine update_backflux
!!$
!!$
!!$  function backflux(p) result(s)
!!$
!!$    use robustmean
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(size(raper)) :: ares,stat
!!$    real(dbl) :: s,a,b,q,r,f
!!$    integer :: i
!!$
!!$    f = p(1)
!!$    b = p(2)
!!$
!!$    s = 0
!!$    do i = nref, naper
!!$       r = raper(i)
!!$       q = f + b*r**2
!!$       s = s + abs(cts_(i) - q) / dcts_(i)
!!$    end do
!!$
!!$  end function
!!$
!!$
!!$  subroutine fackblux(m,np,p,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: m,np,ldfjac
!!$    integer, intent(inout) :: iflag
!!$    real(dbl), dimension(np), intent(in) :: p
!!$    real(dbl), dimension(m), intent(out) :: fvec
!!$    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
!!$    real(dbl), dimension(:), allocatable :: rs,f,df,r
!!$    real(dbl), dimension(2,2) :: dfjac
!!$    integer :: n,i
!!$
!!$    if( iflag == 0 ) then
!!$       write(*,'(6g12.3)') p,fvec
!!$
!!$       if( .true. ) then
!!$          write(*,*) ' jac:',real(fjac(1,:))
!!$          write(*,*) ' jac:',real(fjac(2,:))
!!$          write(*,*) ' jac:',real(fjac(3,:))
!!$
!!$!          call difjac(p(1),p(2),dfjac)
!!$!          write(*,*) 'djac:',dfjac(1,:)
!!$!          write(*,*) 'djac:',dfjac(2,:)
!!$
!!$       end if
!!$
!!$       return
!!$    end if
!!$
!!$    n = naper - nref + 1
!!$    allocate(rs(n),f(n),df(n),r(n))
!!$
!!$    do i = 1,n
!!$       r(i) = raper(i+nref-1)
!!$       rs(i) = (p(1) + p(2)*r(i)**2 - cts_(i+nref-1)) / sig
!!$       f(i) = huber(rs(i))
!!$    end do
!!$!    f = r
!!$
!!$    if( iflag == 1 ) then
!!$
!!$       fvec(1) = sum(f)
!!$       fvec(2) = sum(f*r**2)
!!$       fvec = - fvec / sig
!!$
!!$    else if( iflag == 2 ) then
!!$
!!$       call dhubers(rs,df)
!!$!       df = 1
!!$
!!$       fjac(1,1) = sum(df)
!!$       fjac(1,2) = sum(df*r**2)
!!$       fjac(2,2) = sum(df*r**4)
!!$
!!$       fjac(2,1) = fjac(1,2)
!!$       fjac = fjac / sig**2
!!$
!!$    end if
!!$
!!$    deallocate(rs,f,df,r)
!!$
!!$  end subroutine fackblux
!!$
!!$
!!$
!!$
!!$
!!$
!!$  subroutine asymptote_init5(a,q,beta,da,dq,dbeta)
!!$
!!$    use NelderMead
!!$    use robustmean
!!$    use minpacks
!!$    use rfun
!!$
!!$
!!$    real(dbl), intent(out) :: a,q,beta,da,dq,dbeta
!!$    real(dbl), dimension(3) :: p,dp,p1,dp1
!!$    real(dbl), dimension(3,3) :: fjac,hess
!!$    real(dbl) :: pmin,d,r,fmin,fmax,mad,s,sum2,sum3,rs
!!$    integer :: i,j,ifault,n,info
!!$
!!$    n = naper - nref + 1
!!$    p(1) = 0
!!$    p(2) = 1
!!$    p(3) = 0
!!$    dp(1) = 0.001
!!$    dp(2) = 0.0001
!!$    dp(3) = 1e-4
!!$    n = naper - 1
!!$!    r = raper(naper)
!!$    fmin = grow_(nref)
!!$    fmax = grow_(naper)
!!$    fmin = minval(grow_(nref:naper))
!!$    fmax = maxval(grow_(nref:naper))
!!$    d = fmax - fmin
!!$    grow_ = (grow_ - fmin)  / d
!!$!    write(*,*) fmin,fmax,d
!!$!    stop
!!$
!!$    call rinit(dgrow_(nref:naper)/d,sig,s)
!!$    write(*,*) 'sig:',real(sig)
!!$    lambda = 0
!!$    do i = 1,100
!!$       call nelmin1(fainit5,p,dp,pmin,ifault)
!!$       lambda = lambda + 0.01
!!$       write(*,'(4f13.5,2g12.3,i3)') p,lambda,par_t,par_e,ifault
!!$       if( p(1) > 0 ) then
!!$!             p = p1
!!$!             dp = dp1
!!$          exit
!!$       end if
!!$       p1 = p
!!$       dp1 = dp
!!$    end do
!!$
!!$    i = 0
!!$    call lmder2(funder5,p,epsilon(p),i,info)
!!$    write(*,*) 'lmder:',real(p),info
!!$
!!$    ! osetrit pripad info /= 2
!!$
!!$    i = 2
!!$    call funder5(3,3,p,dp,fjac,3,i)
!!$    call qrinv(fjac,hess)
!!$    sum2 = 0
!!$    sum3 = 0
!!$    do i = 1,naper - nref + 1
!!$       r = radius(i+nref-1)
!!$       rs = (p(2) - p(1)/r + p(3)*r**2 - grow_(i+nref-1)) / sig
!!$       sum2 = sum2 + dhuber(rs)
!!$       sum3 = sum3 + huber(rs)**2
!!$    end do
!!$    s = sig*sqrt(sum3/sum2)
!!$    dp = s
!!$!    do i = 1,3
!!$!       write(*,*) (real(fjac(i,j)),j=1,3)
!!$!    end do
!!$!    do i = 1,3
!!$!       write(*,*) (real(hess(i,j)),j=1,3)
!!$!    end do
!!$!    write(*,*) real(fjac(1,1)),real(fjac(2,2)),real(fjac(3,3))
!!$!    write(*,*) real(hess(1,1)),real(hess(2,2)),real(hess(3,3))
!!$    do i = 1, 3
!!$       if( hess(i,i) > 0 .and. hess(i,i) < 1e3) then ! problem s 10e14
!!$          dp(i) = s *sqrt(hess(i,i))
!!$       end if
!!$    end do
!!$
!!$    if( ifault /= 0 ) write(*,*)  'asymptote_init5:',real(p),real(pmin),ifault
!!$
!!$    r = raper(naper)
!!$    a = d*p(1)*r
!!$    q = fmin + d*p(2)
!!$    beta = d*p(3)/r**2
!!$    da = a*dp(1)/p(1)
!!$    dq = d*p(2)
!!$    dbeta = beta*dp(3)/p(3)
!!$
!!$    write(*,'(5f15.5)') a,q,beta,pmin
!!$!    stop
!!$  end subroutine asymptote_init5
!!$
!!$
!!$  function fainit5(p)
!!$
!!$    use robustmean
!!$
!!$    real(dbl) :: fainit5
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(size(raper)) :: ares,stat
!!$    real(dbl) :: s,a,q,r,x,x1,f,r1,f1,f2,h,beta
!!$    integer :: i,j,n
!!$
!!$    a = p(1)
!!$    q = p(2)
!!$    beta = p(3)
!!$
!!$    h = 0
!!$    s = 0
!!$    f1 = 0
!!$    do i = nref-1, naper
!!$       r = radius(i)
!!$       f = q - a/r + beta*r**2
!!$       f1 = a/r**2 + 2*beta*r
!!$       f2 = -2*a/r**3 + 2*beta
!!$
!!$       ! dgrow_ are not used because initial estimates have no Normal
!!$       ! distribution, errors includes systematic offsets which can
!!$       ! confuse fitting procedure
!!$       s = s + abs(grow_(i) - f)
!!$       h = h + f2 / (1 + f1**2)**1.5
!!$!       write(*,*) real(r),real((grow_(i) - f))
!!$    end do
!!$!    s = s !/ sig
!!$!    r = radius(naper)
!!$    fainit5 = (1-lambda) * s / sig + lambda*(abs(h)*0 + max(0.0,-a))
!!$    par_t = s / sig
!!$    par_e = a
!!$
!!$  end function fainit5
!!$
!!$
!!$  subroutine funder5(m,np,p,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: m,np,ldfjac
!!$    integer, intent(inout) :: iflag
!!$    real(dbl), dimension(np), intent(in) :: p
!!$    real(dbl), dimension(m), intent(out) :: fvec
!!$    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
!!$    real(dbl), dimension(:), allocatable :: rs,f,df,r
!!$    real(dbl), dimension(3,3) :: dfjac
!!$    integer :: n,i,j
!!$    real(dbl) :: q,a,beta
!!$
!!$    if( iflag == 0 ) then
!!$       write(*,'(6g12.3)') p,fvec
!!$
!!$       if( .true. ) then
!!$          write(*,*) ' jac:',real(fjac(1,:))
!!$          write(*,*) ' jac:',real(fjac(2,:))
!!$          write(*,*) ' jac:',real(fjac(3,:))
!!$
!!$!          call difjac(p(1),p(2),dfjac)
!!$!          write(*,*) 'djac:',dfjac(1,:)
!!$!          write(*,*) 'djac:',dfjac(2,:)
!!$
!!$       end if
!!$
!!$       return
!!$    end if
!!$
!!$    n = naper - nref + 1
!!$    allocate(rs(n),f(n),df(n),r(n))
!!$    q = p(2)
!!$    a = p(1)
!!$    beta = p(3)
!!$
!!$    do i = 1,n
!!$       r(i) = radius(i+nref-1)
!!$       rs(i) = (q - a/r(i) + beta*r(i)**2 - grow_(i+nref-1)) / sig
!!$       f(i) = huber(rs(i))
!!$    end do
!!$!    f = r
!!$
!!$    if( iflag == 1 ) then
!!$
!!$       fvec(1) = sum(f)
!!$       fvec(2) = -sum(f/r)
!!$       fvec(3) = sum(f*r**2)
!!$       fvec = fvec / sig
!!$
!!$    else if( iflag == 2 ) then
!!$
!!$       call dhubers(rs,df)
!!$!       df = 1
!!$
!!$       fjac(1,1) = sum(df)
!!$       fjac(1,2) = -sum(df/r)
!!$       fjac(1,3) = sum(df*r**2)
!!$       fjac(2,2) = sum(df/r**2)
!!$       fjac(2,3) = -sum(df*r)
!!$       fjac(3,3) = sum(df*r**4)
!!$
!!$       fjac(2,1) = fjac(1,2)
!!$       fjac(3,1) = fjac(1,3)
!!$       fjac(3,2) = fjac(2,3)
!!$       fjac = fjac / sig**2
!!$
!!$    end if
!!$
!!$    deallocate(rs,f,df,r)
!!$
!!$  end subroutine funder5
!!$
!!$
!!$
!!$
!!$  subroutine asymptote_init3(a,q)
!!$
!!$    use NelderMead
!!$    use robustmean
!!$    use minpacks
!!$
!!$
!!$    real(dbl), intent(out) :: a,q
!!$    real(dbl), dimension(2) :: p,dp,p1,dp1
!!$    real(dbl) :: pmin,d,r,fmin,fmax,mad,s
!!$    integer :: i,ifault,n,info
!!$
!!$    n = naper - nref + 1
!!$    p(1) = 0
!!$    p(2) = 1
!!$    dp(1) = 0.001
!!$    dp(2) = 0.0001
!!$    n = naper - 1
!!$    r = raper(naper)
!!$    fmin = grow_(nref)
!!$    fmax = grow_(naper)
!!$    fmin = minval(grow_(nref:naper))
!!$    fmax = maxval(grow_(nref:naper))
!!$    d = fmax - fmin
!!$    grow_ = (grow_ - fmin)  / d
!!$    write(*,*) fmin,fmax,d
!!$!    stop
!!$
!!$    call rinit(dgrow_(nref:naper)/d,sig,s)
!!$    write(*,*) 'sig:',real(sig)
!!$    lambda = 0
!!$    do i = 1, 6
!!$       call nelmin1(fainit32,p,dp,pmin,ifault)
!!$       lambda = lambda + 0.001
!!$       write(*,'(4f13.5,2g12.3,i3)') p,lambda,pmin,par_t,ifault
!!$       if( .not. (par_t < 1000) ) then
!!$!             p = p1
!!$!             dp = dp1
!!$          exit
!!$       end if
!!$       p1 = p
!!$       dp1 = dp
!!$    end do
!!$
!!$    i = 0
!!$    call lmder2(funder3,p,epsilon(p),i,info)
!!$    write(*,*) 'lmder:',real(p),info
!!$
!!$    if( ifault /= 0 ) write(*,*)  'asymptote_init3:',real(p),real(pmin),ifault
!!$
!!$    p = p*d
!!$    p(2) = fmin + p(2)
!!$    write(*,'(5f15.5)') p
!!$
!!$    r = raper(naper)
!!$    a = p(1)*r
!!$    q = p(2)
!!$
!!$    write(*,'(5f15.5)') a,q,pmin,9999.0
!!$!    stop
!!$  end subroutine asymptote_init3
!!$
!!$  function fainit32(p)
!!$
!!$    use robustmean
!!$
!!$    real(dbl) :: fainit32
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(size(raper)) :: ares,stat
!!$    real(dbl) :: s,a,q,r,x,x1,f,r1,f1,f2,h
!!$    integer :: i,j,n
!!$
!!$    a = p(1)
!!$    q = p(2)
!!$
!!$    h = 0
!!$    s = 0
!!$    f1 = 0
!!$    do i = nref, naper
!!$       r = radius(i)
!!$       f = q - a/r
!!$       f1 = a/r**2
!!$       f2 = -2*a/r**3
!!$!       x = 1/r
!!$!       f =  q - (a + b*x)*x
!!$!       r1 = raper(i+1)
!!$!       f = q - (a + (b + c/r)/r)/r
!!$!       f1 = q - (a + b/r1)/r1
!!$!       x = (grow_(i)/grow_(i+1) - f/f1) !/ (dgrow_(i) / 1.4e-3) !/ raper(i)**2)!*r**2!* x**2
!!$       ! dgrow_ are not used because initial estimates have no Normal
!!$       ! distribution, errors includes systematic offsets which can
!!$       ! confuse fitting procedure
!!$       s = s + abs(grow_(i) - f) ! / 0.07
!!$!       s = s + abs(
!!$       !       if( i > nref ) h = h + abs(abs(x) / abs(x1))
!!$!       h = h + max(0.0,-x)
!!$!       h = h + abs(x)*exp(-abs(x))
!!$       !       f1 = f1 + exp(-abs(x))
!!$!       h = h + f2 / (1 + f1**2)**1.5
!!$!       write(*,*) real(r),real(2*(a + 3*b/r)/r**3),real(1+((a+2*b)/r**2)**2)**1.5
!!$!       write(*,*) i,(1 - (a + (b + c/r)/r)/r) / (1 - (a + (b + c/r1)/r1)/r1)
!!$       write(*,'(5g15.5)') r,grow_(i),f,f1,f2
!!$!       write(*,'(5f15.5)') 1/r,grow_(i),a + 3*b/r + 6*c/r**2
!!$    end do
!!$    s = s / sig
!!$!stop
!!$    r = raper(naper)
!!$    r = radius(naper)
!!$ !   fainit4 = s !+ lambda*(abs(a/r) + abs(b/r**2) + abs(c/r**3)) !(?) + abs(q-1)
!!$ !   fainit4 = s + lambda*abs(1 - (a + (b + c/r)/r)/r + q - 1)
!!$    !    fainit4 = s + lambda*(abs((a + (b + c/r)/r)/r - q) + abs(a) + abs(b) + abs(c))
!!$!    x = 1/raper(naper)
!!$!    fainit4 = s + lambda*(max(0.0,-e) + 0*abs((a + (b + c/r)/r)/r - e) - h)
!!$    fainit32 = s !+ lambda*(max(0.0,-e) + 0*abs((a + b/r)/r - e) - h)
!!$
!!$    !    fainit42 = s + lambda*max(0.0,-h)
!!$    n = (naper - nref + 1)
!!$    h = h / (radius(naper) - radius(nref)) ! rectange integration
!!$    !    fainit42 = s / n + lambda*max(0.0,-h) / n
!!$    f1 = a/r**2
!!$    f2 = - 2*a/r**3
!!$    h = f2 / (1 + f1**2)**1.5
!!$    fainit32 = (1-lambda) * s + lambda*abs(min(0.0,h)) !?-h?
!!$!    fainit42 = (1-lambda) * s + lambda*max(0.0,-h) !?-h?
!!$ !   par_t = max(0.0,-h) / n
!!$    par_t = h
!!$!    write(*,*) real(s),real(h),real(abs(min(0.0,h))),real(par_t)
!!$!    write(*,*) real(p),real(s)
!!$
!!$  end function fainit32
!!$
!!$
!!$
!!$
!!$  subroutine asymptote_init2(a,b,q)
!!$
!!$    use NelderMead
!!$    use robustmean
!!$    use minpacks
!!$
!!$
!!$    real(dbl), intent(out) :: a,b
!!$    real(dbl), intent(out), optional :: q
!!$    real(dbl), dimension(3) :: p,dp,p1,dp1
!!$    real(dbl), dimension(:), allocatable :: res
!!$    real(dbl) :: pmin,d,r,fmin,fmax,mad,s
!!$    integer :: i,ifault,n,info
!!$
!!$    n = naper - nref + 1
!!$!    p(1) = sum((1-grow_(nref:naper))/radius(nref:naper))/sum(1/radius(nref:naper)**2)
!!$!    p(1) = sum((1-grow_(nref:naper))/raper(nref:naper))/sum(1/raper(nref:naper)**2)
!!$!    p(1) = p(1) / 10
!!$!    r = radius(naper)
!!$!        write(*,*) n,p(1),p(1)/r
!!$    p(1) = 0
!!$    p(2) = 0
!!$!    p(3) = 0
!!$    dp(1) = 0.001
!!$    dp(2) = 0.001
!!$!    dp(3) = 0.001
!!$    n = naper - 1
!!$    r = raper(naper)
!!$    fmin = grow_(nref)
!!$    fmax = grow_(naper)
!!$    fmin = minval(grow_(nref:naper))
!!$    fmax = maxval(grow_(nref:naper))
!!$    d = fmax - fmin
!!$    grow_ = (grow_ - fmin)  / d
!!$    write(*,*) fmin,fmax,d
!!$!    stop
!!$
!!$    if( present(q) ) then
!!$       p(3) = 1 !-p(1)/r  !max(-p(1)/r,epsilon(p))
!!$       dp(3) = 0.0001
!!$!       p(4) = 1
!!$       !       dp(4) = 0.1*p(4)
!!$!       write(*,*) real(dgrow_(nref:naper)/d)
!!$       call rinit(dgrow_(nref:naper)/d,sig,s)
!!$       write(*,*) 'sig:',real(sig)
!!$!       stop
!!$       lambda = 0
!!$!       call nelmin1(fainit42,p,dp,pmin,ifault)
!!$!       if( par_t > 0 ) then
!!$!          lambda = 0.1*par_t
!!$       do i = 1, 6
!!$!          dp = 1e-4
!!$!          dp(3) = 1e-4
!!$          call nelmin1(fainit42,p,dp,pmin,ifault)
!!$!       end if
!!$          lambda = lambda + 0.001
!!$          write(*,'(4f13.5,2g10.3,i3)') p,lambda,pmin,par_t,ifault
!!$!          write(*,*) real(p(3) - (p(1) + p(2)/radius(nref:naper))/radius(nref:naper) - grow_(nref:naper))
!!$!          if( .not. (par_t < 0.6745) ) exit
!!$          if( .not. (par_t < 1000) ) then
!!$!             p = p1
!!$!             dp = dp1
!!$             exit
!!$          end if
!!$          p1 = p
!!$          dp1 = dp
!!$       end do
!!$
!!$       ! asi zbytecne ..
!!$       goto 33
!!$       allocate(res(naper))
!!$       do i = nref,naper
!!$          r = radius(i)
!!$          res(i) = p(3) - (p(1) + p(2)/r)/r - grow_(i)
!!$          write(*,*) real(r),real(res(i))
!!$       end do
!!$       call rinit(abs(res(nref:naper))/sig,mad,s)
!!$       sig = mad / 0.6745
!!$       write(*,*) 'mad=',real(mad),real(sig),s
!!$       deallocate(res)
!!$       ! ... az sem
!!$       33 continue
!!$       i = 0
!!$!       sig = 0.05
!!$       call lmder2(funder,p,epsilon(p),i,info)
!!$       write(*,*) 'lmder:',real(p),info
!!$
!!$
!!$!       par_t = 0.99
!!$!       lambda = 1
!!$!       do i = 1, 10
!!$!          call nelmin1(fainit4a,p(1:3),dp(1:3),pmin,ifault,kcount0=1000*nstar,reqmin0=1d-10)
!!$!          lambda = 2*lambda
!!$!          call rmean(cts(1:3,n)/cts(1:3,naper),par_t,d)
!!$!          par_t = 1 - p(1)/r - p(2)/r**2 - p(3)/r**3
!!$!          write(*,*) real(p(1:3)),real(par_t)
!!$!       end do
!!$    else
!!$       call nelmin1(fainit3,p(1:3),dp(1:3),pmin,ifault)
!!$    end if
!!$
!!$    if( ifault /= 0 ) write(*,*)  'asymptote_init:',real(p),real(pmin),ifault
!!$
!!$    p = p*d
!!$    p(3) = fmin + p(3)
!!$    write(*,'(5f15.5)') p
!!$
!!$    r = raper(naper)
!!$    a = p(1)*r
!!$    b = p(2)*r**2
!!$!    c = p(3)*r**3
!!$
!!$    if( present(q) ) q = p(3)
!!$
!!$    write(*,'(5f15.5)') a,b,q,pmin
!!$!    stop
!!$  end subroutine asymptote_init2
!!$
!!$  function fainit42(p)
!!$
!!$    use robustmean
!!$
!!$    real(dbl) :: fainit42
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(size(raper)) :: ares,stat
!!$    real(dbl) :: s,a,b,q,r,x,x1,f,r1,f1,f2,h
!!$    integer :: i,j,n
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    q = p(3)
!!$!    lambda = p(4)
!!$
!!$!    if( q <= 0 ) then
!!$!       fainit4 = 1e5
!!$!       return
!!$!    end if
!!$
!!$    h = 0
!!$    s = 0
!!$    f1 = 0
!!$    do i = nref, naper
!!$       r = radius(i)
!!$       f = q - (a + b/r)/r
!!$!       x = 1/r
!!$!       f =  q - (a + b*x)*x
!!$!       r1 = raper(i+1)
!!$!       f = q - (a + (b + c/r)/r)/r
!!$!       f1 = q - (a + b/r1)/r1
!!$!       x = (grow_(i)/grow_(i+1) - f/f1) !/ (dgrow_(i) / 1.4e-3) !/ raper(i)**2)!*r**2!* x**2
!!$       ! dgrow_ are not used because initial estimates have no Normal
!!$       ! distribution, errors includes systematic offsets which can
!!$       ! confuse fitting procedure
!!$       s = s + abs(grow_(i) - f) ! / 0.07
!!$!       s = s + abs(
!!$       !       if( i > nref ) h = h + abs(abs(x) / abs(x1))
!!$       ares(i) = abs(x)
!!$       x1 = x
!!$!       h = h + max(0.0,-x)
!!$!       h = h + abs(x)*exp(-abs(x))
!!$       !       f1 = f1 + exp(-abs(x))
!!$       f1 = (a + 2*b)/r**2
!!$       f2 = - 2*(a + 3*b/r)/r**3
!!$       h = h + f2 / (1 + f1**2)**1.5
!!$!       write(*,*) real(r),real(2*(a + 3*b/r)/r**3),real(1+((a+2*b)/r**2)**2)**1.5
!!$!       write(*,*) i,(1 - (a + (b + c/r)/r)/r) / (1 - (a + (b + c/r1)/r1)/r1)
!!$!       write(*,'(5g15.5)') r,grow_(i),grow_(i) - f
!!$!       write(*,'(5f15.5)') 1/r,grow_(i),a + 3*b/r + 6*c/r**2
!!$    end do
!!$    s = s / sig
!!$!stop
!!$    r = raper(naper)
!!$    r = radius(naper)
!!$ !   fainit4 = s !+ lambda*(abs(a/r) + abs(b/r**2) + abs(c/r**3)) !(?) + abs(q-1)
!!$ !   fainit4 = s + lambda*abs(1 - (a + (b + c/r)/r)/r + q - 1)
!!$    !    fainit4 = s + lambda*(abs((a + (b + c/r)/r)/r - q) + abs(a) + abs(b) + abs(c))
!!$!    x = 1/raper(naper)
!!$!    fainit4 = s + lambda*(max(0.0,-e) + 0*abs((a + (b + c/r)/r)/r - e) - h)
!!$    fainit42 = s !+ lambda*(max(0.0,-e) + 0*abs((a + b/r)/r - e) - h)
!!$
!!$    !    fainit42 = s + lambda*max(0.0,-h)
!!$    n = (naper - nref + 1)
!!$    h = h / (radius(naper) - radius(nref)) ! rectange integration
!!$    !    fainit42 = s / n + lambda*max(0.0,-h) / n
!!$    h = f2 / (1 + f1**2)**1.5
!!$    fainit42 = (1-lambda) * s + lambda*abs(min(0.0,h)) !?-h?
!!$!    fainit42 = (1-lambda) * s + lambda*max(0.0,-h) !?-h?
!!$ !   par_t = max(0.0,-h) / n
!!$    par_t = h
!!$!    write(*,*) real(s),real(h),real(abs(min(0.0,h))),real(par_t)
!!$    write(*,*) real(p),real(s)
!!$!    fainit42 = s + lambda*(sum(abs(p(1:2))) + abs(p(3)-1))
!!$!    h = maxval(ares(nref:naper))
!!$!    par_t = h / (naper - nref + 1)
!!$
!!$    do i = nref,naper
!!$       stat(i) = 0
!!$       do j = nref,naper
!!$          if( i /= j ) then
!!$             !             write(*,*) res(i),res(j)
!!$             if( ares(i) > epsilon(ares) ) then
!!$                stat(i) = stat(i) + ares(j)/ares(i)
!!$             else
!!$                stat(i) = 1e3
!!$             end if
!!$          end if
!!$       end do
!!$       stat(i) = stat(i) / (naper - nref)
!!$    end do
!!$
!!$
!!$!    call rinit(ares(nref:naper),x1,f1)
!!$!    x1 = ares(nref)
!!$    !    x1 = minval(ares(nref:naper))
!!$!    x1 = 1e-3
!!$!    h = sum(ares(nref:naper)/x1*exp(-ares(nref:naper)/x1)) / (2*x1) / (naper - nref + 1)
!!$!    h = sum(ares(nref:naper)*exp(-ares(nref:naper)))
!!$!    write(*,*) s,h
!!$!    write(*,*) s,log(maxval(stat(nref:naper)) / minval(stat(nref:naper)))
!!$!    fainit42 = s + lambda*(log(maxval(stat(nref:naper)) / minval(stat(nref:naper))))
!!$!    fainit42 = s + lambda*abs(h-0.6745)
!!$!    fainit42 = s + lambda*h
!!$
!!$    ! pridat podminku na zhruba stejne odchylky |x| ... ?
!!$    ! aby se zamezilo problemum s medianem pro malo hodnot, nebo primo odhad medianu?
!!$
!!$!    write(*,*) h/f1
!!$
!!$  end function fainit42
!!$
!!$
!!$  subroutine funder(m,np,p,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: m,np,ldfjac
!!$    integer, intent(inout) :: iflag
!!$    real(dbl), dimension(np), intent(in) :: p
!!$    real(dbl), dimension(m), intent(out) :: fvec
!!$    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
!!$    real(dbl), dimension(:), allocatable :: rs,f,df,r
!!$    real(dbl), dimension(3,3) :: dfjac
!!$    integer :: n,i
!!$    real(dbl) :: q,a,b
!!$
!!$    if( iflag == 0 ) then
!!$       write(*,'(6g12.3)') p,fvec
!!$
!!$       if( .true. ) then
!!$          write(*,*) ' jac:',real(fjac(1,:))
!!$          write(*,*) ' jac:',real(fjac(2,:))
!!$          write(*,*) ' jac:',real(fjac(3,:))
!!$
!!$!          call difjac(p(1),p(2),dfjac)
!!$!          write(*,*) 'djac:',dfjac(1,:)
!!$!          write(*,*) 'djac:',dfjac(2,:)
!!$
!!$       end if
!!$
!!$       return
!!$    end if
!!$
!!$    n = naper - nref + 1
!!$    allocate(rs(n),f(n),df(n),r(n))
!!$    q = p(3)
!!$    a = p(1)
!!$    b = p(2)
!!$
!!$    do i = 1,n
!!$       r(i) = radius(i+nref-1)
!!$       rs(i) = (q - (a + b/r(i))/r(i) - grow_(i+nref-1)) / sig
!!$       f(i) = huber(rs(i))
!!$    end do
!!$    f = r
!!$
!!$    if( iflag == 1 ) then
!!$
!!$       fvec(1) = sum(f)
!!$       fvec(2) = -sum(f/r)
!!$       fvec(3) = -sum(f/r**2)
!!$       fvec = fvec / sig
!!$
!!$    else if( iflag == 2 ) then
!!$
!!$       do i = 1,n
!!$          df(i) = dhuber(rs(i))
!!$       end do
!!$       df = 1
!!$
!!$       fjac(1,1) = sum(df)
!!$       fjac(1,2) = -sum(df/r)
!!$       fjac(1,3) = -sum(df/r**2)
!!$       fjac(2,2) = sum(df/r**2)
!!$       fjac(2,3) = sum(df/r**3)
!!$       fjac(3,3) = sum(df/r**4)
!!$
!!$       fjac(2,1) = fjac(1,2)
!!$       fjac(3,1) = fjac(1,3)
!!$       fjac(3,2) = fjac(2,3)
!!$       fjac = fjac / sig**2
!!$
!!$    end if
!!$
!!$    deallocate(rs,f,df,r)
!!$
!!$  end subroutine funder
!!$
!!$
!!$  subroutine funder3(m,np,p,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: m,np,ldfjac
!!$    integer, intent(inout) :: iflag
!!$    real(dbl), dimension(np), intent(in) :: p
!!$    real(dbl), dimension(m), intent(out) :: fvec
!!$    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
!!$    real(dbl), dimension(:), allocatable :: rs,f,df,r
!!$    real(dbl), dimension(3,3) :: dfjac
!!$    integer :: n,i
!!$    real(dbl) :: q,a,b
!!$
!!$    if( iflag == 0 ) then
!!$       write(*,'(6g12.3)') p,fvec
!!$
!!$       if( .true. ) then
!!$          write(*,*) ' jac:',real(fjac(1,:))
!!$          write(*,*) ' jac:',real(fjac(2,:))
!!$
!!$!          call difjac(p(1),p(2),dfjac)
!!$!          write(*,*) 'djac:',dfjac(1,:)
!!$!          write(*,*) 'djac:',dfjac(2,:)
!!$
!!$       end if
!!$
!!$       return
!!$    end if
!!$
!!$    n = naper - nref + 1
!!$    allocate(rs(n),f(n),df(n),r(n))
!!$    q = p(2)
!!$    a = p(1)
!!$
!!$    do i = 1,n
!!$       r(i) = radius(i+nref-1)
!!$       rs(i) = (q - a/r(i) - grow_(i+nref-1)) / sig
!!$       f(i) = huber(rs(i))
!!$    end do
!!$    f = r
!!$
!!$    if( iflag == 1 ) then
!!$
!!$       fvec(1) = sum(f)
!!$       fvec(2) = -sum(f/r)
!!$       fvec = fvec / sig
!!$
!!$    else if( iflag == 2 ) then
!!$
!!$       do i = 1,n
!!$          df(i) = dhuber(rs(i))
!!$       end do
!!$       df = 1
!!$
!!$       fjac(1,1) = sum(df)
!!$       fjac(1,2) = -sum(df/r)
!!$       fjac(2,2) = sum(df/r**2)
!!$
!!$       fjac(2,1) = fjac(1,2)
!!$       fjac = fjac / sig**2
!!$
!!$    end if
!!$
!!$    deallocate(rs,f,df,r)
!!$
!!$  end subroutine funder3
!!$
!!$
!!$  !--
!!$
!!$
!!$  subroutine asymptote_init(a,b,c,q)
!!$
!!$    use NelderMead
!!$    use robustmean
!!$
!!$    real(dbl), intent(out) :: a,b,c
!!$    real(dbl), intent(out), optional :: q
!!$    real(dbl), dimension(5) :: p,dp
!!$    real(dbl) :: pmin,d,r
!!$    integer :: i,ifault,n
!!$
!!$    n = naper - nref + 1
!!$    p(1) = sum((1-grow_(nref:naper))/radius(nref:naper))/sum(1/radius(nref:naper)**2)
!!$!    p(1) = sum((1-grow_(nref:naper))/raper(nref:naper))/sum(1/raper(nref:naper)**2)
!!$!    p(1) = p(1) / 10
!!$    r = radius(naper)
!!$        write(*,*) n,p(1),p(1)/r
!!$    p(1) = 0
!!$    p(2) = 0
!!$    p(3) = 0
!!$    dp(1) = 0.001
!!$    dp(2) = 0.001
!!$    dp(3) = 0.001
!!$    n = naper - 1
!!$    r = raper(naper)
!!$
!!$!    stop
!!$
!!$    if( present(q) ) then
!!$       p(4) = 1e-5 !-p(1)/r  !max(-p(1)/r,epsilon(p))
!!$       dp(4) = 0.00001
!!$       p(5) = 1e-5
!!$       dp(5) = 0.1*p(5)
!!$       lambda = 1
!!$       do i = 1, 1
!!$          call nelmin1(fainit4,p,dp,pmin,ifault)
!!$          lambda = lambda/10
!!$          write(*,*) real(p),ifault
!!$       end do
!!$!       par_t = 0.99
!!$!       lambda = 1
!!$!       do i = 1, 10
!!$!          call nelmin1(fainit4a,p(1:3),dp(1:3),pmin,ifault,kcount0=1000*nstar,reqmin0=1d-10)
!!$!          lambda = 2*lambda
!!$!          call rmean(cts(1:3,n)/cts(1:3,naper),par_t,d)
!!$!          par_t = 1 - p(1)/r - p(2)/r**2 - p(3)/r**3
!!$!          write(*,*) real(p(1:3)),real(par_t)
!!$!       end do
!!$    else
!!$       call nelmin1(fainit3,p(1:3),dp(1:3),pmin,ifault)
!!$    end if
!!$
!!$    if( ifault /= 0 ) write(*,*)  'asymptote_init:',real(p),real(pmin),ifault
!!$
!!$    r = raper(naper)
!!$    a = p(1)*r
!!$    b = p(2)*r**2
!!$    c = p(3)*r**3
!!$
!!$    if( present(q) ) q = p(4)
!!$
!!$    write(*,'(5f15.5)') a,b,c,p(4),pmin
!!$!    stop
!!$  end subroutine asymptote_init
!!$
!!$  function fainit4(p)
!!$
!!$    real(dbl) :: fainit4
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl) :: s,a,b,c,e,r,x,f,r1,f1,lambda,h
!!$    integer :: i
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    c = p(3)
!!$    e = p(4)
!!$    lambda = p(5)
!!$
!!$!    if( q <= 0 ) then
!!$!       fainit4 = 1e5
!!$!       return
!!$!    end if
!!$
!!$    h = 0
!!$    s = 0
!!$    do i = nref-1, naper
!!$       r = radius(i)
!!$       f = 1 - (a + (b + c/r)/r)/r + e
!!$!       x = 1/r
!!$!       f = 1 - (a + (b + c*x)*x)*x + e
!!$!       r1 = raper(i+1)
!!$!       f = q - (a + (b + c/r)/r)/r
!!$!       f1 = q - (a + (b + c/r1)/r1)/r1
!!$       x = (grow_(i) - f) !/ (dgrow_(i) / raper(i)**2)!*r**2!* x**2
!!$       ! dgrow_ are not used because initial estimates have no Normal
!!$       ! distribution, errors includes systematic offsets which can
!!$       ! confuse fitting procedure
!!$       s = s + abs(x)
!!$       h = h + abs(a + 3*b/r + 6*c/r**2)/r**2
!!$!       write(*,*) i,(1 - (a + (b + c/r)/r)/r) / (1 - (a + (b + c/r1)/r1)/r1)
!!$       write(*,'(5f15.5)') r,f,x,grow_(i),dgrow_(i)
!!$!       write(*,'(5f15.5)') 1/r,grow_(i),a + 3*b/r + 6*c/r**2
!!$    end do
!!$!stop
!!$    r = raper(naper)
!!$    r = radius(naper)
!!$ !   fainit4 = s !+ lambda*(abs(a/r) + abs(b/r**2) + abs(c/r**3)) !(?) + abs(q-1)
!!$ !   fainit4 = s + lambda*abs(1 - (a + (b + c/r)/r)/r + q - 1)
!!$    !    fainit4 = s + lambda*(abs((a + (b + c/r)/r)/r - q) + abs(a) + abs(b) + abs(c))
!!$!    x = 1/raper(naper)
!!$    fainit4 = s + lambda*(max(0.0,-e) + 0*abs((a + (b + c/r)/r)/r - e) - h)
!!$
!!$
!!$  end function fainit4
!!$
!!$  function fainit4a(p)
!!$
!!$    real(dbl) :: fainit4a
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl) :: s,a,b,c,q,r,x,f,r1,f1,t
!!$    integer :: i,j
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    c = p(3)
!!$    !    q = p(4)
!!$    t = par_t
!!$
!!$    s = 0
!!$    do j = 1,1
!!$    do i = nref, naper-1
!!$       r = raper(i)
!!$       r1 = raper(i+1)
!!$       f = 1 - (a + (b + c/r)/r)/r
!!$       f1 = 1 - (a + (b + c/r1)/r1)/r1
!!$!       x = (grow_(i) - f/f1) / dgrow_(i)
!!$       x = (cts(j,i)/cts(j,i+1) - f/f1) !/ dcts(j,i)
!!$       s = s + abs(x)
!!$!       write(*,*) i,f/f1,x,cts(j,i)/cts(j,i+1)
!!$    end do
!!$ end do
!!$!stop
!!$    r = raper(naper)
!!$    fainit4a = s !+ lambda*(abs(a/r) + abs(b/r**2) + abs(c/r**3)) !(?) + abs(q-1)
!!$    fainit4a = s / 3 + lambda*abs(1 - a/r - b/r**2 - c/r**3 - t)
!!$  end function fainit4a
!!$
!!$
!!$  function fainit3(p)
!!$
!!$    real(dbl) :: fainit3
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl) :: s,a,b,c,r,x,r1
!!$    integer :: i
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    c = p(3)
!!$
!!$    s = 0
!!$    do i = nhwhm, naper-1
!!$       r = raper(i)
!!$       r1 = raper(i+1)
!!$       x = (grow_(i) - asymptote(r,a,b,c)/asymptote(r1,a,b,c)) / dgrow_(i)
!!$!       write(*,*) real(r),real(x),real(grow_(i)),real(asymptote(r,a,b,c)),real(dgrow_(i))
!!$       s = s + abs(x)
!!$    end do
!!$    r = raper(naper)
!!$    fainit3 = s !+ 0*lambda*(abs(a/r) + abs(b/r**2) + abs(c/r**3))
!!$!    write(*,*) real(s),real(p)
!!$!    stop 0
!!$
!!$  end function fainit3
!!$
!!$  subroutine anchor(cts,dcts,a,b,c,t,med,back)
!!$
!!$    use fmm
!!$
!!$    real(dbl), dimension(:), intent(in) :: cts,dcts
!!$    real(dbl), intent(in) :: a,b,c,t,med
!!$    real(dbl), intent(in out) :: back
!!$    real(dbl) :: d,tol,bmin,bmax
!!$    integer :: i
!!$
!!$!    par_a = a
!!$!    par_b = b
!!$    !    par_c = c
!!$    cts_ = cts
!!$    dcts_ = dcts
!!$    par_t = t
!!$    do i = nref, naper
!!$       grow_(i) = asymptote(raper(i),a,b,c)
!!$    end do
!!$
!!$    tol = sqrt(epsilon(cts))
!!$    tol = 1e-6
!!$
!!$    if( med > 0 ) then
!!$       d = 5*med*cts(naper)
!!$    else
!!$       d = sqrt(cts(naper))
!!$    end if
!!$    bmin = back - d
!!$    bmax = back + d
!!$!    d = back
!!$    back = fmin(bmin,bmax,ancora,tol)
!!$!    write(*,*) real(back),real(d),real(bmin),real(bmax)
!!$!stop 0
!!$  end subroutine anchor
!!$
!!$
!!$  function ancora(b)
!!$
!!$    use rfun
!!$
!!$    real(dbl), intent(in) :: b
!!$    real(dbl) :: ancora
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow
!!$    real(dbl) :: s,d,f
!!$    integer :: i
!!$    logical :: valid
!!$
!!$    call the_constructor(par_t,b,cts_,dcts_,grow,dgrow,valid)
!!$    if( .not. valid ) then
!!$       ancora = 1e3*(abs(b)+1)
!!$       return
!!$    end if
!!$
!!$    s = 0
!!$    do i = nref, naper
!!$       d = (grow(i) - grow_(i)) / dgrow(i)
!!$       if( type == 0 ) then
!!$          s = s + abs(d)
!!$       else if( type == 1 ) then
!!$          s = s + ihuber(d/sig)
!!$       end if
!!$!       write(*,'(4f10.5)') raper(i),grow(i),grow_(i),dgrow(i)
!!$    end do
!!$!    write(*,*) b,s
!!$    ancora = s
!!$
!!$  end function ancora
!!$
!!$! -- odpad
!!$
!!$
!!$
!!$  subroutine asymptote_update(a,b,c,t)
!!$
!!$    use NelderMead
!!$
!!$    real(dbl), intent(in out) :: a,b,c
!!$    real(dbl), intent(in) :: t
!!$    real(dbl), dimension(3) :: p,dp
!!$    real(dbl) :: pmin
!!$    integer :: i,ifault
!!$    logical :: valid
!!$
!!$    do i = 1,nstar
!!$       call the_constructor(t,back(i),cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$    end do
!!$
!!$    ! estimate parameters
!!$    p(1) = a
!!$    p(2) = b
!!$    p(3) = c
!!$    dp(1) = 0.001
!!$    dp(2) = 0.001
!!$    dp(3) = 0.00001
!!$    call nelmin1(minas,p,dp,pmin,ifault)
!!$ !   write(*,*) ifault
!!$
!!$!    valid = ifault == 0 .and. p(2) < 0.99*maxback
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    c = p(3)
!!$    write(*,*) real(p),real(pmin)
!!$
!!$  end subroutine asymptote_update
!!$
!!$  function minas(p)
!!$
!!$    use rfun
!!$
!!$    real(dbl) :: minas
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$!    real(dbl), dimension(size(raper)) :: grow,dgrow,f
!!$    real(dbl), dimension(size(raper)) :: f
!!$    real(dbl) :: s,d,a,b,c,x,bck,r
!!$    integer :: i,j
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    c = p(3)
!!$    bck = 0
!!$    do i = nref,naper
!!$       f(i) = asymptote(raper(i),a,b,c)
!!$    end do
!!$
!!$    s = 0
!!$    do j = 1,nstar
!!$!       call the_constructor(tx(j),bck,cts(j,:),dcts(j,:),grow,dgrow,valid)
!!$!       if( .not. valid ) then
!!$!          s = s + 1e3
!!$!       else
!!$          do i = nref, naper
!!$             x = (grows(j,i) - f(i)) / dgrows(j,i)
!!$             if( type == 0 ) then
!!$                d = abs(x)
!!$             else if( type == 1 ) then
!!$                d = ihuber(x/sig)
!!$             else
!!$                d = 1e3
!!$             end if
!!$             s = s + d
!!$          end do
!!$!       end if
!!$    end do
!!$
!!$    ! conditions: a,b,c -> 0, a > b/r(max)
!!$    r = raper(naper)
!!$    minas = s / nstar  !+ lambda*(abs(a/r) + abs(b/r**2) + abs(c/r**3))
!!$
!!$  end function minas
!!$
!!$  subroutine asymptote_sig(a,b,c,t,back,sig)
!!$
!!$    use weightedmean
!!$
!!$    real(dbl), intent(in) :: a,b,c,t
!!$    real(dbl), dimension(:), intent(in) :: back
!!$    real(dbl), intent(out) :: sig
!!$
!!$    real(dbl), dimension(size(raper)) :: grow,dgrow,f
!!$    real(dbl), dimension(:), allocatable :: fvec,dvec
!!$    real(dbl) :: x,s
!!$    integer :: i,j,n
!!$    logical :: valid
!!$
!!$    do i = nref,naper
!!$       f(i) = asymptote(raper(i),a,b,c)
!!$    end do
!!$
!!$    n = (naper - nref + 1)*nstar
!!$    allocate(fvec(n),dvec(n))
!!$    n = 0
!!$    do j = 1, nstar
!!$       call the_constructor(t,back(j),cts(j,:),dcts(j,:),grow,dgrow,valid)
!!$       do i = nref,naper
!!$          n = n + 1
!!$          fvec(n) = grow(i) - f(i)
!!$          dvec(n) = dgrow(i)
!!$       end do
!!$    end do
!!$!    write(*,*) n,nstar,real(fvec(1:n))
!!$    !    call rinit(fvec(1:n),d,s)
!!$    call rwmean(fvec,dvec,x,s,sig)
!!$    !    sig = d / 0.6745
!!$    write(*,*) 'sig:',x, s, sig
!!$    deallocate(fvec,dvec)
!!$
!!$  end subroutine asymptote_sig
!!$
!!$
!!$
!!$
!!$  ! laters ************************
!!$
!!$
!!$
!!$
!!$
!!$  subroutine commons1(a,b,c)
!!$
!!$    use NelderMead
!!$
!!$    real(dbl), intent(in out) :: a,b,c
!!$!    integer, intent(in) :: type
!!$    real(dbl), dimension(3) :: p,dp
!!$    real(dbl) :: pmin
!!$    integer :: ifault
!!$
!!$!    par_t = t
!!$
!!$    ! estimate parameters
!!$!    p(1) = t
!!$!    p(2) = a
!!$!    p(3) = b
!!$    p(1) = a
!!$    p(2) = b
!!$    p(3) = c
!!$!    p(4:5) = 1
!!$    dp(1) = 0.001
!!$    dp(2) = 0.001
!!$!    dp(3) = 0.001
!!$    dp(3) = 0.00001
!!$!    dp(3) = 0.001
!!$!    dp(4:5) = 0.1
!!$!    if( type == 0 ) then
!!$!       call nelmin1(mincom1,p,dp,pmin,ifault)
!!$!    else if( type == 1 ) then
!!$!       call nelmin1(mincom1,p,dp,pmin,ifault)
!!$!    end if
!!$!    call nelmin1(mincom3,p,dp,pmin,ifault)
!!$ !   write(*,*) ifault
!!$
!!$!    valid = ifault == 0 .and. p(2) < 0.99*maxback
!!$
!!$!    t = p(1)
!!$!    a = p(2)
!!$!    b = p(3)
!!$!    t = p(1)
!!$    a = p(1)
!!$    b = p(2)
!!$    c = p(3)
!!$!    e = p(3)
!!$!    e = 0
!!$!    t = 1 - e - a/raper(nref) - b/raper(nref)**2
!!$!    t = asymptote(raper(nref),e,a,b)
!!$
!!$    write(*,*) real(p),real(pmin)
!!$
!!$  end subroutine commons1
!!$
!!$
!!$  function mincom0(p)
!!$
!!$!    use robustmean
!!$!    use rfun
!!$
!!$    real(dbl) :: mincom0
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(size(raper)) :: grow,dgrow
!!$    real(dbl) :: s,e,b,d,a,r,g,t,K,q
!!$    integer :: i,j,n
!!$    logical :: valid
!!$
!!$!    t = p(1)
!!$    a = p(1)
!!$    b = p(2)
!!$    e = p(3)
!!$    t = 1 - e - a/raper(nref) - b/raper(nref)**2
!!$    n = (nhwhm + naper) / 2
!!$    n = nref
!!$
!!$    open(1,file='/tmp/7')
!!$    s = 0
!!$    q = 0
!!$    do j = 1,nstar
!!$       call the_constructor(t,back(j),cts(j,:),dcts(j,:),grow,dgrow,valid)
!!$       if( .not. valid ) then
!!$          s = s + 1
!!$       else
!!$          do i = n, naper
!!$             r = raper(i)
!!$             g = 1 - e - (a + b / r) / r
!!$             s = s + abs(grow(i) - g)
!!$             write(1,*) r,grow(i) - g,g,grow(i)
!!$          end do
!!$          do i = nhwhm, naper
!!$             e = sqrt(dgrow(i)**2 + dgrow(i-1)**2)
!!$             q = q + max(-(grow(i) - grow(i-1) - e),0.0)**2
!!$          end do
!!$       end if
!!$    end do
!!$    close(1)
!!$
!!$    d = a - b / raper(naper)
!!$    K = a*(1/raper(naper)**2 - 1/raper(n)**2) + 2*b*(1/raper(naper)**3 - 1/raper(n)**3)
!!$    K = -K
!!$!    d = a + 2*b/raper(n)
!!$!    d = 0
!!$!    do i = n,naper
!!$!       d = d + max(0.0,-(a + 2*b/raper(i)))
!!$!    end do
!!$    ! conditions: e,a,b > 0, a > b/r(max)
!!$!    mincom0 = s + 0*lambda*max(1-t,0.0) + max(0.0,-a) + max(0.0,-b) + max(0.0,-d)
!!$!    mincom0 = s + lambda*(max(1-t,0.0)**2 + max(0.0,-a)**2 + max(0.0,-b)**2 + max(0.0,-d)**2)
!!$!    mincom0 = s + lambda*(a**2 + b**2 + max(0.0,-a)**2 + max(0.0,-b)**2 + max(0.0,-d)**2)
!!$    mincom0 = s + lambda*(max(0.0,-a)**2 + max(0.0,-b)**2 + max(0.0,-e)**2 + 0*max(0.0,-d)**2   &
!!$         + 0*max(0.0,-K)**2 + 0*max(1-t,0.0)**2 + 0*q/nstar)
!!$
!!$  end function mincom0
!!$
!!$  function mincom(p)
!!$
!!$    use rfun
!!$
!!$    real(dbl) :: mincom
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(size(raper)) :: grow,dgrow,f
!!$    real(dbl) :: s,b,d,a,c,t,x,bck,w,e
!!$    integer :: i,j
!!$    logical :: valid
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    c = p(3)
!!$!    e = p(3)
!!$!    e = 0
!!$    bck = 0
!!$!    t = asymptote(raper(nref),e,a,b)
!!$    do i = nref,naper
!!$       f(i) = asymptote(raper(i),a,b,c)
!!$    end do
!!$    !    t = f(nref)
!!$    t = par_t
!!$
!!$!    open(1,file='/tmp/9')
!!$    s = 0
!!$    !    w = e
!!$    w = 0
!!$    do j = 1,nstar
!!$       call the_constructor(t,bck,cts(j,:),dcts(j,:),grow,dgrow,valid)
!!$       if( .not. valid ) then
!!$          s = s + 1e3
!!$       else
!!$          do i = nref, naper
!!$             x = (grow(i) - f(i)) / dgrow(i)
!!$             if( type == 0 ) then
!!$                d = abs(x)
!!$             else if( type == 1 ) then
!!$                d = ihuber(x/sig)
!!$             else
!!$                d = 1e3
!!$             end if
!!$             s = s + d
!!$!             w = w + sqrt(dgrow(i)**2 + dgrow(i-1)**2)
!!$!             write(1,*) raper(i),d,f(i),grow(i)
!!$          end do
!!$       end if
!!$    end do
!!$!    close(1)
!!$
!!$
!!$
!!$
!!$    !        write(*,*) real(e),real(a),real(t),real(s)
!!$    ! conditions: e,a,b > 0, a > b/r(max)
!!$!    write(*,*) sum(abs(t-grow(:,nref)))
!!$!    d = (a + 2*b/raper(naper))/raper(naper)**2
!!$    mincom = s / nstar  &
!!$         + lambda*(max(0.0,-a)**2 + max(0.0,-b)**2 + max(0.0,-e)**2)
!!$    !+ max(0.0,-(a+2*b/raper(naper)))**2) !&
!!$!         + a**2 + b**2 + e**2)
!!$!         + max(0.0,-w)**2*0 + a**2 + b**2 + e**2)
!!$
!!$  end function mincom
!!$
!!$  function asymptote(r,a,b,c)
!!$
!!$    real(dbl), intent(in) :: r,a,b,c
!!$    real(dbl) :: asymptote
!!$
!!$    asymptote = 1 - (a + (b + c/r)/r)/r
!!$
!!$  end function asymptote
!!$


!!$  subroutine arc_init(cts,dcts,flux0,a,b,beta)
!!$
!!$    use NelderMead
!!$
!!$    real(dbl), dimension(:), intent(in) :: cts,dcts
!!$    real(dbl), intent(in) :: flux0
!!$    real(dbl), intent(out) :: a,b
!!$    real(dbl), intent(out) :: beta
!!$    real(dbl), dimension(3) :: p,dp
!!$    real(dbl) :: pmin
!!$    integer :: ifault,i
!!$
!!$    ! here, all dgrows are defined and leaved to no change during computations
!!$    do i = naper-1,nref,-1
!!$       grow_(i) = (cts(i+1) - cts(i)) / flux0
!!$    end do
!!$    dgrow_ = dcts(nhwhm) / flux0
!!$
!!$    p = 0
!!$    dp = 1e-3
!!$    call nelmin1(arcfun,p,dp,pmin,ifault)
!!$
!!$    a = p(1)*raper(naper)
!!$    b = p(2)*raper(naper)**2
!!$    beta = p(3)*flux0/(pi*raper(naper)**2)
!!$
!!$!    write(*,'(3f12.5,2g16.5,i3,a)') a,b,beta,pmin,dgrow_(1),ifault,' ppp'
!!$!    write(*,*) real(p(3)*flux0/pi/raper(naper)**2)
!!$!    write(*,'(5000f10.2)') p(3:)*flux0/pi
!!$!    write(*,*) ifault
!!$!    stop
!!$
!!$  end subroutine arc_init
!!$
!!$
!!$  function arcfun(p) result(s)
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl) :: s,f,f2,a,b,r,r1,beta
!!$    integer :: i,n
!!$    logical :: valid
!!$
!!$    if( p(1) < 0 .or. p(2) < 0 ) then
!!$       s = 1000*nstar
!!$      return
!!$    end if
!!$
!!$    a = p(1)
!!$    b = p(2)
!!$    beta = p(3)
!!$
!!$    s = 0
!!$    f2 = 0
!!$    do i = nref, naper-1
!!$       r = radius(i)
!!$       r1 = radius(i+1)
!!$       f = beta*(r1**2 - r**2) + a*(1/r - 1/r1) + b*(1/r**2 - 1/r1**2)
!!$       f2 = f2 - 2*(a + 3*b/r)/r**3
!!$       s = s + abs(grow_(i) - f) / dgrow_(i)
!!$       !write(*,'(6f12.5)') real(grows(n,i) - f), real(f),real(grows(n,i)),a,b,beta(n)
!!$    end do
!!$   s = s / (naper - nref - 1)
!!$
!!$ end function arcfun
!!$
!!$
!!$
!!$
!!$
!!$  subroutine anchor_init(cts,dcts,flux0,a,b,beta,t,phi)
!!$
!!$    use fmm
!!$
!!$    real(dbl), dimension(:), intent(in) :: cts,dcts
!!$    real(dbl), intent(in) :: flux0,a,b,beta
!!$    real(dbl), intent(out) :: t,phi
!!$    real(dbl) :: tol, tmin, tmax
!!$    integer :: ifault,i
!!$
!!$    ! here, all dgrows are defined and leaved to no change during computations
!!$    do i = nref,naper
!!$       grow_(i) = 1 - cts(i) / flux0 - beta*aper(i)/flux0 + a/raper(i) + b/raper(i)**2
!!$!       write(*,*) real(grow_(i)),real(cts(i) / flux0),real(beta*aper(i)/flux0), &
!!$!            real(a/raper(i)),real(b/raper(i)**2)
!!$    end do
!!$    dgrow_ = dcts(nhwhm) / flux0
!!$
!!$    tol = sqrt(epsilon(cts))
!!$    tol = 1e-6
!!$    tmin = minval(grow_(nref:naper))
!!$    tmax = maxval(grow_(nref:naper))
!!$!    write(*,*) real(tmin), real(tmax)
!!$    phi = fmin(tmin,tmax,anchorfun,tol)
!!$
!!$!    write(*,*) real(t),' anchor_init'
!!$
!!$  end subroutine anchor_init
!!$
!!$
!!$  function anchorfun(t) result(s)
!!$
!!$    real(dbl), intent(in) :: t
!!$
!!$    real(dbl) :: s,f
!!$    integer :: i
!!$
!!$    s = 0
!!$    do i = nref, naper
!!$       s = s + abs(grow_(i) - t) / dgrow_(i)
!!$       !write(*,'(6f12.5)') real(grows(n,i) - f), real(f),real(grows(n,i)),a,b,beta(n)
!!$    end do
!!$   s = s / (naper - nref - 1)
!!$
!!$ end function anchorfun
!!$
!!$



!!$  function ancora(t)
!!$
!!$    use rfun
!!$
!!$    real(dbl), intent(in) :: t
!!$    real(dbl) :: ancora
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow
!!$    real(dbl) :: s,d,b,f
!!$    integer :: i
!!$    logical :: valid
!!$
!!$    b = 0
!!$    call the_constructor(t,b,cts_,dcts_,grow,dgrow,valid)
!!$
!!$    s = 0
!!$    do i = nref, naper
!!$       d = (grow(i) - grow_(i)) / dgrow(i)
!!$       if( type == 0 ) then
!!$          s = s + abs(d)
!!$       else if( type == 1 ) then
!!$          s = s + ihuber(d/sig)
!!$       end if
!!$    end do
!!$
!!$    ancora = s
!!$
!!$  end function ancora


!!$  function asymptote1(r,a,b,c)
!!$
!!$    real(dbl), intent(in) :: r,a,b,c
!!$    real(dbl) :: asymptote1
!!$
!!$    asymptote1 = 1 - (a + (b + c/r)/r)/r
!!$
!!$  end function asymptote1
!!$
!!$  subroutine the_constructor(t,b,cts,dcts,grow,dgrow,valid)
!!$
!!$    real(dbl), intent(in) :: t,b
!!$    real(dbl), dimension(:), intent(in) :: cts, dcts
!!$    real(dbl), dimension(:), intent(out) :: grow,dgrow
!!$    logical, intent(out) :: valid
!!$
!!$    real(dbl), dimension(size(raper)) :: flux,dflux
!!$    integer :: i,n
!!$
!!$    valid = .false.
!!$    n = size(grow)
!!$
!!$    flux = cts - b*aper
!!$    dflux = dcts
!!$
!!$    if( any(flux < epsilon(flux)) ) return
!!$
!!$    ! set up initial point of grow curve
!!$!    grow(nhwhm) = t
!!$!    grow(nref) = t
!!$    grow(naper) = t
!!$
!!$    do i = naper-1,1,-1
!!$       grow(i) = grow(i+1) * (flux(i) / flux(i+1))
!!$    end do
!!$
!!$
!!$    ! compute grow curve
!!$    goto 33
!!$!    do i = nhwhm-1,1,-1
!!$    do i = nref-1,1,-1
!!$       grow(i) = grow(i+1) * (flux(i) / flux(i+1))
!!$    end do
!!$!    do i = nhwhm+1,n
!!$    do i = nref+1,n
!!$       grow(i) = grow(i-1) * (flux(i) / flux(i-1))
!!$    end do
!!$
!!$33 continue
!!$
!!$!    forall( i = 1:n ) dgrow(i) = grow(i)*(dflux(i)/flux(i))
!!$!    dgrow = max(grow * (dflux / flux), 1e-7)
!!$    dgrow = grow * (dflux / flux)
!!$
!!$    valid = .true.
!!$
!!$  end subroutine the_constructor
!!$
!!$  subroutine a_constructor(flux,dflux,grow,dgrow,valid)
!!$
!!$    real(dbl), dimension(:), intent(in) :: flux, dflux
!!$    real(dbl), dimension(:), intent(out) :: grow,dgrow
!!$    logical, intent(out) :: valid
!!$
!!$    integer :: i
!!$
!!$    valid = .false.
!!$
!!$    if( any(flux < epsilon(flux)) ) return
!!$
!!$    ! set up initial point of grow curve
!!$    grow(nhwhm) = 1
!!$
!!$    ! compute grow curve
!!$    do i = nhwhm-1,1,-1
!!$       grow(i) = grow(i+1) * (flux(i) / flux(i+1))
!!$    end do
!!$    do i = nhwhm+1, size(grow)
!!$       grow(i) = grow(i-1) * (flux(i) / flux(i-1))
!!$    end do
!!$
!!$    dgrow = grow * (dflux / flux)
!!$
!!$    valid = .true.
!!$
!!$  end subroutine a_constructor
!!$
!!$
!!$  subroutine part(e,a,b,t)
!!$
!!$    use fmm
!!$
!!$    real(dbl), intent(in) :: e,a,b
!!$    real(dbl), intent(in out) :: t
!!$    real(dbl) :: dt,tol,tmin,tmax
!!$
!!$    par_e = e
!!$    par_a = a
!!$    par_b = b
!!$!    par_c = c
!!$    dt = 1e-2
!!$!    write(*,*) bmax,cts_(nhwhm),par_e,par_a,par_b
!!$!    tol = sqrt(epsilon(cts))*maxval(cts)
!!$    tol = sqrt(epsilon(cts))
!!$    tol = 1e-6
!!$    tmin = 0.5
!!$    tmax = min(1 - tol,asymptote(raper(naper),e,a,b))
!!$!    t = fmin(t-dt,max(t+dt,1.0),tfun,tol)
!!$    t = fmin(tmin,tmax,tfun,tol)
!!$!    if( back > 100 ) then
!!$!       valid = .false.
!!$!    end if
!!$
!!$  end subroutine part
!!$
!!$
!!$  function tfun(t)
!!$
!!$    use rfun
!!$
!!$    real(dbl), intent(in) :: t
!!$    real(dbl) :: tfun
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow, f
!!$    real(dbl) :: s,g,d,x
!!$    integer :: i,j
!!$    logical :: valid
!!$
!!$
!!$!    if( .not. valid ) then
!!$!       tfun = 666
!!$!       return
!!$!    end if
!!$
!!$!    do i = 1,naper
!!$!       write(*,*) raper(i),grow(i)
!!$!    end do
!!$    !    stop 0
!!$
!!$    do i = nref,naper
!!$       f(i) = asymptote(raper(i),par_e,par_a,par_b)
!!$    end do
!!$
!!$    s = 0
!!$    do j = 1,nstar
!!$       call the_constructor(t,back(j),cts(j,:),dcts(j,:),grow,dgrow,valid)
!!$!       if( .not. valid ) s = s + 100
!!$!       d = grow(naper) - f(naper)
!!$!       if( d ) s = s + 1000*()
!!$       do i = nref,naper
!!$!          r = raper(i) !/ hwhm
!!$!          g = 1 - par_e - par_a / r - par_b / r**2
!!$          !          s = s + abs(grow(i) - f(i)) ! / sqrt(dcts**2 + ..?
!!$          x = (grow(i) - f(i)) / dgrow(i)
!!$          if( type == 0 ) then
!!$             d = abs(x)
!!$          else if( type == 1 ) then
!!$             d = ihuber(x/sig)
!!$          else
!!$             d = 1e3
!!$          end if
!!$          s = s + d
!!$!          if( x > 3 ) s = s + 10
!!$       end do
!!$    end do
!!$    tfun = s / nstar
!!$
!!$!    write(*,*) real(t),real(s),t
!!$!    stop 0
!!$
!!$  end function tfun


!!$  subroutine tnoise(t)
!!$
!!$    use fmm
!!$
!!$    real(dbl), intent(out) :: t
!!$    real(dbl) :: tmin,tmax,tol
!!$
!!$    tmin = 0.5
!!$    tmax = 1 - epsilon(tmax)
!!$    tol = 1e-6
!!$    t = fmin(tmin,tmax,noisefun,tol)
!!$!    if( back > 100 ) then
!!$!       valid = .false.
!!$!    end if
!!$
!!$  end subroutine tnoise
!!$
!!$  function noisefun(t)
!!$
!!$    use robustmean
!!$
!!$    real(dbl), intent(in) :: t
!!$    real(dbl) :: noisefun
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow
!!$    real(dbl) :: s,back
!!$    integer :: n,i,j
!!$    logical :: valid
!!$
!!$    back = 0
!!$    n = (nhwhm + naper) / 2
!!$!    t = 1
!!$    do i = 1, nstar
!!$       call the_constructor(t,back,cts(i,:),dcts(i,:),grows(i,:),dgrows(i,:),valid)
!!$    end do
!!$    do i = nref, naper
!!$       call rinit(grows(:,i),grow(i),dgrow(i))
!!$    end do
!!$    s = 0
!!$    do i = nref+1, naper
!!$       !       s = s + abs(dgrow(i) - dgrow(i-1))/(grow(i) + grow(i-1)) !/(raper(i) - raper(i-1))
!!$       s = s + dgrow(i)
!!$    end do
!!$!    write(*,*) s,n
!!$    noisefun = s / (naper - nref)
!!$
!!$  end function noisefun


! ----




!!$  subroutine growmodel1(grows,dgrows,grow,dgrow)
!!$
!!$    use robustmean
!!$
!!$    real(dbl), dimension(:,:), intent(in) :: grows,dgrows
!!$    real(dbl), dimension(:), intent(out) :: grow,dgrow
!!$    real(dbl), dimension(:), allocatable :: f
!!$    integer :: n,i,j,l,naper,nstars
!!$
!!$    nstars = size(grows,1)
!!$    naper = size(grows,2)
!!$
!!$    ! initial estimate of growth curve as mean of all grow-curves
!!$    do i = 1, naper
!!$       call rmean(grows(:,i),grow(i),dgrow(i))
!!$!       write(*,*) i,real(grow(i)),real(dgrow(i))
!!$    end do
!!$
!!$
!!$    allocate(f(naper*nstars))
!!$
!!$    do i = 1,naper
!!$
!!$       n = 0
!!$       do l = 1, nstars
!!$          do j = 1,naper
!!$             n = n + 1
!!$             f(n) = grows(l,i)/grows(l,j) * grow(j)
!!$          end do
!!$       end do
!!$       if( n /= nstars*naper ) stop 'n /= nstars*naper'
!!$
!!$       call rmean(f,grow(i),dgrow(i))
!!$
!!$    end do
!!$
!!$    deallocate(f)
!!$
!!$    ! check for identical curves
!!$    where( abs(dgrow) < epsilon(dgrow) )
!!$       dgrow = 1e-5
!!$    end where
!!$
!!$  end subroutine growmodel1


!!$  function background1(e,a,d)
!!$
!!$    use fmm
!!$
!!$    real(dbl), intent(in) :: e,a,d
!!$    real(dbl) :: b, bmax, tol, background1
!!$
!!$    par_e = e
!!$    par_a = a
!!$    par_b = d
!!$    bmax = 0.1*cts_(nhwhm)
!!$    write(*,*) bmax,cts_(nhwhm),par_e,par_a,par_b
!!$!    tol = sqrt(epsilon(cts))*maxval(cts)
!!$    tol = sqrt(epsilon(cts))
!!$    b = fmin(-bmax,bmax,backfun1,tol)
!!$!    if( b > 100 ) then
!!$!       valid = .false.
!!$!    end if
!!$
!!$    background1 = b
!!$
!!$  end function background1
!!$
!!$  function backfun1(b)
!!$
!!$    real(dbl), intent(in) :: b
!!$    real(dbl) :: backfun1
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow
!!$    real(dbl) :: s,d,x
!!$    integer :: i
!!$    logical :: valid
!!$
!!$    call the_constructor1(par_e,par_a,cts_,dcts_,b,grow,dgrow,valid)
!!$
!!$    if( .not. valid ) then
!!$       backfun1 = 666
!!$       return
!!$    end if
!!$
!!$!    do i = 1,naper
!!$!       write(*,*) raper(i),grow(i)
!!$!    end do
!!$!    stop 0
!!$    s = 0
!!$    do i = (nhwhm+naper)/2,naper
!!$       x = raper(i) !/ hwhm
!!$       d = 1 - par_e - par_a / x - par_b / x**2
!!$       s = s + abs(grow(i) - d) ! / sqrt(dcts**2 + ..?
!!$    end do
!!$    backfun1 = s
!!$
!!$    write(*,*) real(b),real(s)
!!$
!!$  end function backfun1
!!$
!!$    subroutine the_constructor1(e,a,cts,dcts,b,grow,dgrow,valid)
!!$
!!$    real(dbl), intent(in) :: a,e,b
!!$    real(dbl), dimension(:), intent(in) :: cts, dcts
!!$    real(dbl), dimension(:), intent(out) :: grow,dgrow
!!$    logical, intent(out) :: valid
!!$
!!$    real(dbl), dimension(size(raper)) :: flux,dflux
!!$    integer :: i
!!$    real(dbl) :: t
!!$
!!$    valid = .false.
!!$
!!$    flux = cts - b*aper
!!$    dflux = dcts
!!$
!!$    if( any(flux < epsilon(flux)) ) return
!!$
!!$    ! set up initial point of grow curve
!!$    t = 1 - e - a / raper(nhwhm) - 0.13/ raper(nhwhm)**2
!!$    grow(nhwhm) = t
!!$    dgrow(nhwhm) = t * dcts(nhwhm)/cts(nhwhm)
!!$
!!$    ! compute grow curve
!!$    do i = nhwhm-1,1,-1
!!$       grow(i) = grow(i+1)*(flux(i) / flux(i+1))
!!$       dgrow(i) = grow(i)*sqrt((dcts(i)/cts(i))**2 + (dcts(i+1)/cts(i+1))**2)/1.41
!!$    end do
!!$
!!$    do i = nhwhm+1,size(grow)
!!$       grow(i) = grow(i-1)*(flux(i) / flux(i-1))
!!$       dgrow(i) = grow(i-1)*sqrt((dcts(i)/cts(i))**2 + (dcts(i-1)/cts(i-1))**2)/1.41
!!$    end do
!!$
!!$    valid = .true.
!!$
!!$  end subroutine the_constructor1
!!$
!!$
!!$
!!$  function xbackground(e,a,b,t)
!!$
!!$    use fmm
!!$
!!$    real(dbl), intent(in) :: e,a,b,t
!!$    real(dbl) :: back, bmax, tol, xbackground
!!$
!!$!    par_e = e
!!$    par_a = a
!!$    par_b = b
!!$    par_t = 1 - a / raper(nhwhm) - b / raper(nhwhm)**2
!!$    bmax = 0.1*cts_(nhwhm)
!!$!    bmax = 10
!!$!    write(*,*) bmax,cts_(nhwhm),par_e,par_a,par_b
!!$!    tol = sqrt(epsilon(cts))*maxval(cts)
!!$    tol = sqrt(epsilon(cts))
!!$    back = fmin(-bmax,bmax,xbackfun,tol)
!!$!    if( back > 100 ) then
!!$!       valid = .false.
!!$!    end if
!!$
!!$    xbackground = back
!!$
!!$  end function xbackground
!!$
!!$  function xbackfun(b)
!!$
!!$    real(dbl), intent(in) :: b
!!$    real(dbl) :: xbackfun
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow
!!$    real(dbl) :: s,g,r
!!$    integer :: i
!!$    logical :: valid
!!$
!!$    call the_constructor(par_t,b,cts_,dcts_,grow,dgrow,valid)
!!$
!!$    if( .not. valid ) then
!!$       xbackfun = 666
!!$       return
!!$    end if
!!$
!!$!    do i = 1,naper
!!$!       write(*,*) raper(i),grow(i)
!!$!    end do
!!$!    stop 0
!!$    s = 0
!!$    do i = (nhwhm+naper)/2,naper
!!$       r = raper(i) !/ hwhm
!!$!       g = 1 - par_e - par_a / r - par_b / r**2
!!$       g = 1 - par_a / r - par_b / r**2
!!$       s = s + abs(grow(i) - g) ! / sqrt(dcts**2 + ..?
!!$    end do
!!$    xbackfun = s
!!$
!!$!    write(*,*) real(b),real(s)
!!$
!!$  end function xbackfun
!!$
!!$
!!$  function ybackground(e,a,b,t)
!!$
!!$    use fmm
!!$
!!$    real(dbl), intent(in) :: e,a,b,t
!!$    real(dbl) :: back, bmax, tol, ybackground
!!$
!!$    par_a = a
!!$    par_b = b
!!$    par_t = t
!!$    par_t = 1 - a / raper(nhwhm) - b / raper(nhwhm)**2
!!$    bmax = 0.1*cts_(nhwhm)
!!$    bmax = 50
!!$    tol = sqrt(epsilon(cts))
!!$    back = fmin(-bmax,bmax,ybackfun,tol)
!!$!    if( back > 100 ) then
!!$!       valid = .false.
!!$!    end if
!!$
!!$    ybackground = back
!!$
!!$  end function ybackground
!!$
!!$  function ybackfun(b)
!!$
!!$    use rfun
!!$
!!$    real(dbl), intent(in) :: b
!!$    real(dbl) :: ybackfun
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow
!!$    real(dbl) :: s,g,r,x
!!$    integer :: i
!!$    logical :: valid
!!$
!!$    call the_constructor(par_t,b,cts_,dcts_,grow,dgrow,valid)
!!$
!!$    if( .not. valid ) then
!!$       ybackfun = 666
!!$       return
!!$    end if
!!$
!!$    s = 0
!!$    do i = (nhwhm+naper)/2,naper
!!$       r = raper(i) !/ hwhm
!!$!       g = 1 - par_e - par_a / r - par_b / r**2
!!$       g = 1 - par_a / r - par_b / r**2
!!$       x = (grow(i) - g) / (sig*dgrow(i))
!!$       s = s + ihuber(x) ! / sqrt(dcts**2 + ..?
!!$    end do
!!$    ybackfun = s
!!$
!!$!    write(*,*) real(b),real(s)
!!$
!!$  end function ybackfun



!!$  subroutine background2(e,a,b,tx,backs)
!!$
!!$    use NelderMead
!!$    use robustmean
!!$
!!$    real(dbl), intent(in) :: e,a,b
!!$    real(dbl), intent(in out) :: tx
!!$    real(dbl), dimension(:), intent(in out) :: backs
!!$    real(dbl), dimension(:), allocatable :: p,dp
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow, f, r
!!$    real(dbl) :: pmin, d, s
!!$    integer :: n, i, npar, ifault
!!$    logical :: valid
!!$
!!$    par_e = e
!!$    par_a = a
!!$    par_b = b
!!$
!!$    npar = nstar + 1
!!$    allocate(p(npar),dp(npar))
!!$
!!$    p(1) = tx
!!$    p(2:) = backs
!!$    dp(1) = 0.01
!!$    dp(2:) = 0.1
!!$    lambda = 1
!!$    do i = 1, 1
!!$       call nelmin1(backfun2,p,dp,pmin,ifault)
!!$       write(*,*) real(pmin),ifault,real(lambda)
!!$       lambda = 2*lambda
!!$    end do
!!$    tx = p(1)
!!$    backs = p(2:)
!!$!    stop 0
!!$    deallocate(p,dp)
!!$
!!$  end subroutine background2
!!$
!!$  function backfun2(p)
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$    real(dbl) :: backfun2
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow, f
!!$    real(dbl) :: s,g,r
!!$    integer :: n,i,j
!!$    logical :: valid
!!$
!!$    n = (nhwhm + naper) / 2
!!$    call asymptote(n,par_e,par_a,par_b,f)
!!$
!!$    s = 0
!!$    open(1,file='/tmp/c')
!!$    do i = 1, nstar
!!$       call the_constructor(p(1),p(i+1),cts(i,:),dcts(i,:),grow,dgrow,valid)
!!$       if( valid ) then
!!$          s = s + sum(abs(grow(n:) - f(n:)))
!!$          do j = n,naper
!!$             write(1,*) raper(j),abs(grow(j) - f(j)),grow(j),f(j)
!!$          end do
!!$       else
!!$          s = s + naper
!!$       end if
!!$    end do
!!$    close(1)
!!$
!!$    backfun2 = s / (nstar*(naper - n)) !+ lambda*(max(0.0,1-p(1))) !+ sum(abs(p(2:))))
!!$
!!$  end function backfun2
!!$
!!$
!!$
!!$  subroutine background(e,a,b,back,tx)
!!$
!!$    use NelderMead
!!$    use robustmean
!!$
!!$    real(dbl), intent(in) :: e,a,b
!!$    real(dbl), intent(in out) :: back,tx
!!$    real(dbl), dimension(2) :: p,dp
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow, f, r
!!$    real(dbl) :: pmin, d, s
!!$    integer :: n, ifault
!!$    logical :: valid
!!$
!!$    par_e = e
!!$    par_a = a
!!$    par_b = b
!!$
!!$    p(1) = tx
!!$    p(2) = back
!!$    dp(1) = 0.01
!!$    dp(2) = 1
!!$    call nelmin1(backfun,p,dp,pmin,ifault)
!!$    tx = p(1)
!!$    back = p(2)
!!$
!!$    n = (nhwhm + naper) / 2
!!$    call the_constructor(tx,back,cts_,dcts_,grow,dgrow,valid)
!!$    call asymptote(n,e,a,b,f)
!!$    r(n:) = abs(grow(n:) - f(n:)) / dgrow(n:)
!!$    call rinit(r(n:),d,s)
!!$    sig = d
!!$
!!$    call nelmin1(backas,p,dp,pmin,ifault)
!!$    write(*,*) real(p), real(sig)
!!$    tx = p(1)
!!$    back = p(2)
!!$
!!$  end subroutine background
!!$
!!$  function backfun(p)
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$    real(dbl) :: backfun
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow, f
!!$    real(dbl) :: s,g,r
!!$    integer :: n
!!$    logical :: valid
!!$
!!$    call the_constructor(p(1),p(2),cts_,dcts_,grow,dgrow,valid)
!!$
!!$    if( .not. valid ) then
!!$       backfun = 666
!!$       return
!!$    end if
!!$
!!$!    do i = 1,naper
!!$!       write(*,*) raper(i),grow(i)
!!$!    end do
!!$    !    stop 0
!!$
!!$    n = (nhwhm + naper) / 2
!!$    call asymptote(n,par_e,par_a,par_b,f)
!!$    backfun = sum(abs(grow(n:) - f(n:))) / (naper - n)
!!$
!!$    s = 0
!!$    do i = (nhwhm+naper)/2,naper
!!$       r = raper(i) !/ hwhm
!!$       g = 1 - (par_e + (par_a + par_b / r ) / r)
!!$       s = s + abs(grow(i) - g) !/ dgrow(i) ! / sqrt(dcts**2 + ..?
!!$!       write(*,*) real(x),real(d),real(grow(i))
!!$    end do
!!$    backfun = s
!!$
!!$!    write(*,*) real(p),real(s)
!!$
!!$  end function backfun
!!$
!!$  function backas(p)
!!$
!!$    use rfun
!!$
!!$    real(dbl), dimension(:), intent(in) :: p
!!$    real(dbl) :: backas
!!$
!!$    real(dbl), dimension(size(raper)) :: grow, dgrow, f, r, rs, rho
!!$!    real(dbl) :: s,g,r
!!$    integer :: n
!!$    logical :: valid
!!$
!!$    call the_constructor(p(1),p(2),cts_,dcts_,grow,dgrow,valid)
!!$
!!$    if( .not. valid ) then
!!$       backas = 666*naper
!!$       return
!!$    end if
!!$
!!$!    do i = 1,naper
!!$!       write(*,*) raper(i),grow(i)
!!$!    end do
!!$    !    stop 0
!!$
!!$    n = (nhwhm + naper) / 2
!!$    call asymptote(n,par_e,par_a,par_b,f)
!!$    r = (grow - f) / dgrow
!!$    rs = r / sig
!!$    call ihubers(rs,rho)
!!$    backas = sum(rho(n:)) / (naper - n)
!!$
!!$  end function backas
!!$
!!$
!!$  subroutine commons(e,a,b)
!!$
!!$    use NelderMead
!!$
!!$    real(dbl), intent(in out) :: e,a,b
!!$    real(dbl), dimension(3) :: p,dp
!!$    real(dbl) :: pmin
!!$    integer :: ifault
!!$
!!$    ! estimate parameters
!!$    p(1) = e
!!$    p(2) = a
!!$    p(3) = b
!!$    dp(1) = 0.0001
!!$    dp(2) = 0.001
!!$    dp(3) = 0.01
!!$    call nelmin1(mincom,p,dp,pmin,ifault)
!!$!    valid = ifault == 0 .and. p(2) < 0.99*maxback
!!$
!!$    e = p(1)
!!$    a = p(2)
!!$    b = p(3)
!!$
!!$!    write(*,*) real(p),real(pmin)
!!$
!!$  end subroutine commons


!!$  function mincom(p)
!!$
!!$    use robustmean
!!$    use rfun
!!$
!!$    real(dbl) :: mincom
!!$    real(dbl), dimension(:), intent(in) :: p
!!$
!!$    real(dbl), dimension(size(raper)) :: grow,dgrow
!!$    real(dbl) :: s,e,b,d,a,r,g
!!$    integer :: i,j,n
!!$    logical :: valid
!!$
!!$    e = p(1)
!!$    a = p(2)
!!$    b = p(3)
!!$    n = (nhwhm + naper) / 2
!!$
!!$    s = 0
!!$    do j = 1,nstar
!!$       call the_constructor(tx(j),back(j),cts(j,:),dcts(j,:),grow,dgrow,valid)
!!$       if( .not. valid ) s = s + 1
!!$       do i = n, naper
!!$          r = raper(i)
!!$          g = 1 - (e + (a + b / r) / r)
!!$          s = s + abs(grow(i) - g)
!!$!          write(*,*) real(grow(i)),real(dgrow(i)),real((grow(i)-g)/dgrow(i))
!!$!          write(*,*) real(grow(i)),real(dgrow(i)),real((grow(i)-g)/dgrow(i))
!!$!          s = s + huber((grow(i) - g)/dgrow(i))
!!$       end do
!!$    end do
!!$
!!$    d = a - b / raper(naper)
!!$    ! conditions: e,a,b > 0, a > b/r(max)
!!$    mincom = s + lambda*max(0.0,-e) + max(0.0,-a) + max(0.0,-b) + max(0.0,-d)
!!$
!!$
!!$!    if( e < 0 .or. a < 0 .or. b < 0 .or. x < 0 ) mincom = sum(grow_(n:naper))
!!$!    if( a < 0 .or. b < 0 .or. a < b / raper(naper) ) mincom = sum(grow_(n:naper))
!!$

!!$    return
!!$
!!$    n = (nhwhm + naper) / 2
!!$    s = 0
!!$    do k = 1, nstar
!!$
!!$       call the_constructor(e,a,cts(k,:),dcts(k,:),back(k),grow,dgrow,valid)
!!$       if( .not. valid ) s = s + 1e3
!!$       grows(k,:) = grow
!!$       dgrows(k,:) = dgrow
!!$
!!$
!!$       do i = n, naper
!!$!       do i = nhwhm, naper
!!$          x = raper(i) !/ hwhm
!!$          d = 1 + e - a / x
!!$!          d = 1  - a / x
!!$          !d = 1
!!$          s = s + abs(grow(i) - d) !/ dgrow(i)
!!$!          if( grow(i) <= 1 ) then
!!$!             s = s + 1-grow(i)
!!$!          else
!!$!             s = s + 1
!!$          !          end if
!!$          !          s = s + abs(grow(i) - 1) !* raper(i)/raper(naper)
!!$!          if( d < grow(i) .and. grow(i) <= 1 ) then
!!$!             s = s + (1 - grow(i))
!!$!          else if( grow(i) > 1 ) then
!!$!             s = s + 100*(grow(i) - 1)
!!$!          else if( grow(i) <= d ) then
!!$!             s = s + 100*(d - grow(i))
!!$!          end if
!!$       end do
!!$    end do
!!$!    s = s / nstar
!!$
!!$!    mincom = s + lambda(1)*max(0.0,t) + lambda(2)*max(0.0,a)
!!$    mincom = s + abs(e) + a
!!$!    mincom = s !+ 0.01*abs(e)
!!$
!!$!    write(*,*) real(s),real(e),real(a)
!!$
!!$    if( a < 0 ) mincom = 1e3
!!$
!!$    s = 0
!!$    do i = 1,naper
!!$       call rinit(grows(:,i),x,d)
!!$       write(*,*) i,real(x),real(d)
!!$       s = s + d
!!$    end do
!!$    mincom = s + abs(e) + a
!stop 0
!  end function mincom

end module grow_model
