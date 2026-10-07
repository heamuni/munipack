!
!  transformation of an instrumental to standard photometry system
!
!  Copyright © 2013-5, 2017 F.Hroch (hroch@physics.muni.cz)
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
! Additional work:
!
!  * implement Pending function in analytical derivations, derfcn
!  * solve problem of error estimation
!  * provide calibration informations (star identifications,...)
!  * try Augmented Lagrangian method
!  * try non-triagonal matrix
!  * generalize on non-square transformation matrix
!

module fotran

  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter :: ci0 = 10**(-0.4*0.5)
  real(dbl), dimension(:,:), allocatable, private :: cts,dcts
  real(dbl), dimension(:,:), allocatable, private :: pht,dpht
  real(dbl), dimension(:), allocatable, private :: mad
  real(dbl), private :: lambda
  integer, private :: ncts, npht, ndat
  logical, private :: verbose

  private :: res,relerr,robfcn,medtri,derfcn,grouping

contains

  subroutine stdtra(pairs,xpht,xdpht,xcts,xdcts,tra,etra,verb)

    use oakleaf
    use robustplane

    integer, dimension(:,:), intent(in) :: pairs
    real(dbl), dimension(:,:), intent(in) :: xpht,xcts,xdpht,xdcts
    real(dbl), dimension(:,:), intent(out) :: tra,etra
    logical, intent(in) :: verb

    real(dbl),dimension(2) :: p,dp
    real(dbl),dimension(3) :: q,dq
    real(dbl),dimension(:), allocatable :: r,dr,u,du,v,dv
    real(dbl) :: sig
    integer :: i,j,k,l,n,m,ntriplets,ntwins,npairs,lref,info
    integer, dimension(size(pairs,1),2) :: twins
    integer, dimension(size(pairs,1),3) :: triplets
    integer, dimension(2) :: refpair

    npht = size(xpht,2)
    ncts = size(xcts,2)
    verbose = verb

    if( size(xpht,1) /= size(xcts,1) ) stop 'Dimensions pht /= cts.'
    if( size(xdpht,1) /= size(xpht,1) ) stop 'Dimensions pht /= dpht.'
    if( size(xdcts,1) /= size(xcts,1) ) stop 'Dimensions cts /= dcts.'

    ndat = size(xpht,1)
    npairs = size(pairs,1)

    if( .not. (ndat > 0) ) stop 'Data missing.'

    allocate(r(ndat),dr(ndat),pht(ndat,npht),dpht(ndat,npht), &
         cts(ndat,npht),dcts(ndat,npht),u(ndat),du(ndat),v(ndat),dv(ndat))
    pht = xpht
    dpht = xdpht
    cts = xcts
    dcts = xdcts

    if( verbose ) then
       write(*,*) 'Pairs indexes:'
       do k = 1,size(pairs,1)
          write(*,*) pairs(k,1),pairs(k,2)
       end do
    end if

    tra = 0.0_dbl
    etra = 0.0_dbl

    call grouping(pairs,ntriplets,triplets,ntwins,twins,lref,refpair)

    ! precise estimation from triplets
    do n = 1,ntriplets
       i = triplets(n,1)
       l = triplets(n,2)
       j = triplets(n,3)
       call ratio(cts(:,l),pht(:,l),dcts(:,l),dpht(:,l),r,dr)
       call ratio(pht(:,i),pht(:,l),dpht(:,i),dpht(:,l),u,du)
       call ratio(pht(:,j),pht(:,l),dpht(:,j),dpht(:,l),v,dv)
       u = u - ci0
       v = v - ci0
       q = (/ 1.0, 0.0, 0.0 /)
       dq = (/ 0.1, 0.01,0.01 /)
       call rplane(u,du,v,dv,r,dr,q,dq,sig,info)
       if( info == 0 ) then
          tra(l,l) = q(1) - ci0*(q(2) + q(3))
          tra(l,l-1) = q(2)
          tra(l,l+1) = q(3)
          etra(l,l-1:l+1) = dq
          if( verbose ) then
             write(*,*) 'triplet:',triplets(n,:)
             do k = 1,ndat
                write(*,'(5f10.5)') u(k),v(k),r(k),r(k) - (q(1) + q(2)*u(k) + q(3)*v(k))
             end do
             write(*,'(a,5f10.5)') 't-solution:',q
             write(*,'(a,6f10.5)') 't-errors:',dq,sig
          end if
       else
          stop 'No convergence occurred on a plane.'
       end if
    end do

    ! reference filter, if it is not already in some triplet
    if( .not. all( refpair == 0) ) then
       l = lref
       i = refpair(1)
       j = refpair(2)

       call ratio(pht(:,i),pht(:,j),dpht(:,i),dpht(:,j),u,du)
       call ratio(cts(:,l),pht(:,l),dcts(:,l),dpht(:,l),r,dr)
       u = u - ci0
       call rline(u,r,p(1),p(2),dp(1),dp(2),du,dr,sig)
       if( verbose ) then
          write(*,*) "Reference: ",i,j
          do k = 1,ndat
             write(*,'(5f10.5,1x,1pg0.3)') u(k),r(k),du(k),dr(k), &
                  r(k) - (p(1) + p(2)*u(k)),pht(k,l)
          end do
          write(*,'(a,5f10.5)') 'a-solution:',p,dp,sig
       end if
       tra(l,l) = p(1) - p(2)*ci0
       tra(l,i) = p(2)
       etra(l,l) = dp(1)
       etra(l,i) = dp(2)
    end if

    ! colour indexes
    l = pairs(1,1)
    do n = 1, ntwins
       i = twins(n,1)
       j = twins(n,2)

       call ratio(pht(:,i),pht(:,j),dpht(:,i),dpht(:,j),u,du)
       call ratio(cts(:,i),cts(:,j),dcts(:,i),dcts(:,j),r,dr)

       ! correction while previous line is a triplet
       m = 0
       do k = 1, ntriplets
          if( triplets(k,1) == i .and. triplets(k,2) == j ) then
             m = triplets(k,3)
             call ratio(pht(:,m),pht(:,j),dpht(:,m),dpht(:,j),v,dv)
          else if( triplets(k,2) == j .and. triplets(k,3) == i ) then
             m = triplets(k,1)
             call ratio(pht(:,m),pht(:,j),dpht(:,m),dpht(:,j),v,dv)
          end if
       end do

       u = u - ci0
       if( m == 0 ) then
          r = r * (1 + tra(j,i)/tra(j,j) * u + tra(j,i)/tra(j,j) * ci0)
       else
          r = r * (1 + tra(j,i)/tra(j,j) * u + tra(j,m)/tra(j,j) * (v - ci0) + &
               ci0*(tra(j,i) + tra(j,m))/tra(j,j))
       end if
       call rline(u,r,p(1),p(2),dp(1),dp(2),du,dr,sig)
       tra(i,i) = p(2) * tra(j,j)
       tra(i,j) = p(1) * tra(j,j) - tra(i,i)* ci0
       etra(i,i) = dp(2)
       etra(i,j) = dp(1)

       if( verbose ) then
          write(*,*) "Colour index: ",i,j
          call ratio(cts(:,i),cts(:,j),dcts(:,i),dcts(:,j),r,dr)
          do k = 1,ndat
             write(*,'(5f10.5,3(1x,1pg0.5))') u(k),r(k),du(k),dr(k), &
                  r(k) - tra(i,i)/tra(j,j)*(u(k) + tra(i,j)/tra(i,i) + ci0) / &
                  (1 + tra(j,i)/tra(j,j)*u(k) + tra(j,i)/tra(j,j)*ci0),pht(k,i),pht(k,j)
          end do
          write(*,'(a,5f10.5)') 'c-solution:',p,dp,sig
       end if

    end do

   deallocate(r,dr,pht,dpht,cts,dcts,u,du,v,dv)

  end subroutine stdtra


  subroutine fotra(t,pairs,cts,dcts,pht,dpht)

    integer, dimension(:,:), intent(in) :: pairs
    real(dbl), dimension(:,:), intent(in) :: t
    real(dbl), dimension(:), intent(in) :: cts,dcts
    real(dbl), dimension(:), intent(out) :: pht,dpht

    integer, dimension(size(pairs,1),2) :: twins
    integer, dimension(size(pairs,1),3) :: triplets
    integer, dimension(2) :: refpair
    real(dbl) :: ph,c
    integer:: i,j,k,l,n,m,ntriplets,ntwins,lref

    call grouping(pairs,ntriplets,triplets,ntwins,twins,lref,refpair)

    do n = 1,ntriplets
       i = triplets(n,1)
       l = triplets(n,2)
       j = triplets(n,3)
       if( cts(i) > 0 .and. cts(l) > 0 .and. cts(j) > 0 ) then
          ph = t(l,l) + ci0*(t(l,i) + t(l,j)) + &
               t(l,i)*(cts(i)/cts(l) - ci0) + t(l,j)*(cts(j)/cts(l) - ci0)
          pht(l) = ph * cts(l)
       else
          pht(l) = -1
       end if
       if( pht(l) > 0 .and. dcts(l) > 0 .and. cts(l) > 0 ) then
          dpht(l) = pht(l) * (dcts(l) / cts(l))
       else
          dpht(l) = -1
       end if
    end do

    ! reference
    if( .not. all( refpair == 0) ) then
       l = lref
       i = refpair(1)
       j = refpair(2)
       if( cts(i) > 0 .and. cts(j) > 0 ) then
          ph = t(l,l) + t(l,i) * ci0 + t(l,i)*(cts(i)/cts(j) - ci0)
          pht(l) = ph * cts(l)
          dpht(l) = pht(l) * (dcts(l) / cts(l))
       else
          pht(l) = -1
          dpht(l) = -1
       end if
    end if

    do n = 1,ntwins
       i = twins(n,1)
       j = twins(n,2)

       m = 0
       do k = 1, ntriplets
          if( triplets(k,1) == i .and. triplets(k,2) == j ) then
             m = triplets(k,3)
          else if( triplets(k,2) == j .and. triplets(k,3) == i ) then
             m = triplets(k,1)
          end if
       end do

       if( cts(i) > 0 .and. cts(j) > 0 .and. pht(j) > 0 ) then
          c = cts(i)/cts(j) - ci0
          if( m == 0 ) then
             ph = t(i,i) / t(j,j) * (c + t(i,j)/t(i,i) + ci0) / &
                  (1 + t(j,i)/t(j,j)*ci0 + t(j,i)/t(j,j)*c)
          else
             ph = t(i,i) / t(j,j) * (c + t(i,j)/t(i,i) + ci0) / &
                  (1 + t(j,i)/t(j,j) * c + t(j,m)/t(j,j) * (cts(m)/cts(j) - ci0) + &
                  ci0*(t(j,i) + t(j,m))/t(j,j))
          end if
          pht(i) = ph * pht(j)
          dpht(i) = pht(i) * (dcts(j) / cts(j))
       else
          pht(i) = -1
          dpht(i) = -1
       end if
    end do

!    return

!!$    n = size(ctph0)
!!$    allocate(tra(n,n),ctph(n))
!!$    call traex(pairs,tra0,ctph0,tra)
!!$    ctph = 1.0_dbl
!!$
!!$    npht = size(cts)
!!$    allocate(ctc(npht-1))
!!$    do i = 2,npht
!!$       k = pairs(i,1)
!!$       l = pairs(i,2)
!!$       if( cts(l) > epsilon(cts) .and. cts(k) > epsilon(cts) ) then
!!$          ctc(i-1) = (ctph(k)* cts(k)) / (ctph(l) * cts(l))
!!$       else
!!$          ctc(i-1) = 0
!!$       end if
!!$    end do
!!$
!!$    k = pairs(1,1)
!!$    l = pairs(1,2)
!!$    if( cts(k) > epsilon(cts) ) then
!!$       pht(k) = ctph(k) * cts(k) * (tra(1,1) + tra(1,2)*ctc(1))
!!$       dpht(k) = (dcts(k)/cts(k))* pht(k)
!!$    else
!!$       pht = -1
!!$       dpht = -1
!!$       goto 666
!!$    end if
!!$    down = .true.
!!$    do i = 2,npht
!!$       k = pairs(i,1)
!!$       l = pairs(i,2)
!!$       if( k == pairs(1,1) ) down = .false.
!!$       if( down ) then
!!$          if( pht(l) > epsilon(pht) .and. cts(l) > epsilon(cts) ) then
!!$!             write(*,*) pht(l),tra(i,1),tra(i,2),ctc(i-1)
!!$             pht(k)  = pht(l) * (tra(i,1) + tra(i,2)*ctc(i-1))
!!$             dpht(k) = (dcts(l) / cts(l)) * pht(k)
!!$          else
!!$             pht(k) = -1
!!$             dpht(k) = -1
!!$          end if
!!$       else
!!$          if( pht(k) > epsilon(pht) .and. cts(k) > epsilon(cts) .and. &
!!$               abs(tra(i,1) + tra(i,2)*ctc(i-1)) > epsilon(ctc) ) then
!!$             pht(l)  = pht(k) / (tra(i,1) + tra(i,2)*ctc(i-1))
!!$             dpht(l) = (dcts(k) / cts(k)) * pht(l)
!!$          else
!!$             pht(l) = -1
!!$             dpht(l) = -1
!!$          end if
!!$       end if
!!$    end do
!!$
!!$    where( pht < 0 )
!!$       pht = -1
!!$       dpht = -1
!!$    end where
!!$
!!$666 continue
!!$    deallocate(ctc,ctph,tra)

  end subroutine fotra




  subroutine ratio(x,y,dx,dy,r,dr)

    real(dbl), dimension(:), intent(in) :: x,y,dx,dy
    real(dbl), dimension(:), intent(out) :: r,dr

    where( abs(x) > epsilon(dx) .and. abs(y) > epsilon(dy))
       r = x / y
       dr = abs(r)*sqrt((dx/x)**2 + (dy/y)**2)
    elsewhere
       r = 0
       dr = huge(dr)
    end where

  end subroutine ratio


  subroutine grouping(pairs,ntriplets,triplets,ntwins,twins,lref,refpair)

    integer, dimension(:,:), intent(in) :: pairs
    integer, dimension(:,:), intent(out) :: triplets, twins
    integer, dimension(:), intent(out) :: refpair
    integer, intent(out) :: ntriplets, ntwins, lref

    integer :: i,j,n,npairs
    logical :: found, reference

    npairs = size(pairs,1)

    ! form triplets
    ntriplets = 0
    triplets = 0
    do i = 2, npairs
       triplets(ntriplets+1,1:2) = pairs(i,:)
       do j = 2, npairs
          if( triplets(ntriplets+1,2) == pairs(j,1) ) then
             ntriplets = ntriplets + 1
             triplets(ntriplets,3) = pairs(j,2)
          end if
       end do
    end do

    ! form twins
    ntwins = 0
    twins = 0
    do i = 2, npairs
       ntwins = ntwins + 1
       twins(ntwins,:) = pairs(i,:)
       found = .false.
       do j = 2, npairs
          if( pairs(j,1) == twins(ntwins,2) ) found = .true.
       end do
       if( found ) ntwins = ntwins - 1
    end do
    do i = 2, npairs
       ntwins = ntwins + 1
       twins(ntwins,:) = pairs(i,:)
       found = .false.
       do j = 2, npairs
          if( pairs(j,2) == twins(ntwins,1) ) found = .true.
       end do
       do j = 1, ntwins-1
          if( all(twins(j,:) == twins(ntwins,:)) ) found = .true.
       end do
       if( found ) ntwins = ntwins - 1
    end do

    ! reference filter, if it is not already in some triplet
    reference = .true.
    do i = 1,ntriplets
       if( any(pairs(1,1) == triplets(i,:)) ) reference = .false.
    end do
    ! reference = ntriplets > 0

    refpair = 0
    lref = pairs(1,1)
    if( reference ) then
       do n = 2, npairs
          if( pairs(n,1) < lref ) then
             i = pairs(n,1)
             j = pairs(n,2)
          else if( pairs(n,2) > lref ) then
             i = pairs(n,2)
             j = pairs(n,1)
          end if
          if( j == lref  )then
             refpair = (/i,j/)
          end if
          if( i == lref .and. all(refpair == 0) ) then
             refpair = (/i,j/)
          end if
       end do

    end if

    do n = 1, ntwins
       if( twins(n,1) < lref ) then
          i = twins(n,1)
          j = twins(n,2)
       else if( twins(n,2) > lref ) then
          i = twins(n,2)
          j = twins(n,1)
       end if
       twins(n,:) = (/i,j/)
    end do

  end subroutine grouping


  subroutine tratri(xpht,xdpht,xcts,xdcts,tratab,verb)

    use NelderMead
    use minpack
    use minpacks
    use oakleaf, only: qmean, huber, dhuber

    integer :: npar
    real(dbl), dimension(:,:), target, intent(in) :: xpht,xcts,xdpht,xdcts
    real(dbl), dimension(:,:), allocatable, intent(out) :: tratab
    logical, intent(in) :: verb

    real(dbl),dimension(:), allocatable :: p,p0,dp
    real(dbl),dimension(:,:), allocatable :: t,dt,r,relp,fjac,cov
    real(dbl) :: s0,reqmin,sum1,sum2,w,med, rms
    integer :: icount, numres, ifault, info, nprint, it, i, j, n
    character(len=80) :: tfmt,rfmt,resfmt,trfmt


    ! BIG WARNING!
!    pht => xpht
!    cts => xcts
!    dpht => xdpht
!    dcts => xdcts
    npht = size(pht,2)
    ncts = size(cts,2)
    verbose = verb

    if( size(pht,1) /= size(cts,1) ) stop 'Dimensions pht /= cts.'
    if( size(dpht,1) /= size(pht,1) ) stop 'Dimensions pht /= dpht.'
    if( size(dcts,1) /= size(cts,1) ) stop 'Dimensions cts /= dcts.'

    ndat = size(pht,1)

    if( .not. (ndat > 0) ) stop 'Data missing.'

    if( verbose ) then
       write(*,*) "=> Initial data:"
       write(tfmt,'(a,i0,a,i0,a)') '(i5,',npht,'g12.5,',npht,'g12.5)'
       do j = 1,size(pht,1)
          write(*,tfmt) j,pht(j,:),cts(j,:)
       end do
    end if

    write(tfmt,'(a,i0,a)') '(',npht,'f15.6)'
    write(rfmt,'(a,i0,a,i0,a,i0,a)') '(',2*(npht-1),'g12.5,',npht,'f9.5,',npht,'f10.5)'
!    write(*,*) trim(rfmt)

    npar = npht*ncts
    allocate(t(ncts,npht),dt(ncts,npht),r(ndat,ncts),relp(ndat,ncts))

    r = pht / cts

    if( verbose ) then
       write(*,*) "=> Ratios:"
       write(trfmt,'(a,i0,a)') '(i5,',npht,'g12.5)'
       do j = 1,size(pht,1)
          write(*,trfmt) j,r(j,:)
       end do
    end if

    t = 0
    dt = epsilon(dt)
    dt = 0
    do i = 1,ncts
       call qmean(r(:,i),t(i,i),dt(i,i))
       dt(i,i) = dt(i,i) / sqrt(real(ndat))
    end do

    if( verbose ) then
       write(*,*) "=> Initial estimation of diagonal:"
       do i = 1,ncts
          write(*,'(i5,1x,g0.5,a,g0.5)') i,t(i,i)," +- ",dt(i,i)
       end do
    end if

!    do i = 1,ncts
!       do j = 1,npht
!          if( i /= j ) t(i,j) = -(dt(i,i)/t(i,i) + dt(j,j)/t(j,j))*(t(i,i) + t(j,j))/4
!       end do
!    end do


    ! bootstrap for robust method
    allocate(p(npar),p0(npar),dp(npar))
    p = pack(t,.true.)

    if( verbose ) then
       write(*,*) "=> Estimation minimizing of absolute deviations: ndat=",ndat
       write(*,*) " # status      mean deviation    lambda   evaluations"
    end if

    lambda = sum(abs(p))/npar
    reqmin = epsilon(reqmin)
    s0 = huge(s0)
    do it = 1,1!00

       n = 0
       do i = 1,ncts
          do j = 1,npht
             n = n + 1
             dp(n) = max(dt(i,i),dt(j,j),epsilon(1.0))*(-1.0)**(i+j)
          end do
       end do

       p0 = p
       call nelmin(medtri,npar,p0,p,rms,reqmin,dp,1,1000000,icount,numres,ifault)
       t = reshape(p,(/ncts,npht/))

       if( verbose ) then
          write(*,'(i3,i5,g15.5,g10.3,i10)') it,ifault,rms,lambda,icount
          write(*,tfmt) (t(i,:),i=1,ncts)
       end if

       if( abs(rms - s0) < 0.001 ) exit

       lambda = 10*lambda
       s0 = rms

    enddo

    allocate(mad(ncts))

    call res(t,r)
    r = abs(r)
    do j = 1,size(r,2)
       call qmean(r(:,i),mad(j))
       mad(j) = mad(j) / sqrt(real(ndat))
    end do

    if( verbose ) then
       write(*,*) "=> Estimation of scale:"
       write(*,tfmt) mad
    end if

    if( verbose ) then
       write(*,*) "   Residuals       .... relative precision  "
       write(resfmt,'(a,i0,a,i0,a)') '(i5,',npht,'g13.5,',npht,'f10.5)'
       call relerr(t,relp)
       do i = 1,size(r,1)
          write(*,resfmt) i,r(i,:),relp(i,:)
       end do
    end if

    allocate(fjac(npar,npar),cov(npar,npar))

    if( verbose ) then
       nprint = 1
    else
       nprint = 0
    end if

    lambda = sum(abs(p))
    do j = 1,7
       call lmdif2(robfcn,p,epsilon(p),nprint,info)
       lambda = 10*lambda
    end do
    call lmder2(derfcn,p,epsilon(p),nprint,info)

    t = reshape(p,(/ncts,npht/))

    ! residual sum
    ! compute jac !!!
    call qrinv(fjac,cov)

    call res(t,r)
    do j = 1,ncts

       sum1 = 0
       sum2 = 0
       do i = 1, ndat
          w = r(i,j)/mad(j)
          sum1 = sum1 + huber(w)**2
          sum2 = sum2 + dhuber(w)
       end do
       s0 = mad(j)**2*sum1/sum2**2*ndat**2
       rms = sqrt(s0 / (ndat - npar))

       n = 0
       do i = 1, npht
          n = n + 1
          dt(j,i) = sqrt(s0*cov(n,n)/(ndat - npar))
       end do

    enddo

    if( verbose ) then
       write(*,*) "   Residuals       .... relative precision  "
       write(resfmt,'(a,i0,a,i0,a)') '(i5,',npht,'g13.5,',npht,'f10.5)'
       call relerr(t,relp)
       do i = 1,size(r,1)
          write(*,resfmt) i,r(i,:),relp(i,:)
       end do
    end if

    if( verbose ) then
       write(*,*) "=> Robust estimation:"
       write(*,*) ' info=',info
       write(*,*) " Laplace's multiplicator=",lambda
       write(*,*) " Jacobian in minimum:"
       write(*,*) " RMS=",rms

       write(*,*) " Matrix of Transformations:"
       write(*,tfmt) (t(i,:),i=1,ncts)

       write(*,*) " Errors of Matrix of Transformation:"
       write(*,tfmt) (dt(i,:),i=1,ncts)

       if( size(pht,2) > 1 ) then
          write(*,*) " Band ratios (catalogue, instrumental)  Ratios  Relative precision:"
          call relerr(t,relp)
          do i = 1,ndat
             write(*,rfmt) (pht(i,j-1)/pht(i,j),j=2,size(pht,2)), &
                  (cts(i,j-1)/cts(i,j),j=2,size(pht,2)),pht(i,:)/cts(i,:),relp(i,:)
          end do

          do i = 1,ndat
             write(*,*) (pht(i,j-1)/pht(i,j),j=2,size(pht,2)), &
                  (cts(i,j-1)/cts(i,j),j=2,size(pht,2)),&
                  pht(i,:),cts(i,:),cts(i,:)/pht(i,:)
          end do
       end if

    end if

    allocate(tratab(ncts,npht))
    tratab = t

    deallocate(t,dt,p,p0,dp,r,relp,mad,fjac,cov)

  end subroutine tratri

  subroutine res(t,r)

    real(dbl), dimension(:,:), intent(in) :: t
    real(dbl), dimension(:,:), intent(out) :: r
    real(dbl) :: x,d
    integer :: i,n

    do n = 1, ndat
       do i = 1,npht
          x = sum(t(i,:)*cts(n,:))
          d = sqrt(dpht(n,i)**2 + sum(t(i,:)**2 * dcts(n,:)**2))
          r(n,i) = (pht(n,i) - x)/d
       end do
    end do

  end subroutine res

  subroutine relerr(t,r)

    real(dbl), dimension(:,:), intent(in) :: t
    real(dbl), dimension(:,:), intent(out) :: r
    real(dbl) :: x
    integer :: i,n

    do n = 1, ndat
       do i = 1,npht
          x = sum(t(i,:)*cts(n,:))
          r(n,i) = (pht(n,i) - x)/((pht(n,i) + x)/2)
       end do
    end do

  end subroutine relerr

  subroutine robfcn(m,np,p,fvec,iflag)

    use oakleaf

    integer, intent(in) :: m,np
    integer, intent(inout) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(out) :: fvec
    real(dbl), dimension(:,:), allocatable :: t,r,fv
    real(dbl) :: rp,ss,x,d
    integer :: npar,n,i,j

!    if( iflag == 0 .and. verbose ) write(*,*) '#robfcn:',real(p)

    npar = ncts*npht
    allocate(t(ncts,npht),fv(ncts,npht),r(ndat,ncts))
    fv = 0.0_dbl
    t = reshape(p,(/ncts,npht/))

    call res(t,r)

    do n = 1, ndat

       do i = 1,npht
          x = sum(t(i,:)*cts(n,:))
          d = sqrt(sum(t(i,:)**2 * dcts(n,:)**2) + dpht(n,i)**2)
          rp = huber(r(n,i)/mad(i))
          fv(i,:) = fv(i,:) - rp*(cts(n,i)/(d*mad(i)))
       end do
    end do

    ss = 0.0_dbl
    do j = 1,ncts
       do i = 1,npht
          if( abs(i-j) > 1 ) then
             fv(i,j) = fv(i,j) + lambda*t(i,j)**2
          end if
       end do
    end do

    fvec = pack(fv,.true.)

    deallocate(t,fv,r)

  end subroutine robfcn


  subroutine derfcn(m,np,p,fvec,fjac,ldfjac,iflag)

    use oakleaf

    integer, intent(in) :: m,np,ldfjac
    integer, intent(inout) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(out) :: fvec
    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
    real(dbl), dimension(:,:), allocatable :: t,r,fv
    real(dbl) :: rp,x,d,u
    integer :: npar,n,i,j

    npar = ncts*npht
    allocate(t(ncts,npht),r(ndat,ncts))
    t = reshape(p,(/ncts,npht/))
    call res(t,r)

    if( iflag == 1 ) then

       allocate(fv(ncts,npht))
       fv = 0

       do n = 1, ndat

          do i = 1,npht
             x = sum(t(i,:)*cts(n,:))
             d = sqrt(sum(t(i,:)**2 * dcts(n,:)**2) + dpht(n,i)**2)
             rp = huber(r(n,i)/mad(i))
             fv(i,:) = fv(i,:) - rp*(cts(n,i)/(d*mad(i)))
          end do
       end do

       fvec = pack(fv,.true.)
       deallocate(fv)

    else if( iflag == 2 ) then

       fjac = 0
       do n = 1, ndat

          do i = 1,ncts
             x = sum(t(i,:)*cts(n,:))
             d = sqrt(sum(t(i,:)**2 * dcts(n,:)**2) + dpht(n,i)**2)
             rp = dhuber(r(n,i)/mad(i))
             u = rp/(d*mad(i))**2
             do j = 1,npht
                fjac(i,j) = fjac(i,j) + u*cts(n,i)*cts(n,j)
             end do
          end do
       end do


    end if

    deallocate(t,r)

  end subroutine derfcn


  function medtri(p) result(s)

    real(dbl), dimension(:), intent(in) :: p
    real(dbl), dimension(:,:), allocatable :: t,r
    real(dbl) :: s,ss
    integer :: i,j

    allocate(r(ndat,ncts),t(ncts,npht))
    t = reshape(p,(/ncts,npht/))

    call res(t,r)
    s = sum(abs(r))/size(r)

    ss = 0
    do i = 1,ncts
       do j = 1,npht
          if( abs(i-j) > 1 ) then
             ss = ss + abs(t(i,j))
          end if
       end do
    end do

    s = s + lambda *ss

    deallocate(r,t)

  end function medtri


  subroutine trainv(pairs,tra,tra1)

    use minpacks

    integer, dimension(:,:), intent(in) :: pairs
    real(dbl), dimension(:,:), intent(in) :: tra
    real(dbl), dimension(:,:), intent(out) :: tra1
    real(dbl), dimension(size(tra,1),size(tra,1)) :: a,b
    real(dbl) :: det
    integer :: n

!    write(*,*) tra(1,:)
!    write(*,*) tra(2,:)


    call traconv(pairs,tra,a)

!    write(*,*) a(1,:)
!    write(*,*) a(2,:)

    do n = 1,size(tra,1)
       write(*,'(10f10.5)') a(n,:)
    end do

    n = size(a,1)
    if( n == 1 ) then
       b(1,1) = 1.0_dbl / a(1,1)
    else if( n == -2 ) then

       det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
       b(1,:) = (/  a(2,2),-a(1,2) /) / det
       b(2,:) = (/ -a(2,1), a(1,1) /) / det

    else

       call qrinv(a,b)

    end if

    do n = 1,size(tra,1)
       write(*,'(10f10.5)') b(n,:)
    end do


!    write(*,*) b(1,:)
!    write(*,*) b(2,:)

    call traiconv(pairs,b,tra1)

!    write(*,*) tra1(1,:)
!    write(*,*) tra1(2,:)

    do n = 1,size(tra1,1)
       write(*,'(10f10.5)') tra1(n,1:2)
    end do

!    stop

  end subroutine trainv

  subroutine traconv(pairs,tra,a)

    integer, dimension(:,:), intent(in) :: pairs
    real(dbl), dimension(:,:), intent(in) :: tra
    real(dbl), dimension(:,:), intent(out) :: a
    integer :: n,i,k,l
    logical :: down

    n = size(tra,1)
    if( n == 1 ) then
       a(1,1) = tra(1,1)
    else if ( n == -2 ) then
       a(1,1) = tra(1,1)
       a(1,2) = tra(1,2)
       a(2,1) = tra(2,1)*a(1,2)
       a(2,2) = a(1,1)*tra(2,2) + tra(2,1)*a(1,2)
    else
       a = 0.0_dbl
       k = pairs(1,1)
       l = pairs(2,1)
       a(k,k) = tra(1,1)
       a(k,l) = tra(1,2)
       down = .true.
       do i = 2, n
          k = pairs(i,1)
          l = pairs(i,2)
          if( k == pairs(1,1) ) down = .false.
!          write(*,*) k,l,down
          if( down ) then
             a(k,l) = tra(i,1)*a(l,l)
             a(k,k) = a(l,l)*tra(i,2) + a(k,l)*a(l,k)/a(l,l)
          else
             a(k,l) = tra(i,1)*a(k,k)
             a(l,l) = a(k,k)*tra(i,2) + a(k,l)*a(l,k)/a(k,k)
!             a(k,i-1) = tra(i,1)*a(i-1,i-1)
!             a(k,i) = a(i-1,i-1)*tra(i,2) + tra(i,1)*a(i-1,i)
          end if
       end do
    end if

  end subroutine traconv

  subroutine traiconv(pairs,b,tra)

    integer, dimension(:,:), intent(in) :: pairs
    real(dbl), dimension(:,:), intent(in) :: b
    real(dbl), dimension(:,:), intent(out) :: tra
    integer :: n,i,k,l
    logical :: down

    n = size(b,1)
    if( n == 1 ) then
       tra(1,1) = b(1,1)
    else if( n == -2 ) then
       tra(1,1) = b(1,1)
       tra(1,2) = b(1,2)
       tra(2,1) = b(2,1)*tra(1,1)
       tra(2,2) = b(2,2)/tra(1,1) - b(2,1)*tra(1,2)/tra(1,1)**2
    else
       tra = 0.0_dbl
!       tra(1,1:2) = b(1,1:2)
       k = pairs(1,1)
       l = pairs(2,1)
       tra(1,1) = b(k,k)
       tra(1,2) = b(k,l)
       down = .true.
       do i = 2,n
          k = pairs(i,1)
          l = pairs(i,2)
          if( k == pairs(1,1) ) down = .false.
          if( down ) then
             tra(i,1) = b(k,l)/b(l,l)
             tra(i,2) = b(k,k)/b(l,l) - b(k,l)*b(l,k)/b(l,l)**2
          else
             tra(i,1) = b(k,l) / b(k,k)
             tra(i,2) = b(l,l) / b(k,k) - b(k,l)*b(l,k)/b(k,k)**2
          end if
       end do
    end if

  end subroutine traiconv

  subroutine traex(pairs,tratab,t,tra)

    integer, dimension(:,:), intent(in) :: pairs
    real(dbl), dimension(:,:), intent(in) :: tratab
    real(dbl), dimension(:), intent(in) :: t
    real(dbl), dimension(:,:), intent(out) :: tra
    real(dbl), dimension(size(tra,1),size(tra,1)) :: tra1,tra2
    integer :: i

    tra2 = 0.0_dbl
    forall(i=1:size(t))
       tra2(i,i) = t(pairs(i,1))
    end forall
    call traconv(pairs,tratab,tra1)

!    write(*,*) tra1(1,1:2)
!    write(*,*) tra1(2,1:2)

    tra1 = matmul(tra1,tra2)
!    write(*,*) tra1(1,1:2)
!    write(*,*) tra1(2,1:2)
    call traiconv(pairs,tra1,tra)

  end subroutine traex


end module fotran
