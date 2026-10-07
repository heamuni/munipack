!
!  Robust estimation of mean photon rates
!
!  Copyright © 2020 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of OakLeaf.
!
!  OakLeaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  OakLeaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!

module photonrate_init

  use photonrate_common
  use iso_fortran_env

  implicit none
  private

  public :: phinit

  logical :: debug_graph = .true.

  real(REAL64), dimension(:), pointer :: cts, pht
  real(REAL64), dimension(3), parameter :: q = [1, 2, 3] / real(4.0,REAL64)
  real(REAL64), dimension(3) :: tqn
  real(REAL64) :: phi, sx, sy, tauback

  integer, dimension(:), pointer :: freq
  real(REAL64), dimension(:), pointer :: kmean

contains

    subroutine phinit(y,x,kmeans,freqs,sxs,sys,t,tau,dtau,rho,reli,identical)

    ! https://en.wikipedia.org/wiki/Ratio_distribution

    ! This is the very first estimation of the mean rate.

    use fmm
    use medians
    use quantilmean

    real(REAL64), dimension(:), intent(in), target :: y,x,kmeans, t
    integer, dimension(:), intent(in), target :: freqs
    real(REAL64), intent(out) :: tau,dtau,rho
    real(REAL64), intent(in) :: sxs,sys
    logical, intent(out) :: reli, identical

    real(REAL64), dimension(:), allocatable :: xcdf, ycdf
    real(REAL64) :: d,s,tol,rhomin,rhomax,mad, mu, rho0, hromin
    integer :: i, n, stat
    character(len=80) :: errmsg
    character(len=9) :: rhotxt

    reli = .false.
    identical = .false.
    phi = PI / 4
    sx = sxs
    sy = sys

    cts => y
    pht => x
    freq => freqs
    kmean => kmeans

    n = size(t)
    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,'(2a)') 'Error in phinit: ',trim(errmsg)
       return
    end if

    ! empirical CDF
    call ecdf(t,xcdf,ycdf)

    ! The quantile estimate is effective assymptotically for large
    ! data sets. There the limit can be optimal for less than
    ! 10% contamination: in the case, the CDF diverts above/under
    ! to quantiles 0.1 and 0.9. A considerable contamination
    ! diverts CDF more close to the 1/4  and 3/4 quantiles.

    do i = 1, 3
       tqn(i) = quantile(q(i),xcdf,ycdf)
    end do
    mad = (tqn(3) - tqn(1)) / 2

    if( debug_iter ) then
       write(error_unit,*) 'phinit() Quantiles: Q(1/4)=',real(tqn(1)), &
            'Q(1/2)=',real(tqn(2)),'Q(3/4)=',real(tqn(3))
    end if
    deallocate(xcdf,ycdf)

    if( n < 50 ) then ! update for low amount of data

       ! Important! The median is more robust than the quantiles for
       ! short data sets; see the example: -10, 1,2,3,4.

       mad = median(abs(tqn(2) - t))
       tqn(1) = tqn(2) - mad
       tqn(3) = tqn(2) + mad

    end if

    if( debug_iter )  write(error_unit,*) 'phinit() MAD:',real(mad)

    ! identical data
    identical = mad < 10*epsilon(mad)

    if( identical ) then
       if( debug_iter ) write(error_unit,*) 'Phinit detected identical data.'
       reli = .true.
       return
    end if

    ! estimate of mean
    tau = tqn(2)

    ! estimate of standard deviation
    s = mad / Qdev

    ! estimate of standard error
    dtau = s / sqrt(real(n))

    ! estimate of rho
!    d = 2 * s * qmedian(kmean)
!    tol = d / 100
!    write(*,*) 'qmedian:',qmedian(kmean),d


    ! Realiable results are inside the range 0 <= rho < rhomax,
    ! ones aren't so close to the interval the upper endpoint.

!    if( debug_iter ) write(error_unit,*) 'phinit(): rho: ',real(rho)

    ! estimate of rho
    mu = qmedian(kmean)
    d = mu*(mu*s**2 - tau*(tau + 1))
    d = (mu*s)**2
    if( debug_iter ) write(error_unit,*) 'Mean mu:',real(mu),'   d=',real(d)
    if( d > 0 ) then
       rho0 = sqrt(d)
    else
       rho0 = sqrt(mu) ! order estimate, rho value is similar to mu
    end if
    tol = rho0 / 30
    if( debug_iter ) write(error_unit,*) &
         'Initial estimate of rho:',real(rho0),'  tol:',nint(tol)

    rhotxt = ''
    tauback = tau
    rhomin = 0
    hromin = hrofun(rhomin)
    if( hromin > 0 ) then
       do i = 1, 3
          rhomax = 2**i * rho0
          if( debug_iter ) write(error_unit,'(a,i1,a,2(2x,g0),a)') &
               ' Rho bracketing #',i,': 0 <= rho <=',real(rhomax),trim(rhotxt)
          if( hromin*hrofun(rhomax) < 0 ) then
          !       rho = fmin(rhomin,rhomax,hrofun,tol)
             rho = zeroin(rhomin,rhomax,hrofun,tol)
             reli = rho < rhomax - tol !.and. abs(hrofun(rho)) / Qdev < 1e-2
             if( reli ) exit
          end if
          rhotxt = ' (update)'
       end do
    else
       write(error_unit,*) 'Rho(min) < 0, skipping update ...'
    end if

!!$    if( .not. reli ) then
!!$       do i = 1, 3
!!$          rhomin = 0
!!$          rhomax = 2**i * rho0
!!$          if( debug_iter ) then
!!$             write(error_unit,'(a,i1,a,2(2x,g0))') &
!!$                  'Vrho limits min,max: ',i,'-sigma, ',real(rhomin),real(rhomax)
!!$             if( i > 1 ) write(error_unit,'(a,2(2x,g0))') &
!!$                  'updated min, max: ',real(rhomin),real(rhomax)
!!$          end if
!!$          rho = fmin(rhomin,rhomax,Vrho,tol)
!!$          reli = rho < rhomax - tol
!!$          if( reli ) exit
!!$       end do
!!$    end if

    if( .not. reli ) then
       rho = rho0
       reli = .true.
    end if

    ! Realiable results are inside the range 0 <= rho < rhomax,
    ! ones aren't so close to the interval the upper endpoint.

!    if( debug_iter ) write(error_unit,*) 'phinit(): rho: ',real(rho), reli

!    reli = .true.
!    return


    if( debug_graph ) then
       ! graph of hrofun()
       block
         real(REAL64) :: rho  ! it shadows of the main one
         real(REAL64) :: rhomax

         rhomax = 3*rho0
         open(1,file='/tmp/rhomin')
         do i = 0, 100
            rho = i*rhomax/100.0
            write(1,*) rho, hrofun(rho) / Qdev, Vrho(rho)
         end do
         close(1)

!         tauback = 2 !!!!!!!!
         rhomax = 3*rho0
         open(1,file='/tmp/rhomin1')
         do i = 0, 100
            rho = i*rhomax/100.0
            write(1,*) rho, hrofun(rho) / Qdev, Vrho(rho)
         end do
         close(1)
       end block
    end if

  end subroutine phinit

  real(REAL64) function hrofun(rho)

    use medians

    real(REAL64), intent(in) :: rho
    real(REAL64) :: mad, tau, rho2, tt2

    tau = tauback
    tt2 = (tau + 1)/2
    rho2 = rho**2

    mad = sqrt(2.0)*median(abs(cts - tau*pht)/sqrt((cts+tau*pht)*tt2 + rho2))
    hrofun = abs(mad - Qdev)
    hrofun = mad - 1

    if( debug_iter .and. .not. debug_graph ) &
         write(error_unit,*) real(rho),real(mad),real(hrofun)

  end function hrofun


  real(REAL64) function Vrho(rho)

    ! This function returns maximum of the distance between the empirical
    ! and the estimated CDF. The distance is determined at 1/4 and 3/4
    ! quantiles. The horizonal coordinates is determined by the empirical
    ! CDF. The verical distance is best fitted by rho.
    !   The function max(d1,d2) is more efficient than d1 + d2, because
    ! there is only a single global (sharp) minimum -- the sum has more
    ! complicated profile (= more machine operations).
    !   In principle, the method is similar to K-S test: it estimates
    ! the rho parameter by minimising distance between two CDFs.

    use photonrate_normal

    real(REAL64), intent(in) :: rho
    real(REAL64), dimension(2) :: u,v

    u = [ tqn(1), tqn(3) ]
    v = mcdf(kmean,freq,sx,sy,tqn(2),rho,phi,u)
    Vrho = max(abs(v(1) - q(1)), abs(v(2) - q(3)))

    if( debug_iter .and. .not. debug_graph ) &
         write(error_unit,*) 'Vrho:',real(rho),real(Vrho)

  end function Vrho

end module photonrate_init
