!
!  Robust estimation of mean photon rates
!
!  Copyright © 2019-20 F.Hroch (hroch@physics.muni.cz)
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

!-------------------------------------------------------------------------

! This module implements both robust and least squares estimators of
! the ratio between two data sets, the sets of counts and photons.
!
! This is a common astronomical photometry task: to calibrate counts, like
! the integrated stars on CCD frames, against to some expected photons,
! derived by magnitudes of catalogue stars.  The problem is more general
! and can be applied on other various situations when a mixed ratio of
! Normal distributions appears.
!
!
! Data model
!
! Data are modelled by the ratio of sum of two distributions: the main one
! is the asymptotic approximation of Poisson distribution [1] for large rates:
!
!    Po(<n>) ~> N(<n>,sqrt(<n>)),
!
!  and Normal [2] distribution with the zero mean N(0,sigma).
!
! The mean ratio of photons over counts:
!
!           photons
!    tau =  -------
!           counts
!
! is determined by fitting of all parameters (tau, sigpht, sigcts)
! of ratio of distributions
!
!   Po(<n>) + N(0,sigpht)                     N(<n>,sqrt(<n>)) + N(0,sigpht)
!   ---------------------   approximated by   ------------------------------
!   Po(<c>) + N(0,sigcts)                     N(<c>,sqrt(<c>)) + N(0,sigcts)
!
! The exact analytic form of the distribution N(theta1,sig1)/N(theta2,sig2)
! has been computed by Hynkley [3]. This particular implementation is
! suitable for observed photon rates.
!
! Assumptions:
!
! * The approximation Po(n) = N(n,sqrt(n)) for the input data should
!   be valid. No verification is performed. The precision, for n > 666,
!   is high, while n > 100 gives solid results.
! * The final value of the ration should be tau > 1, or at least over 1/10.
!   That's reason for arranging of the ratio in the above form, and
!   no reciprocal one.
! * Additional noise must have Normal distribution with its zero mean,
!   and the Po(n) component should dominate over N(0,sig).
!
! If the Normal component does not dominate over Poisson component,
! this particular approximation is no more valid. A core distribution
! will close to Cauchy distribution [4]. The applicability on wide range
! of fluxes is impossible without the case.
!
! General advice for the usage in photometry calibration: selects
! as much as possible of brighter stars; do not include those stars on
! limit.
!
! Notes
!
! There is a principal difference of the mean ratio of the distributions
! (tau parameter), and the mean of ratios t(i) = photons(i) / counts(i).
! Their mean values E(t) can not be compared in any way. The ratio of mean
! values E(photons) / E(counts) is a good approximation, if all data has
! the same large fluxes, and no contamination is included.
!
! Both photons, and counts, should be real numbers, as product of previous
! computations; all values should be > 0 (zeros are not allowed).
!
! Recommend minimum of data is 5 (3 fitting parameters + 1 free bound
! + 1 safe point); the performance of the algorithm increases for
! large data sets.
!
!
! Parameters
!
! There are three internal fitted parameters (tau, rho, phi) for the whole
! data set. Moreover, data set is sliced on per-intensity basis: the data
! are gruped to clusters with the similar intensity by its statistical
! deviations.
!
! The rho is an effective magnitude of non-Poisson noise components.
! The phi is an effective relative angular measure of the noises component.
! The estimated values are
!
!   sigcts ~= rho*cos(phi),    sigpht = (rho/tau)*sin(phi).
!
! The parameter phi has very weak dependence on the determination of results;
! estimates of both sigcts, sigpht can be rough, especially for a short
! data sets of noisily data.
!
!
! The subroutine
!
! phrate should be called as:
!
!  call phrate(counts,photons, tau,tauerr,sigcts,sigpht,
!              robust, reliable, flag, verbose)
!
! On input:
!   counts - array of measured values
!   photons - array of reference values
!   ctsfac (optional) - the factor to reduce the photon noise of counts,
!                        commonly number of averadged frames
!   phtfac (optional) - the factor to reduce the photon noise of photons,
!                        commonly number of averadged frames
!   robustness (optional) - use/do not use robust method (default .true.)
!   clustering (optional) - do clustering of data (default .false.)
!   verbose(optional) - print an additional processing info
!
! On output are estimated:
!   tau - the mean rate of the distribution
!   tauerr (optional) - the standard error of tau
!   sigcts (optional) - estimate of stdsig for the noise of counts
!   sigpht (optional) - estimate of stdsig for the noise of photons
!   report (optional) - a data structure with a detailed processing report
!   reliable(optional)- indicates reliability of results
!   flag(optional)    - detailed info about reliability of results:
!            0 == success, results are both reliable and accurate
!            1 == consider a rough estimate due convergence fail, try verbose
!            2 == elemental estimates, a low amount of datapoints available
!            4 == any inconsistent input, try verbose
!            5 == failed to allocate memory
!            6 == identical data on input
!
! The results should be interpreted: the true value T of tau of the sample
! photons / counts can be, with 68% probability, found inside the interval
!
!         tau - tauerr  <  T  <  tau + tauerr.
!
!  [1] https://en.wikipedia.org/wiki/Poisson_distribution
!  [2] https://en.wikipedia.org/wiki/Normal_distribution
!  [3] https://doi.org/10.1093/biomet/56.3.635
!  [4] https://en.wikipedia.org/wiki/Ratio_distribution#Derivation


module photonrate

  use photonrate_common
  use iso_fortran_env

  implicit none
  private

  public :: phrate

!  type phrate_report
!     real(REAL64) :: tau, tau_err, rho, rho_err, phi, phi_err
!  end type phrate_report

  ! The lower limit for ratio estimation
  real, parameter :: RLOW_limit = 0.1

  ! The lower limit for Poisson statistics
  integer, parameter :: plimit = 666

  logical :: debug_graph = .true.


!  real(REAL64), dimension(:), allocatable :: res,sig,psi,dpsi!,mju,tt,tmin,tmax,kmean
!  real(REAL64), dimension(4) :: global
!  real(REAL64), dimension(3) :: maxlike_old
!  real(REAL64) :: lambda, mu
!  real(REAL64) :: tmin, tmax, ttiny, thuge, rmax
!  real(REAL64) :: tmin, tmax, rmax
!  real(REAL64) :: rmax!, cdf_min, cdf_max


  character(len=*), parameter :: error_label = 'Error in phrate: '
  character(len=*), parameter :: warning_label = 'Warning in phrate: '


  interface phrate
     module procedure phrate64
  end interface phrate

contains

  subroutine phrate64(counts,photons,&
       tau,tauerr,sigcts,sigpht,ctsfac,phtfac,reliable,flag,robustness, &
       clustering,verbose)

    use medians
    use photonrate_init
    use photonrate_kmean
    use photonrate_opt
    use photonrate_line
!    use photonrate_mood

    real(REAL64), dimension(:), target, intent(in) :: counts, photons
    real(REAL64), intent(out) :: tau
    real(REAL64), intent(out), optional :: tauerr, sigcts, sigpht
    real(REAL64), intent(in), optional :: ctsfac, phtfac
    logical, intent(out), optional :: reliable
    integer, intent(out), optional :: flag
    logical, intent(in), optional :: robustness, clustering, verbose

    ! working arrays
    real(REAL64), dimension(:), pointer :: cts, pht
    real(REAL64), dimension(:), allocatable :: mju,tt,kmean,tmin,tmax
    integer, dimension(:), allocatable :: iset

    integer, dimension(:), allocatable :: freq
    integer :: n,l,stat,nset,iter
!    integer :: n,i,nplus,nsig,nstep,stat,nflag
!    real(REAL64) :: t,dt,t0,dt0,t1,s,s0,r0,c,d,r,w,z,q,tmax,tmin,tol,std,std0,&
    real(REAL64) :: dtau,tau0,dtau0,dphi,phi0,rho,phi,drho,rho0,sigu,sigv,&
         dsigu,dsigv,sx,sy,tauold,s
    logical :: reli, robust, cluster, identical, singular
    character(len=80) :: errmsg

!    nflag = 0
    debug = .false.
    robust = .false.
    cluster = .false.
    if( present(verbose) ) debug = verbose
    if( present(robustness) ) robust = robustness
    if( present(clustering) ) cluster = clustering

    if( present(ctsfac) ) then
       if( ctsfac > 0 ) then
          sx = ctsfac
       else
          if( debug ) write(error_unit,'(2a)') &
               error_label,'assumption ctsfac > 0 failed.'
          if( present(flag) ) flag = 4
          if( present(reliable) ) reliable = .false.
          return
       end if
    else
       sx = 1
    end if
    if( present(phtfac) ) then
       if( phtfac > 0 ) then
          sy = phtfac
       else
          if( debug ) write(error_unit,'(2a)') &
               error_label,'assumption nsumy > 0 failed.'
          if( present(flag) ) flag = 4
          if( present(reliable) ) reliable = .false.
          return
       end if
    else
       sy = 1
    end if

    ! data check
    if( .not. (all(photons > 0) .and. all(counts > 0)) ) then
       write(error_unit,'(2a)') error_label, &
            'some data has non-positive values (assumption: all data > 0).'
       if( present(flag) ) flag = 4
       if( present(reliable) ) reliable = .false.
       return
    end if

    if(  size(counts) /= size(photons) ) then
       write(error_unit,'(2a)') error_label, &
            'The input data-sets has different sizes.'
       if( present(flag) ) flag = 4
       if( present(reliable) ) reliable = .false.
       return
    end if

    if( any(photons < plimit) .or. any(counts < plimit) .and. debug ) then
       write(error_unit,'(2a)') warning_label, &
            'potential accuracy lost: some datapoints under Poisson limit.'
    end if

    n = size(photons)
    pht => photons
    cts => counts

    if( n == 0 ) then ! no data
       tau = 1
       if( present(tauerr) ) tauerr = 0
       if( present(sigcts) ) sigcts = 0
       if( present(sigpht) ) sigpht = 0
       if( present(flag) ) flag = 2
       if( present(reliable) ) reliable = .false.
       return
    endif

    if( n == 1 ) then               ! a single point
       tau = counts(1) / photons(1)
       if( present(tauerr) ) tauerr = 1/sqrt(2*(photons(1) + tau**2*counts(1)))
       if( present(sigcts) ) sigcts = 0
       if( present(sigpht) ) sigpht = 0
       if( present(flag) ) flag = 2
       if( present(reliable) ) reliable = .true.
       return
    endif

    ! global working data buffers
    allocate(tt(n),mju(n),iset(n),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,'(2a)') warning_label, trim(errmsg)
       if( present(reliable) ) reliable = .false.
       if( present(flag) ) flag = 5
       return
    end if

    ! init ratios
    tt = cts / pht

    ! initial estimate of tau
    tau = qmedian(tt)
    if( debug_iter ) write(error_unit,*) 'E(cts/pht): ',real(tau)

    ! the reciprocal case
    if( tau < RLOW_limit ) then
       if( debug ) then
          write(error_unit,*) 'phrate: initial tau=',real(tau)
          write(error_unit,'(2a,f0.1,a/a,a)') error_label, &
            'The very first ratio estimate is under the lower limit <', &
            RLOW_limit,';','The mean ratio > 1 is required.', &
            ' Swap inputs counts and photons, please.'
       end if
       if( present(flag) )     flag = 4
       if( present(reliable) ) reliable = .false.
       goto 666
    end if

    ! initial estimate of mu-s
!    mju = (cts * tau + pht) / (1 + tau**2)
    mju = (cts / tau + pht) / 2

    if( cluster ) then
       ! initial clustering by mju
       call mucluster(mju,iset,nset,kmean,tau)
    else
       nset = size(mju)
       iset = [ (l,l=1,nset) ]
       kmean = mju
    end if

    ! initial estimate of parameters by quantiles
    allocate(freq(nset))
    do l = 1, nset
       freq(l) = count(iset == l)
    end do
    call phinit(cts,pht,kmean,freq,sx,sy,tt,tau,dtau,rho,reli,identical)
    deallocate(freq)

    if( .not. reli ) then
       if( debug ) then
          write(error_unit,'(3a)') error_label, &
               'Init setup by rank estimators failed; results are UNRELIABLE.'
       end if
       if( present(flag) )     flag = 1
       if( present(reliable) ) reliable = .false.
       goto 666
    end if

    if( reli .and. identical ) then
       if( present(tauerr) ) tauerr = 0
       if( present(sigcts) ) sigcts = 0
       if( present(sigpht) ) sigpht = 0
       if( present(reliable) ) reliable = reli
       if( present(flag) ) flag = 6
       goto 666
    end if

    ! listing of initial estimates
    if( debug ) then
       write(error_unit,'(a,g0,a,es8.1,a,g0)') &
            "Initial tau: ",real(tau),' +- ',dtau,'   rho: ',real(rho)
    end if

    if( cluster ) then
       ! update clustering by mju, rho included
       mju = (cts / tau + pht) / 2
       call mucluster(mju,iset,nset,kmean,tau,rho)
    end if

!    stop

!    call moodinit(iset,nset,tau,rho,cts,pht,robust,kmean,phi,singular)
!    write(error_unit,*) 'Moodinit:',real(kmean),real(RAD*phi),singular
!    stop


    !    write(*,*) real(tau),real(rho)
    !    write(*,*) nset,iset
!    write(*,*) real(mju)
!    stop

    sigu = rho / 2
    sigv = rho / 2

    mju = (cts / tau + pht) / 2

    do iter = 1, 1
       tauold = tau

       ! mean mjus
!       if( .not. singular ) &
!            call mood(cts,pht,iset,nset,sx,sy,tau,rho,drho,phi,dphi,mju,kmean,robust)


       !call phline(cts,pht,nset,iset,sx,sy,&
       !     tau,dtau,rho,phi,sigu,dsigu,sigv,dsigv,kmean,mju,robust)


!       sigu = rho*sin(phi)
!       sigv = rho*cos(phi)/tau
       !      write(*,*) '*****',sigu,sigv


    ! Looking for maximal and minimal tmin, tmax (ratio of data).
    if( allocated(tmin) ) deallocate(tmin,tmax)
    allocate(tmin(nset),tmax(nset),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,'(2a)') error_label,trim(errmsg)
       if( present(flag) )     flag = 5
       if( present(reliable) ) reliable = .true.
       goto 666
    end if
    if( robust ) then
       call tlimits(kmean,tau,rho,phi,sx,sy,tmin,tmax,reli)
       if( .not. reli ) then
          if( debug ) then
             write(error_unit,'(3a)') warning_label,'Looking for the ', &
                  'principial interval points of the robust function failed.'
          end if
          goto 666
       end if
    else
       ! this setup effectively switch off robust estimates
       tmin = minval(tt)
       tmin = sign(0.9*abs(tmin),tmin)
       tmax = 1.1*maxval(tt)
    end if

    if( debug_iter ) then
       do l = 1, size(kmean)
!          write(error_unit,*) 'tmin:',real(tmin(l)), 'tmax:',real(tmax(l))
       end do
    end if

    tau0 = tau
    rho0 = rho
    dtau0 = dtau
    phi0 = phi

    call phopt(mju,tt,kmean,iset,tmin,tmax,sx,sy,tau,dtau,rho,drho,phi,dphi,reli)
    if( .not. reli ) then
       if( debug ) write(error_unit,*) error_label, &
            ' failed to get accurate results in robust (or) least squares', &
            ' sense: initial rank estimates are used as the fail-back.'
       tau = tau0
       rho = rho0
       dtau = dtau0
       phi = phi0
       exit
    end if



    if( abs(tau - tauold) < 0.1*dtau .or. singular ) exit

 end do

    ! computation of scatter of the original data sets
    sigu = rho*sin(phi)
    sigv = (rho/tau)*cos(phi)

    if( drho > 0 .and. dphi > 0 ) then

       dsigu = drho**2*sin(phi)**2 + dphi**2*rho**2*cos(phi)**2
       if( dsigu > 0 ) then
          dsigu = sqrt(dsigu)
       else
          dsigu = -1
       end if

       dsigv = drho**2+dtau**2*(rho/tau**2)**2*cos(phi)**2 &
            + dphi**2*(rho/tau)**2*sin(phi)**2
       if( dsigv > 0 ) then
          dsigv = sqrt(dsigv)
       else
          dsigv = -1
       end if
    else
       dsigu = -1
       dsigv = -1
    end if


    ! listing of results
    if( debug ) then
       block
         real(REAL64) :: c,d,q,t,s
         integer :: i, nstep
         logical :: w

         nstep = 10**max(int(log10(real(n))) - 1,0)

         write(error_unit,'(a,a)') '# counts, photons,  rate,', &
              ' mu,  rel.res., poissonity, in winsorise:'
         do i = 1,n,nstep
            c = tau*pht(i)
            d = (cts(i) - c)/((cts(i) + c)/2)
            t = cts(i) / pht(i)
            q = sqrt(mju(i) / (mju(i) + sigv**2))
            w = tmin(iset(i)) < tt(i) .and. tt(i) < tmax(iset(i))
            write(error_unit,101) i,cts(i),pht(i),t,mju(i),d,q,w
101         format(i3,2(en14.3,1x),2x,f0.3,en14.3,1x,f9.5,2x,f0.2,l2)
         enddo
         write(error_unit,'(3(a,g0,a,es7.1:a,2x))') &
              'tau = ',real(tau),' +- ',dtau,',', &
              'rho = ',nint(rho),' +- ',drho,',', &
              'phi[deg] = ',nint(RAD*phi),' +- ',max(RAD*dphi,-real(1.0,REAL64))

         write(error_unit,'(2(a,g0,a,es7.1:a,3x))') &
              'sigcts = ', nint(sigu),' +- ', dsigu,',', &
              'sigpht = ', nint(sigv),' +- ', dsigv

!         call kstest(tt,kmean,iset,tmin,tmax,sx,sy,tau,rho,phi)
         call kstest(tt,kmean,iset,sx,sy,tau,dtau,rho,phi)
         write(error_unit,*) "Data points under/above limits: ",&
              count(tt < tmin(iset)), count(tt > tmax(iset))

         t = sum(tt) / n
         s = sqrt(sum((tt - t)**2) / (n - 1))
         write(error_unit,*) 'E(cst)/E(pht):',real(sum(cts)/sum(pht))
         write(error_unit,*) 'E(cst/pht):',real(t),'+-',real(sqrt(s**2 / n)),&
              '  D(cts/pht):',real(s)

       end block
    end if

    if( present(tauerr) ) tauerr = dtau
    if( present(sigcts) ) sigcts = sigu
    if( present(sigpht) ) sigpht = sigv
    if( present(reliable) ) reliable = reli
    if( present(flag) ) flag = 0

666 continue

    deallocate(tt,mju,iset)
    if( allocated(kmean) ) deallocate(kmean)
    if( allocated(tmin) ) deallocate(tmin,tmax)

  end subroutine phrate64


  subroutine tlimits(kmeans,tau,rho,phi,sx,sy,tmin,tmax,converge)

    use photonrate_normal

    ! setup principial interval end-points for the robust function

    real(REAL64), dimension(:), intent(in) :: kmeans
    real(REAL64), intent(in) :: tau, rho, phi, sx, sy
    real(REAL64), dimension(:), intent(out) :: tmin, tmax
    logical, intent(out) :: converge

    ! If the quantiles are symmetric, the normalised value of z is
    ! asymmetric; only one computation of inverse Normal is need.
!    real(real64), parameter :: qmax = 0.975
    real(real64), parameter :: qmax = 0.99
    real(real64), parameter :: qmin = 1 - qmax

    real(real64) :: zmin, zmax
    integer :: l

    ! Determination of z of N(0,1) for the given quantile.
    zmax = invNormal(qmax)
    zmin = -zmax
    ! This code can be replaced by pre-computed values for small speed up:
    ! qmax   zmax
    ! 0.9    1.28
    ! 0.95   1.64
    ! 0.975  1.96
    ! 0.99   2.32
    ! 0.995  2.58
    ! 0.999  3.09
!    zmin = -1

    do l = 1, size(kmeans)

       ! Solve equations for zmax = (mu/sig)*(tmax - tau), and one for tmin.
       call eq975(kmeans(l),tau,rho,phi,sx,sy,zmax,tmax(l),converge)
       if( .not. converge ) then
          if( debug ) write(error_unit,*) &
               warning_label,'tmax estimate failed',' (tlimits)'
          return
       end if

       call eq975(kmeans(l),tau,rho,phi,sx,sy,zmin,tmin(l),converge)
       if( .not. converge ) then
          if( debug ) write(error_unit,*) &
               warning_label,'tmin estimate failed',' (tlimits)'
          return
       end if
       tmin(l) = max(real(0.1,REAL64),tmin(l))

    end do

  end subroutine tlimits


  subroutine eq975(mu,tau,rho,phi,sx,sy,z,t,converge)

    real(REAL64), intent(in) :: mu,tau, rho, phi, z, sx, sy
    real(REAL64), intent(out) :: t
    logical, intent(out) :: converge

    real(REAL64) :: tol, c2, sig, f, df, d, t0, sig_t, r
    integer :: iter

    converge = .false.
    c2 = cos(phi)**2

!    write(error_unit,*) rho**2, mu*(tau + tau**2), tau, z, rho/mu
    if( rho**2 / 10 > mu*(tau + tau**2) ) then
       t = tau + z*(rho/mu)
       converge = .true.
       return
    end if


    sig = sqrt(mu*(tau + tau**2) + rho**2)
!    write(*,*) z,mu,sig,mu/sig
!    t = tau + z*(sig/mu)
    t0 = tau + z*(sig/mu)
    t = tau
!    write(*,*) t

    tol = 1e-4*(abs(t) + 1)
    do iter = 1, precision(t)
       sig = sqrt(mu*(tau/sy + t**2/sx) + rho**2*(1 + ((t/tau)**2 - 1)*c2))
       sig_t = t*(mu/sig + rho*(rho/sig)/tau**2*c2)
       r = (mu/sig)*(t - tau)
       f =  r - z
       df = mu/sig - r*(sig_t/sig)

       if( .not. (abs(df) > epsilon(df)) ) then
          if( debug ) write(error_unit,'(2a)') &
               warning_label," fail in eq975(): derivative df == 0"
          exit
       end if

       ! estimate update
       d = f / df
       t = t - d

!       if( debug_iter ) write(error_unit,'(a,i2,2x,g0,2es12.2,es14.4)') &
!            "t, incr., f, f': ",iter,real(t),d,f,df

       converge = abs(d) < tol
       if( converge ) exit

    end do

    if( .not. converge ) then
       t = t0
       converge = .true.
    end if

  end subroutine eq975


!  subroutine kstest(tt,kmean,iset,tmin,tmax,sx,sy,tau,rho,phi)
  subroutine kstest(t,kmean,iset,sx,sy,tau,dtau,rho,phi)

    ! K-S test on empirical and computed CDF in q = 0.2 .. 0.8
    ! https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test

    use quantilmean
    use photonrate_normal

    real(REAL64), dimension(:), intent(in) :: t, kmean!, tmin, tmax
    integer, dimension(:), intent(in) :: iset
    real(REAL64), intent(in) :: sx,sy,tau,dtau,rho,phi
    real(REAL64), parameter :: qlo = 0.2
    real(REAL64), parameter :: qhi = 0.8
    integer, parameter :: ncdf = 200
    character(len=80) :: errmsg
    integer, dimension(:), allocatable :: freq
    real(REAL64), dimension(:), allocatable :: xpdf,ypdf,u,v
    real(REAL64) :: dmax, sn, d, tlo, thi, umin, umax, s
    integer :: i,l,n,ns,stat
    logical :: mindef, maxdef

    n = size(t)
    allocate(xpdf(n),ypdf(n),freq(size(kmean)),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       if( debug ) write(error_unit,*) warning_label,trim(errmsg),' (kstest)'
       return
    end if

    do l = 1, size(kmean)
       freq(l) = count(iset == l)
    end do

    ! empirical CDF
    call ecdf(t,xpdf,ypdf)
    ! computed CDF
    v = mcdf(kmean,freq,sx,sy,tau,rho,phi,xpdf)

    ! the interval for the test
    tlo = xpdf(minloc(abs(v - qlo),1))
    thi = xpdf(minloc(abs(v - qhi),1))

    ns = count(tlo <= xpdf .and. xpdf <= thi)
    if( ns > 0 ) then
       dmax = maxval(abs(ypdf - v), tlo <= xpdf .and. xpdf <= thi)
       sn = sqrt(real(ns))
       write(error_unit,'(2(a,f0.2),2(a,f0.3),a,l1)') &
            'K-S test in quantile interval ',qlo,'--',qhi,' (5%, n>80, dmax=',&
            dmax,', crit.=',1.36/sn,'): ',sn*dmax < 1.36
    end if

    if( debug_graph ) then

       ! graph of empirical CDF
       open(1,file='/tmp/ecdf')
       write(1,'(3(a,g0),2(a,f0.3))') '# t = ',real(tau),', rho = ',real(rho), &
            ', phi = ',real(phi),',    dmax = ',dmax,' <? ',1.36/sqrt(real(n))
       do i = 1, n
          write(1,*) xpdf(i),ypdf(i)
       end do
       close(1)

       ! graph of computed CDF
!       umin = tau + 5*(tlo - tau)   !0.95*minval(tmin)
!       umax = tau + 5*(thi - tau) !1.05*maxval(tmax)
       s = dtau*sqrt(real(n))
       umin = tau - 5*s
       umax = tau + 5*s

       mindef = .false.
       maxdef = .false.
       do i = 30,666
          u = tau + 0.1*s*[ -i, i ]
          v = mcdf(kmean,freq,sx,sy,tau,rho,phi,u)
          if( v(1) < 0.05 .and. .not. mindef ) then
             umin = u(1)
             mindef = .true.
          end if
          if( v(2) > 0.95 .and. .not. maxdef ) then
             umax = u(2)
             maxdef = .true.
          end if
          if( maxdef .and. mindef ) exit
       end do
!       write(error_unit,*) umin, umax, tlo, thi
       d = (umax - umin) / (ncdf - 1)
       u = [ (umin + d*(i-1), i = 1,ncdf) ]
       v = mcdf(kmean,freq,sx,sy,tau,rho,phi,u)

       open(1,file='/tmp/fcdf')
       write(1,'(3(a,g0),2(a,f0.3))') '# tau = ',real(tau), &
            ', rho = ',real(rho), ', phi = ',real(phi), &
            ',    dmax = ',dmax,' <? ',1.36/sqrt(real(n))
       do i = 1, size(u)
          write(1,*) u(i),v(i)
       end do
       close(1)

    end if



!!$  subroutine residuals(t,r,res,sig)
!!$
!!$    real(REAL64), intent(in) :: t, r
!!$    real(REAL64), dimension(:), intent(out) :: res, sig
!!$    real(REAL64) :: r2
!!$
!!$    r2 = r**2
!!$!    sig = sqrt(mju*(t + tt**2) + r2*(1 + (tt/t)**2)/2)
!!$    sig = 1
!!$    res = (mju/sig)*(t - tt)
!!$
!!$  end subroutine residuals
!!$
!!$  subroutine residuals1(t,r,phi,res,sig)
!!$
!!$    real(REAL64), intent(in) :: t, r, phi
!!$    real(REAL64), dimension(:), intent(out) :: res, sig
!!$    real(REAL64) :: r2, c2
!!$
!!$    r2 = r**2
!!$    c2 = cos(phi)**2
!!$    res = mju*(t - tt)
!!$    sig = sqrt(mju*(t + tt**2) + r2*(1 + c2*((tt/t)**2-1)))
!!$
!!$  end subroutine residuals1
!!$
!!$
!!$  subroutine reziduals(t,r,z,res,sig)
!!$
!!$    real(REAL64), intent(in) :: t, r, z
!!$    real(REAL64), dimension(:), intent(out) :: res, sig
!!$    real(REAL64) :: r2, t2
!!$
!!$    t2 = t**2
!!$    r2 = r**2
!!$    res = cts - t*pht + z
!!$    sig = sqrt(cts + pht*t2 + r2)
!!$
!!$  end subroutine reziduals


!!$  function mood1(tau,rho,phi,y,x) result(mu)
!!$
!!$    real(REAL64), intent(in) :: tau,rho,phi
!!$    real(REAL64), dimension(:), intent(in) :: x,y
!!$
!!$    real(REAL64), dimension(:), allocatable :: rx, ry
!!$    real(REAL64) :: mu, sigu, sigv, sigx, sigy, f, f2, df, d, tol, dmu
!!$    integer :: n, iter
!!$    logical :: converge
!!$
!!$    converge = .false.
!!$    n = size(y)
!!$    sigu = rho*sin(phi)
!!$    sigv = rho*sin(phi)/tau
!!$
!!$    mu = sum((y*tau + x) / (1 + tau**2)) / n
!!$    mu = sum((y/tau + x) / 2)/n
!!$    ! mu = median((y/tau + x) / 2)
!!$    tol = sqrt(10.0*n)*(0.001*sqrt(mu))
!!$    write(*,*) 'tol:',tol
!!$    do iter = 1, precision(mu)
!!$
!!$       sigx = sqrt(mu + sigv**2)
!!$       sigy = sqrt(mu*tau + sigu**2)
!!$       rx = x - mu
!!$       ry = y - mu*tau
!!$       f = tau*sum(ry/sigy**2) + tau*sum((ry/sigy**2)**2) - 2*tau*n/sigy**2 &
!!$           + sum(rx/sigx**2) + tau*sum((rx/sigx**2)**2) - 2*n/sigx**2
!!$
!!$       df = -tau**2*n/sigy**2 - 4*tau**2*sum((ry/sigy**2)/sigy**2) &
!!$            -4*tau**2*sum((ry/sigy**2)**2/sigy**2) - 4*tau**2*n/sigy**4 &
!!$            -n/sigx**2 - 4*sum((rx/sigx**2)/sigx**2) &
!!$            -4*sum((rx/sigx**2)**2/sigx**2) - 4*n/sigx**4
!!$
!!$
!!$!       f = -(y - mu*tau)/sigy**2 - (x - tau)/sigx**2 &
!!$!            +((y - mu*tau)/sigy**2)**2*tau/2 + ((x - mu)/sigy**2)**2/2 &
!!$!            -tau/sigy**2 - 1/sigx**2
!!$!       df = tau/sigy**2 + (y - mu*tau)*tau/sigy**4 &
!!$!            + 1/sigx**2 + (x - mu)/sigy**4 &
!!$!            + (y - mu*tau)*tau**2/sigy**4 &
!!$!            - 3*((y - mu*tau)/sigy**2)**2*(tau/sigy)**2/4 &
!!$!            -0.5*(y - mu*tau)/sigy**2)**2*tau/sigy &
!!$!            -0.5/sigx**4 - 3*(x - mu)/sigx**5/4 + 0.5*(x - mu)/sigx**5 &
!!$!            +0.5*tau/sigy**5 + 0.5*tau/sigy**3 +0.5/six**4 &
!!$!            +0.5*
!!$
!!$       if( abs(df) < epsilon(df) ) exit
!!$
!!$       d = f / df
!!$       mu = mu - d
!!$
!!$       if( debug_iter ) write(error_unit,'(a,i2,2x,g0,2es12.2,es14.4)') &
!!$            "mu, incr., f, f': ",iter,real(mu),d,f,df
!!$
!!$       converge = abs(d) < tol
!!$
!!$       if( converge ) then
!!$!          f2 = tau**2*sum((ry/sigy**2)**2) + tau**2*sum(((ry/sigy**2)**2)**2) &
!!$!               - (2*tau*n/sigy**2)**2 &
!!$!               + sum((rx/sigx**2)**2) + tau**2*sum(((rx/sigx**2)**2)**2) &
!!$          !               - (2*n/sigx**2)**2
!!$          f2 = sum(ry**2)/sigy**2 + sum(rx**2)/sigx**2
!!$          dmu = sqrt(f2/(-df)/(n-1))
!!$          write(*,*) 'dmu:',dmu
!!$
!!$       end if
!!$
!!$
!!$       if( abs(d) < tol ) exit   ! when a local maximum is reached
!!$       if( df > 0 ) exit         ! when a local minimum is encountered
!!$
!!$
!!$       ! exit of iterations: the absolute error must be at least |d| < tol
!!$       if( converge ) exit
!!$
!!$    end do
!!$
!!$  end function mood1




!!$  subroutine sigmas(mju,t0,r,phi,s1,s2)
!!$
!!$    real(REAL64), parameter :: log2 = -log(2.0_REAL64)
!!$
!!$    real(REAL64), intent(in) :: mju, t0, r, phi
!!$    real(REAL64), intent(out) :: s1, s2
!!$    real(REAL64) :: dmin, d, t, logf
!!$    integer :: i
!!$
!!$    dmin = huge(dmin)
!!$    do i = -1000, 10000
!!$       t = t0 + 0.001*i
!!$!       sig = sqrt(mju*(t0 + t**2) + r**2*(1 + ((t/t0)**2-1)*cos(phi)**2))
!!$!       beta = mju*t0*(t+1) + r**2*(1 + (t/t0)-1)*cos(phi)**2
!!$       !       logf = log(mju*beta/sig**3/sqrtpi) - (mju*(t-t0)/sig)**2/2
!!$       logf = logpdf(t,mju,t0,r,phi)
!!$       d = abs(logf - log2)
!!$!       write(*,*) t,logf,log2,d
!!$       if( d < dmin ) then
!!$          dmin = d
!!$          s2 = t
!!$       end if
!!$!       if( logf < 1.1*log2 ) exit
!!$    end do
!!$
!!$    dmin = huge(dmin)
!!$    do i = -1000, 10000
!!$       t = t0 - 0.001*i
!!$!       if( logf < 1.1*log2 ) exit
!!$       if( t < RLOW_limit ) exit
!!$!       sig = sqrt(mju*(t0 + t**2) + r**2*(1 + ((t/t0)**2-1)*cos(phi)**2))
!!$!       beta = mju*t0*(t+1) + r**2*(1 + (t/t0)-1)*cos(phi)**2
!!$       !       logf = log(mju*beta/sig**3/sqrtpi) - (mju*(t-t0)/sig)**2/2
!!$       logf = logpdf(t,mju,t0,r,phi)
!!$       d = abs(logf - log2)
!!$       if( d < dmin ) then
!!$          dmin = d
!!$          s1 = t
!!$       end if
!!$    end do
!!$
!!$    write(*,*) 'sigmas: ',s1, s2
!!$
!!$  end subroutine sigmas
!!$
!!$  function logpdf(t,mju,t0,r,phi)
!!$
!!$    real(REAL64), parameter :: sqrtpi = 2.5066282746310005_REAL64
!!$
!!$    real(REAL64) :: logpdf
!!$    real(REAL64), intent(in) :: t,mju,t0,r,phi
!!$    real(REAL64) :: sig, beta, c
!!$
!!$    c = cos(phi)
!!$    sig = sqrt(mju*(t0 + t**2) + r**2*(1 + ((t/t0)**2-1)*c**2))
!!$    beta = mju*t0*(t+1) + r**2*(1 + (t/t0)-1)*c**2
!!$    logpdf = log((mju/sig)*(beta/sig**2)/sqrtpi) - (mju*(t-t0)/sig)**2/2
!!$
!!$  end function logpdf
!!$


!!$  subroutine thepeak(t,dt,r,dr,z,dz,s,ds,s0,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, dt, r, s
!!$    real(REAL64), intent(out) :: z, dr, dz, ds, s0
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: npar = 4
!!$    real(REAL64), dimension(npar) :: x, dx
!!$    real(REAL64), dimension(npar,npar) :: cov
!!$    logical, dimension(npar,npar) :: diagonal
!!$    character(len=80) :: errmsg
!!$!    real(REAL64), parameter :: eps2 = sqrt(epsilon(t))
!!$    real(REAL64) :: ttol,rtol,ztol,stol,rmax,rmin, lambda0, f, f0, c, eps
!!$    integer :: iter, info, nprint, stat, n, i
!!$    logical :: converge
!!$    character :: kind
!!$
!!$    diagonal = .false.
!!$    do i = 1, npar
!!$       diagonal(i,i) = .true.
!!$    end do
!!$
!!$    reli = .false.
!!$    nprint = 0
!!$    n = size(res)
!!$    allocate(dpsi(n),stat=stat,errmsg=errmsg)
!!$    ! Computations needs some additional memory,
!!$    ! a fail-back way should be used if it is unavailable.
!!$    if( stat /= 0 ) then
!!$       if( debug ) write(error_unit,*) warning_label,trim(errmsg),' (thepeak)'
!!$!       reli = .false.
!!$       return
!!$    end if
!!$
!!$    eps = 10*n*epsilon(t)*(1 + abs(t))
!!$    rmax = maxval(abs(cts - t*pht))
!!$    if( abs(r) > 10*n*sqrt(epsilon(r)) ) then
!!$       rtol = 1e-5*r
!!$    else
!!$       rtol = 1e-5*rmax
!!$    end if
!!$    ttol = max(1e-3*dt,10*n*sqrt(epsilon(t)))
!!$!    rtol = 0.01*dr
!!$    stol = 0.5e-4
!!$    ztol = ttol
!!$
!!$    dt = -1
!!$    dr = -1
!!$    ds = -1
!!$    dz = -1
!!$    s0 = -1
!!$
!!$    write(*,*) ttol, rtol, stol, ztol, (r / minval(cts))**2
!!$
!!$    z = 0
!!$    dz = (r / minval(cts))**2
!!$    dz = ztol
!!$    call reziduals(t,r,z,res,sig)
!!$    c = s - stukey
!!$!    c = s - 1
!!$    if( abs(c) > 0 ) then
!!$       mu = sum(itukey(res/(s*sig))) / (10*c**2)
!!$    else
!!$       mu = sum(itukey(res/(s*sig))) / 10
!!$    end if
!!$    mu = sum(itukey(res/(s*sig))) / 10
!!$    if( mu < 0.1 ) mu = 1
!!$    lambda = mu*c
!!$
!!$!    write(*,*) mu,lambda,c
!!$!    stop 0
!!$
!!$    block
!!$      t = 2.5
!!$      r = 400
!!$      s = 0.95
!!$      mu = 100
!!$      lambda = -2
!!$      call testders(t,r,s)
!!$      stop 0
!!$    end block
!!$
!!$    converge = .false.
!!$    do iter = 1, 3*precision(r)
!!$
!!$       x = [ t, r, z, s ]
!!$       !       call lmder3(peakder,x,cov,eps,666,1d0,nprint,info)
!!$       info = -1
!!$       kind = 'D'
!!$       !call lmdif2(mountder1,x,eps,nprint,info)
!!$       cov = -1
!!$!       write(*,*) pack(cov,diagonal)
!!$       if( any(x < 0) .or. any( pack(cov,diagonal) < 0) .or. &
!!$            info == 666 .or. info < 0 ) then
!!$          if( debug_iter ) write(error_unit,*) &
!!$               'lmder3 failed at ',real(x),'for',info,' or reached a maximum.'
!!$
!!$          x = [ t, r, z, s ]
!!$          dx = [ 0.1*dt, 0.1*dr, 0.01*dz, 1d-2 ]
!!$          call nelmin1(minpeak,x,dx,f,info)
!!$          kind = 'L'
!!$!          if( debug_iter ) write(error_unit,*) 'nelmin: ',real(x),info
!!$       end if
!!$
!!$       if( debug_iter ) then
!!$          c = s - stukey
!!$          if( s > 0 ) then
!!$             call reziduals(x(1),x(2),x(3),res,sig)
!!$             f = sum(itukey(res/(s*sig)))+sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$          else
!!$             f = 0
!!$          end if
!!$          write(error_unit, &
!!$               '(3a,i0,2x,g0,2x,f0.1,2x,g0,2x,f6.4,1x,f10.5,2x,2(g0,2x),i0)') &
!!$               'peak(',kind,'): #',iter,real(x),lambda,real(mu),real(f),info
!!$       end if
!!$!       stop 0
!!$
!!$       if( (iter > 1 .and. (x(2) > rmax .or. x(4) < 0.1 .or. x(4) > 2 ))  &
!!$            .or. info < 0 ) then
!!$          converge = .false.
!!$          exit
!!$       end if
!!$
!!$!       write(*,*) x(1) - t, x(2) - r, x(3) - stukey
!!$       converge = abs(x(1) - t) < ttol  .and. abs(x(2) - r) < rtol .and. &
!!$            abs(x(3) - z) < ztol .and. abs(x(4) - stukey) < stol
!!$
!!$       t = x(1)
!!$       r = x(2)
!!$       z = x(3)
!!$       s = x(4)
!!$
!!$       call reziduals(t,r,z,res,sig)
!!$       if( s > 0 .and. all(pack(cov,diagonal) > 0) ) then
!!$          s0 = s**2*sum(tukey(res/(s*sig))**2)/ &
!!$               sum(dtukey(res/(s*sig)))*real(n)**2/real(n-1)
!!$!          write(*,*) s0,s,sum(tukey(res/sig/s)**2),sum(dtukey(res/sig/s))
!!$          if( s0 > 0 ) then
!!$             dt = sqrt(s0*cov(1,1)/(n-3))
!!$             dr = sqrt(s0*cov(2,2)/(n-3))
!!$             ds = sqrt(s0*cov(3,3)/(n-3))
!!$             ttol = max(0.001*dt,ttol)
!!$             rtol = max(0.01*dr,rtol)
!!$             ztol = max(0.01*dz,ztol)
!!$!             stol = max(0.01*ds,stol)
!!$             s0 = sqrt(s0 / n)
!!$          end if
!!$       else
!!$          s0 = -1
!!$          dt = -1
!!$          dr = -1
!!$          dz = -1
!!$          ds = -1
!!$       end if
!!$
!!$       if( converge ) exit
!!$
!!$       lambda = lambda - mu*(s - stukey)
!!$       mu = 2*mu
!!$
!!$    end do
!!$
!!$!    if( info == 0 ) then
!!$!       t = x(1)
!!$!       r = x(2)
!!$!       s = x(3)
!!$!    end if
!!$
!!$!    block
!!$!      call testders(t,s,r)
!!$!      stop 0
!!$!    end block
!!$
!!$
!!$
!!$    reli = converge
!!$!    r = abs(r)
!!$    deallocate(dpsi)
!!$
!!$    if( reli .and. debug_iter ) then
!!$       block
!!$         real(REAL64), dimension(npar,npar) :: cor
!!$         real(REAL64) :: q
!!$         integer :: i,j
!!$
!!$         !         if( cov(1,1) > 0 .and. cov(2,2) > 0 .and. cov(3,3) > 0 ) then
!!$         if( any(pack(cov,diagonal) > 0) ) then
!!$            cor = 1
!!$            do i = 1, npar
!!$               do j = i+1, npar
!!$                  q = sqrt(cov(i,i)*cov(j,j))
!!$                  cor(i,j) = cov(i,j) / q
!!$                  cor(j,i) = cor(i,j)
!!$               end do
!!$            end do
!!$         else
!!$            cor = -1
!!$         end if
!!$         write(error_unit,'(a,4(es12.4,1x),t48,a,4(f7.4,1x))') &
!!$              'Cov: ', real(cov(1,:)),'corr: ', cor(1,:)
!!$         do i = 2,npar
!!$            write(error_unit,'(5x,3(es12.4,1x),t48,6x,3(f7.4,1x))') &
!!$                 real(cov(i,:)),cor(i,:)
!!$         end do
!!$       end block
!!$    end if
!!$
!!$
!!$  end subroutine thepeak
!!$
!!$  function minpeak(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64) :: t,r,z,s,c,f
!!$
!!$    t = x(1)
!!$    r = x(2)
!!$    z = x(3)
!!$    s = x(4)
!!$    c = s - stukey
!!$
!!$    if( .not. (0.1 < s .and. s < 2) ) then
!!$       f = mu*c**4
!!$       return
!!$    end if
!!$
!!$    call reziduals(t,r,z,res,sig)
!!$!    f = sum(itukey(res/(s*sig))) + sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$    f = sum((res/(s*sig))**2/2) + sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$    write(error_unit,*) 'minlike:',real(x),real(f)
!!$
!!$  end function minpeak
!!$
!!$
!!$  subroutine thepeak0(t,dt,r,dr,s,ds,s0,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in out) :: t, dt, r, s
!!$    real(REAL64), intent(out) :: dr, ds, s0
!!$    logical, intent(out) :: reli
!!$
!!$    real(REAL64), dimension(3) :: x, dx
!!$    real(REAL64), dimension(3,3) :: cov
!!$    character(len=80) :: errmsg
!!$!    real(REAL64), parameter :: eps2 = sqrt(epsilon(t))
!!$    real(REAL64) :: ttol,rtol,stol,rmax,rmin, lambda0, f, f0, c, eps
!!$    integer :: iter, info, nprint, stat, n
!!$    logical :: converge
!!$    character :: kind
!!$
!!$    reli = .false.
!!$    nprint = 0
!!$    n = size(res)
!!$    allocate(dpsi(n),stat=stat,errmsg=errmsg)
!!$    ! Computations needs some additional memory,
!!$    ! a fail-back way should be used if it is unavailable.
!!$    if( stat /= 0 ) then
!!$       if( debug ) write(error_unit,*) warning_label,trim(errmsg),' (thepeak)'
!!$!       reli = .false.
!!$       return
!!$    end if
!!$
!!$    eps = 10*n*epsilon(t)*(1 + abs(t))
!!$    rmax = maxval(abs(cts - t*pht))
!!$    if( abs(r) > 10*n*sqrt(epsilon(r)) ) then
!!$       rtol = 1e-5*r
!!$    else
!!$       rtol = 1e-5*rmax
!!$    end if
!!$    ttol = max(1e-3*dt,10*n*sqrt(epsilon(t)))
!!$!    rtol = 0.01*dr
!!$    stol = 0.5e-4
!!$
!!$    dt = -1
!!$    dr = -1
!!$    ds = -1
!!$    s0 = -1
!!$
!!$!    write(*,*) ttol, rtol, stol
!!$
!!$    call residuals(t,r,res,sig)
!!$    c = s - stukey
!!$!    c = s - 1
!!$    if( abs(c) > 0 ) then
!!$       mu = sum(itukey(res/(s*sig))) / (10*c**2)
!!$    else
!!$       mu = sum(itukey(res/(s*sig))) / 10
!!$    end if
!!$    mu = sum(itukey(res/(s*sig))) / 10
!!$    if( mu < 0.1 ) mu = 1
!!$    lambda = mu*c
!!$
!!$!    write(*,*) mu,lambda,c
!!$!    stop 0
!!$
!!$    block
!!$      t = 2.5
!!$      r = 400
!!$      s = 0.95
!!$      mu = 100
!!$      lambda = -2
!!$      call testders(t,r,s)
!!$      stop 0
!!$    end block
!!$
!!$    converge = .false.
!!$    do iter = 1, 3*precision(r)
!!$
!!$       x = [ t, r, s ]
!!$       call lmder3(mountder,x,cov,eps,666,1d0,nprint,info)
!!$       kind = 'D'
!!$       !call lmdif2(mountder1,x,eps,nprint,info)
!!$!       cov = 1
!!$       if( any(x < 0) .or. info == 666 .or. info < 0 .or. &
!!$            cov(1,1) < eps .or. cov(2,2) < eps .or. cov(3,3) < eps ) then
!!$          if( debug_iter ) write(error_unit,*) &
!!$               'lmder3 failed at ',real(x),'for',info,' or reached a maximum.'
!!$
!!$          x = [ t, r, s ]
!!$          dx = [ 0.1*dt, 0.1*dr, 1d-2 ]
!!$          call nelmin1(minlike,x,dx,f,info)
!!$          kind = 'L'
!!$!          if( debug_iter ) write(error_unit,*) 'nelmin: ',real(x),info
!!$       end if
!!$
!!$       if( debug_iter ) then
!!$          c = s - stukey
!!$          if( s > 0 ) then
!!$             f = sum(itukey(res/(s*sig)))+sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$          else
!!$             f = 0
!!$          end if
!!$          write(error_unit, &
!!$               '(3a,i0,2x,g0,2x,f0.1,2x,f6.4,1x,f10.5,2x,2(g0,2x),i0)') &
!!$               'peak(',kind,'): #',iter,real(x),lambda,real(mu),real(f),info
!!$       end if
!!$
!!$       if( (iter > 1 .and. (x(2) > rmax .or. x(3) < 0.1 .or. x(3) > 2 ))  &
!!$            .or. info < 0 ) then
!!$          converge = .false.
!!$          exit
!!$       end if
!!$
!!$!       write(*,*) x(1) - t, x(2) - r, x(3) - stukey
!!$       converge = abs(x(1) - t) < ttol  .and. abs(x(2) - r) < rtol .and. &
!!$            abs(x(3) - stukey) < stol
!!$
!!$       t = x(1)
!!$       r = x(2)
!!$       s = x(3)
!!$
!!$       call residuals(t,r,res,sig)
!!$       if( s > 0 .and. cov(1,1) > 0 .and. cov(2,2) > 0 .and. cov(3,3) > 0 ) then
!!$          s0 = s**2*sum(tukey(res/(s*sig))**2)/ &
!!$               sum(dtukey(res/(s*sig)))*real(n)**2/real(n-1)
!!$!          write(*,*) s0,s,sum(tukey(res/sig/s)**2),sum(dtukey(res/sig/s))
!!$          if( s0 > 0 ) then
!!$             dt = sqrt(s0*cov(1,1)/(n-3))
!!$             dr = sqrt(s0*cov(2,2)/(n-3))
!!$             ds = sqrt(s0*cov(3,3)/(n-3))
!!$             ttol = max(0.001*dt,ttol)
!!$             rtol = max(0.01*dr,rtol)
!!$!             stol = max(0.01*ds,stol)
!!$             s0 = sqrt(s0 / n)
!!$          end if
!!$       else
!!$          s0 = -1
!!$          dt = -1
!!$          dr = -1
!!$          ds = -1
!!$       end if
!!$
!!$       if( converge ) exit
!!$
!!$       lambda = lambda - mu*(s - stukey)
!!$       mu = 2*mu
!!$
!!$    end do
!!$
!!$!    if( info == 0 ) then
!!$!       t = x(1)
!!$!       r = x(2)
!!$!       s = x(3)
!!$!    end if
!!$
!!$!    block
!!$!      call testders(t,s,r)
!!$!      stop 0
!!$!    end block
!!$
!!$
!!$
!!$    reli = converge
!!$!    r = abs(r)
!!$    deallocate(dpsi)
!!$
!!$    if( reli .and. debug_iter ) then
!!$       block
!!$         real(REAL64), dimension(3,3) :: cor
!!$         real(REAL64) :: q
!!$         integer :: i,j
!!$
!!$         if( cov(1,1) > 0 .and. cov(2,2) > 0 .and. cov(3,3) > 0 ) then
!!$            cor = 1
!!$            do i = 1, 3
!!$               do j = i+1, 3
!!$                  q = sqrt(cov(i,i)*cov(j,j))
!!$                  cor(i,j) = cov(i,j) / q
!!$                  cor(j,i) = cor(i,j)
!!$               end do
!!$            end do
!!$         else
!!$            cor = -1
!!$         end if
!!$         write(error_unit,'(a,3(es12.4,1x),t48,a,3(f7.4,1x))') &
!!$              'Cov: ', real(cov(1,:)),'corr: ', cor(1,:)
!!$         do i = 2,3
!!$            write(error_unit,'(5x,3(es12.4,1x),t48,6x,3(f7.4,1x))') &
!!$                 real(cov(i,:)),cor(i,:)
!!$         end do
!!$       end block
!!$    end if
!!$
!!$
!!$  end subroutine thepeak0
!!$
!!$
!!$  function minlike(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64) :: t,r,s,c,f
!!$
!!$    t = x(1)
!!$    r = x(2)
!!$    s = x(3)
!!$    c = s - stukey
!!$
!!$    if( .not. (0.1 < s .and. s < 2) ) then
!!$       f = mu*c**4
!!$       return
!!$    end if
!!$
!!$    call residuals(t,r,res,sig)
!!$    f = sum(itukey(res/(s*sig))) + sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$!    f = sum((res/(s*sig))**2/2) + sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$!    write(error_unit,*) 'minlike:',real(x),real(f)
!!$
!!$  end function minlike
!!$
!!$  subroutine mountder(k,l,x,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k, l, ldfjac
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$    real(REAL64), dimension(ldfjac,k), intent(out) :: fjac
!!$
!!$    real(REAL64), dimension(:), allocatable :: tder
!!$    real(REAL64) :: t,r,s,c
!!$    integer :: n, i, j
!!$
!!$    t = x(1)
!!$    r = x(2)
!!$    s = x(3)
!!$    n = size(pht)
!!$    c = s - stukey
!!$
!!$!    if( .not. (s > 0 .and. r >= 0) ) then
!!$    if( .not. (s > 0) ) then
!!$       iflag = -1
!!$       if( debug_iter ) write(error_unit,*) &
!!$            'fder: violation of the condition: s > 0, encountered: s = ',real(s)
!!$       return
!!$    end if
!!$
!!$    if( iflag == 0 .and. debug_iter ) then
!!$       write(error_unit,*) 'x = ',real(x)
!!$       write(error_unit,*) 'fvec = ',real(fvec)
!!$       write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
!!$       write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
!!$       write(error_unit,*) 'fjac(3,:) = ',real(fjac(3,:))
!!$       return
!!$    end if
!!$
!!$!    if( .not. (0.1 < s .and. s < 2) ) then
!!$!       iflag = -1
!!$!       return
!!$!    end if
!!$
!!$    if( iflag == 1 ) then
!!$       !       nev = nev + 1
!!$       call residuals(t,r,res,sig)
!!$!       call tukeys(res/(s*sig),psi,dpsi)
!!$       psi = tukey(res/(s*sig))
!!$!       psi = huber(res/(s*sig))
!!$!       psi = res/(s*sig)
!!$       fvec(1) = -sum(psi*((pht/sig) + t*(res/sig)*(pht/sig)/sig))/s &
!!$            +t*sum((pht/sig)/sig)
!!$       fvec(2) = -sum(psi*(r/sig)*(res/sig)/sig)/s + sum((r/sig)/sig)
!!$       fvec(3) = -sum(psi*(res/sig))/s**2 + n/s - lambda + mu*c
!!$    end if
!!$
!!$    if( iflag == 2 ) then
!!$!       nev = nev + 1
!!$       call residuals(t,r,res,sig)
!!$       call tukeys(res/(s*sig),psi,dpsi)
!!$!       psi = huber(res/(s*sig))
!!$!       dpsi = dhuber(res/(s*sig))
!!$!       psi = res/(s*sig)
!!$!       dpsi = 1
!!$       allocate(tder(n))
!!$       tder = pht/sig + t*(pht/sig)*(res/sig)/sig
!!$       fjac(1,1) = sum(dpsi*tder**2)/s**2 &
!!$            -sum(psi*((res/sig)*(pht/sig)/sig - t*(pht/sig)**2/sig  &
!!$            -t*(pht/sig)**2*(1/sig + 3*t*(res/sig)/sig**2)))/s &
!!$            +sum(pht/sig**2) - 2*t**2*sum((pht/sig)**2/sig**2)
!!$!       fjac(1,2) = +r*sum(dpsi*(pht/sig)*(res/sig)/sig**2)/s**2 &
!!$!            +r*sum(psi*pht/sig**3)/s &
!!$!            +r*sum(dpsi*(cts/sig)*(pht/sig)*(res/sig)/sig**2)/s**2 &
!!$!            +3*t*sum(psi*(cts/sig)*(pht/sig)*(res/sig)/sig**2)/s &
!!$!            -r*t**2*sum(dpsi*(pht/sig)**2*(res/sig)/sig**3)/s**2 &
!!$!            -3*r*t**2*sum(psi*(pht/sig)**2/sig**3)/s
!!$       fjac(1,2) = sum(dpsi*tder*(r/sig)*(res/sig)/sig)/s**2 &
!!$            +sum(psi*((pht/sig)*(r/sig)/sig &
!!$            +3*t*(pht/sig)*(r/sig)*(res/sig)/sig**2))/s &
!!$            -2*t*sum((r/sig)*(pht/sig)/sig**2)
!!$
!!$!       fjac(1,3) = +sum(dpsi*(pht/sig)*(res/sig))/s**3 &
!!$!            +t*sum(dpsi*(cts/sig)*(pht/sig)*(res/sig)/sig)/s**3 &
!!$!            -t**2*sum(dpsi*(pht/sig)**2/sig**2)/s**3 &
!!$!            +sum(psi*pht/sig)/s**2 &
!!$!            +t*sum(psi*(cts/sig)*(pht/sig)/sig)/s**2 &
!!$!            -t**2*sum(psi*(pht/sig)**2/sig)/s**2
!!$       fjac(1,3) = sum(psi*tder)/s**2 +sum(dpsi*tder*(res/sig))/s**3
!!$       fjac(2,2) = -sum(psi*(res/sig)/sig**2)/s + &
!!$            sum(dpsi*((res/sig)**2*(r/sig)**2)/sig**2)/s**2 &
!!$            + 3*sum(psi*((res/sig)*(r/sig)**2/sig**2))/s &
!!$            +sum(1/sig**2) - 2*sum((r/sig**2)**2)
!!$       fjac(2,3) = sum(psi*(res/sig)*(r/sig)/sig)/s**2 &
!!$            +sum(dpsi*(res/sig)**2*(r/sig)/sig)/s**3
!!$       fjac(3,3) = 2*sum(psi*(res/sig))/s**3 + sum(dpsi*(res/sig)**2)/s**4 &
!!$            -n/s**2 + mu
!!$       do i = 2,3
!!$          do j = 1, i-1
!!$             fjac(i,j) = fjac(j,i)
!!$          end do
!!$       end do
!!$    end if
!!$
!!$  end subroutine mountder
!!$
!!$  subroutine mountder1(k,l,x,fvec,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k, l
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$
!!$    real(REAL64) :: t,r,s,c
!!$    integer :: n
!!$
!!$    t = x(1)
!!$    r = x(2)
!!$    s = x(3)
!!$    n = size(pht)
!!$    c = s - stukey
!!$
!!$    if( iflag == 0 .and. debug_iter ) then
!!$       write(error_unit,*) 'x = ',real(x)
!!$       write(error_unit,*) 'fvec = ',real(fvec)
!!$       return
!!$    end if
!!$
!!$    if( .not. (0.5 < s .and. s < 1.5 .and. r >= 0) ) then
!!$       iflag = -1
!!$       return
!!$    end if
!!$
!!$    call residuals(t,r,res,sig)
!!$    psi = tukey(res/(s*sig))
!!$    fvec(1) = -sum(psi*((pht/sig) + t*(res/sig)*(pht/sig)/sig))/s &
!!$         +t*sum(pht/sig**2)
!!$    fvec(2) = -sum(psi*(r/sig)*(res/sig)/sig)/s + sum((r/sig)/sig)
!!$    fvec(3) = -sum(psi*(res/sig))/s**2 + n/s - lambda + mu*c
!!$
!!$  end subroutine mountder1
!!$
!!$
!!$  subroutine testders(t,r,s)
!!$
!!$    real(REAL64), intent(in) :: t,r,s
!!$    real(REAL64), dimension(3) :: fd,fvec,eps,x,x1,x2
!!$    real(REAL64), dimension(3,3) :: fjac, djac
!!$    real(REAL64) :: f0
!!$    integer :: i, j, iflag
!!$
!!$    eps = 1e-7
!!$!    write(*,*) 'qqqq',t,r,s
!!$
!!$    iflag = 1
!!$    x = [ t, r, s ]
!!$    call mountder(3,3,x,fvec,fjac,2,iflag)
!!$    f0 = minlike(x)
!!$    do i = 1, 3
!!$       x = [ t, r, s ]
!!$       x(i) = x(i) + eps(i)
!!$       fd(i) = (minlike(x) - f0) / eps(i)
!!$       write(error_unit,*) i,real(fvec(i)),real(fd(i)),real((fd(i)-fvec(i)))
!!$    end do
!!$
!!$!    return
!!$
!!$    iflag = 2
!!$    x = [ t, r, s ]
!!$    call mountder(3,3,x,fvec,fjac,3,iflag)
!!$    call fdjac1(mdertest,3,x,fvec,djac,3,iflag,3,3,epsilon(x),x1,x2)
!!$    do i = 1,3
!!$       do j = i, 3
!!$          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
!!$               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))!/fjac(i,j))
!!$       end do
!!$    end do
!!$
!!$  end subroutine testders
!!$
!!$  subroutine mdertest(k,x,fvec,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(k), intent(out) :: fvec
!!$
!!$    real(REAL64), dimension(3,3) :: fjac
!!$    integer :: nflag
!!$
!!$    nflag = 1
!!$    call mountder(3,3,x,fvec,fjac,3,nflag)
!!$
!!$  end subroutine mdertest

!!$  function testmnt(x)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64) :: t,r,s,testmnt,c
!!$    integer :: n
!!$
!!$    n = size(res)
!!$    t = x(1)
!!$    r = x(2)
!!$    s = x(3)
!!$    call residuals(t,r,res,sig)
!!$    c = s - stukey
!!$    testmnt = -sum(itukey(res/(s*sig))) - sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$!    testmnt = -sum((res/(s*sig))**2/2) - sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$
!!$  end function testmnt


!!$  subroutine ropt1(t,r,dr,s,ds,reli)
!!$
!!$    use rfun
!!$    use minpacks
!!$    use NelderMead
!!$!    use fmm
!!$!    use likescale
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), intent(in out) :: r,s
!!$    real(REAL64), intent(out) :: dr,ds
!!$!    real, intent(in) :: reps
!!$    logical, intent(out) :: reli
!!$
!!$    character(len=80) :: errmsg
!!$    real(REAL64), parameter :: eps = epsilon(t)
!!$    real(REAL64), dimension(2) :: x, dx
!!$    real(REAL64), dimension(2,2) :: cov
!!$!    integer, parameter :: nprint = 1
!!$    real(REAL64), parameter :: tol = sqrt(epsilon(tol))
!!$    real(REAL64) :: rtol,stol,rmax,rmin, lambda0, f, f0, c, s0
!!$    integer :: iter, info, nprint, stat, n
!!$    logical :: converge
!!$
!!$    n = size(res)
!!$    rmax = maxval(abs(cts - t*pht))
!!$!    rmax = 5*r
!!$!    rmin = minval(cts + t**2*pht)
!!$    rmin = 10*sqrt(n*epsilon(r))
!!$
!!$!    if( debug ) write(error_unit,'(a,2(2x,g0))') &
!!$!         'ropt: rmin, rmax: ',real(rmin),real(rmax)
!!$
!!$!    write(error_unit,*) r,s
!!$!    r = 9600
!!$
!!$!    rtol = reps * r
!!$    if( abs(r) > 10*n*sqrt(epsilon(r)) ) then
!!$       rtol = 1e-3*r
!!$    else
!!$       rtol = 1e-3*rmax
!!$    end if
!!$    stol = 1e-3
!!$    global(1) = t
!!$    nprint = 0
!!$
!!$    allocate(dpsi(n),stat=stat,errmsg=errmsg)
!!$    ! Computations needs some additional memory,
!!$    ! a fail-back way should be used if it is unavailable.
!!$    if( stat /= 0 ) then
!!$       if( debug ) write(error_unit,*) warning_label,trim(errmsg),' (ropt)'
!!$       return
!!$    end if

!!$    block
!!$      r = 400
!!$      s = 0.9
!!$      mu = 1
!!$      lambda = 1
!!$      call testder(t,s,r)
!!$      stop 0
!!$    end block

!!$
!!$    call residuals(t,r,res,sig)
!!$    c = s - stukey
!!$!    write(*,*) sum(itukey(res/(s*sig))), sum(log(s*sig)), c**2/2
!!$!    stop 0
!!$
!!$    if( abs(c) > 0 ) then
!!$       mu = sum(itukey(res/(s*sig))) / (10*c**2)
!!$    else
!!$       mu = sum(itukey(res/(s*sig))) / 10
!!$    end if
!!$    if( mu < 0.1 ) mu = 1
!!$
!!$!    mu = 16
!!$!    lambda = mu*(s - shuber)
!!$    lambda = mu*(s - 1)
!!$    lambda0 = huge(lambda)
!!$    f0 = huge(f)
!!$!    s = stukey
!!$    converge = .false.
!!$    do iter = 1, 3*precision(r)
!!$!       u0 = u
!!$!       x = [ r, s, mu, lambda ]
!!$       x = [ r, s ]
!!$       call lmder3(fder,x,cov,tol,666,1d0,nprint,info)
!!$         !       call lmder2(fder,x,tol,nprint,info)
!!$       !       call lmdif2(fder1,x,tol,nprint,info)
!!$
!!$       if( any(x < 0) .or. cov(1,1) < eps .or. cov(2,2) < eps .or. &
!!$            info == 666 .or. info < 0 ) then
!!$          if( debug_iter ) write(error_unit,*) &
!!$               'lmder3 failed at ',real(x),'for',info,' or reached a minimum.'
!!$          x = [ r, s ]
!!$          dx = [ 10*rtol, 10*stol ]
!!$          call nelmin1(rabsfun,x,dx,f,info)
!!$          if( debug_iter ) write(error_unit,*) 'nelmin: ',real(x),info
!!$          if( info == 2 ) info = -1 ! no convergence
!!$       end if
!!$
!!$       if( debug_iter ) then
!!$
!!$!          block
!!$!            real(REAL64) :: c
!!$!            integer :: n
!!$
!!$!            n = size(res)
!!$          c = s - stukey
!!$
!!$          if( s > 0 ) then
!!$             f = sum(itukey(res/(s*sig)))+sum(log(s*sig)) - lambda*c + mu*c**2/2
!!$             !               f = sum(itukey(res/(s*sig))) + n*log(s) - lambda*c + mu*c**2/2
!!$             !               f = sum(ihuber(res/(s*sig)))+sum(log(s*sig))-lambda*c+mu*c**2/2
!!$          else
!!$             f = 0
!!$          end if
!!$            !         write(*,*) f
!!$!          end block
!!$
!!$          write(error_unit,'(a,i2,2x,f0.1,5x,f6.4,3x,3(g0,3x),i1)') &
!!$               'ropt: ',iter,real(x),real(lambda),real(mu),real(f),info
!!$       end if
!!$
!!$!       write(error_unit,*) x(1),rmax,s,rmin,r
!!$       if( (iter > 1 .and. (x(1) > rmax .or. s < 0.1 .or. s > 2 ))  &
!!$            .or. info < 0 ) then
!!$          converge = .false.
!!$          exit
!!$       end if
!!$
!!$!       if( cov(1,1) < epsilon(cov) .or. cov(2,2) < epsilon(cov) ) then
!!$!          converge = .false.
!!$!          exit
!!$!       end if
!!$
!!$       converge = (abs(x(1) - r) < rtol .and. abs(x(2) - stukey) < stol) .or. &
!!$            abs(f - f0) < 10*n*sqrt(epsilon(f))
!!$!            abs(lambda0 - lambda) < stol) .or. r < rmin .or. &
!!$       r = x(1)
!!$       s = x(2)
!!$
!!$       call residuals(t,r,res,sig)
!!$       if( s > 0 .and. cov(1,1) > 0 .and. cov(2,2) > 0 ) then
!!$          s0 = s**2*sum(tukey(res/sig/s)**2)/ &
!!$               sum(dtukey(res/sig/s)**2)*n**2/(n-1)
!!$          if( s0 > 0 ) then
!!$             dr = sqrt(s0*cov(1,1)/(n-2))
!!$             ds = sqrt(s0*cov(2,2)/(n-2))
!!$             rtol = max(0.01*dr,rtol)
!!$             stol = max(0.01*ds,stol)
!!$          end if
!!$       else
!!$          dr = -1
!!$          ds = -1
!!$       end if



!!$       if( cov(1,1) > 0 .and. cov(2,2) > 0 ) then
!!$       block
!!$!         integer :: flag
!!$         real(REAL64) :: s0, ds, dr
!!$         real(REAL64), dimension(2) :: fvec
!!$         real(REAL64), dimension(2,2) :: fjac
!!$
!!$!         call residuals(t,r,res,sig)
!!$!         write(error_unit,*) s0
!!$!         flag = 2
!!$!         call fder(2,2,x,fvec,fjac,2,flag)
!!$!         call qrinv(fjac,cov)
!!$!         write(error_unit,*) cov(1,:)/sqrt(cov(1,1)*cov(2,2))
!!$!         write(error_unit,*) cov(2,:)/sqrt(cov(1,1)*cov(2,2))
!!$
!!$       end block
!!$    end if

!!$       if( converge ) exit
!!$!       r = x(1)
!!$!       s = x(2)
!!$!       lambda = x(4)
!!$       !       mu = x(3)
!!$       f0 = f
!!$       lambda0 = lambda
!!$       lambda = lambda - mu*(s - stukey)
!!$!       lambda = lambda - mu*(s - shuber)
!!$       mu = 2*mu
!!$    end do
!!$
!!$!    if( r < 0 ) r = 0
!!$
!!$    r = abs(r)
!!$    reli = converge !mu > 0
!!$    deallocate(dpsi)
!!$!    stop 0
!!$
!!$    if( reli .and. debug_iter ) then
!!$       block
!!$         real(REAL64), dimension(2,2) :: cor
!!$         real(REAL64) :: q
!!$
!!$         if( cov(1,1) > 0 .and. cov(2,2) > 0 ) then
!!$            q = sqrt(cov(1,1)*cov(2,2))
!!$            cor = 1
!!$            cor(1,2) = cov(1,2) / q
!!$            cor(2,1) = cor(1,2)
!!$         else
!!$            cor = -1
!!$         end if
!!$         write(error_unit,'(a,2(g0,1x),t45,a,2(f7.4,1x))') &
!!$              'Covariance: ', real(cov(1,:)),'correlation: ', cor(1,:)
!!$         write(error_unit,'(11x,2(g0,1x),t45,13x,2(f7.4,1x))') &
!!$              real(cov(2,:)),cor(2,:)
!!$       end block
!!$    end if
!!$
!!$  end subroutine ropt1
!!$
!!$  subroutine fder(k,l,x,fvec,fjac,ldfjac,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k, l, ldfjac
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$    real(REAL64), dimension(ldfjac,k), intent(out) :: fjac
!!$
!!$!    real(REAL64), dimension(:), allocatable :: dpsi
!!$!    real(REAL64), parameter :: s0 = shuber
!!$!    real(REAL64), parameter :: s0 = stukey
!!$!    real(REAL64) :: t,r,s,mu, lambda, c
!!$    real(REAL64) :: t,r,s,c
!!$    integer :: n
!!$
!!$    t = global(1)
!!$    r = x(1)
!!$    s = x(2)
!!$ !   mu = x(3)
!!$ !   lambda = x(4)
!!$    n = size(pht)
!!$!    c = s - s0
!!$    c = s - stukey
!!$
!!$    if( .not. (s > 0) ) then
!!$       iflag = -1
!!$       if( debug_iter ) write(error_unit,*) &
!!$            'fder: violation of the condition: s > 0, encountered.'
!!$       return
!!$    end if
!!$
!!$    if( iflag == 0 .and. debug_iter ) then
!!$       write(error_unit,*) 'x = ',real(x)
!!$       write(error_unit,*) 'fvec = ',real(fvec)
!!$       write(error_unit,*) 'fjac(1,:) = ',real(fjac(1,:))
!!$       write(error_unit,*) 'fjac(2,:) = ',real(fjac(2,:))
!!$       return
!!$    end if
!!$
!!$!    if( .not. (0.1 < s .and. s < 2) ) then
!!$!       iflag = -1
!!$!       return
!!$!    end if
!!$
!!$    if( iflag == 1 ) then
!!$       !       nev = nev + 1
!!$       call residuals(t,r,res,sig)
!!$!       call tukeys(res/(s*sig),psi,dpsi)
!!$       psi = tukey(res/(s*sig))
!!$!       psi = huber(res/(s*sig))
!!$!       psi = res/(s*sig)
!!$       fvec(1) = -r*sum(psi*(res/sig)/sig**2)/s + sum((r/sig)/sig)
!!$       fvec(2) = -sum(psi*res/sig)/s**2 + n/s - lambda + mu*c
!!$!       fvec(3) = c**2/2
!!$!       fvec(4) = -c
!!$    end if
!!$
!!$    if( iflag == 2 ) then
!!$!       nev = nev + 1
!!$!       c = sum(x(1:2)**2) - 2
!!$       call residuals(t,r,res,sig)
!!$       call tukeys(res/(s*sig),psi,dpsi)
!!$!       psi = huber(res/(s*sig))
!!$!       dpsi = dhuber(res/(s*sig))
!!$!       psi = res/(s*sig)
!!$!       dpsi = 1
!!$       fjac = 0
!!$!       write(*,*) ':'
!!$!       write(*,*) -sum(psi*(res/sig)/sig**2)/s
!!$!       write(*,*) sum(dpsi*(((r*res)/sig)/sig**2)**2)/s
!!$!       write(*,*) 3/s*sum(psi*(res/sig)*(r/sig**2)**2)
!!$!       write(*,*) sum(1/sig**2)
!!$!       write(*,*) 2*sum((r/sig**2)**2)
!!$       fjac(1,1) = -sum(psi*(res/sig)/sig**2)/s + &
!!$            sum(dpsi*((res/sig)**2*(r/sig)**2)/sig**2)/s**2 &
!!$            + 3*sum(psi*((res/sig)*(r/sig)**2/sig**2))/s &
!!$            +sum(1/sig**2) - 2*sum((r/sig**2)**2)
!!$       !+ sum(1/sig**2) - &
!!$!            2*sum((r/sig**2)**2)
!!$!            +3*r**2/s*sum(psi*res/sig**5) + sum(1/sig**2) - 2*r**2*sum(1/sig**4)
!!$       fjac(1,2) = r*sum(psi*(res/sig)/sig**2)/s**2 + &
!!$            r*sum(dpsi*(res/sig)/sig**3)/s**3
!!$       fjac(2,2) = 2/s**3*sum(psi*(res/sig)) + sum(dpsi*(res/sig)**2)/s**4 &
!!$            -n/s**2 + mu
!!$       fjac(2,1) = fjac(1,2)
!!$ !      fjac(2,3) = c
!!$ !      fjac(2,4) = -1
!!$ !      do j = 2,4
!!$ !         fjac(j,1:j-1) = fjac(1:j-1,j)
!!$ !      end do
!!$    end if
!!$
!!$  end subroutine fder
!!$
!!$  subroutine fder1(k,l,x,fvec,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k, l
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(l), intent(out) :: fvec
!!$
!!$!    character(len=80) :: errmsg
!!$!    real(REAL64), dimension(:), allocatable :: dpsi
!!$!    real(REAL64), parameter :: s0 = 1!stukey
!!$    real(REAL64) :: t,r,s,c
!!$    integer :: n!, stat
!!$
!!$    t = global(1)
!!$    r = x(1)
!!$    s = x(2)
!!$!    mu = x(3)
!!$!    lambda = x(4)
!!$    n = size(pht)
!!$    c = s - shuber
!!$
!!$    if( iflag == 0 .and. debug_iter ) then
!!$       write(error_unit,*) 'x = ',real(x)
!!$       write(error_unit,*) 'fvec = ',real(fvec)
!!$       return
!!$    end if
!!$
!!$    if( .not. (0.5 < s .and. s < 1.5 .and. r >= 0) ) then
!!$       iflag = -1
!!$       return
!!$    end if
!!$
!!$
!!$!    allocate(dpsi(n),stat=stat,errmsg=errmsg)
!!$    ! Computations needs some additional memory,
!!$    ! a fail-back way should be used if it is unavailable.
!!$!    if( stat /= 0 ) then
!!$!       if( debug ) write(error_unit,*) warning_label,trim(errmsg)
!!$!       return
!!$!    end if
!!$
!!$    !       nev = nev + 1
!!$    call residuals(t,r,res,sig)
!!$!    call tukeys(res/(s*sig),psi,dpsi)
!!$!    psi = tukeys(res/(s*sig))
!!$    psi = huber(res/(s*sig))
!!$!    psi = res/(s*sig)
!!$    fvec(1) = -r*sum(psi*res/sig**3)/s + r*sum(1/sig**2)
!!$    fvec(2) = -sum(psi*res/sig)/s**2 + n/s - lambda + mu*c
!!$!    fvec(3) = c**2/2
!!$!    fvec(4) = -c
!!$
!!$  end subroutine fder1
!!$
!!$  subroutine fder1test(k,x,fvec,iflag)
!!$
!!$    use rfun
!!$
!!$    integer, intent(in) :: k
!!$    integer, intent(in out) :: iflag
!!$    real(REAL64), dimension(k), intent(in) :: x
!!$    real(REAL64), dimension(k), intent(out) :: fvec
!!$
!!$    real(REAL64) :: t,r,s,c
!!$    integer :: n
!!$
!!$    t = global(1)
!!$    r = x(1)
!!$    s = x(2)
!!$    n = size(pht)
!!$    c = s - shuber
!!$
!!$    if( iflag == 0 .and. debug_iter ) then
!!$       write(error_unit,*) 'x = ',real(x)
!!$       write(error_unit,*) 'fvec = ',real(fvec)
!!$       return
!!$    end if
!!$
!!$    call residuals(t,r,res,sig)
!!$!    call tukeys(res/(s*sig),psi,dpsi)
!!$!    psi = tukeys(res/(s*sig))
!!$    psi = huber(res/(s*sig))
!!$!    psi = res/(s*sig)
!!$    fvec(1) = -r*sum(psi*res/sig**3)/s + r*sum(1/sig**2)
!!$    fvec(2) = -sum(psi*res/sig)/s**2 + n/s - lambda + mu*c
!!$
!!$  end subroutine fder1test
!!$
!!$
!!$  subroutine testder(t,r,s)
!!$
!!$    real(REAL64), intent(in) :: t,r,s
!!$    real(REAL64), dimension(2) :: fd,fvec,eps,x,x1,x2
!!$    real(REAL64), dimension(2,2) :: fjac, djac
!!$    real(REAL64) :: f0
!!$    integer :: i, j, iflag
!!$
!!$    eps = 1e-2
!!$!    eps(1) = 1e-2
!!$    global(1) = t
!!$
!!$    iflag = 1
!!$!    mu = 2
!!$!    lambda = 1
!!$!    x = [ r, s, mu, lambda ]
!!$    x = [ r, s ]
!!$    call fder(2,2,x,fvec,fjac,2,iflag)
!!$    f0 = testropt(x)
!!$    do i = 1, 2
!!$!       x = [ r, s, mu, lambda ]
!!$       x = [ r, s ]
!!$       x(i) = x(i) + eps(i)
!!$       fd(i) = (testropt(x) - f0) / eps(i)
!!$       write(error_unit,*) i,real(fvec(i)),real(fd(i)),real((fd(i)-fvec(i)))
!!$    end do
!!$
!!$!    return
!!$
!!$
!!$    iflag = 2
!!$    x = [ r, s ]
!!$    call fder(2,2,x,fvec,fjac,2,iflag)
!!$    call fdjac1(fder1test,2,x,fvec,djac,2,iflag,2,2,epsilon(x),x1,x2)
!!$    do i = 1,2
!!$       do j = i, 2
!!$!          x = [ r, s, mu, lambda ]
!!$!          x1(i) = x(i) + eps(i)
!!$!          f1 = testropt(x1)
!!$!          x1(i) = x(i) + eps(i)
!!$!          x2(j) = x(j) - eps(j)
!!$!          jd(i,j) = (testropt(x1) + testropt(x2) - 2*f0) / (eps(i)*eps(j))
!!$          write(error_unit,'(2i5,3(2x,g0))') i,j,real(fjac(i,j)), &
!!$               real(djac(i,j)), real((fjac(i,j)-djac(i,j)))!/fjac(i,j))
!!$       end do
!!$    end do
!!$
!!$
!!$  end subroutine testder
!!$
!!$
!!$  function testropt(x)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64) :: t,r,s,testropt,c
!!$    integer :: n
!!$
!!$    n = size(res)
!!$    t = global(1)
!!$    r = x(1)
!!$    s = x(2)
!!$!    mu = x(3)
!!$!    lambda = x(4)
!!$    call residuals(t,r,res,sig)
!!$    c = s - 1!stukey
!!$!    testropt = sum(itukey(res/(s*sig))) + sum(log(s*sig))-lambda*c + mu*c**2/2
!!$!    testropt = sum(itukey(res/(s*sig))) + n*log(s) - lambda*c + mu*c**2/2
!!$    testropt = sum((res/(s*sig))**2)/2 + n*log(s) - lambda*c + mu*c**2/2
!!$
!!$  end function testropt
!!$
!!$
!!$  subroutine ropt2(t,r,s,reps,reli)
!!$
!!$    use NelderMead
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), intent(in out) :: r,s
!!$    real, intent(in) :: reps
!!$    logical, intent(out) :: reli
!!$
!!$    real(REAL64), dimension(2) :: x, dx
!!$    real(REAL64) :: rmin
!!$    integer :: ifault
!!$
!!$    x = [ r, s ]
!!$    dx = [ max(r,1.0d0), 1d-2 ]
!!$    call nelmin1(rabsfun,x,dx,rmin,ifault)
!!$    reli = ifault == 0
!!$
!!$  end subroutine ropt2
!!$
!!$  function rabsfun(x) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), dimension(:), intent(in) :: x
!!$    real(REAL64) :: t,r,s,c,f
!!$
!!$    t = global(1)
!!$    r = x(1)
!!$    s = x(2)
!!$    c = s - stukey
!!$
!!$    if( .not. (0.1 < s .and. s < 2) ) then
!!$       f = mu*c**4
!!$       return
!!$    end if
!!$
!!$    call residuals(t,r,res,sig)
!!$    f = sum(itukey(res/(s*sig))) + sum(log(s*sig)) - lambda*c+mu*c**2/2
!!$!    write(*,*) 'rabsfun:',real(x),real(f)
!!$
!!$  end function rabsfun

!!$  subroutine ropt3(t,r,s,reps,reli)
!!$
!!$    use fmm
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), intent(in out) :: r,s
!!$    real, intent(in) :: reps
!!$    logical, intent(out) :: reli
!!$
!!$    real(REAL64), dimension(2) :: x, dx
!!$    real(REAL64) :: tol,rmax,rmin,rhalf
!!$    integer :: ifault
!!$
!!$    rmin = 0
!!$    rhalf = maxval(abs(cts - t*pht)) / 2 !/ t
!!$!    rmin = max(r-rhalf,0*r)
!!$    rmax = r + 3*rhalf
!!$    tol = max(reps*rmax,10*sqrt(epsilon(r)))
!!$
!!$    global(1) = t
!!$    global(2) = s
!!$    r = fmin(rmin,rmax,roptmin,tol)
!!$
!!$  end subroutine ropt3
!!$
!!$  function roptmin(r) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: r
!!$    real(REAL64) :: t,s,f
!!$
!!$    t = global(1)
!!$    s = global(2)
!!$    call residuals(t,r,res,sig)
!!$    f = sum(itukey(res/(s*sig))) + sum(log(s*sig))
!!$    write(*,*) 'roptmin:',real(r),real(f)
!!$
!!$  end function roptmin



!!$  subroutine ropt(t,r,s,reps,reli)
!!$
!!$    use fmm
!!$    use likescale
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), intent(in out) :: r,s
!!$    real, intent(in) :: reps
!!$    logical, intent(out) :: reli
!!$
!!$    real(REAL64) :: tol,rmax,rmin,rhalf
!!$    integer :: m
!!$
!!$    rmin = max(r / 3, 0*r)
!!$    rhalf = maxval(abs(cts - t*pht)) / 2 !/ t
!!$    rmin = max(r-rhalf,0*r)
!!$    rmax = r + rhalf
!!$!    rmin = 2000
!!$!    rmax = 5000
!!$!    rmin = r / 2
!!$!    rmax = 2 * r
!!$    tol = max(reps*rmax,10*sqrt(epsilon(r)))
!!$!    write(*,*) rmin, rmax, tol
!!$    global(1) = t
!!$!    r = fmin(rmin,rmax,rfitfun,tol)
!!$!    block
!!$!      real(REAL64) :: s
!!$!      call residuals(t,r,res,sig)
!!$!      call iscale(res/sig,s)
!!$!      write(*,*) 's = ', real(s)
!!$!    end block
!!$!    return
!!$
!!$
!!$!    block
!!$!      integer :: i
!!$!      do i = 1,100
!!$!         r = i*rmax/100
!!$
!!$
!!$    do m = 1, 3
!!$       r = zeroin(rmin,rmax,rfitfun1,tol)
!!$       reli = r < rmax - 2*tol .and. global(2) > 0
!!$       if( reli ) exit
!!$       if( r > rmax - 2*tol ) rmax = rmax + rhalf
!!$       if( r < rmin + 2*tol ) rmin = max(rmin - rhalf,0*r)
!!$       if( debug_iter ) write(error_unit,'(a,2(g0,2x))') &
!!$            'updated rmin,rmax: ',real(rmin),real(rmax)
!!$    end do
!!$    reli = reli .and. r < rmax - 2*tol .and.  global(2) > 0
!!$
!!$    if( global(2) > 0 ) s = global(2)
!!$
!!$    if( debug_iter ) write(error_unit,'(a,2(g0,2x),l1)') &
!!$         'ropt: ',real(r),real(s),reli
!!$
!!$!    reli = .false.
!!$
!!$!    if( rfitfun1(rmin)*rfitfun1(rmax) < 0 ) then
!!$    !    else
!!$    if( .not. reli .and. debug ) write(error_unit,'(2a)') &
!!$         warning_label,'ropt: failed to determine optimal r'
!!$!       r = rmax / 2
!!$!    end if
!!$
!!$!    stop 0


!!$    block
!!$      real :: rr
!!$      integer :: i
!!$      do i = 1, size(cts)
!!$         rr = (cts(i) - t*pht(i))
!!$         write(*,*)  rr,  sqrt(cts(i) + t**2*pht(i) + r**2), & !r / sqrt(cts(i) + t**2*pht(i)), &
!!$              sqrt(cts(i) + t**2*pht(i))
!!$!              (r**2-(cts(i) + t**2*pht(i)))
!!$      end do
!!$!      stop 0
!!$    end block


!!$  end subroutine ropt
!!$
!!$  function rfitfun1(r)
!!$
!!$    use rfun
!!$    use likescale
!!$
!!$    real(REAL64), intent(in) :: r
!!$    real(REAL64) :: t, rfitfun1, s, s0, tol
!!$    integer :: n, iter, niter
!!$    logical :: converge
!!$
!!$    t = global(1)
!!$    n = size(res)
!!$    niter = 0

!    call iscale(res/sig,s)
    !    rfitfun1 = sum(psi*res/sig**2) - n*s
!    s = sqrt(sum((psi)**2)/n)
    !    call escale(t,r,s,converge)
!!$    s = 1
!!$    tol = 1e-2
!!$    do iter = 1, 0*2*precision(s)
!!$       s0 = s
!!$       !       s = sum(psi*res/sig/s)/n
!!$       psi = huber(res/sig/s)
!!$       s = sqrt(sum(psi**2*s**2)/n/shuber)
!!$       write(*,*) iter, s, s0
!!$       if( abs(s - s0) < tol) exit
!!$    end do
!!$    rfitfun1 = s - shuber


!    call escale(t,r,s,converge)
!!$    converge = .false.
!!$    if( converge ) then
!!$       continue
!!$!       rfitfun1 = s - stukey
!!$    else ! if( .not. converge ) then
!!$       s = 1!shuber
!!$       s0 = -1
!!$       tol = 1e-3
!!$       converge = .false.
!!$       call residuals(t,r,res,sig)
!!$       do iter = 1, 3*precision(s)
!!$          psi = huber(res/(s*sig))
!!$!          psi = tukey(res/(s*sig))
!!$!          s = sqrt(sum(psi**2*s**2)/n/shuber)
!!$          !       s = sqrt(sum((res/sig)**2)/n)
!!$          s = sum(psi*res/sig)/n
!!$!          s = sqrt(sum(psi**2)/n)*s/stukey
!!$          !          write(*,*) iter, s, s0
!!$          if( .not. (0.1 < s .and. s < 10) ) exit
!!$
!!$          if( abs(s - s0) < -0.01 ) then
!!$             s0 = s
!!$             call escale(t,r,s,converge)
!!$             if( .not. converge ) then
!!$                s = s0
!!$                converge = .true.
!!$             end if
!!$          end if
!!$          converge = abs(s - s0) < tol
!!$          if( converge ) exit
!!$          niter = iter
!!$          s0 = s
!!$       end do
!!$       if( converge ) then
!!$          rfitfun1 = s - shuber
!!$!          rfitfun1 = s - 1!stukey
!!$          global(2) = s
!!$       else
!!$          rfitfun1 = s
!!$          global(2) = -1
!!$       end if
!!$       !       rfitfun1 = s - stukey
!!$!       if( converge ) then
!!$!          rfitfun1 = s - shuber
!!$!       else
!!$!          rfitfun1 = 0.1*huge(s)
!!$!       end if
!!$    end if
!!$
!!$    ! converge test?
!!$!    if( .not. converge ) then
!!$!       s = -1
!!$!       write(*,*) '
!!$!    end if
!!$!    if( converge ) global(2) = s
!!$
!!$    if( debug_iter ) write(error_unit,'(a,3(g0,2x),l2,i3)') &
!!$         'rfitfun1: ',real(r),real(rfitfun1),real(s), converge,niter
!!$    !,sqrt(sum((psi)**2)/n)
!!$!stop 0
!!$  end function rfitfun1
!!$
!!$  subroutine escale(t,r,s,converge)
!!$
!!$    use rfun
!!$    use likescale
!!$
!!$    real(REAL64), intent(in) :: t,r
!!$    real(REAL64), intent(out) :: s
!!$    logical, intent(out) :: converge
!!$
!!$    character(len=80) :: errmsg
!!$    real(REAL64), dimension(:), allocatable :: dpsi, rs
!!$    real(REAL64) :: tol,f,df,d,ss
!!$    integer :: n,iter,stat
!!$
!!$    converge = .false.
!!$    tol = 1e-5
!!$    s = stukey
!!$    n = size(res)
!!$
!!$    allocate(dpsi(n),rs(n),stat=stat,errmsg=errmsg)
!!$    ! Computations needs some additional memory,
!!$    ! a fail-back way should be used if it is unavailable.
!!$    if( stat /= 0 ) then
!!$       if( debug ) write(error_unit,*) warning_label,trim(errmsg)
!!$       return
!!$    end if
!!$
!!$    do iter = 1, precision(r)
!!$
!!$       if( .not. (0.1 < s .and. s < 10) ) return
!!$
!!$       call residuals(t,r,res,sig)
!!$       rs = res/sig
!!$!       call tukeys(res/(s*sig),psi,dpsi)
!!$       call tukeys(rs/s,psi,dpsi)
!!$!       psi = res / (s*sig)
!!$!       dpsi = 1
!!$!       psi = huber(res / (s*sig))
!!$!       dpsi = dhuber(res / (s*sig))
!!$
!!$!       f = sum(psi*res/sig)/s - n
!!$!       df = -(sum(dpsi*res**2/sig**2) + sum(psi*res/sig))/s**2
!!$
!!$!       ss = sum(psi*res/sig)
!!$       ss = sum(psi*rs)
!!$       f = (ss/s - n) / s
!!$!       df = -((ss/s - n) + (sum(dpsi*(res/sig)**2) + ss))/s**2
!!$       df = -((ss/s - n) + (sum(dpsi*rs**2) + ss))/s**2
!!$!       f = f / s
!!$
!!$       block
!!$         real(REAL64), dimension(:), allocatable :: psi1
!!$         real(REAL64) :: e,f1,fd
!!$         allocate(psi1(n))
!!$         e = 1e-7
!!$         psi1 = res / ((s+e)*sig)
!!$         call tukeys(res/((s+e)*sig),psi1,dpsi)
!!$         f1 = (sum(psi1*res/sig)/(s+e) - n) / (s+e)
!!$         fd = (f1 - f)/e
!!$         write(*,*) f,f1,fd,df,(fd-df)/fd
!!$         stop 0
!!$       end block

!       if( .not. (abs(df) < epsilon(df)) ) exit

!!$       ! scale update
!!$       d = f / df
!!$       s = s - d
!!$
!!$       if( debug_iter ) write(error_unit,'(a,i2,1pg20.6,1p3g10.2)') &
!!$            "#, s, inc., f, f' *: ",iter,s,d,-f,-df
!!$
!!$       ! exit of iterations: the absolute error must be at least |d| < tol
!!$       converge = abs(d) < tol
!!$       if( converge ) exit
!!$
!!$    enddo
!!$
!!$!    stop 0
!!$
!!$    ! second derivation of likelog is negative for a stable local minimum
!!$    converge = converge .and. df < 0
!!$
!!$
!!$  end subroutine escale
!!$
!!$
!!$
!!$
!!$
!!$  subroutine rquant(t,reff,reli)
!!$
!!$    use quantilmean
!!$
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), intent(out) :: reff
!!$    logical, intent(out) :: reli
!!$
!!$    real(REAL64), parameter :: qhalf = 0.6745, qmax = 0.87
!!$    real(REAL64), dimension(:), allocatable :: xcdf, ycdf
!!$    character(len=80) :: errmsg
!!$    integer :: n, stat
!!$
!!$    reli = .false.
!!$    ! determine quantile of residuals
!!$    n = size(cts)
!!$    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
!!$    if( stat /= 0 ) then
!!$       write(error_unit,*) warning_label,trim(errmsg),' (rquant)'
!!$       return
!!$    end if
!!$
!!$    call ecdf(abs(cts - t*pht),xcdf,ycdf)
!!$
!!$    block
!!$      real(REAL64) :: q
!!$      integer :: i
!!$      do i = 1,10
!!$         q = 0.1*i
!!$         write(error_unit,*) q, quantile(q,xcdf,ycdf)
!!$      end do
!!$      q = 0.7361
!!$      write(error_unit,*) q, quantile(q,xcdf,ycdf)
!!$
!!$!      open(1,file='/tmp/rcdf')
!!$!      do i = 1, n
!!$!         write(1,*) xcdf(i),ycdf(i)
!!$!      end do
!!$!      close(1)
!!$    end block
!!$
!!$    ! quantile for Q(0.87) for 1.5*Q(0.67), gives upper bound
!!$    if( qmax < ycdf(n) ) then
!!$       reff = quantile(qmax,xcdf,ycdf)
!!$    else ! if( Q(0.87) is out of CDF values ) then
!!$       reff = 5 * quantile(qhalf,xcdf,ycdf)  ! = Q(4/5) * 1.5
!!$    end if
!!$!       write(*,* ) qmax < ycdf(n), reff, quantile(qmax,xcdf,ycdf), quantile(qhalf,xcdf,ycdf)
!!$    deallocate(xcdf, ycdf)
!!$!    else
!!$!       reff = 5 * median(abs(cts - t*pht)) / Qdev
!!$!    end if
!!$    reli = .true.
!!$
!!$  end subroutine rquant

!!$  subroutine rbracket1(t,r,reli)
!!$
!!$    use medians
!!$    use fmm
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), intent(out) :: r
!!$    logical, intent(out) :: reli
!!$
!!$
!!$    real(REAL64) :: tol,rmax,rmin,reff
!!$    integer :: i, n
!!$
!!$    global(1) = t
!!$
!!$    call rquant(t,reff,reli)
!!$    if( .not. reli ) then
!!$       reff = median(abs(cts - t*pht)) / Qdev
!!$    end if
!!$
!!$!    write(*,*) median(abs(cts - t*pht)) / Qdev
!!$    rmin = 0
!!$    tol = reff / 666
!!$    !    rmax = maxval(abs(cts - t*pht))
!!$    !    r = maxval(abs(cts - t*pht))
!!$    reli = .false.
!!$
!!$    do i = 1, 3
!!$       rmax = i*reff
!!$       !    call bisection(rmin,rmax,tol,r,reli)
!!$       if( debug_iter ) then
!!$          write(error_unit,'(a,i1,a,2(2x,g0))') &
!!$               'rbracket rmin, rmax: ',i,'-sigma, ',real(rmin),real(rmax)
!!$          if( i > 3 ) write(error_unit,'(a,2(2x,g0))') &
!!$               'updated rmin, rmax: ',real(rmin),real(rmax)
!!$       end if
!!$
!!$       r = zeroin(rmin,rmax,bifun,tol)
!!$       reli = r < rmax - 2*tol
!!$       if( reli ) exit
!!$    end do
!!$
!!$    if( debug_graph ) then
!!$       block
!!$         use quantilmean
!!$
!!$         real(REAL64), dimension(:), allocatable :: x,p
!!$         real(REAL64) :: rr
!!$         integer :: i,n
!!$
!!$         n = size(res)
!!$         allocate(x(n),p(n))
!!$
!!$         !      r = 6633
!!$         call residuals(t,r,res,sig)
!!$         call ecdf(abs(res)/sig,x,p)
!!$         open(1,file='/tmp/bicdf')
!!$         do i = 1, n
!!$            write(1,*) x(i),p(i)
!!$         end do
!!$         close(1)
!!$
!!$!         open(1,file='/tmp/bifun')
!!$!         do i = 0, 1000
!!$!            rr = rmin + i*(rmax - rmin)/1000.0
!!$!            write(1,*) rr,bifun(rr)
!!$!         end do
!!$!         close(1)
!!$       end block
!!$    end if
!!$
!!$  end subroutine rbracket1

!!$  subroutine bisection(xmin,xmax,tol,x,converge)
!!$
!!$    real(REAL64), intent(in) :: xmin, xmax, tol
!!$    real(REAL64), intent(out) :: x
!!$    logical, intent(out) :: converge
!!$
!!$    real(REAL64) :: fx, fa, fb, a, b
!!$    integer :: n!, nmax
!!$
!!$    converge = .true.
!!$!    tol = 1e-3*xmax
!!$!    nmax = int(abs(log(tol))/log(2.0) + 3)
!!$!    nmax = 32
!!$!    write(*,*) nmax,xmax,xmin,abs(log(tol))/log(2.0)+1
!!$
!!$    a = xmin
!!$    b = xmax
!!$    fa = bifun(a)
!!$    fb = bifun(b)
!!$    do n = 1, 16
!!$
!!$       x = (a + b) / 2
!!$       fx = bifun(x)
!!$       if( fa*fx < 0 ) then
!!$          b = x
!!$          fb = fx
!!$       else !if( fx*fb < 0 ) then
!!$          a = x
!!$          fa = fx
!!$       end if
!!$
!!$       if( abs(b - a) < tol ) return
!!$
!!$    end do
!!$
!!$    converge = .false.
!!$
!!$  end subroutine bisection
!!$
!!$  function bifun(r)
!!$
!!$    use quantilmean
!!$    use medians
!!$
!!$    real(REAL64), intent(in) :: r
!!$    real(REAL64) :: bifun, t
!!$    real(REAL64), dimension(2), parameter :: q = [0.25, 0.75]
!!$!    real(REAL64), dimension(2), parameter :: q = [0.4, ]
!!$    real(REAL64), dimension(2) :: x
!!$    real(REAL64) :: s
!!$    integer :: iflag
!!$
!!$    t = global(1)
!!$    call residuals(t,r,res,sig)
!!$    call quantilefun(res/sig,q,x,iflag)
!!$!    iflag= 1
!!$    if( iflag == 0 ) then
!!$       ! standard deviation is being quantile deviation
!!$!       write(*,*) x(1),x(2)
!!$       s = (x(2) - x(1)) / (2*Qdev)
!!$    else
!!$       s = median(abs(res)/sig) / Qdev
!!$    end if
!!$    bifun = s - 1
!!$
!!$    if( debug_iter ) write(error_unit,*) 'bifun:',real(bifun),real(r),iflag
!!$
!!$  end function bifun

!!$  subroutine rbracket(t,r,reli)
!!$
!!$    use fmm
!!$    use likescale
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), intent(out) :: r
!!$    logical, intent(out) :: reli
!!$!    real, intent(in) :: reps
!!$
!!$    real(REAL64) :: tol,rmax,rmin
!!$    integer :: m
!!$
!!$    rmin = 0
!!$    rmax = maxval(abs(cts - t*pht)) !/ t
!!$    tol = rmax / 100.0
!!$    global(1) = t
!!$
!!$    do m = 1, 2
!!$       r = zeroin(rmin,rmax,bracketfun,tol)
!!$       reli = r < rmax - 2*tol
!!$       if( reli ) exit
!!$       rmax = 10*rmax
!!$       if( debug_iter ) write(error_unit,*) &
!!$            warning_label,' updated rmax=',real(rmax)
!!$    end do
!!$
!!$    if( debug_iter ) write(error_unit,*) 'rbracket: ',real(r),reli
!!$    if( .not. reli .and. debug ) write(error_unit,*) &
!!$         warning_label, 'rbracket: failed to determine optimal r'
!!$
!!$  end subroutine rbracket
!!$
!!$
!!$  function bracketfun(r)
!!$
!!$    use quantilmean
!!$
!!$    real(REAL64), intent(in) :: r
!!$    real(REAL64) :: bracketfun, t
!!$    real(REAL64), dimension(2), parameter :: q = [0.25, 0.75]
!!$    real(REAL64), dimension(2) :: x
!!$    integer :: iflag
!!$
!!$    t = global(1)
!!$    call residuals(t,r,res,sig)
!!$    call quantilefun(res/sig,q,x,iflag)
!!$    if( iflag == 0 ) then
!!$       ! standard deviation is being quantile deviation
!!$       bracketfun = (x(2) - x(1)) / (2*Qdev) - 1
!!$    else
!!$       bracketfun = 0
!!$    end if
!!$
!!$    if( debug_iter ) &
!!$         write(error_unit,*) 'bracketfun:',real(bracketfun),real(r),iflag
!!$
!!$  end function bracketfun

!!$  subroutine walleyofabs(r,d,tmin1,tmax1,tol1,t,reli)
!!$
!!$    use fmm
!!$
!!$    real(REAL64), intent(in) :: r,d,tmin1,tmax1,tol1
!!$    real(REAL64), intent(in out) :: t
!!$    logical, intent(out) :: reli
!!$    real(REAL64) :: h, t0, tmin, tmax, tol
!!$    integer :: i
!!$
!!$    global(1) = r
!!$    tol = d / 100
!!$    h = 2*tol
!!$    t0 = t
!!$
!!$    do i = 1, 5
!!$       tmin = t0 - i*d
!!$       tmax = t0 + i*d
!!$       if( debug_iter ) then
!!$          write(error_unit,'(a,i1,a,2(2x,g0))') &
!!$               'walleyofabs tmin,tmax: ',i,'-sigma, ',real(tmin),real(tmax)
!!$          if( i > 1 ) write(error_unit,'(a,2(2x,g0))') &
!!$               'updated tmin, tmax: ',real(tmin),real(tmax)
!!$       end if
!!$       t = fmin(tmin,tmax,walleyabs,tol)
!!$       reli = tmin + h < t .and. t < tmax - h
!!$       if( reli ) exit
!!$    end do
!!$
!!$    ! Realiable results are inside the range tmin .. tmax,
!!$    ! ones aren't so close to the interval endpoints.
!!$
!!$  end subroutine walleyofabs
!!$
!!$  function walleyabs(t)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64) :: walleyabs, r
!!$    integer :: n
!!$
!!$    r = global(1)
!!$    n = size(res)
!!$!    res = cts - t*pht
!!$    call residuals(t,r,res,sig)
!!$    walleyabs = sum(abs(res)/sig) / n + sum(log(sig)) / n
!!$
!!$    if( debug_iter ) write(error_unit,'(a,g20.10,3g15.5)') &
!!$         'walleyabs',t, walleyabs!, sum(abs(res)/sig)/n, sum(log(sig))/n
!!$
!!$  end function walleyabs

!!$  function rfitfun(r)
!!$
!!$    use rfun
!!$!    use quantilmean
!!$
!!$    real(REAL64), intent(in) :: r
!!$    real(REAL64) :: t, rfitfun, s0, d, w, h11
!!$    integer :: n
!!$
!!$    t = global(1)
!!$    call residuals(t,r,res,sig)
!!$    n = size(res)
!!$    w = real(n)/real(n-1)
!!$    ! fmin:
!!$!    rfitfun = abs(sum(ihuber(res/sig)) / size(res) - 0.6582**2)
!!$    rfitfun = sum(ihuber(res/sig)) + sum(log(sig))
!!$    rfitfun = sum((res/sig)**2)/2 + sum(log(sig))
!!$    rfitfun = rfitfun / n
!!$    write(*,*) r,rfitfun
!!$    return
!!$!    rfitfun = (sqrt(sum(huber(res/sig)**2)/sum(dhuber(res/sig))**2*n)-1)**2
!!$!    rfitfun = (sqrt(sum(tukey(res/sig)**2)/sum(dtukey(res/sig))**2*n)-1)**2
!!$
!!$    ! zeroin:
!!$    rfitfun = sum(huber(res/sig)*(res/sig)) - n*0.6582
!!$    rfitfun = sum((res/sig)**2) - n
!!$    rfitfun = sum(tukey(res/sig)*(res/sig)) - n*0.8061
!!$    write(*,*) r,rfitfun,sum((res/sig)**2),n*0.6582
!!$    return
!!$
!!$!    rfitfun = sqrt(sum((res/sig)**2) / n) - 1
!!$    !    rfitfun = sqrt(sum(huber(res/sig)**2)/sum(dhuber(res/sig))**2*n) - 1
!!$    s0 = sum(tukey(res/sig)**2)
!!$    d = sum(dtukey(res/sig))
!!$    s0 = sum(huber(res/sig)**2)
!!$    d = sum(dhuber(res/sig))
!!$    h11 = sum(dhuber(res/sig)*pht**2/sig**2)
!!$!    s0 = sum((res/sig)**2)
!!$!    d = n
!!$    if( s0 > 0 .and. d > 0 .and. h11 > 0 ) then
!!$       rfitfun = sqrt(s0/(d*h11)*w*n) - 1
!!$!       rfitfun = sqrt(s0/d) - 1
!!$!       rfitfun = (sum((res/sig)**2) -  n) / sqrt(2.0*n)
!!$!       rfitfun = (s0 -  d) / sqrt(2.0*d)
!!$    else
!!$       rfitfun = -1
!!$    end if
!!$
!!$
!!$!   rfitfun = -sum(ihuber(res/sig)) / size(res) - sum(log(sig)) / size(res)
!!$!    rfitfun = abs(sum(ihuber(res/sig)) - n) / sqrt(2.0*n)
!!$!    rfitfun = abs(sum(ihuber(res/sig)) / size(res) - 0.6582**2)
!!$!    rfitfun = abs(sum(itukey(res/sig)) / size(res) - 0.8061**2)
!!$!    rfitfun = abs(sum((res/sig)**2) / size(res) - 1)
!!$!    rfitfun = abs(sum((res/sig)**2) -  n) / sqrt(2.0*n)
!!$
!!$    write(*,*) r,rfitfun,s0,d,h11
!!$
!!$  end function rfitfun
!!$
!!$
!!$!  subroutine zmean(r,tmin1,tmax1,tol,t0,dt0,t,reli)
!!$  subroutine zmean(r,s,tol,t0,dt0,t,reli)
!!$
!!$    use fmm
!!$
!!$!    real(REAL64), intent(in) :: r,tmin1,tmax1,tol,t0,dt0
!!$    real(REAL64), intent(in) :: r,s,tol,t0,dt0
!!$    real(REAL64), intent(out) :: t
!!$    logical, intent(out) :: reli
!!$!    real(REAL64) :: h, tmin, tmax, maxfun, f, zmin, zmax
!!$    real(REAL64) :: h, tmin, tmax, f, zmin, zmax
!!$    integer :: i
!!$
!!$!    write(*,*) r,s,tol,t0,dt0
!!$
!!$    global(1) = r
!!$    global(2) = s
!!$    !    tol = 0.0001 * d
!!$
!!$
!!$!    maxfun = huge(maxfun)
!!$    do i = 1, 5
!!$       tmin = max(t0 - i*dt0,0.01*dt0)
!!$!       tmin = t0 - i*dt0
!!$       tmax = t0 + i*dt0
!!$       if( debug_iter ) then
!!$          write(error_unit,'(a,i1,a,2(2x,g0))') &
!!$               'zmean tmin,tmax: ',i,'-sigma, ',real(tmin),real(tmax)
!!$          if( i > 3 ) write(error_unit,'(a,2(2x,g0))') &
!!$               'updated tmin, tmax: ',real(tmin),real(tmax)
!!$       end if
!!$
!!$       ! Do first steps by bi-section method, the quick way to recognise
!!$       ! that the zero is out of given interval.
!!$       zmin = zfun1(tmin)
!!$       zmax = zfun1(tmax)
!!$       if( zmin*zmax < 0 ) then
!!$          t = (tmax + tmin) / 2
!!$          f = zfun1(t)
!!$          if( zmin*f < 0 ) then
!!$             tmax = t
!!$          else ! if( f*zmax < 0 ) then
!!$             tmin = t
!!$          end if
!!$          t = zeroin(tmin,tmax,zfun1,tol)
!!$          h = 2*tol
!!$          reli = tmin + h < t .and. t < tmax - h
!!$       else
!!$          reli = .false.
!!$       end if
!!$       if( reli ) exit
!!$!       f = abs(zfun1(t))
!!$!       if( f < maxfun ) then
!!$!          maxfun = f
!!$!       else
!!$!          exit
!!$!       end if
!!$    end do
!!$
!!$
!!$    ! a realiable result is inside range tmin .. tmax,
!!$    ! also not too close to the interval endpoints
!!$
!!$  end subroutine zmean

!!$  function zfun(t)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64) :: zfun, r, s
!!$    integer :: n
!!$
!!$    n = size(res)
!!$    r = global
!!$    s = scale
!!$    call residuals(t,r,res,sig)
!!$!    psi = huber(res/sig)
!!$    psi = tukey(res/sig/s)
!!$
!!$    !zfun = sum(psi*pht/sig) / size(res)
!!$    zfun = sum(psi*pht/sig) / (n*s) + &
!!$         sum(psi*res*cts/sig**3) / (n*s) + &
!!$         sum(cts/sig**2) / n
!!$!         sum((psi*res/sig + 1)*cts/sig**2)/2/size(res)
!!$ !   zfun = sum(res*pht/sig**2) / size(res) + &
!!$!         sum(res**2*cts/sig**5)/size(res)/4 - &
!!$!         sum(cts/sig**2) / size(res) / 2
!!$    write(*,'(a,g20.7,4g12.3)') 'zfun ',t, zfun, sum(psi*pht/sig) / n, &
!!$         sum(psi*res*cts/sig**3) / n / 2, sum(cts/sig**2) / n / 2
!!$
!!$  end function zfun

!!$  function zfun1(t)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64) :: zfun1, r, s
!!$!    real(REAL64), dimension(:), allocatable :: wpht
!!$    integer :: n!,i
!!$
!!$    n = size(res)
!!$!    allocate(wpht(n))
!!$    r = global(1)
!!$    s = global(2)
!!$!    s = 1
!!$    call residuals(t,r,res,sig)
!!$    res = cts - t*pht
!!$    sig = sqrt(cts + t**2*pht + 2e3**2 + t**2*2e3**2)
!!$!    psi = huber(res/(s*sig))
!!$    psi = tukey(res/(s*sig))
!!$   psi = res/(sig*s)
!!$
!!$!    zfun1 =   sum(psi*pht/sig) / s &
!!$!          + t*sum(psi*res*pht/sig**3) / s !&
!!$!          - t*sum(pht/sig**2)
!!$
!!$   zfun1 = sum(psi*(pht/sig + t*(res/sig)*(pht/sig)/sig))/s - t*sum(pht/sig**2)
!!$
!!$   zfun1 = sum(psi*(pht/sig + t*(res/sig)*((pht+2e3**2)/sig)/sig))/s - &
!!$         t*sum((pht+2e3**2)/sig**2)
!!$
!!$   zfun1 = sum(psi*(pht/sig + t*(res/sig)*(pht/sig)/sig))/s - t*sum(pht/sig**2)
!!$
!!$   zfun1 = zfun1 / size(cts)
!!$
!!$!    res = cts - t*pht +1.1e2  !+ 1.75e4 !+ t*1e5)
!!$!    sig = sqrt(cts + t**2*pht + 1e5**2 + t**2*1e5**2)
!!$!    psi = res/sig
!!$!    zfun1 = sum(psi*(pht/sig + t*(res/sig)*(pht/sig)/sig))/s - t*sum(pht/sig**2)
!!$
!!$!    return
!!$
!!$!    do i = 1, n
!!$!       if( abs(res(i))/(s*sig(i)) > 1.5 ) then
!!$!          wpht(i) = cts(i)/t
!!$!       end if
!!$!    end do
!!$!    zfun1 =   sum(psi*wpht/sig) / s &
!!$!         + t*sum(psi*res*wpht/sig**3) / s &
!!$!         - t*sum(wpht/sig**2)
!!$!         - t*sum(tukey(1/(s*sig))*pht/itukey(1/(s*sig))/sig**3)/2/s
!!$
!!$
!!$!    write(*,*) sum(psi*pht/sig) / s,t*sum(psi*pht*res/sig**3) / s, t*sum(pht/sig**2)
!!$
!!$!         sum((psi*res/sig + 1)*cts/sig**2)/2/size(res)
!!$ !   zfun = sum(res*pht/sig**2) / size(res) + &
!!$!         sum(res**2*cts/sig**5)/size(res)/4 - &
!!$    !         sum(cts/sig**2) / size(res) / 2
!!$    if( debug_iter ) write(error_unit,'(a,7(g0,2x))') &
!!$         'zfun1 ',real(t),real(zfun1)!,real(r),real(s), &
!!$!         sum(psi*(pht/sig + t*(res/sig)*(pht/sig)/sig))/s, &
!!$!         t*sum(pht/sig**2)
!!$    !, &
!!$!         sum(psi*pht/sig)/s,t*sum(psi*res*pht/sig**3)/s, t*sum(pht/sig**2)
!!$
!!$  end function zfun1
!!$
!!$
!!$  subroutine dtestim(t,r,s,dt,std)
!!$
!!$    use rfun
!!$    use medians
!!$
!!$    real(REAL64), intent(in) :: r,t,s
!!$    real(REAL64), intent(out) :: dt, std
!!$    real(REAL64) :: s0, h11, dp
!!$    integer :: n
!!$
!!$    n = size(res)
!!$!    s = global(2)
!!$    call residuals(t,r,res,sig)
!!$    s0 = sum(tukey(res/sig/s)**2)
!!$    dp = sum(dtukey(res/sig/s))
!!$    h11 = sum(dtukey(res/sig/s)*pht**2/sig**2) / s**2
!!$!    write(*,*) sqrt(s0),h11,&
!!$!         -2*t*sum(tukey(res/sig)*pht**2/sig**3),&
!!$!         -2*t*sum(pht**2/sig**4),sum(pht/sig**2)
!!$    if( s0 > 0 .and. h11 > 0 .and. dp > 0 ) then
!!$       dt = s*sqrt(n*s0/(dp**2*h11))
!!$       std = sqrt(s*s0/dp*n/(n-1))
!!$!       s = sqrt(s0/dp**2*n**2/(n-1))*s
!!$    else
!!$       s0 = sum(huber(res/sig/s)**2)
!!$       dp = sum(dhuber(res/sig/s))
!!$       h11 = sum(dhuber(res/sig)*pht**2/sig**2) / s**2
!!$       if( s0 > 0 .and. h11 > 0 .and. dp > 0 ) then
!!$          dt = s*sqrt(n*s0/(dp**2*h11))
!!$          std = sqrt(s*s0/dp*n/(n-1))
!!$!          s = sqrt(s0/dp**2*n**2/(n-1))*s
!!$       else
!!$          s0 = median(abs(res/sig)) / Qdev
!!$          h11 = median(pht/sig)
!!$          dt = s0 / h11 / sqrt(real(n))
!!$          std = s0
!!$!          s = s0
!!$       end if
!!$    end if
!!$
!!$  end subroutine dtestim
!!$
!!$  subroutine walley(r,tmin,tmax,tol,t,reli)
!!$
!!$    use fmm
!!$
!!$    real(REAL64), intent(in) :: r,tmin,tmax,tol
!!$    real(REAL64), intent(out) :: t
!!$    logical, intent(out) :: reli
!!$    real(REAL64) :: h
!!$
!!$    global(1) = r
!!$    t = fmin(tmin,tmax,walleyhuber,tol)
!!$
!!$    h = 2*tol
!!$    reli = tmin + h < t .and. t < tmax - h
!!$
!!$    ! a realiable result is inside range tmin .. tmax,
!!$    ! also not too close to the interval endpoints
!!$
!!$  end subroutine walley
!!$
!!$  function walleyhuber(t)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64) :: walleyhuber, r
!!$    integer :: n
!!$
!!$    n = size(res)
!!$    r = global(1)
!!$    call residuals(t,r,res,sig)
!!$!    psi = ihuber(res/sig)
!!$!    psi = abs(res)/sig
!!$!    psi = tukey(res/sig)
!!$
!!$    walleyhuber = sum(ihuber(res/sig)) / n + sum(log(sig)) / n
!!$!    walleyhuber = sum(abs(res/sig)) / n + sum(log(sig)) / n
!!$!    walleyhuber = sum(abs(res)) / n
!!$
!!$    if( debug_iter ) &
!!$         write(error_unit,'(a,g20.10,3g15.5)') 'walleyhuber',t, walleyhuber, &
!!$         sum(ihuber(res/sig))/n, sum(log(sig)) / n
!!$
!!$  end function walleyhuber



!  subroutine rnewton(t,dt,sigma,r,converge)
!!$  subroutine rnewton(t,dt,std,r,s,converge)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in out) :: t
!!$    real(REAL64), intent(out) :: dt,std!, sigma
!!$    real(REAL64), intent(in) :: r,s
!!$    logical, intent(out) :: converge
!!$
!!$    character(len=80) :: errmsg
!!$!    real(REAL64), dimension(:), allocatable :: dpsi
!!$!    real(REAL64) :: f,df,d,f2,dp,q,zd,sd,tol,s,sig2,sig3,sig4,cp,p2,t2,t3
!!$    real(REAL64) :: f,df,d,f2,dp,tol!,sig2,sig3,sig4,t2,t3zd,sd,
!!$!    integer :: n, i, iter, stat
!!$    integer :: n, iter, stat
!!$
!!$!    write(*,*) t,dt,r,s
!!$
!!$    converge = .false.
!!$    n = size(res)
!!$    tol = 10*n*epsilon(t)*(abs(t) + 1)
!!$
!!$    allocate(dpsi(n),stat=stat,errmsg=errmsg)
!!$    ! Computations needs some additional memory,
!!$    ! a fail-back way should be used if it is unavailable.
!!$    if( stat /= 0 ) then
!!$       if( debug ) write(error_unit,*) warning_label,trim(errmsg),' (rnewton)'
!!$       return
!!$    end if
!!$
!!$    if( s < epsilon(s) ) then
!!$       if( debug ) &
!!$            write(error_unit,'(2a)') warning_label," fail in rnewton(): s < 0!"
!!$       return
!!$    end if
!!$
!!$    do iter = 1, precision(t)
!!$
!!$       call residuals(t,r,res,sig)
!!$       call tukeys(res/(sig*s),psi,dpsi)
!!$!       psi = huber(res/(s*sig))
!!$!       dpsi = dhuber(res/(s*sig))
!!$       psi = res/(s*sig)
!!$       dpsi = 1
!!$
!!$!       f = sum(psi*pht/sig)/s + t*sum(psi*(res/sig)*(pht/sig)/sig)/s
!!$       f = sum(psi*((pht/sig) + t*(res/sig)*(pht/sig)/sig))/s &
!!$            - t*sum(pht/sig**2)
!!$!       write(*,*) sum(psi*pht/sig)/s,t*sum(psi*res*pht/sig**3)/s,t*sum(pht/sig**2)
!!$
!!$!       df = real(0.0,REAL64)
!!$       tt = t**2
!!$       do i = 1, 0*n
!!$          sig2 = sig(i)**2
!!$          sig3 = sig(i)**3
!!$          sig4 = sig2**2
!!$          cp = cts(i)*pht(i)
!!$          p2 = pht(i)**2
!!$          zd = -pht(i)/sig(i) + t*(t*p2 - cp)/sig3
!!$          sd = t*pht(i)/sig(i)
!!$          df = df + dpsi(i)*pht(i)*zd/sig(i) &
!!$               - psi(i)*pht(i)*sd/sig2 &
!!$               + t*dpsi(i)*cp*zd/sig3 &
!!$               -t3*psi(i)*cp*sd/sig4 &
!!$               + psi(i)*cp/sig3 &
!!$               -tt*dpsi(i)*p2*zd/sig3 &
!!$               +real(3.0,REAL64)*tt*psi(i)*p2*sd/sig4 &
!!$               -t2*psi(i)*p2/sig3 &
!!$               -t2*pht(i)*sd/sig3 &
!!$               + pht(i)/sig2
!!$       end do

!!$       df = real(0.0,REAL64)
!!$       t2 = real(2.0,REAL64)*t
!!$       t3 = real(3.0,REAL64)*t
!!$       do i = 1, 0*n
!!$          sig2 = sig(i)**2
!!$          sig3 = sig(i)**3
!!$          sig4 = sig2**2
!!$          sd = t*pht(i)/sig(i)
!!$          zd = -pht(i)/(s*sig(i)) - t*res(i)*pht(i)/sig3/s
!!$          df = df + dpsi(i)*pht(i)*zd/sig(i)/s &
!!$               -psi(i)*pht(i)*sd/sig2/s &
!!$               +psi(i)*res(i)*pht(i)/sig3/s &
!!$               +t*(dpsi(i)*res(i)*pht(i)*zd/sig3 &
!!$               -3*psi(i)*res(i)*pht(i)*sd/sig4 - psi(i)*pht(i)**2/sig3)/s &
!!$               - 0*pht(i)/sig2 &
!!$               +0*t2*pht(i)*sd/sig3
!!$       end do
!!$!       df = -sum(dpsi*(pht/sig + t*(pht/sig)*(res/sig)/sig)**2)/s**2 &
!!$!            +sum(psi*((pht/sig)*(cts-2*t*pht)/sig**2 &
!!$!            -3*t**2*(res/sig)*(pht/sig)**2/sig**2 - t*(pht/sig)**2/sig))/s**2 &
!!$!            -sum(pht/sig**2) + t**2*sum((pht/sig**2)**2)
!!$!            -3*t**2*res*pht**2/sig**5 - t*pht**2/sig**3))/s**2
!!$!            +sum(psi*(pht*(cts-2*t*pht)/sig**3 &
!!$!            - 3*t**2*(pht/sig)**2*(res/sig)/sig**2 - t*(pht/sig)**2/sig))/s
!!$!            +sum(psi & ! (pht/sig + t*(res/sig)*(pht/sig)/sig)
!!$!            *((pht/sig)*(res/sig)/sig - t*(pht/sig)**2/sig &
!!$!            - 3*t**2*(pht/sig)**2*(res/sig)/sig**2))/s
!!$
!!$!            +sum(psi*(pht/sig)*(res/sig)/sig)/s
!!$!            -t*sum(dpsi*(pht/sig)**2*(res/sig**2))/s &
!!$!            -t**2*sum(dpsi*(res/sig)**2*(pht/sig)**2/sig**2)/s &
!!$!            -3*t**2*sum(psi*(res/sig)*(pht/sig)**2/sig**2)/s
!!$
!!$       df = -sum(dpsi*(pht/sig + t*(pht/sig)*(res/sig)/sig)**2)/s**2 &
!!$            +sum(psi*((res/sig)*(pht/sig)/sig - t*(pht/sig)**2/sig  &
!!$            -t*(pht/sig)**2*(1/sig + 3*t*(res/sig)/sig**2)))/s &
!!$            -sum(pht/sig**2) + 2*t**2*sum((pht/sig)**2/sig**2)
!!$
!!$
!!$       block
!!$         real(REAL64) :: e, f1, f0, fd, g0, g1, dg
!!$         e = 1e-7
!!$         global(1) = r
!!$         global(2) = s
!!$         f1 = zfun1(t+e)
!!$         f0 = zfun1(t)
!!$         fd = (f1 - f0)/e
!!$         g0 = testfun(t)
!!$         g1 = testfun(t+e)
!!$         dg = (g1 - g0)/e
!!$         write(*,*) f,f0,f-f0,fd,df,(fd-df)/df
!!$         write(*,*) ' f(x) = ',f,dg
!!$         write(*,*) "f'(x) = ",df,fd
!!$         stop 0
!!$       end block
!!$
!!$       if( .not. (abs(df) > epsilon(df)) ) then
!!$          if( debug ) write(error_unit,'(2a)') &
!!$               warning_label," fail in rnewton(): derivative df == 0"
!!$          exit
!!$       end if
!!$
!!$       ! estimate update
!!$       d = f / df
!!$       t = t - d
!!$
!!$       if( debug_iter ) write(error_unit,'(a,i2,2x,g0,2es12.2,es14.4)') &
!!$            "mean, incr., f, f': ",iter,real(t),d,f,df
!!$
!!$       f2 = sum(psi**2)  ! mean of psi**2
!!$       dp = sum(dpsi)    ! mean of dpsi
!!$       if( f2 > 0 .and. df < 0 .and. dp > 0 ) then
!!$          dt = s*sqrt(-f2/df/(n-1))*(n/dp) ! df is an element of Jacobian
!!$          std = sqrt(s*f2/dp*n/(n-1))
!!$!          sigma = s*sqrt(f2/dp**2*n**2/(n-1))
!!$!          sigma = dt*sqrt(real(n))*sqrt(-df)
!!$          converge = abs(d) < max(0.001*dt, tol)
!!$!          converge = abs(d) < tol
!!$       else ! f2 == 0 or df > 0, t > |c|
!!$!          dt = 0
!!$!          sigma = 1
!!$          if( debug ) then
!!$             if( df > 0 ) write(error_unit,'(2a)') warning_label,&
!!$                  " fail in rnewton(): convergence to local maximum."
!!$             if( dp < epsilon(dp) ) write(error_unit,'(2a)') warning_label,&
!!$                  " fail in rnewton(): out of range."
!!$          end if
!!$          if( abs(d) < tol ) exit   ! when a local maximum is reached
!!$          if( df > 0 ) exit         ! when a local maximum is encountered
!!$       end if
!!$
!!$       ! exit of iterations: the absolute error must be at least |d| < tol
!!$       if( converge ) exit
!!$
!!$    end do
!!$
!!$    deallocate(dpsi)
!!$    converge = .false.
!!$
!!$  end subroutine rnewton
!!$
!!$  function testfun(t) result(f)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64) :: f,r,s
!!$
!!$    r = global(1)
!!$    s = global(2)
!!$    call residuals(t,r,res,sig)
!!$!    psi = huber(res/sig/s)
!!$!    psi = tukey(res/sig/s)
!!$!    psi = res/sig/s
!!$
!!$!    f = -sum(itukey(res/(s*sig))) !- sum(log(sig*s))
!!$    f = -sum(ihuber(res/(s*sig))) -sum(log(sig*s))
!!$!    f = -sum((res/sig)**2/s**2)/2 !- sum(log(sig*s))
!!$    write(*,*) 'testfun: ',f!sum((res/sig/s)**2)/2!,sum(log(sig*s))
!!$
!!$  end function testfun
!!$
!!$
!!$  subroutine rgrid(t,r,reli)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), intent(out) :: r
!!$    logical, intent(out) :: reli
!!$
!!$    integer, parameter :: ngrid = 23
!!$    real(REAL64), dimension(ngrid+2) :: ord, tab, freq, dis
!!$    real(REAL64) :: h,p,f2,df,s,pmax,rmin, rmax, dmin
!!$    integer :: i,n,cmin,nmin
!!$
!!$    rmin = 0
!!$    call rquant(t,rmax,reli)
!!$    if( .not. reli ) then
!!$       rmax = maxval(abs(cts - t*pht))
!!$    end if
!!$
!!$    write(*,*) 'rmin, rmax: ',rmin, rmax
!!$
!!$    n = size(res)
!!$    h = (rmax - rmin) / (ngrid+2)
!!$
!!$    do i = 1, ngrid+2
!!$       r = rmin + (i-1)*h
!!$       ord(i) = r
!!$       call residuals(t,r,res,sig)
!!$       freq(i) = sum(dtukey(res/sig)) / n
!!$       f2 = sum(tukey(res/sig)**2)
!!$       df = sum(dtukey(res/sig))
!!$       if( f2 > 0 .and. df > 0 ) then
!!$          s = sqrt(f2/df*n/(n-1))
!!$          dis(i) = abs(s - stukey)
!!$       else
!!$          dis(i) = 1
!!$       end if
!!$       tab(i) = -sum(itukey(res/sig)) - sum(log(sig))
!!$ !      write(error_unit,'(4(f15.5))') ord(i),dis(i),freq(i),tab(i)
!!$    end do
!!$
!!$    dmin = 1
!!$    pmax = 1
!!$    cmin = 0
!!$    nmin = 0
!!$    do i = 2, ngrid+1
!!$       if( tab(i-1) < tab(i) .and. tab(i) > tab(i+1) ) then
!!$          cmin = cmin + 1
!!$!          if( abs(freq(i) - 0.683) < abs(pmax -0.683) ) then
!!$!          if( freq(i) > 0.3 .and. nmin == 0 ) then
!!$!             nmin = i
!!$!             pmax = freq(i)
!!$!          end if
!!$          if( dis(i) < dmin ) then
!!$             nmin = i
!!$             dmin = dis(i)
!!$          end if
!!$
!!$          if( debug ) write(error_unit,'(a,i3,2x,f0.1,2x,f7.3,f7.3)') &
!!$               'minloc:',i,ord(i),dis(i),freq(i)
!!$       end if
!!$
!!$    end do
!!$    if( nmin > 0 ) then
!!$       r = ord(nmin)
!!$    end if
!!$!       r = ord(4)
!!$
!!$
!!$    ! found minimum with required data amount
!!$    reli = 2 <= nmin .and. nmin <= ngrid+1
!!$
!!$    if( debug ) write(error_unit,'(3(a,f0.1,1x),2(a,i0,1x),l1)') &
!!$         'Grid search: ',rmin," < ",r," < ",rmax,"  n:",ngrid,' #',cmin,reli
!!$
!!$  end subroutine rgrid
!!$

!!$  subroutine rscale(t,r,converge)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), intent(in out) :: r
!!$    logical, intent(out) :: converge
!!$
!!$    character(len=80) :: errmsg
!!$    real(REAL64), dimension(:), allocatable :: dpsi
!!$    real(REAL64) :: f,df,d,f2,dp,w,tol,dr,fs
!!$    integer :: n, iter, stat
!!$
!!$    converge = .false.
!!$    n = size(res)
!!$    tol = 10*n*epsilon(t)*(abs(t) + 1)
!!$
!!$    allocate(dpsi(n),stat=stat,errmsg=errmsg)
!!$    ! Computations needs some additional memory,
!!$    ! a fail-back way should be used if it is unavailable.
!!$    if( stat /= 0 ) then
!!$       if( debug ) write(error_unit,*) warning_label,trim(errmsg)
!!$       return
!!$    end if
!!$
!!$    do iter = 1, precision(t)
!!$
!!$       call residuals(t,r,res,sig)
!!$       call tukeys(res/sig,psi,dpsi)
!!$!       psi = huber(res/sig)
!!$       !       dpsi = dhuber(res/sig)
!!$       psi = res/sig
!!$       dpsi = 1
!!$
!!$       fs = sum(psi*res/sig**3) - sum(1/sig**2)
!!$       f = r*fs
!!$       df = r**2*sum(dpsi*res**2/sig**6) &
!!$          + sum(psi*res*(sig**2 - 3*r**2)/sig**5) &
!!$          - sum((sig**2 - 2*r**2)/sig**4)
!!$       df = sum(dpsi*res**2/sig**6) - 3*sum(psi*res/sig**5) + 2*sum(1/sig**4)
!!$       df = r**2*df + fs
!!$!       df = df + fs
!!$
!!$
!!$       block
!!$         real(REAL64) :: e, f1, f0, f2, fd, fd2
!!$         e = 1e-3
!!$         global = t
!!$         f1 = -rfitfun(r+e)*n
!!$         f0 = -rfitfun(r)*n
!!$         fd = (f1 - f0)/e
!!$         call residuals(t,r+e,res,sig)
!!$         psi = res/sig
!!$         dpsi = 1
!!$         g1 = (r+e)*(sum(psi*res/sig**3) - sum(1/sig**2))
!!$         call residuals(t,r,res,sig)
!!$         psi = res/sig
!!$         dpsi = 1
!!$         g0 = r*(sum(psi*res/sig**3) - sum(1/sig**2))
!!$!         g0 = -rfitfun(r)*n
!!$!         g1 = -rfitfun(r+e)*n
!!$         f2 = -rfitfun(r-e)*n
!!$         fd2 = (f1 - 2*f0 + f2)/e**2
!!$         write(*,*) f,fd,(f-fd)/f,df,fd2,(df-fd2)/df, &
!!$              sum(dpsi*res**2/sig**6)
!!$!         write(*,*) f,f0,f-f0,fd,df,(fd-df)/df
!!$       end block
!!$       stop 0
!!$
!!$       if( .not. (abs(df) > epsilon(df)) ) exit
!!$
!!$       ! update location
!!$       d = f / df
!!$       r = r - d
!!$
!!$!       if( debug ) write(*,'(a,i2,g0,1p3g12.3)') &
!!$       if( debug ) write(*,'(a,i2,2x,g0,2es12.2,es14.4)') &
!!$            "# r, incr., f, f': ",iter,real(r),d,f,df
!!$
!!$       f2 = sum(psi**2)  ! mean of psi**2
!!$       dp = sum(dpsi)    ! mean of dpsi
!!$       w = -df*dp         ! df is an element of Jacobian
!!$       if( f2 > 0 .and. w > 0 ) then
!!$!          write(*,*) '>',f2,q,dp,n
!!$          dr = sqrt(f2/w)
!!$          converge = abs(d) < 0.001*dr
!!$       else ! f2 == 0 or q == 0, t > |c|
!!$!          dt = 0
!!$!          s = 1
!!$          if( abs(d) < tol ) exit
!!$       end if
!!$
!!$       ! exit of iterations: the absolute errot must be at least |d| < tol
!!$       if( converge ) exit
!!$
!!$    end do
!!$
!!$  end subroutine rscale

!!$  subroutine phmean3d(t1,t2,rmax,t0)
!!$
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: t1,t2,rmax,t0
!!$    real(REAL64) :: t,dt,r,s,stdsig
!!$    integer :: i,j,n
!!$
!!$    !    if( .not. (rmax > 0) ) stop 'Assert (rmax > 0) failed.'
!!$    if( .not. ((t2 - t1) > 0) ) stop 'Assert (t2 - t1 > 0) failed.'
!!$
!!$    n = size(pht)
!!$
!!$    dt = (t2 - t1) / 1000
!!$    s = stukey
!!$    open(1,file='/tmp/ph3d')
!!$    do j = 100, 1000,10
!!$       r = 0.001*j*rmax
!!$       do i = 10,1000,10
!!$          t = t1 + dt*i
!!$          call residuals(t,r,res,sig)
!!$          write(1,*) t,r,-sum(ihuber(res/sig/s))/n - sum(log(s*sig)) / n,&
!!$               -sum(itukey(res/sig/s))/n - sum(log(s*sig)) / n, &
!!$               sum(ihuber(res/sig/s))/n + sum(log(s*sig)) / n, &
!!$               sum(itukey(res/sig/s))/n + sum(log(s*sig)) / n
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)
!!$
!!$!    t = 10.3958244
!!$    !    t = 2
!!$    !t = 1.34372663
!!$    t = t0
!!$    open(1,file='/tmp/rs3d')
!!$    do j = -1000, 5000, 50
!!$       r = j
!!$       do i = 500,1500,50
!!$          s = 0.001*i
!!$          call residuals(t,r,res,sig)
!!$          stdsig = sqrt(s*sum(tukey(res/sig/s)**2)/ &
!!$               max(sum(dtukey(res/sig/s))*n/(n-1),epsilon(s)))
!!$          stdsig = min(stdsig,1d1)
!!$          write(1,*) r,s, &
!!$!               -sum((res/sig/s)**2)/2 - n*log(s), &
!!$!               -sum(ihuber(res/(s*sig))) - n*log(s), &
!!$!               -sum(itukey(res/(s*sig))) - n*log(s)
!!$               -sum((res/sig/s)**2)/2 - sum(log(sig*s)), &
!!$               -sum(ihuber(res/sig/s)) - sum(log(sig*s)), &
!!$               -sum(itukey(res/sig/s)) - sum(log(sig*s)), &
!!$               -sum(itukey(res/sig/s)) - sum(log(sig*s)) - &
!!$               (s-stukey)-50*(s-stukey)**2/2
!!$!               -(sum((res/sig/s)**2)/2 + sum(log(sig*s)))/n -0e7*(s-1)**2/2 &
!!$!               +0e6*log(max(1e4-r,1d-7)), &
!!$!               -(sum(ihuber(res/sig/s)) + sum(log(sig*s)))/n-0e5*(s-0.87)**2/2 &
!!$!               +0e2*log(max(1e4-r,1d-7)), &
!!$!               -(sum(itukey(res/sig/s)) + sum(log(sig*s)))/n-0e5*(s-0.9)**2/2 &
!!$!               +0e2*log(max(1e4-r,1d-7))
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)
!!$
!!$    t = t0
!!$    s = stukey
!!$    open(1,file='/tmp/rd')
!!$    do j = 1, 10000, 10
!!$       r = 0.001*j*rmax
!!$       call residuals(t,r,res,sig)
!!$       write(1,*) r,-sum((res/sig/s)**2)/2 - sum(log(s*sig)), &
!!$            -sum(ihuber(res/sig/s)) - sum(log(sig*s)), &
!!$            -sum(itukey(res/sig/s)) - sum(log(sig*s))!, &
!!$!            sum((res/sig)**2) - n, &
!!$!            sum(huber(res/sig)*(res/sig)) - n*0.6582, &
!!$!            sum(tukey(res/sig)*(res/sig)) - n*0.8061
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine phmean3d

!!$  subroutine signtest(t,r,phi)
!!$
!!$    use quantilmean
!!$    use rfun
!!$
!!$    real(REAL64), intent(in) :: t,r,phi
!!$!    real(REAL64), intent(in) :: mu,t,sigu,sigv
!!$    real(REAL64), dimension(:), allocatable :: u,v,f
!!$    real(REAL64), dimension(2) :: q, tq
!!$!    real(REAL64) :: qmin, qmax, tmin, tmax
!!$    real(REAL64) :: tmin, tmax, dmin, d, mu, s
!!$    integer :: i,n,nplus,nsig,nall,flag
!!$
!!$!    n = size(res)
!!$    n = 1000
!!$    allocate(u(n),v(n),f(n))
!!$
!!$!    call ecdf(tt,x,p)
!!$!    call quantilfun(tt,q,x,flag)
!!$
!!$    q = [ 0.1587, 0.8413 ]
!!$!    q = [ 0.077, 0.77 ]
!!$!    call quantilefun(tt,q,tq,flag)
!!$!    if( flag /= 0 ) return
!!$
!!$!    tmin = tq(1)
!!$!    tmax = tq(2)
!!$
!!$    mu = mju(1)
!!$    s = sqrt(mu + r**2)
!!$!    u(1) = max(tt(1),t - 3*r/mu)
!!$!    u(n) = min(tt(n),t + 5*r/mu)
!!$    u(1) = t - 7*s/mu
!!$    u(n) = t + 7*s/mu
!!$!    write(*,*) tt(1),tt(n),t - 3*s/mu,t + 5*s/mu
!!$    d = (u(n) - u(1))/(n - 1)
!!$    u = [ (u(1) + d*(i-1), i = 1,n) ]
!!$!    call cdf(t,r,phi,u,v)
!!$!    s = 0.245
!!$!    v = NormCDF((u-t)/s)
!!$
!!$!    s  = 0
!!$    f = pdf(u,mu,t,r,phi)
!!$!    f = f / sum(f)
!!$    v(1) = 0
!!$    do i = 2, n
!!$       v(i) = v(i-1) + (f(i) + f(i-1))/2*d
!!$    end do
!!$
!!$    tmin = u(1)
!!$    tmax = u(n)
!!$!    d = 1.5 / n
!!$!    write(*,*) d
!!$!    open(1,file='/tmp/d')
!!$    dmin = 1
!!$    do i = 1, n
!!$       d = abs(v(i) - q(1))
!!$       if( d < dmin ) then
!!$          dmin = d
!!$          tmin = u(i)
!!$          !          write(*,*) '--min',v(i), u(i),dmin
!!$!       else
!!$!          exit
!!$       end if
!!$       if( d > dmin ) exit
!!$    end do
!!$    dmin = 1
!!$    do i = n, 1, -1
!!$       d = abs(v(i) - q(2))
!!$       if( d < dmin ) then
!!$          dmin = d
!!$          tmax = u(i)
!!$          !          write(*,*) '++max',v(i), u(i),dmin
!!$!       else
!!$!          exit
!!$       end if
!!$       if( d > dmin ) exit
!!$!       write(1,*) v(i),u(i),abs(v(i) - q(1)),abs(v(i) - q(2))
!!$    end do
!!$!    close(1)
!!$
!!$    write(*,*) tmin, tmax,s,r
!!$
!!$    open(1,file='/tmp/e666')
!!$    do i = 1, n
!!$       write(1,*) u(i),v(i),f(i)
!!$    end do
!!$    close(1)
!!$
!!$    qmin = 0.1587
!!$    qmax = 0.8413
!!$
!!$    tmin = x(n)
!!$    tmax = x(1)
!!$    do i = 1, n
!!$       if( p(i) < qmin ) tmin = x(i)
!!$       if( p(i) < qmax ) tmax = x(i)
!!$!       if( qmin <= p(i) .and. p(i) <= qmax ) then
!!$!          nsum = nsum + 1
!!$!          if( p(i) >= 0.5 ) nplus = nplus + 1
!!$!       end if
!!$    end do
!    write(*,*) tmin, tmax
!!$
!!$    nplus = 0
!!$    nsig = 0
!!$    nall = size(tt)
!!$    do i = 1, size(tt)
!!$       if( tmin < tt(i) .and. tt(i) < tmax ) then
!!$          nsig = nsig + 1
!!$          if( tt(i) > t ) nplus = nplus + 1
!!$       end if
!!$    end do
!!$!    write(*,*) 'signtest:',nsig,nsig/2,nplus,sqrt(nplus*0.25 + 1e-4)
!!$
!!$    write(error_unit,'(a,g0,a,g0)') &
!!$         'Equivalent 1-sigma interval: ',real(tmin),' .. ',real(tmax)
!!$    write(error_unit,'(a,i0,a,i0,a,f0.1,a)') &
!!$         'Points within the interval: ',&
!!$         nsig,'/',nall,' (',100.*real(nsig)/real(nall),'%)'
!!$    write(error_unit,'(a,i0,"/",f0.1,", ",i0,"+-",f0.1)') &
!!$         'Sign test in the interval (total/expected, +): ',&
!!$         nsig,nsig/2.0,nplus,sqrt(nplus*0.25 + 1e-4)
!!$
!!$  end subroutine signtest


    ! sign test   --- vyhodit? asymetricke?
!!$    block
!!$      integer :: nplus,nsig,nall
!!$
!!$      nplus = 0
!!$      nsig = 0
!!$      nall = size(tt)
!!$      do i = 1, size(tt)
!!$         if( tlo < tt(i) .and. tt(i) < thi ) then
!!$            nsig = nsig + 1
!!$            if( tt(i) > tau ) nplus = nplus + 1
!!$         end if
!!$      end do
!!$      !    write(*,*) 'signtest:',nsig,nsig/2,nplus,sqrt(nplus*0.25 + 1e-4)
!!$
!!$      write(error_unit,'(a,g0,a,g0)') &
!!$           'Equivalent 1-sigma interval: ',real(tlo),' .. ',real(thi)
!!$      write(error_unit,'(a,i0,a,i0,a,f0.1,a)') &
!!$           'Points within the interval: ',&
!!$           nsig,'/',nall,' (',100.*real(nsig)/real(nall),'%)'
!!$      write(error_unit,'(a,i0,"/",f0.1,", ",i0,"+-",f0.1)') &
!!$           'Sign test in the interval (total/expected, +): ',&
!!$           nsig,nsig/2.0,nplus,sqrt(nplus*0.25 + 1e-4)
!!$    end block

  end subroutine kstest


!!$  subroutine cdf1(mu,t,r,phi,u,v)
!!$!  subroutine cdf(t,sigu,sigv,umin,umax,mu,u,v)
!!$
!!$   real(REAL64), intent(in) :: mu,t,r,phi!,umin, umax
!!$ !   real(REAL64), intent(in) :: t,sigu,sigv,umin, umax
!!$!    real(REAL64), dimension(:), intent(in) :: mu
!!$    real(REAL64), dimension(:), intent(in) :: u
!!$    real(REAL64), dimension(:), intent(out) :: v
!!$    real(REAL64), dimension(:), allocatable :: sig,f
!!$    real(REAL64) :: d, s
!!$    integer :: i, m
!!$
!!$    sig = sqrt(mu*t*(1 + u**2)+r**2*(1 + ((u/t)**2 - 1)*cos(phi)**2))
!!$    v = NormCDF((mu/sig)*(u - t))
!!$    return
!!$
!!$    m = size(u)
!!$!    u(1) = umin
!!$!    u(m) = umax
!!$!    d = (umax - umin) / (m + 1)
!!$!    d = 3.0 / m
!!$!    do i = 1, m
!!$!       u(i) = umin + d*i
!!$!       u(i) = 1 + d*i
!!$!    end do
!!$!    u = [(umin+d*i,i=1,m)]
!!$
!!$    s = 0
!!$    do i = 1, m
!!$!       v(i) = pdf(tt(i),mju(i),t,r,phi)
!!$       f(i) = xpdf(u(i),mju,t,r,phi)
!!$!       s = s + v(i)
!!$    end do
!!$    !    write(*,*) 'norm:',sum(f)
!!$    s = sum(f)
!!$    if( s > 0 ) then
!!$       f = f / s
!!$    else
!!$       f = 0 !?
!!$    end if
!!$!    q = v
!!$
!!$
!!$    v(1) = 0
!!$    do i = 2, m
!!$!       v(i) = v(i) + v(i-1)
!!$       v(i) = (sum(f(2:i-1)) + (f(1) + f(i))/2) ! *d
!!$    end do
!!$
!!$!    v = f
!!$
!!$  end subroutine cdf1
!!$
!!$  function xpdf(tt,mju,t,r,phi) result(s)
!!$
!!$    use quicksorter
!!$
!!$    real(REAL64), intent(in) :: tt,t,r,phi
!!$    real(REAL64), dimension(:),intent(in) :: mju
!!$    real(REAL64), dimension(:), allocatable :: xmu
!!$    integer, dimension(:), allocatable :: w
!!$    real(REAL64) :: s, f1, f2, mu0
!!$    integer :: i,j,n,m,l
!!$    logical :: found
!!$
!!$!    write(*,*) 'mju:',real(mju)
!!$
!!$    n = size(cts)
!!$    allocate(xmu(n),w(n))
!!$    !    mu = (pht + cts / t) / 2
!!$!    mu = mju
!!$!    call qsort(mu)
!!$    m = 1
!!$    xmu(1) = mju(1)
!!$    do i = 2, n
!!$       found = .false.
!!$       do j = 1, m
!!$          if( abs(mju(i) - xmu(j)) < 1 ) then
!!$             found = .true.
!!$             goto 10
!!$          end if
!!$       end do
!!$10     continue
!!$       if( .not. found ) then
!!$          m = m + 1
!!$          xmu(m) = mju(i)
!!$       end if
!!$    end do
!!$    call qsort(xmu(1:m))
!!$!    write(*,*) m, real(xmu(1:m))
!!$
!!$    w = 0
!!$    do i = 1, n
!!$       l = 0
!!$       do j = 1, m
!!$          if( abs(mju(i) - xmu(j)) < 1 ) then
!!$             w(j) = w(j) + 1
!!$             exit
!!$          end if
!!$       end do
!!$    end do
!!$!    write(*,*) 'w:',real(w(1:m))
!!$
!!$    block
!!$      integer :: j
!!$      logical :: ex
!!$      real(REAL64) :: u
!!$
!!$      inquire(file='/tmp/f3d',exist=ex)
!!$      if( .not. ex ) then
!!$         open(1,file='/tmp/f3d')
!!$         do i = 1800, 2200, 10
!!$            u = i / 1000.0
!!$            do j = 1, m
!!$               write(1,*) u,xmu(j),w(j)*pdf(u,xmu(j),t,r,phi)
!!$            end do
!!$            write(1,*)
!!$         end do
!!$         close(1)
!!$      end if
!!$    end block
!!$
!!$
!!$!    write(*,*) w(1:m)
!!$
!!$!    m = 100
!!$!    xmu = [ (0.7e3*(i-1),i=1,m) ]
!!$
!!$!    m = 1
!!$!    xmu(1) = 1e4
!!$
!!$!    s = 0
!!$!    f1 = 0
!!$!    mu0 = mju(1)
!!$    f1 = w(1)*pdf(tt,xmu(1),t,r,phi)
!!$    s = f1/2*xmu(1)
!!$!    write(*,*) m,f1,w(1),xmu(1)
!!$!    mu0 = (pht(1) + cts(1) / t) / 2
!!$    do i = 2, m
!!$       !       mu = (pht(i) + cts(i) / t) / 2
!!$       !f2 = pdf0(mu(i),cts(i)/pht(i),t,sigu,sigv)
!!$!       f2 = pdf(cts(i)/pht(i),mju(i),t,r)
!!$       f2 = w(i)*pdf(tt,xmu(i),t,r,phi)
!!$       !       s = s + (f1 + f2)/2!*(mu(i) - mu(i-1))
!!$!       mu = (cts(i)/t + pht(i))/2
!!$       s = s + (f1 + f2)/2*(xmu(i) - xmu(i-1))
!!$!       s = s + f2
!!$       f1 = f2
!!$!       mu0 = mju(i)
!!$    end do
!!$
!!$!    mu0 = 1e4
!!$!    s = pdf(mu0,mu0,t,sigu,sigv)*n
!!$
!!$  end function xpdf

!!$  elemental pure function pdf(t,mu,tau,r,phi)
!!$
!!$    real(REAL64) :: pdf
!!$    real(REAL64), intent(in) :: t,mu,tau,r,phi
!!$
!!$    real(REAL64), parameter :: expmax = log(huge(pdf))
!!$!    real(REAL64), parameter :: sqrtpi = sqrt(PI)
!!$    real(REAL64) :: sig, c2, beta, res, logpdf
!!$
!!$    c2 = cos(phi)**2
!!$    sig = sqrt(mu*(tau + t**2) + r**2*(1+((t/tau)**2 - 1)*c2))
!!$    beta = mu*tau*(t + 1) + r**2*(1 + (t/tau - 1)*c2)
!!$    res = (mu/sig)*(t - t)
!!$
!!$    logpdf = -res**2/2 + log(beta/sig**2) + log(mu/sig) + log(2*PI)/2
!!$    if( abs(logpdf) < expmax ) then
!!$       pdf = exp(logpdf)
!!$    else
!!$       pdf = 0
!!$    end if
!!$
!!$  end function pdf
!!$
!!$  function pdf0(mu,tt,t,sigu,sigv)
!!$
!!$    real(REAL64) :: pdf0
!!$    real(REAL64), intent(in) :: t,tt,mu,sigu,sigv
!!$    real(REAL64), parameter :: sqrtpi = real(2.5066282746310005,REAL64)
!!$    real(REAL64) :: sig, sigu2, sigv2 !sigx, sig2, sigx2
!!$
!!$!    sigx = sqrt(mu + sigv**2)
!!$!    sigx2 = mu + sigv**2
!!$!    sig = sqrt(mu*(t**2 + r) + t**2*sigv**2 + sigu**2)
!!$!    sig2 = mu*(t**2 + r) + t**2*sigv**2 + sigu**2
!!$!    sig = sqrt(sig2)
!!$!    pdf0 = (mu/sig)*((r - t)*t*sigx2/sig2 + 1)*exp(-mu**2/2*(t-r)**2/sig2)
!!$!    pdf = log(mu*((r - t)*t*sigx**2/sig**3 + 1/sig)) - mu**2/2*(t-r)**2/sig**2
!!$!    pdf = exp(pdf)
!!$!    pdf0 = pdf0 / sqrtpi
!!$
!!$    sigu2 = sigu**2
!!$    sigv2 = sigv**2
!!$    sig = sqrt(mu*(tt + t**2) + sigu2 + t**2*sigv2)
!!$    pdf0 = mu*(mu*tt*(1+t)+tt*sigv2 + sigu2)/sig**3*exp(-(mu*(tt - t)/sig)**2/2)
!!$    pdf0 = pdf0 / sqrtpi
!!$
!!$  end function pdf0
!!$
!!$  function pdf2(y,x,t,sigu,sigv)
!!$
!!$    real(REAL64) :: pdf2
!!$    real(REAL64), intent(in) :: y,x,t,sigu,sigv
!!$    real(REAL64), parameter :: sqrtpi = real(2.5066282746310005,REAL64)
!!$    real(REAL64) :: sig, sigu2, sigv2
!!$
!!$    sigu2 = sigu**2
!!$    sigv2 = sigv**2
!!$    sig = sqrt(y + sigu2 + t**2*(x + sigv2))
!!$!    pdf = (y*(x*(1+t) + t*sigv2) + x*sigu2)/sig**3*exp(-((y - t*x)/sig)**2/2)
!!$    pdf2 = (y*t*(x+sigv2) + x*(y+sigu2))/sig**3*exp(-((y - t*x)/sig)**2/2)
!!$    pdf2 = pdf2 / sqrtpi
!!$
!!$  end function pdf2


end module photonrate
