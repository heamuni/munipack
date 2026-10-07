
program testsr1mean_hist

  ! test the estimator

  use toolbox
!  use rfun
  use rfun_REAL64
  use rfun_hampel_REAL64
!  use sr1mean_function
  use rlike_REAL64
!  use rlike
  use iso_fortran_env
  implicit none

  logical, parameter :: winsor = .false.
  character(len=15), parameter :: robfun = 'hampel'
  integer, parameter :: ndata = 4200
  integer, parameter :: nsamples = 10000
  integer, parameter :: k = ndata
  real(REAL64), parameter :: pi = 3.141592653589793_REAL64
  real(REAL64), parameter :: tpar =  0
  real(REAL64), parameter :: sigpar = 1
  real(REAL64), parameter :: beta = 0.8718
  real(REAL64), parameter :: s0 = sqrt(beta)
  real(REAL64), parameter :: gama = 0.5956
  real(REAL64), parameter :: nu = beta
  real(REAL64), parameter :: hnorm = 2.57657
  integer, parameter :: step = 4
  real(REAL64), parameter :: reltol = 1.0 / 42.0 * 1e-3
  real(REAL64), dimension(:), allocatable :: t,s,sig, med, Qn, data, chi2, sig2, &
       res, avgR
  real(REAL64), dimension(:), allocatable, target :: dataset
  real(REAL64), dimension(:,:,:), allocatable :: xhess
  real(REAL64), dimension(2,2) :: hess
  real(REAL64) :: t0, t1, dt, d, mu, logmu, at,asig, fdat, hbin, sqrt2k, chi2med, &
       chi2var, avgs, cscale, avgsig, avgnu, avgvarpi, avgt
  integer :: i, l, nfail, info, niter, nfev, nrst, nout, nidf, nhrd, nbin

  allocate(t(nsamples),s(nsamples),sig(nsamples),dataset(ndata*nsamples),&
       med(nsamples),Qn(nsamples),data(ndata),chi2(nsamples),sig2(nsamples), &
       res(ndata),avgR(nsamples),xhess(nsamples,2,2))

  nfail = 0
  niter = 0
  nfev = 0
  nrst = 0
  nout = 0
  nidf = 0
  nhrd = 0
  hess = 0
  do i = 1, nsamples
     call rmean_new(t(i),s(i),sig(i),med(i),Qn(i),xhess(i,:,:),data,info)
     l = i*ndata
     dataset(1+(l-ndata):l) = data
     at = sum(data)/ndata
     asig = sqrt(sum((data-at)**2)/(ndata-1))
     res = (data-t(i)) / s(i)
     chi2(i) = sig(i)**2*(ndata-1)/sigpar**2
!     chi2(i) = 1/xhess(i,1,1)
     avgR(i) = sum(ipshampel_real64(res)) / size(res)
!     chi2(i) = sig(i)**2*(ndata-1.0)/ndata
!     chi2(i) = 1 / chi2(i)
!     chi2(i) = s(i)**2 * sum(pshampel(res)**2) / sum(dpshampel(res))**2 *(ndata-1)**2
!     chi2(i) = sum(pshampel(res)**2) / sum(dpshampel(res))**2 * ndata**2
!     chi2(i) = sum(pshampel(res)**2)
!     chi2(i) = s0**2 *sum(pshampel(res)**2) / sum(dpshampel(res))**2 *ndata**2
!     chi2(i) = s(i)**2 * sum(ipshampel(res)) / sum(dpshampel(res))**2 *ndata**2
!     chi2(i) = 2*sum(ipshampel(res))*s(i)**2/sig(i)
!     chi2(i) = 2*sum(ipshampel((data-t(i))/sig(i)))
!     chi2(i) = 2*sum(ipshampel(data))
!     chi2(i) = sum(pshampel(data)**2)
!     chi2(i) = (s(i)/sig(i))**2 * sum(res**2)
!     chi2(i) = sum(data**2)
!     chi2(i) = sum((data-at)**2/asig)
!     chi2(i) = sum(data**2)
     if( .not. (info == 0 .or. info == 1) ) nfail = nfail + 1
!     sig(i) = sqrt(chi2(i))
  end do
!  sig2 = sig**2*ndata

  hess = hess / nsamples
  write(*,'("mean H: ",2(2x,g0))') hess
  write(*,*) 'mean t:',sum(t)/nsamples
  write(*,*) 'mean s:',sum(s)/nsamples
  avgsig = sum(sig)/nsamples
  write(*,*) 'mean stdsig:',avgsig

  write(*,*) 'Convergence failed in ',nfail,' case(s).'
  write(*,*) 'Total niter=',niter,' nfev=',nfev,', nres=',nrst, ' nout=',nout, &
       ' nidf=',nidf,' nhrd=',nhrd
  fdat = nsamples
  write(*,'(6(a,f0.2))') &
       'Average  niter=',niter/fdat,' nfev=',nfev/fdat,', nres=',nrst/fdat,&
       ' nout=',nout/fdat,' nidf=',nidf/fdat,' nhrd=',nhrd/fdat
  avgt = sum(t)/nsamples
  avgs = sum(s)/nsamples
  avgnu = hess(1,1)*avgs**2
  avgvarpi = hess(2,2)*avgs**2-1
  write(*,'(5(a,f0.4))') 'ν = ',avgnu,'  √ν = ',sqrt(max(avgnu,1e-3)), &
       "  ϖ = ",avgvarpi,' avgt = ',avgt,' avgs = ',avgs
  cscale = avgs/0.9485
!  stop

  write(error_unit,'(a)',advance="no") "Calculating graphs .. "

  d = sigpar/sqrt(1.0*ndata)
  t0 = tpar - 5*d
  t1 = tpar + 5*d
  dt = d/20
  hbin = 3.5*d/nsamples**(1.0/3.0)
  nbin = (t1 - t0) / hbin
  call histogram(t,t0,t1,nbin,'/tmp/testsr1mean_t',.true.)
  call density_loc(t0,t1,dt,ndata,'/tmp/testsr1mean_density_t')
  d = sigpar/sqrt(2.0*ndata)
  t0 = max(sigpar - 5*d,0.0_REAL64)  ! 0.45
  t1 = sigpar + 5*d  ! 1.45
  d = avgs/sqrt(2.0*ndata)
  t0 = max(avgs - 5*d,0.0_REAL64)
  t1 = avgs + 5*d
  dt = d/20
  call histogram(s,t0,t1,nbin,'/tmp/testsr1mean_s',.true.)
  call density_s(t0,t1,dt,ndata,'/tmp/testsr1mean_density_s')
  logmu = log(2.0)/2 + log_gamma((k+1)/2.0) - log_gamma(k/2.0)
  mu = exp(logmu)
  d = sigpar/sqrt(1.0*ndata)
  t0 = max(sigpar - 5*d,0.0_REAL64)  !3
  t1 = sigpar + 5*d  !10
!  dt = d/20
  call histogram(sig,t0,t1,nbin,'/tmp/testsr1mean_sig',.true.)
  d = sigpar/sqrt(ndata-2.0)
  t0 = max(sqrt(2.0)*sigpar - 5*d,0.0_REAL64)  ! 0.45
!  t0 = sigpar - 5*d  ! 0.45
  t1 = sigpar + 5*d  ! 1.45
!  t0 = 3
!  t1 = 11
!  call histogram(sig*sqrt(ndata-1.0),t0,t1,150,'/tmp/testsr1mean_sig',.true.)
!  call density_sig(t0,t1,dt,ndata,'/tmp/testsr1mean_density_sig')

  t0 = 0
  t1 = 5
  call histogram(avgR,t0,t1,nbin,'/tmp/testsr1mean_avgR',.true.)

  d = sqrt(2.0*ndata)
  t0 = max(k - 5*d,0.0_REAL64)
  t1 = k + 5*d
  dt = d/20
!  call histogram(chi2/sigpar**2,t0,t1,nbin,'/tmp/testsr1mean_chi2',.true.)
  call histogram(chi2,t0,t1,nbin,'/tmp/testsr1mean_chi2',.true.)
!  call histogram(ndata/xhess(:,1,1),t0,t1,nbin,'/tmp/testsr1mean_chi2',.true.)
  call density_chi2(t0,t1,dt,ndata,'/tmp/testsr1mean_density_chi2')

  d = sqrt(2.0*ndata)
  t0 = -0.5
  t1 = 3
  dt = d/20
  call histogram(xhess(:,1,1),t0,t1,nbin,'/tmp/testsr1mean_hess11',.true.)

  d = sqrt(2.0*ndata)
  t0 = -1
  t1 = 1
  dt = d/20
  call histogram(xhess(:,1,2),t0,t1,nbin,'/tmp/testsr1mean_hess12',.true.)

  d = sqrt(2.0*ndata)
  t0 = -1
  t1 = 5
  dt = d/20
  call histogram(xhess(:,2,2)-1,t0,t1,nbin,'/tmp/testsr1mean_hess22',.true.)

  block
    real(REAL64), parameter :: varpi0 = 0.5643
    real(REAL64), parameter :: varpi10 = 0.4817
    real(REAL64), parameter :: varpi25 = 0.3478
    real(REAL64), parameter :: nu25 = 0.5985
    real(REAL64), parameter :: nu10 = 0.7631
    real(REAL64), parameter :: nu0 = 0.8717
    real(REAL64), parameter :: s0 = sqrt(nu0)/1.02 !0.9338163922
    real(REAL64), parameter :: s10 = 0.8790363864
    real(REAL64), parameter :: s25 = 0.784929715
    real(REAL64), dimension(:), allocatable :: std
    real(REAL64) :: nu, varpi, sig, sc
    integer :: N, k

    N = ndata
    k = N
    allocate(std(nsamples))
    if( step == 4 ) then
       nu = nu25
       varpi = varpi25
       sc = s25
    else if( 8 <= step .and. step <= 10 ) then
       nu = nu10
       varpi = varpi10
       sc = s10
    else
       nu = nu0
       varpi = varpi0
       sc = s0
    end if

    d = sqrt(2.0*ndata)
    t0 = tpar - 3
    t1 = tpar + 3
    dt = d/20

    std = sqrt(real(k))*(t-tpar)*cscale/sigpar**2
    sig = sqrt(nu0/nu)*sigpar
    sig = sqrt(nu)*sigpar
    sig = 1 / sqrt(hess(1,1))
    sig = avgsig !avgs / sqrt(nu)
    !    sig = sigpar
    sig = sc / sqrt(nu)
    std = sqrt(real(N))*(t - tpar) / sig
!    std = sqrt(real(N))*(t - tpar) / (sc/sqrt(nu))
    write(*,*) '***',sig,sqrt(nu)*sigpar,avgs
    call histogram(std,t0,t1,nbin,'/tmp/testsr1mean_std_t',.true.)
    open(1,file='/tmp/testsr1mean_std_t_data')
    write(1,*) std
    close(1)

    std = sqrt(real(k))*(s-avgs)/sqrt(1.0/hess(2,2))*cscale/sigpar
    sig = sqrt(nu/nu0*(1 + varpi)/(1 + varpi0))*sigpar
    sig = 1 / sqrt(hess(2,2))!*cscale
    sig = sqrt(nu/(1 + varpi))*sigpar
    sig = sc/sqrt(1 + varpi)*sigpar
!    sig = avgs/sqrt(1 + varpi)
!    sig = avgs/sqrt(1 + avgvarpi)
!    sig = 0.636!0.486
    write(*,*) '****',sig,sqrt(nu)*sigpar,N,sc/sqrt(nu)
!    std = sqrt(real(N))*(s - avgs) / sig
!    std = sqrt(real(N))*(s - sqrt(nu)*sigpar) / sig
    std = sqrt(real(N))*(s - sc*sigpar) / sig
    call histogram(std,t0,t1,nbin,'/tmp/testsr1mean_std_s',.true.)
    open(1,file='/tmp/testsr1mean_std_s_data')
    write(1,*) std
    close(1)

!  chi2 = s**2*pshampel_real64((t-0)/sig)**2/dpshampel_real64((t-0)/sig)
!  chi2 = ((dataset-0)/s)**2
!  do i = 1, nsamples
!     l = i*ndata
!     data = dataset(1+(l-ndata):l)
!     res = (data-t(i)) / s(i)
!     chi2(i) = s(i)**2*sum(pshampel_real64(res)**2)/sum(dpshampel_real64(res))
!     chi2(i) = s(i)**2*k/sum(dpshampel_real64(res))
!     chi2(i) = s(i)**4*sum(pshampel_real64(res)**2)/sum(dpshampel_real64(res))**2*k
!     write(*,*) chi2(i),sig(i)**2*(ndata-1)*avgnu/avgs**2
!     chi2(i) = chi2(i)*hess(1,1)
!     chi2(i) = chi2(i)*avgnu/avgs**2
!     chi2(i) = chi2(i) / (k-1)
!  end do
!  t0 = 0
!  t1 = 100
!  call histogram(chi2,t0,t1,nbin,'/tmp/ddd',.true.)

  !  write(*,*) chi2
!  call histogram((chi2*avgnu/avgs**2-k)/sqrt2k*avgs/sigpar**2,t0,t1,nbin,'/tmp/testsr1mean_std_chi2',.true.)
    !  call histogram((chi2-k)/sqrt2k*avgs/sigpar,t0,t1,nbin,'/tmp/testsr1mean_std_chi2',.true.)
    k = N - 1
    sqrt2k = sqrt(2.0*k)
    std = (chi2*avgnu/avgs**2*sigpar**2-k)/sqrt2k*avgnu
    std = (chi2*sigpar**2/avgsig**2 - k)/sqrt2k*sqrt(nu)!/sigpar
!    std = (chi2 - k)/sqrt2k!*sqrt(nu)!/sigpar
    call histogram(std,t0,t1,nbin,'/tmp/testsr1mean_std_chi2',.true.)
    open(1,file='/tmp/testsr1mean_std_chi2_data')
    write(1,*) std
    close(1)

!  chi2med = sqrt(2.0*k - 1)
!  call histogram(sqrt(2*chi2)-chi2med,t0,t1,nbin,'/tmp/testsr1mean_std_chi2',.true.)
  ! Wilson-Hilferty transformation
!  chi2var = 2.0/(9.0*k)
!  chi2med = 1 - chi2var
!  call histogram(((chi2/k)**(1.0/3.0)-chi2med)/sqrt(chi2var), &
!       t0,t1,nbin,'/tmp/testsr1mean_std_chi2',.true.)

  end block

  block
    real(REAL64), parameter :: tmin =-3e0
    real(REAL64), parameter :: tmax = 3e0
    real(REAL64), parameter :: smin = 0e0
    real(REAL64), parameter :: smax = 3e0
    type(robust_likelihood_hampel), target :: hampelfun
    type(robust_likelihood_tukey), target :: tukeyfun
    type(robust_likelihood_huber), target :: huberfun
    type(nonrobust_likelihood_square), target :: squarefun
    class(robust_likelihood_REAL64), pointer :: fcn


    if( robfun == 'huber' ) then
       fcn => huberfun
    else if( robfun == 'tukey' ) then
       fcn => tukeyfun
    else if( robfun == 'square' ) then
       fcn => squarefun
    else !if( robfun == 'hampel' ) then
       fcn => hampelfun
    end if
    fcn%data => dataset

!    call rland(fcn,tmin,tmax,smin,smax,'/tmp/testsr1mean_hist_rland')

    open(1,file='/tmp/testsr1mean_hist_rland.par')
    write(1,'(a,g0)') 'TMIN = ',tmin
    write(1,'(a,g0)') 'TMAX = ',tmax
    write(1,'(a,g0)') 'SMIN = ',smin
    write(1,'(a,g0)') 'SMAX = ',smax
    close(1)

  end block

  open(1,file='/tmp/testsr1mean_hist')
  do i = 1, nsamples
     write(1,*) t(i),s(i),sig(i),chi2(i)
  end do
  close(1)

  open(1,file='/tmp/testsr1mean_hist.par')
  write(1,'(i0)') ndata
  close(1)

  open(1,file='/tmp/testsr1mean_hist_rank')
  do i = 1, nsamples
     write(1,*) med(i),Qn(i)
  end do
  close(1)

  open(1,file='/tmp/testsr1mean_hist_data')
  write(1,'(a)') '# 6*N(0,1) + 1*N(-10,1), Hampel'
  write(1,*) ndata*nsamples
  write(1,*) dataset
  close(1)

  write(error_unit,*) "done"


contains

  subroutine rmean_new(t,s,sig,med,Qn,xhess,data,info)

!    use trustsr1mean
    use noise_REAL64
!    use quantilmean

    real(REAL64), intent(out) :: t,s,sig,med,Qn
    real(REAL64), dimension(:,:), intent(out) :: xhess
    real(REAL64), dimension(:), intent(out), target :: data
    integer, intent(out) :: info

!    character(len=80) :: msg
!    real(REAL64) :: stderr, stddev
    integer :: i, istep

    call gnoise_fill(tpar,sigpar,data)
    if( step > 0 ) then
       istep = step
       do i = 2, ndata, istep
          data(i) = gnoise(-100.0_REAL64,100.0_REAL64)
          !       write(*,*) data(i)
       end do
    end if

!    open(3,file='/tmp/gnoise')
!    write(3,'(g0)') data
!    close(3)
!    stop

    ! save data
    open(3,file='/tmp/tsr1mean_data',form='unformatted')
    write(3) ndata
    write(3) data
    close(3)


!    goto 33
!    call qmean(data,med,stdsig=Qn)
!    t = qmed(data)
!    s = byQ(data)
!    med = t
!    Qn = s
!    Qn = qmed(abs(data-t)) / 0.6745
!    s = Qn
!    s = 0.8*s
!    call trmean(data,t,s,stderr,stddev,robfun,info,msg,.true.)
!    sig = stddev!*sqrt(real(ndata,REAL64))
!    sig = sqrt(sum(data**2)/(size(data)))
!    sig = sqrt(sum((data-t)**2))
!    sig = sum((data-t)**2)
    !    t = sum(data)/ndata

!    if( info /= 0 ) then
!       write(*,*) info, msg
!       stop
!    end if

!    if( s < 0.1 ) stop 's < 0.1'

!    if( info /= 0 ) stop 'trmean() failed.'
!33  continue

    block
      use robustmean_REAL64
      character(len=*), parameter :: report_file = "/tmp/tsr1mean_hist_report"
      real(REAL64) :: mean,stderr,stdsig,scale
      real(REAL64), dimension(2,2) :: mhess, happrox
      real(REAL64), dimension(size(data)) :: dat
      integer :: iflag,miter,mfev,mrst,mout,midf,mhrd,eq,iostat
      character(len=666) :: line

      if( winsor ) then
         block
           use quantilmean_REAL64
           real(REAL64) :: med, stderr, stdsig, w
           real(REAL64), dimension(size(data)) :: wdata
           call qmean_REAL64(data,med,stderr,stdsig)
           w = 1.5*stdsig
           wdata = data - med
           where( abs(wdata) > 1.5*stdsig )
              dat = sign(w,wdata)
           elsewhere
              dat = data
           end where
         end block
      else
         dat = data
      end if

      mean = tpar
      scale = 0.9*sigpar
      mean = med
      scale = Qn
      mean = sum(data)/size(data)
      scale = sqrt(sum((data-mean)**2)/(size(data)))
      mean = 0
      scale = 1
      call rmean_REAL64(dat,mean,stderr,stdsig,scale,reltol,robfun,sigic=.true.,&
           rankinit=.true.,verbose=.true.,flag=iflag,report=report_file)

      open(13,file=report_file,status='old')
      do
         read(13,'(a)',iostat=iostat) line
         if( iostat == IOSTAT_END ) exit
         eq = index(line,'=') + 1
         if( index(line,'t0') > 0 ) read(line(eq:),*) med
         if( index(line,'s0') > 0 ) read(line(eq:),*) Qn
         if( index(line,'niter') > 0 ) read(line(eq:),*) miter
         if( index(line,'nfev') > 0 ) read(line(eq:),*) mfev
         if( index(line,'nrst') > 0 ) read(line(eq:),*) mrst
         if( index(line,'nout') > 0 ) read(line(eq:),*) mout
         if( index(line,'nidf') > 0 ) read(line(eq:),*) midf
         if( index(line,'nhrd') > 0 ) read(line(eq:),*) mhrd
         if( index(line,'hess') > 0 ) read(line(eq:),*) mhess
         if( index(line,'happrox') > 0 ) read(line(eq:),*) happrox
      end do
      close(13)

      t = mean
      s = scale
      sig = stdsig
      info = iflag
      xhess = mhess
      xhess = happrox
      hess = hess + mhess
      niter = niter + miter
      nfev = nfev + mfev
      nrst = nrst + mrst
      nout = nout + mout
      nidf = nidf + midf
      nhrd = nhrd + mhrd
      if( info >= 5 ) stop 'info >= 5'
    end block

  end subroutine rmean_new

  subroutine rland(fcn,x0,x1,y0,y1,filename)

    class(robust_likelihood_REAL64), intent(in) :: fcn
    real(REAL64), intent(in) :: x0,x1,y0,y1
    character(len=*), intent(in) :: filename
    integer, parameter :: nx = 200
    integer, parameter :: ny = 200

    real(REAL64) :: x,dx,y,dy
    integer :: i,j

    dx = (x1 - x0) / (nx - 1)
    dy = (y1 - y0) / (ny - 1)

    open(1,file=filename)
    do j = 1, ny
       y = y0 + dy*j
       do i = 0, nx
          x = x0 + dx*i
          write(1,*) x,y,fcn%func(x,y)
       end do
       write(1,*)
    end do
    close(1)

  end subroutine rland

  subroutine density_loc(t0,t1,dt,N,somefile)

    real(REAL64), intent(in) :: t0,t1,dt
    integer, intent(in) :: N
    character(len=*), intent(in) :: somefile
    real(REAL64) :: t,sig,f
    integer :: i

    sig = sigpar / sqrt(real(N)) !/ cscale
    open(1,file=somefile)
    do i = 1, int((t1 - t0) / dt)
       t = t0 + i*dt
       f = exp(-(t - tpar)**2/2/sig**2) / (sqrt(2*pi) * sig)
       write(1,*) t,f
    end do
    close(1)

  end subroutine density_loc

  subroutine density_s(s0,s1,ds,N,somefile)

    real(REAL64), intent(in) :: s0,s1,ds
    integer, intent(in) :: N
    character(len=*), intent(in) :: somefile
    real(REAL64) :: s,sig,f,f1,t,sig0,hss
    integer :: i

    t = 0
!    sig = sigpar* beta**2 / sqrt(2.0*real(N))
!    sig = s0 / sqrt((1 + gama)*N)
    sig = sigpar * sqrt(beta/((1 + gama)*N))
    sig0 = sqrt(beta)
    sig0 = sqrt(avgnu)
    sig = avgs * sqrt(nu/((1 + avgvarpi)*N))
    sig = avgs / sqrt((1 + avgvarpi**2)*N) * sqrt(N/(N - 1.0))
    hss = hess(1,1) / (hess(1,1)*hess(2,2) - hess(1,2)**2)
    sig = sqrt(N/(N-1.0)/hess(2,2)/N)
    sig = sqrt(1.0/hess(2,2)/N)/cscale
!    sig = sqrt(N/(N-1.0)/hss/N)
    open(1,file=somefile)
    do i = 1, int((s1 - s0) / ds)
       s = s0 + i*ds
!       f1 = exp(-(t - sig0*sigpar)**2/2/sig**2) / (sqrt(2*pi) * sig)
       f = exp(-(s - sig0*sigpar)**2/2/sig**2) / (sqrt(2*pi) * sig)
       f = exp(-(s - avgs)**2/2/sig**2) / (sqrt(2*pi) * sig)
!       f = exp(-ipshampel((t - sig0*sigpar)/sig)) / (hnorm * sig)
!       f = exp(-ipshampel((s - sig0*sigpar)/sig)) / (hnorm * sig)
       write(1,*) s,f
    end do
    close(1)

  end subroutine density_s

  subroutine density_sig(t0,t1,dt,N,somefile)

    real(REAL64), intent(in) :: t0,t1,dt
    integer, intent(in) :: N
    character(len=*), intent(in) :: somefile
!    integer, parameter :: m = 10000
    real(REAL64), dimension(:), allocatable :: t,f,g,h
    real(REAL64) :: s,sig,logk,u
    integer :: i, j, k, m

    !    m = int((t1 - t0) / dt)
    m = int(t1/dt)
    allocate(t(m),f(m),g(-m:m),h(m))

    k = N !- 1
    logk = (k/2.0 - 1)*log(2.0) + log_gamma(k/2.0)
    do i = 1, size(h)
       t(i) = 0*t0 + i*dt
       u = (k-1)*log(t(i)) - t(i)**2/2 - logk
       h(i) = exp(u)
    end do
    h = h * dt

    sig = sigpar*sqrt(beta/((1+gama)*N))
!    sig = 1.071/sqrt(2.0*k)
    g(0) = 1 / sqrt(2*pi)/sig
    do i = 1, m
       g(i) = exp(-(t(i)-0)**2/2/sig**2)/sqrt(2*pi)/sig
       g(-i) = g(i)
    end do
    g = g * dt

    ! the convolution
    do i = 1, size(f)
       s = 0
       do j = -2*m, 2*m
          if( 1 <= j .and. j <= m .and. -m <= i-j .and. i-j <= m ) then
             s = s + h(j)*g(i-j)
          end if
       end do
       f(i) = s / dt
    end do
    f = g / dt

    open(1,file=somefile)
    do i = int(t0/dt+1), size(t)
       write(1,*) t(i),f(i)
    end do
    close(1)

  end subroutine density_sig


  subroutine density_chi2(t0,t1,dt,N,somefile)

    real(REAL64), intent(in) :: t0,t1,dt
    integer, intent(in) :: N
    character(len=*), intent(in) :: somefile
    real(REAL64), dimension(:), allocatable :: t,f,g,h
    real(REAL64) :: x,s,sig,logk2,u2,chi2
    integer :: i, j, k, m

    k = N
    logk2 = (k/2.0) * log(2.0) + log_gamma(k/2.0)    ! chi2

    open(1,file=somefile)
    do i = 1, int((t1 - t0) / dt)
       x = t0 + i*dt
       u2 = (k/2.0-1.0)*log(x) - x/2 - logk2   ! chi2
       chi2 = exp(u2)
       write(1,*) x,chi2
    end do
    close(1)

    return


!    m = int((t1 - t0) / dt)
    m = int(t1/dt)
    allocate(t(m),f(m),g(-m:m),h(m))

    k = N !- 1
    logk2 = (k/2.0) * log(2.0) + log_gamma(k/2.0)    ! chi2
    do i = 1, size(h)
       t(i) = 0*t0 + i*dt
       u2 = (k/2.0-1.0)*log(t(i)) - t(i)/2 - logk2   ! chi2
       h(i) = exp(u2)
    end do
    h = h *dt

!    sig = beta/sqrt(2*real(N))
    sig = sigpar*sqrt(beta/((1+gama)*N))
!    sig = 1.071/sqrt(2.0*k)
    g(0) = 1 / sqrt(2*pi)/sig
    do i = 1, m
       g(i) = exp(-(t(i)-0)**2/2/sig**2)/sqrt(2*pi)/sig
       g(-i) = g(i)
    end do
    g = g * dt

    ! the convolution
    do i = 1, size(f)
       s = 0
       do j = -2*m, 2*m
          if( 1 <= j .and. j <= m .and. -m <= i-j .and. i-j <= m ) then
             s = s + h(j)*g(i-j)
          end if
       end do
       f(i) = s / dt
    end do

    open(1,file=somefile)
    do i = int(t0/dt+1), size(t)
       write(1,*) t(i),f(i)
    end do
    close(1)

  end subroutine density_chi2

  real(REAL64) function stirlig(x)

    real(REAL64), intent(in) :: x
    real(REAL64), parameter :: pi = 3.141592653589793_REAL64
    real(REAL64), parameter :: e = 2.718281828459045_REAL64

    stirlig = sqrt(2*pi*x)*(x/e)**x

  end function stirlig


end program testsr1mean_hist
