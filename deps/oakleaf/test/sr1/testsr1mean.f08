
! * napsat, ktere metriky pro s nefungovali: KS, von Misses, Fisherova informace,..
!   quantily z cdf,.. uz je jasne i proc
! * prozkoumat, jak windowreizace pripadne Huber davaji vychylena minima, i kdyz jsou
!   jednoznacne
! * downtown (1966) jako odhad scale je uplne nepouzitelny



module gold_ridge

  use rlike_REAL64
  use fminim_REAL64
  use iso_fortran_env
  implicit none

  type, extends(goldfinger_function_REAL64) :: gold_fun
     class(robust_likelihood_REAL64), pointer :: fcn
     integer :: minmax
     real(REAL64), dimension(2) :: q,p
     integer :: nfev
   contains
     procedure :: f => goldfcn
  end type gold_fun

contains

  real(REAL64) function goldfcn(this,x)

    class(gold_fun) :: this
    real(REAL64), intent(in) :: x
    real(REAL64), dimension(2) :: r

    associate( q => this%q, p => this%p, minmax => this%minmax, nfev => this%nfev )

      r = q + x*p
      goldfcn = minmax * this%fcn%func(r(1),r(2))
      nfev = nfev + 1

    end associate

  end function goldfcn

end module gold_ridge

program testsr1mean

  ! test optimization

  !  use symmetric_rank_1
  use rfun_REAL64
  use rfun_hampel_REAL64
  use rlike_REAL64
!  use sr1mean_function
!  use trustsr1mean
!  use quantilmean
!  use zeroin_smin
  use noise_REAL64
  use iso_fortran_env
  implicit none

!  type(tsr1mean_tukey), target :: tukeymean
!  type(tsr1mean_huber), target :: hubermean
!  type(tsr1mean_square), target :: squaremean
!  type(tsr1mean_experimental), target :: bathmean
!  class(tsr1mean_function), pointer :: funmean

!  type(zeroin_data) :: funzero
  integer, parameter :: ndata = 420000
  real(REAL64), parameter :: tpure = 0
  real(REAL64), parameter :: spure = 1e6
  real(REAL64), parameter :: tcont =-10
  real(REAL64), parameter :: scont = 10
  character(len=80), parameter :: robfun = 'hampel'
  logical, parameter :: switch_on_generator = .true.
  logical, parameter :: patolog = .false.
  logical, parameter :: plot_graphs = .true.
  integer, parameter :: step = 0
!  real(REAL64), parameter :: tmin = -3
!  real(REAL64), parameter :: tmax = 3
!  real(REAL64), parameter :: smin = 0e0
!  real(REAL64), parameter :: smax = 3
  real(REAL64), parameter :: tmin = -3
  real(REAL64), parameter :: tmax = 3
  real(REAL64), parameter :: smin = 0e0
  real(REAL64), parameter :: smax = 3
  real(REAL64), parameter :: reltol = 1.0 / 42.0 * 1e-3

  real(REAL64), dimension(:), allocatable, target :: data
  character(len=666), parameter :: datafile = '/tmp/tsr1mean_data'
!  real(REAL64), dimension(3), parameter :: q = [0.25, 0.5, 0.75]
!  real(REAL64), dimension(3), parameter :: qq = [0.4, 0.6745, 0.87]
!  real(REAL64), dimension(3), parameter :: qq = [0.4, 0.5, 0.8]
!  real(REAL64), dimension(size(q)) :: f,g
!  real(REAL64), dimension(2) :: topt,sopt

!  real(REAL64), dimension(2,2) :: cov
!  real(REAL64) :: delta, maxdev,t,s,dt, chi2par, chi2t,df,rmin,rmax, &
!       negtukey_t,fminsig_t, curvature
  !  integer :: maxfev, info, i, ndat, neval, nopt, iopt,istep
  integer :: i, ndat, istep
!  character(len=80) :: msg = ''
!  logical :: verbose, quasi

  open(1,file='/tmp/testsr1mean.par')
  write(1,'(a,g0)') 'TMIN = ',tmin
  write(1,'(a,g0)') 'TMAX = ',tmax
  write(1,'(a,g0)') 'SMIN = ',smin
  write(1,'(a,g0)') 'SMAX = ',smax
  close(1)

  if( switch_on_generator ) then
     ! generate data with N(tpure,spure)
     allocate(data(ndata))
     call gnoise_fill(tpure,spure,data)
!     frac = 0.25
!     step = ndata / max(int(frac*ndata),1)
     !     step = 5
     write(*,*) 'step:',step
     if( step > 0 ) then
        ! contamine the data with N(tcont,scont)
        istep = step
        do i = 2, ndata, istep
           data(i) = gnoise(tcont,scont)
!           write(*,*) 'bad point:',real(data(i))
        end do
!        do i = 3, ndata, step
!           data(i) = gnoise(10.0_REAL64,1.0_REAL64)
!           write(*,*) 'bad point:',real(data(i))
!        end do
     end if
!     write(*,*) 'data:',real(data)

     ! save data
     open(3,file=datafile,form='unformatted')
     write(3) ndata
     write(3) data
     close(3)

  else
     ! saved data
     if( patolog ) then
        write(*,*) 'PATOLOG in ACTION'
        open(3,file='data_20260109',status='old')
        read(3,*)
        read(3,*) ndat
        allocate(data(ndat))
        read(3,*) data
        close(3)
     else
        open(3,file=datafile,form='unformatted')
        read(3) ndat
        allocate(data(ndat))
        read(3) data
        close(3)
     end if
!     write(*,*) 'data:',real(data)
  end if

!  funmean%data => data
!  hubermean%data => data
!  pshubermean%data => data

!  sr1fun%data => data
!  sr1fun%ndata = n
!  allocate(sr1fun%psi(n),sr1fun%dpsi(n),sr1fun%res(n))
!  sr1fun%irob => isquare
!  sr1fun%frob => square
!  sr1fun%drob => dsquare

!!$    if( kind == 'Huber' .or. kind == 'psHuber' ) then
!!$       ifcn => ihuber
!!$       fcn => huber
!!$       dfcn => dhuber
!!$    else if( kind == 'Tukey' ) then
!!$       ifcn => itukey
!!$       fcn => tukey
!!$       dfcn => dtukey
!!$    else if( kind == 'square' ) then
!!$       ifcn => isquare
!!$       fcn => square
!!$       dfcn => dsquare
!!$    else
!!$       stop 'Specified an unknown kind of function.'
!!$    end if

  !  allocate()
!  if( robfun == 'huber' ) then
!     funmean => hubermean
!  else if( robfun == 'hampel' ) then
!     funmean => bathmean
!  else !if( robfun == 'tukey' ) then
!     funmean => tukeymean
!  end if
!  funmean%data => data
!  hubermean%data => data
!  funzero%funmean => funmean


!  call oldway

!!$  goto 99
!!$  block
!!$    use quantilmean
!!$    use kstestscale
!!$
!!$    real(REAL64) :: mean,stderr,stddev,scale
!!$    integer :: info
!!$    character(len=80) :: msg
!!$
!!$    call qmean(data,mean,stdsig=scale)
!!$    write(*,*) 'qmean:',mean,scale
!!$!    call byQ(data,scale)
!!$!    scale = 0.75*scale
!!$
!!$    open(2,file='/tmp/tsr1mean_step')
!!$    call trmean(data,mean,scale,stderr,stddev,robfun,info,msg,.true.)
!!$    write(*,*) 'FINISH trmean(): `',trim(msg),"' flag=",info
!!$    close(2)
!!$  end block
!!$99 continue

  block
    use robustmean_REAL64
!    use sqpfreport
!    type(sqpfmean_diagnosis) :: diagnose
    real(REAL64) :: mean,stderr,stdsig,scale
    real :: tstart,tfinish
    integer :: flag

!    open(17,file="/tmp/tsr1mean_step")
!    diagnose%unit = 17
    mean =-2.5e-10
    scale = 1
    call cpu_time(tstart)
    call rmean_REAL64(data,mean,stderr,stdsig,scale,reltol,robfun,rankinit=.true.,&
         verbose=.true.,flag=flag,report="/tmp/tsr1mean_step")
    call cpu_time(tfinish)
    write(*,'(a,en13.3)') 'Elapsed time [s]: ',tfinish-tstart
!    write(*,*) 'flag=',flag
    !diagnose=diagnose)
!    close(17)

    block
      use rfun
      real(REAL64), dimension(:), allocatable :: res
      integer :: N
!      real(REAL64), parameter :: beta = 0.8712975324409_REAL64
!      real(REAL64), parameter :: gamma = 0.60572050104634123
!      real(REAL64) :: sig,err
      N = size(data)
      allocate(res,source=data)
      res = (data - mean)/scale
      if( robfun == 'tukey' ) then
         write(error_unit,*) "Identities: ", &
              '∑ ψ²  =',sum(tukey(res)**2),&
              "∑ ψ' =",sum(dtukey(res)),&
              "∑ ψ.r =",sum(tukey(res)*res),&
              "∑ ψ =",sum(tukey(res))
      else if( robfun == 'huber' ) then
         write(error_unit,*) "Identities: ", &
              '∑ ψ²  =',sum(pshuber(res)**2),&
              "∑ ψ' =",sum(dpshuber(res)),&
              "∑ ψ.r =",sum(pshuber(res)*res),&
              "∑ ψ =",sum(pshuber(res))
      else
         write(error_unit,*) "Identities: ", &
              '∑ ψ²  =',sum(pshampel(res)**2),&
              "∑ ψ' =",sum(dpshampel(res)),&
              "∑ ψ.r =",sum(pshampel(res)*res),&
              "∑ ψ =",sum(pshampel(res))
         write(error_unit,*) "∑ ψ'. ψ =",sum(dpshampel(res)*pshampel(res)), &
              sum((dpshampel(res)-0.9)**2)/N,sum(dpshampel(res)**2)/N
      end if
!           scale / (sqrt(beta)*spure)
!      sig = scale*N/sum(dpshampel(res))*sqrt(sum(pshampel(res)**2)/(N-1))
!      err = sig / sqrt(1.0*N)
!      write(error_unit,*) 'stdsig=',sig,' stderr=',err
!      write(error_unit,*) 'The asymptotic beta=',sum(pshampel(res)**2)/N
    end block

!!$    write(*,*) 'xxx',sum((data/scale)**2) !/ size(data)
!!$    block
!!$      use rfun_hampel
!!$
!!$      write(*,*) mean,scale
!!$      write(*,*) sum(dhampel((data-mean)/scale)**1)/ndata,&
!!$           sum(hampel((data-mean)/scale)**2)/ndata,&
!!$           sum(((data-mean)/scale)**1)/ndata,&
!!$           '*',sum(dhampel((data-mean)/scale)*(data-mean)/scale)/ndata*sqrt(1.0*ndata),&
!!$           sum(dhampel((data-mean)/scale)*(data-mean)),'<=',&
!!$           sqrt(sum(dhampel((data-mean)/scale)**2)),'.',&
!!$           sqrt(sum(((data-mean)/scale)**2)),&
!!$           sqrt(sum(data**2)/ndata)*sqrt(1.0*ndata)
!!$    end block
  end block

!  stop

  block
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
    fcn%data => data

    if( plot_graphs ) then
       write(error_unit,'(a)',advance="no") "Calculating graphs .. "
       call graph(fcn,tmin,tmax,smin,smax,'/tmp/tsr1mean')
       !  call graph(funmean,290.0,310.0,5.0,15.0,'/tmp/tsr1mean')
       !    call graphd(funmean,data(1)-1e-14,data(1)+1e-14,0.0d0,1d-14,'/tmp/tsr1meanx')
       !    call graph_marginal(funmean,-5.0,5.0,0.0,50.0,'/tmp/tsr1mean_marginal')
       call graph_hess(fcn,tmin,tmax,smin,smax,'/tmp/tsr1mean_hess')
       !  call graph(funmean,-.5,.5,0.0,1.0,'/tmp/tsr1mean')
       call graph_rho(fcn,tmin,tmax,smin,smax,'/tmp/tsr1mean_rho')
       call graph_jac(fcn,tmin,tmax,smin,smax,'/tmp/tsr1mean_jac')
       !    call graph_smin(funzero,-5.0_REAL64,5.0_REAL64,5.0_REAL64,'/tmp/tsr1mean_smin')
       call graph_ridge(fcn,tmin,tmax,smin,smax,'/tmp/tsr1mean_ridge')
       !  call graphreg(funmean,-5.0,5.0,0.0,5.0,'/tmp/tsr1mean_reg')
       !  call graphpole(funmean,-0.2,-0.1,0.0,.01,'/tmp/tsr1mean_pole')
       !  call graphcliffs(funmean,-5.0,5.0,0.0,5.0,'/tmp/tsr1mean_cliffs')
       !  call graph_stringart(funmean,-5.0,5.0,0.0,5.0,'/tmp/tsr1mean_stringart')
       !  call graph_ctn(funmean,-5.0,5.0,0.0,5.0,'/tmp/tsr1mean_ctn')
       !  call graph(funmean,-1.0,0.0,0.0,0.1,'/tmp/tsr1mean')
       !  call graphdev(funmean,-5.0,5.0,0.0,5.0,'/tmp/sr1mean_stddev')
       write(error_unit,*) "done"
    end if
  end block


!  deallocate(sr1fun%psi,sr1fun%dpsi,sr1fun%res)


!  open(2,file='/tmp/data_20260109',status='new')
!  write(2,'(a)') '#  Hampel: 4*N(0,1) + 1*N(0,1e6)'
!  write(2,*) ndat
!  write(2,*) data
!  close(2)

!  block

!    use quantilmean
!    real(REAL64), dimension(:), allocatable :: x,y
!    allocate(x(ndat),y(ndat))
!    call ecdf(abs(data - 0.034),x,y)

!    write(*,*) x
!    write(*,*) y
!    stop

!  end block


contains

!!!!!  subroutine oldway()
!!!!!
!!!!!!  if( patolog ) then
!!!!!  if( .true. ) then
!!!!!     block
!!!!!       use quantilmean
!!!!!
!!!!!       real(REAL64) :: t,s,ss
!!!!!       real(REAL64), dimension(:), allocatable :: w,rho,res
!!!!!       integer :: i,j,n
!!!!!
!!!!!       !       t = 0.83
!!!!!       call qmean(data,t)
!!!!!
!!!!!       n = size(data)
!!!!!       allocate(w(n),rho(n),res(n))
!!!!!
!!!!!       open(1,file='/tmp/patfun')
!!!!!       open(2,file='/tmp/testsr1mean_entropy')
!!!!!       do i = 1, 300
!!!!!          s = 0.05*i
!!!!!          if( robfun == 'huber' ) then
!!!!!!             res = (data-t)/s
!!!!!             ss = sum(ihuber((data-t)/s))/n + log(s)
!!!!!             write(1,*) s,(ihuber((data(j)-t)/s)+log(s),j=1,n),ss
!!!!!             !ss = sum(((data-t)/s)**2)/n/2+ log(s)
!!!!!             !write(1,*) s,(((data(j)-t)/s)**2/2+log(s),j=1,n),ss
!!!!!!             ss = sqrt(n* s**2 * sum(tukey(res)**2) / sum(dtukey(res))**2)
!!!!!!             res = itukey(res)
!!!!!!             rho = 2*itukey(res)
!!!!!             rho = 2*ihuber((data-t)/s)
!!!!!             w = exp(-rho)/s/sqrt(2*3.14)
!!!!!!             ss = - sum(rho*exp(-rho)) / n
!!!!!             ss = sum(w*log(w)) / n
!!!!!             write(2,*) s,(w(j)*log(w(j)),j=1,n),ss
!!!!!!             write(1,*) s,(rho(j)*exp(-rho(j)),j=1,n),ss
!!!!!!             write(1,*) s,(((data(j)-t)/s)**2/2+log(s),j=1,n),ss
!!!!!          else if( robfun == 'hampel' ) then
!!!!!             res = min(max((data-t)/s,-3.0_REAL64),3.0_REAL64)
!!!!!!             ss = sum(ipshampel(res))/n + log(0.9*s)
!!!!!!             ss = sum(pshampel(res)**2)/n - 0.9
!!!!!             ss = sum(ipshampel((data-t)/s))/n + log(s)
!!!!!             !             ss = funmean%val(t,s)
!!!!!
!!!!!!             ss = sum(min(pshampel((data-t)/s)**2,real(9,REAL64))) + log(s)
!!!!!!             write(1,*) s,(min(pshampel((data(j)-t)/s)**2,9.0_REAL64)+log(s),j=1,n),ss
!!!!!             write(1,*) s,(ipshampel((data(j)-t)/s)+log(s),j=1,n),ss
!!!!!!             write(1,*) s,(ipshampel(res(j))+log(0.9*s),j=1,n),ss
!!!!!             !             write(1,*) s,(pshampel(res(j))**2-0.9,j=1,n),ss
!!!!!
!!!!!             rho = ipshampel((data-t)/s)
!!!!!             w = exp(-rho)/s/sqrt(2*3.14)
!!!!!             ss = sum(w*(log(w) - 1)) / n
!!!!!             write(2,*) s,(w(j)*(log(w(j))-1),j=1,n),ss
!!!!!!             ss = sum(s*data)/n + s**4/2 - log(s)
!!!!!!             write(2,*) s,(s*data(j) + s**4/2 - log(s),j=1,n),ss
!!!!!          else
!!!!!             ss = sum(itukey((data-t)/s))/n + log(s)
!!!!!             write(1,*) s,(itukey((data(j)-t)/s)+log(s),j=1,n),ss
!!!!!             rho = 2*itukey((data-t)/s)
!!!!!             w = exp(-rho)/s/sqrt(2*3.14)
!!!!!             ss = sum(w*log(w)) / n
!!!!!             write(2,*) s,(w(j)*log(w(j)),j=1,n),ss
!!!!!          end if
!!!!!!          ss = sum(ipshampel((data-t)/s))/n + log(s)
!!!!!!          write(1,*) s,(ipshampel((data(j)-t)/s)+log(s),j=1,n),ss
!!!!!       end do
!!!!!       close(2)
!!!!!       close(1)
!!!!!     end block
!!!!!  end if
!!!!!
!!!!!!!$  block
!!!!!!!$    use quicksorter
!!!!!!!$    real(REAL64), dimension(:), allocatable :: x
!!!!!!!$    integer :: n
!!!!!!!$    allocate(x(size(data)))
!!!!!!!$    x = data
!!!!!!!$    call qsort(x)
!!!!!!!$    n = size(x) / 2
!!!!!!!$    if( mod(size(x),2) == 0 ) then
!!!!!!!$       t = (x(n) + x(n+1)) / 2
!!!!!!!$    else
!!!!!!!$       t = x(n)
!!!!!!!$    end if
!!!!!!!$    x = abs(data - t)
!!!!!!!$    call qsort(x)
!!!!!!!$    if( mod(size(x),2) == 0 ) then
!!!!!!!$       s = (x(n) + x(n+1)) / 2
!!!!!!!$    else
!!!!!!!$       s = x(n)
!!!!!!!$    end if
!!!!!!!$    write(*,*) 'median, mad:',real(t),real(s)
!!!!!!!$  end block
!!!!!  call quantilefun(data,q,f,info)
!!!!!  t = (f(2) + (f(3) + f(1))/2)/2
!!!!!  s = (f(3) - f(1)) / (2*0.6745)
!!!!!  write(*,*) 'initial t:',real(t),real(s)
!!!!!  t = f(2)
!!!!!
!!!!!!!$  block
!!!!!!!$    real(REAL64), dimension(1), parameter :: q = [0.5]
!!!!!!!$    real(REAL64), dimension(size(q)) :: f
!!!!!!!$    real(REAL64), dimension(:), allocatable :: z,y
!!!!!!!$    integer :: n,i,j,k
!!!!!!!$    n = size(data)
!!!!!!!$    allocate(z(n*n),y(n))
!!!!!!!$    k = 0
!!!!!!!$    do i = 1, n
!!!!!!!$       do j = 1,n
!!!!!!!$          k = k + 1
!!!!!!!$          z(k) = (data(i) + data(j))/2
!!!!!!!$       end do
!!!!!!!$    end do
!!!!!!!$    call quantilefun(z(1:k),q,f,info)
!!!!!!!$    y = 2*f(1) - data
!!!!!!!$    z(1:n) = y - data
!!!!!!!$    call quantilefun(z(1:n),q,f,info)
!!!!!!!$    write(*,*) f
!!!!!!!$    stop
!!!!!!!$  end block
!!!!!
!!!!!!  call quantilefun(abs(data-t),q,g,info)
!!!!!!  s = g(2) / 0.6745
!!!!!!  write(*,*) abs(data-t)
!!!!!
!!!!!  write(*,*) 'f:',real(f)
!!!!!!  write(*,*) 'g:',real(g)
!!!!!
!!!!!
!!!!!
!!!!!  block
!!!!!    use quantilmean
!!!!!    real(REAL64), dimension(:), allocatable :: xcdf,ycdf
!!!!!    real(REAL64) :: qmin,qmax,qmid,eps,x,smax,smin
!!!!!    integer :: n,i
!!!!!
!!!!!    n = size(data)
!!!!!    allocate(xcdf(n),ycdf(n))
!!!!!    call ecdf(data,xcdf,ycdf)
!!!!!
!!!!!    open(1,file='/tmp/testsr1mean_ecdf')
!!!!!    do i = 1, n
!!!!!       write(1,*) xcdf(i),ycdf(i)
!!!!!    end do
!!!!!    close(1)
!!!!!
!!!!!    call ecdf(abs(data-t),xcdf,ycdf)
!!!!!    open(1,file='/tmp/testsr1mean_mad')
!!!!!    do i = 1, n
!!!!!       write(1,*) xcdf(i),ycdf(i)
!!!!!    end do
!!!!!    close(1)
!!!!!
!!!!!!    block
!!!!!!      open(1,file='/tmp/testsr1mean_sbyf')
!!!!!!      do i = 2,n-1
!!!!!!         s = xcdf(i)
!!!!!!         write(1,*) s,funmean%val(t,s)
!!!!!!      end do
!!!!!!      close(1)
!!!!!!    end block
!!!!!
!!!!!
!!!!!!!$    block
!!!!!!!$      real(REAL64) :: s,K
!!!!!!!$      real(REAL64), dimension(n) :: Phi
!!!!!!!$      integer :: i,j,l
!!!!!!!$
!!!!!!!$      open(1,file='/tmp/testsr1mean_sk')
!!!!!!!$      do l = 1, 500
!!!!!!!$         s = 0.01*l
!!!!!!!$         Phi = (1 + erf((data+0.5)/s/sqrt(2.0)))/2
!!!!!!!$         K = 0
!!!!!!!$         do i = 1, n
!!!!!!!$            do j = 1, i-1
!!!!!!!$               K = K + sign(1.0_REAL64,(data(i)-data(j))*(Phi(i) - Phi(j)))
!!!!!!!$            end do
!!!!!!!$         end do
!!!!!!!$         K = K/(3*4*5)
!!!!!!!$         write(1,*) s,K
!!!!!!!$      end do
!!!!!!!$      close(1)
!!!!!!!$    end block
!!!!!
!!!!!
!!!!!!    write(*,*) real(xcdf)
!!!!!!    write(*,*) real(ycdf)
!!!!!
!!!!!
!!!!!    qmax = 0.8
!!!!!    if( qmax < ycdf(n-1) ) then
!!!!!       smax = quantile(qmax,xcdf,ycdf)
!!!!!       s = quantile(real(0.6745,REAL64),xcdf,ycdf)
!!!!!!       qmid = 0.5!0.6745
!!!!!    else
!!!!!!       qmid = 0.5!ycdf(n)*0.6745
!!!!!       smax = xcdf(n-1)
!!!!!       s = quantile(real(0.6745,REAL64)*ycdf(n-1),xcdf,ycdf)
!!!!!!       s = quantile(real(0.5,REAL64),xcdf,ycdf) /0.6745
!!!!!    end if
!!!!!    qmin = 0.4
!!!!!    smin = quantile(qmin,xcdf,ycdf)
!!!!!    qmid = 0.5
!!!!!    s = quantile(qmid,xcdf,ycdf) / 0.6745
!!!!!!    write(*,*) s
!!!!!!    s = quantile(real(0.6745,REAL64),xcdf,ycdf)
!!!!!    smin = 0.5*s
!!!!!    smax = 2*s
!!!!!    write(*,*) 'smin, smax',real(smin),real(smax)
!!!!!
!!!!!    g = [smin, s, smax ]
!!!!!    f = [t-s, t, t+s]
!!!!!!    s = s / 0.6745
!!!!!
!!!!!!    call ecdf(abs(data-t),xcdf,ycdf)
!!!!!!    eps = sqrt(log(2/0.05)/(2*n))
!!!!!!    open(1,file='/tmp/ecdf1')
!!!!!!    do i = 1, n
!!!!!!       write(1,*) xcdf(i),ycdf(i)
!!!!!!    end do
!!!!!!    close(1)
!!!!!!    open(1,file='/tmp/ecdf2')
!!!!!!    do i = -30, 30
!!!!!!       x = 0.1*i
!!!!!!       write(1,*) x,Phi(x,t,s),Phi(x,t,s)-eps,Phi(x,t,s)+eps
!!!!!!    end do
!!!!!!    close(1)
!!!!!
!!!!!  end block
!!!!!
!!!!!  call quantilefun(abs(data-t),qq,g,info)
!!!!!  write(*,*) 'g:',real(g)
!!!!!!  s = g(2) / 0.6745
!!!!!
!!!!!!!$  block
!!!!!!!$    use quantilmean
!!!!!!!$    real(REAL64), dimension(:), allocatable :: x,y
!!!!!!!$    allocate(x(size(data)),y(size(data)))
!!!!!!!$    call ecdf(abs(data-t),x,y)
!!!!!!!$    open(1,file='/tmp/ecdf')
!!!!!!!$    do i = 1, size(data)
!!!!!!!$       write(1,*) x(i),y(i)
!!!!!!!$    end do
!!!!!!!$    close(1)
!!!!!!!$  end block
!!!!!!!$  stop
!!!!!
!!!!!  block
!!!!!    use fmm
!!!!!    use quicksorter
!!!!!    real(REAL64) :: tol,ss,smax,smin
!!!!!    real(REAL64), dimension(size(data)) :: odata
!!!!!!    real(REAL64) :: smin,smax,tol
!!!!!    logical :: found
!!!!!
!!!!!    tol = 0.01*g(2)
!!!!!!    smin = g(1)
!!!!!!    smax = s
!!!!!    chi2t = t
!!!!!    chi2par = 2.5 ! inside the right interval of Tukey
!!!!!    neval = 0
!!!!!    !    simin = zeroin(smin,smax,chi2fun,tol)
!!!!!!    rmin = 0
!!!!!    rmin = sqrt(epsilon(rmin))
!!!!!    rmax = s
!!!!!    smin = zeroin(rmin,rmax,scalemin,tol)
!!!!!
!!!!!    ! Does the minimum falls inside the interval smin < .. < smax ?
!!!!!    found = rmin + tol < smin .and. smin < rmax - tol
!!!!!!    smin = simin
!!!!!
!!!!!    !    if( debug )
!!!!!    write(*,'(a,g20.6,a,g0.2,l2,i5)') &
!!!!!         'Minimal scale limit:',smin," +- ",tol,found,neval
!!!!!
!!!!!
!!!!!    tol = 0.1*g(2)
!!!!!    tol = 0.1*s
!!!!!!    smin = g(1)
!!!!!    smax = g(3)
!!!!!!    odata = data
!!!!!!    call qsort(odata)
!!!!!!    smax = odata(size(odata)-1)
!!!!!!    smax = maxval(abs(data-t))
!!!!!!    smax = 5
!!!!!    tol = (smax - smin) / 1000
!!!!!!    smin =
!!!!!!    smax = 5
!!!!!    chi2par = stukey
!!!!!    chi2par = spshuber
!!!!!    chi2t = t
!!!!!    write(*,*) '*',real(smin),real(smax),real(tol)
!!!!!!    s = fmin(10*smin,smax,chi2min,tol)
!!!!!    s = zeroin(smin,smax,chi2fun,tol)
!!!!!!    s = g(2)
!!!!!
!!!!!    ! Does the minimum falls inside the interval smin < .. < smax ?
!!!!!    found = smin + tol < s .and. s < smax - tol
!!!!!
!!!!!    !    if( debug )
!!!!!    write(*,'(a,g20.6,a,g0.2,l2,i5)') &
!!!!!         'Scale initial estimation:',s," +- ",tol,found
!!!!!
!!!!!!    ss = zeroin(g(1),g(3),chi2fun,tol)
!!!!!!    write(*,*) 'ss=',ss
!!!!!
!!!!!
!!!!!  end block
!!!!!
!!!!!
!!!!!  block
!!!!!    real(REAL64) :: s,smax,ds
!!!!!    ds = maxval(abs(data-t)) / 1000
!!!!!    chi2t = t
!!!!!
!!!!!    open(1,file='/tmp/chi2min')
!!!!!    do i = 1, 1000
!!!!!       s = i*ds
!!!!!       write(1,*) s,chi2min(s),chi2fun1(s)
!!!!!    end do
!!!!!    close(1)
!!!!!
!!!!!  end block
!!!!!
!!!!!
!!!!!  block
!!!!!    integer :: iflag
!!!!!    logical :: reli
!!!!!    call qmean(data,t,dt,s,reliable=reli,flag=iflag)
!!!!!    if( .not. reli ) stop 'not reli'
!!!!!    write(*,*) 'qmean init:',t,s,dt
!!!!!    f(1) = t - s
!!!!!    f(3) = t + s
!!!!!    g(1) = s/2
!!!!!    g(3) = 2*s
!!!!!  end block
!!!!!
!!!!!
!!!!!!!$  ! bi-search
!!!!!!!$  block
!!!!!!!$    real(REAL64) :: smin, smax, fmin, fmax, fzero
!!!!!!!$    integer :: iter
!!!!!!!$    smin = g(1)
!!!!!!!$    fmin = sum(pshuber((data-t)/smin)**2)/size(data) - spshuber
!!!!!!!$!    fmin = sum(itukey((data-t)/smin))/size(data) - stukey
!!!!!!!$    smax = 3*g(3)
!!!!!!!$    fmax = sum(pshuber((data-t)/smax)**2)/size(data) - spshuber
!!!!!!!$!    fmax = sum(itukey((data-t)/smax))/size(data) - stukey
!!!!!!!$    write(*,*) smin, smax, fmin, fmax
!!!!!!!$    do iter = 1, 10
!!!!!!!$       s = (smin + smax) / 2
!!!!!!!$       fzero = sum(pshuber((data-t)/s)**2)/size(data) - spshuber
!!!!!!!$!       fzero = sum(itukey((data-t)/s))/size(data) - stukey
!!!!!!!$       write(*,*) iter,real(s),real(fzero)
!!!!!!!$       if( abs(fzero) < epsilon(fzero) ) exit
!!!!!!!$       if( (smax - smin) < 0.01*g(2) ) exit
!!!!!!!$       if( fzero*fmax < 0 ) then
!!!!!!!$          smin = s
!!!!!!!$          fmin = fzero
!!!!!!!$       else
!!!!!!!$          smax = s
!!!!!!!$          fmax = fzero
!!!!!!!$       end if
!!!!!!!$    end do
!!!!!!!$    write(*,*) 'bi-section search for s=',real(s)
!!!!!!!$!    stop
!!!!!!!$  end block
!!!!!
!!!!!
!!!!!  block
!!!!!    use rfun
!!!!!    use quantilmean
!!!!!    use kstestscale
!!!!!    use fmm, only:zeroin
!!!!!    real(REAL64), dimension(:), allocatable :: xcdf,ycdf,xwins,res,fgrid
!!!!!    real(REAL64), dimension(2,2) :: hess
!!!!!    real(REAL64), dimension(2) :: jac
!!!!!    real(REAL64), dimension(10) :: opt
!!!!!    real(REAL64) :: ss,smin,f,s0,s1,ds,r,frac,h11max,fmin,frac1,tol,tt,sconv,h22,&
!!!!!         dmax,d,chi,ksmin,dmin,chi2min,chi2s,finfo,sc,hormin,klsum,kld
!!!!!    integer :: n,m,i,j,imax,nopt,iopt
!!!!!    logical :: reli
!!!!!
!!!!!    n = size(data)
!!!!!    allocate(xcdf(n),ycdf(n),xwins(n),res(n))
!!!!!    call ecdf(abs(data-t),xcdf,ycdf)
!!!!!    write(*,*) '|data-t|:', xcdf
!!!!!!    call ecdf(abs(tukey(data-t)),xcdf1,ycdf1)
!!!!!!    open(1,file='/tmp/testsr1mean_ecdf')
!!!!!!    do i = 1, size(xcdf)
!!!!!!       write(1,*) xcdf(i),ycdf(i),xcdf1(i),ycdf1(i)
!!!!!!    end do
!!!!!!    close(1)
!!!!!
!!!!!    ! winsorize
!!!!!    do i = 1, n
!!!!!       if( data(i) > 2 ) then
!!!!!          xwins(i) = 2
!!!!!       else if( data(i) < -2 ) then
!!!!!          xwins(i) = -2
!!!!!       else
!!!!!          xwins(i) = data(i)
!!!!!       end if
!!!!!    end do
!!!!!
!!!!!    open(1,file='/tmp/testsr1mean_sbyf')
!!!!!    fmin = huge(fmin)
!!!!!    dmin = huge(dmin)
!!!!!    chi2min = huge(chi2min)
!!!!!    smin = s
!!!!!    chi2t = t
!!!!!    h11max = 0.6
!!!!!    s0 = xcdf(1)
!!!!!    s1 = xcdf(n)
!!!!!    ds = (s1 - s0) / (100*n)
!!!!!    ds = s1 / (10000)
!!!!!    do i = 2,10000
!!!!!!       ss = xcdf(i)
!!!!!       ss = i*ds
!!!!!       f = funmean%val(t,ss)
!!!!!       if( .not. reli ) then
!!!!!          finfo = 0
!!!!!!       else
!!!!!!          finfo = 1 / finfo**2
!!!!!       end if
!!!!!       call funmean%hess(t,ss,hess)
!!!!!       frac = ss**2*hess(1,1)
!!!!!!       frac = funmean%nfrac(t,ss)
!!!!!       frac1 = count(abs((data-t)/ss) < 1) / real(n)
!!!!!       frac1 = ss**2*hess(2,2)
!!!!!       res = (data - t) /ss
!!!!!       h22 = (sum(dpshampel(res)*res**2) + 2*sum(pshampel(res)*res))/n   - ss**3/(ss**2 + 1)**1.5
!!!!!       frac1 = sum((ihuber((xwins-t)*ss)))/n + log(1/ss)
!!!!!       frac1 = sum((itukey((data-t)/ss)))/n + asinh(ss/2)
!!!!!       frac1 = sum((ipshampel((data-t)/ss)))/n + asinh(ss/2)
!!!!!       if( robfun == 'huber' ) then
!!!!!          frac1 = sum((ihuber((data-t)/ss)))/n
!!!!!          finfo = hubermean%stderr(t,ss,reli)
!!!!!       else if( robfun == 'hampel' ) then
!!!!!          frac1 = sum((ipshampel((data-t)/ss)))/n
!!!!!          finfo = bathmean%stderr(t,ss,reli)
!!!!!       else
!!!!!          frac1 = sum((itukey((data-t)/ss)))/n
!!!!!          finfo = tukeymean%stderr(t,ss,reli)
!!!!!       end if
!!!!!       finfo = sqrt(1.0*n)*finfo
!!!!!!       if( ss < exp(1.0/3.0) ) then
!!!!!!          frac1 = frac1 + 0.1226*ss**3
!!!!!!       else
!!!!!!          frac1 = frac1 + log(ss)
!!!!!!       end if
!!!!!       !       frac1 = sum((((data-t)/ss)**2/2))/n + log(ss)
!!!!!       call funmean%grad(t,ss,jac)
!!!!!
!!!!!       dmax = 0
!!!!!       imax = 0
!!!!!       sc = 0
!!!!!       klsum = 0
!!!!!       kld = 0
!!!!!       do j = 2, n
!!!!!          chi = erf(xcdf(j)/ss/sqrt(2.0))
!!!!!          kld = kld + exp(-ipshampel((data(j)-t)/ss))/ss*(xcdf(j) - xcdf(j-1))
!!!!!          kld = exp(-ipshampel((data(j)-t)/ss))/ss*(xcdf(j) - xcdf(j-1))
!!!!!          kld = exp(-ipshampel((data(j)-t)/ss))/ss
!!!!!          klsum = klsum + kld*log(kld)
!!!!!!          klsum = klsum + kld*log(kld/chi)
!!!!!          d = abs(ycdf(j) - chi)
!!!!!          if( d > dmax ) then
!!!!!             dmax = d
!!!!!             imax = j
!!!!!          end if
!!!!!          sc = sc + d!**2
!!!!!       end do
!!!!!!       write(*,*) imax,dmax
!!!!!!       imin = imax
!!!!!       if( dmax < dmin ) then
!!!!!          dmin = dmax
!!!!!          ksmin = ss
!!!!!       end if
!!!!!
!!!!!       if( g(2) < chi2min ) then
!!!!!          chi2min = g(2)
!!!!!          chi2s = ss
!!!!!       end if
!!!!!
!!!!!       write(1,*) ss,f,frac,h22,frac1,dmax,finfo,klsum,jac(2),jac(2)*ss**2, &
!!!!!            lepage(data,t,ss)
!!!!!       if( frac > h11max .and. frac < 0.8) then
!!!!!          if ( f < fmin ) then
!!!!!             smin = ss
!!!!!             fmin = f
!!!!!             h11max = frac
!!!!!          end if
!!!!!!          if( f < fmin ) then
!!!!!!             smin = ss
!!!!!!             fmin = f
!!!!!!             h11max = frac
!!!!!!          end if
!!!!!       end if
!!!!!    end do
!!!!!    close(1)
!!!!!    s = smin
!!!!!    write(*,*) 'bymin (graph):',s
!!!!!
!!!!!    m = n - int((1-0.5) * (n+1)) ! q <= 68%
!!!!!    ss = xcdf(m) / 0.6745
!!!!!    call funmean%hess(t,ss,hess)
!!!!!    ds = ss/sqrt(2*hess(1,1)*n) / 10
!!!!!    ds = funmean%stderr(t,ss,reli)*sqrt(1.0*n) / 5
!!!!!    ds = 1e-1
!!!!!    s0 = xcdf(1)
!!!!!    s1 = xcdf(n)
!!!!!    m = int((s1 - s0)/ds)
!!!!!    write(*,*) ss,hess(1,1), hess(2,2), ds,m,funmean%stderr(t,ss,reli)*sqrt(1.0*n)
!!!!!    allocate(fgrid(m))
!!!!!    do i = 1, m
!!!!!       ss = s0 + ds*i
!!!!!       fgrid(i) = funmean%val(t,ss)
!!!!!    end do
!!!!!    nopt = 0
!!!!!    do i = 2, m-1
!!!!!       if( fgrid(i-1) > fgrid(i) .and. fgrid(i) < fgrid(i+1) ) then
!!!!!          nopt = nopt + 1
!!!!!          opt(nopt) = s0 + ds*i
!!!!!       end if
!!!!!    end do
!!!!!    dmax = huge(dmax)
!!!!!    do i = 1, nopt
!!!!!       call funmean%hess(t,opt(i),hess)
!!!!!       r = opt(i)**2*hess(1,1)
!!!!!       f = funmean%val(t,opt(i))
!!!!!       d = abs(r - 0.68*0.8)
!!!!!       write(*,*) i,opt(i),r,d,hess(2,2),f,opt(i)/sqrt(hess(2,2)*n),hess(1,1)
!!!!!       if( d < dmax .and. hess(2,2) > 0 ) then
!!!!!          dmax = d
!!!!!          iopt = i
!!!!!       end if
!!!!!    end do
!!!!!    if( iopt > 0 ) then
!!!!!       s = opt(iopt)
!!!!!    else
!!!!!       s = -1
!!!!!    end if
!!!!!    write(*,*) 'bymin (brute force):',s
!!!!!
!!!!!
!!!!!    do i = n-1, 2, -1
!!!!!       if( ycdf(i) < 0.68 ) then
!!!!!          s = xcdf(i)
!!!!!          write(*,*) 'by chi:',s
!!!!!          exit
!!!!!       end if
!!!!!    end do
!!!!!
!!!!!
!!!!!    s0 = s
!!!!!    s1 = maxval(abs(data - t))
!!!!!    tol = 0.01*s1
!!!!!    s0 = xcdf(2)
!!!!!    s1 = xcdf(n-1)
!!!!!    write(*,*) s0,s1
!!!!!    s = zeroin(s0,s1,fsuff,tol)
!!!!!    write(*,*) 'byzeroin:',s
!!!!!
!!!!!    call byS(data,s)
!!!!!    write(*,*) 'byS:',s
!!!!!
!!!!!!    call byQ(data,s)
!!!!!!    s = byQ(data)
!!!!!    write(*,*) 'byQ:',s
!!!!!
!!!!!!    stop
!!!!!
!!!!!!    s = ksmin
!!!!!!    write(*,*) 'by K-S test:',s
!!!!!
!!!!!!    s = kscale(data-t)
!!!!!!    write(*,*) 'by K-S test optim:',s
!!!!!
!!!!!!    s = gscale(funmean,t)
!!!!!!    write(*,*) 'by gscale:',s
!!!!!
!!!!!!    s = fscale(funmean,t,xcdf(2),xcdf(n))
!!!!!!    write(*,*) 'by fscale:',s
!!!!!
!!!!!!    s = goldscale(funmean,t)
!!!!!!    write(*,*) 'by goldscale:',s
!!!!!
!!!!!!    s = chi2s
!!!!!!    write(*,*) 'by chi2 test:',s
!!!!!!    s = g(2)
!!!!!
!!!!!    !    s = kscale(data-t)
!!!!!!    s = sum(xcdf(2:4)) /  3
!!!!!!    s = (0.25*xcdf(2) + 0.5*xcdf(3) + 0.25*xcdf(4))
!!!!!!    s = 0.92831
!!!!!
!!!!!!    call fline(funmean,t,xcdf(2),xcdf(5),s)
!!!!!!    write(*,*) 'by fline:',s
!!!!!
!!!!!
!!!!!  end block
!!!!!
!!!!!  block
!!!!!
!!!!!    use biscales
!!!!!    use quantilmean
!!!!!    real(REAL64), dimension(:), allocatable :: xcdf,ycdf
!!!!!    real(REAL64), parameter :: stol = 0.1
!!!!!    real(REAL64) :: s
!!!!!    integer :: n,flag
!!!!!
!!!!!    n = size(data)
!!!!!    allocate(xcdf(n),ycdf(n))
!!!!!!    call ecdf(abs(data-t),xcdf,ycdf)
!!!!!
!!!!!!    call biscale(funmean,t,xcdf(2),xcdf(5),stol,s,flag)
!!!!!
!!!!!  end block
!!!!!
!!!!!!!$  block
!!!!!!!$
!!!!!!!$    use kstestscale
!!!!!!!$    use quantilmean
!!!!!!!$    use trust_global_scale
!!!!!!!$
!!!!!!!$    real(REAL64), dimension(:), allocatable :: xcdf,ycdf
!!!!!!!$    real(REAL64) :: ds,r1,r2,smin,smax,beta,smad,serr,sig,tl,tmin
!!!!!!!$    integer :: i,j,n,m
!!!!!!!$    logical :: converge
!!!!!!!$
!!!!!!!$    n = size(data)
!!!!!!!$    allocate(xcdf(n),ycdf(n))
!!!!!!!$    call ecdf(abs(data-t),xcdf,ycdf)
!!!!!!!$
!!!!!!!$    m = n - int((1-0.5) * (n)) ! q <= 68%
!!!!!!!$    smad = xcdf(m) / 0.6745
!!!!!!!$    write(*,*) 's=',smad,m,xcdf
!!!!!!!$
!!!!!!!$    serr = smad / sqrt(real(n))
!!!!!!!$    write(*,*) 'serr:',serr
!!!!!!!$!    k =
!!!!!!!$
!!!!!!!$    if( n > 13 ) then
!!!!!!!$       m = int(0.9 * (n+1))
!!!!!!!$       smax = xcdf(m)
!!!!!!!$       m = n + 1 - m
!!!!!!!$       smin = xcdf(m)
!!!!!!!$    else
!!!!!!!$       smin = xcdf(2)
!!!!!!!$       smax = xcdf(n)
!!!!!!!$    end if
!!!!!!!$    m = int(0.75 * (n+1))
!!!!!!!$!    smax = xcdf(m)
!!!!!!!$    m = n + 1 - m
!!!!!!!$!    smin = xcdf(m)
!!!!!!!$
!!!!!!!$
!!!!!!!$    call ecdf(data,xcdf,ycdf)
!!!!!!!$
!!!!!!!$!    write(*,*) 'xcdf:',xcdf
!!!!!!!$!    write(*,*) 'ycdf:',ycdf
!!!!!!!$
!!!!!!!$    r1 = sum(abs(ycdf - 0.5)) - n / 2.0
!!!!!!!$
!!!!!!!$    open(1,file='/tmp/testsr1mean_rank')
!!!!!!!$    ds = 0.05
!!!!!!!$    do i = 1,200
!!!!!!!$       s = i*ds
!!!!!!!$
!!!!!!!$!       r1 = 0
!!!!!!!$       r2 = 0
!!!!!!!$       do j = 1, n
!!!!!!!$!          r1 = r1 + abs(ycdf(j) - 0.5) - 0.5
!!!!!!!$          r2 = r2 + abs(Phi(xcdf(j),t,s) - 0.5) - 0.5
!!!!!!!$!          r1 = r1 + (ycdf(j) - 0.5)**2 - 1.0/12.0
!!!!!!!$!          r2 = r2 + (Phi(xcdf(j),t,s) - 0.5)**2 - 1.0/12.0
!!!!!!!$!          r1 = r1 + inorm(ycdf(j))**2 - 1
!!!!!!!$!          r2 = r2 + inorm(Phi(xcdf(j),t,s))**2 - 1
!!!!!!!$!          r2 = r2 + ((xcdf(j) - t)/s)**2 - 1
!!!!!!!$       end do
!!!!!!!$       r2 = sum(abs(Phi(xcdf,t,s) - 0.5)) - n / 2.0
!!!!!!!$
!!!!!!!$       write(1,*) s,(r1-r2)/sqrt(2.0*n),r1,r2
!!!!!!!$
!!!!!!!$    end do
!!!!!!!$    close(1)
!!!!!!!$
!!!!!!!$    beta = 0.8
!!!!!!!$    s = smad
!!!!!!!$    sig = s / sqrt(real(n))
!!!!!!!$
!!!!!!!$    call trgscale(funmean,t,smin,smax,sig,smad,nopt,sopt,converge)
!!!!!!!$!    call trgscale(funmean,t,xcdf,ycdf,smin,smax,s,converge)
!!!!!!!$    write(*,*) 'by trgscale:',sopt(1:nopt),sig,converge
!!!!!!!$    !    stop
!!!!!!!$
!!!!!!!$    if( nopt > 1 ) then
!!!!!!!$       tmin = huge(tmin)
!!!!!!!$       do i = 1, nopt
!!!!!!!$          tl = lepage(data,t,sopt(i))
!!!!!!!$          write(*,*) 'lepage:',sopt(i),tl
!!!!!!!$          if( tl < tmin ) then
!!!!!!!$             s = sopt(i)
!!!!!!!$             tmin = tl
!!!!!!!$          end if
!!!!!!!$       end do
!!!!!!!$    else
!!!!!!!$       s = sopt(1)
!!!!!!!$    end if
!!!!!!!$
!!!!!!!$    s = smad
!!!!!!!$
!!!!!!!$
!!!!!!!$
!!!!!!!$  end block
!!!!!
!!!!!
!!!!!!  block
!!!!!!    use likescale
!!!!!!    real(REAL64) :: smin
!!!!!
!!!!!!    call iscale(abs(data - t),smin,verbose=.true.)
!!!!!!    write(*,*) 'by iscale:',smin
!!!!!!    s = smin
!!!!!!  end block
!!!!!
!!!!!
!!!!!  block
!!!!!
!!!!!    use quicksorter
!!!!!    use kstestscale
!!!!!
!!!!!    real(REAL64), dimension(:), allocatable :: x,d
!!!!!    real(REAL64) :: x1,x2,y1,y2,med,mad,sig,h
!!!!!    integer :: m,n,c1,m1,m2
!!!!!
!!!!!    n = size(data)
!!!!!    allocate(x,source=data)
!!!!!    call qsort(x)
!!!!!    if( mod(n,2) == 0 ) then
!!!!!       m = n / 2
!!!!!       med = (x(m) + x(m+1))/2
!!!!!    else
!!!!!       m = (n + 1)/2
!!!!!       med = x(m)
!!!!!    end if
!!!!!
!!!!!!    m2 = max(int(0.75 * (n+1)),n)
!!!!!!    m1 = n + 1 - m2
!!!!!!    med = (2*med + x(m1) + x(m2))/4
!!!!!
!!!!!    c1 = nint(n/2 - 1.96*sqrt(1.0*n)/2 -0.5)
!!!!!    write(*,*) n/2 - 1.96*sqrt(1.0*n)/2 -0.5,c1
!!!!!
!!!!!    allocate(d(n))
!!!!!    d = abs(x - med)
!!!!!    call qsort(d)
!!!!!    if( mod(n,2) == 0 ) then
!!!!!       m = n / 2
!!!!!       mad = (d(m) + d(m+1))/2
!!!!!    else
!!!!!       m = (n + 1)/2
!!!!!       mad = d(m)
!!!!!    end if
!!!!!    sig = mad / 0.6745
!!!!!    write(*,*) 'sig,mad,med:',sig,mad,med
!!!!!!    stop
!!!!!
!!!!!    h = real(1.0,REAL64) / real(n + 1,REAL64)
!!!!!    m2 = min(int(0.95 * (n+1)),n)
!!!!!    m1 = max(2,n + 1 - m2)
!!!!!    y1 = d(m1)
!!!!!    y2 = d(m2)
!!!!!    call fminsig(t,y1,y2,sig/10,sig)
!!!!!    write(*,*) 'sig,m1,m2:',t,sig,m1,m2,y1,y2
!!!!!    s = sig
!!!!!
!!!!!!    call byQ(data,sig)
!!!!!!    s = byQ(data)
!!!!!    s = 0.81*sig
!!!!!
!!!!!    !    stop
!!!!!
!!!!!    open(3,file='/tmp/tsr1mean_qlim')
!!!!!    if( n > 13 ) then
!!!!!       x1 = x(m-c1)
!!!!!       x2 = x(m+c1)
!!!!!       y1 = d(m-c1)
!!!!!       y2 = d(m+c1)
!!!!!       x1 = med - 2*sig/sqrt(1.0*n)
!!!!!       x2 = med + 2*sig/sqrt(1.0*n)
!!!!!       y1 = sig - 2*sig*sqrt(0.5/n)
!!!!!       y2 = sig + 2*sig*sqrt(0.5/n)
!!!!!    else
!!!!!       x1 = 0.81*x(m-c1)
!!!!!       x2 = 0.81*x(m+c1)
!!!!!!       y1 = d(2)
!!!!!!       y2 = d(n)
!!!!!    end if
!!!!!    write(3,*) x1,y1
!!!!!    write(3,*) x2,y1
!!!!!    write(3,*) x2,y2
!!!!!    write(3,*) x1,y2
!!!!!    write(3,*) x1,y1
!!!!!
!!!!!!    write(3,*) f(1),g(1)
!!!!!!    write(3,*) f(3),g(1)
!!!!!!    write(3,*) f(3),g(3)
!!!!!!    write(3,*) f(1),g(3)
!!!!!!    write(3,*) f(1),g(1)
!!!!!    close(3)
!!!!!  end block
!!!!!
!!!!!!  t = -3
!!!!!!  s = 1
!!!!!!  s = 1
!!!!!!  x = [ 0.0,  1.1]
!!!!!!  x = [ -0.0, 4.1 ]
!!!!!  delta = 1
!!!!!  maxdev = 1e-7
!!!!!  maxfev = 3*precision(t)
!!!!!  verbose = .true.
!!!!!  quasi = .true.
!!!!!  open(2,file='/tmp/tsr1mean_step')
!!!!!
!!!!!  nopt = 1
!!!!!  sopt(1) = s
!!!!!  !  s = 0.7
!!!!!
!!!!!  ! ridge decision
!!!!!  block
!!!!!
!!!!!    use trustsr1mean
!!!!!
!!!!!    logical :: reli
!!!!!    real(REAL64) :: h22, tol, delta, d
!!!!!    real(REAL64), dimension(2,2) :: hess,Q
!!!!!    real(REAL64), dimension(2) :: lam
!!!!!    integer :: ieig,k,l
!!!!!
!!!!!    write(*,*) t,s
!!!!!    call funmean%hess(t,0.75*s,hess)
!!!!!    call q2eig(hess,lam,Q)
!!!!!    write(*,*) 'lam:',lam
!!!!!
!!!!!!    h22 = hess(1,1)*hess(2,2) - hess(1,2)**2
!!!!!!    if( lam(1) > 0 .and. lam(2) > 0 ) then
!!!!!    nopt = 1
!!!!!    topt(1) = t
!!!!!    sopt(1) = 0.75*s
!!!!!    if( lam(1)*lam(2) < 0 ) then
!!!!!!    if( .false. ) then
!!!!!       ieig = 0
!!!!!       do i = 1, 2
!!!!!          !       if( lam(i) < 0 ) pvec = Q(:,i)
!!!!!          if( lam(i) < 0 ) ieig = i
!!!!!       end do
!!!!!       if( ieig == 0 ) stop 'Something wrong.....'
!!!!!!       ieig = 3 - ieig  ! perpedicular vector
!!!!!       tol = 1e-3 * s
!!!!!       delta = s / 4
!!!!!       write(*,*) 'Ridge start with ',t,s,' dir:',real(Q(:,ieig))
!!!!!!       call ridgemax(funmean,t,s,delta,tol,Q(:,ieig),reli)
!!!!!       call smnewton(funmean,Q(:,ieig),t,s,tol,delta,reli)
!!!!!       write(*,*) 'Ridge max at ', t,s,reli
!!!!!       if( reli ) then !stop 'Something wrong...'
!!!!!          call funmean%hess(t,s,hess)
!!!!!          call q2eig(hess,lam,Q)
!!!!!          ieig = 0
!!!!!          do i = 1, 2
!!!!!             if( lam(i) < 0 ) ieig = i
!!!!!          end do
!!!!!          if( ieig == 0 ) stop 'Something wrong.....'
!!!!!          write(*,*) 'Ridge dir ',real(Q(:,ieig))
!!!!!          nopt = 2
!!!!!          d = 0.5*delta
!!!!!          topt(1) = t - d*Q(1,ieig)
!!!!!          topt(2) = t + d*Q(1,ieig)
!!!!!          sopt(1) = s - d*Q(2,ieig)
!!!!!          sopt(2) = s + d*Q(2,ieig)
!!!!!
!!!!!          do k = 1,2*0
!!!!!             d = (-1)**k*delta
!!!!!             do l = 1, 10
!!!!!                topt(k) = t + l*d*Q(1,ieig)
!!!!!                sopt(k) = s + l*d*Q(2,ieig)
!!!!!                write(*,*) k,l,topt(k),sopt(k)
!!!!!                call funmean%hess(topt(k),sopt(k),hess)
!!!!!                if( hess(1,1) > 0 .and. hess(1,1)*hess(2,2) - hess(1,2)**2 > 0 ) exit
!!!!!             end do
!!!!!          end do
!!!!!
!!!!!          write(*,*) topt(1:nopt)
!!!!!          write(*,*) sopt(1:nopt)
!!!!!!          stop
!!!!!       else
!!!!!          write(*,*) 'ridgemax fail...'
!!!!!       end if
!!!!!    end if
!!!!!
!!!!!
!!!!!!       stop 'Ridge in doubts.'
!!!!!!    end if
!!!!!  end block
!!!!!!  stop
!!!!!
!!!!!  curvature = 0
!!!!!  iopt = 0
!!!!!  do i = 1, nopt
!!!!!
!!!!!     !     s = 0.75*sopt(i) ! the correction on s rather than stdsig
!!!!!     t = topt(i)
!!!!!     s = sopt(i) ! the correction on s rather than stdsig
!!!!!     write(*,*) 'start point:',t,s
!!!!!     !  call sr1rob(funmean,x,cov,delta,maxdev,maxfev,quasi,verbose,info,msg)
!!!!!     !  call tsr1mean(data,t,dt,s,info,verbose)
!!!!!     block
!!!!!
!!!!!       use kstestscale
!!!!!
!!!!!       logical :: reli
!!!!!       logical, parameter :: quasi = .true.
!!!!!       integer, parameter :: maxfev = 3*precision(t)
!!!!!       real(REAL64), dimension(2,2) :: hess
!!!!!       real(REAL64) :: tol,delta,K
!!!!!
!!!!! !      s = sopt(i)
!!!!!       tol = 0.01 * s / sqrt(real(size(data)))
!!!!!       delta = 0.3*s
!!!!!
!!!!!       !    call tsr1meansqp(funmean,t,s,hess,delta,tol,maxfev,quasi,verbose,info,msg)
!!!!!       delta = s / sqrt(real(size(data)))
!!!!!       delta = 0.9*s
!!!!!       !    t = -3
!!!!!       !    s = 0.1
!!!!!!           t = 0
!!!!!!           s = 0.01
!!!!!       call trfmeansqp(funmean,t,s,hess,delta,tol,maxfev,verbose,info)
!!!!!       write(*,*) "Lepage's test:",lepage(data,t,s)
!!!!!       call funmean%hess(t,s,hess)
!!!!!       K = hess(1,1)*hess(2,2)-hess(1,2)**2
!!!!!       write(*,*) 'Curvature:',K,hess
!!!!!       if( K > curvature ) then
!!!!!          curvature = K
!!!!!          topt(i) = t
!!!!!          sopt(i) = s
!!!!!          iopt = i
!!!!!          dt = funmean%stderr(t,s,reli)
!!!!!       end if
!!!!!
!!!!!     end block
!!!!!  end do
!!!!!
!!!!!  if( iopt == 0 ) stop 'Solution not found.'
!!!!!  t = topt(iopt)
!!!!!  s = sopt(iopt)
!!!!!
!!!!!  write(*,*) info,trim(msg)
!!!!!  write(*,*) 'result:',real(t),real(s)
!!!!!  write(*,*) 'stderr=',real(dt),' stddev=',real(dt*sqrt(real(size(data))))
!!!!!
!!!!!  open(3,file='/tmp/tsr1mean_result')
!!!!!  write(3,*) t,s
!!!!!  close(3)
!!!!!  close(2)
!!!!!
!!!!!!!$  block
!!!!!!!$    real(REAL64) :: df
!!!!!!!$    integer :: l
!!!!!!!$    df = 0
!!!!!!!$    do l = 1, size(data)
!!!!!!!$       df = df + max(dtukey((data(l)-t)/s),real(0,REAL64))
!!!!!!!$    end do
!!!!!!!$    df = df / size(data)
!!!!!!!$  end block
!!!!!!  df = sum(dtukey((data-t)/s), dtukey((data-t)/s) > 0) / size(data)
!!!!!  df = funmean%nfrac(t,s)
!!!!!  write(*,*) 'df=',real(100*df),'%'
!!!!!
!!!!!  block
!!!!!
!!!!!    integer :: i,j,n
!!!!!    real :: sig
!!!!!    integer, dimension(3) :: cts
!!!!!    n = size(data)
!!!!!    sig = sqrt(1.0*n)*dt
!!!!!    cts = 0
!!!!!    do i = 1, n
!!!!!       do j = 1, 3
!!!!!          if( abs(data(i) - t) < j*sig ) cts(j) = cts(j) + 1
!!!!!       end do
!!!!!    end do
!!!!!    write(*,'(3(i5,a,f6.1,a))') (j,'-sigma:',100*real(cts(j))/n,'%,',j=1,3)
!!!!!  end block
!!!!!
!!!!!  block
!!!!!    use kstestscale
!!!!!
!!!!!    real(REAL64) :: t,s
!!!!!    integer :: n
!!!!!
!!!!!    n = size(data)
!!!!!    call qmean(data,t,stdsig=s)
!!!!!    write(*,*) 'Quatil.mean:',real(t),' dt=',real(s/sqrt(real(n))),' sig=',real(s)
!!!!!!    call byQ(data,s)
!!!!!!    s = byQ(data)
!!!!!    write(*,*) 'Quatil.mean:',real(t),' dt=',real(s/sqrt(real(n))),' sig=',real(s)
!!!!!  end block
!!!!!
!!!!!  block
!!!!!    real(REAL64) :: t,s
!!!!!    integer :: n
!!!!!
!!!!!    n = size(data)
!!!!!    t = sum(data) / n
!!!!!    s = sqrt(sum((data-t)**2)/(n-1))
!!!!!    write(*,*) 'Arith.mean:',real(t),' dt=',real(s/sqrt(real(n))),' sig=',real(s)
!!!!!  end block
!!!!!
!!!!!  block
!!!!!    use rfun
!!!!!    real(REAL64) :: stderr, f2,df
!!!!!    real(REAL64), dimension(2,2) :: hess
!!!!!    integer :: n
!!!!!    logical :: reli
!!!!!
!!!!!    n = size(data)
!!!!!    call funmean%hess(t,s,hess)
!!!!!    stderr = funmean%stderr(t,s,reli)
!!!!!    f2 = sum(tukey((data-t)/s)**2)/n
!!!!!    df = sum(dtukey((data-t)/s))/n
!!!!!    write(*,*) 'hess:',hess
!!!!!    write(*,*) 'stderr:',stderr
!!!!!    write(*,*) 'info:', df/f2, f2,df
!!!!!    write(*,*) 'est.:',sqrt(0.9/hess(1,1)),sqrt(0.9/hess(1,1)/n),sqrt(1.0/hess(1,1)),sqrt(1.0/hess(1,1)/n)
!!!!!    write(*,*) 'est.:',sqrt(hess(2,2)/(hess(1,1)*hess(2,2)-hess(1,2)**2)), &
!!!!!         sqrt(hess(2,2)/(hess(1,1)*hess(2,2)-hess(1,2)**2)/n)
!!!!!  end block
!!!!!
!!!!!
!!!!!  end subroutine oldway
!!!!!

  subroutine graph(fcn,x0,x1,y0,y1,filename)

!    class(tsr1mean_function), intent(in) :: fcn
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

  end subroutine graph

!!$  subroutine graphd(fcn,x0,x1,y0,y1,filename)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL64), intent(in) :: x0,x1,y0,y1
!!$    character(len=*), intent(in) :: filename
!!$    integer, parameter :: nx = 200
!!$    integer, parameter :: ny = 200
!!$    real(REAL64) :: x,dx,y,dy
!!$    integer :: i,j
!!$
!!$    dx = (x1 - x0) / (nx - 1)
!!$    dy = (y1 - y0) / (ny - 1)
!!$
!!$    open(1,file=filename)
!!$    do j = 1, ny
!!$       y = y0 + dy*j
!!$       do i = 0, nx
!!$          x = x0 + dx*i
!!$          write(1,*) x,y,fcn%val(x,y)
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graphd


!!$  subroutine graph_marginal(fcn,x0,x1,y0,y1,filename)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL32), intent(in) :: x0,x1,y0,y1
!!$    character(len=*), intent(in) :: filename
!!$    integer, parameter :: nx = 1000
!!$    integer, parameter :: ny = 1000
!!$    real(REAL64) :: x,dx,y,dy,ss
!!$    integer :: i,j
!!$
!!$    dx = (x1 - x0) / (nx - 1)
!!$    dy = (y1 - y0) / (ny - 1)
!!$
!!$    open(1,file=filename)
!!$    do i = 0, nx
!!$       x = x0 + dx*i
!!$       ss = 0
!!$       do j = 1, ny
!!$          y = y0 + dy*j
!!$          ss = ss + exp(-fcn%val(x,y))
!!$       end do
!!$       write(1,*) x,ss*dy
!!$    end do
!!$    close(1)
!!$
!!$    open(1,file=trim(filename)//'_nuisance')
!!$    do j = 1, ny
!!$       y = y0 + dy*j
!!$       ss = 0
!!$       do i = 0, nx
!!$          x = x0 + dx*i
!!$          ss = ss + exp(-fcn%val(x,y))
!!$       end do
!!$       write(1,*) y,ss*dx
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graph_marginal


  subroutine graph_hess(fcn,x0,x1,y0,y1,filename)

    use modQ2eig_REAL64

    class(robust_likelihood_REAL64), intent(in) :: fcn
!    class(tsr1mean_function), intent(in) :: fcn
    real(REAL64), intent(in) :: x0,x1,y0,y1
    character(len=*), intent(in) :: filename
    real(REAL64), dimension(2,2) :: h,Q
    real(REAL64), dimension(2) :: lam
    real(REAL64) :: x,dx,y,dy
    integer :: nx,ny,i,j,n

    nx = 200
    ny = 200
    dx = (x1 - x0) / (nx - 1)
    dy = (y1 - y0) / (ny - 1)

    open(1,file=filename)
    do j = 1, ny
       y = y0 + dy*j
       do i = 0, nx
          x = x0 + dx*i
          call fcn%hess(x,y,h)
          call q2eig_REAL64(h,lam,Q,n)
          h = y**2*h
          lam = y**2*lam
          write(1,*) x,y,h(1,1),h(1,2),h(2,2),(h(2,2)*h(1,1)-h(1,2)**2),lam,Q
       end do
       write(1,*)
    end do
    close(1)

!!$    open(1,file=trim(filename)//"_Q")
!!$    nx = 40
!!$    ny = 20
!!$    dx = 0.25
!!$    dy = 0.25
!!$    do j = 1, ny
!!$       y = y0 + dy*j
!!$       do i = 0, nx
!!$          x = x0 + dx*i
!!$          call fcn%hess(x,y,h)
!!$          call q2eig(h,lam,Q)
!!$          h = y**2*h
!!$          lam = y**2*lam
!!$          write(1,*) x,y,h(1,1),h(1,2),h(2,2),(h(2,2)*h(1,1)-h(1,2)**2),lam,Q
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)


  end subroutine graph_hess


  subroutine graph_rho(fcn,t0,t1,s0,s1,filename)

    use rfun
!    use rfun_hampel_REAL64

    class(robust_likelihood_REAL64), intent(in) :: fcn
!    class(tsr1mean_function), intent(in) :: fcn
    real(REAL64), intent(in) :: t1,t0,s0,s1
    character(len=*), intent(in) :: filename
    integer, parameter :: n = 200
    real(REAL64), dimension(:), allocatable :: rho,w1,w2
    real(REAL64) :: t,s,ds,dt,chi2,df,df2,df3,ss
    integer :: i,j

    ds = (s1 - s0) / (n - 1)
    dt = (t1 - t0) / (n - 1)

    allocate(rho(size(data)),w1(size(data)),w2(size(data)))

    open(1,file=filename)
    do j = 1, n
       s = s0 + ds*j
       do i = 0, n
          t = t0 + dt*i
          chi2 = sum(itukey((data-t)/s)) / size(data)
          chi2 = fcn%func(t,s) - log(s)
          df = sum(dtukey((data-t)/s), dtukey((data-t)/s) > 0) / size(data)
!          df = sum(dtukey((data-t)/s)) / size(data)
          !          df = sum(dtukey((data-t)/s), dtukey((data-t)/s) > 0) / count(dtukey((data-t)/s) > 0)
          df = sum(dpshuber((data-t)/s), dpshuber((data-t)/s) > 0) / size(data)
!          df = sum(dpshampel((data-t)/s), dpshampel((data-t)/s) > 0) / size(data)
          df2 = sum(dpshampel((data-t)/s)) / size(data)
          df3 = sum(dtukey((data-t)/s)) / size(data)
          rho = ipshampel((data-t)/s)
          rho = itukey((data-t)/s)
          w1 = exp(-rho)/s/sqrt(2*3.14)
!          w1 = exp(-rho)/s
          w2 = exp(-data**2/2)/sqrt(2*3.14)
!          w1 = exp(-((data-t)/s)**2/2)/s/sqrt(2*3.14)
!          ss = sum(w2*log(w2/w1)) / n
          ss = sum(w1*log(w1),w1 > 0) / n
!          ss = sum(tukey((data-t)/s)**2)
!          e1 = sum(ihuber((data-t)/s)*exp(-0.5*ihuber((data-t)/s)))
!          e2 = sum(ihampel((data-t)/s)*exp(-0.5*ihampel((data-t)/s)))
!          e3 = sum(itukey((data-t)/s)*exp(-0.5*itukey((data-t)/s)))
          write(1,*) t,s,chi2,df,df2,df3,ss!,e1,e2,e3
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_rho

  subroutine graph_jac(fcn,x0,x1,y0,y1,filename)

    class(robust_likelihood_REAL64), intent(in) :: fcn
!    class(tsr1mean_function), intent(in) :: fcn
    real(REAL64), intent(in) :: x0,x1,y0,y1
    character(len=*), intent(in) :: filename
    integer, parameter :: nx = 200
    integer, parameter :: ny = 200
    real(REAL64) :: x,dx,y,dy
    real(REAL64), dimension(2) :: g
    integer :: i,j

    dx = (x1 - x0) / (nx - 1)
    dy = (y1 - y0) / (ny - 1)

    open(1,file=filename)
    do j = 1, ny
       y = y0 + dy*j
       do i = 0, nx
          x = x0 + dx*i
          call fcn%grad(x,y,g)
          write(1,*) x,y,y*g,norm2(g),57.3*atan2(g(2),g(1))
       end do
       write(1,*)
    end do
    close(1)

  end subroutine graph_jac


  subroutine graph_ridge(fcn,tmin,tmax,smin,smax,filename)

    use gold_ridge

    class(robust_likelihood_REAL64), intent(in) :: fcn
    real(REAL64), intent(in) :: tmin,tmax,smin,smax
    character(len=*), intent(in) :: filename
    integer, parameter :: n = 200
    real(REAL64), parameter :: rtol = 0.001
    real(REAL64), dimension(2), parameter :: p = [0,1]

    type(gold_fun) :: goldfun
    real(REAL64), dimension(2) :: q
    real(REAL64) :: tol, dt, sridge, a, b, fmax, f, s, t
    integer :: i, j, iflag
    logical :: found

    tol = rtol*smax
    dt = (tmax - tmin) / (n - 1)
    open(1,file=filename)

    do i = 0, n
       t = tmin + i*dt
       q = [t,10*epsilon(s)]
       a = smin
       b = smax
       found = .false.
       fmax = -huge(fmax)
       do j = 1, 50
          goldfun = gold_fun(fcn,-1,q,p,0)
          call goldfinger_REAL64(goldfun,a,b,tol,s,iflag)
!          write(*,*) t,s,goldfun%f(s),iflag,a,b
          if( iflag == 1 ) then
             f = goldfun%f(s)
             if( f > fmax ) then
                fmax = f
                sridge = s
!                write(*,*) t,s,f
                found = .true.
             end if
!             exit
          else if( iflag == 2 ) then
             a = s + 2.01*tol
             goldfun = gold_fun(fcn,1,q,p,0)
             call goldfinger_REAL64(goldfun,a,b,tol,s,iflag)
             a = s
          else if( iflag == 3 ) then
             b = s - 2.01*tol
             goldfun = gold_fun(fcn,1,q,p,0)
             call goldfinger_REAL64(goldfun,a,b,tol,s,iflag)
!             write(*,*) a,b
             b = s
          else
             ! the case of non-redescendent function
             exit
          end if
       end do
       if( found ) write(1,*) t,sridge,0
    end do
    close(1)

  end subroutine graph_ridge

!!!!!  subroutine graph_smin(fcn,tmin,tmax,smax,filename)
!!!!!
!!!!!    use rfun
!!!!!    use zeroinmod
!!!!!
!!!!!    class(zeroin_data), intent(in out) :: fcn
!!!!!    real(REAL64), intent(in) :: tmin,tmax,smax
!!!!!    character(len=*), intent(in) :: filename
!!!!!    integer, parameter :: n = 200
!!!!!    real(REAL64) :: tol = 0.01
!!!!!    real(REAL64), dimension(2,2) :: hess
!!!!!    real(REAL64), dimension(:), allocatable :: r
!!!!!    real(REAL64) :: s,dt,t,smin,d,beta,ex
!!!!!    integer :: i,j,k,iflag,m
!!!!!    logical :: converge
!!!!!
!!!!!    dt = (tmax - tmin) / (n - 1)
!!!!!    smin = sqrt(epsilon(smin))
!!!!!    smin = dt
!!!!!
!!!!!!    write(*,*) (fcn%funmean%data - 0.1)/6
!!!!!!    write(*,*) minval(abs(fcn%funmean%data - 0.1)/6)
!!!!!    !    if( abs(smin) < 1e3*epsilon(smin) ) ! nutne CDF ?
!!!!!
!!!!!    m = size(fcn%funmean%data)
!!!!!    allocate(r(m))
!!!!!
!!!!!    open(1,file=filename)
!!!!!    do i = 0, n
!!!!!       t = tmin + dt*i
!!!!!       fcn%t = t
!!!!!       smin = minval(abs(fcn%funmean%data - t))/(0.9*6)
!!!!!       smin = huge(smin)
!!!!!       do j = 1, size(data)
!!!!!          smin = min(max(abs(data(j) - t),1e3*epsilon(smin)),smin)
!!!!!       end do
!!!!!       smin = smin / 6
!!!!!       s = smin
!!!!!       s = zeroin(fcn,dt,smax,tol,iflag)
!!!!!       if( iflag == 0 ) then
!!!!!          s = max(s,dt)
!!!!!       else
!!!!!          s = smin
!!!!!       end if
!!!!!       !      s = abs(t) / 5.41
!!!!!!       s = smin
!!!!!       tol = 0.1*s
!!!!!
!!!!!       call smnewton(fcn%funmean,real([0,1],REAL64),t,s,tol,0.2*s,converge)
!!!!!       !       call snewton(fcn%funmean,t,tol,s,converge)
!!!!!       r = (fcn%funmean%data - t) / 0.75
!!!!!       beta = sum(tukey(r)**2) / m
!!!!!       ex = sum(r**2) / m
!!!!!       write(1,*) t,s,0,converge,smin,beta,ex,sum(tukey(r)*r)/m,sum(dtukey(r))/m
!!!!!    end do
!!!!!    close(1)
!!!!!
!!!!!!!$    open(1,file='/tmp/rk4')
!!!!!!!$    t = tmin + dt*(n/2)
!!!!!!!$    fcn%t = t
!!!!!!!$    s = smin
!!!!!!!$    tol = smin
!!!!!!!$    call snewton(fcn%funmean,t,tol,s,converge)
!!!!!!!$    do i = n/2,n
!!!!!!!$       t = tmin + dt*i
!!!!!!!$       s = rk4(fcn%funmean,t,s,dt)
!!!!!!!$       write(1,*) t,s
!!!!!!!$    end do
!!!!!!!$    close(1)
!!!!!!    stop
!!!!!
!!!!!    goto 13
!!!!!    open(1,file=filename)
!!!!!    do i = 0, n
!!!!!       t = tmin + dt*i
!!!!!       fcn%t = t
!!!!!       s = smax
!!!!!       do k = 1, 10
!!!!!          s = zeroin(fcn,smin,s,tol,iflag)
!!!!!          if( iflag == 2 ) then
!!!!!             s = smax / 2**k
!!!!!          else if( iflag == 0 ) then
!!!!!             call fcn%funmean%hess(t,s,hess)
!!!!!             d = hess(1,1)*hess(2,2) - hess(1,2)*hess(2,1)
!!!!!             write(*,*) smin,s,t,d,hess(1,1)
!!!!!             if( hess(2,2) < 0 .and. s < 1 ) exit
!!!!!!             exit
!!!!!             s = s / 2
!!!!!          else
!!!!!             write(*,*) t,s,iflag
!!!!!          end if
!!!!!!          if ( s < smin ) exit
!!!!!       end do
!!!!!       write(1,*) t,s,iflag
!!!!!!       stop
!!!!!    end do
!!!!!    close(1)
!!!!!    13 continue
!!!!!
!!!!!!    smin = dt
!!!!!!    open(1,file='/tmp/sssmin')
!!!!!!    do i = 0, n
!!!!!!       t = tmin + dt*i
!!!!!!       fcn%t = t
!!!!!!       s = zeroin(fcn,smin,smax,tol,iflag)
!!!!!!       write(1,*) t,s,iflag
!!!!!!    end do
!!!!!!    close(1)
!!!!!
!!!!!
!!!!!  end subroutine graph_smin
!!!!!
!!!!!
!!!!!  subroutine snewton(fun,t,stol,s,converge)
!!!!!
!!!!!    class(tsr1mean_function), intent(in), target :: fun
!!!!!    real(REAL64), intent(in) :: t, stol
!!!!!    real(REAL64), intent(in out) :: s
!!!!!    logical, intent(out) :: converge
!!!!!
!!!!!!    real(REAL64), parameter :: dlteps = 2*epsilon(s)
!!!!!    real(REAL64) :: fval,fv,fs,h22,d,rho,sp,ared,pred,dlt,one,r2,dtol
!!!!!    integer :: iter
!!!!!
!!!!!    converge = .false.
!!!!!    dlt = 0.5*s
!!!!!    rho = 0
!!!!!
!!!!!!    write(*,*) 'snewton init:',t,s
!!!!!
!!!!!    block
!!!!!      real(REAL64), dimension(2) :: fvec
!!!!!      real(REAL64), dimension(2,2) :: hess
!!!!!      integer :: i
!!!!!      real(REAL64) :: ss
!!!!!
!!!!!      open(17,file='/tmp/y')
!!!!!      do i = 1, 500
!!!!!         ss = 0.01*i
!!!!!         call fun%grad(t,ss,fvec)
!!!!!         call fun%hess(t,ss,hess)
!!!!!         fs = fvec(2)
!!!!!         h22 = hess(2,2)
!!!!!         write(17,*) ss,fun%val(t,ss),fs,h22
!!!!!      end do
!!!!!      close(17)
!!!!!    end block
!!!!!!    stop
!!!!!
!!!!!    do iter = 1, 3*precision(s)
!!!!!
!!!!!       call newstep(fun,t,s,fs,h22)
!!!!!
!!!!!       if( h22 < 0 ) then
!!!!!          d = -fs / h22
!!!!!       else
!!!!!          d = dlt
!!!!!       end if
!!!!!       if( abs(d) >= dlt ) d = sign(dlt,fs)
!!!!!
!!!!!       sp = s + d
!!!!!
!!!!!!       if( debug ) &
!!!!!            !            write(*,'(a,i2,1pg13.5,1p2g10.2,1pg11.3,1pe12.3,3l1)') &
!!!!!!       write(*,*) &
!!!!!!            "snewton,inc.,f,f': ",&
!!!!!!            iter,s,d,fs,h22,dlt
!!!!!
!!!!!       rho = 0
!!!!!       if( sp > 0 ) then
!!!!!          fv = fun%val(t,sp)
!!!!!!       if( debug ) write(17,*) s,fval
!!!!!
!!!!!          ared = fval - fv
!!!!!          pred = -d*(fs + 0.5*h22*d)
!!!!!          if( abs(pred) > 0 ) rho = ared / pred
!!!!!       end if
!!!!!!       write(*,*) ared,pred,rho
!!!!!
!!!!!       if( rho > 0.5e-3 ) then
!!!!!          s = sp
!!!!!          fval = fv
!!!!!
!!!!!          converge = abs(d) < stol/10
!!!!!          if( converge ) exit
!!!!!
!!!!!!          call newstep(fun,t,s,fs,h22)
!!!!!
!!!!!       end if
!!!!!
!!!!!       if( rho > 0.75 ) then
!!!!!          if( abs(d) > 0.9*dlt ) then
!!!!!!             dlt = min(2*dlt,dltmax)
!!!!!             dlt = 2*dlt
!!!!!          end if
!!!!!       else if( rho < 0.25 ) then
!!!!!!          dlt = max(dlt/2,dltmin)
!!!!!          dlt = dlt/2
!!!!!       end if
!!!!!!       dlt = min(dlt,s-smin)
!!!!!
!!!!!
!!!!!    end do
!!!!!
!!!!!!    stop
!!!!!
!!!!!  end subroutine snewton
!!!!!
!!!!!  subroutine newstep(fun,t,s,fs,h22)
!!!!!
!!!!!    class(tsr1mean_function), intent(in), target :: fun
!!!!!    real(REAL64), intent(in) :: t, s
!!!!!    real(REAL64), intent(out) :: fs, h22
!!!!!
!!!!!    real(REAL64), dimension(2) :: fvec
!!!!!    real(REAL64), dimension(2,2) :: hess
!!!!!
!!!!!    call fun%grad(t,s,fvec)
!!!!!    call fun%hess(t,s,hess)
!!!!!
!!!!!    fs = fvec(2)
!!!!!    h22 = hess(2,2)
!!!!!
!!!!!  end subroutine newstep
!!!!!

!!$  subroutine graphreg(fcn,t0,t1,s0,s1,filename)
!!$
!!$    use fmm
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL32), intent(in) :: t0,t1,s0,s1
!!$    character(len=*), intent(in) :: filename
!!$    integer, parameter :: n = 200
!!$    real(REAL64) :: t,dt,ds,smin,smax,tol,f,fmax,chi2,df
!!$    integer :: i,j
!!$
!!$    dt = (t1 - t0) / (n - 1)
!!$    ds = (s1 - s0) / (n - 1)
!!$
!!$    open(1,file=filename)
!!$    do i = 0, n
!!$       t = t0 + dt*i
!!$       fmax = 0
!!$       smax = 0
!!$       do j = 1, n
!!$          s = s0 + j*ds
!!$          f = sum(itukey((data - t)/s)) / size(data) + log(s)
!!$          if( f > fmax ) then
!!$             smin = s - ds
!!$             smax = s + ds
!!$             fmax = f
!!$          end if
!!$       end do
!!$       negtukey_t = t
!!$!       smin = ds
!!$!       smax = 5
!!$       tol = 1e-5
!!$       s = fmin(smin,smax,negtukey,tol)
!!$       chi2 = sum(itukey((data-t)/s)) / size(data)
!!$       df = sum(dtukey((data-t)/s), dtukey((data-t)/s) > 0) / size(data)
!!$       write(1,*) t,s,chi2,df
!!$       write(*,*) t,s,fmax
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graphreg
!!$
!!$  subroutine graphpole(fcn,x0,x1,y0,y1,filename)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL32), intent(in) :: x0,x1,y0,y1
!!$    character(len=*), intent(in) :: filename
!!$    integer, parameter :: nx = 200
!!$    integer, parameter :: ny = 200
!!$    real(REAL64) :: x,dx,y,dy,chi2
!!$    integer :: i,j
!!$
!!$    dx = (x1 - x0) / (nx - 1)
!!$    dy = (y1 - y0) / (ny - 1)
!!$
!!$    open(1,file=filename)
!!$    do j = 1, ny
!!$       y = y0 + dy*j
!!$       do i = 0, nx
!!$          x = x0 + dx*i
!!$          write(1,*) x,y,fcn%val(x,y), itukey((data(12)-x)/y)+log(y), &
!!$               (itukey((data(8)+0*0.03-x)/y)+itukey((data(12)-x)/y))/2+log(y), &
!!$               ((data(12)-x)/y)**2/2+log(y)
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graphpole
!!$
!!$  subroutine graph_ctn(fcn,x0,x1,y0,y1,filename)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL32), intent(in) :: x0,x1,y0,y1
!!$    character(len=*), intent(in) :: filename
!!$    integer, parameter :: nx = 200
!!$    integer, parameter :: ny = 200
!!$    real(REAL64) :: x,dx,y,dy,chi2,f
!!$    integer :: i,j,l,n
!!$
!!$    dx = (x1 - x0) / (nx - 1)
!!$    dy = (y1 - y0) / (ny - 1)
!!$
!!$    open(1,file=filename)
!!$    do j = 1, ny
!!$       y = y0 + dy*j
!!$       do i = 0, nx
!!$          x = x0 + dx*i
!!$          n = count(abs(tukey((data-x)/y)) < 1)
!!$          n = count(abs((data-x)/y) < 1)
!!$          f = sum(dtukey((data-x)/y))/size(data)
!!$          f = 0
!!$          do l = 1, size(data)
!!$             f = f + max(dtukey((data(l)-x)/y),real(0,REAL64))
!!$          end do
!!$          f = f / size(data)
!!$          chi2 = sum(tukey((data-x)/y)**2)/size(data)
!!$          chi2 = sum(itukey((data-x)/y))/size(data)
!!$          write(1,*) x,y,n,chi2,f
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graph_ctn
!!$
!!$
!!$  subroutine graphcliffs(fcn,x0,x1,y0,y1,filename)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL32), intent(in) :: x0,x1,y0,y1
!!$    character(len=*), intent(in) :: filename
!!$!    real(REAL64), dimension(:), pointer :: data
!!$    real(REAL64), dimension(:), allocatable :: t,s
!!$    integer, parameter :: nx = 200
!!$    integer, parameter :: ny = 200
!!$    real(REAL64) :: x,dx,y,dy,smax
!!$    integer :: i,j,l,n
!!$
!!$!    data = fcn%data
!!$    n = size(data)
!!$    allocate(t(2*n*n),s(2*n*n))
!!$
!!$    dx = (x1 - x0) / (nx - 1)
!!$    dy = (y1 - y0) / (ny - 1)
!!$
!!$    open(1,file=filename)
!!$
!!$    ! all possible interscetions
!!$    n = 0
!!$    do i = 1, size(data)
!!$       do j = 1, size(data)
!!$          if( i /= j ) then
!!$             n = n + 1
!!$             t(n) = (data(i) + data(j)) / 2
!!$             s(n) = (data(i) + data(j)) / (2*3.0)
!!$             if( s(n) > 0 ) then
!!$                write(1,*) t(n),s(n)
!!$                write(1,*) -t(n),s(n)
!!$             end if
!!$          end if
!!$!          n = n + 1
!!$!          t(n) = -(data(i) + data(j)) / 2
!!$!          s(n) = (data(i) + data(j)) / (2*3.0)
!!$!          if( s(n) > 0 ) write(1,*) t(n),s(n)
!!$       end do
!!$    end do
!!$
!!$    do j = 1, 0*ny
!!$       y = y0 + dy*j
!!$       do i = 0, nx
!!$          x = x0 + dx*i
!!$
!!$          ! only ones in given range
!!$          smax = -1
!!$          do l = 1, n
!!$             if( x < t(l) .and. t(l) < x + dx ) then
!!$                smax = s(l)
!!$             end if
!!$          end do
!!$          if( smax > 0 ) write(1,*) x,smax
!!$
!!$       end do
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graphcliffs
!!$
!!$  subroutine graph_stringart(fcn,x0,x1,y0,y1,filename)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL32), intent(in) :: x0,x1,y0,y1
!!$    character(len=*), intent(in) :: filename
!!$!    real(REAL64), dimension(:), pointer :: data
!!$    real(REAL64), dimension(:), allocatable :: t,s
!!$    integer, parameter :: nx = 200
!!$    integer, parameter :: ny = 200
!!$    real(REAL64) :: x,dx,y,dy,smax
!!$    integer :: i,j,l,n
!!$
!!$!    data = fcn%data
!!$    n = size(data)
!!$    allocate(t(2*n*n),s(2*n*n))
!!$
!!$    dx = (x1 - x0) / (nx - 1)
!!$    dy = (y1 - y0) / (ny - 1)
!!$
!!$    ! all possible lines
!!$    open(1,file=filename)
!!$    do i = 1, size(data)
!!$       write(1,*) data(i),0,-5*3,5
!!$       write(1,*) data(i),0,+5*3,5
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine graph_stringart



!!$  subroutine graphdev(fcn,x0,x1,y0,y1,filename)
!!$
!!$    class(sr1_mean), intent(in) :: fcn
!!$    real(REAL32), intent(in) :: x0,x1,y0,y1
!!$    character(len=*), intent(in) :: filename
!!$    integer, parameter :: nx = 200
!!$    integer, parameter :: ny = 200
!!$    real(REAL64) :: x,dx,y,dy,d
!!$    integer :: i,j
!!$
!!$    dx = (x1 - x0) / (nx - 1)
!!$    dy = (y1 - y0) / (ny - 1)
!!$
!!$    open(1,file=filename)
!!$    open(2,file='/tmp/c')
!!$    do j = 1, ny
!!$       y = y0 + dy*j
!!$       do i = 1, nx
!!$          x = x0 + dx*i
!!$          d = fcn%stderr([x,y])
!!$          write(1,*) x,y,d
!!$          if( abs(d - 0.9) < 0.01 ) write(2,*) x,y
!!$       end do
!!$       write(1,*)
!!$    end do
!!$    close(1)
!!$
!!$    close(2)
!!$
!!$  end subroutine graphdev


!!$  real(REAL64) function chi2fun(s)
!!$    use rfun
!!$    use rfun_hampel
!!$    real(REAL64), intent(in) :: s
!!$    !    chi2fun = sum(itukey((data-chi2t)/s))/size(data) - chi2par
!!$    chi2fun = sum(tukey((data-chi2t)/s)**2)/size(data) - chi2par
!!$    chi2fun = sum(tukey((data-chi2t)/s)**2)/size(data) - chi2par/s
!!$    chi2fun = sum(tukey((data-chi2t)/s)*(data-chi2t)/s)/size(data) - chi2par
!!$    chi2fun = sum(pshuber((data-chi2t)/s)*(data-chi2t)/s)/size(data) - chi2par
!!$    chi2fun = sum(pshuber((data-chi2t)/s)*(data-chi2t)/s**2)/size(data) - chi2par/s
!!$    chi2fun = sum(pshuber((data-chi2t)/s)**2)/(size(data)-1) - chi2par
!!$    chi2fun = sum(pshampel((data-chi2t)/s)**2)/(size(data)-1) - chi2par
!!$    neval = neval + 1
!!$  end function chi2fun
!!$
!!$  real(REAL64) function chi2fun1(s)
!!$    use rfun
!!$    real(REAL64), intent(in) :: s
!!$    real(REAL64) :: r,ss
!!$    integer :: i
!!$
!!$    ss = 0
!!$    do i = 1,size(data)
!!$       r = (data(i)-chi2t)/s
!!$       if( abs(r) < 5 ) then
!!$          ss = ss + pshuber(r)**2
!!$       end if
!!$    end do
!!$    chi2fun1 = ss/(size(data)-1) - chi2par
!!$
!!$  end function chi2fun1
!!$
!!$
!!$  real(REAL64) function chi2min(s)
!!$    real(REAL64), intent(in) :: s
!!$    chi2min = sum(itukey((data-chi2t)/s))/size(data) + log(s)
!!$    neval = neval + 1
!!$  end function chi2min
!!$
!!$  real(REAL64) function scalemin(s)
!!$    real(REAL64), intent(in) :: s
!!$    scalemin = sum(dtukey((data-t)/s), dtukey((data-t)/s) > 0) / size(data) - 0.3
!!$    neval = neval + 1
!!$  end function scalemin
!!$
!!$  real(REAL64) function negtukey(s)
!!$    real(REAL64), intent(in) :: s
!!$    negtukey = -sum(itukey((data - negtukey_t)/s)) / size(data) - log(s)
!!$  end function negtukey
!!$
!!$  elemental real(REAL64) function Phi(x,t,s)
!!$    real(REAL64), intent(in) :: x,t,s
!!$    Phi = (1 + erf((x-t)/s/sqrt(2.0)))/2
!!$  end function Phi
!!$
!!$  real(REAL64) function fsuff(s)
!!$
!!$    real(REAL64), intent(in) :: s
!!$    real(REAL64), dimension(2,2) :: hess
!!$
!!$    call funmean%hess(t,s,hess)
!!$    fsuff = s**2*hess(1,1) - 0.82
!!$
!!$  end function fsuff

!!$  function fscale(fcn,t,smin,smax) result(s)
!!$
!!$    real(REAL64), intent(in) :: t,smin,smax
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL64) :: a,b,fa,fb,s,f
!!$    integer :: iter
!!$
!!$    a = smin
!!$    b = smax
!!$    fa = fscale_fun(fcn,t,a)
!!$    fb = fscale_fun(fcn,t,b)
!!$    write(*,*) a,b,fa,fb
!!$    do iter = 1, 20
!!$       s = (a + b) / 2
!!$       f = fscale_fun(fcn,t,s)
!!$       write(*,*) s,f
!!$       if( abs(f) < 1e-3 .or. (b-a)/2 < 0.1 ) return
!!$       if( fa*f < 0 ) then
!!$          b = s
!!$          fb = f
!!$       else if( f*fb < 0 ) then
!!$          a = s
!!$          fa = f
!!$       end if
!!$    end do
!!$
!!$  end function fscale
!!$
!!$  real(REAL64) function fscale_fun(fcn,t,s)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64), dimension(2,2) :: hess
!!$
!!$    call fcn%hess(t,s,hess)
!!$    fscale_fun = s**2*hess(1,1) - 0.8
!!$
!!$  end function fscale_fun

!!$  real(REAL64) function inorm(x)
!!$
!!$    real(REAL64), intent(in) :: x
!!$
!!$    if( abs(2*x-1) < 1 ) then
!!$       inorm = sqrt(real(2,REAL64))*inverf(2*x - 1)
!!$    else
!!$       inorm = 0
!!$!       write(*,*) x,2*x-1
!!$!       stop
!!$    end if
!!$
!!$  end function inorm

!!$  subroutine fminsig(t,smin,smax,tol,s)
!!$
!!$    use fmm
!!$    real(REAL64), intent(in) :: t,smin,smax,tol
!!$    real(REAL64), intent(out) :: s
!!$
!!$    fminsig_t = t
!!$    s = fmin(smin,smax,flep,tol)
!!$
!!$  end subroutine fminsig
!!$
!!$
!!$  real(REAL64) function flep(s)
!!$
!!$    use kstestscale
!!$
!!$    real(REAL64), intent(in) :: s
!!$
!!$    flep = lepage(data,fminsig_t,s)
!!$
!!$  end function flep


!!$  subroutine ridgemax(fcn,t,s,delta,tol,pvec,converge)
!!$
!!$!    use trustsr1mean
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL64), intent(in out) :: t, s
!!$    real(REAL64), intent(in) :: delta, tol!, lam
!!$    real(REAL64), dimension(2), intent(in) :: pvec
!!$!    real(REAL64), dimension(2,2), intent(in) :: hess
!!$    logical, intent(out) :: converge
!!$
!!$    real(REAL64), dimension(2) :: x, xp, p, fvec
!!$    real(REAL64), dimension(2,2) :: hess
!!$    real(REAL64) :: fval, fv, g, d, dlt, ared, pred, rho, pHp
!!$    integer :: iter, i
!!$
!!$    write(*,*) 'Ridge start with ',t,s,' dir:',real(pvec)
!!$
!!$!    open(17,file='/tmp/ridgemax1')
!!$!    do i = -100,100
!!$!       d = 0.05*delta*i
!!$!       write(17,*) t+d*pvec(1),s+d*pvec(2),fcn%val(t+d*pvec(1),s+d*pvec(2))
!!$!    end do
!!$!    close(1)
!!$!    stop
!!$
!!$
!!$    converge = .false.
!!$    dlt = delta
!!$    fval = fcn%val(t,s)
!!$    call fcn%hess(t,s,hess)
!!$
!!$!    if( .not. (lam(1)*lam(2) < 0) ) stop 'Something wrong.'
!!$
!!$!    do i = 1, 2
!!$!       if( lam(i) < 0 ) q = Qs(:,i)
!!$!    end do
!!$
!!$    !   if( debug )
!!$    open(17,file='/tmp/ridgemax')
!!$
!!$    x = [ t, s ]
!!$    do iter = 1, 3*precision(s)
!!$
!!$       call fcn%grad(x(1),x(2),fvec)
!!$
!!$       g = dot_product(fvec,pvec)
!!$       pHp = dot_product(pvec,matmul(hess,pvec))
!!$
!!$       if( php < 0 ) then
!!$          ! we're looking for a maximum
!!$!          d = -g / qHq
!!$          d = -g / pHp
!!$       else
!!$          d = dlt
!!$       end if
!!$       if( abs(d) >= dlt ) d = sign(dlt,g)
!!$
!!$       p = d*pvec
!!$       xp = x + p
!!$       write(*,*) x
!!$       write(*,*) xp
!!$
!!$       fv = fcn%val(xp(1),xp(2))
!!$       !       if( debug )
!!$       write(17,*) xp,fval
!!$
!!$       ared = fval - fv
!!$!       pred = -dot_product(g + 0.5*qHq*p,p)
!!$!       pred = -dot_product(fvec + 0.5*lam*p,p)
!!$       pred = -dot_product(fvec + 0.5*matmul(hess,p),p)
!!$       rho = ared / pred
!!$!       rho = 0.5
!!$
!!$       if( rho > 0.5e-3 ) then
!!$          x = xp
!!$          fval = fv
!!$
!!$          converge = abs(d) < tol !.and. qHq < 0
!!$          if( converge ) then
!!$             t = x(1)
!!$             s = x(2)
!!$             exit
!!$          end if
!!$
!!$          call fcn%hess(x(1),x(2),hess)
!!$
!!$       end if
!!$
!!$!       if( debug ) &
!!$            !            write(*,'(a,i2,1pg13.5,1p2g10.2,1pg11.3,1pe12.3,3l1)') &
!!$            write(*,*) &
!!$            "ridgemax,inc.,f,f': ",&
!!$            iter,xp,d,g,pHp,fval,rho,dlt
!!$
!!$
!!$
!!$!       write(*,*) dlt,rho,d*fs, 0.5*h22*d**2
!!$       if( rho > 0.75 ) then
!!$          if( abs(d) > 0.9*dlt ) then
!!$             dlt = min(2*dlt,s)
!!$          end if
!!$       else if( rho < 0.25 ) then
!!$          dlt = dlt/2
!!$       end if
!!$
!!$    end do
!!$
!!$    close(17)
!!$    write(*,*) 'Ridge max at ', t,s,converge
!!$
!!$
!!$  end subroutine ridgemax


!!$  subroutine fline(fcn,t,a,b,s)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL64), intent(in) :: t,a,b
!!$    real(REAL64), intent(in out) :: s
!!$
!!$    real(REAL64) :: f,g,df,p,c,alpha
!!$    real(REAL64), dimension(2) :: fvec
!!$
!!$    s = 0.8
!!$    f = fcn%val(t,s)
!!$    call fcn%grad(t,s,fvec)
!!$    df = fvec(2)
!!$    if( df > 0 ) then
!!$       p = -f/df
!!$    else
!!$       p = 0.1*s
!!$    end if
!!$    write(*,*) s,p,f,df
!!$
!!$    c = 1e-4
!!$    alpha = 1
!!$    do
!!$
!!$       g = fcn%val(t,s + alpha*p)
!!$       write(*,*) g,f + c*alpha*df*p
!!$       if( g <= f + c*alpha*df*p ) exit
!!$       alpha = 0.9*alpha
!!$
!!$    end do
!!$
!!$    s = s + alpha*p
!!$
!!$
!!$
!!$  end subroutine fline

end program testsr1mean


!!$module zeroin_smin
!!$
!!$  use sr1mean_function
!!$  use zeroinmod
!!$  use iso_fortran_env
!!$  implicit none
!!$
!!$  type, extends(zeroin_function) :: zeroin_data
!!$     class(tsr1mean_function), pointer :: funmean
!!$     real(REAL64) :: t
!!$   contains
!!$     procedure :: f => zeroin_fun
!!$  end type zeroin_data
!!$
!!$contains
!!$
!!$  real(REAL64) function zeroin_fun(this,x)
!!$
!!$    class(zeroin_data) :: this
!!$    real(REAL64), intent(in) :: x
!!$    real(REAL64), dimension(2,2) :: hess
!!$    real(REAL64), dimension(2) :: fvec
!!$    real(REAL64), parameter :: zero = real(0,REAL64)
!!$    real(REAL64) :: d,sgn
!!$    integer :: i
!!$
!!$    ! original
!!$    zeroin_fun = this%funmean%nfrac(this%t,x) - 0.25
!!$!    call this%funmean%hess(this%t,x,hess)
!!$!    d = hess(1,1)*hess(2,2) - hess(1,2)*hess(2,1)
!!$!    zeroin_fun = d - 0.01
!!$    return
!!$
!!$    call this%funmean%hess(this%t,x,hess)
!!$    call this%funmean%grad(this%t,x,fvec)
!!$    sgn = 1
!!$    if( d < 0 .or. hess(1,1) < 0 ) sgn = -1
!!$    zeroin_fun = sgn*abs(hess(1,1))*abs(d)
!!$!    write(*,*) this%t,x,x*fvec(2)
!!$    zeroin_fun = x*fvec(2)
!!$    return
!!$
!!$    open(17,file='/tmp/x')
!!$    do i = 1,500
!!$
!!$       call this%funmean%hess(this%t,real(0.01*i,REAL64),hess)
!!$       call this%funmean%grad(this%t,real(0.01*i,REAL64),fvec)
!!$       d = hess(1,1)*hess(2,2) - hess(1,2)*hess(2,1)
!!$       sgn = 1
!!$       if( d < 0 .or. hess(1,1) < 0 ) sgn = -1
!!$       hess = x**2*hess
!!$       write(17,*) 0.01*i,hess(1,1),hess(2,2),fvec
!!$    end do
!!$    close(17)
!!$    stop
!!$  end function zeroin_fun
!!$
!!$end module zeroin_smin
