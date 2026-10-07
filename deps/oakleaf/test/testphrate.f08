
program testphrate

  use oakleaf
  use noise
  use toolbox
  use testsuite
  use iromberg
  use iso_fortran_env

  implicit none

  real(REAL64), parameter :: PI = real(3.14159265358979323844,REAL64)
  real(REAL64), parameter :: RAD = 180 / PI

  integer, parameter :: rp = selected_real_kind(15)
  character(len=15), parameter :: fmt = '(a,l3)'
  integer :: voltage


  call get_environment_variable("VOLTAGE",status=voltage)

  if( .not. (voltage == 0) ) then

     call tester(build)

  else
     call randomnoise
!     call tester(bmuma)
!     call tester(poisson_mix)
!     call tester(poisson_mix_contamined)
!     call tester(poisson_mean)
     call tester(poisson_histo)

!     call lssolve
!     call lsmin
!     call lsmin2d
     !write(*,*) romberg(phfun,0.001_rp,200.0_rp,1d-5)
!     call tabelate
!   call tabelate_params
     !call phtest
     !call maxtab
     !     call marginal
     !     call themean
!     call tabrfun
     !      call InvPhi
     !     call quantilephi
!     call plotqPhi
  end if


contains


  subroutine build

    real(rp), dimension(:), allocatable :: x,y
    type(suite) :: ph
    logical :: cond

    integer :: n,l
    real(rp) :: mean

    mean = 10**4

    x = [ 1e4, 2e4, 3e4, 4e4, 5e4 ]
    y = [ 2.1e4, 1.9e4, 2.1e4, 1.9e4, 2.0e4 ]

!    n = 10
!    allocate(x(n),y(n))

!    do l = 1, n
!       x(l) = pnoise(mean)
!       y(l) = pnoise(mean)
!    end do

!    write(*,*) 'Po(10)  n=',n
    call ph%phstat(x,y)
    cond = abs(ph%t - 1) < 5*ph%dt
    write(*,'(a,l2)') 'TEST build PASS:',cond
    write(*,*)
    deallocate(x,y)

  end subroutine build


  subroutine bmuma

    real(rp), dimension(6) :: x,dx,y,dy
    type(suite) :: std, rob, qt, stdw, robw, poi
    logical :: cond

    x = [1.58682E+07, 4.53635E+05, 4.41066E+05, 4.14017E+05, &
         4.55078E+05, 3.11341E+05]
    dx = [ 3.990E+03, 712., 704., 682., 713., 608. ]
    y = [ 1.51535E+06, 44100., 42143., 39892., 42909., 30737. ]
    dy = [ 1.233E+03, 221., 216., 209., 218., 189. ]

    write(*,*) 'BM UMa'
    call std%stdstat(x/y)
    call rob%robstat(x/y)
    call qt%qtstat(x/y)

    call stdw%stdwstat(x/y,sqrt(dx**2 + dy**2))
    call robw%robwstat(x/y,sqrt(dx**2 + dy**2))
    call poi%poistat(y,dy,x,dx)
    cond = abs(poi%t - 10.463) < 5*poi%dt
    write(*,'(a,l2)') 'TEST bmuma PASS:',cond
    write(*,*)

    block
      real(REAL64) :: t
      call phrate(x,y,t,verbose=.true.)
    end block

  end subroutine bmuma


  subroutine poisson_mix

    use quantilmean

    real(rp), dimension(:), allocatable :: x,dx,y,dy,res,dres,pht,cts,means, &
         mju, sigma, tt, qnt, val, beta, sums
    type(suite) :: std, rob, qt, stdw, robw, poi
    logical :: cond

    integer :: m,n,k,l,i,mg,ns
    real(rp) :: mean, mag, fac, sig, sigx, sigy, zero, q, s, r, e, t, rmax, &
         sigu, sigv, phi, rho

    mean = 1e4
    fac = 2
    zero = 0.0
    sigx = 2e3
    sigy = 2e3
    q = 1

    rho = sqrt(sigx**2*fac**2 + sigy**2)
    write(*,*) 'Expected rho:', real(rho)
    if( sigx > 0 .and. sig > 0 ) then
       phi = atan2(sigy,fac*sigx)
       write(*,*) 'Expected phi:',real(RAD*phi)
    else
       phi = 0
    end if

    n = 13
    allocate(x(n),dx(n),y(n),dy(n),res(n),dres(n),pht(n),cts(n),means(n),sums(2*n))

    n = 0
    do mg = 1000,1600, 100
!    do mg = 1,130000
       mag = mg / 100.0
       mean = 10**((25 - mag)/2.5)
!       mag = 25 - 2.5*log10(mean)
!       m = max(int(exp(mag*1.1)/25.0),1)
  !     m = max(int(exp(mag*1.1)/20e3),1)
!       m = max(int(exp(mag*1.1)/1e3),1)
    m = max(int(exp(mag*1.1)/5e5),1)
!      m = size(x)
!      mean = 1e4
!       m = 1
       write(*,*) mean,mag,m
       ns = 0
       sums = 0
       do l = 1, m
          if( n == size(x) ) exit
          n = n + 1
          means(n) = mean
          pht(n) = pnoise(mean)
          cts(n) = pnoise(fac*mean)
          x(n) = pht(n) + q*gnoise(zero,sigx)
          y(n) = cts(n) + q*gnoise(zero,sigy)
          ns = ns + 2
          sums(ns-1) = x(n)
          sums(ns) = y(n) / fac
!          if( x(n) < 0 .or. y(n) < 0 ) n = n - 1
          e = x(n) + q*gnoise(zero,sigx)**2
          if( e > 0 ) then
             dx(n) = sqrt(e)
          else
             dx(n) = 1e-7
          end if
          e = y(n) + q*fac**2*gnoise(zero,sigy)**2
          if( e > 0 ) then
             dy(n) = sqrt(e)
          else
             dy(n) = 1e-7
          end if
!          if( .not. (x(n) > 0 .and. y(n) > 0) ) n = n - 1
       end do

       mean = sum(sums(1:ns))/ns
       sig = sqrt(sum((sums(1:ns) - mean)**2)/(ns - 1))
!       write(*,*) 'E(mean)=',real(mean),'+-',real(sig/sqrt(1.0*ns)),real(sig)

       if( n == size(x) ) exit
    end do

    if( n /= size(x) ) then
       write(*,*) 'BIG warning: n /= size(x)',n,size(x)
       stop 'BIG'
    end if


    block
      real(REAL64), dimension(:), allocatable :: xpdf, ypdf, res, mju
      real(REAL64) :: tau,s

      tau = fac!+0.001
      allocate(xpdf(n),ypdf(n),mju(n))
      res = (x - means) / (sqrt(means*tau*(1+tau) + rho**2)/tau)
      res = (y - tau*means) / (sqrt(means*tau*(1+tau) + rho**2))*(1 + tau)
      res = ((x - means)**2*tau**2 + (y - tau*means)**2)/sqrt(means*tau*(1+tau) + rho**2)**2
      res = (x - means)*tau/sqrt(means*tau*(1+tau) + rho**2)
      !      res = (y - tau*means)/sqrt(means*tau*(1+tau) + rho**2)
      res = (x - means) / sqrt(means + (rho/tau*cos(phi))**2)
      res = (y - tau*means) / sqrt(means*tau + (rho*sin(phi))**2)
      res = ((x - means) / sqrt(means + sigx**2))**2 &
           +((y - tau*means) / sqrt(means*tau + sigy**2))**2
      res = ((x - means) / sqrt(means + (rho/tau*cos(phi))**2))**2 &
           +((y - tau*means) / sqrt(means*tau + (rho*sin(phi))**2))**2

      do i = 1, n
         mju(i) = means(i) + gnoise(zero,5e3*rho/means(i))
      end do
      mju = (y/tau + x)/2
!      mju = pht
!      res = (x - mju) / sqrt(mju + (rho/tau*cos(phi))**2)
!      res = (y - tau*mju) / sqrt(mju*(tau +2*tau**2*0) + 0*rho**2 + 0*(rho*cos(phi))**2 + 0*mju**2*1e-3**2)
!      res = ((x - mju) / sqrt(mju + (rho/tau*cos(phi))**2))**2 &
!           +((y - tau*mju) / sqrt(mju*tau + (rho*sin(phi))**2))**2

!      res = sqrt(2.0)*(x - mju) / sqrt(mju + (rho/tau*cos(phi))**2)
      res = 2.0*(x - mju) / sqrt(mju*(1+1/tau) + (rho/tau)**2)
!      res = (pht - mju) / sqrt(pht + 0*(rho/tau*cos(phi))**2)
!      res = 2*(y - mju*tau) / sqrt(mju*(tau + tau**2) + rho**2)


      ! test for rho determination
      res = (y - tau*x)/sqrt((y + tau*x)/2*tau + rho**2)
      res = 2*(y - tau*mju)/sqrt(mju*tau*(tau+1) + rho**2)

      call ecdf(res,xpdf,ypdf)

      open(1,file='/tmp/tcdf')
      do i = 1, size(xpdf)
         write(1,*) xpdf(i),ypdf(i)
      end do
      close(1)

      s = sum((x-means)**2+(y/tau-means)**2) / 2
      write(*,*) 'Expected d:',real(sqrt(s/n))

!      write(*,*) real(rho),real(RAD*phi)

!      stop
    end block


    goto 22

    t = 200
    write(*,*) 'E(mu):',sum((x+t*y)/(1+t**2))/n

    ! distributions of Emos
    e = 1e3
    call histogram((x+t*y)/(1+t**2),mean-5*e,mean+100*e,1000,'/tmp/emu')

    ! distributions of weighted Emos
    e = 1

    sigy = sqrt(fac*mean + rho**2*0)
    sigx = sqrt(mean + (rho/fac)**2)

    call histogram((y-fac*mean)/sqrt(fac*mean+sigy**2),-5*e,5*e,400,'/tmp/emuy')
    call histogram(fac*(x-mean)/sqrt(mean+sigx**2)/sqrt(fac),-5*e,5*e,400,'/tmp/emux')
!    call histogram(((x-mean)**2 + (y - fac*mean)**2/fac**2)/sqrt(2*mean + (rho/fac)**2)**2,0*e,10*e,400,'/tmp/emur')
!    call histogram(((x-mean)**2*(t*mean+rho**2*sin(phi)**2) + (y - fac*mean)**2*(mean + (rho*cos(phi))/tau)**2))/(t*mean+rho**2*sin(phi)**2)*(mean + (rho*cos(phi))/tau)**2,-1*e,20*e,100,'/tmp/emur')
    !    call histogram((x-mean)**2 + (y-t*mean)**2/t**2,-1*e,2e6*e,100,'/tmp/emur')

    res = ((x-mean)**2*fac**2+(y-fac*mean)**2)/sqrt(mean*fac*(1+fac)+rho**2)**2
!    res = ((x-mean)**2+(y-fac*mean)**2/fac**2)/sqrt(2*mean+(rho/fac)**2)**2
!    call histogram(res,-1e6*e,1e7*e,100,'/tmp/emur')
    call histogram(res,-e,10*e,100,'/tmp/emur')

    sigy = sqrt(fac*mean + rho**2)
    sigx = sqrt(mean + (rho/fac)**2*0)
    call histogram2d((y-fac*mean)/sqrt(fac*mean+sigy**2), &
         (x-mean)/sqrt(mean+sigx**2), &
         -5*e,5*e,-5*e,5*e,100,'/tmp/emu2d')



    ! distributions of mu
    call histogram(y,5d2,5d6,1000,'/tmp/ydis')
    call histogram(x,5d2,5d6,1000,'/tmp/xdis')
!    e = 1
!    call histogram((y-2e4)/sqrt(2e4+sigy**2),-5*e,+5*e,100,'/tmp/ydis')
!    call histogram((x-1e4)/sqrt(1e4+sigx**2),-5*e,+5*e,100,'/tmp/xdis')

!    allocate(mju(n),sigma(n),tt(n),qnt(n),val(n),beta(n))
    e = 1
    t = 2
    phi = atan(sigu,t*sigv)
    tt = cts / pht
    mju = (x+t*y)/(1+t**2)
    s = sqrt(sigx**2*fac**2 + sigy**2)
    sigma = sqrt(mju*(t+tt**2) + s**2*(1 + ((tt/t)**2 - 1)))
    ! distribution of exponent
    call histogram((mju*(t-tt)/sigma),-1*e,1*e,1000,'/tmp/edis')
    ! distribution of base
    e = 1
    s = 0
    sigma = sqrt(mju*(t+tt**2) + s**2*(1 + ((tt/t)**2 - 1)))
    beta = mju*(t*(tt+1)+s**2*(1+(tt/t-1)*cos(phi)**2))
    call histogram((mju/sigma)*beta/sigma**2,&
         40-6*e,40+6*e,1000,'/tmp/bdis')
    write(*,*) '**',mju(1),tt(1),sigma(1),beta(1)
!    stop 0
    e = 0.1
    call histogram(tt,2-1*e,2+1*e,100,'/tmp/ttdis')

!    block
!      write(*,*) 'minmax:',minval(y/x),maxval(y/x)
!      call ecdf(y/x,val,qnt)
!      open(1,file='/tmp/ttcdf')
!      do i = 1, size(tt)
!         write(1,*) val(i),qnt(i)
!      end do
!      close(1)
!    end block

!    deallocate(mju,sigma,tt,qnt,val)

!    t = 2
!    do i = 1, n
!       write(*,*) real(y(i)),real(x(i)),(x(i)+t*y(i))/(1+t**2)
!    end do
!    stop 0

!    open(1,file='/tmp/yx')
!    do i = 1, n
!       write(1,*) x(i),y(i),pht(i),cts(i),means(i)
!       write(1,*) means(i),y(i)/x(i),y(i),x(i)
!    end do
!    close(1)



    ! distribution of ratios
!    e = 0.1
!    call histogram(y/x,fac-10*e,fac+10*e,100,'/tmp/rq')

!    write(*,*) 'E(y),E(x):',sum(y)/n,sum(x)/n
!    write(*,*) 'E(y)/E(x):',sum(y)/sum(x)


22  continue

    ! estimate
    call phrate(y,x,t,verbose=.true.)

    sigu = sigy
    sigv = sigx

    ! distribution of residuals
!    e = 1
!    res = (y - t*x)/sqrt(y + t**2*x + sigu**2 + t**2*sigv**2)
!    write(*,*) t, res(1:3)
!    call histogram(res,-5*e,5*e,40,'/tmp/rr')

    ! distribution of beta
!    e = 2e2
!    res = (y*x*(1+t) +y*t*sigv**2 + x*sigu**2)/ &
!         (y + t**2*x + sigu**2 + t**2*sigv**2)**1.5
!    call histogram(res,1d0,e,40,'/tmp/beta')

    e = 1
    call histogram(y/x,fac-3*e,fac+3*e,100,'/tmp/fdis')

    ! mean
!    write(*,*) 'mean:',sum(y(:n)/x(:n))/n

!    block
!      real :: t,s
!      t = sum(y/x) / size(tt)
!      s = sqrt(sum((y/x - t)**2)/n)
!      write(*,*) 'mean, stdsig:',real(t),real(s)
!    end block



    stop 0


    block

      real(REAL64), dimension(:), allocatable :: tt,sig,beta,sig_t,beta_t,psi
      real(REAL64) :: c2, t, s, r, f, s1
      integer :: i,j,n

      n = size(y)
      allocate(tt(n),psi(n),sig(n),beta(n),sig_t(n),beta_t(n))
      tt = y / x
      r = sqrt(sigx**2*fac**2 + sigy**2)
      f = atan2(sigy,fac*sigx)
      c2 = cos(f)**2
      write(*,*) 'check:',r,f*57.3,count(tt < 1.5),count(tt > 3.1)
      open(1,file='/tmp/rhynk')
      do i = 1, 500, 1
         t = 0.01*i

         sig = sqrt(mju*(t + tt**2) + r**2*(1 + ((tt/t)**2 - 1)*c2))
         beta = (mju*t*(tt + 1) + r**2*(1 + (tt/t - 1)*c2))
         sig_t = (mju - 2*r**2*(tt/t)**2/t*c2)/(2*sig)
         beta_t = mju*(tt + 1) - r**2*tt/t**2*c2

!         sig = sqrt(mju*(fac + tt**2) + r**2*(1 + ((tt/fac)**2 - 1)*c2))
!         beta = (mju*fac*(tt + 1) + r**2*(1 + (tt/fac - 1)*c2))
!         sig_t = (mju - 2*r**2*(tt/fac)**2/fac*c2)/(2*sig)
!         beta_t = mju*(tt + 1) - r**2*tt/fac**2*c2

         res = mju*(t - tt)/sig
!         if( t < 1.5 ) then
!            psi = 1.8
!         else if( t > 2.4 ) then
!            psi = -3.4
!         else
!            psi = -res*(mju/sig) + res**2*(sig_t/sig) &
!                 + beta_t/beta - real(3.0,REAL64)*(sig_t/sig)
!         end if
!         where( tt < 1.5 )
!            psi = 1.8
!         elsewhere( tt > 2.4 )
!            psi = -3.4
         where( t-tt < -1 )
!         where( res < -1 )
            psi = -1.8*tanh((-1+t)/0.5)
            where( psi > 0 )
               psi = 0
            end where
            !         elsewhere( tt-t < -1 )
!            psi = -1.8
         elsewhere( t-tt > 1 )
!         elsewhere( res > 1 )
            psi = 3.4*tanh((3-t)/0.5)
            where( psi < 0 )
               psi = 0
            end where
 !           psi = 3.4
         elsewhere
!            psi = -res*(mju/sig) + res**2*(sig_t/sig) &
!                 + beta_t/beta - real(3.0,REAL64)*(sig_t/sig)
            psi = res*(mju/sig) - (res**2 - 3)*(sig_t/sig) &
                 - beta_t/beta !- real(3.0,REAL64)*(sig_t/sig)
            where( t < 0.1 .or. tt < 0.1 ) psi = 0
         end where
!         where( tt-t > 1 )
!            psi = 0
!         end where
!         where( tt-t < -1 )
!            psi = 0
!         end where
!         where( tt < 0.1 )
!            psi = 0
!         end where
!         psi = res/sig
         s = -sum(res*(mju/sig)) + sum(res**2*(sig_t/sig)) + &
              sum(beta_t/beta) - real(3.0,REAL64)*sum(sig_t/sig)
         s1 = -sum(res*(mju/sig)) + sum(res**2*(sig_t/sig))
         write(1,*) t, -s/n,-s1/n,sum(res*(mju/sig))/n,sum(psi)/n
      end do
      close(1)

    end block





    res = y - fac*x

    ! the proper ways
!    dres = sqrt(y + x*fac**2 + fac**2*sigx**2 + fac**2*sigy**2)
    dres = sqrt(y + x*fac**2 + fac**2*sigx**2 + sigy**2)
    call histogram(res/dres,-5d0,5d0,1000,'/tmp/q')
    dres = sqrt(y*(1 + fac) + fac**2*sigx**2 + sigy**2)
    call histogram(res/dres,-5d0,5d0,1000,'/tmp/q1')

    ! LS errors
    dres = sqrt(dy**2 + dx**2*fac**2)
    call histogram(res/dres,-5d0,5d0,1000,'/tmp/g')

    ! distribution of errors
!    dres = sqrt(dy**2 + dx**2*fac**2)
    dres = sqrt(y + x*fac**2 + fac**2*sigx**2 + fac**2*sigy**2)
    call histogram(dres,0d3,5.5d3,1000,'/tmp/dd')

    ! distribution of stars by counts
    call histogram(y/dres,0d6,2.1d3,1000,'/tmp/cd')

    ! distribution by estimated parameters
    call histogram(res/sqrt(y+t**2*x+r**2),-5d0,5d0,100,'/tmp/ht')

    block
      real(rp), dimension(:), allocatable :: xcdf, ycdf
      integer :: i
      allocate(xcdf(n),ycdf(n))
      call ecdf(y/x,xcdf,ycdf)
      q = 0.5;  write(*,*) q, quantile(q,xcdf,ycdf)
      q = 0.25;  write(*,*) q, quantile(q,xcdf,ycdf)
      q = 0.75;  write(*,*) q, quantile(q,xcdf,ycdf)
      open(1,file='/tmp/rcdf')
      do i = 1, n
         write(1,*) xcdf(i),ycdf(i)
      end do
      close(1)
    end block

    ! residuals
    block
      integer :: i

      open(1,file='/tmp/rl')
      do i = 1, n
!         write(1,*) x(i), (y(i) - fac*x(i))/dres(i), (y(i) - fac*x(i))
         write(1,*) x(i), (cts(i) - fac*pht(i))/sqrt(cts(i)+fac**2*pht(i)), &
              (y(i) - fac*x(i))/sqrt(y(i)+fac**2*x(i))
      end do
      close(1)
    end block
    stop 0


    ! dispersion
    open(33,file='/tmp/r')
    do l = 1, 100
       r =  5e3 * (l / 100.0) * fac**2
       dres = sqrt(y + x*fac**2 + fac**2*r**2)
       dres = sqrt(y + x*fac**2 + r**2)
       !       write(33,*) r,sqrt(sum(res**2/dres**2)/n)
       write(33,*) r, &
            sqrt(sum(huber(res/dres)**2)/sum(dhuber(res/dres))**2*n),&
            sum(ihuber(res/dres))/n

    end do
    close(33)

    ! estimate of s
    s = sqrt(sum(res**2/(y + x*fac**2)/n))
    r = (s**2 - 1)*sum(y + fac**2*x)/n/fac**2
    if( r > 0 )then
       r = sqrt(r)
    else
       r = 0
    end if
    e = sqrt(sum(res**2/(y + x*fac**2 +1.4**2*fac**2)/n))
    rmax = maxval(abs(res)) / fac
    dres = sqrt(y + x*fac**2 + 1.4**2*fac**2)
    write(*,*) 'estimate s,r,rmax, ropt ',s, r, rmax, e, &
         sqrt(sum((res/dres/e)**2)/n)

    res = y - fac*x
    dres = sqrt(y*(1 + fac) + fac**2*sigx**2 + fac**2*sigy**2)
    write(*,*) 'Study of sums:',sum(huber(res/dres)*x/dres)/n, &
         sum((huber(res/dres)*res/dres + 1)*y/dres**2)/n/2

    open(33,file='/tmp/s')
    do l = 0, 100
       r =  rmax * (l / 100.0)
       dres = sqrt(y + x*fac**2 + fac**2*r**2)
       !       s = sqrt(sum(res**2/dres**2)/n)
!       s = 2.54
       !       write(33,*) r,sqrt(sum(res**2/dres**2)/n)
       write(33,*) r,s,sum(ihuber(res/dres/s))/n+0*sum(log(dres*s))/n, &
            sum((res/dres)**2)/n/4+sum(log(dres))/n,&
            sum((res/dres)*res/dres**2)/n-sum(1/dres)/n

    end do
    close(33)

    ! dependence s on r
!    open(33,file='/tmp/sr')
!    do l = 0, 100
!       r =  rmax * (l / 100.0)
!       dres = sqrt(y + x*fac**2 + fac**2*r**2)
!       write(33,*) r,sum(huber(res/dres))/n-0*sum(log(dres*s))/n, &
!            sum((res/dres)**2)/n/4+sum(log(dres))/n,&
!            sum((res/dres)*res/dres**2)/n-sum(1/dres)/n!
!    end do
!    close(33)



    ! huber's walley
    open(33,file='/tmp/w')
    do k = 9500, 10500, 10
       t = fac * (k / 10000.0)
       res = y - t*x
       m = 0
       do l = 0, 1000, 10
          r =  fac * 5e3 * (l / 1000.0)
          dres = sqrt(y + x*t**2 + r**2)
!          dres = sqrt(y + x*t**2 + t**2*r**2)
!          dres = sqrt(y*(1 + t) + r**2)
!          dres = sqrt(2*y + t**2*r**2)
          s = sqrt(sum(res**2/dres**2)/n)
!          call iscale(res/dres,s)
!          write(33,*) t,s,sqrt(sum((res/dres)**2)/n)
!          m = 1
          e = sum(dhuber(res/dres/s))
!          write(*,*) t,s,e
!          e = 1
          if( e > 0 ) then
             sig = sqrt(sum(huber(res/dres/s)**2)/e**2*n)
             !             s = sqrt(sum(res**2/dres**2)/n)
             if( s < 5 ) then
                write(33,*) t,r,sum(ihuber(res/dres))/n, &
                     sum(itukey(res/dres))/n,sig,s
!                write(33,*) t,s,sum(ihuber(res/dres/s))/n,sig,r, &
!                     sum(itukey(res/dres/s))/n
!                write(33,*) t,s,sum((res/dres/s)**2/2)/n,sig
                m = m + 1
             end if
          else
             s = 0
          end if
       end do
       if( m > 0 ) write(33,*)
    end do
    close(33)

    open(33,file='/tmp/w2')
    do k = 900, 1100, 1
       t = fac * (k / 1000.0)
       res = y - t*x
       m = 0
       do l = 0, 1000, 5
          r =  5e3 * (l / 1000.0)
          dres = sqrt(y + t**2*r**2)
          s = sqrt(sum(res**2/dres**2)/n)
          e = sum(dhuber(res/dres/s))
          if( e > 0 ) then
             sig = sqrt(sum(huber(res/dres/s)**2)/e**2*n)
             if( s < 5 ) then
                write(33,*) t,s,sum(ihuber(res/dres/s))/n,sig,r
                m = m + 1
             end if
          else
             s = 0
          end if
       end do
       if( m > 0 ) write(33,*)
    end do
    close(33)


    ! huber's walley, LS case
    open(33,file='/tmp/w1')
    do k = 900, 1100, 1
       t = fac * (k / 1000.0)
       res = y - t*x
       m = 0
       do l = 5, 500, 5
          s =  (l / 100.0)
          dres = sqrt(dy**2 + dx**2*t**2)
!          s = sqrt(sum(res**2/dres**2)/n)
!          call iscale(res/dres,s)
!          write(33,*) t,s,sqrt(sum((res/dres)**2)/n)
!          m = 1
          e = sum(dhuber(res/dres/s))
!          write(*,*) t,s,e
!          e = 1
          if( e > 0 ) then
             sig = sqrt(sum(huber(res/dres/s)**2)/e**2*n)
!             s = sqrt(sum(res**2/dres**2)/n)
             e = sum(ihuber(res/dres/s))/n
             if( s < 5 .and. e < 3 ) then
                write(33,*) t,s,e,sig
!                write(33,*) t,s,sum((res/dres/s)**2/2)/n,sig
                m = m + 1
             end if
          else
             s = 0
          end if
       end do
       if( m > 0 ) write(33,*)
    end do
    close(33)

  end subroutine poisson_mix

  subroutine poisson_mix_contamined

    real(rp), dimension(:), allocatable :: x,u,y,v,res,dres,cts,pht,sig,psi,tt
    type(suite) :: std, rob, qt, stdw, robw, poi
    logical :: cond

    integer :: m,n,k,l,mg,mm
    real(rp) :: mean, mag, fac, cfac, sigx, sigy, zero, prob, &
         p, q, s, r, e, t, rmax, sigx1, sigy1, sigu, sigv

    mean = 1e4
    fac = 2
    cfac = 2.0
    zero = 0.0
    sigx = 2e3
    sigy = 1e3
    sigx1 = 1e2*sigx
    sigy1 = 1e2*sigy
    q = 1

    write(*,*) 'Expected r:',sqrt(sigx**2*fac**2 + sigy**2)

    n = 13
    prob = 1.9
!    m = int(0.9*n)
!    m = n
    allocate(x(n),u(n),y(n),v(n),res(n),dres(n),cts(n),pht(n),sig(n),psi(n))

    do l = 1, n
       mag = 10 + 0.05*l
       mean = 10**((25 - mag)/2.5)

       do
          call random_number(p)
          pht(l) = pnoise(mean)
          if( p < prob ) then
             x(l) = pht(l) + gnoise(zero,sigx)
          else
             x(l) = pht(l) + gnoise(zero,sigx1)
          end if

          call random_number(p)
          if( p < prob ) then
             cts(l) = pnoise(fac*mean)
             y(l) = cts(l) + gnoise(zero,sigy)
          else
             cts(l) = pnoise(cfac*mean)
             y(l) = cts(l) + gnoise(zero,sigy1)
          end if

          if( x(l) > 0 .and. y(l) > 0 ) exit
       end do

    end do

!!$    n = 0
!!$    do l = 1, 2*size(x)
!!$       mean = 1.1*mean
!!$       if( n == size(x) .or. n == m ) exit
!!$       n = n + 1
!!$       pht(n) = pnoise(mean)
!!$       cts(n) = pnoise(fac*mean)
!!$       u(n) = gnoise(zero,sigy)
!!$       v(n) = gnoise(zero,sigx)
!!$!       x(n) = pnoise(mean) + q*gnoise(zero,sigx)
!!$!       y(n) = pnoise(fac*mean) + q*gnoise(zero,sigy)
!!$       x(n) = pht(n) + v(n)
!!$       y(n) = cts(n) + u(n)
!!$
!!$!       if( .not. (x(n) > 0 .and. y(n) > 0) ) write(*,*) 'under zero',n
!!$!       if( .not. (x(n) > 0 .and. y(n) > 0) ) n = n - 1
!!$    end do
!!$    mean = 1e5
!!$    do l = 1, 2*size(x)
!!$       mean = 1.5*mean
!!$       if( n == size(x) ) exit
!!$       n = n + 1
!!$       pht(n) = pnoise(mean)
!!$       cts(n) = pnoise(cfac*mean)
!!$       u(n) = gnoise(zero,sigy1)
!!$       v(n) = gnoise(zero,sigx1)
!!$       x(n) = pht(n) + v(n)
!!$       y(n) = cts(n) + u(n)
!!$!       x(n) = pnoise(mean) + q*gnoise(zero,sigx1)
!!$!       y(n) = pnoise(cfac*mean) + q*gnoise(zero,sigy1)
!!$!       x(n) = gnoise(mean,sigx1)
!!$!       y(n) = gnoise(fac*mean,sigy1)
!!$!       if( .not. (x(n) > 0 .and. y(n) > 0) ) write(*,*) 'Under zero',n
!!$!       if( .not. (x(n) > 0 .and. y(n) > 0) ) n = n - 1
!!$    end do
!!$    if( n /= size(x) ) stop 'n!'

!    y(1) = 1e6
!    y(2) = 1e-6

    call phrate(y,x,t,robustness=.true.,verbose=.true.)

    block
      character(len=66666) :: line
      open(1,file='/tmp/phl_1')
      open(2,file='/tmp/phr')
      do
         read(1,'(a)',end=99) line
         write(2,'(a)') line
      end do
99    continue
      close(1)
      close(2)
!      call rename('/tmp/phl_1','/tmp/pnr')
    end block
    call phrate(y,x,t,robustness=.false.,verbose=.true.)

    stop 0

!    call histogram(u,-5*sigx,5*sigx,1000,'/tmp/u')
!    call histogram(v,-5*sigy,5*sigy,1000,'/tmp/v')
    e = 1
    call histogram(y/x,fac-3*e,fac+5*e,100,'/tmp/r2')

!    e = 3
!    r = 2*sqrt(sigx**2*fac**2 + sigy**2)
!    r = 2226
!    sig = sqrt(cts + fac**2*pht + u**2 + fac**2*v**2)
!    call histogram((y-fac*x)/sig,-e,e,1000,'/tmp/hh2')
!    call histogram((y-fac*x)/sqrt(sig**2 + r**2),-e,e,1000,'/tmp/ff2')
!    call histogram((cts-fac*pht)/sig,-e,e,1000,'/tmp/gg2')

    ! means
    sig = sqrt(cts + fac**2*pht + u**2 + fac**2*v**2)
    write(*,*) 'rate:',sum(y/x)/n
    write(*,*) 'mean rate:',sum(y*x)/sum(x**2)
!    write(*,*) sum(u/sigx - v/sigy)/n
    write(*,*) sum(u - fac*v)/n
    write(*,*) '1: ',sum((u - fac*v))/n
    write(*,*) sum((cts - fac*pht)/sig)/n
    write(*,*) '<>',sum(v)/n,sum(u)/n,sum(v**2/pht**2)/n

    write(*,*) sum(cts/pht)/n, sum((cts*v)/pht**2)/n, sum(u/pht)/n, &
         sum((v*u)/pht**2)/n/fac,sum((v*u)/x**2)/n/fac

    write(*,*) 'E(y)*E(1/x):',sum(y)/n*sum(1/x)/n
    write(*,*) 'E(cts)*E(1/pht):',sum(cts)/n*sum(1/pht)/n
    write(*,*) 'E(1/pht):',sum(1/cts)/n,sum(sigx**2/cts**2)/n
    t = sum(y/x)/n
    t = fac
    r = sigy**2 + t**2*sigx**2
    res = y - t*x
    psi = res/sqrt(cts + t**2*pht + u**2 + t**2*v**2)
    write(*,*) 'E(z):',sum((y-t*x)/sqrt(cts+t**2*pht+r))/n
    write(*,*) 'dz/dy:',sum(1/sqrt(cts+t**2*pht+r))/n
    write(*,*) "r.r:'",4*sqrt(r)*sigx**2*sum((y-t*x)**2/(cts+t**2*pht+r))/n
    write(*,*) 'mix:',sum(psi*(y*v)/sig**3),sum(psi*(x*v)/sig**3), &
         sum(psi*(res*v)/sig**3),sum(v/sig**2),sum(psi*x/sig)
    write(*,*) 'zfun:',sum(psi*x/sig),sum(psi*res*x/sig**3),sum(x/sig**2)
    write(*,*) 'zfun:',sum(psi*pht/sig),sum(psi*res*pht/sig**3),sum(pht/sig**2)
    res = cts - t*pht
    psi = res/sqrt(cts + t**2*pht + u**2 + t**2*v**2)
    write(*,*) sum(psi*u/sig)/sum(psi*pht/sig)
    write(*,*) sum((y*u)/sig**2),sum((x*u)/sig**2),sum(psi*u/sig)
    write(*,*) sum(psi*pht/sig), sum(psi*x/sig), sum(u**2/sig**2), &
         sum(pht*(u-t*v)/sig**2)

    write(*,*) (sum(cts*pht/sig**2) + sum(u**2/sig**2))/sum(pht**2/sig**2)

    tt = y / x
    block
      real(rp) :: t,mu,tau,r,fi,sig,beta,dt,drho,dsig,dbeta,&
           tmin,tmax,thuge,ttiny,fmin,fmax,x,b,pmin,pmax, &
           psi,dpsi,dpsi1,psi1,psimax,psimin,dpsimax,dpsimin,s,ds,psi0,dpsi0,&
           s0,f, bmin, bmax
      real(real64), parameter :: qmin = 0.01, qmax = 0.99, &
           qtiny = 0.001, qhuge = 0.995
      integer :: i,j,nn
      logical :: flag

      mu = mean
      r = sqrt(sigx**2*fac**2 + sigy**2)
      fi = atan2(sigy/2,sigx)
      write(*,*) 'Params:',real(mu),real(fac),real(r),real(fi*RAD)

      call qbisection(qPhi,mu,fac,r,fi,qmin,0.1_REAL64,2.0_REAL64,tmin,flag)
      if( .not. flag ) stop 1
      call qbisection(qPhi,mu,fac,r,fi,qmax,2.0_REAL64,6.0_REAL64,tmax,flag)
      if( .not. flag ) stop 2
      call qbisection(qPhi,mu,fac,r,fi,qtiny,0.1_REAL64,2.0_REAL64,ttiny,flag)
      if( .not. flag ) stop 3
      call qbisection(qPhi,mu,fac,r,fi,qhuge,2.0_REAL64,9.0_REAL64,thuge,flag)
      if( .not. flag ) stop 4

      call rfunction(tmin,mu,fac,r,fi,fmin,psimin,dpsi)
      call rfunction(tmax,mu,fac,r,fi,fmax,psimax,dpsi)

      b = 5
      pmax = 1/tanh(b/2)
      pmin = 1/tanh(b/2)

      write(*,*) 'tmin = ',real(tmin),' tmax=',real(tmax)
      write(*,*) 'ttiny = ',real(ttiny),' thuge=',real(thuge)
      write(*,*) 'fmin = ',real(fmin),' fmax=',real(fmax)
      write(*,*) 'pmin = ',real(pmin),' pmax=',real(pmax)
      write(*,*) 'psimin = ',real(psimin), ' psimax=',real(psimax)

      tau = fac
      call rfunction_tau(tmin,mu,tau,r,fi,psimin,dpsi)
      call rfunction_tau(tmax,mu,tau,r,fi,psimax,dpsi)
      write(*,*) 'psi(min) = ',real(psimin),' psi(max) = ',real(psimax)


!      tau = fac
!      open(2,file='/tmp/rfun_t')
!      do i = 1, size(tt)
!         call rfunction_tau(tt(i),mu,tau,r,fi,psi,dpsi)
!         write(2,*) tt(i),psi
!      end do
!      close(2)

      open(2,file='/tmp/rfun_tau')
      dt = 1e-3
      psi1 = 0
      do j = 50,6000,10
         tau = j*dt

         s0 = 0
         s = 0
         ds = 0
         do i = 1, size(tt)
            if( tt(i) <= ttiny .or. tt(i) >= thuge) then
               psi = 0
               dpsi = 0
!               psi = psimin
            else if( tt(i) < tmin ) then
               x = b*(ttiny - tt(i))/(ttiny - tmin)/2
               psi = psimin*pmin*tanh(x)
               dpsi = psimin*pmin / cosh(x)**2 * (-b)/(ttiny - tmin)/2
!               psi = psimin
            else if( tt(i) > tmax ) then
               x = b*(thuge - tt(i))/(thuge - tmax)/2
               psi = psimax*pmax*tanh(x)
               dpsi = psimax*pmax / cosh(x)**2 * (-b)/(thuge - tmax)/2
!               psi = psimax
!               nplus = nplus + 1
!               write(*,*) tt(i),psi
            else
               call rfunction_tau(tt(i),mu,tau,r,fi,psi,dpsi)
            end if
            call rfunction_tau(tt(i),mu,tau,r,fi,psi0,dpsi)
            dpsi1 = (psi - psi1)/dt
            s0 = s0 + psi0
            s = s + psi
            ds = ds + dpsi
!            if( 1.9 < tt(i) .and. tt(i) < 2 ) write(*,*) tt(i),psi
         end do
         write(2,*) tau,s/n,ds/n,s0/n

         psi1 = psi
      end do
      close(2)

      tau = fac
      call rfunction_r(tmin,mu,tau,r,fi,psimin,dpsi)
      call rfunction_r(tmax,mu,tau,r,fi,psimax,dpsi)
      write(*,*) 'psi(min) = ',real(psimin),' psi(max) = ',real(psimax)

      open(2,file='/tmp/rfun_r')
      psi1 = 0
      do j = 0,9000,10
         r = j

         s0 = 0
         s = 0
         ds = 0
         do i = 1, size(tt)
            if( tt(i) <= ttiny .or. tt(i) >= thuge) then
               psi = 0
               dpsi = 0
!               psi = psimin
            else if( tt(i) < tmin ) then
               x = b*(ttiny - tt(i))/(ttiny - tmin)/2
               psi = psimin*pmin*tanh(x)
               dpsi = psimin*pmin / cosh(x)**2 * (-b)/(ttiny - tmin)/2
!               psi = psimin
            else if( tt(i) > tmax ) then
               x = b*(thuge - tt(i))/(thuge - tmax)/2
               psi = psimax*pmax*tanh(x)
               dpsi = psimax*pmax / cosh(x)**2 * (-b)/(thuge - tmax)/2
!               psi = psimax
!               nplus = nplus + 1
!               write(*,*) tt(i),psi
            else
               call rfunction_r(tt(i),mu,tau,r,fi,psi,dpsi)
            end if
            call rfunction_r(tt(i),mu,tau,r,fi,psi0,dpsi)
            dpsi1 = (psi - psi1)/1.0
            s0 = s0 + psi0
            s = s + psi
            ds = ds + dpsi
!            if( 1.9 < tt(i) .and. tt(i) < 2 ) write(*,*) tt(i),psi
         end do
         write(2,*) r,s/n,ds/n,s0/n

         psi1 = psi
      end do
      close(2)

      tau = fac
      r = sqrt(sigx**2*fac**2 + sigy**2)

      call rfunction_fi(tmin,mu,tau,r,fi,psimin,dpsi)
      call rfunction_fi(tmax,mu,tau,r,fi,psimax,dpsi)
      write(*,*) 'psi(min) = ',real(psimin),' psi(max) = ',real(psimax)

      open(2,file='/tmp/rfun_fi')
      psi1 = 0
      do j = 0, 90
         fi = j / RAD

         s0 = 0
         s = 0
         ds = 0
         do i = 1, size(tt)
            if( tt(i) <= ttiny .or. tt(i) >= thuge) then
               psi = 0
               dpsi = 0
!               psi = psimin
            else if( tt(i) < tmin ) then
               x = b*(ttiny - tt(i))/(ttiny - tmin)/2
               psi = psimin*pmin*tanh(x)
               dpsi = psimin*pmin / cosh(x)**2 * (-b)/(ttiny - tmin)/2
!               psi = psimin
            else if( tt(i) > tmax ) then
               x = b*(thuge - tt(i))/(thuge - tmax)/2
               psi = psimax*pmax*tanh(x)
               dpsi = psimax*pmax / cosh(x)**2 * (-b)/(thuge - tmax)/2
!               psi = psimax
!               nplus = nplus + 1
!               write(*,*) tt(i),psi
            else
               call rfunction_fi(tt(i),mu,tau,r,fi,psi,dpsi)
            end if
            call rfunction_fi(tt(i),mu,tau,r,fi,psi0,dpsi)
            dpsi1 = (psi - psi1)/1.0
            s0 = s0 + psi0
            s = s + psi
            ds = ds + dpsi
!            if( 1.9 < tt(i) .and. tt(i) < 2 ) write(*,*) tt(i),psi
         end do
         write(2,*) RAD*fi,s/n,ds/n,s0/n

         psi1 = psi
      end do
      close(2)

      open(2,file='/tmp/rfun_x')
      fi = 26.56 /  RAD
      tau = fac
      r = sqrt(sigx**2*fac**2 + sigy**2)
      call rfunction_tau(tmin,mu,tau,r,fi,psimin,dpsimin)
      call rfunction_tau(tmax,mu,tau,r,fi,psimax,dpsimax)
!      call rfunction_tau(tau,mu,tmin,r,fi,psimin,dpsimin)
!      call rfunction_tau(tau,mu,tmax,r,fi,psimax,dpsimax)
      write(*,*) real(tmin),real(tmax)
      write(*,*) real(psimin),real(psimax)
      write(*,*) real(dpsimin),real(dpsimax)

      bmax = abs(dpsimax / psimax)
      bmin = abs(dpsimin / psimin)

      write(*,*) real(bmin), real(bmax)

      psi1 = 0
      do j = 50, 6000, 10
!         tau = j / 1000.0
         t = j / 1e3

         s0 = 0
         s = 0
         ds = 0
         dpsi1 = 0
         nn = 0
         do i = 1,1!size(tt)
!            t = tt(i)
!            call rfunction(t,mu,tau,r,fi,f,psi,dpsi)
            f = 1
            nn = nn + 1
            if( .false. .and. (t <= ttiny .or. t >= thuge) ) then
               psi = 0
               dpsi = 0
!               if( tau <= ttiny ) then
!                  psi = psimin
!               else
!                  psi = psimax
!               end if
            else if( t < tmin ) then
               x = b*(ttiny - t)/(ttiny - tmin)/2
               psi = psimin*pmin*tanh(x)
               x = 3*b*(t - tmin)
               psi = psimin*tanh(x) + psimin
               dpsi = dpsimin/cosh(x)**2
!               dpsi = -psimin*pmin / cosh(x)**2 * (-b)/(ttiny - tmin)/2
!               dpsi = dpsimin/cosh(x)**2 * (-3*b)
               dpsi = -psimin/cosh(x)**2 * (3*b)
               !               psi = psimin
!               x = 3*b*(t-tmin)
!               psi = psimin*exp(x)
               !               dpsi = 3*b*psimin*exp(x)
!               x = 3*b*(t - (ttiny+tmin)/2)
!               psimin = -4.1
!               psi = psimin*tanh(x) + psimin
               !               dpsi = -psimin/cosh(x)**2 * (3*b)
               x = (tmin - t) / (tmin - ttiny)
               psi = psimin*exp(-x**2)
               dpsi = 2*x*psimin*exp(-x**2) / (tmin - ttiny)


               x = bmin * (t - tmin)
               psi = psimin*tanh(x) + psimin
               if( psi <= 0 ) then
                  dpsi = -bmin*psimin/cosh(x)**2
               else
                  psi = 0
                  dpsi = 0
               end if


            else if( t > tmax ) then
               x = b*(thuge - t)/(thuge - tmax)/2
               psi = psimax*pmax*tanh(x)
               x = b*(tmax - t)
  !             psimax = psimax / 2
  !             psi = psimax*tanh(x) + psimax
!               dpsi = -psimax*pmax / cosh(x)**2 * (-b)/(thuge - tmax)/2
!               dpsi = dpsimax/cosh(x)**2
 !              dpsi = -psimax / cosh(x)**2 * (-b)
               !               psi = psimax
!               x = b * (tmax - t)
!               psi = psimax*exp(x)
!               dpsi = b*psimax*exp(x)
!               nplus = nplus + 1
               !               write(*,*) tt(i),psi

 !              x = b*((tmax+thuge)/2 - t)
 !              psimax = 0.75
 !              psi = psimax*tanh(x) + psimax
 !              dpsi = -psimax / cosh(x)**2 * (-b)

!               x = (tmax - t)/(thuge - tmax)
!               psimax = 1.35
!               psi = psimax*exp(-x**2)
!               dpsi = 2*x*psimax*exp(-x**2) /(thuge - tmax)

               x = bmax * (tmax -t)
               psi = psimax*tanh(x) + psimax
               if( psi >= 0 ) then
                  dpsi = -bmax*psimax/cosh(x)**2
               else
                  psi = 0
                  dpsi = 0
               end if
            else
               call rfunction_tau(t,mu,tau,r,fi,psi,dpsi)
!               call rfunction(t,mu,tau,r,fi,f,psi,dpsi)
            end if
            call rfunction_tau(t,mu,tau,r,fi,psi0,dpsi0)
!            dpsi1 = dpsi1 + (psi - psi1)/1e-3
            s0 = s0 + f*psi0
            s = s + f*psi
            ds = ds + dpsi
!            if( 1.9 < tt(i) .and. tt(i) < 2 ) write(*,*) tt(i),psi
         end do
         dpsi1 = (s / nn - psi1)/1e-2
         write(2,*) t,s/nn,ds/nn,s0/nn,dpsi1

         psi1 = s / nn
      end do
      close(2)

    end block


    stop 0

    block
      real(rp), dimension(:), allocatable :: xcdf, ycdf
      integer :: i
      allocate(xcdf(n),ycdf(n))
      call ecdf(y/x,xcdf,ycdf)
      q = 0.5;  write(*,*) q, quantile(q,xcdf,ycdf)
      q = 0.25;  write(*,*) q, quantile(q,xcdf,ycdf)
      q = 0.75;  write(*,*) q, quantile(q,xcdf,ycdf)
      open(1,file='/tmp/rcdf')
      do i = 1, n
         write(1,*) xcdf(i),ycdf(i)
      end do
      close(1)
    end block

!!$    block
!!$      integer :: i
!!$      real(rp) :: dx2, dy2, sig
!!$
!!$      dy2 = sum(y**2)/n
!!$      dx2 = sum(x**2)/n
!!$      sig = sqrt(sqrt(dy2-sigy**2) + fac**2*sqrt(dx2-sigx**2))
!!$      !      write(*,*) dx2,dy2,sig
!!$      r = sigy**2 + fac**2*sigx**2
!!$
!!$      open(1,file='/tmp/rl')
!!$      do i = 1, n
!!$ !        write(1,*) pht(i)/sig(i), cts(i)/sig(i)!  (y(i) - fac*x(i))/sig(i)
!!$         write(1,*) x(i), (cts(i) - fac*pht(i))/sqrt(cts(i)+fac**2*pht(i)), &
!!$              (y(i) - fac*x(i))/sqrt(y(i)**2+fac**2*x(i)**2),&
!!$              (y(i) - fac*x(i))/sqrt(y(i)+fac**2*x(i)),&
!!$              (y(i) - fac*x(i))/sqrt(cts(i)+fac**2*pht(i)+r),&
!!$              (y(i) - fac*x(i))/sig, (cts(i) - fac*pht(i))/sig, y(i)
!!$!         write(1,*) x(i), y(i)/x(i), y(i)/x(i) - 1e4/x(i)
!!$      end do
!!$      close(1)
!!$
!!$      res = y - 1.044920*x
!!$      write(*,*) n/2.0,sqrt(n/1.0),count(res/sig > 0), count((res+0.04*x)/sig > 0)
!!$
!!$    end block


    stop 0

    ! the estimate
    call phrate(y,x,t,verbose=.true.)
!    call phrate(cts,pht,t,verbose=.true.)

    write(*,*) '2: ',-sum((y - fac*x)/sig**2)/sum(1/sig**2)

    write(*,*) sum((y-t*x)**2/sig**2)/n


    sigu = sigy
    sigv = sigx

    ! distribution of residuals
    e = 1
    res = (y - t*x)/sqrt(y + t**2*x + sigu**2 + t**2*sigv**2)
    call histogram(res,-5*e,5*e,40,'/tmp/rr')

    ! distribution of beta
    e = 5e0
    res = (y*x*(1+t) +y*t*sigv**2 + x*sigu**2)/ &
         (y + t**2*x + sigu**2 + t**2*sigv**2)**1.5
    call histogram(res,1d-1,e,40,'/tmp/beta')

    ! distribution function
    e = 1
    call histogram(y/x,t-3*e,t+5*e,40,'/tmp/fdis')


!    write(*,*) sum(u/sig)/n, sum(v/sig)/n

    open(1,file='/tmp/z')
    do l = 990,1010
       t = l / 1000.0
!       sig = sqrt(cts + t**2*pht + u**2 + t**2*v**2)
       r = sqrt(sigx**2*t**2 + sigy**2)
       sig = sqrt(y + t**2*x + r)
       res = y - t*x
       write(1,*) t,sum((res/sig)*x/sig) + sum((res/sig)*t*res*x/sig**3) &
            -t*sum(x/sig**2),sum(res/sig)/n
    end do
    close(1)

!    t = 2
!    res = y - t*x
!    sig = sqrt(y + t**2*x + 0*r**2)
!    call histogram(res/sig,-5d0,5d0,1000,'/tmp/r')


!    write(*,*) sum((y-t*x)/sqrt(y + t**2*x+(sigx**2 + sigy**2)*t**2))/n, &
!         sum(y/sqrt(y + t**2*x+(sigx**2 + sigy**2)*t**2))/n

!    write(*,*) 'add:',t*sigx/sqrt(sigx**2*fac**2 + sigy**2),t*sigy/sqrt(sigx**2*fac**2 + sigy**2)

  end subroutine poisson_mix_contamined

  subroutine poisson_mean

    real(REAL64), dimension(:), allocatable :: x,y,pht,cts
    type(suite) :: std, rob, qt, stdw, robw, poi
    logical :: cond

    integer :: mx,my,n,l,i
    real(rp) :: mu, tau, sigu, sigv, zero, z

    mu = 1e6
    tau = 2
    sigv = 1000
    sigu = 1000
    zero = 0

    write(*,*) 'Expected rho:',real(sqrt(sigv**2*tau**2 + sigu**2))
    if( sigv > 0 .and. sigu > 0 ) &
         write(*,*) 'Expected phi:',real(RAD*atan(sigu/(tau*sigv)))

    mx = 1
    my = 1
    n = 150
    allocate(x(n),y(n),pht(n),cts(n))

    do l = 1, n
!       pht(l) = pnoise(mu)
!       cts(l) = pnoise(tau*mu)
!       x(l) = pht(l) + gnoise(zero,sigv)

       z = 0
       do i = 1, mx
!          z = z + pnoise(mu) + gnoise(zero,sigv)
          z = z + gnoise(mu,sqrt(mu)*25) + gnoise(zero,sigv)*0
       end do
       x(l) = z / mx !+ gnoise(zero,sigv)

       z = 0
       do i = 1, my
          z = z + pnoise(tau*mu) + gnoise(zero,sigu)
       end do
       y(l) = z / my  !+ gnoise(zero,sigu)

    end do

    z = sqrt(mu)/10
    call histogram(x,mu-5*z,mu+5*z,70,'/tmp/zn')


    call phrate(y,x,tau,phtfac=sqrt(real(mx,REAL64)),ctsfac=sqrt(real(my,REAL64)), &
         verbose=.true.)


  end subroutine poisson_mean


  subroutine poisson_histo

    use quantilmean

    real(rp), dimension(:), allocatable :: x,dx,y,dy,res,dres,pht,cts,means, &
         mju, sigma, tt, qnt, val, beta, sums, pmeans, ameans
    type(suite) :: std, rob, qt, stdw, robw, poi
    logical :: cond

    integer :: m, mm,n,k,l,i,mg,ns,ntimes
    real(rp) :: mean, mag, fac, sig, sigx, sigy, zero, q, s, r, e, t, rmax, &
         sigu, sigv, phi, rho

    mean = 1e4
    fac = 2
    zero = 0.0
    sigx = 2e3
    sigy = 2e3
    q = 1

    rho = sqrt(sigx**2*fac**2 + sigy**2)
    write(*,*) 'Expected rho:', real(rho)
    if( sigx > 0 .and. sig > 0 ) then
       phi = atan2(sigy,fac*sigx)
       write(*,*) 'Expected phi:',real(RAD*phi)
    else
       phi = 0
    end if

    ntimes = 6666
    allocate(pmeans(ntimes),ameans(ntimes))

    do mm = 1, ntimes

    n = 13
    allocate(x(n),dx(n),y(n),dy(n),res(n),dres(n),pht(n),cts(n),means(n),sums(2*n))

    n = 0
    do mg = 1000,1600, 100
!    do mg = 1,130000
       mag = mg / 100.0
       mean = 10**((25 - mag)/2.5)
!       mag = 25 - 2.5*log10(mean)
!       m = max(int(exp(mag*1.1)/25.0),1)
  !     m = max(int(exp(mag*1.1)/20e3),1)
!       m = max(int(exp(mag*1.1)/1e3),1)
    m = max(int(exp(mag*1.1)/5e5),1)
!      m = size(x)
!      mean = 1e4
!       m = 1
       write(*,*) mean,mag,m
       ns = 0
       sums = 0
       do l = 1, m
          if( n == size(x) ) exit
          n = n + 1
          means(n) = mean
          pht(n) = pnoise(mean)
          cts(n) = pnoise(fac*mean)
          x(n) = pht(n) + q*gnoise(zero,sigx)
          y(n) = cts(n) + q*gnoise(zero,sigy)
          ns = ns + 2
          sums(ns-1) = x(n)
          sums(ns) = y(n) / fac
!          if( x(n) < 0 .or. y(n) < 0 ) n = n - 1
          e = x(n) + q*gnoise(zero,sigx)**2
          if( e > 0 ) then
             dx(n) = sqrt(e)
          else
             dx(n) = 1e-7
          end if
          e = y(n) + q*fac**2*gnoise(zero,sigy)**2
          if( e > 0 ) then
             dy(n) = sqrt(e)
          else
             dy(n) = 1e-7
          end if
!          if( .not. (x(n) > 0 .and. y(n) > 0) ) n = n - 1
       end do

       mean = sum(sums(1:ns))/ns
       sig = sqrt(sum((sums(1:ns) - mean)**2)/(ns - 1))
!       write(*,*) 'E(mean)=',real(mean),'+-',real(sig/sqrt(1.0*ns)),real(sig)

       if( n == size(x) ) exit
    end do

    if( n /= size(x) ) then
       write(*,*) 'BIG warning: n /= size(x)',n,size(x)
       stop 'BIG'
    end if


    ! estimate
    call phrate(y,x,t,verbose=.true.)

    pmeans(mm) = t
    ameans(mm) = sum(y/x)/n

    deallocate(x,dx,y,dy,res,dres,pht,cts,means,sums)

 end do

 call histogram(pmeans,2-0.05d0,2+0.05d0,40,'/tmp/pm')
 call histogram(ameans,2-0.15d0,2+0.15d0,60,'/tmp/am')



  end subroutine poisson_histo


  subroutine lssolve


    real(rp), dimension(:), allocatable :: x,u,y,v,res,dres,cts,pht,sig,psi

    real(rp) :: mean, fac, sigx, sigy, zero
    real(rp) :: a,b,t,fa,fb,f
    integer :: n, iter

    mean = 1e4
    fac = 2
    zero = 0.0
    sigx = 2e3
    sigy = 2e3

    write(*,*) 'Expected r:',sqrt(sigx**2*fac**2 + sigy**2)

    n = 1300
    allocate(x(n),u(n),y(n),v(n),res(n),dres(n),cts(n),pht(n),sig(n),psi(n))

    do n = 1, size(x)
       pht(n) = pnoise(mean)
       cts(n) = pnoise(fac*mean)
       u(n) = gnoise(zero,sigy)
       v(n) = gnoise(zero,sigx)
       x(n) = pht(n) + v(n)
       y(n) = cts(n) + u(n)
    end do

    a = 0.9*fac
    b = 1.1*fac
    fa = lsfun(a,cts,pht,u,v,x,y)
    fb = lsfun(b,cts,pht,u,v,x,y)

    do iter = 1, 3*precision(t)

       if( abs(b - a) < sqrt(epsilon(t))) exit

       t = (a + b) / 2
       f = lsfun(t,cts,pht,u,v,x,y)
       write(*,*) iter,t,f
       if( f*fa < 0 ) then
          b = t
          fb = f
       else
          a = t
          fa = f
       end if

    end do

  end subroutine lssolve

  function lsfun(t,cts,pht,u,v,x,y)

    real(rp) :: lsfun
    real(rp), intent(in) :: t
    real(rp), dimension(:), intent(in) :: cts,pht,u,v,x,y
    real(rp), dimension(:), allocatable :: res,psi,sig
    integer :: n

    n = size(cts)
    allocate(res(n),psi(n),sig(n))

    res = cts - t*pht
    sig = sqrt(cts + t**2*pht + u**2 + t**2*v**2)
    psi = res/sig
    lsfun = sum(psi*pht/sig) + sum(psi*u/sig) + sum(u**2/sig**2) + &
         sum(pht*(u-t*v)/sig**2)

  end function lsfun

  subroutine lsmin

    real(rp), dimension(:), allocatable :: x,u,y,v,res,cts,pht,sig,sigx2,sigy2,&
         a,b,c,d,sigr, sigx, sigy, tt

    real(rp) :: mean, fac, sigx0, sigy0, zero, mux, muy, tmin, smin
    real(rp) :: t,f,e,ff,r,r0,f0, f1, s, s1, sigma
    integer :: n, i, iter

    mean = 1e4
    fac = 2
    zero = 0.0
    mux = mean
    muy = fac*mean
    sigx0 = 0.2*mean
    sigy0 = 0.2*mean
    r0 = sqrt(sigy0**2 + fac**2*sigx0**2)

    write(*,*) 'Expected r:',r0,' sig:', &
         sqrt(sigy0**2 + fac**2*sigx0**2 + muy + t**2*mux)
    write(*,*) mux,sigx0


    n = 7000
    allocate(x(n),u(n),y(n),v(n),res(n),cts(n),pht(n),sig(n),sigx2(n),sigy2(n),&
         a(n),b(n),c(n),d(n),sigr(n),sigx(n),sigy(n), tt(n))

    do i = 1, size(x)
       pht(i) = pnoise(mux)
       cts(i) = pnoise(muy)
       u(i) = gnoise(zero,sigy0)
       v(i) = gnoise(zero,sigx0)
       x(i) = pht(i) + v(i)
       y(i) = cts(i) + u(i)
    end do
    sigy2 = cts + u**2
    sigx2 = pht + v**2
    sigx = sqrt(sigx2)
    sigy = sqrt(sigy2)
    tt = y / x

!    e = sigy0
!    call histogram(u,-5*e,5*e,100,'/tmp/hu')
!    e = sigx0
!    call histogram(v,-5*e,5*e,100,'/tmp/hv')
    e = sqrt(mux)
    call histogram(pht,mux-5*e,mux+5*e,100,'/tmp/hp')
    e = sqrt(muy)
    call histogram(cts,muy-5*e,muy+5*e,100,'/tmp/hc')
    e = sqrt(mux + sigx0**2)
    call histogram(x,mux-5*e,mux+5*e,100,'/tmp/hx')
    e = sqrt(muy + sigy0**2)
    call histogram(y,muy-5*e,muy+5*e,100,'/tmp/hy')
    e = 1
    call histogram(sigx**2,e,10*mux,100,'/tmp/hsx')
    call histogram(sigy**2,e,10*muy,100,'/tmp/hsy')

    e = 1
    call histogram(y/x,fac-2*e,fac+3*e,100,'/tmp/rh')
    e = 0.01
    call histogram(cts/pht,fac-2*e,fac+3*e,100,'/tmp/rc')
    e = 0.2
    call histogram(log(y/x),0.7-5*e,0.7+5*e,100,'/tmp/logy')
    e = 0.01
    call histogram(log(cts/pht),0.7-5*e,0.7+5*e,100,'/tmp/logc')
    e = 0.2
    call histogram(log(y/x)-log(2.0),-5*e,5*e,100,'/tmp/logd')

    write(*,*) 'E(y), E(c), E(u):', sum(y)/n, sum(cts)/n, sum(u)/n
    write(*,*) 'E(x), E(p), E(v):', sum(x)/n, sum(pht)/n, sum(v)/n
    write(*,*) 'D(y-2*x), D(c-2*n):',sqrt(sum((y-2*x)**2)/n), &
         sqrt(sum((cts-2*pht)**2)/n)
    write(*,*) 'D(u), D(v):', sqrt(sum(u**2)/n), sqrt(sum(v**2)/n)
    write(*,*) 'D(y), D(x):', sqrt(sum(y**2)/n), sqrt(sum(x**2)/n)
    write(*,*) 'D(cts), D(pht):', sqrt(sum(cts**2)/n), sqrt(sum(pht**2)/n)
    !    write(*,*) sqrt(sum((cts-2*pht)**2)/n + sigy0**2 + 2**2*sigx0**2)
    t = fac
    sigma = sqrt(muy + t**2*mux + sigy0**2 + t**2*sigx0**2)
    write(*,*) sum((y*t*sigx0**2+x*sigy0**2)*exp(-((y-t*x)/sigma)**2/2))/sigma**3/n

    write(*,*) sigma, muy*t*sigx0**2+mux*sigy0**2
    write(*,*) sqrt(sum((y-muy)**2)/n),sqrt(sum((x-mux)**2)/n)
    write(*,*) sqrt(sum((cts-muy)**2)/n),sqrt(sum((pht-mux)**2)/n)
    write(*,*) sum((cts-muy)*u)/sum((cts-muy)**2),sum((pht-mux)*v)/sum((pht-mux)**2)

    write(*,*) 'E(cts/pht), E(y/x):',sum(cts/pht)/n, sum(y/x)/n
    write(*,*) sum(u/x)/n, sum(u/cts)/n, sum(v/pht)/n, sum(u*v/pht/cts)/n
    write(*,*) sum(cts*v/pht**2)/n, sum(y*v/x**2)/n
    write(*,*) 'E(y-muy), E(c-muy): ',sum((y-muy)**2)/n,sum((cts-muy)**2)/n

!    write(*,*) sum((y-2*x)**2) -(sigy0**2+t**2*sigx0**2)
!    write(*,*) (sum((y-2*x)**2)/n-sum(cts**2)/n-sigy0**2)/(sigx0**2-sum(pht)/n)


    t = fac
    smin = -huge(smin)
    open(1,file='/tmp/lsapprox')
    open(2,file='/tmp/lst')
    do i = 50,550, 10
       t = i/1e2
       res = y - t*x
       sig = sqrt(sigy2 + t**2*sigx2)
       sigma = sqrt(muy + fac**2*mux + sigy0**2 + fac**2*sigx0**2)/4
!       sigma = sqrt(muy + t**2*mux + sigy0**2 + t**2*sigx0**2)/4
       f = sum(x/sigma*exp(-((y-t*x)/sigma)**2/2))
!       f0 = sum(exp(-((y-t*x)/sigma)**2/2))
       sigma = sqrt(muy + t**2*mux + sigy0**2 + t**2*sigx0**2)!/4
!       f1 = product((x/sigma+t*(y-t*x)*(mux+sigx0**2)/sigma**3)*exp(-((y-t*x)/sigma)**2/2))
!       write(*,*) t,-sum(((y-t*x)/sigma)**2)/2,sum(x/sigma),sum(t*(y-t*x)*(mux+sigx0**2)/sigma**3)
!       s = sum(((y-t*x)/sigma)**2)/2 - sum(log(x/sigma - t*(y-t*x)*(mux+sigx0**2)/sigma**3), &
!            1,x/sigma > t*(y-t*x)*(mux+sigx0**2)/sigma**3)
       s = -sum(((y-t*x)/sigma)**2)/2 + &
            sum(log((y*t*(mux+sigx0**2) + x*(muy+sigy0**2))/sigma**3))
!            sum(log(x/sigma + t*(y-t*x)*(mux+sigx0**2)/sigma**3))
       !            sum(log(max(x/sigma + t*(y-t*x)*(mux+sigx0**2)/sigma**3,1d0)))

       s = -mux**2*sum((t-tt)**2)/2/sigma**2 + &
            sum(log(mux*(t*tt*(mux+sigx0**2)+t*mux+sigy0**2)/sigma**3))

       if( s > smin ) then
          smin = s
          tmin = t
       end if
       s1 = -sum(((cts-t*pht)/sqrt(muy + t**2*mux))**2)/2 + &
!            sum(log(max(pht/sigma - t*(y-t*x)*(mux+0*sigx0**2)/sigma**3,1d0)))
            sum(log((cts*t*mux + pht*muy)/sqrt(muy + t**2*mux)**3))
       write(1,*) t,f,s/n,s1/n
       s = sum((y-t*x)**2/sigma**2)/n
       s1 = sum((cts-t*pht)**2/sqrt(muy + t**2*mux)**2)/n
!       s = sum(log((y*t*(mux+sigx0**2) + x*(muy+sigy0**2))/sigma**3))/n
!       s1 = sum(log((cts*t*mux + pht*muy)/sqrt(muy + t**2*mux)**3))/n
       write(2,*) t,s,s1
    end do
    close(1)
    close(2)
    write(*,*) 't in the extreme:', tmin, smin

    e = 1
    call histogram(x/sigma + t*(y-t*x)*(mux+sigx0**2)/sigma**3,0.8-e,0.8+e,&
         100,'/tmp/xxx')

    open(1,file='/tmp/lsmu')
    do i = 7000, 13000, 100
       mux = i
       sig = sqrt(t*mux + tt**2*mux + sigy0**2 + tt**2*sigx0**2)
       s = -sum((mux*(tt-t)/sig)**2)/2 + &
            sum(log(mux*(mux*t*(tt+1)+t*tt*sigx0**2+sigy0**2)/sig**3))
       write(1,*) mux, s
    end do
    close(1)


    return

    open(1,file='/tmp/lsmin')
!    open(2,file='/tmp/hinkley')
!    do i = 19500,20500,2
    do i = 50,700,5
       t = i/1e2
       res = y - t*x
       sig = sqrt(sigy2 + t**2*sigx2)
!       sigr = sqrt(cts + t**2*pht + r**2)
       !       sigr = sqrt(cts + t**2*pht + sigy0**2 + t**2*sigx0**2)
       r = sqrt(sigy0**2 + t**2*sigx0**2)
       sigr = sqrt(y + t**2*x + r**2)
!       sigr = sqrt(y + t**2*x)

       if( .true. ) then
          a = sqrt(t**2/sigy2 + 1/sigx2)
          b = y*t/sigy2 + x/sigx2
          c = y**2/sigy2 + x**2/sigx2
          d = exp((b/a)**2/2 - c/2)

          a = sqrt(t**2/sigy2 + 1/sigx2)
          b = muy*t/sigy2 + mux/sigx2
          c = muy**2/sigy2 + mux**2/sigx2
          d = exp((b/a)**2/2 - c/2)

          f0 = sum((b**2 - c*a**2)/a**2)/2/n + &
               sum(log(b/(sigx*sigy*a**3*2.506628274631_rp)))/n

          f0 = sum((b**2 - c*a**2)/a**2)/2/n + &
               sum(log(b/(sigx*sigy*a**3*2.506628274631_rp)))/n +&
               sum(log(Phi(b/a) - Phi(-b/a)))/n

          f1 = sum((b**2 - c*a**2)/a**2)/2/n + &
               sum(log(b/(sigx*sigy*a**3*2.506628274631_rp)))/n +&
               sqrt(2/3.14159)*sum(a/b*exp(-(b/a)**2/2))/n


!          f0 = sum((b**2 - c*a**2)/a**2)/2/n + &
!               sum(log(1.25331413731550025_rp*b/a))/n

!               -sum(log(sqrt(sigx2)*sqrt(sigy2)*a**2))/n
!          f0 = sum((b/a)**2 - c)/2 - &
!               sum(log(b/a))
!               sum(log(b/(sqrt(sigx2)*sqrt(sigy2)*a**3))) + &
!               sum(log(Phi(b/a) - Phi(-b/a)))
       end if
!       write(2,*) t, &
!            sum(b*d/sqrt(sigx2)/sqrt(sigy2)/a**3/2.506628275)/n, &
!            sum(Phi(b/a) - Phi(-b/a))/n, &
!            sum(exp(-c/2)/3.14159265358979324/sqrt(sigx2)/sqrt(sigy2)/a**2)/n

!       f = -sum((res/sig)**2)/2 + &
!            sum(log((y*t*sigx2 + x*sigy2)/sig**3/2.506628275)) - &
!            sum(log(sqrt(sigx2)*sqrt(sigy2)/sig**2/3.14159))

       f = -sum((res/sig)**2)/2/n + &
            sum(log((y*t*sigx2 + x*sigy2)/sig**3/2.506628274631_rp))/n !+ &
!            sum(log(0.8*sigx*sigy*sig/(y*t*sigx2+x*sigy2)* &
!            exp(-((y*t*sigx2+x*sigy2)/(sig*sigx*sigy))**2/2)))/n

       f = -mux**2*sum((tt - t)**2/(tt**2*sigx2 + sigy2))/2/n + &
            sum(log(mux*(t*tt*sigx2 + sigy2)/(sigy2 + tt**2*sigx2)**1.5))/n

!       f = -sum((res/sig)**2)/2/n + &
!            0*sum(x)/n +  sum(log(sig))/n

       ff = -sum((res/sigr)**2)/2/n + &
!            sum(log(1.25331413731550025_rp*t*sqrt(y*x)/sigr))/n
!            sum(log(1.25331413731550025_rp*(y*x*(1+t) + x*r**2)/sigr))/n
            sum(log((y*x*(1+t) + x*r**2)/sigr**3/2.506628274631_rp))/n*0

       ff = -mux**2*sum((tt - t)**2/(tt**2*sigx2 + sigy2))/2/n - &
            sum(log(sigy2 + tt**2*sigx2))/n/2 + log(mux)

       ff = -mux*sum((tt - t)**2/(tt**2 + t))/2/n + &
            sum(log((t + t*tt)/(tt**2+t)**1.5))/n + log(mux)/2

       !            sum(log((y*t*sigx2 + x*sigy2)/sig**3))/2.506628275 - &
!            sum(log(sqrt(sigx2)*sqrt(sigy2)/sig**2/3.14159))
!            sum(log(1+1.41*sqrt(sigx2)*sqrt(sigy2)*sig/(y*t*sigx2 + x*sigy2)&
!            *exp((res/sig)**2/2)))

!       f = sum(log((y*t*sigx2 + x*sigy2)/sig**3/2.5*exp(-(res/sig)**2/2) * &
!            (Phi(b/a) - Phi(-b/a)) + &
!            sqrt(sigx2)*sqrt(sigy2)/sig**2/3.14*exp(-c/2)))

!       f = sum(log(b*d/sqrt(sigx2)*sqrt(sigy2)/a**3/2.5*(Phi(b/a) - Phi(-b/a)) &
!            + exp(-c/2)/3.14/sqrt(sigx2)*sqrt(sigy2)/a**2))

!            sum(log(Phi(b/a) - Phi(-b/a))) -&
!            sum(y**2/sigy2 + x**2/sigx2)/2 - &

       write(1,*) t,f, ff, f0, f1, &
            -sum((res/sig)**2)/2/n, &
            sum((b**2 - c*a**2)/a**2)/2/n, sum(log(b/(sqrt(sigx2)*sqrt(sigy2)*a**3*2.506628274631_rp)))/n , &
            -sum((res/sig)**2)/2/n, sum(log((y*t*sigx2 + x*sigy2)/sig**3/2.506628274631_rp))/n, &
            sum(log((y*t*sigx2 + x*sigy2)/sig**3))/n, &
            sum(log((y*x*(1+t) + x*r**2)/sigr**3))/n, &
            sum(log((y*t*sigx2 + x*sigy2)/sig**1))/n, &
            sum(log((y*x*(1+t) + x*r**2)/sigr**1))/n, &
            sum((res/sig)**2)/n-sum((res**2)/sigr**2)/n, &
            sum((3.14*sqrt(sigx2)*sqrt(sigy2)*a**2*exp(-c/2)))/n, &
            sum((b/(sqrt(sigx2)*sqrt(sigy2)*a**3*2.506628274631_rp)))/n, &
            sum((y*t*sigx2 + x*sigy2)/(sqrt(sigx2)*sqrt(sigy2)*a*2.506628274631_rp))/n, &
            -sum(c)/2/n - sum(log(3.14159*sqrt(sigx2)*sqrt(sigy2)*a**2))/n, &
            (-sum(c)/2-sum(log(3.14159*sqrt(sigx2)*sqrt(sigy2)*a**2))) * &
            sum(log(b/(sqrt(sigx2)*sqrt(sigy2)*a**3*2.506628274631_rp))), &
            -mux**2*sum((tt - t)**2/(tt**2*sigx2 + sigy2))/2/n,             &
            sum(log(mux*(t*tt*sigx2 + sigy2)/(sigy2 + tt**2*sigx2)**1.5))/n
!            sum(log(Phi(b/a) - Phi(-b/a)))/n


!            sum(log((y*t*sigx2)/sig**1.5))/n,sum(log((x*sigy2)/sig**1.5))/n
!            sum(log((t*sigx2*(y-t*x)*0 + x*y + x**2*t**2+  x*r**2)/sigr**1.5)),&
!            sum(sig**2)/n,sum(sigr**2)/n
!       write(1,*) t,f/n,sum((res/sig)**2)/2/n, &
!            sum(log((y*t*sigx2 + x*sigy2)/sig**1.5))/n, &
!            ff/n
!       write(1,*) t,f/n,sum(Phi(b/a) - Phi(-b/a)), &
!            sum(exp(-c/2)/3.14/sqrt(sigx2)*sqrt(sigy2)/a**2), &
!            sum(b/sqrt(sigx2)*sqrt(sigy2)/a**3), &
!            sum(d)
    end do
    close(1)
!    close(2)

    open(1,file='/tmp/r')
!    do i = 19500,20500,2
    do i = 50,500
       t = i/1e2
       res = y - t*x
       sig = sqrt(sigy2 + t**2*sigx2)
       r = sqrt(sigy0**2 + t**2*sigx0**2)
       sigr = sqrt(y + t**2*x + r**2)

!       write(1,*) r,-sum((res/sig)**2)/2, -sum((res**2+0*r**2)/sigr**2)/2
!       write(1,*) r,sum(((y*t*sigx2 + x*sigy2)/sig**1.5))/n, &
!            sum(((y*x*(1+t) + x*r**2)/sigr**1.5))/n
       write(1,*) r, &
            sum(((y*t*sigx2 + x*sigy2)))/n, &
            sum(((y*x*(1+t) + x*r**2)))/n, &
            sum(((y*t*sigx2 + x*sigy2)/sig**3))/n, &
            sum(((y*x*(1+t) + x*r**2)/sigr**3))/n

    end do
    close(1)

  end subroutine lsmin

  subroutine lsmin2d

    real(rp), parameter :: rad = 45.0 / atan(1.0)
    real(rp), dimension(:), allocatable :: x,u,y,v,cts,pht,tt,sigx,sigy,sig
    real(rp) :: mean, fac, sigv, sigu, zero, mux, muy, t, mu, f, r, q, f1, phi
    logical, dimension(:), allocatable :: mask
    integer :: n, i, j

    mean = 1e4
    fac = 2
    zero = 0.0
    mux = mean
    muy = fac*mean
    sigu = 0.2*mean
    sigv = 0.2*mean

    n = 70000
    allocate(x(n),u(n),y(n),v(n),cts(n),pht(n),tt(n),sigx(n),sigy(n),sig(n),&
         mask(n))

    do i = 1, size(x)
       pht(i) = pnoise(mux)
       cts(i) = pnoise(muy)
       u(i) = gnoise(zero,sigu)
       v(i) = gnoise(zero,sigv)
       x(i) = pht(i) + v(i)
       y(i) = cts(i) + u(i)
    end do
    tt = y / x

    open(1,file='/tmp/lsmin2d')
    do i = 500, 5000, 10
       t = i / 1e3
       do j = 380, 430
          mu = 10**(j/100.0)
          sigx = sqrt(mu + sigv**2)
          sigy = sqrt(mu*t + sigu**2)
          sig = sqrt(tt**2*sigx**2 + sigy**2)
          f = sum((mu*(t - tt)/sig)**2,tt>0)/2/n - &
               sum(log(mu*(t*tt*sigx**2 + sigy**2)/sig**3),tt>0)/n
          write(1,*) t,mu,f
       end do
       write(1,*)
    end do
    close(1)

    open(1,file='/tmp/lsmin2d_uv')
    t = 2
    mu = 1e4
    do i = -9000, 9000, 100
       sigu = i
       do j = -9000, 9000, 100
          sigv = j
          r = sum(sqrt(tt**2*sigv**2 + sigu**2))/n
          sigx = sqrt(mu + sigv**2)
          sigy = sqrt(mu*t + sigu**2)
          sig = sqrt(tt**2*sigx**2 + sigy**2)
          f = sum((mu*(t - tt)/sig)**2,tt>0)/2/n - &
               sum(log(mu*(t*tt*sigx**2 + sigy**2)/sig**3),tt>0)/n
          write(1,*) sigu,sigv,f
       end do
       write(1,*)
    end do
    close(1)

    open(1,file='/tmp/lsmin2d_tu')
    mu = 1e4
    sigv = 2e3
    do i = 1500, -3500, 10
       t = i/1e3
       do j = 100, 3000, 100
          sigu = j
          r = sum(sqrt(tt**2*sigv**2 + sigu**2))/n
          sigx = sqrt(mu + sigv**2)
          sigy = sqrt(mu*t + sigu**2)
          sig = sqrt(tt**2*sigx**2 + sigy**2)
          f = sum((mu*(t - tt)/sig)**2,tt>0)/2/n - &
               sum(log(mu*(t*tt*sigx**2 + sigy**2)/sig**3),tt>0)/n
          write(1,*) t,sigu,f
       end do
       write(1,*)
    end do
    close(1)

    open(1,file='/tmp/lsmin2d_tv')
    mu = 1e4
    sigu = 2e3
    do i = 1500, -3500, 10
       t = i/1e3
       do j = 100, 3000, 100
          sigv = j
          r = sum(sqrt(tt**2*sigv**2 + sigu**2))/n
          sigx = sqrt(mu + sigv**2)
          sigy = sqrt(mu*t + sigu**2)
          sig = sqrt(tt**2*sigx**2 + sigy**2)
          f = sum((mu*(t - tt)/sig)**2,tt>0)/2/n - &
               sum(log(mu*(t*tt*sigx**2 + sigy**2)/sig**3),tt>0)/n
          write(1,*) t,sigv,f
       end do
       write(1,*)
    end do
    close(1)


    open(1,file='/tmp/lsmin2d_rf')
    t = 2
    x = 1e4
    do i = 2000, 8000, 100
       r = i
       do j = -180, 180, 10
          phi = j / rad
          sigu = r*sin(phi)
          sigv = r*cos(phi)/t
!          sigx = sqrt(mu + s)
          !r = sum(sqrt(tt**2*sigv**2 + sigu**2))/n
          sigx = sqrt(mu + sigv**2)
          sigy = sqrt(mu*t + sigu**2)
          sig = sqrt(tt**2*sigx**2 + sigy**2)
!          sig = sqrt(tt**2*sigx**2 + sigy**2)
          f = sum((mu*(t - tt)/sig)**2,tt>0)/2/n - &
               sum(log(mu*(t*tt*sigx**2 + sigy**2)/sig**3),tt>0)/n
!          phi = rad*atan2(t*sigv,sigu)
          write(1,*) r,rad*phi,f


!!$!          sigu = r*cos(phi)
!!$          !          sigv = r/t*sin(phi)
!!$          r = s
!!$          sig = sqrt(x*(t + tt**2) + r**2)
!!$!          write(*,*) r,q,sig(1),x(1)*t*(tt(1)+1) + r**2*(1+tt(1)*(t-tt(1))*q**2)
!!$          mask = tt > 0 .and. (x*t*(tt+1) + r**2*(t*tt*sigv**2 + sigu**2) ) > 0
!!$          f = sum((x*(t - tt)/sig)**2,mask)/2/n - &
!!$               sum(log(x*(x*t*(tt+1) + r**2*(t*tt*sigv**2 + sigu**2))/sig**3),mask)/n
!!$          f1 = sum((x*(t - tt)/sig)**2,mask)/2/n - &
!!$               sum(log(x*(x*t*(tt+1) + r**2)/sig**3),mask)/n
!!$          write(1,*) r,phi,f,f1
       end do
       write(1,*)
    end do
    close(1)

  end subroutine lsmin2d


  elemental pure function Phi(x)

    real(rp) :: Phi
    real(rp), intent(in) :: x
    real(rp), parameter :: s2 = sqrt(2.0_rp)

    Phi = (1.0_rp + erf(x/s2)) / 2.0_rp

  end function Phi

  subroutine InvPhi

    real(REAL64) :: q,x,f,xmin,xmax,fmin,fmax
    integer :: i

    q = 3.0_REAL64 / 4.0_REAL64
    xmin = 0
    fmin = Phi(xmin) - q
    xmax = 3
    fmax = Phi(xmax) - q
    do i = 1, 5*precision(q)
       x = (xmax + xmin) / 2.0_REAL64
       f = Phi(x) - q
       if( f*fmin < 0 ) then
          xmax = x
          fmax = f
       else if( f*fmax < 0 ) then
          xmin = x
          fmin = f
       else
          write(*,*) '!',x,f
       end if
       if( abs(fmax - fmin) < 10*epsilon(f) ) exit
       write(*,*) i,xmin,xmax
    end do

    write(*,*) q,x,f,Phi(x)

  end subroutine InvPhi



  function phfun(t)

    real(rp) :: phfun
    real(rp), intent(in) :: t
    real(rp) :: mux, muy, sigx, sigy, sig, a,b,c,d
    real(rp), parameter :: pi = 3.14159265358979324_rp
    real(rp), parameter :: pi2 = 2.5066282746310005_rp
    real(rp), parameter :: pii2 = 1.25331413731550025_rp

    mux = 1e4
    muy = 2e4
    sigx = sqrt(0.2e4**2*0 + mux)
    sigy = sqrt(0.2e4**2*0 + muy)

    a = sqrt(t**2/sigy**2 + 1/sigx**2)
    b = muy*t/sigy**2 + mux/sigx**2
    c = muy**2/sigy**2 + mux**2/sigx**2
    d = exp((b**2 - c*a**2)/a**2/2)

!    phfun = b*d/pi2/sigx/sigy/a**3*(Phi(b/a) - Phi(-b/a)) + &
!         exp(-c/2)/pi/sigx/sigy/a**2
!    phfun = b*exp((b/a)**2/2)/pi2/sigx/sigy/a**3*(Phi(b/a) - Phi(-b/a)) + &
!         1/pi/sigx/sigy/a**2
!    phfun = exp(-c/2)*phfun

    phfun = (t*muy*sigx**2 + mux*sigy**2)/(t**2*sigx**2 + sigy**2)**1.5* &
         exp(((b/a)**2-c)/2)/pi2

!    phfun = (muy - t*mux)**2/(t**2*sigx**2 + sigy**2)

!    phfun = (b**2 - c*a**2)/(2*a**2) + log(b/(pi2*sigx*sigy*a**3)

    !    write(*,'(g0,4es16.7)') real(t),a,b,c,d
!    write(*,'(g0,4es16.7)') real(t),&
!         exp(((b/a)**2-c)/2),exp(-c/2), &
!         b/pi2/sigx/sigy/a**3*(Phi(b/a) - Phi(-b/a)), &
!         1/pi/sigx/sigy/a**2

!    sig = sqrt(t**2*sigx**2 + sigy**2)
!    phfun = (muy*t*sigx**2 + mux*sigy**2)/sig**3/pi2 &
!         *exp(-((muy - t*mux)/sig)**2/2) +(sigx*sigy)/sig**2/pi*exp(-c/2)

!    sig = sqrt(t**2*sigx**2 + sigy**2)
!    phfun = (muy*t*sigx**2+mux*sigy**2)/sig**3/pi2*exp(-((muy-t*mux)/sig)**2/2)
!    phfun = t*phfun

!    phfun = ((b/a)**2-c)/2 + log(pii2*(b/a))
!    phfun = exp(phfun)!/17.78979489158902

!    write(*,*) t,(muy-t*mux)/sig,(muy*t*sigx**2+mux*sigy**2)/sig**3

  end function phfun

  subroutine tabelate

    real(rp) :: mux,muy,sigx,sigy,sigu,sigv,sig, a,b,c,d,f,t,g,e, &
         w,f1,eps,beta,dt,df,r,dsig, df1
    real(rp), parameter :: pi = 3.14159265358979324_rp
    real(rp), parameter :: pi2 = 2.5066282746310005_rp
    integer :: i

    mux = 1e4
    muy = 2*mux
    sigu = 0.05*mux
    sigv = 0.05*mux
    sigx = sqrt(mux + sigv**2)
    sigy = sqrt(muy + sigu**2)

    open(2,file='/tmp/hinkley_405')
    do i = 50,5000
       t = i/1e3

       sig = sqrt(t**2*sigx**2 + sigy**2)
       a = sqrt(t**2/sigy**2 + 1/sigx**2)
       b = muy*t/sigy**2 + mux/sigx**2
       c = muy**2/sigy**2 + mux**2/sigx**2
       d = exp((b**2 - c*a**2)/a**2/2)
       f = b*d/pi2/sigx/sigy/a**3*(Phi(b/a) - Phi(-b/a)) + &
            exp(-c/2)/pi/sigx/sigy/a**2
       g = b*d/pi2/sigx/sigy/a**3
!       e = (muy*t*sigx**2 + mux*sigy**2)*exp(-((muy-t*mux)/sig)**2/2)/sig**3/pi2
       e = exp(-((muy-t*mux)/sig)**2/2)
       w = (muy*t*sigx**2+mux*sigy**2)/sig**3/pi2


       write(2,*) t, f, Phi((mux*t-muy)/sig), g, f-g, e, w!, log(f)
!            mux/sig + t*(muy - t*mux)*sigx**2/sig**3
!            b*d/pi2/sigx/sigy/a**3, &
!            Phi(b/a) - Phi(-b/a), &
!            exp(-c/2)/pi/sigx/sigy/a**2, &
!            sqrt(2/pi)*a*exp(-(b/a)**2/2), &
!            (b**2 - c*a**2)/a**2/2, log(b*d/pi2/sigx/sigy/a**3), &
!            sqrt(2/pi)*a*exp(-(b/a)**2/2), &
!            b/a**3/sigx/sigy - &
!            (muy*t*sigx**2+mux*sigy**2)/(t**2*sigx**2+sigy**2)**1.5, &
!            ((b/a)**2-c)+(muy - t*mux)**2/(t**2*sigx**2+sigy**2), &
!            b/a**3/sigx/sigy - &
!            (muy*t*sigx**2 + mux*sigy**2)/sqrt(t**2*sigx**2+sigy**2)**3, &
!            b*d/a, (muy*t*sigx**2 + mux*sigy**2), sqrt(t**2*sigx**2+sigy**2)**3

    end do
    close(2)

    open(2,file='/tmp/dhinkley')
    f1 = 0
    df1 = 0
    dt = 1e-3
    do i = 50,5000
       t = i*dt

       sig = sqrt(t**2*sigx**2 + sigy**2)
       a = sqrt(t**2/sigy**2 + 1/sigx**2)
       b = muy*t/sigy**2 + mux/sigx**2
       c = muy**2/sigy**2 + mux**2/sigx**2
       d = exp((b**2 - c*a**2)/a**2/2)
       f = b*d/pi2/sigx/sigy/a**3*(Phi(b/a) - Phi(-b/a)) + &
            exp(-c/2)/pi/sigx/sigy/a**2
       g = b*d/pi2/sigx/sigy/a**3
!       e = (muy*t*sigx**2 + mux*sigy**2)*exp(-((muy-t*mux)/sig)**2/2)/sig**3/pi2
       w = (muy*t*sigx**2+mux*sigy**2)/sig**3/pi2
       e = ((muy-t*mux)/sig)**2
       beta = muy*t*sigx**2 + mux*sigy**2
       r = muy - t*mux
       dsig = t*sigx**2/sig
       df = (muy*sigx**2 - 3*beta*dsig/sig - &
            beta*(-mux*(muy-t*mux)/sig**2-r**2*dsig/sig**3)) &
            *exp(-r**2/2/sig**2)/pi2/sig**3

!       f = exp(-(t-2)**2/2)
!       df = (t-2)*f

!       df = -mux*beta*(muy-t*mux)*exp(-e/2)/pi2/sig**5

       if( abs(df) > 0 .and. abs(f) > 0 ) then
!       write(*,*) t,df,f,sig,sigx,sigy
          write(2,*) t, -df, -(f-f1)/dt, -df/f, -(df - df1)/dt, &
               muy*sigx**2/beta+ dsig/sig, &
               -mux*(muy-t*mux)/sig**2-r**2*dsig/sig**3
       end if
       f1 = f
       df1 = df

    end do
    close(2)

!    sig = sqrt(muy**2 + 2**2*mux**2)
!    write(*,*) mux**2/sig**2,sig,sigx,sigy

  end subroutine tabelate

  subroutine tabelate_params

    real(rp) :: mu, sigx, sigy, sigu, sigv, sig, r, t
    real(rp), dimension(6) :: f, g
    integer :: i, n, m
!    real(rp), parameter :: sqrtpi = 2.5066282746310005_rp

    ! mu
    r = 2
    open(2,file='/tmp/par_mu')
    open(3,file='/tmp/par_mua')
    do i = -2500,5000, 10
       m = 0
       do n = 38,43,1
          mu = 10.0**(n/10.0)
          sigu = 0.1*mu*10
          sigv = 0.2*mu*10
          sigu = 2e3
          sigv = 2e3
          t = i/1e3
          !          f(n) = pdf(t,r,mu,sigu,sigv)
          m=m+1
          f(m) = pdf2(t,r*mu,mu,sigu,sigv)
          g(m) = pdf3(t,r*mu,mu,sqrt(sigv**2*t**2+sigu**2))
!!$          sigx = sqrt(mu + sigv**2)
!!$          sigy = sqrt(r*mu + sigu**2)
!!$          sig = sqrt(mu*(t**2 + r) + t**2*sigv**2 + sigu**2)
!!$!          sig = sqrt(t**2*sigx**2 + sigy**2)
!!$
!!$          f(n) = mu*((r - t)*t*sigx**2/sig**3 + 1/sig)* &
!!$               exp(-mu**2/2*(t-r)**2/sig**2)
!!$!          f(n) = mu*(t*r*sigx**2 + sigy**2)/sig**3*exp(-mu**2/2*(t-r)**2/sig**2)
!!$          f(n) = f(n) / sqrtpi
       end do
       write(2,*) t, f
       write(3,*) t, g
    end do
    close(2)
    close(3)

    ! r
    mu = 1e4
    sigu = 0.1*mu
    sigv = 0.2*mu
    open(2,file='/tmp/par_r')
    do i = 50,12000, 10
       do n = 1,4
          r = 2.0**(n-1)
          t = i/1e3
          f(n) = pdf(t,r,mu,sigu,sigv)
       end do
       write(2,*) t, f(1:4)
    end do
    close(2)

    ! sigu
    r = 2
    mu = 1e4
    sigv = 0.2*mu
    open(2,file='/tmp/par_sigu')
    do i = 1000,3500,10
       do n = 1,4
          sigu = 0.2*mu/2.0**(2*n-2)
          t = i/1e3
          f(n) = pdf(t,r,mu,sigu,sigv)
       end do
       write(2,*) t, f(1:4)
    end do
    close(2)

    ! sigv
    r = 2
    mu = 1e4
    sigu = 0.1*mu
    open(2,file='/tmp/par_sigv')
    do i = 1000,3500,10
       do n = 1,4
          sigv = 0.2*mu/2.0**(2*n-2)
          t = i/1e3
          f(n) = pdf(t,r,mu,sigu,sigv)
!          f(n) = pdf1(t,r,mu,sqrt(mu+sigv**2),sqrt(r*mu + sigu**2))
       end do
       write(2,*) t, f(1:4)
    end do
    close(2)


    ! 2D for t, mu
    r = 2
    sigu = 2e3
    sigv = 2e3
    open(2,file='/tmp/par_mu2d')
    do i = 1000,3500, 10
       t = i/1e3
       do n = 300,420,2
          mu = 10.0**(n/100.0)
          write(2,*) t, mu, pdf2(t,r*mu,mu,sigu,sigv)
       end do
       write(2,*)
    end do
    close(2)


    ! mu
    t = 2
    r = 2
    sigu = 2e3
    sigv = 2e3
    open(2,file='/tmp/par_mju')
    do i = 5000, 20000, 10
       mu = i
       write(2,*) mu, pdf2(t,r*mu,mu,sigu,sigv)
    end do
    close(2)


  end subroutine tabelate_params

  function pdf(t,r,mu,sigu,sigv)

    real(rp) :: pdf
    real(rp), intent(in) :: t,r,mu,sigu,sigv
    real(rp), parameter :: sqrtpi = 2.5066282746310005_rp
    real(rp) :: sig, sigx

    sigx = sqrt(mu + sigv**2)
    sig = sqrt(mu*(t**2 + r) + t**2*sigv**2 + sigu**2)
    pdf = mu*((r - t)*t*sigx**2/sig**3 + 1/sig)*exp(-mu**2/2*(t-r)**2/sig**2)
!    pdf = log(mu*((r - t)*t*sigx**2/sig**3 + 1/sig)) - mu**2/2*(t-r)**2/sig**2
!    pdf = exp(pdf)
    pdf = pdf / sqrtpi

  end function pdf

  function pdf1(t,r,mu,sigx,sigy)

    real(rp) :: pdf1
    real(rp), intent(in) :: t,r,mu,sigx,sigy
    real(rp), parameter :: sqrtpi = 2.5066282746310005_rp
    real(rp) :: sig

    sig = sqrt(t**2*sigx**2 + sigy**2)
    pdf1 = mu*(t*r*sigx**2 + sigy**2)/sig**3*exp(-mu**2/2*(t-r)**2/sig**2)
    pdf1 = pdf1 / sqrtpi


  end function pdf1

  function pdf2(t,muy,mux,sigu,sigv)

    real(rp) :: pdf2
    real(rp), intent(in) :: t,muy,mux,sigu,sigv
    real(rp), parameter :: sqrtpi = 2.5066282746310005_rp
    real(rp) :: sig, sigx, sigy

    sigx = sqrt(mux + sigv**2)
    sigy = sqrt(muy + sigu**2)
    sig = sqrt(t**2*sigx**2 + sigy**2)
    pdf2 = (muy*t*sigx**2+mux*sigy**2)/sig**3*exp(-((muy-t*mux)/sig)**2/2)/sqrtpi

  end function pdf2

  function pdf3(t,muy,mux,r)

    real(rp) :: pdf3
    real(rp), intent(in) :: t,muy,mux,r
    real(rp), parameter :: sqrtpi = 2.5066282746310005_rp
    real(rp) :: sig

    sig = sqrt(muy + t**2*mux + r**2)
    pdf3 = (mux*(muy*(t+1)+r**2+t*(2-t)*2e3**2))/sig**3*exp(-((muy - t*mux)/sig)**2/2)/sqrtpi

  end function pdf3

  subroutine phtest

    real(rp), dimension(:), allocatable :: x,u,y,v,res,cts,pht,sig,sigx2,sigy2,&
         a,b,c,d,sigr, sigx, sigy, tt

    real(rp) :: mean, fac, sigx0, sigy0, zero, mux, muy
    real(rp) :: t,f,e,ff,r,r0,f0, f1
    integer :: n, i

    mean = 1e4
    fac = 2
    zero = 0.0
    mux = mean
    muy = fac*mean
    sigx0 = 0.2*mean
    sigy0 = 0.4*mean

    n = 1300
    allocate(x(n),u(n),y(n),v(n),cts(n),pht(n))

    do i = 1, size(x)
       pht(i) = pnoise(mux)
       cts(i) = pnoise(muy)
       u(i) = gnoise(zero,sigy0)
       v(i) = gnoise(zero,sigx0)
       x(i) = pht(i) + v(i)
       y(i) = cts(i) + u(i)
    end do

!    write(*,*) sum(pht)/n, sqrt(sum(v**2)/n)

    call phrate(y,x,t,verbose=.true.)

  end subroutine phtest

  subroutine maxtab

    real(rp) :: mu, sigx, sigy, sigu, sigv, sig, r, t, f, fmax, tmax
    logical :: found
    integer :: i, n

    ! mu
    r = 2
    open(2,file='/tmp/maxtab')
    do n = 1,60
       mu = 10.0**(n/10.0)
       sigu = 0.1*mu
       sigv = 0.2*mu

       found = .false.
       fmax = -1
       tmax = -1
       do i = 100,2000
          t = i/1e3
          f = pdf(t,r,mu,sigu,sigv)
          if( f > fmax ) then
             fmax = f
             tmax = t
             found = .true.
          end if
       end do
       if( found ) write(2,*) mu, tmax
    end do
    close(2)

  end subroutine maxtab

  subroutine marginal

    real(rp) :: mu, mu0, sigu, sigv, sig, r, t, f, f1, f2, s
    logical :: found
    integer :: i, n, l

    open(1,file='/tmp/mq')
    open(2,file='/tmp/mf')
    r = 2
    sigu = 2e3
    sigv = 2e3
    do l = 0, 100

!       t = 0.7 + 0.3 * l / 1e2
       t = 1.9 + l / 500.0

       s = 0
       n = 0
       mu0 = 1e3
       do i = -100,200
          mu = 1e3*10.0**(i/100.0)
!      do i = 1, 1000
!          mu = 17e3 + i / 1
          f1 = pdf(t,r,mu,sigu,sigv)
          f2 = (1e3/mu)**1.5
!          f2 = 1
          s = s + f1*f2*(mu - mu0)
          n = n + 1
          mu0 = mu
          !if( l == 100 )
          write(2,*) mu,f1*f2,f1,f2,t
       end do
       write(2,*)
       f = s / n
       write(1,*) t, f

    end do
    close(1)
    close(2)

  end subroutine marginal

  subroutine themean

    use iromberg
    real(rp) :: f0, ft, s, t, w
    integer :: i

!    f0 =  romberg(phfun,0.1_rp,10.0_rp,1d-7)
!    ft =  romberg(meanfun,0.1_rp,10.0_rp,1d-7)
!    write(*,*) 'Mean:',ft/f0,f0,ft

    s = 0
    w = 0
    do i = 100,1000
       t = i / 100.0
       s = s + meanfun(t)
       w = w + phfun(t)
    end do
    write(*,*) 'Mean:',s/w

  end subroutine themean

  function meanfun(t)

    real(rp) :: meanfun
    real(rp), intent(in) :: t
    real(rp) :: mux, muy, sigx, sigy, sig, a,b,c,d
    real(rp), parameter :: pi = 3.14159265358979324_rp
    real(rp), parameter :: pi2 = 2.5066282746310005_rp
    real(rp), parameter :: pii2 = 1.25331413731550025_rp

    mux = 1e4
!    muy = 1.912e4
    muy = 2e4
    sigx = sqrt(0.2e4**2*0 + mux)
    sigy = sqrt(0.2e4**2*0 + muy)

    a = sqrt(t**2/sigy**2 + 1/sigx**2)
    b = muy*t/sigy**2 + mux/sigx**2
    c = muy**2/sigy**2 + mux**2/sigx**2
    d = exp((b**2 - c*a**2)/a**2/2)

!    phfun = b*d/pi2/sigx/sigy/a**3*(Phi(b/a) - Phi(-b/a)) + &
!         exp(-c/2)/pi/sigx/sigy/a**2
!    phfun = b*exp((b/a)**2/2)/pi2/sigx/sigy/a**3*(Phi(b/a) - Phi(-b/a)) + &
!         1/pi/sigx/sigy/a**2
!    phfun = exp(-c/2)*phfun

    sig = sqrt(t**2*sigx**2 + sigy**2)
    meanfun = t*(t*muy*sigx**2 + mux*sigy**2)/sig**3* &
         exp(-((muy - t*mux)/sig)**2/2)/pi2

  end function meanfun

  subroutine tabrfun

    real(rp) :: t,mu,r,phi,fac,sig,beta,dt,drho,dsig,dbeta,&
         tmin,tmax,f,f1,df,df1,thuge,ttiny,fmin,fmax,x,b,pmin,pmax
    real(real64) :: psi,dpsi,dpsi1,psi1,psimax,psimin
    real(real64), parameter :: qmin = 0.05, qmax = 0.95, &
         qtiny = 0.001, qhuge = 0.999
    integer :: i
    logical :: flag

    fac = 2
    mu = 1e4
    r = 4472.1359549
    phi = .4636476 !26.567 / 57.3

    call qbisection(qPhi,mu,fac,r,phi,qmin,0.1_REAL64,2.0_REAL64,tmin,flag)
    if( .not. flag ) stop 1
    call qbisection(qPhi,mu,fac,r,phi,qmax,2.0_REAL64,6.0_REAL64,tmax,flag)
    if( .not. flag ) stop 1
    call qbisection(qPhi,mu,fac,r,phi,qtiny,0.1_REAL64,2.0_REAL64,ttiny,flag)
    if( .not. flag ) stop 1
    call qbisection(qPhi,mu,fac,r,phi,qhuge,2.0_REAL64,9.0_REAL64,thuge,flag)
    if( .not. flag ) stop 1

    call rfunction(tmin,mu,fac,r,phi,fmin,psimin,dpsi)
    call rfunction(tmax,mu,fac,r,phi,fmax,psimax,dpsi)

    b = 5
    pmax = 1/tanh(b/2)
    pmin = 1/tanh(b/2)

    write(*,*) 'tmin = ',real(tmin),' tmax=',real(tmax)
    write(*,*) 'ttiny = ',real(ttiny),' thuge=',real(thuge)
    write(*,*) 'fmin = ',real(fmin),' fmax=',real(fmax)
    write(*,*) 'pmin = ',real(pmin),' pmax=',real(pmax)
    write(*,*) 'psimin = ',real(psimin),' psimax=',real(psimax)

    open(2,file='/tmp/rfun')
    dt = 1e-3
    psi1 = 0
    f1 = -1
    do i = 50,6000
       t = i*dt

       call rfunction(t,mu,fac,r,phi,f,psi,dpsi)
       if( t <= ttiny .or. t >= thuge) then
          f = 0
          psi = 0
          dpsi = 0
       else if( t < tmin ) then
          x = b*(ttiny - t)/(ttiny - tmin)/2
          f = f * (cosh(x)/cosh(b/2))**2
          psi = psimin*pmin*tanh(x)
          dpsi = psimin*pmin / cosh(x)**2 * (-b)/(ttiny - tmin)/2
       else if( t > tmax ) then
          x = b*(thuge - t)/(thuge - tmax)/2
          f = f * (cosh(x)/cosh(b/2))**2
          psi = psimax*pmax*tanh(x)
          dpsi = psimax*pmax / cosh(x)**2 * (-b)/(thuge - tmax)/2
       end if
       dpsi1 = (psi - psi1)/dt
       if( t < ttiny .or. t > thuge) then
          df = 0
          df1 = 0
       else
          df = (f - f1)/dt
          df1 = -df/f
       end if
       if( f1 >= 0 ) &
            write(2,*) t,f,psi,dpsi,df1,dpsi1

       psi1 = psi
       f1 = f
    end do
    close(2)

  end subroutine tabrfun

  subroutine rfunction(t,mu,tau,r,fi,f,psi,dpsi)

    real(REAL64), intent(in) :: t,mu,tau,r,fi
    real(REAL64), intent(out) :: f,psi,dpsi
    real(REAL64) :: rho,drho,beta,dbeta,sig,dsig,d2rho,d2sig
    real(REAL64) :: res,sig_t, sig_r, sig_f, beta_t, beta_r, beta_f
    real(REAL64) :: cf, sf, c2, cs

    sig = sqrt(mu*(tau + t**2) + r**2*(1 + ((t/tau)**2 - 1)*cos(fi)**2))
    beta = mu*tau*(t + 1) + r**2*(1 + (t/tau - 1)*cos(fi)**2)
    rho = (mu*(t - tau)/sig)**2/2
    f = mu*beta/sig**3*exp(-rho) / sqrt(2.0_REAL64*PI)

    dbeta = mu*tau + (r*cos(fi))**2/tau
    dsig = (mu*t + t*(r/tau*cos(fi))**2)/sig
    drho = mu**2*((t - tau) - (t - tau)**2*dsig/sig)/sig**2
    d2sig = (mu + (r/tau*cos(fi))**2)/sig - dsig**2/sig
    d2rho = mu**2*(1/sig**2 - 2*(t - tau)*dsig/sig**3 - &
         (2*(t - tau)*dsig+(t - tau)**2*d2sig)/sig**3 + &
         3*(t - tau)**2*dsig**2/sig**4)

    psi = drho - dbeta/beta + 3*dsig/sig
    dpsi = d2rho + (dbeta/beta)**2 + 3*d2sig/sig - 3*(dsig/sig)**2

  end subroutine rfunction

  subroutine rfunction_tau(t,mu,tau,r,fi,psi,dpsi)

    real(REAL64), intent(in) :: t,mu,tau,r,fi
    real(REAL64), intent(out) :: psi,dpsi
    real(REAL64) :: beta,sig, res, sig_t, beta_t, res_t, &
         sig_tt, beta_tt, res_tt
    real(REAL64) :: cf, c2, cosfi

!    sf = sin(fi)
    cf = cos(fi)
    c2 = cf**2
    cosfi = cos(fi)
!    cs = cf*sf
!    cs = sin(2*fi)/2

    sig = sqrt(mu*(tau + t**2) + r**2*(1 + ((t/tau)**2 - 1)*cos(fi)**2))
    beta = mu*tau*(t + 1) + r**2*(1 + (t/tau - 1)*cos(fi)**2)
    sig_t = (mu - 2*(r*t*cf/tau)**2/tau)/(2*sig)
!    sig_r = r*(1 + ((t/tau)**2 - 1)*c2)/sig
!    sig_f = -r**2*cs*((t/tau)**2 - 1)/sig
    beta_t = mu*(t + 1) - r**2*t/tau**2*c2
!    beta_r = 2*r*(1 + (t/tau - 1)*c2)
!    beta_f = -2*r**2*(t/tau - 1)*cs
    res = mu*(tau - t)/sig

    psi = res*(mu/sig) - (res**2 - 3)*(sig_t/sig) - beta_t/beta
!    psi(2) = (res**2 + 3)*(sig_r/sig) - beta_r/beta
!    psi(3) = (res**2 + 3)*(sig_f/sig) - beta_f/beta
    psi = -psi


    sig = sqrt(mu*(tau + t**2) + r**2*(1 + ((t/tau)**2 - 1)*c2))
    beta = mu*tau*(t + 1) + r**2*(1 + (t/tau - 1)*c2)
    res = (mu/sig)*(tau - t)

    sig_t = mu/(2*sig) - r*((t/tau*cosfi)**2/tau)*(r/sig)
    beta_t = mu*(t + 1) - r**2*t/tau/tau*c2
    psi = -res*(mu/sig) + (res**2 - 3)*(sig_t/sig) + beta_t/beta


    sig = sqrt(mu*(tau + t**2) + r**2*(1 + ((t/tau)**2 - 1)*c2))
    beta = mu*tau*(t + 1) + r**2*(1 + (t/tau - 1)*c2)
    res = (mu/sig)*(tau - t)

    sig_t = mu/(2*sig) - r*((t*cosfi/tau)**2/tau)*(r/sig)
    beta_t = mu*(t + 1) - r**2*t/tau**2*c2
    res_t = mu/sig - res*(sig_t/sig)

    beta_tt = 2*r**2*t/tau**3*c2
    sig_tt = -(mu/sig)*(sig_t/sig)/2  + &
         r*(r/sig)*(t**2/tau**3)*c2*(sig_t/sig + 3/tau)
    res_tt = -(mu/sig)*(sig_t/sig) &
             -res_t*(sig_t/sig) - res*(sig_tt/sig) + res*(sig_t/sig)**2

    dpsi = -res*res_tt - res_t**2 &
         + beta_tt/beta - (beta_t/beta)**2 - 3*sig_tt/sig + 3*(sig_t/sig)**2
    dpsi = -dpsi

  end subroutine rfunction_tau

  subroutine rfunction_r(t,mu,tau,r,fi,psi,dpsi)

    real(REAL64), intent(in) :: t,mu,tau,r,fi
    real(REAL64), intent(out) :: psi,dpsi
    real(REAL64) :: beta,sig, res, sig_r, beta_r
    real(REAL64) :: cf, c2

    cf = cos(fi)
    c2 = cf**2

    sig = sqrt(mu*(tau + t**2) + r**2*(1 + ((t/tau)**2 - 1)*cos(fi)**2))
    beta = mu*tau*(t + 1) + r**2*(1 + (t/tau - 1)*cos(fi)**2)
!    sig_t = (mu - 2*(r*t*cf/tau)**2/tau)/(2*sig)
    sig_r = r*(1 + ((t/tau)**2 - 1)*c2)/sig
!    sig_f = -r**2*cs*((t/tau)**2 - 1)/sig
!    beta_t = mu*(t + 1) - r**2*t/tau**2*c2
    beta_r = 2*r*(1 + (t/tau - 1)*c2)
!    beta_f = -2*r**2*(t/tau - 1)*cs
    res = mu*(tau - t)/sig

!    psi = res*(mu/sig) - (res**2 - 3)*(sig_t/sig) - beta_t/beta
    psi = (res**2 - 3)*(sig_r/sig) + beta_r/beta
!    psi1 = res**2*(sig_r/sig)
!    psi2 = beta_r/beta
!    psi3 = -3*(sig_r/sig)
!    psi(3) = (res**2 + 3)*(sig_f/sig) - beta_f/beta
    psi = -psi

    dpsi = 0

  end subroutine rfunction_r

  subroutine rfunction_fi(t,mu,tau,r,fi,psi,dpsi)

    real(REAL64), intent(in) :: t,mu,tau,r,fi
    real(REAL64), intent(out) :: psi,dpsi
    real(REAL64) :: beta,sig, res, sig_f, beta_f
    real(REAL64) :: cf, c2, cs

    cf = cos(fi)
    c2 = cf**2
    cs = sin(2*fi)/2

    sig = sqrt(mu*(tau + t**2) + r**2*(1 + ((t/tau)**2 - 1)*cos(fi)**2))
    beta = mu*tau*(t + 1) + r**2*(1 + (t/tau - 1)*cos(fi)**2)
!    sig_t = (mu - 2*(r*t*cf/tau)**2/tau)/(2*sig)
!    sig_r = r*(1 + ((t/tau)**2 - 1)*c2)/sig
    sig_f = -r**2*cs*((t/tau)**2 - 1)/sig
!    beta_t = mu*(t + 1) - r**2*t/tau**2*c2
!    beta_r = 2*r*(1 + (t/tau - 1)*c2)
    beta_f = -2*r**2*(t/tau - 1)*cs
    res = mu*(tau - t)/sig

!    psi = res*(mu/sig) - (res**2 - 3)*(sig_t/sig) - beta_t/beta
!    psi = (res**2 - 3)*(sig_r/sig) + beta_r/beta
!    psi1 = res**2*(sig_r/sig)
!    psi2 = beta_r/beta
!    psi3 = -3*(sig_r/sig)
    psi = (res**2 - 3)*(sig_f/sig) + beta_f/beta
    psi = -psi

    dpsi = 0

  end subroutine rfunction_fi

  function rpdf(fac,t,s,phi,mu)

    real(rp) :: rpdf
    real(rp), intent(in) :: fac,t,s,mu,phi
    real(rp), parameter :: sqrtpi = 2.5066282746310005_rp
    real(rp) :: sig, beta

    sig = sqrt(mu*(fac + t**2) + s**2*(1 + ((t/fac)**2 - 1)*cos(phi)**2))
    beta = mu*fac*(t + 1) + s**2*(1 + (t/fac - 1)*cos(phi)**2)
    rpdf = mu * beta * exp(-mu**2/2*(fac-t)**2/sig**2) / sig**3 / sqrtpi

  end function rpdf


  subroutine quantilephi

    real(REAL64), parameter :: q14 = 1.0_REAL64 / 4.0_REAL64
    real(REAL64), parameter :: q34 = 3.0_REAL64 / 4.0_REAL64
    real, parameter :: Qdev = 0.67448975

    real(REAL64) :: mu,t,r,f,a,b,c,t14,t34
    logical :: flag
    integer :: i

    mu = 1e4
    t = 2
    r = 4e3
    do i = 0, 90
       f = i / RAD

!       f = 30.0 / RAD

       call qbisection(qPhi,mu,t,r,f,q14,1.0_REAL64,2.0_REAL64,t14,flag)
       if( .not. flag ) stop 1
!       stop 666
       call qbisection(qPhi,mu,t,r,f,q34,2.0_REAL64,5.0_REAL64,t34,flag)
       if( .not. flag ) stop 2

       a = (mu*(t34 - t)/Qdev)**2 - mu*(t + t34**2)
       b = (mu*(t14 - t)/Qdev)**2 - mu*(t + t14**2)
       c = (t/r)**2 * (a - b) / (t34**2 - t14**2)
       write(*,*) f*RAD,c,acos(sqrt(max(0.0_REAL64,min(c,1.0_REAL64))))-f

    end do

  end subroutine quantilephi

  real(REAL64) function qPhi(x,mu,t,r,f,q)
    real(REAL64), intent(in) :: x,mu,t,r,f,q
    real(REAL64) :: z

    z = mu*(x - t)/sqrt(mu*(t + x**2) + r**2*(1 + ((x/t)**2 - 1)*cos(f)**2))
    qPhi = Phi(z) - q

  end function qPhi

  subroutine plotqPhi

    real(REAL64), parameter :: pi2 = sqrt(2*PI)
    real(REAL64) :: mu, tau, r, fi, q, t, dt, df, sig, dsig, z, f, c2, f1, beta
    integer :: i

    mu = 1e4
    tau = 2
    r = 4.4e3
    fi = 26 / RAD
    c2 = cos(fi)**2
    q = 0.1
    dt = 1e-2
    f1 = 0

    open(1,file='/tmp/qPhi')
    do i = 50, 600
       t = i *dt
       sig = sqrt(mu*(tau + t**2) + r**2*(1 + ((t/tau)**2 - 1)*c2))
       beta = mu*tau*(t + 1) + r**2*(1 + (t/tau - 1)*c2)
       dsig = t*(mu + (r/tau)**2*c2)/sig
       z = mu*(t - tau)/sig
       f = Phi(z) - q
       df = exp(-z**2/2)/(pi2*sig)*(1 - (t-tau)*(dsig/sig))*(mu/sig)
       df = mu*beta*exp(-z**2/2)/sig**3/pi2
       if( i > 50 ) &
            write(1,*) t,f,(f-f1)/dt,df
       f1 = f
    end do
    close(1)

!    dbeta = mu*tau + (r*cos(fi))**2/tau
!    dsig = (mu*t + t*(r/tau*cos(fi))**2)/sig
!    drho = mu**2*((t - tau) - (t - tau)**2*dsig/sig)/sig**2
!    psi = drho - dbeta/beta + 3*dsig/sig


  end subroutine plotqPhi


  subroutine qbisection(fun,mu,t,r,phi,q,xmin,xmax,x,flag)

    interface
       real(REAL64) function fun(x,mu,t,r,phi,q)
         use iso_fortran_env
         real(REAL64), intent(in) :: x,mu,t,r,phi,q
       end function fun
    end interface
    real(REAL64), intent(in) :: xmin, xmax, mu,t,r,phi,q
    real(REAL64), intent(out) :: x
    logical, intent(out) :: flag

    real(REAL64) :: a,b,fa,fb,f
    integer :: i

    flag = .false.
    a = xmin
    fa = fun(a,mu,t,r,phi,q)
    b = xmax
    fb = fun(b,mu,t,r,phi,q)
    do i = 1, 3*precision(x)
       x = (a + b) / 2.0_REAL64
       f = fun(x,mu,t,r,phi,q)
       if( fa*f < 0 ) then
          b = x
          fb = f
       else if( f*fb < 0 ) then
          a = x
          fa = f
!       else
!          write(*,*) '!',x,f
       end if
       if( abs(a - b)/2 < 1e3*epsilon(f) .or. abs(f) < 1e3*epsilon(f) ) then
          flag = .true.
          exit
       end if
!       write(*,*) i,a,b,a-b
    end do

  end subroutine qbisection


end program testphrate
