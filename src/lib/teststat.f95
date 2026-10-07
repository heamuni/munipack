
!
! gfortran -Wall  -fcheck=all  noise.f08 teststat.f95 -L. -L../minpack -lrstat -lsort -llmin -lfmm -lminpacks -lminpack -lm
!

program teststat

  use robustmean
  use weightedmean
  use qmeans
  use rfun
  use noise

  implicit none

  integer, parameter :: rp = selected_real_kind(15)
  real(rp), parameter :: rhomax = log(huge(rp)/2.1)

  real(rp) :: x(7), t(66666), dt(66666), z(7), w(16), psi(66666), rho(66666), dpsi(66666)
  real(rp) :: u,v,s,ss
  integer :: i,n
  real :: t1,t2
  logical :: reli

!  goto 20

!!$  w = (/10035.000000000000,        10001.000000000000,        10030.000000000000 ,    &
!!$       10058.000000000000 ,       10072.000000000000 ,       10012.000000000000  , &
!!$       10043.000000000000  ,      10061.000000000000  ,      10013.000000000000  ,   &
!!$       10047.000000000000    ,    10031.000000000000   ,     10022.000000000000   , &
!!$       10010.000000000000   ,     10072.000000000000  ,      10014.000000000000     ,  &
!!$       10027.000000000000   ,     10033.052157716564 /)
!!$
!!$  w = (/  10049.000000000000   ,     10063.000000000000 ,       10059.000000000000 , &
!!$       9994.0000000000000  ,      10063.000000000000 ,       10016.000000000000 ,&
!!$       10040.000000000000  ,      10014.000000000000 ,       10058.000000000000 ,&
!!$       10012.000000000000  ,      10019.000000000000 ,       10016.000000000000 ,&
!!$       10064.000000000000  ,      10032.000000000000 ,       10070.000000000000 ,&
!!$       10034.000000000000  ,      10038.500000000040 /)

  w = (/ 10109.000000000000 ,       10109.000000000000  ,      10109.000000000000 , &
       10109.000000000000 ,       10109.000000000000  ,      10109.000000000000  , &
       10109.000000000000 ,       10109.000000000000   ,     10109.000000000000  , &
       729415.00000000000 ,       10109.000000000000  ,      10109.000000000000 , &
       10109.000000000000 ,       10109.000000000000  ,      10109.000000000000 , &
       10109.000000000000 /)

  w(1:9) = (/ 0.00000000,       0.00000000,       0.00000000,       0.00000000, &
       0.00000000,       0.00000000,       109.288963,       18.7713432, &
       0.00000000 /)

  call rmean(w(1:9),u,v,s)
  write(*,*) "Rmean (patological):",u, v,s
  u = sum(w(1:9))/9
  v = sqrt(sum((w(1:9)-u)**2)/8)
  write(*,*) "Amean (patological):",u,v/sqrt(9.0),v
  call qmean(w(1:9),u,v)
  write(*,*) "Qmean (patological):",u, v, v/sqrt(9.0)

  ! Another patological
  w(1:4) = (/2.25176118E-05, 4.98100162E-05, 2.31736030E-05,   9.98998730E-05/)
  call qmean(w(1:4),u,v)
  write(*,*) "Qmean (another patological):",u, v, v/sqrt(4.0)

  t(1:10) = (/ 0.998451293, 0.999660850,       1.00406659,       1.00311661, &
       1.00328016, 0.997870922, 0.998634219, 1.00346828, 1.00324047, 1.00227821/)
  dt(1:10) = (/1.56288221E-03,  1.58239249E-03,   1.71578780E-03,   1.52544561E-03, &
       1.66450685E-03,   1.63013267E-03,   1.79429608E-03,   1.97814801E-03,  &
       1.90186617E-03,   2.12184014E-03 /)
!  dt = dt(1)
  call rwmean(t(1:10),dt(1:10),u,v)
  write(*,*) "RWmean (another patological):",u, v


  t(1:14) = (/1.01148117, 1.00847173, 1.00453901, 1.00900149, 1.00946915, 1.00475299, &
       1.01501632, 1.01952410, 1.02912927, 1.01704729, 1.01167893, 1.00378621, &
       1.00812364,       1.01006317/)
  dt(1:14) = (/9.00068413E-03,1.02015464E-02,1.14752613E-02,1.28761921E-02, &
       1.43292071E-02,4.06658743E-03, 4.41102358E-03, 4.76694433E-03,  &
       5.16087562E-03, 5.52221620E-03,5.91059681E-03, 6.32086722E-03, &
       6.86141895E-03, 7.88705982E-03 /)
  call rwmean(t(1:14),dt(1:14),u,v)
  write(*,*) "RWmean (another patological):",u, v

  t(1:5) = (/1.43887237E-05,   1.31587876E-04,   6.08701703E-05,  &
       8.92234530E-05,   7.58275419E-05 /)
  dt(1:5) = (/ 7.17469447E-06,   2.02281808E-05,   1.49832149E-05,   1.97298232E-05, &
       1.98863454E-05 /)
  call rwmean(t(1:5),dt(1:5),u,v)
  write(*,*) "RWmean (another patological):",u, v
  call rmean(t(1:5),u,v,s)
  write(*,*) "Rmean (another patological):",u, v,s

  do i = 0,20
     u = i/20.0
!     v = invnorm(u)
!     write(*,*) u,v,(1+erf(v/sqrt(2.0)))/2-u
  end do
!  stop

  x = (/16, 12, 99,95,18,87,10 /)
  z = (/142, 141, 149, 149, 142, 148,149 /)

  call rmean(z,u,v,s)
  write(*,*) "Rmean (dark):",u, v,s

  u = sum(z)/size(z)
  v = sqrt(sum((z-u)**2)/(size(z)-1))
  write(*,*) "Amean: (dark)",u,v/sqrt(1.0*size(z)),v
  call qmean(z,u,v)
  write(*,*) "Qmean (dark):",u, v, v/sqrt(1.0*size(z))

  !  stop

!  call medmad(x,u,v)
!  write(*,*) u,v

  20 continue

  open(1,file='s')
  do i = 1, size(t)
     t(i) = gdis(0.0_rp,1.0_rp)
     write(1,*) t(i)
  end do
  close(1)
  u = sum(t)/size(t)
  s = sqrt(sum((t-u)**2)/(size(t)-1.0))
  write(*,*) "Amean N(0,1): ",u, s, s/sqrt(real(size(t)-1))
!  stop
!  write(*,*) t
  call rmean(t,u,v,s)
  write(*,*) "Rmean N(0,1): ",u, s, v

!  stop

  s = 0
  open(1,file='huber')
  do i = -50000,50000
     u = i/10000.0
!     write(1,*) u,ihuber(u)
!     s = s + exp(-ihuber(u))
     write(1,*) u,s*0.001/sqrt(2*3.14159)
!     write(1,*) u,exp(-ihuber(u))/sqrt(2*3.14159),exp(-u**2/2)/sqrt(2*3.14159)
!     if( i > 0 .and. i < 100) s = s + exp(-u**2/2)/sqrt(2*3.14159)
!     if( i > 0 .and. i < 100) s = s + exp(-ihuber(u))/sqrt(2*3.14159)
!     if( i > 0 .and. s < 0.25) s = s + 0.01*exp(-ihuber(u))/sqrt(2*3.14159)
     if( i > -100000000 ) then
        if( s < 10.25) then
!           s = s + 0.0001*exp(-u**2/2)/sqrt(2*3.14159)*u**2
!           s = s + 0.0001*exp(-ihuber(u))/sqrt(2*3.14159)*u**2
           s = s + 0.0001*exp(-abs(u))/2*u**2
        else
           write(*,*) u
           exit
        end if
     end if
  end do
  close(1)
!  write(*,*) 'quantil:',2*s
  write(*,*) 'quantil:',s


!  stop

30 continue

  n = int(0.9*size(t))
  do i = 1, n
     t(i) = gdis(0.0_rp,1.0_rp)
  end do
!  n = 0
!  n = 999
  write(*,*)
  do i = n+1, size(t)
     t(i) = gdis(1.0_rp,10.0_rp)
  end do

  call cpu_time(t1)
  u = sum(t)/size(t)
  v = sqrt(sum((t-u)**2)/(size(t)-1.0))
  call cpu_time(t2)
  write(*,*) "Amean 0.9: N(0,1) + 0.1: N(0,100): ",u, v, v/sqrt(size(t)-0.0), &
       'elapsed=',t2-t1

  call cpu_time(t1)
  call rinit(t,u,v)
  call cpu_time(t2)
  write(*,*) "Rinit:",u, v,'elapsed=',t2-t1

!  s = v
!  call cpu_time(t1)
!  call rnewton(t,u,v,s,1d-15,reli)
!  call cpu_time(t2)
!  write(*,*) "Rnewton 0.9: N(0,1) + 0.1: N(0,100): ",u, v,s,reli,'elapsed=',t2-t1

  call rinit(t,u,v)
  s = v
  call cpu_time(t1)
  call rmean1(t,u,v,s,1d-15,i)
  call cpu_time(t2)
  write(*,*) "Rmean1 0.9: N(0,1) + 0.1: N(0,100): ",u, v,s,i,'elapsed=',t2-t1

  call rinit(t,u,v)
  s = v
  call cpu_time(t1)
  call rmean2(t,u,v,s,1d-15,i)
  call cpu_time(t2)
  write(*,*) "Rmean2 0.9: N(0,1) + 0.1: N(0,100): ",u, v,s,i,'elapsed=',t2-t1

  call rinit(t,u,v)
  s = v
  call cpu_time(t1)
  call rmean(t,u,v,s)
  call cpu_time(t2)
  write(*,*) "Rmean 0.9: N(0,1) + 0.1: N(0,100): ",u, v,s,'elapsed=',t2-t1

  call histogram(t,-5.0_rp,5.0_rp,50)

!  stop 0

  open(1,file='info')
!  do i = -200,200
  do i = 10,10
     u = i / 10.0
     !     do n = 20*30,10*666,200
     do n = 1,500
        v = n / 100.0
        call hubers((t - u)/v,psi)
!        s = sum(psi**2)/size(t) !- v**2*0.8
!        psi = psi**2
!        s = - sum(exp(-psi)*psi)
        call ihubers((t - u)/v,rho)
        call dhubers((t - u)/v,dpsi)
!        psi = ((t - u)/v)**2/2
!        psi = 2.1*psi
        ss = - sum(exp(-rho)*rho, rho < rhomax) /count(rho < rhomax)
!        s = - (sum(exp(-rho)*rho) + log(2.66*v)*sum(exp(-rho)))/(2.66*v)
!        s = - sum(exp(-rho)*(rho + log(1.66*v)))/(1.66*v) /size(rho)
!        s = -log(sum(exp(-rho)))
!        s = sum(exp(-rho)*psi**2)
!        s = sum(psi**2) - 0.7*size(psi)
!        s = sum(exp(-rho)*dpsi)**2
!        s = sum(exp(-((t-u)/v)**2)*((t-u)/v)**2)
!        s = - (sum(exp(-psi)*psi) + log(2.5*v)*sum(exp(-psi)))/(2.5*v)
        !        write(1,*) u,v,1/s,ss,sum(exp(-((t - u)/v)**2/2)*((t - u)/v)**2/2)
!        ss = sum(rho) / sqrt(real(size(rho)))
        write(1,*) u,v,ss,s
     end do
  end do
  close(1)

!  stop 0
50 continue

  n = int(0.1*size(t))
  do i = 1, n
     t(i) = pnoise(100000.0_rp)
  end do
!  n = 0
!  n = 999
  write(*,*)
  do i = n+1, size(t)
     t(i) = pnoise(10000.0_rp)
  end do
  call cpu_time(t1)
  u = sum(t)/size(t)
  v = sqrt(sum((t-u)**2)/(size(t)-1.0))
  call cpu_time(t2)
  write(*,*) "Amean 0.9: N(0,1) + 0.1: N(0,100): ",u, v, v/sqrt(size(t)-0.0), &
       'elapsed=',t2-t1

  dt = sqrt(t)
!  where( t >= 1)
!     dt = sqrt(t)
!  elsewhere
!     dt = 1
!  end where
!  where( dt == 0 )
!     dt = huge(dt)
  !  end where
  do i = 1,size(dt)
     if( dt(i) < epsilon(dt) ) write(*,*) t(i),dt(i)
  end do

  u = sum(t/dt**2)/sum(1/dt**2)
  v = sqrt(sum((t-u)**2/dt**2)/size(t))
  v = sqrt(v**2/sum(1/dt**2))
  write(*,*) "WAmean 0.9: N(0,1) + 0.1: N(0,100): ",u, v, sum((t-u)**2/dt**2)/size(t),&
       'elapsed=',t2-t1

  !  dt = 1
!  t = sqrt(t)
!  dt = 1
  call cpu_time(t1)
  call rwmean(t,dt,u,v)
  call cpu_time(t2)
  write(*,*) "rwmean:",u, v,'elapsed=',t2-t1,'s'

!  call rinit(t,u,v)
!  s = v
  call cpu_time(t1)
  call rmean(t,u,v,s)
  call cpu_time(t2)
  write(*,*) "Rmean 0.9: N(0,1) + 0.1: N(0,100): ",u, v,s,'elapsed=',t2-t1

  ! double Normal
  write(*,*)
  do i = 1, size(t)
     call random_number(u)
     if( u > 0.5 ) then
        t(i) = gnoise(1.0_rp,0.1_rp)
     else
        t(i) = gnoise(0.0_rp,0.1_rp)
     end if
  end do
  u = sum(t)/size(t)
  s = sqrt(sum((t-u)**2)/(size(t)-1.0))
  v = s / sqrt(size(t)-0.0)
  write(*,*) "Amean 0.5: N(0,1) + 0.5: N(0,1): ",u,v,s
  call rmean(t,u,v,s)
  write(*,*) "Rmean 0.5: N(0,1) + 0.5: N(0,1): ",u,v,s


  contains

    function gdis(mean, sig)

      ! generate random data with gauss distribution

      real(rp), intent(in) :: mean, sig
      real(rp) :: gdis, x

!      call random_number(x)
!      gdis = invnorm(x)*sig + mean!)!*sig !+ mean !- sig !/ mean
!      gdis = invnorm(x)*sig + mean !- sig/mean
!     gdis = mean +sig - invdist(x)*sig  !- sig/mean
!      gdis = mean - sqrt(2.0)*sig*ierfc(2*x)   ! Numerical Recipes (6.14.3)
!      gdis = mean - 1.618*sqrt(2.0)*sig*ierfc(2*x)   ! Numerical Recipes (6.14.3)
      !      gdis = mean + sqrt(2.0)*sig*ierf(2*x-1)   ! Numerical Recipes (6.14.3)
      gdis = gnoise(mean,sig)

    end function gdis

    subroutine histogram(x, xmin, xmax, nbines)

      real(rp), dimension(:), intent(in) :: x
      real(rp), intent(in) :: xmin, xmax
      integer, intent(in) :: nbines
      integer, dimension(:), allocatable :: hist
      real(rp) :: d,t,s
      integer :: n,l,nbins

      if( nbines < 1 ) then
         nbins = int(log(real(size(x)))/0.7) + 1
      else
         nbins = nbines
      end if
      nbins = nbins / 2

      allocate(hist(-nbins:nbins))

      d = (xmax - xmin) / (2*nbins + 1)
      t = (xmax + xmin) / 2
      hist = 0
      do l = 1, size(x)
         n = nint((x(l) - t) / d)
         if( -nbins <= n .and. n <= nbins ) then
            hist(n) = hist(n) + 1
         end if
      end do

      s = sum(hist) * d
      open(1,file='/tmp/h')
      do n = -nbins,nbins
         write(1,*) (n*d + t), hist(n)/s
      end do
      close(1)

      deallocate(hist)

    end subroutine histogram


!!$
!!$    function ierfc(x)
!!$
!!$      real(rp), intent(in) :: x
!!$      real(rp) :: ierfc
!!$
!!$      ierfc = ierf(x-1)
!!$
!!$    end function ierfc
!!$
!!$    function ierf(x)
!!$
!!$      ! http://en.wikipedia.org/wiki/Error_function
!!$
!!$      real(rp), intent(in) :: x
!!$      real(rp) :: ierf, u, v, a, y, w
!!$
!!$      a = 0.140012288686666
!!$      u = log(1 - x**2)
!!$      v = 2/3.14159265358979/a
!!$      y = v + u/2
!!$
!!$      w = sqrt(y**2 - u/a)
!!$      ierf = sign(1.0d0,x)*sqrt(w - y)
!!$
!!$    end function ierf
!!$
!!$    function invdist(xx)
!!$
!!$      real(rp), intent(in) :: xx
!!$      real(rp) :: invdist
!!$
!!$      ! inverzni fce k distribucni fci Gaussova rozdeleni
!!$      ! s presnosti vetsi jak 0.00045
!!$
!!$      real(rp) :: w,f,x
!!$      logical :: interval
!!$
!!$      x = xx
!!$      if( x < 0.0 ) then
!!$         invdist = 0.0
!!$      elseif( x > 1.0 )then
!!$         invdist = 1.0
!!$      else
!!$         interval = x < 0.5
!!$         if( .not. interval ) x = 1.0 - x + epsilon(1.0)
!!$         w = sqrt(-2.0*log(x));
!!$         f = -w + (2.515517 + w*(0.802853 + w*0.010328))/ &
!!$              (1.0 + w*(1.432788 + w*(0.189269 + w*0.001308)));
!!$         if( interval ) then
!!$            invdist = f
!!$         else
!!$            invdist = -f
!!$         endif
!!$      endif
!!$
!!$    end function invdist
!!$
!!$    ! excelent overwiew of algorithms !
!!$    ! http://home.online.no/~pjacklam/notes/invnorm/
!!$
!!$    function invnorm(p)
!!$
!!$      !ren-raw chen, rutgers business school
!!$      ! normal inverse
!!$      ! translate from
!!$      ! http://home.online.no/~pjacklam/notes/invnorm
!!$      ! a routine written by john herrero
!!$      real*8 invnorm
!!$      real*8 p,p_low,p_high
!!$      real*8 a1,a2,a3,a4,a5,a6
!!$      real*8 b1,b2,b3,b4,b5
!!$      real*8 c1,c2,c3,c4,c5,c6
!!$      real*8 d1,d2,d3,d4
!!$      real*8 z,q,r
!!$      a1=-39.6968302866538
!!$      a2=220.946098424521
!!$      a3=-275.928510446969
!!$      a4=138.357751867269
!!$      a5=-30.6647980661472
!!$      a6=2.50662827745924
!!$      b1=-54.4760987982241
!!$      b2=161.585836858041
!!$      b3=-155.698979859887
!!$      b4=66.8013118877197
!!$      b5=-13.2806815528857
!!$      c1=-0.00778489400243029
!!$      c2=-0.322396458041136
!!$      c3=-2.40075827716184
!!$      c4=-2.54973253934373
!!$      c5=4.37466414146497
!!$      c6=2.93816398269878
!!$      d1=0.00778469570904146
!!$      d2=0.32246712907004
!!$      d3=2.445134137143
!!$      d4=3.75440866190742
!!$      p_low=0.02425
!!$      p_high=1-p_low
!!$      if(p.lt.p_low) goto 201
!!$      if(p.ge.p_low) goto 301
!!$201   q=dsqrt(-2*dlog(p))
!!$      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
!!$      goto 204
!!$301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
!!$      if(p.gt.p_high) goto 302
!!$202   q=p-0.5
!!$      r=q*q
!!$      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
!!$      goto 204
!!$302   if((p.gt.p_high).and.(p.lt.1)) goto 203
!!$203   q=dsqrt(-2*dlog(1-p))
!!$      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
!!$204   invnorm=z
!!$      return
!!$
!!$    end function invnorm
!!$
!!$
!!$    function distgaus(t)
!!$
!!$      real(rp), parameter :: sqrtpi2 = 2.50662827463100050242
!!$      real(rp), intent(in) ::  t
!!$      real(rp) ::  w,f, x, distgaus
!!$      logical :: minus
!!$
!!$      x = t
!!$      minus = x < 0.0
!!$      if( x > 6.0 )then
!!$         distgaus = 1.0
!!$      elseif( x <= -6.0 )then
!!$         distgaus = 0.0
!!$      else
!!$         x = abs(x)
!!$         w = 1.0/(1.0 + 2.316419 - x)
!!$         f = 1.0 - exp(-x**2/2.0)/sqrtpi2*w* &
!!$              (0.3193815+w*(-0.3565638+w*(1.781478+w*(-1.821256+w*1.330274))))
!!$         if( minus )then
!!$            distgaus = 1.0 - f
!!$         else
!!$            distgaus = f
!!$         endif
!!$      endif
!!$
!!$    end function distgaus

end program teststat
