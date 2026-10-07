!
! An implementation of Trust region SQP (sequence quadratic programming)
! minimisation method for robust mean estimation.
!
!  Copyright © 2022-2024 F.Hroch (hroch@physics.muni.cz)
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
!
!  This implementation is due to the algorithms of the second edition
!  of Numerical Optimization by Jorge Nocedal and Stephen J. Wright,
!  modified and specialised Algorithm 18.1 (Local SQP) combined with
!  simplified TR step by More & Sorenson (1981).
!

module trustsr1mean

!  use rfun
  use modQ2eig
  use sr1mean_function
  use zeroinmod
  use iso_fortran_env

  implicit none
  private

  public :: trmean, tsr1mean, tsr1meansqp, trfmeansqp, smnewton

  type, extends(zeroin_function) :: zeroin_data
     class(tsr1mean_function), pointer :: funmean
     real(REAL64) :: t
   contains
     procedure :: f => zeroin_fun
  end type zeroin_data

  type, extends(zeroin_function) :: zeroin_levy
     real(REAL64) :: t, s, xcdf, ycdf
   contains
     procedure :: f => levyfun
  end type zeroin_levy


  logical, parameter :: debug = .true.
  real(REAL64), parameter :: eta = 0.5e-3  ! 0 < eta < 1e-3
  class(tsr1mean_function), pointer :: smin_fcn, iscale_fcn
  real(REAL64) :: smin_t, iscale_t
  integer :: smin_count, iscale_count
  integer, parameter :: recursive_depth_limit = 8
  integer :: recursive_depth

!!$  real(REAL64), dimension(:), pointer :: dat
!!$  real(REAL64), dimension(:), allocatable :: res,psi,dpsi
!!$  integer :: ndat
!!$  procedure(ihuber), pointer :: ifcn
!!$  procedure(huber),  pointer ::  fcn
!!$  procedure(dhuber), pointer :: dfcn

contains

  subroutine trmean(data,mean,scale,stderr,stddev,robfun,info,msg,verbose)

    use quantilmean

    real(REAL64), dimension(:), target, intent(in) :: data
    real(REAL64), intent(out) :: mean, scale
    real(REAL64), intent(out) :: stderr, stddev
    !    real(REAL64), dimesnion(:,:), intent(out) :: hess
    character(len=*), intent(in) :: robfun
    integer, intent(out) :: info
    character(len=*), intent(out) :: msg
    logical, intent(in) :: verbose

    type(tsr1mean_tukey), target :: tukeymean
    type(tsr1mean_huber), target :: hubermean
    type(tsr1mean_experimental), target :: bathmean
    class(tsr1mean_function), pointer :: funmean

    real(REAL64), dimension(2,2) :: hess
    real(REAL64) :: t, s, sig
    integer :: n, iflag
    logical :: reli

    write(*,*) 'START of trmean()'


    n = size(data)

    if( n < 2 ) then
       msg = 'trmean(): at least two data points are required.'
       write(error_unit,*) 'ERROR: ',trim(msg)
       info = 6
       return
    end if

    if( robfun == 'huber' ) then
       funmean => hubermean
    else if( robfun == 'hampel' ) then
       funmean => bathmean
    else
       funmean => tukeymean
    end if
    funmean%data => data

!    t = qmed(data,iflag)

!    if( n < 1 ) then
!       s = downton(data)
!    else
!    s = byQ(data,iflag)
!    end if
!    write(*,*) s
!    s = qmed(abs(data-t),iflag) / 0.6745
    sig = s
!    s = 0.8*s
    write(*,*) 'trmean(): initial rank estimate: t=',t,' sig=',sig
    if( iflag /= 0 ) then
       info = 6
       return
    end if

    ! update s rather than use sig
!    call rscale(funmean,t,sig,s,iflag)
!    write(*,*) 'Rscale: ',real(s),iflag
!    if( iflag /= 0 ) then
!       s = 0.8*sig
!    end if

!    call iscale(data,t,s)

    block
      type(tsr1mean_entropy), target :: etukey
      integer, parameter :: maxfev = 3*precision(t)
      logical, parameter :: quasi = .true.

      real(REAL64), dimension(2,2) :: B
      real(REAL64) :: tol, delta
      etukey%data => data

      tol = 0.1*s / sqrt(real(n))
      delta = 0.5*s

!      call trfmeansqp(etukey,t,s,B,delta,tol,maxfev,quasi,verbose,info)

      write(*,*) 'Entropy:',real(t),real(s),info
    end block


    recursive_depth = 0
    call tropt(funmean,t,s,hess,verbose,.true.,iflag)

    if( iflag == 0 ) then
       msg = ''
       mean = t
       scale = s
       stderr = funmean%stderr(t,s,reli)
       stddev = stderr*sqrt(real(n))
       write(*,*) 'FINAL t:',real(t),' s =',real(s)
       write(*,*) ' stderr:',real(stderr),' stddev =',real(stddev),reli
    else
       msg = 'A jeee...'
       write(*,*) 'Fail for data:',funmean%data
    end if
    info = iflag

  end subroutine trmean


  recursive subroutine tropt(fcn,t,s,hess,verbose,alt,iflag)

    use rfun
    use rfun_hampel

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), intent(in out) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: hess
    logical, intent(in) :: verbose, alt
    integer, intent(out) :: iflag

    integer, parameter :: maxfev = 3*precision(t)
    logical, parameter :: quasi = .true.

    real(REAL64), dimension(2,2) :: B
    real(REAL64), dimension(2) :: topt, sopt, fopt
    real(REAL64) :: tol, delta, K, curvature, frac, nu, H0, Hmin, fmin, fval, &
         KL, entropy, t0, s0, dist, d, stderr, e, p, prob, r, levy, rtol
    integer :: iopt, nopt, i, j,n, info, iflg
    logical :: reli

    t0 = t
    s0 = s
    n = size(fcn%data)

    if( alt ) then
       call ridge(fcn,t,s,nopt,topt,sopt)
    else
       nopt = 1
       topt(1) = t
       sopt(1) = s
    end if

    fmin = huge(fmin)
    Hmin = huge(H0)
    prob = 0
    curvature = 0
    frac = 0
    prob = 0
    dist = huge(dist)
    entropy = huge(entropy)
    stderr = 0.1*huge(stderr)
    iopt = 0
    do i = 1, nopt

       t = topt(i)
       s = sopt(i)

       tol = 0.01 * s / sqrt(real(n))
       rtol = 0.01
       delta = 0.5*s

       call trfmeansqp(fcn,t,s,B,delta,rtol,maxfev,verbose,info)
!       call tsr1meansqp(fcn,t,s,B,delta,tol,maxfev,quasi,verbose,info)
       write(*,*) 'trfmeansqp() info=',info
       topt(i) = t
       sopt(i) = s
       fopt(i) = fcn%val(t,s)

       if( info == 0 ) then
          K = B(1,1)*B(2,2) - B(1,2)**2
          K = s**2*K
          nu = fcn%nfrac(t,s)
          !H0 = lepage(fcn%data,t,s)
          H0 = ksdist(fcn%data,t,s)
          levy = levydist(fcn%data,t,s)
          fval = fcn%val(t,s)
!          fval = (fcn%val(t,s) - log(s))*n
!          fval = chi2(fcn%data,t,s)
          KL = relentropy(fcn%data,t,s)
          d = hypot(t-t0,2*(s-s0))
          e = fcn%stderr(t,s,reli)
!          p = sum(exp(-((fcn%data-t)/s)**2/2.0)/s/2.5)/n
          p = sum(exp(-ipshampel((fcn%data-t)/s))/s/2.5)/n
!          p = exp(-n*fval)
!          p = sum(exp(-abs(fcn%data-t)/s/sqrt(2.0))) / (n*s)
!          if( K >= curvature .and. nu >= frac ) then
!          if( K >= curvature .and. nu > 0 ) then
!          if( fval > 0 .and. fval < fmin .and. nu > 0 ) then
!          if( fval < fmin .and. nu > 0 ) then
          !          if( H0 < Hmin .and. nu > 0 ) then
          !          if( KL < entropy ) then
!          if( d < dist .and. nu > 0 ) then
          !          if( fval < fmin .and. nu > 0  ) then
          if( p > prob ) then
             dist = d
             entropy = KL
             curvature = K
             Hmin = H0
             fmin = 3*min(fval + sqrt(2.0*n),real(n,REAL64))
             frac = nu
             prob = p
             hess = B
             stderr = e
             iopt = i
          end if
          write(*,*) 'PARTICULAR t:',real(t),' s =',real(s), &
               ' K=',real(K),' nu=',real(nu),' KS test=',real(H0), &
               ' LePage=',real(lepage(fcn%data,t,s)), ' fopt=',fval, &
               ' entropy=',real(KL),' R-scale=',real(Rscale_test(fcn%data,t,s)),&
               ' dist=',real(d),' frac=',real(nu),' stderr=',real(e), &
               ' prob=',real(p),'levy=',real(levy)
!          write(*,'(5f10.3)') exp(-ipshampel((fcn%data-t)/s))/2.5/s
          write(*,'(5f10.3)') exp(-itukey((fcn%data-t)/s))/2.5/s
!          write(*,'(5f8.3)') exp(-((fcn%data-t)/s)**2/2)/s/2.5
!          write(*,'(5f8.3)') (fcn%data-t)/s
       else
          iflg = info
          write(error_unit,*) 'WARNING trfmeansqp() failed with iflag=',info
       end if

    end do

!    r = fopt(1) - fopt(2)
!    write(*,*) 'decision:',real(r)

!!$    block
!!$      real(REAL64), dimension(:), allocatable :: w
!!$      real(REAL64) :: qlike
!!$
!!$      allocate(w(n))
!!$
!!$      do j = 0,10
!!$         e = 0.1*j
!!$         do i = 1,n
!!$            p = (1-e)*pdf(fcn%data(i),topt(1),sopt(1)) + e*pdf(fcn%data(i),topt(2),sopt(2))
!!$            w(i) = e*pdf(fcn%data(i),topt(2),sopt(2)) / p
!!$         end do
!!$         qlike = sum((1-w)*fopt(1) + w*fopt(2))
!!$         write(*,'("mix: ",10f8.3)') e,w,qlike,sum(w)/n
!!$      end do
!!$    end block

    if( iopt > 0 ) then
       t = topt(iopt)
       s = sopt(iopt)
       iflag = 0
    else
       iflag = 666
    end if

  end subroutine tropt

  real(REAL64) pure function pdf(x,t,s)

    use rfun_hampel
    real(REAL64), intent(in) :: x,t,s

    pdf = exp(-ipshampel((x-t)/s))/s

  end function pdf


  real(REAL64) function Rscale_test(data,t,s)

    use rfun
    use noise
    use quantilmean
    use quicksorter

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t,s

    real(REAL64), dimension(:), allocatable :: xcdf, Fx, x, Fy, y,z
    integer :: n,i
    real(REAL64) :: R1, R2

 !   if( s > 2 ) return

    n = size(data)
    allocate(xcdf(n),Fy(n),x(n),Fx(n),y(n),z(n))
    call ecdf((data-t)/s,xcdf,Fx)

    z = (data-t)/s
    z = data
    call qsort(z)
    do i = 1, n
      x(i) = 6*real(i - n/2) / real(n+1)
       y(i) = sqrt(2.0)*inverf(real(2*real(i)/real(n+1) - 1,REAL64))
       Fy(i) = Phi(x(i),0.0_REAL64,1.0_REAL64)
    end do

    R1 = sum(Rscale_J(Phi(z,0.0_REAL64,1.0_REAL64)))
    R2 = sum(Rscale_J(Phi(y,0.0_REAL64,1.0_REAL64)))

    Rscale_Test = 1 + (R1 - R2) / sqrt(0.6*n)
!    Rscale_Test = abs(R1)

    if( s < 2 ) then
       open(8,file='/tmp/rs0')
    else
       open(8,file='/tmp/rs')
    end if
    do i = 1, n
       write(8,*) x(i),Fy(i),xcdf(i),Fx(i),&
            Rscale_J(Phi(z(i),0.0_REAL64,1.0_REAL64)),&
            Rscale_J(Phi(y(i),0.0_REAL64,1.0_REAL64)),&
            Phi(z(i),0.0_REAL64,1.0_REAL64)
    end do
    close(8)

  end function Rscale_Test

  elemental real(REAL64) function Rscale_J(x)

    use rfun_hampel
    real(REAL64), intent(in) :: x

    Rscale_J = pshampel(x)*x - 1
!    Rscale_J = huber(x)*x - 1
!    Rscale_J = abs(x - 0.5) - 0.25
!    Rscale_J = x**2 - 1

  end function Rscale_J

  real(REAL64) function relentropy(data,t,s)

    use rfun_hampel
    use quantilmean

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t,s

    real(REAL64), dimension(:), allocatable :: xcdf,ycdf
    integer :: n,i
    real(REAL64) :: D, Fi, r

!    if( s > 2 ) return

    n = size(data)
    allocate(xcdf(n),ycdf(n))
    call ecdf((data-t)/s,xcdf,ycdf)

    D = 0
    if( s < 2 ) then
       open(8,file='/tmp/kl0')
    else
       open(8,file='/tmp/kl')
    end if
    do i = 1,n
!       Fi = Phi(xcdf(i),t,s)
       Fi = Phi(xcdf(i),0.0_REAL64,1.0_REAL64)
       r = Fi/ycdf(i)
       if( r > 0 ) D = D + Fi*log(r)
       if( r > 0 ) write(8,*) xcdf(i),ycdf(i),Fi,Fi*log(r),r
    end do
    close(8)

!    call ecdf(sqrt(ipshampel((data-t)/s)),xcdf,ycdf)
    call ecdf(abs(pshampel((data-t)/s)),xcdf,ycdf)
!    call ecdf(abs((data-t)/s),xcdf,ycdf)
    if( s < 2 ) then
       open(8,file='/tmp/skl0')
    else
       open(8,file='/tmp/skl')
    end if
    do i = 1,n
!       Fi = Phi(xcdf(i),t,s)
       Fi = erf(xcdf(i)/sqrt(2.0))
       r = Fi/ycdf(i)
!       if( r > 0 ) D = D + Fi*log(r)
       write(8,*) xcdf(i),ycdf(i),Fi!,Fi*log(1/r)
    end do
    close(8)


    relentropy = D / n
!    relentropy = sum(Phi(xcdf(2:),t,s)*log(Phi(xcdf(2:),t,s)/ycdf(2:))) / n

  end function relentropy


  real(REAL64) function ksdist(data,t,s)

    use quantilmean

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t,s

    real(REAL64), dimension(:), allocatable :: xcdf,ycdf
    integer :: n

    n = size(data)
    allocate(xcdf(n),ycdf(n))
    call ecdf(data,xcdf,ycdf)

    ksdist = maxval(abs(ycdf - Phi(xcdf,t,s)))

  end function ksdist

  real(REAL64) function levydist(data,t,s)

    use quantilmean

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t,s

    type(zeroin_levy) :: levyzero
    real(REAL64), dimension(:), allocatable :: xcdf,ycdf
    integer :: n, i, iflag
    real(REAL64) :: x,d,tol,levy

    n = size(data)
    allocate(xcdf(n),ycdf(n))
    call ecdf(data,xcdf,ycdf)

    open(8,file='/tmp/levy')
    do i = 1, n
       write(8,*) xcdf(i),ycdf(i)
    end do
    close(8)

    levyzero%t = t
    levyzero%s = s
    tol = 1e-3
    levy = 0
    do i = 2, n-1
       levyzero%xcdf = xcdf(i)
       levyzero%ycdf = ycdf(i)
       x = zeroin(levyzero,xcdf(1),xcdf(n),tol,iflag)
       if( iflag == 0 ) then
          d = sqrt((xcdf(i) - x)**2 + (ycdf(i) - Phi(x,t,s))**2) / sqrt(2.0)
          levy = max(levy,d)
       end if
 !      write(*,*) i,real(x),real(d),iflag
    end do
!    write(*,*) 'Levy dist:',levy
    levydist = levy

  end function levydist

  real(REAL64) function levyfun(this,x)

    class(zeroin_levy) :: this
    real(REAL64), intent(in) :: x
    real(REAL64) :: y

    y = this%ycdf + (x - this%xcdf)
    levyfun = Phi(x,this%t,this%s) - y
!    write(*,*) x,Phi(x,this%t,this%s)-y

  end function levyfun

  elemental real(REAL64) function Phi(x,t,s)
    real(REAL64), intent(in) :: x,t,s
    real(REAL64), parameter :: sqrt2 = sqrt(2.0_REAL64)

    Phi = (1 + erf(((x-t)/s)/sqrt2))/2

  end function Phi

  real(REAL64) function lepage(data,t,s)

    use quicksorter

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), allocatable :: x,y,q
    integer, dimension(:), allocatable :: V
    integer :: i,j,n,m,l,mn,n2
    real(REAL64) :: u,mu1,mu2,sig1,sig2,T1,T2

    m = size(data)
    n = 2*m
    mn = n + m
    allocate(q(mn),x(m),y(n),V(mn))
    do j = 1, n
       u = real(j,REAL64) / real(n+1,REAL64)
       y(j) = inorm(u)
    end do

    x = (data - t) / s
    q(1:m) = x
    q(m+1:) = y

    V(1:m) = 1
    V(m+1:) = 0
!    write(*,*) q
    call qsort_lepage(q,V)

!    write(*,*) V
!    write(*,*) q
!    stop

    n2 = (mn + 1) / 2.0
    T1 = 0
    T2 = m*n2
    do i = 1, mn
       T1 = T1 + i*V(i)
       T2 = T2 - abs(i - n2)*V(i)
    end do

    mu1 = m*n2
    sig1 = m*n*(mn+1)/12.0
    if( sig1 > 0 ) then
       sig1 = sqrt(sig1)
    else
       sig1 = huge(sig1)
    end if
    if( mod(mn,2) == 0 ) then
       ! even
       mu2 = m*(mn + 2)/4.0
       sig2 = m*n*(mn**2 - 4)/(48.0*(mn - 1))
    else
       mu2 = m*(mn + 1)**2/mn/4.0
       sig2 = m*n*(mn + 1)*(mn**2 + 3) / (48.0*mn**2)
    end if
    if( sig2 > 0 ) then
       sig2 = sqrt(sig2)
    else
       sig2 = huge(sig2)
    end if

    lepage = ((T1 - mu1)/sig1)**2 + ((T2 - mu2)/sig2)**2

  end function lepage

  real(REAL64) function inorm(x)

    use noise
    real(REAL64), intent(in) :: x

    !inorm = sqrt(real(2,REAL64))*inverf(real(2,REAL64)*x - real(1,REAL64))
    inorm = sqrt(real(2,REAL64))*inverf(2*x - 1)

  end function inorm


  subroutine trmean_orig(data,mean,scale,stderr,stddev,robfun,info,verbose)

    real(REAL64), dimension(:), target, intent(in) :: data
    real(REAL64), intent(in out) :: mean, scale
    real(REAL64), intent(out) :: stderr, stddev
    !    real(REAL64), dimesnion(:,:), intent(out) :: hess
    character(len=*), intent(in) :: robfun
    integer, intent(out) :: info
    logical, intent(in) :: verbose

    type(tsr1mean_tukey), target :: tukeymean
    type(tsr1mean_huber), target :: hubermean
    type(tsr1mean_experimental), target :: bathmean
    class(tsr1mean_function), pointer :: funmean
    real(REAL64) :: tol, delta, s, t, K, curvature
    real(REAL64), dimension(2,2) :: hess
    real(REAL64), dimension(2) :: topt, sopt
    integer, parameter :: maxfev = 3*precision(mean)
    logical, parameter :: quasi = .true.
    logical :: reli
    character(len=80) :: errmsg
    character(len=666) :: msg
    integer :: iopt, nopt,i,n

    write(*,*) 'start of trmean(): mean=',mean,' scale=',scale

    n = size(data)
    if( n < 2 ) then
       write(error_unit,*) 'ERROR: trmean(): at least two data points are required.'
       info = 6
       return
    end if

    if( .not. (scale > 0) ) then
       write(error_unit,*) 'ERROR: trmean(): Scale parameter should be s > 0.'
       info = 6
       return
    end if

    if( robfun == 'huber' ) then
       funmean => hubermean
    else if( robfun == 'hampel' ) then
       funmean => bathmean
    else
       funmean => tukeymean
    end if
    funmean%data => data

    t = mean
    s = scale
    delta = s
!    call ridge(funmean,t,s,delta,nopt,topt,sopt)

    curvature = 0
    iopt = 0
    do i = 1, nopt

       t = topt(i)
       s = sopt(i)

       tol = 0.01 * s / sqrt(real(n))
       delta = 0.5*s

!       call trfmeansqp(funmean,t,s,hess,K,delta,tol,maxfev,quasi,verbose,info,msg)
       write(*,*) 'trfmeansqp() info=',info,' msg:',trim(msg)

       if( info == 0 ) then
          topt(i) = t
          sopt(i) = s
!          K = hess(1,1)*hess(2,2) - hess(1,2)**2
          if( K > curvature ) then
             curvature = K
             iopt = i
          end if
       else
          write(error_unit,*) 'WARNING trfmeansqp() failed: ',msg
       end if

    end do

    if( iopt > 0 ) then
       t = topt(iopt)
       s = sopt(iopt)
       stderr = funmean%stderr(t,s,reli)
       stddev = stderr*sqrt(real(n))
       write(*,*) 'FINAL t:',real(t),' s =',real(s),' K=',real(curvature)
       write(*,*) ' stderr:',real(stderr),' stddev =',real(stddev)
       mean = t
       scale = s
       if( .not. reli ) stop 'stderr unreliable.'
    end if


  end subroutine trmean_orig



  subroutine tsr1mean(data,mean,stderr,scale,info,verbose)

    real(REAL64), dimension(:), target, intent(in) :: data
    real(REAL64), intent(in out) :: mean, scale
    real(REAL64), intent(out) :: stderr
!    real(REAL64), intent(in) :: smin
    integer, intent(out) :: info
    logical, intent(in) :: verbose

!    type(sr1_mean) :: frob

    type(tsr1mean_tukey), target :: tukeymean
    type(tsr1mean_huber), target :: hubermean
!    type(tsr1mean_square), target :: squaremean
    type(tsr1mean_experimental), target :: bathmean
    class(tsr1mean_function), pointer :: fcn
    real(REAL64) :: sqrtn, tol, delta, s, t, K
    real(REAL64), dimension(2,2) :: hess
    integer, parameter :: maxfev = 3*precision(mean)
    logical, parameter :: quasi = .true.
    logical :: reli
    character(len=80) :: errmsg
    character(len=666) :: msg
    integer :: stat

    if( .not. (size(data) > 1 .and. scale > 0) ) then
!    if( .not. (scale > 0) ) then
       ! improper parameters
       stop 'Parameters improperly defined.'
    end if

    fcn => tukeymean
!    fcn => hubermean
!    fcn => bathmean
    fcn%data => data


!!$    dat = data
!!$    ndat = size(data)
!!$    allocate(res(ndat),psi(ndat),stat=stat,errmsg=errmsg)
!!$    if( stat /= 0 ) then
!!$       write(error_unit,*) 'sr1rmean(): ',trim(errmsg)
!!$       stop 'in tsr1mean'
!!$    end if
!!$
!!$    if( .not. quasi ) then
!!$       allocate(dpsi(ndat),stat=stat,errmsg=errmsg)
!!$       if( stat /= 0 ) then
!!$          write(error_unit,*) 'sr1rmean(): ',trim(errmsg)
!!$          stop 'in tsr1mean'
!!$       end if
!!$    end if

    !    x = [ mean, scale ]
    t = mean
    s = scale
    sqrtn = sqrt(real(size(data)))
    tol = 0.01 * scale / sqrtn
!    tol = 1e-5
!    diag = [ 1, 2 ] / scale**2
!    write(*,*) 'diag:',real(diag)
!    smin = 1e-3 * scale / sqrtn
    !    x(2) = max(smin,x(2))
!
    delta = 0.1*scale!*(1 - 0.01 / sqrtn)
!    delta = 1

!    call trfmeansqp(fcn,t,s,hess,K,delta,tol,maxfev,quasi,verbose,info,msg)

    write(*,*) 'FINAL t:',real(t),' s =',real(s),' res=',real(fcn%val(t,s))
    write(*,*) 'info: ',trim(msg),' info=',info

    mean = t!x(1)
    scale = s! x(2)
    stderr = fcn%stderr(t,s,reli)

    if( .not. reli ) then
       info = 7
    end if

!    deallocate(res,psi)
!    if( .not. quasi ) deallocate(dpsi)

  end subroutine tsr1mean

  recursive subroutine trfmeansqp(fcn,t,s,hess,delta,rtol,maxfev,verb,iflag)

    use rfun
    use quantilmean
      use noise

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), intent(in out) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: hess
!    real(REAL64), intent(out) :: K
    real(REAL64), intent(in) :: delta, rtol
    integer, intent(in) :: maxfev
    logical, intent(in) :: verb
    integer, intent(out) :: iflag
!    character(len=*), intent(out) :: msg

!    real(REAL64), parameter :: rtol = 0.01

    integer, parameter :: maxiter = 3*precision(t)
    real(REAL64), parameter :: nulow = 0.01
    real(REAL64), parameter :: stol = 0.01
    real(REAL64), parameter :: dtol = 0.01

    integer :: i, iter, it, trinfo
    real(REAL64), dimension(2,2) :: B
    real(REAL64), dimension(2) :: p, fgrad, fg, y
    real(REAL64) :: fval, fv, dq, dlt, tp, sp, pnorm, smin, smax, tol, t0, s0, rho
    logical :: feasible, posdef, converge, suffred, intol
    character(len=80) :: fmt

    if( .not. (size(fcn%data) > 1 .and. s > 0) ) then
       iflag = 7
!       msg = 'An improper initial setup: Required N > 0, s > 0.'
       return
    end if

!    t = -2
!    s = 0.2
!    if( .not. (fcn%nfrac(tp,s) > nulow) ) then
!       iflag = 7
!       return
!    end if

    tol = rtol / sqrt(real(size(fcn%data)))
    s0 = s

    fval = huge(fval)
    fval = fcn%val(t,s)
    call fcn%grad(t,s,fgrad)
    if( fcn%break ) then
       iflag = 6
!       msg = 'Interrupt on request.'
       return
    end if
    ! the initial Hessian approximation
!    if( quasi ) then
    B = 0
    forall( i = 1:2 ) B(i,i) = 1
!    else
!       call fcn%hess(t,s,B)
!       if( fcn%break ) then
!          iflag = 6
!          msg = 'Interrupt on request.'
!          return
!       end if
!    end if
    posdef = .false.

    if( debug ) then
       write(*,*) 'trfmeansqp: t,s,tol,delta:',t,s,tol,delta
!       write(fmt,'("(i0,a,",i0,"f8.3,a,2es9.2,a,2es8.1,f5.2i3)")') 2
       fmt = "(i0,a,2f8.3,a,2es9.2,l2,a,2es8.1,f6.3,f5.2,i2,l2)"
       write(*,fmt) 0,':',t,s," >",fval,norm2(fgrad),posdef," |",delta
    end if

!    smax = maxval(abs(fcn%data - t))
    smax = s
    smax = 0.75
    smin = s

    !    call restorescale(fcn,t,smin,converge)
!    if( s > 0 ) then
!       smin = sminlim(fcn%data,t)
!       if( s < 2*smin ) call restorescale(fcn,t,s,converge)
!    end if

    dlt = 0.5*s
    rho = 0.5
    p = 0
    pnorm = 0

    do iter = 1, maxiter

!       dlt = min(delta,max(s,10*tol))
!       dlt = min(delta,0.9*s)
!       dlt = 0.5*s
       ! vyzkouset dlt = s - smin, kde smin je z fcn%nfrac(tp,sp) > 0.25

       if( rho > 0.75 ) then
          if( pnorm > 0.9*dlt ) then
             dlt = 2*dlt
          end if
       else if( 0.25 <= rho .and. rho <= 0.75 ) then
          dlt = dlt
       else ! if ( rho < 0.25 ) then
          dlt = dlt / 2
       end if
       dlt = min(dlt,0.9*s)
!       dlt = 0.9*s


       write(*,*) '>',dlt,tol,s,smin,rho
       do it = 1, 3*precision(t) !nint(3*log(dlt/tol)) + 1

          call trstep_new(B,fgrad,dlt,dtol,p,pnorm,posdef,trinfo)
          ! what about convergence?
          ! info = 0 ok
          ! info = 1 indefinite
          ! info = 6 detected no convergence

          if( verb ) call fcn%print(iter,t,s,fval,fgrad,B,p,dlt,dlt,trinfo)
          if( debug ) write(*,fmt) &
               it,'.',t,s," > ",fval,norm2(fgrad),posdef," |",dlt,pnorm, &
               fcn%nfrac(t+p(1),max(s+p(2),epsilon(s))),0.0,trinfo


          if( trinfo /= 1 ) then
             iflag = 666
 !            msg = 'Trust region step has failed.'
             return
          end if

          if( trinfo /= 6 ) then

             tp = t + p(1)
             sp = s + p(2)

             feasible = sp > 0 !.and. dlt > rtol*s
             if( .not. sp > 0 ) stop 'Trust region too big.'
             if( feasible ) then
                !          smin = rk4(fcn,t,smin,p(1))
             !          call restorescale(fcn,t,smin,converge)
             !             smin = sminlim(fcn%data,tp,sp)
!                write(*,*) tp,sp,smin
                write(*,*) 'sssss',tp,sp,smin,sp>smin,fcn%val(tp,sp), &
                     fcn%nfrac(tp,sp),sum(tukey((fcn%data-tp)/sp)**2), &
                     all(abs(fcn%data-tp) > 10*sp)
                !          if( sp > 2*smin ) then
!                feasible = fcn%nfrac(tp,sp) > 1.1 / size(fcn%data)
                feasible = fcn%nfrac(tp,sp) > nulow
!                feasible = .true.
!                feasible = sum(pshampel((fcn%data-tp)/sp)**2) > 0.1
!                feasible = sum(tukey((fcn%data-tp)/sp)**2) > 0
!                feasible = fcn%nfrac(tp,sp) > 0.01
!                feasible = feasible .and. all(abs(fcn%data-tp) > 10*sp)

!             feasible = fcn%val(t,s) > 0

             ! alternativa (ale vyzaduje poc. odhad mimo indefinite region):
             ! if( fcn%nfrac(tp,sp) < 0.25 ) then
             !    call hess()
             !    feasible = h(1,1)*h(2,2) - h(1,2)**2 > 0
             ! else
             !    feasible = .true.

!             feasible = sp > smin
!             feasible = posdef
!             feasible = .true.
             else
                smin = sminlim(fcn%data,t,s)
!             feasible = .false.
             end if
!          feasible = sp > 2*smin

             write(*,*) 'sol:',tp,sp,feasible

!          if( .not. feasible .and. dlt < rtol*s ) then
             if( .not. feasible ) then

                write(*,*) tp,sp,fcn%nfrac(tp,sp),fcn%val(t,s),recursive_depth
                !             stop 'unfeasible'


             !             call restorescale(fcn,t,sp,s-dlt,s+3*dlt,converge)
!                smin = 0.01
!                sp = qmed(abs(fcn%data - tp),iflag) / 0.6745
                smax = 1.41*sp*inverf(0.8_REAL64)
                smax = maxval(abs(fcn%data - tp))
!                smax = sp*(1 + 3.0/sqrt(real(size(fcn%data))))
!                write(*,*) 'RESTORE begins with:', real(t),real(smin),'..',real(smax)
                write(*,*) 'RESTORE begins with:', real(tp),real(s),real(sp)
                !                call restorescale_newer(fcn,t,smin,smax,sp,converge)
!                sp = s
                call restorescale_line(fcn,real([0,1],REAL64),tp,sp,stol,s0,dlt,converge)
                if( .not. converge ) then
                   iflag = 666
                   !            msg = 'Restore procedure failed.'
                   return
                end if
                write(*,*) 'RESTORE() finished: ',real(tp),real(sp),converge
                if( .not. converge ) then
!                   sp = qmed(abs(fcn%data - tp),iflag) / 0.6745
 !                  sp = 1.41*sp*inverf(0.8_REAL64)
                end if
                smin = sp
!                stop 'restore'
!             call restorescale_new(fcn,t,smin,smax,sp,converge)
             !             if( .not. converge ) stop 'restorescale()'

!                sp = qmed(abs(fcn%data - tp),iflag) / 0.6745
!                sp = 0.8*sp
!                if( iflag /= 0 ) then
!                   iflag = 6
!                   return
!                end if

!                exit
                recursive_depth = recursive_depth + 1
                if( .false. ) then
!                if( recursive_depth < recursive_depth_limit ) then
                   t = tp
                   s = sp
                   write(*,*) 'RESTORE start at:', real(t),real(s)
                   block
                     real(REAL64) :: tol, delta
                     tol = 0.01 * s / sqrt(real(size(fcn%data)))
                     delta = 0.5*s
!                     call trfmeansqp(fcn,t,s,B,delta,tol,maxfev,quasi,verb,iflag)
                   end block
!                   call tropt(fcn,t,s,B,verb,.false.,iflag)
                   write(*,*) 'RESTORE tropt() finished: ',real(t),real(s),iflag
                   hess = B
                else
                   iflag = 6
                end if
!                stop
!                return

!                tp = t
                !             sp = smin
                fv = fcn%val(tp,sp)
                call fcn%grad(tp,sp,fg)
                if( fcn%break ) then
                   iflag = 6
                   return
                end if
                !                call fcn%hess(t,sp,B)
                t = tp
                s = sp
                B = 0
                forall( i = 1:2 ) B(i,i) = 1
                dlt = 0.5*s
                rho = 0
!                write(*,*) 'BBBBBBBBBB',B

                exit

             else

                fv = fcn%val(tp,sp)
                call fcn%grad(tp,sp,fg)
                if( fcn%break ) then
                   iflag = 6
!                msg = 'Interrupt on request.'
                   return
                end if

                dq = -dot_product(fgrad,p) - 0.5*dot_product(p,matmul(B,p))
                rho = -(fval - fv) / (dot_product(fgrad,p) + 0.5*dot_product(p,matmul(B,p)))

                write(*,*) '>>',dq,fval-fv,dot_product(fgrad,p),0.5*dot_product(p,matmul(B,p)),rho

                y = fg - fgrad
                call sr1update(p,y,pnorm,B)

                suffred = dq > 0 .and. fval > fv + 0.1*dq .and. 0.75 < rho .and. rho < 1.25
                !             if( dq > 0 .and. fval - fv > 0.1*dq .and. fv <= fval ) then
                suffred = 0.75 < rho .and. rho < 1.25
                if( suffred ) then
                   exit
                else
                   !                dlt = max(2.1*tol,dlt/2)
                   dlt = dlt / 2
                end if

             end if

          else
             dlt = dlt / 2
          end if

       end do

!       rho = -(fval - fv) / (dot_product(fgrad,p) + 0.5*dot_product(p,matmul(B,p)))

!       y = fg - fgrad
       t = tp
       s = sp
       fval = fv
       fgrad = fg

       if( debug ) write(*,fmt) &
            iter,': ',t,s," > ",fval,norm2(fgrad),posdef," |",dlt,pnorm,fcn%nfrac(t,s),rho,0,s

       intol = abs(p(1)) < tol*s .and. abs(p(2)) < tol*s
       intol = norm2(p) < tol*s
       if( intol .and. posdef .and. iter > 1 ) then
!       if( norm2(p) < tol .and. posdef .and. iter > 1 ) then
!          iflag = 1
          exit
       end if

!       call fcn%grad(t,s,fg)
!       if( fcn%break ) then
!          iflag = 6
!          return
!       end if

!       if( quasi ) then
!          call sr1update(p,y,norm2(p),B)
!       else
!          call fcn%hess(t,s,B)
!          if( fcn%break ) then
!             iflag = 6
!             return
!          end if
!       end if


    end do

    if( debug ) write(*,fmt) &
         666,': ',t,s," > ",fval,norm2(fgrad),posdef," |",dlt,pnorm,1.0,1.0,0,suffred

    write(*,*) 'hess:',real(B)

    hess = B
!    K = (hess(1,1)*hess(2,2) - hess(1,2)**2) / (1 + sum(fgrad**2))**2
!    write(*,*) 'K:',K
    iflag = 0
!    msg = ''

!!$    block
!!$      use rfun
!!$      real(REAL64) :: s0,df
!!$
!!$      s0 = sum(tukey((fcn%data-t)/s)**2) /  size(fcn%data)
!!$      df = sum(dtukey((fcn%data-t)/s)) /  size(fcn%data)
!!$      write(*,*) "stderr:", s*sqrt(s0/B(1,1)**2),s0,df
!!$      write(*,*) "stderr':",s*sqrt(fval/B(1,1)/size(fcn%data)), &
!!$         'stddev:',s*sqrt(fval/B(1,1))
!!$      write(*,*) 'hess:',real(B)
!!$    end block


  end subroutine trfmeansqp

  real(REAL64) function sminlim(data,t,s)

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t,s
    real(REAL64), parameter :: eps = 1e3*epsilon(sminlim)
    real(REAL64) :: smin,a,b,u
    real(REAL64), dimension(2) :: w
    integer :: i,n
    logical :: found

    smin = huge(smin)
    do i = 1, size(data)
       smin = min(max(abs(data(i) - t),eps),smin)
    end do
    sminlim = smin / 6

    smin = s
    a = t - 6*s
    b = t + 6*s
    found = .false.
    n = 0
    do i = 1, size(data)
       if( a < data(i) .and. data(i) < b ) then
          n = n + 1
          w(n) = data(i)
       end if
       if( n == 2 ) then
          found = .true.
          u = (w(1) + w(2))/2
          smin = (maxval(w) - u) / 6
          exit
       end if
    end do

  end function sminlim

  subroutine restorescale_line(fcn,qvec,t,s,rtol,smax,delta,converge)

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), dimension(:), intent(in) :: qvec
    real(REAL64), intent(in) :: t,rtol,delta,smax
    real(REAL64), intent(in out) :: s
    logical, intent(out) :: converge

    type(zeroin_data) :: fraczero
    real(REAL64), dimension(2) :: fg, fgrad, x, p, xp, y
    real(REAL64), dimension(2,2) :: B
    real(REAL64) :: g,pHp,d,ared, pred, dlt, rho, fv, fval, sp, s0
    integer :: i,iter,iflag
    logical :: posdef

    converge = .false.
    dlt = delta
    s0 = s

    goto 33
    fval = fcn%val(t,s)
    call fcn%grad(t,s,fgrad)
    B = 0
    forall( i = 1:2 ) B(i,i) = 1

    x = [ t, s ]
    do iter = 1, 3*precision(x)

       g = dot_product(fgrad,qvec)
       pHp = dot_product(qvec,matmul(B,qvec))

       if( pHp < 0 ) then
          d = -g/pHp
          if( abs(d) > dlt ) d = sign(dlt,g)
       else
          d = sign(dlt,g)
       end if
       p = d*qvec

!       call trstep_new(B,fgrad,dlt,rtol,p,pnorm,posdef,trinfo)
       xp = x + p

       if( xp(2) > 0 ) then
          fv = fcn%val(xp(1),xp(2))
          call fcn%grad(xp(1),xp(2),fg)

          y = fg - fgrad
          ared = fval -  fv
          pred = -(dot_product(fgrad,p) + 0.5*dot_product(p,matmul(B,p)))
!          write(*,*) 'rho:',fv,fval,ared,pred, ared/pred, dot_product(fvec,p),dot_product(p,matmul(B,p))

          if( abs(ared) < epsilon(ared) .and. abs(pred) < epsilon(pred) ) then
             rho = 1
          else if( abs(pred) < epsilon(pred) .and. abs(ared) > 0 ) then
             rho = 0
          else
             rho = ared / pred
          end if
       else
          y = 0
          rho = 0
       end if

       if( rho > eta ) then
          x = xp
          fval = fv
          fgrad = fg
       else
          x = x
       end if

       if( debug ) &
       write(*,'(a,i2,a,2f8.3,":",2f7.3,2es10.2," |",f6.2,f6.2,l2)') &
            "restore ",iter,": ",x(2),fval,d,g,pHp,dlt,rho
!            "restore ",iter,": ",x,fval,d,g,pHp,dlt,rho,posdef

       converge = abs(d) < rtol*x(2) .and. pHp < 0
!       converge = pnorm < rtol*norm2(x) .and. posdef
       if( converge ) then
          s = x(2)
          exit
       end if

       if( rho > 0.75 ) then
          if( abs(d) > 0.9*dlt ) then
             dlt = 2*dlt
          end if
       else if( 0.25 <= rho .and. rho <= 0.75 ) then
          dlt = dlt
       else ! if ( rho < 0.25 ) then
          dlt = dlt / 2
       end if
       dlt = min(dlt,0.9*x(2))

!       write(*,*) p,y,B
       call sr1update(p,y,norm2(p),B)
    end do

!    if( .not. converge ) then
!       converge = .true.
!       s = s0
!    end if

33  continue
    converge = .true.
    s = s0
    if( converge ) then
!       converge = .false.
!!$       goto 10
!!$       do i = 1, int(4*smax/dlt)+1
!!$          sp = s + i*(dlt/4)
!!$!                      write(*,*) s,fcn%nfrac(t,s),fcn%nfrac(tp,s)
!!$          if( fcn%nfrac(t,sp) > 0.1 ) then
!!$             converge = .true.
!!$             s = sp
!!$             exit
!!$          end if
!!$          !                sp = sp + dlt
!!$       end do
!!$10     continue
       !       return

       fraczero%funmean => fcn
       fraczero%t = t
!       s = s0
       write(*,*) 'zeroin on interval: ',s,smax,' tol: ',rtol*smax
       if( .not. s < smax ) stop 's < smax'
       sp = zeroin(fraczero,s,smax,rtol*smax,iflag)
       converge = iflag == 0
 !      write(*,*) 'sp=',sp

!       stop
       !       call zin(fcn,t,s,smax,rtol,sp)
       if( converge ) then
          s = sp
!       write(*,*) iflag
          !       if( iflag /= 0 ) then
       else
          write(error_unit,*) 'ERROR in sr1mean(): zeroin() has failed; flag=',iflag
!          s = smin
       end if
!!$       if ( iflag == 3 ) then
!!$          do i = 1, 0*(int(4*smax/dlt)+1)
!!$             sp = s + i*(dlt/4)
!!$!                      write(*,*) s,fcn%nfrac(t,s),fcn%nfrac(tp,s)
!!$             if( fcn%nfrac(t,sp) > 0.1 ) then
!!$                converge = .true.
!!$                s = sp
!!$                iflag = 0
!!$                exit
!!$             end if
!!$             !                sp = sp + dlt
!!$          end do
!!$       end if
!       converge = iflag == 0
    end if
!    stop

  end subroutine restorescale_line

!!$  subroutine zin(fcn,t,s,smax,rtol,sp)
!!$
!!$    use fmm
!!$
!!$    class(tsr1mean_function), pointer, intent(in) :: fcn
!!$    real(REAL64), intent(in) :: t,s,smax,rtol
!!$    real(REAL64), intent(out) :: sp
!!$
!!$    smin_fcn => fcn
!!$    smin_t = t
!!$    smin_count = 0
!!$    sp = zeroin(s,smax,fraczero1,rtol*smax)
!!$
!!$  end subroutine zin

  real(REAL64) function fraczero(this,x)

    class(zeroin_data) :: this
    real(REAL64), intent(in) :: x

    fraczero = this%funmean%nfrac(this%t,x) - 0.1
!    write(*,*) x,fraczero

  end function fraczero

!!$  real(REAL64) function fraczero1(s)
!!$
!!$    real(REAL64), intent(in) :: s
!!$
!!$!     = this%funmean%nfrac(this%t,x) - 0.1
!!$    fraczero1 = smin_fcn%nfrac(smin_t,s) - 0.1
!!$    write(*,*) s,fraczero1
!!$
!!$  end function fraczero1

  subroutine restorescale_newer(fcn,t,smin,smax,s,converge)

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), intent(in) :: t,smin,smax
    real(REAL64), intent(in out) :: s
    logical, intent(out) :: converge

    real(REAL64), dimension(2), parameter :: q = [0,1]
    real(REAL64), parameter :: z = -1
    real(REAL64) :: d, gtol, u, fmin
    logical :: reli

    d = smax - smin
    gtol = 0.05
    call goldfinger(fcn,z,[t,smin],q,d,gtol,u,fmin,reli)

    converge = 1.5*gtol < u .and. u < 1 - 1.5*gtol .and. reli

    if( converge ) s = smin + u*d

  end subroutine restorescale_newer

  subroutine restorescale_new(fcn,t,smin,smax,s,converge)

    use fmm
    use quantilmean

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), intent(in) :: t,smin,smax
    real(REAL64), intent(in out) :: s
    logical, intent(out) :: converge

    real(REAL64), dimension(:), allocatable :: x,y
    real(REAL64) :: s1,s2,K1,K2,tol,sopt
    real(REAL64), dimension(2,2) :: hess
    integer :: n,i

    converge = .true.

    n = size(fcn%data)
    allocate(x(n),y(n))
    call ecdf(abs(fcn%data - t),x,y)
    s = x(n/2+1) / 0.6745
    s = 0.8*s
    write(*,*) 'restorescale MAD():',s
    write(*,*) real(x)

!      tp = t
      !      call ridge(funmean,tp,s,nopt,topt,sopt)

    call fcn%hess(t,s,hess)
    write(*,*) 'hess:',hess

    if( hess(2,2) > 0 ) return

    ! find minimum in s: sopt(1) ... x(n)
    ! find minimum in s: prvni nenulove x ... sopt(1)

    iscale_t = t
    iscale_fcn => fcn
    tol = 1e-2*s
    sopt = s

    ! IMPLEMENTACE fmin v f08
    s = fmin(sopt,x(n),scaleline,tol)
    if( s > 0 ) then
       call fcn%hess(t,s,hess)
       K1 = hess(1,1)*hess(2,2) - hess(1,2)**2
       if( hess(2,2) <= 0 ) K1 = 0
       s1 = s
       write(*,*) 's1:',s1,K1
    end if

    do i = 1, n
       if( abs(x(i)) > epsilon(x) ) then
          s = fmin(x(i),sopt,scaleline,tol)
          if( s > 0 ) then
             call fcn%hess(t,s,hess)
             K2 = hess(1,1)*hess(2,2) - hess(1,2)**2
             if( hess(2,2) <= 0 ) K2 = 0
             s2 = s
             write(*,*) 's2:',s2,K2
             exit
          end if
       end if
    end do

    if( K2 > K1 ) then
       s = s2
    else
       s = s1
    end if

  end subroutine restorescale_new

  subroutine restorescale(fcn,t,smin,smax,s,converge)

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), intent(in) :: t,smin,smax
    real(REAL64), intent(in out) :: s
    logical, intent(out) :: converge

    type(zeroin_data) :: funzero
    real(REAL64) :: tol
    integer :: iflag

    funzero%funmean => fcn
    funzero%t = t

    tol = 1e-3*smax
    write(*,*) smin,smax,tol

    block
      integer :: i
      real(REAL64) :: dt,f
      dt = (smax - smin) / 100
      open(55,file='/tmp/xyx')
      do i = 1,100
!         f = sum(abs(fcn%data - t)
         write(55,*) i*dt,fcn%nfrac(t,i*dt) - 0.5
      end do
      close(55)
    end block

    s = zeroin(funzero,smin,smax,tol,iflag)
    write(*,*) 'restorescale: zeroin(): ',s,iflag
    if( iflag /= 0 ) then
       write(*,*) iflag
       write(error_unit,*) 'ERROR in sr1mean(): zeroin() has failed.'
       s = smin
    end if
    converge = iflag == 0
!      stop


  end subroutine restorescale

  real(REAL64) function scaleline(x)

    real(REAL64), intent(in) :: x

    scaleline = iscale_fcn%val(iscale_t,x)

  end function scaleline


  real(REAL64) function zeroin_fun(this,x)

    class(zeroin_data) :: this
    real(REAL64), intent(in) :: x

    zeroin_fun = this%funmean%nfrac(this%t,x) - 0.25

  end function zeroin_fun

  subroutine restorescale_old(fcn,t,s,converge)

!    use trust_global_scale
    use quantilmean

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), intent(in) :: t!,smin,smax
    real(REAL64), intent(in out) :: s
    logical, intent(out) :: converge

    real(REAL64), dimension(:), allocatable :: xcdf,ycdf
    real(REAL64), dimension(2), parameter :: qvec = [0, 1]
    real(REAL64), dimension(2,2) :: hess,Q
    real(REAL64), dimension(2) :: lam
    real(REAL64) :: smin,smax,sig,tol,tt,delta,ss,h22
!    real(REAL64) :: sig
    integer :: n,m,i,ieig


    n = size(fcn%data)
    allocate(xcdf(n),ycdf(n))
    call ecdf(abs(fcn%data-t),xcdf,ycdf)

    m = n - int((1-0.5) * n)
    s = xcdf(m) / 0.6745
    write(*,*) 's=',s,m

    if( n > 13 ) then
       m = int(0.95 * (n+1))
       smax = xcdf(m)
       smin = xcdf(n + 1 - m)
    else
       smin = xcdf(2)
       smax = xcdf(n)
    end if

!       smin = minval(abs(fcn%data - t))/6
    do i = 1, n
       if( xcdf(i) > 1e3*epsilon(xcdf) ) then
          smin = xcdf(i) / (0.9*6)
          exit
       end if
    end do

!    smin = huge(smin)
!    do i = 1, n
!       r = abs(fcn%data(i) - t)
!       if( r > 1e3*epsilon(r) ) smin = min(r,smin)
!    end do

    smin = sminlim(fcn%data,t,s)

    sig = s / sqrt(real(n))

!    call trgscale(fcn,t,smin,smax,sig,s,converge)
!    write(*,*) 'by trgscale:',s,sig,converge

    !    s = sminupdate(fcn,t,s)

    tt = t
    s = smin
    tol = 0.1*smin
    delta = s / 2
    call smnewton(fcn,qvec,tt,s,tol,delta,converge)
!    call fcn%hess(tt,s,hess)
!    call q2eig(hess,lam,Q)
!    write(*,*) lam
!    write(*,*) Q
!    ieig = 0
!    do i = 1, 2
!       if( lam(i) < 0 ) ieig = i
!    end do
!    if( ieig == 0 ) stop 'Something wrong.....'
    do i = 0, 16
       tt = t + 0.5*i*s
       ss = s + 0.5*i*s
       call fcn%hess(tt,ss,hess)
       h22 = hess(1,1)*hess(2,2) - hess(1,2)**2
       write(*,*) tt,ss,hess(1,1),h22,hess(1,1)>0 .and. h22 > 0
!       if( hess(1,1) > 0 .and. h22 > 0 ) then
       if( hess(1,1) > 0 ) then
          s = ss
          exit
       end if
    end do
!    s = s + 10*tol
!    s = 1.001*s + delta

!    if( abs(t - tt) > 1e3*epsilon(t)) stop '|t-tt| fail in restorescale.'

!    converge = .true.
    write(*,*) 'Restored s:',s
!    stop 'restored'

  end subroutine restorescale_old



  subroutine tsr1meansqp(fcn,t,s,hess,delta,tol,maxfev,quasi,verb,iflag)

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), intent(in out) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: hess
    real(REAL64), intent(in) :: delta, tol
    integer, intent(in) :: maxfev
    logical, intent(in) :: quasi, verb
    integer, intent(out) :: iflag
!    character(len=*), intent(out) :: msg

    real(REAL64), parameter :: rtol = 0.01
    real(REAL64), parameter :: atol = 0
    real(REAL64), parameter :: dltmin = 1e3*epsilon(delta)

    real(REAL64), dimension(2,2) :: B
    real(REAL64), dimension(2) :: fgrad, fg, y, p
    real(REAL64) :: dlt, fval, fv, ared, pred, rho, pnorm, smin, smin0, tp, sp!, sknot
    integer :: i,iter,trinfo
    logical :: posdef
    character(len=80) :: fmt

    iflag = 666
    smin = sminupdate(fcn,t,s)
    write(*,*) 'smin=',real(smin)

!    s = iscale(fcn,t)

    if( s < 1.1*smin ) then
       write(*,*) 'under smin init: ',real(s),' < ',real(smin)
       s = 1.1*smin
    end if

!    dlt = delta

    ! minimal scale parameter
!    smin = 0.1*x(n)/sqrt(real(666.0))
!    smin = 0.5*x(n)/sqrt(real(n))
!    smin = 0

    dlt = min(delta,0.1*(s - smin))
    write(*,*) delta,s,smin
!    write(*,*) 'ssss',delta,0.5*(s - smin),dlt,smin

    ! the initial Hessian approximation
    if( quasi ) then
       B = 0
       forall( i = 1:2 ) B(i,i) = real(1,REAL64)
!       B(1,1) = 1 / s**2
!       B(2,2) = 2 / s**2
!       call fcn%hessian(x,B)
    else
       call fcn%hess(t,s,B)
    end if

    ! initial setup
    fval = fcn%val(t,s)
    call fcn%grad(t,s,fgrad)
!    call fvalgrad(t,s,fval,fgrad)

    posdef = .false.

    if( verb ) then
       p = real(0,REAL64)
       rho = 0
       call fcn%print(0,t,s,fval,fgrad,B,p,dlt,rho,0)
!       write(2,*) 0,t,s,fval,dlt,p
    end if

    if( debug ) then
!       write(fmt,'("(i0,a,",i0,"f8.3,a,2es9.2,a,2es8.1,f5.2i3)")') 2
       fmt = "(i0,a,2f8.3,a,2es9.2,l2,a,2es8.1,f5.2,i3)"
       write(*,fmt) 0,':',t,s," >",fval,norm2(fgrad),posdef," |",dlt
    end if

    do iter = 1, maxfev

       if( norm2(fgrad) < tol .and. posdef ) then
          iflag = 1
          exit
       end if

!       B(1,:) = [ 5.743E-01,      1.385E+00]
!       B(2,:) = [1.385E+00,      2.471E+00]
       !       fgrad = [-8.668E-01  ,   -4.070E-01]
!       fgrad = 0
!       call trstep(B,fgrad,dlt,rtol,atol,p,posdef,trinfo)
       call trstep_new(B,fgrad,dlt,rtol,p,pnorm,posdef,trinfo)
       ! what about convergnece?
       ! info = 0 ok
       ! info = 1 indefinite
       ! info = 3 detected no convergence

       if( verb ) call fcn%print(iter,t,s,fval,fgrad,B,p,dlt,rho,trinfo)
       if( debug ) write(*,fmt) &
            iter,': ',t,s," > ",fval,norm2(fgrad),posdef," |",dlt,norm2(p),rho,trinfo

       if( trinfo /= 1 ) then
          iflag = 666
!          exit
       end if!stop 'trstep'

       ! attempt the step
       !       xp = x + p
       tp = t + p(1)
       sp = s + p(2)

!!$       block
!!$         integer :: i,j
!!$         real(REAL64) :: s
!!$
!!$         sknot = 0
!!$         do i = 1, size(fcn%data)
!!$            do j = i+1, size(fcn%data)
!!$               if( t-dlt < fcn%data(i) .and. fcn%data(i) < t+dlt ) then
!!$                  s = (fcn%data(i) + fcn%data(j))/(2*3)
!!$                  if( s > sknot ) sknot = s
!!$               end if
!!$            end do
!!$         end do
!!$       end block
!!$       write(*,*) 'smin=',real(sknot)

       ! update smin, once in a while
       !       if( mod(iter,10) == 0 ) then
       if( abs(p(1)) > 0.75*dlt ) then
          smin0 = smin
          smin = sminupdate(fcn,tp,sp)
       end if

       if( sp > smin .and. trinfo == 1 ) then
!       if( sp > sknot .and. trinfo == 1 ) then

!          call fvalgrad(tp,sp,fv,fg)
          fv = fcn%val(tp,sp)
          if( fcn%break ) then
             iflag = 6
             exit
          end if
          call fcn%grad(tp,sp,fg)
          if( fcn%break ) then
             iflag = 6
             exit
          end if

          y = fg - fgrad
          ared = fval - fv
          pred = -(dot_product(fgrad,p) + 0.5*dot_product(p,matmul(B,p)))
!          pred = -(dot_product(fg(?),p) + 0.5*dot_product(p,matmul(B,p)))

          if( abs(ared) < epsilon(ared) .and. abs(pred) < epsilon(pred) ) then
             rho = 1
          else if( abs(pred) < epsilon(pred) .and. abs(ared) > 0 ) then
             rho = 0
          else
             rho = ared / pred
          end if

          write(*,*) 'ared,pred:',real(ared),real(pred),' rho=',real(rho)

       else ! s <= 0

          rho = 0

       end if

       if( rho > eta ) then
          t = tp
          s = sp
          fval = fv
          fgrad = fg
       else
          t = t
          s = s
          smin = smin0
       end if

       pnorm = norm2(p)
       if( rho > 0.75 ) then
          if( pnorm  > 0.9*dlt ) then
             dlt = 2*dlt
          end if
       else if( 0.25 <= rho .and. rho <= 0.75 ) then
          dlt = dlt
       else ! if ( rho < 0.25 ) then
          dlt = dlt / 2;
       end if
!       dlt = min(dlt,s-smin);
       dlt = min(dlt,s/2);
       !       write(*,*) '****',dlt,x(n),pnorm,rho
       if( dlt < dltmin ) then
          write(*,*) 'delta under min:',dlt,dltmin
!       if( dlt < 1e-3*tol ) then
          iflag = 5
          exit
       end if

       if( quasi ) then
          call sr1update(p,y,pnorm,B)
       else
          call fcn%hess(t,s,B)
       end if

       if( iter == maxfev ) iflag = 5
    end do

    hess = B
!    if( iflag == 1 ) then
!       msg = 'The optimum is localised with the required accuracy.'
!    else if( iflag == 5 ) then
!       msg = 'Maximum number of function evaluations reached.'
!    else if( iflag == 6 ) then
!       msg = 'This optimalisation interrupted on request.'
!    else
!       msg = ''
!    end if

  end subroutine tsr1meansqp


!!$  function fvalue(t,s)
!!$
!!$    real(REAL64) :: fval
!!$    real(REAL64), intent(in) :: t,s
!!$
!!$    res = (data - t)/s
!!$    fval = sum(ifcn(res))/n + log(s)
!!$
!!$  end function fvalue
!!$
!!$  subroutine fgradient(t,s,grad)
!!$
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64), dimension(:), intent(out) :: grad
!!$
!!$    res = (data - t)/s
!!$    psi = fcn(res)
!!$
!!$    grad(1) = -sum(psi)/(n*s)
!!$    grad(2) = -sum(psi*res)/(n*s) + real(1.0,REAL64)/s
!!$
!!$  end subroutine fgradient
!!$
!!$  subroutine fvalgrad(t,s,fval,grad)
!!$
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64), intent(out) :: fval
!!$    real(REAL64), dimension(:), intent(out) :: grad
!!$
!!$    res = (data - t)/s
!!$    psi = fcn(res)
!!$
!!$    fval = sum(ifcn(res))/n + log(s)
!!$    grad(1) = -sum(psi)/(n*s)
!!$    grad(2) = -sum(psi*res)/(n*s) + real(1.0,REAL64)/s
!!$
!!$  end subroutine fvalgrad
!!$
!!$  subroutine fhessian(t,s,hess)
!!$
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64), dimension(:,:), intent(out) :: hess
!!$
!!$    res = (data - t)/s
!!$    psi = fcn(res)
!!$    dpsi = dfcn(res)
!!$
!!$    hess(1,1) = sum(dpsi) / s**2 / n
!!$    hess(1,2) = (sum(dpsi*res)/s**2 + sum(psi)/s**2) / n
!!$    hess(2,1) = fhess(1,2)
!!$    hess(2,2) = (sum(dpsi*res**2) + 2*sum(psi*res))/s**2/n - real(1.0,REAL64)/s**2
!!$
!!$  end subroutine fhessian


  subroutine trstep_new(B,g,delta,rtol,p,pnorm,posdef,info)

    real(REAL64), dimension(:,:), intent(in) :: B
    real(REAL64), dimension(:), intent(in) :: g
    real(REAL64), intent(in) :: delta, rtol
    real(REAL64), dimension(:), intent(out) :: p
    real(REAL64), intent(out) :: pnorm
    logical, intent(out) :: posdef
    integer, intent(out) :: info

    real(REAL64), parameter :: atol = 1e-2
    real(REAL64), dimension(2,2) :: Q
    real(REAL64), dimension(2) :: lam,Qg
!    logical, dimension(2) :: Qortho
    real(REAL64) :: pQ, gnorm, p2, t, tau!, atol2
    integer :: iflag

    write(*,*) 'trstep init: delta=',real(delta),' g=',real(g)
    write(*,*) 'B:',real(B)


    info = 666

    ! determine eigen values and vectors
    call q2eig(B,lam,Q)

    Qg = matmul(transpose(Q),g)
!    write(*,*) 'Qg:',Qg
!    write(*,*) 'lam:',lam

    ! Does is B positive definite?
    posdef = all(lam > 0)

    if( posdef ) then

!       par = real(0,REAL64)
       !       call eqsol_new(B,lam,Q,Qg,par,p,pnorm)
       p = -matmul(Q,Qg/lam)
       pnorm = norm2(p)
!       write(*,*) 'sol:',real(p),' norm:',real(pnorm)
       write(*,'(a,2es10.2,a,es15.2)') 'sol:',p,' norm:',pnorm

       ! par = 0: |p| <= delta
       if( pnorm < (1+rtol)*delta ) then
          info = 1
          return
       end if
    end if

!    atol2 = atol**2/2
    gnorm = norm2(g)

    if( gnorm < epsilon(gnorm) ) then
       p = 0
       pnorm = 0
       info = 1
       return
    end if

!    Qortho = abs(1 - abs(Qg)/gnorm) < atol2
!    if( any(abs(1 - abs(Qg)/gnorm) < atol2) ) then
    if( any(abs(Qg) < atol*gnorm) ) then
!    if( abs(1 - abs(dot_product(Q(:,1),g))/gnorm) < atol**2/2 ) then
!    if( abs(1 - abs(Qg(1))/gnorm) < atol**2/2 ) then
       p = -delta*(g/gnorm)
       pnorm = norm2(p)
       write(*,*) 'Steepest descent step:',real(p),':',real(Qg),real(gnorm), &
            abs(1 - abs(Qg)/gnorm) < atol**2/2,abs(Qg)/gnorm,abs(Qg) < atol*gnorm
       if( .not. any(abs(Qg) < atol*gnorm) ) stop
!       stop
       info = 1
!    else if( abs(1 - abs(dot_product(Q(:,2),g))/gnorm) < atol**2/2 ) then
!    else if( abs(1 - abs(Qg(2))/gnorm) < atol**2/2 ) then
!       p = -delta*(g/gnorm)
!       write(*,*) '------------ easy step: 2',real(p),dot_product(Q(:,2),g),gnorm
!       stop
!       info = 1
    else

!!$    if( posdef ) then
!!$       ! |p| > delta, par should be adjusted
!!$       par = max(norm2(g) / delta - norm2(lam),real(0,REAL64))
!!$    else
!!$       ! indefinite
!!$       par = max(norm2(g) / delta - norm2(lam),-1.001*minval(lam))
!!$    end if

       call bpar_new(lam,posdef,Q,Qg,delta,rtol,p,pnorm,iflag)
!       write(*,*) "1:",p
!       p = -delta*(g/gnorm)
!       write(*,*) "2:",-delta*(g/gnorm)
!       write(*,*) 'par=',par, iflag
       if( iflag == 0 ) then
!          call eqsol_new(B,lam,Q,Qg,par,p,pnorm)
          if( posdef ) then
             write(*,*) 'outside sol:',real(p),' norm:',real(pnorm)
          else
             write(*,*) 'indefinite sol:',real(p),' norm:',real(pnorm)
          end if
          info = 1
       else
          if( posdef ) then
             write(*,*) 'outside sol: NO CONVERGENCE'
          else
             write(*,*) 'indefinite sol: NO CONVERGENCE'
          end if
          info = 6
       end if
    end if

!!$    else ! indefinite
!!$
!!$       par = -lam(1)
!!$!       call bpar_new(lam(2:2),Qg(2:2),0.5*delta,rtol,par,iflag)
!!$       !       write(*,*) 'par=',par, iflag
!!$       iflag = 0
!!$       if( iflag == 0 ) then
!!$
!!$          p = -Q(:,2)*Qg(2)/(lam(2) + par)
!!$          p = 0.5*delta*p / norm2(p)
!!$          pQ = dot_product(p,Q(:,1))
!!$
!!$          p2 = Qg(2)/(lam(2) + par)
!!$          write(*,*) 'indef p-part:',real(p2*Q(:,2))
!!$          t = (delta - abs(pQ))*((delta + abs(pQ))/delta)
!!$          !          write(*,*) 'pQ,t,tau:',real(pQ),real(t),real(tau)
!!$          tau = sign(t,pQ)/(abs(pQ) + sqrt(pQ**2 + t/delta))
!!$          !          tau = sign(tau,pz)
!!$          !          p2 = dot_product(p,p)
!!$          if( delta > abs(p2) ) then
!!$!             tau = sign(sqrt(delta**2 - p2**2),Qg(1))
!!$             tau = sqrt(delta**2 - p2**2)
!!$          else
!!$             tau = 0
!!$!             p2 = 0.5*delta
!!$ !            tau = sign(sqrt(delta**2 - p2**2),Qg(1))
!!$          end if
!!$          write(*,*) 'Qg(1),tau:',real(Qg(1)),real(tau),real(tau/delta)
!!$          p = -tau*Q(:,1) - p2*Q(:,2)
!!$          p = sign(delta,-Qg(1))*Q(:,1)
!!$          pnorm = norm2(p)
!!$          write(*,*) 'indef final:',real(p),' norm:',real(pnorm)
!!$!          stop
!!$       else
!!$          write(*,*) 'indefinite sol: NO CONVERGENCE'
!!$          info = 6
!!$       end if
!!$!       stop 'indefinite'
!!$
!!$    end if

!!$    end if
!!$
!!$       if( abs(Qg(1)) < -0.1 ) then
!!$          p = -Q(:,2)*Qg(2)/(lam(2) - lam(1))
!!$          write(*,*) 'hard plain:',real(p),'Qg(1):',real(Qg(1))
!!$          pQ = dot_product(p,Q(:,1))
!!$          t = (delta - abs(pQ))*((delta + abs(pQ))/delta)
!!$
!!$          !          write(*,*) 'pQ,t,tau:',real(pQ),real(t),real(tau)
!!$          tau = sign(t,pQ)/(abs(pQ) + sqrt(pQ**2 + t/delta))
!!$          !          tau = sign(tau,pz)
!!$          p2 = dot_product(p,p)
!!$          tau = sqrt(delta**2 - p2)
!!$          write(*,*) 'pQ,t,tau:',real(pQ),real(t),real(tau),real(tau/delta)
!!$          p = p + tau*Q(:,1)
!!$          pnorm = norm2(p)
!!$          write(*,*) 'hard final:',real(p),' norm:',real(pnorm)
!!$!          stop 'hard'
!!$       else
!!$
!!$
!!$!       Qg(1) = 0
!!$!       Q(:,1) = 0
!!$!       par = max(-lam(1) + 0.1*abs(1+lam(1)),norm2(g)/delta-norm2(lam))
!!$!       par = -lam(1) + 0.001*abs(lam(1))
!!$          par = -1.001*lam(1)
!!$          call bpar_new(lam,Qg,delta,rtol,par,iflag)
!!$          write(*,*) 'par=',par, iflag
!!$          if( iflag == 0 ) then
!!$             call eqsol_new(B,lam,Q,Qg,par,p,pnorm)
!!$             write(*,*) 'indefinite sol:',real(p),' norm:',real(pnorm)
!!$             info = 1
!!$
!!$             block
!!$               use cholesky
!!$               real(REAL64), dimension(2,2) :: B1,R
!!$               real(REAL64), dimension(2) :: p
!!$               integer :: i, indef
!!$
!!$               B1 = B
!!$               forall( i = 1:2 ) B1(i,i) = B(i,i) + par!- 1.001*lam(1)
!!$               call chol_rrfac(B1,R,indef)
!!$               call chol_rrsolve(R,-g,p)
!!$               !         call chol_rrtsolve(R,p,q)
!!$               !         write(*,*) 'chol:',real(norm2(p)),real(norm2(q))
!!$               write(*,*) 'R:',real(R(1,:))
!!$               write(*,*) 'R:',real(R(2,:))
!!$               write(*,*) 'Rp:',norm2(matmul(R,p))**2
!!$               write(*,*) 'Rz:',norm2(matmul(R,Q(:,1)))**2
!!$               write(*,*) par*delta**2
!!$             end block
!!$
!!$       else if( iflag == -1 ) then
!!$          ! the hard case
!!$          par = lam(1)
!!$          call eqsol_new(B,lam,Q,Qg,par,p,pnorm)
!!$          write(*,*) 'p,pnorm:',p,pnorm
!!$          if( pnorm < delta ) then
!!$             pQ = dot_product(p,Q(:,1)) / delta
!!$             t = (delta - pnorm)*((delta + pnorm)/delta)
!!$             !          write(*,*) 'pQ,t,tau:',real(pQ),real(t),real(tau)
!!$             tau = sign(t,pQ)/(abs(pQ) + sqrt(pQ**2 + t/delta))
!!$             write(*,*) 'pQ,t,tau:',real(pQ),real(t),real(tau)
!!$             !          tau = sign(tau,pz)
!!$             p = p + tau*Q(:,1)
!!$             pnorm = norm2(p)
!!$             write(*,*) 'hard final:',real(p),' norm:',real(pnorm)
!!$             info = 1
!!$          else
!!$!             write(*,*) 'indefinite sol: the hard case pnorm > delta'
!!$             write(*,*) 'indefinite sol: NO CONVERGENCE'
!!$             info = 6
!!$          end if
!!$          stop 'indefinite'
!!$!       else
!!$          !          stop 'trstep: The branch in the indefinite case is uncovered.'
!!$!          write(*,*) 'indefinite sol: NO CONVERGENCE'
!!$!          info = 6
!!$       end if
!    end if

  end subroutine trstep_new


  subroutine trstep(B,g,delta,rtol,atol,p,posdef,info)

    real(REAL64), dimension(:,:), intent(in) :: B
    real(REAL64), dimension(:), intent(in) :: g
    real(REAL64), intent(in) :: delta, rtol, atol
    real(REAL64), dimension(:), intent(out) :: p
    logical, intent(out) :: posdef
    integer, intent(out) :: info

    real(REAL64), dimension(2,2) :: Q
    real(REAL64), dimension(2) :: lam,qg,z
    real(REAL64) :: par, pz, pnorm,gnorm,bnorm,t, tau, s2
    logical :: converge!, collapsed
    logical, dimension(size(g)) :: ortho
!    logical, dimensions(2) :: collapses
    integer :: i,k
    logical :: perpendicular

    write(*,*) 'trstep init: delta=',real(delta),' g=',real(g)
    write(*,*) 'B:',real(B)

    info = 666
!    do iter = 1, precision(p)


       !    block
       !      real(REAL64), dimension(2,2) :: BB
       !      BB(1,:) = [2,1]
    !      BB(2,:) = [1,2]
    block
      real(REAL64), dimension(2,2) :: HH
      real(REAL64) :: a,lam1,lam2,c,s
      integer :: i
      HH(1,:) = [1.00063074 ,      2.08462141E-02]
      HH(2,:) = [2.08462141E-02,   1.68902147]
      HH(1,:) = [1.00042152,      -5.55052832E-02]
      HH(2,:) = [-5.55052832E-02,   8.30888748]
      HH(1,:) = [ 0.961768687 ,    -0.161296397]
      HH(2,:) = [-0.161296397 ,     0.740648150]
      HH(1,:) = [0.896569908 ,     0.337577820 ]
      HH(2,:) = [0.337577820 ,     -1.03207862]
      HH(1,:) = [ 1.08179390  ,    0.806767106]
      HH(2,:) = [0.806767106 ,      8.95748234]
      HH(1,:) = [1.11534834   ,   -7.55504519E-02]
      HH(2,:) = [-7.55504519E-02,  0.903266609]

      lam1 = real(2,REAL64)
      lam1 = real(1,REAL64) / real(50,REAL64)
      lam2 = real(1,REAL64)
!      do i = -45, 45, 5
         a = real(i,REAL64)/real(57.29577951308232,REAL64)
         c = cos(a)
         s = sin(a)
         HH(1,:) = [ lam1*c**2+lam2*s**2, (lam1-lam2)*s*c ]
         HH(2,:) = [ (lam1-lam2)*s*c,  lam1*s**2+lam2*c**2]
!         if( i == 18 ) then
!            HH(1,:) = [lam1,real(0,REAL64)]
!            HH(2,:) = [real(0,REAL64),lam2]
!         end if
         call q2eig(B,lam,Q)
!         write(*,*) '== lambda=',i,real(lam),real(norm2(Q(:,1)-[c,-s])), &
!              real(norm2(Q(:,2)-[s,c]))
!         write(*,*) '// Q1=',i,real(Q(:,1))
!         write(*,*) '// Q2=',i,real(Q(:,2))
!         stop
!      end do
!      stop
    end block

    write(*,*) 'jjj',norm2(g)/delta-sum(abs(B)),norm2(g)/delta+sum(abs(B)),norm2(g)/delta

    ! Is B positive definite?
    posdef = all(lam > 0)

!    collapses = abs(dot_product(Q(:,i),g)) < atol
!    collapsed = any(collapsess == .false.)

    if( norm2(g)/delta < 10*epsilon(g) ) then
       ! g == 0 case

       if( all(lam > 0) ) then
          ! in minimum
          p = 0
       else if( lam(1)*lam(2) < 0 ) then
          ! saddle point, looking for min in direction of lam(i) > 0
          do i = 1, 2
             if( lam(i) > 0 ) then
                p = delta*Q(:,i)
             end if
          end do
       else if( all(lam < 0) ) then
          ! maximum
          stop 'trstep in maximum' ! random direction?
       else
          stop 'trstep in doubts while g == 0'
       end if
       write(*,*) 'g==0: p=',real(p),' lam=',real(lam)
       info = 1
       return
    end if

    ! determine whatever any eigen vector and g are orthogonal
    gnorm = norm2(g)
    do i = 1,size(g)
       qg(i) = dot_product(Q(:,i),g)
       ortho(i) = abs(qg(i))/gnorm < sqrt(epsilon(g))
       write(*,*) 'scalar product:',qg(i)
!       if( ortho(i) ) then
!          lam = lam(i)
!       end if
    end do
    perpendicular = any(ortho .eqv. .true.)

    write(error_unit,*) 'trsep initial analysis: positive definite =',posdef, &
         ' perpendicular =',perpendicular

    if( .not. posdef ) then
       ! matrix is not positive definite

       !       if( abs(lam(1) - lam(2)) > 1 ) then
!              Qg(1) = 0
!          call printbpar(lam,Qg,delta)
!          stop
!       end if



       par = -lam(1) + 0.1*abs(lam(1))
       !       par = -1.1*minval(lam)! + 0.1*abs(lam(1))
!       Qg(1) = 0
       call bpar(lam,qg,ortho,delta,par,converge)
       write(*,*) 'par=',par, converge
       if( converge ) then
          call eqsol(B,lam,Q,g,par,p,pnorm)
          write(*,*) 'indefinite sol:',real(p),' norm:',real(pnorm)
          info = 1
          return
       else
          write(*,*) 'indefinite sol:',real(p),' norm:',real(pnorm), 'NON-convergetnt'
          call printbpar(lam,Qg,delta)
          stop
          info = 6
          return
       end if

!!$       if( .false. ) then
!!$          k = 0
!!$          do i = 1, size(ortho)
!!$             if( ortho(i) ) k = i
!!$          end do
!!$
!!$          s2 = 0
!!$          do i = 1, size(qg)
!!$             if( .not. ortho(i) ) then
!!$                s2 = s2 + (qg(i)/(lam(i)+lam(k)))**2
!!$             end if
!!$          end do
!!$          tau = sqrt(delta**2 - s2)
!!$          p = p + tau*Q(:,k)
!!$          pnorm = norm2(p)
!!$          write(*,*) 'the indefinite hard case:',real(tau),real(pnorm),real(s2)
!!$       end if

!       if( .not. collapsed ) then
!       else
!          par = par + lam(1)
!          call eqsol(B,Qvec,g,par,p,pnorm)
!          if( abs(pnorm - delta)/delta < rtol ) then
!             info = 1
!             return
!          else
!            stop 'trstep failed.'
!          end if
!       end if

!!$       if( abs(dot_product(Qvec(:,1),g)) < 100*epsilon(g) .or. &
!!$            abs(dot_product(Qvec(:,2),g)) < 100*epsilon(g) .or. &
!!$            abs(norm2(g)) < 100*epsilon(g) ) then
!!$          write(*,*) 'g:',g
!!$          write(*,*) 'Q:',Qvec
!!$          write(*,*) par
!!$          stop 'Q*g or g is zero.'
!!$       end if

    else
       ! matrix is positive definite


!          par = 0
!       end if
!       write(*,*) 'par=',par, converge

       !    B1 = B
       !    forall( i = 1:2 ) B1(i,i) = B(i,i) + par

!       if( .not. collapsed ) then

!       if( .not. hard ) then

!       par = norm2(g)/delta
       par = 0
       call eqsol(B,lam,Q,g,par,p,pnorm)
!!$    Binv = 0
!!$    forall( i = 1:2 ) Binv(i,i) = 1 / (lam(i) + par)
!!$
!!$    D = matmul(matmul(Qvec,Binv),transpose(Qvec))
!!$!    D = matmul(matmul(transpose(Qvec),Binv),Qvec)
!!$    p = -matmul(D,g)
!!$    pnorm = norm2(p)
       write(*,*) 'sol:',real(p),' norm:',real(pnorm)

!       if(  ) then

          ! inside a hole having parabolic surface
!          if(  pnorm/delta < (1 + rtol) .and. .not. perpendicular ) then
          if(  pnorm/delta < (1 + rtol) ) then
             info = 1
             return
          endif

          ! the current step is outside of TR region,
          ! par is adjusted, step is put on the edge
          !             par = lam(2)

          block
            real(REAL64), dimension(2) :: g1,pp
            real(REAL64) :: a,b,c,d,g2
            g1 = matmul(Q,qg)
            g2 = sum(g1**2)
            write(*,*) real(g1),norm2(g1),sqrt(g2)
            pp = delta*p / pnorm
            a = delta**2
            b = 2*(lam(1)*pp(1)**2 + lam(2)*pp(2)**2)
            c = (lam(1)*pp(1))**2 + (lam(2)*pp(2))**2 - g2
            d = b**2 - 4*a*c
            write(*,*) a,b,c,d
            if( d > 0 ) then
               d = sqrt(d)
               write(*,*) 'q2sol:',(-b-d)/2/a,(-b+d)/2/a
               par = (-b+d)/2/a
            end if
!            stop
          end block

!          par = 0.1*lam(1)
!          par = 0.1*minval(lam)
!          write(*,*) qg(1)/lam(1),qg(2)/lam(2)
          if( abs(qg(1)/lam(1)) < abs(qg(2)/lam(2)) ) then
!             par = abs(qg(1))*delta - lam(1)
          else
!             par = abs(qg(2))*delta - lam(2)
          end if
          bnorm = sum(abs(lam))
!          par = sqrt((gnorm/delta-bnorm)*(gnorm/delta+bnorm))
          !          par = gnorm/delta + sum(abs(lam))
          par = norm2(g)/delta - norm2(lam)
          write(*,*) 'lll',norm2(g)/delta-sum(abs(lam)), &
               norm2(g)/delta+sum(abs(lam)),real(par)
          call bpar(lam,qg,ortho,delta,par,converge)
!          converge = .true.
!          par = (gnorm / delta - sum(abs(lam)))/2
          if( converge ) then
             call eqsol(B,lam,Q,g,par,p,pnorm)
             write(*,*) 'outside sol:',real(p),' norm:',real(pnorm)
          else
             info = 6
             return
          end if
!          stop

          if( abs(pnorm - delta)/delta < rtol ) then
             info = 1
             return
             !       else
             !          info = 2
!          return
          !                stop 'trstep failed.'
          end if
!       end if

!       else
!       ! on the edge
!       if( par > 0  .and. abs(pnorm - delta)/delta < rtol ) then
!          info = 1
!          return
!       end if

          if( pnorm < delta .and. perpendicular ) then
             ! the hard case
             k = 0
             do i = 1, size(ortho)
                if( ortho(i) ) k = i
             end do

             par = lam(k)
             call eqsol(B,lam,Q,g,par,p,pnorm)
             write(*,*) 'p=',real(p),' pnorm=',real(pnorm)

             ! solution only in sub-space lam(i) /= lam(k) ?
             do i = 1, size(qg)
                if( .not. ortho(i) ) then
                   p = qg(i)*Q(:,i) / (lam(i) + lam(k))
                end if
             end do
             write(*,*) 'p=',real(p),' pnorm=',real(norm2(p))

             p = 0
             pnorm = 0
             write(*,*) 'hard sol:',real(p),' norm:',real(pnorm)

             pz = dot_product(p,Q(:,k)) / delta
             t = (delta - pnorm)*((delta + pnorm)/delta)
             tau = t/(abs(pz) + sqrt(pz**2 + t/delta))
             write(*,*) 'pz,t:',real(pz),real(t)
             tau = sign(tau,pz)
!!$       p = p + tau*z

             s2 = 0
             do i = 1, size(qg)
                if( .not. ortho(i) ) then
                   s2 = s2 + (qg(i)/(lam(i)+lam(k)))**2
                end if
             end do
             write(*,*) 's2=',real(s2),' qg=',real(qg)
             !          if( delta**2 - s2 < 0 ) write(*,*) 'sqrt(0)?',delta,s2,delta**2
             !          tau = sqrt(delta**2 - s2)
             p = p + tau*Q(:,k)
             pnorm = norm2(p)
             write(*,*) 'the hard case:',real(tau),real(pnorm),real(s2)
             info = 1
             return
          end if
    !   end if

!          info = 1
!          return
    end if

!!$    parallel = abs(dot_product(Qvec(:,1),g)) < 1e-10 .or. &
!!$         abs(dot_product(Qvec(:,2),g)) < 1e-10
!!$
!!$    if( pnorm > delta .and. .not. parallel ) then
!!$       !       if( abs(par) < epsilon(par) ) par = lam(2)
!!$       par = par + lam(1)
!!$       call bpar(lam(1),lam(2),Qvec(:,1),Qvec(:,2),g,delta,par,converge)
!!$       Binv = 0
!!$       forall( i = 1:2 ) Binv(i,i) = 1 / (lam(i) + par)
!!$       D = matmul(matmul(Qvec,Binv),transpose(Qvec))
!!$       p = -matmul(D,g)
!!$       pnorm = norm2(p)
!!$       write(*,*) 'updated par:',real(par),real(pnorm/delta)
!!$       info = 1
!!$       return
!!$    end if
!!$
!!$    if( pnorm < delta .and. parallel ) then
!!$       ! the hard case
!!$       z = Qvec(:,2)
!!$
!!$       pz = dot_product(p,z) / delta
!!$       t = (delta - pnorm)*((delta + pnorm)/delta)
!!$       tau = t/(abs(pz) + sqrt(pz**2 + t/delta))
!!$       tau = sign(tau,pz)
!!$       p = p + tau*z
!!$
!!$       write(*,*) 'the hard case:',real(tau),real(z)
!!$       stop 'the hard case'
!!$    end if
!!$


    ! check
!    write(*,*) 'check:',real(matmul(D,B))

    ! check
!    write(*,*) 'check:',matmul(B,p)+g

  end subroutine trstep

!!$  subroutine trstep1(B,g,delta,rtol,atol,p,info)
!!$
!!$    use cholesky
!!$
!!$    real(REAL64), dimension(:,:), intent(in) :: B
!!$    real(REAL64), dimension(:), intent(in) :: g
!!$    real(REAL64), intent(in) :: delta, rtol, atol
!!$    real(REAL64), dimension(:), intent(out) :: p
!!$    integer, intent(out) :: info
!!$
!!$    real(REAL64), dimension(2,2) :: B1, R, Qvec
!!$    real(REAL64), dimension(2) :: lam, z, q
!!$    real(REAL64) :: gnorm, bnorm, pnorm, qnorm, rznorm, znorm, rpnorm, &
!!$         par, parc, pars, parl, paru, tau, pz, t
!!$    integer :: i,iter,indef
!!$    character(len=80), parameter :: fmt = "(5x,a,i3,l2,1x,f10.2,2(a,2f8.2,a),f10.2)"
!!$    logical :: converge
!!$
!!$
!!$    par = real(0,REAL64)
!!$    info = 0
!!$    do iter = 1, precision(p)
!!$
!!$       if( debug ) write(*,fmt) "trstep:",iter,indef==0,par,", ",p,", ","",z,", ",pnorm
!!$
!!$       ! update the matrix
!!$       B1 = B
!!$       forall( i = 1:2 ) B1(i,i) = B(i,i) + par
!!$
!!$       if( debug ) then
!!$          do i = 1, 2
!!$             write(*,'(5x,a,i0,a,100es18.6)') 'b1(',i,',:):',b1(i,:)
!!$          end do
!!$       end if
!!$
!!$       ! Cholesky decomposition, and the step estimation, if possible
!!$       call chol_rrfac(B1,R,indef)
!!$       !       write(*,*) 'indef=',indef, lambda
!!$       if( indef == 0 ) then
!!$          call chol_rrsolve(R,-g,p)
!!$          call chol_rrtsolve(R,p,q)
!!$          qnorm = norm2(q)
!!$          pnorm = norm2(p)
!!$
!!$          ! inside a hole having parabolic surface
!!$          if( abs(par) < epsilon(par) .and. pnorm/delta < (1 + rtol) ) then
!!$             info = 1
!!$             exit
!!$          end if
!!$
!!$          ! on the edge
!!$          if( par > 0  .and. abs(pnorm - delta)/delta < rtol ) then
!!$             info = 1
!!$             exit
!!$          end if
!!$
!!$          ! Compute the Newton correction parc to par.
!!$          parc = (pnorm/qnorm)**2*((pnorm - delta)/delta)
!!$
!!$          ! prozkoumat pripady: g == 0, g || Q(:,i)
!!$
!!$       else ! the indefinite case
!!$
!!$          call q2eig(B,lam,Qvec)
!!$          call bpar(lam(1),lam(2),Qvec(:,1),Qvec(:,2),g,delta,par,converge)
!!$          write(*,*) 'par=',par, converge
!!$
!!$          block
!!$            real(REAL64), dimension(2,2) :: E
!!$            integer :: i
!!$            E = 0
!!$            forall(i=1:2) E(i,i) = 1
!!$            do i = 1, 2
!!$               write(*,*) 'ort.test:',real(matmul(B-lam(i)*E,Qvec(:,i)))
!!$            end do
!!$          end block
!!$
!!$       end if
!!$
!!$       par = par + parc
!!$
!!$    end do
!!$
!!$  end subroutine trstep1



  subroutine bpar_new(lam,posdef,Q,Qg,delta,rtol,p,pnorm,iflag)

    real(REAL64), dimension(:,:), intent(in) :: Q
    real(REAL64), dimension(:), intent(in) :: lam,Qg
    logical, intent(in) :: posdef
    real(REAL64), intent(in) :: delta, rtol
    !    real(REAL64), intent(out) :: par
    real(REAL64), dimension(:), intent(out) :: p
    real(REAL64), intent(out) :: pnorm
    integer, intent(out) :: iflag

    real(REAL64), dimension(size(Qg)) :: y,w,u
    real(REAL64) :: d,h,p2,q2,par
    real(REAL64), parameter :: eps = epsilon(par)
    integer :: iter

!    tol = 0.1*rtol!*delta
!    tol = rtol

!    write(*,*) lam
!    write(*,*) Q

    if( posdef ) then
       ! |p| > delta, par should be adjusted
!       write(*,*) norm2(Qg),delta
!       par = max(norm2(Qg) / delta - norm2(lam),real(0.001,REAL64),epsilon(par))
       par = max(norm2(Qg) / delta - norm2(lam),eps)
!       write(*,*) norm2(Qg) / delta - norm2(lam)
    else
       ! indefinite
!       par = max(norm2(Qg) / delta - norm2(lam),-1.001*minval(lam),0.001_REAL64)
       par = max(norm2(Qg) / delta - norm2(lam),-(1+eps)*minval(lam),eps)
!       if( abs(par) < epsilon(par) ) par = epsilon(par)
!       write(*,*) norm2(Qg) / delta - norm2(lam),-1.001*minval(lam)
    end if

    write(*,*) 'bpar:',0,par
    do iter = 1, precision(par)

       if( par < 0 ) exit
 !         iflag = -1
 !         return
 !      end if
       y = lam + par

       !if( any(abs(y) < 1e-12) )
       if( debug .and. any(abs(y) < 10*epsilon(y)) ) &
            write(*,*) 'y==0?',y,lam,par,norm2(Qg),delta

!       write(*,*) Qg

       u = Qg / y
       w = u**2
       p2 = sum(w)
       q2 = sum(w/y)
       ! p2 > 0 is allways, bpar is called only for non-hard case
!       if( p2 > real(0,REAL64) ) then
          pnorm = sqrt(p2)
          h = (pnorm - delta) / delta
          d = (p2 / q2)*h
!       else
!          h =-real(1,REAL64)
!          d = real(0,REAL64)
!       end if
       par = par + d
       write(*,*) 'bpar:',iter,real(par),real(d),real(h)
!       if( abs(d) < tol*(1 + par) ) then
       if( abs(h) < rtol ) then
          p = -matmul(Q,u)
!          pnorm = p
          iflag = 0
!          write(*,*) p,sqrt(q2)
          return
       end if
    end do
    iflag = 666

!!$    block
!!$      real(REAL64) :: f, e2, s
!!$      integer :: i
!!$
!!$       open(9,file='/tmp/bpar')
!!$       do iter = -100,100
!!$          par = 0.01*iter
!!$          y = lam + par
!!$          s = 0
!!$          do i = 1,size(qg)
!!$             if( abs(y(i)) > epsilon(y) ) then
!!$                e2 = (qg(i)/y(i))**2
!!$                s = s + e2
!!$             end if
!!$          end do
!!$          f = 1 / delta - 1 / sqrt(s)
!!$          write(9,*) par,f,delta-sqrt(s)
!!$       end do
!!$       close(9)
!!$!       stop
!!$     end block

  end subroutine bpar_new


  subroutine bpar(lam,qg,ortho,dlt,par,converge)

    real(REAL64), dimension(:), intent(in) :: lam,qg
    logical, dimension(:), intent(in) :: ortho
!    real(REAL64), dimension(:,:), intent(in) :: Q
    real(REAL64), intent(in) :: dlt
    real(REAL64), intent(in out) :: par
    logical, intent(out) :: converge
!    real(REAL64), parameter :: tol = sqrt(epsilon(lam))
    real(REAL64), parameter :: zero = real(0,REAL64)
    real(REAL64), parameter :: one = real(1,REAL64)
    real(REAL64), dimension(size(qg)) :: y
    real(REAL64) :: e2,s,t,f,df,d,p
    integer :: iter,i




    converge = .false.

!    if( norm2(g) < tol ) return

!    a = dot_product(q1,g)
!    a2 = dot_product(q1,g)**2
!    b2 = dot_product(q2,g)**2

!    xlam = 0
!    do i = 1,size(qg)
!       qg(i) = dot_product(Q(:,i),g)
!       ortho(i) = abs(qg(i)) < tol
!       if( ortho(i) ) then
!          xlam = lam(i)
!       end if
!    end do
!    xlam = 0

!    write(*,*) par,xpar,ortho

!    write(*,*) a,sqrt(b2)

!    mu = (lam1 + lam2) / 2
!    mu = -lam1+1e-3
!    mu = 0.9*lam2
!    mu = sqrt(a)/dlt - lam1
!    mu = -lam1 + abs(a)
!    par = par + 0.1*maxval(abs(lam))  !maxval(abs(qg))

    do iter = 1, precision(par)

       y = lam + par
!       u = lam1 + mu
!       u2 = u**2
!       v = lam2 + mu
!       v2 = v**2

!       if( abs(u) < epsilon(u) .or. abs(v) < epsilon(v) ) then
       !          write(*,*) mu, lam1, lam2, a
       if( any(abs(y) < epsilon(y)) ) then
          write(*,*) real(par),real(lam)
          stop 'lambda undefined?'
       end if
!       if( .not. (a2/u2 + b2/v2 > 0) ) then
!          write(*,*) a2,u,b2,v
!          stop 'zeros?'
!       end if

       s = zero
       t = zero
       do i = 1,size(ortho)
!          if( .not. ortho(i) ) then
             e2 = (qg(i)/y(i))**2
             s = s + e2
             t = t + e2/y(i)
!          end if
       end do

!       write(*,*) s,dlt
!       if( s < dlt**2 ) then
!          write(*,*) 'bpar: inside delta.'
!          return
!       end if

       d = sqrt(s)

       f = one / dlt - one / d
       df = t / d**3

!       s = 0
!       do i = 1,2
!          if( .not. ortho(i) ) then
!             s = s + (qg(i)/y(i))**2/y(i)
!          end if
!       end do
!       df = s / d**3

!       f = 1 / dlt - 1 / d
!       df = ((a2/u2)/u + (b2/v2)/v)/d**3
       p = f/df
!       p = (d**2/t)*(1 - d/dlt)
!       p = (s/dlt - 1)*(d/t)
!       write(*,*) p,f/df
       par = par + p

       write(*,*) 'bpar:',iter,real(par),real(p),real(1 - d/dlt)
       !       if( abs(p) < 1e-5 .or. (abs(p) < 1e-5*abs(par)) ) then
!       if( abs(p/(1 + par)) < 1e-3 ) then
       if( abs(p)/(1 + par) < 1e-3 ) then
          converge = .true.
          exit
!          return
       end if


    end do

    if( par < 0 ) then

!!$       open(9,file='/tmp/bpar')
!!$       do iter = -100,100
!!$          par = 0.1*iter
!!$          y = lam + par
!!$          s = 0
!!$          do i = 1,size(ortho)
!!$             if( abs(y(i)) > epsilon(y) ) then
!!$                e2 = (qg(i)/y(i))**2
!!$                s = s + e2
!!$             end if
!!$          end do
!!$          f = 1 / dlt - 1 / sqrt(s)
!!$          write(9,*) par,f,dlt-sqrt(s)
!!$       end do
!!$       close(9)
!!$       stop

       write(*,*) 'bpar: parameter negative:',par,lam(1)
       converge = .false.
       return
    end if


!    stop
  end subroutine bpar

  subroutine printbpar(lam,Qg,dlt)

    real(REAL64), dimension(:), intent(in) :: lam, Qg
    real(REAL64), intent(in) :: dlt
    real(REAL64) :: par,d,d2,f
    integer :: iter

    open(9,file='/tmp/bpar')
    do iter = -200,200
       par = 0.1*iter
!       d2 = sum(Qg / (lam + par)**2, abs(lam+par) > epsilon(d))
!       if( d2 > 0 ) then
!          d = sqrt(d2)
!       else
!          d = 0
!          f = 0
       !       end if
       if( abs(lam(1)+par) > 0 .and. abs(lam(2)+par) > 0 ) then
          d = hypot(Qg(1) / (lam(1) + par),Qg(2) / (lam(2) + par))
       else if( abs(lam(1)+par) > 0 .and. abs(lam(2)+par) < epsilon(lam) ) then
          d = abs(Qg(1) / (lam(1) + par))
       else if( abs(lam(2)+par) > 0 .and. abs(lam(1)+par) < epsilon(lam) ) then
          d = abs(Qg(2) / (lam(2) + par))
       else
          d = 0
       end if
       f = 1 / dlt - 1 / d
       write(9,*) par,f,dlt - d
    end do
    close(9)

  end subroutine printbpar

  subroutine eqsol(B,lam,Q,g,par,p,pnorm)

    real(REAL64), dimension(:,:), intent(in) :: B,Q
    real(REAL64), dimension(:), intent(in) :: lam, g
    real(REAL64), dimension(:), intent(out) :: p
    real(REAL64), intent(in) :: par
    real(REAL64), intent(out) :: pnorm

    real(REAL64), dimension(size(B,1),size(B,2)) :: Binv, D
    integer :: i

    Binv = 0
    forall( i = 1:size(B,1) ) Binv(i,i) = 1 / (lam(i) + par)
    D = matmul(matmul(Q,Binv),transpose(Q))
    p = -matmul(D,g)
    pnorm = norm2(p)
!    write(*,*) 'sol:',real(p),' norm:',real(pnorm)

    ! nasobeni odzadu:
    ! p = matmul(transpose(Q),-g)
    ! p = p / (lam + par)
    ! p = matmul(Q,p)
    p = matmul(Q,matmul(transpose(Q),-g)/(lam + par))
    pnorm = norm2(p)
    if( .not. all(abs(p+matmul(D,g)) < 10*epsilon(p) ) ) stop 'eqsol'
!    write(*,*) 'sol:',real(p),' norm:',real(pnorm),real(p+matmul(D,g))

  end subroutine eqsol

  pure subroutine eqsol_new(B,lam,Q,Qg,par,p,pnorm)

    real(REAL64), dimension(:,:), intent(in) :: B, Q
    real(REAL64), dimension(:), intent(in) :: lam, Qg
    real(REAL64), intent(in) :: par
    real(REAL64), dimension(:), intent(out) :: p
    real(REAL64), intent(out) :: pnorm

    p = -matmul(Q,Qg/(lam + par))
    pnorm = norm2(p)
!    write(*,*) 'sol:',real(p),' norm:',real(pnorm),real(p+matmul(D,g))

  end subroutine eqsol_new


  function sminupdate(fcn,t,s) result(smin)

    use fmm

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), intent(in) :: t,s
    real(REAL64) :: smin,rmin,rmax,tol
    logical :: found

    smin_fcn => fcn
    smin_t = t
    smin_count = 0

    tol = 0.01*s
    rmin = sqrt(epsilon(rmin))
    rmax = s
    smin = zeroin(rmin,rmax,scalemin,tol)

    ! Does the minimum falls inside the interval smin < .. < smax ?
    found = rmin + tol < smin .and. smin < rmax - tol

    !    if( debug )
    write(*,'(a,g20.6,a,g0.2,l2,3x,a,i2)') &
         'Minimal scale limit:',smin," +- ",tol,found,'count:',smin_count

!    if( .not. found ) stop 'smin not found'
!    data => null()

  end function sminupdate

  real(REAL64) function scalemin(s)
!    use rfun
    real(REAL64), intent(in) :: s
!    scalemin = sum(dtukey((data-chi2t)/s), dtukey((data-chi2t)/s) > 0) / size(data) - 0.4
    !    scalemin = sum(dpshuber((data-chi2t)/s), dpshuber((data-chi2t)/s) > 0) / size(data) - 0.4
    scalemin = smin_fcn%nfrac(smin_t,s) - 0.25
    smin_count = smin_count + 1
  end function scalemin


  subroutine smnewton(fcn,qvec,t,s,tol,delta,converge)

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), dimension(:), intent(in) :: qvec
    real(REAL64), intent(in) :: tol, delta
    real(REAL64), intent(in out) :: t,s
    logical, intent(out) :: converge

    real(REAL64), dimension(2) :: fvec, x, p, xp
    real(REAL64), dimension(2,2) :: hess
    real(REAL64) :: fval,fv,g,pHp,d,rho,ared,pred,dlt
    integer :: iter

    converge = .false.
    dlt = delta
    rho = 0

!    write(*,*) 'snewton init:',t,s,qvec,tol

!    block
!      real(REAL64), dimension(2) :: fvec
!      real(REAL64), dimension(2,2) :: hess
!      integer :: i
!      real(REAL64) :: ss
!      open(17,file='/tmp/z')
!      do i = 1, 500
!         ss = 0.01*i
!         call fcn%grad(t,ss,fvec)
!         call fcn%hess(t,ss,hess)
!         g = fvec(2)
!         pHp = hess(2,2)
!         write(17,*) ss,fcn%val(t,ss),g,pHp
!      end do
!      close(17)
!    end block

    fval = fcn%val(t,s)
    x = [ t, s ]
    do iter = 1, precision(x)

!       write(*,*) x
       call fcn%grad(x(1),x(2),fvec)
       call fcn%hess(x(1),x(2),hess)
!       write(*,*) 'grad:',fvec
!       write(*,*) 'hess:',hess

       g = dot_product(fvec,qvec)
       pHp = dot_product(qvec,matmul(hess,qvec))

       if( pHp < 0 ) then
          d = -g / pHp
       else
          d = dlt
!          what = .true.
       end if
       if( abs(d) >= dlt ) d = sign(dlt,g)

!!$          block
!!$            integer :: i
!!$            real(REAL64)  :: d,u,v
!!$            d = 0.1*dlt
!!$            do i = -10, 10
!!$               u = x(1)+i*d
!!$               v = x(2)+i*d
!!$!               if( v > 0 ) write(*,*) u,v,fcn%val(u,v)
!!$            end do
!!$          end block


       p = d*qvec
       xp = x + p

!       if( debug ) &
            !            write(*,'(a,i2,1pg13.5,1p2g10.2,1pg11.3,1pe12.3,3l1)') &
!       write(*,*) "smnewton,inc.,f,f': ",iter,x,d,g,pHp,dlt

       !       write(*,*) 'xp:',xp,xp(2) > 0
       rho = 0
       if( xp(2) > 0 ) then
          fv = fcn%val(xp(1),xp(2))
!       if( debug ) write(17,*) s,fval

          ared = fval - fv
          !       pred = -d*(fs + 0.5*h22*d)
          pred = -dot_product(fvec + 0.5*matmul(hess,p),p)
!          rho = ared / pred
          if( abs(pred) > 0 ) rho = ared / pred
!          rho = ared / pred
!       else
!          rho = 0
       end if
!       write(*,*) 'ared,pred,rho:',ared,pred,rho

       if( rho > 0.5e-3 ) then
          x = xp
          fval = fv

          converge = abs(d) < tol
          if( converge ) then
             t = x(1)
             s = x(2)
             exit
          end if

       end if

       if( rho > 0.75 ) then
          if( abs(d) > 0.9*dlt ) then
!             dlt = min(2*dlt,dltmax)
             dlt = 2*dlt
          end if
       else if( rho < 0.25 ) then
!          dlt = max(dlt/2,dltmin)
          dlt = dlt/2
       end if
!       dlt = min(dlt,s-smin)


    end do

!    if( what ) stop 'what'
!    stop

  end subroutine smnewton

  subroutine ridge_line(fcn,qvec,t,s,rtol,quasi,hess,delta,converge)

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), dimension(:), intent(in) :: qvec
    real(REAL64), intent(in) :: rtol
    real(REAL64), intent(in out) :: t,s,delta
    logical, intent(in) :: quasi
    real(REAL64), dimension(2,2), intent(in) :: hess
    logical, intent(out) :: converge

!    real(REAL64), parameter :: qmin = real(1,REAL64) / sqrt(real(2,REAL64))
    real(REAL64), dimension(2) :: fvec, fgrad, x, p, xp, y!, lam, qvec
    real(REAL64), dimension(2,2) :: B!,Q
    real(REAL64) :: g,pHp,d,ared, pred, dlt, rho, fv, fval, u, gtol, fmax, dmax
    integer :: i,iter
    logical :: reli

    converge = .false.

    ! determine delta
    delta = 0.5*s
!    call fcn%grad(t,s,fvec)
!    call fcn%hess(t,s,hess)
!    g = dot_product(fvec,qvec)
!    pHp = dot_product(qvec,matmul(hess,qvec))
!    if( pHp < 0 ) then
!!       delta = 0.5*abs(g/pHp)
!    end if
!    write(*,*) pHp,delta,g

!    smin = max(s - 5*delta*qvec(2),0.0_REAL64)
!    smax = s + 5*delta*qvec(2)
    write(*,*) 'ridge_line init:',t,s,qvec,rtol,delta

    fval = fcn%val(t,s)
    call fcn%grad(t,s,fvec)
    B = hess
!    if( quasi ) then
!       B = real(0,REAL64)
!       forall( i = 1:2 ) B(i,i) = real(1,REAL64)
!    else
!       B = hess
!       call fcn%hess(t,s,B)
!       if( fcn%break ) then
!          converge = .false.
!          return
!       end if
!    end if

    dlt = delta
    x = [ t, s ]
    d = 2*rtol
    pHp = 1
    do iter = 1, 3*precision(x)

!       if( x(2) < 0 ) exit

!       write(*,*) 'podminka:',d,rtol*x(2),pHp
       converge = abs(d) < rtol*x(2) .and. pHp < 0
       if( converge ) then
          t = x(1)
          s = x(2)
          exit
       end if

!       write(*,*) x
!       call fcn%grad(x(1),x(2),fvec)
!       call fcn%hess(x(1),x(2),hess)
!       write(*,*) 'grad:',fvec
!       write(*,*) 'hess:',hess

!       call q2eig(B,lam,Q)
!       do i = 1, 2
!          if( abs(Q(2,i)) > qmin ) qvec = Q(:,i)
!       end do

       g = dot_product(fvec,qvec)
       pHp = dot_product(qvec,matmul(B,qvec))
       write(*,*) 'g,pHp:',g,pHp
!       if( .not. (pHp < 0) ) exit

       if( pHp < 0 ) then
          d = -g / pHp
          if( abs(d) > dlt ) then
             d = sign(dlt,g)


!!$             gtol = 0.1
!!$            if( g*qvec(2) < 0 .and. abs(qvec(2)) > epsilon(qvec) ) then
!!$               dmax = min(0.9*abs(x(2)/qvec(2)),dlt)
!!$            else
!!$               dmax = dlt
!!$            end if
!!$            dmax = dlt
!!$            write(*,*) g,qvec(2),g*qvec(2),x(2),dmax
!!$            call goldfinger(fcn,x,sign(1.0_REAL64,g)*qvec,dmax,gtol,u,fmax,reli)
!!$            d = u*sign(dlt,g)
!!$            p = d*qvec

          end if
          p = d*qvec
       else
          ! call golden() ?
          block
            real(REAL64), parameter :: z = 1
            real(REAL64) :: dd,umax,f
            real(REAL64), dimension(2) :: xp,pp,pmax

            dd = sign(dlt,g)
            dd = min(0.9*abs(x(2)/qvec(2)),dlt)
            dd = sign(min(0.9*x(2),dlt),g)
            pp = dd*qvec
            xp = x + pp
            umax = 0
            pmax = pp
!            fmax = 0
!            write(*,*) pp,dd
            fmax = fcn%val(xp(1),xp(2))
            do i = 0,10
               u = i / 10.0
               xp = x + u*pp
!               write(*,'(4f10.5)') xp,fcn%val(xp(1),xp(2)), fmax
               f = fcn%val(xp(1),xp(2))
               if( f > fmax ) then
!                  write(*,*) fcn%val(xp(1),xp(2)),fmax,u,xp
                  fmax = f
                  umax = u
                  pmax = u*pp
                  !                  fmax = fcn%val(xp(1),xp(2))
                  !                  stop 'golden ?'
               end if
            end do
!            write(*,*) 'grid:',umax,fmax
            p = pmax
            !            stop

            gtol = 0.05
            dmax = min(0.9*x(2),dlt)
            d = sign(dmax,g)
            write(*,*) 'line search: x=',real(x),' to ',real(x+d*qvec)
            write(*,*) '         ...',real(g),real(qvec),real(x(2)),real(dmax)
            call goldfinger(fcn,z,x,qvec,d,gtol,u,fmax,reli)
            p = u*d*qvec
!            write(*,*) u,x+p
          end block
!          d = sign(dlt,g)
       end if

       xp = x + p


!       if( .not. (smin < x(2) .and. x(2) < smax) ) exit
!       write(*,*) 'abs:',abs(x(2) - s),delta
!       if( abs(x(2) - s) > delta ) exit


       if( xp(2) > 0 ) then

          fv = fcn%val(xp(1),xp(2))
          call fcn%grad(xp(1),xp(2),fvec)
!          call fcn%hess(xp(1),xp(2),hess)

          y = fvec - fgrad
          ared = fval -  fv
          pred = -(dot_product(fvec,p) + 0.5*dot_product(p,matmul(B,p)))
!          write(*,*) 'rho:',fv,fval,ared,pred, ared/pred, dot_product(fvec,p),dot_product(p,matmul(B,p))

          if( abs(ared) < epsilon(ared) .and. abs(pred) < epsilon(pred) ) then
             rho = 1
          else if( abs(pred) < epsilon(pred) .and. abs(ared) > 0 ) then
             rho = 0
          else
             rho = ared / pred
          end if

       else
          rho = 0
       end if

!       if( debug ) &
            !            write(*,'(a,i2,1pg13.5,1p2g10.2,1pg11.3,1pe12.3,3l1)') &
       write(*,'(a,i2,a,2f8.3,":",2f7.3,2es10.2," |",f6.2,f6.2)') &
            "ridge_line ",iter,": ",xp,fv,d,g,pHp,dlt,rho

       if( debug ) write(8,*) xp


       if( rho > eta ) then
          x = xp
          fval = fv
          fgrad = fvec
       else
          x = x
       end if

       if( rho > 0.75 ) then
          if( abs(d) > 0.9*dlt ) then
             dlt = 2*dlt
          end if
       else if( 0.25 <= rho .and. rho <= 0.75 ) then
          dlt = dlt
       else ! if ( rho < 0.25 ) then
          dlt = dlt / 2
       end if
       dlt = min(dlt,s/2)
       dlt = max(rtol*x(2),dlt)

       if( quasi ) then
          call sr1update(p,y,norm2(p),B)
       else
          call fcn%hess(x(1),x(2),B)
          if( fcn%break ) then
             converge = .false.
             return
          end if
       end if


!       x = xp

    end do

    ! the osculation circle radius
!    delta = 1/abs(pHp)
!    fv = sum(fvec**2)
!    d = fvec(1)*hess(2,2) - fvec(2)*hess(1,1)
!    if( fv > 0 .and. d > 0 ) then
!       delta = abs(d)/fv**1.5
!    end if
!!    delta =
!    delta = 0.5
    delta = 0.25*s
!    write(*,*) 'delta determination:',hess,delta,fvec


    block
      real(REAL64) :: u
      open(9,file='/tmp/rl')
      do i = -200, 200
         u = 0.1*i
         xp = [t,s] + u*qvec
         if( xp(2) > 0 ) then
            write(9,*) xp,fcn%val(xp(1),xp(2))
         end if
      end do
      close(9)
    end block

  end subroutine ridge_line


  subroutine ridge(funmean,t0,s0,nopt,topt,sopt)

    class(tsr1mean_function), pointer, intent(in) :: funmean
    real(REAL64), intent(in) :: t0,s0
    integer, intent(out) :: nopt
    real(REAL64), dimension(:), intent(out) :: topt,sopt

    real(REAL64), parameter :: qmin = real(1,REAL64) / sqrt(real(2,REAL64))
    real(REAL64) :: h22, rtol, d, t, s, tp, sp, delta, fmax, u
    real(REAL64), dimension(2,2) :: hess,Q
    real(REAL64), dimension(2) :: lam, fvec
    integer :: ieig,i,k,l
    logical :: reli

    if( size(topt) < 2 .or. size(sopt) < 2 ) stop 'ridge(): dimensions < 2.'

    t = t0
    s = s0

    write(*,*) 'ridge start:',t,s
    nopt = 1
    topt(1) = t
    sopt(1) = s
    return

    call funmean%hess(t,s,hess)
    h22 = hess(1,1)*hess(2,2) - hess(1,2)**2
    write(*,*) hess,h22
!    if( hess(1,1)*h22 > 0 ) return ! we are in safe minimum throat
    !    if( hess(1,1) <= 0 .and. h22 > 0 ) return ! the direction of s is interesting
    if( hess(2,2) > 0 ) then
       write(*,*) 'This ridge seems to be positive definite along s-direction .. all done.'
!       return
    end if

    call q2eig(hess,lam,Q)
    write(*,*) 'lam:',lam
    ieig = 0
    do i = 1, 2
       write(*,*) real(Q(:,i))
       if( abs(dot_product(Q(:,i),real([0,1],REAL64))) > qmin ) ieig = i
!       if( abs(Q(2,i)) > qmin ) ieig = i
    end do
    write(*,*) ieig, dot_product(Q(:,ieig),real([0,1],REAL64)),real(lam(ieig))

    ! the minimalisation along t direction
!    goto 17
    block

      real(REAL64), parameter :: z = -1
      real(REAL64) :: smin, smax, ds, d, gtol, u, fmin, tmin, tmax
      logical :: reli
      integer :: n,i,ie

      n = size(funmean%data)
      ds = s/sqrt(real(n))
      ie = 2 - ieig + 1
      tmin = t - Q(1,ie)*3*ds
      tmax = t + Q(1,ie)*3*ds
      smin = s - Q(2,ie)*3*ds
      smax = s + Q(2,ie)*3*ds
      d = (tmax - tmin) / Q(1,ie)
      gtol = 0.01
      write(*,*) 'Ridge find in s=',real(smin),'..',real(smax),' d=',real(d),&
           ' t=',real(tmin),'..',real(tmax)
      call goldfinger(funmean,z,[tmin,smin],Q(:,ie),d,gtol,u,fmin,reli)
      s = smin + u*d*Q(2,ie)
      t = tmin + u*d*Q(1,ie)
      write(*,*) 'Ridge min for t,s=',real(t),real(s),u,1-u
      if( abs(u) <= 1.5*gtol .or. abs(1-u) <= 1.5*gtol ) then
         write(*,*) 'This ridge seems to be not well formed, skip.'
         t = t0
         s = s0
      end if

    end block
!17 continue

    ! Is the ridge inside expected solution?
    block

      use noise

      real(REAL64), parameter :: z = 1
!      real(REAL64), dimension(2), parameter :: dir = [0,1]
      real(REAL64), dimension(2) :: dir = [0,1]
      real(REAL64) :: smin, smax, ds, dt, gtol, u, fmax, tmin, tmax
      logical :: reli
      integer :: n,i
      n = size(funmean%data)
      ds = s/sqrt(real(n))

      smin = huge(smin)
      do i = 1, n
         d = abs(funmean%data(i) - t)
         if( d < smin .and. d > 10*epsilon(s) ) then
            smin = d
         end if
      end do
      smin = max(s - 0.5*ds,smin)
      smax = s + 3*ds

      smin = 1.41*s*inverf(0.2_REAL64)
      smax = 1.41*s*inverf(0.8_REAL64)
      write(*,*) 'Half-normal quantiles:',real(smin),real(smax)


      call funmean%grad(t,s,dir)
      dir = dir / norm2(dir)
      write(*,*) 'grad=',real(dir)
!      if( d > 0 ) dir = dir / d
!      if( dir(2) < 0 ) dir = matmul(reshape([-1,0,0,-1],[2,2]),dir)
      !      if( abs(dir(2)) < 1e-3 ) stop 'dir in s is zero!'
      dir = Q(:,ieig)
      if( dir(2) < 0 ) dir = matmul(reshape([-1,0,0,-1],[2,2]),dir)
!      d = sign((smax - smin) / abs(dir(2)),dir(2))
      d = (smax - smin) / dir(2)
      dt = d / abs(dir(1))
      tmin = t - dir(1)*(s-smin)
      tmax = t + dir(1)*(smax-s)
      gtol = 0.01
      write(*,*) 'Ridge find in s=',real(smin),'..',real(smax),' d=',real(d),&
           ' t=',real(tmin),'..',real(tmax)
      call goldfinger(funmean,z,[tmin,smin],dir,d,gtol,u,fmax,reli)
      s = smin + u*d*dir(2)
      t = tmin + u*(tmax-tmin)
      write(*,*) 'Ridge max for t,s=',real(t),real(s),u,1-u
      if( abs(u) <= 1.5*gtol .or. abs(1-u) <= 1.5*gtol ) then
         write(*,*) 'This ridge seems to be positive definite along s-direction .. all done.'
         return
      end if
    end block


    ! estimated solution is in an indefinite region
    call funmean%hess(t,s,hess)
    call q2eig(hess,lam,Q)
    write(*,*) 'lam:',lam

    ieig = 0
    do i = 1, 2
       write(*,*) real(Q(:,i))
       if( abs(dot_product(Q(:,i),real([0,1],REAL64))) > qmin ) ieig = i
!       if( abs(Q(2,i)) > qmin ) ieig = i
    end do
    write(*,*) ieig, dot_product(Q(:,ieig),real([0,1],REAL64)),real(lam(ieig))

!    ieig = 0
!    do i = 1, 2
!       if( lam(i) < 0 ) ieig = i
!    end do
    if( ieig == 0 ) &
         stop 'sr1mean:ridge(): an indefinite Hessian without negative lambda? '

    rtol = 1e-2
!    delta = min(3 * s / sqrt(real(size(funmean%data),REAL64)),0.5*s)
!    delta = 0.5*s
    write(*,*) 'Ridge start with ',t,s,' dir:',real(Q(:,ieig))

!    block
!      integer :: i
!      real(REAL64) :: u,v
!      open(17,file='/tmp/l')
!      do i = -20,20
!         u = t + 0.1*i*Q(1,1)
!         v = s + 0.1*i*Q(2,1)
!         if( v > 0 ) write(17,*) u,v,funmean%val(u,v)
!      end do
!      close(17)
!    end block

!    call smnewton(funmean,Q(:,ieig),t,s,tol,delta,reli)
    if( debug ) open(8,file='/tmp/sr1mean_ridge_step')
    call ridge_line(funmean,Q(:,ieig),t,s,rtol,.true.,hess,delta,reli)
    close(8)
    write(*,*) 'Ridge max at ', t,s,reli,delta
!    if( hypot(t - t0,2*(s-s0)) > 0.1*s0 ) return
!    reli = .false.
    if( reli ) then

       goto 44
       call funmean%hess(t,s,hess)
       call q2eig(hess,lam,Q)
       ieig = 0
       do i = 1, 2
          !       if( dot_product(Q(:,i),real([0,1],REAL64)) > qmin ) ieig = i
          if( abs(Q(2,i)) > qmin ) ieig = i
       end do
       write(*,*) ieig, dot_product(Q(:,ieig),real([0,1],REAL64))
       !    do i = 1, 2
       !       if( lam(i) < 0 ) ieig = i
       !    end do
       if( ieig == 0 ) &
            stop 'WARNING at sr1mean:ridge:' // &
            ' Hessian while passsing ridge is positive definite.'
       write(*,*) 'Ridge direction: ',real(Q(:,ieig))

44     continue
       nopt = 2
       d = 0.5*delta
       d = delta
       topt(1) = t - d*Q(1,ieig)
       topt(2) = t + d*Q(1,ieig)
       sopt(1) = s - d*Q(2,ieig)
       sopt(2) = s + d*Q(2,ieig)

!!$       block
!!$         real(REAL64) :: sp, tp, u,f0,f
!!$         real(REAL64), dimension(2) :: fvec,p
!!$         real(REAL64), dimension(2,2) :: B
!!$
!!$         open(8,file='/tmp/ri')
!!$         do k = 2,1,-1
!!$            d = (-1)**k*delta
!!$            do l = 1, 50
!!$               u = 0.1*l*d
!!$               tp = t + u*Q(1,ieig)
!!$               sp = s + u*Q(2,ieig)
!!$               if( sp > 0 ) then
!!$                  call funmean%grad(tp,sp,fvec)
!!$                  call funmean%hess(tp,sp,B)
!!$                  f = funmean%val(tp,sp)
!!$                  p = -[tp - t, sp - s]
!!$                  write(8,*) tp,sp,f,u,norm2(fvec),(f-f0)/(0.5*sum(p*matmul(B,p)))
!!$               end if
!!$            end do
!!$         end do
!!$         close(8)
!!$       end block

    else

       write(error_unit,*) 'WARNING at sr1mean:ridge: newton() did not converge.'
       t = t0
       s = s0
       delta = 0.1*s


       goto 55

       nopt = 2
       do k = 1,2
          d = (-1)**k*delta
          do l = -2, 4
             sp = s + 2**l*d*Q(2,ieig)
             if( sp > 0 ) then
                topt(k) = t + 2**l*d*Q(1,ieig)
                sopt(k) = sp
                call funmean%hess(topt(k),sopt(k),hess)
                call funmean%grad(topt(k),sopt(k),fvec)
                write(*,*) k,l,real(topt(k)),real(sopt(k)),real(dot_product(fvec,Q(:,ieig)))
                if( sopt(k) > 0 ) then
                   if( dot_product(fvec,Q(:,ieig)) < 0 ) exit
!                   if( hess(2,2) > 0 ) exit
!                   if( hess(1,1) > 0 .and. hess(1,1)*hess(2,2) - hess(1,2)**2 > 0 ) exit
                end if
             end if
          end do
       end do

    end if
55 continue

    t = t0
    s = s0
    fmax = funmean%val(t,s)
    ! the search for points where f < fmax, it's necessary
    ! the maximum is localised roughly off-ridge.
    do k = 1,2
       d = (-1)**k*delta
       topt(k) = t + d*Q(1,ieig)
       sopt(k) = s + d*Q(2,ieig)
       do l = 1, 10
          u = 0.5*l*d
          tp = t + u*Q(1,ieig)
          sp = s + u*Q(2,ieig)
          if( sp > 0 ) then
             if( funmean%val(tp,sp) < fmax ) then
                topt(k) = tp
                sopt(k) = sp
                exit
             end if
          end if
       end do
    end do

    write(*,*) 'Ridge final: topt=',topt(1:nopt),'sopt=',sopt(1:nopt)

    ! Note that the subroutine is always successful,
    ! even it is not so successful.

!    stop 'ridge'

  end subroutine ridge


  subroutine goldfinger(fcn,z,r0,qvec,dmax,gtol,alpha,fmax,reli)

    ! A maximalisation of a function: this is a particular variant
    ! of localmin() procedure by Brent (1972). The maximum is searched
    ! along the vector qvec in the interval [0..delta] with starting
    ! point r0. Result is the alpha in [0..1] pointing on a point meaning
    ! r = r0 + alpha*delta*qvec. No parabolic fit, and so checking of a,b
    ! vicinity, is performed.

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), dimension(:), intent(in) :: r0, qvec
    real(REAL64), intent(in) :: z, dmax, gtol
    real(REAL64), intent(out) :: alpha, fmax
    logical, intent(out) :: reli

    real(REAL64), parameter :: Kg = 1.44 != 1/log2((1+sqrt(5))/2)
    real(REAL64), parameter :: c = (3.0_REAL64 - sqrt(5.0_REAL64))/2.0_REAL64
    real(REAL64), parameter :: eps = sqrt(epsilon(c))
    real(REAL64), dimension(2) :: rp
    real(REAL64) :: a,b,e,d,fw,fv,fx,fu,v,w,u,x,m,tol,t2
    integer :: iter, niter
    logical :: cond1, cond2

    niter = int(2*Kg*(log(abs(dmax)/gtol)/log(2.0))**2) + 1
    write(*,*) 'goldfinger() init: niter=',niter,' r0=',real(r0), &
         ' dir.: ',real(qvec),' dmax=',real(dmax)

    a = 0
    b = 1
    x = a + c*(b - a)
    w = x; v = x;
    e = 0
    rp = r0 + x*dmax*qvec
    fx = fcn%val(rp(1),rp(2))
    fv = fx; fw = fx;
    tol = eps*abs(x) + gtol
    t2 = 2*tol

    do iter = 1, niter

       m = (a + b) / real(2,REAL64)

       if( debug ) write(*,'(a,i2,a,2f10.3)') 'goldfinger() ',iter,': ',x,fx

       if( max(x - a, b - x) < t2 )  then
          alpha = x
          fmax = fx
          reli = .true.
          return
       end if

       ! a "golden section" step
       if( x < m ) then
          e = b - x
       else
          e = a - x
       end if
       d = c*e

       u = x + d
       rp = r0 + u*dmax*qvec
       fu = fcn%val(rp(1),rp(2))

       if( z > 0 ) then
          cond1 = fu >= fx
       else
          cond1 = fu <= fx
       end if
!       if( fu >= fx ) then
       if( cond1 ) then
          if( u < x ) then
             b = x
          else
             a = x
          end if
          v = w; w = x; x = u;
          fv = fw; fw = fx; fx = fu
       else
          if( u < x ) then
             a = u
          else
             b = u
          end if
          if( z > 0 ) then
             cond1 = fu >= fw
             cond2 = fu >= fv
          else
             cond1 = fu <= fw
             cond2 = fu <= fv
          end if
!          if( fu >= fw .or. abs(w - x) < eps ) then
          if( cond1 .or. abs(w - x) < eps ) then
             v = w; w = u;
             fv = fw; fw = fu
!          else if( fu >= fv .or. abs(v - c) < eps .or. abs(v - w) < eps ) then
          else if( cond2 .or. abs(v - c) < eps .or. abs(v - w) < eps ) then
             v = u
             fv = fu
          end if
       end if

    end do

    reli = .false.

  end subroutine goldfinger


  subroutine sr1update(p,y,snorm,B)

    real(REAL64), dimension(:), intent(in) :: p,y
    real(REAL64), intent(in) :: snorm
    real(REAL64), dimension(:,:), intent(in out) :: B

    real(REAL64), parameter :: rpar = sqrt(epsilon(rpar))  ! 0 < r < 1

    real(REAL64), dimension(size(p)) :: res
    real(REAL64) :: rnorm, q, r
    integer :: n,i,j

    n = size(p)
    res = y - matmul(B,p)
    rnorm = norm2(res)
    q = dot_product(p,res)
    r = snorm*rnorm
!    write(*,*) real(res)
!    write(*,*) real(y)
!    write(*,*) real(q)

    if( abs(q) >= rpar*r .and. r > 0 ) then
       forall( i = 1:n, j = 1:n ) B(i,j) = B(i,j) + res(i)*res(j)/q
    else
       B = B
    end if

  end subroutine sr1update

  subroutine rscale(fcn,t,sig,s,iflag)

    class(tsr1mean_function), pointer, intent(in) :: fcn
    real(REAL64), intent(in) :: t,sig
    real(REAL64), intent(out) :: s
    integer, intent(out) :: iflag

    real(REAL64) :: a,b,x,fa,fb,f,fzero
    integer :: i
    logical :: reli

    iflag = 666
    fzero = sig / sqrt(real(size(fcn%data)))

    a = 0.5
    fa = fcn%stderr(t,a*sig,reli) - fzero
    if( .not. reli ) return
    b = 1.0
    fb = fcn%stderr(t,b*sig,reli) - fzero
    if( .not. reli ) return

    do i = 1, 2
       x = (a + b) / 2
       f = fcn%stderr(t,x*sig,reli) - fzero
       if( .not. reli ) return
       write(*,*) 'rscale:',real(x),real(f)
       if( fa*f < 0 ) then
          b = x
          fb = f
       else
          a = x
          fa = f
       end if
    end do
    x = (a+b)/2
    s = x*sig
    iflag = 0

    if( debug ) then
       f = fcn%stderr(t,x*sig,reli) - fzero
       write(*,*) 'rscale final parameter:',real(x),real(f),reli

!       open(8,file='/tmp/xxrs')
!       do i = 50, 100
!          x = 0.01*i
!          write(8,*) x,fcn%stderr(t,x*sig,reli) -fzero, fcn%stderr(t,x*sig,reli)
!       end do
!       close(8)
    end if

  end subroutine rscale


  subroutine iscale(data,t,scale)

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t
    real(REAL64), intent(in out) :: scale

    real(REAL64), parameter :: rtol = 1e-2
    real(REAL64) :: s, fval, fgrad, df, d, p, sp, fv, fg, fd, ared, pred, rho, dlt
    integer :: iter

    s = scale
    fval = chi2(data,t,s)
    fgrad = chi2_der(data,t,s)
    df = chi2_der2(data,t,s)
    dlt = 0.5*s

    write(*,*) 'iscale init:',t,s,dlt

    do iter = 1, 3*precision(s)

       if( df > 0 ) then
          d = -fgrad/df
          if( abs(d) > dlt ) d = -sign(dlt,fgrad)
       else
          d = sign(dlt,fgrad)
       end if
       p = d

       write(*,'(a,i2,a,2f10.3,2es10.2,a,2f6.3)') &
            'iscale() ',iter,': ',s,d,df,fval,' | ',rho,dlt

       if( abs(d) < rtol*s .and. df > 0 ) then
          scale = s
          return
       end if

       sp = s + p

       if( sp > 0 ) then

          fv = chi2(data,t,sp)
          fg = chi2_der(data,t,sp)
          fd = chi2_der2(data,t,sp)

          ared = fval - fv
          pred = -(fg*p + 0.5*fd*p**2)
          pred = -(fg*p + 0.5*df*p**2)
          if( abs(ared) < epsilon(ared) .and. abs(pred) < epsilon(pred) ) then
             rho = 1
          else if( abs(pred) < epsilon(pred) .and. abs(ared) > 0 ) then
             rho = 0
          else
             rho = ared / pred
          end if

       else
          rho = 0
       end if

       if( rho > eta ) then
          s = sp
          fval = fv
          fgrad = fg
          df = fd
       end if

       if( rho > 0.75 ) then
          if( abs(d) > 0.9*dlt ) then
             dlt = 2*dlt
          end if
       else if( rho < 0.25 ) then
          dlt = dlt / 2
       end if

    end do


  end subroutine iscale

  real(REAL64) function chi2(data,t,s)

    use rfun_hampel

    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t,s

    chi2 = sum(pshampel((data-t)/s)**2)

  end function chi2

  real(REAL64) function chi2_der(data,t,s)

    use rfun_hampel
    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t,s

    chi2_der = -(2/s)*sum(pshampel((data-t)/s)*dpshampel((data-t)/s)*((data-t)/s))
!    write(*,*) 'c',chi2_der
!    chi2_der = (chi2(data,t,s+1e-6) - chi2(data,t,s)) / 1e-6
!    write(*,*) 'd',chi2_der

  end function chi2_der

  real(REAL64) function chi2_der2(data,t,s)

    use rfun
    real(REAL64), dimension(:), intent(in) :: data
    real(REAL64), intent(in) :: t,s
    real(REAL64), parameter :: d = 1e-6
!    chi2_der2 = 4*sum(dpshampel((data-t)/s)**2 * (data-t)**2)/s**4
!    write(*,*) 'a',chi2_der2
    chi2_der2 = (chi2_der(data,t,s+d) - chi2_der(data,t,s))/d
!    write(*,*) 'b',chi2_der2

  end function chi2_der2



!!$  real(REAL64) function rk4(fcn,t,s,h)
!!$
!!$    class(tsr1mean_function), pointer, intent(in) :: fcn
!!$    real(REAL64), intent(in) :: t,s,h
!!$
!!$!    real(REAL64), dimension(2) :: fvec
!!$    real(REAL64) :: h2,k1,k2,k3,k4
!!$
!!$    h2 = h / 2
!!$
!!$!    call fcn%grad(t,s,fvec)
!!$    !    k1 = -fvec(1) / fvec(2)
!!$    k1 = iR(fcn,t,s)
!!$
!!$!    call fcn%grad(t + h2,s+h2*k1,fvec)
!!$    !    k2 = -fvec(1) / fvec(2)
!!$    k2 = iR(fcn,t + h2,s + h2*k1)
!!$
!!$!    call fcn%grad(t + h2,s+h2*k2,fvec)
!!$    !    k3 = -fvec(1) / fvec(2)
!!$    k3 = iR(fcn,t + h2,s + h2*k2)
!!$
!!$!    call fcn%grad(t+h,s+h2*k3,fvec)
!!$    !    k4 = -fvec(1) / fvec(2)
!!$    k4 = iR(fcn,t + h, s + h*k3)
!!$
!!$    rk4 = s + h*(k1 + 2*(k2 + k3) + k4)/6;
!!$
!!$  end function rk4
!!$
!!$  real(REAL64) function iR(fcn,t,s)
!!$
!!$    class(tsr1mean_function), pointer, intent(in) :: fcn
!!$    real(REAL64), intent(in) :: t,s
!!$    real(REAL64), dimension(2,2) :: hess
!!$
!!$    call fcn%hess(t,s,hess)
!!$    iR = -hess(1,2) / hess(2,2)
!!$
!!$  end function iR


!!$  function iscale(fcn,t) result(s)
!!$
!!$    use fmm
!!$    use quantilmean
!!$    use global_minimization
!!$
!!$    class(tsr1mean_function), pointer, intent(in) :: fcn
!!$    real(REAL64), intent(in) :: t
!!$    real(REAL64), dimension(:), allocatable :: xcdf, ycdf
!!$    real(REAL64), dimension(2,2) :: hess
!!$    real(REAL64) :: s,smin,smax,tol,smad,f,mdiff2
!!$    character(len=80) :: errmsg
!!$    integer :: stat,n,m1,m2,m
!!$    logical :: found
!!$
!!$    iscale_fcn => fcn
!!$    iscale_t = t
!!$    iscale_count = 0
!!$
!!$    n = size(fcn%data)
!!$    allocate(xcdf(n),ycdf(n),stat=stat,errmsg=errmsg)
!!$    if( stat /= 0 ) then
!!$       write(error_unit,*) 'scale: ',trim(errmsg)
!!$       stop
!!$!       if( present(flag) ) flag = 5
!!$!       if( present(reliable) ) reliable = .false.
!!$!       return
!!$    end if
!!$
!!$    ! cummulative distribution function
!!$    call ecdf(abs(fcn%data - t),xcdf,ycdf)
!!$
!!$    if( mod(n,2) == 0 ) then
!!$       m = n / 2
!!$       smad = (xcdf(m) + xcdf(m+1)) / 2
!!$    else
!!$       smad = xcdf((n+1)/2)
!!$    end if
!!$    smad = smad / 0.6745
!!$
!!$    if( smad < epsilon(smad) ) stop 'identical points'
!!$
!!$!    h = real(1.0,REAL64) / real(n + 1,REAL64)
!!$!    m2 = min(nint(0.8 * (n+1)),n)
!!$    m2 = n - int((1 - 0.9) * (n+1)) ! q = 80%
!!$    m1 = max(int(0.4 * (n+1)),2) ! q = 40%
!!$    m = n - int((1-0.68) * (n+1)) ! q = 68%
!!$
!!$    smin = xcdf(m1)
!!$    smax = xcdf(m2)
!!$    smad = xcdf(m)
!!$!    smad = smax
!!$
!!$!    if( smin < 0.1*smad ) then
!!$!       smin = 0.5*smad
!!$!       do i = 1, n
!!$!          if( xcdf(m1) > 0 ) then
!!$!             m1 = i
!!$!             exit
!!$!          end if
!!$!       end do
!!$!    end if
!!$
!!$!    if( smax > 10*smad ) then
!!$!       smax = 1.5*smad
!!$!    end if
!!$
!!$    call fcn%hess(t,smin,hess)
!!$
!!$    write(*,*) 'smin,smax,smad:',smin,smax,smad,m1,m2,m,hess(2,2)
!!$
!!$    tol = 0.2 *(smax + smin) / 2
!!$    tol = 0.1*smad
!!$    s = fmin(smin,smax,scalelike,tol)
!!$    write(*,*) 'by fmin():',s
!!$    mdiff2 = 1.1*abs(hess(2,2))
!!$    mdiff2 = 10
!!$    iscale_count = 0
!!$    f = glomin(smin,smax,smad,mdiff2,tol,tol,scalelike,s)
!!$
!!$    ! Does the minimum falls inside the interval smin < .. < smax ?
!!$    found = smin + tol < s .and. s < smax - tol
!!$!    s = smad
!!$
!!$    !    if( debug )
!!$    write(*,'(a,g20.6,a,g0.2,l2,i5)') &
!!$         'Initial scale estimate:',s," +- ",tol,found,iscale_count
!!$
!!$!    s = xcdf((n+1)/2) / 0.6746
!!$!    write(*,*) 'Used estimate:',s
!!$
!!$  end function iscale
!!$
!!$  function scalelike(s)
!!$
!!$    use rfun
!!$    ! negative of likelihood per data point
!!$
!!$    real(REAL64) :: scalelike
!!$    real(REAL64), intent(in) :: s
!!$    real(REAL64) :: r,q,x
!!$    integer :: i,n
!!$
!!$    if( .not. (s > epsilon(s)) ) then
!!$       scalelike = min(1e3*maxval(abs(iscale_fcn%data - iscale_t))**2,0.1*huge(scalelike))
!!$       return
!!$    end if
!!$
!!$!    scalelike = iscale_fcn%val(iscale_t,s)
!!$
!!$    n = size(iscale_fcn%data)
!!$    q = 0
!!$    do i = 1, n
!!$       r = (iscale_fcn%data(i) - iscale_t) / s
!!$       if( r > 2 ) then
!!$          x = 2
!!$       else if( r < -2 ) then
!!$          x = -2
!!$       else
!!$          x = r
!!$       end if
!!$       q = q + ihuber(x)
!!$!       write(*,*) i, iscale_fcn%data(i), r, x, ihuber(x)
!!$    end do
!!$    scalelike = q / n + log(s)
!!$
!!$    scalelike = sum(ipshampel((iscale_fcn%data - iscale_t) / s))/n !+ asinh(s/2)
!!$!    scalelike = sum(itukey((iscale_fcn%data - iscale_t) / s))/n !+ asinh(s/2)
!!$    if( s < exp(1.0/3.0) ) then
!!$       scalelike = scalelike + 0.1226*s**3
!!$    else
!!$       scalelike = scalelike + log(s)
!!$    end if
!!$
!!$    iscale_count = iscale_count + 1
!!$    write(*,*) s,scalelike,iscale_count
!!$
!!$  end function scalelike



end module trustsr1mean
