!
!  Robust estimate of both scale and background from already known growth-curve.
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
!

module grow_fit

  implicit none

  ! numerical precision of real numbers
  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: pi = 3.14159265358979312_dbl

  ! use analytical Jacobian ?
  logical, parameter, private :: analytic = .true.

  ! verbose (prints intermediate results)
  logical, parameter, private :: verbose = .false.

  ! In minpacks routines, prints intermediate results every nprint-call or be quite (0)
  integer, private :: nprint = 0 ! 0 or >= 1

  integer, parameter :: maxit = 3*precision(0.0_dbl)
  ! Number of iterations is limited by numerical precision.
  ! We belives, that one order is reached by each iteration.

  ! tolerances in numerical precision of iterations
  real(dbl), parameter, private :: tol = 10*epsilon(0.0_dbl)

  real(dbl), dimension(:), allocatable, private :: aper,cts,dcts,curve,curve_err,acts,ects
  real(dbl), private :: rsig, fmin, ctsref
  ! the variance of the mean sky of the star (the sky subtracted from the
  ! apertures), negative when unknown
  real(dbl), private :: skyvar

  private :: growfun, growdif, growder, difjac, res, varpar, &
       snoise, snoise_graph

contains

  subroutine growfit(raper,xcts,xdcts,xcurve,xcurve_err,skyerr,f,df,b,db,sig,info)

    ! skyerr is the error of the mean sky subtracted from the apertures
    ! (SKYERR; <= 0 when unknown); db is then the error of the corrected
    ! sky, sky + b (see varpar)

    use NelderMead
    use minpacks
    use oakleaf

    real, dimension(:), intent(in) :: raper
    real(dbl), dimension(:), intent(in) :: xcts,xdcts,xcurve,xcurve_err
    real, intent(in) :: skyerr
    real(dbl), intent(in out) :: f,df,b,db
    real(dbl), intent(out) :: sig
    integer, intent(out) :: info
    integer :: ifault
    real(dbl), dimension(2) :: p,p1,dp
    real(dbl), dimension(size(raper)) :: r
    real(dbl) :: mad, dsig
    integer :: n, iter, i
    logical :: converge

    info = 3 ! no initialisation

    n = size(raper)
    allocate(aper(n),cts(n),dcts(n),curve(n),curve_err(n),acts(n),ects(n))
    aper = pi*raper**2
    cts = xcts
    dcts = xdcts
    curve = xcurve
    curve_err = xcurve_err
    ctsref = f
    if( skyerr > 0 ) then
       skyvar = real(skyerr,dbl)**2
    else
       skyvar = -1
    end if
    acts = cts / ctsref
    ects = dcts / ctsref

    if( verbose ) then
       do i = 1,n
          write(*,*) i,real(curve(i)),real(curve_err(i)),real(cts(i)/curve(i))
       end do
    end if

!    write(*,*) real(p(1)),real(dp(1))
!    write(*,*) real(cts - p(1)*curve)/aper

    ! setup init values
    p(1) = f
    p(2) = 0
    dp(1) = 1e-4*df + f*curve_err(6)
    dp(2) = 1

    p(1) = 1
    p(2) = 0
    dp = 1e-4

    ! minimal value of p(1) for init
    fmin = max(0.1*median(cts),1.0)
    fmin = minval(acts)

    ! initial estimate of parameters by absolute deviations
    call nelmin1(growfun,p,dp,mad,ifault)
    if( verbose ) &
         write(*,'(a,3g15.5,f8.3,i5)') 'growfit init: ',p*ctsref,mad,mad/p(1),ifault
    if( ifault /= 0 .or. p(1) < 1.01*fmin ) then
       info = 2 ! no convergence
       if( verbose ) write(*,*) &
            "Warning: growfit finished prematurely without nelmin convergence."
       goto 666
    end if
    p = ctsref*p

!    p(2) = 0
!    dp(2)

    ! initial estimate of mean deviation
    r = abs(p(1)*curve + p(2)*aper - cts) / sqrt(dcts**2 + p(1)**2*curve_err**2)
    sig = median(r) / 0.6745
    dsig = median(sqrt(dcts**2 + p(1)**2*curve_err**2)) / 0.6745
    if( verbose ) write(*,'(a,g0.3,5x,2f0.5,a)') 'sig = ',sig,dsig
    rsig = sig

    converge = .false.
    p1 = p
    do iter = 1, maxit

       ! robust estimate
       if( analytic ) then
          call lmder2(growder,p,tol,nprint,info)
       else
          call lmdif2(growdif,p,sqrt(tol),nprint,info)
       end if
       if( verbose ) write(*,'(a,es12.5,f10.5,i5)') 'robust:',p,info
       if( info == 5 ) then
          info = 2 ! no convergence
          if( verbose ) write(*,*) &
               "Warning: growfit finished prematurely without lmdif/lmder convergence."
          if( analytic ) goto 666
       end if

       converge = all(abs(p - p1) < tol)

       if( converge ) exit

       p1 = p

       ! update noise scale
       call snoise(p,sig)
!       rsig = sig

    end do

    if( .not. converge .and. verbose ) &
         write(*,*) "Warning: Grow-curve fitting shows no convergence!"

    if( converge ) then

       ! estimation of uncertainties: the covariance of the fit
       ! propagated from the errors of the (correlated) aperture counts
       ! (both the photon and the sky noise are included in dcts)
       call varpar(p,dp,sig)

       info = 0
       f = p(1)
       b = p(2)
       ! keep initial estimates when the covariance is singular
       if( dp(1) > 0 ) df = dp(1)
       if( dp(2) > 0 ) db = dp(2)
    end if

666 continue

    deallocate(aper,cts,dcts,curve,curve_err,acts,ects)

  end subroutine growfit

  function growfun(p) result(s)

    real(dbl), dimension(:), intent(in) :: p
    real(dbl) :: s,f,b,f2ref

    f = p(1)
    b = p(2)
!    f2ref = cts(6)**2

    if( f < fmin ) then
       s = 100 * size(cts) !* maxval(cts)
    else
       s = sum(abs(f*curve + b*aper - acts)/sqrt(ects**2 + curve_err**2))
!       s = sum(abs(curve + b*aper - cts/f)/sqrt(dcts**2/f2ref + curve_err**2))
!       s = sum(abs(curve + b*aper - cts/f))
    end if
!    write(*,*) real(p),real(s)

  end function growfun


  subroutine res(f,b,r,ds)

    real(dbl), intent(in) :: f,b
    real(dbl), dimension(:), intent(out) :: r,ds

    r = cts - (f*curve + b*aper)
    ds = sqrt(dcts**2 + f**2*curve_err**2 + rsig**2*f**2*0)

  end subroutine res

  subroutine growdif(m,np,p,fvec,iflag)

    use oakleaf

    integer, intent(in) :: m,np
    integer, intent(in out) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(out) :: fvec
    real(dbl), dimension(:), allocatable :: r,f,ds
    integer :: n
    real(dbl) :: f0,b,s

    if( iflag == 0 ) then
       write(*,'(2g13.5)') p
       write(*,'(2g13.5)') fvec
       return
    end if

    n = size(aper)
    allocate(r(n),f(n),ds(n))

    f0 = p(1)
    b = p(2)
    s = rsig

    call res(f0,b,r,ds)
    f = huber(r/ds/s)

    fvec(1) = sum(f*(curve - f0*r*curve_err**2/ds**2)/ds - f0*curve_err**2/ds**2*s)
    fvec(2) = sum(f*aper/ds)

    fvec = - fvec / s

    deallocate(r,f,ds)

  end subroutine growdif

  subroutine growder(m,np,p,fvec,fjac,ldfjac,iflag)

    use oakleaf

    integer, intent(in) :: m,np,ldfjac
    integer, intent(in out) :: iflag
    real(dbl), dimension(np), intent(in) :: p
    real(dbl), dimension(m), intent(out) :: fvec
    real(dbl), dimension(ldfjac,np), intent(out) :: fjac
    real(dbl), dimension(size(p),size(p)) :: dfjac
    real(dbl), dimension(:), allocatable :: r,f,df,ds,dd
    integer :: i,j,n
    real(dbl) :: f0,b,s

    if( iflag == 0 ) then

       write(*,'(6g13.5)') p,fvec

       do i = 1,2
          write(*,'(a,2g15.5)') ' jac:',fjac(i,:)
       end do

       call difjac(p,dfjac)
       do i = 1,2
          write(*,'(a,2g15.5)') 'djac:',dfjac(i,:)
       end do

       return
    end if

    n = size(aper)
    allocate(r(n),f(n),df(n),ds(n),dd(n))

    f0 = p(1)
    b = p(2)
    s = rsig

    call res(f0,b,r,ds)
    f = huber(r/ds/s)
    dd = curve - r*f0*curve_err**2/ds**2

    if( iflag == 1 ) then

       fvec(1) = sum(f*dd/ds - f0*curve_err**2/ds**2*s)
       fvec(2) = sum(f*aper/ds)

    else if( iflag == 2 ) then

       fjac = 0
       df = dhuber(r/ds/s)

       fjac(1,1) = sum(df*(dd/ds)**2 &
            - f*(curve*f0*curve_err**2 + &
            r*curve_err**2*(ds**2-3.0_dbl*f0**2*curve_err**2)/ds**2)/ds**3 &
            - curve_err**2*(ds**2-2.0_dbl*f0**2*curve_err**2)/ds**4)
       fjac(1,2) = sum(df*dd*aper/ds**2 - f*aper*f0*curve_err**2/ds**3)
       fjac(2,2) = sum(df*(aper/ds)**2)

       do j = 1,size(fvec)
          do i = 1,j-1
             fjac(j,i) = fjac(i,j)
          end do
       end do

    end if

    fvec = - fvec / s
    fjac = fjac / s**2

    deallocate(r,f,df,ds,dd)

  end subroutine growder

  subroutine difjac(p,jac)

    ! It has failure near the minimum

    real(dbl), dimension(:), intent(in) :: p
    real(dbl), dimension(:,:), intent(out) :: jac
    real(dbl), parameter :: d = 1e-4
    real(dbl), dimension(size(p)) :: fv1,fv2,q
    integer :: n,i,iflag

    iflag = 1
    n = size(p)

    do i = 1,size(p)
       q = p
       q(i) = p(i) + d
       call growdif(n,n,q,fv1,iflag)

       q = p
       q(i) = p(i) - d
       call growdif(n,n,q,fv2,iflag)

       jac(i,:) = (fv1 - fv2)/(2*d)
    end do

  end subroutine difjac



  subroutine varpar(p,dp,sig)
    ! estimate variability of parameters

!    use minpacks
    use oakleaf

    real(dbl), dimension(:), intent(in) :: p
    real(dbl), dimension(:), intent(out) :: dp
    real(dbl), intent(out) :: sig

    real(dbl), dimension(size(p),size(p)) :: jac,hess,cov
    real(dbl),dimension(:), allocatable :: f,df,r,ds,w
    real(dbl),dimension(:,:), allocatable :: x,g,c
    real(dbl) :: sum2,sum3
    integer :: i,j,k,n

    n = size(aper)
    allocate(r(n),f(n),df(n),ds(n),w(n),x(n,2),g(2,n),c(n,n))
    dp = -1

    call res(p(1),p(2),r,ds)
    f = huber(r/(ds*rsig))
    df = dhuber(r/(ds*rsig))
    sum2 = sum(df)
    sum3 = sum(f**2)

    if( sum2 > 0 .and. sum3 > 0 .and. n > 2 ) then

       ! Near the solution, the robust estimator behaves like the weighted
       ! least squares with weights w = psi'(r/ds) / ds**2 (the matrix jac
       ! is the second derivative of the robust function).
       w = df / ds**2
       x(:,1) = curve
       x(:,2) = aper
       do i = 1,2
          do j = 1,2
             jac(i,j) = sum(w*x(:,i)*x(:,j))
          end do
       end do
       call qrinv(jac,hess)

       ! The apertures are nested: the counts in an aperture include all
       ! the pixels (and their noise) of smaller ones, so the measurements
       ! are strongly correlated. The covariance of the counts in apertures
       ! (sorted by radius) is approximated by cov(i,j) = dcts(min(i,j))**2
       ! (photon noise of the star and the sky, the sky estimate error).
       ! The covariance of the parameters is given by the sandwich formula
       !   cov(p) = H X' W C W X H,  H = (X' W X)**(-1).
       ! With the error of the mean sky known (skyvar), its part of dcts**2,
       ! area**2*skyvar, is replaced by the exact one, area(i)*area(j)*skyvar,
       ! common to all apertures: X' W X H = I, so it is absorbed by b
       ! entirely (b = -(sky error) + the pixel noise) and does not enter f;
       ! the covariance of the pixels alone gives the error of f and the
       ! error of the corrected sky, sky + b, in which the error of the
       ! mean of the ring cancels (it was added to it once more: SKYERR
       ! of gphot was 24-41% too large).
       do i = 1,n
          do j = 1,n
             k = min(i,j)
             if( skyvar >= 0 ) then
                c(i,j) = max(dcts(k)**2 - aper(k)**2*skyvar,0.0_dbl)
             else
                c(i,j) = dcts(k)**2
             end if
          end do
          c(i,i) = c(i,i) + (p(1)*curve_err(i))**2
       end do
       do i = 1,2
          g(i,:) = (hess(i,1)*x(:,1) + hess(i,2)*x(:,2))*w
       end do
       cov = matmul(g,matmul(c,transpose(g)))

       ! Huber (6.6), the scale of normalised residuals (informative)
       sig = rsig*sqrt(sum3/sum2*n/(n-1.0))

       do i = 1, 2
          if( cov(i,i) > 0 ) then
             dp(i) = sqrt(cov(i,i))
          end if
       end do

       if( verbose ) then
          do i = 1,n
             write(*,'(i3,g15.5,3f12.3)') i,r(i),r(i)/(ds(i)*rsig),f(i),df(i)
          end do
          write(*,*) 'jac:',real(jac(1,:))
          write(*,*) 'jac:',real(jac(2,:))
          write(*,*) 'cov:',real(cov(1,:))
          write(*,*) 'cov:',real(cov(2,:))
       end if

    end if

    deallocate(r,f,df,ds,w,x,g,c)

  end subroutine varpar

  subroutine snoise(p,sig)

    use oakleaf

    real(dbl), dimension(:), intent(in) :: p
    real(dbl), intent(in out) :: sig

    integer :: n,iter
    real(dbl) :: s,sum1,sum2,fs,dfs,d
    real(dbl), dimension(:), allocatable :: r,rs,ds,f,df,rho,erho
    logical :: convergent

    n = size(aper)
    allocate(r(n),rs(n),ds(n),rho(n),df(n),f(n),erho(n))

    s = sig
    convergent = .false.

    do iter = 1, maxit

       if( .not. (s > 0) ) exit

       call res(p(1),p(2),r,ds)
       r = r / ds
       rs = r / s
       rho = ihuber(rs)
       f   =  huber(rs)
       df  = dhuber(rs)

!       rho = rs**2/2
!       f = rs
!       df = 1

       erho = exp(-2*rho)
       sum1 = sum(f*r*(1-2*rho)*erho)
       sum2 = sum(((1-2*rho)*(2*f**2-df) + 2*f**2)*r**2*erho)

       fs = -sum1 / s**2
       dfs = 2*sum1 / s**3 - sum2 / s**4

       ! Note. If we have good initial estimate, than sum1 (=fs) is near to zero
       ! and sum1 / sum2 * s**2 is numerically equivalent, but little bit
       ! faster, than fs/dfs. Theirs numerical differences are negligible.

       if( .not. (abs(dfs) > 0) ) exit

       ! Newton's step for scale
       d = fs / dfs
       s = s - d

       !       if( debug )
!       write(*,'(a,i3,4g15.5)') "scale,f,f',incr.",iter,s,fs,dfs,d

       ! exit immediately when required precision is reached
       convergent = abs(d) < 10*epsilon(sig)
       if( convergent ) exit

    end do

    deallocate(r,rs,ds,f,df,rho,erho)

    if( convergent ) sig = s

  end subroutine snoise

  subroutine snoise_graph(p)

    use oakleaf

    real(dbl), dimension(:), intent(in) :: p

    integer :: n,i
    real(dbl) :: s
    real(dbl), dimension(:), allocatable :: r,rho,ds,rs

    n = size(aper)
    allocate(r(n),rho(n),ds(n),rs(n))

    open(1,file='/tmp/e')
    do i = 2,1000
       s = i / 100.0

       call res(p(1),p(2),r,ds)
       rho = ihuber(r/(s*ds))
       write(1,*) s,sum(rho*exp(-2*rho))
    end do
    close(1)

    deallocate(r,rho,ds,rs)

  end subroutine snoise_graph

end module grow_fit
