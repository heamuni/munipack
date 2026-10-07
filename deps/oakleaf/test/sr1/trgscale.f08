
module trust_global_scale

  use sr1mean_function
  use iso_fortran_env
  implicit none

  logical, parameter :: debug = .true.

  public :: trgscale

contains

  subroutine trgscale(fcn,t,smin,smax,sig,s,nopt,sopt,reli)

    class(tsr1mean_function), intent(in) :: fcn
    real(REAL64), intent(in) :: t,smin,smax,sig
    real(REAL64), intent(in out) :: s
    integer, intent(out) :: nopt
    real(REAL64), dimension(:), intent(out) :: sopt
    logical, intent(out) :: reli

    integer, parameter :: maxopt = 13
!    real(REAL64), dimension(maxopt) :: sopt
    integer, dimension(maxopt) :: mopt
    real(REAL64) :: stol, f, fopt, s0
    logical :: converge, valid, ismin
    integer :: n, i, nedge, edge, iter, minmax, dir

    reli = .false.
    n = size(fcn%data)
    minmax = 1 ! minimum
    stol = sig / 10
    dir = 1
    nopt = 0
    nedge = 0
    s0 = s

    write(*,*) 'sig,stol:',sig,stol,size(fcn%data),' smin,smax:',smin,smax

    if( debug ) open(17,file='/tmp/trgscale')

    do iter = 1, maxopt

       call snewton(fcn,t,smin,smax,stol,sig,minmax,stol,s,edge,converge,ismin)
       if( debug ) write(*,*) 'snewton:',s,minmax,edge,converge

       if( converge .and. edge == 0 .and. ismin ) then
          if( nopt < maxopt ) then
             nopt = nopt + 1
             sopt(nopt) = s
             mopt(nopt) = minmax
          end if
       end if

       if( iter == 1 ) dir = sign(real(1,REAL64),s - s0)
!       write(*,*) iter,dir,s0-s

       if( edge /= 0 ) then
          nedge = nedge + 1
          dir = -edge

          if( nedge == 2 ) then

             fopt = huge(fopt)
             do i = 1, nopt
                if( mopt(i) > 0 ) then
                   !                   f = abs(fcn%nfrac(t,sopt(i)) - 0.6)
                   f = fcn%stderr(t,sopt(i),valid)
                   f = fcn%val(t,sopt(i))
                   valid = fcn%nfrac(t,sopt(i)) > 0.25
                   f = abs(fcn%nfrac(t,sopt(i)) - 0.9)
                   write(*,*) '*',f,fopt,sopt(i),fcn%nfrac(t,sopt(i)),valid
                   if( f < fopt .and. valid ) then
                      fopt = f
                      s = sopt(i)
                   end if
                end if
             end do
             reli = fopt < huge(fopt)
             exit
          end if

          do i = 1, nopt
             if( edge == 1 ) then
                if( sopt(i) < s ) then
                   s = sopt(i)
                   minmax = mopt(i)
                end if
             else if( edge == -1 ) then
                if( sopt(i) > s ) then
                   s = sopt(i)
                   minmax = mopt(i)
                end if
             end if
          end do
       end if

       minmax = -minmax
!       write(*,*) s,dir,stol
       s = s + dir * sig

    end do

    if( debug ) close(17)

    write(*,*) nedge, fopt
!    ns = nopt
!    ss(1:ns) = sopt(1:ns)

  end subroutine trgscale



  subroutine snewton(fun,t,smin,smax,dltmin,dltmax,minmax,stol,s,edge,converge,ismin)

    class(tsr1mean_function), intent(in), target :: fun
    real(REAL64), intent(in) :: t, stol, smin,smax, dltmin,dltmax!,r1,rf,best_r
!    real(REAL64), dimension(:), intent(in) :: xcdf
    integer, intent(in) :: minmax
!    logical, intent(in) :: type
    real(REAL64), intent(in out) :: s
    integer, intent(out) :: edge
    logical, intent(out) :: converge, ismin
!   logical, intent(out) ::

    real(REAL64), parameter :: dlteps = 2*epsilon(s)
    real(REAL64) :: fval,fv,fs,h22,d,rho,sp,ared,pred,dlt,one,r2,dtol
    integer :: iter
    logical :: dltlim

!    if( s < smin ) then
!       edge = .true.
!       return
!    end if

!    write(*,*) '.'

    converge = .false.
    ismin = .false.
    dlt = 5*dltmin
!    edge = .false.
    edge = 0
    rho = 0
    if( smin < s .and. s < smax ) then
       fval = fun%val(t,s)
       call newstep(fun,t,s,fs,h22)
       if( debug ) write(17,*) s,fval
    else
       if( s < smin ) edge = -1
       if( s > smax ) edge = 1
       return
    end if


    one = real(minmax,REAL64) ! numerical type conversion
!    if( type ) then
!       ! minimum
!       one = 1
!    else
!       ! maximum
!       one = -1
!    end if

    do iter = 1, 5*precision(s)

       if( one*h22 > 0 ) then
          d = -fs / h22
       else
          d = dlt
       end if

!!$       if( type ) then
!!$          ! minimum
!!$          if( h22 > 0 ) then
!!$             d = -fs / h22
!!$          end if
!!$       else
!!$          ! maximum
!!$          if( h22 < 0 ) then
!!$             d = -fs / h22
!!$          end if
!!$       end if

!!$       if( abs(d) > dlt ) then
!!$          d = -one*sign(dlt,fs)
!!$          if( type ) then
!!$             d = -sign(dlt,fs)
!!$          else
!!$             d = sign(dlt,fs)
!!$          end if
!!$       end if

       if( abs(d) >= dlt ) d = -one*sign(dlt,fs)

       sp = s + d

       if( debug ) &
            !            write(*,'(a,i2,1pg13.5,1p2g10.2,1pg11.3,1pe12.3,3l1)') &
            write(*,*) &
            "trgscale,inc.,f,f': ",&
            iter,s+d,d,fs,h22,rho,abs(abs(d)-dlt)<2*epsilon(dlt),abs(fs/h22)<dlt,&
            smin < sp

!       dltlim = dlt < stol
!       if(  smin + 0*stol > sp .and. dltlim ) then
!          edge = -1
!          exit
!       else if( smax - 0*stol < sp .and. dltlim ) then
!          edge = 1
!          exit
!       end if


       dltlim = abs(dlt - dltmax) < dlteps .or. abs(dlt - dltmin) < dlteps
!       edge = (smin > sp .or. smax < sp ) .and. dltlim
!       write(*,*) sp,smin + stol,dlt,stol
!       edge = (smin + stol > sp .or. smax - stol < sp) .and. dlt < stol +dlteps
!       edge = smin + 1.1*stol > sp .and. dlt < 5.1*stol
!       edge = (sp < smin .or. sp > smax) .and. abs(d) < 3*dltmin
       !       if( abs(d) < 3*dltmin ) then
       dtol = 3*dltmin
       if( abs(d) < dtol .and. (sp < smin+dtol .or. sp > smax-dtol) ) then
          if( sp < smin ) edge = -1
          if( sp > smax ) edge = 1
!          converge = .true.
          exit
       end if
!       if( edge ) then
!          converge = .true.
!          exit
!       end if

!       if( mod(iter,7) == 0 ) then
!          r2 = sum(abs(Phi(xcdf,t,sp) - 0.5)) - size(xcdf) / 2.0
!          if( abs(r1 - r2)/rf > abs(best_r) ) exit
!       end if



!       if( sp > smin ) then
       if( smin < sp .and. sp < smax ) then

          fv = fun%val(t,sp)
          if( debug ) write(17,*) s,fval

          ared = fval - fv
          pred = -d*(fs + 0.5*h22*d)

!          if( abs(pred) > 0 ) then
             rho = ared / pred
!          else
!             rho = 0
!          end if

       else
          rho = 0
       end if

       if( rho > 0.5e-3 ) then
          s = sp
          fval = fv

!          write(*,*) '>',abs(s - smin), abs(s - smax), stol, d

          ismin = h22 > 0
          converge = abs(d) < stol .and. one*h22 > 0
          if( converge ) exit

          call newstep(fun,t,s,fs,h22)

       end if




       write(*,*) dlt,rho,d*fs, 0.5*h22*d**2
       if( rho > 0.75 ) then
          if( abs(d) > 0.9*dlt ) then
!             dlt = min(2*dlt,dltmax)
             dlt = 2*dlt
          end if
       else if( rho < 0.25 ) then
!          dlt = max(dlt/2,dltmin)
          dlt = dlt/2
       end if
       dlt = min(dlt,s-smin,smax-s)

!       call newstep(fun,t,s,fs,h22)


    end do


  end subroutine snewton

  subroutine newstep(fun,t,s,fs,h22)

    class(tsr1mean_function), intent(in), target :: fun
    real(REAL64), intent(in) :: t, s
    real(REAL64), intent(out) :: fs, h22

    real(REAL64), dimension(2) :: fvec
    real(REAL64), dimension(2,2) :: hess

    call fun%grad(t,s,fvec)
    call fun%hess(t,s,hess)

    fs = fvec(2)
    h22 = hess(2,2)

  end subroutine newstep

  elemental real(REAL64) function Phi(x,t,s)
    real(REAL64), intent(in) :: x,t,s
    Phi = (1 + erf((x-t)/s/sqrt(2.0)))/2
  end function Phi

!!$  subroutine trgscale1(fcn,t,xcdf,ycdf,beta,smin,smax,s,reli)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL64), dimension(:), intent(in) :: xcdf,ycdf
!!$    real(REAL64), intent(in) :: t,beta,smin,smax
!!$    real(REAL64), intent(in out) :: s
!!$    logical, intent(out) :: reli
!!$    real(REAL64) :: stol, frac, dir, dltmin, dltmax, sig, best_s, best_frac, best_r,&
!!$         r1,r2,rscale,snew
!!$    logical :: converge, ismin
!!$    integer :: n, iter, edge, minmax
!!$
!!$    reli = .false.
!!$    n = size(fcn%data)
!!$    sig = s / sqrt(real(n))
!!$    stol = sig / 10
!!$    dltmin = 2*stol
!!$    dltmax = 10*sig
!!$    best_s = -1
!!$    best_frac = 0
!!$    best_r = huge(best_r)
!!$
!!$    s = smin + stol
!!$
!!$    r1 = sum(abs(ycdf - 0.5)) - n / 2.0
!!$
!!$
!!$    write(*,*) 'sig,stol:',sig,stol,size(fcn%data)
!!$    write(*,*) 'smin,smax:',smin,smax,dltmin,dltmax
!!$
!!$!    kind = .true. ! minimum
!!$    minmax = 1 ! minimum
!!$!    call snewton(fcn,t,smin,smax,kind,stol,s,converge)
!!$!    frac = fcn%nfrac(t,s)
!!$!    write(*,*) 'frac:',s,frac
!!$
!!$    if( debug ) open(17,file='/tmp/trgscale')
!!$
!!$    do iter = 1, 30
!!$
!!$       snew = s
!!$       call snewton(fcn,t,smin,smax,dltmin,dltmax,minmax,stol,snew,ismin,edge,converge)
!!$       frac = fcn%nfrac(t,s)
!!$
!!$       if( converge ) s = snew
!!$
!!$       if( ismin .and. converge ) then
!!$          r2 = sum(abs(Phi(xcdf,t,s) - 0.5)) - n / 2.0
!!$          rscale = r1 - r2
!!$          write(*,*) 'r1,r2:',r1,r2
!!$
!!$       else
!!$          rscale = huge(rscale)
!!$       end if
!!$
!!$
!!$       write(*,*) 'snew:',s,rscale,ismin,edge,converge
!!$
!!$
!!$!       if( abs(frac - beta) < 1.0 / n .and. ismin .and. converge ) then
!!$       if( abs(rscale) < -1 .and. ismin .and. converge ) then
!!$          reli = .true.
!!$          exit
!!$       end if
!!$
!!$       if( abs(rscale) > abs(best_r) .and. ismin .and. converge ) then
!!$          reli = .true.
!!$          exit
!!$       end if
!!$
!!$
!!$       if( abs(frac - beta) < abs(best_frac - beta) .and. ismin .and. converge ) then
!!$          best_frac = frac
!!$!          best_s = s
!!$       end if
!!$
!!$       if( abs(rscale) < abs(best_r) ) then
!!$          best_r = rscale
!!$          best_s = s
!!$       end if
!!$
!!$       if( frac < beta ) then
!!$          dir = 1
!!$       else
!!$          dir = -1
!!$       end if
!!$
!!$       if( edge /= 0 ) then
!!$          if( best_s > 0 ) exit
!!$          dir = -edge
!!$       end if
!!$
!!$       s = s + dir*stol
!!$       !       kind = .not. kind
!!$       minmax = -minmax
!!$
!!$    end do
!!$
!!$    if( debug ) close(17)
!!$
!!$    if( best_s > 0 ) then
!!$       s = best_s
!!$    end if
!!$
!!$
!!$!    s = s + stol
!!$!    call snewton(fcn,t,smin,smax,.false.,stol,s,converge)
!!$!    s = s + stol
!!$!    call snewton(fcn,t,smin,smax,.true.,stol,s,converge)
!!$!    write(*,*) 'frac:',s,fcn%nfrac(t,s)
!!$
!!$  end subroutine trgscale1
!!$


!!$  subroutine fnewton(fun,t,smin,minmax,sig,stol,s,edge,converge)
!!$
!!$    class(tsr1mean_function), intent(in), target :: fun
!!$    real(REAL64), intent(in) :: t, stol, smin, sig
!!$    integer, intent(in) :: minmax
!!$    real(REAL64), intent(in out) :: s
!!$    logical, intent(out) :: edge, converge
!!$
!!$    real(REAL64), dimension(2,3*precision(s)) :: filter ! (h,f)
!!$    integer :: nf
!!$
!!$    nf = 0
!!$
!!$    do iter = 1, 3*precision(s)
!!$
!!$       do
!!$
!!$          ! step
!!$          call newstep(fun,t,s,fs,h22)
!!$          if( one*h22 > 0 ) then
!!$             d = -fs / h22
!!$          else
!!$             d = dlt
!!$          end if
!!$          if( abs(d) >= dlt ) d = -one*sign(dlt,fs)
!!$          sp = s + d
!!$
!!$          converge = abs(d) < stol .and. one*h22 > 0
!!$          if( converge .and. sp >= smin ) then
!!$             s = sp
!!$             exit
!!$          end if
!!$
!!$          if( sp < smin ) then
!!$             nf = nf + 1
!!$             filter(nf) = [h,f]
!!$             s = smin
!!$             edge = .true.
!!$             return
!!$          else
!!$             h = max(0,sp - smin)
!!$             f = fun%val(t,sp)
!!$             dq = -d*(fs + 0.5*h22*d)
!!$
!!$             if( dq > 0 ) then
!!$                check_dq = (filter(nf-1,2) - f) > 0.5*dq
!!$             else
!!$                check_dq = .false.
!!$             end if
!!$
!!$             beta = 1e-5
!!$             check_filter = f <= f(j) - beta*h(j) .or. h <= h(j)*(1 - beta)
!!$
!!$             if( check_dq .and. check_filter ) then
!!$                exit
!!$             else
!!$                dlt = dlt / 2
!!$             end if
!!$          end if
!!$
!!$       end do
!!$
!!$       if( converge ) exit
!!$
!!$       nf = nf + 1
!!$       filter(nf) = [h,f]
!!$       s = sp
!!$
!!$    end do
!!$
!!$  end subroutine fnewton


  !!$  subroutine trgscale3(fcn,t,xcdf,ycdf,smin,smax,s,reli)
!!$
!!$    class(tsr1mean_function), intent(in) :: fcn
!!$    real(REAL64), dimension(:), intent(in) :: xcdf,ycdf
!!$    real(REAL64), intent(in) :: t,smin,smax
!!$    real(REAL64), intent(in out) :: s
!!$    logical, intent(out) :: reli
!!$    real(REAL64) :: stol, dltmin, dltmax, sig, best_r, r1,r2,rscore,ss,rf
!!$    logical :: converge, edge, best_found, valid
!!$    integer :: n, iter, minmax, dir
!!$
!!$    reli = .false.
!!$    n = size(fcn%data)
!!$    best_r = huge(best_r)
!!$    best_found = .false.
!!$    rf = 1 / sqrt(real(2*n))
!!$    sig = s / sqrt(real(n))
!!$    stol = sig / 10
!!$    dltmin = stol
!!$    dltmax = sig
!!$    write(*,*) 'sig,stol:',sig,stol,size(fcn%data)
!!$
!!$    r1 = sum(abs(ycdf - 0.5)) - n / 2.0
!!$
!!$    write(*,*) 'smin,smax:',smin,smax,dltmin,dltmax
!!$
!!$    minmax = 1 ! minimum
!!$
!!$    if( debug ) open(17,file='/tmp/trgscale')
!!$
!!$    ss = s
!!$    do iter = 1, 13
!!$
!!$       call snewton(fcn,t,smin,smax,dltmin,dltmax,minmax,stol,r1,rf,best_r,xcdf,ss,edge,converge)
!!$!       call fnewton(fcn,t,smin,smax,minmax,sig,stol,ss,edge,converge)
!!$
!!$       if( debug ) write(*,*) 'snewton:',ss,minmax,edge,converge
!!$
!!$       if( converge ) then
!!$
!!$          if( minmax > 0 ) then
!!$             r2 = sum(abs(Phi(xcdf,t,ss) - 0.5)) - n / 2.0
!!$             rscore = (r1 - r2) / rf
!!$             rscore = fcn%stderr(t,ss,valid)
!!$             if( .not. valid ) rscore = huge(rscore)
!!$!             write(*,*) 'r1,r2:',r1,r2
!!$
!!$             if( debug ) write(*,*) 'score:',rscore,best_r
!!$             if( abs(rscore) < abs(best_r) ) then
!!$                best_found = .true.
!!$                best_r = rscore
!!$                s = ss
!!$             end if
!!$
!!$             if( iter == 1 ) then
!!$                if( rscore < 0 ) then
!!$                   dir = 1
!!$                else
!!$                   dir = -1
!!$                end if
!!$             end if
!!$
!!$             ! if |rscale| < 1
!!$!             if( abs(rscore) < -1 .or. abs(rscore) > abs(best_r) ) then
!!$!                reli = .true.
!!$!                exit
!!$!             end if
!!$
!!$          else
!!$             rscore = huge(rscore)
!!$          end if
!!$
!!$!          if( edge ) then
!!$!             reli = best_found
!!$!             exit
!!$!          end if
!!$
!!$       else
!!$
!!$          exit ! we're in unexpected deep doubts.
!!$
!!$       end if
!!$
!!$!       write(*,*) iter,ss,dir,stol
!!$       !ss = ss - sig
!!$       minmax = -minmax
!!$       ss = ss + dir * stol
!!$
!!$       r2 = sum(abs(Phi(xcdf,t,ss) - 0.5)) - n / 2.0
!!$!       if( abs(r1 - r2)/rf > abs(best_r) ) exit
!!$
!!$       if( ss < smin ) then
!!$          reli = best_found
!!$          exit
!!$       end if
!!$
!!$    end do
!!$
!!$    if( debug ) close(17)
!!$
!!$
!!$  end subroutine trgscale3




end module trust_global_scale
