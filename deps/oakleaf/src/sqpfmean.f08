!
!  An implementation of Filter SQP method for the robust mean estimation.
!
!  Copyright © 2023-2026 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Oakleaf.
!
!  Oakleaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Oakleaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!
!  This implementation is due Filter SQP method by Fletcher, Leyffer (2003),
!  Fletcher, Leyffer & Toint (2006).
!  And inspired by all amtter of the second edition of Numerical Optimization
!  by Jorge Nocedal and Stephen J. Wright, modified and specialised
!  Algorithm 18.1 (Local SQP) combined with the simplified TR step
!  by More & Sorenson (1981).
!

module sqpfmean_module_REAL64

  use rlike_REAL64
  use fminim_REAL64
  use iso_fortran_env

  implicit none
  private

  public :: sqpfmean_REAL64

  type, extends(goldfinger_function_REAL64), private :: gold_fun
     class(robust_likelihood_REAL64), pointer :: fcn
     integer :: pm
     real(REAL64), dimension(2) :: q,p
     integer :: nfev
   contains
     procedure :: f => goldfcn
  end type gold_fun

contains

  subroutine sqpfmean_REAL64(fcn,t,s,B,rtol,verbose,report_file,iflag)

    use qtrstep_REAL64

    class(robust_likelihood_REAL64), pointer, intent(in) :: fcn
    real(REAL64), intent(in out) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: B
    real(REAL64), intent(in) :: rtol
    logical, intent(in) :: verbose
    character(len=*), intent(in) :: report_file
    integer, intent(out) :: iflag

    integer, parameter :: maxiter = 5*precision(t)
    real(REAL64), parameter :: macheps = 42*epsilon(1.0_REAL64)
    real(REAL64), parameter :: gold = (1 + sqrt(5.0_REAL64)) / 2
    real(REAL64), parameter :: gold2 = 1 - 1 / gold
    real(REAL64), parameter :: gold7 = 1 / gold**7
    real(REAL64), parameter :: gold5 = 1 / gold**5
    real(REAL64), parameter :: Kg = 1.4404 ! Kg = 1/log2(gold)
    real(REAL64), parameter :: log2 = 0.6931 ! = log(2)
    real(REAL64), parameter :: trtol = 1.0 / 13.0

    integer :: iter, it, trinfo, N, rstflag, nf, maxit, nfev, nrst, nhrd, nout, &
         nidf, niter, runit, iostat
    real(REAL64), dimension(2) :: p, fgrad, fg, y
    real(REAL64) :: fval, fv, dq, dlt, tp, sp, pnorm, xtol, rho, fmax, &
         Ravg, df, ftol, diag, feps, sqrtN
    logical :: feasible, posdef, report, quasi, inside
    character(len=3), dimension(2) :: ch
    character(len=*), parameter :: fmt = &
         "(i0,a,2f9.3,a,f8.3,es9.2,l2,a,2es8.1,2f7.2,i2,l2)"
    character(len=80) :: fmts

    if( .not. (size(fcn%data) > 1 .and. s > tiny(s) .and. rtol > macheps) ) then
       iflag = 6
       if( verbose ) write(error_unit,*) 'sqpfmean(): assumptions unsatisfied.'
       return
    end if

    iflag = 6
    N = size(fcn%data)
    nfev = 0
    nrst = 0
    nhrd = 0
    nout = 0
    nidf = 0
    niter = 0
    sqrtN = sqrt(real(N,REAL64))
    fmax = fcn%fmax()
    ftol = sqrtN*macheps
    maxit = int(2*Kg*(log(sqrtN/max(rtol,macheps))/log2)**2) + 1
    diag = 1 / s
    pnorm = 0
    dlt = 0
    xtol = rtol * s / sqrtN
    call eye(B)
    posdef = .true.
    quasi = .true.

    report = report_file /= ''
    if( report ) then
       open(newunit=runit,file=report_file,iostat=iostat)
       if( iostat /= 0 ) then
          if( verbose ) write(error_unit,*) &
               'sqpfmean(): failed to open the file `',trim(report_file),"'."
          iflag = 6
          close(runit)
          return
       end if
       write(runit,'(a)') '# sqpfmean()'
       write(runit,'(a,g0)') '# t0 = ',t
       write(runit,'(a,g0)') '# s0 = ',s
       write(runit,'(a)') '# iter t s fval dlt p fgrad B flag'
    end if

    call fcn%sqpf(t,s,fval,Ravg,fgrad)
    fgrad = fgrad / diag
    nfev = 3

    if( verbose ) then
       write(error_unit,'(a,i0,2(1x,es10.2),1x,f7.3,1x,2es10.3)') &
            'sqpfmean: N,xtol,ftol,fmax,1/diag,macheps: ',N,xtol,ftol,fmax,1/diag,&
            macheps
       write(error_unit,fmt) 0,':',t,s," >",fval,norm2(fgrad)*diag,posdef," |",&
            dlt/diag,pnorm/diag,Ravg
    end if

    do iter = 1, maxiter

       ! init trust-region
       dlt = s / gold * diag
       feps = (1 + abs(fval))*macheps

       ! find optimal delta
       do it = 1, maxit
          niter = niter + 1

          ! the trust-region step by the eigenvalue factorisation
          call trstep_REAL64(B,fgrad,dlt,trtol,p,pnorm,posdef,trinfo,verbose)
          inside = posdef .and. trinfo == 1
          if( verbose ) then
             select case(trinfo)
             case(3)
                nout = nout + 1
             case(4)
                nidf = nidf + 1
             case(5)
                nhrd = nhrd + 1
             end select
          end if

          if( trinfo == 6 ) then
             ! no convergence in TR step parameter update
             if( verbose ) error stop 'sqpfmean(): no convergence in TR step'
             goto 666
          end if

          tp = t + p(1) / diag
          sp = s + p(2) / diag
          Ravg = -1

          feasible = sp > macheps*max(abs(t),1.0_REAL64)
          if( feasible ) then
             call fcn%sqpf(tp,sp,fv,Ravg,fg)
             nfev = nfev + 3
             fg = fg / diag
             feasible = Ravg < 0.9*fmax ! fmax() can be huge()
          end if

          if( verbose ) then
             write(error_unit,fmt) &
                  it,'.',tp,sp," > ",fv,norm2(fg)*diag,posdef," |",&
                  dlt/diag,pnorm/diag,Ravg,norm2(p)/dlt,trinfo,feasible
          end if

          if( feasible ) then

             dq = -dot_product(fgrad + matmul(B,p)/2,p)
             df = fval - fv
             y = fg - fgrad

             if( verbose ) then
                rho = 99.989
                if( abs(dq) > 0 ) rho = max(min(df / dq, rho),-rho)
                write(error_unit,'(2(a,es9.2),a,es9.2,2(a,f7.3),a,f7.2)') &
                     '  dq=',dq,' df=',df, '  ftol=',ftol,'  fv=',fv, &
                     ' Ravg=',Ravg,' rho=',rho
             end if

             if( (dq > 0 .and. gold2*dq < df .and. df < gold*dq) .or. &
                  (abs(df) < ftol) ) exit

             dlt = dlt / gold
             call sr1update(p,y,pnorm,B)
             ! B is updated to ensure convergence, no matter if quasi == .false.

          else !if( .not. feasible ) then
             call restore(fcn,tp,sp,trtol,dlt/diag,fmax,verbose,nf,rstflag)
             nfev = nfev + nf
             nrst = nrst + 1

             if(verbose) write(error_unit,'(a,2(1x,g0),a,i1)') &
                  'State restored to ',tp,sp,' with iflag= ',rstflag
             if( rstflag /= 1 ) then
                iflag = rstflag
                close(runit)
                return
             end if

             call fcn%sqpf(tp,sp,fv,Ravg,fg)
             nfev = nfev + 3
             fg = fg / diag
             p = 0
             pnorm = dlt
             quasi = .true.
             call eye(B)
             inside = .false.
             trinfo = -1
             exit
          end if

       end do

       if( report ) write(runit,*) &
            iter,t,s,fval,dlt/diag,[tp-t,sp-s],fgrad*diag,B*diag**2,trinfo

       t = tp
       s = sp
       fval = fv
       fgrad = fg
       xtol = rtol * s / sqrtN
       !       quasi = .not. (inside .and. sum(abs(fgrad))/sum(abs(B)) < gold7)
!       write(*,*) abs(fgrad)/s,[B(1,2),B(2,2)],gold5,abs(B(:,2))+macheps
       quasi = .not. (inside .and. all(abs(fgrad)/s < gold5*abs(B(:,2))+macheps))

       if( verbose ) write(error_unit,fmt) &
            iter,': ',t,s," > ",fval,norm2(fgrad)*diag,posdef," |",&
            dlt/diag,pnorm/diag,Ravg,rho,trinfo,quasi

       if( inside ) then
          if( pnorm / diag < xtol ) then
             iflag = 0
             exit
          else if( -feps < df .and. df < ftol ) then
             iflag = 1
             exit
          end if
       end if

       if( quasi ) then
          call sr1update(p,y,pnorm,B)
       else
          call fcn%minhess(t,s,B)
          B = B / diag**2
          nfev = nfev + 3
       end if

    end do

    B = B * diag**2

666 continue

    if( verbose ) then
       if( 0 <= iflag .and. iflag <= 1 ) then
          ch = ''
          ch(iflag+1) = '(*)'
          write(error_unit,'(4(a,es10.2))') &
               'Finish: xtol'//trim(ch(1))//':',pnorm/diag,' < ',xtol, &
               ',   ftol'//trim(ch(2))//':',df,' < ',ftol
       else
          write(error_unit,'(4(a,es10.2))') &
               'Finish fail(!): xtol:',pnorm/diag,' < ',xtol,',  ftol:',df,' < ',ftol
       end if
       write(error_unit,'(a,2(g0,2x),a,i0)') 'sqpfmean: ',t,s,' iflag=',iflag
       fmts = '(a,2(es12.5,2x),a,2(f12.5,2x))'
       write(error_unit,fmts) 'B(1,:)= ',B(1,:),', B(1,:)/s**2= ',B(1,:)/diag**2
       write(error_unit,fmts) 'B(2,:)= ',B(2,:),', B(2,:)/s**2= ',B(2,:)/diag**2
       write(error_unit,'(7(a,i0))') 'niter=',niter, &
            ' nfev=',nfev,' nrst=',nrst,' nout=',nout,' nidf=',nidf,' nhrd=',nhrd,&
            ' flag=',iflag
    end if

    if( report ) then
       write(runit,'(a,g0)') '# tfin = ',t
       write(runit,'(a,g0)') '# sfin = ',s
       write(runit,'(a,g0)') '# niter = ',niter
       write(runit,'(a,g0)') '# nfev = ',nfev
       write(runit,'(a,g0)') '# nrst = ',nrst
       write(runit,'(a,g0)') '# nout = ',nout
       write(runit,'(a,g0)') '# nidf = ',nidf
       write(runit,'(a,g0)') '# nhrd = ',nhrd
       write(runit,'(a,4(g0,1x))') '# happrox = ',B
       block
         real(REAL64), dimension(2,2) :: hess
         call fcn%hess(t,s,hess)
         write(runit,'(a,4(g0,1x))') '# hess = ',hess
       end block
       close(runit)
    end if

  end subroutine sqpfmean_REAL64

  subroutine restore(fcn,t,s,trtol,dlt,fmax,verbose,nfev,lflag)

    class(robust_likelihood_REAL64), pointer, intent(in) :: fcn
    real(REAL64), intent(in) :: trtol,dlt,fmax
    real(REAL64), intent(in out) :: t
    real(REAL64), intent(out) :: s
    logical, intent(in) :: verbose
    integer, intent(out) :: nfev, lflag

    type(gold_fun) :: goldfun
    integer, parameter :: maxiter = 666
    real(REAL64), parameter :: macheps = 42*epsilon(1.0_REAL64)
    real(REAL64), dimension(2), parameter :: p = [0,1]
    real(REAL64), dimension(2) :: q
    real(REAL64) :: tol, alpha, smin, smax, a, b, r
    integer :: iter, iflag
    character(len=*), parameter :: fmt = &
         '(a,i0,": ",2(f10.5,1x),f7.3," |",2es10.2,i3)'

    lflag = 6
    tol = trtol*dlt
    a = 0
    b = dlt
    smin = macheps*(1 + abs(t))
    smax = min(maxval(abs(fcn%data - t)),sqrt(huge(smax)))
    smax = maxval(abs(fcn%data - t))
    q = [t,smin]
    nfev = 0

    if( verbose ) write(error_unit,'(a,2g12.5,3(a,es10.2))') &
         'Restore start point: ',q,' tol: ',tol,' dlt: ',dlt,' smax:',smax

    do iter = 1, maxiter

       if( verbose ) then
          write(error_unit,'(2(a,2(f10.5,1x)))') 'line search: x=',q,' to ',q+b*p
          write(error_unit,'(a,5f10.5)') '         ...',p,t,s,b
       end if

       goldfun = gold_fun(fcn,-1,q,p,0)
       call goldfinger_REAL64(goldfun,a,b,tol,alpha,iflag)
       nfev = nfev + goldfun%nfev
       if( .not. (1 <= iflag .and. iflag <= 3) ) return

       q = q + min(alpha,b-2.01*tol)*p

       if( verbose ) write(error_unit,fmt) &
            "Gold max ",iter,q,fcn%func(q(1),q(2)),b,alpha,iflag

       if( iflag == 1 .or. iflag == 2 ) exit
       if( iter == maxiter ) return

       if( q(2) > smax ) then
          if( abs(fcn%func(q(1),q(2)) - log(q(2))) < macheps ) then
             lflag = 3
             return
          end if
       end if
       b = 2*q(2)
    end do

    tol = trtol*q(2)
    b = 2*q(2)
    do iter = 1, maxiter

       ! skip some wrinkles
       goldfun = gold_fun(fcn,1,q,p,0)
       call goldfinger_REAL64(goldfun,a,b,tol,alpha,iflag)
       nfev = nfev + goldfun%nfev
       if( .not. (1 <= iflag .and. iflag <= 3) ) return

       q = q + min(alpha,b-2.01*tol)*p
       r = fcn%func(q(1),q(2)) - log(q(2))

       if( verbose ) write(error_unit,fmt) &
            "Gold min ",iter,q,fcn%func(q(1),q(2)),b,alpha,iflag

       if( (iflag == 1 .or. iflag == 2) .and. r < fmax) then
          t = q(1)
          s = q(2)
          lflag = 1
          exit
       end if

       b = 2*q(2)
    end do

  end subroutine restore

  real(REAL64) function goldfcn(this,x)

    class(gold_fun) :: this
    real(REAL64), intent(in) :: x
    real(REAL64), dimension(2) :: r

    associate( q => this%q, p => this%p, pm => this%pm, nfev => this%nfev )

      r = q + x*p
      goldfcn = pm * this%fcn%func(r(1),r(2))
      nfev = nfev + 1

    end associate

  end function goldfcn


  pure subroutine sr1update(p,y,pnorm,B)

    real(REAL64), dimension(:), intent(in) :: p,y
    real(REAL64), intent(in) :: pnorm
    real(REAL64), dimension(:,:), intent(in out) :: B

    real(REAL64), parameter :: beta = 1e-8  ! 0 < r < 1
    real(REAL64), parameter :: rmin = 10*epsilon(p)

    real(REAL64), dimension(size(p)) :: r
    real(REAL64) :: rnorm, rp
    integer :: n,i,j

    n = size(p)
    r = y - matmul(B,p)
    rnorm = norm2(r)
    rp = dot_product(r,p)

!    write(error_unit,*) '...',rs,rnorm,snorm

    if( abs(rp) > beta*rnorm*pnorm .and. min(rnorm,pnorm) > rmin ) then

       forall( i = 1:n, j = 1:n ) B(i,j) = B(i,j) + (r(i)*r(j))/rp

       ! Note, there is an alternative way how to write it, perhaps,
       ! a bit quicker. Unfortunately, it has produced non-symmetric
       ! matrices due cancellation of digits. Do not use it.
       ! w = r / rs
       ! forall( i = 1:n ) B(i,:) = B(i,:) + r(i)*w
    end if

  end subroutine sr1update

  pure subroutine eye(B)

    real(REAL64), dimension(:,:), intent(out) :: B

    B(1,:) = [1, 0]
    B(2,:) = [0, 1]

  end subroutine eye


end module sqpfmean_module_REAL64
