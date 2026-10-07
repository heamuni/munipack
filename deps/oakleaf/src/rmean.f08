!
!  Robust mean
!
!  Copyright © 2024-2026 F.Hroch (hroch@physics.muni.cz)
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
!

module robustmean_REAL64

  ! This module provides rmean subroutine for estimation of robust mean.
  !
  ! rmean should be called as:
  !
  !  call rmean(data,mean,stderr,stdsig,scale,reltol,robfun,sigic,
  !                                   rankinit,flag,verbose,report)
  !
  ! on input:
  !   data - array of data values to be estimated
  !   verbose (optional) - print additional informations
  !
  ! on output are estimated:
  !   mean - robust mean
  !   stderr (optional) - standard error
  !   stdsig (optional) - standard deviation
  !   scale (optional) - scale of function normalisation
  !   reltol (optional) - relative tolerance
  !   robfun (optional) - select robust function ('hampel', 'tukey', 'huber')
  !                       or non-robust ('square')
  !   sigic (optional) - .true. to calculate stdsig by influence function
  !   rankinit (optional) - initial estimate by median, MAD (if .false.,
  !                         mean and scale shlould be provided
  !   flag (optional) - indicates reliability of results:
  !          0 == success, results are both reliable and accurate
  !          1 == success, convergence in functional values
  !          2 == basic estimates, few datapoints available
  !          3 == data looks nearly identical
  !          5 == mean or scale unspecified if rankinit == .false.
  !          6 == a general error (try verbose)
  !          8 == failed to allocate memory
  !   verbose (optinal) - verbose prints
  !   report (optional) - a file name to save detailed informations, steps, etc.
  !
  ! The results gives the robust estimate of the true value X
  ! of the sample x (data), with 68% probability, into the interval
  !
  !         mean - stderr  <  X  <  mean + stderr
  !
  ! The distribution is estimated as the Normal distribution with the parameters
  !
  !          N(mean,stdsig).
  !
  ! Robust estimators has been prepared on base of
  !  * Hogg in Launer, Wilkinson: Robustness in Statistics
  !  * Huber: Robust Statistics
  !  * my experiences

  use iso_fortran_env

  implicit none
  private

  public :: rmean_REAL64

contains

  subroutine rmean_REAL64(data,mean,stderr,stdsig,scale,reltol,robfun,sigic, &
       rankinit,flag,verbose,report)

    ! This routine finds the robust mean by searching the maximum of likelihood.

    use sqpfmean_module_REAL64
    use rlike_REAL64
    use quicksort
    use quickmedian

    real(REAL64), dimension(:), target, intent(in) :: data
    real(REAL64), intent(in out) :: mean
    real(REAL64), intent(out), optional :: stderr, stdsig
    real(REAL64), intent(in out), optional :: scale
    real(REAL64), intent(in), optional :: reltol
    character(len=*), intent(in), optional :: robfun
    logical, intent(in), optional :: sigic,rankinit
    integer, intent(out), optional :: flag
    logical, intent(in), optional :: verbose
    character(len=*), intent(in), optional :: report

    type(robust_likelihood_hampel), target :: hampelfun
    type(robust_likelihood_tukey), target :: tukeyfun
    type(robust_likelihood_huber), target :: huberfun
    type(nonrobust_likelihood_square), target :: squarefun
    class(robust_likelihood_REAL64), pointer :: fcn

    integer, parameter :: N1E6 = 2**20
    real(REAL64), parameter :: q34 = 0.6745
    real(REAL64), parameter :: gold = (1 + sqrt(5.0_REAL64)) / 2
    real(REAL64), dimension(:), allocatable :: xcdf
    real(REAL64), dimension(2,2) :: B
    real(REAL64) :: t,dt,s,s0,t0,sig,rtol,cov11,det,Ravg,favg,fmin,fmax,smin
    integer :: N, NK, NN, K, i, aflag, info
    logical :: verb, rinit, converge, influence_curve
    character(len=666) :: modifier, errmsg
    character(len=:), allocatable :: report_file

    verb = .false.
    rinit = .true.
    influence_curve = .false.
    rtol = 1.0 / 42.0
    N = size(data)
    report_file = ''
    if( present(verbose) ) verb = verbose
    if( present(flag) ) flag = 6
    if( present(reltol) ) rtol = reltol
    if( present(sigic) ) influence_curve = sigic
    if( present(report) ) report_file = report

    if( present(rankinit) ) then
       rinit = rankinit
       if( .not. rinit ) then
          if( .not. present(scale) ) then
             if( verb ) write(error_unit,*) &
                  "rmean(): specify mean & scale, if rankinit == .false."
             if( present(flag) ) flag = 6
             return
          end if
       end if
    end if

    if( present(robfun) ) then
       modifier = 'robust'
       if( robfun == 'huber' ) then
          fcn => huberfun
       else if( robfun == 'tukey' ) then
          fcn => tukeyfun
       else if( robfun == 'square' ) then
          fcn => squarefun
          modifier = 'NON-robust'
       else
          fcn => hampelfun
       end if
       if( verb ) write(error_unit,*) &
            'The selected ',trim(modifier),' function: ',trim(robfun)
    else
       fcn => hampelfun
       if( verb ) write(error_unit,*) &
            'The robust function, by default: ',"Hampel's"
    end if
    fcn%data => data


    if( n == 0 ) then

       ! results are undefined
       if( present(flag) ) flag = 6
       if( verb ) write(error_unit,* ) &
            'rmean: finish for no points (all parameters are undefined)'
       return

    else if( n == 1 ) then

       ! dispersion and std.err. are defind as zero
       mean = data(1)
       if( present(stderr) ) stderr = 0
       if( present(stdsig) ) stdsig = 0
       if( present(scale) )   scale = 0
       if( present(flag) )     flag = 2
       if( verb ) write(error_unit,*) 'rmean: ',mean,' stderr:',0,' stdsig:',0, &
            ' scale:',1,' (finish for a single point)'
       return

    else if( n == 2 ) then

       mean = sum(data) / N
       s = norm2(data-mean)
       sig = s
       if( present(stderr) ) stderr = sig / sqrt(real(N,REAL64))
       if( present(stdsig) ) stdsig = sig
       if( present(scale) )   scale = s ! incorrect
       if( present(flag) )     flag = 2
       if( verb ) write(error_unit,*) 'rmean: ',mean,' stderr:',s/sqrt(2.0), &
            ' stdsig:',s,' scale:',s,' (finish for a pair)'
       return

    end if

    !
    ! ... n > 2, below this point ...
    !
    ! values are initialised by quantiles (by median)
    if( rinit ) then

       K = max(N / N1E6, 1)
       NN = N / K
       NK = K*NN
       allocate(xcdf(NN),stat=aflag,errmsg=errmsg)
       if( aflag /= 0 ) then
          if( verb ) write(error_unit,*) 'rmean() allocate: ',trim(errmsg)
          if( present(flag) ) flag = 8
          return
       end if
       if( N < 13 ) then
          xcdf = data
          call qsort(xcdf)
          t0 = (xcdf(N/2+1) + xcdf((N+1)/2)) / 2
          xcdf = abs(data - t0)
          call qsort(xcdf)
          smin = xcdf(N/2) / q34
          fmax = fcn%fmax()
          fmin = fmax + 1
          do i = N/2, N
             s0 = xcdf(i) / q34
             if( s0 > 0 ) then
                favg = fcn%func(t0,s0)
                Ravg = favg - log(s0)
                if( Ravg < fmax / gold .and. favg > fmin ) then
                   s0 = smin
                   exit
                end if
                smin = s0
                fmin = favg
             end if
          end do
       else
          xcdf = data(1:NK:K)
          t0 = qmed(xcdf,NN/2)
          xcdf = abs(data(1:NK:K) - t0)
          s0 = qmed(xcdf,NN/2) / q34
       end if
       if( verb ) write(error_unit,'(a,g0,a,es12.5,3a)') &
            "rob.init: t,s,dt=",t0," (median),",s0," (MAD)"
       deallocate(xcdf)
    else
       t0 = mean
       s0 = scale
       if( verb ) write(error_unit,*) "rob.init by user: t,s=",t0,s0
    end if

    if( abs(s0) < 42*epsilon(s0) ) then
       ! If more than half of data has (nearly) identical values,
       ! the scale estimate is (close to) zero. No application of
       ! robust methods is possible violating s0 > 0.
       !
       ! All the data are supposed to be identical (minority of points
       ! can still deviate), result can be considered as reliable.

       mean = t0
       if( present(stderr) ) stderr = 0
       if( present(stdsig) ) stdsig = 0
       if( present(scale) )  scale = s0
       if( present(flag) )   flag = 3

       if( verb ) write(error_unit,*) 'rmean: ',t0,' stderr:',0,' stdsig:',0, &
            ' scale:',s0,' (finish for identical points)'
       return
    end if

    t = t0
    s = s0
    call sqpfmean_REAL64(fcn,t,s,B,rtol,verb,report_file,info)
    converge = info == 0 .or. info == 1

    if( converge ) then
       if( influence_curve ) then
          dt = fcn%stderr(t,s)
          sig = dt*sqrt(real(N))
       else
          det = B(1,1)*B(2,2) - B(1,2)*B(2,1)
          if( det > 0 .and. B(2,2) > 0 ) then
             cov11 = B(2,2) / det
             sig = sqrt(N/(N - 1.0_REAL64) * cov11)
             dt = sig / sqrt(real(N,REAL64))
          else
             error stop 'rmean: the assumption of pos. def. Hessian unsatisfied.'
             ! this is commonly due SR1 update in vicinity of minimum
          end if
       end if
    else
       t = t0
       s = s0
       dt = s0 / sqrt(real(N,REAL64))
       sig = s
    end if

    if( verb .and. converge ) then
       block
         use rfun

         real(REAL64), dimension(2,2) :: hess
         character(len=80) :: fmt
         integer :: nrob,nplus,i
         integer, dimension(3) :: nsig
         real(REAL64), dimension(:), allocatable :: res
         real(REAL64), parameter :: toleps = 42*epsilon(1.0_REAL64)

         call fcn%hess(t,s,hess)
         fmt = '(a,4es13.5,3x,a,es8.2e2)'
         write(error_unit,fmt) 'B:',B
         write(error_unit,fmt) 'H:',hess,"|B-H|/|H|: ",norm2(B-hess)/norm2(hess)

         fmt = '(a,g0,3(a,es9.2),1x,2l1)'
         write(error_unit,fmt) &
              'rmean: ',t,' stderr:',dt,' stdsig:',sig,' scale:',s,converge
         block
           real(REAL64) :: a,stdsig,stderr

           stderr = fcn%stderr(t,s)
           stdsig = stderr*sqrt(1.0*N)
           write(error_unit,'(a,1x,2(a,es9.2))') "       by influence curve", &
                ' stderr:',stderr,' stdsig:',stdsig

           a = sum(fcn%data)/n
           stdsig = sqrt(sum((fcn%data - a)**2)/(n-1.0))
           stderr = stdsig / sqrt(1.0*n)
           write(error_unit,fmt) 'amean: ',a,' stderr:',stderr,' stdsig:',stdsig

           write(error_unit,*) '1/I / var = ',(1/hess(1,1))/stdsig**2

         end block

         block
           real(REAL64) :: nu, varpi

           nu = hess(1,1)*s**2
           varpi = hess(2,2)*s**2-1
           write(error_unit,'(3(a,f0.5))') &
                'Estimate of nu = ',nu,'  sqrt(nu)= ',sqrt(max(nu,epsilon(nu))),&
                " and varpi= ",varpi
         end block

         ! sign test
         allocate(res,source=data)
         res = (data - t)/sig
         nplus = count(res > 0)
         write(error_unit,'(a,i0,"/",f0.1,", ",i0,"+-",f0.1)') &
              'Sign test (total/expected, +): ', &
              n,n/2.0,nplus, sqrt(nplus*0.25 + 1e-6)
         if( robfun == 'huber' ) then
            res = huber((data - t)/s)
         else if( robfun == 'tukey' ) then
            res = tukey((data - t)/s)
         else if( robfun == 'square' ) then
            res = (data - t)/s
         else
            res = pshampel((data - t)/s)
         end if
         nrob = count(abs(res) > epsilon(res))
         nplus = count(res > 0 .and. abs(res) > epsilon(res))
         write(error_unit,'(a,i0,"/",f0.1,", ",i0,"+-",f0.1)') &
              'Sign test with rob.function (total/expected, +): ', &
              nrob,nrob/2.0,nplus,sqrt(nplus*0.25 + 1e-6)
         res = (data - t)/sig
         do i = 1,3
            nsig(i) = count(abs(res) < i)
         end do
         write(error_unit,'(a,3(i0,"-sig: ",f0.1,"%":",  "))') &
              'Points within the interval:  ', (i,100.0*(nsig(i)/real(n)),i=1,3)
         write(error_unit,'(a,f0.1,"%")') &
              'Points out of the interval 5-sig: ',100*(count(abs(res) > 5)/real(n))
         if( robfun == 'huber' ) then
            res = dhuber((data - t)/s)
         else if( robfun == 'tukey' ) then
            res = dtukey((data - t)/s)
         else if( robfun == 'square' ) then
            res = 1
         else
            res = dpshampel((data - t)/s)
         end if
         write(error_unit,'(a,6(f0.3,a))') "Parts of psi'(x): ",&
              sum(res,res > toleps)/n,'/',real(count(res > toleps))/n,' (+)  ',&
              sum(res,res <-toleps)/n,'/',real(count(res <-toleps))/n,' (-)  _/',&
              real(count(abs(res) <= toleps),REAL64)/n,' (0) [sum/fraction]'
       end block
    end if

    mean = t
    if( present(stderr) ) stderr = dt
    if( present(stdsig) ) stdsig = sig
    if( present(scale) ) scale = s
    if( present(flag) ) flag = info

  end subroutine rmean_REAL64

end module robustmean_REAL64
