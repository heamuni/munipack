
program rmean_newton

  use toolbox
  use iso_fortran_env
  use robnewton
  implicit none

  integer, parameter :: nsamples = 999
  integer, dimension(5), parameter :: ndatas = [ 5, 10, 20, 100, 10000 ]
  real, dimension(2), parameter :: epsc = [ 0.0, 0.2 ]
  integer, dimension(precision(real(1.0,REAL64))) :: riter, triter, ziter
  integer, dimension(3,precision(real(1.0,REAL64))) :: iters
  integer :: i, j, l, rfail, trfail, zfail, nfail, ndata, nbad, zeroiter, trstep
  real(REAL64) :: zmean_scale
  class(tsr1mean_function), pointer :: zmeanfun
  character(len=80) :: filename

  open(1,file='/tmp/rmean_newton.tex')
  write(1,*) '\begin{tabular}{r|ccccc}'
  write(1,*) '{} ',(' & $',ndatas(i),'$',i=1,size(ndatas)),'\\ \hline'

  open(2,file='/tmp/rmean_newton_fails.tex')
  write(2,*) '\begin{tabular}{r|ccccc}'
  write(2,*) '$N=$ ',(' & $',ndatas(i),'$',i=1,size(ndatas)),'\\ \hline'

  do j = 1, size(epsc)
     write(1,'(a,i0,a)') '$\varepsilon=',nint(100*epsc(j)),'\%$'
     write(2,'(a,i0,a)') '$\varepsilon=',nint(100*epsc(j)),'\%$'
     do i = 1, size(ndatas)
        ndata = ndatas(i)
        nbad = nint(epsc(j)*ndata)
        write(*,*) '***',ndata,epsc(j)
        write(1,*) ' & '
        call one_sample(ndata,nbad)
        iters(1,:) = triter
        iters(2,:) = riter
        iters(3,:) = ziter
        write(filename,'(a,i2.2,a,i0)') &
             '/tmp/rmean_newton_iter_',nint(100*epsc(j)),'_',ndata
        open(3,file=filename)
        do l = 1, size(riter)
           write(3,*) l,nint(100.0*iters(:,l)/nsamples)
        end do
        close(3)
!        do l = 1, size(riter)
!           iters(i,3*(l-1)+1) = triter(l)
!           iters(i,3*(l-1)+2) = riter(l)
!           iters(i,3*(l-1)+3) = ziter(l)
!        end do
!        if( .not.(j == size(epsc) .and. i == size(ndatas)) ) then
!           write(1,*) '\\'
!        end if
     end do
     if( .not. (j == size(epsc)) ) then
        write(1,*) '\\'
        write(2,*) '\\'
     end if
  end do

  write(2,*) '\end{tabular}'
  write(1,*) '\end{tabular}'
  close(2)
  close(1)

contains

  subroutine one_sample(ndata,nbad)

    integer, intent(in) :: ndata, nbad
    real(REAL64), dimension(:), allocatable :: t,stderr,stdsig
    real(REAL64) :: f1,f2
    character(len=4) :: eol
    integer :: i

    allocate(t(nsamples),stderr(nsamples),stdsig(nsamples))

    triter = 0
    riter = 0
    ziter = 0
    trstep = 0
    rfail = 0
    trfail = 0
    nfail = 0
    zfail = 0
    do i = 1, nsamples
       call one_round(.false.,ndata,nbad,t(i),stderr(i),stdsig(i))
    end do

    write(*,*) 'Convergence of rnewton failed in ',rfail,' case(s).'
    write(*,*) 'Convergence of trnewton failed in ',trfail,' case(s).'
    write(*,*) 'Convergence of zeroin failed in ',zfail,' case(s).'
    write(*,*) 'Results inconsistent in ',nfail,' case(s).'
    write(*,*) 'Trsteps: ',trstep

    write(1,*) '\begin{tabular}{rcc}'
    do i = 1, size(riter)
       write(*,*) i, triter(i),riter(i), ziter(i)
!       write(3,*) i, nint(100.0*triter(i)/nsamples),nint(100.0*riter(i)/nsamples), nint(100.0*ziter(i)/nsamples)
       f1 = real(riter(i))/real(nsamples)
       f2 = real(triter(i))/real(nsamples)
       if( i < 7 ) then
          eol = ' \\ '
       else
          eol = ''
       end if
       write(1,'(a,3(i0,a),a)') &
            '$',i,'$ & $', nint(100*f1),'$ & $', nint(100*f2),'$',eol
    end do
    write(1,*) '\end{tabular}'

    write(2,'(a,2(i0,a))') &
         ' & $',nint(100.0*rfail/nsamples),'/',nint(100.0*trfail/nsamples),'$'

    deallocate(t,stderr,stdsig)

  end subroutine one_sample

  subroutine one_round(med,ndata,nbad,t,dt,sig)

    use quantilmean
    use sr1mean_function
    use noise

    logical, intent(in) :: med
    integer, intent(in) :: ndata, nbad
    real(REAL64), intent(out) :: t,dt,sig
    !    type(sr1_mean) :: funmean
    type(tsr1mean_tukey), target :: tukeymean
    type(tsr1mean_huber), target :: hubermean
!    type(tsr1mean_square), target :: squaremean
    type(tsr1mean_experimental), target :: bathmean
    class(tsr1mean_function), pointer :: funmean

    real(REAL64), dimension(:), allocatable, target :: data
    integer, dimension(:), allocatable :: ibad
    integer :: i,j,iter,iflag,trs
    real(REAL64) :: s, tmed, tol, t1, t2, tbest, dtbest, x, dtol
    logical :: converge,reli

    allocate(data(ndata))
    call gnoise_fill(0.0_REAL64,1.0_REAL64,data)
    if( nbad > 0 ) then
       allocate(ibad(nbad))
       do i = 1, nbad
          do
             call random_number(x)
             j = int(ndata*x)
             if( findloc(ibad(1:i-1),j,1) == 0 .and. j > 0 ) then
                ibad(i) = j
                exit
             end if
          end do
       end do
       do i = 1, nbad
          data(ibad(i)) = gnoise(-13.0_REAL64,1.0_REAL64)
          !       write(*,*) data(i)
       end do
!       write(*,*) ibad
       deallocate(ibad)
    end if

    call qmean(data,t,dt,s,.true.,reli,iflag)
    tmed = t
!    s = 1
    tol = 1e-2*s/sqrt(real(size(data)))

    block
      type(tsr1mean_experimental) :: hampelmean

      t = tmed
      call hampelmean%init(data,iflag)
      call trnewton(hampelmean,t,dt,s,0.01*tol,converge,iter,trs)
      if( .not. converge ) stop
      tbest = t
      dtbest = dt
      dtol = 1000*dt
      dtol = 1
    end block

    funmean => bathmean
!    funmean => tukeymean
!    funmean => hubermean
!    funmean%data => data
    call funmean%init(data,iflag)
    zmeanfun => funmean

    tol = 1e-2/sqrt(real(size(data)))
    if( med ) then
       t = tmed
    else
       t = 0
       s = 1
    end if
    call rnewton(funmean,t,dt,s,tol,converge,iter)
    t1 = t
    riter(iter) = riter(iter) + 1
    if( .not. converge .or. abs(t - tbest) > dtol  ) then
       rfail = rfail + 1
       write(*,*) 'rnewton fail:',t,tbest,converge
    end if

    if( med ) then
       t = tmed
    else
       t = 0
       s = 1
    end if
    call trnewton(funmean,t,dt,s,tol,converge,iter,trs)
    t2 = t
    triter(iter) = triter(iter) + 1
    trstep = trstep + trs
    if( .not. converge .or. abs(t - tbest) > dtol ) then
       trfail = trfail + 1
       write(*,*) 'trnewton fail:',t,tbest,converge
    end if

    ! zeroin
    call zmean(data,s,t-3*s,t+3*s,t,converge)
    if( zeroiter <= size(ziter) ) ziter(zeroiter) = ziter(zeroiter) + 1
    if( .not. converge .or. abs(t - tbest) > dtol ) then
       zfail = zfail + 1
       write(*,*) 'zeroin fail:',t,tbest,converge
    end if


!    if( converge .and. abs(t1 - t2) > 10*tol ) then
!       nfail = nfail + 1
!       write(*,*) t1,t2
!       trfail = trfail + 1
!       rfail = rfail + 1
!    end if

    sig = sqrt(real(size(data)))*dt

  end subroutine one_round

  subroutine zmean(x,s,tmin,tmax,t,reli)

    ! Zero mean estimates mean of robust function by Brent's method
    ! which is combination of bisection and inverse quadratic interpolation.
    ! It is important that it works without derivatives so there is
    ! no potential risk of divergence due to "re-descending" property
    ! of Tukey's function as in rnewton() with an improper initial.

    use fmm

    real(REAL64), dimension(:),intent(in) :: x
    real(REAL64), intent(in) :: s, tmin, tmax
    real(REAL64), intent(out) :: t
    logical, intent(out) :: reli

    real(REAL64) :: tol, d
    integer :: n

    n = size(x)
    tol = 0.01 * s / sqrt(real(n))
    zmean_scale = s
    zeroiter = 0
    t = zeroin(tmin,tmax,zfun,tol)

    ! tol gives precision of localisation by expected std.err.
    ! Number of calls of zfun will approximately log((tmax-tmin)/d).

    ! result too close to endpoints shows convergence fail
    d = 2*tol
    reli = abs(t - tmin) > d .and. abs(tmax - t) > d

  end subroutine zmean


  function zfun(t)

    use rfun

    real(REAL64) :: zfun
    real(REAL64), intent(in) :: t
    real(REAL64) :: s,g1,h11

    s = zmean_scale
    !psi = tukey((x-t)/s)
    !zfun = sum(psi) / size(psi)
    call zmeanfun%rnewton(t,s,g1,h11)
    zfun = g1
    zeroiter = zeroiter + 1

  end function zfun


end program rmean_newton
