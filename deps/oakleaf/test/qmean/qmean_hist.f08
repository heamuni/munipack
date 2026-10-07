

program qmean_hist

  use toolbox
  use iso_fortran_env
  implicit none

  integer, parameter :: ndata = 100
  integer, parameter :: nsamples = 99999
  integer, parameter :: step = 10
  real(REAL64), dimension(:), allocatable :: t,stderr,stdsig, stdmad
  integer :: i, nfail

  allocate(t(nsamples),stderr(nsamples),stdsig(nsamples),stdmad(nsamples))

  nfail = 0
  do i = 1, nsamples
     call estimator(t(i),stderr(i),stdsig(i),stdmad(i))
  end do

  call histogram(t,real(-0.5,REAL64),real(0.5,REAL64),100,'/tmp/qmean_median',.true.)
  call histogram(stderr,real(0,REAL64),real(3/sqrt(real(ndata)),REAL64),150,'/tmp/qmean_stderr',.true.)
  call histogram(stdsig,real(0,REAL64),real(3,REAL64),100,'/tmp/qmean_stdsig',.true.)
!  call histogram(stdmad,real(0,REAL64),real(3,REAL64),100,'/tmp/qmean_stdmad',.true.)

  open(1,file='/tmp/qmean_hist.par')
  write(1,*) ndata
  close(1)

contains

  subroutine estimator(t,stderr,stdsig,stdmad)

    use noise
    use quantilmean
!    use quicksorter

    real(REAL64), intent(out) :: t,stderr,stdsig,stdmad

    real(REAL64), dimension(:), allocatable, target :: data
    integer :: i,flag
    logical :: reliable

    allocate(data(ndata))
    call gnoise_fill(0.0_REAL64,1.0_REAL64,data)
    do i = 2, ndata, step
       data(i) = gnoise(-15.0_REAL64,1.0_REAL64)
!       write(*,*) data(i)
    end do

!    call qmean(data,t,stderr,stdmad,.true.,reliable,flag)
    call qmean(data,t,stderr,stdsig,'MAD',666,reliable,flag)

    if( .not. reliable .or. flag /= 0  ) then

       ! save data
       open(3,file='/tmp/qmean_datafail',form='unformatted')
       write(3) size(data)
       write(3) data
       close(3)
       stop 'in estimator'
    end if

!!$    data = abs(data - t)
!!$    call qsort(data)
!!$    if( mod(ndata,2) == 0 ) then
!!$       m = ndata / 2
!!$       stdmad = (data(m) + data(m+1)) / real(2.0,REAL64)
!!$    else
!!$       stdmad = data((ndata+1)/2)
!!$    end if
!!$    stdmad = stdmad / 0.6745

  end subroutine estimator

end program qmean_hist
