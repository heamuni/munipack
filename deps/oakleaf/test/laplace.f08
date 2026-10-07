

program testlaplace

  ! https://en.wikipedia.org/wiki/Laplace_distribution

  use noise
  use toolbox
  use quantilmean

  implicit none

  integer, parameter :: n = 6666
  integer, parameter :: kind = selected_real_kind(15)
  character(len=80) :: histname, cdfname
  real(kind), dimension(:), allocatable :: x
  integer :: i

  call get_command_argument(1,histname)
  call get_command_argument(2,cdfname)

  allocate(x(n))
  do i = 1, n
     x(i) = lnoise(0.0_kind,1.0_kind)
  end do

  call histogram(x,-5.0_kind,5.0_kind,100,histname)

  call cdfun(x,cdfname)

contains

  subroutine cdfun(x,filename)

    use oakleaf

    real(kind), dimension(:), intent(in) :: x
    character(len=*), intent(in) :: filename
    real(kind), dimension(:), allocatable :: xcdf, ycdf
    integer :: i,n

    n = size(x)
    allocate(xcdf(n),ycdf(n))

    call ecdf(x,xcdf,ycdf)

    open(1,file=filename)
    do i = 1, size(xcdf)
       write(1,*) xcdf(i),ycdf(i)
    end do
    write(1,*) xcdf(n)+100,1.0
    close(1)

  end subroutine cdfun


end program testlaplace
