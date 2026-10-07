
module toolbox

  use iso_fortran_env
  implicit none

contains

  function arith(x)

    real(REAL64) :: arith
    real(REAL64), dimension(:), intent(in) :: x

    arith = sum(x) / size(x)

  end function arith

  function stddev(r)

    real(REAL64) :: stddev
    real(REAL64), dimension(:), intent(in) :: r
    integer :: n

    n = size(r)
    if( n > 1 ) then
       stddev = sqrt(sum(r**2) / (n - 1))
    else
       stddev = 0
    end if

  end function stddev

  subroutine histogram(x, xmin, xmax, nbines, filename, prob)

    real(REAL64), dimension(:), intent(in) :: x
    real(REAL64), intent(in) :: xmin, xmax
    integer, intent(in) :: nbines
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: prob
    integer, dimension(:), allocatable :: hist
    real(REAL64) :: d,t,s
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

    s = 1
    if( present(prob) ) then
       if( prob ) s = sum(hist) * d
!       if( prob ) s = real(sum(hist)) !/ real(size(x))
    end if

    if( abs(s) < 10*epsilon(s) ) then
       write(*,*) 'Warning: the histogram is full of emptiness.'
       s = 1
    end if

    open(1,file=filename)
    do n = -nbins,nbins
!       write(1,*) (n*d + t), real(hist(n))/ real(size(x))
       write(1,*) (n*d + t), hist(n)/ s
    end do
    close(1)

    deallocate(hist)

  end subroutine histogram


  subroutine histogram2d(x, y, xmin, xmax, ymin, ymax, nbines, filename, prob)

    real(REAL64), dimension(:), intent(in) :: x, y
    real(REAL64), intent(in) :: xmin, xmax, ymin, ymax
    integer, intent(in) :: nbines
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: prob
    integer, dimension(:,:), allocatable :: hist
    real(REAL64) :: dx,tx,dy,ty,s
    integer :: n,m,l,nbins

    if( nbines < 1 ) then
       nbins = int(log(real(size(x)))/0.7) + 1
    else
       nbins = nbines
    end if
    nbins = nbins / 2

    allocate(hist(-nbins:nbins,-nbins:nbins))

    dx = (xmax - xmin) / (2*nbins + 1)
    tx = (xmax + xmin) / 2
    dy = (ymax - ymin) / (2*nbins + 1)
    ty = (ymax + ymin) / 2

    hist = 0
    do l = 1, size(x)
       n = nint((x(l) - tx) / dx)
       m = nint((y(l) - ty) / dy)
       if( -nbins <= n .and. n <= nbins .and. -nbins <= m .and. m <= nbins )then
          hist(n,m) = hist(n,m) + 1
       end if
    end do

    s = 1
    if( present(prob) ) then
       if( prob ) s = sum(hist) * dx * dy
    end if

    if( abs(s) < 10*epsilon(s) ) then
       write(*,*) 'Warning: the histogram is full of emptiness.'
       s = 1
    end if

    open(1,file=filename)
    do n = -nbins,nbins
       do m = -nbins,nbins
          write(1,*) (n*dx + tx), (m*dy + ty), hist(n,m)/s
       end do
       write(1,*)
    end do
    close(1)

    deallocate(hist)

  end subroutine histogram2d


end module toolbox
