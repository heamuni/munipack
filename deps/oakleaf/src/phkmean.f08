!
!  Robust estimation of mean photon rates
!
!  Copyright © 2020 F.Hroch (hroch@physics.muni.cz)
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

module photonrate_kmean

  use photonrate_common
  use iso_fortran_env

  implicit none
  private

  public :: mucluster

  logical :: debug_graph = .false.

contains

  subroutine mucluster(mju,iset,K,kmean,tau,r)


    use quantilmean
    use quicksorter
    use medians

    ! find means of mju by k-mean clustering method
    ! https://en.wikipedia.org/wiki/K-means_clustering
    !
    ! The method is designed by description of k-means in Numerical
    ! recipes. The code is completedlly different, and significantly
    ! simpler; thanks to Fortran array capabilities.
    !
    ! The method depends, of course, on proper initialisation, else
    ! an unwanted local minimum of the clustering is encountered.
    ! As a init strategy, I try all of them: random inits, init by a constant.
    ! All of them fails. Rather, the search of peaks in histogram are used
    ! as starting points. The histogram bins widths are proportional
    ! to sigma of data. The whole init is a heavy routine, but, in this case,
    ! the nexr clustering, by the expectation-minimalisation, takes a few
    ! iterations only.  Moreover, the doubts with choice of the total count
    ! of clusters naturaly disapperar, as the groups are indicated
    ! by the peaks in the histogram.

    real(REAL64), dimension(:), intent(in out) :: mju
    real(REAL64), intent(in) :: tau
    integer, dimension(:), intent(out) :: iset
    integer, intent(out) :: K
    real(REAL64), dimension(:), allocatable, intent(out) :: kmean
    real(REAL64), intent(in), optional :: r

    real(REAL64), dimension(:), allocatable :: x
    integer, dimension(:), allocatable :: kset
    real(REAL64) :: rho
    integer :: i,l,m,iter,jitter
    logical :: change, empty

    ! It's a rought estimate of the basic bin width by scatter of mju
    ! scaled on pht rather than cts. The estimate is simplified sigma
    ! of rate distribution with all phi and 1/tau terms omitted, together
    ! with assumption tau ~ t(i), and divided by tau.

    if( present(r) ) then
       rho = r
    else
       rho = 0
    end if

    ! The estimate of total count of clusters -- K -- by search
    ! for single peaks in a histogram.
!    call histpeak(mju,rho,minval(mju),maxval(mju),K,kset,kmean)
    call sgroup(mju,rho,K,kset,kmean)

    ! Both working arrays kmean, x are allocated without checking
    ! of some available memmory . Is this a proper way?

    if( debug_iter ) then
       block
         real(REAL64) :: mag

         do l = 1, size(kmean)
            if( kmean(l) > 0 ) then
               mag = 25-2.5*log10(kmean(l))
            else
               mag = 99
            end if
            write(error_unit,*) 'Init k-means:',real(mag),real(kmean(l))
         end do
       end block
    end if


    do jitter = 1, 666

       iset = 0

       do iter = 1, 666

          ! E-step
          change = .false.
          do i = 1, size(mju)
!             l = minloc(abs(mju(i) - kmean)/sqrt(mju(i)+rho**2/2)/3,1)
             l = minloc(abs(mju(i) - kmean),1)
             if( l /= iset(i) ) change = .true.
             iset(i) = l
          end do
          !       write(*,*) 'iset:',iset

          if( .not. change ) exit

          ! M-step
          do l = 1, K
             x = pack(mju,iset == l)
             m = size(x)
             if( m > 2 ) then
                kmean(l) = median(x)
             else if( m == 2 .or. m == 1 ) then
                kmean(l) = sum(x) / m
             else if ( m == 0 ) then
                continue
!                if( debug ) write(error_unit,*) 'phmean(): lost class:',l
                !             kmean(l) = -1
             end if
          end do
          !       write(*,*) real(kmean)

       end do

       empty = .false.
       do l = 1, size(kmean)
          if( count(iset == l) == 0 ) then
             ! en empty set

             kmean(l:K-1) = kmean(l+1:K)
             K = K - 1
             empty = .true.
             if( debug ) write(error_unit,*) 'phmean(): class ',l,' is empty.'
          end if
       end do

       if( .not. empty ) exit

       if( empty ) then
          block
            real(REAL64), dimension(:), allocatable :: xmean

            xmean = kmean(1:K)
            kmean = xmean

          end block
       end if

    end do

    ! filter out all empty clusters
!!$    block
!!$      real(REAL64), dimension(:), allocatable :: xmean
!!$      integer, dimension(:), allocatable :: xset
!!$
!!$      xset = k
!!$
!!$      do l = 1, K
!!$         m = count(iset == l)
!!$         if( m == 0 ) then
!!$
!!$
!!$      K = count(kmean > 0)
!!$      xmean = pack(kmean,kmean > 0)
!!$      kmean = xmean
!!$    end block

!    write(*,*) 'Final k-means:',real(kmean)
!    write(*,*) 'Final k-sigs:',real(sig)
!    write(*,*) 'iset:',iset

    do l = 1,K
       where( iset == l )
          mju = kmean(l)
       end where
    end do

!    K = size(mju)
!    iset = [ (l,l=1,K) ]
!    kmean = mju


    if( debug .and. size(mju) < 666 ) then
       !    nset = K
       block
         real(REAL64) :: mag
         integer, dimension(:), allocatable :: idx
         integer :: i,m

         allocate(idx(size(iset)))
         write(error_unit,*) 'K-clusters:'
         do l = 1, K
!            if( kmean(l) > 0 ) then
            mag = real(25 - 2.5*log10(kmean(l)))
!            else
!               mag = 99
!            end if
            m = 0
            do i = 1, size(iset)
               if( iset(i) == l ) then
                  m = m + 1
                  idx(m) = i
               end if
            end do
            write(error_unit,'(f8.3,1x,g0,2x,*(i0,1x))') &
                 mag,real(kmean(l)),idx(1:m)!pack(iset,iset == l)
         end do
       end block
    end if

!    stop 0

  end subroutine mucluster


  subroutine histpeak(mju,rho,mumin,mumax,K,kset,kmean)

    real(REAL64), dimension(:), intent(in) :: mju
    real(REAL64), intent(in) :: mumin, mumax, rho
    integer, intent(out) :: K
    integer, dimension(:), allocatable, intent(out) :: kset
    real(REAL64), dimension(:), allocatable, intent(out) :: kmean

    real(REAL64), dimension(:), allocatable :: u
    integer, dimension(:), allocatable :: hist
    character(len=80) :: errmsg
    real(REAL64) :: c,d,q,qmin!,s, d, smin, smax, q, c
    integer :: i,m, nhist, stat

!    smin = sqrt(mumin + rho**2)
!    smax = sqrt(mumax + rho**2)
!    s = 1 / (1/smin + 1/smax)
!    ds = (smax - smin) / (mumax - mumin)


!    d = 3*s

    c = 3
    d = c*sqrt(mumin + rho**2)
    nhist = 1
    qmin = max(mumin - d/2,0.5*mumin)
    q = qmin
!    write(error_unit,*) 'sss',mumin, mumax,d,q
    do
!       q2 = q**2
!       s = smin + nhist * ds
!       write(error_unit,*) 'bound:',nhist,q,c*sqrt(q)
       if( q > mumax ) exit
       nhist = nhist + 1
       q = q + c*sqrt(q + rho**2)
!       write(error_unit,*) real(q),real(sqrt(q))
    end do
    nhist = nhist + 1
!    write(error_unit,*) 'nhist=',nhist, real(rho)

!    d = (mumax - mumin) / size(mju)
!    nhist = max(nint((mumax - mumin)/d),1)
!    nhist = max(nint(d),1)
!    nhist = nint((mumax - mumin)/d) + 1
!    write(error_unit,*) '>',mumin, mumax,rho,d,nhist

    allocate(hist(0:nhist+1),u(0:nhist+1),stat=stat,errmsg=errmsg)
    if( stat /= 0 ) then
       write(error_unit,*) 'Error in histpeak: `',trim(errmsg),"'"
       return
    end if

    ! non-equidistant grid
    u(0) = qmin !mumin - d/2
!    q = mumin - d/2
    do i = 1, nhist+1
!       q2 = q**2
       u(i) = u(i-1) + c*sqrt(u(i-1) + rho**2)
!       u(i) = u(i-1) + d!sqrt(u(i-1) + d)
    end do

    hist = 0
!    m = 1
    do i = 1, size(mju)
       do m = 1, nhist+1
          if( u(m-1) <= mju(i) .and. mju(i) < u(m) ) then
             hist(m) = hist(m) + 1
 !            exit
          end if
       end do
    end do

    if( debug_graph ) then
       open(1,file='/tmp/hkmean')
       do i = 0, nhist
          write(1,*) u(i), hist(i)
       end do
       close(1)
    end if

    goto 10
    K = 0
    allocate(kset(nhist))
    kset = 0
    do i = 1, nhist
       if( hist(i) >= hist(i-1) .and. hist(i) > hist(i+1) ) then
          ! The asymetry in this condition solves doubts for the equal
          ! frequency hist(i) == hist(i+-1). In the case, it selects
          ! the first (in the leftmost order) coordinate of the peak;
          ! the next bins of the plato are ignored.
          K = K + 1
          kset(K) = i
!          write(*,*) K,u(i)
       end if
    end do
10  continue

    K = 0
    allocate(kset(nhist))
    kset = 0
    do i = 1, nhist
       if( hist(i) > 0 ) then
          K = K + 1
          kset(K) = i
       end if
    end do

    if( K > 0 ) then
       kmean = u(kset(1:K))
    else
       K = 1
       kmean = u(kset(1:1))
    end if

!    deallocate(kmean)
!    K = size(mju) - 1
!    allocate(kmean(K))
!    kmean = sum(mju) / size(mju)


    deallocate(u,hist)

  end subroutine histpeak

  subroutine sgroup(mju,rho,K,kset,kmean)

    use quicksorter

    real(REAL64), dimension(:), intent(in) :: mju
    real(REAL64), intent(in) :: rho
    integer, intent(out) :: K
    integer, dimension(:), allocatable, intent(out) :: kset
    real(REAL64), dimension(:), allocatable, intent(out) :: kmean

    real(REAL64), dimension(:), allocatable :: x
    real(REAL64) :: d,s
    integer :: i,j

    x = mju
    call qsort(x)
!    write(*,*) x
    kmean = x
!    kset = 0

    do j = 1, size(mju)

       K = 1
       kmean(1) = x(1)
!       kset(1) = 1
       do i = 2, size(x)
          d = x(i) - x(i-1)
          s = sqrt(x(i) + 0*rho**2/2)
!          if( d < 3*s ) then
          if( d < -s ) then
             kmean(K) = (kmean(K) + x(i))/2
!             kset(K) = i
          else
             K = K + 1
             kmean(K) = x(i)
!             kset(K) = i
          end if
       end do

       if( K == size(x) ) exit
       x = kmean(1:K)

    end do
    x = kmean(1:K)
    kmean = x

!    write(*,*) K,kmean

  end subroutine sgroup


end module photonrate_kmean
