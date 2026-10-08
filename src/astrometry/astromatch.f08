!
! astromatch - matching engine
!
!
! Copyright © 2011-2026 F.Hroch (hroch@physics.muni.cz)
!
! This file is part of Munipack.
!
! Munipack is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Munipack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
!

module astromatcher

  use iso_fortran_env
  implicit none
  private

  ! maximal number of loop over frame stars
  integer, parameter :: maxloop = 3

  public :: astromatch, uvhisto

contains

  subroutine astromatch(minmatch,maxmatch,type,rec,dec,rflux,xpix,ypix,flux, &
       acen,dcen,xcen,ycen,id1,id2,sig1,sig2,flxmin,flxmax,fluxtest,luckymatch,&
       verbose,status,nskip)

    ! nskip (optional): on input, the number of chunk pairs to skip
    ! (a match found there has been rejected by a previous attempt),
    ! on output, the number of the chunk pair providing the match,
    ! so a next call with it continues by next chunks of stars

    use matcher
    use astrotrafo

    integer, intent(in) :: minmatch, maxmatch, luckymatch
    logical, intent(in) :: verbose, fluxtest
    character(len=*), intent(in) :: type
    real(REAL64), intent(in) :: acen,dcen,xcen,ycen,sig1,sig2,flxmin,flxmax
    real(REAL64), dimension(:), intent(in), target :: rec, dec, &
         xpix, ypix, flux, rflux
    integer, dimension(:), allocatable, intent(out) :: id1,id2
    logical, intent(out) :: status
    integer, intent(in out), optional :: nskip

    real(REAL64), parameter :: rad = 57.29577951308232_REAL64

    real(REAL64), dimension(:), allocatable :: u,v,x,y
    logical, dimension(:), allocatable :: idx1, idx2
    type(AstrometryProjection) :: tp
    integer :: n1,n2,m,n,l,k,ichunk,maxskip

    status = .false.
    ichunk = 0
    maxskip = 0
    if( present(nskip) ) maxskip = nskip

    if(verbose) write(error_unit,*) "=== Matching ==="

    n1 = size(rec)
    n2 = size(xpix)

    allocate(u(n1),v(n1),x(n2),y(n2),idx1(n1),idx2(n2))

    ! projection of spherical coordinates
    call tp % init(type,acen,dcen)
    call tp % proj(rec,dec,u,v)

    ! make it centered
    x = xpix - xcen
    y = ypix - ycen

    ! filter close encounters
    call encounter(rec,dec,cos(dcen/rad),2*sig1,idx1,verbose)
    call encounter(xpix,ypix,1.0_REAL64,2*sig2,idx2,verbose)

    if( verbose ) then
       write(error_unit,'(2(a,i0),a)') " Matching of ",n1, &
            " reference stars against to ",n2," frame stars:"
       write(error_unit,*) &
            "#: {seq1} -> {seq2} | scale, rel.flux: sq.(angle, scale) < Xi2(0.95)"
    end if

    ! the stars are processed in chunks k..k+n2 (l..l+n1) of maxmatch
    ! length at most, the last chunk contains the remaining stars
    k = 1
    do m = 1, min((size(x)-1)/maxmatch + 1, maxloop)
       n2 = min(maxmatch,size(x)-k+1) - 1
       l = 1
       do n = 1, (size(u)-1)/maxmatch + 1
          n1 = min(maxmatch,size(u)-l+1) - 1
          ichunk = ichunk + 1
          if( ichunk <= maxskip ) then
             ! this chunk pair has been already tried, the match rejected
             l = l + maxmatch
             cycle
          end if
          if(verbose) write(error_unit,'(2(2(a,i0),a))') &
               " Reference stars ",l," -- ",l+n1,".", &
               " Frame stars ",k," -- ",k+n2,"."
          call match(u(l:l+n1),v(l:l+n1),rflux(l:l+n1),idx1(l:l+n1), &
               x(k:k+n2),y(k:k+n2),flux(k:k+n2),idx2(k:k+n2), &
               minmatch,maxmatch,sig1,sig2,flxmin,flxmax,fluxtest,id1,id2,&
               luckymatch,verbose,status)
          if( status ) then
             id1 = l + id1 - 1
             id2 = k + id2 - 1
             if( present(nskip) ) nskip = ichunk
             goto 99
          end if
          l = l + maxmatch
       end do
       k = k + maxmatch
    end do
99  continue

  end subroutine astromatch

  subroutine encounter(x,y,cosd,tol,idx,verbose)

    real(REAL64), dimension(:), intent(in) :: x,y
    real(REAL64), intent(in) :: cosd,tol
    logical, dimension(:), intent(out) :: idx
    logical, intent(in) :: verbose

    integer :: n,i,j
    real(REAL64) :: r

    n = size(x)
    idx = .true.
    do i = 1, n
       do j = i+1, n
          r = abs((x(i) - x(j))*cosd) + abs(y(i) - y(j))
          if( r < tol ) then
             idx(j) = .false.
             if( verbose ) &
                  write(error_unit,'(a,2i5,4f11.5,es9.1)') &
                  'Astromatch recognised an encounter:',i,j,x(i),y(i),x(j),y(j),r
          end if
       end do
    end do

  end subroutine encounter

!!$  subroutine matchprint(id1,id2,c,t,a,d,q)
!!$
!!$    use matcher
!!$
!!$    integer, dimension(:), intent(in) :: id1,id2
!!$    real(REAL64), intent(in) :: c,t,a,d,q
!!$    real(REAL64), dimension(size(id1)) :: u,v
!!$    integer :: i,n
!!$    real(REAL64) :: d3
!!$    character(len=80) :: fmt
!!$
!!$    if( size(id1) /= size(id2) ) stop 'Astromatch size(id1) /= size(id2)'
!!$
!!$    n = size(id2)
!!$
!!$    if( verbose ) then
!!$
!!$       write(fmt,"(a,i0,a,i0,a)") "(i0,a,",n,"(1x,i0),a,",n, &
!!$            "(1x,i0),a,f0.1,2x,f0.3,a,2(2x,g0.3),a,f0.1)"
!!$       write(error_unit,fmt) n,":",id1,' -> ',id2,' | ',c,t,': ',a,d,' < ',q
!!$
!!$    end if
!!$
!!$    return
!!$
!!$    if( plog ) then
!!$
!!$       write(fmt,"(a,i0,a)") "(a,i0,",2*n,"(1x,f0.3),f6.3,f8.5)"
!!$       write(*,fmt) "=M> ",n,(xstar(id2(i)),ystar(id2(i)),i=1,n),a,d
!!$
!!$       do i = 3,n
!!$          call triangle(xstar(i-2),ystar(i-2),xstar(i-1), &
!!$               ystar(i-1),xstar(i),ystar(i),d3,u(i),v(i))
!!$       end do
!!$
!!$       write(fmt,"(a,i0,a)") "(a,i0,",2*(n-2),"(1x,f0.3))"
!!$       write(*,fmt) "=MUV> ",n-2,(u(i),v(i),i=3,n)
!!$    end if
!!$
!!$  end subroutine matchprint


!!$  subroutine progress2(ns,ntot,xamin,xdmin,id1,id2)
!!$
!!$    use matcher
!!$    use estimator
!!$
!!$    integer, intent(in) :: ns,ntot
!!$    integer, dimension(:), intent(in), optional :: id1,id2
!!$    real(REAL64), intent(in), optional :: xamin,xdmin
!!$    character(len=80) :: fmt,fm
!!$    real(REAL64), dimension(:), allocatable :: u,v
!!$    integer :: n,m,i,i1,i2,i3
!!$    real(REAL64) :: ac,dc,sc,dsc,rot,drot,rms,amin,dmin,d3,refl,cooerr
!!$    real(REAL64), dimension(:), allocatable :: alpha,delta,x,y
!!$
!!$    if( .not. plog ) return
!!$
!!$    if( present(id2) ) then
!!$
!!$       n = size(id2)
!!$
!!$       allocate(u(n),v(n))
!!$
!!$       ! FIT
!!$       allocate(alpha(n),delta(n),x(n),y(n))
!!$       do i = 1,n
!!$          i1 = id1(i)
!!$          i2 = id2(i)
!!$          alpha(i) = astar(i1)
!!$          delta(i) = dstar(i1)
!!$          x(i) = xstar(i2) - xcen
!!$          y(i) = ystar(i2) - ycen
!!$       end do
!!$
!!$       ac = acen
!!$       dc = dcen
!!$       call inestim(ptype,alpha,delta,x,y,ac,dc,sc,dsc,rot,drot, &
!!$            refl,cooerr,verbose)
!!$       deallocate(alpha,delta,x,y)
!!$
!!$       write(*,'(a,en15.3,en25.15,en25.15,4en25.15)') &
!!$            "=MFIT> ",rms,sc,rot,ac,dc,xcen,ycen
!!$
!!$       if( xdmin < 0.1*huge(xdmin) ) then
!!$          amin = xamin
!!$          dmin = xdmin
!!$       else
!!$          dmin = -1.0
!!$          amin = -1.0
!!$       end if
!!$
!!$       do i = 3,n
!!$          i1 = id2(i)
!!$          i2 = id2(i - 1)
!!$          i3 = id2(i - 2)
!!$          call triangle(xstar(i3),ystar(i3),xstar(i2),ystar(i2), &
!!$               xstar(i1),ystar(i1),d3,u(i),v(i))
!!$       end do
!!$
!!$       write(fmt,"(a,i0,a,i0,a)") "(a,2(1x,g0.7),1x,i0,",2*n,"(1x,f0.3),1x,i0,"
!!$       m = n - 2
!!$       write(fm,'(i0,a)') 2*(n-2),"(1x,f0.3))"
!!$
!!$       fmt = trim(fmt)//fm
!!$       write(*,fmt) "=MPROGRES2> ",amin,dmin,n, &
!!$            (xstar(id2(i)),ystar(id2(i)),i=1,n),m,(u(i),v(i),i=3,n)
!!$       deallocate(u,v)
!!$
!!$    else
!!$       write(*,'(a,2(1x,i0))') "=MPROGRESO> ",ns,ntot
!!$    end if
!!$
!!$
!!$  end subroutine progress2

  subroutine uvhisto(x,y)

    use matcher

    integer, parameter :: nhist = 100

    real(REAL64), dimension(:), intent(in) :: x,y
    integer, allocatable, dimension(:,:) :: hist
    real(REAL64) :: u, v, d3
    integer :: i,j,k,m,n,c
    character(len=80) :: fmt

    allocate(hist(nhist,nhist))
    hist = 0

    ! the number of triangles is limited to enable an acceptable user response
    c = max(size(x)/137,1)

    do i = 1, size(x) - 2,c
       do j = i+1, size(x) - 1,c
          do k = j+1, size(x),c
             call triangle(x(i),y(i),x(j),y(j),x(k),y(k),d3,u,v)
             m = nint(nhist*u)
             n = nint(nhist*v)
             if( 1 <= m .and. m <= nhist .and. 1 <= n .and. n <= nhist ) &
                  hist(m,n) = hist(m,n) + 1
          end do
       end do
    end do

    write(fmt,'(a,i0,a)') "(a,1x,i0,1x,",nhist**2,"(1x,i0))"
    write(*,fmt) "=MTRI> ",nhist,hist

    deallocate(hist)

  end subroutine uvhisto

end module astromatcher
