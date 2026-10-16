!
! Matching on base of looking for nearest stars
!
! Copyright © 2017-2026 F.Hroch (hroch@physics.muni.cz)
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

module nearmatch

  use iso_fortran_env
  implicit none

contains

  subroutine listmatch(type,alpha,delta,xpix,ypix,acen,dcen,xcen,ycen, &
       scale,rot,refl,tol,i1,i2,status)

    use astrotrafo

    character(len=*), intent(in) :: type
    real(REAL64), dimension(:), intent(in) :: alpha,delta,xpix,ypix
    real(REAL64), intent(in) :: acen,dcen,xcen,ycen,scale,rot,refl,tol
    integer, dimension(:), allocatable, intent(out) :: i1,i2
    logical, intent(out) :: status
    real(REAL64), dimension(:), allocatable :: x,y,u,v
    type(AstrometryProjection) :: tp
    type(AstrometryGeometry) :: tg
    integer :: nsize,msize

    status = .false.

    if( size(alpha) /= size(delta) .or. size(xpix) /= size(ypix) ) &
         stop 'Listmach arrays does not conforms.'

    ! select maximum amount of stars by list
    nsize = size(alpha)
    msize = size(xpix)
    allocate(u(nsize),v(nsize),x(msize),y(msize))

    call tp % init(type,acen,dcen)
    call tp % proj(alpha,delta,u,v)

    call tg % init(scale,rot,refl)
    call tg % pix2uv(xpix-xcen,ypix-ycen,x,y)

    call pairmatch(u,v,x,y,tol,i1,i2)
    status = size(i1) > 0

  end subroutine listmatch


  subroutine wcsmatch(type,alpha,delta,rflux,xpix,ypix,flux,acen,dcen,xcen,ycen, &
       scale,rot,refl,frame,maxoff,maxrot,maxscale,bin,i1,i2,ntrial,status,verbose)

    ! The match guided by an approximate calibration (the WCS in the header
    ! of the frame): the centre acen,dcen [deg] at the pixel xcen,ycen,
    ! the scale [pix/deg], the rotation rot [deg], the reflection refl,
    ! known within the offset maxoff [pix], the rotation maxrot [deg]
    ! and the relative scale maxscale. The catalogue stars predicted
    ! on the frame by the calibration are rotated and scaled on a grid
    ! covering the tolerances; for every node, the differences of the
    ! positions of the brightest frame stars and the catalogue stars
    ! vote for the offset in a histogram of the bin size bin [pix].
    ! The node and the offset with the most votes give the transformation
    ! and the pairs of the nearest stars by it (the brighter first) are
    ! the match i1 (catalogue), i2 (frame); ntrial is the number of the
    ! transformations tried (the nodes and the bins of the offsets),
    ! for the significance of the solution. The vote needs three
    ! stars at least, the match five pairs.

    use astrotrafo

    character(len=*), intent(in) :: type
    real(REAL64), dimension(:), intent(in) :: alpha,delta,rflux,xpix,ypix,flux
    real(REAL64), intent(in) :: acen,dcen,xcen,ycen,scale,rot,refl, &
         maxoff,maxrot,maxscale,bin
    integer, dimension(:), intent(in) :: frame
    integer, dimension(:), allocatable, intent(out) :: i1,i2
    real(REAL64), intent(out) :: ntrial
    logical, intent(out) :: status
    logical, intent(in) :: verbose

    real(REAL64), parameter :: rad = 57.295779513082322865_REAL64
    ! the frame stars voting (the brightest), the catalogue stars per
    ! a frame star (scaled by the area of the margin)
    integer, parameter :: maxvote = 60, minvote = 3, minpairs = 5
    real(REAL64), parameter :: catfac = 3
    real(REAL64), dimension(:), allocatable :: u,v,xc,yc,xf,yf,xr,yr
    integer, dimension(:), allocatable :: icat, ifr, jcat, jfr
    integer, dimension(:,:), allocatable :: hist
    logical, dimension(:), allocatable :: used
    type(AstrometryProjection) :: tp
    type(AstrometryGeometry) :: tg
    real(REAL64) :: hw, hh, rmax, margin, drot, dscl, r, s, c, sn, dx, dy, &
         cell, tol, xoff, yoff, smin, d
    integer :: ncat, nfr, nd, nc, nb, nrot, nscl, ir, is, i, j, k, l, m, &
         best, ibest, jbest, rbest, sbest, nsup, npair

    status = .false.
    ntrial = 1
    allocate(i1(0),i2(0))

    ncat = size(alpha)
    nfr = size(xpix)
    if( ncat < minpairs .or. nfr < minpairs ) return

    ! the catalogue predicted on the frame (relative to the centre)
    allocate(u(ncat),v(ncat),xc(ncat),yc(ncat))
    call tp % init(type,acen,dcen)
    call tp % proj(alpha,delta,u,v)
    call tg % init(scale,rot,refl)
    call tg % uv2pix(u,v,xc,yc)

    ! the half-sizes of the frame around the centre, the margin
    hw = max(xcen - 0.5_REAL64, frame(1) + 0.5_REAL64 - xcen)
    hh = max(ycen - 0.5_REAL64, frame(2) + 0.5_REAL64 - ycen)
    rmax = hypot(hw,hh)
    margin = maxoff + rmax*(maxrot/rad + maxscale)

    ! the brightest catalogue stars in the frame with the margin,
    ! the brightest frame stars
    nd = min(nfr,maxvote)
    allocate(ifr(nfr),icat(ncat))
    call fluxorder(flux,ifr)
    call fluxorder(rflux,icat)
    m = 0
    do i = 1, ncat
       k = icat(i)
       if( abs(xc(k)) < hw + margin .and. abs(yc(k)) < hh + margin ) then
          m = m + 1
          icat(m) = k
       end if
    end do
    nc = min(m, max(30, nint(catfac*nd*((hw+margin)*(hh+margin))/(hw*hh))))
    if( nc < minvote ) return
    allocate(xf(nd),yf(nd),xr(nc),yr(nc))
    xf = xpix(ifr(1:nd)) - xcen
    yf = ypix(ifr(1:nd)) - ycen

    ! the grid of rotations and scales: the steps move the stars at
    ! the edge of the frame by a bin at most
    drot = min(maxrot, max(rad*bin/rmax, 0.01_REAL64))
    dscl = min(maxscale, max(bin/rmax, 1e-4_REAL64))
    nrot = 0
    nscl = 0
    if( maxrot > 0 ) nrot = ceiling(maxrot/drot)
    if( maxscale > 0 ) nscl = ceiling(maxscale/dscl)

    ! the histogram of the offsets in the cells of half of the bin,
    ! the votes are counted in the windows of 2 x 2 cells
    cell = bin / 2
    nb = 2*ceiling(maxoff/cell) + 2
    allocate(hist(nb,nb))
    best = 0
    ibest = 0
    jbest = 0
    rbest = 0
    sbest = 0
    do ir = -nrot, nrot
       r = ir*drot / rad
       c = cos(r)
       sn = sin(r)
       do is = -nscl, nscl
          s = 1 + is*dscl
          xr = s*(c*xc(icat(1:nc)) - sn*yc(icat(1:nc)))
          yr = s*(sn*xc(icat(1:nc)) + c*yc(icat(1:nc)))
          hist = 0
          do i = 1, nd
             do j = 1, nc
                dx = xf(i) - xr(j)
                dy = yf(i) - yr(j)
                if( abs(dx) < maxoff .and. abs(dy) < maxoff ) then
                   k = int((dx + maxoff)/cell) + 1
                   l = int((dy + maxoff)/cell) + 1
                   hist(k,l) = hist(k,l) + 1
                end if
             end do
          end do
          do k = 1, nb - 1
             do l = 1, nb - 1
                m = hist(k,l) + hist(k+1,l) + hist(k,l+1) + hist(k+1,l+1)
                if( m > best ) then
                   best = m
                   ibest = k
                   jbest = l
                   rbest = ir
                   sbest = is
                end if
             end do
          end do
       end do
    end do
    ! the trials: the nodes of the grid and the independent bins
    ntrial = real(2*nrot+1,REAL64)*real(2*nscl+1,REAL64)*(nb/2.0_REAL64)**2

    if( verbose ) write(error_unit,'(a,i0,a,i0,a,i0,a,f0.2,a,f0.4,a,es9.2e3,a)') &
         ' WCS-guided match: ',nd,' frame stars, ',nc,' catalogue stars, best vote ', &
         best,' at rotation ',rbest*drot,' deg, scale ',1 + sbest*dscl,' (',ntrial,' trials).'
    if( best < minvote ) return

    ! the offset: the mean of the differences in the best window
    r = rbest*drot / rad
    c = cos(r)
    sn = sin(r)
    s = 1 + sbest*dscl
    xr = s*(c*xc(icat(1:nc)) - sn*yc(icat(1:nc)))
    yr = s*(sn*xc(icat(1:nc)) + c*yc(icat(1:nc)))
    xoff = 0
    yoff = 0
    nsup = 0
    do i = 1, nd
       do j = 1, nc
          dx = xf(i) - xr(j)
          dy = yf(i) - yr(j)
          if( abs(dx) < maxoff .and. abs(dy) < maxoff ) then
             k = int((dx + maxoff)/cell) + 1
             l = int((dy + maxoff)/cell) + 1
             if( (k == ibest .or. k == ibest + 1) .and. (l == jbest .or. l == jbest + 1) ) then
                xoff = xoff + dx
                yoff = yoff + dy
                nsup = nsup + 1
             end if
          end if
       end do
    end do
    xoff = xoff / nsup
    yoff = yoff / nsup

    ! the pairs of the nearest stars by the transformation within the bin
    ! and the error of the grid at the edge, all stars, the brighter first
    deallocate(xr,yr)
    allocate(xr(ncat),yr(ncat),used(ncat),jcat(nfr),jfr(nfr))
    xr = s*(c*xc - sn*yc) + xoff
    yr = s*(sn*xc + c*yc) + yoff
    tol = bin + 0.5*rmax*(drot/rad + dscl)
    used = .false.
    npair = 0
    do i = 1, nfr
       k = ifr(i)
       l = 0
       smin = tol
       do j = 1, ncat
          if( used(j) ) cycle
          d = hypot(xpix(k) - xcen - xr(j), ypix(k) - ycen - yr(j))
          if( d < smin ) then
             smin = d
             l = j
          end if
       end do
       if( l > 0 ) then
          used(l) = .true.
          npair = npair + 1
          jcat(npair) = l
          jfr(npair) = k
       end if
    end do

    if( verbose ) write(error_unit,'(a,2f9.2,a,i0,a,f0.2,a)') &
         ' WCS-guided match: offset ',xoff,yoff,' pix, ',npair, &
         ' pairs within ',tol,' pix.'

    if( npair < minpairs ) return
    deallocate(i1,i2)
    allocate(i1(npair),i2(npair))
    i1 = jcat(1:npair)
    i2 = jfr(1:npair)
    status = .true.

  end subroutine wcsmatch


  subroutine fluxorder(f,idx)

    ! the indices of f in the descending order (the brightest first),
    ! the equal values in their original order (heapsort by the pair
    ! of the value and the index)

    real(REAL64), dimension(:), intent(in) :: f
    integer, dimension(:), intent(out) :: idx
    integer :: n, i, j, k, t, root, child

    n = size(f)
    idx = [(i, i = 1, n)]
    ! a heap with the last star of the order on the top
    ! (the faintest), moved to the end one by one
    do k = n/2, 1, -1
       call sift(k,n)
    end do
    do j = n, 2, -1
       t = idx(1)
       idx(1) = idx(j)
       idx(j) = t
       call sift(1,j-1)
    end do

  contains

    logical function after(a,b)
      ! the star a is in the order after b: fainter, or as bright later
      integer, intent(in) :: a,b
      after = f(a) < f(b) .or. (f(a) == f(b) .and. a > b)
    end function after

    subroutine sift(start,last)
      integer, intent(in) :: start,last
      root = start
      do
         child = 2*root
         if( child > last ) exit
         if( child < last ) then
            if( after(idx(child+1),idx(child)) ) child = child + 1
         end if
         if( after(idx(child),idx(root)) ) then
            t = idx(root)
            idx(root) = idx(child)
            idx(child) = t
            root = child
         else
            exit
         end if
      end do
    end subroutine sift

  end subroutine fluxorder


  subroutine pairmatch(u,v,x,y,tol,i1,i2)

    ! The pairs of the nearest reference (u,v) and frame (x,y) stars
    ! closer than tol, in the same units. Both lists are supposed sorted
    ! by flux: every reference star takes its nearest frame star,
    ! a frame star is taken by the brightest reference star only.

    real(REAL64), dimension(:), intent(in) :: u,v,x,y
    real(REAL64), intent(in) :: tol
    integer, dimension(:), allocatable, intent(out) :: i1,i2
    integer, dimension(:), allocatable :: id1,id2
    real(REAL64) :: smin,smin0,s
    integer :: nsize,msize,n,i,j,idx
    logical :: found

    nsize = size(u)
    msize = size(x)
    allocate(id1(nsize),id2(nsize))

    smin0 = tol

    n = 0
    do i = 1, nsize
       smin = smin0
!       write(*,'(i5,5f15.5)') i,u(i),v(i),x(i),y(i),smin
!       write(*,*) alpha(i),delta(i),smin
       idx = 0
       do j = 1, msize
          s = hypot(u(i) - x(j),v(i) - y(j))
          if( s < smin ) then
             ! the nearest and the brighest star (a flux-sorted list is assumed)
             smin = s
             idx = j
          end if
       end do
       if( idx > 0 ) then

          ! check whatever the star is already presented,
          ! it will be ignored, because brighter ones are prefered
          ! (and we are supposing flux-sorted lists).
          found = .false.
          do j = 1, n
             if( idx == id2(j) ) then
                found = .true.
                exit
             end if
          end do

          if( .not. found ) then
             n = n + 1
             id1(n) = i
             id2(n) = idx
          end if

          ! write(*,'(f10.3,2i5,2f15.2)') smin,i,idx,u(i)-x(idx),v(i)-y(idx)

          if( n == size(id1) ) exit
       end if
    end do

    allocate(i1(n),i2(n))
    i1 = id1(1:n)
    i2 = id2(1:n)

  end subroutine pairmatch


end module nearmatch
