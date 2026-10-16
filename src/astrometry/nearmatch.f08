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
       scale,rot,refl,frame,maxoff0,maxrot,maxscale,bin,i1,i2,ntrial,status,verbose, &
       reflboth)

    ! The match guided by an approximate calibration (the WCS in the header
    ! of the frame, or the pointing and the scale): the centre acen,dcen [deg]
    ! at the pixel xcen,ycen, the scale [pix/deg], the rotation rot [deg],
    ! the reflection refl, known within the offset maxoff0 [pix; at most
    ! the diagonal of the frame: beyond, the pointing is unknown], the rotation
    ! maxrot [deg] and the relative scale maxscale. The catalogue stars
    ! predicted on the frame by the calibration are rotated and scaled on
    ! a grid covering the tolerances; for every node, the differences of the
    ! positions of the brightest frame stars and the catalogue stars
    ! vote for the offset in a histogram of the bin size bin [pix].
    ! The node and the offset with the most votes give the transformation
    ! and the pairs of the nearest stars by it (the brighter first) are
    ! the match i1 (catalogue), i2 (frame); ntrial is the number of the
    ! transformations tried (the nodes and the bins of the offsets),
    ! for the significance of the solution. The vote needs three
    ! stars at least, the match five pairs.
    !
    ! With maxrot >= 180, the rotation is unknown (any roll): the grid
    ! of the whole circle would need thousands of nodes (the step moves
    ! the stars at the edge by a bin), so the candidates of the rotation
    ! and the offset are found first by the pairs of the brightest stars
    ! (rollvote: invariant to the rotation and the offset), and the vote
    ! is done around each of them; ntrial is the number of the nodes and
    ! the bins of the whole search (the circle, the scales, the offsets).
    ! reflboth (optional): both reflections are tried (the trials doubled).

    use astrotrafo

    character(len=*), intent(in) :: type
    real(REAL64), dimension(:), intent(in) :: alpha,delta,rflux,xpix,ypix,flux
    real(REAL64), intent(in) :: acen,dcen,xcen,ycen,scale,rot,refl, &
         maxoff0,maxrot,maxscale,bin
    integer, dimension(:), intent(in) :: frame
    integer, dimension(:), allocatable, intent(out) :: i1,i2
    real(REAL64), intent(out) :: ntrial
    logical, intent(out) :: status
    logical, intent(in) :: verbose
    logical, intent(in), optional :: reflboth

    real(REAL64), parameter :: rad = 57.295779513082322865_REAL64
    real(REAL64), parameter :: pi = 3.14159265358979323846_REAL64
    ! the frame stars voting (the brightest), the catalogue stars per
    ! a frame star (scaled by the area of the margin)
    integer, parameter :: maxvote = 60, minvote = 3, minpairs = 5
    real(REAL64), parameter :: catfac = 3
    ! any rotation: the stars of the pairs (the brightest frame stars,
    ! the catalogue stars at most), the candidates per reflection
    integer, parameter :: maxroll = 30, maxrollcat = 200, maxcand = 8
    real(REAL64), dimension(:), allocatable :: u,v,xc,yc,xf,yf,xs,ys,xr,yr, &
         crot,cxoff,cyoff
    integer, dimension(:), allocatable :: icat, ifr, jcat, jfr, ksel, cvote
    logical, dimension(:), allocatable :: used
    type(AstrometryProjection) :: tp
    type(AstrometryGeometry) :: tg
    real(REAL64) :: hw, hh, rmax, margin, drot, dscl, r, s, c, sn, d, &
         cell, tol, xoff, yoff, smin, rcat, ratio, lmin, dth, dt, hrot, hoff, &
         prot, pscl, pxoff, pyoff, qrot, qscl, qxoff, qyoff, pr, brefl, maxoff
    integer :: ncat, nfr, nd, nc, nb, nrot, nscl, i, j, k, l, m, ip, np, &
         best, pbest, qbest, npair, nh, nch, ncand, ic, nr
    integer(INT64) :: t0, t1, t2, trate
    logical :: anyrot, near

    status = .false.
    ntrial = 1
    allocate(i1(0),i2(0))

    ncat = size(alpha)
    nfr = size(xpix)
    if( ncat < minpairs .or. nfr < minpairs ) return

    anyrot = maxrot >= 180
    np = 1
    if( present(reflboth) ) then
       if( reflboth ) np = 2
    end if

    ! the catalogue projected (the prediction on the frame by a reflection
    ! below, relative to the centre)
    allocate(u(ncat),v(ncat),xc(ncat),yc(ncat))
    call tp % init(type,acen,dcen)
    call tp % proj(alpha,delta,u,v)

    ! the half-sizes of the frame around the centre, the margin
    hw = max(xcen - 0.5_REAL64, frame(1) + 0.5_REAL64 - xcen)
    hh = max(ycen - 0.5_REAL64, frame(2) + 0.5_REAL64 - ycen)
    rmax = hypot(hw,hh)
    maxoff = min(maxoff0, 2*rmax)

    ! the brightest frame stars, the catalogue by brightness
    nd = min(nfr,maxvote)
    allocate(ifr(nfr),icat(ncat),ksel(ncat))
    call fluxorder(flux,ifr)
    call fluxorder(rflux,icat)
    allocate(xf(nd),yf(nd))
    xf = xpix(ifr(1:nd)) - xcen
    yf = ypix(ifr(1:nd)) - ycen

    ! the grid of rotations and scales: the steps move the stars at
    ! the edge of the frame by a bin at most
    drot = min(maxrot, max(rad*bin/rmax, 0.01_REAL64))
    dscl = min(maxscale, max(bin/rmax, 1e-4_REAL64))
    nrot = 0
    nscl = 0
    if( anyrot ) then
       ! the whole circle
       nrot = ceiling(180/drot)
    else if( maxrot > 0 ) then
       nrot = ceiling(maxrot/drot)
    end if
    if( maxscale > 0 ) nscl = ceiling(maxscale/dscl)

    ! the histogram of the offsets in the cells of half of the bin,
    ! the votes are counted in the windows of 2 x 2 cells
    cell = bin / 2
    nb = 2*ceiling(maxoff/cell) + 2


    if( anyrot ) then
       ! the catalogue stars within the circle of the frame (any rotation),
       ! the offset and the scale
       margin = 0
       rcat = rmax*(1 + maxscale) + maxoff
       ratio = pi*rcat**2/(4*hw*hh)
       ! the pairs (rollvote): the shortest of 15% of the frame (the precision
       ! of the angle), the cells of the rotation and the offset by it
       nh = min(nfr,maxroll)
       lmin = max(8*bin, 0.15_REAL64*2*min(hw,hh))
       dth = max(0.5_REAL64/rad, 2*bin/lmin)
       dt = max(2*bin, dth*rcat/2)
       ! the vote around a candidate: its window and the steps of the grid
       hrot = rad*dth + drot
       hoff = dt + bin
    else
       margin = maxoff + rmax*(maxrot/rad + maxscale)
       ratio = ((hw+margin)*(hh+margin))/(hw*hh)
    end if

    best = 0
    brefl = refl
    prot = 0
    pscl = 1
    pxoff = 0
    pyoff = 0
    do ip = 1, np

       pr = refl
       if( ip == 2 ) pr = -refl
       call tg % init(scale,rot,pr)
       call tg % uv2pix(u,v,xc,yc)

       ! the brightest catalogue stars in the frame with the margin
       m = 0
       do i = 1, ncat
          k = icat(i)
          if( anyrot ) then
             near = hypot(xc(k),yc(k)) < rcat
          else
             near = abs(xc(k)) < hw + margin .and. abs(yc(k)) < hh + margin
          end if
          if( near ) then
             m = m + 1
             ksel(m) = k
          end if
       end do
       nc = min(m, max(30, nint(catfac*nd*ratio)))
       if( nc < minvote ) cycle
       ! the trials: the nodes of the grid and the independent bins
       ! (of the whole search, the reflections)
       ntrial = real(2*nrot+1,REAL64)*real(2*nscl+1,REAL64)*(nb/2.0_REAL64)**2*np
       xs = xc(ksel(1:nc))
       ys = yc(ksel(1:nc))

       if( .not. anyrot ) then

          call offvote(xf,yf,xs,ys,0.0_REAL64,nrot,drot,nscl,dscl, &
               0.0_REAL64,0.0_REAL64,maxoff,cell,pbest,prot,pscl,pxoff,pyoff)

       else

          ! the candidates by the pairs, the vote around each of them
          nch = min(nc, max(20, nint(catfac*nh*ratio)), maxrollcat)
          call system_clock(t0,trate)
          call rollvote(xf(1:nh),yf(1:nh),xs(1:nch),ys(1:nch),lmin,2*rmax, &
               maxscale,bin,maxoff,dth,dt,maxcand,ncand,crot,cxoff,cyoff,cvote)
          call system_clock(t1)
          pbest = 0
          do ic = 1, ncand
             nr = ceiling(hrot/drot)
             call offvote(xf,yf,xs,ys,crot(ic),nr,drot,nscl,dscl, &
                  cxoff(ic),cyoff(ic),hoff,cell,qbest,qrot,qscl,qxoff,qyoff)
             if( verbose ) write(error_unit,'(a,i0,a,f0.1,a,f7.1,f7.1,a,i0,a,i0,a,f0.2,a,f0.4)') &
                  ' Any rotation: candidate #',ic,' rotation ',crot(ic),' deg, offset ', &
                  cxoff(ic),cyoff(ic),' pix (',cvote(ic),' pairs): vote ',qbest, &
                  ' at ',modulo(qrot,360.0_REAL64),' deg, scale ',qscl
             if( qbest > pbest ) then
                pbest = qbest
                prot = qrot
                pscl = qscl
                pxoff = qxoff
                pyoff = qyoff
             end if
          end do
          call system_clock(t2)
          if( verbose ) write(error_unit,'(a,i0,a,i0,a,i0,a,f0.2,a,f0.1,a,f0.1,a,f0.0,a,f0.3,a,f0.3,a)') &
               ' Any rotation (reflection ',nint(pr),'): ',nh,' frame stars, ',nch, &
               ' catalogue stars in pairs of ',lmin,' pix at least, cells ', &
               rad*dth,' deg, ',dt,' pix within ',maxoff,' pix; the pairs ', &
               real(t1-t0,REAL64)/trate,' s, the votes ',real(t2-t1,REAL64)/trate,' s.'

       end if

       if( pbest > best ) then
          best = pbest
          r = prot
          s = pscl
          xoff = pxoff
          yoff = pyoff
          brefl = pr
       end if
    end do

    if( verbose ) write(error_unit,'(a,i0,a,i0,a,i0,a,f0.2,a,f0.4,a,es9.2e3,a)') &
         ' WCS-guided match: ',nd,' frame stars, ',nc,' catalogue stars, best vote ', &
         best,' at rotation ',r,' deg, scale ',s,' (',ntrial,' trials).'
    if( best < minvote ) return

    ! the pairs of the nearest stars by the transformation within the bin
    ! and the error of the grid at the edge, all stars, the brighter first
    if( np > 1 ) then
       call tg % init(scale,rot,brefl)
       call tg % uv2pix(u,v,xc,yc)
    end if
    r = r / rad
    c = cos(r)
    sn = sin(r)
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


  subroutine offvote(xf,yf,xc,yc,rot0,nrot,drot,nscl,dscl,x0,y0,maxoff,cell, &
       best,rbest,sbest,xoff,yoff)

    ! The vote for the offset: the catalogue stars xc,yc (predicted on
    ! the frame, relative to its centre) rotated by rot0 + ir*drot [deg]
    ! (ir = -nrot .. nrot) and scaled by 1 + is*dscl (is = -nscl .. nscl);
    ! for every node, the differences of the frame stars xf,yf and the
    ! catalogue stars within maxoff of x0,y0 [pix] are counted in the cells
    ! of the size cell [pix], the votes in the windows of 2 x 2 cells.
    ! On output, the most votes best (the first node and window of them),
    ! the rotation rbest [deg], the scale sbest and the offset xoff, yoff
    ! [pix]: the mean of the differences in the window. Only the cells
    ! voted are visited (and cleared), the windows of the most votes of
    ! a node are the first ones in the order of the cells (k, then l),
    ! as by the scan of the whole histogram.

    real(REAL64), dimension(:), intent(in) :: xf,yf,xc,yc
    real(REAL64), intent(in) :: rot0,drot,dscl,x0,y0,maxoff,cell
    integer, intent(in) :: nrot,nscl
    integer, intent(out) :: best
    real(REAL64), intent(out) :: rbest,sbest,xoff,yoff

    real(REAL64), parameter :: rad = 57.295779513082322865_REAL64
    real(REAL64), dimension(size(xc)) :: xr,yr
    integer, dimension(:,:), allocatable :: hist
    integer, dimension(:), allocatable :: tk, tl
    real(REAL64) :: r,c,sn,s,dx,dy
    integer :: nb,ir,is,i,j,k,l,m,ibest,jbest,irb,isb,nsup,nt,it,k0,l0, &
         nbest,kb,lb

    nb = 2*ceiling(maxoff/cell) + 2
    allocate(hist(nb,nb),tk(size(xf)*size(xc)),tl(size(xf)*size(xc)))
    hist = 0
    best = 0
    ibest = 0
    jbest = 0
    irb = 0
    isb = 0
    do ir = -nrot, nrot
       r = (rot0 + ir*drot) / rad
       c = cos(r)
       sn = sin(r)
       do is = -nscl, nscl
          s = 1 + is*dscl
          xr = s*(c*xc - sn*yc)
          yr = s*(sn*xc + c*yc)
          nt = 0
          do i = 1, size(xf)
             do j = 1, size(xc)
                dx = xf(i) - xr(j) - x0
                dy = yf(i) - yr(j) - y0
                if( abs(dx) < maxoff .and. abs(dy) < maxoff ) then
                   k = int((dx + maxoff)/cell) + 1
                   l = int((dy + maxoff)/cell) + 1
                   if( hist(k,l) == 0 ) then
                      nt = nt + 1
                      tk(nt) = k
                      tl(nt) = l
                   end if
                   hist(k,l) = hist(k,l) + 1
                end if
             end do
          end do
          ! the windows of 2 x 2 cells with a voted cell
          nbest = 0
          kb = 0
          lb = 0
          do it = 1, nt
             do k0 = max(tk(it) - 1, 1), min(tk(it), nb - 1)
                do l0 = max(tl(it) - 1, 1), min(tl(it), nb - 1)
                   m = hist(k0,l0) + hist(k0+1,l0) + hist(k0,l0+1) + hist(k0+1,l0+1)
                   if( m > nbest .or. (m == nbest .and. (k0 < kb .or. &
                        (k0 == kb .and. l0 < lb))) ) then
                      nbest = m
                      kb = k0
                      lb = l0
                   end if
                end do
             end do
          end do
          do it = 1, nt
             hist(tk(it),tl(it)) = 0
          end do
          if( nbest > best ) then
             best = nbest
             ibest = kb
             jbest = lb
             irb = ir
             isb = is
          end if
       end do
    end do
    rbest = rot0 + irb*drot
    sbest = 1 + isb*dscl
    xoff = x0
    yoff = y0
    if( best == 0 ) return

    ! the offset: the mean of the differences in the best window
    r = rbest / rad
    c = cos(r)
    sn = sin(r)
    s = sbest
    xr = s*(c*xc - sn*yc)
    yr = s*(sn*xc + c*yc)
    xoff = 0
    yoff = 0
    nsup = 0
    do i = 1, size(xf)
       do j = 1, size(xc)
          dx = xf(i) - xr(j) - x0
          dy = yf(i) - yr(j) - y0
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
    xoff = x0 + xoff / nsup
    yoff = y0 + yoff / nsup

  end subroutine offvote


  subroutine rollvote(xf,yf,xc,yc,lmin,lmax,maxscale,bin,maxoff,dth,dt,maxcand, &
       ncand,crot,cxoff,cyoff,cvote)

    ! The candidates of the rotation (any) and the offset of the catalogue
    ! stars xc,yc against the frame stars xf,yf (both relative to the centre
    ! of the frame, the catalogue predicted on it) [pix]. Every pair of the
    ! frame stars of the length lmin .. lmax and every pair of the catalogue
    ! stars of the same length (within the relative scale maxscale and bin)
    ! give the rotation (both orientations of the catalogue pair) and the
    ! offset (by their midpoints and the ratio of the lengths); the votes
    ! within maxoff are counted in the cells of dth [rad] and dt [pix],
    ! in the windows of 2 x 2 x 2 cells. The true pairs of stars vote for
    ! one transformation, the others are spread over the circle and the
    ! offsets. The windows of the most votes (maxcand at most, not
    ! overlapping, three votes at least) are the candidates: the rotation
    ! crot [deg], the offset cxoff, cyoff [pix] (the centres of the windows)
    ! and the votes cvote.

    real(REAL64), dimension(:), intent(in) :: xf,yf,xc,yc
    real(REAL64), intent(in) :: lmin,lmax,maxscale,bin,maxoff,dth,dt
    integer, intent(in) :: maxcand
    integer, intent(out) :: ncand
    real(REAL64), dimension(:), allocatable, intent(out) :: crot,cxoff,cyoff
    integer, dimension(:), allocatable, intent(out) :: cvote

    real(REAL64), parameter :: rad = 57.295779513082322865_REAL64
    real(REAL64), parameter :: pi = 3.14159265358979323846_REAL64
    integer, parameter :: minvote = 3
    real(REAL64), dimension(:), allocatable :: pd,pa,pmx,pmy,w
    integer, dimension(:), allocatable :: idx, ucnt, wsum
    integer(INT64), dimension(:), allocatable :: key, ukey, tkey
    logical, dimension(:), allocatable :: free
    integer :: nf,nc,npc,i,j,k,l,q,io,nth,nbx,ith,ix,iy,nv,nu,a,b,e,da,db,de, &
         ibest,kbest,n,lo,hi,mid
    integer(INT64) :: kk
    real(REAL64) :: dxf,dyf,df,af,mxf,myf,dlo,dhi,s,th,c,sn,tx,ty,dthc

    nf = size(xf)
    nc = size(xc)
    ncand = 0
    allocate(crot(maxcand),cxoff(maxcand),cyoff(maxcand),cvote(maxcand))

    ! the catalogue pairs (of the lengths possible) sorted by the length
    npc = nc*(nc - 1)/2
    allocate(pd(npc),pa(npc),pmx(npc),pmy(npc),w(npc),idx(npc))
    npc = 0
    do k = 1, nc - 1
       do l = k + 1, nc
          df = hypot(xc(l) - xc(k), yc(l) - yc(k))
          if( df*(1 + maxscale) + bin >= lmin .and. df*(1 - maxscale) - bin <= lmax ) then
             npc = npc + 1
             w(npc) = df
             pa(npc) = atan2(yc(l) - yc(k), xc(l) - xc(k))
             pmx(npc) = (xc(k) + xc(l)) / 2
             pmy(npc) = (yc(k) + yc(l)) / 2
          end if
       end do
    end do
    if( npc == 0 ) return
    call realorder(w(1:npc),idx(1:npc))
    pd(1:npc) = w(idx(1:npc))
    w(1:npc) = pa(idx(1:npc))
    pa(1:npc) = w(1:npc)
    w(1:npc) = pmx(idx(1:npc))
    pmx(1:npc) = w(1:npc)
    w(1:npc) = pmy(idx(1:npc))
    pmy(1:npc) = w(1:npc)

    ! the cells: the rotation over the circle, the offsets within maxoff
    nth = max(1,nint(2*pi/dth))
    dthc = 2*pi / nth
    nbx = max(1,ceiling(2*maxoff/dt))
    allocate(key(1024))
    nv = 0
    do i = 1, nf - 1
       do j = i + 1, nf
          dxf = xf(j) - xf(i)
          dyf = yf(j) - yf(i)
          df = hypot(dxf,dyf)
          if( df < lmin .or. df > lmax ) cycle
          af = atan2(dyf,dxf)
          mxf = (xf(i) + xf(j)) / 2
          myf = (yf(i) + yf(j)) / 2
          dlo = (df - bin) / (1 + maxscale)
          dhi = (df + bin) / (1 - maxscale)
          ! the first catalogue pair not shorter than dlo (bisection)
          lo = 1
          hi = npc + 1
          do while( lo < hi )
             mid = (lo + hi) / 2
             if( pd(mid) < dlo ) then
                lo = mid + 1
             else
                hi = mid
             end if
          end do
          do q = lo, npc
             if( pd(q) > dhi ) exit
             s = df / pd(q)
             do io = 1, 2
                ! the orientations of the catalogue pair: k -> l, l -> k
                th = af - pa(q)
                if( io == 2 ) th = th + pi
                c = cos(th)
                sn = sin(th)
                tx = mxf - s*(c*pmx(q) - sn*pmy(q))
                ty = myf - s*(sn*pmx(q) + c*pmy(q))
                if( abs(tx) < maxoff .and. abs(ty) < maxoff ) then
                   ith = modulo(floor(th / dthc), nth)
                   ix = min(int((tx + maxoff)/dt), nbx - 1)
                   iy = min(int((ty + maxoff)/dt), nbx - 1)
                   nv = nv + 1
                   if( nv > size(key) ) then
                      allocate(tkey(2*size(key)))
                      tkey(1:size(key)) = key
                      call move_alloc(tkey,key)
                   end if
                   key(nv) = (int(ith,INT64)*nbx + ix)*nbx + iy
                end if
             end do
          end do
       end do
    end do
    if( nv == 0 ) return

    ! the cells voted and their votes
    call int64sort(key(1:nv))
    allocate(ukey(nv),ucnt(nv))
    nu = 0
    do i = 1, nv
       if( nu > 0 ) then
          if( key(i) == ukey(nu) ) then
             ucnt(nu) = ucnt(nu) + 1
             cycle
          end if
       end if
       nu = nu + 1
       ukey(nu) = key(i)
       ucnt(nu) = 1
    end do

    ! the windows of 2 x 2 x 2 cells starting by the voted cells
    ! (the rotation is periodic)
    allocate(wsum(nu),free(nu))
    do i = 1, nu
       call decode(ukey(i),a,b,e)
       n = 0
       do da = 0, 1
          do db = 0, 1
             do de = 0, 1
                if( b + db >= nbx .or. e + de >= nbx ) cycle
                kk = (int(modulo(a + da,nth),INT64)*nbx + (b + db))*nbx + (e + de)
                n = n + votes(kk)
             end do
          end do
       end do
       wsum(i) = n
    end do

    ! the candidates: the windows of the most votes, not overlapping
    ! (the first one of the same votes: the lowest cell)
    free = .true.
    do ncand = 0, maxcand - 1
       ibest = 0
       kbest = minvote - 1
       do i = 1, nu
          if( free(i) .and. wsum(i) > kbest ) then
             kbest = wsum(i)
             ibest = i
          end if
       end do
       if( ibest == 0 ) exit
       call decode(ukey(ibest),a,b,e)
       crot(ncand+1) = rad*(a + 1)*dthc
       cxoff(ncand+1) = -maxoff + (b + 1)*dt
       cyoff(ncand+1) = -maxoff + (e + 1)*dt
       cvote(ncand+1) = kbest
       ! the windows overlapping it are excluded
       do i = 1, nu
          if( .not. free(i) ) cycle
          call decode(ukey(i),j,k,l)
          if( min(modulo(j - a,nth),modulo(a - j,nth)) <= 1 .and. &
               abs(k - b) <= 1 .and. abs(l - e) <= 1 ) free(i) = .false.
       end do
    end do

  contains

    subroutine decode(kk,a,b,e)
      integer(INT64), intent(in) :: kk
      integer, intent(out) :: a,b,e
      e = int(modulo(kk,int(nbx,INT64)))
      b = int(modulo(kk / nbx,int(nbx,INT64)))
      a = int(kk / (int(nbx,INT64)*nbx))
    end subroutine decode

    integer function votes(kk)
      ! the votes of the cell kk (bisection in the sorted cells)
      integer(INT64), intent(in) :: kk
      integer :: lo,hi,mid
      votes = 0
      lo = 1
      hi = nu
      do while( lo <= hi )
         mid = (lo + hi) / 2
         if( ukey(mid) == kk ) then
            votes = ucnt(mid)
            return
         else if( ukey(mid) < kk ) then
            lo = mid + 1
         else
            hi = mid - 1
         end if
      end do
    end function votes

  end subroutine rollvote


  subroutine realorder(f,idx)

    ! the indices of f in the ascending order, the equal values in their
    ! original order (heapsort by the pair of the value and the index)

    real(REAL64), dimension(:), intent(in) :: f
    integer, dimension(:), intent(out) :: idx
    integer :: n, i, j, k, t, root, child

    n = size(f)
    idx = [(i, i = 1, n)]
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
      integer, intent(in) :: a,b
      after = f(a) > f(b) .or. (f(a) == f(b) .and. a > b)
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

  end subroutine realorder


  subroutine int64sort(a)

    ! the ascending order of a (heapsort in place)

    integer(INT64), dimension(:), intent(in out) :: a
    integer(INT64) :: t
    integer :: n, j, k, root, child

    n = size(a)
    do k = n/2, 1, -1
       call sift(k,n)
    end do
    do j = n, 2, -1
       t = a(1)
       a(1) = a(j)
       a(j) = t
       call sift(1,j-1)
    end do

  contains

    subroutine sift(start,last)
      integer, intent(in) :: start,last
      root = start
      do
         child = 2*root
         if( child > last ) exit
         if( child < last ) then
            if( a(child+1) > a(child) ) child = child + 1
         end if
         if( a(child) > a(root) ) then
            t = a(root)
            a(root) = a(child)
            a(child) = t
            root = child
         else
            exit
         end if
      end do
    end subroutine sift

  end subroutine int64sort


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
