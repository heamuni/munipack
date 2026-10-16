!
!  match    finding of cross identification between two lists
!           the core of the matching algorithm
!
!  Copyright © 2010 - 2026 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Munipack.
!
!  Munipack is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Munipack is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
!
!
!  http://en.wikipedia.org/wiki/Backtracking
!
!  Key ideas of matching by backtracking has origin in book:
!    N.Wirth: Data Structures + Algorithms = Programs
!

module matcher

  use iso_fortran_env
  implicit none
  private

  type, private :: dataseq
     real(REAL64), dimension(:), pointer :: x,y,f
     real(REAL64), dimension(:,:,:,:), allocatable, private :: tris
     real(REAL64), dimension(:,:), allocatable, private :: dist
     logical, dimension(:,:,:), allocatable, private :: tris_def
     logical, dimension(:,:), allocatable, private :: dist_def
   contains
     procedure :: init => dataseq_init
     procedure :: enorm => dataseq_enorm
     procedure :: triangle => dataseq_triangle
     final :: dataseq_clean
  end type dataseq


  type, private :: sequences
     real(REAL64) :: sig1 = 1.0/3600.0 ! spherical coordinates limit [deg]
     real(REAL64) :: sig2 = 1.0        ! on chip coordinates limit [pix]
     real(REAL64) :: flxmin = 0.5    ! minimal flux ratio
     real(REAL64) :: flxmax = 2.0    ! maximal flux ratio
     real(REAL64) :: amin, dmin, tscale, tflux
     real(REAL64) :: dtol1,dtol2
     integer :: minmatch = 5
     integer :: maxmatch = 42
     integer :: minseq = 5 + 3
     logical :: reciprocal             ! is the flux ratio reciprocal?
     logical :: fluxtest = .true.  ! use = .true., or don't use flux
     integer :: nseq = 0  ! current sequence length
     integer :: nid = 0   ! current lenght of ?
     logical :: found = .false.
     integer, dimension(:), allocatable :: id1, id2, seq1, seq2
     logical, dimension(:), allocatable :: idx1, idx2
     real(REAL64), dimension(:), allocatable :: qchi2
     real(REAL64), dimension(:), allocatable :: rflux, scale
   contains
     procedure :: init => sequences_init
     final :: sequences_clean
  end type sequences

  logical, parameter :: debug = .false.
  real(REAL64), parameter :: xsig = 2.0 ! Normal distribution limit at 95%
  real(REAL64), parameter :: rsig = 3.8 ! Rayleigh distribution limit at 95%

  ! the parallel match of the chunk pairs (astromatch): the first chunk
  ! pair which succeeded, the matches of the later ones are abandoned
  integer, public :: matchfound = huge(1)

  public :: match, triangle

contains

  subroutine match(x1,y1,f1,index1,x2,y2,f2,index2,minmatch,maxmatch,sigma1,sigma2, &
       flxmin,flxmax,fluxtest,id1,id2,luckymatch,verbose,status,scale0,scaletol, &
       ntrial,ichunk)

    ! scale0 (optional, > 0) is an expected scale (the frame units per
    ! the reference units: pix/deg) known within the relative tolerance
    ! scaletol: the starting pairs of the sequences whose ratio of
    ! lengths is far from it are skipped (the scale prior);
    ! ntrial (optional) is the number of the starting pairs tried
    ! (the transformations tested, for the significance of the match);
    ! ichunk (optional) is the index of the chunk pair of the parallel
    ! match: it stops (without a match) when a chunk pair before it
    ! has succeeded (matchfound)

    implicit none

!!$  interface
!!$     subroutine matchprint(id1,id2,c,t,p,d,q)
!!$       import :: REAL64
!!$       integer, dimension(:), intent(in) :: id1,id2
!!$       real(REAL64), intent(in) :: c,t,p,d,q
!!$     end subroutine matchprint
!!$
!!$     subroutine progress2(ns,ntot,amin,dmin,id1,id2)
!!$       import :: REAL64
!!$       integer, intent(in) :: ns,ntot
!!$       integer, dimension(:), intent(in), optional :: id1,id2
!!$       real(REAL64), intent(in), optional :: amin,dmin
!!$     end subroutine progress2
!!$
!!$  end interface

    real(REAL64), dimension(:), intent(in), target :: x1,y1,x2,y2,f1,f2
    logical, dimension(:), intent(in) :: index1, index2
    integer, dimension(:), allocatable, intent(out) :: id1, id2
    logical, intent(in) :: fluxtest
    integer, intent(in) :: minmatch, maxmatch, luckymatch
    real(REAL64), intent(in) :: sigma1, sigma2, flxmin, flxmax
    logical, intent(in) :: verbose
    logical, intent(out) :: status
    real(REAL64), intent(in), optional :: scale0, scaletol
    real(REAL64), intent(out), optional :: ntrial
    integer, intent(in), optional :: ichunk

    integer :: i1,i2,j1,j2,ns,ntot,n1,n2,jfound
    real(REAL64) :: r,d1,d2,s0,stol,sabs,ntried
    type(sequences) :: seq
    type(dataseq) :: ref, lst

    status = .false.

    if( size(x1) /= size(y1) .or. size(x2) /= size(y2) ) &
         stop 'Match: Array dimensions does not matches.'

    if( fluxtest ) then
       if( (.not. all(f1 > 0)) .or. (.not. all(f2 > 0)) ) &
            stop 'Match: All fluxes must be > 0.'
    end if

    if( .not. (minmatch > 0 .and. maxmatch > minmatch .and. &
         sigma1 > 0 .and. sigma2 > 0) ) &
         stop 'Match improper parameters? - minmatch, maxmatch, sigma1, sigma2'

    seq%minmatch = minmatch
    seq%maxmatch = maxmatch
    seq%sig1 = sigma1
    seq%sig2 = sigma2
    seq%fluxtest = fluxtest
    seq%flxmin = flxmin
    seq%flxmax = flxmax

    call ref%init(x1,y1,f1,maxmatch)
    call lst%init(x2,y2,f2,maxmatch)
    call seq%init(maxmatch)

    seq%nseq = 0
    seq%nid = 0
    seq%amin = huge(seq%amin)
    seq%dmin = huge(seq%dmin)
    seq%dtol1 = xsig * seq%sig1
    seq%dtol2 = xsig * seq%sig2
    seq%minseq = minmatch + luckymatch

    ! the scale prior: the lengths of a starting pair (frame d2,
    ! reference d1) must agree with scale0 within the relative tolerance
    ! and the errors of both ends (Normal limit of a difference)
    s0 = -1
    stol = 0
    sabs = 0
    if( present(scale0) .and. present(scaletol) ) then
       if( scale0 > 0 .and. scaletol > 0 ) then
          s0 = scale0
          stol = scaletol
          sabs = xsig*sqrt(2.0_REAL64)*max(sigma2,scale0*sigma1)
       end if
    end if

    ns = 0
    ntried = 0
    n1 = size(x1)
    n2 = size(x2)
    ntot = (n1 + 1) * (n1 / 2)


!!$    block
!!$      integer :: i1,i2,j1,j2
!!$      real :: r1,r2
!!$      open(1,file='/var/tmp/xxx1')
!!$      do i1 = 1, size(f1)
!!$         do i2 = i1,size(f1)
!!$            do j1 = 1, size(f2)
!!$               do j2 = j1,size(f2)
!!$                  r1 = f1(i1)/f2(j1)
!!$                  r2 = f1(i2)/f2(j2)
!!$                  if( 0.5 < r1/r2 .and. r1/r2 < 1.5 ) then
!!$                     write(1,*) f1(i1)/f2(j1),f1(i2)/f2(j2)
!!$                  end if
!!$               end do
!!$            end do
!!$         end do
!!$      end do
!!$      close(1)
!!$      stop
!!$    end block


    do i1 = 1, n1 - 1
       do i2 = i1 + 1, n1

          ns = ns + 1
          !  call progress2(ns,ntot)
          !  if( debug .and. mod(ns,100) == 0 ) &
          !     write(error_unit,'(a,2i7,f7.1,a)') 'Pair ',i1,i2,(1e2*ns)/ntot,'%'
          !  end if

          if( index1(i1) .and. index1(i2) ) then
             seq%id1(1:2) = [ i1, i2 ]

             d1 = ref%enorm(i1,i2)
             seq%idx1 = index1
             seq%idx1(seq%id1(1:2)) = .false.

             do j1 = 1, n2 - 1
                do j2 = j1+1, n2

                   if( index2(j1) .and. index2(j2) ) then
                      seq%id2(1:2) = [ j1, j2 ]

                      seq%nid = 2
                      d2 = lst%enorm(j1,j2)
                      if( s0 > 0 ) then
                         if( abs(d2 - s0*d1) > stol*s0*d1 + sabs ) cycle
                      end if
                      seq%idx2 = index2
                      seq%idx2(seq%id2(1:2)) = .false.

                      if( seq%fluxtest ) then
                         seq%rflux(1:2) = ref%f(seq%id1(1:2)) / lst%f(seq%id2(1:2))
                         seq%reciprocal = sum(seq%rflux(1:2)) < 2
                         if( seq%reciprocal ) seq%rflux(1:2) = 1 / seq%rflux(1:2)
                         r = seq%rflux(2) / seq%rflux(1)
!                         if( .not. (seq%flxmin < r .and. r < seq%flxmax) ) write(*,*) seq%flxmin,r,seq%flxmax
                         if( .not. (seq%flxmin < r .and. r < seq%flxmax) ) cycle
                      end if

                      seq%scale(1:2) = d2 / d1

                      ! pre-select stars by similar flux
!!$                      if( seq%fluxtest ) then
!!$                         block
!!$                           real(REAL64) :: rf, flxmin, flxmax,rflux
!!$                           integer :: i,j
!!$                           logical :: xflux
!!$                           open(1,file='/var/tmp/xxx')
!!$                           rf = sum(seq%rflux(1:2)) / 2
!!$                           flxmin = seq%flxmin * rf
!!$                           flxmax = seq%flxmax * rf
!!$                           do j = 1, size(seq%idx2)
!!$                              xflux = .false.
!!$                              do i = 1, size(seq%idx1)
!!$                                 if( seq%reciprocal ) then
!!$                                    rflux = lst%f(j) / ref%f(i)
!!$                                 else
!!$                                    rflux = ref%f(i) / lst%f(j)
!!$                                 end if
!!$                                 write(1,*) i,j, seq%rflux(1:2) ,rflux / rf,ref%f(i),lst%f(j), flxmin, flxmax
!!$                                 if( flxmin <= rflux .and. rflux <= flxmax ) then
!!$                                    xflux = .true.
!!$ !                                   write(*,*) i,j,rflux
!!$                                    exit
!!$                                 end if
!!$                              end do
!!$                              if( xflux ) then
!!$                                 seq%idx2(j) = .true.
!!$                              end if
!!$                           end do
!!$                           close(1)
!!$!                           write(*,*) count(seq%idx2)
!!$!                           stop
!!$                         end block
!!$                      end if

                      ! a chunk pair before this one has the match
                      if( present(ichunk) ) then
                         !$omp atomic read
                         jfound = matchfound
                         if( jfound < ichunk ) then
                            seq%nseq = 0
                            goto 666
                         end if
                      end if

                      ! start of recursive search
                      ntried = ntried + 1
                      seq%found = .false.
                      call aqualung(seq,ref,lst)

                      ! show the longest sequence
                      ! if( seq%nid >= max(seq%minmatch,seq%nseq) ) &
                      !     call progress2(ns,ntot,seq%amin,seq%dmin, &
                      !     seq%seq1(1:seq%nseq),seq%seq2(1:seq%nseq))
                      if( seq%found .and. verbose .and. seq%nseq > 3 ) &
                           call progress2(seq%seq1(1:seq%nseq),seq%seq2(1:seq%nseq),&
                           seq%tscale,seq%tflux,seq%amin,seq%dmin,seq%qchi2(seq%nseq))

                      ! finish on first success attempt (luckymatch == 0),
                      ! or finish when results looks more reliable (>0)
                      if( seq%nseq >= seq%minseq ) goto 666

                   endif
                end do
             end do
          end if
       end do
    end do

666 continue

    if( present(ntrial) ) ntrial = ntried

    if( seq%nseq >= minmatch ) then

       allocate(id1(seq%nseq),id2(seq%nseq))
       id1 = seq%seq1(1:seq%nseq)
       id2 = seq%seq2(1:seq%nseq)
       status = .true.

    end if

    if( debug ) then
       block
         integer :: i,j,k,i0,j0
         real(REAL64) :: u1,v1,u2,v2,stol1,stol2,d1,d2

         open(1,file='/var/tmp/uv1')
         do i = 1, maxmatch
            do j = 1, maxmatch
               do k = 1, maxmatch
                  if( ref%tris_def(i,j,k) ) write(1,*) ref%tris(i,j,k,1:2)
               end do
            end do
         end do
         close(1)

         open(1,file='/var/tmp/uv2')
         do i = 1, maxmatch
            do j = 1, maxmatch
               do k = 1, maxmatch
                  if( lst%tris_def(i,j,k) ) write(1,*) lst%tris(i,j,k,1:2)
               end do
            end do
         end do
         close(1)

         open(1,file='/var/tmp/duv')
         if( seq%nseq > 0 ) then
            i0 = seq%seq1(seq%nseq)
            j0 = seq%seq2(seq%nseq)
         end if
         do i = 1, seq%nseq-1
            j = seq%seq1(i)
            k = seq%seq1(i+1)
            if( ref%tris_def(i0,j,k) ) then
               u1 = ref%tris(i0,j,k,1)
               v1 = ref%tris(i0,j,k,2)
               d1 = ref%tris(i0,j,k,3)
               stol1 = seq%sig1 / ref%tris(i0,j,k,3)
            else
               u1 = -1
            end if
            j = seq%seq2(i)
            k = seq%seq2(i+1)
            if( lst%tris_def(i0,j,k) ) then
               u2 = lst%tris(j0,j,k,1)
               v2 = lst%tris(j0,j,k,2)
               d2 = lst%tris(j0,j,k,3)
               stol2 = seq%sig2 / lst%tris(j0,j,k,3)
            else
               u2 = -1
            end if
            if( u1 > -0.5 .and. u2 > -0.5 ) &
                 write(1,*) u1-u2, v1 - v2, stol1, stol2,seq%sig1/d1,seq%sig2/d2
            i0 = seq%seq1(i)
            j0 = seq%seq2(i)
         end do
         close(1)

         open(1,file='/var/tmp/id2')
         do i = 1, seq%nseq
            j = seq%seq2(i)
            write(1,*) lst%x(j), lst%y(j), j
         end do
         close(1)

         block
           real(REAL64) :: r,r0,s
           open(1,file='/var/tmp/fr')
           s = 0
           do i = 1, seq%nseq
              r = ref%f(seq%seq1(i)) / lst%f(seq%seq2(i))
              r = 1 / r
              if( i == 1 ) r0 = r
              s = s + r/r0
              write(1,*) r,r/r0
              r0 = r
           end do
           close(1)
           if( seq%nseq > 0  ) write(*,*) 'mean r = ',s / seq%nseq
         end block

         block
           real(REAL64), dimension(seq%nseq) :: d1,d2

           call distances(ref%x,ref%y,seq%seq1(1:seq%nseq),d1)
           call distances(lst%x,lst%y,seq%seq2(1:seq%nseq),d2)

           open(1,file='/var/tmp/scale')
           do i = 1, seq%nseq
              write(1,*) d2(i) / d1(i), d2(i) - seq%scale(seq%nseq)*d1(i)
           end do
           close(1)
         end block

       end block
    end if

  end subroutine match


  recursive subroutine aqualung(seq,ref,lst)

    type(sequences), intent(in out) :: seq
    type(dataseq), intent(in out) :: ref, lst

!    integer, dimension(:), allocatable :: idx2
    integer :: i,j,i1,i2,j1,j2,i0,j0,n,nref,nlst,nid!,nidx2,jidx2
    real(REAL64) :: u1, v1, u2, v2, d1, d2, c1, c2, tol, &
         tol1, tol2, a, d, xlim, rflux, tflux, &
         rf,rsig1, rsig2, scale, sc, flxmin, flxmax, dtol
!    logical :: xflux

    nref = size(ref%x)
    nlst = size(lst%x)
    nid = seq%nid
    n = nid - 1
    i1 = seq%id1(n)
    i2 = seq%id1(nid)
    j1 = seq%id2(n)
    j2 = seq%id2(nid)
    rsig1 = xsig*seq%sig1
    rsig2 = xsig*seq%sig2

    xlim = seq%qchi2(nid)

    rf = sum(seq%rflux(2:nid)) / (nid - 1)
    flxmin = seq%flxmin * rf
    flxmax = seq%flxmax * rf

    scale = sum(seq%scale(2:nid))  / (nid - 1)
    dtol = rsig*max(seq%sig2,scale*seq%sig1)

!!$    allocate(idx2(size(seq%idx2)))
!!$    nidx2 = 0
!!$    do j = 1, size(seq%idx2)
!!$
!!$       if( seq%fluxtest ) then
!!$          xflux = .false.
!!$          do i = 1, nref
!!$             if( .not. seq%idx1(i) ) cycle
!!$             if( seq%reciprocal ) then
!!$                rflux = lst%f(j) / ref%f(i)
!!$             else
!!$                rflux = ref%f(i) / lst%f(j)
!!$             end if
!!$             if( flxmin <= rflux .and. rflux <= flxmax ) then
!!$                xflux = .true.
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          xflux = .true.
!!$       end if
!!$
!!$       if( seq%idx2(j) .and. xflux ) then
!!$!       if( seq%idx2(j) ) then
!!$          nidx2 = nidx2 + 1
!!$          idx2(nidx2) = j
!!$       end if
!!$    end do
!!$
!!$    write(*,*) nidx2

    do i = 1, nref
       if( .not. seq%idx1(i) ) cycle

       d1 = ref%enorm(i,i2)

       ! reject close stars including identical ones
       if( d1 < rsig1 ) cycle

       ! In front of launching of a long-duration search,
       ! we step over the stars already presented
       ! in previously formed sequences and those which
       ! distances are under coordinate uncertainities.

       call ref%triangle(i1,i2,i,u1,v1,c1)
       tol1 = rsig1 / c1

       ! rejects triangles with too short sides
       if( v1 < tol1 ) cycle

!       do jidx2 = 1, nidx2
!          j = idx2(jidx2)
       do j = 1, nlst
          if( .not. seq%idx2(j) ) cycle

          if( seq%fluxtest ) then
             if( seq%reciprocal ) then
                rflux = lst%f(j) / ref%f(i)
             else
                rflux = ref%f(i) / lst%f(j)
             end if
!             if( rflux < flxmin .or. rflux > flxmax ) write(*,*) rflux,flxmin,flxmax
             if( rflux < flxmin .or. rflux > flxmax ) cycle
          end if

          d2 = lst%enorm(j,j2)
          if( abs(d2 - scale*d1) < dtol ) then
             call lst%triangle(j1,j2,j,u2,v2,c2)
             tol2 = rsig2 / c2
          else
             cycle
          end if

          if( v2 < tol2 ) cycle

          tol = max(tol1,tol2)
          if( abs(u1 - u2) < tol .and. abs(v1 - v2) < tol ) then

             ! check head and tail of the sequence on similarity of triangles
             i0 = seq%id1(1)
             j0 = seq%id2(1)
             call ref%triangle(i0,i2,i,u1,v1,c1)
             call lst%triangle(j0,j2,j,u2,v2,c2)
             tol1 = rsig1 / c1
             tol2 = rsig2 / c2
             tol = max(tol1, tol2)

             seq%nid = seq%nid + 1
             seq%id1(seq%nid) = i
             seq%id2(seq%nid) = j
             seq%scale(seq%nid) = d2 / d1
             if( seq%fluxtest ) seq%rflux(seq%nid) = rflux

             ! final checks for realiable match:
             ! dispersions of distances and angles
             sc = sum(seq%scale(2:seq%nid)) / (seq%nid - 1)
             d = sdisp(sc,ref%x,ref%y,seq%id1(1:seq%nid),seq%sig1,&
                  lst%x,lst%y,seq%id2(1:seq%nid),seq%sig2)
             a = angels(ref%x,ref%y,seq%id1(1:seq%nid),seq%sig1,&
                  lst%x,lst%y,seq%id2(1:seq%nid),seq%sig2)
             if ( seq%fluxtest ) then
                tflux = sum(seq%rflux(2:seq%nid)) / (seq%nid - 1)
             else
                tflux = 1
             end if

             if( abs(u1 - u2) < tol .and. abs(v1 - v2) < tol &
                  .and. v2 > tol2 .and. v1 > tol1 &
                  .and. a < xlim .and. d < xlim ) then

                if( seq%nid < seq%maxmatch .and. seq%nseq <= seq%minseq ) then
!                if( seq%nid < seq%maxmatch ) then

                   seq%idx1(i) = .false.
                   seq%idx2(j) = .false.

                   ! sink "with deep-sea-diver sounds"
                   call aqualung(seq,ref,lst)

                   seq%idx1(i) = .true.
                   seq%idx2(j) = .true.
                end if

                if ( seq%nid > seq%nseq .or. &
                     (seq%nid == seq%nseq .and. seq%amin > a .and. seq%dmin > d) ) then
                   seq%nseq = seq%nid
                   seq%amin = a
                   seq%dmin = d
                   seq%tscale = sc
                   seq%tflux = tflux
                   seq%seq1(1:seq%nseq) = seq%id1(1:seq%nid)
                   seq%seq2(1:seq%nseq) = seq%id2(1:seq%nid)
                   seq%found = .true.
                end if
             end if

             seq%nid = seq%nid - 1

          end if

       end do ! j
    end do ! i

!    deallocate(idx2)

  end subroutine aqualung


  subroutine triangle(x1,y1,x2,y2,x3,y3,d3,u,v)

    ! computes coordinates in uv space

    real(REAL64), parameter :: tol = 1000*epsilon(0.0)

    real(REAL64), intent(in) :: x1,y1,x2,y2,x3,y3
    real(REAL64), intent(out) :: d3,u,v

    real(REAL64), dimension(0:3) :: d
    real(REAL64) :: x
    integer :: i,j

    ! compute lengths of triangle's sides
    d(1) = hypot(x1 - x2, y1 - y2)
    d(2) = hypot(x1 - x3, y1 - y3)
    d(3) = hypot(x2 - x3, y2 - y3)

    ! sorting, an implementation of the insert sort algorithm by N.Wirth
    do i = 2, 3
       x = d(i)
       d(0) = x
       j = i - 1
       do while ( x < d(j) )
          d(j+1) = d(j)
          j = j - 1
       end do
       d(j+1) = x
    end do

    if( d(3) > tol ) then
       ! compute triangle coordinates
       u = d(1)/d(3)
       v = d(2)/d(3)
       ! the longest side
       d3 = d(3)
    else
       ! a triangle is collapsing to a point
       u = 0
       v = 0
       d3 = 0
    end if

  end subroutine triangle


  subroutine distances(x,y,id,d)

    real(REAL64), dimension(:), intent(in) :: x,y
    integer, dimension(:), intent(in) :: id
    real(REAL64), dimension(:), intent(out) :: d
    integer :: i,j,k,l,n

    n = size(id)
    l = n  ! connecting head and tail
    do k = 1, n
       i = id(k)
       j = id(l)
       d(k) = hypot(x(i)-x(j),y(i)-y(j))
       l = k
    end do

  end subroutine distances

  real(REAL64) function sdisp(s,x1,y1,id1,sig1,x2,y2,id2,sig2)

    ! Approximated standard residuals of edge lenghts
    ! on frame and catalogue in pixels (lenghts in degrees
    ! are scaled):
    !            (d2 - s*d1)**2
    !            --------------
    !                sig**2
    ! where sig = 1.41*sqrt((s*sig1)**2 + sig2**2).

    real(REAL64), intent(in) :: s, sig1, sig2
    integer, dimension(:), intent(in) :: id1, id2
    real(REAL64), dimension(:), intent(in) :: x1,y1,x2,y2
    real(REAL64), dimension(size(id1)) :: d1, d2, r2

    call distances(x1,y1,id1,d1)
    call distances(x2,y2,id2,d2)
    r2 = (d2 - s*d1)**2 / ((s*sig1)**2 + sig2**2)
    sdisp = sum(r2)

  end function sdisp


  subroutine angles(x,y,id,a,da)

    ! computes angles in polygon given by the sequence
    ! the angles and theirs errors are in radians

    real(REAL64), dimension(:), intent(in) :: x,y
    integer, dimension(:), intent(in) :: id
    real(REAL64), dimension(:), intent(out) :: a,da
    real(REAL64), dimension(2) :: u,v
    real(REAL64) :: uu,vv,t
    integer :: i,l,k,n

    n = size(id)

    l = n
    do i = 1, n
       k = mod(i,n) + 1
       u = [ x(id(i)) - x(id(l)), y(id(i)) - y(id(l)) ]
       v = [ x(id(k)) - x(id(i)), y(id(k)) - y(id(i)) ]
       uu = norm2(u)
       vv = norm2(v)
       if( uu > 0 .and. vv > 0 ) then
          t = sum(u*v)/(uu*vv)
          t = max(-1.0_REAL64,min(1.0_REAL64,t))
          a(i) = acos(t)
          da(i) = abs(t)*sqrt(1/uu**2 + 1/vv**2)
       else
          a(i) = acos(1.0_REAL64)
          da(i) = acos(1.0_REAL64)
       end if
       l = i
    end do

  end subroutine angles


  function angels(x1,y1,id1,sig1,x2,y2,id2,sig2)

    ! returns sum of squares suitable for test by Xi2 distribution
    ! of residuals of angles:
    !
    !            |a1 - a2|**2
    !            ---------
    !              sig**2
    !
    ! where sig=sqrt(da1**2 + da2**2)

    real(REAL64) :: angels
    integer, dimension(:), intent(in) :: id1, id2
    real(REAL64), dimension(:), intent(in) :: x1,y1,x2,y2
    real(REAL64), intent(in) :: sig1, sig2
    real(REAL64), dimension(size(id1)) :: a1,a2,da1,da2,r

    call angles(x1,y1,id1,a1,da1)
    call angles(x2,y2,id2,a2,da2)

    !  write(*,*) 'da1',57.3*da1
    !  write(*,*) 'da2',57.3*da2
    !  write(*,*) 'a1-a2',57.3*(a1-a2)

    if( all(da1 > 0) .and. all(da2 > 0) ) then
       r = (a1 - a2)/sqrt((sig1*da1)**2 + (sig2*da2)**2)
       angels = sum(r**2)
    else
       angels = 7*size(r)
    end if

  end function angels

  real function qchi2_approx(n)

    ! approximation of quantiles of chi2 distribution for choosed probability
    ! accuracy for n < 30 is poor (few thents only)
    ! u are quantiles of N(0,1):
    !    P     u
    !   0.9   1.282
    !   0.95  1.645
    !   0.99  2.326
    !   0.999 3.090
    !

    real, parameter :: u = 1.645
    integer, intent(in) :: n

    qchi2_approx = (sqrt(2.0*n - 1.0) + u)**2 / 2.0

  end function qchi2_approx

  subroutine sequences_init(this, maxmatch)

    class(sequences) :: this
    integer, intent(in) :: maxmatch
    integer :: n

    this%reciprocal = .false.
    allocate(this%id1(maxmatch),this%id2(maxmatch), &
         this%seq1(maxmatch),this%seq2(maxmatch), &
         this%idx1(maxmatch),this%idx2(maxmatch), &
         this%rflux(maxmatch),this%scale(maxmatch),this%qchi2(maxmatch))

    do n = 1, maxmatch
       this%qchi2(n) = qchi2_approx(n)
    end do

  end subroutine sequences_init

  subroutine sequences_clean(this)

    type(sequences) :: this

    deallocate(this%id1,this%id2,this%idx1,this%idx2,this%seq1,this%seq2, &
         this%rflux,this%qchi2,this%scale)

  end subroutine sequences_clean


  real(REAL64) function dataseq_enorm(this,i,j)

    class(dataseq), intent(in out) :: this
    integer, intent(in) :: i,j
    real(REAL64) :: d

    if( this%dist_def(i,j) ) then
       dataseq_enorm = this%dist(i,j)
       return
    end if

    d = hypot(this%x(i) - this%x(j), this%y(i) - this%y(j))
    this%dist(i,j) = d
    this%dist(j,i) = d
    this%dist_def(i,j) = .true.
    dataseq_enorm = d

  end function dataseq_enorm

  subroutine dataseq_triangle(this,i1,i2,i3,u,v,c)

    ! computes coordinates in uv space

    class(dataseq), intent(in out) :: this
    integer, intent(in) :: i1,i2,i3
    real(REAL64), intent(out) :: u,v,c

    real(REAL64), dimension(0:3) :: d
    real(REAL64) :: x
    integer :: i,j

    if( this%tris_def(i1,i2,i3) ) then

       u = this%tris(i1,i2,i3,1)
       v = this%tris(i1,i2,i3,2)
       c = this%tris(i1,i2,i3,3)

       return
    end if

    ! compute lengths of triangle's sides
    d(1) = this%enorm(i1,i2)
    d(2) = this%enorm(i1,i3)
    d(3) = this%enorm(i2,i3)

    ! sorting, an implementation of the insert sort algorithm by N.Wirth
    do i = 2, 3
       x = d(i)
       d(0) = x
       j = i - 1
       do while ( x < d(j) )
          d(j+1) = d(j)
          j = j - 1
       end do
       d(j+1) = x
    end do

    ! compute coordinates in the triangle space
    u = d(1)/d(3)
    v = d(2)/d(3)

    ! the longest side
    c = d(3)

    this%tris(i1,i2,i3,:) = [u,v,c]
    this%tris_def(i1,i2,i3) = .true.

  end subroutine dataseq_triangle


  subroutine dataseq_init(this, x,y,f, m)

    class(dataseq) :: this
    real(REAL64), dimension(:), pointer, intent(in) :: x,y,f
    integer, intent(in) :: m

    this%x => x
    this%y => y
    this%f => f

    allocate(this%tris(m,m,m,3),this%dist(m,m),this%tris_def(m,m,m), &
         this%dist_def(m,m))

    ! the logical masks seems to be slightly faster, perhaps in conditions
    this%tris_def = .false.
    this%dist_def = .false.

  end subroutine dataseq_init

  subroutine dataseq_clean(this)

    type(dataseq) :: this

    deallocate(this%tris,this%dist,this%tris_def,this%dist_def)

  end subroutine dataseq_clean

  subroutine progress2(id1,id2,c,t,a,d,q)

    integer, dimension(:), intent(in) :: id1,id2
    real(REAL64), intent(in) :: c,t,a,d,q
    integer :: n
    character(len=80) :: fmt

    n = size(id2)
    write(fmt,"(a,i0,a,i0,a)") "(i0,a,",n,"(1x,i0),a,",n, &
         "(1x,i0),a,f0.1,2x,f0.3,a,2(2x,g0.3),a,f0.1)"
    write(error_unit,fmt) n,":",id1,' -> ',id2,' | ',c,t,': ',a,d,' < ',q

  end subroutine progress2


end module matcher
