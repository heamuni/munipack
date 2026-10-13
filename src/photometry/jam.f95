!
!  Merging of photometric catalogues
!
!  Copyright © 2013-7 F.Hroch (hroch@physics.muni.cz)
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

module jamming


  implicit none

  integer, parameter, private :: dbl = selected_real_kind(15)
  real(dbl), parameter, private :: rad = 57.295779513082322865_dbl

contains


  subroutine cooid(refra,refdec,ra,dec,tol,idx)

    real(dbl), dimension(:), intent(in) :: ra, dec, refra, refdec
    real(dbl), intent(in) :: tol
    integer, dimension(:), intent(out) :: idx
    integer :: i,j,jmin
    real(dbl) :: r,rmin

    idx = 0

    do i = 1,size(ra)
       rmin = tol
       jmin = 0
       do j = 1,size(refra)
          r = spmetr(refra(j),refdec(j),ra(i),dec(i))
!          write(*,*) refra(j),refdec(j),ra(i),dec(i),r
          if( r < rmin )then
             rmin = r
             jmin = j
          end if
       end do

       if( jmin > 0 ) then
          idx(jmin) = i
       end if

    end do

  end subroutine cooid

  function spmetr(a1,d1,a2,d2)

    ! angular distance in degrees

    real(dbl), parameter :: rad = 57.295779513082322865_dbl
    real(dbl) :: spmetr
    real(dbl), intent(in) :: a1,d1,a2,d2
    real(dbl) :: cosr

    cosr = sin(d1/rad)*sin(d2/rad) + cos(d1/rad)*cos(d2/rad)*cos((a1 - a2)/rad)
    spmetr = rad*acos(min(cosr,1.0_dbl))
!    spmetr = sqrt((d1 - d2)**2 + cos((a1+a2)/2/rad)**2*(a1 - a2)**2)

  end function spmetr


  subroutine jamcat(tol, refra, refdec, mag, dmag, ra,dec,dn,ddn, cts, dcts )

    real(dbl), dimension(:), intent(in) :: tol
    real(dbl), dimension(:), intent(in out) :: refra, refdec
    real(dbl), dimension(:,:), intent(in out) :: mag, dmag
    real(dbl), dimension(:,:), intent(in) :: dn,ddn, ra, dec
    real(dbl), dimension(:,:), allocatable, intent(out) :: cts,dcts

    real(dbl), dimension(:,:), allocatable :: xcts,xdcts
    integer, dimension(:), allocatable :: id
    integer :: ncat,ncols,j,n,l


    ncat = size(mag,1)
    ncols = size(dn,2)

    allocate(id(ncat))
    allocate(xcts(ncat,ncols),xdcts(ncat,ncols))
    xcts = -1
    xdcts = 0

    do n = 1, ncols

       call cooid(refra,refdec,ra(:,n),dec(:,n),tol(n),id)

       do l = 1,ncat
          if( id(l) > 0 ) then
             j = id(l)
!             write(*,*) dn(j,n),mag(l,n)
             if(  dn(j,n) > 0 .and. mag(l,n) < 99 ) then
                xcts(l,n) = dn(j,n)
                xdcts(l,n) = ddn(j,n)
             end if
          end if
       end do

    end do


    ! select valid records
    ncat = 0
    do j = 1,size(mag,1)
       if( all(xcts(j,:) > 0) ) then
!          if(  all(sqrt(xcts(j,:))/xcts(j,:) < 0.005) ) then
             ncat = ncat + 1
             refra(ncat) = refra(j)
             refdec(ncat) = refdec(j)
             mag(ncat,:) = mag(j,:)
             dmag(ncat,:) = dmag(j,:)
             xcts(ncat,:) = xcts(j,:)
             xdcts(ncat,:) = xdcts(j,:)
!          end if
       end if
    end do

    deallocate(id)

    allocate(cts(ncat,ncols),dcts(ncat,ncols))
    cts(1:ncat,:) = xcts(1:ncat,:)
    dcts(1:ncat,:) = xdcts(1:ncat,:)
    deallocate(xcts,xdcts)


  end subroutine jamcat


  subroutine jamcatx(cat,phsystable,keys,col_ra,col_dec,col_mag,col_magerr,utol, &
       area, init_area, filters, filename, catid, filter, ra, dec, airmass,ph,dph, &
       cts,dcts, photosys_ref, photosys_instr, pairs)

    use phsysfits
    use photoconv
    use mfits

    character(len=*), intent(in) :: cat,phsystable,col_ra,col_dec,photosys_ref
    character(len=*), intent(out) :: photosys_instr
    character(len=*), dimension(:), intent(in) :: keys,col_mag,col_magerr,filters,filename
    real(dbl), intent(in) :: utol, area
    logical, intent(in) :: init_area
    real(dbl), dimension(:), allocatable, intent(out) :: ra, dec, airmass
    real(dbl), dimension(:,:), allocatable, intent(out) :: ph,dph,cts,dcts
    character(len=FLEN_VALUE), dimension(:), allocatable, intent(out) :: filter
    character(len=FLEN_VALUE), intent(out) :: catid
    integer, dimension(:,:), allocatable, intent(out) :: pairs

    real(dbl), allocatable, dimension(:) :: refra,refdec,xra,xdec,adu,dadu,exptime,areas
    real(dbl), allocatable, dimension(:,:) :: mag,dmag, xcts,xdcts
    real(dbl) :: ftol,tol,w, aper
    integer, dimension(:), allocatable :: id
    type(type_phsys) :: phsys
    integer :: status,nfiles,ncat,ncols,j,n,l

    ! setup to prove
    real(dbl), parameter :: waterline = 5
    real(dbl), parameter :: maxerr = 0.1

    call readcat(cat,(/col_ra,col_dec/),col_mag,col_magerr, &
         refra,refdec,mag,dmag,catid,status)
    if( status /= 0 ) stop 'Failed to read a catalogue.'

    ncat = size(mag,1)
    ncols = size(mag,2)
    nfiles = size(filename)

    if( ncols /= nfiles ) stop 'Count of color bands and files differs.'

    allocate(xcts(ncat,nfiles),xdcts(ncat,nfiles))
    xcts = -1
    xdcts = 0
    allocate(exptime(nfiles),areas(nfiles),filter(nfiles),airmass(nfiles))

    do n = 1, nfiles
       call readframe(filename(n),keys,waterline,maxerr,xra,xdec,adu,dadu,&
            ftol,exptime(n),areas(n),photosys_instr,filter(n), &
            init_area,0,aper,.false.,status,airmass(n))
       ! SAPER=0 is set at any case, only apperture corrected data are used !

       if( status /= 0 ) then
          write(*,*) "Failed to read `",trim(filename(n)),"'."
          stop 666
       end if

       if( utol > 0 ) then
          tol = utol
       else
          tol = ftol
       end if

       if( init_area ) areas = area

       if( filter(n) /= filters(n) ) write(*,*) "Warning: Different filters (?): `",&
            trim(filters(n)),"'(catalogue), `",trim(filter(n)),"'(frame)."

       allocate(id(ncat))
       call cooid(refra,refdec,xra,xdec,tol,id)

       do l = 1,ncat
          if( id(l) > 0 ) then
             j = id(l)
             if( adu(j) > 0 ) then
                xcts(l,n) = adu(j)
                xdcts(l,n) = dadu(j)
             end if
          end if
       end do

       deallocate(id,xra,xdec,adu,dadu)
    end do

    ! select just valid records
    ncat = 0
    do j = 1,size(mag,1)
       if( all(xcts(j,:) > epsilon(xcts)) .and. all(mag(j,:) < 99) ) then
          ncat = ncat + 1
          refra(ncat) = refra(j)
          refdec(ncat) = refdec(j)
          mag(ncat,:) = mag(j,:)
          dmag(ncat,:) = dmag(j,:)
          xcts(ncat,:) = xcts(j,:)
          xdcts(ncat,:) = xdcts(j,:)
       end if
    end do

    allocate(ph(ncat,ncols),dph(ncat,ncols),cts(ncat,nfiles),dcts(ncat,nfiles),&
         ra(ncat),dec(ncat))

    ra = refra(1:ncat)
    dec = refdec(1:ncat)

    call phselect(phsystable,photosys_ref,phsys)
    call phsyspairs(phsys,filters,pairs)
    call phsysmagph(phsys,filters,pairs,mag(1:ncat,:),dmag(1:ncat,:),ph,dph)
!    call phsysmagflux(phsys,filters,pairs,mag(1:ncat,:),dmag(1:ncat,:),flux,dflux)

!    write(*,'(1p,5g15.6)') ph(1,:)
!    write(*,'(1p,5g15.5)') flux(1,:)

!    do j = 1,ncat
!       write(*,'(5f10.3)') mag(j,:)
!       write(*,'(i5,10g10.3)') i,flux(j,:),dflux(j,:)
!    end do

!    call phsysphmag(phsys,filters,pairs,ph,dph,mag(1:ncat,:),dmag(1:ncat,:))
!    call phsysphmag(phsys,filters,pairs,ph,dph,flux,dmag(1:ncat,:))
!    do j = 1,ncat
!       write(*,'(10f10.3)') mag(j,:),flux(j,:),mag(j,:)-flux(j,:)
!    end do



!!$!    mag(1,:) = 0
!!$    call phsysmagflux(phsys,filters,mag(1:ncat,:),dmag(1:ncat,:),flux,dflux)
!!$    write(*,'(5f10.3)') mag(1,:)
!!$    write(*,'(5e10.3)') flux(1,:)
!!$    do j = 1,ncat
!!$!       flux(j,1) = flux(j,2)*10**(-0.4*(mag(j,1)-mag(j,2)))
!!$ !      flux(j,2) = flux(j,3)*10**(-0.4*(mag(j,2)-mag(j,3)))
!!$ !      flux(j,1) = flux(j,2)*10**(-0.4*(mag(j,1)-mag(j,2)))
!!$ !      flux(j,4) = flux(j,3)*10**(0.4*(mag(j,3)-mag(j,4)))
!!$ !      flux(j,5) = flux(j,4)*10**(0.4*(mag(j,4)-mag(j,5)))
!!$    end do
!!$!    write(*,'(5e10.3)') flux(1,:)
!!$    k = pairs(1,1)
!!$    do n = 2,size(pairs,1)
!!$       i = pairs(n,1)
!!$       j = pairs(n,2)
!!$       if( pairs(n,1) < k ) then
!!$!          i = pairs(n,1)
!!$!          j = pairs(n,2)
!!$          flux(:,i) = flux(:,j)*10**(-0.4*(mag(1:ncat,i)-mag(1:ncat,j)))
!!$       else if( pairs(n,2) > k ) then
!!$!          i = pairs(n,1)
!!$!          j = pairs(n,2)
!!$          flux(:,j) = flux(:,i)*10**(0.4*(mag(1:ncat,i)-mag(1:ncat,j)))
!!$       end if
!!$       write(*,*) i,j,mag(1,i)-mag(1,j)
!!$!       do l = 1,ncat
!!$!          flux(l,i) = flux(l,j)*10**(-0.4*(mag(l,i)-mag(l,j)))
!!$!       end do
!!$    end do
!!$
!!$    write(*,'(5e10.3)') flux(1,:)
!!$    write(*,'(5a10)') filters
!!$
!!$    write(*,'(5e10.3)') ph(1,:)
!!$    !       ph(j,2) = ph(j,2)*(1+0.5*(75*0.7)**2/(553.0*(553.0-444.0))*(1.96*flux(j,1)/flux(j,2)-1))
!!$    !    write(*,*) flux(j,1)/flux(j,2),0.5*(75*0.7)**2/(553.0*(553.0-444.0))*(1.96*flux(j,1)/flux(j,2)-1)
!!$    k = pairs(1,1)
!!$    do n = 2, size(pairs,1)
!!$       i = pairs(n,1)
!!$       j = pairs(n,2)
!!$       if( pairs(n,1) < k ) then
!!$          i = pairs(n,1)
!!$          j = pairs(n,2)
!!$          w = phsys%lam_eff(j) / phsys%lam_eff(i)
!!$          ph(:,i) = ph(:,j)*(flux(:,i)/flux(:,j))*w
!!$       else if( pairs(n,2) > k ) then
!!$          i = pairs(n,1)
!!$          j = pairs(n,2)
!!$          w = phsys%lam_eff(i) / phsys%lam_eff(j)
!!$          ph(:,j) = ph(:,i)*(flux(:,j)/flux(:,i))*w
!!$       end if
!!$       write(*,*) i,j,w
!!$       do l = 1,ncat
!!$!          ph(l,i) = ph(l,j)*(flux(l,i)/flux(l,j))*w!*1.25
!!$       end do
!!$    end do
!!$!    write(*,*) phsys%lam_eff(2)/phsys%lam_eff(1)
!!$!    ph(1,1) = ph(1,2)*(flux(1,1)/flux(1,2))*1.26
!!$!    ph(1,4) = ph(1,3)*(flux(1,4)/flux(1,3))*0.80
!!$!    ph(1,5) = ph(1,4)*(flux(1,5)/flux(1,4))*0.80
!!$    write(*,'(5e10.3)') ph(1,:)

!!$    do l = 1,ncat
!!$       do n = 1,size(ph,2)
!!$          write(*,'(1pg10.3)',advance="no") &
!!$               !               mag(l,n)-mag(l,n+1) - 2.5*log10(flux(l,n+1)/flux(l,n))
!!$               flux(l,n)
!!$       end do
!!$       write(*,*)
!!$    end do
!!$
!!$    do l = 1,ncat
!!$       do n = 1,size(ph,2)
!!$          write(*,'(1pg10.3)',advance="no") &
!!$               !               mag(l,n)-mag(l,n+1) - 2.5*log10(ph(l,n+1)/ph(l,n))
!!$               ph(l,n)
!!$       end do
!!$       write(*,*)
!!$    end do

!    stop

    ! rates per second and square metre (the errors of the counts
    ! scale as the counts; the square root was the scale of the variance)
    do n = 1,nfiles
       w = exptime(n) * areas(n)
       cts(:,n) = xcts(1:ncat,n) / w
       dcts(:,n) = xdcts(1:ncat,n) / w
    end do

    deallocate(mag,dmag,xcts,xdcts,exptime,areas)

  end subroutine jamcatx


  subroutine jamframes(tol, ndat, ra, dec, id )

    use phsysfits
    use photoconv
    use mfits

    real(dbl), dimension(:), intent(in) :: tol
    integer, dimension(:), intent(in) :: ndat
    real(dbl), dimension(:,:), intent(in) :: ra, dec
    integer, dimension(:,:), intent(out) :: id
    integer :: n

    id = 0

    forall(n=1:ndat(1))
       id(n,1) = n
    end forall

    do n = 2, size(ra,2)
       call cooid(ra(1:ndat(1),1),dec(1:ndat(1),1), &
            ra(1:ndat(n),n),dec(1:ndat(n),n),tol(n),id(:,n))
    end do

  end subroutine jamframes

  subroutine jamref(tol, refra, refdec, refph, drefph, ra,dec, cts, dcts)

    real(dbl), dimension(:), intent(in) :: tol
    real(dbl), dimension(:), intent(in out) :: refra, refdec
    real(dbl), dimension(:,:), allocatable, intent(in out) :: ra,dec, &
         refph, drefph, cts, dcts
!    real(dbl), dimension(:,:), intent(in) :: dn,ddn, ra, dec
!    real(dbl), dimension(:,:), allocatable, intent(out) :: cts,dcts

    real(dbl), dimension(:,:), allocatable :: xcts,xdcts
    integer, dimension(:), allocatable :: id
    integer :: nref,ncols,j,n,l


!    ncat = size(mag,1)
!    ncols = size(dn,2)

    nref= size(refra)
    ncols = size(refph,2)

    allocate(id(nref))
    allocate(xcts(nref,ncols),xdcts(nref,ncols))
    xcts = -1
    xdcts = 0

    do n = 1, ncols

       call cooid(refra,refdec,ra(:,n),dec(:,n),tol(n),id)

       do l = 1,nref
          if( id(l) > 0 ) then
             j = id(l)
             if( cts(j,n) > 0 .and. refph(l,n) > 0 ) then
                xcts(l,n) = cts(j,n)
                xdcts(l,n) = dcts(j,n)
             end if
          end if
       end do

    end do


    ! select just valid records
    nref = 0
    do j = 1,size(refph,1)
       if( all(xcts(j,:) > 0) .and. all(refph(j,:) > 0) ) then
!          if(  all(sqrt(xcts(j,:))/xcts(j,:) < 0.01) ) then
             nref = nref + 1
             refra(nref) = refra(j)
             refdec(nref) = refdec(j)
             refph(nref,:) = refph(j,:)
             drefph(nref,:) = drefph(j,:)
             xcts(nref,:) = xcts(j,:)
             xdcts(nref,:) = xdcts(j,:)
!          end if
       end if
    end do

    deallocate(id)

    deallocate(cts,dcts)
    allocate(cts(nref,ncols),dcts(nref,ncols))
    cts = xcts(1:nref,:)
    dcts = xdcts(1:nref,:)
    deallocate(xcts,xdcts)


  end subroutine jamref

end module jamming
