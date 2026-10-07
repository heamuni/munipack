!
!  Testing utility for WCSLIB F2008+ interface
!
!  cd lib/
!  make && gfortran -Wall -fcheck=all -g -p cfitsio.f08 fitsio.f08 wcslib.f08 test_wcslib.f08 -L. -ltitsio -lcfitsio -lwcs  -lm && ./a.out
!
! Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
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

program test_wcslib

  use wcslib
  use fitsio
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  implicit none

!  character(len=FLEN_CARD), dimension(:), allocatable :: header
  character(len=:), allocatable :: header
  character(len=FLEN_CARD), dimension(1) :: exclist
  integer :: status, nreject, nwcs, &
       WCSHDO_all, nkeyrec, nkeys
  type(fitsfiles) :: fits
!  integer, dimension(27) :: alts
  character(len=FLEN_FILENAME) :: filename
!  type(wcsprm) :: wcs
  type(wcspar) :: wcs

  exclist(1) = 'PSF'

!  goto 555

  call get_command_argument(1,filename)
  ! filename = '/var/tmp/artificial.fits'

  status = 0
  call fits_open_image(fits,filename,FITS_READONLY,status)
!  call fits_get_hdrspace(fits,keysexist,morekeys,status)
!  write(*,*) keysexist,morekeys,status

!!$  keynum = 1
!!$  do !i = 1, keysexist
!!$     call fits_read_record(fits,keynum,rec,status)
!!$     if( status /= 0 ) exit
!!$!     write(*,*) trim(rec)
!!$     keynum = keynum + 1
!!$  end do
!!$  write(*,*) keynum

!!  call wcspih(header,
!  status = 0
  call fits_hdr2str(fits,1,exclist,nkeys,header,status)

  call fits_close_file(fits,status)
  write(error_unit,*) 'fits status = ',status

!  write(error_unit,*) header
!  do i = 1, size(header)
!     write(error_unit,*) trim(header(i))
!  end do

!  rec = "x'0x000FFFFF'"
!  read(rec,*) WCSHDR_all
  !  WCSHDR_all = int(Z'0x000FFFFF')
!  WCSHDR_all = 1048575
  status = wcs_image_header(header,nkeys,WCSHDR_all, -3, nreject, nwcs,wcs)
  write(*,*) nreject, nwcs, status, c_associated(wcs % ptr)
  deallocate(header)

!  write(*,*) wcs % flag, wcs % flag

!  status = wcs_index(nwcs,wcs, alts)
!  write(*,*) status, alts

    block
      type(c_ptr) :: val!,iptr,jptr,wptr
!      integer, target :: xxx
      integer(kind=C_INT), target :: what,i,j
      real(REAL64), dimension(:), pointer :: crpix
      real(REAL64), target :: crpix1
      real(REAL64), dimension(2), target :: crpix0
      integer, pointer :: naxis
      what = 102
      crpix0 = 299
      crpix1 = 300
      val = C_LOC(crpix0(1))
      i = 1
      j = 0
!      wptr = C_LOC(what)
!      iptr = C_LOC(i)
!      jptr = C_LOC(j)
      status = wcsput(wcs % ptr,what,val,i,j)
      status = wcsget(wcs % ptr,what,val)
      call c_f_pointer(val,crpix,[2])
      write(*,*) 'crpix',crpix
      what = 101
      status = wcsget(wcs % ptr,what,val)
      call c_f_pointer(val,naxis)
      write(*,*) 'naxis',naxis
    end block
!    stop

  status = wcs_set(wcs)
  write(*,*) status


  status = wcsprt(wcs % ptr)
  write(*,*) 'prt:',status

  block
    ! distortions
    real(REAL64), dimension(:,:), allocatable :: pv
    real(REAL64), dimension(:), allocatable :: dpv
    integer :: naxis, npv, n, i, j, k

    ! it look as the wcsget() way does not works properly
    ! as npv is set on zero, and the internal loop in wcs_f.c is not executed
    status = wcs_get(wcs,WCS_NAXIS,naxis)
    status = wcs_get(wcs,WCS_NPV,npv)
    write(*,*) 'naxis, npv: ',naxis, npv, status

    allocate(dpv(3*npv))
!    do j = 1, npv
    status = wcs_get(wcs,WCS_PV,dpv)
!    write(*,*) status
!    write(*,*) dpv
!    end do

    n = npv/naxis
    allocate(pv(naxis,n))
    k = 1
    do i = 1, naxis
       do j = 1, n
          pv(i,j) = dpv(k)
          k = k + 3
       end do
    end do

    ! the proper way
    call wcs_get_dis(wcs,pv)

    write(*,*) 'pv(1,:): ',pv(1,:)
    write(*,*) 'pv(2,:): ',pv(2,:)
  end block

!!$  block
!!$    type(AstroWCS) :: wcs
!!$    real(REAL64), dimension(2) :: crval, crpix
!!$    real(REAL64) :: sc, rot, refl
!!$    real(REAL64), dimension(:,:), allocatable :: pv
!!$
!!$    crval = [ 333.333, 22.222 ]
!!$    crpix = [ 500, 500 ]
!!$    sc = 1 / 3600.0
!!$    rot = -20
!!$    refl = 1
!!$
!!$    call wcs % init('TAN',crval,crpix,sc,rot,refl)
!!$
!!$    status = wcsprt(wcs % ptr)
!!$    write(*,*) 'prt:',status
!!$
!!$!    call wcs_ini_dis(wcs)
!!$!    call wcs_get_dis(wcs,pv)
!!$
!!$!    write(*,*) 'pv(1,:): ',pv(1,:)
!!$!    write(*,*) 'pv(2,:): ',pv(2,:)
!!$  end block

!  stop

 ! write(*,*) '*'
!  write(*,*) header
!  write(*,*) '**'

  WCSHDO_all = 255
  status = wcs_header_output(WCSHDO_all,wcs,nkeyrec,header)
  write(*,*) status
  write(*,*) header
  deallocate(header)

!  do i = 1, size(header)
!     write(error_unit,*) trim(header(i))
!  end do

!  status = wcsset(wcs)

  block
    real(REAL64), dimension(:,:), allocatable :: pix,img,world
    real(REAL64), dimension(:), allocatable :: phi,theta
    integer, dimension(:), allocatable :: stat

    allocate(pix(3,2),img(3,2),world(3,2),phi(3),theta(3),stat(3))

    pix(1,:) = [ 100,100 ]
    pix(2,:) = [ 200,200 ]
    pix(3,:) = [ 300,300 ]
    !  img = 0
    !  world = 0
    !  phi = 0
    !  theta = 0
    !  stat = 0
    status = wcs_pixel2world(wcs,pix,img,phi,theta,world,stat)
    !  write(*,*) status,pix,img,phi,theta,world,stat
    block
      integer :: i
      do i = 1, size(pix,1)
         write(*,*) i,pix(i,:),world(i,:), img(i,:),stat(i),phi(i),theta(i)
         ! phi, theta reports uninitialised jump by valgrind
      end do
    end block

    world(1,:) = [ 56.6647, 24.1036 ]
    world(2,:) = [ 56.5, 24.1 ]
    world(3,:) = [ 56.6, 24.5 ]
    status = wcs_world2pixel(wcs,world,phi,theta,img,pix,stat)
    !  write(*,*) status,pix,img,phi,theta,world,stat
    block
      integer :: i
      do i = 1, size(pix,1)
         write(*,*) i,pix(i,:),world(i,:), img(i,:),stat(i),phi(i),theta(i)
         ! phi, theta reports uninitialised jump by valgrind
      end do
    end block

    call wcs_free(wcs)
    write(*,*) status

    deallocate(pix,img,world,phi,theta,stat)
!    deallocate(header)
  end block


  block
    real(REAL64), dimension(2) :: crpix, crval
    real(REAL64), dimension(2,2) :: cd
    real(REAL64), dimension(2,7) :: pv
    character(len=88), dimension(2) :: ctype
    character(len=88) :: radesys
    integer :: i, j

    wcs = wcs_alloc()
    status = wcs_ini(2,wcs)

    crpix = [ 335, 335 ]
    do i = 1, 2
       status = wcs_put(wcs,WCS_CRPIX,crpix(i),i,0)
    end do
    crval = [ 50, 25 ]
    do i = 1, 2
       status = wcs_put(wcs,WCS_CRVAL,crval(i),i,0)
    end do
    status = wcs_put(wcs,WCS_CTYPE,'RA---TAN',1,0)
    status = wcs_put(wcs,WCS_CTYPE,'DEC--TAN',2,0)

    cd(1,:) = [ -2.4649788212210852E-04, 2.5684664959419055E-06 ]
    cd(2,:) = [  2.5684664959419055E-06, 2.4649788212210852E-04 ]
    do i = 1, 2
       do j = 1,2
          status = wcs_put(wcs,WCS_CD,cd(i,j),i,j)
       end do
    end do

    status = wcs_set(wcs)
    write(*,*) 'wcs_set: ', status

    status = wcs_header_output(WCSHDO_all,wcs,nkeyrec,header)
    write(*,*) status
    write(*,*) header

    status = wcs_get(wcs,WCS_CRVAL,crval)
    write(*,*) 'crval = ',crval
    status = wcs_get(wcs,WCS_CRPIX,crpix)
    write(*,*) 'crpix = ',crpix
    status = wcs_get(wcs,WCS_CTYPE,ctype)
    write(*,*) 'ctype = ',trim(ctype(1)),' | ',trim(ctype(2))
    status = wcs_get(wcs,WCS_CD,cd)
    write(*,*) 'cd = ',real(cd(1,:))
    write(*,*) 'cd = ',real(cd(2,:))
    status = wcs_get(wcs,WCS_RADESYS,radesys)
    write(*,*) "radesys = '",trim(radesys),"'"

    pv = 0
    pv(1,2) = 1
    pv(2,3) = 1
!    ctype = 'TPV'
    status = wcs_put(wcs,WCS_CTYPE,'RA---TPV',1,0)
    status = wcs_put(wcs,WCS_CTYPE,'DEC--TPV',2,0)

    do i = 1, size(pv,1)
       do j = 1, size(pv,2)
          status = wcs_put(wcs,WCS_PV, pv(i,j), i, j)
!          write(*,*) i,j,status
       end do
    end do

!    call wcs_put_dis(wcs,ctype,pv)
    status = wcs_set(wcs)
    write(*,*) status

    status = wcsprt(wcs % ptr)
    call wcs_perr(wcs)
    call wcs_free(wcs)
    write(*,*) status
!    stop
  end block


  block
    type(prjpar) :: prj
    real(REAL64), dimension(2) :: crval
    real(REAL64), dimension(3) :: alpha, delta, u, v
    real(REAL64) :: R0 = 57.3
    integer, dimension(3) :: stat
    character(len=3) :: ctype, a
    integer :: i, status

    status = wcserr_enable(1);
    write(*,*) 'prj:'

    crval = [ 30, 60 ]
    crval = [  299.93283007737159_REAL64, 22.772932147394929_REAL64 ]
    ctype = "TAN"

    status = prj_init(prj)
    write(*,*) 'prjini:',status,c_associated(prj%ptr),associated(prj%prj)
    status = prj_perr(prj,"Error:")
    status = prj_put(prj,PRJ_FLAG,0,1)
    write(*,*) 'prjflag',status
    status = prj_put(prj,PRJ_CODE,ctype,1)
    status = prj_put(prj,PRJ_R0,R0,1)
    status = prj_put(prj,PRJ_THETA0,crval(2),1)
    status = prj_put(prj,PRJ_PHI0,crval(1),1)
!    status = prj_put(prj,PRJ_X0,0.0_REAL64,1)
!    status = prj_put(prj,PRJ_Y0,0.0_REAL64,1)
    write(*,*) 'prjphi0',status
    status = prj_set(prj)
    write(*,*) 'prjset',status,c_associated(prj%prj%prjs2x)

    status = prj_get(prj,PRJ_CODE,a)
    write(*,*) status, a

    alpha = crval(1) + [ -1, 0, 1 ]
    delta = crval(2) + [ -1, 0, 1 ]
    alpha(1) = 299.27718766610883_REAL64
    delta(1) = 22.357468143545056_REAL64
    ! ma byt 0.60639345349398233      -0.41415148831720683
    status = prj_cel2cart(prj,alpha,delta,u,v,stat)
    do i = 1, size(alpha)
       write(*,*) alpha(i),delta(i),u(i),v(i),stat(i)
    end do


    status = prj_free(prj)
    write(*,*) status

  end block

  555 continue


  block
    type(wcspar) :: wcs
    real(REAL64), dimension(2) :: crval
    real(REAL64), dimension(3) :: alpha, delta
    real(REAL64), dimension(3,2) :: pix,img,world,world1
    real(REAL64), dimension(3) :: phi,theta
    integer, dimension(3) :: stat
    character(len=3) :: ctype
    integer :: i, status

    status = wcserr_enable(1);

    crval = [  299.93283007737159_REAL64, 22.772932147394929_REAL64 ]
    ctype = "TAN"

    wcs = wcs_alloc()
    status = wcs_ini(2,wcs)
    do i = 1, 2
       status = wcs_put(wcs,WCS_CRVAL,crval(i),i,0)
    end do
    status = wcs_put(wcs,WCS_CTYPE,'RA---TAN',1,0)
    status = wcs_put(wcs,WCS_CTYPE,'DEC--TAN',2,0)
!    status = wcs_put(wcs,WCS_CD,-1.0_REAL64,1,1)
!    status = wcs_put(wcs,WCS_CDELT,-1.0_REAL64,1,0)
!    status = wcs_put(wcs,WCS_PC,-1.0_REAL64,1,1)

    status = wcs_set(wcs)
    write(*,*) 'wcs_set: ', status

    status = wcs_header_output(WCSHDO_all,wcs,nkeyrec,header)
    write(*,*) header
    deallocate(header)


    alpha(1) = 299.27718766610883_REAL64
    delta(1) = 22.357468143545056_REAL64
    alpha(2) = 299.29328739856180_REAL64
    delta(2) = 22.495856679832567_REAL64
    alpha(3) = 299.53273685594661_REAL64 + 0.5
    delta(3) = 22.850505560821176_REAL64

    do i = 1, 3
       world(i,:) = [ alpha(i), delta(i) ]
    end do!
!    phi(1) = 0
!    theta(1) = 90
    status = wcs_world2pixel(wcs,world,phi,theta,img,pix,stat)
    !  write(*,*) status,pix,img,phi,theta,world,stat
    do i = 1, size(pix,1)
       write(*,'(i3,8f9.3,i2)') i,world(i,:), phi(i),theta(i),img(i,:),pix(i,:),stat(i)
    end do

    status = wcs_pixel2world(wcs,pix,img,phi,theta,world1,stat)
    do i = 1, size(pix,1)
       write(*,*) i,world1(i,:)-world(i,:)
    end do

    call wcs_free(wcs)

  end block


end program test_wcslib
