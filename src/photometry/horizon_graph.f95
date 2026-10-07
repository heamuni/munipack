
! plot of likelihood and realted functions

! make install && gfortran -Wall -I../lib -o horizon fmin.f95 calibre.f95 horizon_graph.f95 -L../lib -L../minpack -lrstat -lsort -llmin -lminpacks -lminpack -llmin -lm && ./horizon

program horizon_graph

  use etacalibre

  integer, parameter :: dbl = selected_real_kind(15)
  real(dbl), dimension(36) :: cts,dcts,pht,dpht
  real(dbl), dimension(2) :: fv
  real(dbl) :: ctph,sig
  integer :: i,j

  
  ! E101_std_000002_sci.log, 20130402
  pht = (/ &
       4.69611E+07, 6.81299E+06, 3.92770E+06, 2.60938E+06, 1.97214E+06, 1.66471E+06, &
       2.05560E+06, 1.47822E+06, 1.04457E+06, 1.14430E+06, 8.76875E+05, 7.75065E+05, &
       8.65641E+05, 7.63727E+05, 8.90713E+05, 6.03862E+05, 6.04976E+05, 4.88583E+05, &
       4.64877E+05, 4.16618E+05, 4.36654E+05, 3.02369E+05, 3.06294E+05, 3.17203E+05, &
       2.84013E+05, 88907., 1.96851E+05, 2.01436E+05, 1.52524E+05, 1.20709E+05, &
       1.63885E+05, 1.19822E+05, 1.59565E+05, 2.05560E+05, 97395., 74223. /)

  dpht = (/ &
       2.16271E+06, 1.88256E+05, 1.08530E+05, 72102., 54494., 61332., 56800., 54461., &
       38485., 42159., 40383., 28555., 63785., 35172., 32816., 22248., 27861., 36001.,&
       29973., 23024., 32175., 16710., 22569., 20452., 20928., 4913.4, 10879., 20409.,&
       8429.1, 8894.4, 7547.4, 20969., 13227., 1893.3, 3588.3, 6836.4 /)

  cts = (/ &
       15226., 1.67221E+06, 9.62377E+05, 6.43891E+05, 4.88831E+05, 4.14225E+05, &
       5.07297E+05, 3.61345E+05, 2.55680E+05, 2.80821E+05, 2.13939E+05, 1.89409E+05, &
       2.13363E+05, 1.90943E+05, 2.16388E+05, 1.51464E+05, 1.51199E+05, 1.21140E+05, &
       1.12006E+05, 1.01042E+05, 1.07607E+05, 74211., 75318., 77542., 69805., 21779.,&
       47927., 48658., 37025., 29353., 39961., 30006., 37721., 32405., 24649., 17591. /)

  dcts = (/ &
       123.61, 1293.1, 981.01, 802.43, 699.17, 643.61, 712.26, 601.13, 505.66, 529.94, &
       462.55, 435.22, 461.92, 436.98, 465.19, 389.20, 388.86, 348.07, 334.69, 317.89, &
       328.05, 272.44, 274.46, 278.48, 264.22, 147.61, 218.94, 220.61, 192.44, 171.36, &
       199.93, 173.26, 194.24, 180.04, 157.03, 132.67 /)

  cts(1) = pht(1)/4.069
  dcts(1) = dpht(1)/2
  cts(34) = pht(34)/4.069
  dcts(34) = dpht(34)/2
  
  open(1,file='grad')
  do i = -80,80
     ctph = 4.0685 + i/10000.0
     do j = 2,77
        sig = j/100.0
        call graph(pht,dpht,cts,dcts,ctph,sig,"grad",fv)
        write(1,*) ctph,sig,sum(fv**2)
     end do
     write(1,*)
  end do
  close(1)

!  open(1,file='grad1')
!  do i = -80,80,10
!     ctph = 4.0685 + i/10000.0
!     j = 1
!     do while ( j <= 78 )
!        sig = j/100.0
!        call graph(pht,dpht,cts,dcts,ctph,sig,"grad",fv)
!        write(1,*) ctph,sig,sum(fv**2)
!        if( j < 20 ) then
!           j = j + 1
!        else
!           j = j + 7
!        end if
!     end do
!     write(1,*)
!  end do
!  close(1)
!
  open(2,file='like')
  do i = -80,80
     ctph = 4.0685 + i/1000.0
     do j = 2,77
        sig = j/100.0
        call graph(pht,dpht,cts,dcts,ctph,sig,"like",fv)
        write(2,*) ctph,sig,fv(1)
     end do
     write(2,*)
  end do
  close(2)



end program horizon_graph
