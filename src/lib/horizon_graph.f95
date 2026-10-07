
! plot of likelihood and realted functions

! make install && gfortran -Wall -o horizon_graph horizon_graph.f95 -L../lib -L../minpack -lrstat -lsort -llmin -lminpacks -lminpack -llmin -lm && ./horizon_graph

program horizon_graph

  use robustmean

  integer, parameter :: dbl = selected_real_kind(15)
  real(dbl), dimension(666) :: x
  real(dbl), dimension(2) :: fv
  real(dbl) :: t,sig
  integer :: i,j,n


  do i = 1, int(0.95*size(x))
     x(i) = gdis(0.0_dbl,1.0_dbl)
     n = i
  end do
  do i = n+1, size(x)
     x(i) = gdis(0.0_dbl,5.0_dbl)
  end do
  
  
  do i = -100,100,10
     t = 0.0 + i/100.0
     do j = 1,20
        sig = j/10.0
        call graph(x,t,sig,"grad",fv)
        write(*,*) t,sig,sum(fv**2)
     end do
     write(*,*)
  end do


contains

  function gdis(mean, sig)

    ! generate random data with gauss distribution

    real(dbl), intent(in) :: mean, sig
    real(dbl) :: gdis, x
    
    call random_number(x)
    gdis = invnorm(x)*sig + mean

  end function gdis


    ! excelent overwiew of algorithms !
    ! http://home.online.no/~pjacklam/notes/invnorm/

    function invnorm(p)

      !ren-raw chen, rutgers business school
      ! normal inverse
      ! translate from
      ! http://home.online.no/~pjacklam/notes/invnorm
      ! a routine written by john herrero
      real*8 invnorm
      real*8 p,p_low,p_high
      real*8 a1,a2,a3,a4,a5,a6
      real*8 b1,b2,b3,b4,b5
      real*8 c1,c2,c3,c4,c5,c6
      real*8 d1,d2,d3,d4
      real*8 z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   invnorm=z
      return

    end function invnorm


end program horizon_graph
