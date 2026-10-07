!
!  Fortran 90+ interface for Minpack
!
!  Copyright © 2010, 2013 F.Hroch (hroch@physics.muni.cz)
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
!  along with OakLeaf.  If not, see <http://www.gnu.org/licenses/>.

module minpack

  implicit none

  ! precision for double real
  integer, parameter, private :: dbl = selected_real_kind(15)

  interface

     subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
       integer, parameter :: dbl = selected_real_kind(15)
       integer, intent(in) :: m,n,ldfjac,mode
       real(dbl) :: x(n),fvec(m),fjac(ldfjac,n),xp(n),fvecp(m),err(m)
     end subroutine chkder

     subroutine covar(n,r,ldr,ipvt,tol,wa)
       ! it should be not available in all distributions (Debian)
       integer, parameter :: dbl = selected_real_kind(15)
       integer, intent(in) ::  n,ldr
       integer::  ipvt(n)
       real(dbl) :: tol
       real(dbl) :: r(ldr,n),wa(n)
     end subroutine covar

     subroutine dmchar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp, &
          maxexp,eps,epsneg,xmin,xmax)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: i,ibeta,iexp,irnd,it,iz,j,k,machep,maxexp,minexp, &
            mx,negep,ngrd
       real(dbl) :: a,b,beta,betain,betam1,eps,epsneg,one,xmax, &
            xmin,y,z,zero
     end subroutine dmchar

     subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,lr
       real(dbl) :: delta
       real(dbl) :: r(lr),diag(n),qtb(n),x(n),wa1(n),wa2(n)
     end subroutine dogleg

     function dpmpar(i)
       integer, parameter :: dbl = selected_real_kind(15)
       real(dbl) :: dpmpar
       integer, intent(in) :: i
     end function dpmpar

     function enorm(n,x)
       integer, parameter :: dbl = selected_real_kind(15)
       integer, intent(in) :: n
       real(dbl) :: enorm,x(n)
     end function enorm

     subroutine errjac(n,x,fjac,ldfjac,nprob)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,ldfjac,nprob
       real(dbl) :: x(n),fjac(ldfjac,n)
     end subroutine errjac

     subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,ldfjac,iflag,ml,mu
       real(dbl) :: epsfcn
       real(dbl) :: x(n),fvec(n),fjac(ldfjac,n),wa1(n),wa2(n)
       interface
          subroutine fcn(n,x,fvec,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(n), intent(out) :: fvec
          end subroutine fcn
       end interface
     end subroutine fdjac1

     subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,ldfjac,iflag
       real(dbl) :: epsfcn
       real(dbl) :: x(n),fvec(m),fjac(ldfjac,n),wa(m)
       interface
          subroutine fcn(m,n,x,fvec,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: m,n
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(m), intent(out) :: fvec
          end subroutine fcn
       end interface
     end subroutine fdjac2

     subroutine grdfcn(n,x,g,nprob)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,nprob
       real(dbl) :: x(n),g(n)
     end subroutine grdfcn

     subroutine hesfcn(n,x,h,ldh,nprob)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,ldh,nprob
       real(dbl) :: x(n),h(ldh,n)
     end subroutine hesfcn

     subroutine initpt(n,x,nprob,factor)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,nprob
       real(dbl) :: factor
       real(dbl) :: x(n)
     end subroutine initpt

     subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
          mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr, &
          qtf,wa1,wa2,wa3,wa4)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
       real(dbl) :: xtol,epsfcn,factor
       real(dbl) :: x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr), &
            qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
       interface
          subroutine fcn(n,x,fvec,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(n), intent(out) :: fvec
          end subroutine fcn
       end interface
     end subroutine hybrd

     subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,info,lwa
       real(dbl) :: tol
       real(dbl) :: x(n),fvec(n),wa(lwa)
       interface
          subroutine fcn(n,x,fvec,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(n), intent(out) :: fvec
          end subroutine fcn
       end interface
     end subroutine hybrd1

     subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode, &
          factor,nprint,info,nfev,njev,r,lr,qtf,wa1,wa2, &
          wa3,wa4)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,ldfjac,maxfev,mode,nprint,info,nfev,njev,lr
       real(dbl) :: xtol,factor
       real(dbl) :: x(n),fvec(n),fjac(ldfjac,n),diag(n),r(lr), &
            qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
       interface
          subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n, ldfjac
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(n), intent(out) :: fvec
            real(dbl), dimension(ldfjac,n), intent(out) :: fjac
          end subroutine fcn
       end interface
     end subroutine hybrj

     subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,ldfjac,info,lwa
       real(dbl) :: tol
       real(dbl) :: x(n),fvec(n),fjac(ldfjac,n),wa(lwa)
       interface
          subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n, ldfjac
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(n), intent(out) :: fvec
            real(dbl), dimension(ldfjac,n), intent(out) :: fjac
          end subroutine fcn
       end interface
     end subroutine hybrj1

     subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, &
          maxfev,diag,mode,factor,nprint,info,nfev,njev, &
          ipvt,qtf,wa1,wa2,wa3,wa4)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
       integer :: ipvt(n)
       real(dbl) :: ftol,xtol,gtol,factor
       real(dbl) :: x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n), &
            wa1(n),wa2(n),wa3(n),wa4(m)
       interface
          subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n, m, ldfjac
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(m), intent(out) :: fvec
            real(dbl), dimension(ldfjac,n), intent(out) :: fjac
          end subroutine fcn
       end interface
     end subroutine lmder

     subroutine lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,ldfjac,info,lwa
       integer :: ipvt(n)
       real(dbl) :: tol
       real(dbl) :: x(n),fvec(m),fjac(ldfjac,n),wa(lwa)
       interface
          subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n, m, ldfjac
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(m), intent(out) :: fvec
            real(dbl), dimension(ldfjac,n), intent(out) :: fjac
          end subroutine fcn
       end interface
     end subroutine lmder1

     subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn, &
          diag,mode,factor,nprint,info,nfev,fjac,ldfjac, &
          ipvt,qtf,wa1,wa2,wa3,wa4)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,maxfev,mode,nprint,info,nfev,ldfjac
       integer :: ipvt(n)
       real(dbl) :: ftol,xtol,gtol,epsfcn,factor
       real(dbl) :: x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n), &
            wa1(n),wa2(n),wa3(n),wa4(m)
       interface
          subroutine fcn(m,n,x,fvec,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n, m
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(m), intent(out) :: fvec
          end subroutine fcn
       end interface
     end subroutine lmdif

     subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,info,lwa
       integer :: iwa(n)
       real(dbl) :: tol
       real(dbl) :: x(n),fvec(m),wa(lwa)
       interface
          subroutine fcn(m,n,x,fvec,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n, m
            integer, intent(in out) :: iflag
            real(dbl), dimension(n), intent(in) :: x
            real(dbl), dimension(m), intent(out) :: fvec
          end subroutine fcn
       end interface
     end subroutine lmdif1

     subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,ldr
       integer :: ipvt(n)
       real(dbl) :: delta,par
       real(dbl) :: r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n),wa2(n)
     end subroutine lmpar

     subroutine lmstr(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, &
          maxfev,diag,mode,factor,nprint,info,nfev,njev, &
          ipvt,qtf,wa1,wa2,wa3,wa4)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
       integer :: ipvt(n)
       logical :: sing
       real(dbl) :: ftol,xtol,gtol,factor
       real(dbl) :: x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n), &
            wa1(n),wa2(n),wa3(n),wa4(m)
       interface
          subroutine fcn(m,n,x,fvec,fjrow,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n, m
            integer, intent(in out) :: iflag
            real(dbl), dimension(:), intent(in) :: x
            real(dbl), dimension(:), intent(out) :: fvec, fjrow
          end subroutine fcn
       end interface
     end subroutine lmstr

     subroutine lmstr1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,ldfjac,info,lwa
       integer :: ipvt(n)
       real(dbl) :: tol
       real(dbl) :: x(n),fvec(m),fjac(ldfjac,n),wa(lwa)
       interface
          subroutine fcn(m,n,x,fvec,fjrow,iflag)
            integer, parameter :: dbl = selected_real_kind(15)
            integer, intent(in) :: n, m
            integer, intent(in out) :: iflag
            real(dbl), dimension(:), intent(in) :: x
            real(dbl), dimension(:), intent(out) :: fvec, fjrow
          end subroutine fcn
       end interface
     end subroutine lmstr1

     subroutine objfcn(n,x,f,nprob)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,nprob
       real(dbl) :: f
       real(dbl) :: x(n)
     end subroutine objfcn

     subroutine qform(m,n,q,ldq,wa)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,ldq
       real(dbl) :: q(ldq,m),wa(m)
     end subroutine qform

     subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,lda,lipvt
       integer :: ipvt(lipvt)
       logical :: pivot
       real(dbl) :: a(lda,n),rdiag(n),acnorm(n),wa(n)
     end subroutine qrfac

     subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,ldr
       integer :: ipvt(n)
       real(dbl) :: r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
     end subroutine qrsolv

     subroutine r1mpyq(m,n,a,lda,v,w)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,lda
       real(dbl) :: a(lda,n),v(n),w(n)
     end subroutine r1mpyq

     subroutine r1updt(m,n,s,ls,u,v,w,sing)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,ls
       logical :: sing
       real(dbl) :: s(ls),u(m),v(n),w(m)
     end subroutine r1updt

     subroutine rwupdt(n,r,ldr,w,b,alpha,cos,sin)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,ldr
       real(dbl) :: alpha
       real(dbl) :: r(ldr,n),w(n),b(n),cos(n),sin(n)
     end subroutine rwupdt

     subroutine ssqfcn(m,n,x,fvec,nprob)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,nprob
       real(dbl) :: x(n),fvec(m)
     end subroutine ssqfcn

     subroutine ssqjac(m,n,x,fjac,ldfjac,nprob)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: m,n,ldfjac,nprob
       real(dbl) :: x(n),fjac(ldfjac,n)
     end subroutine ssqjac

     subroutine vecfcn(n,x,fvec,nprob)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,nprob
       real(dbl) :: x(n),fvec(n)
     end subroutine vecfcn

     subroutine vecjac(n,x,fjac,ldfjac,nprob)
       integer, parameter :: dbl = selected_real_kind(15)
       integer :: n,ldfjac,nprob
       real(dbl) :: x(n),fjac(ldfjac,n)
     end subroutine vecjac

  end interface

end module minpack
