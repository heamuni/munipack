!
!  Likelihood, and related, with robust minimization functions
!
!  This is do-loop version. There is also a vectorised version.
!
!  Copyright © 2025-26 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Oakleaf.
!
!  Oakleaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Oakleaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!


module rlike_REAL64

  use rfun
  use iso_fortran_env
  use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN
  implicit none
  private

  type, abstract, public :: robust_likelihood_REAL64
     real(REAL64), dimension(:), pointer :: data
   contains
     procedure(sqpf), public, deferred :: sqpf
     procedure(func), public, deferred :: func
     procedure(grad), public, deferred :: grad
     procedure(hess), public, deferred :: hess
     procedure(minhess), public, deferred :: minhess
     procedure(stderr), public, deferred :: stderr
     procedure(fmax), public, deferred, nopass :: fmax
  end type robust_likelihood_REAL64

  abstract interface

     subroutine sqpf(this,t,s,f,r,g)
       import :: REAL64, robust_likelihood_REAL64
       class(robust_likelihood_REAL64) :: this
       real(REAL64), intent(in) :: t,s
       real(REAL64), intent(out) :: f,r
       real(REAL64), dimension(:), intent(out) :: g
     end subroutine sqpf

     real(REAL64) function func(this,t,s)
       import :: REAL64, robust_likelihood_REAL64
       class(robust_likelihood_REAL64) :: this
       real(REAL64), intent(in) :: t,s
     end function func

     subroutine grad(this,t,s,g)
       import :: REAL64, robust_likelihood_REAL64
       class(robust_likelihood_REAL64) :: this
       real(REAL64), intent(in) :: t,s
       real(REAL64), dimension(:), intent(out) :: g
     end subroutine grad

     subroutine hess(this,t,s,H)
       import :: REAL64, robust_likelihood_REAL64
       class(robust_likelihood_REAL64) :: this
       real(REAL64), intent(in) :: t,s
       real(REAL64), dimension(:,:), intent(out) :: H
     end subroutine hess

     subroutine minhess(this,t,s,H)
       import :: REAL64, robust_likelihood_REAL64
       class(robust_likelihood_REAL64) :: this
       real(REAL64), intent(in) :: t,s
       real(REAL64), dimension(:,:), intent(out) :: H
     end subroutine minhess

     real(REAL64) function stderr(this,t,s)
       import :: REAL64, robust_likelihood_REAL64
       class(robust_likelihood_REAL64) :: this
       real(REAL64), intent(in) :: t,s
     end function stderr

     real(REAL64) function fmax()
       import :: REAL64
     end function fmax

  end interface

  type, extends(robust_likelihood_REAL64), public :: robust_likelihood_hampel
   contains
     procedure :: sqpf => hampel_sqpf
     procedure :: func => hampel_func
     procedure :: grad => hampel_grad
     procedure :: hess => hampel_hess
     procedure :: minhess => hampel_minhess
     procedure :: stderr => hampel_stderr
     procedure, nopass :: fmax => hampel_fmax
  end type robust_likelihood_hampel

  type, extends(robust_likelihood_REAL64), public :: robust_likelihood_tukey
   contains
     procedure :: sqpf => tukey_sqpf
     procedure :: func => tukey_func
     procedure :: grad => tukey_grad
     procedure :: hess => tukey_hess
     procedure :: minhess => tukey_minhess
     procedure :: stderr => tukey_stderr
     procedure, nopass :: fmax => tukey_fmax
  end type robust_likelihood_tukey

  type, extends(robust_likelihood_REAL64), public :: robust_likelihood_huber
   contains
     procedure :: sqpf => huber_sqpf
     procedure :: func => huber_func
     procedure :: grad => huber_grad
     procedure :: hess => huber_hess
     procedure :: minhess => huber_minhess
     procedure :: stderr => huber_stderr
     procedure, nopass :: fmax => huber_fmax
  end type robust_likelihood_huber

  type, extends(robust_likelihood_REAL64), public :: nonrobust_likelihood_square
   contains
     procedure :: sqpf => square_sqpf
     procedure :: func => square_func
     procedure :: grad => square_grad
     procedure :: hess => square_hess
     procedure :: minhess => square_minhess
     procedure :: stderr => square_stderr
     procedure, nopass :: fmax => square_fmax
  end type nonrobust_likelihood_square

!  logical, parameter :: kahan = .false.

contains

  ! Smooth Hampel's function

  subroutine hampel_sqpf(this,t,s,f,r,g)

    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: f,r
    real(REAL64), dimension(:), intent(out) :: g
    real(REAL64) :: res, psi, rho, s0, gt, gs, s1
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      s0 = 0
      gt = 0
      gs = 0
      s1 = 1 / s
      do i = 1, n
         res = (data(i) - t) * s1
         call hampelsqpf(res,rho,psi)
         s0 = s0 + rho
         if( .not. (psi == 0) ) then
            ! This test is essential for proper handling of infinities and NaNs.
            ! 0*infinity is undefined, but robust function properties guided us
            ! to be zero.
            gt = gt + psi
            gs = gs + psi*res
         end if
      end do
      r = s0 / n
      f = r + log(s)
      g = [ -gt, n - gs ] / (n*s)
    end associate

!!$    if( kahan ) then
!!$       ! https://en.wikipedia.org/wiki/Kahan_summation_algorithm
!!$       block
!!$         real(REAL64) :: s2,tt,y
!!$         real(REAL64), dimension(:), allocatable :: v
!!$         associate( data => this%data )
!!$           n = size(data)
!!$           allocate(v(n))
!!$           s0 = 0
!!$           s1 = 1 / s
!!$           s2 = 0
!!$           do i = 1, n
!!$              res = (data(i) - t)*s1
!!$              y = ishampel(res)
!!$              v(i) = y
!!$              s2 = s2 + y
!!$              tt = s0 + s2
!!$              s2 = (s0 - tt) + s2
!!$              s0 = tt
!!$           end do
!!$           write(error_unit,*) 'Hampel sum by Kahan:',s0/n,r,r - s0/n, &
!!$                sum(v)/n,sum(v)/n-s0/n,sum(v)/n-r
!!$         end associate
!!$       end block
!!$    end if

  end subroutine hampel_sqpf

  function hampel_func(this,t,s) result(func)

    real(REAL64) :: func
    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64) :: s0, s1
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      s0 = 0
      s1 = 1 / s
      do i = 1, n
         s0 = s0 + ipshampel((data(i) - t)*s1)
      end do
      func = s0 / n + log(s)
    end associate

  end function hampel_func

  subroutine hampel_grad(this,t,s,g)

    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: g
    real(REAL64) :: res,psi,gt,gs,s1
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      gt = 0
      gs = 0
      s1 = 1 / s
      do i = 1, n
         res = (data(i) - t) * s1
         psi = pshampel(res)
         if( .not. (psi == 0) ) then
            gt = gt + psi
            gs = gs + psi*res
         end if
      end do
      g = [ -gt/n, -gs/n + 1 ] / s
    end associate

  end subroutine hampel_grad

  subroutine hampel_hess(this,t,s,H)

    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H
    real(REAL64) :: res, psi, dpsi, htt, hts, hss1, hss2, w, dpsires, s1
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      htt = 0
      hts = 0
      hss1 = 0
      hss2 = 0
      s1 = 1 / s
      do i = 1, n
         res = (data(i) - t) * s1
         psi = pshampel(res)
         dpsi = dpshampel(res)
         if( psi /= 0 .and. dpsi /= 0 ) then
            dpsires = dpsi*res
            htt = htt + dpsi
            hts = hts + dpsires + psi
            hss1= hss1+ dpsires*res
            hss2= hss2+ psi*res
         end if
      end do
      w = n*s**2
      H(1,1) = htt / w
      H(1,2) = hts / w
      H(2,1) = H(1,2)
      H(2,2) = (hss1 + 2*hss2) / w - 1/s**2
    end associate

  end subroutine hampel_hess

  subroutine hampel_minhess(this,t,s,H)

    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H
    real(REAL64) :: res, dpsi, htt, hts, hss, w, dpsires, s1
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      htt = 0
      hts = 0
      hss = 0
      s1 = 1 / s
      do i = 1, n
         res = (data(i) - t) * s1
         dpsi = dpshampel(res)
         if( .not. (dpsi == 0) ) then
            dpsires = dpsi*res
            htt = htt + dpsi
            hts = hts + dpsires
            hss = hss + dpsires*res
         end if
      end do
      w = 1 / (n*s**2)
      H(1,1) = htt * w
      H(1,2) = hts * w
      H(2,1) = H(1,2)
      H(2,2) = (n + hss) * w
    end associate

  end subroutine hampel_minhess


  function hampel_stderr(this,t,s) result(stderr)


    real(REAL64) :: stderr
    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64) :: f2,df,f,res
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      f2 = 0
      df = 0
      do i = 1, n
         res = (data(i) - t) / s
         f2 = f2 + pshampel(res)**2
         df = df + dpshampel(res)
      end do
    end associate

    if( f2 > 0 .and. df > 0 ) then
       f = sqrt(n/(n-1.0) * f2)
       stderr = s*(f/df)
    else if( df > 0 ) then
       stderr = 0
    else
       stderr = IEEE_VALUE(stderr, IEEE_QUIET_NAN)
    end if

  end function hampel_stderr

  function hampel_fmax()

    real(REAL64) :: hampel_fmax

    hampel_fmax = ipshampel(huge(1.0_REAL64))

  end function hampel_fmax

  ! Tukey's function

  subroutine tukey_sqpf(this,t,s,f,r,g)

    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: f,r
    real(REAL64), dimension(:), intent(out) :: g
    real(REAL64) :: res, rho, psi, s0, gt, gs, s1
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      s0 = 0
      gt = 0
      gs = 0
      s1 = 1 / s
      do i = 1, n
         res = (data(i) - t) * s1
         call tukeysqpf(res,rho,psi)
         s0 = s0 + rho
         if( .not. (psi == 0) ) then
            gt = gt + psi
            gs = gs + psi*res
         end if
      end do
      r = s0 / n
      f = r + log(s)
      g = [ -gt/n, -gs/n + 1 ] / s
    end associate

  end subroutine tukey_sqpf


  function tukey_func(this,t,s) result(func)

    real(REAL64) :: func
    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64) :: s0,s1
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      s0 = 0
      s1 = 1 / s
      do i = 1, n
         s0 = s0 + itukey((data(i) - t)*s1)
      end do
      func = s0 / n + log(s)
    end associate

  end function tukey_func

  subroutine tukey_grad(this,t,s,g)

    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: g
    real(REAL64) :: res, psi, gt, gs, s1
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      gt = 0
      gs = 0
      s1 = 1 / s
      do i = 1, n
         res = (data(i) - t) * s1
         psi = tukey(res)
         if( .not. (psi == 0) ) then
            gt = gt + psi
            gs = gs + psi*res
         end if
      end do
      g = [ -gt/n, -gs/n + 1 ] / s
    end associate

  end subroutine tukey_grad

  subroutine tukey_hess(this,t,s,H)

    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H
    real(REAL64) :: res, psi, dpsi, htt, hts, hss1, hss2, w, dpsires
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      htt = 0
      hts = 0
      hss1 = 0
      hss2 = 0
      do i = 1, n
         res = (data(i) - t) / s
         psi = tukey(res)
         dpsi = dtukey(res)
         if( psi /= 0 .and. dpsi /= 0 ) then
            dpsires = dpsi*res
            htt = htt + dpsi
            hts = hts + dpsires + psi
            hss1= hss1+ dpsires*res
            hss2= hss2+ psi*res
         end if
      end do
      w = n*s**2
      H(1,1) = htt / w
      H(1,2) = hts / w
      H(2,1) = H(1,2)
      H(2,2) = (hss1 + 2*hss2) / w - 1/s**2
    end associate

  end subroutine tukey_hess

  subroutine tukey_minhess(this,t,s,H)

    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H
    real(REAL64) :: res, dpsi, htt, hts, hss, w, dpsires, s1
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      htt = 0
      hts = 0
      hss = 0
      s1 = 1 / s
      do i = 1, n
         res = (data(i) - t) * s1
         dpsi = dtukey(res)
         if( .not. (dpsi == 0) ) then
            dpsires = dpsi*res
            htt = htt + dpsi
            hts = hts + dpsires
            hss = hss + dpsires*res
         end if
      end do
      w = n*s**2
      H(1,1) = htt / w
      H(1,2) = hts / w
      H(2,1) = H(1,2)
      H(2,2) = (n + hss) / w
    end associate

  end subroutine tukey_minhess

  function tukey_stderr(this,t,s) result(stderr)

    real(REAL64) :: stderr
    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64) :: f2,df,f,res
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      f2 = 0
      df = 0
      do i = 1, n
         res = (data(i) - t) / s
         f2 = f2 + tukey(res)**2
         df = df + dtukey(res)
      end do
    end associate

    if( f2 > 0 .and. df > 0 ) then
       f = sqrt((n*f2)/(n-1))
       stderr = s*(f/df)
    else if( df > 0 ) then
       stderr = 0
    else
       stderr = IEEE_VALUE(stderr, IEEE_QUIET_NAN)
    end if


  end function tukey_stderr

  function tukey_fmax() result(fmax)

    real(REAL64) :: fmax

    fmax = itukey(huge(1.0_REAL64))

  end function tukey_fmax



  ! Smooth Huber's function

  subroutine huber_sqpf(this,t,s,f,r,g)

    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: f,r
    real(REAL64), dimension(:), intent(out) :: g
    real(REAL64) :: res, rho, psi, s0, gt, gs, s1
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      s0 = 0
      gt = 0
      gs = 0
      s1 = 1 / s
      do i = 1, n
         res = (data(i) - t)*s1
         call hubersqpf(res,rho,psi)
         s0 = s0 + rho
         gt = gt + psi
         gs = gs + psi*res
      end do
      r = s0 / n
      f = r + log(s)
      g = [ -gt/n, -gs/n + 1 ] / s
    end associate

  end subroutine huber_sqpf

  function huber_func(this,t,s) result(func)

    real(REAL64) :: func
    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64) :: s0, s1
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      s0 = 0
      s1 = 1 / s
      do i = 1, n
         s0 = s0 + ipshuber((data(i) - t)*s1)
      end do
      func = s0 / n + log(s)
    end associate

  end function huber_func

  subroutine huber_grad(this,t,s,g)

    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: g
    real(REAL64) :: res,psi,gt,gs
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      gt = 0
      gs = 0
      do i = 1, n
         res = (data(i) - t) / s
         psi = pshuber(res)
         gt = gt + psi
         gs = gs + psi*res
      end do
      g = [ -gt/n, -gs/n + 1 ] / s
    end associate

  end subroutine huber_grad

  subroutine huber_hess(this,t,s,H)

    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H
    real(REAL64) :: res, psi, dpsi, htt, hts, hss1, hss2, w, dpsires
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      htt = 0
      hts = 0
      hss1 = 0
      hss2 = 0
      do i = 1, n
         res = (data(i) - t) / s
         psi = pshuber(res)
         dpsi = dpshuber(res)
         dpsires = dpsi*res
         htt = htt + dpsi
         hts = hts + dpsires + psi
         hss1= hss1+ dpsires*res
         hss2= hss2+ psi*res
      end do
      w = n*s**2
      H(1,1) = htt / w
      H(1,2) = hts / w
      H(2,1) = H(1,2)
      H(2,2) = (hss1 + 2*hss2) / w - 1/s**2
    end associate

  end subroutine huber_hess

  subroutine huber_minhess(this,t,s,H)

    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H
    real(REAL64) :: res, dpsi, htt, hts, hss, w, dpsires
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      htt = 0
      hts = 0
      hss = 0
      do i = 1, n
         res = (data(i) - t) / s
         dpsi = dpshuber(res)
         dpsires = dpsi*res
         htt = htt + dpsi
         hts = hts + dpsires
         hss = hss + dpsires*res
      end do
      w = n*s**2
      H(1,1) = htt / w
      H(1,2) = hts / w
      H(2,1) = H(1,2)
      H(2,2) = (n + hss) / w
    end associate

  end subroutine huber_minhess

  function huber_stderr(this,t,s) result(stderr)

    real(REAL64) :: stderr
    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64) :: f2,df,f,res
    integer :: n, i

    associate( data => this%data)
      n = size(data)
      f2 = 0
      df = 0
      do i = 1, n
         res = (data(i) - t) / s
         f2 = f2 + pshuber(res)**2
         df = df + dpshuber(res)
      end do
    end associate

    if( f2 > 0 .and. df > 0 ) then
       f = sqrt((n*f2)/(n-1))
       stderr = s*(f/df)
    else if( df > 0 ) then
       stderr = 0
    else
       stderr = IEEE_VALUE(stderr, IEEE_QUIET_NAN)
    end if

  end function huber_stderr

  function huber_fmax() result(fmax)

    real(REAL64) :: fmax

    fmax = huge(fmax)

  end function huber_fmax

  ! Squares, non-robust

  subroutine square_sqpf(this,t,s,f,r,g)

    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: f, r
    real(REAL64), dimension(:), intent(out) :: g
    real(REAL64) :: res, s0, gt, s1
    integer :: n,i

    associate( data => this%data )
      n = size(data)
      s0 = 0
      gt = 0
      s1 = 1 / s
      do i = 1, n
         res = (data(i) - t)*s1
         s0 = s0 + res**2
         gt = gt + res
      end do
      r = s0 / (2*n)
      f = r + log(s)
      g = [ -gt/n, -s0/n + 1 ] / s
    end associate

  end subroutine square_sqpf


  function square_func(this,t,s) result(func)

    real(REAL64) :: func
    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64) :: s0,s1
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      s0 = 0
      s1 = 1 / s
      do i = 1, n
         s0 = s0 + ((data(i) - t)*s1)**2
      end do
      func = s0 / (2*n) + log(s)
    end associate

  end function square_func

  subroutine square_grad(this,t,s,g)

    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: g
    real(REAL64) :: res,gt,gs
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      gt = 0
      gs = 0
      do i = 1, n
         res = (data(i) - t) / s
         gt = gt + res
         gs = gs + res**2
      end do
      g = [ -gt/n, -gs/n + 1 ] / s
    end associate

  end subroutine square_grad

  subroutine square_hess(this,t,s,H)

    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H
    real(REAL64) :: res, hts, hss, w
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      hts = 0
      hss = 0
      do i = 1, n
         res = (data(i) - t) / s
         hts = hts + res
         hss = hss + res**2
      end do
      w = n*s**2
      H(1,1) = 1 / s**2
      H(1,2) = 2*hts / w
      H(2,1) = H(1,2)
      H(2,2) = 3*hss / w - 1/s**2
    end associate

  end subroutine square_hess

  subroutine square_minhess(this,t,s,H)

    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H

    associate( data => this%data )
      H(1,1) = 1 / s**2
      H(1,2) = 0
      H(2,1) = 0
      H(2,2) = 2 / s**2
    end associate

  end subroutine square_minhess

  function square_stderr(this,t,s) result(stderr)

    real(REAL64) :: stderr
    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64) :: f2,df,f
    integer :: n, i

    associate( data => this%data )
      n = size(data)
      f2 = 0
      do i = 1, n
         f2 = f2 + ((data(i) - t) / s)**2
      end do
    end associate

    df = n
    if( f2 > 0 ) then
       f = sqrt((n*f2)/(n-1))
       stderr = s*(f/df)
    else
       stderr = 0
    end if

  end function square_stderr

  function square_fmax() result(fmax)

    real(REAL64) :: fmax

    fmax = huge(fmax)

  end function square_fmax

end module rlike_REAL64
