!
!  Likelihood, and related, with robust minimization functions
!
!  This is vector, well readable, variant. Unfotunatelly,
!  it is also significantly slower and memory consuming than
!  the non-vercor one.
!
!  Copyright © 2024 F.Hroch (hroch@physics.muni.cz)
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
!
!  There is no extra memory charge except hess() functions.


module rlike

  use iso_fortran_env
  implicit none
  private

  type, abstract, public :: robust_likelihood
     real(REAL64), dimension(:), pointer, private :: data
     real(REAL64), dimension(:), allocatable, private :: res, psi
   contains
     procedure, public :: init => initialication
     procedure(sqpf), public, deferred :: sqpf
     procedure(func), public, deferred :: func
     procedure(grad), public, deferred :: grad
     procedure(hess), public, deferred :: hess
     procedure(stderr), public, deferred :: stderr
     procedure(frac), public, deferred :: frac
  end type robust_likelihood

  abstract interface

     subroutine sqpf(this,t,s,f,g)
       import :: REAL64, robust_likelihood
       class(robust_likelihood) :: this
       real(REAL64), intent(in) :: t,s
       real(REAL64), intent(out) :: f
       real(REAL64), dimension(:), intent(out) :: g
     end subroutine sqpf

     real(REAL64) function func(this,t,s)
       import :: REAL64, robust_likelihood
       class(robust_likelihood) :: this
       real(REAL64), intent(in) :: t,s
     end function func

     subroutine grad(this,t,s,g)
       import :: REAL64, robust_likelihood
       class(robust_likelihood) :: this
       real(REAL64), intent(in) :: t,s
       real(REAL64), dimension(:), intent(out) :: g
     end subroutine grad

     subroutine hess(this,t,s,H)
       import :: REAL64, robust_likelihood
       class(robust_likelihood) :: this
       real(REAL64), intent(in) :: t,s
       real(REAL64), dimension(:,:), intent(out) :: H
     end subroutine hess

     real(REAL64) function stderr(this,t,s)
       import :: REAL64, robust_likelihood
       class(robust_likelihood) :: this
       real(REAL64), intent(in) :: t,s
     end function stderr

     real(REAL64) function frac(this,t,s)
       import :: REAL64, robust_likelihood
       class(robust_likelihood) :: this
       real(REAL64), intent(in) :: t,s
     end function frac

  end interface

  type, extends(robust_likelihood), public :: robust_likelihood_hampel
   contains
     procedure :: sqpf => hampel_sqpf
     procedure :: func => hampel_func
     procedure :: grad => hampel_grad
     procedure :: hess => hampel_hess
     procedure :: stderr => hampel_stderr
     procedure :: frac => hampel_frac
     final :: hampel_final
  end type robust_likelihood_hampel

  type, extends(robust_likelihood), public :: robust_likelihood_tukey
   contains
     procedure :: sqpf => tukey_sqpf
     procedure :: func => tukey_func
     procedure :: grad => tukey_grad
     procedure :: hess => tukey_hess
     procedure :: stderr => tukey_stderr
     procedure :: frac => tukey_frac
     final :: tukey_final
  end type robust_likelihood_tukey

  type, extends(robust_likelihood), public :: robust_likelihood_huber
   contains
     procedure :: sqpf => huber_sqpf
     procedure :: func => huber_func
     procedure :: grad => huber_grad
     procedure :: hess => huber_hess
     procedure :: stderr => huber_stderr
     procedure :: frac => huber_frac
     final :: huber_final
  end type robust_likelihood_huber

  type, extends(robust_likelihood), public :: nonrobust_likelihood_square
   contains
     procedure :: sqpf => square_sqpf
     procedure :: func => square_func
     procedure :: grad => square_grad
     procedure :: hess => square_hess
     procedure :: stderr => square_stderr
     procedure :: frac => square_frac
     final :: square_final
  end type nonrobust_likelihood_square


contains

  ! Smooth Hampel's function

  subroutine hampel_sqpf(this,t,s,f,g)

    use rfun

    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: f
    real(REAL64), dimension(:), intent(out) :: g
    integer :: n,i

    associate( data => this%data, res => this%res, psi => this%psi )
      n = size(data)
      res = (data - t)/s
      f = sum(ipshampel(res))/n + log(s)
      psi = pshampel(res)
      g = [ -sum(psi)/n, -sum(psi*res)/n + 1 ] / s
    end associate

    ! This vector version is significantly slower than
    ! the same procedure in rlike.f08 implemented via loops with -O2.
    ! I've no idea why.

  end subroutine hampel_sqpf

  function hampel_func(this,t,s) result(func)

    use rfun

    real(REAL64) :: func
    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    integer :: n

    associate( data => this%data )
      n = size(data)
      func = sum(ipshampel((data - t)/s))/n + log(s)
    end associate

  end function hampel_func

  subroutine hampel_grad(this,t,s,g)

    use rfun

    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: g

    integer :: n

    associate( data => this%data, res => this%res, psi => this%psi )
      n = size(res)
      res = (data - t) / s
      psi = pshampel(res)
      g = [ -sum(psi)/n, -sum(psi*res)/n + 1 ] / s
    end associate

  end subroutine hampel_grad

  subroutine hampel_hess(this,t,s,H)

    use rfun

    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H

    real(REAL64), dimension(:), allocatable :: dpsi
    integer :: n

    associate( data => this%data, res => this%res, psi => this%psi )

      n = size(res)
      allocate(dpsi(n))
      res = (data - t) / s
      psi = pshampel(res)
      dpsi = dpshampel(res)

      H(1,1) = sum(dpsi) / s**2 / n
      H(1,2) = (sum(dpsi*res)/s + sum(psi)/s**2) / n
      H(2,1) = H(1,2)
      H(2,2) = (sum(dpsi*res**2) + 2*sum(psi*res))/s**2/n - 1/s**2
    end associate

  end subroutine hampel_hess

  function hampel_stderr(this,t,s) result(stderr)

    use rfun

    real(REAL64) :: stderr
    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s

    real(REAL64) :: f2,df,f
    integer :: n

    associate( data => this%data, res => this%res )

      n = size(res)
      res = (data - t) / s
      f2 = sum(pshampel(res)**2)
      df = sum(dpshampel(res))

    end associate

    if( f2 > 0 .and. df > 0 ) then
       f = sqrt((n*f2)/(n-1))
       stderr = s*(f/df)
    else if( df > 0 ) then
       stderr = 0
    else
       error stop 'hampel_stderr(): NOT (f2 > 0 .and. df > 0)'
    end if

  end function hampel_stderr

  function hampel_frac(this,t,s) result(frac)

    use rfun

    real(REAL64) :: frac
    class(robust_likelihood_hampel) :: this
    real(REAL64), intent(in) :: t,s
    integer :: n

    associate( data => this%data )
      n = size(data)
      frac = sum(dpshampel((data - t)/s))/n
    end associate

  end function hampel_frac


  ! Tukey's function

  subroutine tukey_sqpf(this,t,s,f,g)

    use rfun

    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: f
    real(REAL64), dimension(:), intent(out) :: g
    integer :: n

    associate( data => this%data, res => this%res, psi => this%psi )

      n = size(data)
      res = (data - t)/s
      f = sum(itukey(res))/n + log(s)
      psi = tukey(res)
      g = [ -sum(psi)/n, -sum(psi*res)/n + 1 ] / s

    end associate

  end subroutine tukey_sqpf


  function tukey_func(this,t,s) result(func)

    use rfun

    real(REAL64) :: func
    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    integer :: n

    associate( data => this%data )
      n = size(data)
      func = sum(itukey((data - t)/s))/n + log(s)
    end associate

  end function tukey_func


  subroutine tukey_grad(this,t,s,g)

    use rfun

    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: g

    integer :: n

    associate( data => this%data, res => this%res, psi => this%psi )

      n = size(res)
      res = (data - t) / s
      psi = tukey(res)
      g = [ -sum(psi)/n, -sum(psi*res)/n + 1 ] / s

    end associate

  end subroutine tukey_grad

  subroutine tukey_hess(this,t,s,H)

    use rfun

    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H

    real(REAL64), dimension(:), allocatable :: dpsi
    integer :: n

    associate( data => this%data, res => this%res, psi => this%psi )

      n = size(res)
      allocate(dpsi(n))
      res = (data - t) / s
      psi = tukey(res)
      dpsi = dtukey(res)

      H(1,1) = sum(dpsi) / s**2 / n
      H(1,2) = (sum(dpsi*res)/s + sum(psi)/s**2) / n
      H(2,1) = H(1,2)
      H(2,2) = (sum(dpsi*res**2) + 2*sum(psi*res))/s**2/n - 1/s**2

    end associate

  end subroutine tukey_hess

  function tukey_stderr(this,t,s) result(stderr)

    use rfun

    real(REAL64) :: stderr
    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s

    real(REAL64) :: f2,df,f
    integer :: n

    associate( data => this%data, res => this%res )

      n = size(res)
      res = (data - t) / s

      f2 = sum(tukey(res)**2)
      df = sum(dtukey(res))

    end associate

    if( f2 > 0 .and. df > 0 ) then
       f = sqrt((n*f2)/(n-1))
       stderr = s*(f/df)
    else if( df > 0 ) then
       stderr = 0
    else
       error stop 'tukey_stderr(): NOT (f2 > 0 .and. df > 0)'
    end if


  end function tukey_stderr

  function tukey_frac(this,t,s) result(frac)

    use rfun

    real(REAL64) :: frac
    class(robust_likelihood_tukey) :: this
    real(REAL64), intent(in) :: t,s
    integer :: n

    associate( data => this%data )
      n = size(data)
      frac = sum(dtukey((data - t)/s))/n
    end associate

  end function tukey_frac


  ! Smooth Huber's function

  subroutine huber_sqpf(this,t,s,f,g)

    use rfun

    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: f
    real(REAL64), dimension(:), intent(out) :: g
    integer :: n

    associate( data => this%data, res => this%res, psi => this%psi )

      n = size(data)
      res = (data - t)/s
      f = sum(ipshuber(res))/n + log(s)
      psi = pshuber(res)
      g = [ -sum(psi)/n, -sum(psi*res)/n + 1 ] / s

    end associate

  end subroutine huber_sqpf

  function huber_func(this,t,s) result(func)

    use rfun

    real(REAL64) :: func
    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    integer :: n

    associate( data => this%data )
      n = size(data)
      func = sum(ipshuber((data - t)/s))/n + log(s)
    end associate

  end function huber_func

  subroutine huber_grad(this,t,s,g)

    use rfun

    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: g

    integer :: n

    associate( data => this%data, res => this%res, psi => this%psi )

      n = size(res)
      res = (data - t) / s
      psi = pshuber(res)
      g = [ -sum(psi)/n, -sum(psi*res)/n + 1 ] / s
    end associate

  end subroutine huber_grad

  subroutine huber_hess(this,t,s,H)

    use rfun

    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H

    real(REAL64), dimension(:), allocatable :: dpsi
    integer :: n

    associate( data => this%data, res => this%res, psi => this%psi )

      n = size(res)
      allocate(dpsi(n))
      res = (data - t) / s
      psi = pshuber(res)
      dpsi = dpshuber(res)

      H(1,1) = sum(dpsi) / s**2 / n
      H(1,2) = (sum(dpsi*res)/s + sum(psi)/s**2) / n
      H(2,1) = H(1,2)
      H(2,2) = (sum(dpsi*res**2) + 2*sum(psi*res))/s**2/n - 1/s**2
    end associate

  end subroutine huber_hess

  function huber_stderr(this,t,s) result(stderr)

    use rfun

    real(REAL64) :: stderr
    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s

    real(REAL64) :: f2,df,f
    integer :: n

    associate( data => this%data, res => this%res)

      n = size(data)
      res = (data - t) / s
      f2 = sum(pshuber(res)**2)
      df = sum(dpshuber(res))
    end associate

!    if( f2 > 0 .and. df > 0 ) then
    if( f2 > 0 ) then
       f = sqrt((n*f2)/(n-1))
       stderr = s*(f/df)
    else
       stderr = 0
!       error stop 'huber_stderr(): NOT (f2 > 0 .and. df > 0)'
    end if

  end function huber_stderr

  function huber_frac(this,t,s) result(frac)

    use rfun

    real(REAL64) :: frac
    class(robust_likelihood_huber) :: this
    real(REAL64), intent(in) :: t,s
    integer :: n

    associate( data => this%data )
      n = size(data)
      frac = sum(dpshuber((data - t)/s))/n
    end associate

  end function huber_frac


  ! Squares, non-robust

  subroutine square_sqpf(this,t,s,f,g)

    use rfun

    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), intent(out) :: f
    real(REAL64), dimension(:), intent(out) :: g
    integer :: n

    associate( data => this%data, res => this%res )
      n = size(data)
      res = (data - t)/s
      f = sum(res**2)/(2*n) + log(s)
      g = [ -sum(res)/n, -sum(res**2)/n + 1 ] / s
    end associate

  end subroutine square_sqpf


  function square_func(this,t,s) result(func)

    real(REAL64) :: func
    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    integer :: n

    associate( data => this%data, res => this%res )
      n = size(data)
      res = (data - t)/s
      func = sum(res**2)/(2*n) + log(s)
    end associate

  end function square_func

  subroutine square_grad(this,t,s,g)

    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:), intent(out) :: g

    integer :: n

    associate( data => this%data, res => this%res )

      n = size(res)
      res = (data - t) / s
      g = [ -sum(res)/n, -sum(res**2)/n + 1 ] / s
    end associate

  end subroutine square_grad

  subroutine square_hess(this,t,s,H)

    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s
    real(REAL64), dimension(:,:), intent(out) :: H

    real(REAL64), dimension(:), allocatable :: dpsi
    integer :: n

    associate( data => this%data, res => this%res, psi => this%psi )

      n = size(res)
      allocate(dpsi(n))
      res = (data - t) / s
      psi = res
      dpsi = 1

      H(1,1) = sum(dpsi) / s**2 / n
      H(1,2) = (sum(dpsi*res)/s + sum(psi)/s**2) / n
      H(2,1) = H(1,2)
      H(2,2) = (sum(dpsi*res**2) + 2*sum(psi*res))/s**2/n - 1/s**2
    end associate

  end subroutine square_hess

  function square_stderr(this,t,s) result(stderr)

    real(REAL64) :: stderr
    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s

    real(REAL64) :: f2,df,f
    integer :: n

    associate( data => this%data, res => this%res )

      n = size(res)
      res = (data - t) / s
      f2 = sum(res**2)
      df = n

    end associate

    if( f2 > 0 ) then
       f = sqrt((n*f2)/(n-1))
       stderr = s*(f/df)
    else
       stderr = 0
    end if

  end function square_stderr

  function square_frac(this,t,s) result(frac)

    real(REAL64) :: frac
    class(nonrobust_likelihood_square) :: this
    real(REAL64), intent(in) :: t,s

    frac = 1

  end function square_frac

  subroutine initialication(this,data,iflag)

    class(robust_likelihood) :: this
    real(REAL64), dimension(:), target, intent(in) :: data
    integer, intent(out) :: iflag

    integer :: n
    character(len=81) :: errmsg

    this%data => data
    n = size(data)
    if( .not. (n > 0) ) error stop 'rlike$init(): NOT (n > 0)'

    allocate(this%res(n),this%psi(n),stat=iflag,errmsg=errmsg)
    if( iflag /= 0 ) then
       write(error_unit,*) 'rlike%init(): ',trim(errmsg)
       return
    end if

  end subroutine initialication

  subroutine hampel_final(this)

    type(robust_likelihood_hampel) :: this

!    write(*,*) ' Hampel ............................... Finished.',allocated(this%res), allocated(this%psi)
    if( allocated(this%res) ) deallocate(this%res)
    if( allocated(this%psi) ) deallocate(this%psi)

  end subroutine hampel_final

  subroutine tukey_final(this)

    type(robust_likelihood_tukey) :: this

!    write(*,*) ' Tukey ............................... Finished.',allocated(this%res), allocated(this%psi)
    if( allocated(this%res) ) deallocate(this%res)
    if( allocated(this%psi) ) deallocate(this%psi)

  end subroutine tukey_final

  subroutine huber_final(this)

    type(robust_likelihood_huber) :: this

!    write(*,*) ' Huber ............................... Finished.',allocated(this%res), allocated(this%psi)
    if( allocated(this%res) ) deallocate(this%res)
    if( allocated(this%psi) ) deallocate(this%psi)

  end subroutine huber_final

  subroutine square_final(this)

    type(nonrobust_likelihood_square) :: this

!    write(*,*) ' Square .............................. Finished.',allocated(this%res), allocated(this%psi)
    if( allocated(this%res) ) deallocate(this%res)
    if( allocated(this%psi) ) deallocate(this%psi)

  end subroutine square_final


end module rlike
