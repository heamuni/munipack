!
! This module is a adopted version of algorithm described below.
!
! Download:    http://people.sc.fsu.edu/~burkardt/f_src/asa047/asa047.html
!
! Thx to J. Burkardt http://people.sc.fsu.edu/~burkardt/index.html
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


Module NelderMead

  implicit none
  private

  public :: nelmin, nelmin1

contains

  ! a simple wrapper to the original subroutine

  subroutine nelmin1(fn,xmin,step,ynewlo,ifault, &
       reqmin0,konvge0,kcount0,icount0,numres0)

    interface
       function fn(t)
         real(kind=8) :: fn
         real(kind=8), dimension(:), intent(in) :: t
       end function fn
    end interface

    real(kind=8), dimension(:), intent(out) :: xmin
    real(kind=8), dimension(:), intent(in) :: step
    real(kind=8), intent(out) :: ynewlo
    integer, intent(out) :: ifault
    real(kind=8), optional, intent(in) :: reqmin0
    integer, optional, intent(in) :: konvge0, kcount0
    integer, optional, intent(out) :: icount0,numres0

    real(kind=8), dimension(size(xmin)) :: start
    integer :: n,konvge,kcount,icount,numres
    real(kind=8) :: reqmin

    n = size(xmin)
    start = xmin

    if( present(reqmin0) ) then
       reqmin = reqmin0
    else
       reqmin = epsilon(ynewlo)
    end if

    if( present(konvge0) ) then
       konvge = konvge0
    else
       konvge = 1
    end if

    if( present(kcount0) ) then
       kcount = kcount0
    else
       kcount = 1000*n
    end if

    call nelmin(fn,n,start,xmin,ynewlo,reqmin,step,konvge,kcount,icount,numres,ifault)

    if( present(icount0) ) icount0 = icount
    if( present(numres0) ) numres0 = numres

  end subroutine nelmin1


  subroutine nelmin ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
  icount, numres, ifault )

!*****************************************************************************80
!
!! NELMIN minimizes a function using the Nelder-Mead algorithm.
!
!  Discussion:
!
!    This routine seeks the minimum value of a user-specified function.
!
!    Simplex function minimisation procedure due to Nelder and Mead (1965),
!    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
!    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
!    25, 97) and Hill(1978, 27, 380-2)
!
!    The function to be minimized must be defined by a function of
!    the form
!
!      function fn ( x, f )
!      real ( kind = 8 ) fn
!      real ( kind = 8 ) x(*)
!
!    and the name of this subroutine must be declared EXTERNAL in the
!    calling routine and passed as the argument FN.
!
!    This routine does not include a termination test using the
!    fitting of a quadratic surface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by R ONeill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Nelder, Roger Mead,
!    A simplex method for function minimization,
!    Computer Journal,
!    Volume 7, 1965, pages 308-313.
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, external FN, the name of the function which evaluates
!    the function to be minimized.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    0 < N is required.
!
!    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
!    for the iteration.  On output, this data may have been overwritten.
!
!    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
!    is estimated to minimize the function.
!
!    Output, real ( kind = 8 ) YNEWLO, the minimum value of the function.
!
!    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
!    of the function values.  0 < REQMIN is required.
!
!    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
!    initial simplex.  The relative magnitudes of its elements should reflect
!    the units of the variables.
!
!    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
!    every KONVGE iterations. 0 < KONVGE is required.
!
!    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
!    evaluations.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
!    used.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no errors detected.
!    1, REQMIN, N, or KONVGE has an illegal value.
!    2, iteration terminated because KCOUNT was exceeded without convergence.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: ccoeff = 0.5D+00
  real ( kind = 8 ) del
  real ( kind = 8 ), parameter :: ecoeff = 2.0D+00
  real ( kind = 8 ), parameter :: eps = 0.001D+00
!  real ( kind = 8 ), external :: fn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcount
  integer ( kind = 4 ) kcount
  integer ( kind = 4 ) konvge
  integer ( kind = 4 ) l
  integer ( kind = 4 ) numres
  real ( kind = 8 ) p(n,n+1)
  real ( kind = 8 ) p2star(n)
  real ( kind = 8 ) pbar(n)
  real ( kind = 8 ) pstar(n)
  real ( kind = 8 ), parameter :: rcoeff = 1.0D+00
  real ( kind = 8 ) reqmin
  real ( kind = 8 ) rq
  real ( kind = 8 ) start(n)
  real ( kind = 8 ) step(n)
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin(n)
  real ( kind = 8 ) y(n+1)
  real ( kind = 8 ) y2star
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ynewlo
  real ( kind = 8 ) ystar
  real ( kind = 8 ) z

  interface
     function fn(t)
       real(kind=8) :: fn
       real(kind=8), dimension(:), intent(in) :: t
     end function fn
  end interface


!
!  Check the input parameters.
!
  if ( reqmin <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  if ( n < 1 ) then
    ifault = 1
    return
  end if

  if ( konvge < 1 ) then
    ifault = 1
    return
  end if
!
!  Initialization.
!
  icount = 0
  numres = 0
  jcount = konvge
  del = 1.0D+00
  rq = reqmin * real ( n, kind = 8 )
!
!  Initial or restarted loop.
!
  do

    p(1:n,n+1) = start(1:n)
    y(n+1) = fn ( start )
    icount = icount + 1
!
!  Define the initial simplex.
!
    do j = 1, n
      x = start(j)
      start(j) = start(j) + step(j) * del
      p(1:n,j) = start(1:n)
      y(j) = fn ( start )
      icount = icount + 1
      start(j) = x
    end do
!
!  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
!  the vertex of the simplex to be replaced.
!
    ilo = minloc ( y(1:n+1), 1 )
    ylo = y(ilo)
!
!  Inner loop.
!
    do while ( icount < kcount )
!
!  YNEWLO is, of course, the HIGHEST value???
!
      ihi = maxloc ( y(1:n+1), 1 )
      ynewlo = y(ihi)
!
!  Calculate PBAR, the centroid of the simplex vertices
!  excepting the vertex with Y value YNEWLO.
!
      do i = 1, n
        pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = 8 )
      end do
!
!  Reflection through the centroid.
!
      pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
      ystar = fn ( pstar )
      icount = icount + 1
!
!  Successful reflection, so extension.
!
      if ( ystar < ylo ) then

        p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
        y2star = fn ( p2star )
        icount = icount + 1
!
!  Retain extension or contraction.
!
        if ( ystar < y2star ) then
          p(1:n,ihi) = pstar(1:n)
          y(ihi) = ystar
        else
          p(1:n,ihi) = p2star(1:n)
          y(ihi) = y2star
        end if
!
!  No extension.
!
      else

        l = 0
        do i = 1, n + 1
          if ( ystar < y(i) ) then
            l = l + 1
          end if
        end do

        if ( 1 < l ) then

          p(1:n,ihi) = pstar(1:n)
          y(ihi) = ystar
!
!  Contraction on the Y(IHI) side of the centroid.
!
        else if ( l == 0 ) then

          p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
          y2star = fn ( p2star )
          icount = icount + 1
!
!  Contract the whole simplex.
!
          if ( y(ihi) < y2star ) then

            do j = 1, n + 1
              p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
              xmin(1:n) = p(1:n,j)
              y(j) = fn ( xmin )
              icount = icount + 1
            end do

            ilo = minloc ( y(1:n+1), 1 )
            ylo = y(ilo)

            cycle
!
!  Retain contraction.
!
          else
            p(1:n,ihi) = p2star(1:n)
            y(ihi) = y2star
          end if
!
!  Contraction on the reflection side of the centroid.
!
        else if ( l == 1 ) then

          p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
          y2star = fn ( p2star )
          icount = icount + 1
!
!  Retain reflection?
!
          if ( y2star <= ystar ) then
            p(1:n,ihi) = p2star(1:n)
            y(ihi) = y2star
          else
            p(1:n,ihi) = pstar(1:n)
            y(ihi) = ystar
          end if

        end if

      end if
!
!  Check if YLO improved.
!
      if ( y(ihi) < ylo ) then
        ylo = y(ihi)
        ilo = ihi
      end if

      jcount = jcount - 1

      if ( 0 < jcount ) then
        cycle
      end if
!
!  Check to see if minimum reached.
!
      if ( icount <= kcount ) then

        jcount = konvge

        x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
        z = sum ( ( y(1:n+1) - x )**2 )

        if ( z <= rq ) then
          exit
        end if

      end if

    end do
!
!  Factorial tests to check that YNEWLO is a local minimum.
!
    xmin(1:n) = p(1:n,ilo)
    ynewlo = y(ilo)

    if ( kcount < icount ) then
      ifault = 2
      exit
    end if

    ifault = 0

    do i = 1, n
      del = step(i) * eps
      xmin(i) = xmin(i) + del
      z = fn ( xmin )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      end if
      xmin(i) = xmin(i) - del - del
      z = fn ( xmin )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      end if
      xmin(i) = xmin(i) + del
    end do

    if ( ifault == 0 ) then
      exit
    end if
!
!  Restart the procedure.
!
    start(1:n) = xmin(1:n)
    del = eps
    numres = numres + 1

  end do

  return
end subroutine nelmin

end Module NelderMead
