!
!  An implementation of the decomposition of a 2x2 matrix on
!  eigenvalues and eigenvectors by Schur complement method.
!
!  Copyright © 2023-2026 F.Hroch (hroch@physics.muni.cz)
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
!  The implementation is inspired by Golub & and Van Loan book:
!  Matrix Computations, 2013 edition, Sect. 8.5.2.

module modQ2eig_REAL64

  use iso_fortran_env

  implicit none
  private

  logical, parameter :: debug = .false. ! = .false. for the release !
  logical, parameter :: debug_print = .false.

  public :: q2eig_REAL64

contains

    subroutine q2eig_REAL64(H,lam,Q,nroots)

    ! Determination of eigenvalues by solution of the quadratic equation:
    !
    ! lam**2 - [H(1,1) + H(2,2)]*lam + H(1,1)*H(2,2) - H(1,2)*H(2,1) = 0
    !
    ! and corresponding eigenvectors by Schur complemet method.

    real(REAL64), dimension(:,:), intent(in) :: H
    real(REAL64), dimension(:), intent(out) :: lam
    real(REAL64), dimension(:,:), intent(out) :: Q
    integer, intent(out) :: nroots
    real(REAL64) :: b,c,d,p,r,s,t,x,y,u,v,l,absy

    if( debug_print ) then
       write(error_unit,*) 'q2eig H(1,:):',H(1,:)
       write(error_unit,*) 'q2eig H(2,:):',H(2,:)
    end if

    if( debug ) then
       if( abs(H(1,2) - H(2,1)) > epsilon(H) ) then
          write(error_unit,*) 'H:',H, 'res:',H(1,2)-H(2,1)
          error stop 'Symmetric matrix H(1,2) == H(2,1) assumption violated.'
       end if
    end if

    x = H(2,2) - H(1,1)
    y = H(1,2) + H(2,1)
    b = -(H(1,1) + H(2,2))
    d = hypot(x,y)
    if( d > 0 ) then
       c = H(1,1)*H(2,2) - H(1,2)*H(2,1)

       ! eigenvalues
       lam(2) = -(b + sign(d,b)) / 2
       lam(1) = c / lam(2)

       if( debug_print ) write(error_unit,'(a,2(1x,g0),a,2es10.2)') &
            'quadra eigen.vals:',lam,' res:',lam**2+b*lam+c

       ! eigenvectors
       l = lam(2) - lam(1)
       absy = abs(y)
       if( absy > 0 ) then
          r = x / y
          p = d / absy
          u = x / l
          v = y / l
          if( u >= 0 ) then
             t = r + sign(p,v)
          else
             t = sign(1 / (abs(r) + p), v)
          end if
          s = sign(1 / hypot(1.0_REAL64,t), v)
          c = t*s
       else
          if( x / l > 0  ) then
             c = 1
             s = 0
          else
             c = 0
             s = 1
          end if
       end if
       nroots = 2
    else !if ( abs(d) < epsilon(d) ) then
       ! double root
       lam = -b / 2
       c = 1
       s = 0
       nroots = 1
!    else ! if( d < -macheps ) then
       ! complex roots, weird for H symmetric
!       error stop 'q2eig() requires a symmetric real matrix.'
    end if
    Q(:,1) = [c,-s]
    Q(:,2) = [s, c]

    ! results check
    if( debug ) then
       block
         integer :: i
         real(REAL64) :: c
         c = H(1,1)*H(2,2) - H(1,2)*H(2,1)
         if( debug_print ) then
            write(error_unit,*) 'Q',real(Q(1,:))
            write(error_unit,*) 'Q',real(Q(2,:))
            write(error_unit,'(a,2(1x,es10.2))') 'quadra eigen vectors residuals:',&
                 (norm2(matmul(H,Q(:,i))-lam(i)*Q(:,i)),i=1,2)
         end if
         do i = 1,2
            r = lam(i)**2 + b*lam(i) + c
            if( abs(r) > 1e3*(1+sum(lam**2))*epsilon(lam) ) then
               write(error_unit,*) i,lam(i),b,c,abs(r),1e3*(1+lam(i)**2)*epsilon(lam)
               error stop 'Eigen values check failed.'
            end if
         end do

         do i = 1,2
            r = norm2(matmul(H,Q(:,i))-lam(i)*Q(:,i))
            if( r > 1e3*(1+maxval(abs(lam)))*epsilon(Q) ) then
               write(error_unit,*) i,lam,':',Q,' res:',r
               error stop 'Eigen vectors check failed.'
            end if
         end do
       end block

       if( abs(s) > 1 .or. abs(c) > 1 ) then
          write(error_unit,*) 'cos = ',c,' sin = ',s
          error stop 'Assumption |cos| <= 1 & |sin| <= 1 violated.'
       end if

    end if

  end subroutine q2eig_REAL64

end module modQ2eig_REAL64
