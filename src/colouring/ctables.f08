!
! ctables manipulations
!
! Copyright © 2018 F.Hroch (hroch@physics.muni.cz)
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

module ctables

  use iso_fortran_env

  implicit none

contains

  subroutine ctable_load(ctable,cspace,outspace,cmatrix,verbose)

    character(len=*), intent(in) :: ctable,cspace,outspace
    real, dimension(:,:), allocatable, intent(out) :: cmatrix
    logical, intent(in) :: verbose

    character(len=80) :: label1, label2, msg
    integer :: n,stat,ndim,mdim

    open(1,file=ctable,status='old',iostat=stat,iomsg=msg)
    if( stat > 0 ) then
       write(error_unit,*) "Error: ",trim(msg)
       stop 'A color table file not found.'
    end if

    do
       read(1,*,end=90,err=666,iomsg=msg) label1, label2, ndim, mdim
       allocate(cmatrix(ndim,mdim))
       do n = 1, ndim
          read(1,*,end=90,err=666,iomsg=msg) cmatrix(n,:)
       end do
       if( label1 == cspace .and. outspace == label2 ) exit
       deallocate(cmatrix)
    end do
90  close(1)

    if( .not. allocated(cmatrix) ) stop 'Required an unknown colour space.'

    if( verbose ) then
       write(error_unit,*) 'Info for colour transformation matrix: ', &
            trim(label1),' -> ',trim(label2)
       do n = 1, ndim
          write(error_unit,'(a,5f10.4)') ' Info: ',cmatrix(n,:)
       end do
    end if

    return

666 continue
    write(error_unit,*) "Error: ",trim(msg)
    stop 'A color table file read error.'

  end subroutine ctable_load

  subroutine ctable_list(ctable)

    character(len=*), intent(in) :: ctable

    character(len=80) :: msg,cspace,outspace
    integer :: stat,ndim,mdim
    real, dimension(:,:),allocatable :: cmatrix

    open(1,file=ctable,status='old',iostat=stat,iomsg=msg)
    if( stat > 0 ) then
       write(error_unit,*) "Error: ",trim(msg)
       stop 'A color table file not found.'
    end if

    write(*,*) 'Available colour spaces:'
    write(*,'(a20,1x,a,1x,a)') 'Input','>','Output'

    do
       read(1,*,end=90,err=666,iomsg=msg) cspace,outspace, ndim, mdim
       write(*,'(a20,3x,a)') trim(cspace), trim(outspace)
       allocate(cmatrix(ndim,mdim))
       read(1,*,end=90,err=666,iomsg=msg) cmatrix
       deallocate(cmatrix)
    end do
90  close(1)
    return

666 continue
    write(error_unit,*) "Error: ",trim(msg)
    stop 'A color table file read error.'

  end subroutine ctable_list


end module ctables
