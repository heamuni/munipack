!
!  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
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
!

module fitskeys

  implicit none

contains

  subroutine include_keywords(fitskeyfile,keys_include,cards_include)

    use titsio
    use iso_fortran_env

    character(len=*), intent(in) :: fitskeyfile
    character(len=*), dimension(:), allocatable, intent(out) :: keys_include, &
         cards_include

    character(len=FLEN_CARD) :: line, card, iomsg
    character(len=FLEN_KEYWORD) :: key
    character(len=FLEN_CARD), dimension(666) :: cards
    character(len=FLEN_KEYWORD), dimension(666) :: keys
    integer :: keytype, keylen, iostat, status, n
    logical :: report_head

    open(9,file=fitskeyfile,action='read',status='old',iostat=iostat,iomsg=iomsg)
    if( iostat /= 0 ) then
       write(error_unit,*) 'Error: ',trim(iomsg)
       return
    end if

    n = 0
    report_head = .true.
    do
       read(9,'(a)',iostat=iostat,iomsg=iomsg) line
       if( iostat == iostat_end ) exit
       if( iostat /= 0 ) then
          write(error_unit,*) 'Error: ',trim(iomsg),' for the line: `',trim(line),"'"
          exit
       end if

       status = 0
       call fits_parse_template(line,card,keytype,status)
       call fits_get_keyname(card,key,keylen,status)

       if( status == 0 ) then
          n = n + 1
          cards(n) = card
          keys(n) = key
          if( n == size(cards) ) then
             write(error_unit,*) "Warning: A limit of fits keys reached."
             exit
          end if
       else
          if( report_head ) then
             report_head = .false.
             write(error_unit,*) 'Error in `',trim(fitskeyfile),"':"
          end if
          write(error_unit,*) "Failed to parse `",trim(line),"'."
          call fits_report_error(error_unit,status)
       end if

    end do
    close(9)

    if( n > 0 ) then
       allocate(keys_include(n),cards_include(n))
       keys_include = keys(1:n)
       cards_include = cards(1:n)
    end if

  end subroutine include_keywords

end module fitskeys
