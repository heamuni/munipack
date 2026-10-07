!
!  profile photometry
!
!  Copyright © 2013, 2016 F.Hroch (hroch@physics.muni.cz)
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


program pphot

  use titsio
  use fitspphot

  implicit none

  character(len=4*FLEN_FILENAME) :: record,key,val
  character(len=FLEN_FILENAME) :: outname, file,backup, output = ''
  character(len=FLEN_KEYWORD), dimension(1) :: fkeys
  logical :: verbose = .false., plog = .false.
  integer :: eq,status

  fkeys(1) = FITS_KEY_GAIN

  do
     read(*,'(a)',end=20) record

     eq = index(record,'=')
     if( eq == 0 ) stop 'Malformed input record.'
     key = record(:eq-1)
     val = record(eq+1:)

     if( key == 'VERBOSE' ) then

        read(val,*) verbose

     else if( key == 'PIPELOG' ) then

        read(val,*) plog

     else if( key == 'FITS_KEY_GAIN' ) then

           read(val,*) fkeys(1)

     else if( key == 'FILE' ) then

        read(val,*) file, backup, output

        status = 0
!        call fitsback(file,backup,output,.false.,outname,status)

        if( verbose ) then
           write(*,*)
           write(*,'(a)') "========  Processing file: "//trim(outname)
           write(*,*)
        end if

     end if

  end do

20 continue

  stop 0

end program pphot
