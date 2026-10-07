!
!  QUICKSORT - interfaces for QuickSort algorithm
!
!  Copyright © 2025 F.Hroch (hroch@physics.muni.cz)
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


module quicksort

  use quicksort_REAL32
  use quicksort_REAL64
#ifdef HAVE_REAL128
  use quicksort_REAL128
#endif

  interface qsort
     procedure qsort_REAL32
     procedure qsort_REAL64
#ifdef HAVE_REAL128
     procedure qsort_REAL128
#endif
  end interface qsort

end module quicksort
