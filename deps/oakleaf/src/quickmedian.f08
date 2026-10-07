!
!  quick median - interfaces to quickmedian
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
!

module quickmedian

  use quickmedian_REAL32
  use quickmedian_REAL64
#ifdef HAVE_REAL128
  use quickmedian_REAL128
#endif

  interface qmed
     procedure qmed_REAL32
     procedure qmed_REAL64
#ifdef HAVE_REAL128
     procedure qmed_REAL128
#endif
  end interface qmed

end module quickmedian
