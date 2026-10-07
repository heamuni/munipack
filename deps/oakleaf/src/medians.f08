!
!  Medians - the interface
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

module medians

  use medians_REAL32
  use medians_REAL64
#ifdef HAVE_REAL128
  use medians_REAL128
#endif

  interface median
     procedure median_REAL32
     procedure median_REAL64
#ifdef HAVE_REAL128
     procedure median_REAL128
#endif
  end interface median

  interface qmedian
     procedure qmedian_REAL32
     procedure qmedian_REAL64
#ifdef HAVE_REAL128
     procedure qmedian_REAL128
#endif
  end interface qmedian

end module medians
