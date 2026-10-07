!
!  Quantile mean
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

module quantilmean

  use quantilmean_REAL32
  use quantilmean_REAL64
#ifdef HAVE_REAL128
  use quantilmean_REAL128
#endif

  interface qmean
     procedure qmean_REAL32
     procedure qmean_REAL64
#ifdef HAVE_REAL128
     procedure qmean_REAL128
#endif
  end interface qmean

  interface byQ
     procedure byQ_REAL32
     procedure byQ_REAL64
#ifdef HAVE_REAL128
     procedure byQ_REAL128
#endif
  end interface byQ

  interface ecdf
     procedure ecdf_REAL32
     procedure ecdf_REAL64
#ifdef HAVE_REAL128
     procedure ecdf_REAL128
#endif
  end interface ecdf

  interface quantile
     procedure quantile_REAL32
     procedure quantile_REAL64
#ifdef HAVE_REAL128
     procedure quantile_REAL128
#endif
  end interface quantile

  ! obsolete
  interface quantilefun
     procedure quantilefun_REAL32
     procedure quantilefun_REAL64
#ifdef HAVE_REAL128
     procedure quantilefun_REAL128
#endif
  end interface quantilefun

  ! obsolete
  interface qbracket
     procedure qbracket_REAL32
     procedure qbracket_REAL64
#ifdef HAVE_REAL128
     procedure qbracket_REAL128
#endif
  end interface qbracket

end module quantilmean
