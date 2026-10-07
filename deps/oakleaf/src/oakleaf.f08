!
!  Oakleaf
!
!  Copyright © 2001-2026 F.Hroch (hroch@physics.muni.cz)
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
!

module oakleaf

  use rfun
  use likescale

!  use madmean
  use quantilmean
  use quicksort
  use medians

  use bfgs_method
  use huber_curvefit

  use qrhouse

  use fluxmean
!  use photonrate
  use robustline

  use robustmean_REAL32
  use robustmean_REAL64
#ifdef HAVE_REAL128
  use robustmean_REAL128
#endif

  interface rmean
     procedure rmean_REAL32
     procedure rmean_REAL64
#ifdef HAVE_REAL128
     procedure rmean_REAL128
#endif
  end interface rmean

end module oakleaf
