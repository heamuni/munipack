/*

  FITS utility for C/C++

  Copyright © 2012, 2018 F.Hroch (hroch@physics.muni.cz)

  This file is part of Munipack.

  Munipack is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Munipack is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Munipack.  If not, see <http://www.gnu.org/licenses/>.


*/

#include <fitsio.h>
#include <string>

int mfitsio_copy(const std::string&, const std::string&);
int mfitsio_unlink(const std::string&);
/*
std::string mfitsio_scratch_unique();
int fitsio_scratch_init(fitsfile **, std::string *, int *);
int fitsio_scratch_open(fitsfile **, const std::string&, std::string *, int *);
int fitsio_scratch_keep(const std::string&, const std::string&,
			const std::string&, int *);
*/
