/*

  FITS cat utility

  Copyright © 2018 F.Hroch (hroch@physics.muni.cz)

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

#include "fits.h"
#include "mfitsio.h"
#include <string>
#include <fitsio.h>

#include <cassert>

using namespace std;


int cat(const string& filename, const string& output)
{
  fitsfile *f, *o;
  int status;

  status = 0;
  if( fits_open_file(&f, filename.c_str(), READONLY, &status) == 0 ) {

    fits_create_file(&o, output.c_str(), &status);
    if( status == FILE_NOT_CREATED ) {
      mfitsio_unlink(output);
      status = 0;
      fits_create_file(&o, output.c_str(), &status);
    }
    if( status == 0)
      fits_copy_file(f,o, 1, 1, 1, &status);

    fits_close_file(o, &status);
  }
  fits_close_file(f, &status);

  fits_report_error(stderr, status);

  return status;
}
