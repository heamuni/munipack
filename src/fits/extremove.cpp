/*

  FITS related utility


  Copyright © 2016, 2018 F.Hroch (hroch@physics.muni.cz)

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
#include <vector>
#include <fitsio.h>

#include <cassert>

using namespace std;


int ext_remove(const string& filename, const vector<string>& extnames,
	       const string& output)
{
  fitsfile *f = 0;
  int status, hdutype;
  const int extver = 0;
  string name;

  if( filename == output )
    name = filename;
  else {
    mfitsio_unlink(output);
    mfitsio_copy(filename,output);
    name = output;
  }

  status = 0;
  if( fits_open_file(&f, name.c_str(), READWRITE, &status) == 0 ) {

    vector<string>::const_iterator extname;
    for( extname = extnames.begin(); extname < extnames.end(); extname++ ) {

      fits_movabs_hdu(f,1,&hdutype,&status);
      while( status == 0 ) {

	fits_movnam_hdu(f,ANY_HDU,const_cast<char*>(extname->c_str()),
			extver,&status);
	if( status == BAD_HDU_NUM ) {
	  break;
	}
	fits_delete_hdu(f,&hdutype,&status);
      }
    }
  }

  fits_close_file(f, &status);
  fits_report_error(stderr, status);

  return status;
}
