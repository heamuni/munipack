/*

  FITS related utility


  Copyright © 2011-2, 2018 F.Hroch (hroch@physics.muni.cz)

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
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <fitsio.h>

using namespace std;


int table(const string& filename)
{
  fitsfile *f;
  int status = 0;

  status = 0;
  if( fits_open_table(&f, filename.c_str(), READONLY, &status) )
    fits_report_error(stderr, status);

  long nrows, ncols;
  int nc;
  fits_get_num_rows(f,&nrows,&status);
  fits_get_num_cols(f,&nc,&status);
  ncols = nc;

  if (status != 0 ) {
    fits_report_error(stderr, status);
    return status;
  }

  int widths[ncols];
  for(int k = 0; k < ncols && status == 0; k++)
    fits_get_col_display_width(f,k+1,&widths[k],&status);

  if (status != 0 ) {
    fits_report_error(stderr, status);
    return status;
  }

  char **table = new char*[ncols*nrows];
  for(int k = 0; k < ncols; k++)
    for(int i = 0; i < nrows; i++)
      table[i+k*nrows] = new char[widths[k]];


  long frow = 1, felem = 1; char *nullval = 0; int dummy;
  for(int k = 0; k < ncols && status == 0; k++ )
    fits_read_col(f, TSTRING, k+1, frow, felem, nrows, &nullval,
		  table + k*nrows,&dummy,&status);

  fits_close_file(f, &status);

  if( status == 0 )

    for(int i = 0; i < nrows; i++) {

      for(int k = 0; k < ncols; k++)
	cout << table[i+nrows*k] << " ";

      cout << endl;
    }

  else

    fits_report_error(stderr, status);

  for(int i = 0; i < ncols*nrows; i++)
    delete[] table[i];
  delete[] table;

  return status;
}
