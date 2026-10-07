/*

  FITS dump utility


  Copyright © 2012, 2022 F.Hroch (hroch@physics.muni.cz)

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
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <fitsio.h>

#include <cassert>

using namespace std;

int dump(const string& filename, const string& outname)
{
  fitsfile *f;
  int status = 0;
  int nhdu = 0;

  ofstream fout(outname.c_str());
  if( ! fout.good() ) {
    cerr << "Failed to write to `" << filename << "'" << endl;
    return 1;
  }

  status = 0;
  fits_open_file(&f, filename.c_str(), READONLY, &status);
  fits_get_num_hdus(f,&nhdu,&status);

  if( status ) goto hell;

  for(int k = 0; k < nhdu; k++) {

    fout << "# BEGIN HDU " << k << endl;

    int htype, dummy;

    fits_movabs_hdu(f,k+1,&htype,&status);
    if( status ) goto hell;

    int nhead;
    char record[FLEN_CARD];
    fits_get_hdrspace(f,&nhead,&dummy,&status);
    for(int n = 0; status == 0 && n < nhead; n++) {
      if( fits_read_record(f,n+1,record,&status) == 0 )
	fout << record << endl;
    }
    if( status ) goto hell;

    fout << "END" << endl;

    if( htype == IMAGE_HDU ) {

      int bitpix, naxis;

      fits_get_img_type(f,&bitpix,&status);
      fits_get_img_dim(f,&naxis,&status);

      if( naxis > 0 ) {

	long *naxes = new long[naxis];
	fits_get_img_size(f,naxis,naxes,&status);
	if( status ) { delete[] naxes; goto hell; }

	long ndata = 1; for(int i = 1; i < naxis; i++ ) ndata = ndata*naxes[i];
	assert(ndata > 0);

	if( bitpix > 0 ) {

	  int nullval = 0;
	  int *data = new int[ndata];
	  long fpixel = 1;
	  for(int j = 0; j < naxes[0] && status == 0; j++) {
	    fits_read_img(f,TINT,fpixel,ndata,&nullval,data,&dummy,&status);
	    if( status == 0 ) {
	      fpixel = fpixel + ndata;
	      for(int i = 0; i < ndata; i++)
		fout << " " << data[i];
	      fout << endl;
	    }
	  }
	  delete[] data;
	  if( status ) goto hell;
	}
	else {

	  double nullval = 0;
	  double *data = new double[ndata];
	  long fpixel = 1;
	  for(int j = 0; j < naxes[0] && status == 0; j++) {
	    fits_read_img(f,TDOUBLE,fpixel,ndata,&nullval,data,&dummy,&status);
	    if( status == 0 ) {
	      fpixel = fpixel + ndata;
	      for(int i = 0; i < ndata; i++)
		fout << " " << data[i];
	      fout << endl;
	    }
	  }
	  delete[] data;
	  if( status ) goto hell;

	}


	delete[] naxes;
      }

    }
    else if( htype == ASCII_TBL || htype == BINARY_TBL ) {

      long nrows;
      int ncols;
      fits_get_num_rows(f,&nrows,&status);
      fits_get_num_cols(f,&ncols,&status);

      if( status ) goto hell;

      int *widths = new int[ncols];
      int *types = new int[ncols];
      long repeat, widths2;
      for(int k = 0; k < ncols && status == 0; k++) {
	fits_get_coltype(f,k+1,&types[k],&repeat,&widths2,&status);
	fits_get_col_display_width(f,k+1,&widths[k],&status);
      }


      if( status ) { delete[] widths; delete[] types; goto hell;}

      char **table = new char*[ncols*nrows];
      for(int k = 0; k < ncols; k++)
	for(int i = 0; i < nrows; i++)
	  table[i+k*nrows] = new char[widths[k]];


      long frow = 1, felem = 1; char *nullval = 0; int dummy;
      for(int k = 0; k < ncols && status == 0; k++ )
	fits_read_col(f, TSTRING, k+1, frow, felem, nrows, &nullval,
		      table + k*nrows,&dummy,&status);

      if( status == 0 ) {

	for(int i = 0; i < nrows; i++) {
	  for(int k = 0; k < ncols; k++) {
	    int b = widths[k];
	    char *cell = table[i+nrows*k];
	    char c[2*b+3];
	    if( types[k] == TSTRING ) {
	      c[0] = '\'';
	      int l = 0;
	      for(l = 0; l != b && cell[l] != '\0'; l++) {
		if( cell[l] == '\'' ) {
		  c[l++] = '\'';
		  c[l+1] = '\'';
		}
		else
		  c[l+1] = cell[l];
	      }
	      c[l+1] = '\'';
	      c[l+2] = '\0';
	    }
	    else {
	      strncpy(c,cell,b);
	      c[b] = '\0';
	    }
	    fout << c << " ";

	  }
	  fout << endl;
	}
      }

      for(int i = 0; i < ncols*nrows; i++)
	delete[] table[i];
      delete[] table;
      delete[] widths;
      delete[] types;

    }

    fout << "# END HDU " << k << endl;
  }

 hell:

  fits_report_error(stderr, status);
  fits_close_file(f, &status);
  return status;
}
