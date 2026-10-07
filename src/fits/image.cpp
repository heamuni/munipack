/*

  FITS related utility


  Copyright © 2011-2 F.Hroch (hroch@physics.muni.cz)

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
#include <string.h>
#include <fitsio.h>

using namespace std;


int image(const string& filename)
{
  fitsfile *f;
  int status = 0;
  int fpixel = 1;
  int dummy,naxis;
  double nullval = 0.0;

  status = 0;
  if( fits_open_image(&f, filename.c_str(), READONLY, &status) )
    fits_report_error(stderr, status);

  fits_get_img_dim(f,&naxis,&status);

  if( naxis > 0 && status == 0 ) {

    long naxes[naxis];
    fits_get_img_size(f,naxis,naxes,&status);
    long ndata = 1; for(int i = 0; i < naxis; i++ ) ndata = ndata*naxes[i];
    double *image = new double[ndata];
    fits_read_img(f,TDOUBLE,fpixel,ndata,&nullval,image,&dummy,&status);

    for(int idx = 0; idx < ndata; idx++) {

      long q[naxis];
      long s =  idx;
      for(int l = naxis-1; 0 < l; l--) {
	q[l] = s/naxes[l];
	s = s - q[l]*naxes[l];
      }
      q[0] = s;

      for(int i = 0; i < naxis; i++)
	cout << q[i]+1 << " ";
      cout << image[idx] << endl;

    }

    /*
    for(int n = 0; n < naxis; n++) q[n] = 0;
    long idx = 0;
    for(int n = 0; n < naxis; n++) {
      for(int l = 0; l < naxes[n]; l++) {
	cout << naxes[n] << endl;

	q[n] = l;
	long s = q[0];
	for(int i = 1; i < naxis; i++)
	  s = q[i]*naxes[i];
	idx = s;

	for(int i = 0; i < naxis; i++)
	  cout << q[i]+1 << " ";
	cout << idx << " " << image[idx] << endl;
      }
    }
    */
    delete[] image;

  }

  fits_close_file(f, &status);

  if( status != 0 )
    fits_report_error(stderr, status);


  return status;
}
