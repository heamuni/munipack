/*

  FITS related utility


  Copyright © 2011 F.Hroch (hroch@physics.muni.cz)

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
#include <fitsio.h>

using namespace std;


int structure(const string& filename)
{
  fitsfile *f;
  int status = 0;
  int nhdu = 0;
  long *naxes = 0;
  int hdutype,bitpix,naxis,ncols;
  long nrows;
  char extname[FLEN_CARD],comment[FLEN_CARD];

  cout << "#";
  cout.width(11);
  cout << "EXTNAME";
  cout.width(11);
  cout << "TYPE";
  cout.width(7);
  cout << "BITPIX";
  cout << " SIZE";
  cout << endl;

  status = 0;
  if( fits_open_file(&f, filename.c_str(), READONLY, &status) )
    fits_report_error(stderr, status);

  if( fits_get_num_hdus(f,&nhdu,&status) )
    fits_report_error(stderr, status);

  for(int k = 0; k < nhdu; k++) {

    if( fits_movabs_hdu(f,k+1,&hdutype,&status) )
      fits_report_error(stderr, status);

    fits_read_keyword(f,"EXTNAME",extname,comment,&status);
    if( status != 0 ) {
      status = 0;
      fits_read_keyword(f,"HDUNAME",extname,comment,&status);
      status = 0;
    }

    if( hdutype == IMAGE_HDU ) {

      fits_get_img_type(f,&bitpix,&status);
      fits_get_img_dim(f,&naxis,&status);

      if( naxis > 0 ) {
	naxes = new long[naxis];
	fits_get_img_size(f,naxis,naxes,&status);
      }
    }
    else if( hdutype == ASCII_TBL || hdutype == BINARY_TBL ) {

      fits_get_num_rows(f,&nrows,&status);
      fits_get_num_cols(f,&ncols,&status);

    }

    /*
    fits_read_key(f,TINT,"BITPIX",&bitpix,comment,&status);
    fits_read_key(f,TINT,"NAXIS",&naxis,comment,&status);

    if( status == 0 ) {

    }
    */

    string htype /*= "**********"*/;
    switch (hdutype) {
    case IMAGE_HDU:  htype = "IMAGE"; break;
    case ASCII_TBL:  htype = "ASCII_TBL"; break;
    case BINARY_TBL: htype = "BINARY_TBL"; break;
    }

    if( naxis == 0 && hdutype == IMAGE_HDU )
      htype = "DUMMY";

    cout << k + 1;
    cout.width(11);
    cout << extname;
    cout.width(11);
    cout << htype;
    cout.width(7);

    if( hdutype == IMAGE_HDU )
      cout << bitpix;
    else
      cout << "";

    cout << " ";

    if( hdutype == IMAGE_HDU && naxes) {

      for(size_t i = 0; i < (size_t) naxis; i++) {
	//	if( i != 0 ) cout << "x";
	cout << (i != 0 ? "x" : "") << naxes[i];
      }

    }
    else if( hdutype == ASCII_TBL || hdutype == BINARY_TBL ) {

      cout << nrows << "x" << ncols;

    }

    cout << endl;

    delete[] naxes;
    naxes = 0;
  }

  fits_close_file(f, &status);

  if( status != 0 )
    fits_report_error(stderr, status);

  return status;
}
