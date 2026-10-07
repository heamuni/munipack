/*

  FITS restore utility


  Copyright © 2012-6, 2022 F.Hroch (hroch@physics.muni.cz)

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
#include "fortranio.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <fitsio.h>

#include <cassert>

using namespace std;



bool ParseRecord(const string& record, string& key,string& value, string& comment)
{
  int status = 0;
  const char *card = record.c_str();
  char keyname[FLEN_CARD],cvalue[FLEN_CARD],com[FLEN_CARD];
  int keylen;

  fits_get_keyname((char*)card,keyname,&keylen,&status);
  fits_parse_value((char*)card,cvalue,com,&status);

  if( status == 0 ) {
    key = keyname;
    comment = com;
    value = cvalue;

    vector<string> vs = strftok(cvalue);
    if( vs.size() == 1 )
      value = vs[0];

    return true;
  }
  else
    return false;
}

template<class T>
bool getnum(const string& a, T& t)
{
  istringstream i(a);
  i >> t;
  return(!i.fail());
}

template<class T>
T *parse_line(const string& line, int& ndata, T *data)
{
  vector<T> d;
  istringstream dataline(line);
  while( ! dataline.eof() ) {
    T x;
    dataline >> x;
    d.push_back(x);
  }
  ndata = d.size();
  data = new T[ndata];
  for(size_t i = 0; i < d.size(); i++)
    data[i] = d[i];
  return data;
}


int restore(const string& filename, const string& output)
{
  fitsfile *f;
  int status = 0;
  int nhdu = 0;

  ifstream fin(filename.c_str(),ifstream::in);
  if( ! fin.good() ) {
    cerr << "Failed to open `" << filename << "'" << endl;
    return 1;
  }

  status = 0;
  fits_create_file(&f, output.c_str(), &status);
  if( status == FILE_NOT_CREATED ) {
    mfitsio_unlink(output);
    status = 0;
    fits_create_file(&f, output.c_str(), &status);
  }

  if( status ) goto hell;

  while( fin.good() ) {

    string line;
    getline(fin,line);

    int hdutype = ANY_HDU;

    if( line.substr(0,11) == "# BEGIN HDU" ) {

      // parsing header
      int naxis = 0;
      int bitpix = 0;
      vector<int> naxes;
      int tfields = 0;
      vector<string> tform;
      vector<int> tform_type;
      string extension;

      vector<string> head;
      while( fin.good() ) {

	string record;
	getline(fin,record);
	if( record == "END" ) break;
	head.push_back(record);

	string key,value,comment;
	if( ParseRecord(record,key,value,comment) ) {

	  if( key == "NAXIS" ) {
	    int n;
	    if( getnum(value,n) )
	      naxis = n;
	  }

	  if( key == "BITPIX" ) {
	    int n;
	    if( getnum(value,n) )
	      bitpix = n;
	  }

	  if( key.substr(0,5) == "NAXIS" && key.size() > 5 ) {
	    int n;
	    if( getnum(value,n) )
	      naxes.push_back(n);
	  }

	  if( key == "TFIELDS" ) {
	    int n;
	    if( getnum(value,n) )
	      tfields = n;
	  }

	  if( key.substr(0,5) == "TFORM" ) {
	    tform.push_back(value);

	    if( value.find("E") != string::npos )
	      tform_type.push_back(TFLOAT);
	    else if( value.find("D") != string::npos )
	      tform_type.push_back(TDOUBLE);
	    else if( value.find("I") != string::npos )
	      tform_type.push_back(TSHORT);
	    else if( value.find("B") != string::npos )
	      tform_type.push_back(TBYTE);
	    else if( value.find("A") != string::npos )
	      tform_type.push_back(TSTRING);
	    else {
	      cerr << "Unimplemented TFORM value `" << value << "'." << endl;
	      abort();
	    }
	  }

	  if( key == "XTENSION" )
	    extension = value;
	}
      }

      assert((size_t)naxis == naxes.size());

      if( naxis > 0 && tfields == 0 )
	hdutype = IMAGE_HDU;
      else if( naxis > 0 && tfields > 0 ) {
	assert(!extension.empty());
	if( extension == "CHARTABLE" )
	  hdutype = ASCII_TBL;
	else
	  hdutype = BINARY_TBL;
      }

      nhdu = nhdu + 1;
      fits_create_hdu(f,&status);

      for(size_t i = 0; i < head.size() && status == 0; i++)
	fits_write_record(f,head[i].c_str(),&status);

      if( status ) goto hell;

      // write data
      if( hdutype == IMAGE_HDU ) {

	int fpixel = 1;
	int ndata;

	while( fin.good() ) {
	  string line;
	  getline(fin,line);

	  if( line.substr(0,9) == "# END HDU" ) break;

	  if( bitpix == BYTE_IMG ) {
	    unsigned char *data = 0;
	    data = parse_line(line,ndata,data);
	    fits_write_img(f,TBYTE,fpixel,ndata,data,&status);
	    delete[] data;
	  }
	  else if( bitpix == SHORT_IMG ) {
	    short *data = 0;
	    data = parse_line(line,ndata,data);
	    fits_write_img(f,TSHORT,fpixel,ndata,data,&status);
	    delete[] data;
	  }
	  else if( bitpix == LONG_IMG ) {
	    long *data = 0;
	    data = parse_line(line,ndata,data);
	    fits_write_img(f,TLONG,fpixel,ndata,data,&status);
	    delete[] data;
	  }
	  else if( bitpix == FLOAT_IMG ) {
	    float *data = 0;
	    data = parse_line(line,ndata,data);
	    fits_write_img(f,TFLOAT,fpixel,ndata,data,&status);
	    delete[] data;
	    }
	  else if( bitpix == DOUBLE_IMG ) {
	    double *data = 0;
	    data = parse_line(line,ndata,data);
	    fits_write_img(f,TDOUBLE,fpixel,ndata,data,&status);
	    delete[] data;
	  }

	  fpixel = fpixel + ndata;
	  if( status ) goto hell;
	}
      }
      else if( hdutype == ASCII_TBL || hdutype == BINARY_TBL ) {

	int nrows = naxes.size() >= 2 ? naxes[1] : 0;
	int ncols = tform.size();
	assert(ncols > 0 && nrows > 0);

	void **columns = new void*[ncols];
	for(int k = 0; k < ncols; k++) {
	  if( tform_type[k] == TBYTE ) {
	    char *d = new char[nrows];
	    columns[k] = static_cast<void *>(d);
	  }
	  else if( tform_type[k] == TSHORT ) {
	    short *d = new short[nrows];
	    columns[k] = static_cast<void *>(d);
	  }
	  else if( tform_type[k] == TSTRING ) {
	    char **d = new char*[nrows];
	    columns[k] = static_cast<void *>(d);
	  }
	  else if( tform_type[k] == TFLOAT ) {
	    float *d = new float[nrows];
	    columns[k] = static_cast<void *>(d);
	  }
	  else if( tform_type[k] == TDOUBLE ) {
	    double *d = new double[nrows];
	    columns[k] = static_cast<void *>(d);
	  }
	}

	int nc = 0;
	while( fin.good() && nc != nrows ) {
	  string line;
	  getline(fin,line);
	  if( line.substr(0,9) == "# END HDU" ) break;

	  vector<string> items(strftok(line));
	  if( items.size() != tform.size() ) {
	    cerr << "Failed to parse table line: `" << line << "'" << endl;
	    status = 666;
	    goto hell;
	  }

	  for(size_t k = 0; k < items.size(); k++) {

	    if( tform_type[k] == TBYTE ) {
	      string s(items[k]);
	      char d = s.size() == 1 ? s[0] : ' ';
	      char *data = static_cast<char *>(columns[k]);
	      data[nc] = d;
	    }
	    else if( tform_type[k] == TSHORT ) {
	      short d;
	      short *data = static_cast<short *>(columns[k]);
	      if( getnum(items[k],d) )
		data[nc] = d;
	      else
		data[nc] = 0;
	    }
	    else if( tform_type[k] == TSTRING ) {
	      string d(items[k]);
	      char **data = static_cast<char **>(columns[k]);
	      data[nc] = strdup(d.c_str());
	    }
	    else if( tform_type[k] == TFLOAT ) {
	      float d;
	      float *data = static_cast<float *>(columns[k]);
	      if( getnum(items[k],d) )
		data[nc] = d;
	      else
		data[nc] = 0.0;
	    }
	    else if( tform_type[k] == TDOUBLE ) {
	      double d;
	      double *data = static_cast<double *>(columns[k]);
	      if( getnum(items[k],d) )
		data[nc] = d;
	      else
		data[nc] = 0;
	    }
	  }
	  nc++;
	}

	if( nc != nrows ) {
	  cerr << "Count of table rows exceeds definition." << endl;
	  status = 666;
	  goto hell;
	}

	for(int k = 0; k < ncols && status == 0; k++) {

	  int c = k + 1;
	  long firstrow = 1, firstelem = 1;

	  if( tform_type[k] == TBYTE ) {
	    char *d = static_cast<char *>(columns[k]);
	    fits_write_col(f,TBYTE,c,firstrow,firstelem,nrows,d,&status);
	    delete[] d;
	  }
	  else if( tform_type[k] == TSHORT ) {
	    short *d = static_cast<short *>(columns[k]);
	    fits_write_col(f,TSHORT,c,firstrow,firstelem,nrows,d,&status);
	    delete[] d;
	  }
	  else if( tform_type[k] == TSTRING ) {
	    char **d = static_cast<char **>(columns[k]);
	    fits_write_col(f,TSTRING,c,firstrow,firstelem,nrows,d,&status);
	    for(int l = 0; l < nrows; l++)
	      free(d[l]);
	    delete[] d;
	  }
	  else if( tform_type[k] == TFLOAT ) {
	    float *d = static_cast<float *>(columns[k]);
	    fits_write_col(f,TFLOAT,c,firstrow,firstelem,nrows,d,&status);
	    delete[] d;
	  }
	  else if( tform_type[k] == TDOUBLE ) {
	    double *d = static_cast<double *>(columns[k]);
	    fits_write_col(f,TDOUBLE,c,firstrow,firstelem,nrows,d,&status);
	    delete[] d;
	  }
	}
	delete[] columns;

	// end parsing data
      }

      // end hdu
    }
  }

 hell:

  fits_close_file(f,&status);
  fits_report_error(stderr, status);

  return status;
}
