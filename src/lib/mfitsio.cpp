/*

  C++ versions of soubroutines in fitsio.f95

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

#include "mfitsio.h"
#include <string>
#include <fitsio.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/*
struct {
  fitsfile *fptr;
  string filename, scratch, output;
} FITSFILE;

static FITSFILE fitsfiles[50];
*/

int mfitsio_copy(const string& ifile, const string& ofile)
{
  fitsfile *fin,*fout;
  int status;

  status = 0;
  fits_open_file(&fin,ifile.c_str(),READONLY,&status);
  fits_create_file(&fout,ofile.c_str(),&status);
  fits_copy_file(fin,fout,1,1,1,&status);
  fits_close_file(fin,&status);
  fits_close_file(fout,&status);
  fits_report_error(stderr,status);
  return status;
}

int mfitsio_unlink(const string& filename)
{
  fitsfile *f;
  int status = 0;

  fits_open_file(&f,filename.c_str(),READWRITE,&status);
  fits_delete_file(f,&status);
  fits_report_error(stderr,status);
  return status;
}
