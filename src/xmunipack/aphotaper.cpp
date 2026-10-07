/*

  Apertures for the display of the perture photometry tool

  Copyright © 2019 F.Hroch (hroch@physics.muni.cz)

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

#include "aphot.h"
#include <wx/wx.h>
#include <cmath>
#include <cfloat>

using namespace std;


MuniApertures::MuniApertures(int n, int rmax, bool s)
{
  wxASSERT(n > 0);

  naper = n;
  aper = new double[naper];
  spiral = s;

  if( spiral ) {
    const double pi = 3.1415926536;
    double pitch = log(rmax) / pi;
    for(int i = 0; i < naper; i++) {
      double theta = i / (naper - 1.0) * pi;
      aper[i] = exp(pitch*theta);
    }
  }
  else {
    double dr = (rmax - 1.0) / (naper - 1.0);
    for(int i = 0; i < naper; i++) {
      aper[i] = 1 + i*dr;
    }
  }
}

// copy constructor
MuniApertures::MuniApertures(const MuniApertures& orig):
  naper(orig.naper), aper(new double[naper]), spiral(orig.spiral)
{
  copy(orig.aper,orig.aper+naper,aper);
}

// assign constructor
MuniApertures& MuniApertures::operator=(const MuniApertures& other)
{
  if( this != &other ) {
    int n = other.naper;
    double *tmp = new double[n];
    copy(other.aper,other.aper+n,tmp);
    if( naper > 0 )
      delete[] aper;
    naper = n;
    aper = tmp;
    spiral = other.spiral;
  }
  return *this;
}


MuniApertures::~MuniApertures()
{
  delete[] aper;
}

double MuniApertures::GetAper(int i) const
{
  wxASSERT(1 <= i && i <= naper);
  return aper[i-1];
}

bool MuniApertures::FindAper(double r, int tol, int *n) const
{
  bool found = false;
  double rmin = 2*tol;

  for(int i = 0; i < naper; i++) {
    double d = abs(r - aper[i]);
    if( d < rmin ) {
      found = true;
      rmin = d;
      *n = i + 1;
    }
  }

  return found;
}
