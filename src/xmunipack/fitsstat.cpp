/*

  xmunipack - fits implementation of statistics

  Copyright © 2018-2021 F.Hroch (hroch@physics.muni.cz)

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
#include <cfloat>
#include <wx/wx.h>
#include <limits>

using namespace std;

FitsArrayStat::FitsArrayStat(const FitsArray& a, int nskip):
  FitsArray(a),minmax_initialised(false),medmad_initialised(false),
  skip(nskip > 0 ? nskip : max(int(Npixels())/32768,1)),
  med(0.0),mad(0.0),xmin(0.0),xmax(0.0)
{
  wxASSERT(PixelData() && Naxis() > 0);
}

float FitsArrayStat::GetMin()
{
  if( ! minmax_initialised ) {
    Setup_minimax();
    minmax_initialised = true;
  }
  return xmin;
}
float FitsArrayStat::GetMax()
{
  if( ! minmax_initialised ) {
    Setup_minimax();
    minmax_initialised = true;
  }
  return xmax;
}

void FitsArrayStat::Setup_minimax()
{
  // maximum, minimum
  xmin = numeric_limits<double>::max();
  xmax = numeric_limits<double>::min();
  const float *d = PixelData();
  for(long i = 0; i < Npixels(); i++) {
    float x = d[i];
    if( x > xmax ) xmax = x;
    if( x < xmin ) xmin = x;
  }
  wxLogDebug("FitsArrayStat::Setup_minimax() %f %f",xmin,xmax);
}


float FitsArrayStat::GetMed()
{
  if( !medmad_initialised ) {
    Setup_medmad();
    medmad_initialised = true;
  }
  return med;
}

float FitsArrayStat::GetMad()
{
  if( !medmad_initialised ) {
    Setup_medmad();
    medmad_initialised = true;
  }
  return mad;
}

void FitsArrayStat::Setup_medmad()
{
  long nd = Npixels() / skip;
  float *d = new float[nd+1];

  // median
  for(long i = 0, j = 0; j < nd; i += skip, j++) {
    float x = Pixel(i);
    d[j] = x;
  }
  med = QMed(nd,d,nd/2+1);

  // mad = median of absolute deviations
  long n = 0;
  int imax = Npixels() - skip;
  for(int i = 0; n < nd && i < imax; i += skip) {
    float r = Pixel(i) - med;
    if( r > 0.0 )
      d[n++] = r;
  }
  mad = QMed(n,d,n/2+1);

  delete[] d;

  wxLogDebug("Setup_medmad(): %f %f %ld %f %d",mad,med,n,Pixel(1),skip);
}


float FitsArrayStat::QMed(long n, float *a, int k)
{
  float w,x;
  long l,r,i,j;

  x = 0.0;

  l = 0;
  r = n - 1;

  while( l < r ) {
    x = a[k]; i = l; j = r;

    do {

      while( a[i] < x && i < r) i++;
      while( x < a[j] && j > l) j--;

      if( i <= j ) { w = a[i]; a[i] = a[j]; a[j] = w; i++; j--; }

    } while ( i <= j );

    if( j < k ) l = i;
    if( k < i ) r = j;
  }

  return(x);
}
