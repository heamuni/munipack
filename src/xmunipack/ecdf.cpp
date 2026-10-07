/*

  xmunipack - empirical cummulative distribution function

  Copyright © 2018-2025 F.Hroch (hroch@physics.muni.cz)

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
#include <cstdio>
#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;

// cumulative distribution function of data
EmpiricalCDF::EmpiricalCDF(long n, const float *data): ncdf(0),xcdf(0),ycdf(0)
{
  const double macheps = std::numeric_limits<float>::epsilon();
  wxASSERT(n > 0 && data);
  //wxLogDebug("EmpiricalCDF::EmpiricalCDF(long n, const float *data):");

  float *d = new float[n];
  copy(data,data+n,d);
  sort(d,d+n);

  // remove duplicite points + Nans + Infs
  long m = 0;
  long k = 0;
  while( !isfinite(d[k]) && k < n ) k++;
  d[0] = d[k];
  for(long i = k; i < n; i++) {
    if( fabsf(d[i] - d[m]) > 2*macheps && isfinite(d[i]) )
      d[m++] = d[i];
  }

  ncdf = m;
  xcdf = new float[m];
  ycdf = new float[m];
  copy(d,d+m,xcdf);
  float h = 1.0 / float(m + 1);
  for(long i = 0; i < m; i++)
    ycdf[i] = i*h;

  /*
  FILE *f = fopen("/tmp/cdf","w");
  for(int i = 0; i < m; i++)
    fprintf(f,"%f %f\n",xcdf[i],ycdf[i]);
  fclose(f);
  */
  //  wxLogFatalError("ss");


  delete[] d;
}

// copy constructor
EmpiricalCDF::EmpiricalCDF(const EmpiricalCDF& cdf):
  ncdf(cdf.ncdf), xcdf(new float[ncdf]), ycdf(new float[ncdf])
{
  //  wxLogDebug("EmpiricalCDF::EmpiricalCDF(const EmpiricalCDF& cdf)");
  copy(cdf.xcdf,cdf.xcdf+ncdf,xcdf);
  copy(cdf.ycdf,cdf.ycdf+ncdf,ycdf);
}

// assign constructor
EmpiricalCDF& EmpiricalCDF::operator=(const EmpiricalCDF& other)
{
  //  wxLogDebug("EmpiricalCDF& EmpiricalCDF::operator=(const EmpiricalCDF&)");
  if( this != &other ) {
    long n = other.ncdf;
    float *x = new float[n];
    float *y = new float[n];

    copy(other.xcdf,other.xcdf+n,x);
    copy(other.ycdf,other.ycdf+n,y);

    if( ncdf > 0 ) {
      delete[] xcdf;
      delete[] ycdf;
    }
    ncdf = n;
    xcdf = x;
    ycdf = y;
  }
  return *this;
}


EmpiricalCDF::~EmpiricalCDF()
{
  //  wxLogDebug("EmpiricalCDF::~EmpiricalCDF() %ld",ncdf);
  ncdf = 0;
  delete[] xcdf;
  delete[] ycdf;
}

bool EmpiricalCDF::IsOk() const
{
  return ycdf && xcdf;
}

float EmpiricalCDF::GetQuantile(float q) const
{
  const double macheps = std::numeric_limits<float>::epsilon();
  wxASSERT(xcdf && ycdf);

  long n = ncdf;

  if( n == 0 )
    return 0.0;
  else if( n == 1 )
    return xcdf[0];

  if( q < ycdf[0] )
    return xcdf[0];
  else if( q > ycdf[n-1] )
    return xcdf[n-1];
  else {
    float h = 1.0 / float(n + 1);
    float r = q / h;
    int m = round(r);
    if( fabsf(m - r) < 10*macheps )
      return xcdf[m];
    else {
      int low = int(r);
      int high = low + 1;
      //      wxLogDebug("%d %d %ld %f %f",low,high,ncdf,q,h);

      float dy = ycdf[high] - ycdf[low];
      if( fabsf(dy) > 10*macheps )
	// inverse by linear interpolation
	return (xcdf[high] - xcdf[low])/dy*(q - ycdf[low]) + xcdf[low];
      else {
	// nearly singular
	return (xcdf[high] + xcdf[low]) / 2;
      }
    }
  }
}

void EmpiricalCDF::GetPositiveQuantile(float& qblack, float& black) const
{
  wxASSERT(xcdf && ycdf);

  for(long i = 0; i < ncdf; i++)
    if( xcdf[i] > 0 ) {
      black = xcdf[i];
      qblack = ycdf[i];
      return;
    }
}

float EmpiricalCDF::GetInverse(float x) const
{
  wxASSERT(xcdf && ycdf);

  long n = ncdf;

  if( n == 0 )
    return 0.0;
  else if( n == 1 )
    return ycdf[0];

  if( x < xcdf[0] )
    return ycdf[0];
  else if( x >= xcdf[n-1] )
    return ycdf[n-1];
  else {
    for(long i = 1; i < ncdf; i++)
      if( xcdf[i-1] <= x && x < xcdf[i] ) {
	float a = (ycdf[i] - ycdf[i-1]) / (xcdf[i] - xcdf[i-1]);
	return a*(x - xcdf[i-1]) + ycdf[i-1];
      }
  }

  return 0.5;
}
