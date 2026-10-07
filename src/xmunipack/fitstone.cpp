/*

  xmunipack - pre-scaling

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
#include <algorithm>
#include <limits>

#define QBLACK 0.25
#define QSENSE 0.95

using namespace std;

FitsTone::FitsTone(): initialised(false), black(0.0),qblack(QBLACK),
		      sense(1.0),rsense(1.0),refsense(1.0),refblack(0.0),
		      refqblack(QBLACK) {}

FitsTone::FitsTone(const FitsArray& array):
  initialised(false), black(0.0),qblack(QBLACK),sense(1.0),rsense(1.0),refsense(1.0),
  refblack(0.0),refqblack(QBLACK)
{
  wxASSERT(array.IsOk());

  bool iscolour = array.GetKey("CSPACE").Matches("*XYZ*") && array.GetDepth() == 3;

  wxLogDebug("FitsTone::FitsTone(const FitsArray&): %d %d",
	     int(array.GetDepth()),iscolour);

  int k = iscolour ? 1 : 0;
  long npix = array.GetWidth()*array.GetHeight();
  const float *data = array.PixelData() + k*npix;
  Setup(npix,data,iscolour);
}

FitsTone::FitsTone(const FitsTone& tone):
  initialised(tone.initialised), black(tone.black), qblack(tone.qblack),
  sense(tone.sense),rsense(tone.rsense),refsense(tone.refsense),
  refblack(tone.refblack),refqblack(tone.refqblack), cdf_back(tone.cdf_back) {}

FitsTone& FitsTone::operator = (const FitsTone& other)
{
  if( this != & other ) {
    initialised = other.initialised;
    black = other.black;
    qblack = other.qblack;
    sense = other.sense;
    rsense = other.rsense;
    refsense = other.refsense;
    refblack = other.refblack;
    refqblack = other.refqblack;
    cdf_back = other.cdf_back;
  }
  return *this;
}

FitsTone::FitsTone(long npix, const float *data):
  initialised(false), black(0.0),qblack(QBLACK),sense(1.0),rsense(1.0),refsense(1.0),
  refblack(0.0), refqblack(QBLACK)
{
  Setup(npix,data,false);
}

bool FitsTone::IsOk() const
{
  return initialised && cdf_back.IsOk();
}

void FitsTone::Setup(long npix, const float *data, bool iscolour)
{
  const double macheps = std::numeric_limits<float>::epsilon();
  wxASSERT(npix > 0 && data);

  const long nmax = 32768;
  const long skip = max(npix / nmax, long(1));
  const long nd = npix / skip;
  wxLogDebug("FitsTone::Setup %ld %ld %f",nd,skip,qblack);

  float *d = new float[nd+1];

  // Part I.
  // CDF over image grid,
  // if sparse field is analysed, it's CDF of sky and: Q(0.5) = sigma^2
  long m = 0;
  for(long i = 0; i < npix - skip && m < nd; i += skip)
    d[m++] = data[i];

  cdf_back = EmpiricalCDF(m,d);

  qblack = QBLACK;
  black = cdf_back.GetQuantile(qblack);
  if( iscolour )
    cdf_back.GetPositiveQuantile(qblack,black);

  float med = cdf_back.GetQuantile(0.5);
  float mad = (cdf_back.GetQuantile(0.75) - cdf_back.GetQuantile(0.25)) / 2;
  // MAD, computed this way, can be treacherous.

  // Part II
  // To estimate starlight, only pixels above the threshold are accepted,
  // the pixels are collected over a grid, eventually decreasing threshold,
  // to reach sufficient amount of data; the limit 16 for the side corresponds
  // with default value of nmax 32768.
  long n = 0;
  for(int kappa = 30; kappa >= 0; kappa=kappa-5 ) {
    n = 0;
    float thresh = (0.1*kappa) * mad / 0.6745;
    long side = 1;
    while( n < (nmax / 10) && side <= 16) {
      side = 2*side;
      long skip2 = max(skip / side,long(1));
      long imax = npix - skip2;
      for(long i = 0; n < nd && i < imax; i += skip2) {
	float r = data[i] - med;
	if( r > thresh )
	  d[n++] = r;
      }
    }
    wxLogDebug("kappa=%f thresh=%f n=%ld",0.1*kappa,thresh,n);
    if( n > 42 ) break;
  }

  if( n > 0 ) {
    EmpiricalCDF cdf_light(n,d);
    refsense = cdf_light.GetQuantile(QSENSE);
  }
  else if( cdf_back.IsOk() && med > 100*macheps ) {
    refsense = 3*(cdf_back.GetQuantile(QSENSE) - cdf_back.GetQuantile(1-QSENSE));
  }
  else {
    // uniform intensity
    refsense = black > 0 ? sqrt(black) : 1;
  }
  delete[] d;

  // last resort
  if( fabs(refsense) < 10*macheps ) refsense = 1;

    //    wxLogDebug("Estim: %f %f %f %f %ld %ld",black,med,sig,refsense,m,n);
  wxLogDebug("Median=%f MAD=%f Black=%f Qlight=%f n=%ld",med,mad,black,refsense,n);

  refblack = black;
  refqblack = refqblack;
  sense = refsense;
  rsense = 1;

  initialised = true;
}



void FitsTone::SetBlack(float b)
{
  black = b;
  qblack = cdf_back.GetInverse(black);
}

void FitsTone::SetSense(float s)
{
  wxASSERT(s > 0);
  sense = s;
  rsense = refsense / sense;
}

void FitsTone::SetQblack(float q)
{
  qblack = q;
  black = cdf_back.GetQuantile(qblack);
}

void FitsTone::SetRsense(float r)
{
  wxASSERT(r > 0);
  rsense = r;
  sense = refsense / rsense;
}


float *FitsTone::Scale(long n, const float *a) const
{
  wxASSERT(initialised);

  float *f = new float[n];
  for(long i = 0; i < n; i++) {
    float t = (a[i] - black) / sense;
    f[i] = max(t, 0.0f);
  }
  return f;
}

void FitsTone::Reset()
{
  qblack = refqblack;
  black = refblack;
  sense = refsense;
  rsense = 1.0;
}
