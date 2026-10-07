/*

  xmunipack - fits image histogram

  Copyright © 2018, 2021-25 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/wx.h>
#include <algorithm>
#include <limits>

#if wxDEBUG_LEVEL > 1
#include <wx/debug.h>
#include <wx/txtstrm.h>
#include <wx/wfstream.h>
#endif


// ----------  FitsHisto

class FitsHistoData : public wxObjectRefData
{
public:
  FitsHistoData();
  FitsHistoData(int, int *, double *);
  FitsHistoData(const FitsHistoData&);
  FitsHistoData& operator = (const FitsHistoData&);
  virtual ~FitsHistoData();

  int nbin;
  int *count;  // counts of data in every bin
  double *bin;  // left range (minimal) coordinate of bins
  // counts[1] is amount of counts of data in the range bin[1] ... bin[2]
};


FitsHistoData::FitsHistoData(): nbin(0),count(0),bin(0) {}

FitsHistoData::FitsHistoData(int n, int *h, double *t): wxObjectRefData(),
  nbin(n),count(h),bin(t)
{
  wxASSERT(count && bin);
}


FitsHistoData::FitsHistoData(const FitsHistoData& other): wxObjectRefData()
{
  nbin = other.nbin;
  count = new int[nbin];
  bin = new double[nbin];
  std::copy(other.count,other.count+nbin,count);
  std::copy(other.bin,other.bin+nbin,bin);
}

FitsHistoData& FitsHistoData::operator = (const FitsHistoData& other)
{
  wxFAIL_MSG("*** FitsHistoData: WE ARE REALLY NEED ASSIGNMENT CONSTRUCTOR");
  if( this != & other ) {
    ;
  }
  return *this;
}


FitsHistoData::~FitsHistoData()
{
  delete[] count;
  delete[] bin;
}


FitsHisto::FitsHisto(): wbin(0.0),xmin(0.0),xmax(0.0),cmax(0),cmin(0),ntotal(0) {}

FitsHisto::FitsHisto(long n, const float *d, double q, int nb)
{
  const double macheps = std::numeric_limits<double>::epsilon();
  wxASSERT(n > 0 && d && -macheps < q && q < 1.0+macheps);

  // http://en.wikipedia.org/wiki/Histogram
  // By default, 128 bins is used, that corresponds,
  // by 1/3-rule to 2**18 elements

  int nbin = 128;
  if( nb > 0 )
    nbin = nb;
  else
    nbin = int(2.0*pow(float(n),1.0/3.0)) + 1;

  EmpiricalCDF ecdf(n,d);

  // quantiles are suitable for well distributed data
  xmin = ecdf.GetQuantile(1.0-q);
  xmax = ecdf.GetQuantile(q);

  // fall-back for a few datapoints
  if( xmax - xmin < macheps ) {
    xmin = std::numeric_limits<double>::max();
    xmax = std::numeric_limits<double>::lowest();
    for(int i = 0; i < n; i++) {
      if( isfinite(d[i]) ) {
	if( xmin > d[i] ) xmin = d[i];
	if( xmax < d[i] ) xmax = d[i];
      }
    }
  }

  // equal values
  if( xmax - xmin < macheps ) {
    float x = d[0];
    if( fabs(x) > macheps ) {
      xmin = x / 2;
      xmax = 2 * x;
    }
    else {
      xmin = 0;
      xmax = 1;
    }
  }

  wbin = (xmax - xmin) / nbin;
  Create(nbin,n,d);

  wxLogDebug("FitsHisto(%ld,%p,%lf,%d): %lf %lf %d %lf %f %d",
	     n,d,q,nb,xmin,xmax,int(nbin),wbin,
	     float(ecdf.GetQuantile(0.5)),int(n));
}

FitsHisto::FitsHisto(const FitsHisto& a):
  wbin(a.wbin),xmin(a.xmin),xmax(a.xmax),cmax(a.cmax),cmin(a.cmin),ntotal(a.ntotal)
{
  wxASSERT(a.IsOk());
  Ref(a);
}

FitsHisto& FitsHisto::operator = (const FitsHisto& other)
{
  if( this != & other ) {
    wbin = other.wbin;
    xmin = other.xmin;
    xmax = other.xmax;
    cmin = other.cmin;
    cmax = other.cmax;
    ntotal = other.ntotal;
    Ref(other);
  }
  return *this;
}

wxObjectRefData *FitsHisto::CreateRefData() const
{
  return new FitsHistoData;
}


wxObjectRefData *FitsHisto::CloneRefData(const wxObjectRefData *other) const
{
  return new FitsHistoData(* (FitsHistoData *) other);
}


bool FitsHisto::IsOk() const
{
  const FitsHistoData *data = static_cast<const FitsHistoData *>(m_refData);
  return data && data->nbin > 0 && data->count && data->bin;
}

void FitsHisto::Create(int nbin, int n, const float *d)
{
  wxASSERT(wbin > std::numeric_limits<float>::epsilon() && d && n > 0);

  int *count = new int[nbin];
  double *bin = new double[nbin];

  for(int i = 0; i < nbin; i++) count[i] = 0;
  for(int i = 0; i < nbin; i++) bin[i] = xmin + i*wbin;

  ntotal = 0;
  for(long i = 0; i < n; i++ ) {
    int m = int((d[i] - xmin)/wbin);
    if( 0 <= m  && m < nbin ) {
      count[m]++;
      ntotal++;
    }
  }
  SetRefData(new FitsHistoData(nbin,count,bin));

  cmax = 0;
  cmin = n;
  for(int i = 0; i < nbin; i++) {
    if( count[i] > cmax ) cmax = count[i];
    if( count[i] < cmin && count[i] > 0 ) cmin = count[i];
  }

#if wxDEBUG_LEVEL > 1
  wxFFileOutputStream file("/tmp/fitshisto","w+");
  wxTextOutputStream cout(file);
  for(int i = 0; i < nbin; i++)
    cout << i << " " << bin[i] << " " << count[i] << endl;
#endif
}

int FitsHisto::Count(int n) const
{
  const FitsHistoData *data = static_cast<const FitsHistoData *>(m_refData);
  wxASSERT(data && 0 <= n && n < data->nbin && data->count);
  return data->count[n];
}

float FitsHisto::Freq(int n) const
{
  const FitsHistoData *data = static_cast<const FitsHistoData *>(m_refData);
  wxASSERT(data && 0 <= n && n < data->nbin && data->count);

  return float(data->count[n]) / float(ntotal);
}

float FitsHisto::Bin(int n) const
{
  const FitsHistoData *data = static_cast<const FitsHistoData *>(m_refData);
  wxASSERT(data && 0 <= n && n < data->nbin && data->bin);
  return data->bin[n];
}

float FitsHisto::BinMin() const
{
  const FitsHistoData *data = static_cast<const FitsHistoData *>(m_refData);
  wxASSERT(data && data->nbin > 0 && data->bin && data->count);

  for(int i = 0; i < data->nbin-1; i++)
    if( data->count[i] > 0 )
      return data->bin[i];

  return data->bin[0];
}

float FitsHisto::BinMax() const
{
  const FitsHistoData *data = static_cast<const FitsHistoData *>(m_refData);
  wxASSERT(data && data->nbin > 0 && data->bin && data->count);

  for(int i = data->nbin-1; i > 0; i--)
    if( data->count[i] > 0 )
      return data->bin[i];

  return data->bin[data->nbin-1];
}

float FitsHisto::BinWidth() const
{
  return wbin;
}


int FitsHisto::NBins() const
{
  const FitsHistoData *data = static_cast<const FitsHistoData *>(m_refData);
  if( data )
    return data->nbin;
  else
    return 0;
}

int FitsHisto::MaxCount() const
{
  return cmax;
}

int FitsHisto::MinCount() const
{
  return cmin;
}

float FitsHisto::MaxFreq() const
{
  return float(cmax) / float(ntotal);
}

float FitsHisto::MinFreq() const
{
  return float(cmin) / float(ntotal);
}

float FitsHisto::MinVal() const
{
  return xmin;
}

float FitsHisto::MaxVal() const
{
  return xmax;
}

long FitsHisto::Ntotal() const
{
  return ntotal;
}
