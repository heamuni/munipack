/*

  xmunipack - fits array

  Copyright © 2021-22 F.Hroch (hroch@physics.muni.cz)

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
#include <algorithm>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#include <wx/txtstrm.h>
#include <wx/wfstream.h>
#endif


using namespace std;


// ------------   FitsArrayData

class FitsArrayData : public wxObjectRefData
{

public:
  FitsArrayData();
  FitsArrayData(int, long *, float *);
  FitsArrayData(const FitsArrayData&);
  virtual ~FitsArrayData();

  FitsArrayData& operator = (const FitsArrayData&);
  bool operator == (const FitsArrayData&) const;

  int naxis;
  long *naxes;
  float *array;
  // int *array;
  // double *array;

};

FitsArrayData::FitsArrayData(): naxis(0),naxes(0),array(0)
{
  //  wxLogDebug("FitsArrayData::FitsArrayData() %p",array);
}

FitsArrayData::FitsArrayData(int n, long *ns, float *a): wxObjectRefData(),
  naxis(n),naxes(ns),array(a)
{
  wxASSERT(naxes && array);
  //  wxLogDebug("FitsArrayData being CREATED... %p %d",array,int(naxes[0]));
}


FitsArrayData::FitsArrayData(const FitsArrayData& other): wxObjectRefData()
{
  // COPY constructor
  long n = 1;
  for(long i = 0; i < other.naxis; i++)
    n = n * other.naxes[i];
  naxis = other.naxis;
  naxes = new long[naxis];
  array = new float[n];
  copy(other.naxes,other.naxes+naxis,naxes);
  copy(other.array,other.array+n,array);
  //  wxLogDebug("FitsArrayData COPY %p",array);
}

FitsArrayData& FitsArrayData::operator = (const FitsArrayData& other)
{
 wxFAIL_MSG("*** FitsArrayData: WE ARE REALLY NEED ASSIGNMENT CONSTRUCTOR");
 if( this != & other ) {
   ;
 }
 return *this;
}

bool FitsArrayData::operator == (const FitsArrayData& other) const
{
  if( naxis != other.naxis )
    return false;

  int isum = 0;
  for(int i = 0; i < naxis; i++)
    isum = isum + abs(naxes[i] - other.naxes[i]);

  if( isum > 0 )
    return false;

  long n = 1;
  for(int i = 0; i < naxis; i++)
    n = n * naxes[i];

  for(int i = 0; i < n; i++) {
    if( fabs(array[i] - other.array[i]) > FLT_EPSILON )
      return false;
  }

  return true;
  //return naxis == other.naxis && naxes == other.naxes && array = other.array;
}

FitsArrayData::~FitsArrayData()
{
  //  wxLogDebug("FitsArrayData being DELETED... %p %d",array,int(naxes[0]));
  delete[] naxes;
  delete[] array;
}


// FitsArray() ---------------------------------------------------------

FitsArray::FitsArray(const FitsHeader& h, long width, long height):
  FitsHdu(h,IMAGE_HDU,HDU_IMAGE,HDU_IMAGE_FRAME)
{
  wxASSERT(width > 0 && height > 0 && m_refData == 0);

  int naxis = 2;
  long *naxes = new long[naxis];
  naxes[0] = width;
  naxes[1] = height;
  float *array = new float[width*height];
  //  fill(array,array+width*height,0);
  SetRefData(new FitsArrayData(naxis,naxes,array));

  /* IMPLEMENT update of header array strings:
   * copy on write the full header
   * modify image dimensions
   * add sub-window coordinates
   */
}

FitsArray::FitsArray(const FitsHeader& h, long width, long height, long depth):
  FitsHdu(h,IMAGE_HDU,HDU_IMAGE,HDU_IMAGE_COLOUR)
{
  wxASSERT(width > 0 && height > 0 && depth > 0 && m_refData == 0);

  int naxis = 3;
  long *naxes = new long[naxis];
  naxes[0] = width;
  naxes[1] = height;
  naxes[2] = depth;
  float *array = new float[width*height*depth];
  //  fill(array,array+width*height*depth,0);
  SetRefData(new FitsArrayData(naxis,naxes,array));
}

FitsArray::FitsArray(const FitsHeader& h, int t, int n, long *ns, float *a):
  FitsHdu(h,t)
{
  wxASSERT(n >= 0 && ns && a && m_refData == 0);
  SetRefData(new FitsArrayData(n,ns,a));
}

FitsArray::FitsArray(const FitsArray& a):
  FitsHdu(a.GetHeader(),a.GetHduType())
{
  //  wxASSERT(a.IsOk());
  if( a.IsOk() )
    Ref(a);
  //  wxLogDebug("FitsArray::FitsArray(const FitsArray&) %p %d",m_refData,int(m_refData->GetRefCount()));
}

FitsArray::FitsArray(const FitsHdu& h):
  FitsHdu(h.GetHeader(),h.GetHduType())
{
  //  wxLogDebug("FitsArray::FitsArray HDU");
  wxASSERT(h.Type() == HDU_IMAGE && h.IsOk());
  Ref(h);
  //  wxLogDebug("FitsArray::FitsArray(const FitsHdu&) %p %d",m_refData,int(m_refData->GetRefCount()));
}

wxObjectRefData *FitsArray::CreateRefData() const
{
  wxLogDebug("FitsArray::CreateRefData()");
  return new FitsArrayData;
}


wxObjectRefData *FitsArray::CloneRefData(const wxObjectRefData *other) const
{
  wxLogDebug("FitsArray::CloneRefData");
  //  const FitsArrayData *other = static_cast<const FitsArrayData *>(that);
  //  wxASSERT(other);
  return new FitsArrayData(* (FitsArrayData *) other);
}


FitsArray::~FitsArray()
{
  /*
  const FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  if( data )
    wxLogDebug("FitsArray::~FitsArray() %d: %ld %ld %d "+GetKey("DATE-OBS"),
	       int(data->GetRefCount()),GetWidth(),GetHeight(),int(Npixels()));
  */
}


FitsArray& FitsArray::operator = (const FitsArray& other)
{
  //  wxFAIL_MSG("*** FitsArray: WE ARE REALLY NEED ASSIGNMENT CONSTRUCTOR");
  if( this != & other ) {
    header = other.header;
    htype = other.htype;
    type = other.type;
    flavour = other.flavour;
    modified = other.modified;
    Ref(other);
  }
  return *this;
}


bool FitsArray::operator == (const FitsArray& other) const
{
  if (m_refData == other.m_refData)
    return true;

  if (!m_refData || !other.m_refData)
    return false;

  return ( *(FitsArrayData*)m_refData == *(FitsArrayData*)other.m_refData );
}

bool FitsArray::operator != (const FitsArray& other) const
{
  return  !(*this == other);
}

bool FitsArray::IsOk() const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  return m_refData && data && data->naxis > 0 && data->naxes && data->array;
}


int FitsArray::Naxis() const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  //  wxASSERT(data);
  if( data )
    return data->naxis;
  else
    return 0;
}

long FitsArray::Naxes(int n) const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  //  wxASSERT(data);

  if( data && 0 <= n && n < data->naxis )
    return data->naxes[n];
  else
    return 0;
}

long* FitsArray::Naxes() const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  //  wxASSERT(data);

  if( data ) {
    long *ns = new long[data->naxis];
    for(int n = 0; n < data->naxis; n++)
      ns[n] = data->naxes[n];
    return ns;
  }
  else
    return 0;
}

long FitsArray::Npixels() const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  wxASSERT(data);

  long npixels = 1;
  for(int k = 0; k < data->naxis; k++)
    npixels = npixels*data->naxes[k];
  return npixels;
}

float FitsArray::Pixel(int x) const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  wxASSERT(data && data->array && data->naxes && data->naxis == 1);
  wxASSERT(0 <= x && x < data->naxes[0]);
  return *(data->array + x);
  /*
  float *array = PixelData();
  wxASSERT(array && 0 <= i && i < Npixels());
  return *(array + i);
  */
}

float FitsArray::Pixel(int x, int y) const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  wxASSERT(data && data->array && data->naxes && data->naxis == 2);
  wxASSERT(0 <= x && x < data->naxes[0] && 0 <= y && y < data->naxes[1]);
  return *(data->array+y*data->naxes[0] + x);
}

float FitsArray::Pixel(int x, int y, int z) const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  wxASSERT(data && data->array && data->naxes && data->naxis == 3);
  wxASSERT(0 <= x && x < data->naxes[0] && 0 <= y && y < data->naxes[1] &&
	   0 <= z && z < data->naxes[2]);
  return *(data->array + (z*data->naxes[1] + y)*data->naxes[0] + x);
}

const float *FitsArray::PixelData() const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  wxASSERT(data);
  const float *array = data->array;
  wxASSERT(array);
  return array;
}

float *FitsArray::GetData() const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  wxASSERT(data);
  float *array = data->array;
  wxASSERT(array);
  return array;
}

FitsArray FitsArray::GetGrid(int n, int stride) const
{
  FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  wxASSERT(data && (Naxis() == 2 || Naxis() == 3));
  long npixels = Naxes(0) * Naxes(1);
  const float *array = data->array + n*npixels;

  int np = 2;
  long *ns = new long[np];
  ns[0] = max(Naxes(0) / stride, 1L);
  ns[1] = max(Naxes(1) / stride, 1L);
  long npix = ns[0]*ns[1];
  float *a = new float[npix];
  for(  long j = 0, k = 0; k < ns[1]; j = j + stride, k++)
    for(long i = 0, l = 0; l < ns[0]; i = i + stride, l++) {
      wxASSERT(k*ns[0] + l < npix && j*Naxes(0) + i < npixels);
      a[k*ns[0] + l] = array[j*Naxes(0) + i];
      //      wxLogDebug("%ld %ld %ld %ld %f",k,l,i,j,a[k*ns[0] + l]);
    }
  return FitsArray(header,GetHduType(),np,ns,a);
}


FitsArray FitsArray::Plane(int n) const
{
  wxLogDebug("FitsArray::Plane");
  wxASSERT(0 <= n && n < 3);

  int np = 2;
  long *ns = new long[np];
  ns[0] = Naxes(0);
  ns[1] = Naxes(1);
  long npixels = ns[0]*ns[1];

  const FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);
  wxASSERT(data && data->array);
  const float *array = data->array + n*npixels;

  float *a = new float[npixels];
  copy(array,array+npixels,a);
  return FitsArray(header,GetHduType(),np,ns,a);
}


void FitsArray::GetSubImage(int x, int y, FitsArray& subwin) const
{
  int depth = GetDepth();
  int width = GetWidth();
  int w = subwin.GetWidth();
  int h = subwin.GetHeight();
  int npixels = w*h;
  int mpixels = GetWidth()*GetHeight();
  const FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);

  for(int k = 0; k < depth; k++) {
    float *src = data->array + x + k*mpixels;
    float *dest = subwin.GetData() + k*npixels;
    for(int j = 0; j < h; j++) {
      float *s = src + (y+j)*width;
      float *d = dest + j*w;
      copy(s,s+w,d);
    }
  }
}

void FitsArray::SetSubImage(int x, int y, const FitsArray& subwin)
{
  int depth = GetDepth();
  int width = GetWidth();
  int w = subwin.GetWidth();
  int h = subwin.GetHeight();
  int npixels = w*h;
  int mpixels = GetWidth()*GetHeight();
  const FitsArrayData *data = static_cast<FitsArrayData *>(m_refData);

  for(int k = 0; k < depth; k++) {
    float *src = subwin.GetData() + k*npixels;
    float *dest = data->array + x + k*mpixels;
    for(int j = 0; j < h; j++) {
      float *s = src + j*w;
      float *d = dest + (y+j)*width;
      copy(s,s+w,d);
    }
  }
}

bool FitsArray::HasWCS() const
{
  return IsOk() &&
    GetKey("CTYPE1").Find("TAN") != wxNOT_FOUND &&
    GetKey("CTYPE2").Find("TAN") != wxNOT_FOUND &&
    GetKey("CRPIX1") != ""  && GetKey("CRPIX1") != "" &&
    GetKey("CRVAL1") != ""  && GetKey("CRVAL1") != "" &&
    GetKey("CD1_1") != "";
}
