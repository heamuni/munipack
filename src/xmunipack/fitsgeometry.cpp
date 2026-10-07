/*

  xmunipack - fits image geometry

  Copyright © 2018-2022 F.Hroch (hroch@physics.muni.cz)

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


// ------------   FitsGeometry


FitsGeometry::FitsGeometry(const FitsArray& a): FitsArray(a)
{
  wxASSERT(IsOk());
}

void FitsGeometry::ShrinkSubwin(int shrink, int x, int y, int w, int h,
				int iwidth, int iheight, float *idata) const
{
  //wxLogDebug("FitsGeometry::ShrinkSubwin");

  wxASSERT(shrink >= 1 && iwidth > 0 && iheight > 0 && x >= 0 && y >= 0 &&
	   w > 0 && h > 0);

  const float *array = PixelData();
  int width = GetWidth();
  int height = GetHeight();
  size_t npixels = width*height;
  size_t ipixels = iwidth*iheight;

  const int imin = x / shrink;
  const int jmin = y / shrink;
  const int imax = min(imin + w / shrink,iwidth);
  const int jmax = min(jmin + h / shrink,iheight);

  for(int k = 0; k < GetDepth(); k++) {
    const float *src = array + k*npixels;
    float *dst = idata + k*ipixels;
    for(int j = jmin; j < jmax; j++) {
      int j0 = shrink*j;
      float *drow = dst + j*iwidth;
      for(int i = imin; i < imax; i++) {
	int i0 = shrink*i;
	drow[i] = MeanSquare(src,width,height,i0,j0,shrink);
      }
    }
  }
}



FitsArray FitsGeometry::GetSubArray(int xoff, int yoff, int w, int h)
{
  long npixels = w*h;
  long *ns = new long[2];
  float *a = new float[npixels];

  ns[0] = w;
  ns[1] = h;

  const float *array = PixelData();
  wxASSERT(array);
  int width = GetWidth();

  switch(Naxis()) {

  case 2:

    /*
    wxLogDebug("FitsGeometry::GetSubArray %d %d %d %d %d %d",xoff,yoff,w,h,
	       int(data->naxes[0]),int(data->naxes[1]));
    */
    wxASSERT(xoff >= 0 && yoff >= 0 &&
	     xoff + w <= GetWidth() && yoff + h <= GetHeight());

    for(int j = 0; j < h; j++) {
      const float *src = array + (yoff + j)*width + xoff;
      copy(src,src+w,a+j*w);
    }

    break;

  default:

    wxFAIL_MSG("FitsGeom::SubArray arbitrary dimension isn't implemented yet.");
    break;
  }

  return FitsArray(header,GetHduType(),2,ns,a);
}

void FitsGeometry::SetSubArray(int xoff, int yoff, const FitsArray& sub)
{
  // JUST ONLY *TWO* DIMENSIONS ARE IMPLEMENTED

  wxASSERT(sub.Naxis() == 2 && Naxis() == 2);
  float *array = GetData();
  wxASSERT(array);
  int width = GetWidth();

  /*
  wxLogDebug("FitsGeometry::SetSubArray X: %d %d %d",
	     (int)xoff,(int)sub.Naxes(0),int(data->naxes[0]));
  wxLogDebug("FitsGeometry::SetSubArray Y: %d %d %d",
	     (int)yoff,(int)sub.Naxes(1),int(data->naxes[1]));
  */
  wxASSERT(0 <= xoff && (xoff + sub.Naxes(0)) <= GetWidth());
  wxASSERT(0 <= yoff && (yoff + sub.Naxes(1)) <= GetHeight());

  int w = sub.GetWidth();
  int h = sub.GetHeight();
  const float *a = sub.PixelData();
  wxASSERT(a);

  /*
  wxLogDebug("%d %d %d %d %d %d %d %d",(int)h,(int)sub.Naxes(0),(int)w,
	     (int)sub.Naxes(1),(int)data->naxes[0],(int)xoff,
	     (int)data->naxes[1],(int)yoff);
  */

  for(int j = 0; j < h; j++) {
    const float *src = a + j*w;
    float *dest = array + (yoff + j)*width + xoff;
    copy(src,src+w,dest);
  }
}


float FitsGeometry::MeanLine(int x, int w) const
{
  const float *array = PixelData();
  wxASSERT(array && Naxis() == 1);

  int i0 = x;
  int dx = w / 2;

  int i1 = max(i0 - dx,0);
  int i2 = min(long(i0 + dx),Naxes(0));

  float s = 0;
  float n = 0;
  for(int i = i1; i <= i2; i++) {
    // s = s + Pixel(i);
    // wxASSERT(0 <= i && i < data->npixels);
    //    s = s + *(data->array + i);
    s = s + *(array + i);
    n = n + 1;
  }

  if( n > 0 ) return s/n; else return 0.0;
}


float FitsGeometry::MeanRect(int x, int y, int w, int h) const
{
  return MeanRect(PixelData(),Naxes(0),Naxes(1),x,y,w,h);
}

float FitsGeometry::MeanRect(const float *data, long width, long height,
			     int x, int y, int w, int h) const
{
  int i0 = x;
  int j0 = y;
  int dx = w / 2;
  int dy = h / 2;

  int i1 = max(i0 - dx,0);
  int j1 = max(j0 - dy,0);
  int i2 = min(i0 + dx,int(width)-1);
  int j2 = min(j0 + dy,int(height)-1);

  float s = 0.0;
  int n = (j2 - j1 + 1)*(i2 - i1 + 1);
  const float *ai = data + i1;
  for(int j = j1; j <= j2; j++) {
    const float *a = ai + j*width;
    for(int i = i1; i <= i2; i++)
      s += *a++;
  }
  return s / n;
}

float FitsGeometry::MeanSquare(const float *pic, int width, int height, int x, int y, int d) const
{
  wxASSERT(pic && x >= 0 && y >= 0 && d >= 1);

  const int i1 = x;
  const int j1 = y;
  const int i2 = min(x+d,width);
  const int j2 = min(y+d,height);

  float s = 0.0;
  for(int j = j1; j < j2; j++) {
    const float *row = pic + j*width;
    for(int i = i1; i < i2; i++)
      s += row[i];
  }
  int n = (j2 - j1)*(i2 - i1);
  return n > 0 ? s / n : 0;
}


float FitsGeometry::MeanRect_debug(int x, int y, int w, int h) const
{
  const float *array = PixelData();
  wxASSERT(array && Naxis() == 2);
  int width = Naxes(0);
  //  int height = Naxes(1);

  /*
  const int dmax = max(w,h)/2;

  int i0 = int(x + 0.5);
  int j0 = int(y + 0.5);
  int dx = min(int(w/2.0 - 0.5),dmax);
  int dy = min(int(h/2.0 - 0.5),dmax);

  wxASSERT(dx >= 0 && dy >= 0);

  int i1 = max(i0 - dx,0);
  int j1 = max(j0 - dy,0);
  int i2 = min(long(i0 + dx),data->naxes[0]-1);
  int j2 = min(long(j0 + dy),data->naxes[1]-1);
  */

  int i0 = x;
  int j0 = y;
  int dx = w / 2;
  int dy = h / 2;

  wxASSERT(dx >= 0 && dy >= 0);

  int i1 = max(i0 - dx,0);
  int j1 = max(j0 - dy,0);
  int i2 = min(long(i0 + dx),Naxes(0)-1);
  int j2 = min(long(j0 + dy),Naxes(1)-1);
  //  int i2 = min(long(i0 + dx),data->naxes[0]-1);
  //  int j2 = min(long(j0 + dy),data->naxes[1]-1);


  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //  if( dx > 0 && dy > 0 )
  //    wxLogDebug(_("%d %d"),dx,dy);

//   int ww = i2 - i1 + 1;
//   int hh = j2 - j1 + 1;
//   ww = dx;
//   hh = dy;
//   if(  ww > 0 && hh > 0 ) {
//     float a[ww*hh];
//     const float *array = data->array;
//     wxASSERT(array);
//     for(int j = 0; j < hh; j++)
//       memcpy(a+j*ww,array+(j1+j)*data->naxes[0]+i1,ww*sizeof(float));

//     //  return a[0];
//     return FitsArrayStat::QMed(ww*hh,a,max((ww*hh)/2,0));
//   }
//   else
//     return Pixel(i1,j1);


  ///!!!!!!!!!!!!!!!!!!!!!!!!!!!

  float s = 0.0;

#ifdef __WXDEBUG__
  int n = 0;
#else
  int n = (j2 - j1 + 1)*(i2 - i1 + 1);
#endif

  //  float *ai = data->array + i1;
  const float *ai = array + i1;
  //  int idim = data->naxes[0];
  int idim = width;
  for(int j = j1; j <= j2; j++) {
    //    float *a = data->array + j*data->naxes[0] + i1;
    const float *a = ai + j*idim;
    for(int i = i1; i <= i2; i++) {
      //      s = s + Pixel(i,j);
      //      wxASSERT(0 <= i && i < data->npixels);
      s += *a++;
      //      s = s + *(a + i);
      //      wxLogDebug("%f",(float)*(a + i));
#ifdef __WXDEBUG__
      n = n + 1;
#endif
    }
  }
#ifdef __WXDEBUG__
  int nx = j2 - j1;
  int ny = i2 - i1;
  int nn = (nx + 1)*(ny + 1);
  if( n != nn )
    wxLogDebug("FitsGeometry::MeanRect: %d %d %d %d",n,nn,nx,ny);
#endif
  wxASSERT(n == (j2 - j1 + 1)*(i2 - i1 +1) );
  if( n > 0 ) return s/n; else return 0.0;
}
