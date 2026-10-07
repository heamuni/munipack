/*

  xmunipack - fits image display

  Copyright © 2018-22 F.Hroch (hroch@physics.muni.cz)

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

#include "fitsdisplay.h"
#include <cfloat>
#include <wx/wx.h>
#include <cmath>
#include <algorithm>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

#include "Luv_Gamut.h"

#define RAD 57.29577951308232


// ---- FitsDisplay  ------


FitsDisplay::FitsDisplay(bool is): is_colour(is), shrink(1),
				   ngamma(0), gstep(1.0), gamma(0)
{
  if( is_colour )
    MakeGamma();
}

FitsDisplay::FitsDisplay(const FitsDisplay& display)
{
  wxFAIL_MSG("FitsDisplay ----- WE ARE REALY NEED COPY CONSTRUCTOR -----");
}

FitsDisplay& FitsDisplay::operator= (const FitsDisplay& other)
{
  wxFAIL_MSG("FitsDisplay ----- WE ARE REALY NEED ASSIGN CONSTRUCTOR -----");
  if( this != &other ) {
    ;
  }
  return *this;
}

FitsDisplay::~FitsDisplay()
{
  delete[] gamma;
}

void FitsDisplay::SetShrink(int s)
{
  shrink = s;
}


void FitsDisplay::SetItt(const FitsItt& i)
{
  itt = i;
}

void FitsDisplay::SetTone(const FitsTone& t)
{
  tone = t;
}

void FitsDisplay::SetPalette(const FitsPalette& p)
{
  pal = p;
}

void FitsDisplay::SetColour(const FitsColour& c)
{
  colour = c;
}

void FitsDisplay::GetRGB(int width, int height, int depth, const float *array,
			 int x, int y, int w, int h, int xoff, int yoff,
			 int iwidth, int iheight, unsigned char *rgb) const
{
  if( is_colour )
    GetColour(width,height,depth,array,x,y,w,h,xoff,yoff,iwidth,iheight,rgb);
  else
    GetGrey(width,height,depth,array,x,y,w,h,xoff,yoff,iwidth,iheight,rgb);
}


void FitsDisplay::GetGrey(int fwidth, int fheight, int depth, const float *fimage,
			     int x, int y, int w, int h, int xoff, int yoff,
			     int iwidth, int iheight, unsigned char *idata) const
{
  wxASSERT(idata && fimage && fwidth > 0 && fheight > 0 && depth == 1 &&
	   xoff >= 0 && yoff >= 0 && 0 <= x && x < fwidth && 0 <= y && y < fheight);

  const int j1 = y;
  const int j2 = y + h;
  const int i1 = x;
  const int i2 = x + w;
  wxASSERT(i2 <= fwidth && j2 <= fheight);
  /*
  wxLogDebug("FITS x:%d y:%d w:%d h:%d, BMP x:%d ymin:%d ymax:%d",
	     x,y,w,h,x-xoff,iheight - (j2-yoff),iheight-1 -(j1-yoff));
  */

  int istride = 3*iwidth;
  int y1 = iheight-1 + yoff;

  for(int j = j1; j < j2; j++ ) {
    int jj = y1 - j;
    wxASSERT(0<= jj && jj < iheight);

    unsigned char *iline = idata + jj*istride;
    const float *fline = fimage + j*fwidth;

    for(int i = i1; i < i2; i++ ) {
      int ii = i - xoff;
      wxASSERT(00 <= ii && ii < iwidth);
      unsigned char *p = iline + 3*ii;

      float f;
      f = tone.Scale(fline[i]);
      f = itt.Scale(f);

      float g = GammaFun(f);

      pal.RGB(g,p);
    }
  }

}



  // This part does colour transformation from an CIE 1931 XYZ colour space
  // to a display colour space. Optionaly, colour parameters are tuned.
  //
  // Digital cameras usually produces colours in sRGB or AdobeRGB.
  //
  // Astronomical filter series usually stores data in their proper systems.
  // Plain usage leads to a deformation of colours.
  // Documents
  //         http://www.fho-emden.de/~hoffmann/ciexyz29082000.pdf
  //     and http://www.fho-emden.de/~hoffmann/cielab03022003.pdf
  // detaily describes the general transformation.
  //
  // Therefore, we are converting the input to the XYZ colour space.
  // (XYZ should be obtained from any RAW digital photography with dcraw).
  // (XYZ should be obtained from Landolt's BVR by a custom transformation.)

  // http://en.wikipedia.org/wiki/RGB_color_spaces
  // http://en.wikipedia.org/wiki/CIE_1931_color_space
  // http://en.wikipedia.org/wiki/D65

  // Luminosity scaled XYZ is transformed to some RGB.
  // (Adobe RGB (Mac) or sRGB is not directly supported by wxWidgets.)

void FitsDisplay::GetColour(int width, int height, int depth,
			       const float *pixeldata, int x, int y, int w, int h,
			       int xoff, int yoff, int iwidth, int iheight,
			       unsigned char *rgb) const
{
  wxASSERT(rgb && pixeldata && width > 0 && height > 0 && depth == 3);

  float S = colour.GetSaturation();
  bool nite = colour.GetNiteVision();

  int j1 = y;
  int j2 = y + h;
  int i1 = x;
  int i2 = x + w;
  wxASSERT(i2 <= width && j2 <= height);

  int npixels = width*height;
  int n2pixels = 2*npixels;
  int istride = 3*iwidth;
  int y1 = iheight-1 + yoff;
  for(int j = j1; j < j2; j++ ) {
    int jj = y1 - j;
    wxASSERT(0<= jj && jj < iheight);
    unsigned char *iline = rgb + jj*istride;
    const float *pixeline = pixeldata + j*width;
    for(int i = i1; i < i2; i++ ) {
      int ii = i - xoff;
      wxASSERT(0 <= ii && ii < iwidth);
      unsigned char *p = iline + 3*ii;

      const float *f = pixeline + i;
      float Z = *f;
      float Y = *(f+npixels);
      float X = *(f+n2pixels);

      if( Z < 0.0f ) Z = 0.0f;
      if( X < 0.0f ) X = 0.0f;
      if( Y < 0.0f ) Y = 0.0f;

      if( Y > 0.0f ) {

	/*
	  To save colours during non-linear scaling,
	  XYZ coordinates are transformed to xyY.
	  Only Y has scaled intensity.
	*/

	float W = X + Y + Z;
	if( W > 0.0f ) {

	  // to xyY
	  float x = X / W;
	  float y = Y / W;

	  if( nite )
	    colour.NiteVision(X,Z,x,y,Y);

	  // pre-scale
	  Y = 100*tone.Scale(Y);
	  if( Y < 0.0f ) Y = 0.0f;

	  // back to XYZ
	  float Q = Y / y;
	  X = Q*x;
	  Z = Q*(1.0f - x - y);
	  if( X < 0.0f ) X = 0.0f;
	  if( Z < 0.0f ) Z = 0.0f;
	}

	if( Y > 0 ) {
	  float L,u,v;
	  colour.XYZ_Luv(X,Y,Z,L,u,v);

	  float c = hypotf(u,v);
	  float coshue = 0;
	  float sinhue = 0;
	  if( c > 0 ) {
	    coshue = u / c;
	    sinhue = v / c;
	  }

	  // keep colour saturation
	  float s = c / L;

	  L = 100.0f*itt.Scale(L/100.0f);

	  // recover colour saturation
	  c = s*L;

	  // user saturation
	  c = S*c;

	  // limit chroma to the valid gamut
	  int lum = round(L);
	  if( 0 <= lum && lum <= 100 && c > 0.0f ) {
	    float hue = atan2f(v,u);
	    int deg = round(RAD*hue);
	    wxASSERT(0 <= deg+180 && deg+180 < 361 && 0 <= lum && lum <= 100);
	    if( Luv_Gamut[180+deg][lum] < c ) {
	      // wxLogDebug("shot %f %f",Luv_Gamut[180+deg][lum],c);
	      //	      c = std::min(Luv_Gamut[180+deg][lum],c);
	      c = Luv_Gamut[180+deg][lum];
	    }
	  }
	  else
	    c = 0;


	  u = c*coshue;
	  v = c*sinhue;
	  colour.Luv_XYZ(L,u,v,X,Y,Z);
	}

	XYZ_RGB(X,Y,Z,p);
      }
      else {
	//	R = G = B = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
      }

    } // i
  } // j

}


void FitsDisplay::MakeGamma()
{
  ngamma = 2048;
  gamma = new unsigned char[ngamma+1];
  gstep = 1.0 / float(ngamma - 1);
  for(int i = 0; i < ngamma; i++) {
    float x = i*gstep;
    float g = GammaFun(x);
    gamma[i] = std::min(std::max(int(255*g+0.5),0),255);
    //    printf("%f %d %d\n",x,gtab[i],(int)(blim*Gamma(x)+0.5));
  }
  gamma[ngamma] = gamma[ngamma-1];
}


wxImage FitsDisplay::MakeIcon(const FitsArray& array, int maxwidth, int maxheight)
  const
{
  wxASSERT(array.IsOk() && array.GetWidth() > 0 && array.GetHeight() > 0 &&
	   maxwidth > 0 && maxheight > 0);

  int depth = array.GetDepth();
  int shrink = wxMax(array.GetWidth()/maxwidth,array.GetHeight()/maxheight);
  if( shrink == 0 ) shrink = 1;
  int w = wxMax(array.GetWidth() / shrink,1);
  int h = wxMax(array.GetHeight()/ shrink,1);
  wxASSERT(w > 0 && h > 0 && shrink > 0 && depth > 0);

  float *cdata = new float[w*h*depth];
  wxASSERT(cdata);

  FitsGeometry geometry(array);
  geometry.ShrinkSubwin(shrink,0,0,array.GetWidth(),array.GetHeight(),w,h,cdata);

  unsigned char *idata = (unsigned char *) malloc(3*w*h);
  wxASSERT(idata);

  GetRGB(w,h,depth,cdata,0,0,w,h,0,0,w,h,idata);
  delete[] cdata;

  return wxImage(w,h,idata);
}


wxImage FitsDisplay::ConvertTowxImage(const FitsArray& array) const
{
  wxASSERT(array.IsOk() && array.GetWidth() > 0 && array.GetHeight() > 0 &&
	   shrink > 0);

  int depth = array.GetDepth();
  int w = array.GetWidth() / shrink;
  int h = array.GetHeight() / shrink;
  wxASSERT(w > 0 && h > 0 && shrink > 0 && depth > 0);

  float *cdata = new float[w*h*depth];
  wxASSERT(cdata);

  unsigned char *idata = (unsigned char *) malloc(3*w*h);
  wxASSERT(idata);

  if( shrink > 1 ) {
    FitsGeometry geometry(array);
    geometry.ShrinkSubwin(shrink,0,0,array.GetWidth(),array.GetHeight(),w,h,cdata);
    GetRGB(w,h,depth,cdata,0,0,w,h,0,0,w,h,idata);
  }
  else {
    wxASSERT(w == array.GetWidth() && h == array.GetHeight());
    GetRGB(w,h,depth,array.PixelData(),0,0,w,h,0,0,w,h,idata);
  }
  delete[] cdata;

  return wxImage(w,h,idata);
}
