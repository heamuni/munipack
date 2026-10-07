/*

  XMunipack - FITS display

  Copyright © 2022 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _FITSDISPLAY_H_
#define _FITSDISPLAY_H_


#include "fits.h"
#include <fitsio.h>
#include <wx/wx.h>
#include <vector>
#include <algorithm>
#include <cmath>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

#ifdef __WXOSX__
#define __ADOBERGB__
#else
#define __sRGB__
#endif

class FitsDisplay
{
  bool is_colour;
  int shrink;
  FitsTone tone;
  FitsItt itt;
  FitsPalette pal;
  FitsColour colour;

  int ngamma;
  float gstep;
  unsigned char *gamma;

  void MakeGamma();
  inline unsigned char Gamma(float) const;
  inline void XYZ_RGB(float,float,float, unsigned char *) const;

  void GetGrey(int,int,int,const float*,int,int,int,int,int,int,int,int,
	       unsigned char *) const;
  void GetColour(int,int,int,const float*,int,int,int,int,int,int,int,int,
		 unsigned char *) const;

public:
  FitsDisplay(bool);
  FitsDisplay(const FitsDisplay&);
  FitsDisplay& operator = (const FitsDisplay&);
  virtual ~FitsDisplay();

  virtual void SetPalette(const FitsPalette&);
  virtual void SetTone(const FitsTone&);
  virtual void SetItt(const FitsItt& );
  virtual void SetColour(const FitsColour& );
  virtual void SetShrink(int);
  virtual void GetRGB(int,int,int,const float*,int,int,int,int,int,int,int,int,
		      unsigned char *) const;
  virtual wxImage MakeIcon(const FitsArray&, int, int) const;
  virtual wxImage ConvertTowxImage(const FitsArray&) const;

  FitsTone GetTone() const { return tone; }
  FitsItt GetItt() const { return itt; }
  FitsPalette GetPalette() const { return pal; }
  FitsColour GetColour() const { return colour; }
  int GetShrink() const { return shrink; }

  inline float GammaFun(float) const;

};

inline unsigned char FitsDisplay::Gamma(float r) const
{
  int n = r / gstep;
  return n < 0 ? gamma[0] : (n >= ngamma ? gamma[ngamma] : gamma[n]);
  /*
  if( n < 0 )
    return gamma[0];
  else if( n >= ngamma )
    return gamma[ngamma];
  else {
    wxASSERT(0 <= n && n < ngamma);
    return gamma[n];
  }
  */
}

#ifdef __sRGB__

inline float FitsDisplay::GammaFun(float r) const
{
  const float q = 1.0 / 2.4;
  const float a = 0.055;
  const float a1 = 1.0 + a;

  return r < 0.0031308f ? 12.92f*r : a1*powf(r,q) - a;
}

inline void FitsDisplay::XYZ_RGB(float X,float Y, float Z, unsigned char *rgb) const
{
  // http://en.wikipedia.org/wiki/SRGB
  rgb[0] = Gamma( 3.2406e-2f*X - 1.5372e-2f*Y - 0.4986e-2f*Z);
  rgb[1] = Gamma(-0.9689e-2f*X + 1.8758e-2f*Y + 0.0415e-2f*Z);
  rgb[2] = Gamma( 0.0557e-2f*X - 0.2040e-2f*Y + 1.0570e-2f*Z);
}

#endif

#ifdef __ADOBERGB__

inline float FitsDisplay::GammaFun(float r) const
{
  const float q = 1.0/2.19921875;
  return powf(r,q);
}

inline void FitsDisplay::XYZ_RGB(float X, float Y, float Z, unsigned char *rgb) const
{
  // http://en.wikipedia.org/wiki/Adobe_RGB_color_space
  // http://www.adobe.com/digitalimag/pdfs/AdobeRGB1998.pdf
  rgb[0] = Gamma( 2.04159e-2f*X - 0.56501e-2f*Y - 0.34473e-2f*Z);
  rgb[1] = Gamma(-0.96924e-2f*X + 1.87597e-2f*Y + 0.04156e-2f*Z);
  rgb[2] = Gamma( 0.01344e-2f*X - 0.11836e-2f*Y + 1.01517e-2f*Z);
}
#endif

#endif
