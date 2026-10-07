/*

  xmunipack - palettes

  Copyright © 1997-2011, 2018-22 F.Hroch (hroch@physics.muni.cz)

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

using namespace std;


// ---- FitsPalette

FitsPalette::FitsPalette(int n, bool i):
  pal(n),inverse(i),npal(0),npal1(-1),rpal(0),gpal(0),bpal(0),pals(0)
{
  CreatePalette();
}

FitsPalette::FitsPalette(const FitsPalette& p)
{
  pal = p.pal;
  inverse = p.inverse;
  npal = p.npal;
  npal1 = p.npal1;
  pals = new unsigned char[3*npal];
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];
  std::copy(p.pals,p.pals+3*npal,pals);
  std::copy(p.rpal,p.rpal+npal,rpal);
  std::copy(p.gpal,p.gpal+npal,gpal);
  std::copy(p.bpal,p.bpal+npal,bpal);
}

FitsPalette& FitsPalette::operator = (const FitsPalette& other)
{
 if( this != & other ) {

   delete[] rpal;
   delete[] gpal;
   delete[] bpal;
   delete[] pals;

   pal = other.pal;
   inverse = other.inverse;
   npal = other.npal;
   npal1 = other.npal1;
   pals = new unsigned char[3*npal];
   rpal = new unsigned char[npal];
   gpal = new unsigned char[npal];
   bpal = new unsigned char[npal];

   std::copy(other.pals,other.pals+3*npal,pals);
   std::copy(other.rpal,other.rpal+npal,rpal);
   std::copy(other.gpal,other.gpal+npal,gpal);
   std::copy(other.bpal,other.bpal+npal,bpal);
 }
 return *this;
}


FitsPalette::~FitsPalette()
{
  delete[] rpal;
  delete[] gpal;
  delete[] bpal;
  delete[] pals;
}

void FitsPalette::SetInverse(bool n)
{
  inverse = n;
  CreatePalette();
}

bool FitsPalette::GetInverse() const
{
  return inverse;
}

int FitsPalette::GetPalette() const
{
  return pal;
}

wxString FitsPalette::GetPalette_str() const
{
  return Type_str(pal);
}

int FitsPalette::GetColours() const
{
  return npal;
}

wxString FitsPalette::Type_str(int n)
{
  switch(n) {
  case PAL_GREY:   return "Grey";
  case PAL_SEPIA:  return "Sepia";
  case PAL_ROYAL:  return "Royal";
  case PAL_COLOUR: return "Colour";
  case PAL_HIGHLIGHT: return "Highlight";
  case PAL_RAINBOW:   return "Rainbow";
  case PAL_MADNESS:   return "Madness";
  case PAL_VGA:    return "VGA";
  case PAL_AIPS0:  return "AIPS0";
  case PAL_FIRST:
  case PAL_LAST:
  default:
    wxFAIL_MSG("Undefined type of palette.");
    return "Unknown";
  }
}

wxArrayString FitsPalette::Type_str()
{
  wxArrayString a;

  for(int i = PAL_FIRST+1; i < PAL_LAST; i++)
    a.Add(Type_str(i));

  return a;
}


void FitsPalette::SetPalette(int l)
{
  wxASSERT(PAL_FIRST < l && l < PAL_LAST);
  pal = l;
  CreatePalette();
}


unsigned char FitsPalette::R(int i) const
{
  wxASSERT(rpal && 0 <= i && i < npal);
  return rpal[i];
}

unsigned char FitsPalette::G(int i) const
{
  wxASSERT(gpal && 0 <= i && i < npal);
  return gpal[i];
}

unsigned char FitsPalette::B(int i) const
{
  wxASSERT(bpal && 0 <= i && i < npal);
  return bpal[i];
}


unsigned char *FitsPalette::RGB(long n, float *f)
{
  unsigned char *rgb = (unsigned char *) malloc(3*n);
  unsigned char *p = rgb;
  int m = npal - 1;

  for(long i = 0; i < n; i++) {
    int l = max(min(int(255.0f * f[i] + 0.5f),m),0);
    *p++ = rpal[l];
    *p++ = gpal[l];
    *p++ = bpal[l];
  }

  return rgb;
}


void FitsPalette::CreatePalette()
{
  switch (pal) {
  case PAL_SEPIA:  Create_Sepia(); break;
  case PAL_ROYAL:   Create_Royal(); break;
  case PAL_COLOUR:  Create_Colour(); break;
  case PAL_HIGHLIGHT:  Create_Highlight(); break;
  case PAL_RAINBOW:    Create_Rainbow();break;
  case PAL_MADNESS:    Create_Madness();break;
  case PAL_VGA:    Create_VGA();   break;
  case PAL_AIPS0:  Create_AIPS0(); break;
  case PAL_FIRST:
  case PAL_LAST:
  default:  Create_Grey();
  }

  if( inverse ) {
    for(int i = 0; i < npal/2; i++) {
      unsigned char x;
      int l = npal - i -1;
      x = rpal[i]; rpal[i] = rpal[l]; rpal[l] = x;
      x = gpal[i]; gpal[i] = gpal[l]; gpal[l] = x;
      x = bpal[i]; bpal[i] = bpal[l]; bpal[l] = x;
    }
  }

  pals = new unsigned char[3*npal];
  for(int i = 0; i < npal; i++) {
    int k = 3*i;
    pals[k] = rpal[i];
    pals[k+1] = gpal[i];
    pals[k+2] = bpal[i];
  }
  npal1 = npal-1;
}

void FitsPalette::SetPalette(const wxString& a)
{
  for(int i = PAL_FIRST+1; i < PAL_LAST; i++)
    if( a == Type_str(i) ) {
      SetPalette(i);
      return;
    }
  wxFAIL_MSG("FitsPalette::SetPalette(const wxString& a): unknown palette.");
}

void FitsPalette::Create_Grey()
{
  npal = 256;
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];

  for(int i = 0; i < npal; i++ )
    rpal[i] = gpal[i] = bpal[i] = i;
}

void FitsPalette::Create_Sepia()
{
  // https://en.wikipedia.org/wiki/Sepia_(color)

  npal = 363; // 256*sqrt(2)
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];

  for(int i = 0; i < npal; i++ ) {
    float L = float(i) / float(npal-1);
    HSL_RGB(30.0,1.0,L,rpal[i],gpal[i],bpal[i]);
  }
}

void FitsPalette::Create_Royal()
{
  // Royal blue
  npal = 363;
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];

  for(int i = 0; i < npal; i++ ) {
    float L = float(i) / float(npal-1);
    HSL_RGB(225.0,1.0,L,rpal[i],gpal[i],bpal[i]);
  }
}

void FitsPalette::Create_VGA()
{
  // VGA pallete (http://en.wikipedia.org/wiki/Web_colors)
  const int vga[16][3] = {
    {0,0,0},       // black #000000
    {0,0,128},     // navy, #000080
    {0,128,0},     // green, #008000
    {128,0,0},     // maroon, #800000
    {128,128,0},   // olive, #808000
    {0,128,128},   // teal, #008080
    {128,0,128},   // purple, #800080
    {128,128,128}, // grey, #808080
    {0,0,255},     // blue, #0000FF
    {0,255,255},   // aqua, #00FFFF
    {255,0,255},   // fuchsia, #FF00FF
    {0,255,0},     // lime, #00FF00
    {255,255,0},   // yelow, #FFFF00
    {255,0,0},     // red, #FF0000
    {192,192,192}, // silver , #C0C0C0
    {255,255,255}};// white, #FFFFFF

  npal = 256;
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];

  int n = 0;
  for(int i = 0; i < 16; i++ ) {
    for(int j = 0; j < 16; j++ ) {
      rpal[n] = vga[i][0];
      gpal[n] = vga[i][1];
      bpal[n] = vga[i][2];
      n++;
    }
  }
  wxASSERT(n==npal);
}


void FitsPalette::Create_AIPS0()
{
  // by ds9
  const int aips0[9][3] = {
    {49,49,49},
    {121,0,155},
    {0,0,200},
    {95,167,235},
    {0,151,0},
    {0,246,0},
    {255,255,0},
    {255,176,0},
    {255,0,0}};

  npal = 252;
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];

  int n = 0;
  for(int i = 0; i < 9; i++ ) {
    for(int j = 0; j < 28; j++ ) {
      rpal[n] = aips0[i][0];
      gpal[n] = aips0[i][1];
      bpal[n] = aips0[i][2];
      n++;
    }
  }
  wxASSERT(n==npal);
}

void FitsPalette::Create_Madness()
{
  npal = 3*256;
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];

  float l = 256;
  for(int i = 0; i < npal; i++) {
    float x = 3.33*(i - 256);
    rpal[i] = int(255.0*cos(3.14*(l+x)/l));
    gpal[i] = int(255.0*sin(3.14*x/l));
    bpal[i] = int(255.0*cos(3.14*x/l));
  }
}

void FitsPalette::Create_Rainbow()
{
  npal = 2*240;
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];

  for(int i = 0; i < npal; i++) {
    float H = 240.0 - i / 2.0;
    HSL_RGB(H,1.0,0.5,rpal[i],gpal[i],bpal[i]);
  }
}

void FitsPalette::Create_Highlight()
{
  npal = 256 + 50 + 241;
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];

  for(int i = 0; i < 256; i++)
    rpal[i] = gpal[i] = bpal[i] = i;

  int n = 256;
  for(int i = 100; i > 50; i-- ) {
    float L = 0.01*i;
    HSL_RGB(240.0,0.7,L,rpal[n],gpal[n],bpal[n]);
    n++;
  }

  for(int i = 240; i >= 0; i--) {
    float H = i;
    HSL_RGB(H,0.99,0.5,rpal[n],gpal[n],bpal[n]);
    n++;
  }

  wxASSERT(n==npal);
}

void FitsPalette::Create_Colour()
{
  npal = 32 + 2*240;
  rpal = new unsigned char[npal];
  gpal = new unsigned char[npal];
  bpal = new unsigned char[npal];

  for(int i = 0; i < 32; i++) {
    float L = i / 62.0;
    HSL_RGB(240.0,0.75,L,rpal[i],gpal[i],bpal[i]);
  }

  int n = 32;
  for(int i = 0; i < 2*240; i++) {
    float H = 240.0 - i / 2.0;
    HSL_RGB(H,1.0,0.5,rpal[n],gpal[n],bpal[n]);
    n++;
  }

  wxASSERT(n==npal);
}


/*

  Functions below implements formula "HSL to RGB alternative" on the
  wiki page: https://en.wikipedia.org/wiki/HSL_and_HSV

*/

void FitsPalette::HSL_RGB(float H, float S, float L,
			  unsigned char& r, unsigned char & g, unsigned char& b)
{
  r = int(255*f(0,H,S,L));
  g = int(255*f(8,H,S,L));
  b = int(255*f(4,H,S,L));
}

float FitsPalette::f(int n, float H, float S, float L) const
{
  float k = fmod((n + H / 30.0),12.0);
  float a = S*min(L,float(1.0-L));
  return L - a*max(-1.0,min(min(k-3.0,9.0-k),1.0));
}
