/*

  xmunipack - colours

  Copyright © 2009-2014, 2019-22 F.Hroch (hroch@physics.muni.cz)

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


  Color images:

  rawtran -o IMG_5807.fits -X "-q 3 -n 500" IMG_5807.CR2

 * important: -q selects adequate interpolation method,
              -n 500 selects threshold for wavelets


*/


#include "fits.h"
#include <wx/wx.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>
#include <wx/tokenzr.h>

/* CIE XYZ, D65 white point */
#define D65_Xn 95.047
#define D65_Yn 100.000
#define D65_Zn 108.883
#define D65_WHITEPOINT_X 0.31271
#define D65_WHITEPOINT_Y 0.32902

/* white-point for Luv, D65 */
#define D65_uw 0.19783
#define D65_vw 0.46832


// --- FitsColour


FitsColour::FitsColour():
  ispace(COLOURSPACE_XYZ),saturation(1.0),
  nitevision(false),meso_level(0.0),meso_width(1.0),
  ncolours(3),nbands(3),trafo(0),level(0),weight(0),
  Sn(Scotopic(D65_Xn,D65_Yn,D65_Zn)/100), uw(D65_uw),vw(D65_vw)
{
  SetTrans(3,3);
}

FitsColour::FitsColour(const wxString& cdatafile, const FitsArray& array):
  ispace(COLOURSPACE_XYZ),saturation(1.0),
  nitevision(false),meso_level(0.0),meso_width(1.0),
  ncolours(0),nbands(0),trafo(0),level(0),weight(0),
  Sn(Scotopic(D65_Xn,D65_Yn,D65_Zn)/100), uw(D65_uw),vw(D65_vw)
{
  Init(cdatafile,array);
}

FitsColour::FitsColour(const FitsColour& c):
  ispace(c.ispace),saturation(c.saturation),
  nitevision(c.nitevision),meso_level(c.meso_level),meso_width(c.meso_width),
  ncolours(c.ncolours),nbands(c.nbands),trafo(0),level(0),weight(0),
  Sn(c.Sn),uw(c.uw),vw(c.vw)
{
  int n = ncolours*nbands;
  trafo = new float[n];
  std::copy(c.trafo,c.trafo+n,trafo);
  level = new float[nbands];
  weight = new float[nbands];
  std::copy(c.level,c.level+nbands,level);
  std::copy(c.weight,c.weight+nbands,weight);
}

FitsColour& FitsColour::FitsColour::operator=(const FitsColour& other)
{
  if( this != &other ) {
    ispace = other.ispace;
    saturation = other.saturation;
    nitevision = other.nitevision;
    meso_level = other.meso_level;
    meso_width = other.meso_width;
    ncolours = other.ncolours;
    nbands = other.nbands;
    Sn = other.Sn;
    uw = other.uw;
    vw = other.vw;

    delete[] trafo;
    delete[] level;
    delete[] weight;

    int n = ncolours*nbands;
    trafo = new float[n];
    level = new float[nbands];
    weight = new float[nbands];
    std::copy(other.trafo,other.trafo+n,trafo);
    std::copy(other.level,other.level+nbands,level);
    std::copy(other.weight,other.weight+nbands,weight);
  }
  return *this;
}


FitsColour::~FitsColour()
{
  delete[] trafo;
  delete[] level;
  delete[] weight;
}

void FitsColour::Init(const wxString& cdatafile, const FitsArray& array)
{
  wxString cs = array.GetKey(FITS_KEY_CSPACE);
  if( cs.IsEmpty() )
    cs = "RGB";
  // RGB is supposed due to compatibility with other software

  SetTrans(cs,cdatafile);
}

float FitsColour::GetWeight(int n) const
{
  wxASSERT(weight && 0 <= n && n < nbands);
  return weight[n];
}

float FitsColour::GetLevel(int n) const
{
  wxASSERT(level && 0 <= n && n < nbands);
  return level[n];
}


void FitsColour::SetWeight(int n, float w)
{
  wxASSERT(weight && 0 <= n && n < nbands);
  weight[n] = w;
}

void FitsColour::SetLevel(int n, float x)
{
  wxASSERT(level && 0 <= n && n < nbands);
  level[n] = x;
}

void FitsColour::SetTrans(int n, int m, float x)
{
  wxASSERT(trafo && 0 <= n && n < ncolours && 0 <= m && m < nbands);
  trafo[m*nbands + n] = x;
}


void FitsColour::SetTrans(int n, int m)
{
  delete[] weight;
  delete[] level;
  delete[] trafo;

  ncolours = n;
  nbands = m;
  weight = new float[n];
  level = new float[n];
  trafo = new float[n*m];
  for(int i = 0; i < n; i++) {
    weight[i] = 1.0;
    level[i] = 0.0;
    for(int j = 0; j < m; j++)
      trafo[j*n + i] = i == j ? 1.0 : 0.0;
  }
}


void FitsColour::SetTrans(const wxString& cs)
{
  cspace = cs;

  if( cspace.Find("XYZ") != wxNOT_FOUND ) {
    size_t n = 3;
    SetTrans(n,n);
    for(size_t i = 0; i < n; i++)
      for(size_t j = 0; j < n; j++)
	SetTrans(i,j,0.0);

    for(size_t i = 0; i < n; i++) {
      SetTrans(i,i,1.0);
      SetLevel(i,0.0);
      SetWeight(i,1.0);
    }
  }
}

void FitsColour::SetTrans(const wxString& cs, const wxString& filename)
{
  cspace = cs;

  wxFileInputStream input(filename);
  wxTextInputStream text(input," ,\t");

  while(input.IsOk() && ! input.Eof()) {
    wxString line,ilabel,olabel;
    int n,m;

    line = text.ReadLine();
    line.Trim();
    if( line.IsEmpty() ) continue;

    wxArrayString a;
    wxStringTokenizer t(line,"'");
    int i = 0;
    while ( t.HasMoreTokens() ) {
      wxString x = t.GetNextToken();
      x.Trim();
      if( ! x.IsEmpty() ) {
	a.Add(x);
      }
      i++;
    }

    if( a.GetCount() == 2 ) {
      ilabel = a[0];
      olabel = a[1];
    }

    if( input.Eof() ) break;
    n = text.Read32();
    if( input.Eof() ) break;
    m = text.Read32();
    if( input.Eof() ) break;

    float *cmatrix = new float[n*m];
    for(int i = 0; i < n*m; i++) {
      if( input.Eof() ) break;
      cmatrix[i]= text.ReadDouble();
    }

    //   wxLogDebug(ilabel + " > " +olabel+ " , " +cspace+ " , " + Type_str(COLORSPACE_XYZ));

    if( ilabel == cspace && olabel == Type_str(COLOURSPACE_XYZ) ) {

      SetTrans(n,m);
      return;
    }
    delete[] cmatrix;
  }


  // proper colorspace data not found
  SetTrans(3,3);
}


int FitsColour::GetColours() const
{
  return ncolours;
}

int FitsColour::GetBands() const
{
  return nbands;
}

float FitsColour::GetTrans(int n, int m) const
{
  wxASSERT(trafo && 0 <= n && n < ncolours && 0 <= m && m < nbands);
  return *(trafo+m*nbands + n);
}

wxString FitsColour::GetColourspace() const
{
  return cspace;
}


float FitsColour::GetWhitePointX() const
{
  float s = 6*uw - 16*vw + 12;
  return 9*uw/s;
}

float FitsColour::GetWhitePointY() const
{
  float s = 6*uw - 16*vw + 12;
  return 4*vw/s;
}


void FitsColour::SetWhitePoint(float x, float y)
{
  wxASSERT(0 < x && x < 1 && 0 < y && y < 1);

  float s = -2*x + 12*y + 3;
  uw = 4*x/s;
  vw = 9*y/s;
}


void FitsColour::SetSaturation(float x)
{
  saturation = x;
}

void FitsColour::SetMesoLevel(float x)
{
  meso_level = x;
}

void FitsColour::SetMesoWidth(float x)
{
  meso_width = x;
}

void FitsColour::SetNiteVision(bool t)
{
  nitevision = t;
}

wxString FitsColour::Type_str(int n)
{
  switch(n){
  case COLOURSPACE_XYZ:      return "XYZ";
  default:        return wxEmptyString;
  }
}


wxArrayString FitsColour::Type_str()
{
  wxArrayString a;

  for(int i = COLOURSPACE_XYZ+1; i < COLOURSPACE_LAST; i++)
    a.Add(Type_str(i));

  return a;
}


float FitsColour::Scotopic(float X, float Y, float Z) const
{
  return 0.36169f*Z + 1.18214f*Y - 0.80498f*X;
}


void FitsColour::Instr_XYZ(long npix, size_t nband, const float **d,
			   float *Z, float *Y, float *X)
{
  wxASSERT(nband == (size_t) nbands);

  long nbytes = npix*sizeof(float);
  memcpy(X,d[0],nbytes);
  memcpy(Y,d[1],nbytes);
  memcpy(Z,d[2],nbytes);
  return;

  if( GetColourspace().Find("XYZ") != wxNOT_FOUND ) {
    wxASSERT(ncolours == nbands);

    long nbytes = npix*sizeof(float);
    memcpy(X,d[0],nbytes);
    memcpy(Y,d[1],nbytes);
    memcpy(Z,d[2],nbytes);
  }
  else {
    wxASSERT(weight && level && trafo);

    // allocate all on heap?

    float cb[ncolours][nbands];
    for(int i = 0; i < ncolours; i++)
      for(int j = 0; j < nbands; j++) {
	cb[i][j] = GetTrans(i,j);
      }

    float b[nbands];
    float c[ncolours];
    for(long i = 0; i < npix; i++) {

      for(int j = 0; j < nbands; j++) {
	b[j] = weight[j]*(d[j][i] - level[j]);
	if( b[j] < 0.0 ) b[j] = 0.0;
      }

      for(int l = 0; l < ncolours; l++) {
	float s = 0.0;
	for(int j = 0; j < nbands; j++)
	  s = s + cb[l][j]*b[j];
	c[l] = s;
      }

      X[i] = c[0];
      Y[i] = c[1];
      Z[i] = c[2];
    }
  }

}


void FitsColour::NiteVision(float X, float Z, float& x, float& y, float& Y) const
{
  const float xw = D65_WHITEPOINT_X;
  const float yw = D65_WHITEPOINT_Y;
  const float xyw = xw / yw;

  // simplified transformation for photo-,meso-, and scotopic
  float S = Scotopic(X,Y,Z) / Sn;
  float w = (1.0f + tanh((Y - meso_level)/meso_width)) / 2.0f;

  float Yw = Y*w;
  float Sw = S*(1 - w);
  float D = Yw/y + Sw/yw;

  Y = Yw + Sw;
  x = ((x/y)*Yw  + xyw*Sw) / D;
  y = Y / D;
}
