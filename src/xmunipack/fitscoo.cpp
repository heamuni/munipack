/*
   xmunipack - coordinates

   Copyright © 2013 - 2026 F.Hroch (hroch@physics.muni.cz)

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
#include <wcslib.h>
#include <wx/wx.h>
#include <math.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

#define RAD 57.295779513082322865


// -----------------------------------------------------------------

FitsProjection::FitsProjection():
  acen(0.0),dcen(0.0),xcen(0.0),ycen(0.0),scale(1.0),angle(0.0),reflex(1.0) {}

FitsProjection::FitsProjection(const wxString& t, double a,double d,double x,
			       double y,double s,double r,double z):
  type(t.Upper()),acen(a/RAD),dcen(d/RAD),xcen(x),ycen(y),scale(s*RAD),
  angle(r/RAD), reflex(z) {}

FitsProjection::FitsProjection(const wxString& t, double a,double d, double x,
			       double y, double cd11, double cd12, double cd21,
			       double cd22):
  type(t.Upper()),acen(a/RAD),dcen(d/RAD),xcen(x),ycen(y),scale(1.0/RAD),
  angle(0.0), reflex(1.0)
{
  reflex = -cd11*cd22 > 0 ? 1.0 : -1.0;

  if( fabs(cd12) > 0 || fabs(cd22) > 0 )
    angle = atan2(reflex*cd12,cd22);
  else
    type = "";

  double q = fabs(cd11*cd22) + fabs(cd12*cd21);
  if( q > 0 )
    scale = RAD/sqrt(q);
  else
    type = "";

}

bool FitsProjection::IsOk() const { return type != "";  }

double FitsProjection::GetScale() const { return scale; }

void FitsProjection::ad2xy(double a,double d,double& x,double& y) const
{
  if( ! IsOk() ) {
    x = 0;
    y = 0;
    return;
  }

  double aa = cos(angle)*scale;
  double bb = sin(angle)*scale;

  double u,v;

  gnomon(a/RAD,d/RAD,u,v);
  x = xcen +  (aa*u + bb*v)*reflex;
  y = ycen + (-bb*u + aa*v);
}

void FitsProjection::xy2ad(double x,double y,double& a,double& d) const
{
  if( ! IsOk() ) {
    a = 0;
    d = 0;
    return;
  }

  double aa = cos(angle)/scale;
  double bb = sin(angle)/scale;

  double u,v,xx,yy;

  u = x - xcen;
  v = y - ycen;
  xx = (aa*u - bb*v)*reflex;
  yy =  bb*u + aa*v;

  ignomon(xx,yy,a,d);
  a = RAD*a;
  d = RAD*d;
}


void FitsProjection::gnomon(double a,double d,double& x,double& y) const
{
  double c,p,q,r,v,w,s;

  c = cos(d);
  p = sin(d);
  q = c*sin(a - acen);
  r = c*cos(a - acen);
  v = sin(dcen);
  w = cos(dcen);
  s = p*v + r*w;

  x = -q/s;
  y = (w*p - v*r)/s;
}

void FitsProjection::ignomon(double x,double y,double& a,double& d) const
{
  double p,q,r,v,w,t;

  v = sin(dcen);
  w = cos(dcen);
  t = sqrt(1.0 + x*x + y*y);
  p = (v + w*y)/t;
  q = -x/t;
  r = (w - v*y)/t;

  d = asin(p);
  a = atan2(q,r) + acen;
}

double FitsProjection::xy2dist(double x1,double y1,double x2,double y2) const
{
  // https://en.wikipedia.org/wiki/Haversine_formula
  double r1,d1,r2,d2;
  xy2ad(x1,y1,r1,d1);
  xy2ad(x2,y2,r2,d2);

  double sind = sin((d2 - d1)/2/RAD);
  double sinr = sin((r2 - r1)/2/RAD);
  double h = sind*sind + cos(d1/RAD)*cos(d2/RAD)*sinr*sinr;
  return h > 0 ? 2*RAD*asin(sqrt(h)) : 0;
}

// -----------------------------------------------------------------


FitsCoo::FitsCoo(): FitsArray(),wcs(0),nwcs(0),type(COO_PIXEL),haswcs(false),digits(4) {}

FitsCoo::FitsCoo(const FitsArray& a):
  FitsArray(a),wcs(0),nwcs(0),type(COO_PIXEL),haswcs(false),digits(4)
{
  // decode WCS
  char *blockhead;
  int nkeyrec, nreject;
  int status = 0;
  status = hdr2str(1,NULL,0,&blockhead,&nkeyrec,&status);

  if ((status = wcspih(blockhead, nkeyrec, WCSHDR_all, 0, &nreject, &nwcs,&wcs)))
    fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);

  wxASSERT(nwcs == 1);
  // if more coordinate representations are included, it will crash else some memory leak will appear

 if ((status = wcsset(wcs)))
    fprintf(stderr, "wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);

  free(blockhead);

  haswcs = false;
  if( status == 0 && wcs ) {
    wxString ctype1 = wcs->ctype[0];
    haswcs = (wcs->naxis == 2 || wcs->naxis == 3) && !ctype1.IsEmpty();
  }

  if( haswcs ) {
    double *pc = wcs->pc;
    double cd11 = *pc;
    double cd12 = *(pc+1);
    double cd21 = *(pc+2);
    double cd22 = *(pc+3);

    double q = fabs(cd11*cd22) + fabs(cd12*cd21);
    if( q > 0 ) {
      double scale = 1 / sqrt(q);
      // significant digits for the print format are determined by the scale
      digits = wxMax(int(log10(scale) + 1),0);
    }
  }

  wxLogDebug("FitsCoo WCS decode status: %s",haswcs ? "T" : "F");
}

FitsCoo::FitsCoo(const FitsCoo& a)
{
  wxFAIL_MSG("*** FitsArray: WE ARE REALLY NEED COPY CONSTRUCTOR");
}

FitsCoo::~FitsCoo()
{
  if( wcs ) {
    int status;
    status = wcsvfree(&nwcs, &wcs);
    if( status )
      fprintf(stderr, "wcsvfree ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
    free(wcs);
  }
}

FitsCoo& FitsCoo::operator = (const FitsCoo& other)
{
  //  wxLogDebug("*** FitsCoo: WE ARE REALLY NEED ASSIGNMENT CONSTRUCTOR");

  if( this != & other ) {

    header = other.header;
    htype = other.htype;
    type = other.type;
    flavour = other.flavour;
    modified = other.modified;
    Ref(other);

    int status;
    wcs = (wcsprm *) malloc(sizeof(struct wcsprm));
    wcs->flag = -1;
    status = wcscopy(1,other.wcs,wcs);
    if( status )
      fprintf(stderr, "wcscopy ERROR %d: %s.\n", status, wcs_errmsg[status]);

    haswcs = status == 0 && other.haswcs;

    nwcs = other.nwcs;
    digits = other.digits;
  }
  return *this;
}

bool FitsCoo::HasWCS() const
{
  return haswcs;
}

void FitsCoo::SetType(int t)
{
  type = static_cast<coords_type>(t);
}

void FitsCoo::SetType(const wxString& a)
{
  for(int i = COO_FIRST+1; i < COO_LAST; i++)
    if( a == Label_str(i) ) {
      SetType(i);
      return;
    }
}

coords_type FitsCoo::GetType() const
{
  return type;
}

void FitsCoo::GetEq(int xcen, int ycen, double& alpha, double& delta) const
{
  double imgcrd[3], phi, pixcrd[3], theta, world[3];
  int stat[1];

  alpha = 0;
  delta = 0;

  if( haswcs ) {
    pixcrd[0] = xcen;
    pixcrd[1] = ycen;
    pixcrd[2] = 0;
    int status = wcsp2s(wcs, 1, wcs->naxis, pixcrd, imgcrd, &phi, &theta, world, stat);
    if( status == 0 && stat[0] == 0 ) {
      alpha = world[0];
      delta = world[1];
    }
    else
      wxLogDebug("wcsp2s ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
  }
}


void FitsCoo::RaSix(double alpha, int& h, int& min, double& sec) const
{
  double x = alpha/15.0;
  h = int(x);
  min = int(60.0*(x - h));
  sec = 3600.0*(x - (h + min/60.0));
}

void FitsCoo::DecSix(double delta, char& sign, int& d, int& m, double& s) const
{
  sign = delta > 0.0 ? '+' : '-';

  double y = fabs(delta);
  d = int(y);
  m = int(60.0*(y - d));
  s = 3600.0*(y - (d + m/60.0));
}


void FitsCoo::GetPix(int i, int j, wxString& a, wxString& b) const
{
  if( 0 <= i && i < Width() && 0 <= j && j < Height() ) {
    a.Printf("%d",i+1); // C-style indexing should be considered very funny!
    b.Printf("%d",j+1);
  }
  else {
    a.Clear();
    b.Clear();
  }
}

void FitsCoo::GetCoo(int i, int j, wxString& a, wxString& b) const
{
  if( haswcs && 0 <= i && i < Width() && 0 <= j && j < Height() ) {

    int h,min,d,m;
    double x,y, sec,s;
    char sg;

    GetEq(i,j,x,y);

    if( type == COO_EQDEG ) {
      a.Printf(L"%.*f°",digits,x);
      b.Printf(L"%.*f°",digits,y);
    }

    else if( type == COO_EQSIX ) {
      int n = wxMax(digits-4,0);
      RaSix(x,h,min,sec);
      DecSix(y,sg,d,m,s);
      a.Printf(L"%02dh%02dm%04.*fs",h,min,n+1,sec);
      b.Printf(L"%c%02d°%02d'%02.*f\"",sg,d,m,n,s);
    }
  }
  else {
    a.Clear();
    b.Clear();
  }
}

wxString FitsCoo::Label_str(int n)
{
  switch(n){
  case COO_PIXEL: return  "pixels";
  case COO_EQDEG: return L"Equatorial (°)";
  case COO_EQSIX: return  "Equatorial (h:m:s)";
  default:        return  "";
  }
}

wxArrayString FitsCoo::Label_str()
{
  wxArrayString a;

  for(int i = COO_FIRST+1; i < COO_LAST; i++)
    a.Add(Label_str(i));

  return a;
}
