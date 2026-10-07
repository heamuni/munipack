/*

  xmunipack - find plot

  Copyright © 2019 F.Hroch (hroch@physics.muni.cz)

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

#include "xmunipack.h"
#include "mathplot.h"
#include "plot.h"
#include <wx/wx.h>
#include <vector>
#include <algorithm>
#include <cstdlib>

#define SQR(x) (x)*(x)

const double hwhm = 1.1774; // HWHM = sqrt(2*log(2))

// Gauss
double gauss(double x, double y, double g0, double s, double x0, double y0)
{
  double ss = 2*SQR(s);
  double r2 = SQR(x - x0) + SQR(y - y0);
  return g0*exp(-r2 / ss) / (2.50663 * s);
}

// Random number in range -0.5 .. 0.5
double rnd()
{
  const double rmax = RAND_MAX;
  return rand() / rmax - 0.5;
}

using namespace std;

// -- Slice

class mpGauss: public mpProfile
{
  double g0, s;
  double x0, y0;

public:
  mpGauss(double g, double d, int i0, int j0):
    mpProfile("Template"),g0(g), s(d), x0(i0), y0(j0) {}

  double GetY(double x) {
    return gauss(x,0.0,g0,s,x0,y0);
  }

};


// -- MuniPlotFind  --------------------------------------------

MuniPlotFind::MuniPlotFind(wxWindow *w, const FitsArray& a):
  mpWindow(w,wxID_ANY), array(a), i0(-1), j0(-1), side(1), fwhm(1.0), back(0),
  bestsize(wxSize(300,300)), update(false)
{
  AddLayer(new mpText(L"⌖  Pick a star ...",10,30));
  Bind(wxEVT_IDLE,&MuniPlotFind::OnIdle,this);

  // initial illustration
  vector<double> xslice, yslice, xpoints, ypoints;
  double dq = rnd();
  for(int j = -10; j < 10; j++) {
    double g = gauss(0.0,double(j),100.0,1.0,0.0,5.0);
    double qx = g + 10*rnd();
    double qy = g + 10*rnd();
    double x = j - 0.5;
    double y = x + dq;
    xpoints.push_back(x);
    xpoints.push_back(x+1);
    yslice.push_back(qy);
    yslice.push_back(qy);
    ypoints.push_back(y);
    ypoints.push_back(y+1);
    xslice.push_back(qx);
    xslice.push_back(qx);
  }

  mpFXYVector *xprofile = new mpFXYVector();
  xprofile->ShowName(false);
  wxPen xpen(*wxRED,1);
  xprofile->SetPen(xpen);
  xprofile->SetData(xpoints,xslice);
  xprofile->SetContinuity(true);
  AddLayer(xprofile);

  mpFXYVector *yprofile = new mpFXYVector();
  yprofile->ShowName(false);
  wxPen ypen(*wxBLUE,1);
  yprofile->SetPen(ypen);
  yprofile->SetData(ypoints,yslice);
  yprofile->SetContinuity(true);
  AddLayer(yprofile);
  Fit();


  /* m27_01R.fits
  i0 = 576;
  j0 = 372;
  SetPoint(576,372);
  SetFwhm(2.0);
  */
  //  Update();
}

void MuniPlotFind::OnIdle(wxIdleEvent&)
{
  if( update ) {
    update = false;
    Refresh();
  }
}


void MuniPlotFind::SetPoint(int i, int j)
{
  i0 = i;
  j0 = j;
  int d = wxMax(side / 2,1);

  double qmax = array.Pixel(i0,j0);
  for(int k = i - d; k < i + d; k++)
    for(int l = j - d; l < j + d; l++) {
      if( 0 <= k && k < array.GetWidth() && 0 <= l && l < array.GetHeight() ) {
	if( array.Pixel(k,l) > qmax ) {
	  qmax = array.Pixel(k,l);
	  i0 = k;
	  j0 = l;
	}
      }
    }
}

void MuniPlotFind::SetFwhm(double f)
{
  if( fwhm > 0.01 ) {
    fwhm = f;
    side = int(3*fwhm + 0.5);
  }
  else {
    fwhm = 0.01;
    side = 1;
  }
}



void MuniPlotFind::Update()
{
  update = true;
}

void MuniPlotFind::Refresh()
{
  if( i0 < 0 )
    return;

  DelAllLayers(true,false);

  // sub-window
  int imin = wxMax(i0-side,0);
  int imax = wxMin(i0+side,array.GetWidth()-1);
  int jmin = wxMax(j0-side,0);
  int jmax = wxMin(j0+side,array.GetHeight()-1);
  //  wxLogDebug("bouding box: %d %d %d %d", imin, imax,jmin,jmax);
  FitsGeometry garray(array);
  FitsArray sub(garray.GetSubArray(imin,jmin,imax - imin, jmax - jmin));

  // background
  vector<double> b;
  b.push_back(array.Pixel(imin,jmin));
  b.push_back(array.Pixel(imax,jmin));
  b.push_back(array.Pixel(imin,jmax));
  b.push_back(array.Pixel(imax,jmax));
  sort(b.begin(),b.end());
  back = (b[1] + b[2]) / 2;
  //  wxLogDebug("%f %f %f %f back= %f",b[0],b[1],b[2],b[3],back);

  // X slice
  FitsArray xarray(garray.GetSubArray(imin,j0,imax - imin,1));
  vector<double> xslice;
  vector<double> points;
  for(int i = 0; i < xarray.GetWidth(); i++) {
    double q = xarray.Pixel(i,0) - back;
    points.push_back(i-0.5);
    xslice.push_back(q);
    points.push_back(i+0.5);
    xslice.push_back(q);
    //    wxLogDebug("x: %f %f",points[i],xslice[i]);
  }
  mpFXYVector *xprofile = new mpFXYVector("Horizontal");
  xprofile->ShowName(false);
  wxPen xpen(*wxRED,1);
  xprofile->SetPen(xpen);
  xprofile->SetData(points,xslice);
  xprofile->SetContinuity(true);
  AddLayer(xprofile);

  // Y slice
  FitsArray yarray(garray.GetSubArray(i0,jmin,1,jmax - jmin));
  vector<double> yslice;
  points.clear();
  for(int j = 0; j < yarray.GetHeight(); j++) {
    double q = yarray.Pixel(0,j) - back;
    points.push_back(j-0.4);
    yslice.push_back(q);
    points.push_back(j+0.6);
    yslice.push_back(q);
    //    wxLogDebug("y: %f %f",points[j],yslice[j]);
  }
  mpFXYVector *yprofile = new mpFXYVector("Vertical");
  yprofile->ShowName(false);
  wxPen ypen(*wxBLUE,1);
  yprofile->SetPen(ypen);
  yprofile->SetData(points,yslice);
  yprofile->SetContinuity(true);
  AddLayer(yprofile);

  // Gaussian
  double s = fwhm / (2*hwhm);
  double g0 = gnorm(imin,imax,jmin,jmax,s);
  double x0 = ((i0 - imin) + (j0 - jmin)) / 2;
  //  wxLogDebug("g0=%f s=%f",g0,s);
  mpGauss *f = new mpGauss(g0,s,x0,0.0);
  f->ShowName(false);
  wxPen gpen(*wxGREY_PEN);
  gpen.SetWidth(2);
  f->SetPen(gpen);
  AddLayer(f);

  // FWHM
  vector<double> marker;
  points.clear();
  points.push_back(x0 - s*hwhm);
  marker.push_back(gauss(s*hwhm,0.0,g0,s,0.0,0.0));
  points.push_back(x0 + s*hwhm);
  marker.push_back(marker[0]);
  mpFXYVector *sline = new mpFXYVector("FWHM");
  sline->ShowName(false);
  wxPen spen(*wxBLACK_PEN);
  spen.SetWidth(5);
  sline->SetPen(spen);
  sline->SetData(points,marker);
  sline->SetContinuity(true);
  AddLayer(sline);

  mpInfoLegend *legend = new mpInfoLegend();
  wxFont fn(*wxSMALL_FONT);
  legend->SetFont(fn);
  legend->SetPen(wxNullPen);
  AddLayer(legend);

  Fit();
}

double MuniPlotFind::gnorm(int imin, int imax, int jmin, int jmax,
			   double s) const
{
  double s1 = 0;
  double s2 = 0;
  double x0 = i0;
  double y0 = j0;

  for(int i = imin; i < imax; i++) {
    double x = i;
    for(int j = jmin; j < jmax; j++) {
      double y = j;
      double gij = gauss(x,y,1.0,s,x0,y0);
      s1 = s1 + gij*(array.Pixel(i,j) - back);
      s2 = s2 + SQR(gij);
    }
  }
  return s1 / s2;
}
