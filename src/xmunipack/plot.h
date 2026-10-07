
/*

  XMunipack - plotting


  Copyright © 2012-14, 2019-21  F.Hroch (hroch@physics.muni.cz)

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

#ifndef _XMUNIPACK_PLOT_H_
#define _XMUNIPACK_PLOT_H_

#include "event.h"
#include "mathplot.h"
#include "fits.h"
#include <plplot.h>
#include <wx/wx.h>


/*
class MuniPlotTable
{
public:
  MuniPlotTable(const std::vector<wxRealPoint>&,const wxColour& =*wxBLACK);

  std::vector<wxRealPoint> points;
  wxColour colour;

};
*/

// class MuniPlot: public mpWindow
// 		//class MuniPlot: public wxPLplotwindow
// {
// public:
//   MuniPlot(wxWindow *);
//   void AddData(const MuniPlotTable&);

// private:

//   std::vector<MuniPlotTable> tables;
//   double xmin,xmax,ymin,ymax;

//   //  void OnSize(wxSizeEvent&);
//   void Draw();

// };


/*

  A general wrapper for plplot with output to a memory with Cairo driver.
  The user should implement MakeGraph() to provide a plotting function.

 */
class PlPlotCairo: public wxPanel
{
  wxImage graph;
  bool ready;

  void InitGraph();
  void Plotter(bool=false);
  void OnPaint(wxPaintEvent&);

public:
  PlPlotCairo(wxWindow *, const wxSize&);
  void Update();
  virtual void MakeGraph() = 0;

};


class MuniPlotHistoBase: public PlPlotCairo
{

protected:
  FitsHisto hist;
  PLFLT xmin,xmax,ymin,ymax;
  void StrokeHistoProfile(int, int, PLINT&, PLFLT*, PLFLT*);

public:
  MuniPlotHistoBase(wxWindow *, const wxSize&);
  void SetHisto(const FitsHisto&);
};


class MuniPlotHisto: public MuniPlotHistoBase
{
  float black,sense;

  void MakeGraph();

public:
  MuniPlotHisto(wxWindow *, const wxSize&, double, double);
  void SetScale(double,double);
};


class MuniPlotNite: public MuniPlotHistoBase
{
  float level, width;
  bool enabled;

  void MakeGraph();
  void MakePalette(PLFLT, PLINT, PLINT *, PLINT *, PLINT *);

public:
  MuniPlotNite(wxWindow *, const wxSize&, float, float, bool);
  void SetMeso(float, float);
  void SetEnabled(bool);

};


class MuniPlotItt: public PlPlotCairo
{
  FitsItt itt;

  void MakeGraph();

public:
  MuniPlotItt(wxWindow *, const wxSize&, const FitsItt&);
  void SetItt(int);
};


class MuniPlotPalette: public PlPlotCairo
{
  FitsPalette pal;
  FitsItt itt;
  float black, sense;

  void MakeGraph();

public:
  MuniPlotPalette(wxWindow *, const wxSize&, int, bool, int, float, float);
  void SetPalette(int);
  void SetInverse(bool);
  void SetItt(int);
  void SetScale(float,float);
};



class MuniPlotSaturation: public PlPlotCairo
{
  float saturation;
  float white_point[2];

  void MakeGraph();
  void MakeTable(PLFLT, PLFLT, PLINT, PLINT *, PLINT *, PLINT *);
  void plfbox(PLFLT,PLFLT,PLFLT,PLFLT);
  void plftri(PLFLT,PLFLT,PLFLT,PLFLT);

public:
  MuniPlotSaturation(wxWindow *, const wxSize&, float, float, float);
  void SetSaturation(float);
  void SetWhitePoint(float,float);

};




/*
class MuniPlotUV: public wxPLplotwindow
{
public:
  MuniPlotUV(wxWindow *);
  virtual ~MuniPlotUV();
  wxSize DoGetBestSize() const;
  void DrawTri(const std::vector<double>&, const std::vector<double>&);
  void Clear();

private:

  int nuv;
  PLFLT *u,*v;
  void Draw();

};
*/

class MuniPlotFind: public mpWindow
{
  FitsArray array;
  int i0, j0, side;
  double fwhm, back;
  wxSize bestsize;
  bool update;

  double gnorm(int,int,int,int, double) const;

  void OnIdle(wxIdleEvent&);
  void Refresh();

 public:
  MuniPlotFind(wxWindow *, const FitsArray&);
  wxSize DoGetBestSize() const { return bestsize; }
  void SetPoint(int,int);
  void SetFwhm(double);
  void Update();
};

#endif
