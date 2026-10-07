/*

  xmunipack - display's headers

  Copyright © 2012-13, 2018, 2021 F.Hroch (hroch@physics.muni.cz)

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
#include "enum.h"
#include "event.h"
#include "mconfig.h"
#include "types.h"
#include <wx/wx.h>
#include <list>


class DisplayShrinkRender: public wxThread
{
  wxEvtHandler *handler;
  FitsImage image;
  FitsTone tone;
  int shrink;
  int thisid;
  ExitCode Entry();

 public:
  DisplayShrinkRender(wxEvtHandler *, const FitsImage&, const FitsTone&,
		      int, int);
  virtual ~DisplayShrinkRender();
};


class DisplayTuneRender: public wxThread
{
  wxEvtHandler *handler;
  FitsImage image;
  FitsTone tone;
  FitsItt itt;
  FitsPalette pal;
  FitsColour colour;
  int thisid;
  ExitCode Entry();

 public:
  DisplayTuneRender(wxEvtHandler *, const FitsImage&, const FitsTone&,
		    const FitsItt&, const FitsPalette&, const FitsColour&, int);
  virtual ~DisplayTuneRender();
};


class MuniDisplayCanvas: public wxWindow
{
public:
  MuniDisplayCanvas(wxWindow *w, MuniConfig *c);
  virtual ~MuniDisplayCanvas();

  void SetHdu(const FitsArray&, const wxImage&);
  void SetInitShrink(int);
  void InvokeRendering();
  void StopRendering();
  FitsTone GetTone() const { return tone; };
  FitsItt GetItt() const { return itt; };
  FitsPalette GetPalette() const { return pal; }
  FitsColour GetColour() const { return color; }
  wxImage GetImage() const;
  wxPoint GetCartesianPosition(const wxPoint&) const;

  void OnTuneFine(MuniTuneEvent&);
  void OnZoom(MuniZoomEvent&);
  void OnValueType(wxCommandEvent&);
  void OnCooType(wxCommandEvent&);
  void OnAstrometry(MuniAstrometryEvent&);
  void OnPhotometry(MuniPhotometryEvent&);
  void OnDraw(MuniDrawEvent&);
  void AddLayer(const MuniLayer&);
  void RemoveLayers(int);

private:

  friend class DisplayShrinkRender;
  friend class DisplayTuneRender;
  wxCriticalSection renderCS;
  wxThread *render;
  int canvasid;


  wxWindow *topwin;
  MuniConfig *config;
  wxEvtHandler *handler;

  FitsArray array;
  wxImage icon, canvas;
  int xoff,yoff,xsubwin,ysubwin,xdrag0,ydrag0,vbX,vbY;
  FitsImage fitsimage, scaled;
  FitsPalette pal;
  FitsTone tone;
  FitsItt itt;
  FitsColour color;
  int value_type, coo_type;
  int scaleicon;
  int shrink, zoom;
  bool zooming, update_zoom, shrinking, tunning, completed;
  bool dragging,invoking,rendering,finished;
  bool astrometry;
  double xcen, ycen, acen, dcen, ascale, aangle;
  wxPoint crosshair;
  std::list<MuniLayer> layers;

  void OnClose(wxCloseEvent&);
  void OnSize(wxSizeEvent&);
  void OnPaint(wxPaintEvent&);
  void OnMouseMotion(wxMouseEvent&);
  void OnMouseEnter(wxMouseEvent&);
  void OnMouseLeave(wxMouseEvent&);
  void OnMouseWheel(wxMouseEvent&);
  void OnClick(wxMouseEvent&);
  void OnLeftUp(wxMouseEvent&);
  void OnRenderFinish(MuniRenderEvent&);
  void OnSubRender(MuniRenderEvent&);
  void Reset();
  //  void OverlayBitmap(const vector<MuniDrawBase *>&);
  //  void OverlayGrid(const vector<wxGraphicsPath>&);
  void OnClipValue(wxCommandEvent&);
  void OnClipCoo(wxCommandEvent&);
  void OnKeyDown(wxKeyEvent&);
  void OnMenu(wxMouseEvent& event);
  void OnIdle(wxIdleEvent&);
  //  void DrawSVG(wxPaintDC&, const wxXmlDocument&);
  void DrawLayer(wxPaintDC&, const MuniLayer&);
  void ConvertCoo(double,double,double&,double&);

  void InitCanvas();
  void UpdateCanvas(int,int);
  void UpdateCanvas();
  void UpdateTune();
  void UpdateScale(int);
  void UpdateScale();
  //  void UpdateScroll();
  void StartRendering();

};
