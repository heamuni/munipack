/*

  xmunipack - zoomer

  Copyright © 2021 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _ZOOMER_H
#define _ZOOMER_H

#include "render.h"
#include "mconfig.h"
#include <wx/wx.h>


class MuniZoomSet
{
  int nelements, index;
  double *table;

  void Init();

public:
  MuniZoomSet();
  MuniZoomSet(double);
  MuniZoomSet(int);
  virtual ~MuniZoomSet();
  MuniZoomSet(const MuniZoomSet&);
  MuniZoomSet& operator=(const MuniZoomSet&);

  bool SetZoom(int);
  bool SetZoom(double);
  double GetZoom();
  int GetIndex() const;

  double SetBestFitZoom(const wxSize&, const wxSize&);
  double SetNormalZoom();
  double Increase();
  double Decrease();

};

class MuniZoomWindow: public wxPanel
{
  const MuniConfig *config;
  bool dragging, init, init_subwin, subwin_focused, zoom_changed, tone_changed,
    render_start;
  int dx, dy, shrink, iwidth, iheight, rid;
  double display_zoom, side_ratio;
  unsigned char *idata;
  wxPoint preview_offset;
  wxPoint subwin_center;
  wxRect subwin;
  std::vector<wxRect> refresh;
  wxColour bgcolour;
  wxSize display_size, image_size, canvas_size, preview_size;
  wxBitmap canvas, overlay;
  MuniDisplayRender *render;
  FitsArray fitsimage;
  FitsTone tone;
  FitsItt itt;
  FitsPalette pal;
  FitsColour colour;

  void OnIdle(wxIdleEvent&);
  void OnPaint(wxPaintEvent&);
  void OnSize(wxSizeEvent&);
  void OnZoom(MuniZoomEvent&);
  void OnMouseMotion(wxMouseEvent&);
  void OnMouseUp(wxMouseEvent&);
  void OnMouseDown(wxMouseEvent&);
  void OnRender(MuniRenderEvent&);
  void OnTune(MuniTuneEvent&);
  wxString GetLabel() const;
  wxSize ScaledSize() const;
  bool IsInitialised() const;
  bool ShowSubwin() const;
  void InitCanvas();
  void RefreshCanvas(int,int,int,int);
  void InitScale();
  void InitSubwin();
  void DrawSubwin(wxGraphicsContext *);
  void DrawZoom(wxGraphicsContext *);

public:
  MuniZoomWindow(wxWindow *, const MuniConfig *c);
  virtual ~MuniZoomWindow();
  void SetDisplayZoom(double);
  void SetArray(const FitsArray&);
  void SetImageSize(const wxSize&);
  void SetDisplaySize(const wxSize&);
  void SetTone(const FitsTone&);
  void SetItt(const FitsItt&);
  void SetPalette(const FitsPalette&);
  void SetColour(const FitsColour&);
  void SetTune(const FitsTone&,const FitsItt&, const FitsPalette&, const FitsColour&);
  void StartRendering();

};

class MuniZoomSlider: public wxPanel
{
  wxSlider *slider;

  void OnScroll(wxScrollEvent&);

public:
  MuniZoomSlider(wxWindow *, wxWindowID);
  void SetZoom(double);
  double GetZoom() const;
  double GetZoom(int) const;

};

#endif
