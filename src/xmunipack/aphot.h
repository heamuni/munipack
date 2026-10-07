/*

  xmunipack - aphot headers

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


#include "mtool.h"
#include "mprocess.h"
#include "types.h"
#include "event.h"
#include <wx/wx.h>

class MuniApertures
{
  int naper;
  double *aper;
  bool spiral;

public:
  MuniApertures(): naper(0), aper(0), spiral(1) {}
  MuniApertures(int, int, bool);
  MuniApertures(const MuniApertures&);
  MuniApertures& operator=(const MuniApertures&);
  virtual ~MuniApertures();
  int GetNaper() const { return naper; }
  double GetAper(int) const;
  bool FindAper(double, int, int *) const;

};

class MuniAphotDisplay: public wxWindow
{
  int naper, saper, zoom, i0, j0, width, height;
  double rmin, rmax, fwhm;
  bool spiral;
  wxImage image;
  wxBitmap canvas, subimage;
  MuniApertures aper;
  void OnSize(wxSizeEvent&);
  void OnPaint(wxPaintEvent&);
  void Draw(wxBitmap&);
  void OnClick(MuniClickEvent&);
  void OnKlick(wxMouseEvent&);
  void OnZoom(wxCommandEvent&);
  int EstimSaper(double) const;
  void Redraw();
  void Rebase();

public:
  MuniAphotDisplay(wxWindow *, const wxImage&, int,int,int,int,int,bool);
  wxSize DoGetBestSize() const { return wxSize(400,400); }
  void SetNaper(int n) { naper = n; saper = wxMin(saper,n); Redraw(); }
  void SetSpiral(bool s) { spiral = s; Redraw(); }
  void SetRingMin(double r) { rmin = r; Redraw(); }
  void SetRingMax(double r) { rmax = r; Redraw(); }
  void SetPosition(int,int);
  void SetPosition(double x, double y) { i0 = int(x+0.5); j0 = int(y+0.5); }
  void SetFwhm(double f) { fwhm = f; saper = EstimSaper(f); Redraw(); }
  int GetZoom() const { return zoom; }
  int GetSaper() const { return saper; }
  double GetRaper() const { return aper.GetAper(saper); }
  MuniApertures GetAper() const { return aper; }

};


class MuniAphot: public MuniTool
{
  wxString fitsname, tmpfits;
  wxImage image;
  MuniProcess *aproc;
  MuniAphotDisplay *canvas;
  int naper, saper;
  double zoom, ring_min, ring_max;
  bool snap, spiral;

  void OnIdle(wxIdleEvent&);
  void OnTimer(wxTimerEvent&);
  void OnFinish(wxProcessEvent&);
  void OnInput(MuniProcess *);
  void OnOutput(const wxArrayString&);
  void OnSave(wxCommandEvent&);
  void OnClick(MuniClickEvent&);
  void OnSnap(wxCommandEvent&);
  void OnZoom(wxCommandEvent&);
  void OnSpiral(wxCommandEvent&);
  void OnNaper(wxSpinEvent&);
  void OnRingMin(wxSpinEvent&);
  void OnRingMax(wxSpinEvent&);
  void OnUpdateRing(wxUpdateUIEvent&);
  void ByHand(int,int);
  void OnHandFinish(wxProcessEvent&);
  bool StarParser(const wxString&, double *, double *, double *) const;
  bool HandParser(const wxString&, double *, double *, double *,
		  double *, double *, double *) const;
  virtual void CleanDraw() const;
  virtual void DrawStars(const wxArrayString&) const;

public:
  MuniAphot(wxWindow *, MuniConfig *, const FitsFile&, const wxImage&);
  virtual ~MuniAphot();

  void SetPoint(int,int) const;

};
