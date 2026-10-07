/*

  Display for Aperture photometry tool

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

#include "aphot.h"
#include <wx/wx.h>
#include <wx/graphics.h>
#include <cstdlib>
#include <cmath>
#include <cfloat>

using namespace std;

MuniAphotDisplay::MuniAphotDisplay(wxWindow *w, const wxImage& i,
				   int n, int m, int z, int rin, int rout,
				   bool s):
  wxWindow(w,wxID_ANY), naper(n), saper(m), zoom(z), rmin(rin), rmax(rout),
  fwhm(-1), spiral(s), image(i)
{
  i0 = image.GetWidth() / 2;
  j0 = image.GetHeight() / 2;
  width = -1;
  height = -1;
  aper = MuniApertures(naper,rmin,spiral);
  if( ! (0 < saper && saper < naper) )
    saper = naper / 2 + 1;
  Bind(wxEVT_SIZE,&MuniAphotDisplay::OnSize,this);
  Bind(wxEVT_PAINT,&MuniAphotDisplay::OnPaint,this);
  Bind(EVT_CLICK,&MuniAphotDisplay::OnClick,this);
  Bind(wxEVT_LEFT_DOWN,&MuniAphotDisplay::OnKlick,this);
  Bind(wxEVT_SLIDER,&MuniAphotDisplay::OnZoom,this);
}

int MuniAphotDisplay::EstimSaper(double fwhm) const
{
  for(int i = 1; i <= aper.GetNaper(); i++)
    if( aper.GetAper(i) > fwhm )
      return i;
  return aper.GetNaper() / 2 + 1;
}

void MuniAphotDisplay::Redraw()
{
  aper = MuniApertures(naper,rmin,spiral);
  canvas = subimage;
  Draw(canvas);
  Refresh();
}

void MuniAphotDisplay::Draw(wxBitmap& bmp)
{
  int s = zoom;
  int i0 = width / 2;
  int j0 = height / 2;
  double r;

  wxMemoryDC dc;
  dc.SelectObject(bmp);

  wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
  if( gc ) {
    gc->SetAntialiasMode(wxANTIALIAS_NONE);
    gc->BeginLayer(1.0);

    // filled ring
    wxGraphicsPath gpath = gc->CreatePath();
    gpath.AddCircle(i0,j0,s*rmin);
    gpath.AddCircle(i0,j0,s*rmax);

    // default parameters
    wxColour cring(138,184,230,96);
    if( rmax < rmin )
      cring = wxColour(178,34,34,96);

    gc->SetPen(wxPen(cring,1));
    wxBrush bring(cring);
    gc->SetBrush(bring);
    gc->SetAntialiasMode(wxANTIALIAS_DEFAULT); //switch-off in fast-update mode?
    gc->FillPath(gpath);

    // apertures in blue
    gc->SetBrush(*wxTRANSPARENT_BRUSH);
    wxColour caper(138,184,230,196);
    gc->SetPen(wxPen(caper,1));
    for(int i = 1; i <= aper.GetNaper(); i++) {
      r = s*aper.GetAper(i);
      gc->DrawEllipse(i0-r,j0-r,2*r,2*r);
    }

    // the selected aperture
    gc->SetPen(wxPen(caper,3));
    r = s*aper.GetAper(saper);
    gc->DrawEllipse(i0-r,j0-r,2*r,2*r);

    // apertures in black
    wxColour black(65,105,225,128);
    gc->SetPen(wxPen(black,1));
    for(int i = 1; i <= aper.GetNaper()-1; i++) {
      r = s*aper.GetAper(i);
      gc->DrawEllipse(i0-r,j0-r,2*r,2*r);
    }

    // sky ring
    r = s*rmin;
    gc->DrawEllipse(i0-r,j0-r,2*r,2*r);
    r = s*rmax;
    gc->DrawEllipse(i0-r,j0-r,2*r,2*r);

    gc->EndLayer();
    delete gc;
  }
}

void MuniAphotDisplay::OnPaint(wxPaintEvent& e)
{
  wxPaintDC dc(this);
  if( canvas.IsOk() )
    dc.DrawBitmap(canvas,0,0);
}

void MuniAphotDisplay::OnSize(wxSizeEvent& event)
{
  wxSize size(event.GetSize());
  width = size.x;
  height = size.y;
  if( image.IsOk() && width > 17 && height > 17 )
    Rebase();

  event.Skip();
}

void MuniAphotDisplay::OnClick(MuniClickEvent& event)
{
  SetPosition(event.x,event.y);
  Rebase();
}

void MuniAphotDisplay::OnKlick(wxMouseEvent& event)
{
  wxClientDC dc(this);
  wxPoint p = event.GetLogicalPosition(dc);

  double x = p.x - width / 2;
  double y = p.y - height / 2;
  double r = sqrt(x*x + y*y) / zoom;

  int n;
  if( aper.FindAper(r,3,&n) ) {
    saper = n;
    Redraw();
  }

  MuniClickEvent ev(EVT_CLICK,GetId());
  ev.r = aper.GetAper(saper);
  wxQueueEvent(GetParent(),ev.Clone());
}

void MuniAphotDisplay::Rebase()
{
  wxASSERT(zoom > 0 && image.IsOk());

  int w = width / zoom;
  int h = height / zoom;

  int i1 = i0 - w/2;
  int j1 = image.GetHeight() - (j0+h/2);
  wxRect rect(i1, j1, w, h);
  wxRect irect(0,0,image.GetWidth(),image.GetHeight());
  wxRect r = rect.Intersect(irect);
  int xoff = i1 < 0 ? w - r.width : 0;
  int yoff = j1 < 0 ? h - r.height : 0;

  wxBitmap sub(image.GetSubImage(rect));
  wxMemoryDC mdc;
  mdc.SelectObject(sub);

  wxMemoryDC dc;
  wxBitmap bmp(width,height);
  dc.SelectObject(bmp);
  dc.SetBackground(*wxGREY_BRUSH);
  dc.Clear();
  dc.StretchBlit(zoom*xoff,zoom*yoff,zoom*r.width,zoom*r.height,&mdc,
		 0,0,r.width,r.height);
  dc.SelectObject(wxNullBitmap);
  subimage = bmp;
  canvas = bmp;

  aper = MuniApertures(naper,rmin,spiral);
  Draw(canvas);
  Refresh();
}

void MuniAphotDisplay::OnZoom(wxCommandEvent& event)
{
  zoom = event.GetInt();
  Rebase();
}

void MuniAphotDisplay::SetPosition(int i, int j)
{
  i0 = i;
  j0 = j;
  Rebase();
}
