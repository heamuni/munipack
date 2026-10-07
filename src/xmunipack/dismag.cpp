/*

  xmunipack - magnifier

  Copyright © 2019-21 F.Hroch (hroch@physics.muni.cz)

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

  * There're unsolved mystery in sizer, layer algorithm which
    prevent setup of an initial size of this window to value
    of previously adjusted by user.

*/

#include "dismag.h"
#include "fits.h"
#include <wx/wx.h>
#include <wx/graphics.h>

MuniMagnifier::MuniMagnifier(wxWindow *w, MuniConfig *config):
  wxWindow(w,wxID_ANY), scale(config->magnifier_scale), zoom(1.0), inside(false)
{
  SetBackgroundStyle(wxBG_STYLE_PAINT);
  bgcolour = wxColour(config->display_bgcolour);

  Bind(wxEVT_PAINT,&MuniMagnifier::OnPaint,this);
  Bind(wxEVT_UPDATE_UI,&MuniMagnifier::OnUpdate,this);
  Bind(EVT_SLEW,&MuniMagnifier::OnMouseMotion,this);
}

void MuniMagnifier::SetImage(const wxImage& i)
{
  wxASSERT(i.IsOk());
  //  wxLogDebug("MuniMagnifierWindow::SetImage %d %d %d",i.IsOk(),i.GetWidth(),i.GetHeight());

  image = i;
}


void MuniMagnifier::OnUpdate(wxUpdateUIEvent& event)
{
  /*
  wxLogDebug("MuniMagnifierWindow::OnUpdate %d %d %d %d",crosshair.x,crosshair.y,
	     int(image.IsOk()),int(canvas.IsOk()));
  */

  event.SetUpdateInterval(40); // = 1/25 s
  Refresh();
}

void MuniMagnifier::OnPaint(wxPaintEvent& event)
{
  wxSize size(GetSize());
  int w = size.GetWidth() / scale;
  int h = size.GetHeight() / scale;
  int i = size.GetWidth() / (2*scale);
  int j = size.GetHeight() / (2*scale);
  int x1 = crosshair.x - i;
  int y1 = crosshair.y - j;


  wxBrush brush(bgcolour);
  wxPen pen(bgcolour);

  wxPaintDC dc(this);
  wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
  if( gc ) {

    gc->SetBrush(brush);
    gc->SetPen(pen);
    gc->DrawRectangle(0,0,size.GetWidth(),size.GetHeight());

    if( image.IsOk() && inside) {

      wxRect rsub(x1,y1,w,h);
      wxRect rimg(0,0,image.GetWidth(),image.GetHeight());
      wxRect r = rimg.Intersect(rsub);

      if( ! r.IsEmpty() ) {
	wxImage sub(image.GetSubImage(r));
	sub.Rescale(sub.GetWidth()*scale,sub.GetHeight()*scale,
		    wxIMAGE_QUALITY_NEAREST);

	// offset
	int xoff = x1 < 0 ? wxMax(size.GetWidth() - sub.GetWidth(),0) : 0;
	int yoff = y1 < 0 ? wxMax(size.GetHeight()-sub.GetHeight(),0) : 0;

	wxGraphicsBitmap bmp = gc->CreateBitmapFromImage(sub);
	gc->DrawBitmap(bmp,xoff,yoff,sub.GetWidth(),sub.GetHeight());
      }
    }

    int alpha = 200;
    int s = scale * (zoom > 1.0 ? int(zoom + 0.5) : 1);
    int ss = wxMax(int(scale * zoom),1);
    int i1 = size.GetWidth() / 2 - ss / 2;
    int j1 = size.GetHeight()/ 2 - ss / 2;
    gc->SetBrush(*wxTRANSPARENT_BRUSH);

    wxColour whitea(255,255,255,alpha);
    gc->SetPen(wxPen(whitea));
    gc->DrawRectangle(i1,j1,s,s);

    wxColour blacka(0,0,0,alpha);
    gc->SetPen(wxPen(blacka));
    gc->DrawRectangle(i1-1,j1-1,s+2,s+2);

    delete gc;
  }
}


void MuniMagnifier::OnMouseMotion(MuniSlewEvent& event)
{
  //  wxLogDebug("MuniMagnifierWindow::OnMouseMotion %d %d",event.xsub,event.ysub);
  inside = event.inside;
  crosshair = wxPoint(event.xsub,event.ysub);
  zoom = event.zoom;
}
