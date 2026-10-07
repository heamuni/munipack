/*

  xmunipack - preview

  Copyright © 2021-2022 F.Hroch (hroch@physics.muni.cz)

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

#include "dispreview.h"
#include "zoomer.h"
#include "types.h"
#include "event.h"
#include <wx/wx.h>
#include <vector>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

MuniDisplayPreview::MuniDisplayPreview(wxWindow *w, const MuniConfig *config):
  wxWindow(w,wxID_ANY), shrink(0),xoff(0), yoff(0), width(0), height(0),
  display_size(wxDefaultSize)
{
  SetBackgroundColour(config->display_bgcolour);
  Bind(wxEVT_PAINT,&MuniDisplayPreview::OnPaint,this);
  Bind(wxEVT_SIZE,&MuniDisplayPreview::OnSize,this);
}

void MuniDisplayPreview::OnSize(wxSizeEvent& event)
{
  wxSize size(event.GetSize());
  if( size != display_size ) {
    wxLogDebug("MuniDisplayPreview::OnSize changed %dx%d",size.GetWidth(),size.GetHeight());
    display_size = size;

    MuniSizeChangedEvent event(EVT_SIZE_CHANGED);
    event.size = size;
    wxQueueEvent(GetParent(),event.Clone());
  }
  event.Skip();
}

void MuniDisplayPreview::OnPaint(wxPaintEvent& event)
{
  if( canvas.IsOk() ) {
    wxPaintDC dc(this);
    dc.DrawBitmap(wxBitmap(canvas),xoff,yoff);
  }
}

void MuniDisplayPreview::InitCanvas()
{
  const wxSize size(GetSize());
  wxASSERT(width > 0 && height > 0 && size.GetWidth() > 0 && size.GetHeight() > 0 &&
	   shrink > 0);

  MuniZoomSet zset;
  double z = zset.SetBestFitZoom(size,wxSize(width,height));
  if( z < 1.0 )
    shrink = int(1.0 / z + 0.5);
  else
    shrink = 1;

  int w = wxMax(width / shrink,1);
  int h = wxMax(height / shrink,1);

  wxASSERT(w > 0 && h > 0);
  canvas = wxImage(w,h);
  xoff = wxMax(size.GetWidth() - w,0)/2;
  yoff = wxMax(size.GetHeight() - h,0)/2;

  //  wxLogDebug("%d %d %d %d %d %d %d %d %d",shrink,w,h,width,height,xoff,yoff,
  //	     size.GetWidth(),size.GetHeight());
  wxASSERT(xoff >= 0 && yoff >= 0);
  // what about removing of the negative poisitons from the image?
}

void MuniDisplayPreview::RefreshCanvas(int crow, const wxImage& subwin)
{
  if( ! (subwin.GetWidth() > 0 && subwin.GetHeight() > 0) )
    return;

  wxASSERT(canvas.GetWidth() > 0 && canvas.GetHeight() > 0);
  int y = canvas.GetHeight() - crow / shrink;
  canvas.Paste(subwin,0,y);
  RefreshRect(wxRect(xoff+0,yoff+y,subwin.GetWidth(),subwin.GetHeight()));
}

wxImage MuniDisplayPreview::GetCanvas() const
{
  return canvas;
}

void MuniDisplayPreview::SetZoom(double zoom)
{
  wxASSERT(zoom > 0);

  if( zoom < 1.01 )
    shrink = int(1.0/zoom + 0.5);
  else
    shrink = 1;
}


void MuniDisplayPreview::SetImageSize(const wxSize& size)
{
  wxASSERT(size.IsFullySpecified());
  width = size.GetWidth();
  height = size.GetHeight();
  InitCanvas();
}

void MuniDisplayPreview::SetZoom(double zoom, const wxSize& size)
{
  if( shrink > 0 ) return;
  SetZoom(zoom);
  SetImageSize(size);
}
