/*

  xmunipack - Thumbnails canvas

  Copyright © 2019-2022 F.Hroch (hroch@physics.muni.cz)

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

#include "thumbcanvas.h"
#include <wx/wx.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

// ------   MuniThumbCanvas


MuniThumbCanvas::MuniThumbCanvas(wxWindow *w, const wxBitmap& b):
  wxWindow(w,wxID_ANY,wxDefaultPosition,wxSize(b.GetWidth()+1,b.GetHeight()+1)),
  icon(b)
{
  SetBackgroundStyle(wxBG_STYLE_PAINT);
  Bind(wxEVT_PAINT,&MuniThumbCanvas::OnPaint,this);
}

void MuniThumbCanvas::OnPaint(wxPaintEvent& WXUNUSED(event))
{
  if( icon.IsOk() ) {
    wxSize s = GetSize();
    int xoff = (s.GetWidth() - icon.GetWidth())/2;
    int yoff = (s.GetHeight() - icon.GetHeight())/2;
    wxPaintDC dc(this);
    dc.DrawBitmap(icon,xoff,yoff,true);
  }
}

void MuniThumbCanvas::SetIcon(const wxBitmap& i)
{
  icon = i;
  Refresh();
}

void MuniThumbCanvas::SetIcon(const wxImage& i)
{
  SetIcon(wxBitmap(i));
}
