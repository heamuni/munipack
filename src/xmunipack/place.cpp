/*

  xmunipack - base of places


  Copyright © 2009-2011 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/wx.h>
#include <wx/dcbuffer.h>
#include <wx/graphics.h>


// MuniPlaceBase::MuniPlaceBase(wxWindow *w, const wxPoint& pos, 
// 			     const wxSize& size, MuniConfig *c):
//   wxScrolledWindow(w,wxID_ANY,pos,size),
//   topwin(GetParent()),config(c),hdusel(-2)
// {
//   wxASSERT(topwin);
//   //  SetBackgroundStyle(wxBG_STYLE_CUSTOM);
//   //  SetBackgroundStyle(wxBG_STYLE_PAINT);

//   //Connect(wxEVT_ERASE_BACKGROUND,wxEraseEventHandler(MuniPlaceBase::OnErase));
//   Connect(wxEVT_PAINT,wxPaintEventHandler(MuniPlaceBase::OnPaint));
// }

// void MuniPlaceBase::Assign(const FitsFile& f, int sel, const wxImage& i)
// {
//   wxASSERT(f.IsOk());
//   fits = f;
//   hdusel = sel;
// }

// void MuniPlaceBase::Assign(const FitsFile& f)
// {
//  //  wxASSERT(f.IsOk() && meta.IsOk() && meta.GetFullPath() == f.GetFullPath());
//   // test all header content for equivalence ?
// #ifdef __WXDEBUG__
//   if( f.IsOk() && meta.IsOk() ) {
//     wxASSERT(meta.GetFullPath() == f.GetFullPath());
//   }
// #endif
//   fits = f;
// }

// void MuniPlaceBase::Assign(const FitsMeta& m)
// {
//   meta = m;
// }

// void MuniPlaceBase::SetHdu(int i)
// {
//   wxASSERT(meta.IsOk());
//   hdusel = i;
// }

// int MuniPlaceBase::GetHdu() const
// {
//   return hdusel;
// }

// FitsFile MuniPlaceBase::GetFits() const
// {
//   return fits;
// }

// FitsMeta MuniPlaceBase::GetMeta() const
// {
//   return meta;
// }

// void MuniPlaceBase::OnErase(wxEraseEvent& event) 
// {
//   wxLogDebug("MuniPlaceBase::OnErase");

//   wxBrush brush(*wxLIGHT_GREY,wxCROSSDIAG_HATCH);

//   wxClientDC *dc = (wxClientDC *) event.GetDC();
//   if( dc ) {
//     dc->SetBackground(brush);
//     dc->Clear();
//   }
//   else {
//     dc = new wxClientDC(this);
//     dc->SetBackground(brush);
//     dc->Clear();
//     delete dc;
//   }
// }

// void MuniPlaceBase::OnPaint(wxPaintEvent& event) 
// {
//   //  wxLogDebug("MuniPlaceBase::OnPaint");

//   wxBitmap bitmap(config->munipack_icon);
//   if( bitmap.IsOk() ) {
//     wxSize canvas = GetClientSize();
//     int xoff = (canvas.GetWidth() - bitmap.GetWidth())/2;
//     int yoff = (canvas.GetHeight() - bitmap.GetHeight())/2;
    
//     if( xoff < 0 ) xoff = 0;
//     if( yoff < 0 ) yoff = 0;

//     wxBrush brush(*wxLIGHT_GREY,wxCROSSDIAG_HATCH);

//     wxPaintDC dc(this);
//     dc.SetBackground(brush);
//     dc.Clear();
//     dc.DrawBitmap(bitmap,xoff,yoff,false);
//   }

//   //  event.Skip();
// }

// void MuniPlaceBase::OnMenu(wxMouseEvent& WXUNUSED(event)) {}
// void MuniPlaceBase::OnMouseMotion(wxMouseEvent& event) {}
// void MuniPlaceBase::OnKeyDown(wxKeyEvent& event) {}
// void MuniPlaceBase::Render() {}
// bool MuniPlaceBase::IsRendering() const { return false; }
// void MuniPlaceBase::OnLeaveFullscreen(wxCommandEvent& WXUNUSED(event)) {}
// void MuniPlaceBase::OnEnterWin(wxMouseEvent& event) {}
// void MuniPlaceBase::OnLeaveWin(wxMouseEvent& event) {}
// void MuniPlaceBase::OnConfigUpdated(wxCommandEvent& event) {}
// wxSize MuniPlaceBase::DoGetBestSize() const { return wxSize(162,100); }
