/*

  Preview -- a panel which shows an image during loading

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

#ifndef _XMUNIPACK_DISPREVIEW_H_
#define _XMUNIPACK_DISPREVIEW_H_

#include "fits.h"
#include "event.h"
#include "mconfig.h"
#include <wx/wx.h>
#include <vector>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


class MuniDisplayPreview: public wxWindow
{
  int shrink, xoff, yoff, width, height;
  wxSize display_size;
  wxImage canvas;

  void OnSize(wxSizeEvent&);
  void OnPaint(wxPaintEvent&);
  void InitCanvas();

 public:
  MuniDisplayPreview(wxWindow *, const MuniConfig *);
  wxImage GetCanvas() const;
  void SetZoom(double,const wxSize&);
  void RefreshCanvas(int, const wxImage&);

  void SetZoom(double);
  void SetImageSize(const wxSize&);

};

#endif
