/*

  xmunipack - display's headers

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

#ifndef _XMUNIPACK_DISMAG_H
#define _XMUNIPACK_DISMAG_H

#include "event.h"
#include "mconfig.h"
#include <wx/wx.h>


class MuniMagnifier: public wxWindow
{
  int scale;
  double zoom;
  bool inside;
  wxColour bgcolour;
  wxImage image;
  wxPoint crosshair;

  void OnPaint(wxPaintEvent&);
  void OnUpdate(wxUpdateUIEvent&);
  void OnMouseMotion(MuniSlewEvent&);

public:
  MuniMagnifier(wxWindow *, MuniConfig *);
  void SetImage(const wxImage&);

};

#endif
