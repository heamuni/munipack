/*

  xmunipack - image captions

  Copyright © 2021-22 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _XMUNIPACK_CAPTION_H_
#define _XMUNIPACK_CAPTION_H_

#include "fits.h"
#include "event.h"
#include "mconfig.h"
#include <wx/wx.h>


class MuniDisplayCaptionInfo: public wxPanel
{
  MuniConfig *config;
  bool init;
  wxStaticText *object, *colour, *exptime, *date, *time,
    *label_object, *label_colour;
  FitsArray array;

  void OnIdle(wxIdleEvent&);

public:
  MuniDisplayCaptionInfo(wxWindow *, MuniConfig *);
  void SetArray(const FitsArray&);
  void ConfigUpdate();

};

class MuniDisplayCaptionMotion: public wxPanel
{
  MuniConfig *config;
  bool init, update, hascal, show_wcs;
  int x,y;
  wxStaticText *quantity, *intensity, *magnitude, *label_quantity,
    *alpha, *delta, *xpix, *ypix;
  wxFlexGridSizer *quantitysizer, *qsizer, *pixsizer, *coosizer, *coloursizer;
  wxStaticText *colours_XYZ[3];
  FitsValue value;
  FitsCoo coords;
  FitsArray array;

  void OnIdle(wxIdleEvent&);
  void OnMouseMotion(MuniSlewEvent&);
  void InitArray();

public:
  MuniDisplayCaptionMotion(wxWindow *, MuniConfig *);
  void SetArray(const FitsArray&);
  void ConfigUpdate();

};


class MuniDisplayCaption: public wxPanel
{
  bool inside;

  MuniDisplayCaptionInfo *info;
  MuniDisplayCaptionMotion *motion;

  void OnMouseMotion(MuniSlewEvent&);

public:
  MuniDisplayCaption(wxWindow *, MuniConfig *);
  void SetArray(const FitsArray&);
  void ConfigUpdate();

};

#endif
