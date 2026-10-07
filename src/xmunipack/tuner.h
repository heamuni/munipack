/*

  xmunipack - tunner


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


#ifndef _XMUNIPACK_TUNER_H
#define _XMUNIPACK_TUNER_H

#include "fits.h"
#include "event.h"
#include <wx/wx.h>

class MuniTunerAdjuster: public wxWindow
{
  wxSlider *slider;
  int irange;
  double rmin, rmax;

  double GetValue(int) const;
  void OnScroll(wxScrollEvent&);

public:
  MuniTunerAdjuster(wxWindow *,wxWindowID,double, double, int,
		    const wxString& ="", const wxString& ="");

  void SetValue(double);
  double GetValue() const;
};

class MuniTunerLogjuster: public wxWindow
{
  wxSlider *slider;
  double scale, rmin, rmax;

  double GetValue(int) const;
  void OnScroll(wxScrollEvent&);

public:
  MuniTunerLogjuster(wxWindow *,wxWindowID,double, double, int,
		     const wxString& ="", const wxString& ="");

  void SetValue(double);
  double GetValue() const;
};


class MuniTuner: public wxPanel
{
  MuniTunerAdjuster *adjblack;
  MuniTunerLogjuster *adjsense;

public:
  MuniTuner(wxWindow *, wxWindowID);
  void SetTone(double, double);
};


#endif
