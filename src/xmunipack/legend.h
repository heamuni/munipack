/*

  xmunipack - display legend of the display

  Copyright © 2022 F.Hroch (hroch@physics.muni.cz)

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

#include "fits.h"
#include "mconfig.h"
#include <wx/wx.h>
#include <wx/graphics.h>
#include <vector>

class MuniDisplayLegend
{
  wxGraphicsContext *gc;
  const MuniConfig *config;
  const wxRect exposed_area;
  const double xoff, yoff, zoom;
  wxDouble height, width;
  std::vector<double> atics, dtics;
  wxColour fgcolour;

  double X(double) const;
  double Y(double) const;
  double W(double) const;
  double H(double) const;

  void DrawRectangle(wxDouble,wxDouble,wxDouble,wxDouble);
  void DrawLine(wxDouble,wxDouble,wxDouble,wxDouble);
  void DrawText(const wxString&,wxDouble,wxDouble,bool=true);
  void GetTextDim(const wxString&, wxDouble&, wxDouble&) const;
  void DrawAxisLabel(const wxString&,wxDouble,wxDouble,wxDouble=0.0,
		     double=0.0,double=0.0);
  bool DrawGrid(const FitsArray&);
  bool DrawRuler(const FitsArray&);
  bool DrawArrow(const FitsArray&);

  void TicsArray(double, double, const std::vector<double>&, std::vector<double>&,
		 std::vector<double>&, double&) const;
  void DrawMeridians(const FitsProjection&, double, double, const std::vector<double>&, double);
  void DrawParallels(const FitsProjection&, double, double, const std::vector<double>&, double);
  void DrawTicsTop(const FitsProjection&, double, double, double, const std::vector<double>&, double, double);
  void DrawTicsBottom(const FitsProjection&, double, double, double, const std::vector<double>&, double, double);
  void DrawTicsLeft(const FitsProjection&, double, double, double, const std::vector<double>&, double, double);
  void DrawTicsRight(const FitsProjection&, double, double, double, const std::vector<double>&, double, double);
  wxString RA(double, double) const;
  wxString Dec(double, double) const;
  void DrawLabelsBottom(const FitsProjection&,double,double,const std::vector<double>&,double,double,double);
  void DrawLabelsLeft(const FitsProjection&,double,double,const std::vector<double>&,double,double,double);
  void DrawRa(const wxString&, const wxString&, const wxString&, double, double);


 public:
  MuniDisplayLegend(wxGraphicsContext *, const MuniConfig *, const wxRect&,
		    double, double, double);
  bool Draw(const FitsArray&);

};
