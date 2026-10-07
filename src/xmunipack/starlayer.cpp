/*

  Star layer

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

#include "xmunipack.h"
#include <wx/wx.h>
#include <vector>

using namespace std;

MuniStarLayer::MuniStarLayer():
  gold(wxColour(255,215,0,240)),DarkOrange2(wxColour(238,118,0)),hwhm(1) {}

MuniLayer MuniStarLayer::GetLayer() const
{
  MuniLayer layer(ID_PHOTOMETRY,objects);
  return layer;
}


void MuniStarLayer::DrawObjects(const vector<double>& xcoo,
				const vector<double>& ycoo,
				const vector<double>& flux)
{
  wxASSERT(xcoo.size() == ycoo.size() && flux.size() == xcoo.size());

  objects.push_back(new MuniDrawFont(*wxNORMAL_FONT,gold));
  objects.push_back(new MuniDrawPen(wxPen(gold,1.8)));
  objects.push_back(new MuniDrawBrush(wxColour(90,90,255,190)));

  double fmax = 0.0;
  for(size_t i = 0; i < flux.size(); i++)
    if( flux[i] > 0.0 && flux[i] > fmax )
      fmax = flux[i];

  double fmin = fmax;
  for(size_t i = 0; i < flux.size(); i++)
    if( flux[i] > 0.0 && flux[i] < fmin )
      fmin = flux[i];

  const double size1 = hwhm > 1 ? 3*hwhm : 3;
  const double size2 = 2.0;
  const double sizec = 3.0;
  double r1 = (size1 - size2)/(fmax - fmin);

  for(size_t i = 0; i < xcoo.size(); i++) {

    if( flux[i] > 0.0 ) {
      double r = r1*(flux[i] - fmin) + size2;
      objects.push_back(new MuniDrawCircle(xcoo[i],ycoo[i],r));
    }
    else
      objects.push_back(new MuniDrawCross(xcoo[i],ycoo[i],sizec));
  }
}
