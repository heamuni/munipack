/*

  VOTable to SVG convertor

  Copyright © 2018-2026 F.Hroch (hroch@physics.muni.cz)

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

#include "votable.h"
#include <wx/wx.h>
#include <wx/txtstrm.h>
#include <cmath>
#include <vector>

using namespace std;

SVGcanvas::SVGcanvas(const wxString& filename): VOTable(filename)
{
  proj_alpha = 0;
  proj_delta = 0;
  proj_scale = 1000;
  mag_limit = 15;
  canvas_width = 500;
  canvas_height = 500;
}


void SVGcanvas::SetProjection(const wxString& a) { proj_type = a;}
void SVGcanvas::SetProjectionCenter(double a,double d)
{
  proj_alpha = a;
  proj_delta = d;
}
void SVGcanvas::SetCanvasSize(long w, long h)
{
  canvas_width = w;
  canvas_height = h;
}
void SVGcanvas::SetScale(double c) { proj_scale = c; }
void SVGcanvas::SetMaglim(double x) { mag_limit = x; }
void SVGcanvas::SetMagkey(const wxString& a) { mag_key = a; }
void SVGcanvas::SetAlphakey(const wxString& a) { alpha_key = a; }
void SVGcanvas::SetDeltakey(const wxString& a) { delta_key = a; }

bool SVGcanvas::SaveSVG(wxOutputStream& output) const
{
  if( ! IsOk() ) return false;
  if( HasError() ) return false;

  if( resources.size() == 0 )
    return false;

  // We are loosers, only the first table is converted.
  const VOTableTable table(resources[0].table);

  int nra = -1, ndec = -1, nmag = -1;
  vector<wxString> fs = table.GetFields();
  for(vector<wxString>::const_iterator i = fs.begin(); i != fs.end(); ++i) {
    if( *i == alpha_key ) {
      nra = i - fs.begin();
    }
    else if( *i == delta_key ) {
      ndec = i - fs.begin();
    }
    else if( *i == mag_key ) {
      nmag = i - fs.begin();
    }
  }

  if( nra == -1 || ndec == -1 || nmag == -1 ) {
    wxLogError("RA, Dec of magnitude column has no data.");
    return false;
  }

  vector<wxString> xra = table.GetColumn(nra);
  vector<wxString> xdec = table.GetColumn(ndec);
  vector<wxString> xmag = table.GetColumn(nmag);

  wxTextOutputStream cout(output);
  cout << "<svg xmlns=\"http://www.w3.org/2000/svg\""
       << " width=\"" << int(canvas_width)
       << "\" height=\"" << int(canvas_height) << "\">" << endl;

  // WARNING: only simple gnonomical projection is implemented

  double c = xdec.size() > 0 ? cos(ToDouble(xdec[0])/57.3) : 1.0;
  for(size_t i = 0; i < xra.size(); i++) {
    double ra = ToDouble(xra[i]);
    double dec = ToDouble(xdec[i]);
    double mag = ToDouble(xmag[i]);
    double x = canvas_width/2 - proj_scale*(ra - proj_alpha)*c;
    double y = canvas_height/2 - proj_scale*(dec - proj_delta);//reverse y-axis!
    double r = 3.0*pow(10.0,0.11*(mag_limit - mag));

    if( 0 <= x && x <= canvas_width && 0 <= y && y <= canvas_height )
      cout << "<circle cx=\"" << x << "\" cy=\"" << y
	   << "\" r=\"" << r << "\"/>" << endl;
  }

  cout << "</svg>" << endl;
  return true;
}

double SVGcanvas::ToDouble(const wxString& a) const
{
  double x;
  a.ToDouble(&x);
  return x;
}
