/*

  Aperture Photometry

  Copyright © 2010 - 2026 F.Hroch (hroch@physics.muni.cz)

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

#include "munipack.h"
#include <wx/wx.h>
#include <wx/tokenzr.h>
#include <vector>

using namespace std;


MuniAphot::MuniAphot(wxEvtHandler *h, const MuniCmdLineParser& cmd):
    MuniAction(h,"aphot")
{
  // determine character representing the decimal point
  const wxString decimal = wxLocale::GetInfo(wxLOCALE_DECIMAL_POINT,
					     wxLOCALE_CAT_NUMBER);
  const wxString sep = decimal == "." ? ",;" : ";";

  CommonOutputMultiple(cmd);

  wxString a;
  long n;
  double x;

  if( cmd.Found("apertures",&a) ) {

    n = 0;
    wxStringTokenizer t(a,sep);
    while ( t.HasMoreTokens() ) {
      t.GetNextToken();
      n++;
    }

    Write("NAPER = %ld",n);
    Write("APER = "+a);
  }

  if( cmd.Found("ring",&a) )
    Write("RING = "+a);

  if( cmd.Found("enable-ellipticity") )
    Write("ELLIPTICITY = T");

  if( cmd.Found("eccentricity",&x) )
    Write("ECCENTRICITY = %lf",x);

  if( cmd.Found("inclination",&x) )
    Write("INCLINATION = %lf",x);

  if( cmd.Found("coo",&a) ) {

    vector<wxString> coo;
    wxStringTokenizer objtokenizer(a," |");
    while( objtokenizer.HasMoreTokens() ) {
      wxString l(objtokenizer.GetNextToken());
      coo.push_back(l);
    }

    vector<double> x,y;
    for(size_t n = 0; n < coo.size(); n++ ) { // over objects
      double q;
      wxString token;
      wxStringTokenizer tokenizer(coo[n],sep);

      token = tokenizer.GetNextToken();
      if( ! token.IsEmpty() && token.ToDouble(&q) )
	x.push_back(q);
      else
	wxLogFatalError("Failed to interpret `%s' as a coordinate.",token);

      token = tokenizer.GetNextToken();
      if( ! token.IsEmpty() && token.ToDouble(&q) )
	y.push_back(q);
      else
	wxLogFatalError("Failed to interpret `%s' as a coordinate.",token);
    }

    Write("NSTARS = %d",x.size());
    for(size_t i = 0; i < x.size(); i++)
      Write("STAR = %.15f %.15f",x[i],y[i]);
  }

  WriteFilesMultiple(cmd.GetFiles());
}
