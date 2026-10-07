/*

  timeseries (light curves)

  Copyright © 2012 - 2026 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/intl.h>
#include <vector>

using namespace std;

MuniTimeseries::MuniTimeseries(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"timeseries")
{
  wxString s;
  double x;

  // determine which character is used as the decimal point
  const wxString decimal = wxLocale::GetInfo(wxLOCALE_DECIMAL_POINT,
					     wxLOCALE_CAT_NUMBER);
  const wxString sep = decimal == "." ? ",;" : ";";

  CommonOutputSingle(cmd);

  if( cmd.Found("c",&s) || cmd.Found("coo",&s) ) {

    vector<wxString> coo;
    wxStringTokenizer objtokenizer(s," |");
    while( objtokenizer.HasMoreTokens() ) {
      wxString a(objtokenizer.GetNextToken());
      coo.push_back(a);
    }

    vector<double> a,d;
    for(size_t n = 0; n < coo.size(); n++ ) { // over objects
      double q;
      wxString token;
      wxStringTokenizer tokenizer(coo[n],sep);

      // RA
      token = tokenizer.GetNextToken();
      if( ! token.IsEmpty() && token.ToDouble(&q) )
	a.push_back(q);
      else
	wxLogFatalError("Failed to interpret `%s' as Right Ascension.",token);

      // Dec
      token = tokenizer.GetNextToken();
      if( ! token.IsEmpty() && token.ToDouble(&q) )
	d.push_back(q);
      else
	wxLogFatalError("Failed to interpret `%s' as Declination.",token);
    }

    Write("NCOO = %d",a.size());
    for(size_t i = 0; i < a.size(); i++)
      Write("COO = %.15f %.15f",a[i],d[i]);
  }

  if( cmd.Found("l",&s) || cmd.Found("col",&s) ) {

    wxString l;
    long n;
    apstr(s,n,l);
    Write("NCOL = %d",n);
    Write("COL = " + l);
  }

  if( cmd.Found("K",&s) || cmd.Found("keys",&s) ) {

    wxString l;
    long n;
    apstr(s,n,l);
    Write("NKEY = %d",n);
    Write("KEY = " + l);
  }

  if( cmd.Found("coo-col",&s) ) {

    wxString l;
    long n;
    apstr(s,n,l);
    Write("NCOOCOL = %d",n);
    Write("COOCOL = " + l);
  }

  if( cmd.Found("T",&s) || cmd.Found("time-type",&s) )
    Write("TIMETYPE = '" + s.MakeUpper() + "'");
  else if( cmd.Found("lc-epoch",&x) || cmd.Found("lc-period",&x) )
    Write("TIMETYPE = 'PHASE'");

  if( cmd.Found("time-stamp",&s) )
    Write("TIMESTAMP = '" + s.MakeUpper() + "'");

  if( cmd.Found("cat",&s) )
    Write("CATALOGUE = '" + s + "'");

  if( cmd.Found("lc-epoch",&x) )
    Write("LC_EPOCH = %e",x);

  if( cmd.Found("lc-period",&x) )
    Write("LC_PERIOD = %e",x);

  if( cmd.Found("tol",&x) )
    Write("TOL = %e",x);

  if( cmd.Found("coo-type",&s) )
    Write("COOTYPE = '" + s.MakeUpper() + "'");

  if( cmd.Found("extname",&s) )
    Write("EXTNAME = '" + s + "'");

  if( cmd.Found("stdout") )
    Write("STDOUT = T");

  if( cmd.Found("enable-filename") )
    Write("PRINTNAME = T");

  if( cmd.Found("disable-timetype") )
    Write("PRINTTIME = F");

  if( cmd.Found("enable-horizon") )
    Write("PRINTHORIZON = T");

  if( cmd.Found("enable-airmass") )
    Write("PRINTAIRMASS = T");

  WriteFilesSingle(cmd.GetFiles());
}
