/*

  Growth-Curve Photometry

  Copyright © 2016 - 2026 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/log.h>
#include <wx/filename.h>
#include <wx/tokenzr.h>

using namespace std;


MuniGphot::MuniGphot(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"gphot")
{
  // determine character representing the decimal point
  const wxString decimal = wxLocale::GetInfo(wxLOCALE_DECIMAL_POINT,
					     wxLOCALE_CAT_NUMBER);
  const wxString sep = decimal == "." ? ",;" : ";";

  CommonOutputMultiple(cmd);

  double x;
  long n;
  wxString a;

  bool legacy = cmd.Found("legacy");
  if( legacy )
    Write("METHOD = 'LEGACY'");

  if( cmd.Found("th",&x) || cmd.Found("threshold",&x) ) {
    if( ! legacy )
      MuniFatalError("--threshold applies to --legacy only.");
    Write("THRESHOLD = %.17g",x);
  }

  if( cmd.Found("var",&n) ) {
    if( n < 0 || n > 2 )
      MuniFatalError("--var is 0, 1 or 2.");
    if( legacy )
      MuniFatalError("--var and --legacy are exclusive.");
    Write("EEVAR = %ld",n);
  }

  if( cmd.Found("radial") ) {
    if( legacy )
      MuniFatalError("--radial and --legacy are exclusive.");
    Write("EERADIAL = T");
  }

  if( cmd.Found("centre",&a) ) {
    if( ! cmd.Found("radial") )
      MuniFatalError("--centre requires --radial.");
    wxStringTokenizer t(a,sep);
    double c[2];
    int i = 0;
    while( t.HasMoreTokens() && i < 2 ) {
      wxString token(t.GetNextToken());
      if( ! token.ToDouble(&c[i]) )
	MuniFatalError("Failed to interpret `%s' as a coordinate.",token);
      i++;
    }
    if( i != 2 || t.HasMoreTokens() )
      MuniFatalError("--centre requires x,y.");
    Write("EECENTRE = %.17g %.17g",c[0],c[1]);
  }

  WriteFilesMultiple(cmd.GetFiles());
}
