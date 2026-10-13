/*

  Color images creator

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
#include <wx/filename.h>
#include <wx/tokenzr.h>

using namespace std;

MuniColouring::MuniColouring(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"colouring")
{
  wxString a;
  double x,y,r;

  // determine which character is used as the decimal point
  const wxString decimal = wxLocale::GetInfo(wxLOCALE_DECIMAL_POINT,
					     wxLOCALE_CAT_NUMBER);
  const wxString sep = decimal == "." ? ",;" : ";";

  // common options
  CommonOutputSingle(cmd);
  CommonOptionsBitpix(cmd);

  if( cmd.Found("ctable",&a) )
    Write("CTABLE = '" + a + "'");
  else if( FindDataDir("ctable.dat",&a) )
    Write("CTABLE = '" + a + "'");

  if( cmd.Found("phsystab",&a) )
    Write("PHSYSTABLE = '" + a + "'");
  else if( FindDataDir("photosystems.fits",&a) )
    Write("PHSYSTABLE = '" + a + "'");

  if( cmd.Found("list") ) {
    Write("LIST = T");
    return;
  }

  wxArrayString files = cmd.GetFiles();

  if( cmd.Found("c",&a) || cmd.Found("cspace-input",&a) )
    Write("COLOURSPACE = '" + a + "'");

  if( cmd.Found("cspace-output",&a) )
    Write("OUTERSPACE = '" + a + "'");

  if( cmd.Found("disable-back") )
    Write("ESTIMBACKS = F");

  if( cmd.Found("white-radius",&r) )
    Write("WRADIUS = %f",r);

  bool wspot = cmd.Found("white-spot",&a);
  bool wstar = cmd.Found("white-star",&a);
  if( wspot || wstar ) {
    wxString b,c;
    wxStringTokenizer tok(a,sep);
    b = tok.GetNextToken();
    c = tok.GetNextToken();
    if( ! b.IsEmpty() && ! c.IsEmpty() && b.ToDouble(&x) && c.ToDouble(&y) ) {
      wxString label(wspot ? "WSPOT" : "WSTAR");
      Write(label + " = %f %f",x,y);
    }
    else {
      MuniFatalError("Failed to interpred `"+a+"' as cartesian coordinates.");
      return;
    }
  }

  if( cmd.Found("w",&a) || cmd.Found("weights",&a) ) {
    wxString line;
    size_t n;
    dblstr(a,sep,n,line);
    if( n != files.GetCount() ) {
      MuniFatalError("Count of weight items differs to passed filenames.");
      return;
    }

    Write("NWEIGHTS = %d",int(n));
    Write("WEIGHTS = "+line);
  }

  if( cmd.Found("q",&a) || cmd.Found("ctphs",&a) ) {
    wxString line;
    size_t n;
    dblstr(a,sep,n,line);
    if( n != files.GetCount() ) {
      MuniFatalError("Count of ctph items differs to passed filenames.");
      return;
    }

    Write("NCTPHS = %d",int(n));
    Write("CTPHS = "+line);
  }

  if( cmd.Found("b",&a) || cmd.Found("backs",&a) ) {
    wxString line;
    size_t n;
    dblstr(a,sep,n,line);
    if( n != files.GetCount() ) {
      MuniFatalError("Count of background items differs to filenames.");
      return;
    }

    Write("NBACKS = %d",int(n));
    Write("BACKS = "+line);
  }

  Write("NBANDS = %d",int(files.GetCount()));
  for(size_t i = 0; i < files.GetCount(); i++)
    Write("BAND = '"+files[i]+"'");

}
