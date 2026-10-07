/*

  Cross-match

  Copyright © 2013 - 2026 F.Hroch (hroch@physics.muni.cz)

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

using namespace std;


MuniCross::MuniCross(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"cross")
{
  CommonOutputSingle(cmd);

  wxString a,c/*,vocat*/;
  long n;
  double x;

  /*
  if( cmd.Found("vocat",&a) )
    vocat = a;
  else if( FindDataDir("VOcat_conf.xml",&a) )
    vocat = a;

  VOCatConf catalogs(vocat);
  VOCatResources cat(catalogs.GetCatFits(cone));
  */

  if( cmd.GetParamCount() != 2 )
    wxLogFatalError("Cross-match needs exactly two tables.");

  if( cmd.Found("tol",&x) )
    Write("TOL = %e",x);

  if( cmd.Found("ftol",&x) )
    Write("FTOL = %e",x);

  /*
  if( cmd.Found("col-ra",&a) )
    Write("COL_RA = '" + a + "'");
  else if( cat.IsOk() )
    Write("COL_RA = '" + cat.GetLabel("POS_EQ_RA") + "'");

  if( cmd.Found("col-dec",&a) )
    Write("COL_DEC = '" + a + "'");
  else if( cat.IsOk() )
    Write("COL_DEC = '" + cat.GetLabel("POS_EQ_DEC") + "'");

  if( cmd.Found("col-pm-ra",&a) )
    Write("COL_PMRA = '" + a + "'");
  else if( cat.IsOk() )
    Write("COL_PMRA = '" + cat.GetLabel("POS_EQ_PMRA") + "'");

  if( cmd.Found("col-pm-dec",&a) )
    Write("COL_PMDEC = '" + a + "'");
  else if( cat.IsOk() )
    Write("COL_PMDEC = '" + cat.GetLabel("POS_EQ_PMDEC") + "'");
*/


  if( cmd.Found("col-ra",&a) ) {
    apstr(a,n,c);
    Write("COL_RA = " + c);
  }

  if( cmd.Found("col-dec",&a) ) {
    apstr(a,n,c);
    Write("COL_DEC = " + c);
  }

  if( cmd.Found("col-pm-ra",&a) ) {
    apstr(a,n,c);
    Write("COL_PMRA = " + c);
  }

  if( cmd.Found("col-pm-dec",&a) ) {
    apstr(a,n,c);
    Write("COL_PMDEC = " + c);
  }


    /*
  if( cmd.Found("col-epoch-ra",&a) )
    Write("COL_EPRA = '" + a + "'");
  else if( cat.IsOk() )
    Write("COL_EPRA = '" + cat.GetLabel("TIME_EPOCH_RA") + "'");

  if( cmd.Found("col-epoch-dec",&a) )
    Write("COL_EPDEC = '" + a + "'");
  else if( cat.IsOk() )
    Write("COL_EPDEC = '" + cat.GetLabel("TIME_EPOCH_DEC") + "'");
    */

  if( cmd.Found("col-mag",&a) ) {
    apstr(a,n,c);
    Write("COL_MAG = " + c);
  }

  /*
  if( cmd.Found("col-mag",&a) )
    Write("COL_MAG = '" + a + "'");
  else if( cat.IsOk() )
    Write("COL_MAG = '" + cat.GetLabel("PHOT_MAG") + "'");
  */
  WriteFilesSingle(cmd.GetFiles());
}
