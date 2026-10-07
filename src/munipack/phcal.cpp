/*

  Photometric calibration

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
#include "vocatconf.h"
#include <wx/wx.h>
#include <wx/filename.h>


MuniPhcal::MuniPhcal(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"phcal")
{
  wxString vocat;

  CommonOutputMultiple(cmd);

  wxString a;
  if( !  cmd.Found("O") &&
      ! (cmd.Found("target-directory",&a) || cmd.Found("t",&a)) ) {
    SetAdvanced(true);
    SetMask("\\1_cal.\\2");
  }

  if( cmd.Found("phsystab",&a) )
    Write("PHSYSTABLE = '" + a + "'");
  else if( FindDataDir("photosystems.fits",&a) )
    Write("PHSYSTABLE = '" + a + "'");

  if( cmd.Found("vocat",&a) )
    vocat = a;
  else if( FindDataDir("VOcat_conf.xml",&a) )
    vocat = a;

  if( cmd.Found("list") )

    Write("LIST = T");

  else {

    long n;
    double x;
    wxString a,c;

    if( cmd.Found("advanced") )
      Write("ADVANCED = T");

    if( cmd.Found("f",&a) || cmd.Found("filters",&a) ) {
      long i;
      apstr(a,i,c);
      Write("NFILTERS = %ld",i);
      Write("FILTERS = " + c);
    }

    if( cmd.Found("col-mag",&a) ) {
      long i;
      apstr(a,i,c);
      Write("COL_NMAG = %ld",i);
      Write("COL_MAG = " + c);
    }

    if( cmd.Found("col-magerr",&a) ) {
      long i;
      apstr(a,i,c);
      Write("COL_NMAGERR = %ld",i);
      Write("COL_MAGERR = " + c);
    }

    if( cmd.Found("q",&a) || cmd.Found("quantity",&a) ) {
      long i;
      apstr(a,i,c);
      Write("NQUANTITIES = %ld",i);
      Write("QUANTITIES = " + c);
    }

    if( cmd.Found("tratab",&a) )
      Write("TRATABLE = '" + a + "'");

    if( cmd.Found("photsys-instr",&a) )
      Write("PHOTSYS_INSTR = '"+a+"'");

    if( cmd.Found("photsys-ref",&a) )
      Write("PHOTSYS_REF = '"+a+"'");

    if( cmd.Found("tol",&x) )
      Write("TOL = %e",x);

    if( cmd.Found("area",&x) )
      Write("AREA = %e",x);

    if( cmd.Found("saper",&n) )
      Write("SAPER = %ld",n);

    if( cmd.Found("th",&x) || cmd.Found("threshold",&x))
      Write("THRESHOLD = %e",x);

    if( cmd.Found("e",&x) || cmd.Found("maxerr",&x))
      Write("MAXERR = %e",x);

    if( cmd.Found("apcorr",&x))
      Write("APCORR = %e",x);

    if( cmd.Found("C",&a) || cmd.Found("cal",&a) ) {
      long n;
      apstr(a,n,c);
      Write("NCTPH = %ld",n);
      Write("CTPH = " + a);
    }

    if( cmd.Found("r",&a) || cmd.Found("ref",&a) )
      Write("REF = '" + a + "'");

    if( cmd.Found("c",&c) || cmd.Found("cat",&c) ) {

      VOCatConf catalogs(vocat);
      VOCatResources cat(catalogs.GetCatFits(c));

      if( cmd.Found("col-ra",&a) )
	Write("COL_RA = '" + a + "'");
      else if( cat.IsOk() )
	Write("COL_RA = '" + cat.GetLabel("POS_EQ_RA") + "'");

      if( cmd.Found("col-dec",&a) )
	Write("COL_DEC = '" + a + "'");
      else if( cat.IsOk() )
	Write("COL_DEC = '" + cat.GetLabel("POS_EQ_DEC") + "'");

      Write("CAT = '" + c + "'");
    }

    WriteFilesMultiple(cmd.GetFiles());
  }
}
