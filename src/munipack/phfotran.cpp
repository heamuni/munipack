/*

  Photometric system transformation

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
#include "vocatconf.h"
#include <wx/wx.h>
#include <wx/filename.h>


MuniPhfotran::MuniPhfotran(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"phfotran")
{
  CommonOutputSingle(cmd);

  double x;
  wxString a,c,vocat;

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

    if( cmd.Found("f",&a) || cmd.Found("filters",&a) ) {
      long i;
      apstr(a,i,c);
      Write("NFILTERS = %ld",i);
      Write("FILTERS = " + c);
    }

    if( cmd.Found("C",&a) || cmd.Found("cal",&a) ) {
      long n;
      apstr(a,n,c);
      Write("NCTPH = %ld",n);
      Write("CTPH = " + a);
    }

    if( cmd.Found("E",&a) || cmd.Found("extin",&a) ) {
      long n;
      apstr(a,n,c);
      Write("NEXTIN = %ld",n);
      Write("EXTIN = " + a);
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

    // the model: mag (colour terms of magnitudes, default) or flux
    // (the linear mixing of photon rates)
    if( cmd.Found("model",&a) ) {
      wxString m = a.Upper();
      if( ! (m == "MAG" || m == "FLUX") )
	MuniFatalError("Unknown --model `%s' (mag, flux).",a);
      Write("MODEL = '" + m + "'");
    }

    long n;
    if( cmd.Found("color-order",&n) ) {
      if( ! (n == 1 || n == 2) )
	MuniFatalError("--color-order is 1 (linear) or 2 (quadratic).");
      Write("COLOR_ORDER = %ld",n);
    }

    if( cmd.Found("photsys-instr",&a) )
      Write("PHOTSYS_INSTR = '"+a+"'");

    if( cmd.Found("photsys-ref",&a) )
      Write("PHOTSYS_REF = '"+a+"'");

    if( cmd.Found("tol",&x) )
      Write("TOL = %.17g",x);

    if( cmd.Found("area",&x) )
      Write("AREA = %.17g",x);

    if( cmd.Found("c",&c) || cmd.Found("cat",&c) ) {

      VOCatConf catalogs(vocat);
      if( ! catalogs.IsOk() )
	MuniFatalError("Failed to read the definition file of catalogues `"
		       + vocat + "'.");
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

    WriteFilesSingle(cmd.GetFiles());
  }
}
