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
#include <wx/tokenzr.h>


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

    // the colour terms: the catalogue columns A,B of the colour A - B,
    // optionally with the columns of their errors (A,B,errA,errB); the option
    // is repeatable, every occurrence adds a colour (at most three)
    long ncterm = 0;
    for(wxCmdLineArgs::const_iterator arg = cmd.GetArguments().begin();
	arg != cmd.GetArguments().end(); ++arg) {
      if( arg->GetKind() == wxCMD_LINE_OPTION &&
	  arg->GetLongName() == "color-term" ) {
	long i;
	a = arg->GetStrVal();
	apstr(a,i,c);
	if( ! (i == 2 || i == 4) )
	  MuniFatalError("--color-term `%s': two columns A,B (the colour A - B), "
			 "or four A,B,errA,errB.",a);
	if( ++ncterm > 3 )
	  MuniFatalError("At most three colour terms (--color-term).");
	Write("COL_NCOLOR = %ld",i);
	Write("COL_COLOR = " + c);
      }
    }

    // the orders and the reference colours: one for all the colours,
    // or one for every colour (a list in the order of --color-term)
    if( cmd.Found("color-order",&a) ) {
      wxStringTokenizer tok(a,",");
      wxString l;
      long i = 0;
      while( tok.HasMoreTokens() ) {
	wxString b = tok.GetNextToken();
	if( ! (b.ToLong(&n) && (n == 1 || n == 2)) )
	  MuniFatalError("--color-order `%s': 1 (linear) or 2 (quadratic), "
			 "one for all the colours or one for every --color-term.",a);
	l += wxString::Format(" %ld",n);
	i++;
      }
      if( ! (i == 1 || i == ncterm) )
	MuniFatalError("--color-order `%s': one value, or one for every "
		       "--color-term (%ld).",a,ncterm);
      if( i > 1 )
	Write("COLOR_NORDER = %ld",i);
      Write("COLOR_ORDER =" + l);
    }

    if( cmd.Found("color-ref",&a) ) {
      wxStringTokenizer tok(a,",");
      wxString l;
      long i = 0;
      while( tok.HasMoreTokens() ) {
	wxString b = tok.GetNextToken();
	if( b.Upper() == "MEAN" )
	  l += " 1.7976931348623157e308";
	else if( b.ToDouble(&x) )
	  l += wxString::Format(" %.17g",x);
	else
	  MuniFatalError("Unknown --color-ref `%s' (mean, or a colour in mag; "
			 "one for all the colours or one for every --color-term).",a);
	i++;
      }
      if( ! (i == 1 || i == ncterm) )
	MuniFatalError("--color-ref `%s': one value, or one for every "
		       "--color-term (%ld).",a,ncterm);
      if( i > 1 )
	Write("COLOR_NREF = %ld",i);
      if( ! (i == 1 && a.Upper() == "MEAN") )
	Write("COLOR_REF =" + l);
    }

    if( cmd.Found("photsys-instr",&a) )
      Write("PHOTSYS_INSTR = '"+a+"'");

    if( cmd.Found("photsys-ref",&a) )
      Write("PHOTSYS_REF = '"+a+"'");

    if( cmd.Found("tol",&x) )
      Write("TOL = %.17g",x);

    if( cmd.Found("area",&x) )
      Write("AREA = %.17g",x);

    if( cmd.Found("saper",&n) )
      Write("SAPER = %ld",n);

    if( cmd.Found("photometry",&a) ) {
      wxString ph = a.Upper();
      if( ! (ph == "APER" || ph == "GROW" || ph == "PSF") )
	MuniFatalError("Unknown photometry `%s' (aper, grow, psf).",a);
      Write("PHOTOMETRY = '" + ph + "'");
    }

    if( cmd.Found("th",&x) || cmd.Found("threshold",&x))
      Write("THRESHOLD = %.17g",x);

    if( cmd.Found("e",&x) || cmd.Found("maxerr",&x))
      Write("MAXERR = %.17g",x);

    // aperture correction: the enclosed energy APEE (auto, by default),
    // none, or a constant factor (total / aperture flux)
    if( cmd.Found("apcorr",&a) ) {
      if( a.Upper() == "AUTO" )
	;
      else if( a.Upper() == "NONE" )
	Write("APEE = F");
      else if( a.ToDouble(&x) && x > 0 ) {
	Write("APEE = F");
	Write("APCORR = %.17g",x);
      }
      else
	MuniFatalError("Unknown --apcorr `%s' (auto, none, or a positive factor).",a);
    }

    // systematic error: auto (by residuals, default) or a value [mag]
    if( cmd.Found("syserr",&a) ) {
      if( a.Upper() == "AUTO" )
	Write("SYSERR = -1");
      else if( a.ToDouble(&x) && x >= 0 )
	Write("SYSERR = %.17g",x);
      else
	MuniFatalError("Unknown --syserr `%s' (auto, or a non-negative value in mag).",a);
    }

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

    WriteFilesMultiple(cmd.GetFiles());
  }
}
