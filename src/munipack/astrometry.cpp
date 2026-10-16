/*

  Astrometry

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

  1. create list to match:
     * projection from catalogue
     * photometry of ref
  2. for every image:
     * photometry
     * match
     * astrometry


   TODO:

     * print default parameters (in help ?)

*/

#include "munipack.h"
#include "vocatconf.h"
#include <wx/wx.h>
#include <wx/tokenzr.h>

MuniAstrometry::MuniAstrometry(wxEvtHandler *h, const MuniCmdLineParser& cmd):
    MuniAction(h,"astrometry")
{
#ifdef HAVE_CAF
  SetCAF(cmd);
#endif

  CommonOutputMultiple(cmd);

  if( cmd.Found("show-defaults") ) {
    Write("MODE = 'DEFAULTS'");
    return;
  }

  if( cmd.Found("remove") ) {
    Write("MODE = 'REMOVE'");
  }

  else {

    wxString a, cone;
    double x,y;
    long i;
    size_t nseq = 0;
    bool mancal = false;
    bool seqcal = false;
    bool crosscal = false;

    if( cmd.Found("m",&a) || cmd.Found("mode",&a) ) {
      a = a.Upper();
      Write("MODE = '" + a + "'");
      mancal = a == "MANUAL";
      seqcal = a == "SEQUENCE";
      crosscal = a == "CROSS";
    }

    // PROJ ** must ** precede both REF and CAT
    if( cmd.Found("p",&a) || cmd.Found("projection",&a) )
      Write("PROJECTION = '" + a.Upper() + "'");

    if( cmd.Found("fit",&a) )
      Write("FIT = '" + a.Upper() + "'");

    if( cmd.Found("report") )
      Write("REPORT = T");

    if( cmd.Found("xcen",&x) && cmd.Found("ycen",&y) )
      Write("CRPIX = %20.15f %20.15f",x,y);

    if( cmd.Found("rcen",&x) && cmd.Found("dcen",&y) )
      Write("CRVAL = %20.15f %20.15f",x,y);

    if( cmd.Found("scale",&x) )
      Write("SCALE = %25.15e",x);

    if( cmd.Found("angle",&x) )
      Write("ANGLE = %20.15f",x);

    if( cmd.Found("scale-tol",&x) )
      Write("SCALETOL = %e",x);

    if( cmd.Found("scale-free") )
      Write("SCALETOL = 0");

    if( cmd.Found("blind") )
      Write("WCSGUIDE = F");

    if( cmd.Found("wcs-offset",&a) ) {
      // a number with optional units: pix (default), arcsec ("), arcmin ('),
      // deg, % (of the larger side of the frame)
      wxString num = a.Strip(wxString::both), unit;
      size_t n = 0;
      while( n < num.Len() && (wxIsdigit(num[n]) || num[n] == '.' || num[n] == '+' ||
			       num[n] == 'e' || num[n] == 'E' ||
			       (n > 0 && num[n] == '-' && (num[n-1] == 'e' || num[n-1] == 'E'))) )
	n++;
      unit = num.Mid(n).Strip(wxString::both).Lower();
      num = num.Left(n);
      if( unit == "" || unit == "pix" || unit == "px" )
	unit = "PIX";
      else if( unit == "arcsec" || unit == "as" || unit == "\"" )
	unit = "ARCSEC";
      else if( unit == "arcmin" || unit == "am" || unit == "'" )
	unit = "ARCMIN";
      else if( unit == "deg" || unit == "d" )
	unit = "DEG";
      else if( unit == "%" )
	unit = "PERCENT";
      else
	MuniFatalError("--wcs-offset " + a + ": the units must be pix, arcsec, arcmin, deg or %.");
      if( ! num.ToCDouble(&x) || ! (x > 0) )
	MuniFatalError("--wcs-offset " + a + ": a positive number (with units) is expected.");
      Write("WCSOFFSET = %e",x);
      Write("WCSOFFUNIT = '" + unit + "'");
    }

    if( cmd.Found("wcs-rot",&x) )
      Write("WCSROT = %e",x);

    if( cmd.Found("roll-unknown") && cmd.Found("roll-known") )
      MuniFatalError("--roll-unknown and --roll-known exclude each other.");

    if( cmd.Found("roll-unknown") )
      Write("ROLL = 'UNKNOWN'");

    if( cmd.Found("roll-known") )
      Write("ROLL = 'KNOWN'");

    if( cmd.Found("reflex-unknown") )
      Write("REFLBOTH = T");

    if( cmd.Found("wcs-scale",&x) )
      Write("WCSSCALE = %e",x);

    if( cmd.Found("reflex") )
      Write("REFLEX = T");

    if( cmd.Found("radesys",&a) )
      Write("RADESYS = '" + a + "'");

    if( cmd.Found("epoch",&x) )
      Write("EPOCH = %e",x);

    if( cmd.Found("sig",&x) )
      Write("SIG = %e",x);

    if( cmd.Found("sigcat",&x) )
      Write("SIGCAT = %e",x);

    if( cmd.Found("flxmin",&x) )
      Write("FLXMIN = %e",x);

    if( cmd.Found("flxmax",&x) )
      Write("FLXMAX = %e",x);

    if( cmd.Found("units",&a) && ( a == "uas" || a == "mas" || a == "arcsec" ||
				   a == "arcmin" || a == "deg") )
      Write("AUNITS = '" + a.Lower() + "'");

    if( cmd.Found("rcen",&x) ||  cmd.Found("dcen",&x) ||
	cmd.Found("scale",&x) || cmd.Found("angle",&x) )
      Write("INITPAR = F");

    if( cmd.Found("minmatch",&i) ) {
      // the solution needs five stars (four parameters and their errors)
      if( i < 5 )
	MuniFatalError("--minmatch %ld: the match sequence needs five stars at least.",i);
      Write("MINMATCH = %ld",i);
    }

    if( cmd.Found("maxmatch",&i) )
      Write("MAXMATCH = %ld",i);

    if( cmd.Found("maxstars",&i) )
      Write("MAXSTARS = %ld",i);

    if( cmd.Found("luckymatch",&i) )
      Write("LUCKYMATCH = %ld",i);

    if( cmd.Found("disable-lucky-match") )
      Write("LUCKYMATCH = 0");

    if( cmd.Found("disable-flux-check") )
      Write("FLUXCHECK = F");

    if( cmd.Found("disable-rms-check") )
      Write("RMSCHECK = F");

    if( cmd.Found("disable-save") )
      Write("WCSSAVE = F");

    if( cmd.Found("r",&a) || cmd.Found("ref",&a) ) {
      Write("REF = '" + a + "'");
      cone = "";
    }

    if( cmd.Found("R",&a) || cmd.Found("rel",&a) ) {
      Write("REL = '" + a + "'");
      cone = "";
    }

    if( cmd.Found("c",&a) || cmd.Found("cat",&a) ) {
      cone = a;
    }

    if( cone != "" && ! (mancal && seqcal && crosscal) ) {

      wxString a,vocat;

      if( cmd.Found("vocat",&a) )
	vocat = a;
      else if( FindDataDir("VOcat_conf.xml",&a) )
	vocat = a;

      VOCatConf catalogs(vocat);
      if( ! catalogs.IsOk() )
	MuniFatalError("Failed to read the definition file of catalogues `"
		       + vocat + "'.");
      VOCatResources cat(catalogs.GetCatFits(cone));

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

      if( cmd.Found("col-epoch-ra",&a) )
	Write("COL_EPRA = '" + a + "'");
      else if( cat.IsOk() )
	Write("COL_EPRA = '" + cat.GetLabel("TIME_EPOCH_RA") + "'");

      if( cmd.Found("col-epoch-dec",&a) )
	Write("COL_EPDEC = '" + a + "'");
      else if( cat.IsOk() )
	Write("COL_EPDEC = '" + cat.GetLabel("TIME_EPOCH_DEC") + "'");

      if( cmd.Found("col-mag",&a) )
	Write("COL_MAG = '" + a + "'");
      else if( cat.IsOk() )
	Write("COL_MAG = '" + cat.GetLabel("PHOT_MAG") + "'");
      //	Write("COL_MAG = '" + cat.GetSort() + "'");

      Write("CAT = '" + cone + "'");
    }

    if( cmd.Found("seq1",&a) ) {
      wxStringTokenizer t(a,",");
      if( nseq == 0 ) {
	nseq = t.CountTokens();
	Write("NSEQ = %ld",nseq);
      }
      if( nseq != t.CountTokens() )
	MuniFatalError("Specified sequences are unequal in length.");

      Write("SEQ1 = '" + a + "'");
    }

    if( cmd.Found("seq2",&a) ) {
      wxStringTokenizer t(a,",");
      if( nseq == 0 ) {
	nseq = t.CountTokens();
	Write("NSEQ = %ld",nseq);
      }
      if( nseq != t.CountTokens() )
	MuniFatalError("Specified sequences are unequal in length.");

      Write("SEQ2 = '" + a + "'");
    }

  }

  WriteFilesMultiple(cmd.GetFiles());


  // catalogue or image ?
  // if catalogue
  //     ra,dec -> vo (UCAC)
  // photometry
  // astrometry



  // catalogue or image ?
  // if catalogue
  //    ra,dec -> vo (UCAC) -> projection -> table
  // if image
  //    list table

  // for every image
  //    aperture photometry
  //    match to reference table
  //    if catalogue
  //       astrometry
}
