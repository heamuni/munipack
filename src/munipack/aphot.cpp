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
    Write("ECCENTRICITY = %.17g",x);

  if( cmd.Found("inclination",&x) )
    Write("INCLINATION = %.17g",x);

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
	MuniFatalError("Failed to interpret `%s' as a coordinate.",token);

      token = tokenizer.GetNextToken();
      if( ! token.IsEmpty() && token.ToDouble(&q) )
	y.push_back(q);
      else
	MuniFatalError("Failed to interpret `%s' as a coordinate.",token);
    }

    Write("NSTARS = %d",x.size());
    for(size_t i = 0; i < x.size(); i++)
      Write("STAR = %.15f %.15f",x[i],y[i]);
  }

  // the enclosed energy of apertures
  wxString psf;
  if( cmd.Found("psf",&a) ) {
    psf = a.Upper();
    if( ! (psf == "GAUSS" || psf == "MOFFAT" || psf == "FILE") )
      MuniFatalError("Unknown PSF `%s' (GAUSS, MOFFAT, FILE).",a);
  }

  if( cmd.Found("psf-file",&a) ) {
    if( ! (psf == "" || psf == "FILE") )
      MuniFatalError("--psf-file requires --psf=FILE.");
    psf = "FILE";
    Write("PSFFILE = '" + a + "'");
  }
  else if( psf == "FILE" )
    MuniFatalError("--psf=FILE requires --psf-file.");

  if( psf != "" )
    Write("PSF = '" + psf + "'");

  if( cmd.Found("fwhm",&x) ) {
    if( ! (x > 0) )
      MuniFatalError("--fwhm must be positive.");
    Write("FWHM = %.17g",x);
  }

  if( cmd.Found("beta",&x) ) {
    if( ! (x > 1) )
      MuniFatalError("--beta must be greater than 1.");
    Write("BETA = %.17g",x);
  }

  if( cmd.Found("scale",&x) )
    Write("SCALE = %.17g",x);

  if( cmd.Found("psf-centre",&a) ) {
    wxStringTokenizer t(a,sep);
    double c[2];
    int i = 0;
    while( t.HasMoreTokens() && i < 2 ) {
      wxString token(t.GetNextToken());
      if( ! token.ToDouble(&c[i]) )
	MuniFatalError("Failed to interpret `%s' as a coordinate.",token);
      i++;
    }
    if( i != 2 )
      MuniFatalError("--psf-centre requires x,y.");
    Write("PSFCENTRE = %.17g %.17g",c[0],c[1]);
  }

  if( cmd.Found("apcorr",&a) ) {
    if( a.Upper() != "AUTO" )
      MuniFatalError("Unknown --apcorr `%s' (auto).",a);
    if( psf != "" )
      MuniFatalError("--apcorr=auto and --psf, --psf-file are exclusive.");
    Write("APCORR = 'AUTO'");
  }

  if( cmd.Found("apcorr-var",&n) ) {
    if( n < 0 || n > 2 )
      MuniFatalError("--apcorr-var is 0, 1 or 2.");
    if( ! cmd.Found("apcorr") )
      MuniFatalError("--apcorr-var requires --apcorr=auto.");
    Write("APCORRVAR = %ld",n);
  }

  if( cmd.Found("apcorr-radial") ) {
    if( ! cmd.Found("apcorr") )
      MuniFatalError("--apcorr-radial requires --apcorr=auto.");
    Write("APCORRRADIAL = T");
  }

  if( cmd.Found("coo") && (psf != "" || cmd.Found("apcorr")) )
    MuniFatalError("--coo does not support the enclosed energy (--psf, --apcorr).");

  WriteFilesMultiple(cmd.GetFiles());
}
