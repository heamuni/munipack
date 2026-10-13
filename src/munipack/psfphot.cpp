/*

  PSF Photometry

  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)

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

using namespace std;


MuniPsfphot::MuniPsfphot(wxEvtHandler *h, const MuniCmdLineParser& cmd):
    MuniAction(h,"psfphot")
{
  // determine character representing the decimal point
  const wxString decimal = wxLocale::GetInfo(wxLOCALE_DECIMAL_POINT,
					     wxLOCALE_CAT_NUMBER);
  const wxString sep = decimal == "." ? ",;" : ";";

  CommonOutputMultiple(cmd);

  wxString a, psf;
  double x;
  long n;

  if( cmd.Found("psf",&a) ) {
    psf = a.Upper();
    if( ! (psf == "GAUSS" || psf == "MOFFAT" || psf == "FILE" ||
	   psf == "EMPIRICAL") )
      MuniFatalError("Unknown PSF `%s' (GAUSS, MOFFAT, FILE, EMPIRICAL).",a);
  }

  if( cmd.Found("psf-file",&a) ) {
    if( psf == "" )
      psf = "FILE";
    else if( ! (psf == "FILE" || psf == "EMPIRICAL") )
      MuniFatalError("--psf-file requires --psf=FILE or EMPIRICAL.");
    Write("PSFFILE = '" + a + "'");
  }
  else if( psf == "FILE" )
    MuniFatalError("--psf=FILE requires --psf-file.");

  if( psf != "" )
    Write("PSF = '" + psf + "'");

  // the empirical PSF
  const char *eopts[] = { "psf-core", "psf-zoom", "psf-radius", "psf-var",
			  "psf-stars", "psf-list", "psf-save", "psf-smooth" };
  for(size_t i = 0; i < sizeof(eopts)/sizeof(eopts[0]); i++)
    if( cmd.Found(eopts[i]) && psf != "EMPIRICAL" )
      MuniFatalError("--%s requires --psf=EMPIRICAL.",eopts[i]);

  if( cmd.Found("psf-core",&a) ) {
    wxString core = a.Upper();
    if( ! (core == "GAUSS" || core == "MOFFAT") )
      MuniFatalError("Unknown core `%s' of the empirical PSF (GAUSS, MOFFAT).",a);
    Write("PSFCORE = '" + core + "'");
  }

  if( cmd.Found("psf-zoom",&n) ) {
    if( n < 1 )
      MuniFatalError("--psf-zoom must be positive.");
    Write("PSFZOOM = %ld",n);
  }

  if( cmd.Found("psf-radius",&x) )
    Write("PSFRADIUS = %.17g",x);

  if( cmd.Found("psf-var",&n) ) {
    if( n < 0 || n > 2 )
      MuniFatalError("--psf-var is 0, 1 or 2.");
    Write("PSFVAR = %ld",n);
  }

  if( cmd.Found("psf-stars",&n) ) {
    if( n < 1 )
      MuniFatalError("--psf-stars must be positive.");
    Write("PSFSTARS = %ld",n);
  }

  if( cmd.Found("psf-smooth",&x) ) {
    if( x < 0 )
      MuniFatalError("--psf-smooth must be non-negative.");
    Write("PSFSMOOTH = %.17g",x);
  }

  if( cmd.Found("psf-list",&a) )
    Write("PSFLIST = '" + a + "'");

  if( cmd.Found("psf-save",&a) ) {
    if( cmd.GetFiles().GetCount() > 1 )
      MuniFatalError("--psf-save is available for a single file only.");
    Write("PSFSAVE = '" + a + "'");
  }

  if( cmd.Found("group",&x) ) {
    if( x < 0 )
      MuniFatalError("--group must be non-negative.");
    Write("GROUP = %.17g",x);
  }

  if( cmd.Found("residual-detect",&n) ) {
    if( n < 0 )
      MuniFatalError("--residual-detect must be non-negative.");
    Write("RESDETECT = %ld",n);
  }

  if( cmd.Found("residual-threshold",&x) ) {
    if( ! (x > 0) )
      wxLogFatalError("--residual-threshold must be positive.");
    Write("RESTHRESH = %.17g",x);
  }

  if( cmd.Found("fwhm",&x) )
    Write("FWHM = %.17g",x);

  if( cmd.Found("beta",&x) )
    Write("BETA = %.17g",x);

  if( cmd.Found("elliptic") )
    Write("ELLIPTIC = T");

  if( cmd.Found("fitrad",&x) )
    Write("FITRAD = %.17g",x);

  if( cmd.Found("sky",&a) ) {
    wxString sky = a.Upper();
    if( ! (sky == "FIXED" || sky == "LOCAL") )
      MuniFatalError("Unknown sky mode `%s' (fixed, local).",a);
    Write("SKY = '" + sky + "'");
  }

  if( cmd.Found("maxiter",&n) )
    Write("MAXITER = %ld",n);

  if( cmd.Found("npass",&n) )
    Write("NPASS = %ld",n);

  if( cmd.Found("saturate",&x) )
    Write("SATURATE = %.17g",x);

  if( cmd.Found("gain",&x) )
    Write("GAIN = %.17g",x);

  if( cmd.Found("read-noise",&x) )
    Write("READNOISE = %.17g",x);

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

  if( cmd.Found("residual",&a) ) {
    if( cmd.GetFiles().GetCount() > 1 )
      MuniFatalError("--residual is available for a single file only.");
    Write("RESIDUAL = '" + a + "'");
  }

  WriteFilesMultiple(cmd.GetFiles());
}
