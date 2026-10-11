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
    if( ! (psf == "GAUSS" || psf == "MOFFAT" || psf == "FILE") )
      wxLogFatalError("Unknown PSF `%s' (GAUSS, MOFFAT, FILE).",a);
  }

  if( cmd.Found("psf-file",&a) ) {
    if( psf == "" )
      psf = "FILE";
    else if( psf != "FILE" )
      wxLogFatalError("--psf-file requires --psf=FILE.");
    Write("PSFFILE = '" + a + "'");
  }
  else if( psf == "FILE" )
    wxLogFatalError("--psf=FILE requires --psf-file.");

  if( psf != "" )
    Write("PSF = '" + psf + "'");

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
      wxLogFatalError("Unknown sky mode `%s' (fixed, local).",a);
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
	wxLogFatalError("Failed to interpret `%s' as a coordinate.",token);
      i++;
    }
    if( i != 2 )
      wxLogFatalError("--psf-centre requires x,y.");
    Write("PSFCENTRE = %.17g %.17g",c[0],c[1]);
  }

  if( cmd.Found("residual",&a) ) {
    if( cmd.GetFiles().GetCount() > 1 )
      wxLogFatalError("--residual is available for a single file only.");
    Write("RESIDUAL = '" + a + "'");
  }

  WriteFilesMultiple(cmd.GetFiles());
}
