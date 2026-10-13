/*

  Artificial frames

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
#include "vocatconf.h"
#include <wx/wx.h>
#include <wx/filename.h>

using namespace std;

// both options of a pair (like --rcen, --dcen) must be given together
static void CheckPair(const MuniCmdLineParser& cmd, const wxString& a,
		      const wxString& b)
{
  if( cmd.Found(a) != cmd.Found(b) )
    wxLogFatalError("Both --"+a+" and --"+b+" must be given together.");
}

MuniArtificial::MuniArtificial(wxEvtHandler *han, const MuniCmdLineParser& cmd):
  MuniAction(han,"artificial")
{
#ifdef HAVE_CAF
  long np;
  if( cmd.Found("mproc",&np) && np <= 0 )
    wxLogFatalError("The number of parallel processes (--mproc) must be positive.");
  SetCAF(cmd);
#endif

  CommonOptionsBitpix(cmd);

  double x,y;
  long w,h,n;
  wxString a, cname, filter;

  if( cmd.Found("verbose") )
    Write("VERBOSE = T");

  if( cmd.Found("phsystab",&a) )
    Write("PHSYSTABLE = '" + a + "'");
  else if( FindDataDir("photosystems.fits",&a) )
    Write("PHSYSTABLE = '" + a + "'");

  if( cmd.Found("psf",&a) )
    Write("PROFILE = '" + a.Upper() + "'");

  if( cmd.Found("psf-file",&a) )
    Write("PSFFILE = '" + a + "'");

  if( cmd.Found("fits-keys-file",&a) )
    Write("FITSKEYFILE = '" + a + "'");

  if( cmd.Found("flat-model-file",&a) )
    Write("FLATMODELFILE = '" + a + "'");

  if( cmd.Found("spread",&a) )
    Write("SPREAD = '" + a.Upper() + "'");

  if( cmd.Found("hwhm",&x) )
    Write("HWHM = %.17g",x);

  if( cmd.Found("eccentricity",&x) )
    Write("ECCENTRICITY = %.17g",x);

  if( cmd.Found("inclination",&x) )
    Write("INCLINATION = %.17g",x);

  if( cmd.Found("airy",&x) )
    Write("AIRY = %.17g",x);

  if( cmd.Found("beta",&x) )
    Write("BETA = %.17g",x);

  if( cmd.Found("maglim",&w) )
    Write("MAGLIM = %ld",w);

  if( cmd.Found("lc-model",&a) )
    Write("LCMODEL = '" + a.Upper() + "'");

  if( cmd.Found("lc-table",&a) )
    Write("LCTABLE = '" + a + "'");

  if( cmd.Found("lc-fourier",&a) )
    Write("LCFOURIER = '" + a + "'");

  if( cmd.Found("lc-mag",&x) )
    Write("LCMAG = %.17g",x);

  if( cmd.Found("lc-amp",&x) )
    Write("LCAMP = %.17g",x);

  if( cmd.Found("lc-jd0",&x) )
    Write("LCJD0 = %.17g",x);

  if( cmd.Found("lc-per",&x) )
    Write("LCPER = %.17g",x);

  CheckPair(cmd,"lc-ra","lc-dec");
  if( cmd.Found("lc-ra",&x) && cmd.Found("lc-dec",&y) )
    Write("LCCOO = %.17g %.17g",x,y);

  if( cmd.Found("sky-mag",&x) )
    Write("SKYMAG = %.17g",x);

  if( cmd.Found("sky-grad-x",&x) )
    Write("SKYGRADX = %.17g",x);

  if( cmd.Found("sky-grad-y",&x) )
    Write("SKYGRADY = %.17g",x);

  if( cmd.Found("area",&x) )
    Write("AREA = %.17g",x);

  if( cmd.Found("diameter",&x) )
    Write("DIAMETER = %.17g",x);

  if( cmd.Found("exptime",&x) )
    Write("EXPTIME = %.17g",x);

  if( cmd.Found("qeff",&x) )
    Write("QEFF = %.17g",x);

  if( cmd.Found("gain",&x) )
    Write("GAIN = %.17g",x);

  if( cmd.Found("saturate",&x) )
    Write("SATURATE = %.17g",x);

  if( cmd.Found("cosmics-freq",&x) )
    Write("COSMICS_FREQ = %.17g",x);

  if( cmd.Found("cosmics-cmax",&x) )
    Write("COSMICS_CMAX = %.17g",x);

  if( cmd.Found("cosmics-len",&x) )
    Write("COSMICS_LEN = %.17g",x);

  if( cmd.Found("gnoise",&x) )
    Write("GNOISE = %.17g",x);

  if( cmd.Found("dark",&x) ) {
    if( x < 0 )
      wxLogFatalError("Dark current (--dark) must be non-negative.");
    Write("DARK = %.17g",x);
  }

  if( cmd.Found("readnoise",&x) ) {
    if( x < 0 )
      wxLogFatalError("Read-out noise (--readnoise) must be non-negative.");
    Write("READNOISE = %.17g",x);
  }

  if( cmd.Found("bias",&x) )
    Write("BIAS = %.17g",x);

  if( cmd.Found("seed",&n) )
    Write("SEED = %ld",n);

  if( cmd.Found("extk",&x) )
    Write("EXTINK = %.17g",x);

  if( cmd.Found("timestep",&x) )
    Write("TIMESTEP = %.17g",x);

  if( cmd.Found("count",&n) )
    Write("NFILES = %ld",n);

  if( cmd.Found("date",&a) ) {
    a.Replace("-"," ");
    Write("DATE = " + a);
  }

  if( cmd.Found("time",&a) ) {
    a.Replace(":"," ");
    Write("TIME = " + a);
  }

  if( cmd.Found("fov",&x) )
    Write("FOV = %.17g",x);

  if( cmd.Found("p",&a) || cmd.Found("projection",&a) )
    Write("PROJECTION = '" + a.Upper() + "'");

  if( cmd.Found("f",&filter) || cmd.Found("filter",&filter) )
    Write("FILTER = '" + filter + "'");

  CheckPair(cmd,"width","height");
  if( cmd.Found("width",&w) && cmd.Found("height",&h) )
    Write("NAXES = %ld %ld",w,h);

  CheckPair(cmd,"rcen","dcen");
  if( cmd.Found("rcen",&x) && cmd.Found("dcen",&y) )
    Write("CRVAL = %.17g %.17g",x,y);

  if( cmd.Found("scale",&x) )
    Write("SCALE = %.17g",x);

  if( cmd.Found("angle",&x) )
    Write("ANGLE = %.17g",x);

  CheckPair(cmd,"long","lat");
  if( cmd.Found("long",&x) && cmd.Found("lat",&y) )
    Write("GEOGRAPHIC = %.17g %.17g",x,y);

  if( cmd.Found("photsys",&a) )
    Write("PHOTSYS = '"+a+"'");

  if( cmd.Found("atmosphere") )
    Write("ATMOSPHERE = T");

  if( cmd.Found("mask",&a) )
    Write("OUTPUT = '"+a+"'");

  if( cmd.Found("disable-noise") )
    Write("NOISE = F");

  if( cmd.Found("no-watermark") )
    Write("WATERMARK = F");

  if( cmd.Found("c",&cname) || cmd.Found("cat",&cname) ) {

    wxString vocat;
    if( cmd.Found("vocat",&a) )
      vocat = a;
    else if( FindDataDir("VOcat_conf.xml",&a) )
      vocat = a;

    VOCatConf catalogs(vocat);
    VOCatResources cat(catalogs.GetCatFits(cname));

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
    else if( filter != "" )
      Write("COL_MAG = '" + filter + "mag'");
    else if( cat.IsOk() )
      Write("COL_MAG = '" + cat.GetSort() + "'");

    Write("CAT = '" + cname + "'");
  }

  if( cmd.Found("galcat",&cname) )
    Write("GALCAT = '" + cname + "'");
}
