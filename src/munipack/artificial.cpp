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

MuniArtificial::MuniArtificial(wxEvtHandler *han, const MuniCmdLineParser& cmd):
  MuniAction(han,"artificial")
{
#ifdef HAVE_CAF
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
    Write("HWHM = %lf",x);

  if( cmd.Found("eccentricity",&x) )
    Write("ECCENTRICITY = %lf",x);

  if( cmd.Found("inclination",&x) )
    Write("INCLINATION = %lf",x);

  if( cmd.Found("airy",&x) )
    Write("AIRY = %lf",x);

  if( cmd.Found("beta",&x) )
    Write("BETA = %lf",x);

  if( cmd.Found("maglim",&w) )
    Write("MAGLIM = %ld",w);

  if( cmd.Found("lc-model",&a) )
    Write("LCMODEL = '" + a.Upper() + "'");

  if( cmd.Found("lc-table",&a) )
    Write("LCTABLE = '" + a + "'");

  if( cmd.Found("lc-fourier",&a) )
    Write("LCFOURIER = '" + a + "'");

  if( cmd.Found("lc-mag",&x) )
    Write("LCMAG = %lf",x);

  if( cmd.Found("lc-amp",&x) )
    Write("LCAMP = %lf",x);

  if( cmd.Found("lc-jd0",&x) )
    Write("LCJD0 = %lf",x);

  if( cmd.Found("lc-per",&x) )
    Write("LCPER = %lf",x);

  if( cmd.Found("lc-ra",&x) && cmd.Found("lc-dec",&y) )
    Write("LCCOO = %20.15f %20.15f",x,y);

  if( cmd.Found("sky-mag",&x) )
    Write("SKYMAG = %lf",x);

  if( cmd.Found("sky-grad-x",&x) )
    Write("SKYGRADX = %lf",x);

  if( cmd.Found("sky-grad-y",&x) )
    Write("SKYGRADY = %lf",x);

  if( cmd.Found("area",&x) )
    Write("AREA = %lf",x);

  if( cmd.Found("diameter",&x) )
    Write("DIAMETER = %lf",x);

  if( cmd.Found("exptime",&x) )
    Write("EXPTIME = %lf",x);

  if( cmd.Found("qeff",&x) )
    Write("QEFF = %lf",x);

  if( cmd.Found("gain",&x) )
    Write("GAIN = %lf",x);

  if( cmd.Found("saturate",&x) )
    Write("SATURATE = %lf",x);

  if( cmd.Found("cosmics-freq",&x) )
    Write("COSMICS_FREQ = %lf",x);

  if( cmd.Found("cosmics-cmax",&x) )
    Write("COSMICS_CMAX = %lf",x);

  if( cmd.Found("cosmics-len",&x) )
    Write("COSMICS_LEN = %lf",x);

  if( cmd.Found("gnoise",&x) )
    Write("GNOISE = %lf",x);

  if( cmd.Found("extk",&x) )
    Write("EXTINK = %lf",x);

  if( cmd.Found("timestep",&x) )
    Write("TIMESTEP = %lf",x);

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
    Write("FOV = %lf",x);

  if( cmd.Found("p",&a) || cmd.Found("projection",&a) )
    Write("PROJECTION = '" + a.Upper() + "'");

  if( cmd.Found("f",&filter) || cmd.Found("filter",&filter) )
    Write("FILTER = '" + filter + "'");

  if( cmd.Found("width",&w) && cmd.Found("height",&h) )
    Write("NAXES = %ld %ld",w,h);

  if( cmd.Found("rcen",&x) && cmd.Found("dcen",&y) )
    Write("CRVAL = %20.15f %20.15f",x,y);

  if( cmd.Found("scale",&x) )
    Write("SCALE = %25.15e",x);

  if( cmd.Found("angle",&x) )
    Write("ANGLE = %20.15f",x);

  if( cmd.Found("long",&x) && cmd.Found("lat",&y) )
    Write("GEOGRAPHIC = %20.15f %20.15f",x,y);

  if( cmd.Found("photsys",&a) )
    Write("PHOTSYS = '"+a+"'");

  if( cmd.Found("atmosphere") )
    Write("ATMOSPHERE = T");

  if( cmd.Found("mask",&a) )
    Write("OUTPUT = '"+a+"'");

  if( cmd.Found("disable-noise") )
    Write("NOISE = F");

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
