/*

  Kombine

  Copyright © 2011 - 2026 F.Hroch (hroch@physics.muni.cz)

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


MuniKombine::MuniKombine(wxEvtHandler *han, const MuniCmdLineParser& cmd):
  MuniAction(han,"kombine")
{
  double x,y;
  long w,h;
  wxString a;

  CommonOutputSingle(cmd);
  CommonOptionsBitpix(cmd);

  if( cmd.Found("p",&a) || cmd.Found("projection",&a) )
    Write("PROJECTION = '" + a.Upper() + "'");

  if( cmd.Found("width",&w) && cmd.Found("height",&h) )
    Write("NAXES = %ld %ld",w,h);

  if( cmd.Found("xcen",&x) && cmd.Found("ycen",&y) )
    Write("CRPIX = %20.15f %20.15f",x,y);

  if( cmd.Found("rcen",&x) && cmd.Found("dcen",&y) )
    Write("CRVAL = %20.15f %20.15f",x,y);

  if( cmd.Found("pm-ra",&x) && cmd.Found("pm-dec",&y) )
    Write("CRMOV = %20.15f %20.15f",x,y);

  if( cmd.Found("pm-jdref",&x) )
    Write("JDREF = %20.15f",x);

  if( cmd.Found("scale",&x) )
    Write("SCALE = %25.15e",x);

  if( cmd.Found("angle",&x) )
    Write("ANGLE = %20.15f",x);

  if( cmd.Found("reflex",&a) ) {
    a.MakeUpper();
    wxString reflex = a.Find("Y") != wxNOT_FOUND ? "T" : "F";
    Write("REFLEX = " + reflex);
  }

  if( cmd.Found("i",&a) || cmd.Found("interpol",&a) )
    Write("INTERPOL = '" + a.Upper() + "'");

  if( cmd.Found("arith") )
    Write("ROBUST = F");

  if( cmd.Found("disable-back") )
    Write("BACKGROUND = F");

  if( cmd.Found("flux-conserve") )
    Write("FLUXCONSERVE = T");

  WriteFilesSingle(cmd.GetFiles());
}
