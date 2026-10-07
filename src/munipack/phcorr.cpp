/*

  Photometric corrections

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

MuniPhcorr::MuniPhcorr(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"phcorr")
{
  wxString a;

  CommonOptionsBitpix(cmd);
  CommonOutputMultiple(cmd);
  CommonOptionsEnableOverwrite(cmd);
  EnvironmentPhCorr();

  if( cmd.Found("normalise-flat") )
    Write("NORMALISE = T");

  if( cmd.Found("flat",&a) )
    Write("FLAT = '"+a+"'");

  if( cmd.Found("bias",&a) )
    Write("BIAS = '"+a+"'");

  if( cmd.Found("dark",&a) )
    Write("DARK = '"+a+"'");

  if( cmd.Found("bitmask",&a) )
    Write("MASK = '"+a+"'");

  if( cmd.Found("xbitmask",&a) )
    Write("XMASK = '"+a.Upper()+"'");

  long i;
  if( cmd.Found("box",&i) )
    Write("BOX = %ld",i);

  double x;
  if( cmd.Found("xdark",&x) )
    Write("XDARK = %lf",x);

  if( cmd.Found("gain",&x) )
    Write("GAIN = %lf",x);

  WriteFilesMultiple(cmd.GetFiles());
}
