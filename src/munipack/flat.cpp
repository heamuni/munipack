/*

  Average of flats.

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

MuniFlat::MuniFlat(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"flat")
{
  CommonOutputSingle(cmd);
  CommonOptionsBitpix(cmd);
  CommonOptionsPhCorr(cmd);
  EnvironmentPhCorr();

  double x;
  wxString a;

  if( cmd.Found("bias",&a) )
    Write("BIAS = '"+a+"'");

  if( cmd.Found("dark",&a) )
    Write("DARK = '"+a+"'");

  if( cmd.Found("xdark",&x) )
    Write("XDARK = %lf",x);

  if( cmd.Found("gain",&x) )
    Write("GAIN = %lf",x);

  if( cmd.Found("approximation",&a) )
    Write("APPROXIMATION = '" + a.Upper() + "'");

  WriteFilesSingle(cmd.GetFiles());
}
