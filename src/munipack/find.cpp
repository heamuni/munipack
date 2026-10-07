/*

  Find Stars

  Copyright © 2013 - 2026 F.Hroch (hroch@physics.muni.cz)

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

MuniFind::MuniFind(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"find")
{
  CommonOutputMultiple(cmd);

  double x;

  if( cmd.Found("f",&x) || cmd.Found("fwhm",&x) )
    Write("FWHM = %lf",x);

  if( cmd.Found("th",&x) || cmd.Found("threshold",&x) )
    Write("THRESHOLD = %lf",x);

  if( cmd.Found("read-noise",&x) )
    Write("READNOISE = %lf",x);

  if( cmd.Found("saturate",&x) )
    Write("SATURATE = %lf",x);

  if( cmd.Found("lothresh",&x) )
    Write("LOWER_THRESHOLD = %lf",x);

  if( cmd.Found("shrplo",&x) )
    Write("SHARP_LOWER = %lf",x);

  if( cmd.Found("shrphi",&x) )
    Write("SHARP_HIGHER = %lf",x);

  if( cmd.Found("rndlo",&x) )
    Write("ROUND_LOWER = %lf",x);

  if( cmd.Found("rndhi",&x) )
    Write("ROUND_HIGHER = %lf",x);

  WriteFilesMultiple(cmd.GetFiles());
}
