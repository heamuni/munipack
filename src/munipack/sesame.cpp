/*

  Sesame name resolver

  Copyright © 2019 - 2026 F.Hroch (hroch@physics.muni.cz)

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

MuniSesame::MuniSesame(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"sesame")
{
  if( cmd.Found("verbose") )
    Write("VERBOSE = T");

  if( cmd.GetParamCount() > 0 )
    Write("OBJECT = '"+cmd.GetParam(0)+"'");
  else
    wxLogFatalError("Object name not given.");
}
