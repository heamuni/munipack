/*

  VO Table manipulations

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
#include <wx/filename.h>

using namespace std;


MuniVotable::MuniVotable(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"votable")
{
  CommonOutputSingle(cmd);

  wxString a, type;
  double x;
  long l;

  if( cmd.Found("o",&a) || cmd.Found("output",&a) ) {
    wxFileName fn(a);
    if( fn.IsOk() )
      type = fn.GetExt();
  }

  if( cmd.Found("type",&a) )
    type = a;

  if( ! type.IsEmpty() )
    Write("TYPE = '"+type.Upper()+"'");

  if( cmd.Found("pt",&a) )
    Write("PROJ TYPE = '" + a.MakeUpper() + "'");

  if( cmd.Found("pa",&x) )
    Write("PROJ ALPHA = %lf",x);

  if( cmd.Found("pd",&x) )
    Write("PROJ DELTA = %lf",x);

  if( cmd.Found("pw",&l) )
    Write("PROJ WIDTH = %ld",l);

  if( cmd.Found("ph",&l) )
    Write("PROJ HEIGHT = %ld",l);

  if( cmd.Found("ps",&x) )
    Write("PROJ SCALE = %lf",x);

  if( cmd.Found("ml",&x) )
    Write("MAG LIMIT = %lf",x);

  if( cmd.Found("col-mag",&a) )
    Write("COL_MAG = '" + a + "'");

  if( cmd.Found("col-ra",&a) )
    Write("COL_RA = '" + a + "'");

  if( cmd.Found("col-dec",&a) )
    Write("COL_DEC = '" + a + "'");

  WriteFilesSingle(cmd.GetFiles());
}
