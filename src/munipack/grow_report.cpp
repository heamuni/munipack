/*

  Growth-Curve Photometry Report

  Copyright © 2016, 2018 F.Hroch (hroch@physics.muni.cz)

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


bool Munipack::grow_report(MuniProcess *action, MuniCmdLineParser& cmd)
{
  for(size_t i = 0; i < cmd.GetFilesCount(); i++)
    action->Write("FILE = '"+cmd.GetFile(i)+"'");

  return true;
}
