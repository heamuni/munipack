/*

  Implements FITS file name utilities (taylored by wxFileName):

  * GNU backup-file conventions

  * file-existing


  Copyright © 2013, 2018 F.Hroch (hroch@physics.muni.cz)

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

#include <wx/wx.h>

class FITSFileName
{
public:
  FITSFileName(const wxString&);
  wxString GetFullName() const;
  wxString GetType() const { return filetype; }
  wxString GetFullPath() const { return fullpath; }
  bool IsOk() const { return status == 0; }

private:

  int status;
  wxString filename, filetype, fullpath, tmpname, suffix;

};
