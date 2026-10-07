/*

  XMunipack -- export of FITS HDU

  Copyright © 2022 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _XMUNIPACK_EXPORT_H
#define _XMUNIPACK_EXPORT_H

#include "fits.h"
#include "view.h"
#include <wx/wx.h>
#include <wx/thread.h>


class MuniExportDialog: public wxFileDialog
{
public:
  MuniExportDialog(wxWindow *, const wxString&,const wxString&,const wxString&,
		   const wxString&, long);
  bool Save(const FitsArray&, const MuniDisplay*);

};

#endif
