/*

  Art Icons

  Copyright © 2019-22 F.Hroch (hroch@physics.muni.cz)

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

#include "mconfig.h"
#include <wx/wx.h>
#include <wx/artprov.h>


MuniArtIcons::MuniArtIcons(const wxArtClient& c, const wxSize& s):
  client(c), size(s) {}

wxBitmap MuniArtIcons::Icon(const wxArtID& id) const
{
  wxBitmap b = wxArtProvider::GetBitmap(id,client,size);
  if( b.IsOk() )
    return b;
  else
    return wxArtProvider::GetBitmap(wxART_MISSING_IMAGE,client,size);
}
