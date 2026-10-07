/*

  xmunipack - config

  Copyright © 2012, 2018-22 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _XMUNIPACK_ICON_H_
#define _XMUNIPACK_ICON_H_

#include "fits.h"
#include "mconfig.h"


class MuniIcon
{
  const MuniConfig *config;
  const wxColour colour;

  wxImage Padding(const wxImage&,int) const;

public:
  MuniIcon(const MuniConfig*);

  void GetIcons(const FitsFile&, wxImage&, std::vector<wxImage>&) const;
  wxImage GetCover(const wxImage&, const wxString& ="", int =-1) const;
  std::vector<wxImage> GetList(const FitsMeta&) const;
  std::vector<wxImage> GetList(const FitsFile&) const;
  std::vector<wxImage> MergeList(const FitsFile&, const std::vector<wxImage>&) const;
  wxImage ImageIcon(const wxImage&) const;
  wxBitmap Utf8Icon(double, const wxString&,const wxString& ="") const;

  // both obsolete
  static wxImage ListIcon(const wxImage&,int,
			  const wxColour& =wxColour(255,255,255));
  static wxImage BulletIcon(const wxSize&, const wxColour&);

};

#endif
