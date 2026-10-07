/*

  XMunipack


  Copyright © 2009-2015, 2017-22 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _XMUNIPACK_H_
#define _XMUNIPACK_H_

#include "../config.h"
#include "mconfig.h"
#include "view.h"
#include "browser.h"
#include <wx/wx.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


class XMunipack: public wxApp
{
  bool OnInit();
  int OnExit();
  void OnEventLoopEnter(wxEventLoopBase*);

  MuniConfig *config;
  MuniView *view;
  MuniBrowser *browser;

};

#endif
