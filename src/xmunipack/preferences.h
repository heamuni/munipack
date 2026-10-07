/*

  xmunipack preferences

  Copyright © 2021-22 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _PREFERENCES_H_
#define _PREFERENCES_H_

#include "../config.h"
#include "config.h"
#include "mconfig.h"
#include <wx/wx.h>
#include <wx/propdlg.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


class MuniPreferences: public wxPropertySheetDialog
{
  MuniConfig *config;

  void FeedView(wxPanel *);
  void FeedKeywords(wxPanel *);

  void OnClose(wxCloseEvent&);
  void OnIdle(wxIdleEvent&);
  void OnCooType(wxCommandEvent&);
  void OnQPhType(wxCommandEvent&);
  void OnOk(wxCommandEvent&);
  void OnCancel(wxCommandEvent&);

public:
  MuniPreferences(wxWindow *, MuniConfig *);
  void SelectPage(int);

};

#endif
