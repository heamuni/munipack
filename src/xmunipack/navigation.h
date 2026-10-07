/*

  XMunipack -- Navigation bar for Browser

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

#include "../config.h"
#include "mconfig.h"
#include "event.h"
#include <wx/wx.h>
#include <vector>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


/*
class MuniTextCtrl: public wxTextCtrl
{

  void OnChar(wxKeyEvent&);
  void OnText(wxCommandEvent&);
  wxArrayString GetSubdirs(const wxString&) const;

public:
  MuniTextCtrl(wxWindow *,const wxString&,long);
};
*/

class MuniNavigation: public wxPanel
{
  MuniConfig *config;
  wxTextCtrl *cmd;
  wxStaticBitmap *moon;
  wxBitmapButton *bopen;
  double phase;
  wxString pwd,mask;

  void OnUpdateMoon(wxUpdateUIEvent& event);
  void OnEnter(wxCommandEvent&);
  void GetPath(const wxString&, wxString&, wxString&) const;

public:
  MuniNavigation(wxWindow *, MuniConfig *);
  void NewMoon();
  void HideMoon();
  void PhaseMoon(double);
  void ChangeDir(const wxString&, const wxString&);

};
