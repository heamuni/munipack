/*

  XMunipack


  Copyright © 2021 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _CONSOLE_H_
#define _CONSOLE_H_

#include "../config.h"
#include "config.h"
#include "mconfig.h"
#include <wx/wx.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


class MuniConsoleItem
{
public:
  MuniConsoleItem(const wxString& t, int i): text(t),iindex(i) {}
  wxString GetText() const { return text; }
  int GetImageIndex() const { return iindex; }

private:

  wxString text;
  int iindex;

};


class MuniConsole: public wxFrame
{
public:
  MuniConsole(wxWindow *, MuniConfig *);
  virtual ~MuniConsole();
  void AppendOutput(const wxArrayString&);
  void AppendError(const wxArrayString&);
  void AppendOutput(const wxString&);
  void AppendError(const wxString&);

private:

  MuniConfig *config;
  //  wxTextCtrl *output,*error;
  wxListCtrl *list;
  wxMenu *menuView;
  wxToolBar *tsel;
  wxSearchCtrl *search;
  wxFont fixed;
  std::vector<MuniConsoleItem> listitem;

  int logs;
  bool wrap;

  void OnClose(wxCloseEvent&);
  void OnSize(wxSizeEvent&);
  void FileClose(wxCommandEvent&);
  void SaveOutput(wxCommandEvent&);
  void ClearOutput(wxCommandEvent&);
  void SelectOutput(wxCommandEvent&);
  void AppendItems(size_t);
  void SetLogs(int);
  void SetColumnWidth();
  void Search(wxCommandEvent&);
  void SearchFinish(wxCommandEvent&);
  void OnWrap(wxCommandEvent&);

  DECLARE_EVENT_TABLE()
};

#endif
