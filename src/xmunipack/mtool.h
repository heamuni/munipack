/*

  xmunipack - base for tools

  Copyright © 2019 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _XMUNIPACK_MTOOL_H_
#define _XMUNIPACK_MTOOL_H_

#include "mconfig.h"
#include "mprocess.h"
#include <wx/wx.h>

class MuniTool: public wxFrame
{
  const wxString action;
  MuniProcess *mproc;
  wxTimer timer;
  wxAnimationCtrl *anim;
  wxPanel *panel;
  wxMenu *menu_file, *menu_action;
  wxToolBar *tbar;
  wxToolBarToolBase *tool_exec, *tool_stop, *tool_warn, *tool_info;
  wxStaticText *status;
  wxString manpage;
  long index;
  int exitcode;
  bool interrupted;
  wxArrayString in, out, err;

  void OnHelp(wxCommandEvent&);
  void OnAbout(wxCommandEvent&);
  void OnClose(wxCommandEvent&);
  void OnCloseWin(wxCloseEvent&);
  void OnSave(wxCommandEvent&);
  void OnExec(wxCommandEvent&);
  void OnStop(wxCommandEvent&);
  void OnInfo(wxCommandEvent&);
  void OnTimer(wxTimerEvent&);
  void OnFinish(wxProcessEvent&);
  wxMenuBar *CreateMenuBar();
  wxSizer *CreateTools();

 protected:

  MuniConfig *config;
  void SetManPage(const wxString&);
  void SetStatus(const wxString&, ...);
  void SetStatusDisplay(const wxString&, ...);
  void SetStatusError(const wxString&, ...);
  void SetStatusWarning(const wxString&, ...);
  wxString Parser(const wxString&, const wxString&) const;
  wxArrayString GetOutput() const;
  wxArrayString GetLastOutput();
  virtual void OnInput(MuniProcess *) = 0;
  virtual void OnOutput(const wxArrayString&) {}
  virtual void CleanDraw() const {}
  void SetPanel(wxPanel *, const wxSizerFlags&);

  MuniTool(wxWindow *, MuniConfig *, const wxString&, const wxString&);

};

#endif
