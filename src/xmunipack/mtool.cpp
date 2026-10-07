/*

  Base support for tools

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

#include "xmunipack.h"
#include "mtool.h"
#include "mconfig.h"
#include "mprocess.h"
#include "event.h"
#include "help.h"
#include <wx/wx.h>
#include <wx/regex.h>
#include <wx/toolbar.h>

MuniTool::MuniTool(wxWindow *w, MuniConfig *c, const wxString& a,
		       const wxString& t):
  wxFrame(w,wxID_ANY,t,wxDefaultPosition,wxDefaultSize, wxDEFAULT_FRAME_STYLE|
	  wxFRAME_TOOL_WINDOW|wxFRAME_FLOAT_ON_PARENT),
  action(a), mproc(0), timer(this), anim(0), panel(0), index(-1),
  exitcode(-1), interrupted(false), config(c)
{
  SetIcon(config->munipack_icon);
  manpage = "man_" + action + ".html";

  SetMenuBar(CreateMenuBar());

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(CreateTools(),wxSizerFlags().Expand());
  SetSizer(topsizer);
  CentreOnParent();
};

void MuniTool::OnStop(wxCommandEvent& event)
{
  wxASSERT(mproc);
  mproc->Kill();
  interrupted = true;
}


void MuniTool::OnExec(wxCommandEvent& event)
{
  wxLogDebug("MuniTool::OnExec of "+action);
  wxASSERT(mproc == 0);
  wxASSERT(action != "");
  wxASSERT(tbar);

  // setup GUI elements
  EnableCloseButton(false);

  tool_exec = tbar->RemoveTool(wxID_EXECUTE);
  tbar->InsertTool(1,tool_stop);
  if( tbar->GetToolPos(ID_WARNING) != wxNOT_FOUND ) {
    tool_warn = tbar->RemoveTool(ID_WARNING);
    tbar->InsertTool(2,tool_info);
  }
  tbar->EnableTool(wxID_SAVE,false);
  tbar->EnableTool(wxID_INFO,false);
  tbar->Realize();

  menu_file->Enable(wxID_SAVE,false);
  menu_file->Enable(wxID_CLOSE,false);
  menu_action->Enable(wxID_EXECUTE,false);
  menu_action->Enable(wxID_STOP,true);
  menu_action->Enable(wxID_INFO,false);

  anim->Show(true);
  anim->Play();
  CleanDraw();

  SetStatus("Launching "+action+" ...");

  // start a new process
  index = 0;
  interrupted = false;
  exitcode = -1;
  out.Empty();
  timer.Start(250);
  mproc = new MuniProcess(this,action);
  mproc->SetEcho(false);
  mproc->Write("PIPELOG = T");
  OnInput(mproc);
  mproc->OnStart();

  Bind(wxEVT_END_PROCESS,&MuniTool::OnFinish,this);
  Bind(wxEVT_TIMER,&MuniTool::OnTimer,this);

  Layout();
  event.Skip();
}

void MuniTool::OnFinish(wxProcessEvent& event)
{
  wxASSERT(mproc && tbar);
  wxLogDebug("MuniTool::OnFinish of "+action+": %d",mproc->GetExitCode());

  exitcode = mproc->GetExitCode();
  bool warn = false;

  if( exitcode == 0 ) {
    // success
    out = mproc->GetOutput();
    OnOutput(out);

  }
  else {
    // fail
    if( interrupted )
      SetStatus("Interrupted.");
    else {
      warn = true;
      SetStatusError("Command "+action+" failed.");
      err = mproc->GetErrors();
      out = mproc->GetOutput();
      in = mproc->GetInput();

    }
  }

  delete mproc;
  mproc = 0;
  timer.Stop();

  // update GUI
  anim->Stop();
  anim->Show(false);

  EnableCloseButton(true);
  tool_stop = tbar->RemoveTool(wxID_STOP);
  tbar->EnableTool(wxID_SAVE,exitcode == 0);
  tbar->InsertTool(1,tool_exec);
  tbar->EnableTool(wxID_INFO,true);
  if( warn ) {
    tool_info = tbar->RemoveTool(wxID_INFO);
    tbar->InsertTool(2,tool_warn);
  }
  tbar->Realize();

  menu_file->Enable(wxID_SAVE,exitcode == 0);
  menu_file->Enable(wxID_CLOSE,true);
  menu_action->Enable(wxID_EXECUTE,true);
  menu_action->Enable(wxID_STOP,false);
  menu_action->Enable(wxID_INFO,true);

  Layout();

  Unbind(wxEVT_END_PROCESS,&MuniTool::OnFinish,this);
  Unbind(wxEVT_TIMER,&MuniTool::OnTimer,this);

  event.Skip();
}

void MuniTool::OnTimer(wxTimerEvent& event)
{
  out = mproc->GetOutput();
  event.Skip();
}


void MuniTool::OnHelp(wxCommandEvent& event)
{
  MuniHelp(config->munipack_html_dir,manpage);
}

void MuniTool::OnInfo(wxCommandEvent& event)
{
  wxLogMessage("Command: "+action);
  for(size_t i = 0; i < in.GetCount(); i++)
    wxLogMessage("Input:"+in[i]);
  for(size_t i = 0; i < out.GetCount(); i++)
    wxLogMessage("Output:"+out[i]);
  for(size_t i = 0; i < err.GetCount(); i++)
    wxLogError("Error:"+err[i]);

  if( exitcode == 0 )
    wxLogMessage("Command `"+action+"' successful.");
  else
    wxLogError("Command `"+action+"' failed.");
}

void MuniTool::SetManPage(const wxString& h) { manpage = h; }


wxSizer *MuniTool::CreateTools()
{
  MuniArtIcons ico(wxART_TOOLBAR,wxSize(22,22));

  tbar = new wxToolBar(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
		       wxTB_HORIZONTAL|wxTB_TEXT);
  tbar->AddTool(wxID_SAVE,"Save",ico.Icon("document-save"));
  tbar->AddTool(wxID_EXECUTE,"Exec",ico.Icon("system-run"));
  tbar->AddTool(wxID_STOP,"Stop",ico.Icon("process-stop"));
  tbar->AddTool(wxID_INFO,"Info",ico.Icon(wxART_INFORMATION));
  tbar->AddTool(ID_WARNING,"Info",ico.Icon(wxART_WARNING));
  tbar->AddSeparator();

  tool_stop = tbar->RemoveTool(wxID_STOP);
  tool_warn = tbar->RemoveTool(ID_WARNING);
  tbar->EnableTool(wxID_SAVE,false);
  tbar->EnableTool(wxID_INFO,false);
  tbar->Realize();

  status = new wxStaticText(this,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
			    wxALIGN_LEFT | wxST_ELLIPSIZE_END);

  anim = new wxAnimationCtrl(this,wxID_ANY,config->throbber,
    			     wxDefaultPosition,wxSize(32,32));
  anim->Show(false);

  wxBoxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);
  sizer->Add(tbar,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL));
  sizer->Add(status,wxSizerFlags(1).Align(wxALIGN_CENTER_VERTICAL).Border());
  sizer->Add(anim,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).Border());

  return sizer;
}


wxMenuBar *MuniTool::CreateMenuBar()
{
  wxMenuBar *bar = new wxMenuBar();

  menu_file = new wxMenu();
  menu_file->Append(wxID_SAVE);
  menu_file->Append(wxID_CLOSE);

  menu_action = new wxMenu();
  menu_action->Append(wxID_EXECUTE);
  menu_action->Append(wxID_STOP);
  menu_action->Append(wxID_INFO);

  wxMenu *menu_help = new wxMenu();
  menu_help->Append(wxID_HELP);
  menu_help->Append(wxID_ABOUT);

  bar->Append(menu_file,"&File");
  bar->Append(menu_action,"&Action");
  bar->Append(menu_help,"&Help");

  Bind(wxEVT_TOOL,&MuniTool::OnExec,this,wxID_EXECUTE);
  Bind(wxEVT_TOOL,&MuniTool::OnStop,this,wxID_STOP);
  Bind(wxEVT_TOOL,&MuniTool::OnSave,this,wxID_SAVE);
  Bind(wxEVT_TOOL,&MuniTool::OnInfo,this,wxID_INFO);
  Bind(wxEVT_TOOL,&MuniTool::OnInfo,this,ID_WARNING);
  Bind(wxEVT_MENU,&MuniTool::OnHelp,this,wxID_HELP);
  Bind(wxEVT_MENU,&MuniTool::OnAbout,this,wxID_ABOUT);
  Bind(wxEVT_MENU,&MuniTool::OnClose,this,wxID_CLOSE);
  Bind(wxEVT_CLOSE_WINDOW,&MuniTool::OnCloseWin,this);

  menu_file->Enable(wxID_SAVE,false);
  menu_action->Enable(wxID_STOP,false);
  menu_action->Enable(wxID_INFO,false);

  return bar;
}

void MuniTool::SetPanel(wxPanel *p, const wxSizerFlags& flags)
{
  panel = p;
  wxSizer *topsizer = GetSizer();
  topsizer->Add(panel,flags);
  Layout();
}


wxString MuniTool::Parser(const wxString& out, const wxString& keyword) const
{
  wxRegEx re("^[ ]?=(.*)> (.+)");
  wxASSERT(re.IsValid());

  if( re.Matches(out) ) {

    wxString key(re.GetMatch(out,1));
    wxString value(re.GetMatch(out,2));

    if( key == keyword )
      return value;
  }

  return "";
}

wxArrayString MuniTool::GetOutput() const
{
  return out;
}


wxArrayString MuniTool::GetLastOutput()
{
  wxArrayString last;
  for(size_t i = index; i < out.GetCount(); i++)
    last.Add(out[i]);

  index = out.GetCount();

  return last;
}

void MuniTool::SetStatus(const wxString& fmt, ...)
{
  wxASSERT(status);

  wxString label;

  va_list par;
  va_start(par, fmt);
  label.PrintfV(fmt,par);
  va_end(par);

  status->SetLabelMarkup("<small>"+label+"</small>");
  Layout();
}

void MuniTool::SetStatusDisplay(const wxString& fmt, ...)
{
  wxASSERT(status);

  wxString label;

  va_list par;
  va_start(par, fmt);
  label.PrintfV(fmt,par);
  va_end(par);

  status->SetLabelMarkup(label);
  Layout();
}

void MuniTool::SetStatusWarning(const wxString& fmt, ...)
{
  wxASSERT(status);

  wxString label;

  va_list par;
  va_start(par, fmt);
  label.PrintfV(fmt,par);
  va_end(par);

  status->SetLabelMarkup("<b>"+label+"</b>");
  Layout();
}

void MuniTool::SetStatusError(const wxString& fmt, ...)
{
  wxASSERT(status);

  wxString label;

  va_list par;
  va_start(par, fmt);
  label.PrintfV(fmt,par);
  va_end(par);

  status->SetLabelMarkup("<span fgcolor=\"#ff0000\"><b>"+label+"</b></span>");
  Layout();
}

void MuniTool::OnSave(wxCommandEvent& event)
{
  wxASSERT(tbar);

  wxLogDebug("MuniTool::OnSave");
  menu_file->Enable(wxID_SAVE,false);

  tbar->EnableTool(wxID_SAVE,false);
  tbar->Realize();
}

void MuniTool::OnClose(wxCommandEvent& WXUNUSED(event))
{
  wxLogDebug("MuniTool::OnClose");
  Close();
}

void MuniTool::OnCloseWin(wxCloseEvent& event)
{
  wxLogDebug("MuniTool::OnCloseWin");

  /*
  if( event.CanVeto() ) {
    wxMessageDialog dialog(this,"Unsaved work: Do you want exit anyway? ",
			   GetTitle(),
			   wxICON_EXCLAMATION|wxYES_NO|wxCANCEL|wxNO_DEFAULT);
    dialog.SetExtendedMessage("Your changes will be lost if you don't save them.");
    if( dialog.ShowModal() == wxID_NO )
      return;
  }
  */

  /*
  if( mproc && event.CanVeto() ) {
    wxMessageDialog dialog(this,"Work in progress ...", GetTitle(),
			   wxICON_EXCLAMATION|wxYES_NO|wxCANCEL|wxNO_DEFAULT);
    dialog.SetExtendedMessage("Do you want interrput it anyway?");
    if( dialog.ShowModal() == wxID_NO )
      return;
  }
  */

  // clear detected sources
  CleanDraw();

  // notify parent
  wxCommandEvent e(EVT_TOOL_FINISH,this->GetId());
  e.SetInt(exitcode);
  wxQueueEvent(GetParent(),e.Clone());
}

void MuniTool::OnAbout(wxCommandEvent& WXUNUSED(event))
{
  MuniAbout(config->munipack_icon);
}
