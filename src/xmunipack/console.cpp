/*

  xmunipack - console

  * beautiful replacement for wxLogWindow 


  **** CURRENTLY UNUSED ***


  Copyright © 2009-2011 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/wx.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>
#include <wx/textctrl.h>
#include <wx/srchctrl.h>
#include <wx/listctrl.h>

BEGIN_EVENT_TABLE(MuniConsole, wxFrame)
EVT_CLOSE(MuniConsole::OnClose)
EVT_SIZE(MuniConsole::OnSize)
EVT_MENU(wxID_CLOSE, MuniConsole::FileClose)
EVT_MENU(wxID_SAVE, MuniConsole::SaveOutput)
EVT_MENU(wxID_CLEAR, MuniConsole::ClearOutput)
EVT_MENU(ID_LOG_WRAP, MuniConsole::OnWrap)
EVT_MENU_RANGE(ID_LOG_ALL,ID_LOG_ERRORS, MuniConsole::SelectOutput)
EVT_TEXT_ENTER(wxID_ANY,MuniConsole::Search)
EVT_SEARCHCTRL_SEARCH_BTN(wxID_ANY,MuniConsole::Search)
EVT_SEARCHCTRL_CANCEL_BTN(wxID_ANY,MuniConsole::SearchFinish)
END_EVENT_TABLE()



MuniConsole::MuniConsole(wxWindow *w, MuniConfig *c):
wxFrame(w,wxID_ANY,"Console",wxDefaultPosition,c->console_size),config(c),
  logs(ID_LOG_ALL),wrap(c->console_wrap)
{
  SetIcon(config->munipack_icon);
  MuniArtIcons iconses(wxART_MENU,wxSize(16,16));

  // menus
  wxMenu *menuFile = new wxMenu;
  menuFile->Append(wxID_SAVE);
#ifndef __WXMAC__
  menuFile->AppendSeparator();
  menuFile->Append(wxID_CLOSE);
#endif

  menuView = new wxMenu;
  menuView->AppendRadioItem(ID_LOG_ALL,"All");
  menuView->AppendRadioItem(ID_LOG_LOGS,"Logs");
  menuView->AppendRadioItem(ID_LOG_ERRORS,"Errors");
  menuView->AppendSeparator();
  menuView->AppendCheckItem(ID_LOG_WRAP,"Wrap Lines");
  menuView->Append(wxID_CLEAR);
  menuView->Check(ID_LOG_WRAP,wrap);

  wxMenu *menuHelp = new wxMenu;
  menuHelp->Append(wxID_HELP);
  menuHelp->Append(ID_BUG,"Report a Bug ...");
  menuHelp->Append(wxID_ABOUT);

  wxMenuBar *menuBar = new wxMenuBar();
  menuBar->Append(menuFile,"&File");
  menuBar->Append(menuView,"&View");
  menuBar->Append(menuHelp,"&Help");
  SetMenuBar(menuBar);

  // toolbar
  wxToolBar *tbar = new wxToolBar(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
				  wxTB_TOP | wxTB_HORIZONTAL | wxTB_TEXT);
  tbar->SetFont(*wxSMALL_FONT);
  MuniArtIcons icons(wxART_TOOLBAR,wxSize(22,22));
  tbar->AddTool(wxID_SAVE,"Save",icons.Icon(wxART_FILE_SAVE),"Save to a file.");
  tbar->AddStretchableSpace();
  wxSize s = GetSize();
  search = new wxSearchCtrl(tbar,wxID_ANY,wxEmptyString,wxDefaultPosition,
			    wxSize((20*s.GetWidth())/100,wxDefaultCoord),
			    wxTE_PROCESS_ENTER);
  tbar->AddControl(search);
  tbar->Realize();

  tsel = new wxToolBar(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
		       wxTB_BOTTOM|wxTB_HORIZONTAL|wxTB_TEXT|wxTB_HORZ_TEXT);

  // create fake image
  wxImage img(1,1);
  img.InitAlpha();
  img.SetAlpha(0,0,0);
  wxBitmap bmp(img);
  tsel->AddRadioTool(ID_LOG_ALL,"All",bmp,bmp);
  tsel->AddRadioTool(ID_LOG_LOGS,"Logs",iconses.Icon(wxART_INFORMATION),
		     iconses.Icon(wxART_INFORMATION));
  tsel->AddRadioTool(ID_LOG_ERRORS,"Errors",iconses.Icon(wxART_ERROR),
		     iconses.Icon(wxART_ERROR));
  tsel->Realize();
  
  wxToolBar *tbot = new wxToolBar(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
				  wxTB_BOTTOM | wxTB_HORIZONTAL);
  tbot->AddTool(wxID_CLEAR,wxEmptyString,icons.Icon(wxART_CUT),"Clear");
  tbot->Realize();
  

  // font
  wxFont sf = wxSystemSettings::GetFont(wxSYS_ANSI_FIXED_FONT);
  wxFont fixed(sf.GetPointSize(),wxFONTFAMILY_MODERN,wxFONTSTYLE_NORMAL,
	       wxFONTWEIGHT_NORMAL);

  wxImageList *ilist = new wxImageList();
  ilist->Add(iconses.Icon(wxART_INFORMATION));
  ilist->Add(iconses.Icon(wxART_ERROR));

  list = new wxListCtrl(this,wxID_ANY,wxDefaultPosition, wxDefaultSize,
			wxLC_REPORT|wxLC_NO_HEADER);
  list->SetFont(fixed);
  list->AssignImageList(ilist,wxIMAGE_LIST_SMALL);
  list->InsertColumn(0,wxEmptyString);

  // sizer
  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(tbar,wxSizerFlags().Expand());
  topsizer->Add(tsel,wxSizerFlags().Expand());
  topsizer->Add(list,wxSizerFlags(1).Expand());
  topsizer->Add(tbot,wxSizerFlags().Expand());
  SetSizer(topsizer);
}

MuniConsole::~MuniConsole()
{
  config->console_size = GetSize();
  config->console_wrap = wrap;
}

void MuniConsole::AppendOutput(const wxString& a)
{
  MuniConsoleItem item(a,0);
  listitem.push_back(item);
  AppendItems(listitem.size() - 1);
}

void MuniConsole::AppendOutput(const wxArrayString& a)
{
  size_t i0 = listitem.size();

  for(size_t i = 0; i < a.GetCount(); i++) {
    MuniConsoleItem item(a.Item(i),0);
    listitem.push_back(item);
  }

  AppendItems(i0);
}

void MuniConsole::AppendError(const wxString& a)
{
  MuniConsoleItem item(a,1);
  listitem.push_back(item);
  AppendItems(listitem.size() - 1);
}

void MuniConsole::AppendError(const wxArrayString& a)
{
  size_t i0 = listitem.size();

  for(size_t i = 0; i < a.GetCount(); i++) {
    MuniConsoleItem item(a.Item(i),1);
    listitem.push_back(item);
  }

  AppendItems(i0);
}

void MuniConsole::AppendItems(size_t i0)
{
  long k = 0;
  for(size_t i = i0; i < listitem.size(); i++) {
    MuniConsoleItem item(listitem[i]);

    wxString a = item.GetText();
    wxArrayString lines;
    if( wrap ) {
      const size_t step = 80;
      lines.Add(a.Mid(0,step));
      for(size_t l = step; l < a.Len(); l+=step) {
	lines[lines.GetCount()-1] += " &";
	lines.Add("  "+a.Mid(l,step));
      }
    }
    else
      lines.Add(a);

    for(size_t l = 0; l < lines.GetCount(); l++) {

      if( logs == ID_LOG_ALL )
	list->InsertItem(k+i0,lines[l],item.GetImageIndex());
      else if( logs == ID_LOG_LOGS && item.GetImageIndex() == 0 )
	list->InsertItem(k+i0,lines[l],item.GetImageIndex());
      else if( logs == ID_LOG_ERRORS && item.GetImageIndex() == 1 )
	list->InsertItem(k+i0,lines[l],item.GetImageIndex());

      k++;
    }
  }
}

void MuniConsole::SetLogs(int logs) 
{
  list->ClearAll();
  list->SetFont(fixed);
  list->InsertColumn(0,wxEmptyString);
  SetColumnWidth();

  AppendItems(0);
}


void MuniConsole::SetColumnWidth()
{
  wxSize s = GetClientSize();
  list->SetColumnWidth(0,s.GetWidth());
}

void MuniConsole::OnSize(wxSizeEvent& event)
{
  SetColumnWidth();
  event.Skip();
}


void MuniConsole::OnClose(wxCloseEvent& event)
{
  //  wxPostEvent(GetParent(),event);
  wxQueueEvent(GetParent(),event.Clone());
}

void MuniConsole::FileClose(wxCommandEvent& event)
{
  //  wxCloseEvent e(wxEVT_CLOSE_WINDOW,GetId());
  //  wxPostEvent(GetParent(),e);
  wxQueueEvent(GetParent(),new wxCloseEvent(wxEVT_CLOSE_WINDOW,GetId()));
}


void MuniConsole::SaveOutput(wxCommandEvent& event)
{
  wxFileDialog select(this,"Choose a text file",wxEmptyString,wxEmptyString,
		      "Text files (*.txt)|*.txt| All files (*)|*",
		      wxFD_SAVE|wxFD_OVERWRITE_PROMPT|wxFD_CHANGE_DIR);
  if (select.ShowModal() == wxID_OK ) {

    wxFileOutputStream fout(select.GetPath());
    wxTextOutputStream tout(fout);

    for(size_t i = 0; i < listitem.size(); i++) {
	MuniConsoleItem item(listitem[i]);
	wxString label = 
	  item.GetImageIndex() == 0 ? "Output:" : "Error:";
	tout << label+item.GetText() + "\n";
      }
  }
}

void MuniConsole::ClearOutput(wxCommandEvent& event)
{
  list->ClearAll();
  list->InsertColumn(0,wxString('*',88));
  list->SetColumnWidth(0,wxLIST_AUTOSIZE_USEHEADER);
  listitem.clear();

  logs = ID_LOG_ALL;
  menuView->Check(logs,true);
  tsel->ToggleTool(logs,true);  
}

void MuniConsole::SelectOutput(wxCommandEvent& event)
{
  logs = event.GetId();
  menuView->Check(logs,true);
  tsel->ToggleTool(logs,true);
  SetLogs(logs);
}

void MuniConsole::Search(wxCommandEvent& event)
{
  wxLogDebug(event.GetString());
  wxLogDebug(search->GetValue());

  wxString x = search->GetValue();

  for(size_t i = 0; i < (size_t) list->GetItemCount(); i++){
    wxListItem item;
    item.SetId(i);
    list->GetItem(item);
    wxString line = item.GetText();
    if( line.Find(x) != wxNOT_FOUND ) {
      long s = item.GetState();
      item.SetState(s|wxLIST_STATE_SELECTED);
    }
    list->SetItem(item);

  }

}

void MuniConsole::SearchFinish(wxCommandEvent& event)
{
  search->ShowCancelButton(false);

  for(size_t i = 0; i < (size_t) list->GetItemCount(); i++){
    wxListItem item;
    item.SetId(i);
    list->GetItem(item);
    long s = item.GetState();
    item.SetState(s & ~wxLIST_STATE_SELECTED);
    list->SetItem(item);
  }

}

void MuniConsole::OnWrap(wxCommandEvent& event)
{
  wrap = event.IsChecked();
  SetLogs(logs);
}

