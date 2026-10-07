/*

  xmunipack - file browser window


  Copyright © 2009-2012, 2019-22 F.Hroch (hroch@physics.muni.cz)

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


#include "browser.h"
#include "view.h"
#include "help.h"
#include <wx/wx.h>
#include <wx/dirdlg.h>
#include <wx/listctrl.h>
#include <wx/filename.h>
#include <wx/aboutdlg.h>
#include <wx/dirctrl.h>
#include <wx/tooltip.h>
#include <wx/srchctrl.h>
#include <wx/wupdlock.h>
#include <wx/utils.h>
#include <wx/busyinfo.h>
#include <vector>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

#define MENU_LABEL_ARRANGE "Arrange by"
#define MENU_LABEL_LABELS "Labels as"

using namespace std;


MuniBrowser::MuniBrowser(wxWindow *w, MuniConfig *c):
  wxFrame(w, wxID_ANY,"Munipack: Browser",wxDefaultPosition,
	  c->browser_size, wxDEFAULT_FRAME_STYLE),
  config(c),
  menuView(new wxMenu), //menuAct(new wxMenu),
  menuArrange(new wxMenu), menuLabels(new wxMenu), list(0),preferences(NULL),
  tbar(0), topsizer(0),
  pwd(wxFileName::GetCwd()), mask(config->browser_fitsmask),
  metaload(false), shutdown(false)
{
  SetIcon(config->munipack_icon);

#ifdef __WXMAC__
  // all main windows?
  SetExitOnFrameDelete(true);
  wxMenuBar *menubar = new wxMenuBar;
  wxMenuBar::MacSetCommonMenuBar(menubar);
#endif

  // toolbars
  tbar = new wxToolBar(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
		       wxTB_TOP | wxTB_HORIZONTAL | wxTB_TEXT);
  tbar->SetFont(*wxSMALL_FONT);
  MuniArtIcons icons(wxART_TOOLBAR,wxSize(22,22));
  tbar->AddTool(wxID_OPEN,"Open",icons.Icon(wxART_FILE_OPEN),"Change directory");
  SetToolBar(tbar);
  tbar->Show(config->browser_toolbar);

  // menus
  /*
  wxMenu *menuWizards = new wxMenu;
  menuWizards->Append(ID_LC,"Light curve");
  menuWizards->Append(ID_CMD,"Color diagram");
  */

  menuFile = new wxMenu;
  menuFile->Append(wxID_NEW);
  menuFile->Append(wxID_OPEN);
  menuFile->Append(ID_NEW_VIEW,"View...");

  //  menuFile->AppendSubMenu(menuExport,wxT("Export"));
  //  menuFile->Append(ID_EXPORT,"Export As...");
  //  menuFile->Append(wxID_SAVEAS,wxT("Export Files"));
  //  menuFile->AppendSeparator();
  //  menuFile->Append(wxID_ADD,wxT("New Archive..."));
  //  menuFile->Append(wxID_REMOVE,wxT("Exclude Archive"));
#ifdef __WXMAC__
  menuFile->Append(wxID_CLOSE);
#endif
  menuFile->AppendSeparator();
  menuFile->Append(wxID_PROPERTIES);
#ifndef __WXMAC__
  menuFile->AppendSeparator();
  menuFile->Append(wxID_CLOSE,wxEmptyString,
		 "Close this window (or finish xmunipack).");
#endif
#ifdef __WXMAC__
  menuFile->AppendSeparator();
  //  menuFile->Append(wxID_EXIT);
#endif

  wxMenu *menuEdit = new wxMenu;
  /*
  menuEdit->Append(wxID_CUT);
  menuEdit->Append(wxID_COPY);
  menuEdit->Append(wxID_PASTE);
  menuEdit->AppendSeparator();
  menuEdit->Append(wxID_SELECTALL);
  menuEdit->AppendSeparator();
  */
  menuEdit->Append(wxID_PREFERENCES);
  menuEdit->AppendSeparator();
  menuEdit->Append(wxID_REFRESH);
  menuEdit->Append(wxID_STOP);
  // edit header?

  menuArrange->AppendRadioItem(ID_SORT_UNSORT,"Unsorted");
  menuArrange->AppendRadioItem(ID_SORT_FILENAME,"Filename");
  menuArrange->AppendRadioItem(ID_SORT_OBJECT,"Object");
  // menuArrange->AppendRadioItem(wxID_ANY,wxT("By Type"));
  //  menuArrange->AppendRadioItem(wxID_ANY,wxT("By Dimensions"));
  menuArrange->AppendRadioItem(ID_SORT_DATEOBS,"Time");
  menuArrange->AppendRadioItem(ID_SORT_FILTER,"Filter");
  menuArrange->AppendRadioItem(ID_SORT_EXPOSURE,"Exposure Time");
  menuArrange->AppendRadioItem(ID_SORT_SIZE,"Filesize");
  //  menuArrange->AppendRadioItem(ID_SORT_KEY,"Predefined Keyword");
  menuArrange->AppendSeparator();
  menuArrange->AppendCheckItem(ID_SORT_REVERSE,"Reverse");

  menuLabels->AppendRadioItem(ID_LABEL_NOLABEL,"No Labels");
  menuLabels->AppendRadioItem(ID_LABEL_FILENAME,"Filename");
  menuLabels->AppendRadioItem(ID_LABEL_OBJECT,"Object");
  menuLabels->AppendRadioItem(ID_LABEL_DATE,"Date");
  menuLabels->AppendRadioItem(ID_LABEL_TIME,"Time");
  menuLabels->AppendRadioItem(ID_LABEL_FILTER,"Filter");
  menuLabels->AppendRadioItem(ID_LABEL_EXPOSURE,"Exposure Time");
  //  menuLabels->AppendRadioItem(ID_LABEL_KEY,"Predefined Keyword");

  menuView->AppendRadioItem(wxID_VIEW_LARGEICONS,"View as Icons");
  menuView->AppendRadioItem(wxID_VIEW_LIST,"View as List");
  menuView->AppendSubMenu(menuArrange,MENU_LABEL_ARRANGE);
  menuView->AppendSubMenu(menuLabels,MENU_LABEL_LABELS);
  /*
  menuView->AppendCheckItem(ID_VIEW_MARK,"Mark by Type",
		   "Mark scientific, flats, dark, .. exposures.");
  */
  menuView->AppendSeparator();
  menuView->AppendCheckItem(ID_BROWSER_TOOLBAR,"Show Toolbar",
    "Change visibility of toolbar (shorthand buttons with icons on top)");
  menuView->Check(ID_BROWSER_TOOLBAR,config->browser_toolbar);

  //  menuAct->Append(ID_COLORING,wxT("Coloring..."));
  //  menuAct->Append(ID_FIND,wxT("Find Stars..."),
  //		  wxT("Detection of stars on images."));
  //  menuAct->Append(ID_PHOT,"Photometry...","Aperture photometry.");
  //  menuAct->Append(ID_PROFILE,wxT("Profile Photometry..."),
  //		  wxT("Profile photometry."));
  //  menuAct->Append(ID_MATCH,wxT("Matching of Images..."),
  //		  wxT("By genetic algorithm or Fourier transformation."));
  /*
  menuAct->Append(ID_ASTROMET,"Astrometry...","Astrometry calibration");
  menuAct->Append(ID_KOMBINE,"Compose Images...","Composition of images.");
  menuAct->Append(ID_DEKON,"Deconvolution...","Deconvolution of images.");
  menuAct->Append(ID_DARKBAT,"Correction");
  menuAct->Append(ID_AVERAGE,"Average");
  menuAct->Append(ID_ARITH,"Arithmetic...","Basic arithmetical operations with images.");
  */
  //  menuAct->AppendSeparator();
  //  menuAct->AppendCheckItem(ID_LOG, wxT("&Log..."));

  menuHelp = new wxMenu;
  //  menuHelp->AppendCheckItem(ID_LOG,"&Log...");
  //  menuHelp->Append(ID_LOG,"Log");
  menuHelp->Append(wxID_ABOUT);

  wxMenuBar *menuBar = new wxMenuBar;
  menuBar->Append(menuFile,"&File");
  menuBar->Append(menuEdit,"&Edit");
  menuBar->Append(menuView,"&View");
  // operations?
  /*
  menuBar->Append(menuAct,"&Tools");
  menuBar->Append(menuWizards,"&Assistants");
  */
  menuBar->Append(menuHelp,"&Help");

  // gray temporary unaccesible tools
  //  menuFile->Enable(wxID_SAVEAS,false);
  //  menuFile->Enable(ID_EXPORT,false);

  /*
  menuArrange->Enable(ID_SORT_KEY,false);
  menuLabels->Enable(ID_LABEL_KEY,false);
  */
  //  menuView->Enable(ID_VIEW_MARK,false);
  menuEdit->Enable(wxID_STOP,false);

  //  menuAct->Enable(ID_PHOT,false);
  //  menuAct->Enable(ID_PROFILE,false);
  //  menuAct->Enable(ID_MATCH,false);
  //  menuAct->Enable(ID_ASTROMET,false);
  /*
  menuAct->Enable(ID_KOMBINE,false);
  menuAct->Enable(ID_DEKON,false);
  menuAct->Enable(ID_ARITH,false);
  menuWizards->Enable(ID_LC,false);
  menuWizards->Enable(ID_CMD,false);
  */
  //  menuBar->Enable(menuBar->FindMenu(wxT("Assistants")),false);

  SetMenuBar(menuBar);

  // navigation bar
  navbar = new MuniNavigation(this,config);

  // icon space
  if( config->browser_iconlist == wxLC_ICON )
    list = new MuniListIcon(this,ID_LIST,config);
  else if( config->browser_iconlist == wxLC_REPORT )
    list = new MuniListList(this,ID_LIST,config);
  else
    wxFAIL_MSG("Error: Unsupported list type.");

  // sizers
  topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(navbar,wxSizerFlags().Expand());
  topsizer->Add(list,wxSizerFlags(1).Expand().Border());
  SetSizer(topsizer);
  SetSizeHints(wxSize(320,240));

  menuLabels->Check(config->browser_labeltype,true);
  menuArrange->Check(config->browser_sorttype,true);
  menuArrange->Check(ID_SORT_REVERSE,config->browser_reverse);
  switch (config->browser_iconlist) {
  case wxLC_ICON:   menuView->Check(wxID_VIEW_LARGEICONS,true); break;
  case wxLC_REPORT: menuView->Check(wxID_VIEW_LIST,true); break;
  }
  menuView->Enable(menuView->FindItem(MENU_LABEL_LABELS),
		   config->browser_iconlist==wxLC_ICON);

  Bind(wxEVT_CLOSE_WINDOW,&MuniBrowser::OnClose,this);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnNewBrowser,this,wxID_NEW);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnNewView,this,ID_NEW_VIEW);
  //  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::FileOpen,this,wxID_OPEN);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnDirOpen,this,wxID_OPEN);
  Bind(wxEVT_BUTTON,&MuniBrowser::OnDirOpen,this,wxID_OPEN);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnFileClose,this,wxID_CLOSE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnPreferences,this,wxID_PREFERENCES);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnProperties,this,wxID_PROPERTIES);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnRefresh,this,wxID_REFRESH);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnStop,this,wxID_STOP);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnSelectAll,this,wxID_SELECTALL);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnCut,this,wxID_CUT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnCopy,this,wxID_COPY);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnPaste,this,wxID_PASTE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::SelectItem,this,wxID_BACKWARD);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::SelectItem,this,wxID_FORWARD);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnLabel,this,ID_LABEL_FILENAME,ID_LABEL_NOLABEL);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnSort,this,ID_SORT_UNSORT,ID_SORT_KEY);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnReverse,this,ID_SORT_REVERSE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnIconList,this,wxID_VIEW_LARGEICONS);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnIconList,this,wxID_VIEW_LIST);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnShowToolbar,this,ID_BROWSER_TOOLBAR);
  //  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::ViewLog,this,ID_LOG);
  //  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::ShowLog,this,ID_WARNING);
  //  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::FindStars,this,ID_FIND);
  //  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::AperturePhot,this,ID_PHOT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::ProfilePhot,this,ID_PROFILE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::Matching,this,ID_MATCH);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::Astrometry,this,ID_ASTROMET);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::Stacking,this,ID_KOMBINE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::Deconvolution,this,ID_DEKON);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnAverage,this,ID_AVERAGE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::OnDarkbat,this,ID_DARKBAT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowser::HelpAbout,this,wxID_ABOUT);
  Bind(EVT_LIST_ACTIVATED,&MuniBrowser::OnActivated,this,ID_VIEW_FILE);
  Bind(EVT_META_PROGRESS,&MuniBrowser::OnMetaProgress,this);

  Bind(wxEVT_UPDATE_UI,&MuniBrowser::OnUpdateCut,this,wxID_CUT);
  Bind(wxEVT_UPDATE_UI,&MuniBrowser::OnUpdatePaste,this,wxID_PASTE);
  Bind(wxEVT_UPDATE_UI,&MuniBrowser::OnUpdateStop,this,wxID_STOP);
  Bind(wxEVT_UPDATE_UI,&MuniBrowser::OnUpdateOpen,this,wxID_OPEN);

  Bind(EVT_NAVIGATION,&MuniBrowser::OnEnterNavigator,this);

  list->ChangeDir(pwd,mask);
  navbar->ChangeDir(pwd,mask);
  SetTitle(pwd);

  wxAcceleratorEntry entries[2];
  entries[0].Set(wxACCEL_CTRL, (int) 'S', wxID_STOP);
  entries[1].Set(wxACCEL_CTRL, (int) 'R', wxID_REFRESH);
  wxAcceleratorTable accel(2, entries);
  SetAcceleratorTable(accel);

}

void MuniBrowser::OnClose(wxCloseEvent& event)
{
  if( metaload ) {
    wxASSERT(list);
    list->Stop();
    event.Veto();
    shutdown = true;
    return;
  }

  config->browser_size = GetSize();
  config->browser_fitsmask = mask;
  config->browser_toolbar = tbar->IsShown();
  config->browser_reverse = menuArrange->IsChecked(ID_SORT_REVERSE);

  Destroy();
}



void MuniBrowser::OnUpdateCut(wxUpdateUIEvent& event)
{
  //  event.Enable(list->GetSelectedMeta().size() > 0);
}

void MuniBrowser::OnUpdatePaste(wxUpdateUIEvent& event)
{
  //  event.Enable(list->GetClipboard().size() > 0);
}

void MuniBrowser::OnUpdateStop(wxUpdateUIEvent& event)
{
  event.Enable(metaload);
}

void MuniBrowser::OnUpdateOpen(wxUpdateUIEvent& event)
{
  event.Enable(metaload == false);
}

void MuniBrowser::OnDirOpen(wxCommandEvent& event)
{
  wxDirDialog select(this,"Choose a directory","",
		     wxDD_DEFAULT_STYLE|wxDD_DIR_MUST_EXIST);
  if (select.ShowModal() == wxID_OK ) {
    pwd = select.GetPath();
    list->ChangeDir(pwd,mask);
    navbar->ChangeDir(pwd,mask);
    SetTitle(pwd);
  }
}

// void MuniBrowser::FileOpen(wxCommandEvent& WXUNUSED(event))
// {
//   wxFileDialog select(this,"Choose a file",wxEmptyString,wxEmptyString,
// 		      "FITS files ("+config->dirmask+")|"+config->dirmask+
// 		      "|RAW files ("+config->rawmask+")|"+config->rawmask+
// 		      "|All files (*)|*",
// 		      wxFD_FILE_MUST_EXIST|wxFD_CHANGE_DIR|wxFD_MULTIPLE);

//   if (select.ShowModal() == wxID_OK ) {

//     wxArrayString files;
//     select.GetPaths(files);

//     if( select.GetFilterIndex() == 1 ) {
//       // RAW files

//       MuniImportRawOptions opt(this,config);
//       if( opt.ShowModal() == wxID_OK ) {
// 	MuniImportRaw *iraw = new MuniImportRaw(this,opt,files);
// 	iraw->Show();
// 	//	ImportRaw(opt,files);
//       }
//     }
//     else
//       FileLoad(files);
//   }
// }

// void MuniBrowser::ImportRaw(const MuniImportRawOptions& opt, const wxArrayString& files)
// {
//   MuniImportRaw *iraw = new MuniImportRaw(this,opt,files);
//   iraw->Show();
// }


void MuniBrowser::OnMetaProgress(MetaProgressEvent& event)
{
  //  wxLogDebug("MuniBrowser::OnMetaProgress()");
  wxASSERT(navbar);

  MuniArtIcons ico(wxART_TOOLBAR,wxSize(24,24));

  if( event.GetId() == ID_METALOAD_START ) {
    metaload = true;
    navbar->NewMoon();
    tbar->AddTool(wxID_STOP,"Stop",ico.Icon("process-stop"),"Stop loading");
    tbar->Realize();
  }
  else if( event.GetId() == ID_METALOAD_FINISH ) {
    metaload = false;
    navbar->HideMoon();
    tbar->RemoveTool(wxID_STOP);
    tbar->Realize();
    if( shutdown )
      Close();
  }
  else if( event.GetId() == ID_METALOAD_PHASE )
    navbar->PhaseMoon(event.phase);

}


void MuniBrowser::OnProperties(wxCommandEvent& event)
{
  MuniDirProperties *w = new MuniDirProperties(this,config,list->GetMetaList());
  w->Show();
}


void MuniBrowser::OnPreferences(wxCommandEvent& event)
{
  wxASSERT(!preferences);
  preferences = new MuniPreferences(this,config);
  Bind(wxEVT_CLOSE_WINDOW,&MuniBrowser::OnClosePreferences,this,preferences->GetId());
  preferences->SelectPage(1);
  preferences->Show();
}

void MuniBrowser::OnClosePreferences(wxCloseEvent& event)
{
  wxLogDebug("MuniBrowser::OnClosePreferences "+config->fits_key_filter);
  wxASSERT(preferences);
  Unbind(wxEVT_CLOSE_WINDOW,&MuniBrowser::OnClosePreferences,this,
	 preferences->GetId());
  preferences->Destroy();
  preferences = 0;

  wxASSERT(list);
  list->ChangeDir(pwd,mask);
}


void MuniBrowser::OnRefresh(wxCommandEvent& event)
{
  wxASSERT(list);
  list->ChangeDir(pwd,mask);
}

void MuniBrowser::OnStop(wxCommandEvent& event)
{
  wxASSERT(list);
  list->Stop();
}


void MuniBrowser::OnFileClose(wxCommandEvent& event)
{
  Close();
}

void MuniBrowser::OnSelectAll(wxCommandEvent& WXUNUSED(event))
{
  list->SelectAll();
}

void MuniBrowser::OnCut(wxCommandEvent& event)
{
  list->Cut();
}

void MuniBrowser::OnCopy(wxCommandEvent& event)
{
  list->Copy();
}

void MuniBrowser::OnPaste(wxCommandEvent& event)
{
  list->Paste();
}


void MuniBrowser::SelectItem(wxCommandEvent& event)
{
  int dir = 0;
  switch(event.GetId()){
  case wxID_BACKWARD: dir = -1; break;
  case wxID_FORWARD: dir = 1;   break;
  }

  list->SelectItemRelative(dir);
}

void MuniBrowser::OnNewBrowser(wxCommandEvent& WXUNUSED(event))
{
  MuniBrowser *newb = new MuniBrowser(NULL,config);
  newb->Show(true);
}


void MuniBrowser::OnNewView(wxCommandEvent& WXUNUSED(event))
{
  MuniView *view = new MuniView(NULL,config);
  view->Show();
}


void MuniBrowser::FindStars(wxCommandEvent& WXUNUSED(event))
{
  //  SetStatusText(wxT("FindStars is not implemenetd yet."));
}

void MuniBrowser::AperturePhot(wxCommandEvent& WXUNUSED(event))
{
  //  SetStatusText(wxT("AperturePhot is not implemenetd yet."));
}

void MuniBrowser::ProfilePhot(wxCommandEvent& WXUNUSED(event))
{
  //  SetStatusText(wxT("ProfilePhot is not implemenetd yet."));
}

void MuniBrowser::Matching(wxCommandEvent& WXUNUSED(event))
{
  //  SetStatusText(wxT("Matching is not implemenetd yet."));
}

void MuniBrowser::Astrometry(wxCommandEvent& WXUNUSED(event))
{
  //  SetStatusText(wxT("Astrometry is not implemenetd yet."));
  //  wxLogDebug("Selected files count: %d.",(int)list->GetSelectedMeta().size());
  /*
  MuniAstrometer ameter(this,config,list->GetSelectedMeta());
  ameter.ShowModal();
  */
}

void MuniBrowser::Stacking(wxCommandEvent& WXUNUSED(event))
{
  //  SetStatusText(wxT("Stacking is not implemenetd yet."));
}

void MuniBrowser::Deconvolution(wxCommandEvent& WXUNUSED(event))
{
  //  SetStatusText(wxT("Deconvolution is not implemenetd yet."));
}

void MuniBrowser::OnDarkbat(wxCommandEvent& WXUNUSED(event))
{
  //  wxLogDebug("Running darkbat...");
}


void MuniBrowser::OnAverage(wxCommandEvent& WXUNUSED(event))
{
  //  wxLogDebug("Preparing darks...");
}

void MuniBrowser::HelpAbout(wxCommandEvent& WXUNUSED(event))
{
  MuniAbout(config->munipack_icon);
}

void MuniBrowser::OnLabel(wxCommandEvent& event)
{
  wxASSERT(list);
  config->browser_labeltype = event.GetId();
  list->SelectLabel();
}

void MuniBrowser::OnSort(wxCommandEvent& event)
{
  wxASSERT(list);
  config->browser_sorttype = event.GetId();
  list->Sort();
}

void MuniBrowser::OnReverse(wxCommandEvent& event)
{
  wxASSERT(list);
  config->browser_reverse = event.IsChecked();
  list->Sort();
}

void MuniBrowser::OnActivated(wxNotifyEvent& event)
{
  OpenView(event.GetString());
}

void MuniBrowser::OpenView(const wxString& filename)
{
  MuniView *view = new MuniView(NULL,config);
  view->Show();
  view->LoadFile(filename);
}

void MuniBrowser::OnIconList(wxCommandEvent& event)
{
  MuniListCtrl *old_list = list;
  if( event.GetId() == wxID_VIEW_LARGEICONS ) {
    list = new MuniListIcon(this,ID_LIST,config);
    config->browser_iconlist = wxLC_ICON;
    menuView->Enable(menuView->FindItem(MENU_LABEL_LABELS),true);
  }
  else if( event.GetId() == wxID_VIEW_LIST ) {
    list = new MuniListList(this,ID_LIST,config);
    config->browser_iconlist = wxLC_REPORT;
    menuView->Enable(menuView->FindItem(MENU_LABEL_LABELS),false);
  }
  list->ChangeDir(pwd,mask);
  topsizer->Replace(old_list,list);
  Layout();
}

void MuniBrowser::OnEnterNavigator(NavigationEvent& event)
{
  //  wxLogDebug("OnEnter: "+event.pwd+" "+event.mask);
  pwd = event.pwd;
  mask = event.mask;
  list->ChangeDir(pwd,mask);
  SetTitle(pwd);
}


void MuniBrowser::FileLoad(const wxArrayString& files)
{
  list->FitsLoad(files);
}


void MuniBrowser::OnShowToolbar(wxCommandEvent& event)
{
  tbar->Show(event.IsChecked());
  Layout();
}
