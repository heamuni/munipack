/*

  xmunipack - FITS view


  Copyright © 2009-2014, 2017-22 F.Hroch (hroch@physics.muni.cz)

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


#include "help.h"
#include "tune.h"
#include "tuner.h"
#include "dispreview.h"
#include "caption.h"
#include "zoomer.h"
#include "export.h"
#include "view.h"
#include "fileprop.h"
#include <wx/wx.h>
#include <wx/toolbar.h>
#include <wx/txtstrm.h>
#include <wx/wfstream.h>
#include <wx/progdlg.h>
#include <wx/filename.h>
#include <wx/filefn.h>

#if wxUSE_FSWATCHER
#include <wx/fswatcher.h>
#endif


#define MENU_IMAGE "&Image"
#define MENU_HEAD "&Head"
#define MENU_TABLE "&Table"
#define MENU_TOOLS "&Tools"
#define MENU_EXTENSION "E&xtension"


using namespace std;


MuniView::MuniView(wxWindow *w, MuniConfig *c):
  wxFrame(w, wxID_ANY,"Untitled", wxDefaultPosition, c->view_size),
  config(c), loadfile(false), hdusel(0), zoom_init(false),
  switch_hdu(false),display_update(false),
  sloader(0), /* coloring(0),*/
  tune(0), preferences(0), preview(NULL), aphot(0), display_size(wxDefaultSize),
  image_size(wxDefaultSize)
#if wxUSE_FSWATCHER
  ,fswatch(0)
#endif
{
  SetIcon(config->munipack_icon);

  // menus
  menuFile = new wxMenu;
  menuFile->Append(wxID_NEW);
  menuFile->Append(wxID_OPEN);
  /*
  wxMenu *menuFileCreate = new wxMenu;
  menuFileCreate->Append(ID_COLOURING,"Colour image...");
  menuFileCreate->Append(wxID_ANY,"Image server ...");
  menuFileCreate->Append(wxID_ANY,"Artifical image...");
  menuFile->AppendSubMenu(menuFileCreate,"Create");
  */
  /*
  wxMenu *menuFileVO = new wxMenu;
  menuFileVO->Append(ID_CONE,"Cone search");
  menuFileVO->Append(ID_CONE,"Image search");
  menuFile->AppendSubMenu(menuFileVO,"Virtual observatory");
  */
  menuFile->AppendSeparator();
  menuFile->Append(wxID_SAVE);
  menuFile->Append(ID_EXPORT,"Export as...");
#ifdef __WXMAC__
  menuFile->Append(wxID_CLOSE);
#endif
  menuFile->AppendSeparator();
  menuFile->Append(wxID_PROPERTIES);
  /*
  menuFile->AppendSeparator();
  menuFile->Append(ID_PAGE_SETUP,"Page setup");
  menuFile->Append(wxID_PREVIEW);
  menuFile->Append(wxID_PRINT);
  */
#ifndef __WXMAC__
  menuFile->AppendSeparator();
  menuFile->Append(wxID_CLOSE);
#endif
#ifdef __WXMAC__
  menuFile->AppendSeparator();
  menuFile->Append(wxID_EXIT);
#endif

  wxMenu *menuEdit =  new wxMenu;
  /*
  menuEdit->Append(wxID_UNDO);
  menuEdit->Append(wxID_CUT);
  menuEdit->Append(wxID_COPY);
  menuEdit->Append(wxID_PASTE);
  */
#ifndef __WXMAC__
  menuEdit->AppendSeparator();
  menuEdit->Append(wxID_PREFERENCES);
#endif
  menuEdit->AppendSeparator();
  menuEdit->Append(wxID_STOP);

  menuView = new wxMenu;
  menuView->AppendCheckItem(ID_VIEW_TOOLBAR,"Show Toolbar",
    "Change visibility of toolbar (shorthand buttons with icons on top)");
  menuView->Check(ID_VIEW_TOOLBAR,config->view_toolbar);
  menuView->AppendCheckItem(ID_VIEW_CONTROLS,"Show control panels");
  menuView->Check(ID_VIEW_CONTROLS,config->view_controls);
  /*
  menuView->AppendSeparator();
  menuView->Append(ID_FULLSCREEN,"Fu&llscreen",
		   "Enable fullscreen image display");
  */

  wxMenu *menuHelp = new wxMenu;
  //  menuHelp->AppendCheckItem(ID_LOG,"&Log");
  menuHelp->Append(wxID_ABOUT);

  // to be determined by actual FITS file
  menuExt = new wxMenu();
  menuExt->AppendSeparator();

  wxMenuBar *menuBar = new wxMenuBar;
  menuBar->Append(menuFile,"&File");
  menuBar->Append(menuEdit,"&Edit");
  menuBar->Append(menuView,"&View");
  menuBar->Append(menuExt,MENU_EXTENSION);
  menuBar->Append(menuHelp,"&Help");
  SetMenuBar(menuBar);

  // toolbars
  wxToolBar *tbar = CreateToolBar(wxTB_HORIZONTAL|wxTB_TEXT);
  SetupToolbar_init();
  tbar->Show(config->view_toolbar);
  SetToolBar(tbar);

  // Extension list
  extlist = new MuniExtensionList(this,config);

  // magnifier
  magnifier = new MuniMagnifier(this,config);
  magnifier->Show(false);

  // zoom panel
  zoomview = new MuniZoomWindow(this,config);
  zoomview->Show(false);

  // zoom slider
  zoomer = new MuniZoomSlider(this,ID_ZOOMER);
  zoomer->Show(false);

  // tune mini-panel
  tuner = new MuniTuner(this,ID_TUNER);
  tuner->Show(false);

  // panel settings
  /*
  wxButton *poptions = new wxButton(this,wxID_PREFERENCES,L"",wxDefaultPosition,
				    wxDefaultSize,wxBU_EXACTFIT|wxBORDER_NONE);
  */

  panelsizer = new wxBoxSizer(wxVERTICAL);
  panelsizer->SetMinSize(wxSize(18*GetCharWidth(),2*GetCharHeight()));
  //panelsizer->SetMinSize(wxSize(config->icon_size+2,-1));
  panelsizer->Add(extlist,wxSizerFlags(1).Expand());
  panelsizer->Add(zoomview,wxSizerFlags().Centre().Expand().Shaped());
  panelsizer->Add(zoomer,wxSizerFlags().Expand());
  panelsizer->Add(tuner,wxSizerFlags().Expand());

  // the workplace
  // workplace filling
  splash = new MuniSplashing(this,config);

  // caption, display only
  caption = new MuniDisplayCaption(this,config);
  caption->Show(config->view_controls);

  placesizer = new wxBoxSizer(wxVERTICAL);
  placesizer->Add(splash,wxSizerFlags(1).Expand());
  placesizer->Add(caption,wxSizerFlags().Expand());

  wxBoxSizer *topsizer = new wxBoxSizer(wxHORIZONTAL);
  topsizer->Add(panelsizer,wxSizerFlags().Expand().Border(wxLEFT,1));
  topsizer->Add(placesizer,wxSizerFlags(1).Expand().Border(wxLEFT|wxRIGHT|wxBOTTOM,1));
  topsizer->Show(panelsizer,config->view_controls);

  SetSizer(topsizer);
  SetSizeHints(wxSize(600,500));

  menuFile->Enable(wxID_SAVE,false);
  menuFile->Enable(ID_EXPORT,false);
  /*
  menuFile->Enable(ID_PAGE_SETUP,false);
  menuFile->Enable(wxID_PREVIEW,false);
  menuFile->Enable(wxID_PRINT,false);
  */

  /*
  menuEdit->Enable(wxID_UNDO,false);
  menuEdit->Enable(wxID_CUT,false);
  menuEdit->Enable(wxID_COPY,false);
  menuEdit->Enable(wxID_PASTE,false);
  */

  //  menuView->Enable(ID_FULLSCREEN,false);
  menuFile->Enable(wxID_PROPERTIES,false);
  //  menuFileCreate->Enable(ID_COLOURING,false);

  Bind(wxEVT_CLOSE_WINDOW,&MuniView::OnClose,this);
  Bind(wxEVT_IDLE,&MuniView::OnIdle,this);
  Bind(EVT_SIZE_CHANGED,&MuniView::OnDisplaySize,this);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::FileClose,this,wxID_CLOSE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::FileClose,this,wxID_EXIT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::NewView,this,wxID_NEW);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::FileOpen,this,wxID_OPEN);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::FileSave,this,wxID_SAVE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::FileExport,this,ID_EXPORT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnFileProperties,this,wxID_PROPERTIES);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnPreferences,this,wxID_PREFERENCES);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnConeSearch,this,ID_CONE);
  Bind(EVT_FINISH_DIALOG,&MuniView::OnCloseCone,this,ID_CONE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnShowToolbar,this,ID_VIEW_TOOLBAR);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnShowControls,this,ID_VIEW_CONTROLS);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnFullScreen,this,ID_FULLSCREEN);
  //  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::Coloring,this,ID_COLORING);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::HelpAbout,this,wxID_ABOUT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnMenuExt,this,wxID_BACKWARD);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnMenuExt,this,wxID_FORWARD);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnStop,this,wxID_STOP);

  Bind(EVT_DRAW,&MuniView::OnDraw,this);
  Bind(EVT_CLICK,&MuniView::OnClick,this);
  Bind(EVT_SLEW,&MuniView::OnMouseMotion,this);
  Bind(EVT_TUNE,&MuniView::OnTuner,this,ID_TONE_QBLACK);
  Bind(EVT_TUNE,&MuniView::OnTuner,this,ID_TONE_RSENSE);
  Bind(EVT_ZOOM,&MuniView::OnZoomer,this,ID_ZOOMER);
  Bind(EVT_ZOOM,&MuniView::OnZoomview,this,ID_ZOOMVIEW);
  Bind(EVT_RENDER,&MuniView::OnRender,this);

  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnExtChanged,this,ID_EXTENSION_SELECT);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateStop,this,wxID_STOP);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateBackward,this,wxID_BACKWARD);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateForward,this,wxID_FORWARD);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateShowLegend,this,ID_LEGEND);
  //  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateShowSources,this,ID_SOURCES);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateShowTune,this,ID_TUNE);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateZoomer,this,zoomer->GetId());
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateFullScreen,this,ID_FULLSCREEN);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateZoom100,this,wxID_ZOOM_100);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateZoomIn,this,wxID_ZOOM_IN);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateZoomOut,this,wxID_ZOOM_OUT);
  Bind(wxEVT_UPDATE_UI,&MuniView::OnUpdateExportAs,this,ID_EXPORT);

  // Binds common to all extensions
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnHeader,this,ID_INFO);

  // Binds for Image
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnMenuZoom,this,wxID_ZOOM_100);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnMenuZoom,this,wxID_ZOOM_IN);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnMenuZoom,this,wxID_ZOOM_OUT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnTune,this,ID_TUNE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnShowLegend,this,ID_LEGEND);
  //  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnShowSources,this,ID_SOURCES);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnFind,this,ID_FIND);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnAphot,this,ID_APHOT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnAstrometry,this,ID_ASTROMETRY);
  Bind(EVT_FITS_STREAM,&MuniView::OnLoadFits,this);

  wxAcceleratorEntry entries[3];
  // The keys on regular keyboard are supporred directly.
  entries[0].Set(wxACCEL_CTRL, WXK_NUMPAD0, wxID_ZOOM_100);
  entries[1].Set(wxACCEL_CTRL, WXK_NUMPAD_ADD, wxID_ZOOM_IN);
  entries[2].Set(wxACCEL_CTRL, WXK_NUMPAD_SUBTRACT, wxID_ZOOM_OUT);
  wxAcceleratorTable accel(3, entries);
  SetAcceleratorTable(accel);
}

MuniView::~MuniView()
{
  wxASSERT(sloader == 0);

#if wxUSE_FSWATCHER
  if( fswatch ) {
    fswatch->SetOwner(0);
    Unbind(wxEVT_FSWATCHER, &MuniView::OnFileSystemEvent, this);
    delete fswatch;
  }
#endif

  wxLogDebug("MuniView::~MuniView()");
}

void MuniView::OnIdle(wxIdleEvent& event)
{
  /*
   wxLogDebug("MuniView::OnIdle %d %d %d %d, %d %d %p %ld",display_size.GetWidth(),
	      display_size.GetHeight(),display_size.IsFullySpecified(),zoom_init,
	      event.GetId(),GetId(),event.GetEventObject(),event.GetTimestamp());
  */

  if( display_update ) {
    display_update = false;
    SetupZoom();
  }

  if( switch_hdu ) {
    wxLogDebug("MuniView::OnIdle() .. switch_hdu");
    switch_hdu = false;

    if( fits.Hdu(hdusel).IsDisplayImplemented() ) {
      wxASSERT(zoomview && zoomer && GetDisplay() && GetDisplay()->GetZoom() > 0);
      double zoom = GetDisplay()->GetZoom();
      zoomer->SetZoom(zoom);
      zoomview->SetDisplayZoom(zoom);
    }
  }


#if wxUSE_FSWATCHER
  // File monitor -- re-load
  if( loadfile ) {
    wxString filename(fits.GetFullPath());
    double dt = MonitorDelay(filename);
    if( dt < 0 ) return;
    if( dt < 0.1 ) {
      event.RequestMore();
      return;
    }

    loadfile = false;
    wxLogDebug("Re-loading "+filename);
    LoadFile(filename);
  }
#endif
}

void MuniView::SetupZoom()
{
  if( zoom_init == false ) return;
  zoom_init = false;

  wxASSERT(display_size.IsFullySpecified() && image_size.IsFullySpecified());
  float zoom = GetBestFitZoom(display_size,image_size);

  wxLogDebug("MuniView::SetupZoom %f %dx%d -> %dx%d",zoom,
	     image_size.GetWidth(),image_size.GetHeight(),
	     display_size.GetWidth(),display_size.GetHeight());

  if( preview )
    preview->SetZoom(zoom,image_size);

  if( sloader )
    sloader->SetZoom(zoom);

  wxASSERT(zoomer && zoomview);
  zoomer->SetZoom(zoom);
  zoomview->SetDisplayZoom(zoom);

  MuniDisplay *display = GetDisplay();
  if( display )
    display->SetZoom(zoom);
}


void MuniView::OnDisplaySize(MuniSizeChangedEvent& event)
{
  wxLogDebug("MuniView::OnDisplaySize");

  display_update = true;
  display_size = event.size;

  wxASSERT(zoomview);
  zoomview->SetDisplaySize(display_size);
}



void MuniView::OnClick(MuniClickEvent& event)
{
  //  wxLogDebug("MuniView::OnClick");
  if( aphot )
    wxQueueEvent(aphot,event.Clone());
}

void MuniView::OnClose(wxCloseEvent& event)
{
  wxLogDebug("MuniView::OnClose");

  if( event.CanVeto() && IsModified() ) {
    if( Unsaved(fits.GetName()) == wxID_CANCEL ) {
      event.Veto();
      return;
    }
  }

  // if( ! backup.IsEmpty() ) {
  //   wxBusyCursor wait;
  //   for(size_t i = 0; i < backup.GetCount(); i++)
  //     wxRemoveFile(backup[i]);
  // }

  if( sloader ) {
    sloader->StopAndWait();
    sloader->Wait();
    delete sloader;
    sloader = 0;
  }

  config->view_toolbar = GetToolBar()->IsShown();
  config->view_size = GetSize();

  event.Skip();
}



void MuniView::CreateFSWatch()
{
  wxLogDebug("MuniView::CreateFSWatch");

#if wxUSE_FSWATCHER

  if( fswatch ) return;

  fswatch = new wxFileSystemWatcher();
  fswatch->SetOwner(this);
  Bind(wxEVT_FSWATCHER, &MuniView::OnFileSystemEvent, this);
#endif
}


void MuniView::MonitorFile(const wxString& filename)
{
#if wxUSE_FSWATCHER

  const int flags = wxPATH_GET_VOLUME | wxPATH_GET_SEPARATOR;

  if( !fswatch ) return;

  wxFileName file = wxFileName::FileName(filename);
  wxString path(file.GetPath(flags));
  if( path.IsEmpty() ) {
    file.AssignCwd();
    path = file.GetPath(flags);
  }

  fswatch->RemoveAll();
  fswatch->Add(path);
  wxLogDebug("Monitoring: "+filename+" in "+path);
#endif
}

#if wxUSE_FSWATCHER
void MuniView::OnFileSystemEvent(wxFileSystemWatcherEvent& event)
{

  wxString fitsname(fits.GetFullPath());
  wxString eventname(event.GetPath().GetFullPath());

  // We're monitoring MODIFY events on our file. Multiply events can be
  // repored, perhaps for every file-block update; we are launch trigger,
  // which calls LoadFile in Idle time, some time later.
  if( fitsname == eventname && event.GetChangeType() == wxFSW_EVENT_MODIFY ) {
    loadfile = true;
    wxWakeUpIdle();
  }
  //wxLogDebug("MuniView::OnFileSystemEvent: "+eventname+": "+event.ToString());

  // malfunction, sometimes reported with the empty string description
  if( event.GetChangeType() == wxFSW_EVENT_WARNING &&
      event.GetErrorDescription() != "" )
    wxLogWarning(event.GetErrorDescription());

  if( event.GetChangeType() == wxFSW_EVENT_ERROR &&
      event.GetErrorDescription() != "" ) {
    wxLogError(event.GetErrorDescription());
    fswatch->RemoveAll();
  }
}
#endif

double MuniView::MonitorDelay(const wxString& filename) const
{
  // return delay in seconds since last modification of the filename
  time_t tt = wxFileModificationTime(filename);
  if( tt < 0 )
    return -1.0;

  wxDateTime tmod(tt);
  wxDateTime t = wxDateTime::Now();
  //  wxLogDebug("MuniView::MonitorDelay %f",86400*(t.GetMJD()-tmod.GetMJD()));
  return 86400.0*(t.GetMJD() - tmod.GetMJD());
}



bool MuniView::IsModified() const
{
  return ! backup.IsEmpty();
}

void MuniView::RemoveBackup()
{
  for(size_t i = 0; i < backup.GetCount(); i++) {
    if( wxFileExists(backup[i]) )
      wxRemoveFile(backup[i]);
  }
  backup.Clear();
}


void MuniView::OnStop(wxCommandEvent& event)
{
  wxASSERT(sloader);
  sloader->Stop();
}


void MuniView::LoadFile(const wxString& filename)
{
  //  wxASSERT(backup.IsEmpty());
  //  RemoveBackup();

  if( IsModified() ) {
    if( Unsaved(fits.GetName()) == wxID_CANCEL )
      return;
  }

  // cleanup GUI
  DeletePlaces();
  //  backup.Clear();

  // setup toolbar
  ClearToolbar();
  SetupToolbar_load();

  extlist->Clear();
  wxASSERT(preview == NULL);
  splash->Play();
  zoom_init = true;
  image_size = wxDefaultSize;

  wxASSERT(sloader == 0);
  sloader = new FitsStream(this,filename,config->icon_size);
  if( sloader ) {
    wxThreadError code = sloader->Create();
    wxASSERT(code == wxTHREAD_NO_ERROR);
    sloader->Run();
  }
  else
    wxLogError("MuniView::LoadFile FITS load thread failed to start.");


#if wxUSE_FSWATCHER
  if( fswatch )
    fswatch->RemoveAll();
#endif
}

void MuniView::LoadFileBackup(const wxString& file, const wxString& b)
{
  // NEEDS CHECK
  wxString a;
  a.Printf(b+"_%d",(int)backup.GetCount());
  backup.Add(a);

  wxFAIL_MSG("MuniView::LoadFileBackup SCHEDULED FOR RE-IMPLEMENT");
}

int MuniView::Unsaved(const wxString& filename)
{
  wxMessageDialog dialog(this,"Do you want to save changes you made in \""+
			 filename+"\" ?",filename,wxICON_EXCLAMATION |
			 wxYES_NO | wxCANCEL | wxNO_DEFAULT);
  dialog.SetExtendedMessage("Your changes will be lost if you don't save them.");
  dialog.SetYesNoLabels("Save","Don't Save");

  int code = dialog.ShowModal();

  if( code == wxID_YES ) {
    wxCommandEvent e;
    FileSave(e);
  }
  else if( code == wxID_NO ) {
    if( IsModified() ) {
      wxBusyCursor wait;
      RemoveBackup();
    }
  }
  else if( code == wxID_CANCEL )
    ;

  return code;
}


void MuniView::OnLoadFits(FitsStreamEvent& event)
{
  //wxLogDebug("MuniView::OnLoadFits %d %d",int(event.GetId()),int(ID_FITSLOAD_OPEN));

  if( event.GetId() == ID_FITSLOAD_OPEN ) {
    wxLogDebug("MuniView::OnLoadFits = ID_FITSLOAD_OPEN");
    wxFileName name(event.filename);
    SetTitle(name.GetName()+" ...");
    extlist->Init();
    splash->Stop();
    splash->Show(false);
    Layout();
  }

  else if( event.GetId() == ID_FITSLOAD_HEAD ) {
    wxLogDebug("ID_FITSLOAD_HEAD extname=`"+event.extname+"' type: "+
	       FitsHeader::GetType_str(event.hdutype)+" ("+
	       FitsHeader::GetFlavour_str(event.flavour)+")");

    if( preview == 0 && event.chdu == 1 && event.hdutype == HDU_IMAGE &&
	(event.flavour == HDU_IMAGE_FRAME || event.flavour == HDU_IMAGE_COLOUR)	) {
      preview = new MuniDisplayPreview(this,config);
      placesizer->Prepend(preview,wxSizerFlags(1).Expand());
    }

    // hide live-zoom
    magnifier->Show(false);

    wxString label(event.extname);
    extlist->Append(label,event.hdutype);
    Layout();
  }

  else if( event.GetId() == ID_FITSLOAD_IMGHEAD ) {
    vector<long> naxes = event.naxes;
    wxLogDebug("IMGHEAD: BITPIX=%d, NAXIS=%d, NAXES=%ldx%ld, cHDU=%d",
	       event.bitpix,event.naxis,naxes[0],naxes[1],event.chdu);
    if( event.chdu == 1 && preview ) {
      wxASSERT(preview);
      preview_size = wxSize(int(naxes[0]),int(naxes[1]));
      image_size = wxSize(int(naxes[0]),int(naxes[1]));
      display_update = true;
    }
  }

  else if( event.GetId() == ID_FITSLOAD_IMAGE ) {
    //    wxLogDebug("LOAD_IMAGE: cHDU:%d",event.chdu);
    extlist->UpdateProgress(event.chdu,event.progress);
    //    wxQueueEvent(extlist,event.Clone());
  }

  else if( event.GetId() == ID_FITSLOAD_IMGPREVIEW ) {
    //    wxLogDebug("LOAD_IMAGE PREVIEW: cHDU:%d",event.chdu);

    if( event.chdu == 1 && preview )
      preview->RefreshCanvas(event.crow,event.preview);
  }


  else if( event.GetId() == ID_FITSLOAD_TBLHEAD ) {

    wxLogDebug("TBLHEAD: %ldx%d",event.nrows,event.ncols);
  }

  else if( event.GetId() == ID_FITSLOAD_TABLE ) {
    int n = event.crow - 100*(event.crow / 100);
    if( n == 0 )
      wxLogDebug("MuniView::OnLoadFits = ID_FITSLOAD_TABLE, rows=%ld",event.crow);
    // display current state
    extlist->UpdateProgress(event.chdu,event.progress);
  }

  else if( event.GetId() == ID_FITSLOAD_HDUFIN ) {
    extlist->FinishProgress(event.chdu);
  }

  else if( event.GetId() == ID_FITSLOAD_FINISH ) {
    wxLogDebug("MuniView::OnLoadFits = ID_FITSLOAD_FINISH");

    sloader->Wait();
    delete sloader;
    sloader = 0;

    wxImage preimage;
    if( preview ) {
      preimage = preview->GetCanvas();
      placesizer->Detach(preview);
      preview->Destroy();
      preview = 0;
    }
    extlist->Finish();

    fits = event.fitsfile;
    icons = event.icons;

    if( fits.Status() ) {

      wxASSERT(fits.HduCount() > 0);
      SetFits(fits,event.tones);
      hdusel = 0;
      SwitchHdu(hdusel);
      MonitorFile(fits.GetFullPath());
      SetRepresentedFilename(fits.GetFullPath());
      MuniDisplay *display = GetDisplay();
      if( display && preimage.IsOk() ) {
	display->SetPreview(preimage);
      }
    }
    else {

      // corrupted file
      wxArrayString es = fits.GetErrorMessage();
      for(size_t i = 0; i < es.size(); i++)
	wxLogError(es[i]);
      wxLogError("Loading failed on `"+fits.GetErrorDescription()+"'");

      // CLEANUP
    }

    menuFile->Enable(wxID_PROPERTIES,fits.IsOk());
  }

  else if( event.GetId() == ID_FITSLOAD_CRASH ||
	   event.GetId() == ID_FITSLOAD_INTERRUPT ) {
    wxLogDebug("MuniView::OnLoadFits = ID_FITSLOAD_CRASH or INTERRUPT");

    sloader->Wait();
    delete sloader;
    sloader = 0;

    for(size_t i = 0; i < event.errlog.GetCount(); i++)
      wxLogError(event.errlog[i]);
    wxLogError(event.errmsg);

    // CLEANUP
    extlist->Clear();
    splash->Stop();
    ClearToolbar();
    SetupToolbar_init();
  }

}


void MuniView::SetFits(const FitsFile& fits, const vector<FitsTone>& tones)
{
  wxASSERT(fits.IsOk() &&  places.size() == 0);

  wxSizerFlags sizer_flags(1);
  sizer_flags.Expand();

  for(size_t i = 0; i < fits.HduCount(); i++) {

    if( fits.Hdu(i).Type() == HDU_IMAGE ) {
      if( fits.Hdu(i).IsDisplayImplemented() ) {
	MuniDisplay *display = new MuniDisplay(this,config);
	display->SetArray(FitsArray(fits.Hdu(i)),tones[i]);
	places.push_back(display);
	placesizer->Prepend(display,sizer_flags);
	display->Show(false);
      }
      else {
	MuniDummy *dummy = new MuniDummy(this,config);
	places.push_back(dummy);
	placesizer->Prepend(dummy,sizer_flags);
	dummy->Show(false);
      }
    }
    else if( fits.Hdu(i).Type() == HDU_TABLE ) {
      MuniGrid *grid = new MuniGrid(this,config);
      grid->SetHdu(fits.Hdu(i));
      places.push_back(grid);
      placesizer->Prepend(grid,sizer_flags);
      grid->Show(false);
    }
    else if( fits.Hdu(i).Type() == HDU_HEAD ) {
      MuniHead *head = new MuniHead(this,config);
      head->SetHdu(fits.Hdu(i));
      places.push_back(head);
      placesizer->Prepend(head,sizer_flags);
      head->Show(false);
    }
  }

  // menu
  ReplaceExtensionMenu(fits);

  // title
  SetupTitle();

}


void MuniView::ReplaceExtensionMenu(const FitsFile& fits)
{
  wxASSERT(fits.IsOk() && menuExt);

  menuExt = new wxMenu;

  wxMenuBar *menuBar = GetMenuBar();
  wxMenu *oldmenu = menuBar->Replace(menuBar->FindMenu(MENU_EXTENSION),
				     menuExt,MENU_EXTENSION);
  wxASSERT(oldmenu);

  while( oldmenu->GetMenuItemCount() > 0 ) {
    int i = oldmenu->GetMenuItemCount() - 1;
    wxMenuItem *item = oldmenu->FindItemByPosition(i);
    wxASSERT(item);
    int id = item->GetId();
    Unbind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnMenuExt,this,id);
    oldmenu->Destroy(id);
  }
  delete oldmenu;
  viewid.clear();

  // create the Extension menu
  for(size_t k = 0; k < fits.HduCount(); k++ ) {

    wxString label = fits.Hdu(k).GetControlLabel();

    if( k < 9 ) {
      wxString l;
      l.Printf("\tAlt+%d",(int) k+1);
      label += l;
    }
    wxMenuItem* item = menuExt->AppendRadioItem(wxID_ANY,label);
    wxASSERT(item);
    int id = item->GetId();
    Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnMenuExt,this,id);
    viewid.push_back(id);
  }

  menuExt->AppendSeparator();
  menuExt->Append(wxID_BACKWARD,"Previous\tAlt+Left");
  menuExt->Append(wxID_FORWARD,"Next\tAlt+Right");
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnMenuExt,this,wxID_BACKWARD);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniView::OnMenuExt,this,wxID_FORWARD);

  menuFile->Enable(wxID_SAVE,false);
  menuFile->Enable(ID_EXPORT,true);
}

void MuniView::DeletePlaces()
{
  for(size_t i = 0; i < places.size(); i++) {
    wxWindow *win = places[i];
    wxASSERT(win);
    bool b = placesizer->Detach(win);
    wxASSERT(b);
    win->Destroy();
  }
  places.clear();
}


void MuniView::OnShowToolbar(wxCommandEvent& event)
{
  bool show = event.IsChecked();
  GetToolBar()->Show(show);
  config->view_toolbar = show;
}

void MuniView::OnShowControls(wxCommandEvent& event)
{
  bool show = event.IsChecked();
  GetSizer()->Show(panelsizer,show,true);
  magnifier->Show(false);
  caption->Show(show);
  config->view_controls = show;
  Layout();
}


void MuniView::FileClose(wxCommandEvent& WXUNUSED(event))
{
  Close();
}

void MuniView::FileOpen(wxCommandEvent& WXUNUSED(event))
{
  if( IsModified() ) {
    if( Unsaved(fits.GetName()) == wxID_CANCEL )
      return;
  }

  wxFileDialog select(this,"Choose a file",wxEmptyString,wxEmptyString,
		      "FITS files "+config->dirmask+")|"+
		      config->dirmask+"| All files (*)|*",
		      wxFD_FILE_MUST_EXIST|wxFD_CHANGE_DIR);
  if (select.ShowModal() == wxID_OK )
    LoadFile(select.GetPath());
}

void MuniView::SaveFile(const wxString& savename)
{
  bool s = fits.Save("!"+savename);
  if( s ) {
    SetupTitle();
    RemoveBackup();
  }

}

void MuniView::FileSave(wxCommandEvent& WXUNUSED(event))
{
  wxFileDialog select(this,"Choose a file",fits.GetPath(),fits.GetFullName(),
		      "FITS files "+config->dirmask+")|"+
		      config->dirmask+"| All files (*)|*",
		      wxFD_SAVE|wxFD_OVERWRITE_PROMPT|wxFD_CHANGE_DIR);

  if (select.ShowModal() == wxID_OK ) {
    wxBusyCursor wait;
    SaveFile(select.GetPath());
  }
}



void MuniView::FileExport(wxCommandEvent& WXUNUSED(event))
{
  wxASSERT(fits.IsOk());

  // exporting images
  if( fits.Hdu(hdusel).IsDisplayImplemented() ) {
    wxASSERT(GetDisplay());

    MuniExportDialog eselect(this,"Export As","",fits.GetName(),
			    "PNG files (*.png)|*.png|JPEG files (*.jpg)|*.jpg|"
			    "TIFF files (*.tif)|*.tif|PNM files (*.pnm)|*.pnm",
			    wxFD_SAVE|wxFD_OVERWRITE_PROMPT|wxFD_CHANGE_DIR);
    wxBusyCursor wait;
    if( eselect.ShowModal() == wxID_OK )
      eselect.Save(FitsArray(fits.Hdu(hdusel)),GetDisplay());
  }

  // exporting tables
  else if( GetHduType() == HDU_TABLE ) {

    wxFileDialog select(this,"Export As",wxEmptyString,fits.GetName(),
			"TEXT files (*.txt)|*.txt",
			wxFD_SAVE|wxFD_OVERWRITE_PROMPT|wxFD_CHANGE_DIR);

    wxBusyCursor wait;
    if (select.ShowModal() == wxID_OK ) {

      /*
      wxProgressDialog dialog(_("Export of Table"),_("Exporting ")+
			      select.GetFilename()+_(" ... "),100,this,
			      wxPD_APP_MODAL|wxPD_AUTO_HIDE);
      dialog.Pulse();
      */


      wxFileOutputStream output(select.GetPath());
      wxTextOutputStream cout(output);

      const FitsTable table(fits.Hdu(/*HduSel()*/hdusel));
      cout << "#";
      for(int i = 0; i < table.Width(); i++) {
	wxString key;
	key.Printf("TTYPE%d",(int) i+1);
	cout << " " << table.GetKey(key);
      }
      cout << endl;

      for(int j = 0; j < table.Height(); j++) {
	cout << j;
	for(int i = 0; i < table.Width(); i++)
	  cout << " " << table.Pixel_str(i,j) ;
	cout << endl;
	//	if( j % 100 == 0) dialog.Pulse();
      }
    }
  }

  // exporting headers
  else if( GetHduType() == HDU_HEAD ) {

    wxFileDialog select(this,"Export As",wxEmptyString,fits.GetName(),
			"TEXT files (*.txt)|*.txt",
			wxFD_SAVE|wxFD_OVERWRITE_PROMPT|wxFD_CHANGE_DIR);

    wxBusyCursor wait;
    if (select.ShowModal() == wxID_OK ) {

      /*
      wxProgressDialog dialog(_("Export of Header"),_("Exporting ")+
			      select.GetFilename()+_(" ... "),100,this,
			      wxPD_APP_MODAL|wxPD_AUTO_HIDE);
      dialog.Pulse();
      */

      wxFileOutputStream output(select.GetPath());
      wxTextOutputStream cout(output);

      const FitsHdu head = fits.Hdu(hdusel);
      for(size_t i = 0; i < head.GetCount(); i++)
	cout << head.Item(i) << endl;
    }
  }

  else
    wxFAIL_MSG("----- WARNING: Unreachable code.");

}


void MuniView::OnFileProperties(wxCommandEvent& WXUNUSED(event))
{
  MuniIcon micon(config);
  std::vector<wxImage> ilist(micon.MergeList(fits,icons));

  wxASSERT(ilist.size() > 0);
  wxImage cover = micon.GetCover(ilist[0]);
  FitsMeta meta(fits,cover,ilist);

  MuniFileProperties *w = new MuniFileProperties(this,config,meta,cover,ilist);
  w->Show();
}

void MuniView::OnPreferences(wxCommandEvent& WXUNUSED(event))
{
  wxASSERT(!preferences);
  preferences = new MuniPreferences(this,config);
  Bind(wxEVT_CLOSE_WINDOW,&MuniView::OnClosePreferences,this,preferences->GetId());
  Bind(EVT_CONFIG_UPDATED,&MuniView::OnPreferencesEvent,this);
  preferences->Show();
}

void MuniView::OnClosePreferences(wxCloseEvent& event)
{
  wxLogDebug("MuniView::OnClosePreferences "+config->fits_key_filter);
  wxASSERT(preferences);
  Unbind(wxEVT_CLOSE_WINDOW,&MuniView::OnClosePreferences,this,preferences->GetId());
  Unbind(EVT_CONFIG_UPDATED,&MuniView::OnPreferencesEvent,this);
  preferences->Destroy();
  preferences = 0;

  wxASSERT(caption);
  MuniDisplay *display = GetDisplay();
  if( display ) {
    caption->ConfigUpdate();
    display->ConfigUpdate();
  }
}


void MuniView::OnPreferencesEvent(MuniConfigEvent& event)
{
  if( event.GetId() == ID_PREFERENCES_COOTYPE ||
      event.GetId() == ID_PREFERENCES_QPHTYPE ) {

    wxASSERT(caption);
    MuniDisplay *display = GetDisplay();
    if( display ) {
      caption->ConfigUpdate();
      display->ConfigUpdate();
    }
  }
}

double MuniView::GetBestFitZoom(const wxSize& display_size, const wxSize& image_size) const
{
  wxASSERT(display_size.IsFullySpecified() && image_size.IsFullySpecified());
  MuniZoomSet set;
  double z = set.SetBestFitZoom(display_size,image_size);
  return z > 1 ? 1.0 : z;
}


void MuniView::OnMenuZoom(wxCommandEvent& event)
{
  wxASSERT(zoomview && zoomer && GetDisplay());

  double zoom = GetDisplay()->GetZoom();
  MuniZoomSet set(zoom);

  switch(event.GetId()) {
  case wxID_ZOOM_100: zoom = 1.0; break;
  case wxID_ZOOM_IN:  zoom = set.Increase(); break;
  case wxID_ZOOM_OUT: zoom = set.Decrease(); break;
  }
  zoomview->SetDisplayZoom(zoom);
  zoomer->SetZoom(zoom);
  GetDisplay()->SetZoom(zoom);
}

void MuniView::OnZoomer(MuniZoomEvent& event)
{
  wxASSERT(zoomview && GetDisplay() );

  // redirect to ZoomView
  wxQueueEvent(zoomview,event.Clone());

  // redirect to display
  wxQueueEvent(GetDisplay(),event.Clone());
}

void MuniView::OnZoomview(MuniZoomEvent& event)
{
  MuniDisplay *display = GetDisplay();
  if( display )
    wxQueueEvent(display,event.Clone());
}

void MuniView::NewView(wxCommandEvent& WXUNUSED(event))
{
  MuniView *w = new MuniView(0,config);
  w->Show();
}

void MuniView::OnFullScreen(wxCommandEvent& WXUNUSED(event))
{
  wxASSERT(GetDisplay());

  if( IsFullScreen() ) {
    bool show = config->view_controls;
    caption->Show(show);
    extlist->Show(show);
    zoomview->Show(show);
    magnifier->Show(show);
    zoomer->Show(show);
    tuner->Show(show);
    GetToolBar()->Show(config->view_toolbar);
    GetSizer()->Show(panelsizer,config->view_controls && fits.HduCount() > 1,true);
    ShowFullScreen(false);

  }
  else {

    GetToolBar()->Show(false);
    caption->Show(false);
    extlist->Show(false);
    zoomview->Show(false);
    magnifier->Show(false);
    zoomer->Show(false);
    tuner->Show(false);
    GetSizer()->Show(panelsizer,false,true);
    ShowFullScreen(true);
  }

  Layout();
}


void MuniView::OnTune(wxCommandEvent& event)
{
  wxASSERT(!tune && GetDisplay());

  wxPoint point = GetScreenRect().GetTopLeft();

  FitsArray array(fits.Hdu(hdusel));
  FitsTone tone = GetDisplay()->GetTone();
  FitsItt itt = GetDisplay()->GetItt();
  if( array.IsColour() )
    tune = new MuniTune(this,wxID_ANY,point,wxDefaultSize,
			config->icon_size,array,tone,itt,
			GetDisplay()->GetColour());
  else
    tune = new MuniTune(this,wxID_ANY,point,wxDefaultSize,
			config->icon_size,array,tone,itt,
			GetDisplay()->GetPalette());

  Bind(wxEVT_CLOSE_WINDOW,&MuniView::OnCloseTune,this,tune->GetId());
  Bind(EVT_TUNE_SCALE,&MuniView::OnTuneScale,this);
  Bind(EVT_TUNE_ITT,&MuniView::OnTuneItt,this);
  Bind(EVT_TUNE_PAL,&MuniView::OnTunePal,this);
  Bind(EVT_TUNE_COLOUR,&MuniView::OnTuneColour,this);
  Bind(EVT_TUNE_NITE,&MuniView::OnTuneNite,this);
  tune->Show(true);
}

void MuniView::OnCloseTune(wxCloseEvent& event)
{
  wxLogDebug("MuniView::OnCloseTune");
  wxASSERT(tune);
  Unbind(EVT_TUNE_NITE,&MuniView::OnTuneNite,this);
  Unbind(EVT_TUNE_COLOUR,&MuniView::OnTuneColour,this);
  Unbind(EVT_TUNE_PAL,&MuniView::OnTunePal,this);
  Unbind(EVT_TUNE_ITT,&MuniView::OnTuneItt,this);
  Unbind(EVT_TUNE_SCALE,&MuniView::OnTuneScale,this);
  Unbind(wxEVT_CLOSE_WINDOW,&MuniView::OnCloseTune,this,tune->GetId());
  tune->Destroy();
  tune = 0;
}


void MuniView::OnHeader(wxCommandEvent& event)
{
  MuniHeader *header = new MuniHeader(this,config);
  header->SetHdu(fits.Hdu(hdusel));
  header->Show();
}


void MuniView::OnAstrometry(wxCommandEvent& event)
{
  //  FitsTable table;

  for(size_t i = 0; i < fits.HduCount(); i++) {
    const FitsHdu h(fits.Hdu(i));

    if( h.IsOk() && h.GetExtname().Find(APEREXTNAME) != wxNOT_FOUND ) {
      //      table = FitsTable(fits.Hdu(i));
      GetDisplay()->Astrometry(fits.GetFullPath(),FitsTable(fits.Hdu(i)));
      break;
    }
  }

  //  GetDisplay()->Astrometry(fits.GetFullPath(),table);
}

void MuniView::OnFind(wxCommandEvent& event)
{
  GetDisplay()->OnFind(fits.GetFullPath(),fits.Hdu(hdusel));
}

void MuniView::OnShowLegend(wxCommandEvent& event)
{
  wxASSERT(GetDisplay());
  config->display_legend = event.IsChecked();
  GetDisplay()->ConfigUpdate();
}

void MuniView::OnShowSources(wxCommandEvent& event)
{
  wxASSERT(GetDisplay());
  config->display_sources = event.IsChecked();
  // ConfigUpdate ?
  GetDisplay()->ShowSources(event.IsChecked(),fits);
}

void MuniView::HelpAbout(wxCommandEvent& WXUNUSED(event))
{
  MuniAbout(config->munipack_icon);
}

void MuniView::OnMenuExt(wxCommandEvent& event)
{
  wxLogDebug("MuniView::OnMenuExt ");

  size_t hdu = hdusel;

  if( event.GetId() == wxID_BACKWARD )
    hdu = hdusel - 1;

  else if( event.GetId() == wxID_FORWARD )
    hdu = hdusel + 1;

  else {
    for(size_t i = 0; i < viewid.size(); i++)
      if( viewid[i] == event.GetId() )
	hdu = i;
  }

  wxASSERT(0 <= hdu && hdu < fits.HduCount());
  if( hdu == hdusel )
    return;

  SwitchHdu(hdu);
}


void MuniView::OnExtChanged(wxCommandEvent& event)
{
  wxLogDebug("MuniView::OnExtChanged %d",int(event.GetInt()));
  wxASSERT( !places.empty() );

  int n = event.GetInt();
  SwitchHdu(n);
}

void MuniView::SwitchHdu(int hdu)
{
  //  wxLogDebug("MuniView::SwitchHdu");
  switch_hdu = true;
  wxASSERT(hdu < int(places.size()) && places.size() > 0);
  Clear();

  places[hdusel]->Show(false);
  hdusel = hdu;
  places[hdusel]->Show(true);

  int type = GetHduType(hdusel);
  SetupMenubar(type);
  SetupToolbar(type);
  SetupPanel(type);
  SetupCaption(type);

  if( fits.Hdu(hdusel).IsDisplayImplemented() ) {
    if( GetDisplay()->GetZoom() < 0 ) {
      FitsArray array(fits.Hdu(hdusel));
      wxSize isize = wxSize(array.GetWidth(),array.GetHeight());
      wxLogDebug("MuniView::SwitchHdu zoom: %d %d %d %d",display_size.GetWidth(),
		 display_size.GetHeight(),isize.GetWidth(),isize.GetHeight());
      display_update = true;
      image_size = wxSize(array.GetWidth(),array.GetHeight());
      zoom_init = true;
    }
    GetDisplay()->StartRendering();
  }

  menuExt->Check(viewid[hdusel],true);
  extlist->ChangeSelection(hdusel);
  SetupTitle();
  Layout();
}


void MuniView::ClearMenubar()
{
  // menu
  wxMenuBar *menuBar = GetMenuBar();
  wxASSERT(menuBar);
  const char *names[] = { MENU_IMAGE, MENU_TABLE, MENU_HEAD, /*MENU_TOOLS,*/ 0};
  for(size_t i = 0; names[i] != 0; i++) {
    int n = menuBar->FindMenu(names[i]);
    if( n != wxNOT_FOUND ) {
      wxMenu *menu = menuBar->Remove(n);
      delete menu;
    }
  }
  SetMenuBar(menuBar);
}

void MuniView::ClearToolbar()
{
  wxToolBar *tbar = GetToolBar();
  wxASSERT(tbar);
  tbar->ClearTools();
  tbar->Realize();
}



void MuniView::ClearPanel()
{
  magnifier->Show(false);
  tuner->Show(false);
  zoomview->Show(false);
  zoomer->Show(false);
}

void MuniView::ClearCaption()
{
  caption->Show(false);
}


void MuniView::Clear()
{
  ClearMenubar();
  ClearToolbar();
  ClearPanel();
  ClearCaption();
}


void MuniView::SetupMenubar(int type)
{
  switch(type) {
  case HDU_IMAGE: SetupMenubar_image(); break;
  case HDU_TABLE: SetupMenubar_table(); break;
  case HDU_HEAD:  SetupMenubar_head(); break;
  }
}

void MuniView::SetupMenubar_image()
{
  wxMenu *menuImage = new wxMenu;
  menuImage->Append(ID_INFO,"Header...");
  menuImage->Append(ID_TUNE,"Tune...");
  menuImage->AppendSeparator();
  //  menuImage->Append(wxID_ZOOM_FIT,"Best fit\tCtrl+*");
  menuImage->Append(wxID_ZOOM_100,"Normal size\tCtrl+0");
  menuImage->Append(wxID_ZOOM_IN,"Zoom\tCtrl++");
  menuImage->Append(wxID_ZOOM_OUT,"Shrink\tCtrl+-");
  menuImage->AppendSeparator();
  menuImage->AppendCheckItem(ID_LEGEND,"Show Legend");
  //  menuImage->AppendCheckItem(ID_SOURCES,"Show detected sources\tCtrl+G");
  //  menuImage->AppendSeparator();
  //  menuImage->AppendCheckItem(ID_DETAIL,"Show Detail panel");
  //  menuImage->Check(ID_DETAIL,config->detail_show);
  //  menuImage->Check(ID_SOURCES,config->display_sources);
  menuImage->Check(ID_LEGEND,config->display_legend);

  wxMenu *menuTools = new wxMenu;
  menuTools->Append(ID_FIND,"Find stars...");
  menuTools->Append(ID_APHOT,"Photometry...");
  menuTools->Append(ID_ASTROMETRY,"Astrometry...");

  GetMenuBar()->Insert(3,menuImage,MENU_IMAGE);
  //  GetMenuBar()->Insert(4,menuTools,MENU_TOOLS);

  //  menuView->Enable(ID_FULLSCREEN,true);
}

void MuniView::SetupMenubar_table()
{
  wxMenu *menuTable = new wxMenu;
  menuTable->Append(ID_INFO,"Header...");
  GetMenuBar()->Insert(3,menuTable,MENU_TABLE);
}

void MuniView::SetupMenubar_head()
{
  wxMenu *menuHead = new wxMenu;
  menuHead->Append(ID_INFO,"Header...");
  GetMenuBar()->Insert(3,menuHead,MENU_HEAD);
}


void MuniView::SetupToolbar(int type)
{
  switch(type) {
  case HDU_IMAGE: SetupToolbar_image(); break;
  case HDU_TABLE:
  case HDU_HEAD:  SetupToolbar_table(); break; // just for a moment
  }
}

void MuniView::SetupToolbar_image()
{
  MuniArtIcons ico(wxART_TOOLBAR,wxSize(22,22));
  wxToolBar *tbar = GetToolBar();
  wxASSERT(tbar);

  tbar->AddTool(ID_INFO,"Header",ico.Icon(wxART_INFORMATION),"Show header");
  tbar->AddSeparator();

  // tbar->AddTool(wxID_ZOOM_FIT,"Best",ico.Icon("zoom-fit-best"),"Fit to size");
  // tbar->AddTool(wxID_ZOOM_100,"Normal",ico.Icon("zoom-original"),"Normal size");

  tbar->EnableTool(ID_INFO,true);
  tbar->Realize();
}


void MuniView::SetupToolbar_table()
{
  MuniArtIcons ico(wxART_TOOLBAR,wxSize(22,22));
  wxToolBar *tbar = GetToolBar();
  wxASSERT(tbar);
  tbar->AddTool(ID_INFO,"Header",ico.Icon(wxART_INFORMATION),"Show header");
  tbar->EnableTool(ID_INFO,true);
  tbar->Realize();
}


void MuniView::SetupPanel(int type)
{
  if( type == HDU_IMAGE && GetDisplay()) {
    wxASSERT(GetDisplay() && tuner && zoomview && zoomer);
    FitsTone tone = GetDisplay()->GetTone();
    tuner->SetTone(tone.GetQblack(),tone.GetRsense());

    if( fits.Hdu(hdusel).IsDisplayImplemented() ) {
      FitsArray array(fits.Hdu(hdusel));
      zoomview->SetArray(FitsArray(fits.Hdu(hdusel)));

      wxSize isize(array.GetWidth(),array.GetHeight());
      zoomview->SetImageSize(isize);
      zoomview->SetTune(tone,GetDisplay()->GetItt(),
			GetDisplay()->GetPalette(),GetDisplay()->GetColour());
      zoomview->Show(true);
    }
    zoomer->Show(true);
    tuner->Show(true);
    magnifier->Show(false);
  }
}

void MuniView::SetupCaption(int type)
{
  if( type == HDU_IMAGE) {
    wxASSERT(caption);
    FitsArray array = static_cast<FitsArray>(fits.Hdu(hdusel));
    caption->SetArray(array);
    caption->Show(config->view_controls);
  }
}


void MuniView::OnDraw(MuniDrawEvent& event)
{
  // The repeater, redirects all draw events (by toolboxes, etc) to Display

  if( fits.IsOk() && GetHduType() == HDU_IMAGE && aphot && GetDisplay() )
    wxQueueEvent(GetDisplay(),event.Clone());
}


// void MuniView::Coloring(wxCommandEvent& WXUNUSED(event))
// {
//   if( coloring ) return;

//   Bind(EVT_FILELOAD,&MuniView::OnColoringFinish,this);

//   coloring = new MuniColoring(this,config);
//   coloring->Show(true);
// }

// void MuniView::OnColoringFinish(wxCommandEvent& event)
// {
//   wxLogDebug("MuniView::OnColoringFinish");

//   wxASSERT(coloring);
//   coloring->Destroy();
//   coloring = 0;

//   Unbind(EVT_FILELOAD,&MuniView::OnColoringFinish,this);

//   // CHECK!!! (probably unproper handling of backups !!!!

//   wxString file(event.GetString());
//   if( ! file.IsEmpty() )
//     LoadFile(file);
// }

void MuniView::OnConeSearch(wxCommandEvent& event)
{
  // fill the object entry by the appropriate value by header
  wxString object;
  if( fits.IsOk() ) {
    const FitsHdu head = fits.Hdu(hdusel);
    object = head.GetKey(config->fits_key_object);
  }

  MuniCone *cone = new MuniCone(this,config,object);
  cone->Show();
}

void MuniView::OnCloseCone(wxCommandEvent& event)
{
  wxString cone(event.GetString());
  if( ! cone.IsEmpty() ) {
    MuniView *w = new MuniView(0,config);
    w->Show();
    w->LoadFile(cone);
  }
}

void MuniView::OnAphot(wxCommandEvent& event)
{
  wxASSERT(aphot == 0 && GetDisplay());
  aphot = new MuniAphot(this,config,fits,GetDisplay()->GetImage());
  aphot->Show();
  Bind(EVT_TOOL_FINISH,&MuniView::OnAphotFinish,this);
}

void MuniView::OnAphotFinish(wxCommandEvent& event)
{
  wxLogDebug("MuniView::OnAphotFinish %d",event.GetInt());

  Unbind(EVT_TOOL_FINISH,&MuniView::OnAphotFinish,this);
  aphot->Destroy();
  aphot = 0;
}


void MuniView::OnUpdateBackward(wxUpdateUIEvent& event)
{
  wxASSERT(fits.IsOk());
  event.Enable(fits.HduCount() > 1  && hdusel > 0);
}

void MuniView::OnUpdateForward(wxUpdateUIEvent& event)
{
  wxASSERT(fits.IsOk());
  event.Enable(fits.HduCount() > 1  && hdusel < fits.HduCount() - 1);
}


void MuniView::OnUpdateShowLegend(wxUpdateUIEvent& event)
{
  wxASSERT(0 <= hdusel && hdusel < fits.HduCount());

  bool enable = false;
  if( fits.IsOk() && fits.Hdu(hdusel).IsDisplayImplemented() ) {
    FitsArray a(fits.Hdu(hdusel));
    enable = a.HasWCS();
  }
  event.Enable(enable);
}

void MuniView::OnUpdateShowSources(wxUpdateUIEvent& event)
{
  event.Enable(fits.HasPhotometry() || fits.HasFind());
}

void MuniView::OnUpdateFullScreen(wxUpdateUIEvent& event)
{
  if( fits.IsOk() )
    event.Enable(fits.IsOk() && fits.Hdu(hdusel).IsDisplayImplemented());
}

void MuniView::OnUpdateShowTune(wxUpdateUIEvent& event)
{
  event.Enable(!tune && fits.IsOk() && fits.Hdu(hdusel).IsDisplayImplemented());
}

void MuniView::OnUpdateZoomer(wxUpdateUIEvent& event)
{
  event.Enable(display_size.IsFullySpecified() && !magnifier->IsShown());
}

void MuniView::OnUpdateStop(wxUpdateUIEvent& event)
{
  event.Enable(sloader);
}

void MuniView::OnUpdateZoom100(wxUpdateUIEvent& event)
{
  MuniDisplay *display = GetDisplay();
  if( display && fits.IsOk() && fits.Hdu(hdusel).IsDisplayImplemented() )
    event.Enable(fabs(display->GetZoom() - 1.0) > 1e-3);
  else
    event.Enable(false);
}

void MuniView::OnUpdateZoomIn(wxUpdateUIEvent& event)
{
  MuniDisplay *display = GetDisplay();
  if( fits.IsOk() && fits.Hdu(hdusel).IsDisplayImplemented() && display ) {
    double zoom = display->GetZoom();

    MuniZoomSet set(zoom);
    event.Enable(fabs(set.Increase() - zoom) > 0.001);
  }
  else
    event.Enable(false);
}

void MuniView::OnUpdateZoomOut(wxUpdateUIEvent& event)
{
  MuniDisplay *display = GetDisplay();
  if( fits.IsOk() && fits.Hdu(hdusel).IsDisplayImplemented() && display ) {
    double zoom = display->GetZoom();
    MuniZoomSet set(zoom);
    event.Enable(fabs(set.Decrease() - zoom) > 0.001);
  }
  else
    event.Enable(false);
}

void MuniView::OnUpdateExportAs(wxUpdateUIEvent& event)
{
  event.Enable(fits.IsOk() && fits.Hdu(hdusel).IsDisplayImplemented());
}

int MuniView::GetHduType() const
{
  wxASSERT(fits.IsOk() && 0 <= hdusel && hdusel < fits.HduCount());
  return fits.Hdu(hdusel).Type();
}


int MuniView::GetHduType(int hdu) const
{
  wxASSERT(fits.IsOk() && 0 <= hdu && hdu < int(fits.HduCount()));
  return fits.Hdu(hdu).Type();
}


MuniDisplay *MuniView::GetDisplay() const
{
  if(0 <= hdusel && hdusel < places.size() && places.size() > 0) {
    MuniDisplay *display = dynamic_cast<MuniDisplay *>(places[hdusel]);
    return display;
  }
  else
    return NULL;
}

void MuniView::SetupTitle()
{
  if( fits.IsOk() && 0 <= hdusel && hdusel < fits.HduCount()) {

    wxString extname = fits.Hdu(hdusel).GetControlLabel().Capitalize();
    wxString title = fits.GetName() + (IsModified() ? "*" : "") +
      (fits.HduCount() > 1 ? ": " + extname : "");
    SetTitle(title);
  }
}


void MuniView::SetupToolbar_load()
{
  MuniArtIcons ico(wxART_TOOLBAR,wxSize(22,22));

  wxToolBar *tbar = GetToolBar();
  tbar->AddTool(wxID_STOP,"Stop",ico.Icon("process-stop"),"Stop loading");
  tbar->Realize();
}

void MuniView::SetupToolbar_init()
{
  MuniArtIcons ico(wxART_TOOLBAR,wxSize(22,22));

  wxToolBar *tbar = GetToolBar();
  tbar->AddTool(wxID_OPEN,"Open",ico.Icon(wxART_FILE_OPEN),"Open FITS");
  tbar->Realize();
}


void MuniView::OnMouseMotion(MuniSlewEvent& event)
{
  //wxLogDebug("MuniView::OnMouseMotion %d %d",event.x,event.y);
  /*
  if( event.entering ) wxLogDebug("MuniView::OnMouseMotion --- entering");
  if( event.leaving ) wxLogDebug("MuniView::OnMouseMotion --- leaving");
  */

  if( config->view_controls ) {

    if( event.entering ) {
      zoomview->Show(false);
      magnifier->Show(true);
      panelsizer->Replace(zoomview,magnifier);
    }
    else if( event.leaving ) {
      zoomview->Show(true);
      magnifier->Show(false);
      panelsizer->Replace(magnifier,zoomview);
    }
    Layout();

    wxASSERT(magnifier && caption);
    wxQueueEvent(magnifier,event.Clone());
    wxQueueEvent(caption,event.Clone());
  }

}

void MuniView::OnRender(MuniRenderEvent& event)
{
  //  wxLogDebug("MuniView::OnRender() ..... events due MuniDisplay() ");
  if( event.GetId() == ID_RENDER_FINISH &&
      fits.Hdu(hdusel).IsDisplayImplemented() ) {
    wxASSERT(zoomview && GetDisplay());
    zoomview->StartRendering();
    magnifier->SetImage(GetDisplay()->GetImage());
  }
}


void MuniView::OnTuner(MuniTuneEvent& event)
{
  wxASSERT(zoomview && GetDisplay());
  wxLogDebug("MuniView::OnTuner %d %d %d",int(event.GetId()),int(ID_TONE_BLACK),
	     int(tuner->GetId()));

  FitsTone tone = GetDisplay()->GetTone();
  switch(event.GetId()) {
  case ID_TONE_QBLACK:   tone.SetQblack(event.x); break;
  case ID_TONE_RSENSE:   tone.SetRsense(event.x);break;
  }
  GetDisplay()->SetTone(tone);

  zoomview->SetTone(tone);
  if( tune )
    tune->SetScale(tone.GetBlack(),tone.GetSense());
}


void MuniView::OnTuneScale(MuniTuneScaleEvent& event)
{
  wxASSERT(tuner && tune && GetDisplay());
  wxLogDebug("MuniView::OnTuneScale");

  FitsTone tone = GetDisplay()->GetTone();
  if( event.reset ) {
    tone.Reset();
    tune->SetScale(tone.GetBlack(),tone.GetSense());
  }
  else {
    tone.SetBlack(event.black);
    tone.SetSense(event.sense);
  }
  GetDisplay()->SetTone(tone);

  zoomview->SetTone(tone);
  tuner->SetTone(tone.GetQblack(),tone.GetRsense());
}

void MuniView::OnTuneItt(MuniTuneIttEvent& event)
{
  wxASSERT(zoomview && GetDisplay());
  wxLogDebug("MuniView::OnTuneItt");

  FitsItt itt = GetDisplay()->GetItt();
  itt.SetItt(event.itt);
  GetDisplay()->SetItt(itt);

  zoomview->SetItt(itt);
}

void MuniView::OnTunePal(MuniTunePalEvent& event)
{
  wxASSERT(zoomview && GetDisplay());
  wxLogDebug("MuniView::OnTunePal");

  FitsPalette pal = GetDisplay()->GetPalette();
  pal.SetInverse(event.inverse);
  pal.SetPalette(event.palette);
  GetDisplay()->SetPalette(pal);
  zoomview->SetPalette(pal);
}

void MuniView::OnTuneColour(MuniTuneColourEvent& event)
{
  wxASSERT(zoomview && GetDisplay());
  wxLogDebug("MuniView::OnTuneColour");

  FitsColour c = GetDisplay()->GetColour();
  c.SetSaturation(event.saturation);
  c.SetWhitePoint(event.x_white,event.y_white);
  GetDisplay()->SetColour(c);
  zoomview->SetColour(c);
}

void MuniView::OnTuneNite(MuniTuneNiteEvent& event)
{
  wxASSERT(zoomview && GetDisplay());
  wxLogDebug("MuniView::OnTuneNite");

  FitsColour c = GetDisplay()->GetColour();
  c.SetNiteVision(event.nite);
  c.SetMesoLevel(event.level);
  c.SetMesoWidth(event.width);
  GetDisplay()->SetColour(c);
  zoomview->SetColour(c);
}
