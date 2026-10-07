/*

  XMunipack


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

#ifndef _XMUNIPACK_VIEW_H_
#define _XMUNIPACK_VIEW_H_

#include "../config.h"
#include "config.h"
#include "fits.h"
#include "enum.h"
#include "types.h"
#include "event.h"
#include "mconfig.h"
#include "vocatconf.h"
#include "aphot.h"
#include "enum.h"
#include "fits.h"
#include "icon.h"
#include "render.h"
#include "dispreview.h"
#include "caption.h"
#include "zoomer.h"
#include "tune.h"
#include "tuner.h"
#include "dismag.h"
#include "preferences.h"
#include <wx/wx.h>
#include <wx/listctrl.h>
#include <wx/dataview.h>
#include <wx/thread.h>
#include <wx/timer.h>
#include <wx/filepicker.h>
#include <wx/spinctrl.h>
#include <wx/grid.h>
#include <wx/artprov.h>
#include <wx/graphics.h>
#include <queue>
#include <vector>
#include <list>

#if wxUSE_FSWATCHER
#include <wx/fswatcher.h>
#endif

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif



class FitsStream: public wxThread
{
  wxEvtHandler *handler;
  wxString filename;
  FitsTone tone;
  int shrink,icon_size;
  bool stop, stop_wait;

  ExitCode Entry();
  wxImage RenderPreview(int,int,const float*);
  void make_icons(const FitsFile&, const std::vector<FitsTone>&,
		  std::vector<wxImage>&);

public:
  FitsStream(wxEvtHandler *,const wxString&, int);
  void Stop();
  void StopAndWait();
  void SetZoom(float);

};





/*
class MuniLog: public wxLogWindow
{
public:
  MuniLog(wxWindow *);
  virtual ~MuniLog();

  virtual bool OnFrameClose(wxFrame *);
  bool Visibility();

private:

  wxWindow *topwin;
  bool visible;

};
*/




// class MuniDialog: public wxDialog
// {
//   MuniConfig *config;
//   MuniProcess *mproc;
//   wxTimer timer;
//   wxAnimationCtrl *anim;
//   wxButton *execute, *stop;
//   wxStaticText *status;
//   wxString action, help;
//   long index;
//   wxArrayString out;

//   void OnHelp(wxCommandEvent&);
//   void OnExecute(wxCommandEvent&);
//   void OnStop(wxCommandEvent&);
//   void OnTimer(wxTimerEvent&);
//   void OnFinish(wxProcessEvent&);
//   void OnStdButton(wxCommandEvent&);

//  protected:

//   wxString fitsname, tmpfits;

//   void SetAction(const wxString&);
//   void SetHelp(const wxString&);
//   void SetStatus(const wxString&);
//   wxString Parser(const wxString&, const wxString&) const;
//   //  void Write(const wxString&, ...);
//   wxSizer *CreateStatusSizer(const wxString&);
//   wxSizer *CreateButtonSizer(long);
//   wxArrayString GetOutput() const;
//   wxArrayString GetLastOutput();
//   virtual void OnInput(MuniProcess *) = 0;
//   virtual void OnOutput(const wxArrayString&);
//   virtual void CleanDraw() const {}
//   //  virtual void DrawStars(const wxArrayString&) const;

//   // public:
//   MuniDialog(wxWindow *, MuniConfig *, const wxString&, const wxString&);
//   virtual ~MuniDialog();

// };


/*
class MuniPhotometry: public MuniBaseDialog
{
public:
  MuniPhotometry(wxWindow *,MuniConfig *, const wxString&);
  virtual ~MuniPhotometry();

  wxString GetBackup() const;

private:

  wxGauge *progress;
  wxButton *calbutt;
  double fwhm,thresh,saturation,readnoise,phpadu;
  wxString backupfile;
  bool erase,apply,showtooltip;
  std::list<int> ids;

  bool ParseProcessing(const wxArrayString&);
  void EraseCanvas();
  void OnTimer(wxTimerEvent&);

  void OnUpdateUI(wxUpdateUIEvent&);
  void OnIdle(wxIdleEvent&);
  void OnFwhm(wxSpinDoubleEvent&);
  void OnThresh(wxSpinDoubleEvent&);
  void OnSaturation(wxSpinDoubleEvent&);
  void OnReadNoise(wxSpinDoubleEvent&);
  void OnPhpADU(wxSpinDoubleEvent&);
  void OnApply(wxCommandEvent&);
  void OnCancel(wxCommandEvent&);
  void OnFinish(wxProcessEvent&);
  void OnPhotometry(wxCommandEvent&);

};
*/

class MuniFind: public wxDialog
{
  MuniConfig *config;
  MuniProcess *mproc;
  wxTimer timer;
  wxAnimationCtrl *anim;
  wxWindow *plot;
  wxString fitsname, tmpfits;
  const FitsArray array;
  double fwhm, thresh, satur;
  long index;
  wxButton *start, *stop;
  wxStaticText *status;

  double InitSatur() const;
  wxString Parser(const wxString&) const;
  bool StarParser(const wxString&, long *, double *, double *) const;
  void DrawStars(long, long, const wxArrayString&) const;
  void CleanDraw() const;
  long LastStar(const wxString&) const;

  void OnHelp(wxCommandEvent&);
  void OnFind(wxCommandEvent&);
  void OnStop(wxCommandEvent&);
  void OnTimer(wxTimerEvent&);
  void OnFindFinish(wxProcessEvent&);

  void OnFwhm(wxSpinDoubleEvent&);
  void OnThresh(wxSpinDoubleEvent&);
  void OnSatur(wxSpinDoubleEvent&);
  void OnStdButton(wxCommandEvent&);

 public:
  MuniFind(wxWindow *, MuniConfig *, const wxString&, const FitsHdu&);
  virtual ~MuniFind();

  void SetPoint(int,int) const;

};


class MuniCone: public wxDialog
{
public:
  MuniCone(wxWindow *,MuniConfig *,const wxString& =wxEmptyString);
  virtual ~MuniCone();
  wxString GetPath() const { return tmpfile; }

private:

  MuniConfig *config;
  MuniProcess *mproc;
  wxTimer timer;
  VOCatConf catconf;
  wxAnimationCtrl *throbber;
  wxStaticText *status, *johnson_label;
  wxTextCtrl *objentry, *alpha, *delta;
  wxButton *search, *stop;
  wxCheckBox *johnson;
  bool apply;

  wxSpinCtrlDouble *radius, *magmin, *magmax;
  wxString tmpfile, catfits, object;
  double ra,dec;
  long index;
  std::list<int> ids;

  void OnUpdateUI(wxUpdateUIEvent&);
  void OnOk(wxCommandEvent&);
  void OnClose(wxCloseEvent&);
  void OnHelp(wxCommandEvent&);
  void OnObjectName(wxCommandEvent&);
  void OnObjectEnter(wxCommandEvent&);
  void OnGetCoo(wxCommandEvent&);
  void OnRightAscension(wxCommandEvent&);
  void OnDeclination(wxCommandEvent&);
  void OnSearch(wxCommandEvent&);
  void OnFinish(wxProcessEvent&);
  void OnStop(wxCommandEvent&);
  void OnTimer(wxTimerEvent&);
  void OnService(wxCommandEvent&);
  wxString ParseOutput(const wxArrayString&);
  void Resolve();
  void ResolveFinish(wxProcessEvent&);
  double deg(const wxString&, double =1.0);

};


class MuniSelectSource: public wxDialog
{
public:
  MuniSelectSource(wxWindow *, MuniConfig *, bool =false);
  virtual ~MuniSelectSource();
  wxString GetPath() const;
  wxString GetId() const;
  int GetType() const;
  bool GetRelative() const;
  bool IsTemporary() const;
  wxString GetLabelRA() const;
  wxString GetLabelDec() const;
  wxString GetLabelPMRA() const;
  wxString GetLabelPMDec() const;
  wxString GetLabelMag() const;

private:

  MuniConfig *config;
  int page;
  wxChoice *choice_ra, *choice_dec, *choice_pmra, *choice_pmdec, *choice_mag;
  wxString reffile,catfile,tmpcatfile,idlabel,label_ra,label_dec,
    label_pmra,label_pmdec,label_mag;
  bool xframe, astrorel;

  void OnUpdateUI(wxUpdateUIEvent&);
  void CreateControls();
  void OnCheckRel(wxCommandEvent&);
  void OnRefFile(wxFileDirPickerEvent&);
  void OnCatFile(wxFileDirPickerEvent&);
  void OnCatVO(wxCommandEvent&);
  bool CheckCatalogue(const wxString&, wxArrayString&);
  void OnChoice(wxCommandEvent&);
  void OnBookChange(wxBookCtrlEvent&);
  void EraseTemp();
  void SetLabels(const wxArrayString&);

};


class MuniAstrometryOptions: public wxPanel
{
public:
  MuniAstrometryOptions(wxWindow *,MuniConfig *);
  ~MuniAstrometryOptions();

  int GetMatchType() const;
  int GetMinMatch() const;
  int GetMaxMatch() const;
  double GetSig() const;
  double GetFSig() const;
  bool GetFullMatch() const;
  wxString GetOutputUnits() const;

private:

  MuniConfig *config;

  double sig, fsig;
  int minmatch, maxmatch, matchtype;
  bool full_match;
  wxString output_units;
  std::list<int> ids;

  void OnUpdateUI(wxUpdateUIEvent&);
  void OnSpinSig(wxSpinDoubleEvent&);
  void OnSpinFSig(wxSpinDoubleEvent&);
  void OnSpinMinMatch(wxSpinEvent&);
  void OnSpinMaxMatch(wxSpinEvent&);
  void OnMatchType(wxCommandEvent&);
  void OnFullMatch(wxCommandEvent&);
  void OnChoiceUnits(wxCommandEvent&);

};

class MuniAstrometry: public wxDialog
{
public:
  MuniAstrometry(wxWindow *, MuniConfig *);
  virtual ~MuniAstrometry();
  void SetFile(const wxString&, const FitsArray&);
  void SetDetectedSources(const FitsTable&);
  wxString GetBackup() const;

private:

  MuniConfig *config;
  wxSpinCtrlDouble *wscale, *wangle, *acenter, *dcenter;
  wxButton *sunit, *calbutt, *stopbutt, *savebutt, *rembutt;
  wxGauge *progress;
  wxStaticText *overlayid, *refcatid, *info, *infolabel, *proginfo, *proglabel;
  wxBoxSizer *autosizer;
  wxCheckBox *overlay_check, *reflex_checkbox;
  FitsArray array;
  FitsTable catalogue, stars, coverlay;
  wxTimer timer;
  wxDateTime start;
  MuniPipe pipe;
  MuniAstrometryOptions *astropt;

  int nstars,nhist;
  double maglim, s0, rms;
  double xoff,yoff,alpha, delta, scale, reflex, angle, amin, dmin, fmin;
  wxString proj,label_ra,label_dec,label_mag, label_pmra, label_pmdec;
  bool init_par, init_file, init_ref, running, parsing, unsaved,
    draw_overlay, edited, relative, tmpcat,readonly,showtooltip;
  wxString file,catfile,reffile,backupfile,workingfile;
  int output_index, page;
  std::vector<int> hist;
  std::vector<double> x,y,u,v;
  std::list<int> ids;

  void OnUpdateUI(wxUpdateUIEvent&);
  void OnIdle(wxIdleEvent&);
  void OnStdButton(wxCommandEvent&);
  void OnCalibrate(wxCommandEvent&);
  void OnCalibrateStop(wxCommandEvent&);
  void OnCalibrateFinish(wxProcessEvent&);
  void OnOverlay(wxCommandEvent&);
  void OnReflex(wxCommandEvent&);
  void OnDrawOverlay(wxCommandEvent&);
  void OnReference(wxCommandEvent&);
  void OnTimer(wxTimerEvent&);
  void CreateControls();
  void RunProcessing();
  void DrawOverlay(const FitsTable&);
  void OnPopScaleUnit(wxCommandEvent&);
  void OnScaleUnit(wxCommandEvent&);
  void OnInitRef(wxCommandEvent&);
  void OnSpinDouble(wxSpinDoubleEvent&);
  void OnTextDouble(wxCommandEvent&);
  FitsTable LoadCatalogue(const wxString&);
  void InitByCatalogue(const FitsTable&);
  void ParseProcessing(const wxArrayString&);
  void Reset();
  void OnSaveWCS(wxCommandEvent&);
  void OnSaveWCSFinish(wxProcessEvent&);
  void OnChoiceProj(wxCommandEvent&);
  void RemoveWorkingFile();
  void EraseTemp();
  bool FinishClean(wxCommandEvent&);
  double Scale(double) const;
  double Period(double) const;
  void OnRemove(wxCommandEvent&);
  void OnRemoveFinish(wxProcessEvent&);
  void OnBookChange(wxBookCtrlEvent&);

};

class MuniAstrolog
{
  wxDateTime datetime;

 public:
  MuniAstrolog();
  MuniAstrolog(const wxDateTime&);

  wxString GetSign() const;
};

class MuniAstrometer: public wxDialog
{
public:
  MuniAstrometer(wxWindow *, MuniConfig *, const std::vector<FitsMeta> &);
  virtual ~MuniAstrometer();

private:

  MuniConfig *config;
  const std::vector<FitsMeta> list;
  wxDataViewListCtrl *mtable;
  wxFilePickerCtrl *fpick;
  wxStaticText *refcatid;
  wxBoxSizer *topsizer;
  wxGauge *gstat;
  wxButton *butt, *catbutt;
  wxTimer timer;
  MuniPipe pipe;
  MuniAstrometryOptions *astropt;

  int findex,lastrow;
  wxString proj, reffile,catfile,label_ra,label_dec,label_pmra,label_pmdec,
    label_mag;
  bool running,astrorel,relative,tmpcat;
  std::list<int> ids;

  void OnClose(wxCloseEvent&);
  void OnUpdateUI(wxUpdateUIEvent&);
  void CreateControls();
  void CreateProcess();
  void SetTable();
  void ParseOutput();
  void OnProcess(wxCommandEvent&);
  void OnTimer(wxTimerEvent&);
  void OnFinish(wxProcessEvent&);
  void OnReference(wxCommandEvent&);
  void OnChoiceProj(wxCommandEvent&);
  void EraseTemp();

};


class MuniCalibrate: public wxDialog
{
public:
  MuniCalibrate(wxWindow *,MuniConfig *, const wxString&);
  wxString GetResult() const;

private:

  MuniConfig *config;
  wxAnimationCtrl *throbber;
  double fwhm,thresh, alpha, delta, radius;
  wxString file, output, catalog, projection, coutput;
  MuniPipe pipe;

  void OnFwhm(wxSpinDoubleEvent&);
  void OnThresh(wxSpinDoubleEvent&);
  void OnCatalog(wxCommandEvent&);
  void OnProjection(wxCommandEvent&);
  void OnAlpha(wxCommandEvent&);
  void OnDelta(wxCommandEvent&);
  void OnRadius(wxCommandEvent&);
  void OnApply(wxCommandEvent&);
  void OnCancel(wxCommandEvent&);
  void OnFinish(wxProcessEvent&);

};

class MuniColoring: public wxDialog
{
public:
  MuniColoring(wxWindow *,MuniConfig *);
  bool SetDropMeta(int, int, const std::vector<FitsMeta>&);

private:

  const int iSize;

  MuniConfig *config;
  wxChoice *cspace;
  wxString filename,dirname,colorspace;
  wxArrayString cchoices, opt;
  std::vector<FitsMeta> metalist;
  std::vector<wxString> param_lines;
  wxListView *list;
  wxImageList *icons;
  int index;
  MuniPipe pipe;
  wxAnimationCtrl *throbber;

  void Init();
  void CreateControls();
  void InitList(const wxString&);
  void SetMeta(int, const FitsMeta&, double =0.0, double =-1.0);

  void OnBandfile(wxFileDirPickerEvent&);
  void OnFilename(wxCommandEvent&);
  void OnDirname(wxFileDirPickerEvent&);
  void OnApply(wxCommandEvent&);
  void OnCancel(wxCommandEvent&);
  void OnFinish(wxProcessEvent&);
  void OnColorspace(wxCommandEvent&);
  void OnListSelected(wxListEvent&);
  void OnUpdateBandfile(wxUpdateUIEvent&);
  void OnUpdateOk(wxUpdateUIEvent&);

};

class MuniDisplay: public wxPanel
{
  MuniConfig *config;
  MuniDisplayRender *render;
  FitsArray array;
  wxBitmap canvas;
  unsigned char *idata;
  int iwidth, iheight, rid;
  double xcen, ycen, zoom;
  bool render_start;
  wxSize display_size;
  wxRect exposed;
  wxPoint offset;
  wxColour bgcolour;
  MuniFind *find;
  MuniAstrometry *astrometry;

  FitsTone tone;
  FitsItt itt;
  FitsPalette pal;
  FitsColour colour;
  wxImage coogrid;

  void OnCloseAstrometry(wxCommandEvent&);
  void OnCloseFind(wxCommandEvent&);
  void OnLeaveFullscreen(wxCommandEvent&);
  void OnClick(MuniClickEvent&);
  void OnMouseMotion(wxMouseEvent&);
  void OnMouseEnter(wxMouseEvent&);
  void OnMouseLeave(wxMouseEvent&);
  void OnKeyDown(wxKeyEvent&);
  void OnZoomer(MuniZoomEvent&);

  void OnSize(wxSizeEvent&);
  void OnIdle(wxIdleEvent&);
  void OnRender(MuniRenderEvent&);
  void OnPaint(wxPaintEvent&);
  wxRect GetArrayRect() const;
  wxRect GetExposed() const;
  wxPoint GetOffset(const wxRect&) const;
  void DrawCanvas();
  void Render();

public:
  MuniDisplay(wxWindow *, MuniConfig *);
  virtual ~MuniDisplay();

  void SetPreview(const wxImage&);
  void SetArray(const FitsArray&, const FitsTone&);
  void StartRendering();

  void SetZoom(double);
  void SetTone(const FitsTone&);
  void SetItt(const FitsItt&);
  void SetColour(const FitsColour&);
  void SetPalette(const FitsPalette&);
  //  virtual void SetOverlay(wxInputStream&);
  virtual void SetStars(FitsTable&);
  double GetZoom() const;
  FitsTone GetTone() const;
  FitsItt GetItt() const;
  FitsPalette GetPalette() const;
  FitsColour GetColour() const;
  wxImage GetImage() const;
  void ShowSources(bool, const FitsFile&);
  void OnFind(const wxString&, const FitsHdu&);
  void Photometry(const wxString&);
  void Astrometry(const wxString&,const FitsTable&);
  void Calibrate(const wxString&);
  void ConfigUpdate();

};



class MuniGrid: public wxGrid
{
public:
  MuniGrid(wxWindow *, MuniConfig *);

  bool SetHdu(const FitsHdu&);
  void InvokeRendering() {}

private:

  FitsTable table;
  int rows_filled;

  void OnIdle(wxIdleEvent&);
};


class MuniHead: public wxTextCtrl
{
  const MuniConfig *config;
  FitsHdu head;

public:
  MuniHead(wxWindow *, const MuniConfig *);
  bool SetHdu(const FitsHdu&);

};

class MuniHeader: public wxFrame
{
  MuniConfig *config;
  MuniHead *header;
  FitsHdu head;
  wxMenu *menuFile/*, *menuEdit*/;

  void OnClose(wxCloseEvent&);
  void FileClose(wxCommandEvent&);
  void FileExport(wxCommandEvent&);
  void EditCut(wxCommandEvent&);
  void EditCopy(wxCommandEvent&);
  void EditPaste(wxCommandEvent&);

public:
  MuniHeader(wxWindow *, MuniConfig *);
  void SetHdu(const FitsHdu&);

};

class MuniSplashing: public wxWindow
{
private:

  wxBitmap logo;
  wxAnimationCtrl *anim;

public:
  MuniSplashing(wxWindow *, const MuniConfig *);
  void Play();
  void Stop();

};



class MuniDummy: public wxPanel
{
  const MuniConfig *config;

public:
  MuniDummy(wxWindow *, const MuniConfig *);
};


class MuniExtensionList: public wxDataViewListCtrl
{
  const MuniConfig *config;

  std::map<int,int> idxtype;
  std::vector<wxIcon> icons, phases;
  wxStopWatch timer;

  std::vector<wxIcon> MoonPhases(int) const;
  void OnExtChanged(wxDataViewEvent&);
  wxIcon DrawIcon(int, int) const;

public:
  MuniExtensionList(wxWindow*, const MuniConfig *);
  void Init();
  void Finish();
  void UpdateProgress(int,double);
  void FinishProgress(int);
  void Clear();
  bool IsOk() const;
  void ChangeSelection(int);
  void Append(const wxString&, int);

};



class MuniView: public wxFrame
{
  MuniConfig *config;
  bool loadfile;  // FSWatch trigger to launch LoadFile in idle time
  FitsFile fits;
  size_t hdusel;
  bool zoom_init, switch_hdu, display_update;
  std::vector<long> viewid;
  wxArrayString backup;
  FitsStream *sloader;

  wxMenu *menuFile,*menuView,*menuExt;
  wxSizer *panelsizer, *placesizer;
  MuniExtensionList *extlist;
  MuniMagnifier *magnifier;
  MuniZoomSlider *zoomer;
  MuniZoomWindow *zoomview;
  MuniTune *tune;
  MuniTuner *tuner;
  MuniPreferences *preferences;
  MuniDisplayCaption *caption;
  MuniDisplayPreview *preview;
  MuniSplashing *splash;
  std::vector<wxWindow *> places;
  std::vector<wxImage> icons;
//  MuniColoring *coloring;
  //  wxLogWindow *console;
  MuniAphot *aphot;
  wxSize display_size, preview_size, image_size;

#if wxUSE_FSWATCHER
  wxFileSystemWatcher *fswatch;
  void OnFileSystemEvent(wxFileSystemWatcherEvent&);
#endif

  void SetFits(const FitsFile&, const std::vector<FitsTone>&);
  void SwitchHdu(int);
  void ReplaceExtensionMenu(const FitsFile&);
  void DeletePlaces();
  void ClearToolbar();
  void ClearMenubar();
  void ClearPanel();
  void ClearCaption();
  void Clear();
  void SetupMenubar(int);
  void SetupMenubar_image();
  void SetupMenubar_table();
  void SetupMenubar_head();
  void SetupToolbar(int);
  void SetupToolbar_image();
  void SetupToolbar_table();
  void SetupToolbar_head();
  void SetupToolbar_init();
  void SetupToolbar_load();
  void SetupPanel(int);
  void SetupCaption(int);
  void SetupZoom();

  int GetHduType() const;
  int GetHduType(int) const;
  double GetBestFitZoom(const wxSize&, const wxSize&) const;
  void SetupTitle();
  MuniDisplay *GetDisplay() const;

  void OnClose(wxCloseEvent& event);
  void OnIdle(wxIdleEvent& event);
  void OnDisplaySize(MuniSizeChangedEvent&);
  void OnClick(MuniClickEvent&);
  void FileClose(wxCommandEvent& WXUNUSED(event));
  void FileOpen(wxCommandEvent& WXUNUSED(event));
  void OnStop(wxCommandEvent& WXUNUSED(event));
  void NewView(wxCommandEvent& WXUNUSED(event));
  void OnConeSearch(wxCommandEvent&);
  void OnCloseCone(wxCommandEvent&);
  void OnFullScreen(wxCommandEvent& WXUNUSED(event));
  void OnTune(wxCommandEvent& WXUNUSED(event));
  void OnCloseTune(wxCloseEvent&);
  void OnPreferences(wxCommandEvent& WXUNUSED(event));
  void OnClosePreferences(wxCloseEvent&);

  void OnFind(wxCommandEvent& WXUNUSED(event));
  void OnAphot(wxCommandEvent& WXUNUSED(event));
  void OnAphotFinish(wxCommandEvent&);
  void OnAstrometry(wxCommandEvent& WXUNUSED(event));
  void OnShowLegend(wxCommandEvent&);
  void OnShowSources(wxCommandEvent&);
  void OnHeader(wxCommandEvent& WXUNUSED(event));
  void HelpAbout(wxCommandEvent& WXUNUSED(event));
  //  void ViewLog(wxCommandEvent&);
  void OnFileProperties(wxCommandEvent& WXUNUSED(event));
  void OnShowToolbar(wxCommandEvent&);
  void OnShowControls(wxCommandEvent&);
  void OnMenuExt(wxCommandEvent&);
  void OnMenuZoom(wxCommandEvent&);
  void OnZoomer(MuniZoomEvent&);
  void OnZoomview(MuniZoomEvent&);
  void FileSave(wxCommandEvent& WXUNUSED(event));
  void FileExport(wxCommandEvent& WXUNUSED(event));
  void OnExtChanged(wxCommandEvent&);
  void OnUpdateStop(wxUpdateUIEvent&);
  void OnUpdateBackward(wxUpdateUIEvent&);
  void OnUpdateForward(wxUpdateUIEvent&);
  void OnUpdateShowLegend(wxUpdateUIEvent&);
  void OnUpdateShowSources(wxUpdateUIEvent&);
  void OnUpdateShowTune(wxUpdateUIEvent&);
  void OnUpdateZoomer(wxUpdateUIEvent&);
  void OnUpdateZoom100(wxUpdateUIEvent&);
  void OnUpdateZoomIn(wxUpdateUIEvent&);
  void OnUpdateZoomOut(wxUpdateUIEvent&);
  void OnUpdateExportAs(wxUpdateUIEvent&);
  void OnUpdateFullScreen(wxUpdateUIEvent&);
  /*
  void Coloring(wxCommandEvent& WXUNUSED(event));
  void OnColoringFinish(wxCommandEvent&);
  */

  void OnRender(MuniRenderEvent&);
  void OnMouseMotion(MuniSlewEvent&);
  void OnDraw(MuniDrawEvent&);
  void OnTuneScale(MuniTuneScaleEvent&);
  void OnTuneItt(MuniTuneIttEvent&);
  void OnTunePal(MuniTunePalEvent&);
  void OnTuneColour(MuniTuneColourEvent&);
  void OnTuneNite(MuniTuneNiteEvent&);
  void OnTuner(MuniTuneEvent&);
  void OnPreferencesEvent(MuniConfigEvent&);

  void OnLoadFits(FitsStreamEvent&);

  int Unsaved(const wxString&);
  void RemoveBackup();

  // File monitor
  void MonitorFile(const wxString&);
  double MonitorDelay(const wxString&) const;

public:
  MuniView(wxWindow *,MuniConfig *);
  virtual ~MuniView();
  void LoadFile(const wxString&);
  void SaveFile(const wxString&);
  void LoadFileBackup(const wxString&, const wxString&);
  bool IsModified() const;
  void CreateFSWatch();

};

#endif
