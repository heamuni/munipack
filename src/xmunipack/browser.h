/*

  XMunipack

  Copyright © 2021-26 F.Hroch (hroch@physics.muni.cz)

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

  There remains some obsolete routines for re-implementation:

  * search
  * average
  * darkbat
  * clipboard-oriented
  * thumbnails ?

*/

#include "../config.h"
#include "mconfig.h"
#include "fits.h"
#include "event.h"
#include "mprocess.h"
#include "list.h"
#include "navigation.h"
#include "preferences.h"
#include "fileprop.h"
#include <wx/wx.h>
#include <wx/filepicker.h>
#include <vector>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif



class MuniAnim
{
  // ve stylu totemu nebo vlc pod mac: dole ovladaci prvky,
  // info o obrazku (text) do statusbaru

};

/*
class MuniThumbnail
{
public:
  MuniThumbnail(const FitsMeta&);
  MuniThumbnail(const wxString&);
  MuniThumbnail(wxInputStream&);

  void Load(const wxString&);
  void Load(wxInputStream&);
  void Save(const wxString&);
  void Save(wxOutputStream&);
  bool IsOk() const;


  FitsMeta GetMeta() const;
  void SetMeta(const FitsMeta&);

  wxString GetURL() const;

private:

  FitsMeta meta;

  wxString icon;
  wxArrayString icons;

  wxXmlDocument CreateXML(const FitsMeta&,const wxString&);
  FitsMeta ParseXML(const wxXmlDocument&);
  wxString CreateIconame(const wxString&,const wxString& =wxEmptyString);

};
*/



/*
class MuniDataObjectMeta: public wxDataObjectSimple
{
public:
  MuniDataObjectMeta();
  MuniDataObjectMeta(const std::vector<FitsMeta>& mlist);
  virtual ~MuniDataObjectMeta();

  size_t GetDataSize() const;
  bool GetDataHere(void *) const;
  bool SetData(size_t, const void *);
  std::vector<FitsMeta> GetMetafitses() const;

private:

  size_t len;
  char *data;
};



class MuniAverage: public MuniListWindow
{
public:
  MuniAverage(wxWindow *,wxWindowID, long =0, MuniConfig * = 0);

private:

  MuniConfig *config;
  MuniPipe pipe;
  wxTimer timer;

  wxTextCtrl *flabel;
  wxDirPickerCtrl *dirpic;
  wxRadioButton *btype[3];
  wxGauge *gauge;
  wxStaticText *label;

  wxString dirname,bitpix,level,filename[3];
  bool robust;
  int xtype;

  void OnDirname(wxFileDirPickerEvent&);
  void OnOptions(wxCommandEvent&);
  void OnCreate(wxCommandEvent&);
  void OnFlabel(wxCommandEvent&);
  void OnBtype(wxCommandEvent&);
  void OnUpdateButt(wxUpdateUIEvent&);
  void OnFinish(wxProcessEvent&);
  void OnUpdate(wxTimerEvent&);

};

class MuniDarkbat: public MuniListWindow
{
public:
  MuniDarkbat(wxWindow *,wxWindowID, long =0, MuniConfig * = 0);

private:

  MuniConfig *config;
  wxFilePickerCtrl *fpic,*dpic,*bpic;
  wxDirPickerCtrl *dirpic;
  wxRadioButton *r0, *r1, *r2;
  wxButton *bcre;
  wxGauge *gauge;
  wxStaticText *label;
  wxString ffilename,dfilename,bfilename,dirname,bitpix,suffix;
  wxArrayString results;
  int mode;
  MuniPipe pipe;
  wxTimer timer;

  void OnFlatname(wxFileDirPickerEvent&);
  void OnDarkname(wxFileDirPickerEvent&);
  void OnBiasname(wxFileDirPickerEvent&);
  void OnDirname(wxFileDirPickerEvent&);
  void OnResult(wxCommandEvent&);
  void OnOptions(wxCommandEvent&);
  void OnCreate(wxCommandEvent&);
  void OnUpdateDirpic(wxUpdateUIEvent&);
  void OnUpdateButt(wxUpdateUIEvent&);
  void OnClearBias(wxCommandEvent&);
  void OnClearDark(wxCommandEvent&);
  void OnClearFlat(wxCommandEvent&);
  wxString CreateResult(const wxString&) const;
  void OnFinish(wxProcessEvent&);
  void OnUpdate(wxTimerEvent&);

};
*/




class MuniImportRawOptions: public wxDialog
{
public:
  MuniImportRawOptions(wxWindow *, MuniConfig *);

  wxString GetType() const;
  wxString GetBand() const;
  wxString GetBitpix() const;
  wxString GetOverWrite() const;
  wxString GetInterpol() const;
  wxString GetDarkframe() const;
  wxString GetDcoptions() const;
  wxString GetDirs() const;

private:

  MuniConfig *config;
  wxChoice *filters,*interpols;
  wxCheckBox *overs;
  wxRadioButton *type0, *type1, *bitpix0, *bitpix1;
  wxFilePickerCtrl *darks;
  wxDirPickerCtrl *dirs;
  wxTextCtrl *dcopts;

  bool over,type_colour,type_grey,bitpix_16bit,bitpix_float;
  int filter,interpol;
  wxString dark,dcopt,dir;
  wxArrayString fchoices,fopt,ichoices,iopt;

  void Init();
  void CreateControls();
  void OnUpdateUI(wxUpdateUIEvent&);
  void OnDarks(wxFileDirPickerEvent&);
  void OnDirs(wxFileDirPickerEvent&);

};


/*
class MuniBrowserSearch: public wxSearchCtrl
{
public:
  MuniBrowserSearch(wxWindow *, wxWindowID, const wxString& =wxEmptyString,
		    const wxPoint& =wxDefaultPosition,
		    const wxSize& =wxDefaultSize, long =0);

  std::vector<long> Find(const std::vector<FitsMeta>&) const;

private:

  int type;
  wxString muster;

  std::vector<long> FindByName(const std::vector<FitsMeta>&) const;
  std::vector<long> FindByKey(const std::vector<FitsMeta>&) const;
  std::vector<long> FindByAdv(const std::vector<FitsMeta>&) const;

  void OnSearchEnter(wxCommandEvent&);
  void OnSearchButton(wxCommandEvent&);
  void OnSearchFinish(wxCommandEvent&);
  void OnSearchMenu(wxCommandEvent&);
  void OnSearchUpdate(wxCommandEvent&);
  void OnUpdateUI(wxUpdateUIEvent&);

};
*/

class MuniImportRaw: public wxDialog
{
private:

  int fcount;
  wxStaticText *label;
  wxGauge *gauge;
  MuniPipe pipe;
  wxString dcdark;

  void OnFinish(wxProcessEvent&);
  void OnClose(wxCloseEvent&);
  void OnCancel(wxCommandEvent&);
  void CreateControls();
  void CreatePipe(const MuniImportRawOptions&, const wxArrayString&);

public:

  MuniImportRaw(wxWindow *, const MuniImportRawOptions&, const wxArrayString&);
  virtual ~MuniImportRaw();

  void UpdateGauge(int);
  void LoadFile(const wxString&);

};



class MuniBrowser: public wxFrame
{
  MuniConfig *config;
  wxMenu *menuFile, *menuView, /*menuAct,*/ *menuArrange, *menuLabels, *menuHelp;
  MuniListCtrl *list;
  MuniNavigation *navbar;
  MuniPreferences *preferences;

  wxToolBar *tbar;
  wxBoxSizer *topsizer;
  wxArrayString errmsg;
  wxString pwd, mask;
  bool metaload, shutdown;

  void OnClose(wxCloseEvent&);
  void OnDirOpen(wxCommandEvent&);
  //  void FileOpen(wxCommandEvent& WXUNUSED(event));
  void OnPreferences(wxCommandEvent&);
  void OnClosePreferences(wxCloseEvent&);
  void OnStop(wxCommandEvent&);
  void OnRefresh(wxCommandEvent&);
  void OnProperties(wxCommandEvent&);
  void OnFileClose(wxCommandEvent&);
  void OnSelectAll(wxCommandEvent&);
  void OnCut(wxCommandEvent&);
  void OnCopy(wxCommandEvent&);
  void OnPaste(wxCommandEvent&);
  void SelectItem(wxCommandEvent&);
  void OnIconList(wxCommandEvent&);
  void OnNewBrowser(wxCommandEvent&);
  void OnNewView(wxCommandEvent&);
  void OnActivated(wxNotifyEvent&);
  void OnShowToolbar(wxCommandEvent&);
  void FindStars(wxCommandEvent& WXUNUSED(event));
  void AperturePhot(wxCommandEvent& WXUNUSED(event));
  void ProfilePhot(wxCommandEvent& WXUNUSED(event));
  void Matching(wxCommandEvent& WXUNUSED(event));
  void Astrometry(wxCommandEvent& WXUNUSED(event));
  void Stacking(wxCommandEvent& WXUNUSED(event));
  void Deconvolution(wxCommandEvent& WXUNUSED(event));
  void OnAverage(wxCommandEvent& WXUNUSED(event));
  void OnDarkbat(wxCommandEvent& WXUNUSED(event));
  void HelpAbout(wxCommandEvent& WXUNUSED(event));
  void OnMetaProgress(MetaProgressEvent&);
  void OnLabel(wxCommandEvent&);
  void OnSort(wxCommandEvent&);
  void OnReverse(wxCommandEvent&);

  void OnEnterNavigator(NavigationEvent&);


  //  void ImportRaw(const MuniImportRawOptions&,const wxArrayString&);

  void OnUpdateCut(wxUpdateUIEvent&);
  void OnUpdatePaste(wxUpdateUIEvent&);
  void OnUpdateStop(wxUpdateUIEvent&);
  void OnUpdateOpen(wxUpdateUIEvent&);

  void OpenView(const wxString&);

public:
  MuniBrowser(wxWindow *, MuniConfig *);

  void FileLoad(const wxString&);
  void FileLoad(const wxArrayString&);
};
