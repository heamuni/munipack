/*

  XMunipack -- list for Browser

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
#include "fits.h"
#include "enum.h"
#include "event.h"
#include "dirscan.h"
#include "icon.h"
#include <wx/wx.h>
#include <wx/listctrl.h>
#include <wx/imaglist.h>
#include <vector>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif



class MetaRender: public wxThread
{
  wxEvtHandler *handler;
  const size_t id;
  const wxArrayString files;
  const wxSize size;
  std::vector<FitsMetaHdu> hdu;
  FitsArray array;

  ExitCode Entry();
  bool MagicFile(const wxString&) const;
  bool FitsOpen(const wxString&, bool pick_array);

public:
  MetaRender(wxEvtHandler *, const size_t, const wxArrayString&, const wxSize&);

};


class MuniListCtrl: public wxListCtrl
{
  DirScan *dirscan;
  size_t DirScanId;
  long selected_item;

protected:

  size_t metaId;
  MuniConfig *config;
  wxImageList thumbs;
  std::vector<FitsMeta> metalist;
  wxThread *metarender;
  wxBitmap default_icon;

  void OnRightClick(wxListEvent&);
  void OnClick(wxListEvent&);
  void OnView(wxCommandEvent& WXUNUSED(event));
  void OnProperties(wxCommandEvent& WXUNUSED(event));
  void SortItems();
  void OnMouse(wxMouseEvent&);
  void OnCut(wxCommandEvent&);
  void OnCopy(wxCommandEvent&);
  void OnPaste(wxCommandEvent&);
  void OnSelall(wxCommandEvent&);
  virtual void OnMetaRender(MetaOpenEvent&) = 0;
  void OnDirScan(DirScanEvent&);
  void AddItem(const wxString&);
  void StartMetarender(const wxArrayString&);
  void SendActivate(long);

public:
  MuniListCtrl(wxWindow *,wxWindowID, long, MuniConfig *);

  virtual void Stop();
  virtual void ChangeDir(const wxString&, const wxString&);

  virtual void FitsLoad(const wxArrayString&);
  virtual void SelectAll();
  virtual void DeSelectAll();

  virtual void PasteMeta(const FitsMeta&);
  virtual void PasteMeta(const std::vector<FitsMeta>&);

  //  std::vector<FitsMeta> GetClipboard() const;
  wxArrayString GetClipboard() const;
  void SetClipboard(const wxArrayString&);
  //  void SetClipboard(const std::vector<FitsMeta>&);
  std::vector<FitsMeta> GetMetaList() const { return metalist; }

  std::vector<unsigned int> GetSelectedIndex() const;
  void SelectItemLast();
  void SelectItem(long);
  void SelectItemRelative(long);

  virtual void Cut();
  virtual void Copy();
  virtual void Paste();
  virtual void SelectLabel() = 0;
  virtual void Sort() = 0;

};

class MuniListIcon: public MuniListCtrl
{
  wxString LabelFits(const FitsMeta&, int) const;
  void OnMetaRender(MetaOpenEvent&);
  void OnMetaRenderFinish(MetaOpenEvent&);

public:

  MuniListIcon(wxWindow *, wxWindowID, MuniConfig *);
  virtual void SelectLabel();
  virtual void Sort();

};

class MuniListList: public MuniListCtrl
{

  void OnMetaRender(MetaOpenEvent&);
  void OnMetaRenderFinish(MetaOpenEvent&);
  void OnSize(wxSizeEvent&);
  void SetItemMeta(long,const FitsMeta&);

public:
  MuniListList(wxWindow *, wxWindowID, MuniConfig *);
  virtual void SelectLabel();
  virtual void Sort();

};
