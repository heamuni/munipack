/*

  xmunipack - preview icon list


  Copyright © 2009-2013, 2019-2022 F.Hroch (hroch@physics.muni.cz)

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

#include "list.h"
#include "fileprop.h"
#include <vector>
#include <wx/wx.h>
#include <wx/utils.h>
#include <wx/thread.h>
#include <wx/filename.h>
#include <wx/uri.h>
#include <wx/filesys.h>
#include <wx/graphics.h>
#include <wx/clipbrd.h>
#include <wx/dnd.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

using namespace std;

// keys for sorting
static wxString sort_key, sort_dateobs, sort_object, sort_filter,
  sort_exposure;


static bool CmpName(const FitsMeta& meta1, const FitsMeta& meta2)
{
  wxString name1 = meta1.GetName();
  wxString name2 = meta2.GetName();
  return wxStrcmp(name1,name2) < 0;
}

static bool CmpObject(const FitsMeta& meta1, const FitsMeta& meta2)
{
  wxString key1 = meta1.GetKeys(sort_object);
  wxString key2 = meta2.GetKeys(sort_object);
  return wxStrcmp(key1,key2) < 0;
}

static bool CmpFilter(const FitsMeta& meta1, const FitsMeta& meta2)
{
  wxString key1 = meta1.GetFilter(sort_filter);
  wxString key2 = meta2.GetFilter(sort_filter);
  return wxStrcmp(key1,key2) < 0;
}

static bool CmpKey(const FitsMeta& meta1, const FitsMeta& meta2)
{
  wxString key1 = meta1.GetKeys(sort_key);
  wxString key2 = meta2.GetKeys(sort_key);
  return wxStrcmp(key1,key2) < 0;
}

static bool CmpExposure(const FitsMeta& meta1, const FitsMeta& meta2)
{
  wxString key1 = meta1.GetKeys(sort_exposure);
  wxString key2 = meta2.GetKeys(sort_exposure);
  double e1,e2;
  if( key1.ToDouble(&e1) && key2.ToDouble(&e2) )
    return e1 < e2;
  else
    return false;
}

static bool CmpSize(const FitsMeta& meta1, const FitsMeta& meta2)
{
  wxULongLong s1 = meta1.GetSize();
  wxULongLong s2 = meta2.GetSize();
  return s1 < s2;
}

static bool CmpJulday(const FitsMeta& meta1, const FitsMeta& meta2)
{
  FitsTime t1(meta1.GetDateobs(sort_dateobs));
  FitsTime t2(meta2.GetDateobs(sort_dateobs));
  return t1.GetJd() < t2.GetJd();
}



// -- DropClass

// add files drop

/*
class xDropTarget: public wxDropTarget
{
public:
  xDropTarget(wxWindow *w): target(w) {
    SetDataObject(new MuniDataObjectMeta);
  }

  wxDragResult OnData(wxCoord x, wxCoord y, wxDragResult def)
  {
    if ( !GetData() )
      return wxDragNone;

    MuniDataObjectMeta *dobj = (MuniDataObjectMeta *) GetDataObject();
    wxASSERT(dobj);
    vector<FitsMeta> slist = dobj->GetMetafitses();
    if( slist.size() > 0 )
      static_cast<MuniListCtrl *>(target)->PasteMeta(slist);

    return wxDragCopy;
  }

private:

  wxWindow *target;
};
*/


// ---- MuniListIcon

MuniListIcon::MuniListIcon(wxWindow *w, wxWindowID id, MuniConfig *c):
  MuniListCtrl(w,id,wxLC_ICON,c)
{
  int size = config->icon_size;
  thumbs.Create(size,size);
  MuniIcon iengine(config);
  default_icon = iengine.Utf8Icon(2.0,L"❦");
  thumbs.Add(wxBitmap(default_icon));
  SetImageList(&thumbs,wxIMAGE_LIST_NORMAL);

  Bind(EVT_META_OPEN,&MuniListIcon::OnMetaRender,this,ID_METALOAD_PHASE);
  Bind(EVT_META_OPEN,&MuniListIcon::OnMetaRenderFinish,this,ID_METALOAD_FINISH);
}

void MuniListIcon::OnMetaRender(MetaOpenEvent& event)
{
  //  wxLogDebug("MuniListIcon::OnMetaRender()");

  if( event.id != metaId ) return;

  size_t i = event.index;

  MuniIcon iengine(config);
  if( event.meta.IsOk() ) {
    metalist[i] = event.meta;
    if( event.icon.IsOk() )
      metalist[i].SetIcon(event.icon);
    else
      metalist[i].SetIcon(iengine.Utf8Icon(2.0,config->head_symbol).ConvertToImage());
    wxString label = LabelFits(metalist[i],config->browser_labeltype);
    wxBitmap icon = iengine.GetCover(metalist[i].GetIcon(),label,
				     event.meta.HduCount());
    thumbs.Replace(i,icon);
  }
  else {
    wxBitmap icon = iengine.Utf8Icon(2.0,L"☙",
				     wxURI::Unescape(metalist[i].GetFullName()));
    thumbs.Replace(i,icon);
  }
  SetItemImage(i,i);

  MetaProgressEvent ev(EVT_META_PROGRESS,ID_METALOAD_PHASE);
  ev.phase = double(i) / double(metalist.size());
  wxQueueEvent(this,ev.Clone());
}

void MuniListIcon::OnMetaRenderFinish(MetaOpenEvent& event)
{
  //  wxLogDebug("MuniListIcon::OnMetaRenderFinish()");

  if( event.id != metaId ) return;

  metarender = 0;

  MetaProgressEvent ev(EVT_META_PROGRESS,ID_METALOAD_FINISH);
  wxQueueEvent(this,ev.Clone());

  Sort();
  return;
}


wxString MuniListIcon::LabelFits(const FitsMeta& f, int type) const
{
  switch (type) {
  case ID_LABEL_NOLABEL:  return "";
  case ID_LABEL_FILENAME: return wxURI::Unescape(f.GetFullName());
  case ID_LABEL_OBJECT:   return f.GetKeys(config->fits_key_object);
  case ID_LABEL_DATE:     return f.GetDate(config->fits_key_dateobs);
  case ID_LABEL_TIME:     return f.GetTime(config->fits_key_dateobs);
  case ID_LABEL_FILTER:   return f.GetFilter(config->fits_key_filter);
  case ID_LABEL_EXPOSURE: return f.GetExposure(config->fits_key_exptime);
  case ID_LABEL_KEY:      return f.GetKeys(config->browser_labelkey);
  default:                return "...";
  }
}

void MuniListIcon::SelectLabel()
{
  for(size_t i = 0; i < metalist.size(); i++) {
    if( metalist[i].IsOk() ) {
      wxString label = LabelFits(metalist[i],config->browser_labeltype);
      MuniIcon iengine(config);
      wxBitmap icon = iengine.GetCover(metalist[i].GetIcon(),label,
				       metalist[i].HduCount());
      thumbs.Replace(i,icon);
      SetItemImage(i,i);
    }
  }
}


void MuniListIcon::Sort()
{
  DeleteAllItems();
  thumbs.RemoveAll();

  SortItems();

  for(size_t i = 0; i < metalist.size(); i++) {

    wxBitmap icon;
    MuniIcon iengine(config);
    if( metalist[i].GetIcon().IsOk() ) {
      wxString label = LabelFits(metalist[i],config->browser_labeltype);
      icon = iengine.GetCover(metalist[i].GetIcon(),label,metalist[i].HduCount());
    }
    else
      icon = iengine.Utf8Icon(2.0,L"☙",wxURI::Unescape(metalist[i].GetFullName()));
    thumbs.Add(icon);

    wxListItem item;
    item.SetId(i);
    item.SetData(i);
    item.SetImage(i);
    item.SetMask(wxLIST_MASK_IMAGE|wxLIST_MASK_DATA);
    InsertItem(item);
  }
}


// ---- MuniListList

MuniListList::MuniListList(wxWindow *w, wxWindowID id, MuniConfig *c):
  MuniListCtrl(w,id,wxLC_REPORT|wxLC_VRULES,c)
{
  InsertColumn(0,"File name",wxLIST_FORMAT_LEFT,wxLIST_AUTOSIZE);
  InsertColumn(1,"Object",wxLIST_FORMAT_LEFT,wxLIST_AUTOSIZE_USEHEADER);
  InsertColumn(2,"Filter",wxLIST_FORMAT_LEFT,wxLIST_AUTOSIZE_USEHEADER);
  InsertColumn(3,"Date",wxLIST_FORMAT_RIGHT,wxLIST_AUTOSIZE_USEHEADER);
  InsertColumn(4,"Time",wxLIST_FORMAT_RIGHT,wxLIST_AUTOSIZE_USEHEADER);
  InsertColumn(5,"Exposure",wxLIST_FORMAT_RIGHT,wxLIST_AUTOSIZE_USEHEADER);
  InsertColumn(6,"Size",wxLIST_FORMAT_RIGHT,wxLIST_AUTOSIZE_USEHEADER);

  int wch = GetCharWidth();
  SetColumnWidth(3,11*wch);
  SetColumnWidth(4,12*wch);

  Bind(wxEVT_SIZE,&MuniListList::OnSize,this);
  Bind(EVT_META_OPEN,&MuniListList::OnMetaRender,this,ID_METALOAD_PHASE);
  Bind(EVT_META_OPEN,&MuniListList::OnMetaRenderFinish,this,ID_METALOAD_FINISH);
}

void MuniListList::SelectLabel()
{
  wxFAIL_MSG("MuniListIcon::Label(): Unreachable code.");
}

void MuniListList::OnSize(wxSizeEvent& event)
{
  int w = 0;
  for(int i = 1; i < GetColumnCount(); i++)
    w += GetColumnWidth(i);

  wxSize size(GetClientSize());
  SetColumnWidth(0,size.GetWidth()-w);

  event.Skip();
}


void MuniListList::OnMetaRender(MetaOpenEvent& event)
{
  if( event.id != metaId ) return;

  size_t i = event.index;

  metalist[i] = event.meta;
  if( event.meta.IsOk() ) {
    SetItemMeta(i,metalist[i]);
  }
  else {
    wxColour c = wxSystemSettings::GetColour(wxSYS_COLOUR_GRAYTEXT);
    SetItemTextColour(i,c);
  }

  MetaProgressEvent ev(EVT_META_PROGRESS,ID_METALOAD_PHASE);
  ev.phase = double(i) / double(metalist.size());
  wxQueueEvent(this,ev.Clone());

}


void MuniListList::OnMetaRenderFinish(MetaOpenEvent& event)
{
  wxASSERT(event.GetId() == ID_METALOAD_FINISH);
  if( event.id != metaId ) return;

  metarender = 0;
  SetColumnWidth(6,wxLIST_AUTOSIZE);

  int w = 0;
  for(int i = 1; i < GetColumnCount(); i++)
    w += GetColumnWidth(i);

  wxSize size(GetClientSize());
  SetColumnWidth(0,size.GetWidth()-w);

  MetaProgressEvent ev(EVT_META_PROGRESS,ID_METALOAD_FINISH);
  wxQueueEvent(this,ev.Clone());

  Sort();
  return;
}

void MuniListList::Sort()
{
  if( config->browser_sorttype == ID_SORT_UNSORT ) return;

  SortItems();

  DeleteAllItems();

  for(size_t i = 0; i < metalist.size(); i++) {
    //    wxLogDebug("sort:"+metalist[i].GetURL()+" "+metalist[i].GetFullName());

    wxListItem item;
    item.SetId(i);
    item.SetData(i);
    item.SetText(wxURI::Unescape(metalist[i].GetFullName()));
    item.SetMask(wxLIST_MASK_TEXT|wxLIST_MASK_DATA);
    InsertItem(item);

    if( metalist[i].IsOk() )
      SetItemMeta(i,metalist[i]);
    else {
      wxColour c = wxSystemSettings::GetColour(wxSYS_COLOUR_GRAYTEXT);
      SetItemTextColour(i,c);
    }
  }
}

void MuniListList::SetItemMeta(long i, const FitsMeta& meta)
{
  SetItem(i,0,wxURI::Unescape(meta.GetFullName()));
  SetItem(i,1,meta.GetKeys(config->fits_key_object));
  SetItem(i,2,meta.GetFilter(config->fits_key_filter));
  SetItem(i,3,meta.GetDate(config->fits_key_dateobs));
  SetItem(i,4,meta.GetTime(config->fits_key_dateobs));
  SetItem(i,5,meta.GetExposure(config->fits_key_exptime));
  SetItem(i,6,meta.GetHumanReadableSize());
}

// ---- MuniListCtrl


MuniListCtrl::MuniListCtrl(wxWindow *w, wxWindowID id, long style, MuniConfig *c):
  wxListCtrl(w,id,wxDefaultPosition,wxDefaultSize,style),dirscan(0),DirScanId(007),
  selected_item(-1),metaId(007),
  config(c), metarender(0)
{
  //  SetDropTarget(new xDropTarget(this));

  Bind(wxEVT_MOTION,&MuniListCtrl::OnMouse,this);
  Bind(wxEVT_COMMAND_LIST_ITEM_RIGHT_CLICK,&MuniListCtrl::OnRightClick,this);
  Bind(wxEVT_COMMAND_LIST_ITEM_ACTIVATED,&MuniListCtrl::OnClick,this);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniListCtrl::OnView,this,ID_VIEW);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniListCtrl::OnCut,this,wxID_CUT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniListCtrl::OnCopy,this,wxID_COPY);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniListCtrl::OnPaste,this,wxID_PASTE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniListCtrl::OnProperties,this,wxID_PROPERTIES);
  Bind(EVT_DIR_SCAN,&MuniListCtrl::OnDirScan,this);
}

void MuniListCtrl::Stop()
{
  if( dirscan )
    dirscan->Delete();
  if( metarender )
    metarender->Delete();
}

// vector<FitsMeta> MuniListCtrl::GetClipboard() const
// {
//   vector<FitsMeta> slist;

//   wxDataFormat dfm;
//   dfm.SetId("MUNIPACK_METAFITS");

//   if( wxTheClipboard->IsOpened() )
//     return slist;
//   // this test is required by UI updating in browser

//   if( wxTheClipboard->Open() ) {

//     if( wxTheClipboard->IsSupported(dfm) ) {

//       MuniDataObjectMeta data;
//       if( wxTheClipboard->GetData(data) )
// 	slist = data.GetMetafitses();
//       else
// 	wxLogDebug("MuniList::EditPaste: failed to paste ");

//     }
//     wxTheClipboard->Close();
//   }

wxArrayString MuniListCtrl::GetClipboard() const
{
  wxArrayString files;
  if (wxTheClipboard->Open() ) {
    if( wxTheClipboard->IsSupported(wxDF_FILENAME) ) {

      wxFileDataObject data;
      if( wxTheClipboard->GetData(data) ) {
	files = data.GetFilenames();
	for(size_t i = 0; i < files.GetCount(); i++)
	  wxLogDebug(files[i]);
      }
    }
    wxTheClipboard->Close();
  }
  return files;
}

void MuniListCtrl::SetClipboard(const wxArrayString& files)
{
  if( files.IsEmpty() ) return;

  if( wxTheClipboard->Open() ) {
    //  MuniDataObjectMeta *fs = new MuniDataObjectMeta(slist);
    //    wxTheClipboard->SetData(files);
    wxTheClipboard->Close();
  }
}

/*
void MuniListCtrl::SetClipboard(const vector<FitsMeta>& slist)
{
  if( ! (slist.size() > 0) ) return;

  if( wxTheClipboard->Open() ) {
    MuniDataObjectMeta *fs = new MuniDataObjectMeta(slist);
    wxTheClipboard->SetData(fs);
    wxTheClipboard->Close();
  }
}
*/

void MuniListCtrl::FitsLoad(const wxArrayString& fitses)
{
  if( fitses.GetCount() == 0 ) return;

  for(size_t i = 0; i < fitses.GetCount(); i++)
    AddItem(fitses[i]);

  // start metarender
  StartMetarender(fitses);
}


void MuniListCtrl::OnRightClick(wxListEvent& event)
{
  if( event.GetIndex() >= 0 ) {
    selected_item = event.GetIndex();

    wxMenu *menuMark = new wxMenu;
    menuMark->AppendRadioItem(ID_MARK_SCI,"Scientific (light)");
    menuMark->AppendRadioItem(ID_MARK_FLAT,"Flat field");
    menuMark->AppendRadioItem(ID_MARK_DARK,"Dark frame");
    menuMark->AppendRadioItem(ID_MARK_BIAS,"Bias frame");

    wxMenu popup;
    popup.Append(ID_VIEW,"&View");
    /*
    popup.AppendSeparator();
    popup.Append(wxID_CUT);
    popup.Append(wxID_COPY);
    */
    //    popup.AppendSeparator();
    //    popup.AppendSubMenu(menuMark,"&Mark as");
    //    popup.AppendSeparator();
    popup.Append(wxID_PROPERTIES);
    /*
    menuMark->Enable(ID_MARK_SCI,false);
    menuMark->Enable(ID_MARK_FLAT,false);
    menuMark->Enable(ID_MARK_DARK,false);
    menuMark->Enable(ID_MARK_BIAS,false);
    */
    bool isok = metalist[selected_item].IsOk();
    popup.Enable(ID_VIEW,GetSelectedItemCount() == 1 && isok);
    popup.Enable(wxID_PROPERTIES,GetSelectedItemCount() == 1 && isok);

    PopupMenu(&popup);
  }
  else {

    wxMenu popup;
    popup.Append(wxID_PASTE);
    popup.Enable(wxID_PASTE,GetClipboard().GetCount()/*size()*/ > 0);
    PopupMenu(&popup);

  }

}

void MuniListCtrl::OnClick(wxListEvent& event)
{
  SendActivate(event.GetIndex());
}

void MuniListCtrl::OnView(wxCommandEvent& event)
{
  wxASSERT(0 <= selected_item && selected_item < long(metalist.size()));
  SendActivate(selected_item);
}

void MuniListCtrl::SendActivate(long idx)
{
  wxASSERT(0 < idx && idx < (long) metalist.size());
  if( metalist[idx].IsOk() ) {
    wxNotifyEvent event(EVT_LIST_ACTIVATED,ID_VIEW_FILE);
    event.SetString(wxURI::Unescape(metalist[idx].GetFullPath()));
    wxQueueEvent(this,event.Clone());
  }
}

void MuniListCtrl::OnProperties(wxCommandEvent& WXUNUSED(event))
{
  wxASSERT(0 <= selected_item && selected_item < long(metalist.size()));
  FitsMeta meta(metalist[selected_item]);

  MuniIcon iengine(config);
  wxImage cover;
  std::vector<wxImage> ilist(iengine.GetList(meta));

  if( InReportView() )
    // no icons already generated
    cover = iengine.Utf8Icon(2.0,L"❦").ConvertToImage();
  else
    cover = iengine.GetCover(meta.GetIcon());

  for(size_t i = 0; i < meta.HduCount(); i++)
    if( meta.Hdu(i).Type() == HDU_IMAGE ) {
      ilist[i] = cover;
      break;
    }

  MuniFileProperties *w = new MuniFileProperties(this,config,meta,cover,ilist);
  w->Show();
}

void MuniListCtrl::SelectAll()
{
  for(long i = 0; i < GetItemCount(); i++)
    SetItemState(i,wxLIST_STATE_SELECTED,wxLIST_STATE_SELECTED);
}

void MuniListCtrl::SelectItemLast()
{
  long i = GetItemCount();
  if( i > 0 )
    SetItemState(i-1,wxLIST_STATE_SELECTED,wxLIST_STATE_SELECTED);
}

void MuniListCtrl::SelectItem(long i)
{
  //  wxASSERT(0 <= i && i < GetItemCount() && GetItemCount() ==long(flist.size()));
  SetItemState(i,wxLIST_STATE_SELECTED,wxLIST_STATE_SELECTED);
}

void MuniListCtrl::SelectItemRelative(long dir)
{
  long n = GetNextItem(-1,wxLIST_NEXT_ALL,wxLIST_STATE_SELECTED);
  if( n > -1 ) {
    long i = n + dir;
    if( 0 <= i && i < GetItemCount() ) {
      DeSelectAll();
      SetItemState(i,wxLIST_STATE_SELECTED,wxLIST_STATE_SELECTED);
      wxQueueEvent(this,new wxListEvent(wxEVT_COMMAND_LIST_ITEM_ACTIVATED));
    }
  }
}


void MuniListCtrl::DeSelectAll()
{
  for(long i = 0; i < GetItemCount(); i++)
    SetItemState(i,0,wxLIST_STATE_SELECTED);
}


void MuniListCtrl::SortItems()
{
  sort_key = config->browser_labelkey;
  sort_object = config->fits_key_object;
  sort_dateobs = config->fits_key_dateobs;
  sort_filter = config->fits_key_filter;
  sort_exposure = config->fits_key_exptime;

  switch (config->browser_sorttype) {
  case ID_SORT_FILENAME: sort(metalist.begin(),metalist.end(),CmpName); break;
  case ID_SORT_OBJECT:   sort(metalist.begin(),metalist.end(),CmpObject); break;
  case ID_SORT_DATEOBS:  sort(metalist.begin(),metalist.end(),CmpJulday); break;
  case ID_SORT_FILTER:   sort(metalist.begin(),metalist.end(),CmpFilter); break;
  case ID_SORT_EXPOSURE: sort(metalist.begin(),metalist.end(),CmpExposure); break;
  case ID_SORT_SIZE:     sort(metalist.begin(),metalist.end(),CmpSize);  break;
  case ID_SORT_KEY:      sort(metalist.begin(),metalist.end(),CmpKey);  break;
  }

  if( config->browser_reverse )
    reverse(metalist.begin(),metalist.end());
}




void MuniListCtrl::OnMouse(wxMouseEvent& event)
{
  if( event.Dragging() ) {
    //    wxLogDebug(_("MuniListCtrls::OnMouse dragging"));

    // vector<FitsMeta> slist = GetSelectedMeta();
    // if( slist.size() > 0 ) {
    //   //      MuniDataObjectMeta fs(slist);

    //   /*
    //   FitsMeta f(slist[0]);
    //   wxImage ico(f.GetIcon());
    //   wxIconOrCursor ic = wxDROP_ICON(
    //   */
    //   /*
    //   wxDropSource dragSource(fs,this);
    //   dragSource.DoDragDrop();
    //   */
    // }
  }

  event.Skip();
}



void MuniListCtrl::PasteMeta(const FitsMeta& meta)
{
}

void MuniListCtrl::PasteMeta(const vector<FitsMeta>& slist)
{
  for(vector<FitsMeta>::const_iterator i = slist.begin(); i != slist.end(); ++i)
    PasteMeta(*i);
}

void MuniListCtrl::OnCopy(wxCommandEvent& event)
{
  Copy();
}

void MuniListCtrl::OnCut(wxCommandEvent& event)
{
  Cut();
}

void MuniListCtrl::OnPaste(wxCommandEvent& event)
{
  Paste();
}

void MuniListCtrl::Cut()
{

  /*
  SetClipboard(slist);
  */
}

void MuniListCtrl::Copy()
{
  /*
  SetClipboard(slist);
  */
}

void MuniListCtrl::Paste()
{
  //  PasteMeta(GetClipboard());
}



void MuniListCtrl::ChangeDir(const wxString& path, const wxString& mask)
{
  DeleteAllItems();
  thumbs.RemoveAll();
  metalist.clear();

  if( dirscan && dirscan->IsRunning()  )
    dirscan->Delete();
  if( metarender )
    metarender->Delete();

  DirScanId++;
  dirscan = new DirScan(this,path,mask,DirScanId);
  if( dirscan->Create() == wxTHREAD_NO_ERROR ) {
    wxThreadError code = dirscan->Run();
    if( !(code == wxTHREAD_NO_ERROR) )
      wxLogError("MuniListCtrl: Can't run a file scan thread.");
  }
  else {
    delete dirscan;
    dirscan = 0;
  }
}

void MuniListCtrl::OnDirScan(DirScanEvent& event)
{
  //  wxLogDebug("MuniListCtrl::OnDirScan() %d %d",event.finish,int(event.id));

  if( event.id != DirScanId ) return;

  wxASSERT((size_t)GetItemCount() == metalist.size());

  // add icon
  for(size_t i = 0; i < event.files.GetCount(); i++)
    AddItem(event.files[i]);

  // finish
  if( event.finish ) {
    dirscan = 0;

    // start metarender
    StartMetarender(event.files);
  }
}

void MuniListCtrl::AddItem(const wxString& file)
{
  int n = GetItemCount();

  if( ! InReportView() )
    thumbs.Add(wxBitmap(default_icon));

  FitsMeta meta(file);
  metalist.push_back(meta);
  wxListItem item;
  item.SetId(n);

  wxFileName filename(file);
  wxString label(wxURI::Unescape(filename.GetFullName()));

  item.SetData(n);
  if( InReportView() ) {
    item.SetText(label);
    item.SetMask(wxLIST_MASK_TEXT|wxLIST_MASK_DATA);
  }
  else {// icons, no labels by wxListCrtl
    item.SetImage(n);
    item.SetMask(wxLIST_MASK_IMAGE|wxLIST_MASK_DATA);
  }

  InsertItem(item);
}

void MuniListCtrl::StartMetarender(const wxArrayString& files)
{
  if( metarender ) {
    metarender->Delete();
    MetaProgressEvent ev(EVT_META_PROGRESS,ID_METALOAD_FINISH);
    wxQueueEvent(this,ev.Clone());
  }

  metaId++;
  int n = InReportView() ? 0 : config->icon_size;
  metarender = new MetaRender(this,metaId,files,wxSize(n,n));
  if( metarender->Create() == wxTHREAD_NO_ERROR ) {
    wxThreadError code = metarender->Run();
    if( !(code == wxTHREAD_NO_ERROR) )
      wxLogError("MuniListCtrl: Can't run a metarender thread.");
  }
  else {
    delete metarender;
    metarender = 0;
  }

  MetaProgressEvent ev(EVT_META_PROGRESS,ID_METALOAD_START);
  wxQueueEvent(this,ev.Clone());

}
