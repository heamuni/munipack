/*

  xmunipack - side panel

  Copyright © 1997-2013 F.Hroch (hroch@physics.muni.cz)

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
//#include "plot.h"
#include <cfloat>
#include <vector>
#include <wx/wx.h>
#include <wx/dcbuffer.h>
#include <wx/graphics.h>

using namespace std;

//#define SIZE 22

/*
BEGIN_EVENT_TABLE(MuniSide, wxSplitterWindow)
  EVT_TREE_SEL_CHANGED(wxID_ANY, MuniSide::OnTreeActivated)
  EVT_SIDE(ID_SIDE_HDU, MuniSide::OnHduSelect)
END_EVENT_TABLE()
*/

// the tree node storage:
//   *  index is index in FitsFile
//   *  type is for derived data, color images etc.
//   *  the tree hubs has no data defined

class xTreeItemData: public wxTreeItemData
{
public:
  xTreeItemData(const int n, const int t=0): index(n), type(t) {}

  int index, type;
};

// class xPanel : public wxPanel
// {
// public:
//   xPanel(wxWindow *w);
//   void AssignHistos(const std::vector<FitsHisto> &);

// private:

//   std::vector<FitsHisto> hlist;
//   int hmax;
//   wxSize last;

//   void OnSize(wxSizeEvent&);
//   void OnPaint(wxPaintEvent&);

//   DECLARE_EVENT_TABLE()
// };

// BEGIN_EVENT_TABLE(xPanel, wxPanel)
// EVT_PAINT(xPanel::OnPaint)
// EVT_SIZE(xPanel::OnSize)
// END_EVENT_TABLE()


// --- MuniSide

MuniSide::MuniSide():wxTreeCtrl(),topwin(0),config(0) {}

MuniSide::MuniSide(wxWindow *w, const wxPoint& pos, const wxSize& size, 
		   MuniConfig *c):
  wxTreeCtrl(w,wxID_ANY,pos,size,(wxTR_DEFAULT_STYLE | wxTR_HAS_BUTTONS | wxTR_HIDE_ROOT | wxTR_FULL_ROW_HIGHLIGHT) & ~wxTR_NO_LINES),
  topwin(w->GetParent()), config(c)
{
  wxASSERT(topwin);
  Bind(wxEVT_COMMAND_TREE_SEL_CHANGED,&MuniSide::OnTreeActivated,this);
  Bind(xEVT_SIDE,&MuniSide::OnHduSelect,this,ID_SIDE_HDU);
}

bool MuniSide::Create(wxWindow *w, const wxPoint& pos, const wxSize& size, 
		   MuniConfig *c)
{
  if( wxTreeCtrl::Create(w,wxID_ANY,pos,size,(wxTR_DEFAULT_STYLE | wxTR_HAS_BUTTONS | wxTR_HIDE_ROOT | wxTR_FULL_ROW_HIGHLIGHT) & ~wxTR_NO_LINES) ) {

    topwin = w->GetParent();
    config = c;

    Bind(wxEVT_COMMAND_TREE_SEL_CHANGED,&MuniSide::OnTreeActivated,this);
    Bind(xEVT_SIDE,&MuniSide::OnHduSelect,this,ID_SIDE_HDU);

    return true;
  }
  else
    return false;
}



// MuniSide::MuniSide(wxWindow *w, const wxPoint& pos, const wxSize& size, 
// 		   MuniConfig *c):
// wxSplitterWindow(w,wxID_ANY,pos,size), topwin(w->GetParent()), config(c)
// {
//   wxASSERT(topwin);

//   SetSashGravity(1.0);

//   // tree
//   fstruct = new wxTreeCtrl(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
//    (wxTR_DEFAULT_STYLE | wxTR_HAS_BUTTONS | wxTR_HIDE_ROOT | wxTR_FULL_ROW_HIGHLIGHT) & ~wxTR_NO_LINES);

//   // infopage
//   infopage = new wxPanel(this);
//   infopage->Show(config->sideinfo_show);

//   wxASSERT(fstruct && infopage);

//   if( config->sideinfo_show ) {
//     wxSize s(infopage->GetBestSize());
//     SplitHorizontally(fstruct,infopage,-s.GetHeight());
//   }
//   else
//     Initialize(fstruct);
  
//   // icons
//   /*
//   int ssize = config->info_symbol.GetWidth();
//   wxImageList *icons = new wxImageList(ssize, ssize, true);
//   icons->Add(wxBitmap(config->info_symbol));
//   icons->Add(wxBitmap(config->picture_symbol));
//   icons->Add(wxBitmap(config->table_symbol));
//   icons->Add(wxBitmap(config->color_symbol));
//   fstruct->AssignImageList(icons);
//   */
// }

void MuniSide::Clear()
{
  //  fits.Clear();
  meta.Clear();
  //  hists.clear();
  /*  fstruct->*/DeleteAllItems();
  //  infopage->DestroyChildren();
}

// void MuniSide::UpdateInfoPage(int hdusel)
// {
//   //  wxLogDebug(_("%d"),hdusel);
//   //  return;
//   if( ! IsSplit() ) return;

//   wxASSERT(infopage);

//   infopage->DestroyChildren();

//   wxFont sf = wxFont(*wxSMALL_FONT);

//   wxBoxSizer *infosizer = 0;

//   /*
//   wxTreeItemId itemId = fstruct->GetSelection();
//   wxASSERT(itemId.IsOk());
//   xTreeItemData item = static_cast<xTreeItemData>(itemId);
//   */

//   if( hdusel == -1 && fits.Type() == FITS_COLOR )
//     //  if( hdusel == -1 && meta.Type() == FITS_COLOR )
//     infosizer = UpdateInfoColor(infopage);
//   else if( hdusel != -1 ) {
//     const FitsHdu hdu = fits.Hdu(hdusel);
//     //const FitsHdu hdu = meta.Hdu(hdusel);
//     if( hdu.Type() == HDU_IMAGE )
//       infosizer = UpdateInfoImage(infopage,hdu,hists[hdusel]);
//     else if( hdu.Type() == HDU_TABLE )
//       infosizer = UpdateInfoTable(infopage,hdu);
//   }
  
//   if( ! infosizer ) 
//     infosizer = UpdateInfoEmpty(infopage);

//   infopage->SetSizerAndFit(infosizer);

//   wxSize s(infopage->GetBestSize());
//   SetSashPosition(-s.GetHeight());
// }


// wxBoxSizer *MuniSide::UpdateInfoEmpty(wxPanel *infopage)
// {
//   wxBoxSizer *infosizer = new wxBoxSizer(wxVERTICAL);
//   infosizer->Add(new wxStaticText(infopage,wxID_ANY,wxT("")));
//   return infosizer;
// }

// wxBoxSizer *MuniSide::UpdateInfoTable(wxPanel *infopage, const FitsHdu& hdu)
// {
//   wxFont sf = wxFont(*wxSMALL_FONT);

//   const FitsTable table(hdu);

//   wxBoxSizer *infosizer = new wxBoxSizer(wxVERTICAL);

//   wxString a;
//   a.Printf(wxT("%d rows"),(int) table.Height());
//   wxStaticText *trow = new wxStaticText(infopage,wxID_ANY,a);
//   trow->SetFont(sf);
//   infosizer->Add(trow,wxSizerFlags().Center());
//   a.Printf(wxT("%d columns"),(int) table.Width());
//   wxStaticText *tcols = new wxStaticText(infopage,wxID_ANY,a);
//   tcols->SetFont(sf);
//   infosizer->Add(tcols,wxSizerFlags().Center());

//   return infosizer;
// }

// wxBoxSizer *MuniSide::UpdateInfoImage(wxPanel *infopage, const FitsHdu& hdu,
// 				      const FitsHisto& h)
// {
//   wxFont sf = wxFont(*wxSMALL_FONT);
//   wxString a;

//   wxBoxSizer *infosizer = new wxBoxSizer(wxVERTICAL);

//   //  xPanel *p = new xPanel(infopage);
//   MuniPlotHisto *p = new MuniPlotHisto(infopage,false);
//   //  p->SetAxes(false);
//   //  vector<FitsHisto> histo1;
//   //  histo1.push_back(h);
//   //  p->AssignHistos(histo1);
//   p->SetHisto(h);
//   infosizer->Add(p,wxSizerFlags(4).Expand());

//   const FitsArray array(hdu);
//   //const FitsHdu array(hdu);

//   wxFlexGridSizer *sgrid = new wxFlexGridSizer(3);
 
//   // date & time
//   FitsTime t(array.GetKey(config->fits_dateobs));
//   wxStaticText *tdate = new wxStaticText(infopage,wxID_ANY,t.Date());
//   tdate->SetFont(sf);
//   wxStaticText *ttime = new wxStaticText(infopage,wxID_ANY,t.Time());
//   ttime->SetFont(sf);

//   // exp
//   a = array.Exposure_str(config->fits_exposure) + wxT(" ") 
//     + array.GetUnit(config->fits_exposure);
//   wxStaticText *texp = new wxStaticText(infopage,wxID_ANY,a);
//   texp->SetFont(sf);
  
//   // filter
//   wxStaticText *tfilt = new wxStaticText(infopage,wxID_ANY,
// 					 array.GetKey(config->fits_filter));
//   tfilt->SetFont(sf);

//   sgrid->Add(tdate,wxSizerFlags());
//   sgrid->AddSpacer(sf.GetPointSize());
//   sgrid->Add(tfilt,wxSizerFlags());
//   sgrid->Add(ttime,wxSizerFlags());
//   sgrid->AddSpacer(sf.GetPointSize());
//   sgrid->Add(texp,wxSizerFlags());

//   infosizer->Add(sgrid,wxSizerFlags().Center());

//   // dimensions
//   a.Printf(wxT("%d"),(int) array.Naxes(0));
//   for(int k = 1; k < array.Naxis(); k++) {
//     wxString b;
//     b.Printf(wxT(" × %d"),(int) array.Naxes(k));
//     a += b;
//   }
//   a += wxT(" pixels");

//   wxStaticText *tdim = new wxStaticText(infopage,wxID_ANY,a);
//   tdim->SetFont(sf);
//   infosizer->Add(tdim,wxSizerFlags(1).Center());

//   return infosizer;
// }

// wxBoxSizer *MuniSide::UpdateInfoColor(wxPanel *infopage)
// {
//   wxFont sf = wxFont(*wxSMALL_FONT);
//   wxString a;

//   wxBoxSizer *infosizer = new wxBoxSizer(wxVERTICAL);

//   //  xPanel *p = new xPanel(infopage);
//   MuniPlotHisto *p = new MuniPlotHisto(infopage,false);
//   //  p->SetAxes(false);
//   p->SetHisto(hists);
//   //  p->AssignHistos(hists);
//   infosizer->Add(p,wxSizerFlags(4).Expand());

//   const FitsArray array(fits.Hdu(1));
//   //  const FitsArray array(meta.Hdu(1));
//   const FitsHdu head(array);
//   //  const FitsHdu head(meta.Hdu(1));

//   wxFlexGridSizer *sgrid = new wxFlexGridSizer(3);
 
//   // date & time
//   FitsTime t(head.GetKey(config->fits_dateobs));
//   wxStaticText *tdate = new wxStaticText(infopage,wxID_ANY,t.Date());
//   tdate->SetFont(sf);
//   wxStaticText *ttime = new wxStaticText(infopage,wxID_ANY,t.Time());
//   ttime->SetFont(sf);

//   // exp
//   a = head.Exposure_str(config->fits_exposure) + wxT(" ") 
//     + head.GetUnit(config->fits_exposure);
//   wxStaticText *texp = new wxStaticText(infopage,wxID_ANY,a);
//   texp->SetFont(sf);
  
//   // filter
//   wxStaticText *tfilt = new wxStaticText(infopage,wxID_ANY,wxT("RGB"));
//   tfilt->SetFont(sf);

//   sgrid->Add(tdate,wxSizerFlags());
//   sgrid->AddSpacer(sf.GetPointSize());
//   sgrid->Add(tfilt,wxSizerFlags());
//   sgrid->Add(ttime,wxSizerFlags());
//   sgrid->AddSpacer(sf.GetPointSize());
//   sgrid->Add(texp,wxSizerFlags());

//   infosizer->Add(sgrid,wxSizerFlags().Center());

//   // dimensions
//   a.Printf(wxT("%d"),(int)array.Naxes(0)/*head.Naxes(0)*/);
//   for(int k = 1; k < array.Naxis()/*head.Naxis()*/; k++) {
//     wxString b;
//     b.Printf(wxT(" × %d"),(int)array.Naxes(k)/*head.Naxes(k)*/);
//     a += b;
//   }
//   a += wxT(" pixels");
//   wxStaticText *tdim = new wxStaticText(infopage,wxID_ANY,a);
//   tdim->SetFont(sf);
//   infosizer->Add(tdim,wxSizerFlags(1).Center());

//   return infosizer;
// }


// void MuniSide::Assign(const FitsFile& f, const std::vector<FitsHisto>& hs,
// 		      const wxArrayString& mitems)
// {
//   FitsFile fits = f;
//   hists = hs;
//   treelabels = mitems;
//   wxASSERT(fits.Status() && fits.HduCount() > 0);
  
//   fstruct->DeleteAllItems();

//   // fits structure
//   wxTreeItemId root = fstruct->AddRoot(fits.GetName());

//   if( fits.Type() == FITS_COLOR )
//     fstruct->AppendItem(root,wxT("Color"),3,-1,new xTreeItemData(-1,HDU_COLOR));

//   wxASSERT(fits.HduCount() == mitems.GetCount());
//   for(size_t n = 0; n < fits.HduCount(); n++ ) {
//     const FitsHdu hdu = fits.Hdu(n);

//     wxTreeItemId hdus = fstruct->AppendItem(root,mitems.Item(n),-1,-1,
// 					    new xTreeItemData(n));
//     fstruct->AppendItem(hdus,wxT("Header"),0,-1,new xTreeItemData(n,HDU_HEAD));

//     if( hdu.Type() == HDU_IMAGE || hdu.Type() == HDU_TABLE ) 
//       fstruct->AppendItem(hdus,hdu.Type_str(),hdu.Type() == HDU_IMAGE?1:2,-1,
// 			  new xTreeItemData(n,hdu.Type()));
//     else if( hdu.Type() == HDU_HEAD )
//       ;
//     else
//       fstruct->AppendItem(hdus,hdu.Type_str(),-1,-1,
// 			  new xTreeItemData(n,HDU_UNKNOWN));
      
//   }
//   fstruct->ExpandAll();
// }

void MuniSide::Assign(const FitsMeta& m)
{
  meta = m;
  wxASSERT(meta.IsOk());

  //  treelabels.clear();
  /*  fstruct->*/DeleteAllItems();

  // icons
  wxFont fn(*wxNORMAL_FONT);
  int ssize = 2*(fn.GetPointSize()+3);
  //  int ssize = 2*config->info_symbol.GetWidth();
  wxImageList *icons = new wxImageList(ssize, ssize, true);
  for(size_t n = 0; n < meta.HduCount(); n++ ) {
    const FitsMetaHdu hdu = meta.Hdu(n);
    //    wxImage ico(hdu.GetIcon());
    //    double r = double(ico.GetWidth()) / double(ico.GetHeight());
    //    icons->Add(wxBitmap(ico.Scale(ssize,int(ssize/r))));
    //   icons->Add(wxBitmap(MuniIcon::BrowserIcon(hdu.GetIcon(),ssize,ssize)));
    icons->Add(wxBitmap(MuniIcon::ListIcon(hdu.GetIcon(),ssize)));
  }
  //  wxImage ico(meta.GetIcon());
  //  double r = double(ico.GetWidth()) / double(ico.GetHeight());
  //  icons->Add(wxBitmap(ico.Scale(ssize,int(ssize/r))));
  //  icons->Add(wxBitmap(MuniIcon::BrowserIcon(meta.GetIcon(),ssize,ssize)));
  icons->Add(wxBitmap(MuniIcon::ListIcon(meta.GetIcon(),ssize)));
  /*  fstruct->*/AssignImageList(icons);

  // fits structure
  //  wxTreeItemId root = fstruct->AddRoot(meta.GetName());
  wxTreeItemId root = /*fstruct->*/AddRoot(wxEmptyString);

  if( meta.Type() == FITS_COLOR )
    /*   fstruct->*/AppendItem(root,"Color",meta.HduCount(),-1,new xTreeItemData(-1,HDU_COLOR));

  wxTreeItemId hub = /*fstruct->*/AppendItem(root,meta.GetName());

  //  wxASSERT(fits.HduCount() == mitems.GetCount());
  for(size_t n = 0; n < meta.HduCount(); n++ ) {
    const FitsMetaHdu hdu = meta.Hdu(n);

    /*
    wxString label = hdu.GetKey(wxT("EXTNAME"));
    if( label.IsEmpty() )
      label.Printf(hdu.Type_str()+wxT(" %d"),n);
    */
    //    treelabels.push_back(label);

    wxString label = hdu.GetControlLabel();

    wxTreeItemId hdus = /*fstruct->*/AppendItem(hub,label,n,-1,
					    new xTreeItemData(n));
    /*
    fstruct->AppendItem(hdus,wxT("Header"),0,-1,new xTreeItemData(n,HDU_HEAD));

    if( hdu.Type() == HDU_IMAGE || hdu.Type() == HDU_TABLE ) 
      fstruct->AppendItem(hdus,hdu.Type_str(),hdu.Type() == HDU_IMAGE?1:2,-1,
			  new xTreeItemData(n,hdu.Type()));
    else if( hdu.Type() == HDU_HEAD )
      ;
    else
      fstruct->AppendItem(hdus,hdu.Type_str(),-1,-1,
			  new xTreeItemData(n,HDU_UNKNOWN));
    */
  }
  /*  fstruct->*/Expand(hub);
}

/*
void MuniSide::Assign(const FitsFile& f, const std::vector<FitsHisto>& hs)
{
  fits = f;
  hists = hs;
  wxASSERT( fits.IsOk() );
}
*/

void MuniSide::OnTreeActivated(wxTreeEvent& event) 
{
  wxTreeItemId itemId = event.GetItem();

  xTreeItemData *data = 
    dynamic_cast<xTreeItemData *>(/*fstruct->*/GetItemData(itemId));

  if( ! data ) {

    if( /*fstruct->*/IsExpanded(itemId) )
      /*     fstruct->*/CollapseAllChildren(itemId);
    else
      /*      fstruct->*/ExpandAllChildren(itemId);

  }
  else {

    MuniSideEvent ev(xEVT_SIDE,wxID_ANY);
    ev.hdu = data->index;
    ev.type = data->type;
    //    if( data->index >= 0 )
    //      ev.SetString(treelabels.Item(data->index));

    //wxPostEvent(topwin,ev);
    wxQueueEvent(topwin,ev.Clone());

    //    UpdateInfoPage(data->index);

  }

//   return;

//   xTreeItemData *item = 
//     dynamic_cast<xTreeItemData *>(fstruct->GetItemData(itemId));
//   wxASSERT(item);
//   //  wxLogDebug(_("MuniSide::OnTreeActivated %d %d"),item->index,item->type);

//   int type = 0;
//   if( item )
//     type = item->type;

//   if( type == 0 ) {

//     if( fstruct->IsExpanded(itemId) )
//       fstruct->CollapseAllChildren(itemId);
//     else
//       fstruct->ExpandAllChildren(itemId);

//   }
//   else {
      
//     wxASSERT(item);

//     MuniSideEvent ev(xEVT_SIDE,wxID_ANY);
//     ev.hdu = item->index;
//     ev.type = item->type;
//     if( item->index >= 0 )
//       ev.SetString(treelabels.Item(item->index));

//     wxPostEvent(topwin,ev);
//     UpdateInfoPage(item->index);
//   }
}

void MuniSide::OnHduSelect(MuniSideEvent& event)
{
  HduSelect(event.hdu,/*fstruct->*/GetRootItem());
}

// void MuniSide::HduSelect(int n/*, bool h*/)
// {

//   wxTreeItemId root = fstruct->GetRootItem();
//   HduSelect(n,root);

// }

void MuniSide::HduSelect(int n, const wxTreeItemId& root)
{
  wxTreeItemIdValue cookie;

  wxTreeItemId i = /*fstruct->*/GetFirstChild(root,cookie);
  while( i.IsOk() ) {

    xTreeItemData *d = dynamic_cast<xTreeItemData *>(/*fstruct->*/GetItemData(i));
    if( d && d->index == n ) {

      /*fstruct->*/SelectItem(i);
      return;
      
    }
    else {

      HduSelect(n,i);

//       wxTreeItemId i = fstruct->GetFirstChild(hub,cookie);
//       while( i.IsOk() ) {

// 	xTreeItemData *d = dynamic_cast<xTreeItemData *>(fstruct->GetItemData(i));
// 	if( d && d->index == n ) {
// 	  fstruct->SelectItem(i);
// 	  return;
// 	}
// 	i = fstruct->GetNextChild(hub,cookie);
//       }
    }  
    i = /*fstruct->*/GetNextChild(root,cookie);
  }


//   while( i.IsOk() ) {

//     if( n == -1 ) {
//       fstruct->SelectItem(i);
//       return;
//     }

//     xTreeItemData *k = dynamic_cast<xTreeItemData *>(fstruct->GetItemData(i));

//     //    wxLogDebug(_("MuniSide::HduSelect %d %d %d"),n,k->index,k->type);

//     if( k && k->index == n ){

//       wxTreeItemId x;
//       if( h ) {
// 	wxTreeItemIdValue cook;
// 	x = fstruct->GetFirstChild(i,cook);
//       }
//       else
// 	x = fstruct->GetLastChild(i);
      
//       if( x.IsOk() ) {
// 	fstruct->SelectItem(x);
// 	return;
//       }
//     }
//     i = fstruct->GetNextChild(root,cookie);
//   }
}

// void MuniSide::ShowInfo(bool v)
// {
//   config->sideinfo_show = v;
//   if( config->sideinfo_show && ! IsSplit() ) {
//     wxSize s(infopage->GetBestSize());
//     SplitHorizontally(fstruct,infopage,-s.GetHeight());

//     wxTreeItemId i = fstruct->GetSelection();
//     wxASSERT(i.IsOk());
//     xTreeItemData *d = dynamic_cast<xTreeItemData *>(fstruct->GetItemData(i));
//     //    xTreeItemData *d = static_cast<xTreeItemData>(itemId);
//     wxASSERT(d);
//     UpdateInfoPage(d->index);//item.index);
//   }
//   else
//     Unsplit();
// }



// ---- xPanel


// xPanel::xPanel(wxWindow *w): wxPanel(w,wxID_ANY),hmax(0) 
// {
//   SetBackgroundStyle(wxBG_STYLE_CUSTOM);
// }

// void xPanel::AssignHistos(const std::vector<FitsHisto>& hl)
// {
//   hlist = hl;

//   // search for maximum in histogram
//   hmax = 0;
//   for(std::vector<FitsHisto>::iterator k = hlist.begin(); k != hlist.end();++k){
//     FitsHisto hist = *k;  
//     for(int i = 0; i < hist.NBins(); i++ ) {
//       if( hist.Hist(i) > hmax )
// 	hmax = hist.Hist(i);
//     }
//   }
// }

// void xPanel::OnSize(wxSizeEvent& ev) 
// { 
//   wxSize size = GetClientSize(); 
//   if( size != last ) { 
//     last = size;
//     wxPaintEvent e; 
//     wxPostEvent(this,e);
//   }
// }


// void xPanel::OnPaint(wxPaintEvent& event)
// {
//   wxSize size = GetClientSize();

//   wxImage temp(size.GetWidth(),size.GetHeight(),false);
//   wxBitmap bmp(temp);
//   wxMemoryDC mdc(bmp);

//   // draw
//   wxGraphicsContext *gc = wxGraphicsContext::Create(mdc);
//   if( gc ) {

//     // clear
//     gc->SetBrush(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOW));
//     wxGraphicsPath gp(gc->CreatePath());
//     gp.AddRectangle(0,0,size.GetWidth(),size.GetHeight());
//     gc->FillPath(gp);

//     const unsigned char opacity = hlist.size() > 1 ? 128 : 255;
//     wxColour colour(128,128,128,opacity);
//     int pit = 0;
//     for(std::vector<FitsHisto>::iterator k = hlist.begin(); k!=hlist.end();++k){
//       FitsHisto hist = *k;

//       if( ! hist.IsOk() ) continue;
      
//       double a = hist.Cents(0);
//       double b = double(size.GetWidth())/(hist.Cents(hist.NBins()-1) -
// 					  hist.Cents(0));
//       double hh = double(size.GetHeight())/double(hmax);
      
//       if( hlist.size() > 1 ) {
// 	switch (pit) {
// 	case 0: colour.Set(255,0,0,opacity); break;
// 	case 1: colour.Set(0,255,0,opacity); break;
// 	case 2: colour.Set(0,0,255,opacity); break;
// 	default: colour.Set(64,64,64,opacity); break;
// 	}
//       }
//       gc->SetBrush(gc->CreateBrush(wxBrush(colour)));
//       gc->SetPen(gc->CreatePen(wxPen(colour,0,wxSOLID)));
	
//       for(int l = 0; l < hist.NBins(); l++ ) {
// 	double h = double(hist.Hist(l))*hh;
// 	gc->DrawRectangle((hist.Cents(l)-a)*b,size.GetHeight()-h,1.0,h);
//       }

//       pit++;
//       if( pit == 3 ) pit = 0;
//     }
//     delete gc;
//   }

//   wxAutoBufferedPaintDC dc(this);
//   dc.Blit(0,0,bmp.GetWidth(),bmp.GetHeight(),&mdc,0,0);
// }




