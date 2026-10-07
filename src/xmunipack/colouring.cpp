/*

  xmunipack - colouring


  Copyright © 2010-2013, 2018 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/statline.h>
#include <wx/filesys.h>
#include <vector>

using namespace std;



// -- colourization target ---

class ColourDropTarget: public wxDropTarget
{
public:
  ColourDropTarget(wxWindow *);
  wxDragResult OnData(wxCoord, wxCoord, wxDragResult);

private:
  wxWindow *target;
};

ColourDropTarget::ColourDropTarget(wxWindow *w): target(w)
{
  SetDataObject(new MuniDataObjectMeta);
}

wxDragResult ColourDropTarget::OnData(wxCoord x, wxCoord y, wxDragResult def)
{
  if ( !GetData() )
    return wxDragNone;

  wxDragResult result = wxDragNone;

  MuniDataObjectMeta *dobj = (MuniDataObjectMeta *) GetDataObject();
  wxASSERT(dobj);
  vector<FitsMeta> slist = dobj->GetMetafitses();
  if( slist.size() == 1 && slist[0].Type() == FITS_GREY ) {
    if( static_cast<MuniColoring *>(target)->SetDropMeta(x,y,slist) )
      result = wxDragCopy;
  }

  return result;
}


// - Coloring ------------------------------------------------


MuniColoring::MuniColoring(wxWindow *w, MuniConfig *c):
  wxDialog(w,wxID_ANY,"Colouring"),iSize(48),config(c),
  filename("colouring.fits"),pipe(this),
  throbber(new wxAnimationCtrl(this,wxID_ANY,c->throbber))
{
  SetIcon(config->munipack_icon);
  EnableCloseButton(false);
  SetDropTarget(new ColourDropTarget(this));

  Init();
  CreateControls();
}

void MuniColoring::Init()
{
  //  cchoices.Add("RGB");              opt.Add("RGB");
  cchoices.Add("CIE 1931 XYZ");     opt.Add("XYZ");
  cchoices.Add("Johnson BVR");      opt.Add("Johnson BVR");

  colorspace = opt[0];

  index = -1;
}

void MuniColoring::CreateControls()
{
  throbber->Show(false);

  wxFont bf(*wxNORMAL_FONT);
  bf.SetWeight(wxFONTWEIGHT_BOLD);

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);

  wxBoxSizer *ysizer = new wxBoxSizer(wxVERTICAL);

  // top icon

  wxBoxSizer *osizer = new wxBoxSizer(wxHORIZONTAL);
  wxImage coloring_icon = MuniConfig::LoadImage("google-gadgets-gtk.png");
  MuniThumbCanvas *th = new MuniThumbCanvas(this,coloring_icon);
  osizer->Add(th,wxSizerFlags().Center().DoubleBorder());

  wxStaticText *title = new wxStaticText(this,wxID_ANY,"Colour picture");
  title->SetFont(bf);
  osizer->Add(title,wxSizerFlags().DoubleBorder().Align(wxALIGN_CENTER_VERTICAL));
  ysizer->Add(osizer,wxSizerFlags().Center());

  // colourspace
  wxBoxSizer *ctrls = new wxBoxSizer(wxHORIZONTAL);
  ctrls->Add(new wxStaticText(this,wxID_ANY,"Colourspace:"),
	  wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).Border());

  cspace = new wxChoice(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
			 cchoices);
  cspace->SetSelection(0);
  cspace->SetToolTip("Select input colospace");
  ctrls->Add(cspace,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).Border());

  ysizer->Add(ctrls,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).Border());

  // imagelist
  list = new wxListView(this,wxID_ANY,wxDefaultPosition,wxSize(-1,200),
			wxLC_REPORT|wxLC_NO_HEADER|wxLC_SINGLE_SEL);
  list->SetToolTip("Drag fits files here.");
  ysizer->Add(list,wxSizerFlags().Expand().DoubleBorder(wxLEFT|wxRIGHT));

  // image file load
  wxBoxSizer *fsizer = new wxBoxSizer(wxHORIZONTAL);
  fsizer->Add(new wxStaticText(this,wxID_ANY,"Load file:"),
	  wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).Border());

  wxFilePickerCtrl *fpic =
    new wxFilePickerCtrl(this,wxID_ANY,wxEmptyString,"Choose a file",
			 "FITS files ("+config->dirmask+")|"+config->dirmask+
			 "|All files (*)|*");
  fsizer->Add(fpic,wxSizerFlags(1).Expand().Border());
  ysizer->Add(fsizer,wxSizerFlags().Right().Border());

  // separator
  ysizer->Add(new wxStaticLine(this,wxID_ANY),wxSizerFlags().Expand().DoubleBorder());

  // output filename
  wxFlexGridSizer *gsizer = new wxFlexGridSizer(2);
  gsizer->AddGrowableCol(1);

  gsizer->Add(new wxStaticText(this,wxID_ANY,"Filename:"),
	      wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL));

  wxTextCtrl *fname = new wxTextCtrl(this,wxID_ANY,filename);
  gsizer->Add(fname,wxSizerFlags(1).Expand().Border());

  gsizer->Add(new wxStaticText(this,wxID_ANY,"Location:"),
	      wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL));


  wxDirPickerCtrl *dpic = new wxDirPickerCtrl(this,wxID_ANY);
  gsizer->Add(dpic,wxSizerFlags(1).Expand().Border(wxLEFT|wxRIGHT));
  ysizer->Add(gsizer,wxSizerFlags().Expand().Border());


  topsizer->Add(ysizer,wxSizerFlags().Border());

  wxBoxSizer *bot = new wxBoxSizer(wxHORIZONTAL);
  bot->Add(throbber,wxSizerFlags().Border().Align(wxALIGN_LEFT));
  bot->AddStretchSpacer(1);
  wxSizer *buttons = CreateButtonSizer(wxOK|wxCANCEL);
  if( buttons )
    bot->Add(buttons,wxSizerFlags().Border().Align(wxALIGN_RIGHT));
  topsizer->Add(bot,wxSizerFlags().Expand());
  SetSizerAndFit(topsizer);

  FindWindowById(wxID_OK,this)->Enable(false);

  // actions

  Bind(wxEVT_COMMAND_FILEPICKER_CHANGED,&MuniColoring::OnBandfile,this,fpic->GetId());
  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniColoring::OnFilename,this,fname->GetId());
  Bind(wxEVT_COMMAND_DIRPICKER_CHANGED,&MuniColoring::OnDirname,this,dpic->GetId());
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniColoring::OnColorspace,this,cspace->GetId());
  Bind(wxEVT_COMMAND_LIST_ITEM_ACTIVATED,&MuniColoring::OnListSelected,this,list->GetId());
  Bind(wxEVT_COMMAND_LIST_ITEM_SELECTED,&MuniColoring::OnListSelected,this,list->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniColoring::OnApply,this,wxID_OK);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniColoring::OnCancel,this,wxID_CANCEL);
  Bind(wxEVT_UPDATE_UI,&MuniColoring::OnUpdateBandfile,this,fpic->GetId());
  Bind(wxEVT_UPDATE_UI,&MuniColoring::OnUpdateOk,this,wxID_OK);

  InitList(cspace->GetStringSelection());
}

bool MuniColoring::SetDropMeta(int x, int y, const std::vector<FitsMeta>& l)
{
  wxASSERT(list);

  if( l.size() != 1 ) return false;

  wxPoint p = list->GetPosition();

  for(int i = 0; i < list->GetItemCount(); i++) {
    wxRect r;
    list->GetItemRect(i,r);
    if( r.Contains(x-p.x,y-p.y) ) {
      SetMeta(i,l[0]);
      return true;
    }
  }

  return false;
}

void MuniColoring::OnUpdateOk(wxUpdateUIEvent& event)
{

  bool ok = metalist.size() > 0;
  for(size_t i = 0; i < metalist.size(); i++)
    if( ! metalist[i].IsOk() ){
      ok = false;
      return;
    }

  wxFileName f(dirname,filename);
  event.Enable(ok && f.IsOk());
}

void MuniColoring::OnUpdateBandfile(wxUpdateUIEvent& event)
{
  event.Enable(list->GetSelectedItemCount() == 1);
}

void MuniColoring::OnBandfile(wxFileDirPickerEvent& event)
{
  wxASSERT(index > -1);

  FitsFile fits(event.GetPath());
  if( fits.IsOk() && fits.Type() == FITS_GREY ) {
    /*
    MuniIcon micon(fits,config->display_colourspace,config->cdatafile,
		   config->default_icon,config->table_icon,config->head_icon,
		   config->icon_size);
    FitsMeta meta(fits,micon.GetIcon(),micon.GetList());
    FitsArrayStat array(fits.Hdu(0));
    */
    //    SetMeta(index,meta,array.Med(),1.0);
  }
}

void MuniColoring::OnFilename(wxCommandEvent& event)
{
  filename = event.GetString();
}

void MuniColoring::OnDirname(wxFileDirPickerEvent& event)
{
  dirname = event.GetPath();
}

void MuniColoring::OnApply(wxCommandEvent& ev)
{
  Bind(wxEVT_END_PROCESS,&MuniColoring::OnFinish,this);

  MuniProcess *action = new MuniProcess(&pipe,"colouring");
  //MuniProcess *action = new MuniProcess(&pipe,"cat");
  pipe.push(action);

  wxFileName f(dirname,filename);

  action->Write("CTABLE = '" + config->cdatafile + "'");
  action->Write("PHSYSTABLE = '"+config->phsystemfile+"'");
  action->Write("OUTPUT = '" + f.GetFullPath() + "'");

  if( ! colorspace.IsEmpty() )
    action->Write("COLOURSPACE = '" + colorspace + "'");

  action->Write("NBANDS = %d",int(param_lines.size()));

  for(size_t i = 0; i < param_lines.size(); i++)
    action->Write("BAND = '"+param_lines[param_lines.size()-i-1]+"'");

  pipe.Start();

  FindWindowById(wxID_OK,this)->Enable(false);

  throbber->Show(true);
  throbber->Play();
  Layout();
}

void MuniColoring::OnCancel(wxCommandEvent& ev)
{
  Unbind(wxEVT_END_PROCESS,&MuniColoring::OnFinish,this);

  wxCommandEvent event(EVT_FILELOAD,GetId());
  wxQueueEvent(GetParent(),event.Clone());
}

void MuniColoring::OnFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniColoring::OnFinish");

  Unbind(wxEVT_END_PROCESS,&MuniColoring::OnFinish,this);

  throbber->Stop();
  throbber->Show(false);

  FindWindowById(wxID_OK,this)->Enable(true);

  Layout();

  if( event.GetExitCode() == 0 ) {

    wxFileName f(dirname,filename);

    wxCommandEvent e(EVT_FILELOAD,GetId());
    e.SetString(f.GetFullPath());
    wxQueueEvent(GetParent(),e.Clone());
  }
  else
    wxLogError("Failed creation of a color image.");
}

void MuniColoring::OnListSelected(wxListEvent& ev)
{
  index = ev.GetIndex();
}

void MuniColoring::OnColorspace(wxCommandEvent& ev)
{
  colorspace = opt[ev.GetSelection()];
  InitList(ev.GetString());
}

void MuniColoring::InitList(const wxString& cspace)
{
  wxASSERT(list);

  wxArrayString labels;
  vector<wxColour> colors;

  // if( cspace == cchoices[0] ) {

  //   // RGB
  //   labels.Add("R");
  //   labels.Add("G");
  //   labels.Add("B");

  //   colors.push_back(wxColour(*wxRED));
  //   colors.push_back(wxColour("FOREST GREEN"));
  //   colors.push_back(wxColour(*wxBLUE));
  // }
  // else
  if( cspace == cchoices[0] ) {

    // XYZ
    labels.Add("X");
    labels.Add("Y");
    labels.Add("Z");

    colors.push_back(wxColour(*wxRED));
    colors.push_back(wxColour("FOREST GREEN"));
    colors.push_back(wxColour(*wxBLUE));

  }
  else if( cspace == cchoices[1] ) {

    // BVR
    labels.Add("Johnson R");
    labels.Add("Johnson V");
    labels.Add("Johnson B");

    colors.push_back(wxColour(*wxRED));
    colors.push_back(wxColour("FOREST GREEN"));
    colors.push_back(wxColour(*wxBLUE));

  }

  wxASSERT(colors.size() == labels.GetCount() && ! colors.empty());


  // clear lists
  list->ClearAll();
  metalist.clear();
  param_lines.clear();

  // fill list
  wxSize is(iSize,iSize);
  icons = new wxImageList(iSize, iSize, true);
  for(size_t i = 0; i < colors.size(); i++)
    icons->Add(MuniIcon::BulletIcon(is,colors[i]));
  list->AssignImageList(icons,wxIMAGE_LIST_SMALL);

  // fill list
  list->InsertColumn(0,wxEmptyString);
  for(size_t i = 0; i < labels.GetCount(); i++)
    list->InsertItem(i,labels[i],i);

  // set item width
  wxSize s = list->GetClientSize();
  list->SetColumnWidth(0,s.GetWidth());

  // fill auxliary arrays
  for(size_t i = 0; i < colors.size(); i++) {
    metalist.push_back(FitsMeta());
    param_lines.push_back(wxEmptyString);
  }
}

void MuniColoring::SetMeta(int idx, const FitsMeta& meta,
			   double black, double weight)
{
  unsigned char r = 255, g = 255, b = 255;

  if( idx == 0 ) {
    r = 255; g = 96; b = 96;
  }
  else if( idx == 1 ) {
    r = 96; g = 192; b = 96;
  }
  else if( idx == 2 ) {
    r = 96; g = 96; b = 255;
  }

  wxColour colour(r,g,b);
  wxImage img(MuniIcon::ListIcon(meta.GetIcon(),iSize,colour));


  icons->Replace(idx,wxBitmap(img));
  list->SetItemText(idx,meta.GetName());

  wxSize s = list->GetClientSize();
  list->SetColumnWidth(0,s.GetWidth());


  // statistics
  if( weight < 0.0 ) {
    FitsFile fits(meta.GetURL());
    if( fits.IsOk() && fits.Type() == FITS_GREY) {
      FitsArrayStat array(fits.Hdu(0));
      //      black = array.Med();
      weight = 1.0;
    }
  }

  wxString filter = meta.GetFilter(config->fits_key_filter);
  wxFileName name(wxFileSystem::URLToFileName(meta.GetURL()));
  //  wxString a;
  //  a.Printf("'" + name.GetFullPath() + "' '" + filter + "'");

  param_lines[idx] = name.GetFullPath();
  metalist[idx] = meta;
}
