/*

  xmunipack - select astrometry source

  Copyright © 2012-3 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/regex.h>
#include <wx/utils.h>
#include <wx/filefn.h> 
#include <wx/listbook.h>
#include <cfloat>

using namespace std;



MuniSelectSource::MuniSelectSource(wxWindow *w, MuniConfig *c, bool t):
  wxDialog(w,wxID_ANY,"Coordinate Source",wxDefaultPosition,wxDefaultSize,
	   wxDEFAULT_DIALOG_STYLE|wxRESIZE_BORDER),
  config(c), page(0), xframe(t), astrorel(false)
{
  SetIcon(config->munipack_icon);
  CreateControls();  
}

MuniSelectSource::~MuniSelectSource()
{
  EraseTemp();
}

void MuniSelectSource::CreateControls()
{
  wxSizerFlags lf;
  lf.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT);

  wxFont sf(*wxSMALL_FONT);

  const wxSize imageSize(24, 24);
  wxImageList *ilist = new wxImageList(imageSize.GetWidth(), imageSize.GetHeight());
  ilist->Add(wxArtProvider::GetIcon(wxART_GO_DOWN,wxART_TOOLBAR,imageSize));
  ilist->Add(wxArtProvider::GetIcon(wxART_FILE_OPEN,wxART_TOOLBAR,imageSize));

  // title
  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);

  wxListbook *book = new wxListbook(this,wxID_ANY);
  book->AssignImageList(ilist);

  
  wxPanel *cpanel = new wxPanel(book);
  wxBoxSizer *csizer = new wxBoxSizer(wxVERTICAL);

  wxFlexGridSizer *grid = new wxFlexGridSizer(2);
  
  wxStaticText *vlabel = new wxStaticText(cpanel,wxID_ANY,"Virtual Observatory:");
  grid->Add(vlabel,lf);
  
  wxButton *vobutt = new wxButton(cpanel,wxID_ANY,"Search Catalogues");
  grid->Add(vobutt,wxSizerFlags().Border().Left());
  
  wxStaticText *flabel = new wxStaticText(cpanel,wxID_ANY,"File:");
  grid->Add(flabel,lf);
  wxFilePickerCtrl *cpick = new wxFilePickerCtrl(cpanel,wxID_ANY,"","Select A File",
			       "FITS files "+config->dirmask+")|"+
			       config->dirmask+"| All files (*)|*");
  grid->Add(cpick,wxSizerFlags().Expand());
  csizer->Add(grid,wxSizerFlags().Border().Center());

  csizer->Add(new wxStaticLine(cpanel,wxID_ANY),wxSizerFlags().Border().Expand());

  wxFlexGridSizer *s = new wxFlexGridSizer(4);
  s->AddGrowableCol(1);
  s->AddGrowableCol(3);

  s->Add(new wxStaticText(cpanel,wxID_ANY,L"Label α:"),lf);
  choice_ra = new wxChoice(cpanel,ID_ASTRO_CH_RA);
  s->Add(choice_ra,wxSizerFlags().DoubleBorder(wxRIGHT));
  s->Add(new wxStaticText(cpanel,wxID_ANY,L"Label δ:"),lf);
  choice_dec = new wxChoice(cpanel,ID_ASTRO_CH_DEC);
  s->Add(choice_dec);
  s->Add(new wxStaticText(cpanel,wxID_ANY,L"Label μ(α):"),lf);
  choice_pmra = new wxChoice(cpanel,ID_ASTRO_CH_PMRA);
  s->Add(choice_pmra,wxSizerFlags().DoubleBorder(wxRIGHT));
  s->Add(new wxStaticText(cpanel,wxID_ANY,L"Label μ(δ):"),lf);
  choice_pmdec = new wxChoice(cpanel,ID_ASTRO_CH_PMDEC);
  s->Add(choice_pmdec);

  s->Add(new wxStaticText(cpanel,wxID_ANY,"Label mag:"),lf);
  choice_mag = new wxChoice(cpanel,ID_ASTRO_CH_MAG);
  s->Add(choice_mag,wxSizerFlags().Border(wxTOP));

  csizer->Add(s,wxSizerFlags().Border().Center());
  cpanel->SetSizer(csizer);
  book->AddPage(cpanel,"Catalogue",false,0);

  Bind(wxEVT_COMMAND_FILEPICKER_CHANGED,&MuniSelectSource::OnCatFile,this,cpick->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniSelectSource::OnCatVO,this,vobutt->GetId());
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniSelectSource::OnChoice,this,ID_ASTRO_CH_RA);
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniSelectSource::OnChoice,this,ID_ASTRO_CH_DEC);
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniSelectSource::OnChoice,this,ID_ASTRO_CH_PMRA);
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniSelectSource::OnChoice,this,ID_ASTRO_CH_PMDEC);
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniSelectSource::OnChoice,this,ID_ASTRO_CH_MAG);


  if( xframe ) {

    wxPanel *rpanel = new wxPanel(book);
    wxBoxSizer *rsizer = new wxBoxSizer(wxVERTICAL);
    
    wxBoxSizer *ffsizer = new wxBoxSizer(wxHORIZONTAL);
    flabel = new wxStaticText(rpanel,wxID_ANY,"File:");
    ffsizer->Add(flabel,lf);

    wxFilePickerCtrl *rpick = new wxFilePickerCtrl(rpanel,wxID_ANY,"","Select A File",
						   "FITS files "+config->dirmask+")|"+
						   config->dirmask+"| All files (*)|*");
    ffsizer->Add(rpick,wxSizerFlags(1));
    rsizer->Add(ffsizer,wxSizerFlags().Expand().DoubleBorder());
    
    wxCheckBox *check_rel = new wxCheckBox(rpanel,wxID_ANY,"Relative");
    rsizer->Add(check_rel,wxSizerFlags().Border().Center());

    rpanel->SetSizer(rsizer);

    book->AddPage(rpanel,"Frame",false,1);

    Bind(wxEVT_COMMAND_FILEPICKER_CHANGED,&MuniSelectSource::OnRefFile,this,rpick->GetId());
    Bind(wxEVT_COMMAND_CHECKBOX_CLICKED,&MuniSelectSource::OnCheckRel,this,check_rel->GetId());

  }

  topsizer->Add(book,wxSizerFlags(1).Expand().Border());

  wxStdDialogButtonSizer *buttsize = 
    static_cast<wxStdDialogButtonSizer *>(CreateButtonSizer(wxOK|wxCANCEL));
  if( buttsize )
    topsizer->Add(buttsize,wxSizerFlags().Expand().Border());

  SetSizerAndFit(topsizer);

  Bind(wxEVT_UPDATE_UI,&MuniSelectSource::OnUpdateUI,this);
  Bind(wxEVT_COMMAND_LISTBOOK_PAGE_CHANGED,&MuniSelectSource::OnBookChange,this,book->GetId());

}

void MuniSelectSource::EraseTemp()
{
  if( !tmpcatfile.IsEmpty() && wxFileExists(tmpcatfile) ) 
    wxRemoveFile(tmpcatfile);
  tmpcatfile.Clear();
}

void MuniSelectSource::OnUpdateUI(wxUpdateUIEvent& event)
{
  bool def = ! catfile.IsEmpty() || ! reffile.IsEmpty() || ! tmpcatfile.IsEmpty();
      
  wxASSERT(FindWindow(wxID_OK));
  FindWindow(wxID_OK)->Enable(def);
}

wxString MuniSelectSource::GetPath() const
{
  if( page == 0 ) {
    if( !catfile.IsEmpty() )
      return catfile;
    else if( !tmpcatfile.IsEmpty() )
      return tmpcatfile;
  }
  else if( page == 1 )
    return reffile;

  return wxEmptyString;
}

wxString MuniSelectSource::GetId() const
{
  return idlabel;
}

int MuniSelectSource::GetType() const
{
  switch (page) {
  case 0: return ID_ASTRO_CAT;
  case 1: return ID_ASTRO_REF;
  }
  return 0;
}

bool MuniSelectSource::IsTemporary() const
{
  return !tmpcatfile.IsEmpty();
}

bool MuniSelectSource::GetRelative() const
{
  return astrorel;
}

wxString MuniSelectSource::GetLabelRA() const
{
  return label_ra;
}

wxString MuniSelectSource::GetLabelDec() const
{
  return label_dec;
}

wxString MuniSelectSource::GetLabelPMRA() const
{
  return label_pmra;
}

wxString MuniSelectSource::GetLabelPMDec() const
{
  return label_pmdec;
}

wxString MuniSelectSource::GetLabelMag() const
{
  return label_mag;
}

void MuniSelectSource::OnCheckRel(wxCommandEvent& event)
{
  astrorel = event.IsChecked();
}

void MuniSelectSource::OnRefFile(wxFileDirPickerEvent& event)
{
  bool ext = false;

  FitsFile fits(event.GetPath());
  for(size_t i = 0; i < fits.HduCount(); i++) {
    if( fits.Hdu(i).GetExtname() == APEREXTNAME )
      ext = true;
  }
   
  if( ext ) {
    EraseTemp();
    idlabel = fits.GetName();
    Layout();
    Fit();
    catfile.Clear();
    reffile = event.GetPath();
  }
  else {
    wxLogMessage("MUNIPACK extension with detected sources not found.");
    wxASSERT(FindWindow(event.GetId()));
    static_cast<wxFilePickerCtrl *>(FindWindow(event.GetId()))->SetPath("");
    idlabel.Clear();
  }
}

void MuniSelectSource::OnCatFile(wxFileDirPickerEvent& event)
{
  wxArrayString cols;

  if( CheckCatalogue(event.GetPath(),cols) ) {
    
    catfile = event.GetPath();
    reffile.Clear();
    EraseTemp();

    SetLabels(cols);
  }
  Layout();
  Fit();
}

void MuniSelectSource::OnCatVO(wxCommandEvent& event)
{

  MuniCone cone(this,config);
  if( cone.ShowModal() == wxID_OK ) {

    wxArrayString cols;
    if( CheckCatalogue(cone.GetPath(),cols) ) {

      wxString xfile = wxFileName::CreateTempFileName("xmunipack-astrometer-search_");

      if( wxCopyFile(cone.GetPath(),xfile) ) {

	catfile.Clear();
        reffile.Clear();
	EraseTemp();	    

	tmpcatfile = xfile;

	SetLabels(cols);
	return;
      }
    }
  }
  idlabel.Clear();
}

void MuniSelectSource::SetLabels(const wxArrayString& cols)
{
  choice_ra->Set(cols);
  choice_dec->Set(cols);
  choice_pmra->Set(cols);
  choice_pmdec->Set(cols);
  choice_mag->Set(cols);

  VOCatConf catalogs;
  VOCatResources cat(catalogs.GetCat());

  if( cat.IsOk() ) {
  
    label_ra = cat.GetLabel("POS_EQ_RA");
    label_dec = cat.GetLabel("POS_EQ_DEC");
    label_pmra = cat.GetLabel("POS_EQ_PMRA");
    label_pmdec = cat.GetLabel("POS_EQ_PMDEC");
    label_mag = cat.GetSort();

    choice_ra->SetSelection(choice_ra->FindString(label_ra));
    choice_dec->SetSelection(choice_dec->FindString(label_dec));
    choice_pmra->SetSelection(choice_ra->FindString(label_pmra));
    choice_pmdec->SetSelection(choice_dec->FindString(label_pmdec));
    choice_mag->SetSelection(choice_mag->FindString(label_mag));
  }
}


bool MuniSelectSource::CheckCatalogue(const wxString& filename, wxArrayString& cols)
{
  FitsFile fits(filename);

  if( fits.IsOk() ) {
    for(size_t i = 0; i < fits.HduCount(); i++) {
      if( fits.Hdu(i).Type() == HDU_TABLE ) {
	cols = static_cast<FitsTable>(fits.Hdu(i)).GetColLabels();
	idlabel = fits.Hdu(i).GetExtname();
	return true;
      }
    }
  }
  return false;
}

void MuniSelectSource::OnChoice(wxCommandEvent& event)
{
  if( event.GetId() == ID_ASTRO_CH_RA ) 
    label_ra = event.GetString();
  if( event.GetId() == ID_ASTRO_CH_DEC ) 
    label_dec = event.GetString();
  if( event.GetId() == ID_ASTRO_CH_PMRA ) 
    label_pmra = event.GetString();
  if( event.GetId() == ID_ASTRO_CH_PMDEC ) 
    label_pmdec = event.GetString();
  if( event.GetId() == ID_ASTRO_CH_MAG )
    label_mag = event.GetString();
}

void MuniSelectSource::OnBookChange(wxBookCtrlEvent& event)
{
  page = event.GetSelection();
}

