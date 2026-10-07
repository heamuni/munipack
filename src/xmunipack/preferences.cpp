/*

  xmunipack - preferences


  Copyright © 2010-2013, 2018-22 F.Hroch (hroch@physics.muni.cz)

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


  Notes.

  The derivation from wxPropertySheetDialog has stranges:
     - default constructor must not be called
     - validation and style style set work just only when
       Create is called after SetExtraStyle and SetSheetStyle
     - Style shrinkof works horribly for (almost) emtpy tab

*/

#include "xmunipack.h"
#include "fits.h"
#include <wx/wx.h>
#include <wx/valgen.h>

MuniPreferences::MuniPreferences(wxWindow *w, MuniConfig *c): config(c)
{
  SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);

  // SetSheetStyle(wxPROPSHEET_TOOLBOOK/*|wxPROPSHEET_SHRINKTOFIT*/);
  Create(w,wxID_ANY,wxEmptyString);
  CreateButtons();

  wxBookCtrlBase* book = GetBookCtrl();

  // View
  wxPanel *vpanel = new wxPanel(book);
  FeedView(vpanel);
  book->AddPage(vpanel,"View",true);

  // FITS keywords
  wxPanel *kpanel = new wxPanel(book);
  FeedKeywords(kpanel);
  book->AddPage(kpanel,"FITS",false);

  LayoutDialog();

  Bind(wxEVT_IDLE,&MuniPreferences::OnIdle,this);
  Bind(wxEVT_CLOSE_WINDOW,&MuniPreferences::OnClose,this,GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniPreferences::OnOk,this,wxID_OK);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniPreferences::OnCancel,this,wxID_CANCEL);

  SetAffirmativeId(wxID_OK);
  SetEscapeId(wxID_CANCEL);
}

void MuniPreferences::OnClose(wxCloseEvent& event)
{
  wxQueueEvent(GetParent(),event.Clone());
}

void MuniPreferences::OnOk(wxCommandEvent& event)
{
  Validate() && TransferDataFromWindow();
  Close();
}

void MuniPreferences::OnCancel(wxCommandEvent& event)
{
  Close();
}

void MuniPreferences::FeedView(wxPanel *panel)
{
  // types of coordinates
  wxArrayString cootypes = FitsCoo::Label_str();
  wxArrayString labels;
  for(size_t i = 1; i < cootypes.GetCount(); i++)
    labels.Add(cootypes[i]);
  wxChoice *ch = new wxChoice(panel,wxID_ANY,wxDefaultPosition,wxDefaultSize,labels);
  ch->SetSelection(config->display_cootype - (COO_FIRST+2));

  // types of photometric quantities
  wxArrayString qphtypes = FitsValue::Label_str();
  qphtypes.RemoveAt(0);
  wxChoice *chph =
    new wxChoice(panel,wxID_ANY,wxDefaultPosition,wxDefaultSize,qphtypes);
  chph->SetSelection(config->display_qphtype - (PHQUANTITY_FIRST+2));

  wxSizerFlags label_flags;
  label_flags.Align(wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT);
  wxSizerFlags choice_flags;
  choice_flags.Align(wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT|wxTOP);

  wxBoxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);
  sizer->Add(new wxStaticText(panel,wxID_ANY,"Format:"),label_flags);
  sizer->Add(ch,choice_flags);

  wxBoxSizer *pizzer = new wxBoxSizer(wxHORIZONTAL);
  pizzer->Add(new wxStaticText(panel,wxID_ANY,"Quantity:"),label_flags);
  pizzer->Add(chph,choice_flags);

  wxStaticBoxSizer *cs = new wxStaticBoxSizer(wxVERTICAL,panel,"Coordinates");
  cs->Add(sizer,wxSizerFlags().Border());

  wxStaticBoxSizer *ps = new wxStaticBoxSizer(wxVERTICAL,panel,"Surface photometry");
  ps->Add(pizzer,wxSizerFlags().Border());

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(cs,wxSizerFlags().Border().Expand());
  topsizer->Add(ps,wxSizerFlags().Border().Expand());

  panel->SetSizer(topsizer);

  Bind(wxEVT_CHOICE,&MuniPreferences::OnCooType,this,ch->GetId());
  Bind(wxEVT_CHOICE,&MuniPreferences::OnQPhType,this,chph->GetId());
}


void MuniPreferences::FeedKeywords(wxPanel *panel)
{
  wxFont bf(*wxNORMAL_FONT);
  bf.SetWeight(wxFONTWEIGHT_BOLD);

  wxSizerFlags label_sizer;
  label_sizer.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT);

  wxSizerFlags control_sizer(1);
  control_sizer.Align(wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT).Expand();

  /*
  wxStaticText *label = new wxStaticText(panel,wxID_ANY,"FITS Keywords:");
  label->SetFont(bf);
  */

  wxFlexGridSizer *sgrid = new wxFlexGridSizer(2);
  sgrid->AddGrowableCol(1);

  sgrid->Add(new wxStaticText(panel,wxID_ANY,"Object:"),
	     label_sizer);
  wxTextCtrl *tobject = new wxTextCtrl(panel,wxID_ANY);
  sgrid->Add(tobject,control_sizer);

  sgrid->Add(new wxStaticText(panel,wxID_ANY,"Observer:"),
	     label_sizer);
  wxTextCtrl *tobser = new wxTextCtrl(panel,wxID_ANY);
  sgrid->Add(tobser,control_sizer);

  sgrid->Add(new wxStaticText(panel,wxID_ANY,"Date-obs:"),
	     label_sizer);
  wxTextCtrl *tdate = new wxTextCtrl(panel,wxID_ANY);
  sgrid->Add(tdate,control_sizer);

  sgrid->Add(new wxStaticText(panel,wxID_ANY,"Exposure:"),
	     label_sizer);
  wxTextCtrl *texp = new wxTextCtrl(panel,wxID_ANY);
  sgrid->Add(texp,control_sizer);

  sgrid->Add(new wxStaticText(panel,wxID_ANY,"Filter:"),
	     label_sizer);
  wxTextCtrl *tfilter = new wxTextCtrl(panel,wxID_ANY);
  sgrid->Add(tfilter,control_sizer);

  sgrid->Add(new wxStaticText(panel,wxID_ANY,"Area:"),
	     label_sizer);
  wxTextCtrl *tarea = new wxTextCtrl(panel,wxID_ANY);
  sgrid->Add(tarea,control_sizer);

  sgrid->Add(new wxStaticText(panel,wxID_ANY,"Gain:"),
	     label_sizer);
  wxTextCtrl *tgain = new wxTextCtrl(panel,wxID_ANY);
  sgrid->Add(tgain,control_sizer);

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  //  topsizer->Add(label,wxSizerFlags().Align(wxALIGN_CENTER).Border());
  topsizer->Add(sgrid,wxSizerFlags().Border().Expand());
  panel->SetSizer(topsizer);
  topsizer->Fit(panel);

  tobject->SetValidator(wxTextValidator(wxFILTER_NONE,&config->fits_key_object));
  tobser->SetValidator(wxTextValidator(wxFILTER_NONE,&config->fits_key_observer));
  tdate->SetValidator(wxTextValidator(wxFILTER_NONE,&config->fits_key_dateobs));
  texp->SetValidator(wxTextValidator(wxFILTER_NONE,&config->fits_key_exptime));
  tfilter->SetValidator(wxTextValidator(wxFILTER_NONE,&config->fits_key_filter));
  tgain->SetValidator(wxTextValidator(wxFILTER_NONE,&config->fits_key_gain));
  tarea->SetValidator(wxTextValidator(wxFILTER_NONE,&config->fits_key_area));

}


void MuniPreferences::OnCooType(wxCommandEvent& event)
{
  config->display_cootype = event.GetSelection() + COO_FIRST+2;
  wxQueueEvent(this,new MuniConfigEvent(EVT_CONFIG_UPDATED,ID_PREFERENCES_COOTYPE));
}

void MuniPreferences::OnQPhType(wxCommandEvent& event)
{
  config->display_qphtype = event.GetSelection() + PHQUANTITY_FIRST+2;
  wxQueueEvent(this,new MuniConfigEvent(EVT_CONFIG_UPDATED,ID_PREFERENCES_QPHTYPE));
}


void MuniPreferences::OnIdle(wxIdleEvent& event)
{
  int n = GetBookCtrl()->GetSelection();
  if( n >= 0 )
    SetTitle(GetBookCtrl()->GetPageText(n));
}

void MuniPreferences::SelectPage(int n)
{
  GetBookCtrl()->SetSelection(n);
}
