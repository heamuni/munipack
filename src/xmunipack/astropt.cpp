/*

  xmunipack - astrometry options

  Copyright © 2013-5 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/collpane.h>
#include <wx/statline.h>
#include <wx/regex.h>
#include <wx/utils.h>
#include <wx/sstream.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>
#include <wx/filefn.h>
#include <wx/richtooltip.h>
#include <list>
#include <cfloat>


#define ARCUAS L"μas"
#define ARCMAS "mas"
#define ARCSEC "arcsec"
#define ARCMIN "arcmin"
#define ARCDEG "deg"

using namespace std;


MuniAstrometryOptions::MuniAstrometryOptions(wxWindow *w, MuniConfig *c):
  wxPanel(w),
  config(c),
  sig(config->astrometry_sig),
  fsig(config->astrometry_fsig),
  minmatch(config->astrometry_minmatch),
  maxmatch(config->astrometry_maxmatch),
  matchtype(ID_ASTRO_MATCH),
  full_match(config->astrometry_fullmatch),
  output_units(config->astrometry_units)
{
  wxSizerFlags sf, rl;
  sf.Align(wxALIGN_CENTER_VERTICAL);
  rl.Align(wxALIGN_CENTER_VERTICAL).Border(wxRIGHT);

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);

  wxBoxSizer *ts = new wxBoxSizer(wxHORIZONTAL);

  ts->Add(new wxStaticText(this,wxID_ANY,"Matching method:"),
		rl.DoubleBorder(wxRIGHT));

  wxRadioButton *mmet = new wxRadioButton(this,ID_ASTRO_MATCH,"Back-tracking",
					  wxDefaultPosition,wxDefaultSize,wxRB_GROUP);
  wxRadioButton *nmet = new wxRadioButton(this,ID_ASTRO_NEAR,"Nearly");
  mmet->SetToolTip("Back-tracking (limited by some heuristics) explores all permutations of stars to get mutual correlation between star lists.");
  nmet->SetToolTip("Star lists are correlated by looking for nearly stars. One needs an valuable initial estimate of parameters.");
  ts->Add(mmet,rl);
  ts->Add(nmet,rl);
  topsizer->Add(ts,wxSizerFlags().Center().Border());

  wxFlexGridSizer *gsizer = new wxFlexGridSizer(4);
  gsizer->AddGrowableCol(0);
  gsizer->AddGrowableCol(2);

  gsizer->Add(new wxStaticText(this,wxID_ANY,"Minimal Sequence:"),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));
  wxSpinCtrl *snmatch = new wxSpinCtrl(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
				       wxDefaultSize,wxSP_ARROW_KEYS,3,100,minmatch);
  gsizer->Add(snmatch,sf);

  gsizer->Add(new wxStaticText(this,wxID_ANY,"Position Errors:"),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));

  wxSpinCtrlDouble *ssig =
    new wxSpinCtrlDouble(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
			 wxDefaultSize,wxSP_ARROW_KEYS,0.1,100.0,sig,1.0);
  ssig->SetToolTip("Position errors of centroids of objects (about half or less of FWHM) in pixels.");
  gsizer->Add(ssig,sf);

  gsizer->Add(new wxStaticText(this,wxID_ANY,"Maximal Sequence:"),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));
  wxSpinCtrl *smmatch = new wxSpinCtrl(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
				       wxDefaultSize,wxSP_ARROW_KEYS,5,137,maxmatch);
  gsizer->Add(smmatch,sf);

  snmatch->SetToolTip("Minimal length of a matching sequence. Values under 5 leads to an uncertain matching, over 20 may be slow.");
  smmatch->SetToolTip("Maximum length of a matching sequence. Crowded fields requires higher values over 30.");

  gsizer->Add(new wxStaticText(this,wxID_ANY,"Flux Errors:"),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));

  wxSpinCtrlDouble *fxsig =
    new wxSpinCtrlDouble(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
			 wxDefaultSize,wxSP_ARROW_KEYS,0.1,100.0,fsig,1.0);
  fxsig->SetToolTip("Relative errors of fluxes of objects.");
  gsizer->Add(fxsig,sf);

  gsizer->Add(new wxStaticText(this,wxID_ANY,""),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));

  wxCheckBox *full = new wxCheckBox(this,wxID_ANY,"Full Match");
  full->SetToolTip("Scan throughout all stars. Normally, one finish when the first successful match has occurred.");
  full->SetValue(full_match);
  gsizer->Add(full,wxSizerFlags().Border(wxTOP|wxBOTTOM).Align(wxALIGN_LEFT));

  gsizer->Add(new wxStaticText(this,wxID_ANY,"Output Units:"),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));
  wxArrayString xunits;
  xunits.Add(ARCUAS);
  xunits.Add(ARCMAS);
  xunits.Add(ARCSEC);
  xunits.Add(ARCMIN);
  xunits.Add(ARCDEG);
  wxChoice *ounits = new wxChoice(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,xunits);
  ounits->SetToolTip("Units used for storing both the list of residuals and the scale in FITS headers.");
  gsizer->Add(ounits,sf);
  ounits->SetSelection(ounits->FindString(output_units));
  topsizer->Add(gsizer,wxSizerFlags().Center().Border());

  SetSizer(topsizer);

  Bind(wxEVT_UPDATE_UI,&MuniAstrometryOptions::OnUpdateUI,this);
  Bind(wxEVT_COMMAND_CHECKBOX_CLICKED,&MuniAstrometryOptions::OnFullMatch,
       this,full->GetId());
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniAstrometryOptions::OnChoiceUnits,
       this,ounits->GetId());
  Bind(wxEVT_COMMAND_SPINCTRL_UPDATED,&MuniAstrometryOptions::OnSpinMinMatch,
       this,snmatch->GetId());
  Bind(wxEVT_COMMAND_SPINCTRL_UPDATED,&MuniAstrometryOptions::OnSpinMaxMatch,
       this,smmatch->GetId());
  Bind(wxEVT_COMMAND_RADIOBUTTON_SELECTED,&MuniAstrometryOptions::OnMatchType,
       this,ID_ASTRO_MATCH);
  Bind(wxEVT_COMMAND_RADIOBUTTON_SELECTED,&MuniAstrometryOptions::OnMatchType,
       this,ID_ASTRO_NEAR);
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniAstrometryOptions::OnSpinSig,
       this,ssig->GetId());
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniAstrometryOptions::OnSpinFSig,
       this,fxsig->GetId());

  ids.push_back(snmatch->GetId());
  ids.push_back(smmatch->GetId());
  ids.push_back(full->GetId());
}

MuniAstrometryOptions::~MuniAstrometryOptions()
{
  config->astrometry_sig = sig;
  config->astrometry_fsig = fsig;
  config->astrometry_units = output_units;
  config->astrometry_minmatch = minmatch;
  config->astrometry_maxmatch = maxmatch;
  config->astrometry_fullmatch = full_match;
}

void MuniAstrometryOptions::OnChoiceUnits(wxCommandEvent& event)
{
  output_units = event.GetString();
}

void MuniAstrometryOptions::OnFullMatch(wxCommandEvent& event)
{
  full_match = event.IsChecked();
}

void MuniAstrometryOptions::OnSpinMinMatch(wxSpinEvent& event)
{
  minmatch = event.GetPosition();
}

void MuniAstrometryOptions::OnSpinMaxMatch(wxSpinEvent& event)
{
  maxmatch = event.GetPosition();
}

void MuniAstrometryOptions::OnMatchType(wxCommandEvent& event)
{
  matchtype = event.GetId();
}

void MuniAstrometryOptions::OnSpinSig(wxSpinDoubleEvent& event)
{
  sig = event.GetValue();
}

void MuniAstrometryOptions::OnSpinFSig(wxSpinDoubleEvent& event)
{
  fsig = event.GetValue();
}

void MuniAstrometryOptions::OnUpdateUI(wxUpdateUIEvent& e)
{
  bool enable = matchtype == ID_ASTRO_MATCH;
  for(list<int>::const_iterator i = ids.begin(); i != ids.end(); ++i)
    FindWindow(*i)->Enable(enable);
}

int MuniAstrometryOptions::GetMatchType() const { return matchtype; }
int MuniAstrometryOptions::GetMinMatch() const { return minmatch; }
int MuniAstrometryOptions::GetMaxMatch() const { return maxmatch; }
double MuniAstrometryOptions::GetSig() const { return sig; }
double MuniAstrometryOptions::GetFSig() const { return fsig; }
bool MuniAstrometryOptions::GetFullMatch() const { return full_match; }
wxString MuniAstrometryOptions::GetOutputUnits() const { return output_units; }
