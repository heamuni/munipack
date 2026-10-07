/*

  xmunipack - photometric and astrometric image

  Copyright © 2011-3 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/animate.h>
#include <wx/collpane.h>
#include <cfloat>

using namespace std;


// --- MuniCalibrate


MuniCalibrate::MuniCalibrate(wxWindow *w, MuniConfig *c, const wxString& f):
  wxDialog(w,wxID_ANY,"Calibrate"),config(c),
  throbber(new wxAnimationCtrl(this,wxID_ANY,c->throbber)),
  fwhm(2.0),thresh(5.0),radius(0.01),file(f),pipe(this)
{
  SetIcon(config->munipack_icon);
  EnableCloseButton(false);

  wxFont bf(*wxNORMAL_FONT);
  bf.SetWeight(wxFONTWEIGHT_BOLD);

  wxSpinCtrlDouble *sfwhm = new wxSpinCtrlDouble(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
						 wxDefaultSize,wxSP_ARROW_KEYS,0.0,666.0,fwhm,1.0);
  sfwhm->SetDigits(1);

  wxSpinCtrlDouble *sthresh = new wxSpinCtrlDouble(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
						   wxDefaultSize,wxSP_ARROW_KEYS,0.0,666.0,thresh,1.0);
  sthresh->SetDigits(1);

  wxGridSizer *dtsizer = new wxGridSizer(2);

  wxStaticText *tfwhm = new wxStaticText(this,wxID_ANY,"Fwhm:");
  tfwhm->SetFont(bf);
  dtsizer->Add(tfwhm,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  dtsizer->Add(sfwhm,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT));

  wxStaticText *ttresh = new wxStaticText(this,wxID_ANY,L"Threshold (in σ):");
  ttresh->SetFont(bf);
  dtsizer->Add(ttresh,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  dtsizer->Add(sthresh,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT));

  wxStaticBoxSizer *dsizer = new wxStaticBoxSizer(wxVERTICAL,this,"Detection");
  dsizer->Add(dtsizer,wxSizerFlags().Align(wxALIGN_CENTER));

  wxStaticBoxSizer *asizer = new wxStaticBoxSizer(wxVERTICAL,this,"Astrometry");

  wxGridSizer *gsizer = new wxGridSizer(2);

  wxStaticText *tcat = new wxStaticText(this,wxID_ANY,"Catalog:");
  tcat->SetFont(bf);
  gsizer->Add(tcat,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  wxArrayString catalogs;
  catalogs.Add("UCAC 3");
  wxChoice *scatalog = new wxChoice(this,wxID_ANY,wxDefaultPosition,
				    wxDefaultSize,catalogs);
  scatalog->SetSelection(0);
  gsizer->Add(scatalog,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT));

  wxStaticText *tproj = new wxStaticText(this,wxID_ANY,"Projection:");
  tproj->SetFont(bf);
  gsizer->Add(tproj,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  wxArrayString projections;
  projections.Add("Gnomonic");
  wxChoice *proj = new wxChoice(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
				projections);
  proj->SetSelection(0);
  gsizer->Add(proj,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT));
  asizer->Add(gsizer);

  wxStaticText *talpha = new wxStaticText(this,wxID_ANY,L"α [°]:");
  talpha->SetFont(bf);
  wxStaticText *tdelta = new wxStaticText(this,wxID_ANY,L"δ [°]:");
  tdelta->SetFont(bf);

  wxStaticText *tradius = new wxStaticText(this,wxID_ANY,"Radius [']:");
  tradius->SetFont(bf);

  wxTextCtrl *salpha = new wxTextCtrl(this,wxID_ANY);
  wxTextCtrl *sdelta = new wxTextCtrl(this,wxID_ANY);
  wxTextCtrl *sradius = new wxTextCtrl(this,wxID_ANY);

  wxBoxSizer *coosizer = new wxBoxSizer(wxHORIZONTAL);
  coosizer->Add(talpha,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL));
  coosizer->Add(salpha,wxSizerFlags(1));
  coosizer->Add(tdelta,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL));
  coosizer->Add(sdelta,wxSizerFlags(1));

  asizer->Add(coosizer,wxSizerFlags().Expand());

  wxBoxSizer *rsizer = new wxBoxSizer(wxHORIZONTAL);
  rsizer->Add(tradius,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL));
  rsizer->Add(sradius,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT));
  asizer->Add(rsizer,wxSizerFlags().Left());

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(dsizer,wxSizerFlags().Expand().Border());
  topsizer->Add(asizer,wxSizerFlags().Expand().Border());

  wxCollapsiblePane *collPane = new wxCollapsiblePane(this,wxID_ANY,"Details");
  wxWindow *win = collPane->GetPane();
  wxBoxSizer *lsizer = new wxBoxSizer(wxHORIZONTAL);
  wxTextCtrl *logwin = new wxTextCtrl(win,wxID_ANY);
  lsizer->Add(logwin,wxSizerFlags(1).Expand());
  win->SetSizer(lsizer);
  lsizer->SetSizeHints(win);
  topsizer->Add(collPane,wxSizerFlags().Border(wxLEFT|wxRIGHT|wxBOTTOM));

  wxBoxSizer *bot = new wxBoxSizer(wxHORIZONTAL);
  bot->Add(throbber,wxSizerFlags().Border().Align(wxALIGN_LEFT));
  bot->AddStretchSpacer(1);
  wxSizer *buttons = CreateButtonSizer(wxOK|wxCANCEL);
  if( buttons )
    bot->Add(buttons,wxSizerFlags().Border().Align(wxALIGN_RIGHT));
  topsizer->Add(bot,wxSizerFlags().Expand());
  SetSizerAndFit(topsizer);

  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniCalibrate::OnFwhm,this,sfwhm->GetId());
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniCalibrate::OnThresh,this,sthresh->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniCalibrate::OnApply,this,wxID_OK);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniCalibrate::OnCancel,this,wxID_CANCEL);
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniCalibrate::OnProjection,this,proj->GetId());
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniCalibrate::OnCatalog,this,scatalog->GetId());
  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniCalibrate::OnAlpha,this,salpha->GetId());
  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniCalibrate::OnDelta,this,sdelta->GetId());
  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniCalibrate::OnRadius,this,sradius->GetId());

  throbber->Show(false);
  catalog = "http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=UCAC-3&";
}

wxString MuniCalibrate::GetResult() const { return output; }

void MuniCalibrate::OnFwhm(wxSpinDoubleEvent& event)
{
  fwhm = event.GetValue();
}

void MuniCalibrate::OnThresh(wxSpinDoubleEvent& event)
{
  thresh = event.GetValue();
}

void MuniCalibrate::OnCancel(wxCommandEvent& ev)
{
  wxLogDebug("Leaving Photometry...");

  EndModal(wxID_CANCEL);
}

void MuniCalibrate::OnProjection(wxCommandEvent& ev)
{
  wxLogDebug("Projection..." + ev.GetString());
}

void MuniCalibrate::OnCatalog(wxCommandEvent& ev)
{
  /*
  map<wxString,wxString> catalogs;
  catalogs["Simbad"] = "http://simbad.u-strasbg.fr/simbad/sim-cone?";
  catalogs["UCAC 3"] = "http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=UCAC-3&";
  */

  wxLogDebug("catalog..." + ev.GetString());


  /*
  map<wxString,wxString>::const_iterator i = catalogs.find(ev.GetString());
  wxASSERT(i != catalogs.end());
  catalog = i->second;
  */
}

void MuniCalibrate::OnAlpha(wxCommandEvent& ev)
{
  //  wxLogDebug("catalog..." + ev.GetString());
  wxString a = ev.GetString();
  double x;
  if( a.ToDouble(&x) )
    alpha = x;
}

void MuniCalibrate::OnDelta(wxCommandEvent& ev)
{
  //  wxLogDebug("catalog..." + ev.GetString());
  wxString a = ev.GetString();
  double x;
  if( a.ToDouble(&x) )
    delta = x;
}

void MuniCalibrate::OnRadius(wxCommandEvent& ev)
{
  //  wxLogDebug("catalog..." + ev.GetString());
  wxString a = ev.GetString();
  double x;
  if( a.ToDouble(&x) )
    radius = x/60.0;
}

void MuniCalibrate::OnApply(wxCommandEvent& ev)
{
  wxLogDebug("Running OnCalibrate...");

  output = wxFileName::CreateTempFileName("xmunipack-photometry");
  wxRemoveFile(output);
  //  FitsCopyFile(file,output);
  //  wxCopyFile(file,output);
  wxCopyFile(file,file+"~");


  Bind(wxEVT_END_PROCESS,&MuniCalibrate::OnFinish,this);

  //  MuniProcess *action = new MuniProcess(&pipe,"cat");
  MuniProcess *action = new MuniProcess(&pipe,"photometry");
  pipe.push(action);

  wxString a;

  a.Printf("FWHM = %lf",fwhm);
  action->Write(a);

  a.Printf("THRESHOLD = %lf",thresh);
  action->Write(a);

  action->Write("OUTPUT = '" + output + "'");
  action->Write("FILE = '" + file + "'");

  MuniProcess *cone = new MuniProcess(&pipe,"cone");
  pipe.push(cone);

  coutput = wxFileName::CreateTempFileName("xmunipack-catalogue");
  wxRemoveFile(coutput);

  cone->Write("SORT = 'Va'"); // !!!! just only for UCAC-3 !!!!
  cone->Write("TYPE = 'FITS'");
  cone->Write("OUTPUT = '" + coutput + "'");
  cone->Write("URL = '" + catalog + "'");

  a.Printf("SR = '%lf'",radius);
  cone->Write(a);
  a.Printf("RA = '%lf'",alpha);
  cone->Write(a);
  a.Printf("DEC = '%lf'",delta);
  cone->Write(a);

  MuniProcess *astrometry = new MuniProcess(&pipe,"astrofit");
  pipe.push(astrometry);

  FitsFile fits(file);
  FitsArray array(fits.Hdu(0));
  double crpix1 = array.Naxes(0)/2.0;
  double crpix2 = array.Naxes(1)/2.0;

  wxLogDebug("%lf %lf",crpix1,crpix2);

  a.Printf("CRPIX1 = %lf",crpix1);
  astrometry->Write(a);
  a.Printf("CRPIX2 = %lf",crpix2);
  astrometry->Write(a);

  astrometry->Write("OUTPUT = '" + file + "'");
  astrometry->Write("CAT = '" + coutput + "'");
  astrometry->Write("FILE = '" + output + "'");


  throbber->Show(true);
  throbber->Play();

  //  wxWindow *bok = FindWindowById(wxID_APPLY,this);
  //  bok->Enable(false);
  FindWindowById(wxID_OK,this)->Enable(false);

  Layout();

  pipe.Start();
}

void MuniCalibrate::OnFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniCalibrate::OnFinish");

  Unbind(wxEVT_END_PROCESS,&MuniCalibrate::OnFinish,this);

  throbber->Stop();
  throbber->Show(false);

  //  wxWindow *bok = FindWindowById(wxID_APPLY,this);
  //  bok->Enable(true);
  FindWindowById(wxID_OK,this)->Enable(true);

  Layout();




  if( event.GetExitCode() != 0 )
    wxLogDebug("Failed with exit code %d",event.GetExitCode());
  else {

    FitsFile fits(output);

    if( fits.IsOk() )
      for(size_t i = 0; i < fits.HduCount(); i++)
	if( fits.Hdu(i).GetKey("EXTNAME").Find("MUNIPACK") != wxNOT_FOUND ) {
	  FitsTable t(fits.Hdu(i));
	  dynamic_cast<MuniDisplay *>(GetParent())->SetStars(t);
	}

    EndModal(wxID_OK);

  }
}
