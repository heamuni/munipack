/*

  xmunipack - dark, flat, .. corrections


  Copyright © 2010-2013 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/filename.h>
#include <wx/valgen.h>
#include <wx/msgdlg.h>
#include <wx/filesys.h>
#include <vector>


using namespace std;

enum {
  ID_RESULT_DEFAULT, ID_RESULT_OVERWRITE, ID_RESULT_DIR
};



// --- Options

class MuniCorrectOptions: public wxDialog
{
public:
  MuniCorrectOptions(wxWindow *, MuniConfig *, const wxString&);
  wxString GetBitpix() const;
  wxString GetMask() const;
  wxString GetSuffix() const;

private:

  MuniConfig *config;
  wxRadioButton *bitpix0, *bitpix1;
  wxTextCtrl *xmask;

  bool bitpix_16bit,bitpix_float;
  wxString suff;

  void Init(const wxString&);
  void CreateControls();

};


MuniCorrectOptions::MuniCorrectOptions(wxWindow *w, MuniConfig *c,
				       const wxString& suf):
  wxDialog(w,wxID_ANY,"Corrections Options"),config(c)
{
  SetIcon(config->munipack_icon);
  EnableCloseButton(false);

  Init(suf);
  CreateControls();
}

void MuniCorrectOptions::Init(const wxString& s)
{
  bitpix_16bit = false;
  bitpix_float = true;
  suff = s;
}

void MuniCorrectOptions::CreateControls()
{
  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);

  wxFlexGridSizer *gsizer = new wxFlexGridSizer(2);
  gsizer->AddGrowableCol(1);

  // output file
  gsizer->Add(new wxStaticText(this,wxID_ANY,"Output FITS:"),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));

  wxBoxSizer *psizer = new wxBoxSizer(wxVERTICAL);

  bitpix0 = new wxRadioButton(this,wxID_ANY,"16-bit",wxDefaultPosition,
			    wxDefaultSize,wxRB_GROUP);
  bitpix1 = new wxRadioButton(this,wxID_ANY,"float");
  bitpix0->SetToolTip("Select representation of numbers in output image. The float numbers are intended for general usage. 16-bit numbers are less-precise but saves some space.");
  bitpix1->SetToolTip("Select representation of numbers in output image. The float numbers are intended for general usage. 16-bit numbers are less-precise but saves some space.");
  psizer->Add(bitpix0,wxSizerFlags());
  psizer->Add(bitpix1,wxSizerFlags());
  gsizer->Add(psizer,wxSizerFlags().Border());

  xmask = new wxTextCtrl(this,wxID_ANY);
  xmask->SetToolTip("Sets a suffix for output filenames.");
  gsizer->Add(new wxStaticText(this,wxID_ANY,"Suffix:"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  gsizer->Add(xmask,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).Border());
  topsizer->Add(gsizer,wxSizerFlags().Border().Expand());

  topsizer->Add(CreateButtonSizer(wxOK|wxCANCEL),wxSizerFlags().Expand().Border());
  SetSizerAndFit(topsizer);

  // data
  bitpix0->SetValidator(wxGenericValidator(&bitpix_16bit));
  bitpix1->SetValidator(wxGenericValidator(&bitpix_float));
  xmask->SetValidator(wxTextValidator(wxFILTER_NONE,&suff));
}


wxString MuniCorrectOptions::GetBitpix() const
{
  if( bitpix_16bit )
    return "16";

  if( bitpix_float )
    return "-32";

  return wxEmptyString;
}


wxString MuniCorrectOptions::GetMask() const
{
  return suff;
}

wxString MuniCorrectOptions::GetSuffix() const
{
  return suff;
}


// -- MuniDarkbat

MuniDarkbat::MuniDarkbat(wxWindow *w, wxWindowID id, long t, MuniConfig *c):
  MuniListWindow(w,id,t,c), config(c), suffix("_final"),
  mode(ID_RESULT_DEFAULT),pipe(this), timer(this)
{
  wxFont bold(wxSystemSettings::GetFont(wxSYS_DEFAULT_GUI_FONT));
  bold.SetWeight(wxFONTWEIGHT_BOLD);

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);

  wxBoxSizer *hsizer = new wxBoxSizer(wxHORIZONTAL);

  wxWindow *window = new wxPanel(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
				 wxTAB_TRAVERSAL|wxBORDER_THEME);

  wxBoxSizer *sl = new wxBoxSizer(wxVERTICAL);

  wxFlexGridSizer *gsizer = new wxFlexGridSizer(3);
  gsizer->AddGrowableCol(1);

  gsizer->Add(new wxStaticText(window,wxID_ANY,"Bias:"),
  	      wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  bpic = new wxFilePickerCtrl(window,wxID_ANY,wxEmptyString,
			      "Select Bias",config->dirmask);
  gsizer->Add(bpic,wxSizerFlags(1).Expand());

  wxButton *bb = new wxButton(window,wxID_ANY,L"↻",wxDefaultPosition,wxDefaultSize,
			      wxBU_EXACTFIT);
  bb->SetToolTip("Clear");
  gsizer->Add(bb,wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  gsizer->Add(new wxStaticText(window,wxID_ANY,"Dark:"),
  	      wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  dpic = new wxFilePickerCtrl(window,wxID_ANY,wxEmptyString,
			      "Select Dark",config->dirmask);
  gsizer->Add(dpic,wxSizerFlags(1).Expand().Border(wxTOP|wxBOTTOM));

  wxButton *bd = new wxButton(window,wxID_ANY,L"↻",wxDefaultPosition,wxDefaultSize,
			      wxBU_EXACTFIT);
  bd->SetToolTip("Clear");
  gsizer->Add(bd,wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  gsizer->Add(new wxStaticText(window,wxID_ANY,"Flat:"),
	      wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  fpic = new wxFilePickerCtrl(window,wxID_ANY,wxEmptyString,
			      "Select Flat-field",config->dirmask);
  gsizer->Add(fpic,wxSizerFlags(1).Expand());

  wxButton *bf = new wxButton(window,wxID_ANY,L"↻",wxDefaultPosition,wxDefaultSize,
			      wxBU_EXACTFIT);
  bf->SetToolTip("Clear");
  gsizer->Add(bf,wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  sl->Add(gsizer,wxSizerFlags().Expand().Align(wxALIGN_CENTER_VERTICAL));
  hsizer->Add(sl,wxSizerFlags(1).Border());

  hsizer->Add(new wxStaticLine(window,wxID_ANY,wxDefaultPosition,wxDefaultSize,wxLI_VERTICAL),wxSizerFlags().Expand().DoubleBorder());

  wxBoxSizer *sr = new wxBoxSizer(wxVERTICAL);

  r0 = new wxRadioButton(window,ID_RESULT_DEFAULT,"Modify filenames",
			 wxDefaultPosition,wxDefaultSize,wxRB_GROUP);
  r0->SetToolTip("Adds a suffix to the original filename (x.fits to x_suffix.fits). See Options... to setup the suffix.");
  r1 = new wxRadioButton(window,ID_RESULT_OVERWRITE,"In place");
  r1->SetToolTip("Overwrite original files. BE VERY CAREFUL !");
  r2 = new wxRadioButton(window,ID_RESULT_DIR,"Destination:");
  r2->SetToolTip("Selects another directory to store resuls. Keeps original filenames.");
  dirpic = new wxDirPickerCtrl(window,wxID_ANY);

  wxFlexGridSizer *rsizer = new wxFlexGridSizer(2);
  rsizer->AddGrowableCol(1);
  int space = bold.GetPointSize();
  rsizer->Add(new wxStaticText(window,wxID_ANY,"Results:"/*"Filename:"*/),
	      wxSizerFlags().Border(wxRIGHT));
  rsizer->Add(r0,wxSizerFlags());
  rsizer->AddSpacer(space);
  rsizer->Add(r1,wxSizerFlags());
  rsizer->AddSpacer(space);
  rsizer->Add(r2,wxSizerFlags());

  rsizer->AddSpacer(space);
  wxBoxSizer *xs = new wxBoxSizer(wxHORIZONTAL);
  xs->AddSpacer(space);
  xs->Add(dirpic,wxSizerFlags());
  rsizer->Add(xs,wxSizerFlags());
  sr->Add(rsizer,wxSizerFlags(1).Border());
  hsizer->Add(sr,wxSizerFlags().Border());
  topsizer->Add(hsizer,wxSizerFlags().Expand());

  bcre = new wxButton(window,ID_DARKBAT_CREATE,"Create");
  wxButton *bdet = new wxButton(window,ID_OPTIONS,"Options...");
  wxSizer *butt = new wxBoxSizer(wxHORIZONTAL);
  gauge = new wxGauge(window,wxID_ANY,1);
  label = new wxStaticText(window,wxID_ANY,"Label");

  butt->Add(bcre,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL));
  butt->Add(bdet,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL));
  butt->Add(gauge,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL));
  butt->Add(label,wxSizerFlags(1).Border().Align(wxALIGN_CENTER_VERTICAL));
  topsizer->Add(butt,wxSizerFlags().Border().Expand());

  window->SetSizer(topsizer);
  AddWindow(window);

  dirpic->Enable(false);
  gauge->Show(false);
  label->Show(false);

  Bind(wxEVT_UPDATE_UI,&MuniDarkbat::OnUpdateButt,this,bcre->GetId());
  Bind(wxEVT_UPDATE_UI,&MuniDarkbat::OnUpdateDirpic,this,dirpic->GetId());
  Bind(wxEVT_COMMAND_DIRPICKER_CHANGED,&MuniDarkbat::OnDirname,this,dirpic->GetId());
  Bind(wxEVT_COMMAND_FILEPICKER_CHANGED,&MuniDarkbat::OnFlatname,this,fpic->GetId());
  Bind(wxEVT_COMMAND_FILEPICKER_CHANGED,&MuniDarkbat::OnDarkname,this,dpic->GetId());
  Bind(wxEVT_COMMAND_FILEPICKER_CHANGED,&MuniDarkbat::OnBiasname,this,bpic->GetId());
  Bind(wxEVT_COMMAND_RADIOBUTTON_SELECTED,&MuniDarkbat::OnResult,this,r1->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniDarkbat::OnClearBias,this,bb->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniDarkbat::OnClearDark,this,bd->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniDarkbat::OnClearFlat,this,bf->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniDarkbat::OnOptions,this,bdet->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniDarkbat::OnCreate,this,bcre->GetId());
}

void MuniDarkbat::OnCreate(wxCommandEvent& event)
{
  results.Clear();

  wxString bpix;

  if( ! bitpix.IsEmpty() ) {
    if( bitpix == "float" )
      bpix = "BITPIX = -32";
    else if( bitpix == "16-bit" )
      bpix = "BITPIX = 16";
    else
      bpix = "BITPIX = -32";
  }
  else
    bpix = "BITPIX = -32";


  wxString cmd,cline;

  if( ! bfilename.IsEmpty() ) {
    cmd = "dark";
    cline = "Darkframe = '"+bfilename+"'";
  }

  if( ! dfilename.IsEmpty() ) {
    cmd = "dark";
    cline = "Darkframe = '"+dfilename+"'";
  }

  if( ! ffilename.IsEmpty() ) {
    cmd = "flat";
    cline = "Flatframe = '"+ffilename+"'";
  }

  MuniProcess *c = new MuniProcess(&pipe,cmd);
  pipe.push(c);
  c->Write(cline);
  c->Write(bpix);

  vector<FitsMeta> metalist = GetAllMeta();
  vector<FitsMeta>::const_iterator m;
  for(m = metalist.begin(); m != metalist.end(); ++m){
    wxASSERT(m->IsOk());
    wxFileName fname(wxFileSystem::URLToFileName(m->GetURL()));
    wxString orig = fname.GetFullPath();
    c->Write("FILE = '"+orig+"'");
    results.Add(orig);
  }


  /*

    The processing can be sterted only when there is at least
    one MuniProcess. The condition is checked by OnUpdateBut.

  */

  timer.Start(500);
  gauge->SetRange(metalist.size());
  gauge->Show(true);
  label->Show(true);
  label->SetLabel("Correcting ...");

  Bind(wxEVT_END_PROCESS,&MuniDarkbat::OnFinish,this);
  Bind(wxEVT_TIMER,&MuniDarkbat::OnUpdate,this);

  pipe.Start();

  Layout();
}

void MuniDarkbat::OnFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniDarkbat::OnFinish");

  timer.Stop();
  gauge->Show(false);
  label->Show(false);

  Unbind(wxEVT_END_PROCESS,&MuniDarkbat::OnFinish,this);
  Unbind(wxEVT_TIMER,&MuniDarkbat::OnUpdate,this);

  if( event.GetExitCode() == 0 ) {

    if( mode == ID_RESULT_OVERWRITE )
      DeleteAllMeta();

    DeSelectAll();
    AddFits(results);
    Layout();
  }
}

void MuniDarkbat::OnUpdate(wxTimerEvent& event)
{
  wxLogDebug("MuniDarkbat::OnUpdate ");
  wxArrayString out(pipe.GetOutput());

  int n = out.GetCount();
  if( n <= 0 || n > gauge->GetRange() )
    gauge->Pulse();
  else
    gauge->SetValue(n);
}

wxString MuniDarkbat::CreateResult(const wxString& orig) const
{
  if( mode == ID_RESULT_DEFAULT ) {
    wxFileName name(orig);
    wxFileName result(name.GetPath(),name.GetName()+suffix,name.GetExt());
    return result.GetFullPath();
  }
  else if( mode == ID_RESULT_OVERWRITE ) {
    return orig;
  }
  else if( mode == ID_RESULT_DIR ) {
    wxFileName name(orig);
    wxFileName result(dirname,name.GetFullName());
    return result.GetFullPath();
  }
  wxFAIL_MSG("REACHED !!! ");
  return wxEmptyString;
}

void MuniDarkbat::OnFlatname(wxFileDirPickerEvent& event)
{
  ffilename = event.GetPath();
}

void MuniDarkbat::OnDarkname(wxFileDirPickerEvent& event)
{
  dfilename = event.GetPath();
}

void MuniDarkbat::OnBiasname(wxFileDirPickerEvent& event)
{
  bfilename = event.GetPath();
}

void MuniDarkbat::OnClearBias(wxCommandEvent& event)
{
  bpic->SetPath(wxEmptyString);
  bfilename = wxEmptyString;
}

void MuniDarkbat::OnClearDark(wxCommandEvent& event)
{
  dpic->SetPath(wxEmptyString);
  dfilename = wxEmptyString;
}

void MuniDarkbat::OnClearFlat(wxCommandEvent& event)
{
  fpic->SetPath(wxEmptyString);
  ffilename = wxEmptyString;
}

void MuniDarkbat::OnDirname(wxFileDirPickerEvent& event)
{
  dirname = event.GetPath();
}

void MuniDarkbat::OnResult(wxCommandEvent& event)
{
  if( event.GetId() == ID_RESULT_OVERWRITE ) {
    wxMessageDialog dialog(this,
			   "Are you sure you want permanently replace all originals ?",
			   "Really Overwrite ?",wxOK|wxCENTRE|wxICON_EXCLAMATION);
    dialog.SetExtendedMessage("Select another choice to leave originals untouched.");
    dialog.ShowModal();
  }
}

void MuniDarkbat::OnOptions(wxCommandEvent& event)
{
  MuniCorrectOptions opt(this,config,suffix);
  if( opt.ShowModal() == wxID_OK ) {
    bitpix = opt.GetBitpix();
    suffix = opt.GetSuffix();
  }
}


void MuniDarkbat::OnUpdateDirpic(wxUpdateUIEvent& event)
{
  event.Enable(r2->GetValue());
}


void MuniDarkbat::OnUpdateButt(wxUpdateUIEvent& event)
{
  bool c = ! bfilename.IsEmpty() || ! dfilename.IsEmpty() || ! ffilename.IsEmpty();

  event.Enable(GetItemCount() > 0 && c);

  if( r0->GetValue() )
    mode = ID_RESULT_DEFAULT;
  else if( r1->GetValue() )
    mode = ID_RESULT_OVERWRITE;
  else if( r2->GetValue() )
    mode = ID_RESULT_DIR;
}
