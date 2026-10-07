/*

  xmunipack - find stars

  Copyright © 2019-22 F.Hroch (hroch@physics.muni.cz)

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
#include "plot.h"
#include "help.h"
#include <wx/wx.h>
#include <wx/dialog.h>
#include <wx/regex.h>
#include <wx/tokenzr.h>
#include <wx/filefn.h>
#include <vector>
#include <cfloat>

using namespace std;

MuniFind::MuniFind(wxWindow *w, MuniConfig *c, const wxString& s,
		   const FitsHdu& h):
  wxDialog(w,wxID_ANY,"Find stars",wxDefaultPosition,wxDefaultSize,
	   wxDEFAULT_DIALOG_STYLE|wxRESIZE_BORDER), config(c),
  mproc(0), timer(this), anim(0), fitsname(s), array(FitsArray(h)),
  fwhm(config->find_fwhm), thresh(config->find_thresh)
{
  SetIcon(config->munipack_icon);
  satur = InitSatur();

  wxSizerFlags label_flag, entry_flag;
  label_flag.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).Border(wxRIGHT);
  entry_flag.Align(wxALIGN_LEFT);

  wxSpinCtrlDouble *spin_fwhm =
    new wxSpinCtrlDouble(this,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
			 wxSP_ARROW_KEYS,0.01,666.0,fwhm,0.5);
  spin_fwhm->SetDigits(1);
  spin_fwhm->SetToolTip("An effective spread of stars in pixels.\n"
			"The spread is measured as the full width\n"
			"at half of maximum (FWHM).");

  wxSpinCtrlDouble *spin_thresh =
    new wxSpinCtrlDouble(this,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
			 wxSP_ARROW_KEYS,1.0,666.0,thresh,1.0);
  spin_thresh->SetDigits(1);
  spin_thresh->SetToolTip("Selects stars with maximum over the threshold.\n"
			  "Threshold is to relative to averadged standard\n"
			  "deviation of background.");

  wxSpinCtrlDouble *spin_satur =
    new wxSpinCtrlDouble(this,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
			 wxSP_ARROW_KEYS,0.0,1e10,satur,65e3);
  spin_satur->SetDigits(2);

  wxFlexGridSizer *fsizer = new wxFlexGridSizer(2);

  fsizer->Add(new wxStaticText(this,wxID_ANY,"FWHM"),label_flag);
  fsizer->Add(spin_fwhm,entry_flag);

  wxStaticBoxSizer *psizer =
    new wxStaticBoxSizer(wxVERTICAL,this," Spread ");
  psizer->Add(fsizer,wxSizerFlags().Centre().Border());

  plot = new MuniPlotFind(this,array);
  static_cast<MuniPlotFind *>(plot)->SetFwhm(fwhm);
  psizer->Add(plot,wxSizerFlags(1).Expand().Border());

  wxFlexGridSizer *fxsizer = new wxFlexGridSizer(2);

  fxsizer->Add(new wxStaticText(this,wxID_ANY,"Threshold"),label_flag);
  fxsizer->Add(spin_thresh,entry_flag);

  fxsizer->Add(new wxStaticText(this,wxID_ANY,"Saturation"),label_flag);
  fxsizer->Add(spin_satur,entry_flag);

  wxStaticBoxSizer *lsizer = new wxStaticBoxSizer(wxVERTICAL,this," Limits ");
  lsizer->Add(fxsizer,wxSizerFlags().Centre().Border());

  start = new wxButton(this,wxID_ANY,"Find");
  stop = new wxButton(this,wxID_ANY,"Stop");
  stop->Show(false);
  stop->SetLabelMarkup("<span fgcolor=\"#ff0000\"><b>Stop</b></span>");
  status = new wxStaticText(this,wxID_ANY,"");

  wxBoxSizer *statsizer = new wxBoxSizer(wxHORIZONTAL);
  wxSizerFlags button_flag;
  button_flag.Align(wxALIGN_CENTRE_VERTICAL|wxALIGN_RIGHT).Border(wxRIGHT);
  statsizer->Add(status,wxSizerFlags(1).DoubleBorder().Align(wxALIGN_CENTRE_VERTICAL).ReserveSpaceEvenIfHidden());
  statsizer->Add(start,button_flag);
  statsizer->Add(stop,button_flag);

  wxBoxSizer *footsizer = new wxBoxSizer(wxHORIZONTAL);
  anim = new wxAnimationCtrl(this,wxID_ANY,config->throbber,
			     wxDefaultPosition,wxSize(32,32));
  anim->Show(false);
  footsizer->Add(anim,wxSizerFlags().Align(wxALIGN_CENTRE_VERTICAL).Border(wxLEFT).ReserveSpaceEvenIfHidden());
  footsizer->AddStretchSpacer(1);
  wxStdDialogButtonSizer *bsizer = static_cast<wxStdDialogButtonSizer *>
    (CreateButtonSizer(wxCANCEL|wxAPPLY|wxHELP));
  if( bsizer )
    footsizer->Add(bsizer,wxSizerFlags().Right());

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(psizer,wxSizerFlags(1).Expand().Border());
  topsizer->Add(lsizer,wxSizerFlags().Expand().Border());
  topsizer->Add(statsizer,wxSizerFlags().Border().Expand());
  topsizer->Add(footsizer,wxSizerFlags().Expand().Border(wxTOP|wxBOTTOM));

  SetSizerAndFit(topsizer);

  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniFind::OnHelp,this,wxID_HELP);
  Bind(wxEVT_BUTTON,&MuniFind::OnFind,this,start->GetId());
  Bind(wxEVT_BUTTON,&MuniFind::OnStop,this,stop->GetId());
  Bind(wxEVT_SPINCTRLDOUBLE,&MuniFind::OnFwhm,this,spin_fwhm->GetId());
  Bind(wxEVT_SPINCTRLDOUBLE,&MuniFind::OnThresh,this,spin_thresh->GetId());
  Bind(wxEVT_SPINCTRLDOUBLE,&MuniFind::OnSatur,this,spin_satur->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniFind::OnStdButton,this,wxID_APPLY);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniFind::OnStdButton,this,wxID_CANCEL);

  SetAffirmativeId(wxID_APPLY);
  SetEscapeId(wxID_CANCEL);

  FindWindow(wxID_APPLY)->Enable(false);
}

MuniFind::~MuniFind()
{
  config->find_fwhm = fwhm;
  config->find_thresh = thresh;

  if( wxFileExists(tmpfits) )
    wxRemoveFile(tmpfits);
}

void MuniFind::OnStdButton(wxCommandEvent& event)
{
  if( event.GetId() == wxID_APPLY ) {
    if( ! FitsCopyHdu(tmpfits,fitsname,"FIND") )
      wxLogError("Failed to save results to " + fitsname +
		 " (space, access, ... ?).");
  }

  // stop a potentially working subprocess
  if( mproc )
    mproc->Kill();

  // clear detected sources
  CleanDraw();

  // notify parent
  wxCommandEvent e(EVT_FINISH_DIALOG,this->GetId());
  e.SetString(fitsname);
  e.SetInt(event.GetId());
  wxQueueEvent(GetParent(),e.Clone());

  // good-bye
  SetReturnCode(event.GetId());
  event.Skip();
}

void MuniFind::OnHelp(wxCommandEvent& event)
{
  MuniHelp(config->munipack_html_dir,"man_find.html");
}

double MuniFind::InitSatur() const
{
  //  if( saturate keyword

  if( array.Bitpix() > 0 ) {
    int n = 1;
    for(int i = 0; i < array.Bitpix(); i++)
      n = 2*n;
    return n - 1;
  }
  else {
    FitsArrayStat s(array);
    return s.GetMax()*(1.0 - FLT_EPSILON);
  }
}

void MuniFind::OnFwhm(wxSpinDoubleEvent& event)
{
  fwhm = event.GetValue();
  static_cast<MuniPlotFind *>(plot)->SetFwhm(fwhm);
  static_cast<MuniPlotFind *>(plot)->Update();
}

void MuniFind::SetPoint(int i0, int j0) const
{
  static_cast<MuniPlotFind *>(plot)->SetPoint(i0,j0);
  static_cast<MuniPlotFind *>(plot)->Update();
}

void MuniFind::OnThresh(wxSpinDoubleEvent& event)
{
  thresh = event.GetValue();
}

void MuniFind::OnSatur(wxSpinDoubleEvent& event)
{
  satur = event.GetValue();
}

void MuniFind::OnStop(wxCommandEvent& event)
{
  wxASSERT(mproc);
  mproc->Kill();
}

void MuniFind::OnFind(wxCommandEvent& event)
{
  wxLogDebug("MuniFind::OnFind");

  // setup GUI elements
  start->Show(false);
  stop->Show(true);
  anim->Show(true);
  anim->Play();
  CleanDraw();

  // work on the copy
  status->SetLabel("Making a copy of data...");
  Layout();
  if( ! tmpfits.IsEmpty() )
    wxRemoveFile(tmpfits);
  tmpfits = wxFileName::CreateTempFileName("xmunipack-find_");
  wxCopyFile(fitsname,tmpfits,true);

  status->SetLabelMarkup("<small>Starting find...</small>");
  Layout();

  index = 0;
  timer.Start(250);
  mproc = new MuniProcess(this,"find");
  mproc->SetEcho(false);

  mproc->Write("PIPELOG = T");
  mproc->Write("FWHM = %e",fwhm);
  mproc->Write("SATURATE = %e",satur);
  mproc->Write("THRESHOLD = %e",thresh);
  mproc->Write("NFILES = 1");
  mproc->Write("FILE = '"+tmpfits+"' '"+tmpfits+"'");

  mproc->OnStart();

  Bind(wxEVT_END_PROCESS,&MuniFind::OnFindFinish,this);
  Bind(wxEVT_TIMER,&MuniFind::OnTimer,this);

  FindWindow(wxID_APPLY)->Enable(false);

  Layout();
}

void MuniFind::OnTimer(wxTimerEvent& event)
{
  wxArrayString out(mproc->GetOutput());
  if( ! (out.GetCount() - index > 0) )
    return;

  //  wxLogDebug("OnTimer %d %d",int(index),int(out.GetCount()));
  DrawStars(index,out.GetCount(),out);
  index = out.GetCount();

  if( out.Last().IsEmpty() )
    return;

  long n = LastStar(out.Last());
  if( n > 0 ) {
    wxString row;
    row.Printf("<small>Finished %ld rows.</small>",n);
    status->SetLabelMarkup(row);
  }
  else
    status->SetLabel("?");

  Layout();
}

void MuniFind::OnFindFinish(wxProcessEvent& event)
{
  wxASSERT(mproc);
  wxLogDebug("MuniFind::OnFindFinish %d",mproc->GetExitCode());

  if( mproc->GetExitCode() == 0 ) {
    wxArrayString out(mproc->GetOutput());
    DrawStars(index,out.GetCount(),out);

    if( out.GetCount() == 0 )
      status->SetLabelMarkup("<small>No stars found.</small>");
    else {
      long n = LastStar(out.Last());
      if( n > 0 ) {
	wxString row;
	row.Printf("<small>Found %ld stars.</small>",n);
	status->SetLabelMarkup(row);
      }
    }

    FindWindow(wxID_APPLY)->Enable(true);
  }
  else {
    status->SetLabelMarkup("<span fgcolor=\"#ff0000\"><b>Search failed.</b></span>");
    wxLogError("Command: find");
    wxArrayString err(mproc->GetErrors());
    wxArrayString out(mproc->GetOutput());
    wxArrayString in(mproc->GetInput());
    for(size_t i = 0; i < in.GetCount(); i++)
      wxLogError("Input:"+in[i]);
    for(size_t i = 0; i < out.GetCount(); i++)
      wxLogError("Output:"+out[i]);
    for(size_t i = 0; i < err.GetCount(); i++)
      wxLogError("Error:"+err[i]);
    wxLogError("Search failed.");
  }

  delete mproc;
  mproc = 0;

  timer.Stop();
  anim->Stop();
  anim->Show(false);
  start->Show(true);
  stop->Show(false);
  Layout();

  Unbind(wxEVT_END_PROCESS,&MuniFind::OnFindFinish,this);
  Unbind(wxEVT_TIMER,&MuniFind::OnTimer,this);
}

wxString MuniFind::Parser(const wxString& out) const
{
  wxRegEx re("^=(.*)> (.+)");
  wxASSERT(re.IsValid());

  if( re.Matches(out) ) {

    wxString key(re.GetMatch(out,1));
    wxString value(re.GetMatch(out,2));

    if( key == "FIND" )
      return value;
  }

  return "";
}

bool MuniFind::StarParser(const wxString& line, long *n, double *x, double *y) const
{
  if( line.IsEmpty() )
    return false;

  wxStringTokenizer tk(line);
  int m = 0;
  while( tk.HasMoreTokens() ) {
    wxString a = tk.GetNextToken();
    m++;
    if( m == 1 && ! a.ToLong(n) )
      return false;
    if( m == 2 && ! a.ToDouble(x) )
      return false;
    if( m == 3 && ! a.ToDouble(y) )
      return false;
  }
  return true;
}

long MuniFind::LastStar(const wxString& line) const
{
  long n;
  double x,y;
  if( StarParser(Parser(line),&n,&x,&y) )
    return n;

  return -1;
}


void MuniFind::DrawStars(long n, long m, const wxArrayString& out) const
{
  const double csize = wxMax(fwhm,1.0);
  vector<wxObject *> objects;

  wxColour b(255,255,0,120);
  wxPen pen(b,2);
  objects.push_back(new MuniDrawPen(pen));

  for(long i = n; i < m; i++) {
    long n;
    double x,y;
    if( StarParser(Parser(out[i]),&n,&x,&y) )
      objects.push_back(new MuniDrawCross(x,y,csize));
  }

  if( ! objects.empty() ) {
    MuniLayer layer(ID_FIND,objects);
    MuniDrawEvent ev(EVT_DRAW);
    ev.layer = layer;
    wxQueueEvent(GetParent(),ev.Clone());
  }
}

void MuniFind::CleanDraw() const
{
  MuniLayer layer(ID_FIND);
  MuniDrawEvent ev(EVT_DRAW);
  ev.layer = layer;
  wxQueueEvent(GetParent(),ev.Clone());
}
