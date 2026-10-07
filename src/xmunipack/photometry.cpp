/*

  xmunipack - photometry

  Copyright © 2011-2 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/regex.h>
#include <wx/sstream.h>
#include <wx/txtstrm.h>
#include <wx/utils.h>
#include <wx/richtooltip.h>
#include <cfloat>
#include <vector>
#include <list>

#define BUTTLABEL "Start  Processing"
#define MAXGRANGE 100

using namespace std;

// --- BaseDialog

MuniBaseDialog::MuniBaseDialog(wxWindow *w, const wxString& t, MuniConfig *c, const wxString& f):
  wxDialog(w,wxID_ANY,t),config(c),pipe(this),timer(this),
  running(false),readonly(false),exitcode(0),file(f)
{
  SetIcon(config->munipack_icon);

  readonly = ! wxFileName::IsFileWritable(file);

  workingfile = wxFileName::CreateTempFileName("xmunipack-"+t+"_");
  wxRemoveFile(workingfile);
  FitsCopyFile(file,workingfile);
};

MuniBaseDialog::~MuniBaseDialog()
{
  if( ! workingfile.IsEmpty() && wxFileExists(workingfile) )
    wxRemoveFile(workingfile);
}

bool MuniBaseDialog::RunProcess(MuniProcess *p)
{
  wxASSERT(pipe.empty());

  pipe.push(p);

  Bind(wxEVT_END_PROCESS,&MuniBaseDialog::OnFinish,this);
  Bind(wxEVT_TIMER,&MuniBaseDialog::OnTimer,this);

  wxBeginBusyCursor();

  running = true;
  timer.Start(100);
  pipe.Start();

  return true; //?????
}

void MuniBaseDialog::OnTimer(wxTimerEvent& event) { event.Skip(); }

void MuniBaseDialog::OnFinish(wxProcessEvent& event)
{
  event.Skip();

  timer.Stop();
  running = false;
  exitcode = event.GetExitCode();
  wxEndBusyCursor();

  Unbind(wxEVT_END_PROCESS,&MuniBaseDialog::OnFinish,this);
  Unbind(wxEVT_TIMER,&MuniBaseDialog::OnTimer,this);
}

// --- MuniPhotometryLayer



MuniLayer MuniPhotometryLayer::GetLayer() const
{
  MuniLayer layer(ID_PHOTOMETRY,objects);
  return layer;
}


void MuniPhotometryLayer::DrawObjects(const vector<double>& xcoo, const vector<double>& ycoo, const vector<double>& flux)
{
  wxASSERT(xcoo.size() == ycoo.size() && flux.size() == xcoo.size());

  objects.push_back(new MuniDrawFont(*wxNORMAL_FONT,gold));
  objects.push_back(new MuniDrawPen(wxPen(gold,1.8)));
  objects.push_back(new MuniDrawBrush(wxColour(90,90,255,190)));

  double fmax = 0.0;
  for(size_t i = 0; i < flux.size(); i++)
    if( flux[i] > 0.0 && flux[i] > fmax )
      fmax = flux[i];

  double fmin = fmax;
  for(size_t i = 0; i < flux.size(); i++)
    if( flux[i] > 0.0 && flux[i] < fmin )
      fmin = flux[i];

  const double size1 = hwhm > 1 ? 3*hwhm : 3;
  const double size2 = 2.0;
  const double sizec = 3.0;
  double r1 = (size1 - size2)/(fmax - fmin);

  for(size_t i = 0; i < xcoo.size(); i++) {

    if( flux[i] > 0.0 ) {
      double r = r1*(flux[i] - fmin) + size2;
      objects.push_back(new MuniDrawCircle(xcoo[i],ycoo[i],r));
    }
    else
      objects.push_back(new MuniDrawCross(xcoo[i],ycoo[i],sizec));
  }
}


// --- MuniPhotometry


MuniPhotometry::MuniPhotometry(wxWindow *w, MuniConfig *c, const wxString& f):
  MuniBaseDialog(w,"Photometry",c,f),
  /*fwhm(c->photo_fwhm),thresh(c->photo_thresh),
  saturation(c->photo_satur),readnoise(c->photo_readns),phpadu(c->photo_phpadu),
  */
  erase(false),apply(false),showtooltip(true)
{
  SetIcon(config->munipack_icon);
  //  EnableCloseButton(true);

  wxSpinCtrlDouble *sfwhm = new wxSpinCtrlDouble(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
						 wxDefaultSize,wxSP_ARROW_KEYS,0.0,666.0,fwhm,1.0);
  sfwhm->SetDigits(1);

  wxSpinCtrlDouble *sthresh = new wxSpinCtrlDouble(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
						   wxDefaultSize,wxSP_ARROW_KEYS,0.0,666.0,thresh,0.5);
  sthresh->SetDigits(1);

  wxSpinCtrlDouble *satur = new wxSpinCtrlDouble(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
						   wxDefaultSize,wxSP_ARROW_KEYS,0.0,FLT_MAX,saturation,1e5);
  satur->SetDigits(2);

  wxSpinCtrlDouble *sreadns = new wxSpinCtrlDouble(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
						   wxDefaultSize,wxSP_ARROW_KEYS,0.0,666.0,readnoise,0.1);
  sreadns->SetDigits(1);

  wxSpinCtrlDouble *sphpadu = new wxSpinCtrlDouble(this,wxID_ANY,wxEmptyString,wxDefaultPosition,
						   wxDefaultSize,wxSP_ARROW_KEYS,0.0,666.0,phpadu,0.1);
  sphpadu->SetDigits(1);

  // parameters
  wxFlexGridSizer *fsizer = new wxFlexGridSizer(2);
  fsizer->AddGrowableCol(0);

  fsizer->Add(new wxStaticText(this,wxID_ANY,"FWHM:"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  fsizer->Add(sfwhm,wxSizerFlags().Align(wxALIGN_CENTER));

  fsizer->Add(new wxStaticText(this,wxID_ANY,"Threshold:"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  fsizer->Add(sthresh,wxSizerFlags().Align(wxALIGN_CENTER));

  fsizer->Add(new wxStaticText(this,wxID_ANY,"Saturation:"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  fsizer->Add(satur,wxSizerFlags().Align(wxALIGN_CENTER));

  fsizer->Add(new wxStaticText(this,wxID_ANY,"Read Noise:"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  fsizer->Add(sreadns,wxSizerFlags().Align(wxALIGN_CENTER));

  fsizer->Add(new wxStaticText(this,wxID_ANY,"Photo-e- per ADU:"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  fsizer->Add(sphpadu,wxSizerFlags().Align(wxALIGN_CENTER));

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(fsizer,wxSizerFlags().Center().Border());

  // run processing
  wxBoxSizer *prsizer = new wxBoxSizer(wxHORIZONTAL);
  progress = new wxGauge(this,wxID_ANY,MAXGRANGE);
  prsizer->Add(progress,wxSizerFlags(1).Expand().Border());

  calbutt = new wxButton(this,wxID_ANY,BUTTLABEL);
  prsizer->Add(calbutt,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  topsizer->Add(prsizer,wxSizerFlags().Expand().Border());

  wxSizer *buttsize = CreateSeparatedButtonSizer(wxAPPLY|wxCANCEL);

  if( buttsize ) {
    topsizer->Add(buttsize,wxSizerFlags().Right().Border());
    SetAffirmativeId(wxID_APPLY);
  }
  SetSizerAndFit(topsizer);

  Bind(wxEVT_UPDATE_UI,&MuniPhotometry::OnUpdateUI,this);
  Bind(wxEVT_IDLE,&MuniPhotometry::OnIdle,this);
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniPhotometry::OnFwhm,this,sfwhm->GetId());
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniPhotometry::OnThresh,this,sthresh->GetId());
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniPhotometry::OnSaturation,this,satur->GetId());
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniPhotometry::OnReadNoise,this,sreadns->GetId());
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniPhotometry::OnPhpADU,this,sphpadu->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniPhotometry::OnApply,this,wxID_APPLY);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniPhotometry::OnCancel,this,wxID_CANCEL);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniPhotometry::OnPhotometry,this,calbutt->GetId());

  ids.push_back(sfwhm->GetId());
  ids.push_back(sthresh->GetId());
  ids.push_back(satur->GetId());
  ids.push_back(sreadns->GetId());
  ids.push_back(sphpadu->GetId());

}

MuniPhotometry::~MuniPhotometry()
{
  wxLogDebug("MuniPhotometry::~MuniPhotometry");

  /*
  config->photo_fwhm = fwhm;
  config->photo_thresh = thresh;
  config->photo_satur = saturation;
  */
}

wxString MuniPhotometry::GetBackup() const
{
  return backupfile;
}


void MuniPhotometry::OnUpdateUI(wxUpdateUIEvent&)
{
  wxASSERT(FindWindow(wxID_APPLY));
  FindWindow(wxID_APPLY)->Enable((apply && !running) && !readonly);

  wxASSERT(FindWindow(wxID_CANCEL));
  FindWindow(wxID_CANCEL)->Enable(!running);

  //  EnableCloseButton(!running);

  for(list<int>::const_iterator i = ids.begin(); i != ids.end(); ++i)
    FindWindow(*i)->Enable(!running);
}

void MuniPhotometry::OnIdle(wxIdleEvent&)
{
  if( showtooltip ) {
    showtooltip = false;

    wxASSERT(FindWindow(wxID_APPLY));
    if( readonly ) {
      wxRichToolTip tip("Read Only Data",
			"This file is not writable.\n"
			"Any results will not saved.");
      tip.SetIcon(wxICON_WARNING);
      tip.ShowFor(FindWindow(wxID_APPLY));
    }
  }

}


void MuniPhotometry::OnFwhm(wxSpinDoubleEvent& event)
{
  fwhm = event.GetValue();
}

void MuniPhotometry::OnThresh(wxSpinDoubleEvent& event)
{
  thresh = event.GetValue();
}

void MuniPhotometry::OnSaturation(wxSpinDoubleEvent& event)
{
  saturation = event.GetValue();
}

void MuniPhotometry::OnReadNoise(wxSpinDoubleEvent& event)
{
  readnoise = event.GetValue();
}

void MuniPhotometry::OnPhpADU(wxSpinDoubleEvent& event)
{
  phpadu = event.GetValue();
}

void MuniPhotometry::OnCancel(wxCommandEvent& ev)
{
  wxLogDebug("Leaving Photometry...");
  EraseCanvas();
  EndModal(wxID_CANCEL);
}

void MuniPhotometry::OnApply(wxCommandEvent& ev)
{
  wxBusyCursor wait;
  wxLogDebug("Applying Photometry...");

  backupfile = file + config->backup_suffix;
  wxCopyFile(file,backupfile);
  wxCopyFile(workingfile,file);

  EraseCanvas();
  EndModal(wxID_OK);
}

void MuniPhotometry::OnPhotometry(wxCommandEvent& event)
{
  wxLogDebug("Running OnPhotometry...");

  Bind(wxEVT_END_PROCESS,&MuniPhotometry::OnFinish,this);
  Bind(wxEVT_TIMER,&MuniPhotometry::OnTimer,this);
  calbutt->SetLabel("STOP Processing");

  MuniProcess *action = new MuniProcess(&pipe,"aphot");

  action->Write("PIPELOG = T");
  action->Write("FWHM = %lf",fwhm);
  action->Write("THRESHOLD = %lf",thresh);
  action->Write("FILE = '"+workingfile+"' ''");

  wxWindow *bapply = FindWindowById(wxID_APPLY,this);
  bapply->Enable(false);

  Layout();

  apply = false;
  EraseCanvas();
  RunProcess(action);
}

void MuniPhotometry::OnFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniPhotometry::OnFinish");

  Unbind(wxEVT_END_PROCESS,&MuniPhotometry::OnFinish,this);
  Unbind(wxEVT_TIMER,&MuniPhotometry::OnTimer,this);
  calbutt->SetLabel(BUTTLABEL);

  wxWindow *bapply = FindWindowById(wxID_APPLY,this);
  bapply->Enable(true);

  Layout();

  if( event.GetExitCode() != 0 )
    wxLogDebug("Failed with exit code %d",event.GetExitCode());
  else {
    // LoadFile(result);
    //wxCommandEvent e(wxEVT_COMMAND_MENU_SELECTED,wxID_OPEN);
    //    static_cast<MuniView *>(GetGrandParent())->LoadFile(file);

    /*
    wxCommandEvent e(EVT_FILELOAD,GetId());
    e.SetEventObject(this);
    e.SetString(file);
    wxQueueEvent(GetParent(),e.Clone());
    */

    apply = ParseProcessing(pipe.GetOutput());


  }
}


// parse output of astrometry fit and fill the data table
void MuniPhotometry::OnTimer(wxTimerEvent& event)
{
  progress->Pulse();
  ParseProcessing(pipe.GetOutput());
}

bool MuniPhotometry::ParseProcessing(const wxArrayString& out)
{
  vector<double> xcoo,ycoo,flux;

  wxRegEx re("^=(.*)> (.+)");
  wxASSERT(re.IsValid());

  for(size_t i = 0; i < out.GetCount(); i++) {

    if( re.Matches(out[i]) ) {

      wxString key(re.GetMatch(out[i],1));
      wxString value(re.GetMatch(out[i],2));

      if( key == "APHOT" ) {

	wxStringInputStream ss(value);
	wxTextInputStream t(ss);

	double x,y,f;
	t >> x >> y >> f;
	xcoo.push_back(x);
	ycoo.push_back(y);
	flux.push_back(f);
      }
      else if( key == "FIND" ) {

	wxStringInputStream ss(value);
	wxTextInputStream t(ss);

	double x,y;
	t >> x >> y;
	xcoo.push_back(x);
	ycoo.push_back(y);
	flux.push_back(-1.0);
      }

    }
  }

  if( ! xcoo.empty() ) {
    MuniPhotometryLayer layer;
    layer.DrawObjects(xcoo,ycoo,flux);

    MuniPhotometryEvent ev(EVT_PHOTOMETRY);
    ev.layer = layer.GetLayer();
    wxQueueEvent(GetParent(),ev.Clone());
    erase = true;

    return true;
  }

  return false;
}

void MuniPhotometry::EraseCanvas()
{
  if( erase ) {
    MuniPhotometryEvent e(EVT_PHOTOMETRY);
    e.erase = true;
    wxQueueEvent(GetParent(),e.Clone());
  }
}
