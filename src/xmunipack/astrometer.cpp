 /*

  xmunipack - batch astrometry calibration

  Copyright © 2012-5 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/statline.h>
#include <wx/dataview.h>
#include <wx/richtooltip.h>
#include <wx/sstream.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>
#include <wx/regex.h>
#include <list>
#include <fitsio.h>
#include <string.h>

#define BUTTLABEL "Process All Frames"

#define MAS    "mas"
#define ARCSEC "arcsec"
#define ARCMIN "arcmin"
#define ARCDEG "deg"

#define PROJ_IDENT " "
#define PROJ_GNO "Gnomonic"

#define Horigin L"תשע״ב"

using namespace std;



MuniAstrometer::MuniAstrometer(wxWindow *w, MuniConfig *c, const vector<FitsMeta>& l):
  wxDialog(w,wxID_ANY,"Astrometry",wxDefaultPosition,wxDefaultSize,wxDEFAULT_DIALOG_STYLE|wxRESIZE_BORDER),
  config(c),
  list(l),fpick(0),catbutt(0),
  timer(this),pipe(this),findex(-1),lastrow(0),
  proj(config->astrometry_proj),
  running(false),relative(false), tmpcat(false)
{
  SetIcon(config->munipack_icon);

  CreateControls();

  SetTable();
}

MuniAstrometer::~MuniAstrometer()
{
  config->astrometry_proj = proj;
  EraseTemp();
}

void MuniAstrometer::CreateControls()
{
  wxSizerFlags lf, cf, sl, rl;
  lf.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT);
  cf.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT);
  sl.Border().Expand();
  rl.Align(wxALIGN_CENTER_VERTICAL).Border(wxRIGHT);

  wxStaticText *label;

  wxFont bf(*wxNORMAL_FONT);
  bf.SetWeight(wxFONTWEIGHT_BOLD);

  wxFont sf(*wxSMALL_FONT);

  // title
  topsizer = new wxBoxSizer(wxVERTICAL);


  label = new wxStaticText(this,wxID_ANY,"Reference:");
  label->SetFont(sf);
  topsizer->Add(label,wxSizerFlags().DoubleBorder(wxRIGHT|wxLEFT|wxTOP));

  wxBoxSizer *rhsizer = new wxBoxSizer(wxHORIZONTAL);

  refcatid = new wxStaticText(this,wxID_ANY,Horigin);
  rhsizer->Add(refcatid,wxSizerFlags(1).Border().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT));

  wxButton *refbutt = new wxButton(this,wxID_ANY,"Reference...");
  rhsizer->Add(refbutt,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  topsizer->Add(rhsizer,wxSizerFlags().Expand().Border());

  topsizer->Add(new wxStaticLine(this,wxID_ANY),sl);
  label = new wxStaticText(this,wxID_ANY,"Parameters:");
  label->SetFont(sf);
  topsizer->Add(label,wxSizerFlags().DoubleBorder(wxLEFT));

  wxFlexGridSizer *grid = new wxFlexGridSizer(2);
  grid->AddGrowableCol(1);

  label = new wxStaticText(this,wxID_ANY,"Projection:");
  grid->Add(label,lf);
  wxArrayString projections;
  projections.Add(PROJ_IDENT);
  projections.Add(PROJ_GNO);
  wxChoice *chproj = new wxChoice(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,projections);
  chproj->SetSelection(chproj->FindString(proj));
  grid->Add(chproj,cf);
  topsizer->Add(grid,wxSizerFlags().Center());

  topsizer->Add(new wxStaticLine(this,wxID_ANY),wxSizerFlags().Expand().Border(wxLEFT|wxRIGHT|wxTOP));

  mtable = new wxDataViewListCtrl(this,wxID_ANY);
  mtable->AppendTextColumn("Files");
  mtable->AppendTextColumn("Status",wxDATAVIEW_CELL_ACTIVATABLE,100);
  topsizer->Add(mtable,wxSizerFlags(1).Border().Expand());

  gstat = new wxGauge(this,wxID_ANY,100,wxDefaultPosition,
		      wxDefaultSize,wxGA_HORIZONTAL|wxGA_SMOOTH);
  topsizer->Add(gstat,wxSizerFlags().Expand().Border(wxRIGHT|wxLEFT|wxBOTTOM));

  astropt = new MuniAstrometryOptions(this,config);
  topsizer->Add(astropt,wxSizerFlags().Border().Expand());

  // show/hide via CollapsiblePane

  butt = new wxButton(this,wxID_ANY,BUTTLABEL);
  topsizer->Add(butt,wxSizerFlags().Border().Expand());

  SetSizerAndFit(topsizer);


  Bind(wxEVT_CLOSE_WINDOW,&MuniAstrometer::OnClose,this);
  Bind(wxEVT_UPDATE_UI,&MuniAstrometer::OnUpdateUI,this);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometer::OnReference,this,refbutt->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometer::OnProcess,this,butt->GetId());
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniAstrometer::OnChoiceProj,this,chproj->GetId());

  ids.push_back(chproj->GetId());
  ids.push_back(refbutt->GetId());
}

void MuniAstrometer::OnClose(wxCloseEvent& event)
{
  wxLogDebug("MuniAstrometer::OnClose %s",running ? "T" : "F");

  // Arranging of the code is really important
  // A use of Close() doesn't work because the event-system have no break
  // to execute other events (finishing subprocesses).
  // Perhaps, emmiting of a new signal starts next event-loop, whilst
  // Close or event.Skip() just passed signal to other routines, so
  // dialog can be exit before the event arrived.
  if( running ) {
    pipe.Stop();
    wxQueueEvent(this,event.Clone());
  }
  else
    // Skipping here passes processing of Close event to default handler.
    event.Skip();
}


void MuniAstrometer::OnUpdateUI(wxUpdateUIEvent& e)
{
  for(std::list<int>::const_iterator i = ids.begin(); i != ids.end(); ++i)
    FindWindow(*i)->Enable(!running);

  wxASSERT(astropt);
  astropt->Enable(!running);

  butt->Enable(!(catfile == "" && reffile == ""));
}

void MuniAstrometer::SetTable()
{
  for(size_t i = 0; i < list.size(); i++) {
    wxVector<wxVariant> line;
    line.push_back(list[i].GetName());
    line.push_back("");
    mtable->AppendItem(line);
  }
}

void MuniAstrometer::OnChoiceProj(wxCommandEvent& event)
{
  proj = event.GetString();
}

void MuniAstrometer::OnReference(wxCommandEvent& event)
{
  MuniSelectSource ms(this,config,true);
  if( ms.ShowModal() == wxID_OK ) {

    EraseTemp();

    if( ms.GetType() == ID_ASTRO_CAT ) {

      tmpcat = ms.IsTemporary();
      if( tmpcat ) {
	catfile = wxFileName::CreateTempFileName("xmunipack-astrometer_");
	wxCopyFile(ms.GetPath(),catfile);
      }
      else
	catfile = ms.GetPath();
      label_ra = ms.GetLabelRA();
      label_dec = ms.GetLabelDec();
      label_pmra = ms.GetLabelPMRA();
      label_pmdec = ms.GetLabelPMDec();
      label_mag = ms.GetLabelMag();
      refcatid->SetLabel(ms.GetId());
    }
    else if( ms.GetType() == ID_ASTRO_REF ) {
      reffile = ms.GetPath();
      relative = ms.GetRelative();
      label_ra.Clear();
      label_dec.Clear();
      label_pmra.Clear();
      label_pmdec.Clear();
      label_mag.Clear();
      refcatid->SetLabel(ms.GetId());
    }
    Layout();
  }
}

void MuniAstrometer::EraseTemp()
{
  if( tmpcat && !catfile.IsEmpty() && wxFileExists(catfile) )
    wxRemoveFile(catfile);
}

void MuniAstrometer::OnProcess(wxCommandEvent& event)
{
  // check input
  if( reffile.IsEmpty() && catfile.IsEmpty() ) {
    wxRichToolTip tip("Undefined Reference",
		      "Please define an astrometric\n"
                      "catalogue or an reference frame.");
    tip.SetIcon(wxICON_WARNING);
    if( reffile.IsEmpty() && fpick != 0 )
      tip.ShowFor(fpick);
    if( catfile.IsEmpty() && catbutt != 0)
      tip.ShowFor(catbutt);
    return;
  }

  if( list.empty() ) {
    wxRichToolTip tip("Nothing to Calibrate",
		      "Please select some files\n"
                      "to astrometry calibration.");
    tip.SetIcon(wxICON_WARNING);
    tip.ShowFor(mtable);
    return;
  }


  if( running ) {

    pipe.Stop();
    if( findex >= 0 )
      mtable->SetTextValue("Interrupted.",findex,1);
  }

  else {

    Bind(wxEVT_END_PROCESS,&MuniAstrometer::OnFinish,this);
    Bind(wxEVT_TIMER,&MuniAstrometer::OnTimer,this);

    running = true;
    butt->SetLabel("STOP Processing");
    CreateProcess();
    timer.Start(100);
    pipe.Start();
    lastrow = 0;
    findex = 0;
  }
}


void MuniAstrometer::CreateProcess()
{
  wxASSERT(astropt);

  MuniProcess *c = new MuniProcess(&pipe,"astrometry");
  pipe.push(c);

  c->Write("PIPELOG = T");
  c->Write("PROJECTION = '"+proj.Upper()+"'");
  c->Write("WCSSAVE = T");
  c->Write("AUNITS = '"+astropt->GetOutputUnits()+"'");
  c->Write("SIG = %e",astropt->GetSig());
  c->Write("FSIG = %e",astropt->GetFSig());


  if( label_ra != "" && label_dec != "" ) {
    c->Write("LABEL_RA = '"+label_ra+"'");
    c->Write("LABEL_DEC = '"+label_dec+"'");
    if( label_pmra != "" )  c->Write("LABEL_PMRA = '" + label_pmra + "'");
    if( label_pmdec != "" ) c->Write("LABEL_PMDEC = '" + label_pmdec + "'");
  }

  c->Write("MATCH = 'BACKTRACKING'");
  c->Write("MINMATCH = %ld",long(astropt->GetMinMatch()));
  c->Write("MAXMATCH = %ld",long(astropt->GetMaxMatch()));

  if( ! reffile.IsEmpty() ) {
    if( astrorel )
      c->Write("REL = '" + reffile + "'");
    else
      c->Write("REF = '" + reffile + "'");
  }

  else if( ! catfile.IsEmpty() )
    c->Write("CAT = '" + catfile + "'");

  for(size_t i = 0; i < list.size(); i++)
    c->Write("FILE = '"+list[i].GetFullPath()+"' ''");

}

void MuniAstrometer::OnFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniAstrometer::OnFinish");

  running = false;
  butt->SetLabel(BUTTLABEL);
  timer.Stop();

  Unbind(wxEVT_END_PROCESS,&MuniAstrometer::OnFinish,this);
  Unbind(wxEVT_TIMER,&MuniAstrometer::OnTimer,this);

  ParseOutput();
}

void MuniAstrometer::OnTimer(wxTimerEvent& event)
{
  ParseOutput();
}

void MuniAstrometer::ParseOutput()
{
  int x = gstat->GetRange();

  wxArrayString out(pipe.GetOutput());

  wxRegEx re("^=(.*)> (.+)");
  wxASSERT(re.IsValid());

  for(size_t j = lastrow; j < out.GetCount(); j++, lastrow++) {

    if( re.Matches(out[j]) ) {

      wxString key(re.GetMatch(out[j],1));
      wxString value(re.GetMatch(out[j],2));

      if( key == "ASTROMETRY" ) {

	if( value.Find("Start") != wxNOT_FOUND ) {
	  gstat->SetValue(0);
	  mtable->SetTextValue(value,findex,1);
	}
	else if( value.Find("Finish") != wxNOT_FOUND ) {
	  gstat->SetValue(0);
	  mtable->SetTextValue(value,findex,1);
	  findex++;
	}

      }
      else if( key == "MPROGRES2" ) {

	wxStringInputStream ss(value);
	wxTextInputStream t(ss);

	double ns, ntot;
	t >> ns >> ntot;

	gstat->SetValue(wxMin(int(x*(ns/ntot)),x));
      }

    }
  }
}
