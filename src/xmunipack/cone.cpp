/*

  cone search

  Copyright © 2012-4, 2019-22 F.Hroch (hroch@physics.muni.cz)

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
#include "help.h"
#include <wx/wx.h>
#include <wx/tokenzr.h>
#include <wx/animate.h>
#include <wx/statline.h>
#include <wx/regex.h>
#include <wx/richtooltip.h>
#include <list>
#include <map>


using namespace std;


// --- MuniCone
MuniCone::MuniCone(wxWindow *w, MuniConfig *c, const wxString& o):
  wxDialog(w,wxID_ANY,"Cone search"),config(c),mproc(0),timer(this),
  throbber(0),apply(false),object(o),ra(-999.9),dec(-99.9),index(0)
{
  SetIcon(config->munipack_icon);

  wxSizerFlags lf, cf;
  lf.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).Border(wxRIGHT);
  cf.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT);

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);

  // Name resolver
  wxBoxSizer *namesizer = new wxBoxSizer(wxHORIZONTAL);
  namesizer->Add(new wxStaticText(this,wxID_ANY,"Object"),
		 wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL));
  objentry = new wxTextCtrl(this,wxID_ANY,object,wxDefaultPosition,
			    wxDefaultSize,wxTE_PROCESS_ENTER);
  namesizer->Add(objentry,wxSizerFlags(1).Border(wxRIGHT|wxLEFT));

  wxButton *getcoo = new wxButton(this,wxID_ANY,"Get Coordinates");
  namesizer->Add(getcoo);
  topsizer->Add(namesizer,wxSizerFlags().Expand().DoubleBorder());

  topsizer->Add(new wxStaticLine(this,wxID_ANY),
		wxSizerFlags().Expand().DoubleBorder());

  // Coordinates
  wxBoxSizer *coosizer = new wxBoxSizer(wxHORIZONTAL);
  wxString tip(" enter a value in degrees with decimal point"
	       " or use the sexadecimal format");
  wxSizerFlags coo_label;
  coo_label.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT);

  // RA
  coosizer->Add(new wxStaticText(this,wxID_ANY,L"α [°]"),coo_label);
  alpha = new wxTextCtrl(this,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
			 wxTE_PROCESS_ENTER);
  alpha->SetToolTip("Right Ascension:"+tip);
  coosizer->Add(alpha,wxSizerFlags(1).DoubleBorder(wxRIGHT));

  // Dec
  coosizer->Add(new wxStaticText(this,wxID_ANY,L"δ [°]"),coo_label);
  delta = new wxTextCtrl(this,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
			 wxTE_PROCESS_ENTER);
  delta->SetToolTip("Declination:"+tip);
  coosizer->Add(delta,wxSizerFlags(1));
  topsizer->Add(coosizer,wxSizerFlags().Expand().DoubleBorder());

  // Options
  wxFlexGridSizer *grid = new wxFlexGridSizer(2);

  // radius
  grid->Add(new wxStaticText(this,wxID_ANY,L"Cone radius [°]"),lf);

  radius=new wxSpinCtrlDouble(this,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
		    wxSP_ARROW_KEYS|wxTE_PROCESS_ENTER,0.0,100.0,0.1,0.1);
  grid->Add(radius,cf.Border());

  // magnitudes
  magmin=new wxSpinCtrlDouble(this,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
		      wxSP_ARROW_KEYS|wxTE_PROCESS_ENTER,-30.0,30.0,0.0,0.5);
  magmax=new wxSpinCtrlDouble(this,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
	              wxSP_ARROW_KEYS|wxTE_PROCESS_ENTER,-30.0,30.0,15.0,0.5);

  grid->Add(new wxStaticText(this,wxID_ANY,"Magnitude (brightest)"),lf);
  grid->Add(magmin,cf.Border(wxLEFT|wxTOP));
  grid->Add(new wxStaticText(this,wxID_ANY,"Magnitude (faintest)"),lf);
  grid->Add(magmax,cf.Border(wxLEFT|wxBOTTOM));

  grid->Add(new wxStaticText(this,wxID_ANY,"Catalogue"),lf);
  wxArrayString catalogs;
  vector<VOCatResources> cats(catconf.GetCatalogues());
  int idefault = 0;
  for(size_t i = 0; i < cats.size(); i++) {
    if( cats[i].GetName() == catconf.GetName() )
      idefault = i;
    catalogs.Add(cats[i].GetName());
  }
  wxChoice *cat =
    new wxChoice(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,catalogs);
  cat->SetSelection(idefault);
  grid->Add(cat,cf.Border().Align(wxALIGN_LEFT));

  // Johnson's patch
  johnson_label = new wxStaticText(this,wxID_ANY,"Johnson's patch");
  grid->Add(johnson_label,lf);
  johnson = new wxCheckBox(this,wxID_ANY,"");
  grid->Add(johnson,cf.Border());
  johnson->SetValue(true);
  johnson->SetToolTip("Gunn's ri magnitudes are converted to Johnson RI (UCAC4 only)");

  topsizer->Add(grid,wxSizerFlags().Centre().DoubleBorder());

  topsizer->Add(new wxStaticLine(this,wxID_ANY),
		wxSizerFlags().Expand().Border());

  wxBoxSizer *sb = new wxBoxSizer(wxHORIZONTAL);

  status = new wxStaticText(this,wxID_ANY,"Some description");
  status->SetFont(*wxSMALL_FONT);
  status->Show(false);
  sb->Add(status,wxSizerFlags(1).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).DoubleBorder().ReserveSpaceEvenIfHidden());

  search = new wxButton(this,wxID_ANY,"Search");
  sb->Add(search,wxSizerFlags().DoubleBorder().Align(wxALIGN_RIGHT));
  stop = new wxButton(this,wxID_ANY,"Stop");
  sb->Add(stop,wxSizerFlags().DoubleBorder().Align(wxALIGN_RIGHT));
  stop->Show(false);
  stop->SetLabelMarkup("<span fgcolor=\"#ff0000\"><b>Stop</b></span>");

  topsizer->Add(sb,wxSizerFlags().Expand().Border());

  wxBoxSizer *footsizer = new wxBoxSizer(wxHORIZONTAL);
  throbber = new wxAnimationCtrl(this,wxID_ANY,c->throbber,wxDefaultPosition,
				 wxSize(32,32));
  footsizer->Add(throbber,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).DoubleBorder(wxLEFT|wxRIGHT).ReserveSpaceEvenIfHidden());
  throbber->Show(false);
  footsizer->AddStretchSpacer(1);
  wxSizer *buttsize = CreateButtonSizer(wxOK|wxCANCEL|wxHELP);
  if( buttsize )
    footsizer->Add(buttsize,wxSizerFlags().Right().Border());
  topsizer->Add(footsizer,wxSizerFlags().Expand());

  SetSizerAndFit(topsizer);

  Bind(wxEVT_UPDATE_UI,&MuniCone::OnUpdateUI,this);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniCone::OnSearch,this,search->GetId());
  Bind(wxEVT_TEXT,&MuniCone::OnObjectName,this,objentry->GetId());
  Bind(wxEVT_TEXT_ENTER,&MuniCone::OnObjectEnter,this,objentry->GetId());
  Bind(wxEVT_BUTTON,&MuniCone::OnGetCoo,this,getcoo->GetId());
  Bind(wxEVT_TEXT,&MuniCone::OnRightAscension,this,alpha->GetId());
  Bind(wxEVT_TEXT,&MuniCone::OnDeclination,this,delta->GetId());
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniCone::OnService,this,cat->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniCone::OnOk,this,wxID_OK);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniCone::OnHelp,this,wxID_HELP);
  Bind(wxEVT_BUTTON,&MuniCone::OnStop,this,stop->GetId());

  ids.push_back(cat->GetId());
  ids.push_back(alpha->GetId());
  ids.push_back(delta->GetId());
  ids.push_back(radius->GetId());
  ids.push_back(magmax->GetId());
  ids.push_back(magmin->GetId());
  ids.push_back(getcoo->GetId());
  ids.push_back(johnson->GetId());

  tmpfile = wxFileName::CreateTempFileName("xmunipack-cone_");
}

MuniCone::~MuniCone()
{
  if( ! tmpfile.IsEmpty() && wxFileExists(tmpfile) )
    wxRemoveFile(tmpfile);
}

void MuniCone::OnClose(wxCloseEvent& event)
{
  //  wxLogDebug("MuniCone::OnClose");

  // stop a potentially working subprocess
  if( mproc )
    mproc->Kill();

  // Skipping here passes processing of Close event to default handler.
  event.Skip();
}

void MuniCone::OnOk(wxCommandEvent& event)
{
  wxLogDebug("Applying Cone...");
  event.Skip();

  wxString cone = wxFileName::CreateTempFileName("ConeSearch_") + ".fits";
  FitsCopyFile(tmpfile,cone);

  wxCommandEvent ev(EVT_FINISH_DIALOG,ID_CONE);
  ev.SetString(cone);
  wxQueueEvent(GetParent(),ev.Clone());
}

void MuniCone::OnHelp(wxCommandEvent& event)
{
  MuniHelp(config->munipack_html_dir,"man_cone.html");
}

void MuniCone::OnUpdateUI(wxUpdateUIEvent& event)
{
  wxASSERT(search);

  wxASSERT(FindWindow(wxID_OK));
  FindWindow(wxID_OK)->Enable(apply && !mproc);

  wxASSERT(FindWindow(wxID_CANCEL));
  FindWindow(wxID_CANCEL)->Enable(!mproc);

  search->Enable(ra > -999 && dec > -99 && catconf.IsOk());
  bool morgan = catconf.GetName() == "UCAC4";
  johnson->Enable(morgan);
  johnson_label->Enable(morgan);

  for(list<int>::const_iterator i = ids.begin(); i != ids.end(); ++i)
    FindWindow(*i)->Enable(!mproc);
}

void MuniCone::OnObjectName(wxCommandEvent& event)
{
  object = event.GetString();
}

void MuniCone::OnObjectEnter(wxCommandEvent& event)
{
  object = event.GetString();
  Resolve();
}

void MuniCone::OnGetCoo(wxCommandEvent& event)
{
  Resolve();
}

void MuniCone::OnRightAscension(wxCommandEvent& event)
{
  double r = deg(event.GetString(),15.0);
  //  wxLogDebug("MuniCone::OnRightAscension %f",r);

  if( 0.0 <= r && r <= 360.0 ) {
    ra = r;
    alpha->SetForegroundColour(wxNullColour);
  }
  else {
    ra = -999.9;
    alpha->SetForegroundColour(*wxRED);
  }
}

void MuniCone::OnDeclination(wxCommandEvent& event)
{
  double d = deg(event.GetString());
  //  wxLogDebug("MuniCone::OnDeclination %f",d);

  if( -90.0 <= d && d <= 90.0 ) {
    dec = d;
    delta->SetForegroundColour(wxNullColour);
  }
  else {
    dec = -99.9;
    delta->SetForegroundColour(*wxRED);
  }
}


void MuniCone::Resolve()
{
  Bind(wxEVT_END_PROCESS,&MuniCone::ResolveFinish,this);

  mproc = new MuniProcess(this,"sesame");
  mproc->SetEcho(false);
  mproc->Write("OBJECT = '"+object+"'");
  mproc->OnStart();

  wxBeginBusyCursor();
  throbber->Show(true);
  throbber->Play();
  Layout();
}

void MuniCone::ResolveFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniCone::ResolveFinish");

  if( event.GetExitCode() != 0 ) {
    objentry->SetSelection(-1,-1);
    wxRichToolTip tip("Unknown object",
		      "The object has been not found\nby Sesame resolver.");
    tip.SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWFRAME ));
    tip.SetIcon(wxICON_ERROR);
    tip.ShowFor(objentry);
  }
  else {
    wxString ra,dec;
    wxString line(mproc->GetOutput().Last());
    wxLogDebug(line);
    wxStringTokenizer tk(line);
    for(int n = 0; tk.HasMoreTokens(); n++) {
      wxString item = tk.GetNextToken();
      double q;
      if( n == 0 && item.ToDouble(&q) )
	ra = item;
      if( n == 1 && item.ToDouble(&q) )
	dec = item;
    }
    if( ! ra.IsEmpty() && ! dec.IsEmpty() ) {
      alpha->SetValue(ra);
      delta->SetValue(dec);
    }
  }

  Unbind(wxEVT_END_PROCESS,&MuniCone::ResolveFinish,this);

  wxEndBusyCursor();
  throbber->Stop();
  throbber->Show(false);
  delete mproc;
  mproc = 0;
  Layout();
  UpdateWindowUI();
}


void MuniCone::OnService(wxCommandEvent& event)
{
  catconf.SetCat(event.GetString());
}

void MuniCone::OnStop(wxCommandEvent& event)
{
  wxASSERT(mproc);
  mproc->Kill();
}

void MuniCone::OnSearch(wxCommandEvent& event)
{
  if( wxFileExists(tmpfile) )
    wxRemoveFile(tmpfile);

  Bind(wxEVT_END_PROCESS,&MuniCone::OnFinish,this);
  Bind(wxEVT_TIMER,&MuniCone::OnTimer,this);

  mproc = new MuniProcess(this,"cone");
  mproc->SetEcho(false);

  map<wxString,wxString> pars;
  pars["SR"] = wxString::FromCDouble(radius->GetValue());
  pars["RA"] = wxString::FromCDouble(ra);
  pars["DEC"] = wxString::FromCDouble(dec);

  mproc->Write("PIPELOG = T");
  //  mproc->Write("VERBOSE = T");
  mproc->Write("SORT = '"+catconf.GetSort()+"'");
  mproc->Write("MAGMIN = "+wxString::FromCDouble(magmin->GetValue()));
  mproc->Write("MAGMAX = "+wxString::FromCDouble(magmax->GetValue()));
  mproc->Write("TYPE = 'FITS'");
  mproc->Write("URL = '"+catconf.GetUrl(pars)+"'");
  mproc->Write("OUTPUT = '"+tmpfile+"'");
  if( johnson->GetValue() == 1 && catconf.GetName() == "UCAC4" )
    mproc->Write("PATCH = T");

  timer.Start(250);
  throbber->Show(true);
  throbber->Play();
  status->Show(true);
  status->SetLabel("Connecting Virtual observatory...");
  search->Show(false);
  stop->Show(true);
  apply = false;
  index = 0;
  Layout();
  UpdateWindowUI();

  mproc->OnStart();
}

void MuniCone::OnFinish(wxProcessEvent& event)
{
  //  wxLogDebug("MuniCone::OnFinish");

  if( event.GetExitCode() != 0 ) {
    apply = false;
    wxLogDebug("Failed with exit code %d",event.GetExitCode());
    status->SetLabelMarkup("<span fgcolor=\"#ff0000\"><b>Cone search failed.</b></span>");
    wxLogError("Command: cone");
    wxArrayString err(mproc->GetErrors());
    wxArrayString out(mproc->GetOutput());
    wxArrayString in(mproc->GetInput());
    for(size_t i = 0; i < in.GetCount(); i++)
      wxLogError("Input:"+in[i]);
    for(size_t i = 0; i < out.GetCount(); i++)
      wxLogError("Output:"+out[i]);
    for(size_t i = 0; i < err.GetCount(); i++)
      wxLogError("Error:"+err[i]);
    wxLogError("Cone search failed.");
  }
  else {
    wxString msg = ParseOutput(mproc->GetOutput());
    status->SetLabel(msg);
    apply = msg.Find("No objects found") == wxNOT_FOUND;
  }

  Unbind(wxEVT_END_PROCESS,&MuniCone::OnFinish,this);
  Unbind(wxEVT_TIMER,&MuniCone::OnTimer,this);

  timer.Stop();
  throbber->Stop();
  throbber->Show(false);
  search->Show(true);
  stop->Show(false);

  delete mproc;
  mproc = 0;

  Layout();
  UpdateWindowUI();
}

void MuniCone::OnTimer(wxTimerEvent& event)
{
  //  wxLogDebug("OnTimer");
  wxArrayString out(mproc->GetOutput());

  if( out.GetCount() == 0 )
    return;

  wxString msg = ParseOutput(out);
  status->SetLabel(msg);
  status->Layout();
}


wxString MuniCone::ParseOutput(const wxArrayString& out)
{
  wxString msg;

  wxRegEx re("^=(.*)> (.+)");
  wxASSERT(re.IsValid());

  for(size_t i = index; i < out.GetCount(); i++) {

    if( re.Matches(out[i]) ) {

      wxString key(re.GetMatch(out[i],1));
      wxString value(re.GetMatch(out[i],2));

      if( key == "CONE" )
	msg = value;

    }
  }
  index = out.GetCount() - 1;

  return msg;
}

double MuniCone::deg(const wxString& line, double sexadecimal)
{
  /* This routine is adapted Fortran code distributed along
     with USNO 1.0 and 2.0 catalogues at mid '90. No reference
     is available recently.
  */

  if( line.IsEmpty() )
    return -999.9;

  double dp = 1.0;
  double piece[4] = { 0.0, 0.0, 0.0, 0.0}; // p[0] is left unused
  size_t n = 1;
  bool dpfind = false;
  double sgn = 1.0;

  for(wxString::const_iterator i = line.begin(); i != line.end(); i++) {
    wxChar c = *i;

    if( c == '-' )
      sgn = -1.0;
    else if( c == ':' || c == ' ' ) {
      dpfind = false;
      n = n + 1;
      if( n > 3 ) break;
    }
    else if( c == '.' )
      dpfind = true;
    else if( '0' <= c && c <= '9' ) {
      double z = int(c) - int('0');
      if( dpfind ) {
	dp = dp / 10.0;
	piece[n] = piece[n] + dp*z;
      }
      else
	piece[n] = 10*piece[n] + z;
    }
    else
      return -999.9;
  }

  for(size_t i = 2; i <= 3; i++) {
    if( piece[i] < 0.0 || piece[i] > 60.0 )
      return -999.9;
  }

  double six = n > 2 ? sexadecimal : 1.0;
  return six * sgn * (piece[1] + (piece[2] + piece[3] / 60.0) / 60.0);
}
