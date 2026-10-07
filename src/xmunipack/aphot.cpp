/*

  Aperture photometry tool

  Copyright © 2019 F.Hroch (hroch@physics.muni.cz)

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

#include "aphot.h"
#include "mtool.h"
#include "mconfig.h"
#include "types.h"
#include "event.h"
#include "fits.h"
#include <wx/wx.h>
#include <wx/regex.h>
#include <wx/tokenzr.h>
#include <wx/graphics.h>
#include <wx/spinctrl.h>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cfloat>

using namespace std;



MuniAphot::MuniAphot(wxWindow *w, MuniConfig *c, const FitsFile& fits,
		     const wxImage& i):
  MuniTool(w,c,"aphot","Aperture photometry"), fitsname(fits.GetFullPath()),
  image(i), aproc(0), naper(12), zoom(4), ring_min(20), ring_max(30),
  snap(true), spiral(true)
{
  naper = config->aphot_naper;
  ring_min = config->aphot_rmin;
  ring_max = config->aphot_rmax;
  spiral = config->aphot_spiral;
  snap = config->aphot_snap;
  zoom = config->aphot_zoom;
  saper = config->aphot_saper;

  wxPanel *panel = new wxPanel(this,wxID_ANY);

  wxSpinCtrl *naper_spin, *rmin_spin, *rmax_spin;
  wxCheckBox *snap_check;
  wxSlider *zoom_slider;

  canvas = new MuniAphotDisplay(panel,image,naper,saper,zoom,ring_min,
				ring_max,spiral);
  naper_spin = new wxSpinCtrl(panel,wxID_ANY,"",wxDefaultPosition,
			      wxDefaultSize,wxSP_ARROW_KEYS,2,666,naper);
  wxRadioButton *spiral_radio =
    new wxRadioButton(panel,ID_APHOT_SPIRAL,"Spiral",
		      wxDefaultPosition,wxDefaultSize,wxRB_GROUP);
  wxRadioButton *equidistant_radio =
    new wxRadioButton(panel,ID_APHOT_EQL,"Equidistant");
  if( spiral )
    spiral_radio->SetValue(true);
  else
    equidistant_radio->SetValue(true);

  rmin_spin = new wxSpinCtrl(panel,wxID_ANY,"",wxDefaultPosition,
			    wxDefaultSize,wxSP_ARROW_KEYS,1,666,ring_min);
  rmax_spin = new wxSpinCtrl(panel,wxID_ANY,"",wxDefaultPosition,
			     wxDefaultSize,wxSP_ARROW_KEYS,1,666,ring_max);

  snap_check = new wxCheckBox(panel,wxID_ANY,"To Snap");
  snap_check->SetToolTip("If checked, snap to a local maximum brightness.");
  snap_check->SetValue(snap);

  zoom_slider = new wxSlider(panel,wxID_ANY,zoom,1,10,wxDefaultPosition,
			     wxDefaultSize,wxSL_HORIZONTAL|wxSL_BOTTOM);

  wxSizerFlags label_flag, entry_flag, label_radio, label_slide;
  label_flag.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).Border(wxRIGHT);
  label_radio.Align(wxALIGN_TOP|wxALIGN_RIGHT).Border(wxTOP|wxRIGHT);
  entry_flag.Align(wxALIGN_LEFT);
  label_slide.Align(wxALIGN_CENTER_VERTICAL);


  wxFlexGridSizer *asizer = new wxFlexGridSizer(2);

  asizer->Add(new wxStaticText(panel,wxID_ANY,"Count"),
	      label_flag.DoubleBorder(wxLEFT));
  asizer->Add(naper_spin,entry_flag);

  wxBoxSizer *spysizer = new wxBoxSizer(wxVERTICAL);
  spysizer->Add(spiral_radio);
  spysizer->Add(equidistant_radio);

  asizer->Add(new wxStaticText(panel,wxID_ANY,"Radii"),label_radio);
  asizer->Add(spysizer,entry_flag.Border(wxTOP));


  wxStaticBoxSizer *aper_sizer = new wxStaticBoxSizer(wxVERTICAL,panel,
						      " Apertures ");
  aper_sizer->Add(asizer,wxSizerFlags().Border().Center());

  wxFlexGridSizer *rsizer = new wxFlexGridSizer(2);

  rsizer->Add(new wxStaticText(panel,wxID_ANY,"Inner"),label_flag);
  rsizer->Add(rmin_spin,entry_flag);
  rsizer->Add(new wxStaticText(panel,wxID_ANY,"Outer"),label_flag);
  rsizer->Add(rmax_spin,entry_flag);

  wxStaticBoxSizer *ring_sizer = new wxStaticBoxSizer(wxVERTICAL,panel,
						      " Background ring ");
  ring_sizer->Add(rsizer,wxSizerFlags().Border().Center());

  wxBoxSizer *snap_sizer = new wxBoxSizer(wxHORIZONTAL);
  snap_sizer->Add(snap_check);

  wxBoxSizer *zoom_sizer = new wxBoxSizer(wxHORIZONTAL);
  zoom_sizer->Add(new wxStaticText(panel,wxID_ANY,L"⊝"),label_slide);
  zoom_sizer->Add(zoom_slider,wxSizerFlags(1));
  zoom_sizer->Add(new wxStaticText(panel,wxID_ANY,L"⊕"),label_slide);

  wxBoxSizer *csizer = new wxBoxSizer(wxVERTICAL);
  csizer->Add(aper_sizer,wxSizerFlags().Expand());
  csizer->Add(ring_sizer,wxSizerFlags().Expand().DoubleBorder(wxTOP|wxBOTTOM));
  csizer->Add(snap_sizer,wxSizerFlags().Centre());
  csizer->AddStretchSpacer();
  csizer->Add(zoom_sizer,wxSizerFlags().Border(wxLEFT|wxRIGHT).Expand());


  wxBoxSizer *topsizer = new wxBoxSizer(wxHORIZONTAL);
  topsizer->Add(canvas,wxSizerFlags(1).Expand().Border(wxRIGHT));
  topsizer->Add(csizer,wxSizerFlags().Expand());


  panel->SetSizer(topsizer);
  topsizer->SetSizeHints(panel);
  SetPanel(panel,wxSizerFlags(1).Expand().Border());
  Fit();

  Bind(wxEVT_TOOL,&MuniAphot::OnSave,this,wxID_SAVE);
  Bind(wxEVT_TIMER,&MuniAphot::OnTimer,this);
  Bind(wxEVT_IDLE,&MuniAphot::OnIdle,this);
  Bind(wxEVT_END_PROCESS,&MuniAphot::OnFinish,this);
  Bind(EVT_CLICK,&MuniAphot::OnClick,this);

  Bind(wxEVT_CHECKBOX,&MuniAphot::OnSnap,this,snap_check->GetId());
  Bind(wxEVT_SLIDER,&MuniAphot::OnZoom,this,zoom_slider->GetId());
  Bind(wxEVT_SPINCTRL,&MuniAphot::OnNaper,this,naper_spin->GetId());
  Bind(wxEVT_SPINCTRL,&MuniAphot::OnRingMin,this,rmin_spin->GetId());
  Bind(wxEVT_SPINCTRL,&MuniAphot::OnRingMax,this,rmax_spin->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniAphot::OnSpiral,this,ID_APHOT_SPIRAL);
  Bind(wxEVT_RADIOBUTTON,&MuniAphot::OnSpiral,this,ID_APHOT_EQL);
  Bind(wxEVT_UPDATE_UI,&MuniAphot::OnUpdateRing,this,rmin_spin->GetId());
  Bind(wxEVT_UPDATE_UI,&MuniAphot::OnUpdateRing,this,rmax_spin->GetId());

  FitsHdu hdu = fits.FindHdu(FINDEXTNAME);
  if( hdu.IsOk() ) {
    FitsTable table(hdu);
    int i = int(table.Nrows() * float(rand()) / float(RAND_MAX));
    canvas->SetPosition(table.Cell(i,0),table.Cell(i,1));
    canvas->SetFwhm(hdu.GetKeyDouble("FWHM"));
  }
  else
    SetStatusWarning("Stars find not done yet.");
}

MuniAphot::~MuniAphot()
{
  if( ! tmpfits.IsEmpty() && wxFileExists(tmpfits) )
    wxRemoveFile(tmpfits);

  config->aphot_naper = naper;
  config->aphot_rmin = ring_min;
  config->aphot_rmax = ring_max;
  config->aphot_zoom = zoom;
  config->aphot_spiral = spiral;
  config->aphot_snap = snap;
  config->aphot_saper = canvas->GetSaper();
  config->aphot_size = GetSize();
}

void MuniAphot::OnNaper(wxSpinEvent& event)
{
  naper = event.GetInt();
  canvas->SetNaper(naper);
  CleanDraw();
}

void MuniAphot::OnRingMin(wxSpinEvent& event)
{
  ring_min = event.GetInt();
  canvas->SetRingMin(ring_min);
  CleanDraw();
}

void MuniAphot::OnRingMax(wxSpinEvent& event)
{
  ring_max = event.GetInt();
  canvas->SetRingMax(ring_max);
  CleanDraw();
}

void MuniAphot::OnUpdateRing(wxUpdateUIEvent& event)
{
  wxWindow *w =  FindWindow(event.GetId());
  wxASSERT(w);

  if( ring_max < ring_min )
      w->SetForegroundColour(*wxRED);
  else
      w->SetForegroundColour(wxNullColour);
}

void MuniAphot::OnSpiral(wxCommandEvent& event)
{
  spiral = event.GetId() == ID_APHOT_SPIRAL;
  canvas->SetSpiral(spiral);
  CleanDraw();
}


void MuniAphot::OnIdle(wxIdleEvent& event)
{
  if( tmpfits.IsEmpty() ) {
    tmpfits = wxFileName::CreateTempFileName("xmunipack-aphot_");
    wxLogDebug("Making a copy of data..."+fitsname+" to "+tmpfits);
    wxCopyFile(fitsname,tmpfits,true);
  }
}

void MuniAphot::OnTimer(wxTimerEvent& event)
{
  wxArrayString out(GetLastOutput());
  if( ! (out.GetCount() > 0) )
    return;

  DrawStars(out);
  wxString line = Parser(out.Last(),"APHOT");
  wxLogDebug(line);
  SetStatus(line);
}


void MuniAphot::OnFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniAphot::OnFinish");

  if( event.GetExitCode() == 0 ) {

    wxArrayString out(GetLastOutput());
    if( out.GetCount() > 0 ) {
      DrawStars(out);
      wxString line = Parser(out.Last(),"APHOT");
      wxLogDebug(line);
      SetStatus(line);
    }
    SetStatus("");
  }
  else {
    wxLogDebug("Removing "+tmpfits);
    if( ! tmpfits.IsEmpty() && wxFileExists(tmpfits) ) {
      wxRemoveFile(tmpfits);
      tmpfits = "";
    }
  }

}

void MuniAphot::OnOutput(const wxArrayString& out)
{
  if( out.GetCount() > 0 ) {
    DrawStars(out);
    wxString line = Parser(out.Last(),"APHOT");
    wxLogDebug(line);
    SetStatus(line);
  }
}

void MuniAphot::OnInput(MuniProcess *mproc)
{
  mproc->Write("NFILES = 1");
  mproc->Write("FILE = '"+tmpfits+"' '"+tmpfits+"'");
}

void MuniAphot::OnSave(wxCommandEvent& event)
{
  wxLogDebug("MuniAphot::OnSave");

  if( ! FitsCopyHdu(tmpfits,fitsname,APEREXTNAME) ) {
    wxLogError("Failed to save results to " + fitsname +
	       " (space, access, ... ?).");
    return;
  }

  event.Skip(); // it does to propagate this event to MuniTool
}

void MuniAphot::OnClick(MuniClickEvent& event)
{
  wxLogDebug("MuniAphot::OnClick %d %d",event.x, event.y);

  if( event.GetId() == canvas->GetId() ) {
    SetStatus("Aperture radius is %.1f pixels.",event.r);
    return;
  }

  if( ! (0 <= event.x && event.x < image.GetWidth() &&
	 0 <= event.y && event.y < image.GetHeight()) ) return;

  SetStatus("Aperture centre at %d,%d.",event.x,event.y);

  ByHand(event.x,event.y);

  wxQueueEvent(canvas,event.Clone());

  // Prepare Draw for main display
  vector<wxObject *> objects;
  MuniApertures aper = canvas->GetAper();

  wxColour cring(138,184,230,96);
  objects.push_back(new MuniDrawPen(wxPen(cring,1)));
  objects.push_back(new MuniDrawBrush(cring));
  objects.push_back(new MuniDrawRing(event.x,event.y,ring_min,ring_max));
  objects.push_back(new MuniDrawBrush(*wxTRANSPARENT_BRUSH));
  wxColour caper(138,184,230,196);
  objects.push_back(new MuniDrawPen(wxPen(caper,1)));
  for(int i = 1; i <= aper.GetNaper(); i++)
    objects.push_back(new MuniDrawCircle(event.x,event.y,aper.GetAper(i)));

  MuniLayer layer(ID_PHOTOMETRY,objects);
  MuniDrawEvent ev(EVT_DRAW);
  ev.layer = layer;
  wxQueueEvent(GetParent(),ev.Clone());
}

void MuniAphot::ByHand(int x, int y)
{
  if( aproc ) return;

  wxLogDebug("MuniAphot::ByHand %d %d",x,y);
  wxASSERT(aproc == 0);

  aproc = new MuniProcess(this,"aphot");
  aproc->SetEcho(false);
  aproc->Write("PIPELOG = T");
  aproc->Write("NAPER = 1");
  aproc->Write("APER = %f",canvas->GetRaper());
  aproc->Write("NSTARS = 1");
  aproc->Write("STAR = %d %d",x,y);
  aproc->Write("SNAP = %c",snap ? 'T' : 'F');
  aproc->Write("FILE = '"+tmpfits+"' '"+tmpfits+"'");
  aproc->OnStart();

  Bind(wxEVT_END_PROCESS,&MuniAphot::OnHandFinish,this);
}

void MuniAphot::OnHandFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniAphot::OnHandFinish");

  if( aproc->GetExitCode() == 0 ) {
    wxArrayString out(aproc->GetOutput());
    for(size_t i = 0; i < out.GetCount(); i++) {
      wxString line = Parser(out[i],"BYHAND");
      double x,y,b,berr,f,ferr;
      if( !line.IsEmpty() && HandParser(line,&x,&y,&b,&berr,&f,&ferr)  ) {
	//wxLogDebug("Counts %f %f; sky %f %f",f,ferr,b,berr);
	if( f >= 0 )
	  SetStatusDisplay("Counts <b>"+HumanFormat(f)+"</b> (%.3f); Sky <b>"+
		   HumanFormat(b)+"</b>",ferr / f);
	else
	  SetStatusDisplay("Counts indeterminable; Sky <b>"
			   +HumanFormat(b)+"</b>");
	if( snap )
	  canvas->SetPosition(int(x),int(y));
      }
    }
  }
  else {
    wxArrayString out(aproc->GetOutput());
    wxArrayString err(aproc->GetErrors());
    wxLogDebug("ByHand Finish failed:");
    for(size_t i = 0; i < out.GetCount(); i++)
      wxLogDebug(out[i]);
    for(size_t i = 0; i < err.GetCount(); i++)
      wxLogDebug(err[i]);
  }

  delete aproc;
  aproc = 0;
  Unbind(wxEVT_END_PROCESS,&MuniAphot::OnHandFinish,this);
}


bool MuniAphot::StarParser(const wxString& line, double *x, double *y,
			   double *f) const
{
  if( line.IsEmpty() )
    return false;

  wxStringTokenizer tk(line);
  int m = 0;
  while( tk.HasMoreTokens() ) {
    wxString a = tk.GetNextToken();
    m++;
    if( m == 1 && ! a.ToDouble(x) )
      return false;
    if( m == 2 && ! a.ToDouble(y) )
      return false;
    if( m == 3 && ! a.ToDouble(f) )
      return false;
  }
  return true;
}

bool MuniAphot::HandParser(const wxString& line, double *x, double *y,
			   double *b,double *berr,double *f,double *ferr) const
{
  if( line.IsEmpty() )
    return false;

  wxStringTokenizer tk(line);
  int m = 0;
  while( tk.HasMoreTokens() ) {
    wxString a = tk.GetNextToken();
    m++;
    if( m == 2 && ! a.ToDouble(x) )
      return false;
    if( m == 3 && ! a.ToDouble(y) )
      return false;
    if( m == 4 && ! a.ToDouble(b) )
      return false;
    if( m == 5 && ! a.ToDouble(berr) )
      return false;
    if( m == 6 && ! a.ToDouble(f) )
      return false;
    if( m == 7 && ! a.ToDouble(ferr) )
      return false;
  }
  return true;
}


void MuniAphot::DrawStars(const wxArrayString& out) const
{
  vector<double> xcoo,ycoo,flux;

  for(size_t i = 0; i < out.GetCount(); i++) {

    wxString line = Parser(out[i],"APHOT");
    double x,y,f;
    if( StarParser(line,&x,&y,&f) ) {
      xcoo.push_back(x);
      ycoo.push_back(y);
      flux.push_back(f);
    }
  }

  MuniStarLayer layer;
  layer.DrawObjects(xcoo,ycoo,flux);

  MuniDrawEvent ev(EVT_DRAW);
  ev.layer = layer.GetLayer();
  wxQueueEvent(GetParent(),ev.Clone());
}

void MuniAphot::CleanDraw() const
{
  MuniLayer layer(ID_PHOTOMETRY);
  MuniDrawEvent ev(EVT_DRAW);
  ev.layer = layer;
  wxQueueEvent(GetParent(),ev.Clone());
}

void MuniAphot::OnSnap(wxCommandEvent& event)
{
  snap = event.IsChecked();
}

void MuniAphot::OnZoom(wxCommandEvent& event)
{
  zoom = event.GetInt();
  wxQueueEvent(canvas,event.Clone());
  SetStatus(L"Zoom: %d×",int(zoom+0.5));
}
