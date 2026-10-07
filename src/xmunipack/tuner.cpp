/*

  xmunipack - mini-panel for fast intensity scale tunnings

  Copyright © 2021 F.Hroch (hroch@physics.muni.cz)

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

#include "tuner.h"
#include <wx/wx.h>

using namespace std;


MuniTuner::MuniTuner(wxWindow *w, wxWindowID id): wxPanel(w,id)
{
  double rmin, rmax;
  int steps;

  steps = 11;
  rmin = 1.0;
  rmax = 0.0;
  adjblack = new MuniTunerAdjuster(this,ID_TONE_QBLACK,rmin,rmax,steps,L"◑",L"○");
  steps = 25;
  rmin = 1e-3;
  rmax = 1e3;
  adjsense = new MuniTunerLogjuster(this,ID_TONE_RSENSE,rmin,rmax,steps,L"☾",L"☼");

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(adjsense,wxSizerFlags().Border(wxLEFT|wxRIGHT).Expand());
  topsizer->Add(adjblack,wxSizerFlags().Border(wxLEFT|wxRIGHT).Expand());
  SetSizer(topsizer);

  adjblack->Enable(false);
  adjsense->Enable(false);
}

void MuniTuner::SetTone(double qblack, double rsense)
{
  adjblack->Enable(true);
  adjsense->Enable(true);
  adjblack->SetValue(qblack);
  adjsense->SetValue(rsense);
}

// --------------------------------------------------------------------

MuniTunerAdjuster::MuniTunerAdjuster(wxWindow *w, wxWindowID id,
				     double xmin, double xmax, int steps,
				     const wxString& lmin, const wxString& lmax):
  wxWindow(w,id), irange(steps), rmin(xmin), rmax(xmax)
{
  slider = new wxSlider(this,wxID_ANY,0,0,irange);

  wxSizerFlags flaglabel; flaglabel.Align(wxALIGN_CENTER_VERTICAL);
  wxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);
  sizer->Add(new wxStaticText(this,wxID_ANY,lmin),flaglabel);
  sizer->Add(slider,wxSizerFlags(1).Align(wxALIGN_CENTER_VERTICAL));
  sizer->Add(new wxStaticText(this,wxID_ANY,lmax),flaglabel);
  SetSizer(sizer);

  int sid = slider->GetId();

  Bind(wxEVT_SCROLL_THUMBTRACK,&MuniTunerAdjuster::OnScroll,this,sid);
  Bind(wxEVT_SCROLL_LINEUP,&MuniTunerAdjuster::OnScroll,this,sid);
  Bind(wxEVT_SCROLL_LINEDOWN,&MuniTunerAdjuster::OnScroll,this,sid);
  Bind(wxEVT_SCROLL_PAGEUP,&MuniTunerAdjuster::OnScroll,this,sid);
  Bind(wxEVT_SCROLL_PAGEDOWN,&MuniTunerAdjuster::OnScroll,this,sid);
}

void MuniTunerAdjuster::SetValue(double r)
{
  double k = (r - rmin) / (rmax - rmin);
  slider->SetValue(int(irange*k+0.5));
}

double MuniTunerAdjuster::GetValue() const
{
  return GetValue(slider->GetValue());
}

double MuniTunerAdjuster::GetValue(int n) const
{
  return rmin + (double(n)/double(irange))*(rmax - rmin);
}

void MuniTunerAdjuster::OnScroll(wxScrollEvent& event)
{
  //  wxLogDebug("MuniTunerAdjuster::OnScroll %d %f",event.GetPosition(),
  //  	     GetValue(event.GetPosition()));

  MuniTuneEvent e(EVT_TUNE,GetId());
  e.SetEventObject(GetParent());
  e.x = GetValue(event.GetPosition());
  wxQueueEvent(GetGrandParent(),e.Clone());
}

MuniTunerLogjuster::MuniTunerLogjuster(wxWindow *w, wxWindowID id,
				       double minValue, double maxValue, int steps,
				       const wxString& cmin, const wxString& cmax):
  wxWindow(w,id), scale(1.0), rmin(minValue), rmax(maxValue)
{
  wxASSERT(rmin > 0 && rmax > 0 && rmax > rmin && steps > 0);

  scale = steps / log(rmax/rmin);
  slider = new wxSlider(this,wxID_ANY,0,0,steps);

  wxSizerFlags flaglabel; flaglabel.Align(wxALIGN_CENTER_VERTICAL);
  wxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);
  sizer->Add(new wxStaticText(this,wxID_ANY,cmin),flaglabel);
  sizer->Add(slider,wxSizerFlags(1).Align(wxALIGN_CENTER_VERTICAL));
  sizer->Add(new wxStaticText(this,wxID_ANY,cmax),flaglabel);
  SetSizer(sizer);

  int sid = slider->GetId();

  Bind(wxEVT_SCROLL_THUMBTRACK,&MuniTunerLogjuster::OnScroll,this,sid);
  Bind(wxEVT_SCROLL_LINEUP,&MuniTunerLogjuster::OnScroll,this,sid);
  Bind(wxEVT_SCROLL_LINEDOWN,&MuniTunerLogjuster::OnScroll,this,sid);
  Bind(wxEVT_SCROLL_PAGEUP,&MuniTunerLogjuster::OnScroll,this,sid);
  Bind(wxEVT_SCROLL_PAGEDOWN,&MuniTunerLogjuster::OnScroll,this,sid);
}


void MuniTunerLogjuster::SetValue(double r)
{
  double k = scale*log(r/rmin);
  slider->SetValue(int(k+0.5));
}

double MuniTunerLogjuster::GetValue() const
{
  return GetValue(slider->GetValue());
}

double MuniTunerLogjuster::GetValue(int n) const
{
  double x = n / scale;
  return rmin*exp(x);
}


void MuniTunerLogjuster::OnScroll(wxScrollEvent& event)
{
  //  wxLogDebug("MuniTunerLogjusterBase::OnScroll %d %f",int(event.GetPosition()),
  //  	     GetValue(event.GetPosition()));

  MuniTuneEvent e(EVT_TUNE,GetId());
  e.x = GetValue(event.GetPosition());
  wxQueueEvent(GetParent(),e.Clone());
}
