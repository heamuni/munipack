/*

  xmunipack - figure caption

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

#include "caption.h"
#include <wx/wx.h>

using namespace std;

#define OBJECT "Object"
#define FILTER "Filter"
#define CSPACE "Colours"
#define DATE "Date"
#define TIME "Time"
#define EXPTIME "Exposure"
#define COUNTS "Counts"
#define INTENSITY "I"
#define MAGNITUDE L"μ"
#define ALPHA L"α"
#define DELTA L"δ"

static const char *cspace_label_XYZ[] = { "X","Y","Z" };

static wxString LABEL(wxString a)
{
  if( a != "" )
    return "<span fgcolor=\"gray\">"  + a +  "</span>";
  else
    return "";
}

MuniDisplayCaptionInfo::MuniDisplayCaptionInfo(wxWindow *w, MuniConfig *c):
  wxPanel(w,wxID_ANY), config(c), init(false)
{
  wxSizerFlags label_flags, value_flags;
  label_flags.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).DoubleBorder(wxRIGHT);
  value_flags.Align(wxALIGN_CENTER_VERTICAL).Expand();

  wxStaticText *l;

  // Object
  label_object = new wxStaticText(this, wxID_ANY,"");
  label_object->SetLabelMarkup(LABEL(OBJECT));
  object = new wxStaticText(this, wxID_ANY, "Object");

  wxFlexGridSizer *objectsizer = new wxFlexGridSizer(2);
  objectsizer->Add(label_object,label_flags);
  objectsizer->Add(object,value_flags);

  // Exposure time
  l = new wxStaticText(this, wxID_ANY,"");
  l->SetLabelMarkup(LABEL(EXPTIME));
  exptime = new wxStaticText(this, wxID_ANY, "Exptime [s]");

  wxFlexGridSizer *filtersizer = new wxFlexGridSizer(2);
  filtersizer->Add(l,label_flags);
  filtersizer->Add(exptime,value_flags);

  // Filter; the filter is
  label_colour = new wxStaticText(this, wxID_ANY,"");
  label_colour->SetLabelMarkup(LABEL(FILTER));
  colour = new wxStaticText(this, wxID_ANY, "Filter");
  filtersizer->Add(label_colour,label_flags);
  filtersizer->Add(colour,value_flags);

  wxFlexGridSizer *datesizer = new wxFlexGridSizer(2);
  l = new wxStaticText(this, wxID_ANY,"");
  l->SetLabelMarkup(LABEL(DATE));
  date = new wxStaticText(this, wxID_ANY, "Date");
  datesizer->Add(l,label_flags);
  datesizer->Add(date,value_flags);

  l = new wxStaticText(this, wxID_ANY,"");
  l->SetLabelMarkup(LABEL(TIME));
  time = new wxStaticText(this, wxID_ANY, "Time");
  datesizer->Add(l,label_flags);
  datesizer->Add(time,value_flags);

  wxBoxSizer *topsizer = new wxBoxSizer(wxHORIZONTAL);
  topsizer->AddStretchSpacer(13);
  topsizer->Add(objectsizer,wxSizerFlags().Border());
  topsizer->AddStretchSpacer();
  topsizer->Add(filtersizer,wxSizerFlags().Border());
  topsizer->AddStretchSpacer();
  topsizer->Add(datesizer,wxSizerFlags().Border());
  topsizer->AddStretchSpacer(13);
  SetSizer(topsizer);

  Bind(wxEVT_IDLE,&MuniDisplayCaptionInfo::OnIdle,this);
}

void MuniDisplayCaptionInfo::SetArray(const FitsArray& a)
{
  wxASSERT(a.IsOk());
  array = a;
  init = true;
}

void MuniDisplayCaptionInfo::OnIdle(wxIdleEvent& event)
{
  /*
    All the info is filled during the idle time.
    If FITS keyword for object is missing, the label is hidden.
    All other labels (not values) are visible, even if they are
    missing to inform user by subliminal way.
  */

  if( init ) {
    init = false;

    wxString obj(array.GetKey(config->fits_key_object));

    // to format of exposure time
    wxString expline(array.GetKey(config->fits_key_exptime));
    double e;
    if( expline.ToDouble(&e) ) {
      wxString line;
      if( e > 0.1 ) {
	if( abs(round(e) - e) < 0.05 )
	  line.Printf("%.0f s",e);
	else
	  line.Printf("%.1f s",e);
      }
      else
	line.Printf("1/%.0f s",1.0/e);
      expline = line;
    }

    // DATE-OBS split
    FitsTime ft(array.GetKey(config->fits_key_dateobs));

    // Filter or Colour-space
    wxString clabel;
    if( array.IsColour() ) {
      label_colour->SetLabelMarkup(LABEL(CSPACE));
      clabel = "CIE 1931 XYZ";
    }
    else
      clabel = array.GetKey(config->fits_key_filter);

    label_object->SetLabelMarkup(obj == "" ? "" : LABEL(OBJECT));
    object->SetLabel(obj);
    exptime->SetLabel(expline);
    colour->SetLabel(clabel);
    date->SetLabel(ft.GetDate());
    time->SetLabel(ft.GetTime());

    Layout();
  }
}

void MuniDisplayCaptionInfo::ConfigUpdate()
{
  init = true;
  Refresh();
}

MuniDisplayCaptionMotion::MuniDisplayCaptionMotion(wxWindow *w, MuniConfig *c):
  wxPanel(w,wxID_ANY), config(c), init(false), update(false), hascal(false),
  show_wcs(false),x(0),y(0)
{
  wxSizerFlags label_flags, value_flags, number_flags, units_flags;
  label_flags.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).DoubleBorder(wxRIGHT);
  value_flags.Align(wxALIGN_CENTER_VERTICAL).Expand().FixedMinSize();
  number_flags.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).FixedMinSize();
  units_flags.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT);

  wxStaticText *l;

  // counts or photons
  label_quantity = new wxStaticText(this, wxID_ANY,"");
  label_quantity->SetLabelMarkup(LABEL(COUNTS));
  quantity = new wxStaticText(this, wxID_ANY, "_9.999M_");

  quantitysizer = new wxFlexGridSizer(2);
  quantitysizer->Add(label_quantity,label_flags);
  quantitysizer->Add(quantity,value_flags);
  quantitysizer->AddStretchSpacer();

  // intensity
  wxStaticText *label_intensity = new wxStaticText(this, wxID_ANY,"");
  label_intensity->SetLabelMarkup(LABEL(INTENSITY));
  intensity = new wxStaticText(this, wxID_ANY, "9.999M");

  // magnitude
  wxStaticText *label_magnitude = new wxStaticText(this, wxID_ANY,"");
  label_magnitude->SetLabelMarkup(LABEL(MAGNITUDE));
  magnitude = new wxStaticText(this, wxID_ANY,"99.99");

  qsizer = new wxFlexGridSizer(3);
  qsizer->Add(label_intensity,label_flags);
  qsizer->Add(intensity,value_flags);
  l = new wxStaticText(this, wxID_ANY,"");
  l->SetLabelMarkup(LABEL(L"eV/s/m²/\"²"));
  qsizer->Add(l,units_flags);
  qsizer->Add(label_magnitude,label_flags);
  qsizer->Add(magnitude,value_flags);
  l = new wxStaticText(this, wxID_ANY,"");
  l->SetLabelMarkup(LABEL(L"mag/\"²"));
  qsizer->Add(l,units_flags);

  // pixels coordinates
  pixsizer = new wxFlexGridSizer(2);

  l = new wxStaticText(this, wxID_ANY,"");
  l->SetLabelMarkup(LABEL("x"));
  xpix = new wxStaticText(this, wxID_ANY, "99999");
  pixsizer->Add(l,label_flags);
  pixsizer->Add(xpix,number_flags);

  l = new wxStaticText(this, wxID_ANY,"");
  l->SetLabelMarkup(LABEL("y"));
  ypix = new wxStaticText(this, wxID_ANY, "99999");
  pixsizer->Add(l,label_flags);
  pixsizer->Add(ypix,number_flags);

  // spherical coordinates
  coosizer = new wxFlexGridSizer(2);

  l = new wxStaticText(this, wxID_ANY,"");
  l->SetLabelMarkup(LABEL(ALPHA));
  alpha = new wxStaticText(this, wxID_ANY, "_11:22:333.999");
  coosizer->Add(l,label_flags);
  coosizer->Add(alpha,number_flags);

  l = new wxStaticText(this, wxID_ANY,"");
  l->SetLabelMarkup(LABEL(DELTA));
  delta = new wxStaticText(this, wxID_ANY, "_11:22:333.999");
  coosizer->Add(l,label_flags);
  coosizer->Add(delta,number_flags);

  // XYZ colours
  wxStaticText *labels_XYZ[3];
  for(size_t i = 0; i < 3; i++) {
    labels_XYZ[i] = new wxStaticText(this, wxID_ANY,"");
    labels_XYZ[i]->SetLabelMarkup(LABEL(cspace_label_XYZ[i]));
    colours_XYZ[i] = new wxStaticText(this, wxID_ANY, "__99999.9");
  }
  colours_XYZ[0]->SetForegroundColour(wxColour(92,0,0));
  colours_XYZ[1]->SetForegroundColour(wxColour(0,92,0));
  colours_XYZ[2]->SetForegroundColour(wxColour(0,0,2*92));

  coloursizer = new wxFlexGridSizer(6);
  for(size_t i = 0; i < 3; i++) {
    coloursizer->Add(labels_XYZ[i],label_flags);
    coloursizer->Add(colours_XYZ[i],number_flags);
  }

  // final packing
  wxBoxSizer *topsizer = new wxBoxSizer(wxHORIZONTAL);
  topsizer->AddStretchSpacer(13);
  topsizer->Add(quantitysizer,wxSizerFlags().Border());
  topsizer->Add(qsizer,wxSizerFlags().Border());
  topsizer->Add(coloursizer,wxSizerFlags().Border());
  topsizer->AddStretchSpacer(13);
  topsizer->Add(pixsizer,wxSizerFlags().Border());
  topsizer->AddStretchSpacer(1);
  topsizer->Add(coosizer,wxSizerFlags().Border());
  topsizer->AddStretchSpacer(13);
  SetSizer(topsizer);

  Bind(wxEVT_IDLE,&MuniDisplayCaptionMotion::OnIdle,this);
  Bind(EVT_SLEW,&MuniDisplayCaptionMotion::OnMouseMotion,this);
}

void MuniDisplayCaptionMotion::SetArray(const FitsArray& a)
{
  wxASSERT(a.IsOk());
  array = a;
  x = array.GetWidth() / 2;
  y = array.GetHeight() / 2;
  init = true;

  InitArray();
}

void MuniDisplayCaptionMotion::InitArray()
{
  value = FitsValue(array,config->phsystemfile, config->fits_key_area,
		    config->fits_key_exptime, config->fits_key_filter);
  coords = FitsCoo(array);

  // spherical coordinates are visible only when the FITS is calibrated
  show_wcs = coords.HasWCS();
  coords.SetType(show_wcs ? config->display_cootype : COO_PIXEL);

  // photometry calibration
  hascal = value.HasPhcal();

  // select label for counts or photons
  if( hascal )
    label_quantity->SetLabelMarkup(LABEL(value.GetName()));
}


void MuniDisplayCaptionMotion::OnIdle(wxIdleEvent& event)
{
  if( init ) {
    init = false;

    coosizer->Show(show_wcs);
    qsizer->Show(hascal);

    bool iscolour = array.IsColour();
    quantitysizer->Show(!iscolour);
    coloursizer->Show(iscolour);
    Layout();
  }

  if( update ) {
    update = false;

    if( array.IsColour() ) {

      vector<wxString> labels = value.Get_str(x,y,{0,1,2});
      wxASSERT(labels.size() == 3);

      for(size_t i = 0; i < 3; i++)
	colours_XYZ[i]->SetLabel(labels[i]);

    }
    else {
      if( hascal ) {
	value.SetType(UNIT_PHOTON);
	label_quantity->SetLabelMarkup(LABEL(value.GetName()));
	quantity->SetLabel(value.Get_str(x,y));
	value.SetType(UNIT_INTENSITY);
	intensity->SetLabelMarkup(value.Get_str(x,y));
	value.SetType(UNIT_MAG);
	magnitude->SetLabelMarkup(value.Get_str(x,y));
      }
      else {
	label_quantity->SetLabelMarkup(LABEL(value.GetName()));
	value.SetType(UNIT_COUNT);
	quantity->SetLabel(value.Get_str(x,y));
      }
    }

    wxString xstr,ystr;
    coords.GetPix(x,y,xstr,ystr);
    xpix->SetLabel(xstr);
    ypix->SetLabel(ystr);

    if( show_wcs ) {
      wxString a,d;
      coords.GetCoo(x,y,a,d);
      alpha->SetLabel(a);
      delta->SetLabel(d);
    }

    Layout();
  }
}

void MuniDisplayCaptionMotion::OnMouseMotion(MuniSlewEvent& event)
{
  update = true;
  x = event.x;
  y = event.y;
}

void MuniDisplayCaptionMotion::ConfigUpdate()
{
  InitArray();
  update = true;
  Refresh();
}

//---- MuniDisplayCaption ----------------------------------------

MuniDisplayCaption::MuniDisplayCaption(wxWindow *w, MuniConfig *config):
  wxPanel(w,wxID_ANY), inside(false)
{
  info = new MuniDisplayCaptionInfo(this,config);
  motion = new MuniDisplayCaptionMotion(this,config);

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(info,wxSizerFlags().Expand());
  SetSizer(topsizer);
  motion->Show(false);

  Bind(EVT_SLEW,&MuniDisplayCaption::OnMouseMotion,this);
}

void MuniDisplayCaption::SetArray(const FitsArray& array)
{
  wxASSERT(array.IsOk());
  info->SetArray(array);
  motion->SetArray(array);
}


void MuniDisplayCaption::OnMouseMotion(MuniSlewEvent& event)
{
  //  wxLogDebug("MuniDisplayCaption::OnMouseMotion %d %d",event.x,event.y);

  if( event.entering ) {
    inside = true;
    //    wxLogDebug("Hidding spanish inquisition...");
    info->Show(false);
    motion->Show(true);
    GetSizer()->Replace(info,motion);
    Layout();
  }

  if( event.leaving ) {
    inside = false;
    //    wxLogDebug("Revealing spanish inquisition...");
    info->Show(true);
    motion->Show(false);
    GetSizer()->Replace(motion,info);
    Layout();
  }

  if( inside )
    wxQueueEvent(motion,event.Clone());
}

void MuniDisplayCaption::ConfigUpdate()
{
  info->ConfigUpdate();
  motion->ConfigUpdate();
  Refresh();
}
