/*

  XMunipack - derived events

  Copyright © 2009-2012, 2017-22 F.Hroch (hroch@physics.muni.cz)

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


#ifndef _XMUNIPACK_EVENT_H_
#define _XMUNIPACK_EVENT_H_

#include "types.h"
#include "fits.h"

/*

  Definition of new events. Adds some data structures to standard ones.

*/


wxDECLARE_EVENT(EVT_FILELOAD, wxCommandEvent);
wxDECLARE_EVENT(EVT_FINISH_DIALOG, wxCommandEvent);
wxDECLARE_EVENT(EVT_TOOL_FINISH, wxCommandEvent);
wxDECLARE_EVENT(EVT_LIST_ACTIVATED, wxNotifyEvent);


// -- Slew ---
// positions of zoomed areas

class MuniSlewEvent: public wxEvent
{
 public:
 MuniSlewEvent(wxEventType eventType =wxEVT_NULL, int id =wxID_ANY):
   wxEvent(id, eventType),x(0),y(0),xsub(0),ysub(0),inside(false),leaving(false),
   entering(false), zoom(0.0) {}
  virtual wxEvent *Clone(void) const { return new MuniSlewEvent(*this); }
  int x,y,xsub,ysub;
  bool inside, leaving, entering;
  double zoom;
};

wxDECLARE_EVENT(EVT_SLEW,MuniSlewEvent);

#define MuniSlewEventHandler(func) (&func)

// -- Mouse click (replaces mouse events, ones are not propagated.)

class MuniClickEvent: public wxCommandEvent
{
 public:
 MuniClickEvent(wxEventType eventType =wxEVT_NULL, int id =wxID_ANY):
  wxCommandEvent(eventType,id),x(0),y(0),r(0.0) {}
  // MuniClickEvent(wxCommandEvent& e): wxCommandEvent(e) { }
  virtual wxEvent *Clone(void) const { return new MuniClickEvent(*this); }
  int x,y;
  double r;
};

wxDECLARE_EVENT(EVT_CLICK,MuniClickEvent);

#define MuniClickEventHandler(func) (&func)

// -- Zoom

class MuniZoomEvent: public wxCommandEvent
{
 public:
  MuniZoomEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(e,id),x(0),y(0),zoom(1.0),render(false) {}
  MuniZoomEvent(wxCommandEvent& e):
    wxCommandEvent(e),x(0),y(0),zoom(1.0),render(false) {}
  wxEvent *Clone(void) const { return new MuniZoomEvent(*this); }
  double x,y;
  double zoom;
  bool render;
};

wxDECLARE_EVENT(EVT_ZOOM,MuniZoomEvent);

#define MuniZoomEventHandler(func) (&func)



// -- Tune

class MuniTuneEvent: public wxCommandEvent
{
 public:
  MuniTuneEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(e,id),x(0.0),n(0),index(0) {}
  MuniTuneEvent(wxCommandEvent& e): wxCommandEvent(e),x(0.0),n(0),index(0) {}
  wxEvent *Clone(void) const { return new MuniTuneEvent(*this); }
  double x;
  int n;
  int index;
};

wxDECLARE_EVENT(EVT_TUNE,MuniTuneEvent);

#define MuniTuneEventHandler(func) (&func)

// -- Tune scale
class MuniTuneScaleEvent: public wxCommandEvent
{
 public:
  MuniTuneScaleEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(e,id),black(0.0),sense(1.0),reset(false) {}
  MuniTuneScaleEvent(wxCommandEvent& e):
    wxCommandEvent(e),black(0.0),sense(1.0),reset(false) {}
  wxEvent *Clone(void) const { return new MuniTuneScaleEvent(*this); }
  double black, sense;
  bool reset;
};

wxDECLARE_EVENT(EVT_TUNE_SCALE,MuniTuneScaleEvent);

#define MuniTuneScaleEventHandler(func) (&func)

// -- Tune itt
class MuniTuneIttEvent: public wxCommandEvent
{
 public:
  MuniTuneIttEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(e,id), itt(ITT_LINE) {}
  MuniTuneIttEvent(wxCommandEvent& e): wxCommandEvent(e), itt(ITT_LINE) {}
  wxEvent *Clone(void) const { return new MuniTuneIttEvent(*this); }
  int itt;
};

wxDECLARE_EVENT(EVT_TUNE_ITT,MuniTuneIttEvent);

#define MuniTuneIttEventHandler(func) (&func)


// -- Tune Palette
class MuniTunePalEvent: public wxCommandEvent
{
 public:
  MuniTunePalEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(e,id), palette(PAL_GREY), inverse(false) {}
  MuniTunePalEvent(wxCommandEvent& e):
    wxCommandEvent(e), palette(PAL_GREY), inverse(false) {}
  wxEvent *Clone(void) const { return new MuniTunePalEvent(*this); }
  int palette;
  bool inverse;
};

wxDECLARE_EVENT(EVT_TUNE_PAL,MuniTunePalEvent);

#define MuniTunePalEventHandler(func) (&func)

// -- Tune Colour
class MuniTuneColourEvent: public wxCommandEvent
{
 public:
  MuniTuneColourEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(e,id), saturation(1.0), x_white(1.0/3.0), y_white(1.0/3.0) {}
  MuniTuneColourEvent(wxCommandEvent& e):
    wxCommandEvent(e), saturation(1.0), x_white(1.0/3.0), y_white(1.0/3.0) {}
  wxEvent *Clone(void) const { return new MuniTuneColourEvent(*this); }
  double saturation, x_white, y_white;
};

wxDECLARE_EVENT(EVT_TUNE_COLOUR,MuniTuneColourEvent);

#define MuniTuneColourEventHandler(func) (&func)


// -- Tune Nite
class MuniTuneNiteEvent: public wxCommandEvent
{
 public:
  MuniTuneNiteEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(e,id), nite(false), level(0.0), width(1.0) {}
  MuniTuneNiteEvent(wxCommandEvent& e):
    wxCommandEvent(e), nite(false), level(0.0), width(1.0) {}
  wxEvent *Clone(void) const { return new MuniTuneNiteEvent(*this); }
  bool nite;
  double level, width;
};

wxDECLARE_EVENT(EVT_TUNE_NITE,MuniTuneNiteEvent);

#define MuniTuneNiteEventHandler(func) (&func)


// -- Astrometry

class MuniAstrometryEvent: public wxCommandEvent
{
 public:
  MuniAstrometryEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
  wxCommandEvent(e,id), scale(0.0),reflex(1.0),angle(0.0),xcen(0.0),
    ycen(0.0),acen(0.0),dcen(0.0),astrometry(false) {}
 MuniAstrometryEvent(wxCommandEvent& e):
  wxCommandEvent(e), scale(0.0), reflex(1.0),angle(0.0),xcen(0.0),ycen(0.0),
    acen(0.0),dcen(0.0),astrometry(false) {}
  wxEvent *Clone(void) const { return new MuniAstrometryEvent(*this); }
  wxString proj;
  double scale,reflex,angle,xcen,ycen,acen,dcen;
  bool astrometry;
  MuniLayer layer;
};

wxDECLARE_EVENT(EVT_ASTROMETRY,MuniAstrometryEvent);

//#define MuniAstrometryEventHandler(func) (&func)


// -- Photometry

class MuniPhotometryEvent: public wxCommandEvent
{
 public:
 MuniPhotometryEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY): wxCommandEvent(e,id),erase(false) {}
  MuniPhotometryEvent(wxCommandEvent& e): wxCommandEvent(e),erase(false) {}
  wxEvent *Clone(void) const { return new MuniPhotometryEvent(*this); }
  bool erase;
  MuniLayer layer;
};

wxDECLARE_EVENT(EVT_PHOTOMETRY,MuniPhotometryEvent);

// -- Thread

// pssing of rendered image in thread

class MuniRenderEvent: public wxThreadEvent
{
 public:
 MuniRenderEvent(WXTYPE commandEventType =wxEVT_NULL, int i =wxID_ANY):
  wxThreadEvent(commandEventType,i),id(0),completed(false),x(0),y(0) {}
 MuniRenderEvent(wxThreadEvent& e):
  wxThreadEvent(e),id(0),completed(false),x(0),y(0) {}
  wxEvent *Clone(void) const { return new MuniRenderEvent(*this); }
  int id, rid;
  bool completed;
  int x,y,w,h;
};

wxDECLARE_EVENT(EVT_RENDER,MuniRenderEvent);



class FitsStreamEvent: public wxThreadEvent
{
 public:
 FitsStreamEvent(WXTYPE commandEventType =wxEVT_NULL, int id =wxID_ANY):
   wxThreadEvent(commandEventType,id),progress(0.0),nrows(0),crow(0),cube(0),
   naxis(0), bitpix(0), ncols(0), chdu(0) {}
  FitsStreamEvent(wxThreadEvent& e): wxThreadEvent(e) {}
  wxEvent *Clone(void) const { return new FitsStreamEvent(*this); }
  wxString filename;
  FitsFile fitsfile;
  FitsHeader head;
  wxString extname;
  double progress;
  long nrows, crow, cube;
  int naxis, bitpix, ncols, chdu;
  hdu_type hdutype;
  hdu_flavour flavour;
  wxImage preview;
  std::vector<long> naxes;
  std::vector<wxImage> icons;
  std::vector<FitsTone> tones;
  wxArrayString errlog;
  wxString errmsg;
};

wxDECLARE_EVENT(EVT_FITS_STREAM,FitsStreamEvent);

#define MuniFitsStreamEventHandler(func) (&func)


// -- meta open

// passing of loaded Meta

class MetaOpenEvent: public wxThreadEvent
{
public:
  MetaOpenEvent(WXTYPE commandEventType =wxEVT_NULL, int id =wxID_ANY):
    wxThreadEvent(commandEventType,id) {}
  MetaOpenEvent(wxThreadEvent& e): wxThreadEvent(e) {}
  wxEvent *Clone(void) const { return new MetaOpenEvent(*this); }
  FitsMeta meta;
  wxImage icon;
  size_t id, index;
};

wxDECLARE_EVENT(EVT_META_OPEN,MetaOpenEvent);
#define MuniMetaOpenEventHandler(func) (&func)

// -- FitsMetaOpen progress

class MetaProgressEvent: public wxCommandEvent
{
public:
  MetaProgressEvent(WXTYPE commandEventType =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(commandEventType,id) {}
  MetaProgressEvent(wxCommandEvent& e): wxCommandEvent(e) {}
  wxEvent *Clone(void) const { return new MetaProgressEvent(*this); }
  double phase;
};

wxDECLARE_EVENT(EVT_META_PROGRESS,MetaProgressEvent);
#define MuniMetaProgressEventHandler(func) (&func)

// -- file scan

class DirScanEvent: public wxThreadEvent
{
public:
  DirScanEvent(WXTYPE commandEventType =wxEVT_NULL, int id =wxID_ANY):
    wxThreadEvent(commandEventType,id) {}
  DirScanEvent(wxThreadEvent& e): wxThreadEvent(e) {}
  wxEvent *Clone(void) const { return new DirScanEvent(*this); }
  bool finish;
  size_t id;
  wxArrayString files;
};

wxDECLARE_EVENT(EVT_DIR_SCAN,DirScanEvent);
#define DirScanEventHandler(func) (&func)

// -- Navigation event

class NavigationEvent: public wxCommandEvent
{
public:
  NavigationEvent(WXTYPE commandEventType =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(commandEventType,id) {}
  NavigationEvent(wxCommandEvent& e): wxCommandEvent(e) {}
  wxEvent *Clone(void) const { return new NavigationEvent(*this); }
  wxString pwd,mask;
};

wxDECLARE_EVENT(EVT_NAVIGATION,NavigationEvent);
#define NavigationEventHandler(func) (&func)



// -- general draw layer event

class MuniDrawEvent: public wxCommandEvent
{
 public:
 MuniDrawEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY): wxCommandEvent(e,id) {}
 MuniDrawEvent(wxCommandEvent& e): wxCommandEvent(e) {}
  wxEvent *Clone(void) const { return new MuniDrawEvent(*this); }
  MuniLayer layer;
};

wxDECLARE_EVENT(EVT_DRAW,MuniDrawEvent);


// ---- Histogram event

class MuniHistoEvent: public wxThreadEvent
{
 public:
  MuniHistoEvent(WXTYPE commandEventType =wxEVT_NULL, int id =wxID_ANY):
    wxThreadEvent(commandEventType,id) {}
  MuniHistoEvent(wxThreadEvent& e): wxThreadEvent(e) {}
  wxEvent *Clone(void) const { return new MuniHistoEvent(*this); }
  FitsHisto hist;
};

wxDECLARE_EVENT(EVT_HISTO,MuniHistoEvent);
#define MuniHistoEventHandler(func) (&func)


// -- Config event
class MuniConfigEvent: public wxCommandEvent
{
 public:
  MuniConfigEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
    wxCommandEvent(e,id) {}
  MuniConfigEvent(wxCommandEvent& e): wxCommandEvent(e) {}
  wxEvent *Clone(void) const { return new MuniConfigEvent(*this); }
};
wxDECLARE_EVENT(EVT_CONFIG_UPDATED, MuniConfigEvent);
#define MuniConfigEventHandler(func) (&func)

// -- Size notify event
class MuniSizeChangedEvent: public wxNotifyEvent
{
 public:
  MuniSizeChangedEvent(WXTYPE e =wxEVT_NULL, int id =wxID_ANY):
    wxNotifyEvent(e,id),size(wxDefaultSize) {}
  MuniSizeChangedEvent(wxNotifyEvent& e): wxNotifyEvent(e),size(wxDefaultSize) {}
  wxEvent *Clone(void) const { return new MuniSizeChangedEvent(*this); }
  wxSize size;
};
wxDECLARE_EVENT(EVT_SIZE_CHANGED, MuniSizeChangedEvent);
#define MuniSizeChangedEventHandler(func) (&func)

// last card
#endif
