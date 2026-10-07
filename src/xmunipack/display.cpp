/*

  xmunipack - display image

  Copyright © 2009-2013, 2017-2024 F.Hroch (hroch@physics.muni.cz)

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
#include "render.h"
#include "legend.h"
#include <wx/wx.h>
#include <wx/graphics.h>
#include <wx/mstream.h>
#include <cfloat>
#include <cmath>

using namespace std;



// -- MuniDisplay

MuniDisplay::MuniDisplay(wxWindow *w, MuniConfig *c):
  wxPanel(w,wxID_ANY,wxDefaultPosition,wxSize(100,100),
	  wxTAB_TRAVERSAL|wxWANTS_CHARS),
  config(c),render(0),idata(0),iwidth(0),iheight(0),
  rid(007),xcen(0.5),ycen(0.5),
  zoom(-1.0),render_start(false), display_size(wxDefaultSize),
  find(0),astrometry(0)
{
  wxASSERT(config);

  bgcolour = wxColour(config->display_bgcolour);
  SetBackgroundStyle(wxBG_STYLE_SYSTEM);
  SetBackgroundColour(bgcolour);

  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniDisplay::OnLeaveFullscreen,this,
       ID_MENU_FULLSCREEN);
  Bind(wxEVT_MOTION,&MuniDisplay::OnMouseMotion,this);
  Bind(wxEVT_ENTER_WINDOW,&MuniDisplay::OnMouseEnter,this);
  Bind(wxEVT_LEAVE_WINDOW,&MuniDisplay::OnMouseLeave,this);
  Bind(wxEVT_KEY_DOWN,&MuniDisplay::OnKeyDown,this);
  Bind(EVT_ZOOM,&MuniDisplay::OnZoomer,this);
  Bind(wxEVT_IDLE,&MuniDisplay::OnIdle,this);
  Bind(wxEVT_SIZE,&MuniDisplay::OnSize,this);
  Bind(EVT_RENDER,&MuniDisplay::OnRender,this);
  Bind(wxEVT_PAINT,&MuniDisplay::OnPaint,this);
}

MuniDisplay::~MuniDisplay()
{
  if( render && render->IsRunning() )
    render->StopAndWait();
  delete render;
  free(idata);

  wxLogDebug("MuniDisplay::~MuniDisplay()");
}

void MuniDisplay::OnSize(wxSizeEvent& event)
{
  wxSize size(event.GetSize());

  if( size != display_size ) {
    wxLogDebug("MuniDisplay::OnSize changed %dx%d",size.GetWidth(),size.GetHeight());
    display_size = size;
    render_start = true;

    MuniSizeChangedEvent event(EVT_SIZE_CHANGED);
    event.size = size;
    wxQueueEvent(GetParent(),event.Clone());
  }
  event.Skip();
}


void MuniDisplay::OnIdle(wxIdleEvent& event)
{
  //wxLogDebug("MuniDisplay::OnIdle");
  if( array.IsOk() && render_start ) {
    wxASSERT(render);
    Render();
  }
}

void MuniDisplay::SetArray(const FitsArray& a, const FitsTone& t)
{
  wxASSERT(a.IsOk() && t.IsOk());
  array = a;
  xcen = 0.5;
  ycen = 0.5;
  tone = t;
  zoom = -1;
  display_size = wxDefaultSize;

  delete render;
  render = new MuniDisplayRender(this,array);
}

wxRect MuniDisplay::GetArrayRect() const
{
  wxASSERT(array.IsOk() && zoom > 0);

  wxSize size(GetSize());
  int w = int(size.GetWidth() / zoom);
  int h = int(size.GetHeight()/ zoom);
  int x1 = xcen * array.GetWidth() - w/2;
  int y1 = ycen * array.GetHeight() - h/2;

  /*
  wxLogDebug("MuniDisplay::GetCurrentRect() %d %d %d %d, %ld %ld %d %d, %f %f %f",
	     x1,y1,w,h,array.GetWidth(),array.GetHeight(),size.GetWidth(),
	     size.GetHeight(),zoom,xcen,ycen);
  */
  return wxRect(x1,y1,w,h);
}

wxRect MuniDisplay::GetExposed() const
{
  wxASSERT(array.IsOk());

  wxRect drect(GetArrayRect());
  wxRect irect(0,0,array.GetWidth(),array.GetHeight());
  wxRect exposed = irect.Intersect(drect);

  /*
  wxLogDebug("MuniDisplay::GetExposed() %d %d %d %d",
	     exposed.x,exposed.y,exposed.width,exposed.height);
  wxLogDebug("%d %d %d %d",drect.width,drect.height,irect.width,irect.height);
  */
  return exposed;
}


wxPoint MuniDisplay::GetOffset(const wxRect& exposed) const
{
  wxASSERT(array.IsOk() && exposed.GetWidth() > 0 && exposed.GetHeight() > 0);

  wxSize size(GetSize());
  wxRect drect(GetArrayRect());

  // offset
  int zwidth =  int(zoom*exposed.GetWidth());
  int zheight = int(zoom*exposed.GetHeight());
  int rx = zwidth > size.GetWidth() ? 1 : 2;
  int ry = zheight > size.GetHeight() ? 1 : 2;
  int xoff = drect.x < 0 ? wxMax(size.GetWidth() - zwidth,0) / rx : 0;
  int yoff = drect.y < 0 ? wxMax(size.GetHeight()- zheight,0)/ ry : 0;

  /*
    wxLogDebug("MuniDisplay::GetOffset %d %d %d %d, off: %d %d",
    rx,ry,size.GetWidth(),size.GetHeight(),xoff,yoff);
  */

  return wxPoint(xoff,yoff);
}


void MuniDisplay::StartRendering()
{
  render_start = true;
}

void MuniDisplay::Render()
{
  rid = rid + 1;
  render_start = false;

  if( render->IsRunning() ) {
    render->StopAndWait();
    ::wxEndBusyCursor();
  }

  exposed = GetExposed();
  offset = GetOffset(exposed);
  if( exposed.GetWidth() == 0 || exposed.GetHeight() == 0 ) return;
  wxRect iexposed(exposed);
  iexposed.SetY(array.GetHeight() - (exposed.y + exposed.height));

  int wold = iwidth;
  int hold = iheight;
  iwidth = wxMax(int(zoom*exposed.GetWidth()),1);
  iheight = wxMax(int(zoom*exposed.GetHeight()),1);
  wxLogDebug("MuniDisplay::OnIdle  iwidth=%d iheight=%d",iwidth,iheight);

  if( iwidth != wold || iheight != hold ) {
    unsigned char *old = idata;
    size_t nbytes = 3*iwidth*iheight;
    idata = (unsigned char *) malloc(nbytes);
    if( old ) {
      wxImage img(wold,hold,old);
      img.Rescale(iwidth,iheight);
      memcpy(idata,img.GetData(),nbytes);
      // free(old) is called by wxImage destructor.
    }
    else
      memset(idata,bgcolour.Green(),nbytes);
  }

  DrawCanvas();

  render->SetImageData(iwidth,iheight,idata);
  render->SetZoom(zoom);
  render->SetExposed(iexposed);
  render->SetTone(tone,itt,pal,colour);
  render->Run(rid);
  ::wxBeginBusyCursor();

}

void MuniDisplay::OnRender(MuniRenderEvent& event)
{
  //  wxLogDebug("MuniDisplay::OnRender");

  if( ! (rid == event.rid) ) return;

  if( event.GetId() == ID_SUBRENDER && render_start == false ) {

    wxASSERT(exposed.GetWidth() > 0 && exposed.GetHeight() > 0);
    //wxLogDebug("MuniDisplay::OnRender %d %d %d %d",event.x,exposed.x,event.y,exposed.y);
    int j = event.y;
    int x = event.x + offset.x;
    int y = event.y + offset.y;
    int w = event.w;
    int h = event.h;

    wxASSERT(0 <= j && j+h <= iheight && w == iwidth);

    // draw the current sub-window
    wxClientDC dc(this);
    wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
    if( gc ) {
      wxRect rect(x,y,w,h);
      gc->Clip(wxRegion(rect));
      unsigned char *data = idata + 3*w*j;
      wxImage subimg(w,h,data,true);
      wxGraphicsBitmap bmp = gc->CreateBitmapFromImage(subimg);
      gc->DrawBitmap(bmp,x,y,w,h);
      gc->ResetClip();
      delete gc;
    }

  }
  else if( event.GetId() == ID_RENDER_FINISH ) {

    wxLogDebug("MuniDisplay::OnRender FINISH");
    render->Wait();
    ::wxEndBusyCursor();
    DrawCanvas();
    Refresh();

    // notify parent
    wxQueueEvent(GetParent(),event.Clone());
  }

  else if( event.GetId() == ID_RENDER_INTERRUPT ) {

    wxLogDebug("MuniDisplay::OnRender INTERRUPT");
    render->Wait();
    ::wxEndBusyCursor();
    Refresh();
  }
}


void MuniDisplay::DrawCanvas()
{
  wxASSERT(idata);

  wxSize size(GetSize());
  const wxRect exposed(GetExposed());
  int xoff = offset.x;
  int yoff = offset.y;
  int width = size.GetWidth();
  int height = size.GetHeight();
  canvas = wxBitmap(width,height);

  //  wxLogDebug("DrawCanvas(): %d %d",display_size.GetWidth(),display_size.GetHeight());
  wxMemoryDC mdc(canvas);
  wxGraphicsContext *gc = wxGraphicsContext::Create(mdc);
  if( gc ) {

    wxBrush brush(bgcolour);
    gc->SetBrush(brush);
    gc->DrawRectangle(0,0,width,height);

    wxImage image(iwidth,iheight,idata,true);
    wxGraphicsBitmap bmp = gc->CreateBitmapFromImage(image);
    gc->DrawBitmap(bmp,xoff,yoff,iwidth,iheight);

    if( config->display_legend && array.HasWCS() ) {
      MuniDisplayLegend *legend =
	new MuniDisplayLegend(gc,config,exposed,xoff,yoff,zoom);
      legend->Draw(array);
      delete legend;
    }
    delete gc;
  }
}


void MuniDisplay::OnPaint(wxPaintEvent& event)
{
  //  wxLogDebug("MuniDisplay::OnPaint");

  wxPaintDC dc(this);
  wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
  if( gc && canvas.IsOk() ) {
    wxRegionIterator damaged(GetUpdateRegion());
    while( damaged ) {
      gc->Clip(wxRegion(damaged.GetRect()));
      gc->DrawBitmap(canvas,0,0,canvas.GetWidth(),canvas.GetHeight());
      gc->ResetClip();
      damaged++;
    }
    delete gc;
  }
}


void MuniDisplay::SetPreview(const wxImage& img)
{
  wxLogDebug("MuniDisplay::SetPreview(%p)",&img);

  if( img.IsOk() ) {
    wxASSERT(idata == 0);
    iwidth = img.GetWidth();
    iheight = img.GetHeight();
    size_t npix = 3*iwidth*iheight;
    idata = (unsigned char *) malloc(npix);
    memcpy(idata,img.GetData(),npix);
  }
}

void MuniDisplay::SetZoom(double z)
{
  if( abs(zoom - z) < 1e-3 ) return;
  zoom = z;
  render_start = true;
}

double MuniDisplay::GetZoom() const
{
  return zoom;
}

void MuniDisplay::SetTone(const FitsTone& t)
{
  tone = t;
  render_start = true;
}

FitsTone MuniDisplay::GetTone() const
{
  return tone;
}

void MuniDisplay::SetItt(const FitsItt& it)
{
  itt = it;
  render_start = true;
}

FitsItt MuniDisplay::GetItt() const
{
  return itt;
}

void MuniDisplay::SetColour(const FitsColour& c)
{
  colour = c;
  render_start = true;
}

FitsColour MuniDisplay::GetColour() const
{
  return colour;
}

FitsPalette MuniDisplay::GetPalette() const
{
  return pal;
}

void MuniDisplay::SetPalette(const FitsPalette& p)
{
  pal = p;
  render_start = true;
}

wxImage MuniDisplay::GetImage() const
{
  int nbytes = 3*iwidth*iheight;
  unsigned char *rgb = (unsigned char *) malloc(nbytes);
  memcpy(rgb,idata,nbytes);
  return wxImage(iwidth,iheight,rgb);
}


void MuniDisplay::OnMouseMotion(wxMouseEvent& event)
{
  if( ! (array.GetWidth() > 0 && array.GetHeight() > 0) ) return;

  wxClientDC dc(this);
  wxPoint coo = event.GetLogicalPosition(dc) - offset;
  //  wxLogDebug("MuniDisplay::OnMouseMotion %d %d %d %d %f %d",coo.x,coo.y);

  wxRect r = GetExposed();
  int x = int(coo.x / zoom) + r.x;
  int y = array.GetHeight() - (int(coo.y / zoom) + r.y) - 1;

  MuniSlewEvent evt(EVT_SLEW);
  evt.x = x;
  evt.y = y;
  evt.xsub = coo.x;
  evt.ysub = coo.y;
  evt.zoom = zoom;
  evt.inside = x >= 0 && y >= 0 && x < array.GetWidth() && y < array.GetHeight();
  wxQueueEvent(GetParent(),evt.Clone());
}

void MuniDisplay::OnMouseEnter(wxMouseEvent& event)
{
  MuniSlewEvent evt(EVT_SLEW);
  evt.entering = true;
  wxQueueEvent(GetParent(),evt.Clone());
}

 void MuniDisplay::OnMouseLeave(wxMouseEvent& event)
{
  MuniSlewEvent evt(EVT_SLEW);
  evt.leaving = true;
  wxQueueEvent(GetParent(),evt.Clone());
}


void MuniDisplay::OnClick(MuniClickEvent& event)
{
  if( find ) {
    find->SetPoint(event.x,event.y);
  }

  wxQueueEvent(GetParent(),event.Clone());
}

void MuniDisplay::OnLeaveFullscreen(wxCommandEvent& event)
{
  wxCommandEvent e(wxEVT_COMMAND_MENU_SELECTED,ID_FULLSCREEN);
  wxQueueEvent(GetParent(),e.Clone());
}

/*
void MuniDisplay::SetOverlay(wxInputStream& istream)
{
  wxASSERT(canvas);

  MuniSvg svg(istream);
  static_cast<MuniDisplayCanvas *>(canvas)->SetOverlay(svg.GetDrawObjects());
}
*/

void MuniDisplay::SetStars(FitsTable& table)
{
  //  wxASSERT(canvas);

  wxMemoryOutputStream ostream;
  table.GetStarChart(ostream);

  wxMemoryInputStream istream(ostream);
  //  MuniSvg svg(istream);
  //  static_cast<MuniDisplayCanvas *>(canvas)->AddOverlay(svg.GetDrawObjects());
}

/*
vector<MuniDrawBase *> MuniSvg::GetDrawObjects() const
{
  vector<MuniDrawBase *> objects;

  const wxXmlNode *node = GetRoot()->GetChildren();
  while(node) {

    if( node->GetName() == "circle" ) {

      double x;
      float cx=0.0,cy=0.0,r=0.0;
      wxString a;
      wxXmlAttribute *prop = node->GetAttributes();
      while(prop) {

	if( prop->GetName() == "cx" ) {
	  a = prop->GetValue();
	  if( a.ToDouble(&x) )
	    cx = x;
	}
	else if( prop->GetName() == "cy" ) {
	  a = prop->GetValue();
	  if( a.ToDouble(&x) )
	    cy = x;
	}
	else if( prop->GetName() == "r" ) {
	  a = prop->GetValue();
	  if( a.ToDouble(&x) )
	    r = x;
	}

	prop = prop->GetNext();
      }

      objects.push_back(new MuniDrawCircle(cx,cy,r));
    }

    node = node->GetNext();
  }

  return objects;
}
*/


void MuniDisplay::ShowSources(bool show, const FitsFile& fits)
{
  //  wxASSERT(canvas);

  if( show ) {

    bool apext = false;
    for(size_t i = 0; i < fits.HduCount(); i++) {
      if( fits.Hdu(i).GetExtname() == APEREXTNAME )
	apext = true;
    }

    for(size_t i = 0; i < fits.HduCount(); i++) {
      if( fits.Hdu(i).GetExtname() == APEREXTNAME ||
	  fits.Hdu(i).GetExtname() == PHOTOEXTNAME ||
         (fits.Hdu(i).GetExtname() == FINDEXTNAME && !apext) ) {
	FitsTable table(fits.Hdu(i));
	const double *x = table.GetColumn(FITS_LABEL_X).GetCol_double();
	const double *y = table.GetColumn(FITS_LABEL_Y).GetCol_double();
	double *f = 0;
	if( fits.Hdu(i).GetExtname() == APEREXTNAME )
	  f = (double *)table.GetColumn(FITS_LABEL_APCOUNT "1").GetCol_double();
	else if( fits.Hdu(i).GetExtname() == PHOTOEXTNAME )
	  f = (double *) table.GetColumn(FITS_COL_PHOTONS).GetCol_double();

	vector<double> xcoo,ycoo,flux;
	for(long i = 0; i < table.Nrows(); i++) {
	  xcoo.push_back(x[i]);
	  ycoo.push_back(y[i]);
	  if( f )
	    flux.push_back(f[i]);
	  else
	    flux.push_back(-1.0);
	}

	MuniStarLayer layer;
	if( fits.Hdu(i).GetExtname() == APEREXTNAME )
	  layer.SetHWHM(fits.Hdu(i).GetKeyDouble("HWHM"));
	layer.DrawObjects(xcoo,ycoo,flux);
	//	static_cast<MuniDisplayCanvas *>(canvas)->AddLayer(layer.GetLayer());
      }
    }
  }
  else
    ;
    //    static_cast<MuniDisplayCanvas *>(canvas)->RemoveLayers(ID_PHOTOMETRY);
}


void MuniDisplay::OnFind(const wxString& fitsname, const FitsHdu& hdu)
{
  wxASSERT(find == 0);

  find = new MuniFind(this,config,fitsname,hdu);
  find->Show();
  Bind(EVT_FINISH_DIALOG,&MuniDisplay::OnCloseFind,this,find->GetId());
}

void MuniDisplay::OnCloseFind(wxCommandEvent& event)
{
  wxLogDebug("MuniDisplay::OnCloseFind");
  wxASSERT(find && static_cast<MuniView *>(GetParent()));

  /*
  int code = event.GetInt();
  wxString file = event.GetString();
  */

  Unbind(EVT_FINISH_DIALOG,&MuniDisplay::OnCloseFind,this,find->GetId());
  find = 0;
}


void MuniDisplay::Calibrate(const wxString& file)
{
  // check availability of write access
  wxFileName fn(file);
  if( ! fn.IsFileWritable() ) {
    wxLogWarning("File `"+file+"' is not writtable. Results of calibration should not be saved.");
  }

  MuniCalibrate calibrate(this,config,file);
  if( calibrate.ShowModal() == wxID_OK ) {

    /*
    wxString cfile = calibrate.GetResult();

    FitsFile fits(cfile);

    if( fits.IsOk() )
      for(size_t i = 0; i < fits.HduCount(); i++)
	if( fits.Hdu(i).GetKey("EXTNAME").Find("MUNIPACK") != wxNOT_FOUND ) {
	  FitsTable t(fits.Hdu(i));
	  SetStars(t);
	}

    wxRemoveFile(cfile);
    */

    wxCommandEvent e(EVT_FILELOAD,GetId());
    e.SetEventObject(GetParent());
    e.SetString(file);
    wxQueueEvent(GetParent(),e.Clone());


  }
}

void MuniDisplay::Astrometry(const wxString& file, const FitsTable& t)
{
  if( astrometry ) return;

  astrometry = new MuniAstrometry(this,config);
  astrometry->SetFile(file,array);
  astrometry->SetDetectedSources(t);
  astrometry->Show();

  Bind(EVT_FINISH_DIALOG,&MuniDisplay::OnCloseAstrometry,this,astrometry->GetId());

}

void MuniDisplay::OnCloseAstrometry(wxCommandEvent& event)
{
  wxASSERT(astrometry && static_cast<MuniView *>(GetParent()));

  int code = event.GetInt();
  wxString file = event.GetString();

  //wxLogDebug("Return code= %d %d %d %d %d",code,wxID_APPLY,code==wxID_APPLY,astrometry->GetAffirmativeId(),wxID_CANCEL);

  Unbind(EVT_FINISH_DIALOG,&MuniDisplay::OnCloseAstrometry,this,astrometry->GetId());
  //  astrometry->Destroy();
  astrometry = 0;

  if( code == wxID_APPLY )
    static_cast<MuniView *>(GetParent())->LoadFileBackup(file,file+"~");

}


void MuniDisplay::OnZoomer(MuniZoomEvent& e)
{
  if( e.GetId() == ID_ZOOMER ) {
    // zoom

    zoom = e.zoom;
    render_start = true;

  }
  else if( e.GetId() == ID_ZOOMVIEW ) {
    // subwin offset; e.x, e.y are relative coordinates of centre

    if( array.IsOk() ) {
      xcen = e.x;
      ycen = e.y;
      render_start = e.render;
    }
  }
}


void MuniDisplay::OnKeyDown(wxKeyEvent& event)
{
  //  wxLogDebug("MuniDisplay::OnKeyDown() %d %d",int(event.GetX()),int(event.GetY()));

  long keycode = event.GetKeyCode();

  // ESCAPE full-screen
  wxTopLevelWindow *twin = static_cast<wxTopLevelWindow *>(GetParent());
  if( twin && twin->IsFullScreen() && keycode == WXK_ESCAPE ){
    wxQueueEvent(GetParent(),
		 new wxCommandEvent(wxEVT_COMMAND_MENU_SELECTED,ID_FULLSCREEN));
    event.Skip();
    return;
  }

  long x = event.GetX();
  long y = event.GetY();

  switch(keycode) {
  case WXK_LEFT:
    x--; break;
  case WXK_UP:
    y--; break;
  case WXK_RIGHT:
    x++; break;
  case WXK_DOWN:
    y++; break;
  }
  WarpPointer(x,y);

  event.Skip();
}


void MuniDisplay::ConfigUpdate()
{
  DrawCanvas();
  Refresh();
}

/*

   HOW to create a throbbler

   - download glunarclock source package:
        http://glunarclock.sourceforge.net/

   - separate the long image sequence and resize images to single exposures:

       for A in `seq 0 55`; do
          P=$((A*48));
          X=$(printf "%02d" $A);
          convert moon_56frames.png -crop 48x48+$P+1 -resize 32x32 m_$X.png;
       done

   - animate:

       convert -delay 10  +repage m_*.png throbbler.gif


*/
