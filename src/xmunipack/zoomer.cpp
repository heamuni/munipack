/*

  xmunipack - zoomer

  Copyright © 2021-22 F.Hroch (hroch@physics.muni.cz)

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

#include "tune.h"
#include "event.h"
#include "zoomer.h"
#include <wx/wx.h>
#include <wx/graphics.h>
#include <cmath>


MuniZoomWindow::MuniZoomWindow(wxWindow *w, const MuniConfig *c):
  wxPanel(w), config(c), dragging(false), init(true), init_subwin(true),
  subwin_focused(false), zoom_changed(true), tone_changed(true), render_start(false),
  dx(0), dy(0), shrink(-1), iwidth(0), iheight(0), rid(007),
  display_zoom(-1.0), side_ratio(-1), idata(0), render(0)
{
  SetBackgroundStyle(wxBG_STYLE_PAINT);
  bgcolour = wxColour(config->display_bgcolour);

  Bind(wxEVT_IDLE,&MuniZoomWindow::OnIdle,this);
  Bind(wxEVT_SIZE,&MuniZoomWindow::OnSize,this);
  Bind(wxEVT_PAINT,&MuniZoomWindow::OnPaint,this);
  Bind(EVT_ZOOM,&MuniZoomWindow::OnZoom,this);
  Bind(wxEVT_LEFT_UP,&MuniZoomWindow::OnMouseUp,this);
  Bind(wxEVT_LEFT_DOWN,&MuniZoomWindow::OnMouseDown,this);
  Bind(wxEVT_MOTION,&MuniZoomWindow::OnMouseMotion,this);
  Bind(EVT_RENDER,&MuniZoomWindow::OnRender,this);
  Bind(EVT_TUNE,&MuniZoomWindow::OnTune,this);
}

MuniZoomWindow::~MuniZoomWindow()
{
  if( render && render->IsRunning() )
    render->StopAndWait();
  delete render;

  free(idata);
}



void MuniZoomWindow::OnSize(wxSizeEvent& event)
{
  init = true;
  event.Skip();
}

void MuniZoomWindow::OnIdle(wxIdleEvent& event)
{
  //  wxLogDebug("MuniZoomWindow::OnIdle %d %d",render_start,shrink);

  wxSize size(GetSize());
  if( canvas_size != size && size.GetWidth() > 0 && size.GetHeight() ) {
    canvas_size = size;
    InitCanvas();
    init = true;

  }

  if( init && IsInitialised() ) {
    init = false;
    InitScale();
  }

  if( init_subwin && IsInitialised() ) {
    init_subwin = false;
    InitSubwin();
  }


  if( render_start && display_zoom > 0 ) {

    if( render && render->IsRunning() )
      render->StopAndWait();

    wxASSERT(render && ! render->IsRunning());

    render_start = false;

    int w = iwidth;
    int h = iheight;
    iwidth = preview_size.GetWidth();
    iheight = preview_size.GetHeight();

    if( w != iwidth || h != iheight ) {
      free(idata);
      int nbytes = 3*iwidth*iheight;
      idata = (unsigned char *) malloc(nbytes);
      memset(idata,bgcolour.Green(),nbytes);
    }

    wxLogDebug("MuniZoomWindow::OnIdle %d %d %d",iwidth,iheight,shrink);
    render->SetImageData(iwidth,iheight,idata);
    render->SetZoom(1.0/double(shrink));
    render->SetTone(tone,itt,pal,colour);
    rid = rid + 1;
    render->Run(rid);
  }

}

void MuniZoomWindow::OnRender(MuniRenderEvent& event)
{
  //wxLogDebug("MuniZoomWindow::OnRender");

  if( ! (event.rid == rid) ) return;

  if( event.GetId() == ID_SUBRENDER && event.rid == rid ) {

    wxASSERT(canvas.IsOk() && canvas.GetWidth() > 2 && canvas.GetHeight() > 2 &&
	     idata && iwidth > 0 && iheight > 0);

    RefreshCanvas(event.x,event.y,event.w,event.h);
    Refresh();
  }

  else if( event.GetId() == ID_RENDER_FINISH ) {

    wxLogDebug("MuniZoomWindow::OnRender  FINISH");
    render->Wait();

    RefreshCanvas(0,0,canvas.GetWidth(),canvas.GetHeight());
    Refresh();
    zoom_changed = false;
    tone_changed = false;

  }

  else if( event.GetId() == ID_RENDER_INTERRUPT ) {

    wxLogDebug("MuniZoomWindow::OnRender  INTERRUPT");
    render->Wait();
  }

}

void MuniZoomWindow::OnPaint(wxPaintEvent& event)
{
  //  wxLogDebug("MuniZoomWindow::OnPaint %d",canvas.IsOk());

  wxPaintDC dc(this);
  wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
  if( gc ) {
    gc->SetAntialiasMode(wxANTIALIAS_NONE);

    if( canvas.IsOk() ) {

      wxRegionIterator upd(GetUpdateRegion());
      while (upd) {
	int x = upd.GetX();
	int y = upd.GetY();
	int w = upd.GetWidth();
	int h = upd.GetHeight();
	gc->Clip(wxRegion(x,y,w,h));
	gc->DrawBitmap(canvas,0,0,canvas.GetWidth(),canvas.GetHeight());
	upd++;
      }

      // zoom label
      if( display_zoom > 0 )
	DrawZoom(gc);

      // subwindow
      if( ShowSubwin() )
	DrawSubwin(gc);

    }
    else {
      wxSize size(GetSize());
      wxBrush brush(bgcolour);
      gc->SetBrush(brush);
      gc->DrawRectangle(0,0,size.GetWidth(),size.GetHeight());
    }

    delete gc;
  }

}

wxSize MuniZoomWindow::ScaledSize() const
{
  // scaled size of image in the display window
  wxSize scaled_size = display_zoom * image_size;
  return scaled_size;
}

bool MuniZoomWindow::ShowSubwin() const
{
  wxSize scaled_size = ScaledSize();

  // flag whatever the currently scaled image is larger than the display window
  bool o_x = scaled_size.GetWidth()  > display_size.GetWidth();
  bool o_y = scaled_size.GetHeight() > display_size.GetHeight();
  return o_x || o_y;
}

void MuniZoomWindow::InitScale()
{
  wxASSERT(canvas.IsOk() &&
	   canvas.GetWidth() > 0 && canvas.GetHeight() > 0 &&
	   image_size.GetWidth() > 0 && image_size.GetHeight() > 0);

  // width:height
  side_ratio = double(image_size.GetWidth()) / double(image_size.GetHeight());

  // image preview size
  double s1 = std::ceil(double(image_size.GetWidth()) / double(canvas_size.GetWidth()));
  double s2 = std::ceil(double(image_size.GetHeight()) / double(canvas_size.GetHeight()));
  if( s1 > 1 || s2 > 1 )
    shrink = wxMax(int(s1+0.5),int(s2+0.5));
  else
    shrink = 1;
  int w = wxMax(image_size.GetWidth() / shrink,1);
  int h = wxMax(image_size.GetHeight()/ shrink,1);
  preview_size = wxSize(w,h);

  //wxLogDebug("MuniZoomWindow::InitScale preview %d %d %d",w,h,shrink);

  // offset of preview in the canvas
  int xoff = wxMax(canvas_size.GetWidth() - w,0) / 2;
  int yoff = wxMax(canvas_size.GetHeight() - h,0) / 2;
  preview_offset = wxPoint(xoff,yoff);
  //  wxLogDebug("MuniZoomWindow::InitScale preview offsets %d %d",xoff,yoff);

  // center of preview
  subwin_center = wxPoint(preview_size.GetWidth()/2,preview_size.GetHeight()/2);

}

void MuniZoomWindow::InitSubwin()
{
  //  wxLogDebug("MuniZoomWindow::InitSubwin() %d",init_subwin);

  wxSize scaled_size = ScaledSize();

  // subwindow size
  double sx = double(display_size.GetWidth()) / double(scaled_size.GetWidth());
  double sy = double(display_size.GetHeight()) / double(scaled_size.GetHeight());
  sx = wxMin(sx,1.0);
  sy = wxMin(sy,1.0);
  subwin.width = wxMax(int(sx*preview_size.GetWidth() + 0.5),1);
  subwin.height= wxMax(int(sy*preview_size.GetHeight() + 0.5),1);
  wxASSERT(subwin.width > 0 && subwin.height > 0);

  subwin.x = preview_offset.x + subwin_center.x - subwin.width / 2;
  subwin.y = preview_offset.y + subwin_center.y - subwin.height / 2;
  wxLogDebug("MuniZoomWindow::InitSubwin: %d %d %d %d %d %d",
	     subwin.x,subwin.y,subwin.width,subwin.height,
	     subwin_center.x,subwin_center.y);

  MuniZoomEvent event(EVT_ZOOM,ID_ZOOMVIEW);
  event.x = double(subwin_center.x) / double(preview_size.GetWidth());
  event.y = double(subwin_center.y) / double(preview_size.GetHeight());
  event.render = false;
  wxQueueEvent(GetParent(),event.Clone());

}

void MuniZoomWindow::InitCanvas()
{
  int width = canvas_size.GetWidth();
  int height = canvas_size.GetHeight();
  wxASSERT(width > 0 && height > 0);
  wxLogDebug("MuniZoomWindow::InitCanvas() %d %d",width,height);

  canvas = wxBitmap(width,height);
  wxMemoryDC mdc(canvas);
  wxBrush brush(bgcolour);
  wxGraphicsContext *gc = wxGraphicsContext::Create(mdc);
  if( gc ) {
    gc->SetAntialiasMode(wxANTIALIAS_NONE);
    gc->SetBrush(brush);
    gc->DrawRectangle(0,0,width,height);
    delete gc;
  }
}

void MuniZoomWindow::RefreshCanvas(int x, int y, int w, int h)
{
  //  wxLogDebug("MuniZoomWindow::RefreshCanvas() ");
  wxASSERT(canvas.IsOk() && canvas.GetWidth() > 0 && canvas.GetHeight() > 0);
  wxASSERT(iwidth > 0 && iheight > 0 && idata);
  //  wxLogDebug("MuniZoomWindow::RefreshCanvas %d %d %d %d",x,y,w,h);

  wxMemoryDC dc(canvas);
  if( idata ) {
    wxImage display(iwidth,iheight,idata,true);
    wxBitmap bmp(display);
    wxMemoryDC mdc(bmp);
    dc.SetClippingRegion(1,1,canvas.GetWidth()-2,canvas.GetHeight()-2);
    dc.Blit(x+preview_offset.x,y+preview_offset.y,w,h,&mdc,x,y);
  }
}


void MuniZoomWindow::DrawSubwin(wxGraphicsContext *gc)
{
  //  wxLogDebug("MuniZoomWindow::DrawSubwin()");

  wxPenStyle ps = subwin_focused ? wxPENSTYLE_DOT : wxPENSTYLE_SOLID;
  wxColour cpen(wxSystemSettings::GetColour(wxSYS_COLOUR_INACTIVEBORDER));

  // white inner rectangle
  int x1 = subwin.x;
  int y1 = subwin.y;
  int x2 = subwin.x + subwin.width;
  int y2 = subwin.y + subwin.height;
  gc->SetPen(wxPen(cpen,1,ps));
  wxColour cg = subwin_focused ? wxColour(196,196,196) : wxColour(128,128,128,196);
  gc->SetPen(wxPen(*wxWHITE,1,ps));
  wxGraphicsPath path = gc->CreatePath();
  path.MoveToPoint(x1,y1);
  path.AddLineToPoint(x1,y2);
  path.AddLineToPoint(x2,y2);
  path.AddLineToPoint(x2,y1);
  path.AddLineToPoint(x1,y1);
  path.CloseSubpath();
  gc->StrokePath(path);

  // black outer rectangle
  gc->SetPen(wxPen(*wxBLACK,1,ps));
  int u1 = x1 - 1;
  int v1 = y1 - 1;
  int u2 = x2 + 1;
  int v2 = y2 + 1;
  path = gc->CreatePath();
  path.MoveToPoint(u1,v1);
  path.AddLineToPoint(u1,v2);
  path.AddLineToPoint(u2,v2);
  path.AddLineToPoint(u2,v1);
  path.AddLineToPoint(u1,v1);
  path.CloseSubpath();
  gc->StrokePath(path);

  // keep for redraw
  refresh.push_back(wxRect(x1-2,y1-2,subwin.width+5,1+5));
  refresh.push_back(wxRect(x1-2,y2-2,subwin.width+5,1+5));
  refresh.push_back(wxRect(x1-2,y1-2,1+5,subwin.height+5));
  refresh.push_back(wxRect(x2-2,y1-2,1+5,subwin.height+5));

}


void MuniZoomWindow::DrawZoom(wxGraphicsContext *gc)
{
  int width = canvas_size.GetWidth();
  int height = canvas_size.GetHeight();
  wxString label = GetLabel();
  double cw,ch,descent,leading;
  gc->SetFont(*wxNORMAL_FONT,*wxLIGHT_GREY);
  gc->GetTextExtent(label,&cw,&ch,&descent,&leading);
  wxColour whitea(255,255,255,50);
  wxBrush bt(whitea);
  gc->SetBrush(bt);
  gc->SetPen(wxPen(whitea));

  double u = (width-cw)/2;
  double v = height-5-ch;

  gc->DrawRoundedRectangle(u-4,v,cw+12,ch,2);
  gc->DrawText(label,u+2,v+1);
}


wxString MuniZoomWindow::GetLabel() const
{
  wxString label;
  if( display_zoom > 0.99 )
    label.Printf("%d:1",int(display_zoom + 0.5));
  else
    label.Printf("1:%d",int(1.0/display_zoom + 0.5));
  return label;
}


void MuniZoomWindow::OnMouseMotion(wxMouseEvent& event)
{
  //wxLogDebug("%d %d",event.Dragging(),event.LeftUp());

  wxClientDC dc(this);
  wxPoint crosshair = event.GetLogicalPosition(dc);

  if( canvas.IsOk() && ShowSubwin() ) {
    subwin_focused = subwin.Contains(crosshair);
  }
  //    wxLogDebug("%d %d %d",crosshair.x,crosshair.y,subwin_focused);


  if( dragging && event.Dragging() && canvas.IsOk() && ShowSubwin() ) {

    int x = wxMin(wxMax(crosshair.x - dx,0),canvas.GetWidth()-subwin.GetWidth());
    int y = wxMin(wxMax(crosshair.y - dy,0),canvas.GetHeight()-subwin.GetHeight());
    wxRect sub(x,y,subwin.GetWidth(),subwin.GetHeight());
    wxRect full(0,0,canvas.GetWidth(),canvas.GetHeight());

    if( full.Contains(sub) ) {

      //	wxLogDebug("ref regs %d",int(refresh.size()));

      for(size_t i = 0; i < refresh.size(); i++)
	RefreshRect(refresh[i]);
      refresh.clear();

      int xcen = x - preview_offset.x + subwin.GetWidth() / 2;
      int ycen = y - preview_offset.y + subwin.GetHeight() / 2;
      subwin_center = wxPoint(xcen,ycen);
      subwin.x = x;
      subwin.y = y;
    }

    // send up the relative (0.0..1.0) coordinates of left bottom corner
    // Intentionally commented out, it's slow.
    /*
    MuniZoomEvent event(EVT_ZOOM,ID_ZOOMVIEW);
    event.x = double(subwin_center.x) / double(preview_size.GetWidth());
    event.y = double(subwin_center.y) / double(preview_size.GetHeight());
    wxQueueEvent(GetParent(),event.Clone());
    */
  }
}

void MuniZoomWindow::OnMouseUp(wxMouseEvent& event)
{
  if( event.LeftUp() && canvas.IsOk() && ShowSubwin() ) {

    dragging = false;
    wxClientDC dc(this);
    wxPoint crosshair = event.GetLogicalPosition(dc);
    wxLogDebug("MuniZoomWindow::OnMouseMotion UP %d %d",crosshair.x,crosshair.y);
    if( subwin.Contains(crosshair) ) {

      // send up the relative (0.0..1.0) coordinates of the centre
      MuniZoomEvent event(EVT_ZOOM,ID_ZOOMVIEW);
      event.x = double(subwin_center.x) / double(preview_size.GetWidth());
      event.y = double(subwin_center.y) / double(preview_size.GetHeight());
      event.render = true;
      wxQueueEvent(GetParent(),event.Clone());
    }
    wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
    if( gc ) {
      DrawSubwin(gc);
      delete gc;
    }
    Refresh();
  }
}

void MuniZoomWindow::OnMouseDown(wxMouseEvent& event)
{
  if( event.LeftDown() && canvas.IsOk() && ShowSubwin() ) {

    dragging = true;
    wxClientDC dc(this);
    wxPoint crosshair = event.GetLogicalPosition(dc);

    wxLogDebug("MuniZoomWindow::OnMouse DOWN %d %d %d %d %d %d %d",
	       crosshair.x,crosshair.y,subwin.width,subwin.height,
	       subwin.Contains(crosshair),subwin.x,subwin.y);

    if( subwin.Contains(crosshair) ) {
      dx = crosshair.x - subwin.x;
      dy = crosshair.y - subwin.y;
    }

    wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
    if( gc ) {
      DrawSubwin(gc);
      delete gc;
    }
  }
}

void MuniZoomWindow::OnZoom(MuniZoomEvent& event)
{
  wxLogDebug("MuniZoomWindow::OnZoom %f",event.zoom);
  SetDisplayZoom(event.zoom);
}


void MuniZoomWindow::SetDisplayZoom(double z)
{
  display_zoom = z;
  init_subwin = true;
  zoom_changed = true;
  Refresh();
}

void MuniZoomWindow::SetArray(const FitsArray& a)
{
  //  wxLogDebug("MuniZoomWindow::SetArray");
  fitsimage = a;
  display_zoom = -1;
  free(idata);
  idata = 0;
  iwidth = 0;
  iheight = 0;

  delete render;
  render = new MuniDisplayRender(this,fitsimage);
  zoom_changed = true;
  tone_changed = true;
  init_subwin = true;
}


void MuniZoomWindow::SetImageSize(const wxSize& i)
{
  image_size = i;
  init = true;
  wxLogDebug("MuniZoomWindow::SetImageSizes %d %d ",image_size.GetWidth(),
	     image_size.GetHeight());
}

void MuniZoomWindow::SetDisplaySize(const wxSize& d)
{
  display_size = d;
  init = true;
  init_subwin = true;
  zoom_changed = true;
  wxLogDebug("MuniZoomWindow::SetDisplaySize %d %d ",display_size.GetWidth(),
	     display_size.GetHeight());
}

void MuniZoomWindow::SetTone(const FitsTone& t)
{
  wxLogDebug("MuniZoomWindow::SetTone");
  tone = t;
  tone_changed = true;
}

void MuniZoomWindow::SetItt(const FitsItt& it)
{
  wxLogDebug("MuniZoomWindow::SetItt");
  itt = it;
  tone_changed = true;
}

void MuniZoomWindow::SetPalette(const FitsPalette& p)
{
  wxLogDebug("MuniZoomWindow::SetPal");
  pal = p;
  tone_changed = true;
}

void MuniZoomWindow::SetColour(const FitsColour& c)
{
  wxLogDebug("MuniZoomWindow::SetColour");
  colour = c;
  tone_changed = true;
}

void MuniZoomWindow::SetTune(const FitsTone& t, const FitsItt& i,
			     const FitsPalette& p, const FitsColour& c)
{
  tone = t;
  itt = i;
  pal = p;
  colour = c;
  tone_changed = true;
}

void MuniZoomWindow::OnTune(MuniTuneEvent& e)
{
  //  wxLogDebug("MuniZoomWindow::OnTune");

  bool found = true;

  switch(e.GetId()) {

    // tone scaling
  case ID_TONE_BLACK:  tone.SetBlack(e.x);  break;
  case ID_TONE_SENSE:  tone.SetSense(e.x);  break;
  case ID_TONE_QBLACK: tone.SetQblack(e.x); break;
  case ID_TONE_RSENSE: tone.SetRsense(e.x); break;

    // fine tune
  case ID_ITT_TYPE:   itt.SetItt(e.GetString()); break;

  // colours
  case ID_COLOUR_NITEVISION: colour.SetNiteVision(e.GetInt()); break;
  case ID_COLOUR_SATUR:   colour.SetSaturation(e.x);           break;
  case ID_COLOUR_MESO_LEVEL:  colour.SetMesoLevel(e.x);        break;
  case ID_COLOUR_MESO_WIDTH:  colour.SetMesoWidth(e.x);        break;
  case ID_COLOUR_MEAN:    colour.SetLevel(e.index,e.x);        break;
  case ID_COLOUR_WEIGHT:  colour.SetWeight(e.index,e.x);       break;

    // palette
  case ID_PALETTE_TYPE:    pal.SetPalette(e.GetString()); break;
  case ID_PALETTE_INVERSE: pal.SetInverse(e.GetInt());    break;

  default: found = false;
  }

  tone_changed = found;
}


bool MuniZoomWindow::IsInitialised() const
{
  return display_size.GetWidth() > 0 && display_size.GetHeight() &&
    image_size.GetWidth() > 0 && image_size.GetHeight() && display_zoom > 0;
}

void MuniZoomWindow::StartRendering()
{
  wxASSERT(IsInitialised());
  render_start = init || tone_changed || zoom_changed;
}


// -----------------


MuniZoomSlider::MuniZoomSlider(wxWindow *w, wxWindowID id): wxPanel(w)
{
  // ranges must exactly match table in MuniZoomSet
  slider = new wxSlider(this,id,5,0,10);

  wxSizerFlags flaglabel; flaglabel.Align(wxALIGN_CENTER_VERTICAL);
  wxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);
  sizer->Add(new wxStaticText(this,wxID_ANY,L"⭑"),flaglabel.Border(wxLEFT));
  sizer->Add(slider,wxSizerFlags(1).Align(wxALIGN_CENTER_VERTICAL));
  sizer->Add(new wxStaticText(this,wxID_ANY,L"★"),flaglabel.Border(wxRIGHT));
  SetSizer(sizer);

  Bind(wxEVT_SCROLL_THUMBTRACK,&MuniZoomSlider::OnScroll,this,id);
  Bind(wxEVT_SCROLL_LINEUP,&MuniZoomSlider::OnScroll,this,id);
  Bind(wxEVT_SCROLL_LINEDOWN,&MuniZoomSlider::OnScroll,this,id);
  Bind(wxEVT_SCROLL_PAGEUP,&MuniZoomSlider::OnScroll,this,id);
  Bind(wxEVT_SCROLL_PAGEDOWN,&MuniZoomSlider::OnScroll,this,id);
}

void MuniZoomSlider::OnScroll(wxScrollEvent& event)
{
  MuniZoomSet zset(event.GetPosition());
  MuniZoomEvent e(EVT_ZOOM,slider->GetId());
  e.zoom = zset.GetZoom();
  wxQueueEvent(GetParent(),e.Clone());
  wxLogDebug("MuniZoomSlider::OnScroll %d %f",event.GetPosition(),e.zoom);
}


void MuniZoomSlider::SetZoom(double zoom)
{
  wxLogDebug("MuniZoomSlider::SetZoom %f",zoom);
  MuniZoomSet zset(zoom);
  slider->SetValue(zset.GetIndex());
}

double MuniZoomSlider::GetZoom() const
{
  return GetZoom(slider->GetValue());
}

double MuniZoomSlider::GetZoom(int z) const
{
  MuniZoomSet set(z);
  return set.GetZoom();
}
