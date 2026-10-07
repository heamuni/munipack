/*

  xmunipack - display canvas

  Copyright © 2012-14, 2018-2021 F.Hroch (hroch@physics.muni.cz)

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
#include "display.h"
#include <wx/wx.h>
#include <wx/scrolwin.h>
#include <wx/brush.h>
#include <wx/graphics.h>
#include <wx/renderer.h>
#include <wx/thread.h>
#include <wx/clipbrd.h>
#include <list>

using namespace std;


// -- MuniDisplayCanvas

MuniDisplayCanvas::MuniDisplayCanvas(wxWindow *w, MuniConfig *c):
  wxWindow(w,wxID_ANY),render(0),canvasid(0),
  topwin(GetGrandParent()),config(c),handler(0),
  xoff(0), yoff(0), xsubwin(0), ysubwin(0),
  value_type(UNIT_COUNT), coo_type(COO_PIXEL),
  scaleicon(false), shrink(1),zoom(1),
  zooming(false), shrinking(false),tunning(false), completed(false),
  dragging(false),invoking(false), rendering(false), finished(false),
  astrometry(false)
{
  wxLogDebug("MuniDisplayCanvas constructor..");

  SetBackgroundColour(wxColour("black"));

  pal.SetPalette(config->display_pal);
  pal.SetInverse(config->display_palinv);

  Bind(wxEVT_CLOSE_WINDOW,&MuniDisplayCanvas::OnClose,this);
  //  Bind(wxEVT_SIZE,&MuniDisplayCanvas::OnSize,this);
  Bind(wxEVT_PAINT,&MuniDisplayCanvas::OnPaint,this);
  Bind(wxEVT_IDLE,&MuniDisplayCanvas::OnIdle,this);
  Bind(EVT_RENDER,&MuniDisplayCanvas::OnRenderFinish,this);
  Bind(EVT_RENDER,&MuniDisplayCanvas::OnSubRender,this,ID_SUBRENDER);
  Bind(EVT_DRAW,&MuniDisplayCanvas::OnDraw,this);
  Bind(wxEVT_KEY_DOWN,&MuniDisplayCanvas::OnKeyDown,this);
  Bind(wxEVT_RIGHT_DOWN,&MuniDisplayCanvas::OnMenu,this);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniDisplayCanvas::OnClipValue,this,
       ID_CLIP_VALUE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniDisplayCanvas::OnClipCoo,this,
       ID_CLIP_COO);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniDisplayCanvas::OnValueType,this,
       ID_VALTYPE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniDisplayCanvas::OnCooType,this,
       ID_COOTYPE);

  // SetHandler clone
  Bind(wxEVT_MOTION,&MuniDisplayCanvas::OnMouseMotion,this);
  Bind(wxEVT_ENTER_WINDOW,&MuniDisplayCanvas::OnMouseEnter,this);
  Bind(wxEVT_LEAVE_WINDOW,&MuniDisplayCanvas::OnMouseLeave,this);
  Bind(wxEVT_LEFT_UP,&MuniDisplayCanvas::OnLeftUp,this);
  Bind(wxEVT_LEFT_DOWN,&MuniDisplayCanvas::OnClick,this);
  Bind(EVT_TUNE,&MuniDisplayCanvas::OnTuneFine,this);
  Bind(EVT_ZOOM,&MuniDisplayCanvas::OnZoom,this);

}

MuniDisplayCanvas::~MuniDisplayCanvas()
{
  /*
  vector<MuniDrawBase *>::const_iterator i;
  for(i = drawobj.begin(); i != drawobj.end(); ++i)
    delete *i;
  */

}

void MuniDisplayCanvas::OnClose(wxCloseEvent& event)
{
  //  wxLogDebug("MuniDisplayCanvas::OnClose");

  if( render ) {
    StopRendering();
    //    render->Delete();
    event.Veto();
    wxQueueEvent(this,new wxCloseEvent());
  }
  else
    event.Skip();
}

void MuniDisplayCanvas::OnIdle(wxIdleEvent& event)
{
  //  wxLogDebug("MuniDisplayCanvas::OnIdle");

  /*
  wxSize s = GetSize();
  wxLogDebug("%d %d",s.GetWidth(),s.GetHeight());
  */

  // Rendering is started by Idle event FOLLOWING all resize(!) events
  // The canvas with a right size is prepared for drawing.

  if( invoking ) {

    if( render ) {
      StopRendering();
      //      event.RequestMore();
      //      return;
    }

    invoking = false;
    finished = false;
    ::wxBeginBusyCursor();
    StartRendering();
  }

  if( finished ) {
    ::wxEndBusyCursor();
    rendering = false;
  }

}

void MuniDisplayCanvas::OnSize(wxSizeEvent& event)
{
  /*
  wxSize s = event.GetSize();
  wxLogDebug("%d %d",s.GetWidth(),s.GetHeight());
  */
  //  UpdateScroll();
  event.Skip();
}

void MuniDisplayCanvas::OnPaint(wxPaintEvent& event)
{
  //  wxLogDebug("MuniDisplayCanvas::OnPaint");

  // display size
  wxSize size = GetSize();
  int width = size.GetWidth();
  int height = size.GetHeight();

  wxPaintDC dc(this);
  //dc.SetLogicalScale(1.0/zoom,1.0/zoom);
  if( pal.GetInverse() )
    dc.SetBackground(*wxWHITE_BRUSH);
  else
    //    dc.SetBackground(*wxBLACK_BRUSH);
    dc.SetBackground(wxColour("grey15"));

  if( canvas.IsOk() ) {
    wxBitmap bmp(canvas);
    wxMemoryDC mdc(bmp);

    // image size
    int iwidth = fitsimage.GetWidth();
    int iheight = fitsimage.GetHeight();

    double sx = double(width) / double(zoom*iwidth);
    double sy = double(height) / double(zoom*iheight);
    sx = wxMin(sx,1.0);
    sy = wxMin(sy,1.0);

    int w = zoom*iwidth;
    int h = zoom*iheight;
    int swidth = int(iwidth*sx + 0.5);
    int sheight= int(iheight*sy + 0.5);
    xoff = wxMax((width-w)/2,0);
    yoff = wxMax((height-h)/2,0);

    /*
    wxLogDebug("MuniDisplayCanvas::OnPaint %d %d %d %d",
	       iwidth,iheight,swidth,sheight);
    */

    /*
    if( zoom > 1 )
      dc.StretchBlit(xoff,yoff,w,h,&mdc,0,0,iwidth,iheight);
    else
      dc.Blit(xoff,yoff,iwidth,iheight,&mdc,0,0);

    dc.SetBrush(*wxTRANSPARENT_BRUSH);
    dc.SetPen(*wxYELLOW_PEN);

    if( shrink > 1 )
      dc.DrawRectangle(xoff+xsubwin,yoff+ysubwin,iwidth/shrink,iheight/shrink);
    else if( zoom > 1 )
      dc.DrawRectangle(xoff+xsubwin,yoff+ysubwin,width,height);
    else
      dc.DrawRectangle(xoff+xsubwin,yoff+ysubwin,iwidth,iheight);

    return;
    */


    if( zoom > 1 ) {
      //      wxLogDebug("%d %d %d %d %d",width,height,iwidth,iheight,zoom);
      w = wxMin(w,width);
      h = wxMin(h,height);
      double r = double(width) / double(height);
      swidth = int(width/zoom + 0.5) + 1;
      sheight = int(swidth/r + 0.5);
      dc.StretchBlit(xoff,yoff,w,h,&mdc,xsubwin,ysubwin,swidth,sheight);
    }
    else {
      w = iwidth/shrink;
      h = iheight/shrink;
      xoff = wxMax((width-w)/2,0);
      yoff = wxMax((height-h)/2,0);
      dc.Blit(xoff,yoff,w,h,&mdc,xsubwin,ysubwin);
    }
  }

  //  dc.Clear();

  /*
  if( canvas.IsOk() )
    dc.DrawBitmap(wxBitmap(canvas),xoff,yoff);
  */
  // Layers
  /*
  list<MuniLayer>::const_iterator i;
  for(i = layers.begin(); i != layers.end(); ++i)
    DrawLayer(dc,*i);
  */

  return;

  // it's kept for later development

  int x0,y0,wc,hc;
  /*
  GetViewStart(&x0,&y0);
  GetClientSize(&wc,&hc);
  */

  wxRegionIterator upd = GetUpdateRegion();
  while(upd) {
    wxRect r = upd.GetRect();

    int xs = r.GetX();
    int ys = r.GetY();
    int ws = r.GetWidth();
    int hs = r.GetHeight();
    wxLogDebug("damaged: %d %d %d %d",xs,ys,ws,hs);

    int xw = xoff > 0 ? xs - xoff : x0 + xs;
    int yw = yoff > 0 ? ys - yoff : y0 + ys;
    int ww = ws - 2*xoff;
    int hw = hs - 2*yoff;

    if( false && canvas.IsOk() ) {

      int xb = wxMax(xw,0);
      int yb = wxMax(yw,0);
      int wb = wxMin(ww,canvas.GetWidth());
      int hb = wxMin(hw,canvas.GetHeight());

      bool wt = 0 < wb && wb <= canvas.GetWidth();
      bool ht = 0 < hb && hb <= canvas.GetHeight();
      bool xt = xb >= 0 && xb + wb <= canvas.GetWidth();
      bool yt = yb >= 0 && yb + hb <= canvas.GetHeight();

      //      wxLogDebug("canvas: %d %d %d %d %d %d %d %d %f %d %d %d %d %d %d %d",
      //		 xs,ys,ws,hs,xb,yb,wb,hb,zoom,canvas.GetWidth(),canvas.GetHeight(),x0,y0,xoff,yoff,
      //		 (int) wt && ht && xt && yt);

      if( wt && ht && xt && yt ) {

	wxImage i(canvas.GetSubImage(wxRect(xb,yb,wb,hb)));
	if( zoom > 1 ) {
	  i.Rescale(zoom*wb,zoom*hb);
	}
	dc.DrawBitmap(wxBitmap(i),xs+xoff,ys+yoff);
      }
      else
	wxLogDebug("failed canvas: %d %d %d %d",xb,yb,wb,hb);
    }

    upd++;
  }

}

void MuniDisplayCanvas::SetHdu(const FitsArray& a, const wxImage& i)
{
  wxLogDebug("MuniDisplayCanvas::SetHdu EXTNAME="+a.GetExtname()+" %p",this);
  array = a;
  icon = i;
  if( icon.HasAlpha() )
    icon.ClearAlpha();
  wxASSERT(array.IsOk() && icon.IsOk());

  fitsimage = FitsImage(array);
  wxASSERT(fitsimage.IsOk());

  /*
  FitsArray stat(array.IsColour() ? array.Plane(1) : array);
  tone.Setup(stat.Npixels(),stat.PixelData());
  */
  tone = FitsTone(array);

  // Init canvas
  //  shrink = 1;
  InitCanvas();
  //  rendering = true;
  shrinking = true;
  tunning = true;
}

void MuniDisplayCanvas::SetInitShrink(int s)
{
  shrink = s;
  InitCanvas();
}

void MuniDisplayCanvas::InitCanvas()
{
  canvasid = 007;
  //  shrink = 1;
  //  zoom = 1;
  canvas = icon.Scale(fitsimage.GetWidth()/shrink,fitsimage.GetHeight()/shrink);
  //  SetVirtualSize(canvas.GetWidth(),canvas.GetHeight());
  //  UpdateScroll();
}

void MuniDisplayCanvas::UpdateCanvas(int shrink, int zoom)
{
  if( canvas.IsOk() ) {
    int w = fitsimage.GetWidth();
    int h = fitsimage.GetHeight();
    if( shrink >= 1 )
      canvas.Rescale(w/shrink,h/shrink);
  }
  //  SetVirtualSize(zoom*canvas.GetWidth(),zoom*canvas.GetHeight());

  canvasid = canvasid + 1;
}

void MuniDisplayCanvas::UpdateCanvas()
{
  canvasid = canvasid + 1;
}


void MuniDisplayCanvas::Reset()
{
  StopRendering();

  //  itt.Reset();
  //  pal.Reset();
  //  color.Reset();
  tone.Reset();
  pal.SetPalette(config->display_pal);

  InitCanvas();
  //  rendering = true;
  invoking = true;
  shrinking = true;
  tunning = true;
}

void MuniDisplayCanvas::StartRendering()
{
  wxASSERT(render == 0);

  //  ::wxBeginBusyCursor();

  if( shrinking ) {
    completed = false;
    render = new DisplayShrinkRender(this,fitsimage,tone,shrink,canvasid);
    wxThreadError code = render->Create();
    wxASSERT(code == wxTHREAD_NO_ERROR);
    render->Run();
    shrinking = false;
  }
  else if( tunning ) {
    wxASSERT(scaled.IsOk());
    render = new DisplayTuneRender(this,scaled,tone,itt,pal,color,canvasid);
    wxThreadError code = render->Create();
    wxASSERT(code == wxTHREAD_NO_ERROR);
    render->Run();
    tunning = false;
  }
}


void MuniDisplayCanvas::OnSubRender(MuniRenderEvent& event)
{
  if( event.id != canvasid || ! event.picture.IsOk() )
    return;

  wxImage img(event.picture.GetWidth(),event.picture.GetHeight(),
  	      event.picture.NewTopsyTurvyRGB());

  if( img.IsOk() && canvas.IsOk() ) {
    /*
    wxLogDebug("MuniDisplayCanvas::OnSubRender %d %d %d %d",
	       (int)event.x,(int)event.y,
	       (int)event.x + img.GetWidth(),
	       (int)event.y + img.GetHeight());
    wxASSERT(event.x >= 0 &&  event.y >= 0 &&
	     event.x + img.GetWidth() <= canvas.GetWidth() &&
	     event.y + img.GetHeight() <= canvas.GetHeight() );
    */
    if( event.x >= 0 &&  event.y >= 0 &&
	     event.x + img.GetWidth() <= canvas.GetWidth() &&
	event.y + img.GetHeight() <= canvas.GetHeight() )
      canvas.Paste(img,event.x,event.y);
  }

  Refresh();
}

void MuniDisplayCanvas::OnRenderFinish(MuniRenderEvent& event)
{
  wxLogDebug("MuniDisplayCanvas::OnRenderFinish %d %d %p",event.id,event.GetId(),this);

  finished = true;

  if( event.image.IsOk() && event.GetId() == ID_RENDER_SHRINK &&
      event.completed ) {
    scaled = event.image;
    completed = event.completed;
  }

  if( event.GetId() == ID_RENDER_SHRINK && scaled.IsOk() ) {
    rendering = true;
    tunning = true;
  }

  if( event.completed )
    wxQueueEvent(GetParent(),event.Clone());

}

void MuniDisplayCanvas::InvokeRendering()
{
  invoking = true;
}

void MuniDisplayCanvas::StopRendering()
{
  wxLogDebug("MuniDisplayCanvas::StopRendering()");
  {
    wxCriticalSectionLocker enter(renderCS);
    if( render )
      render->Delete();
  }

  while(true) {
    {
      wxCriticalSectionLocker enter(renderCS);
      if( ! render ) break;
    }
    ::wxMilliSleep(1);
  }

  rendering = false;
  finished = true;
}

void MuniDisplayCanvas::OnDraw(MuniDrawEvent& event)
{
  //  wxLogDebug("MuniDisplayCanvas::OnDraw");
  MuniLayer layer(event.layer);

  // remove all layers identified by MuniLayer::id
  list<MuniLayer>::iterator i = layers.begin();
  while( i != layers.end() ) {
    if( i->GetId() == event.layer.GetId() )
      i = layers.erase(i);
    else
      ++i;
  }
  layers.push_back(event.layer);

  Refresh();
}

void MuniDisplayCanvas::ConvertCoo(double u, double v, double& x, double& y)
{
  if( zoom > 1 ) {
    x = xoff + zoom*u;
    y = yoff + (canvas.GetHeight() - v)*zoom;
    x = x - vbX;
    y = y - vbY;
  }
  else {
    x = xoff + u / shrink;
    y = yoff - v / shrink + canvas.GetHeight();
    x = x - vbX;
    y = y - vbY;
  }
}

void MuniDisplayCanvas::DrawLayer(wxPaintDC& dc, const MuniLayer& layer)
{
  //  wxLogDebug("DrawLayer");
  //  wxStopWatch sw;
  //DoPrepareDC(dc);

  double height = canvas.GetHeight();

  // Find Out where the window is scrolled to
  //  vbX,vbY;                     // Top left corner of client
  //  GetViewStart(&vbX,&vbY);
  /*
  vbX = 0;
  vbY = 0;
  */
  //  wxLogDebug("%d",(int)layer.objects.size());

  //  wxColour blue(wxBLUE->Red(),wxBLUE->Green(),wxBLUE->Blue(),50);

  /*
  dc.SetBrush(wxBrush(blue));
  dc.SetPen(wxPen(*wxYELLOW,1));

  for(size_t i = 0; i < layer.objects.size(); i++) {

    vector<wxObject *>::const_iterator i;
    for(i = layer.objects.begin(); i != layer.objects.end(); ++i) {

      MuniDrawCircle *circle = dynamic_cast<MuniDrawCircle *>(*i);
      if( circle ) {
	float d = circle->r*zoom;
	float x = xoff + zoom*circle->x;
	float y = yoff + (height - zoom*circle->y);
	dc.DrawCircle(x,y,d);
      }

      MuniDrawLine *line = dynamic_cast<MuniDrawLine *>(*i);
      if( line ) {
	wxPoint2DDouble points[2];
	points[0] = wxPoint2DDouble(zoom*line->x1,height-zoom*line->y1);
	points[1] = wxPoint2DDouble(zoom*line->x2,height-zoom*line->y2);
	dc.DrawLine(zoom*line->x1,height-zoom*line->y1,
		    zoom*line->x2,height-zoom*line->y2);
      }

      MuniDrawText *text = dynamic_cast<MuniDrawText *>(*i);
      if( text ) {
	dc.DrawText(text->text,zoom*text->x,height-zoom*text->y);
      }
    }

  }
  */

  //  wxLogDebug("%f %f %f %f %f %f ",(float)xoff,(float)yoff,(float)zoom,(float)vbX,(float)vbY,height);

  wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
  if( gc ) {
    gc->SetAntialiasMode(wxANTIALIAS_NONE);
    gc->BeginLayer(1.0);

    // default parameters
    gc->SetFont(*wxNORMAL_FONT,*wxWHITE);
    gc->SetPen(wxPen(*wxWHITE,1));
    gc->SetBrush(*wxBLUE_BRUSH);

    // switch-off in fast-update mode?
    gc->SetAntialiasMode(wxANTIALIAS_DEFAULT);

    //    gc->SetBrush(wxBrush(wxColour(90,90,255,130)));
    //    gc->SetBrush(*wxBLUE_BRUSH);
    //    gc->SetPen(wxPen(*wxWHITE,1));
    //    gc->DrawEllipse(100.0,100.0,100.0,100.0);


    vector<wxObject *>::const_iterator i;
    vector<wxObject *> objects = layer.GetObjects();
    for(i = objects.begin(); i != objects.end(); ++i) {

      MuniDrawPen *pen = dynamic_cast<MuniDrawPen *>(*i);
      if( pen )
	gc->SetPen(pen->pen);
	//	gc->SetPen(wxPen(pen->colour,pen->width));

      MuniDrawBrush *brush = dynamic_cast<MuniDrawBrush *>(*i);
      if( brush )
	gc->SetBrush(brush->brush);

      MuniDrawCircle *circle = dynamic_cast<MuniDrawCircle *>(*i);
      if( circle ) {
	double r = wxMax(circle->r * zoom,2);
	double d = 2*r;

	double x,y;
	ConvertCoo(circle->x,circle->y,x,y);
	//wxLogDebug("%f %f %f %f",x,y,circle->r,d);
	x = x - r;
	y = y - r;
	gc->DrawEllipse(x,y,d,d);
      }

      MuniDrawEllipse *ellipse = dynamic_cast<MuniDrawEllipse *>(*i);
      if( ellipse ) {
	double a = ellipse->a;
	double e = ellipse->e;
	double b = a*sqrt(1 - e*e);
	double x,y;
	//	wxGraphicsMatrix m = gc->GetTransform();
	ConvertCoo(ellipse->x,ellipse->y,x,y);
	gc->Translate(x,y);
	gc->Rotate(ellipse->i / 57.3);
	gc->DrawEllipse(-a,-b,2*a,2*b);
	gc->Rotate(-ellipse->i / 57.3);
	gc->Translate(-x,-y);
      }

      MuniDrawRing *ring = dynamic_cast<MuniDrawRing *>(*i);
      if( ring ) {
	double rin = ring->rin;
	double rout = ring->rout;
	double x,y;
	ConvertCoo(ring->x,ring->y,x,y);
	wxGraphicsPath gpath = gc->CreatePath();
	gpath.AddEllipse(x-rin,y-rin,2*rin,2*rin);
	gpath.AddEllipse(x-rout,y-rout,2*rout,2*rout);
	gc->FillPath(gpath);

	/*
	double e = ring->e;
	double q = sqrt(1 - e*e);
	double x,y,a,b;
	ConvertCoo(ring->x,ring->y,x,y);
	gc->Translate(x,y);
	gc->Rotate(-ring->i / 57.3);
	wxGraphicsPath gpath = gc->CreatePath();
	a = rin;
	b = q*rin;
	gpath.AddEllipse(-a,-b,2*a,2*b);
	a = rout;
	b = q*rout;
	gpath.AddEllipse(-a,-b,2*a,2*b);
	gc->FillPath(gpath);
	gc->Rotate(ring->i / 57.3);
	gc->Translate(-x,-y);
	*/
      }


      MuniDrawLine *line = dynamic_cast<MuniDrawLine *>(*i);
      if( line ) {
	wxPoint2DDouble points[2];
	double x,y;
	ConvertCoo(line->x1,line->y1,x,y);
	points[0] = wxPoint2DDouble(x,y);
	ConvertCoo(line->x2,line->y2,x,y);
	points[1] = wxPoint2DDouble(x,y);
	//	wxLogDebug("%f %f %f %f %f %f",line->x1,line->y1,line->x2,line->y2,x,y);
	/*
	points[0] = wxPoint2DDouble(xoff+zoom*line->x1,yoff+height-zoom*line->y1);
	points[1] = wxPoint2DDouble(xoff+zoom*line->x2,yoff+height-zoom*line->y2);
	*/
	gc->DrawLines(2,points);
      }

      MuniDrawCross *cross = dynamic_cast<MuniDrawCross *>(*i);
      if( cross ) {
	wxPoint2DDouble points[2];
	double x,y;
	ConvertCoo(cross->x,cross->y,x,y);
	points[0] = wxPoint2DDouble(x-cross->r,y);
	points[1] = wxPoint2DDouble(x+cross->r,y);
	gc->DrawLines(2,points);
	points[0] = wxPoint2DDouble(x,y-cross->r);
	points[1] = wxPoint2DDouble(x,y+cross->r);
	gc->DrawLines(2,points);
	/*
	points[0] = wxPoint2DDouble(xoff+zoom*(cross->x-cross->r),yoff+height-zoom*cross->y);
	points[1] = wxPoint2DDouble(xoff+zoom*(cross->x+cross->r),yoff+height-zoom*cross->y);
	gc->DrawLines(2,points);
	points[0] = wxPoint2DDouble(xoff+zoom*cross->x,yoff+height-zoom*(cross->y-cross->r));
	points[1] = wxPoint2DDouble(xoff+zoom*cross->x,yoff+height-zoom*(cross->y+cross->r));
	gc->DrawLines(2,points);
	gc->DrawLines(2,points);
	*/
      }

      MuniDrawFont *font = dynamic_cast<MuniDrawFont *>(*i);
      if( font )
	gc->SetFont(font->font,font->colour);

      MuniDrawText *text = dynamic_cast<MuniDrawText *>(*i);
      if( text ) {
	double w,h,d,l;
	double x = text->x > 0.0 ? text->x : text->x;
	double y = text->y > 0.0 ? text->y : text->y;
	gc->GetTextExtent(text->text,&w,&h,&d,&l);
	if( text->angle < 0.01 )
	  gc->DrawText(text->text,xoff+zoom*(x-w/2.0),yoff+height-zoom*y);
	else if( fabs(text->angle - 1.5708) < 0.01 )
	  gc->DrawText(text->text,xoff+zoom*(x-h),yoff+height-zoom*(y - w/2.0),text->angle);
      }

      MuniDrawRectangle *rect = dynamic_cast<MuniDrawRectangle *>(*i);
      if( rect ) {
	double x,y;
	ConvertCoo(rect->x,rect->y-rect->h,x,y);
	gc->DrawRectangle(x,y,rect->w,rect->h);
	//	gc->DrawRectangle(xoff+zoom*rect->x,yoff+zoom*(height-rect->y-rect->h),rect->w,rect->h);
      }

      MuniDrawBitmap *bitmap = dynamic_cast<MuniDrawBitmap *>(*i);
      if( bitmap ) {
	double x,y;
	ConvertCoo(bitmap->x,bitmap->y-bitmap->h,x,y);
	gc->DrawBitmap(bitmap->bitmap,x,y,bitmap->w,bitmap->h);

	//	gc->DrawBitmap(bitmap->bitmap,xoff+bitmap->x,yoff+(height-bitmap->y-bitmap->h),bitmap->w,bitmap->h);
      }

    }

    gc->EndLayer();
    delete gc;
  }


  //  wxLogDebug("Drawing took: %f sec",sw.Time()/1000.0);
}

// void MuniDisplayCanvas::DrawSVG(wxPaintDC& dc, const wxXmlDocument& svg)
// {
//   if( !(svg.IsOk() && svg.GetRoot()->GetName() == "svg") ) return;

//   wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
//   if( gc ) {

//     // default colors
//     gc->SetBrush(wxBrush(*wxBLUE));
//     gc->SetPen(wxPen(*wxYELLOW,1));
//     gc->SetFont(*wxNORMAL_FONT,*wxWHITE);

//     wxXmlNode *child = svg.GetRoot()->GetChildren();
//     while (child) {

//       if (child->GetName() == "line") {

// 	double x1,y1,x2,y2;

// 	wxXmlAttribute *prop = child->GetAttributes();
// 	while(prop) {

// 	  wxString a = prop->GetValue();
// 	  double x;
// 	  a.ToDouble(&x);

// 	  if( prop->GetName() == "x1" )
// 	    x1 = x;
// 	  else if( prop->GetName() == "y1" )
// 	    y1 = x;
// 	  else if( prop->GetName() == "x2" )
// 	    x2 = x;
// 	  else if( prop->GetName() == "y2" )
// 	    y2 = x;

// 	  prop = prop->GetNext();
// 	}
// 	gc->StrokeLine(x1,y1,x2,y2);

//       }

//       else if( child->GetName() == "text") {

// 	double x,y;
// 	wxString text(child->GetNodeContent());

// 	wxXmlAttribute *prop = child->GetAttributes();
// 	while(prop) {

// 	  wxString a = prop->GetValue();
// 	  double t;
// 	  a.ToDouble(&t);

// 	  if( prop->GetName() == "x" )
// 	    x = t;
// 	  else if( prop->GetName() == "y" )
// 	    y = t;

// 	  prop = prop->GetNext();
// 	}

// 	gc->DrawText(text,x,y);
//       }

//       else if( child->GetName() == "circle") {

// 	double cx,cy,r;

// 	wxXmlAttribute *prop = child->GetAttributes();
// 	while(prop) {

// 	  wxString a = prop->GetValue();
// 	  double x;
// 	  a.ToDouble(&x);

// 	  if( prop->GetName() == "cx" )
// 	    cx = x;
// 	  else if( prop->GetName() == "cy" )
// 	    cy = x;
// 	  else if( prop->GetName() == "r" )
// 	    r = x;

// 	  prop = prop->GetNext();
// 	}

// 	gc->DrawEllipse(cx-r/2.0,cy-r/2.0,r,r);
//       }
//       child = child->GetNext();
//     }
//     delete gc;
//   }

// }

wxPoint MuniDisplayCanvas::GetCartesianPosition(const wxPoint& point) const
{
  int height = fitsimage.GetHeight();

  wxPoint cursor = point - wxPoint(xoff,yoff);

  int x = cursor.x;
  int y = cursor.y;
  if( zoom > 1 ) {
    x = x / zoom;
    y = height - y / zoom;
  }
  else if( shrink >= 1 ) {
    x = shrink * x;
    y = height - shrink * y;
  }

  return wxPoint(x,y);
}

void MuniDisplayCanvas::OnMouseMotion(wxMouseEvent& event)
{
  if( ! canvas.IsOk()) return;
  if( astrometry ) return;

  wxClientDC dc(this);
  //  DoPrepareDC(dc);

  if( event.Moving() ) {

    MuniSlewEvent evt(EVT_SLEW);
    crosshair = GetCartesianPosition(event.GetLogicalPosition(dc));
    evt.x = crosshair.x;
    evt.y = crosshair.y;

    wxRect rect(0,0,fitsimage.GetWidth(),fitsimage.GetHeight());
    evt.inside = rect.Contains(crosshair);

    wxQueueEvent(GetParent(),evt.Clone());
  }
  else if( event.Dragging() ) {


    if( ! dragging ) {

      SetCursor(wxCursor(wxCURSOR_HAND));
      xdrag0 = event.GetLogicalPosition(dc).x;
      ydrag0 = event.GetLogicalPosition(dc).y;
      dragging = true;

    }
    else if( dragging ) {

      int x0,y0,dx,dy;
      //      GetViewStart(&x0,&y0);
      //      GetScrollPixelsPerUnit(&dx,&dy);
      if( dx > 0  && dy > 0 ) {
	int i = x0 - (event.GetLogicalPosition(dc).x - xdrag0) / dx;
	int j = y0 - (event.GetLogicalPosition(dc).y - ydrag0) / dy;
	//	Scroll(i,j);
      }
    }

  }

}

void MuniDisplayCanvas::OnMouseEnter(wxMouseEvent& event)
{
  MuniSlewEvent evt(EVT_SLEW);
  evt.entering = true;
  wxQueueEvent(GetParent(),evt.Clone());
}

 void MuniDisplayCanvas::OnMouseLeave(wxMouseEvent& event)
{
  MuniSlewEvent evt(EVT_SLEW);
  evt.leaving = true;
  wxQueueEvent(GetParent(),evt.Clone());
}

void MuniDisplayCanvas::OnMouseWheel(wxMouseEvent& event)
{
  /*
  int id = event.GetWheelRotation() > 0 ? wxID_ZOOM_IN : wxID_ZOOM_OUT ;
  wxQueueEvent(topwin,new wxCommandEvent(wxEVT_COMMAND_MENU_SELECTED,id));
  */
}

void MuniDisplayCanvas::OnClick(wxMouseEvent& event)
{
  wxClientDC dc(this);
  wxPoint p = GetCartesianPosition(event.GetLogicalPosition(dc));
  MuniClickEvent e(EVT_CLICK);
  e.x = p.x;
  e.y = p.y;
  wxQueueEvent(GetParent(),e.Clone());
}

void MuniDisplayCanvas::OnLeftUp(wxMouseEvent& event)
{
  if( dragging ) {
    SetCursor(*wxSTANDARD_CURSOR);
    dragging = false;
  }
}

void MuniDisplayCanvas::OnMenu(wxMouseEvent& event)
{
  wxTopLevelWindow *twin = static_cast<wxTopLevelWindow *>(topwin);
  if( !twin ) return;

  wxMenu popup;

  if( twin->IsFullScreen() ) {
    popup.Append(ID_MENU_FULLSCREEN,"Leave Fullscreen");
  }
  else /*if( ! twin->IsFullScreen() )*/ {
    popup.Append(ID_CLIP_VALUE,"Copy value");
    popup.Append(ID_CLIP_COO,"Copy coordinates");
  }

  PopupMenu(&popup);
}

void MuniDisplayCanvas::OnClipValue(wxCommandEvent& event)
{
  wxPoint p(crosshair);
  FitsValue value(array,config->phsystemfile, config->fits_key_area,
		    config->fits_key_exptime, config->fits_key_filter);
  value.SetType(value_type);

  wxString val;
  if( array.IsColour() ) {
    /*
    val =       "X = " + value.Get_str(p.x,p.y,2);
    val = val + " Y = " + value.Get_str(p.x,p.y,1);
    val = val + " Z = " + value.Get_str(p.x,p.y,0);
    */
  }
  else
    val = value.Get_str(p.x,p.y);

  if(wxTheClipboard->Open() ) {
    wxTheClipboard->SetData(new wxTextDataObject(val));
    wxTheClipboard->Close();
  }
}

void MuniDisplayCanvas::OnClipCoo(wxCommandEvent& event)
{
  wxPoint p(crosshair);
  FitsCoo coords(array);
  coords.SetType(coo_type);

  wxString xcoo_str,ycoo_str;
  if( coo_type == COO_PIXEL )
    coords.GetPix(p.x,p.y,xcoo_str,ycoo_str);
  else
    coords.GetCoo(p.x,p.y,xcoo_str,ycoo_str);


  if(wxTheClipboard->Open() ) {
    wxTheClipboard->SetData(new wxTextDataObject(xcoo_str+" "+ycoo_str));
    wxTheClipboard->Close();
  }
}

void MuniDisplayCanvas::OnValueType(wxCommandEvent& event)
{
  wxASSERT(event.GetId() == ID_VALTYPE);
  value_type = event.GetInt();
}

void MuniDisplayCanvas::OnCooType(wxCommandEvent& event)
{
  wxASSERT(event.GetId() == ID_COOTYPE);
  coo_type = event.GetInt();
}

void MuniDisplayCanvas::OnKeyDown(wxKeyEvent& event)
{
  long keycode = event.GetKeyCode();

  wxTopLevelWindow *twin = static_cast<wxTopLevelWindow *>(topwin);

  if( twin && twin->IsFullScreen() && keycode == WXK_ESCAPE ){
    wxQueueEvent(topwin,new wxCommandEvent(wxEVT_COMMAND_MENU_SELECTED,
					   ID_FULLSCREEN));
    return;
  }

  if( keycode == WXK_ESCAPE ) {
    StopRendering();
    return;
  }

  if( event.GetModifiers() & wxMOD_CMD ) {
    int id = wxID_ANY;
    if( keycode == WXK_NUMPAD_MULTIPLY )
      ;
    else {
      switch(keycode) {
      case WXK_LEFT:            id = wxID_BACKWARD; break;
      case WXK_RIGHT:           id = wxID_FORWARD;  break;
      }
    }

    if( id != wxID_ANY ) {
      wxQueueEvent(topwin,new wxCommandEvent(wxEVT_COMMAND_MENU_SELECTED,id));
      return;
    }
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

void MuniDisplayCanvas::UpdateTune()
{
  //rendering = true;
  invoking = true;
  tunning = true;
  UpdateCanvas();
}

void MuniDisplayCanvas::UpdateScale()
{
  //  rendering = true;
  invoking = true;
  //  scaling = false;
  shrinking = true;
  tunning = true;
  UpdateCanvas();
}

void MuniDisplayCanvas::UpdateScale(int s)
{
  shrinking = true;
  int oldshrink = shrink;

  if( s > 0 ) {
    shrink = s;
    zoom = 1;
    wxLogDebug("MuniDisplayCanvas::UpdateScale shrink: 1:%d",shrink);
  }
  else {
    shrink = 1;
    zoom = abs(s);
    wxLogDebug("MuniDisplayCanvas::UpdateScale zoom: %d:1",zoom);
  }

  //  rendering = oldshrink != 1;
  invoking = oldshrink != 1;

  if( invoking ) {
    shrinking = true;
    tunning = true;
  }

  wxLogDebug("MuniDisplayCanvas::UpdateScale: %d %d %d %d %d",
	     oldshrink,zoom,completed,rendering,zooming);

  UpdateCanvas(shrink,zoom);
  //  UpdateScroll();
  Refresh();
}

/*
void MuniDisplayCanvas::UpdateScroll()
{
  if( ! canvas.IsOk() ) return;

  int wc,hc;
  GetClientSize(&wc,&hc);
  int w = zoom*canvas.GetWidth();
  int h = zoom*canvas.GetHeight();

  int dw = -1;
  int dh = -1;
  if( w > wc )
    dw = zoom;
  if( h > hc )
    dh = zoom;

  if( dw > 0 || dh > 0 )
    SetScrollRate(dw,dh);

  xoff = max((wc - w)/2,0);
  yoff = max((hc - h)/2,0);

  Refresh();
}
*/

void MuniDisplayCanvas::OnTuneFine(MuniTuneEvent& e)
{
  //  wxLogDebug("MuniDisplayCanvas::OnTuneFine %d",e.GetId()==ID_COLOR_NITEVISION);

  bool found = false;

  switch(e.GetId()) {

    // fine tune
  case ID_ITT_TYPE:   itt.SetItt(e.GetString()); found = true; break;
    //  case ID_ITT_AMP:    itt.SetAmp(e.x); found = true; break;
    //  case ID_ITT_ZERO:   itt.SetZero(e.x); found = true; break;

  case ID_COLOUR_NITEVISION: color.SetNiteVision(e.GetInt()); found = true;break;
  case ID_COLOUR_SATUR:   color.SetSaturation(e.x); found = true; break;
  case ID_COLOUR_HUE:     color.SetHue(e.x); found = true; break;
    //  case ID_COLOUR_NITETHRESH: color.SetNiteThresh(e.x); found = true; break;
    //  case ID_COLOUR_NITEWIDTH:  color.SetNiteWidth(e.x); found = true; break;
  case ID_COLOUR_MEAN:    color.SetLevel(e.index,e.x); found = true; break;
  case ID_COLOUR_WEIGHT:  color.SetWeight(e.index,e.x); found = true; break;

  case ID_PALETTE_TYPE:     pal.SetPalette(e.GetString()); found = true; break;
  case ID_PALETTE_INVERSE: pal.SetInverse(e.GetInt()); found = true; break;
  }
  if( found ) {
    UpdateTune();
    return;
  }

  found = false;
  switch(e.GetId()) {

    // scaling
  case ID_TONE_BLACK: tone.SetBlack(e.x); found = true; break;
  case ID_TONE_SENSE:  tone.SetSense(e.x); found = true; break;
  case ID_TONE_QBLACK: tone.SetQblack(e.x); found = true; break;
  case ID_TONE_RSENSE:  tone.SetRsense(e.x); found = true; break;
  }
  if( found ) {
    UpdateScale();
    return;
  }

  // shrink or zoom
  if( e.GetId() == ID_ZOOM_SCALE ) {
    UpdateScale(e.n);
    return;
  }

  if( e.GetId() == ID_RESET ) {

    Reset();
    return;
  }

  //  if( e.GetId() == ID_ZOOM ) // ignored
  //    return;

  wxLogDebug("WARNING: MuniDisplay::OnTuneFine - reached unknown ID");
}

void MuniDisplayCanvas::OnAstrometry(MuniAstrometryEvent& e)
{
  //  wxLogDebug("MuniDisplayCanvas::OnAstrometry");

  RemoveLayers(ID_ASTROMETRY);

  /*
  list<MuniLayer>::iterator i;
  for( i = layers.begin(); i != layers.end(); ++i)
    if( i->GetId() == ID_ASTROMETRY ) {
      layers.erase(i);
      break;
    }
  */

  astrometry = e.astrometry;
  layers.clear();

  if( astrometry ) {

    //wxLogDebug(e.proj+" %f %f %f",e.xcen,e.ycen,e.scale);

    aangle = e.angle;
    ascale = e.scale;
    acen = e.acen;
    dcen = e.dcen;
    xcen = e.xcen;
    ycen = e.ycen;
    //  layers.push_back(e.svg);
    layers.push_back(e.layer);
//    wxLogDebug("astrometryevent: %d",(size_t)e.layer.objects.size());
//    wxLogDebug("astrometryevent: %d",(size_t)e.layer.IsEmpty());
  }

  Refresh();
}

void MuniDisplayCanvas::OnZoom(MuniZoomEvent& e)
{
  if( e.GetId() == ID_ZOOMER ) {

    double zoom = e.zoom;
    if( abs(zoom - 1) < 0.01 )
      UpdateScale(1);
    else if( zoom > 1.01 )
      UpdateScale(-int(zoom+0.5));
    else // if zoom < 0.99
      UpdateScale(int(1.0/zoom + 0.5));
  }
  else if ( e.GetId() == ID_ZOOMVIEW ) {
    double width  = fitsimage.GetWidth();
    double height = fitsimage.GetHeight();

    xsubwin = int(e.x*width/shrink + 0.5);
    ysubwin = int(e.y*height/shrink+ 0.5);

    Refresh();
  }


}

void MuniDisplayCanvas::OnPhotometry(MuniPhotometryEvent& e)
{
  wxLogDebug("MuniDisplayCanvas::OnPhotometry %d",(int)e.erase);

  if( e.erase )
    RemoveLayers(ID_PHOTOMETRY);
  else
    layers.push_back(e.layer);

  Refresh();
}

void MuniDisplayCanvas::AddLayer(const MuniLayer& ml)
{
  layers.push_back(ml);
  Refresh();
}

void MuniDisplayCanvas::RemoveLayers(int id)
{
  list<MuniLayer>::iterator i = layers.begin();
  while( i != layers.end() ) {
    if( i->GetId() == id )
      i = layers.erase(i);
    else
      ++i;
  }

  Refresh();
}


wxImage MuniDisplayCanvas::GetImage() const
{
  /*
    Returns non-scaled or non-shrink size of frame.
    That's a placeholder, the proper rendering is necessary.

  */

  if( shrink == 1 )
    return canvas;
  else
    return canvas.Scale(array.GetWidth(),array.GetHeight());
}



/*
void MuniDisplayCanvas::OverlayBitmap(const vector<MuniDrawBase *>& drawobj)
{
  wxASSERT(bitmap.IsOk());

  wxBitmap layer = wxBitmap(bitmap.GetWidth(),bitmap.GetHeight());

  wxMemoryDC mdc(layer);
  mdc.SetBackground(*wxBLACK_BRUSH);
  mdc.Clear();
  int height = layer.GetHeight();

  wxGraphicsContext *gc = wxGraphicsContext::Create(mdc);
  if( gc ) {

    gc->SetBrush(wxBrush(*wxBLUE));
    gc->SetPen(wxPen(*wxYELLOW,1));

    vector<MuniDrawBase *>::const_iterator i;
    for(i = drawobj.begin(); i != drawobj.end(); ++i) {

      MuniDrawCircle *e = static_cast<MuniDrawCircle *>(*i);
      if( e ) {
	float r = e->r;
	float x = e->cx - r/2.0;
	float y = height - (e->cy + r/2.0);
	gc->DrawEllipse(x,y,r,r);
      }

    }
    delete gc;
  }
  mdc.SelectObject(wxNullBitmap);

  wxImage i = layer.ConvertToImage();
  i.SetMaskColour(0,0,0);

  overlay.push_back(wxBitmap(i));

}
*/

 //void MuniDisplayCanvas::OverlayGrid(const vector<wxGraphicsPath>& paths)
 //{
  // wxASSERT(bitmap.IsOk());

  // wxBitmap layer = wxBitmap(bitmap.GetWidth(),bitmap.GetHeight());

  // wxMemoryDC mdc(layer);
  // mdc.SetBackground(*wxBLACK_BRUSH);
  // mdc.Clear();
  // int height = layer.GetHeight();

  // wxGraphicsContext *gc = wxGraphicsContext::Create(mdc);
  // if( gc ) {

  //   gc->SetBrush(wxBrush(*wxBLUE));
  //   gc->SetPen(wxPen(*wxYELLOW,1,wxPENSTYLE_SOLID));

  //   wxGraphicsPath path = gc->CreatePath();

  //   path.MoveToPoint(333,333);
  //   path.AddLineToPoint(500,500);

  //   /*
  //   vector<wxGraphicsPath>::const_iterator i;
  //   for(i = paths.begin(); i != paths.end(); ++i) {
  //   */

  //   gc->StrokePath(path);

  // }
  // delete gc;
  // mdc.SelectObject(wxNullBitmap);

  // wxImage i = layer.ConvertToImage();
  // i.SetMaskColour(0,0,0);

  // overlay.push_back(wxBitmap(i));

 //}
