/*

  xmunipack - image rendering in threads

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

#include "render.h"
#include "fits.h"
#include <wx/wx.h>
#include <wx/thread.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;


// ---   MuniThreadRender

MuniThreadRender::MuniThreadRender(wxEvtHandler *eh, wxMutex *mx, wxCondition *c,
				   int *cnt, const FitsArray& a,
				   const FitsGeometry *g, const FitsDisplay *d,
				   int id):
  wxThread(wxTHREAD_JOINABLE), // important: one'll be destroyed in std::vector
  handler(eh), kerberos(mx), condition(c), count(cnt), stop(false),
  cdata_valid(false),
  zoom(-1.0), shrink(-1), iwidth(0), iheight(0), cwidth(0), cheight(0), cdepth(0),
  rid(id), ntiles(0), fwidth(int(a.GetWidth())), fheight(int(a.GetHeight())),
  fdepth(int(a.GetDepth())), fdata(a.PixelData()), idata(0), cdata(0),
  geometry(g),display(d)
{
  wxASSERT(handler);
  //  wxLogDebug("MuniThreadRender::MuniThreadRender");
}

wxThread::ExitCode MuniThreadRender::Entry()
{
  wxASSERT(idata && fdata && shrink > 0);

  int z = zoom > 1.01 ? int(zoom + 0.5) : 1;

  for(size_t n = 0; n < rects.size() && stop == false; n++) {

    /*
    wxLogDebug("MuniThreadRender::Entry() rect[%ld]: %d %d %d %d, %d %.3f, %d %d",
	       n,rects[n].x,rects[n].y,rects[n].width,rects[n].height,shrink,zoom,
	       int(exposed.x),int(exposed.y));
    */
    int x = rects[n].x;
    int y = rects[n].y;
    int w = rects[n].width;
    int h = rects[n].height;
    int xoff = exposed.x;
    int yoff = exposed.y;

    if( w >= shrink && h >= shrink ) {

      int xs = x / shrink;
      int ys = y / shrink;
      int ws = w / shrink;
      int hs = h / shrink;
      int xoffs = xoff / shrink;
      int yoffs = yoff / shrink;

      if( shrink > 1 ) {

	// we begin with shrinking into the temporary array
	if( cdata && cdata_valid == false )
	  geometry->ShrinkSubwin(shrink,x,y,w,h,cwidth,cheight,cdata);
	wxASSERT(cdata && cwidth >= iwidth && cheight >= iheight);

	// the shrinked image is rendered
	display->GetRGB(cwidth,cheight,cdepth,cdata,xs,ys,ws,hs,xoffs,yoffs,
			iwidth,iheight,idata);
      }
      else if( z == 1 )
	display->GetRGB(fwidth,fheight,fdepth,fdata,x,y,w,h,xoff,yoff,
			iwidth,iheight,idata);

      else if( z > 1 ) {
	int swidth = w;
	int sheight = h;
	unsigned char *sdata = new unsigned char[3*w*h];
	display->GetRGB(fwidth,fheight,fdepth,fdata,x,y,w,h,x,y,
			swidth,sheight,sdata);
	Zoom(swidth,sheight,sdata,z,x-xoff,y-yoff,w,h,iwidth,iheight,idata);
	delete[] sdata;

      }

      // all computations are finished, notify calling thread
      MuniRenderEvent ev(EVT_RENDER,ID_SUBRENDER);
      // coordinates relative to idata, the topsy turvy image
      ev.rid = rid;
      ev.x = int((x - exposed.x)*zoom);
      ev.y = iheight - int((h + y - exposed.y)*zoom);
      ev.w = int(w*zoom);
      ev.h = int(h*zoom);
      wxQueueEvent(handler,ev.Clone());
      //      Sleep(250);

      ntiles++;
    }
  }

  //  wxLogDebug("MuniThreadRender::Entry() at finish %d %d",*count,ntiles);

  kerberos->Lock();
  *count = *count - 1;
  if( *count == 0 ) {
    condition->Signal();
    wxLogDebug("condition->Signal()");
  }
  kerberos->Unlock();

  return (wxThread::ExitCode) stop;
}

void MuniThreadRender::Zoom(int swidth, int sheight, const unsigned char *sdata,
			    int zoom, int x,int y,int w, int h,
			    int iwidth, int iheight, unsigned char *idata)
{
  wxASSERT(swidth > 0 && sheight > 0 && sdata &&
	   iwidth > 0 && iheight > 0 && idata);

  int sstride = 3*swidth;
  int istride = 3*iwidth;
  int y1 = iheight - zoom*(h+y);
  for(int j = 0; j < sheight; j++) {
    int jj = zoom*j;
    const unsigned char *sline = sdata +  j*sstride;
    unsigned char *iline = idata + (y1 + jj)*istride;

    // top line
    for(int i = 0; i < swidth; i++) {
      int ii = zoom*(i + x);
      const unsigned char *spix = sline + 3*i;
      unsigned char *ipix = iline + 3*ii;
      for(int l = 0; l < zoom; l++) {
	unsigned char *lpix = ipix + 3*l;
	memcpy(lpix,spix,3);
      }
    }

    // rest of lines
    for(int l = 1; l < zoom; l++)
      memcpy(iline + l*istride,iline,istride);
  }
}


int MuniThreadRender::GetRenderedTiles() const
{
  return ntiles;
}

void MuniThreadRender::Stop()
{
  stop = true;
}

void MuniThreadRender::SetRectangles(const wxRect& e, const vector<wxRect>& r)
{
  if( r.size() == 0 )
    wxLogDebug("MuniThreadRender::SetRectangles() *** WARNING: Empty area. ****");

  exposed = e;
  rects = r;
}

void MuniThreadRender::SetImageData(int w, int h, unsigned char *d)
{
  wxASSERT(w > 0 && h > 0 && d);
  iwidth = w;
  iheight = h;
  idata = d;
}

void MuniThreadRender::SetZoom(double z)
{
  zoom = z;
  shrink = zoom < 0.99 ? int(1.0/zoom + 0.5) : 1;
}

void MuniThreadRender::SetShrinkArray(int w, int h, int d, float *c, bool a)
{
  wxASSERT(w > 0 && h > 0 && d > 0 && c);
  cwidth = w;
  cheight = h;
  cdepth = d;
  cdata = c;
  cdata_valid = a;
}


// -------   MuniDisplayRender

MuniDisplayRender::MuniDisplayRender(wxEvtHandler *eh, const FitsArray& a):
  handler(eh), condition(kerberos), stop(false), cdata_valid(false), zoom(-1.0),
  ncpu(1), count(0), iwidth(0), iheight(0), cwidth(0), cheight(0), cdepth(0),
  idata(0), cdata(0), fitsimage(a),geometry(a), display(fitsimage.IsColour()),
  signalling(0)
{
  ncpu = std::max(wxThread::GetCPUCount(),1);
  //ncpu = 1; // testing
}


MuniDisplayRender::~MuniDisplayRender()
{
  wxASSERT(renders.size() == 0 && signalling == 0);
  free(cdata);
}

MuniDisplayRender::MuniDisplayRender(const MuniDisplayRender& r):
  condition(kerberos),geometry(fitsimage),display(0)
{
  wxFAIL_MSG("MuniDisplayRender: A COPY CONSTRUCTOR IS INTENTIONALLY LEFT UNIMPLEMENTED");
}

MuniDisplayRender& MuniDisplayRender::operator = (const MuniDisplayRender& r)
{
  wxFAIL_MSG("MuniDisplayRender: AN ASSIGNMENT CONSTRUCTOR IS INTENTIONALLY LEFT UNIMPLEMENTED ");
  return *this;
}

wxThreadError MuniDisplayRender::Run(int id)
{
  wxASSERT(zoom > 0);
  wxASSERT(fitsimage.IsOk() && idata && iwidth > 0 && iheight > 0);
  wxASSERT(renders.size() == 0 && signalling == 0); // all resources should be free

  stop = false;
  kerberos.Lock();

  // init the control counter
  count = ncpu;

  // signalling thread
  signalling = new MuniSignalThread(handler,&condition,id);
  wxThreadError code = signalling->Create();
  if( code != wxTHREAD_NO_ERROR ) {
    wxLogError("Can't create a signalling thread.");
    return code;
  }

  // working threads
  for(int i = 0; i < ncpu; i++) {

    MuniThreadRender *render =
      new MuniThreadRender(handler,&kerberos,&condition,&count,fitsimage,
			   &geometry,&display,id);
    wxThreadError code = render->Create();
    if( code == wxTHREAD_NO_ERROR ) {
      renders.push_back(render);
    }
    else {
      wxLogError("Can't create a thread.");
      delete render;
      for(size_t i = 0; i < renders.size(); i++)
	delete renders[i];
      renders.clear();
      return code;
    }
  }

  // init shrink
  if( zoom < 1.0 && cdata == 0 ) {
    cdata_valid = false;
    cwidth = wxMax(int(zoom*fitsimage.GetWidth()),1);
    cheight = wxMax(int(zoom*fitsimage.GetHeight()),1);
    cdepth = fitsimage.GetDepth();
    cdata = (float *) malloc(cwidth*cheight*cdepth*sizeof(float));
  }
  else if( zoom > 0.99 ) {
    free(cdata);
    cdata = 0;
  }

  // line block
  const vector<wxRect> tiles(Lines());
  ntiles = tiles.size();

  // setup parameters
  for(int i = 0; i < ncpu; i++) {

    // every ncpu tile is used
    vector<wxRect> tales;
    for(int n = i; n < ntiles; n += ncpu )
      tales.push_back(tiles[n]);

    //    wxLogDebug("MuniDisplayRender::Run() i=%d %d",int(i),int(tales.size()));
    renders[i]->SetRectangles(exposed,tales);
    renders[i]->SetImageData(iwidth,iheight,idata);
    renders[i]->SetZoom(zoom);
    if( cdata )
      renders[i]->SetShrinkArray(cwidth,cheight,cdepth,cdata,cdata_valid);
  }

  // signalling
  code = signalling->Run();
  if( code != wxTHREAD_NO_ERROR ) {
    wxLogError("Can't run a signalling thread.");
    return code;
  }

  // Run() is postponed until all threads has set parameters
  for(size_t i = 0; i < renders.size(); i++) {
    wxThreadError code = renders[i]->Run();
    if( code != wxTHREAD_NO_ERROR ) {
      wxLogError("Can't run the thread: #%d.",int(i));
      return code;
    }
  }

  return wxTHREAD_NO_ERROR;
}


void MuniDisplayRender::Wait()
{
  wxLogDebug("MuniDisplayRender::Wait()");

  int ctiles = 0;
  for(size_t i = 0; i < renders.size(); i++) {
    ctiles += renders[i]->GetRenderedTiles();
    renders[i]->Wait();
    delete renders[i];
  }
  renders.clear();

  wxASSERT(signalling);
  signalling->Wait();
  delete signalling;
  signalling = 0;

  if( cdata_valid == false && ! stop )
    cdata_valid = ctiles == ntiles;

  stop = false;
  kerberos.Unlock(); // the last
}

void MuniDisplayRender::Stop()
{
  stop = true;
  signalling->Stop();
  for(size_t i = 0; i < renders.size(); i++)
    renders[i]->Stop();
}

void MuniDisplayRender::StopAndWait()
{
  wxASSERT(signalling);

  signalling->StopAndWait();
  Stop();
  Wait();
}

bool MuniDisplayRender::IsRunning() const
{
  if( signalling ) return true;

  for(size_t i = 0; i < renders.size(); i++) {
    if( renders[i]->IsRunning() )
      return true;
  }

  return false;
}


vector<wxRect> MuniDisplayRender::Lines() const
{
  wxASSERT(fitsimage.IsOk());

  int x = 0;
  int y = 0;
  int width = fitsimage.GetWidth();
  int height = fitsimage.GetHeight();

  if( exposed.GetWidth() > 0 &&  exposed.GetHeight() > 0 ) {
    x = exposed.GetX();
    y = exposed.GetY();
    width = exposed.GetWidth();
    height = exposed.GetHeight();
  }

  int shrink = zoom < 0.99 ? int(1/zoom + 0.5) : 1;
  int n = height / shrink > ncpu ? height / (shrink*ncpu) : 1;
  for(size_t i = 0; i < 10; i++ ) {
    int m = (width*n) / shrink;
    if( m < 8192 ) break;
    n = n / 2;
  }
  int d = shrink*wxMax(n,1);

  vector<wxRect> rects;
  for(int j = 0; j < height; j += d ) {
    int h = std::min(d,height-j);
    if( h > 0 ) {
      wxRect r(x,y+j,width,h);
      rects.push_back(r);
      //      wxLogDebug("%d %d %d %d",x,y+j,width,h);
    }
  }
  wxLogDebug("DisplayRenders::Lines() %d %d %d %d, %d %d %d %d",x,y,width,height,
	     d,int(rects.size()),n,shrink);

  return rects;
}

void MuniDisplayRender::SetImageData(int w, int h, unsigned char *d)
{
  wxASSERT(w > 0 && h > 0 && d);
  iwidth = w;
  iheight = h;
  idata = d;
}

void MuniDisplayRender::SetZoom(double z)
{
  if( fabs(zoom - z) > 0.001 ) {
    free(cdata);
    cdata = 0;
    cdata_valid = false;
    zoom = z;
  }
}


void MuniDisplayRender::SetTone(const FitsTone& tone, const FitsItt& itt,
				const FitsPalette& pal, const FitsColour& colour)
{
  display.SetTone(tone);
  display.SetItt(itt);
  display.SetPalette(pal);
  display.SetColour(colour);
}

void MuniDisplayRender::SetExposed(const wxRect& r)
{
  if( exposed.x != r.x || exposed.y != r.y ||
      exposed.width != r.width || exposed.height != r.height ) {
    cdata_valid = false;
    exposed = r;
  }
}

// ---- MuniSignalThread

MuniSignalThread::MuniSignalThread(wxEvtHandler *eh, wxCondition *c, int id):
  wxThread(wxTHREAD_JOINABLE), handler(eh), condition(c), rid(id),
  stop(false), stop_wait(false)
{
  sw.Start();
}

MuniSignalThread::~MuniSignalThread()
{
  wxLogDebug("Rendering took: %.3f ms",sw.TimeInMicro().ToDouble()/1000.0);
}

void MuniSignalThread::Stop()
{
  stop = true;
}

void MuniSignalThread::StopAndWait()
{
  stop_wait = true;
}

wxThread::ExitCode MuniSignalThread::Entry()
{
  condition->Wait();

  if( stop_wait )
    // no event
    ;
  else if( stop ) {
    // interrupted
    MuniRenderEvent ev(EVT_RENDER,ID_RENDER_INTERRUPT);
    ev.rid = rid;
    wxQueueEvent(handler,ev.Clone());
  }
  else {
    // finish
    MuniRenderEvent ev(EVT_RENDER,ID_RENDER_FINISH);
    ev.rid = rid;
    wxQueueEvent(handler,ev.Clone());
  }
  return (wxThread::ExitCode) 0;
}
