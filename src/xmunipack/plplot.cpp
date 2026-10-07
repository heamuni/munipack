/*

  xmunipack - plplot backend

  Copyright © 2010-2014, 2018-21 F.Hroch (hroch@physics.muni.cz)

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
#include "plot.h"
#include <plplot.h>
#include <wx/graphics.h>
#include <wx/wx.h>
#include <cstring>
#include <cstdlib>

using namespace std;

#define WhiteSmoke 245


PlPlotCairo::PlPlotCairo(wxWindow *w, const wxSize& s):
  wxPanel(w,wxID_ANY,wxDefaultPosition,s), ready(false)
{
  SetBackgroundStyle(wxBG_STYLE_PAINT);
  SetBackgroundColour(*wxWHITE);

  Plotter(true);

  Bind(wxEVT_PAINT,&PlPlotCairo::OnPaint,this);
}

void PlPlotCairo::Update()
{
  Plotter();
  Refresh();
}

void PlPlotCairo::OnPaint(wxPaintEvent& event)
{
  wxASSERT(graph.IsOk());

  wxPaintDC dc(this);
  wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
  if( gc ) {
    wxGraphicsBitmap bmp = gc->CreateBitmapFromImage(graph);
    gc->DrawBitmap(bmp,0,0,graph.GetWidth(),graph.GetHeight());
    delete gc;
  }
}

void PlPlotCairo::Plotter(bool initialise)
{
  ready = true;

  wxSize s = GetSize();
  PLINT width = s.GetWidth();
  PLINT height = s.GetHeight();
  size_t npixels = width*height;
  PLPointer *plotmem = (PLPointer *) malloc(4*npixels);
  unsigned char *smoke = (unsigned char *) plotmem;
  for(size_t n = 0; n < npixels; smoke = smoke + 4, n++) {
    memset(smoke,WhiteSmoke,3);
    memset(smoke+3,255,1);
  }

  plsdev("memcairo");
  plsmema(width,height,plotmem);

  if( initialise )
    InitGraph();
  else
    MakeGraph();

  unsigned char *data = (unsigned char *) malloc(3*npixels);
  unsigned char *alpha= (unsigned char *) malloc(npixels);

  unsigned char *d = data;
  unsigned char *a = alpha;
  unsigned char *p = (unsigned char *) plotmem;
  for(size_t n = 0; n < npixels; p = p + 4, d = d + 3, a = a + 1, n++) {
    memcpy(d,p,3);
    memcpy(a,p+3,1);
  }
  graph = wxImage(width,height,data,alpha);
  free(plotmem);
}

void PlPlotCairo::InitGraph()
{
  wxColour ct = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWTEXT);

  plinit();
  pladv(0);
  plvsta();
  plsfont(PL_FCI_SERIF,0,0);
  plscol0a(1,ct.Red(),ct.Green(),ct.Blue(),0.5);
  plcol0(1);
  plwind(0.0,1.0,0.0,1.0);
  plbox( "bctn", 0.0, 0, "bctn", 0.0, 0 );
  plptex(0.5,0.5,0.0,0.0,0.5,"Initialising...");
  plend();
}

// -- MuniPlotHistoBase  --------------------------------------------

MuniPlotHistoBase::MuniPlotHistoBase(wxWindow *w, const wxSize& sz):
  PlPlotCairo(w,sz), xmin(0.0),xmax(1.0),ymin(0.0),ymax(1.0) {}

void MuniPlotHistoBase::SetHisto(const FitsHisto& h)
{
  wxASSERT(h.IsOk());

  hist = h;

  float bw = hist.BinWidth();
  float bw2 = bw / 2.0;
  float fmin = hist.MinFreq();
  float fmax = hist.MaxFreq();

  xmin = hist.MinVal() - bw2;
  xmax = hist.MaxVal() + bw2;
  ymin = round(log10(fmin)-0.5);
  ymax = round(log10(fmax)+0.5);
  Update();
}

void MuniPlotHistoBase::StrokeHistoProfile(int n1, int n2,
					   PLINT& n, PLFLT *x, PLFLT *y)
{
  wxASSERT(hist.IsOk());

  x[0] = hist.Bin(n1);
  y[0] = ymin;
  int i = 0;
  for(int j = n1; j < n2; j++) {
    float f = hist.Freq(j) > 0 ? log10(hist.Freq(j)) : ymin;
    x[2*i+1] = hist.Bin(j);
    y[2*i+1] = f;
    x[2*i+2] = hist.Bin(j+1);
    y[2*i+2] = f;
    i++;
  }
  x[2*i+1] = hist.Bin(n2);
  y[2*i+1] = ymin;
  n = 2*i + 2;
}


// -- MuniPlotHisto  --------------------------------------------

MuniPlotHisto::MuniPlotHisto(wxWindow *w, const wxSize& sz,double b, double s):
  MuniPlotHistoBase(w,sz), black(b),sense(s) {}

void MuniPlotHisto::SetScale(double b, double s)
{
  black = b;
  sense = s;
  Update();
}

void MuniPlotHisto::MakeGraph()
{
  wxASSERT(hist.IsOk());

  wxColour cb = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWFRAME);
  wxColour ct = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWTEXT);

  float bw2 = hist.BinWidth() / 2.0;

  plinit();
  pladv(0);
  plvsta();
  plscol0a(0,cb.Red(),cb.Green(),cb.Blue(),1.0);
  plscol0a(1,ct.Red(),ct.Green(),ct.Blue(),1.0);
  plcol0(1);
  plsfont(PL_FCI_SERIF,0,0);
  plwind(xmin,xmax,ymin,ymax);
  plbox( "bctn", 0.0, 0, "bncvtl", 0.0, 0 );
  pllab( "Intensity", "Frequency", "" );

  const PLINT nline = 2*hist.NBins() + 2;
  PLFLT x[nline],y[nline];
  PLINT n;
  StrokeHistoProfile(0,hist.NBins()-1,n,x,y);
  plline(n,x,y);

  // gradient palette
  PLINT npal = hist.NBins();
  PLINT rpal[npal],gpal[npal],bpal[npal];
  for(int i = 0; i < hist.NBins(); i++) {
    float f = (hist.Bin(i) + bw2 - black) / sense;
    rpal[i] = gpal[i] = bpal[i] = min(max(int(255*f + 0.5),0),255);
  }
  plscmap1(rpal,gpal,bpal,npal);
  plgradient(n,x,y,0.0);

  plend();
}


//  PlotNite   -------------
MuniPlotNite::MuniPlotNite(wxWindow *w, const wxSize& sz,float t, float u, bool e):
  MuniPlotHistoBase(w,sz), level(t),width(u), enabled(e) { }

void MuniPlotNite::SetMeso(float t, float w)
{
  level = t;
  width = w;
  Update();
}

void MuniPlotNite::SetEnabled(bool e)
{
  enabled = e;
  Update();
}


void MuniPlotNite::MakeGraph()
{
  wxASSERT(hist.IsOk());

  wxColour cb = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWFRAME);
  wxColour ct = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWTEXT);

  // pallete
  PLINT npal = hist.NBins();
  PLINT rpal[npal], gpal[npal], bpal[npal];

  plinit();

  pladv(0);
  plvsta();
  plscol0a(0,cb.Red(),cb.Green(),cb.Blue(),1.0);
  plscol0a(1,ct.Red(),ct.Green(),ct.Blue(),1.0);
  plcol0( 1 );
  plsfont(PL_FCI_SERIF,0,0);

  plwind(xmin,xmax,ymin,ymax);
  plbox( "bctn", 0.0, 0, "bncvtl", 0.0, 0 );
  pllab( "Luminance Y component in CIE 1931 XYZ", "Frequency", "" );

  // the histogram silouhette
  PLINT nline = 2*hist.NBins()+2;
  PLFLT x[nline],y[nline];
  PLINT n;
  StrokeHistoProfile(0,hist.NBins()-1,n,x,y);
  plline(n,x,y);

  if( enabled ) {

    int n1, n2;

    // scotopic, nearly black background
    n1 = 0;
    n2 = hist.NBins() - 1;
    for(int i = 0; i < hist.NBins(); i++) {
      if( hist.Bin(i) < level - width / 2 )
	n2 = i;
    }
    if( n2 > n1 ) {
      StrokeHistoProfile(n1,n2,n,x,y);
      plscol0(0,64,64,64);
      plcol0(0);
      plfill(n,x,y);
    }

    // mescopic, having low saturated palette
    n1 = 0;
    n2 = hist.NBins()-1;
    for(int i = 0; i < hist.NBins(); i++) {
      if( hist.Bin(i) < level - width / 2 )
	n1 = i;
      if( hist.Bin(i) > level - width / 2 )
	n2 = i - 1;
    }
    n1 = min(n1,hist.NBins()-1);
    n2 = max(n2,0);
    if( n2 > n1 ) {
      StrokeHistoProfile(n1,n2,n,x,y);
      MakePalette(0.5,npal,rpal,gpal,bpal);
      plscmap1(rpal,gpal,bpal,npal);
      plgradient(n,x,y,0.0);
    }

    // photopic, with full colour palette
    n1 = 0;
    n2 = hist.NBins() - 1;
    for(int i = 0; i < hist.NBins(); i++) {
      if( hist.Bin(i) < level + width / 2 )
	n1 = i + 1;
    }
    n2 = min(n2,hist.NBins() - 1);
    if( n2 > n1 ) {
      StrokeHistoProfile(n1,n2,n,x,y);
      MakePalette(1.0,npal,rpal,gpal,bpal);
      plscmap1(rpal,gpal,bpal,npal);
      plgradient(n,x,y,0.0);
    }

  }
  else {

    // full colour gradient in range yellow to blue in HSL is demonstrated
    MakePalette(1.0,npal,rpal,gpal,bpal);
    plscmap1(rpal,gpal,bpal,npal);
    plgradient(n,x,y,0.0);
  }

  plend();
}


  void MuniPlotNite::MakePalette(PLFLT sat, PLINT npal, PLINT *rpal, PLINT *gpal, PLINT *bpal)
{
  float d = (240.0 - 60.0) / float(npal); // yellow to blue, counter-clockwise
  for(int n = 0; n < npal; n++) {
    PLFLT H = 240.0 - n*d;
    PLFLT r,g,b;
    plhlsrgb(H,0.5,sat,&r,&g,&b);
    rpal[n] = int(255.0*r);
    gpal[n] = int(255.0*g);
    bpal[n] = int(255.0*b);
  }
}



// ----------------------------------------------------------------

MuniPlotItt::MuniPlotItt(wxWindow *w, const wxSize& s, const FitsItt& it):
  PlPlotCairo(w,s),itt(it)
{
  Update();
}

void MuniPlotItt::SetItt(int i)
{
  itt.SetItt(i);
  Update();
}

void MuniPlotItt::MakeGraph()
{
  wxLogDebug("MuniPlotItt::MakeGraph()");

  wxColour cb = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWFRAME);
  wxColour ct = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWTEXT);
  wxSize s = GetSize();
  PLINT width = s.GetWidth();

  plinit();

  pladv(0);
  plvsta();
  plscol0a(0,cb.Red(),cb.Green(),cb.Blue(),1.0);
  plscol0a(1,ct.Red(),ct.Green(),ct.Blue(),1.0);
  plcol0( 1 );
  plsfont(PL_FCI_SERIF,0,0);

  plwind(0.0,3.0,0.0,1.0);
  plbox( "bctn", 0.0, 0, "bncvt", 0.0, 0 );
  pllab( "Input", "Output", "" );

  PLINT nitt = width/3;
  float dx = 3.0 / nitt; // step
  PLFLT x[nitt],y[nitt];
  for(int l = ITT_FIRST+1; l < ITT_LAST; l++) {
    FitsItt it(l);
    for(int i = 0; i < nitt; i++) {
      x[i] = dx*i;
      y[i] = it.Scale(x[i]);
    }
    int wpen = it.GetItt() == itt.GetItt() ? 3 : 1;
    plwidth(wpen);
    pllsty(l - ITT_FIRST);
    plline(nitt,x,y);
  }
  plend();
}


// -------------------------------------------------------

MuniPlotPalette::MuniPlotPalette(wxWindow *w, const wxSize& sz, int palid, bool inv,
				 int i, float b, float s):
  PlPlotCairo(w,sz), pal(palid,inv), itt(i), black(b), sense(s)
{
  Update();
}

void MuniPlotPalette::SetPalette(int p)
{
  pal.SetPalette(p);
  Update();
}

void MuniPlotPalette::SetItt(int i)
{
  itt.SetItt(i);
  Update();
}

void MuniPlotPalette::SetInverse(bool inv)
{
  pal.SetInverse(inv);
  Update();
}

void MuniPlotPalette::SetScale(float b, float s)
{
  black = b;
  sense = s;
  Update();
}

void MuniPlotPalette::MakeGraph()
{
  float white = black + sense*itt.InvScale(pal.GetColours()/256.0);
  bool invert = pal.GetInverse();

  wxColour cb = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWFRAME);
  wxColour ct = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWTEXT);

  plinit();
  pladv(0);
  plvsta();
  plscol0a(0,cb.Red(),cb.Green(),cb.Blue(),1.0);
  plscol0a(1,ct.Red(),ct.Green(),ct.Blue(),1.0);
  plcol0( 1 );
  plsfont(PL_FCI_SERIF,0,0);

  PLINT npal = pal.GetColours();
  PLINT rpal[npal], gpal[npal], bpal[npal];
  for(int i = 0; i < npal; i++) {
    rpal[i] = pal.R(i);
    gpal[i] = pal.G(i);
    bpal[i] = pal.B(i);
  }
  plscmap1(rpal,gpal,bpal,npal);

  plsfont(PL_FCI_SERIF,0,0);
  plvpor(0.1,0.9,0.2,0.9);

  float xmin = invert ? white : black;
  float xmax = invert ? black : white;

  plwind(xmin,xmax,0.0,1.0);

  PLFLT xbox[4] = { xmin, xmax, xmax, xmin };
  PLFLT ybox[4] = { 0.0, 0.0, 1.0, 1.0 };
  PLFLT angle = invert ? 180.0 : 0.0;
  plgradient(4,xbox,ybox,angle);

  plbox( "bctn", 0.0, 0, "bc", 0.0, 0 );
  pllab( "Intensity", "", "" );

  plend();
}

//  PlotSaturation
MuniPlotSaturation::MuniPlotSaturation(wxWindow *w, const wxSize& sz, float s,
				       float x, float y):
  PlPlotCairo(w,sz),saturation(s)
{
  white_point[0] = x;
  white_point[1] = y;
  Update();
}

void MuniPlotSaturation::SetSaturation(float s)
{
  saturation = s;
  Update();
}

void MuniPlotSaturation::SetWhitePoint(float x, float y)
{
  white_point[0] = x;
  white_point[1] = y;
  Update();
}

void MuniPlotSaturation::MakeTable(PLFLT hue, PLFLT l,
				   PLINT n, PLINT *red, PLINT *green, PLINT *blue)
{
  PLFLT r,g,b;
  for(int i = 0; i < n; i++) {
    PLFLT s = 0.1*i;
    plhlsrgb(hue,l,s,&r,&g,&b);
    red[i] = max(min(int(255.0*r),255),0);
    green[i] = max(min(int(255.0*g),255),0);
    blue[i] = max(min(int(255.0*b),255),0);
  }
}

void MuniPlotSaturation::MakeGraph()
{
  wxColour cb = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWFRAME);
  wxColour ct = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWTEXT);

  const float w = 0.1;
  const float h = 0.1;
  const float x0 = 0;
  const float y0 = 0;
  const PLFLT lum = 0.35;

  plinit();
  pladv(0);
  plvsta();
  plvpor(0.1,0.9,0.2,0.9);
  plwind( -w/2.0, 2.0+w/2.0, 0.0, 0.3);

  PLINT nsat = 21;
  PLINT red[nsat], green[nsat], blue[nsat];
  PLFLT opaq[nsat];
  for(int i = 0; i < nsat; i++) opaq[i] = 1.0;

  for(int k = 0; k < 3; k++) {
    PLFLT hue = 120*k;
    MakeTable(hue,lum,nsat,red,green,blue);
    plscmap0a(red, green, blue, opaq, nsat );
    for(int i = 0; i < nsat; i++) {
      plcol0(i);
      plfbox(x0+i*w - w/2.0,y0+k*h,w,h);
    }
  }

  plscol0a(0,cb.Red(),cb.Green(),cb.Blue(),1.0);
  plscol0a(1,ct.Red(),ct.Green(),ct.Blue(),1.0);
  plcol0( 1 );
  plsfont(PL_FCI_SERIF,0,0);
  plbox( "bctn", 0.0, 0, "bc", 0.0, 0 );
  pllab( "Saturation", "", "" );

  plscol0a(1,ct.Red(),ct.Green(),ct.Blue(),0.7);
  plcol0( 1 );
  plftri(x0 + 10*saturation*w,0*3*h,w/2.0,h/4.0);

  plend();
}

void MuniPlotSaturation::plfbox(PLFLT x0, PLFLT y0, PLFLT w, PLFLT h)
{
  PLFLT x[4], y[4];

  x[0] = x0;
  y[0] = y0;
  x[1] = x0 + w;
  y[1] = y0;
  x[2] = x0 + w;
  y[2] = y0 + h;
  x[3] = x0;
  y[3] = y0 + h;
  plfill( 4, x, y );
}

void MuniPlotSaturation::plftri(PLFLT x0, PLFLT y0, PLFLT w, PLFLT h)
{
  PLFLT x[3], y[3];

  // bottom
  x[0] = x0;
  y[0] = y0;
  x[1] = x0 - w / 2.0;
  y[1] = y0 + h / 2.0;
  x[2] = x0 + w / 2.0;
  y[2] = y0 + h / 2.0;
  plfill( 3, x, y );

  // top
  float height = 0.3; // must match plwind() last parameter
  x[0] = x0;
  y[0] = height;
  x[1] = x0 - w / 2.0;
  y[1] = height - h/2.0;
  x[2] = x0 + w / 2.0;
  y[2] = height - h/2.0;
  plfill( 3, x, y );

  // connecting line
  x[0] = x0;
  y[0] = 0;
  x[1] = x0;
  y[1] = height;
  plline(2,x,y);
}


// --- MuniPlotTable -------------------------------------

/*
MuniPlotTable::MuniPlotTable(const vector<wxRealPoint>& p,const wxColour& c):
  points(p),colour(c)
{
}
*/

// --- MuniPlot ------------------------------------------


// MuniPlot::MuniPlot(wxWindow *w):
//   wxPLplotwindow(w,wxID_ANY,wxDefaultPosition,wxDefaultSize,wxWANTS_CHARS,
// 		 PLPLOT_OPTIONS),
//   xmin(numeric_limits<double>::max()),xmax(numeric_limits<double>::min()),
//   ymin(numeric_limits<double>::max()),ymax(numeric_limits<double>::min())
// {
//   Draw();
// }

// void MuniPlot::AddData(const MuniPlotTable& t)
// {
//   for(size_t i = 0; i < t.points.size(); i++) {
//     wxRealPoint p(t.points[i]);
//     if( p.x < xmin ) xmin = p.x;
//     if( p.x > xmax ) xmax = p.x;
//     if( p.y < ymin ) ymin = p.y;
//     if( p.y > ymax ) ymax = p.y;
//   }

//   tables.push_back(t);
//   Draw();
// }



// void MuniPlot::Draw()
// {
//   wxColour cb = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWFRAME);

//   wxPLplotstream* pls = GetStream();
//   pls->adv( 0 );
//   pls->schr(0.0,3.0);
//   //  pls->scol0a(0,255,255,255,1.0);
//   pls->scol0a(0,cb.Red(),cb.Green(),cb.Blue(),1.0);
//   pls->scol0a(1,0,0,0,1.0);

//   if( tables.empty() ) {
//     pls->vpor(0.01, 0.99, 0.2, 0.99);
//     pls->wind(-0.5, 10.5, 0.0, 1.0);
//     pls->box("bcnt",0.0,0.0,"bct",0.0,0.0);
//   }
//   else {

//     pls->scol0a(2,255,0,0,0.2);
//     pls->scol0a(3,0,255,0,0.2);
//     pls->scol0a(4,0,0,255,0.2);
//     pls->scol0a(5,170,170,170,0.2);

//     pls->scol0a(12,255,0,0,1.0);
//     pls->scol0a(13,0,255,0,1.0);
//     pls->scol0a(14,0,0,255,1.0);
//     pls->scol0a(15,170,170,170,1.0);

//     pls->vpor(0.01, 0.99, 0.2, 0.99);
//     pls->wind(xmin, xmax, ymin, ymax);
//     pls->box("bcnt",0.0,0.0,"bct",0.0,0.0);

//   }
//   //    const size_t nh = 50;
// //   PLFLT         hx[nh],hr[nh],hg[nh],hb[nh];
// //     for ( size_t i = 0; i < nh; i++ ) {
// //       double t = (i - 25.0)/10.0;
// //       hx[i] = t;
// //       hr[i] = exp(-(t-1)*(t-1)/2.0);
// //       hg[i] = exp(-(t-0)*(t-0)/2.0);
// //       hb[i] = exp(-(t+1)*(t+1)/2.0);
// //       //      wxLogDebug(_("%f %f"),hx[i],hy[i]);
// //     }


//   RenewPlot();
// }




// MuniPlotUV::MuniPlotUV(wxWindow *w):
//   wxPLplotwindow(w,wxID_ANY,wxDefaultPosition,wxDefaultSize,wxWANTS_CHARS,PLPLOT_OPTIONS),
//   nuv(0),u(0),v(0)
// {
//   Draw();
// }

// MuniPlotUV::~MuniPlotUV()
// {
//   delete[] u;
//   delete[] v;
// }

// wxSize MuniPlotUV::DoGetBestSize() const
// {
//   return wxSize(200,200);
// }


// void MuniPlotUV::Clear()
// {
//   nuv = 0;
//   delete[] u;
//   delete[] v;
//   Draw();
// }

// void MuniPlotUV::DrawTri(const vector<double>& uu, const vector<double>& vv)
// {
//   wxASSERT(uu.size() == vv.size());

//   delete[] u;
//   delete[] v;

//   nuv = uu.size();
//   u = new PLFLT[nuv];
//   v = new PLFLT[nuv];

//   for(int i = 0; i < nuv; i++) {
//     u[i] = uu[i];
//     v[i] = vv[i];
//   }

//   Draw();
// }


// void MuniPlotUV::Draw()
// {
//   wxColour cb = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWFRAME);

//   wxPLplotstream* pls = GetStream();
//   pls->adv(0);

//   pls->scol0a(0,cb.Red(),cb.Green(),cb.Blue(),1.0);
//   pls->scol0a(1,0,0,0,1.0);

//   //  pls->schr(0.0,3.0);

//   //  pls->vpor(0.0, 1.0, 0.0, 1.0);
//   //  pls->wind(0.0, 1.0, 0.0, 1.0);
//   //  pls->box("bcnt",0.1,0.0,"bcnt",0.1,0.0);
//   pls->env(0.0,1.0,0.0,1.0,2.0,0.0);
//   //  pls->lab("u","v","Triangles in uv space");

//   /*
//   const int npoints = 3;
//   PLFLT x[npoints],y[npoints];
//   x[0] = 0.9; y[0] = 0.9;
//   x[1] = 0.1; y[1] = 0.7;
//   x[2] = 0.7; y[2] = 0.1;
//   */
//   if( nuv > 0 )
//     pls->line(nuv,u,v);

//   RenewPlot();
// }
