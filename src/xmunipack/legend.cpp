/*

  xmunipack - draws legend on the display

  Copyright © 2022 F.Hroch (hroch@physics.muni.cz)

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

#include "legend.h"
#include "mconfig.h"
#include <wx/wx.h>
#include <wx/graphics.h>
#include <vector>

#define RAD 57.29577951308232
#define MAJOR_TICS 7
#define MINOR_TICS 3

using namespace std;


MuniDisplayLegend::MuniDisplayLegend(wxGraphicsContext *wgc, const MuniConfig *c,
				     const wxRect& r, double x, double y, double z):
  gc(wgc),config(c),exposed_area(r),xoff(x),yoff(y),zoom(z),fgcolour(65,105,255)
{
  wxASSERT(gc);
  wxLogDebug("MuniDisplayLegend() %lf %lf",xoff,yoff);
  gc->GetSize(&width,&height);


  if( config->display_cootype == COO_EQDEG ) {
    const int ntics = 25;
    const double tics[25] = { 90.0, 60.0, 30.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1,
    0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 5e-5, 2e-5,
    1e-5, 5e-6, 2e-6, 1e-6};
    dtics.assign(tics,tics+ntics);
    atics.assign(tics,tics+ntics);
  }
  else if( config->display_cootype == COO_EQSIX ) {
    const int ntics = 15;
    const double tics[15] = { 90.0, 30.0, 10.0, 5.0, 1.0,
      30.0 / 60.0, 10.0 / 60.0, 5.0 / 60.0, 1.0 / 60.0,
      30.0 / 3600.0, 10.0 / 3600.0, 5.0 / 3600.0, 1.0 / 3600.0,
      0.5 / 3600.0, 0.1 / 3600.0 };
    const int mtics = 13;
    const double rtics[mtics] = { 6.0, 3.0, 1.0,
      30.0 / 60.0, 10.0 / 60.0, 5.0 / 60.0, 1.0 / 60.0,
      30.0 / 3600.0, 10.0 / 3600.0, 5.0 / 3600.0, 1.0 / 3600.0,
      0.5 / 3600.0, 0.1 / 3600.0 };
    atics.assign(rtics,rtics+mtics);
    dtics.assign(tics,tics+ntics);
  }
}


double MuniDisplayLegend::Y(double y) const
{
  return height-1 - zoom*y - yoff;
}

double MuniDisplayLegend::X(double x) const
{
  return xoff + x*zoom;
}

double MuniDisplayLegend::W(double w) const
{
  return zoom*w;
}

double MuniDisplayLegend::H(double h) const
{
  return zoom*h;
}


bool MuniDisplayLegend::Draw(const FitsArray& array)
{
  wxASSERT(array.IsOk());

  wxColour gridcolour(fgcolour.Red(),fgcolour.Green(),fgcolour.Blue(),64);
  wxPen pen(gridcolour,1);
  gc->SetPen(pen);
  gc->SetBrush(*wxTRANSPARENT_BRUSH);
  gc->SetFont(*wxSMALL_FONT,fgcolour);
  gc->DrawRectangle(xoff,yoff,zoom*array.GetWidth()-1,zoom*array.GetHeight()-1);

  return DrawGrid(array) /*&& DrawRuler(array) && DrawArrow(array)*/;
}

bool MuniDisplayLegend::DrawGrid(const FitsArray& array)
{

  double xcen, ycen, rcen, dcen, scale, phi, reflex;
  if( ! array.GetWCS(xcen,ycen,rcen,dcen,scale,phi,reflex) && fabs(phi) > 60 ) {
    wxLogDebug("MuniDisplayLegend::DrawGrid(): "
	       "*** WARNING: HIGH declination > 60 deg ****");
    // The singular points of are poles of my fear.
    return false;
  }

  double xmin = exposed_area.x;
  double ymin = exposed_area.y;
  double xmax = exposed_area.x + exposed_area.width;
  double ymax = exposed_area.y + exposed_area.height;
  wxLogDebug("DrawGrid() %lf %f %f %f %lf",xmin,ymin,xmax,ymax,zoom);

  // declination
  FitsProjection p("GNOMONIC",rcen,dcen,xcen,ycen,scale,phi,reflex);
  wxLogDebug("FitsProjection: %f %f %f %f %f %f",rcen,dcen,xcen,ycen,scale,phi);

  double rmin, rmax, dmin, dmax;
  p.xy2ad(xmin,ymin,rmin,dmin); // left-bottom
  p.xy2ad(xmax,ymax,rmax,dmax); // right-top
  //  wxLogDebug(": %lf %lf %lf %lf",rmin,rmax,dmin,dmax);

  // tics
  double dfov = fabs(dmax - dmin);
  vector<double> dmajor, dminor;
  double dtic;
  TicsArray(dmin,dmax,dtics,dmajor,dminor,dtic);

  double rfov = fabs(rmax - rmin);
  vector<double> rmajor, rminor;
  double rtic;
  double rfac = config->display_cootype == COO_EQSIX ? 15 : 1;
  TicsArray(rmax/rfac,rmin/rfac,atics,rmajor,rminor,rtic);
  if( config->display_cootype == COO_EQSIX ) {
    rtic = 15*rtic;
    for(size_t i = 0; i < rmajor.size(); i++)
      rmajor[i] = 15*rmajor[i];
    for(size_t i = 0; i < rminor.size(); i++)
      rminor[i] = 15*rminor[i];
  }

  double iwidth = array.GetWidth();
  double iheight = array.GetHeight();

  gc->Clip(xoff,yoff,zoom*iwidth,zoom*iheight);
  DrawMeridians(p,rcen,dcen,rmajor,dfov);
  DrawParallels(p,rcen,dcen,dmajor,rfov);
  gc->ResetClip();

  wxPen pen2(fgcolour,1);
  gc->SetPen(pen2);

  DrawTicsBottom(p,iwidth,rcen,dcen,rmajor,dfov,MAJOR_TICS);
  DrawTicsBottom(p,iwidth,rcen,dcen,rminor,dfov,MINOR_TICS);
  DrawTicsTop(p,iheight,rcen,dcen,rmajor,dfov,MAJOR_TICS);
  DrawTicsTop(p,iheight,rcen,dcen,rminor,dfov,MINOR_TICS);
  DrawTicsLeft(p,iheight,rcen,dcen,dmajor,rfov,MAJOR_TICS);
  DrawTicsLeft(p,iheight,rcen,dcen,dminor,rfov,MINOR_TICS);
  DrawTicsRight(p,iwidth,rcen,dcen,dmajor,rfov,MAJOR_TICS);
  DrawTicsRight(p,iwidth,rcen,dcen,dminor,rfov,MINOR_TICS);

  DrawLabelsBottom(p,rcen,dcen,rmajor,dfov,MAJOR_TICS,rtic);
  DrawLabelsLeft(p,rcen,dcen,dmajor,rfov,MAJOR_TICS,dtic);

  gc->ResetClip();

  return true;
}


void MuniDisplayLegend::DrawRectangle(wxDouble x, wxDouble y, wxDouble w, wxDouble h)
{
  wxLogDebug("DrawRectangle(0): %lf %lf %lf %lf",x,y,w,h);
  gc->DrawRectangle(X(x),Y(y)-H(h),W(w),H(h));
  wxLogDebug("DrawRectangle(1): %lf %lf %lf %lf",X(x),Y(y),W(w),H(h));
}

void MuniDisplayLegend::DrawLine(double x1, double y1, double x2, double y2)
{
  wxGraphicsPath path = gc->CreatePath();
  path.MoveToPoint(X(x1),Y(y1));
  path.AddLineToPoint(X(x2),Y(y2));
  path.CloseSubpath();
  gc->StrokePath(path);
}

void MuniDisplayLegend::DrawText(const wxString& a, double x, double y, bool horizon)
{
  wxDouble w,h,ld,ll;
  gc->GetTextExtent(a,&w,&h,&ld,&ll);
  if( horizon )
    gc->DrawText(a,X(x)-w/2,Y(y)+h/2,0.0);
  else // vertical
    gc->DrawText(a,X(x)-h,Y(y)+w/2,90.0/RAD);
}

void MuniDisplayLegend::DrawAxisLabel(const wxString& a, double x, double y,
				      double phi, double xinc, double yinc)
{
  wxDouble w,h,td,tl;
  gc->GetTextExtent(a,&w,&h,&td,&tl);
  double dx = copysign(fabs(xinc) + 0.5,xinc);
  double dy = copysign(fabs(yinc) + 0.5,yinc);
  gc->DrawText(a,X(x)-dx*w,Y(y)-dy*h,phi/RAD);
  //  wxLogDebug(" %lf %lf",double(X(x)-w/2),double(Y(y)-dy*h));
}

void MuniDisplayLegend::DrawRa(const wxString& hour, const wxString& min,
			       const wxString& sec, double x, double y)
{
  wxDouble w1,h1,w2,h2,ld,ll,w3,w4,w5,w6;
  gc->GetTextExtent(hour,&w1,&h1,&ld,&ll);
  gc->GetTextExtent(min,&w2,&h2,&ld,&ll);

  wxFont font(*wxSMALL_FONT);
  font.MakeSmaller();
  gc->SetFont(font,fgcolour);
  gc->GetTextExtent("h",&w3,&h2,&ld,&ll);
  gc->GetTextExtent("m",&w4,&h2,&ld,&ll);
  gc->GetTextExtent("s",&w6,&h2,&ld,&ll);

  double h = h1;
  double w = w1 + w2 + w3 + w4;
  gc->SetFont(*wxSMALL_FONT,fgcolour);
  gc->DrawText(hour,X(x)-w/2,Y(y)+h/2);
  gc->SetFont(font,fgcolour);
  gc->DrawText("h",X(x)-w/2+w1,Y(y)+h/2);
  gc->SetFont(*wxSMALL_FONT,fgcolour);
  gc->DrawText(min,X(x)-w/2+w1+w3,Y(y)+h/2);
  gc->SetFont(font,fgcolour);
  gc->DrawText("m",X(x)-w/2+w1+w2+w3,Y(y)+h/2);
  gc->SetFont(*wxSMALL_FONT,fgcolour);

  if( sec != "" ) {
    gc->GetTextExtent(sec,&w5,&h2,&ld,&ll);
    gc->DrawText(sec,X(x)-w/2+w1+w2+w3+w4,Y(y)+h/2);
    gc->SetFont(font,fgcolour);
    gc->DrawText("s",X(x)-w/2+w1+w2+w3+w4+w5,Y(y)+h/2);
    gc->SetFont(*wxSMALL_FONT,fgcolour);
  }
}


void MuniDisplayLegend::TicsArray(double dmin, double dmax,
				  const std::vector<double>& ticpool,
				  std::vector<double>& major,
				  std::vector<double>& minor, double& tic) const
{
  minor.clear();

  for(size_t i = 0; i < ticpool.size(); i++) {
    double t = ticpool[i];
    double d1 = round(dmin/t)*t;
    double d2 = round(dmax/t)*t;
    int c = round((d2 - d1) / t) + 1;

    major.clear();
    for(int l = 0; l < c; l++) {
      double d = d1 + l*t;
      if( dmin <= d && d <= dmax )
	major.push_back(d);
    }

    if( major.size() > 1) {

      double t = ticpool[min(i+1,ticpool.size()-1)];
      double d1 = round(dmin/t)*t;
      double d2 = round(dmax/t)*t;
      int c = round((d2 - d1) / t) + 1;
      for(int l = 0; l < c; l++) {
	double d = d1 + l*t;
	if( dmin <= d && d <= dmax )
	  minor.push_back(d);
      }

      tic = t;
      return;
    }
  }

  major.clear();
}


void MuniDisplayLegend::DrawMeridians(const FitsProjection& p, double rcen, double dcen, const std::vector<double>& rmajor, double dfov)
{
  const int m = 10;
  const double d = dfov / m;

  // draw Declination circles
  for(size_t i = 0; i < rmajor.size(); i++) {
    double ra = rmajor[i];
    for( int j = -2*m; j <= 2*m; j++) {
      double dec = dcen + j*d;
      double x1,y1,x2,y2;
      p.ad2xy(ra,dec,x1,y1);
      p.ad2xy(ra,dec+d,x2,y2);
      DrawLine(x1,y1,x2,y2);
    }
  }
}

void MuniDisplayLegend::DrawParallels(const FitsProjection& p, double rcen, double dcen, const std::vector<double>& major, double rfov)
{
  const int m = 10;
  const double d = rfov / m;

  for(size_t j = 0; j < major.size(); j++) {
    double dec = major[j];

    for( int i = -m; i < m; i++) {
      double ra = rcen + i*d;
      double x1,y1,x2,y2;
      p.ad2xy(ra,dec,x1,y1);
      p.ad2xy(ra+d,dec,x2,y2);
      DrawLine(x1,y1,x2,y2);
    }
  }
}



void MuniDisplayLegend::DrawTicsBottom(const FitsProjection& p, double iwidth, double rcen, double dcen, const std::vector<double>& tics, double dfov, double h)
{
  const int m = 10;
  const double d = dfov / m;

  for(size_t i = 0; i < tics.size(); i++) {
    double ra = tics[i];
    for( int j = -2*m; j < 2*m; j++) {
      double dec = dcen + j*d;
      double x1,y1,x2,y2;
      p.ad2xy(ra,dec,x1,y1);
      p.ad2xy(ra,dec+d,x2,y2);

      if( y1*y2 < 0 ) {
	double x = (x1 + x2) / 2;
	double q = atan2(y2-y1,x2-x1);
	y2 = -h*sin(q);
	x2 = x - h*cos(q);
	if( 0 <= x && x < iwidth )
	  DrawLine(x,0.0,x2,y2);
	break;
      }
    }
  }

}

void MuniDisplayLegend::DrawTicsTop(const FitsProjection& p, double h,
				    double rcen, double dcen,
				    const std::vector<double>& tics,
				    double dfov, double t)
{
  const int m = 10;
  const double d = dfov / m;

  for(size_t i = 0; i < tics.size(); i++) {
    double ra = tics[i];
    for( int j = -2*m; j < 2*m; j++) {
      double dec = dcen + j*d;
      double x1,y1,x2,y2;
      p.ad2xy(ra,dec,x1,y1);
      p.ad2xy(ra,dec+d,x2,y2);

      if( (y1-h)*(y2-h) < 0 ) {
	double x = (x1 + x2) / 2;
	double q = atan2(y2-y1,x2-x1);
	y2 = h + t*sin(q);
	x2 = x + t*cos(q);
	if( 0 <= x && x < h )
	  DrawLine(x,h,x2,y2);
	break;
      }
    }
  }

}

void MuniDisplayLegend::DrawTicsLeft(const FitsProjection& p, double iheight, double rcen, double dcen, const std::vector<double>& tics, double rfov, double h)
{
  const int m = 10;
  const double d = rfov / m;

  for(size_t j = 0; j < tics.size(); j++) {
    double dec = tics[j];
    for( int i = -m; i < m; i++) {
      double ra = rcen + i*d;
      double x1,y1,x2,y2;
      p.ad2xy(ra,dec,x1,y1);
      p.ad2xy(ra+d,dec,x2,y2);

      if( x1*x2 < 0 ) {
	double y = (y1 + y2) / 2;
	double q = atan2(y2-y1,x2-x1);
	y2 = y + h*sin(q);
	x2 = h*cos(q);
	if( 0 <= y && y < iheight )
	  DrawLine(0.0,y,x2,y2);
	break;
      }
    }
  }
}



void MuniDisplayLegend::DrawTicsRight(const FitsProjection& p, double w,
				      double rcen, double dcen,
				      const std::vector<double>& tics,
				      double rfov, double h)
{
  const int m = 10;
  const double d = rfov / m;

  for(size_t j = 0; j < tics.size(); j++) {
    double dec = tics[j];
    for( int i = -2*m; i < 2*m; i++) {
      double ra = rcen + i*d;
      double x1,y1,x2,y2;
      p.ad2xy(ra,dec,x1,y1);
      p.ad2xy(ra+d,dec,x2,y2);

      if( (x1-w)*(x2-w) < 0 ) {
	double y = (y1 + y2) / 2;
	double q = atan2(y2-y1,x2-x1);
	y2 = y - h*sin(q);
	x2 = w - h*cos(q);
	if( 0 <= y && y < w )
	  DrawLine(w,y,x2,y2);
	break;
      }
    }
  }
}

void MuniDisplayLegend::DrawLabelsBottom(const FitsProjection& p,
					 double rcen, double dcen,
					 const std::vector<double>& rmajor,
					 double dfov, double lentic, double rtic)
{
  const int m = 10;
  const double d = dfov / m;

  for(size_t i = 0; i < rmajor.size(); i++) {
    double ra = rmajor[i];
    for( int j = -2*m; j < 2*m; j++) {
      double dec = dcen + j*d;
      double x1,y1,x2,y2;
      p.ad2xy(ra,dec,x1,y1);
      p.ad2xy(ra,dec+d,x2,y2);

      if( y1*y2 < 0 ) {
	double x = (x1 + x2) / 2;

	wxString a = RA(ra,rtic);

	if( config->display_cootype == COO_EQSIX ) {
	  wxString h = a.BeforeFirst('h');
	  wxString min = a.Mid(3).BeforeFirst('m');
	  wxString sec = a.Mid(6).BeforeFirst('s');
	  DrawRa(h,min,sec,x,0);
	}
	else
	  DrawText(a,x,0);
	break;
      }
    }
  }

}

void MuniDisplayLegend::DrawLabelsLeft(const FitsProjection& p,
				       double rcen, double dcen,
				       const std::vector<double>& dmajor,
				       double rfov, double lentic, double dtic)
{
  const int m = 10;
  const double d = rfov / m;

  for(size_t j = 0; j < dmajor.size(); j++ ) {
    double dec = dmajor[j];
    for( int i = -2*m; i < 2*m; i++) {
      double ra = rcen + i*d;
      double x1,y1,x2,y2;
      p.ad2xy(ra,dec,x1,y1);
      p.ad2xy(ra+d,dec,x2,y2);

      if( x1*x2 < 0 ) {
	double y = (y1 + y2) / 2;
	wxString a = Dec(dec,dtic);
	DrawText(a,-lentic,y,false);
	break;
      }
    }
  }

}

wxString MuniDisplayLegend::RA(double r, double tic) const
{
  int digits = wxMax(int(-log10(tic) + 0.5),0);

  wxString a;
  if( config->display_cootype == COO_EQDEG ) {
    a.Printf(L"%.*f°",digits,r);
    return a;
  }
  else if( config->display_cootype == COO_EQSIX ) {
    FitsCoo coo;
    int h,min;
    double s,alpha;
    alpha = r / 15.0;
    h = int(alpha);
    if( tic > 0.015 ) {
      min = int(60*(alpha-h)+0.5);
      a.Printf(L"%02dh%02dm",h,min);
    }
    else {
      min = int(60*(alpha-h));
      s = 3600.0*(alpha - (h + min/60.0));
      a.Printf(L"%02dh%02dm%02.0fs",h,min,round(s));
    }
  }
  return a;
}

wxString MuniDisplayLegend::Dec(double d, double tic) const
{
  int digits = wxMax(int(-log10(tic) + 0.5),0);
  wxString a;

  if( config->display_cootype == COO_EQDEG ) {
    a.Printf(L"%+.*f°",digits,d);
  }
  else if( config->display_cootype == COO_EQSIX ) {
    FitsCoo coo;
    int deg,min;
    char c = d >= 0.0 ? '+' : '-';
    double delta = fabs(d);
    deg = int(delta);
    min = int(60*(delta - deg)+0.5);
    //    coo.DecSix(d,c,deg,min);
      //    if( c == ' ' ) c = '+';
    a.Printf(L"%c%02d°%02d'",c,deg,min);
  }
  return a;
}

bool MuniDisplayLegend::DrawRuler(const FitsArray& array)
{
  if( width < 666 || height < 333 ) return false;

  double x0, y0, a0, d0, scale, phi, reflex;
  if( ! array.GetWCS(x0,y0,a0,d0,scale,phi,reflex) )
    return false;

  FitsProjection p("GNOMONIC",a0,d0,x0,y0,scale,phi,reflex);
  double fov = p.xy2dist(0,0,array.GetWidth(),array.GetHeight());
  //  wxLogDebug("FOV: %f",fov);

  //  double h = fov / 10;
  double rlen, ulen;
  if( fov > 1 )
    // degree
    ;
  else if( fov > 0.01 ) {
    // minutes
    rlen = 5.0 / 60.0;
    ulen = 1.0 / 60.0;
  }
  else if( fov > 1e-4 )
    // seconds
    ;

  double d = rlen * scale;
  double u = ulen * scale;

  wxColour colour(fgcolour.Red(),fgcolour.Green(),fgcolour.Blue(),128);
  wxPen pen(colour,1);
  gc->SetPen(pen);
  gc->SetBrush(*wxTRANSPARENT_BRUSH);

  wxDouble tw,th;
  GetTextDim("6'",tw,th);

  int thickness = wxMax(height / 100,10);
  int xoff = 3*thickness;
  int yoff = 2*thickness + 1.5*th;

  gc->SetAntialiasMode(wxANTIALIAS_NONE);
  DrawRectangle(xoff,yoff,d,thickness);
  gc->SetAntialiasMode(wxANTIALIAS_DEFAULT);

  wxBrush brush(fgcolour);
  gc->SetBrush(brush);
  DrawRectangle(xoff+u,yoff,u,thickness);

  double baseline = yoff + 0.5*th;
  DrawText("1'",xoff+u/2,baseline);
  DrawText("5'",xoff+d,baseline);

  return true;
}

void MuniDisplayLegend::GetTextDim(const wxString& text, wxDouble &w, wxDouble &h) const
{
  wxDouble descent,eleading;
  gc->GetTextExtent(text,&w,&h,&descent,&eleading);
}


bool MuniDisplayLegend::DrawArrow(const FitsArray& array)
{
  if( width < 666 || height < 333 ) return false;

  wxPen pen(fgcolour,3);
  gc->SetPen(pen);

  DrawLine(800,40,790,90);

  return true;
}
