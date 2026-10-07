/*

  xmunipack - icon


  Copyright © 2009-2012, 2018-2022 F.Hroch (hroch@physics.muni.cz)

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



  IMPORTANT

  ***  All functions MUST BE called in main thread ! ***

  The wxBitmap.ConvertToImage() calls under GTK+ function assigned
  directly to a window system (to a GC).

  .. else the execution will be locked by wxImage::ConvertToImage();

*/

#include "icon.h"
#include <cstring>
#include <wx/wx.h>
#include <wx/graphics.h>

using namespace std;


MuniIcon::MuniIcon(const MuniConfig *c): config(c),colour(*wxLIGHT_GREY) {}


void MuniIcon::GetIcons(const FitsFile& fits, wxImage& cover,
			std::vector<wxImage>& list) const
{
  // all HDU icons
  list = GetList(fits);

  // cover icon
  switch(fits.Type()) {
  case FITS_GREY:
  case FITS_COLOUR:
  case FITS_MULTI:
    wxASSERT(list.size() > 0);
    cover = GetCover(list[0],"",1);
    break;

  default:
    cover = Utf8Icon(2.0,config->default_symbol).ConvertToImage();
    break;
  }

}

std::vector<wxImage> MuniIcon::MergeList(const FitsFile& fits,
					 const std::vector<wxImage>& icons) const
{
  int s = config->icon_size;
  wxSize size(s,s);
  std::vector<wxImage> list;

  // create icons for all HDUs
  for(size_t k = 0; k < fits.HduCount(); k++) {

    if( icons[k].IsOk() )

      list.push_back(GetCover(icons[k]));

    else {

      switch(fits.Hdu(k).Type()) {
      case HDU_IMAGE:
	list.push_back(Utf8Icon(2.0,config->image_symbol).ConvertToImage());
	break;

      case HDU_HEAD:
	list.push_back(Utf8Icon(2.0,config->head_symbol).ConvertToImage());
	break;

      case HDU_TABLE:
	list.push_back(Utf8Icon(2.0,config->table_symbol).ConvertToImage());
	break;

      default:
	list.push_back(Utf8Icon(2.0,config->default_symbol).ConvertToImage());
	break;
      }
    }
  }

  return list;
}


std::vector<wxImage> MuniIcon::GetList(const FitsFile& fits) const
{
  wxSize size(config->icon_size,config->icon_size);
  std::vector<wxImage> list;

  // create icons for all HDUs
  for(size_t k = 0; k < fits.HduCount(); k++) {

    switch(fits.Hdu(k).Type()) {
    case HDU_IMAGE:
      list.push_back(Utf8Icon(2.0,config->image_symbol).ConvertToImage());
      break;

    case HDU_HEAD:
      list.push_back(Utf8Icon(2.0,config->head_symbol).ConvertToImage());
      break;

    case HDU_TABLE:
      list.push_back(Utf8Icon(2.0,config->table_symbol).ConvertToImage());
      break;

    default:
      list.push_back(Utf8Icon(2.0,config->default_symbol).ConvertToImage());
      break;
    }
  }

  return list;
}

std::vector<wxImage> MuniIcon::GetList(const FitsMeta& meta) const
{
  wxSize size(config->icon_size,config->icon_size);
  std::vector<wxImage> list;

  // create icons for all HDUs
  for(size_t k = 0; k < meta.HduCount(); k++) {

    switch(meta.Hdu(k).Type()) {
    case HDU_IMAGE:
      list.push_back(Utf8Icon(2.0,config->image_symbol).ConvertToImage());
      break;

    case HDU_HEAD:
      list.push_back(Utf8Icon(2.0,config->head_symbol).ConvertToImage());
      break;

    case HDU_TABLE:
      list.push_back(Utf8Icon(2.0,config->table_symbol).ConvertToImage());
      break;

    default:
      list.push_back(Utf8Icon(2.0,config->default_symbol).ConvertToImage());
      break;
    }
  }

  return list;
}


wxImage MuniIcon::ImageIcon(const wxImage& icon) const
{
  wxASSERT(icon.IsOk());
  return Padding(icon,config->icon_size);
}


wxImage MuniIcon::Padding(const wxImage& icon, int size) const
{
  int w,h;
  float r = float(icon.GetWidth()) / float(icon.GetHeight());
  if( r > 1.0  ) {
    w = size;
    h = int(w / r);
  }
  else {
    h = size;
    w = int(h * r);
  }
  wxImage img(icon.Scale(w,h));
  int i0 = 0;
  int j0 = 0;
  if( size > w )
    i0 = wxMax((size - w)/2, 0);
  if( size > h )
    j0 = wxMax((size - h)/2, 0);

  wxImage i(size,size);
  wxColour bgcolour(config->display_bgcolour);
  i.SetRGB(wxRect(0,0,size,size),bgcolour.Red(),bgcolour.Green(),bgcolour.Blue());
  i.Paste(icon,i0,j0);
  return i;
}


wxImage MuniIcon::GetCover(const wxImage& icon, const wxString& label, int nhdu) const
{
  wxSize size(config->icon_size,config->icon_size);
  wxBitmap bmp(size);
  wxMemoryDC mdc(bmp);
  wxGraphicsContext *gc = wxGraphicsContext::Create(mdc);
  if( gc ) {

    // clear canvas
    gc->SetBrush(wxBrush(wxColour(config->display_bgcolour)));
    gc->DrawRectangle(0,0,size.GetWidth(),size.GetHeight());

    // frame
    wxSize s(icon.GetSize());
    double xoff = (size.GetWidth() - s.GetWidth()) / 2;
    double yoff = (size.GetHeight() - s.GetHeight()) / 2;
    wxGraphicsBitmap b(gc->CreateBitmapFromImage(icon));
    gc->DrawBitmap(b,xoff,yoff,s.GetWidth(),s.GetHeight());

    // frame border
    gc->SetAntialiasMode(wxANTIALIAS_NONE);
    gc->SetBrush(*wxTRANSPARENT_BRUSH);
    wxPen bpen(colour);
    bpen.SetWidth(1);
    gc->SetPen(bpen);
    gc->DrawRectangle(0,0,size.GetWidth()-1,size.GetHeight()-1);

    // text
    gc->SetAntialiasMode(wxANTIALIAS_DEFAULT);

    // hardcoded font + inverse background
    wxFont font = wxSystemSettings::GetFont(wxSYS_DEFAULT_GUI_FONT);
    font.MakeSmaller();

    font.SetFamily(wxFONTFAMILY_ROMAN);
    const int f = 245;
    gc->SetFont(font,wxColour(f,f,f,222));

    const int g = 64;
    wxColour c(wxColour(g,g,g,96));
    wxBrush brush(c);
    gc->SetBrush(brush);
    wxPen pen(c,1);
    gc->SetPen(pen);

    if( nhdu > 1 ) {
      wxString label;
      label.Printf("%d+",int(nhdu));
      wxDouble w,h,ld,ll;
      gc->GetTextExtent(label,&w,&h,&ld,&ll);
      double x = size.GetWidth() - w;
      double y = 2;
      gc->DrawRectangle(x-2,y,w+4,h+2);
      gc->DrawText(label,x,y+1,0.0);
    }

    if( label != "" ) {
      wxDouble w,h,ld,ll;
      gc->GetTextExtent(label,&w,&h,&ld,&ll);
      double xoff = (size.GetWidth() - w) / 2;
      double yoff = size.GetHeight() - (h+2);

      wxString a(label);
      if( w > size.GetWidth() ) {
	double q = double(w) / double(a.Len());
	int m = size.GetWidth() / q + 0.5;
	a.Truncate(m-3);
	a += L"…";
	xoff = q/2+1;
      }
      gc->DrawRectangle(0,yoff,size.GetWidth(),h+2);
      gc->DrawText(a,xoff,yoff,0.0);
    }

    delete gc;
  }
  return bmp.ConvertToImage();
}

wxBitmap MuniIcon::Utf8Icon(double scale, const wxString& c, const wxString& label) const
{
  wxSize size(config->icon_size,config->icon_size);
  wxBitmap bmp(size);
  wxMemoryDC mdc(bmp);
  wxGraphicsContext *gc = wxGraphicsContext::Create(mdc);
  if( gc ) {

    // border
    gc->SetBrush(wxBrush(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOW)));
    gc->DrawRectangle(0,0,size.GetWidth(),size.GetHeight());
    double width = size.GetWidth();
    double height = size.GetHeight();
    double x = (size.GetWidth() - width) / 2;
    double y = size.GetHeight() - height - 2;
    gc->SetBrush(*wxTRANSPARENT_BRUSH);
    gc->SetPen(wxPen(colour));
    gc->DrawRectangle(1,1,width-2,height-2);

    if( label != "" ) {
      wxFont font(*wxNORMAL_FONT);
      font.MakeSmaller();
      gc->SetFont(font,wxSystemSettings::GetColour(wxSYS_COLOUR_GRAYTEXT));

      wxDouble w,h,ld,ll;
      gc->GetTextExtent(label,&w,&h,&ld,&ll);
      double xoff = (size.GetWidth() - w) / 2;
      double yoff = size.GetHeight() - (h+2);

      wxString a(label);
      if( w > size.GetWidth() ) {
	double q = double(w) / double(a.Len());
	int m = size.GetWidth() / q + 0.5;
	a.Truncate(m-3);
	a += L"…";
	xoff = q/2+1;
      }
      gc->DrawText(a,xoff,yoff,0.0);
    }

    wxFont font(*wxNORMAL_FONT);
    font.MakeBold();
    font.Scale(scale);
    gc->SetFont(font,colour);
    wxDouble w,h,ld,ll;
    gc->GetTextExtent(c,&w,&h,&ld,&ll);
    wxDouble xoff = (size.GetWidth() - w) / 2;
    wxDouble yoff = (size.GetHeight() - h) / 2;
    xoff = (width - w) / 2;
    yoff = (height - h) / 2;
    gc->DrawText(c,x+xoff,y+yoff,0.0);
    delete gc;
  }
  return bmp;
}



// obsolete

wxImage MuniIcon::BulletIcon(const wxSize& size, const wxColour& c)
{
  int x0 = size.GetWidth()/2;
  int y0 = size.GetHeight()/2;
  int radius = 20;

  wxColour colour(c.Red(),c.Green(),c.Blue(),255);
  int npix = size.GetWidth()*size.GetHeight();
  unsigned char *rgb = (unsigned char *) malloc(3*npix);

  wxColour cb = wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOW);
  wxImage image(size.GetWidth(),size.GetHeight(),rgb);
  for(int i = 0; i < image.GetWidth(); i++)
    for(int j = 0; j < image.GetHeight(); j++)
      image.SetRGB(i,j,cb.Red(),cb.Green(),cb.Blue());

  wxBitmap bitmap(image);
  wxMemoryDC dc(bitmap);
  wxGraphicsContext *gc = wxGraphicsContext::Create(dc);
  if( gc ) {
    gc->SetPen(wxPen(colour));
    gc->SetBrush(wxBrush(colour));
    wxGraphicsBrush b =
      gc->CreateRadialGradientBrush(x0-radius/8,y0-radius/8,x0+radius/8,y0+radius/8,
				    3*radius/5,*wxWHITE,colour);
    gc->SetBrush(b);
    gc->DrawEllipse(x0-radius/2,y0-radius/2,radius,radius);
    delete gc;
  }
  dc.SelectObjectAsSource(wxNullBitmap);
  return bitmap.ConvertToImage();
}


wxImage MuniIcon::ListIcon(const wxImage& icon, int size, const wxColour& c)
{
  wxASSERT(icon.IsOk() && icon.GetHeight() > 0);
  wxSize s;
  wxPoint p;
  double r = double(icon.GetWidth()) / double(icon.GetHeight());
  if( icon.GetWidth() > icon.GetHeight() ) {
    float x = r > 0.0 ? size/r : 1;
    if( x < 1 ) x = 1;
    s = wxSize(size,int(x));
    p = wxPoint(0,(size-s.GetHeight())/2);
  }
  else {
    float x = size*r;
    if( x < 1 ) x = 1;
    s = wxSize(int(x),size);
    p = wxPoint((size-s.GetWidth())/2,0);
  }
  wxImage img = icon.Scale(s.GetWidth(),s.GetHeight(),wxIMAGE_QUALITY_HIGH);

  for(int i = 0; i < img.GetWidth(); i++)
    for(int j = 0; j < img.GetHeight(); j++) {
      unsigned char r = (c.Red()*img.GetRed(i,j))/255;
      unsigned char g = (c.Green()*img.GetGreen(i,j))/255;
      unsigned char b = (c.Blue()*img.GetBlue(i,j))/255;
      img.SetRGB(i,j,r,g,b);
    }

  return img.Size(wxSize(size,size),p);
}
