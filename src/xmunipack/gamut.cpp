/*

  This is helper for xmunipack.

  Creates color gamut diagram.

  Copyright © 1997-2010 F.Hroch (hroch@physics.muni.cz)

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

  

  Usage:

  * compile as

     g++ `wx-config --cxxflags` gamut.cpp color.cpp `wx-config --libs`

  * get cccie64.cvs from http://cvrl.ioo.ucl.ac.uk/index.htm
    (look for: Chromaticity coordinates 
               CIE 1964 10-deg chromaticity coordinates)

  * drop comas: tr ',' ' ' < cccie64.csv > cccie64.dat

  * run this code:
      
     a.out < cccie64.dat

  * output is stored as gamut.png

  * creates Luv chromaticity diagram:
      http://en.wikipedia.org/wiki/CIELUV_color_space

 */

#include "fits.h"
#include <wx/wx.h>
#include <iostream>
#include <algorithm>

using namespace std;


class DerivedApp : public wxApp
{
public:
  virtual bool OnInit();
};

IMPLEMENT_APP(DerivedApp)

bool DerivedApp::OnInit()
{
  // set options
  int width = 500;
  int height = 500;
  float xr = 1.7;
  float yr = 1.6;

  // input gamut
  vector<double> wave,xg,yg;
  while( ! cin.eof() ) {
    double x,y,z,t;
    cin >> x >> y >> z >> t;
    wave.push_back(x);
    xg.push_back(y);
    yg.push_back(z);
    //    wxLogDebug(_("%f %f %f"),x,y,z);
  }
  
  wxPoint points[wave.size()];
  for(size_t i = 0; i < wave.size(); i++) {
    float x = xg[i];
    float y = yg[i];
    float X = 0, Y = 0, Z = 0;
    if( y != 0.0 ) { // don't divide by zero!
      X = x/y;
      Y = 1.0;
      Z = (1.0 - x - y)/y;  
    }
    float t = X + 15.0*Y + 3.0*Z;
    float u = 4.0*X/t;
    float v = 9.0*Y/t;
    points[i] = wxPoint(xr*width*u,height*(1.0-yr*v));
  }

  wxRegion region(wave.size(),points);

  wxInitAllImageHandlers();


  wxImage image(width,height);
  //  image.SetAlpha();
  wxBitmap bitmap(image);
  wxMemoryDC mdc(bitmap);

  mdc.SetBrush(*wxBLACK_BRUSH);
  mdc.Clear();

  mdc.DestroyClippingRegion();
  mdc.SetClippingRegion(region);

  //  int xoff = size.GetWidth()/2;
  //  int yoff = size.GetHeight()/2;

  //  wxLogDebug(_("%d %d %d %d"),xoff,yoff,dx,dy);

  FitsColor color;
  for(int j = 0; j < height; j++)
    for(int i = 0; i < width; i++) {
      //      int j = height - jj - 1;

      float u = float(i)/float(width)/xr;
      float v = float(height - j - 1)/float(height)/yr;
      //      y = y*yr;

      float X = 0, Y = 0, Z = 0;
//       if( y != 0.0 ) { // don't divide by zero!
// 	X = x/y;
// 	Y = 1.0;
// 	Z = (1.0 - x - y)/y;  
//       }
//       float t = X + 15.0*Y + 3.0*Z;
//       float u = 4.0*X/t;
//       float v = 9.0*Y/t;

      float L = 70;
      //      u = 13*L*(u - 0.2);
      //      v = 13*L*(v - 0.45);
      //      float X,Y,Z;
      
      Y = powf((L+16.0)/116.0,3.0);
      X = Y*9.0*u/4.0/v;
      Z = Y*(12.0 - 3.0*u - 20.0*v)/4.0/v;
      
      //      color.Luv_XYZ(L,500*u,500*v,1.0,X,Y,Z);
      //      wxLogDebug(_("%d %d %f %f %f %f %f"),i,j,x,y,X,Y,Z);
      //      if( X > 0.0 && Z > 0.0 && X <= 1.1 && Z < 1.8 ) {

      unsigned char R,G,B;
      color.XYZ_sRGB(X,Y,Z,R,G,B);
      /*
      R = X;
      G = Y;
      B = Z;
      */
      //      if( R >= 0 && G >= 0 && B >= 0 && R < 256 && G < 256 && B < 256 ) {
      mdc.SetPen(wxColor(R,G,B));
      //      mdc.DrawRectangle(i,j,1,1);
      mdc.DrawPoint(i,j);
	//      }
    }
  mdc.SetBrush(*wxTRANSPARENT_BRUSH);
  mdc.SetPen(wxPen(*wxLIGHT_GREY,5));
  mdc.DrawPolygon(wave.size(),points);
  //  mdc.SetClippingRegion(0,0,width,height);
  //  mdc.DrawRectangle(0,0,width,height);
  mdc.SelectObjectAsSource(wxNullBitmap);

  wxImage img(bitmap.ConvertToImage());
  img.SetAlpha();
  for(int i = 0; i < img.GetWidth(); i++)
    for(int j = 0; j < img.GetHeight(); j++)
      if( img.GetRed(i,j) == 255 && img.GetGreen(i,j) == 255  && img.GetBlue(i,j) == 255 )
	img.SetAlpha(i,j,0);
      else
	img.SetAlpha(i,j,255);
  img.SaveFile(wxT("gamut.png"));

  return 0;
}
