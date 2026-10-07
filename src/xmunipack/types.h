/*

  xmunipack - basic types for drawing

  Copyright © 2012, 2019 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _XMUNIPACK_TYPES_H_
#define _XMUNIPACK_TYPES_H_


#include "enum.h"
#include <wx/wx.h>
#include <vector>


class MuniDrawFont: public wxObject
{
public:
  MuniDrawFont(const wxFont &f, const wxColour &c): font(f), colour(c) {}
  wxFont font;
  wxColour colour;
};

class MuniDrawPen: public wxObject
{
public:
  MuniDrawPen(const wxPen& p): pen(p) {}
  wxPen pen;
};

class MuniDrawBrush: public wxObject
{
public:
  MuniDrawBrush(const wxColour& c): brush(c) {}
  MuniDrawBrush(const wxBrush& b): brush(b) {}
  wxBrush brush;
};

class MuniDrawCircle: public wxObject
{
public:
  MuniDrawCircle(double xx, double yy, double rr):x(xx),y(yy),r(rr) {}
  double x,y,r;
};

class MuniDrawEllipse: public wxObject
{
public:
 MuniDrawEllipse(double xx, double yy, double aa, double ee, double ii):
  x(xx),y(yy),a(aa),e(ee),i(ii) {}
  double x,y,a,e,i;
};

class MuniDrawRing: public wxObject
{
public:
 MuniDrawRing(double xx, double yy, double r0, double r1):
  x(xx),y(yy),rin(r0),rout(r1) {}
  double x,y,rin,rout;
};


class MuniDrawLine: public wxObject
{
public:
  MuniDrawLine(double xx1, double yy1, double xx2, double yy2): x1(xx1),y1(yy1),x2(xx2),y2(yy2) {}
  double x1,y1,x2,y2;
};

class MuniDrawCross: public wxObject
{
public:
  MuniDrawCross(double xx, double yy, double rr):x(xx),y(yy),r(rr) {}
  double x,y,r;
};

class MuniDrawText: public wxObject
{
public:
  MuniDrawText(double xx, double yy, const wxString& t): x(xx),y(yy),angle(0.0),text(t) {}
  MuniDrawText(double xx, double yy, double a, const wxString& t): x(xx),y(yy),angle(a),text(t) {}
  double x,y,angle;
  wxString text;
};

class MuniDrawRectangle: public wxObject
{
public:
 MuniDrawRectangle(double xx, double yy, double ww, double hh): x(xx),y(yy),w(ww),h(hh) {}
  double x,y,w,h;
};

class MuniDrawBitmap: public wxObject
{
public:
 MuniDrawBitmap(const wxBitmap& bmp,double xx, double yy, double ww, double hh): bitmap(bmp),x(xx),y(yy),w(ww),h(hh) {}
  wxBitmap bitmap;
  double x,y,w,h;
};

/*
class MuniSvg: public wxXmlDocument
{
public:
  MuniSvg(const wxXmlDocument& t): wxXmlDocument(t) {}
  vector<MuniDrawBase *> GetDrawObjects() const;


};
*/

class MuniLayer {
public:
  MuniLayer(): id(ID_NULL) {}
  MuniLayer(int i): id(i) {}
  MuniLayer(int i, const std::vector<wxObject *> o): id(i),objects(o) {}

  bool IsOk() { return id != ID_NULL && objects.size() > 0; }
  int GetId() const { return id; }
  size_t IsEmpty() const { return objects.empty(); } // flag: to remove layer
  std::vector<wxObject *> GetObjects() const { return objects; }

private:

  int id;
  std::vector<wxObject *> objects;

};

class MuniStarLayer
{
  const wxColour gold,DarkOrange2;
  std::vector<wxObject *> objects;
  double hwhm;

 public:
  MuniStarLayer();
  MuniLayer GetLayer() const;
  bool IsEmpty() const { return objects.empty(); }
  void DrawObjects(const std::vector<double>&,const std::vector<double>&,
		   const std::vector<double>&);
  void SetHWHM(double x) { hwhm = x; }

};


#endif
