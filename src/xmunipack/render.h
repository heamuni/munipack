/*
  xmunipack - image render

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

#ifndef _XMUNIPACK_RENDER_H
#define _XMUNIPACK_RENDER_H

#include "enum.h"
#include "fits.h"
#include "fitsdisplay.h"
#include "event.h"

#include <wx/wx.h>
#include <wx/thread.h>
#include <vector>

class MuniThreadRender: public wxThread
{
  wxEvtHandler *handler;
  wxMutex *kerberos;
  wxCondition *condition;
  int *count;

  bool stop, cdata_valid;
  double zoom;
  int shrink, iwidth, iheight, cwidth, cheight, cdepth, rid, ntiles;
  const int fwidth, fheight, fdepth;
  const float *fdata;
  unsigned char *idata;
  float *cdata;
  const FitsGeometry *geometry;
  const FitsDisplay *display;
  wxRect exposed;
  std::vector<wxRect> rects;

  ExitCode Entry();
  void Zoom(int,int,const unsigned char *sdata,int,int,int,int,int,
	    int, int, unsigned char *);

 public:
  MuniThreadRender(wxEvtHandler *, wxMutex *, wxCondition *,
		   int *, const FitsArray&, const FitsGeometry *,
		   const FitsDisplay *, int);

  int GetRenderedTiles() const;
  void SetRectangles(const wxRect&, const std::vector<wxRect>&);
  void SetImageData(int,int,unsigned char *);
  void SetZoom(double);
  void SetShrinkArray(int,int,int,float*,bool);

  void Stop();

};


class MuniSignalThread: public wxThread
{
  wxEvtHandler *handler;
  wxCondition *condition;
  int rid;
  bool stop, stop_wait;
  wxStopWatch sw;
  ExitCode Entry();

public:
  MuniSignalThread(wxEvtHandler *, wxCondition *, int);
  virtual ~MuniSignalThread();

  void Stop();
  void StopAndWait();

};

class MuniDisplayRender
{
  wxEvtHandler *handler;
  wxMutex kerberos;
  wxCondition condition;

  bool stop, cdata_valid;
  double zoom;
  int ncpu, count, iwidth, iheight, cwidth, cheight, cdepth, ntiles;
  unsigned char *idata;
  float *cdata;
  const FitsArray fitsimage;
  const FitsGeometry geometry;
  FitsDisplay display;
  wxRect exposed;
  MuniSignalThread *signalling;
  std::vector<MuniThreadRender *> renders;

  std::vector<wxRect> Lines() const;

public:
  MuniDisplayRender(wxEvtHandler *, const FitsArray&);
  MuniDisplayRender(const MuniDisplayRender&);
  MuniDisplayRender& operator = (const MuniDisplayRender&);
  virtual ~MuniDisplayRender();

  void SetImageData(int,int,unsigned char *);
  void SetZoom(double);
  void SetTone(const FitsTone&,const FitsItt&, const FitsPalette&, const FitsColour&);
  void SetExposed(const wxRect&);

  wxThreadError Run(int);
  void Stop();
  void Wait();
  void StopAndWait();
  bool IsRunning() const;

};

#endif
