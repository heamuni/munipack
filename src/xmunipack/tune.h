/*

  xmunipack - tune window

  Copyright © 2021-2022 F.Hroch (hroch@physics.muni.cz)

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


#ifndef _XMUNIPACK_TUNE_H
#define _XMUNIPACK_TUNE_H

#include "plot.h"
#include "fits.h"
#include "event.h"
#include <wx/wx.h>
#include <wx/notebook.h>



class MuniHistoThread: public wxThread
{
  wxEvtHandler *handler;
  FitsArray array;
  double qrange;

  ExitCode Entry();

public:
  MuniHistoThread(wxEvtHandler *, const FitsArray&, double);

};

class MuniTune: public wxFrame
{
  const FitsArray array;
  double black, sense;
  int itt;
  int palette;
  bool inverse;
  double saturation, x_white, y_white, meso_level, meso_width, qrange;
  bool undo;
  bool nitevision;

  MuniPlotHisto *plthist;
  MuniPlotItt *pltitt;
  MuniPlotPalette *pltpal;
  MuniPlotSaturation *pltsat;
  MuniPlotNite *pltnite;
  wxRadioButton *itt_line, *itt_asinh, *itt_tanh, *itt_snlike, *itt_sqr, *itt_photo;
  wxRadioButton *pal_grey, *pal_sepia, *pal_royal, *pal_colour, *pal_highlight,
    *pal_rainbow, *pal_madness, *pal_vga, *pal_aips0;
  wxRadioButton *wp_d50, *wp_d65, *wp_tungsten, *wp_flour, *wp_led, *wp_equal;
  wxCheckBox *invcheck, *nitecheck;
  wxRadioButton *rb99, *rb999, *rbfull;
  wxSlider *slider_sat;
  wxTextCtrl *entry_sense, *entry_black, *entry_meso_level, *entry_meso_width;
  MuniHistoThread *thread_histo;
  bool shutdown;

  void ScaleUpdated();
  void NiteUpdated();
  wxString ToString(const char *, double x) const;
  void SetIttRadio();
  void SetWhitePointRadio();
  void RenderHisto();

  void OnEscape(wxKeyEvent&);
  void OnApplyNite(wxCommandEvent&);
  void OnEnterMesoLevel(wxCommandEvent&);
  void OnEnterMesoWidth(wxCommandEvent&);
  void OnSat(wxCommandEvent&);
  void OnPal(wxCommandEvent&);
  void OnCheckInverse(wxCommandEvent&);
  void OnCheckNite(wxCommandEvent&);
  void OnClose(wxCloseEvent&);
  void OnItt(wxCommandEvent&);
  void OnWhitePoint(wxCommandEvent&);
  void OnQrange(wxCommandEvent&);
  void OnApply(wxCommandEvent&);
  void OnUndo(wxCommandEvent&);
  void OnEnterBlack(wxCommandEvent&);
  void OnEnterSense(wxCommandEvent&);
  void OnUpdateNite(wxUpdateUIEvent&);
  void OnUpdateNiteApply(wxUpdateUIEvent&);
  void OnUpdateUndo(wxUpdateUIEvent&);
  void OnUpdateApply(wxUpdateUIEvent&);
  void OnHistoFinish(MuniHistoEvent&);
  void FeedScaleTab(wxPanel *);
  void FeedIttTab(wxPanel *);
  void FeedPaletteTab(wxPanel *);
  void FeedColourTab(wxPanel *);
  void FeedNiteTab(wxPanel *);

public:
  MuniTune(wxWindow *, wxWindowID, const wxPoint&, const wxSize&, int,
		const FitsArray&, const FitsTone&, const FitsItt&,
		const FitsPalette&);
  MuniTune(wxWindow *, wxWindowID, const wxPoint&, const wxSize&, int,
		const FitsArray&, const FitsTone&, const FitsItt&,
		const FitsColour&);
  virtual ~MuniTune();

  void SetScale(double, double);

};

#endif
