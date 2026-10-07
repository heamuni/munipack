/*

  xmunipack - config

  Copyright © 2012-5, 2019-22 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _XMUNIPACK_CONFIG_H_
#define _XMUNIPACK_CONFIG_H_

#include <wx/config.h>
#include <wx/animate.h>
#include <wx/artprov.h>


class MuniArtIcons
{
  const wxArtClient client;
  const wxSize size;

 public:
  MuniArtIcons(const wxArtClient& =wxART_OTHER, const wxSize& =wxDefaultSize);
  wxBitmap Icon(const wxArtID&) const;

};


class MuniConfig: public wxConfig
{
public:
  MuniConfig(const wxString&);
  virtual ~MuniConfig();

  wxSize browser_size, view_size, /*console_size,*/colourization_size,
    header_size, aphot_size;
  int icon_size;
  int view_controls,
    /*detail_show,*/ detail_zoom,detail_scale;
  int display_pal,display_cootype,display_qphtype;
  bool display_palinv;
  int browser_labeltype,browser_sorttype, browser_iconlist;
  bool astrometry_fullmatch;
  int astrometry_minmatch, astrometry_maxmatch;
  double astrometry_sig, astrometry_fsig;
  double find_fwhm, find_thresh;
  int aphot_naper, aphot_rmin, aphot_rmax, aphot_zoom, aphot_saper;
  bool aphot_spiral, aphot_snap;
  wxString astrometry_proj, astrometry_units;
  bool browser_reverse;
  int icon_menu_width,icon_menu_height;
  bool browser_toolbar, view_toolbar/*, console_wrap*/;
  bool display_legend, display_sources;
  wxString dirmask, rawmask, browser_fitsmask,
    fits_key_object, fits_key_dateobs,fits_key_exptime,
    fits_key_filter, fits_key_observer, fits_key_gain, fits_key_area,
    browser_labelkey, browser_sortkey,display_colourspace,display_bgcolour,
    cdatafile, phsystemfile, munipack_html_dir;
  wxString backup_suffix;
  wxString default_symbol, head_symbol, image_symbol, table_symbol;
  wxIcon munipack_icon;
  wxImage moon_56frames;
  wxAnimation throbber;

  // Magnifier
  int magnifier_scale;

  static wxString FindIconPath(const wxString&);
  static wxIcon LoadIcon(const wxString&);
  static wxImage LoadImage(const wxString&);
  static wxString FindDataDir(const wxString&);
  static wxString FindHtmlDir(const wxString&);

};

#endif
