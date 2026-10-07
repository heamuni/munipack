/*

  xmunipack - config

  Copyright © 2012-2025 F.Hroch (hroch@physics.muni.cz)

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

#include "mconfig.h"
#include "enum.h"
#include "icon.h"
#include <wx/wx.h>
#include <wx/config.h>
#include <wx/stdpaths.h>
#include <wx/artprov.h>
#include <wx/filename.h>
#include <wx/listctrl.h>

// config
#define MUNIBROWSE_WIDTH     "MuniBrowse_width"
#define MUNIBROWSE_HEIGHT    "MuniBrowse_height"
#define MUNIBROWSE_ICONLIST  "MuniBrowse_iconlist"
#define MUNIBROWSE_LABELTYPE "MuniBrowse_labeltype"
#define MUNIBROWSE_LABELKEY  "MuniBrowse_labelkey"
#define MUNIBROWSE_SORTTYPE  "MuniBrowse_sorttype"
#define MUNIBROWSE_SORTKEY   "MuniBrowse_sortkey"
#define MUNIBROWSE_REVERSE   "MuniBrowse_reverse"
#define MUNIBROWSE_FITSMASK  "MuniBrowse_fitsmask"
#define MUNIBROWSE_TOOLBAR   "MuniBrowse_toolbar"
#define MUNIVIEW_WIDTH       "MuniView_width"
#define MUNIVIEW_HEIGHT      "MuniView_height"
#define MUNIVIEW_TOOLBAR     "MuniView_toolbar"
#define MUNIVIEW_CONTROLS    "MuniView_controls"
#define MUNIMAGNIFIER_SCALE  "MuniMagnifier_scale"
/*
#define MUNICONSOLE_WIDTH    "MuniConsole_width"
#define MUNICONSOLE_HEIGHT   "MuniConsole_height"
#define MUNICONSOLE_WRAP     "MuniConsole_wrap"
*/
#define MUNIHEADER_WIDTH     "MuniHeader_width"
#define MUNIHEADER_HEIGHT    "MuniHeader_height"
#define ASTROMETRY_FULLMATCH "Astrometry_fullmatch"
#define ASTROMETRY_MINMATCH  "Astrometry_minmatch"
#define ASTROMETRY_MAXMATCH  "Astrometry_maxmatch"
#define ASTROMETRY_UNITS     "Astrometry_units"
#define ASTROMETRY_PROJ      "Astrometry_proj"
#define ASTROMETRY_SIG       "Astrometry_sig"
#define ASTROMETRY_FSIG      "Astrometry_fsig"
#define FIND_FWHM            "Find_fwhm"
#define FIND_THRESH          "Find_thresh"
#define APHOT_NAPER          "Aphot_naper"
#define APHOT_RMIN           "Aphot_rmin"
#define APHOT_RMAX           "Aphot_rmax"
#define APHOT_ZOOM           "Aphot_zoom"
#define APHOT_SAPER          "Aphot_saper"
#define APHOT_SPIRAL         "Aphot_spiral"
#define APHOT_SNAP           "Aphot_snap"
#define APHOT_WIDTH          "Aphot_width"
#define APHOT_HEIGHT         "Aphot_height"
#define ICON_SIZE            "Icon_size"
#define DISPLAY_PAL          "Display_pal"
#define DISPLAY_PALINV       "Display_palinv"
#define DISPLAY_COOTYPE      "Display_cootype"
#define DISPLAY_QPHTYPE      "Display_qphtype"
#define DISPLAY_LEGEND       "Display_legend"
#define DISPLAY_SOURCES      "Display_sources"
#define DISPLAY_BGCOLOUR     "Display_bgcolour"
#define FITS_KEY_OBJECT      "Fits_key_object"
#define FITS_KEY_OBSERVER    "Fits_key_observer"
#define FITS_KEY_EXPTIME     "Fits_key_exptime"
#define FITS_KEY_GAIN        "Fits_key_gain"
#define FITS_KEY_AREA        "Fits_key_area"
#define FITS_KEY_FILTER      "Fits_key_filter"
#define FITS_KEY_DATEOBS     "Fits_key_dateobs"
#define COLOUR_DISPLAY       "Display_colourspace"
#define CDATAFILE            "Colourspace_data"
#define PHSYSTEMFILE         "Photometric_systems"


MuniConfig::MuniConfig(const wxString& config_file):
  wxConfig("Xmunipack","",config_file,"",wxCONFIG_USE_LOCAL_FILE),
  dirmask("*.fits*;*.fit*;*.fts*;*.FITS*;*.FTS*;*.FIT*"),
  rawmask("*.cr2;*.crw;*.CR2;*.CRW;*.RAW"),
  backup_suffix("~")
{
  // backup suffix
  wxString simple_backup_suffix;
  if( wxGetEnv("SIMPLE_BACKUP_SUFFIX",&simple_backup_suffix) )
    backup_suffix = simple_backup_suffix;

  // initial browser window size
  int width,height,w,h,k;
  wxSize display_size = wxGetDisplaySize();
  width = display_size.GetWidth();
  height = display_size.GetHeight();

  wxString val;

  Read(MUNIBROWSE_WIDTH,&w,3*width/5);
  Read(MUNIBROWSE_HEIGHT,&h,3*height/5);
  browser_size = wxSize(w,h);
  Read(MUNIBROWSE_FITSMASK,&browser_fitsmask,"*.fit*");
  Read(MUNIBROWSE_TOOLBAR,&browser_toolbar,0);

  Read(MUNIVIEW_WIDTH,&w,4*width/5);
  Read(MUNIVIEW_HEIGHT,&h,4*height/5);
  Read(MUNIVIEW_TOOLBAR,&view_toolbar,1);
  Read(MUNIVIEW_CONTROLS,&view_controls,1);
  view_size = wxSize(w,h);

  Read(MUNIMAGNIFIER_SCALE,&magnifier_scale,4);


  /*
  Read(MUNICONSOLE_WIDTH,&w,width/2);
  Read(MUNICONSOLE_HEIGHT,&h,height/2);
  console_size = wxSize(w,h);
  Read(MUNICONSOLE_WRAP,&console_wrap,0);
  */

  wxFont font = wxSystemSettings::GetFont(wxSYS_ANSI_FIXED_FONT);
  wxSize sfont = font.GetPixelSize();
  Read(MUNIHEADER_WIDTH,&w,85*sfont.GetWidth());
  Read(MUNIHEADER_HEIGHT,&h,45*sfont.GetHeight());
  header_size = wxSize(w,h);

  Read(ASTROMETRY_FULLMATCH,&astrometry_fullmatch,false);
  Read(ASTROMETRY_MINMATCH,&astrometry_minmatch,5);
  Read(ASTROMETRY_MAXMATCH,&astrometry_maxmatch,30);
  Read(ASTROMETRY_PROJ,&astrometry_proj,"GNOMONIC");
  Read(ASTROMETRY_SIG,&astrometry_sig,1.0);
  Read(ASTROMETRY_FSIG,&astrometry_fsig,5.0);
  Read(ASTROMETRY_UNITS,&astrometry_units,"arcsec");

  Read(FIND_FWHM,&find_fwhm,3.0);
  Read(FIND_THRESH,&find_thresh,7.0);

  Read(APHOT_NAPER,&aphot_naper,12);
  Read(APHOT_RMIN,&aphot_rmin,20);
  Read(APHOT_RMAX,&aphot_rmax,30);
  Read(APHOT_ZOOM,&aphot_zoom,4);
  Read(APHOT_SAPER,&aphot_saper,0);
  Read(APHOT_SPIRAL,&aphot_spiral,1);
  Read(APHOT_SNAP,&aphot_snap,1);
  Read(APHOT_WIDTH,&w,-1);
  Read(APHOT_HEIGHT,&h,-1);
  aphot_size = wxSize(w,h);

  Read(ICON_SIZE,&icon_size,150);

  Read(MUNIBROWSE_ICONLIST,&k,0);
  switch(k){
  case 0:  browser_iconlist = wxLC_ICON; break;
  case 1:  browser_iconlist = wxLC_REPORT; break;
  default: browser_iconlist = wxLC_ICON;
  }

  Read(MUNIBROWSE_LABELTYPE,&k,0);
  switch(k){
  case 0:  browser_labeltype = ID_LABEL_NOLABEL;  break;
  case 1:  browser_labeltype = ID_LABEL_FILENAME; break;
  case 2:  browser_labeltype = ID_LABEL_OBJECT;   break;
  case 3:  browser_labeltype = ID_LABEL_DATE;     break;
  case 4:  browser_labeltype = ID_LABEL_TIME;     break;
  case 5:  browser_labeltype = ID_LABEL_FILTER;   break;
  case 6:  browser_labeltype = ID_LABEL_EXPOSURE; break;
  case 7:  browser_labeltype = ID_LABEL_KEY;      break;
  default: browser_labeltype = ID_LABEL_NOLABEL;
  }

  Read(MUNIBROWSE_LABELKEY,&browser_labelkey,"");

  Read(MUNIBROWSE_SORTTYPE,&k,0);
  switch(k){
  case 0:  browser_sorttype = ID_SORT_UNSORT;   break;
  case 1:  browser_sorttype = ID_SORT_FILENAME; break;
  case 2:  browser_sorttype = ID_SORT_OBJECT;   break;
  case 3:  browser_sorttype = ID_SORT_DATEOBS;  break;
  case 4:  browser_sorttype = ID_SORT_FILTER;   break;
  case 5:  browser_sorttype = ID_SORT_EXPOSURE; break;
  case 6:  browser_sorttype = ID_SORT_SIZE;     break;
  case 7:  browser_sorttype = ID_SORT_KEY;      break;
  default: browser_sorttype = ID_SORT_FILENAME;
  }

  Read(MUNIBROWSE_SORTKEY,&browser_sortkey,"");
  Read(MUNIBROWSE_REVERSE,&browser_reverse,0);

  Read(DISPLAY_PAL,&val,wxEmptyString);
  display_pal = PAL_GREY;
  for(int i = PAL_FIRST+1; i < PAL_LAST; i++) {
    if( val == FitsPalette::Type_str(i) )
      display_pal = i;
  }

  Read(DISPLAY_QPHTYPE,&val,wxEmptyString);
  display_qphtype = UNIT_COUNT;
  for(int i = UNIT_FIRST+1; i < UNIT_LAST; i++) {
    if( val == FitsValue::Label_str(i) )
      display_qphtype = i;
  }

  Read(DISPLAY_COOTYPE,&val,wxEmptyString);
  display_cootype = COO_EQDEG;
  for(int i = COO_FIRST+1; i < COO_LAST; i++) {
    if( val == FitsCoo::Label_str(i) )
      display_cootype = i;
  }

  Read(DISPLAY_PALINV,&display_palinv, 0);
  Read(DISPLAY_LEGEND, &display_legend, 0);
  Read(DISPLAY_SOURCES, &display_sources, 0);

  Read(FITS_KEY_OBJECT, &fits_key_object, "OBJECT");
  Read(FITS_KEY_OBSERVER, &fits_key_observer, "OBSERVER");
  Read(FITS_KEY_DATEOBS, &fits_key_dateobs, "DATE-OBS");
  Read(FITS_KEY_EXPTIME, &fits_key_exptime, "EXPTIME");
  Read(FITS_KEY_FILTER, &fits_key_filter, "FILTER");
  Read(FITS_KEY_GAIN, &fits_key_gain, "GAIN");
  Read(FITS_KEY_AREA, &fits_key_area, "AREA");

  Read(COLOUR_DISPLAY, &display_colourspace, "sRGB");
  Read(DISPLAY_BGCOLOUR, &display_bgcolour, "grey15");

  cdatafile = FindDataDir("ctable.dat");
  phsystemfile = FindDataDir("photosystems.fits");
  munipack_html_dir = FindHtmlDir("munipack.html");

  default_symbol = "?";
  head_symbol = L"🦉";
  image_symbol = L"❦";
  table_symbol = "#";

  munipack_icon = LoadIcon("munipack_icon.png");
  moon_56frames = LoadImage("moon_56frames.png");
  throbber = wxAnimation(FindIconPath("throbber.gif"));

}

MuniConfig::~MuniConfig()
{
  int k;

  Write(MUNIBROWSE_WIDTH,browser_size.GetWidth());
  Write(MUNIBROWSE_HEIGHT,browser_size.GetHeight());
  Write(MUNIBROWSE_FITSMASK,browser_fitsmask);
  Write(MUNIBROWSE_TOOLBAR,browser_toolbar);
  Write(MUNIVIEW_WIDTH,view_size.GetWidth());
  Write(MUNIVIEW_HEIGHT,view_size.GetHeight());
  Write(MUNIVIEW_TOOLBAR,view_toolbar);
  Write(MUNIVIEW_CONTROLS,view_controls);
  Write(MUNIMAGNIFIER_SCALE,magnifier_scale);
  /*
  Write(MUNICONSOLE_WIDTH,console_size.GetWidth());
  Write(MUNICONSOLE_HEIGHT,console_size.GetHeight());
  Write(MUNICONSOLE_WRAP,console_wrap);
  */
  Write(MUNIHEADER_WIDTH,header_size.GetWidth());
  Write(MUNIHEADER_HEIGHT,header_size.GetHeight());
  Write(ICON_SIZE,icon_size);

  switch(browser_iconlist){
  case wxLC_ICON:   k = 0; break;
  case wxLC_REPORT: k = 1; break;
  default: k = 0;
  }
  Write(MUNIBROWSE_ICONLIST,k);

  switch(browser_labeltype){
  case ID_LABEL_NOLABEL:  k = 0; break;
  case ID_LABEL_FILENAME: k = 1; break;
  case ID_LABEL_OBJECT:   k = 2; break;
  case ID_LABEL_DATE:     k = 3; break;
  case ID_LABEL_TIME:     k = 4; break;
  case ID_LABEL_FILTER:   k = 5; break;
  case ID_LABEL_EXPOSURE: k = 6; break;
  case ID_LABEL_KEY:      k = 7; break;
  default: k = 0;
  }
  Write(MUNIBROWSE_LABELTYPE,k);
  Write(MUNIBROWSE_LABELKEY,browser_labelkey);

  switch(browser_sorttype){
  case ID_SORT_UNSORT:   k = 0; break;
  case ID_SORT_FILENAME: k = 1; break;
  case ID_SORT_OBJECT:   k = 2; break;
  case ID_SORT_DATEOBS:  k = 3; break;
  case ID_SORT_FILTER:   k = 4; break;
  case ID_SORT_EXPOSURE: k = 5; break;
  case ID_SORT_SIZE:     k = 6; break;
  case ID_SORT_KEY:      k = 7; break;
  default: k = 0;
  }
  Write(MUNIBROWSE_SORTTYPE,k);

  Write(MUNIBROWSE_SORTKEY,browser_sortkey);
  Write(MUNIBROWSE_REVERSE,browser_reverse);
  Write(ASTROMETRY_FULLMATCH,astrometry_fullmatch);
  Write(ASTROMETRY_MINMATCH,astrometry_minmatch);
  Write(ASTROMETRY_MAXMATCH,astrometry_maxmatch);
  Write(ASTROMETRY_PROJ,astrometry_proj);
  Write(ASTROMETRY_SIG,astrometry_sig);
  Write(ASTROMETRY_FSIG,astrometry_fsig);
  Write(ASTROMETRY_UNITS,astrometry_units);
  Write(FIND_FWHM,find_fwhm);
  Write(FIND_THRESH,find_thresh);

  Write(APHOT_NAPER,aphot_naper);
  Write(APHOT_RMIN,aphot_rmin);
  Write(APHOT_RMAX,aphot_rmax);
  Write(APHOT_ZOOM,aphot_zoom);
  Write(APHOT_SAPER,aphot_saper);
  Write(APHOT_SPIRAL,aphot_spiral);
  Write(APHOT_SNAP,aphot_snap);
  Write(APHOT_WIDTH,aphot_size.GetWidth());
  Write(APHOT_HEIGHT,aphot_size.GetHeight());

  Write(DISPLAY_PAL,FitsPalette::Type_str(display_pal));
  Write(DISPLAY_PALINV,display_palinv);
  Write(DISPLAY_LEGEND,display_legend);
  Write(DISPLAY_SOURCES,display_sources);
  Write(DISPLAY_COOTYPE,FitsCoo::Label_str(display_cootype));
  Write(DISPLAY_QPHTYPE,FitsValue::Label_str(display_qphtype));
  Write(FITS_KEY_OBJECT,fits_key_object);
  Write(FITS_KEY_OBSERVER,fits_key_observer);
  Write(FITS_KEY_DATEOBS,fits_key_dateobs);
  Write(FITS_KEY_EXPTIME,fits_key_exptime);
  Write(FITS_KEY_FILTER,fits_key_filter);
  Write(FITS_KEY_GAIN,fits_key_gain);
  Write(FITS_KEY_AREA,fits_key_area);
  Write(DISPLAY_BGCOLOUR,display_bgcolour);
  Write(COLOUR_DISPLAY,display_colourspace);

  wxLogDebug("MuniConfig::~MuniConfig()");
}

wxString MuniConfig::FindDataDir(const wxString& file)
{
  wxString dir;

  if( wxGetEnv("MUNIPACK_DATA_DIR",&dir) == false ) {
#ifdef MUNIPACK_DATA_DIR
    dir = wxString(MUNIPACK_DATA_DIR);
#endif
  }

  wxFileName path(dir,file);
  if( path.FileExists() )
    return path.GetFullPath();

  // fail-back
  wxPrintf("%s: an improper setup. Please, set MUNIPACK_DATA_DIR.\n",__FILE__);
  return "";
}

wxString MuniConfig::FindHtmlDir(const wxString& file)
{
  wxString dir;

  if( wxGetEnv("MUNIPACK_HTML_DIR",&dir) == false ) {
#ifdef MUNIPACK_HTML_DIR
    dir = wxString(MUNIPACK_HTML_DIR);
#endif
  }

  wxFileName path(dir,file);
  if( path.FileExists() )
    return path.GetFullPath();

  // fail-back
  /* The help is inactive currently. */
  //  wxPrintf("%s: an improper setup. Please, set MUNIPACK_HTML_DIR.\n",__FILE__);
  return "";
}


wxString MuniConfig::FindIconPath(const wxString& name)
{
  wxArrayString p;

#ifdef __WXDEBUG__
  p.Add("icons");
#endif

  wxPathList paths(p);
  wxString dir;

  if( wxGetEnv("MUNIPACK_ICON_DIR",&dir) )
    paths.Add(dir);

#ifdef MUNIPACK_ICON_DIR
  paths.Add(MUNIPACK_ICON_DIR);
#endif

  wxString fullpath = paths.FindValidPath(name);
  return fullpath;
}

wxIcon MuniConfig::LoadIcon(const wxString& name)
{
  wxString fullpath = FindIconPath(name);

  wxIcon icon;
  if( ! fullpath.IsEmpty() ) {

    //    wxLogDebug(fullpath);
    wxImage icona(fullpath);
    if( icona.Ok() ) {
      wxBitmap bitmap(icona);
      icon.CopyFromBitmap(bitmap);
    }
  }

  if( ! icon.IsOk() )
    icon = wxArtProvider::GetIcon(wxART_MISSING_IMAGE,wxART_OTHER,wxDefaultSize);

  return icon;
}

wxImage MuniConfig::LoadImage(const wxString& name)
{
  wxString fullpath = FindIconPath(name);

  wxImage image;
  if( ! fullpath.IsEmpty() && image.LoadFile(fullpath) ) {
    //    wxLogDebug(fullpath);
    return image;
  }
  else {
    wxBitmap b(wxArtProvider::GetBitmap(wxART_MISSING_IMAGE,wxART_OTHER,
					wxDefaultSize));
    return b.ConvertToImage();
  }
}
