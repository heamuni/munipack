/*

  xmunipack - definition of events


  Copyright © 2009-2012, 2019-22 F.Hroch (hroch@physics.muni.cz)

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

wxDEFINE_EVENT(EVT_SLEW, MuniSlewEvent);
wxDEFINE_EVENT(EVT_CLICK, MuniClickEvent);
wxDEFINE_EVENT(xEVT_RENDER, MuniRenderEvent);
wxDEFINE_EVENT(EVT_ZOOM, MuniZoomEvent);
wxDEFINE_EVENT(EVT_TUNE, MuniTuneEvent);
wxDEFINE_EVENT(EVT_TUNE_SCALE, MuniTuneScaleEvent);
wxDEFINE_EVENT(EVT_TUNE_ITT, MuniTuneIttEvent);
wxDEFINE_EVENT(EVT_TUNE_PAL, MuniTunePalEvent);
wxDEFINE_EVENT(EVT_TUNE_COLOUR, MuniTuneColourEvent);
wxDEFINE_EVENT(EVT_TUNE_NITE, MuniTuneNiteEvent);
wxDEFINE_EVENT(EVT_ASTROMETRY, MuniAstrometryEvent);
wxDEFINE_EVENT(EVT_PHOTOMETRY, MuniPhotometryEvent);
wxDEFINE_EVENT(EVT_CONFIG_UPDATED, MuniConfigEvent);
wxDEFINE_EVENT(EVT_SIZE_CHANGED, MuniSizeChangedEvent);
wxDEFINE_EVENT(EVT_FITS_STREAM, FitsStreamEvent);
wxDEFINE_EVENT(EVT_META_OPEN, MetaOpenEvent);
wxDEFINE_EVENT(EVT_META_PROGRESS, MetaProgressEvent);
wxDEFINE_EVENT(EVT_DIR_SCAN, DirScanEvent);
wxDEFINE_EVENT(EVT_NAVIGATION, NavigationEvent);
wxDEFINE_EVENT(EVT_RENDER, MuniRenderEvent);
wxDEFINE_EVENT(EVT_FILELOAD, wxCommandEvent);
wxDEFINE_EVENT(EVT_FINISH_DIALOG, wxCommandEvent);
wxDEFINE_EVENT(EVT_TOOL_FINISH, wxCommandEvent);
wxDEFINE_EVENT(EVT_LIST_ACTIVATED, wxNotifyEvent);
wxDEFINE_EVENT(EVT_DRAW, MuniDrawEvent);
wxDEFINE_EVENT(EVT_HISTO, MuniHistoEvent);
