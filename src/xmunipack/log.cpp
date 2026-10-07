
/*

  xmunipack - log


  Copyright © 2009-2011 F.Hroch (hroch@physics.muni.cz)

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

  $Id$

  Note.

   The wxLogWindow doesn't supports set of window's size.

*/


#include "xmunipack.h"
#include <wx/wx.h>
#include <wx/log.h>
#include <wx/event.h>

MuniLog::MuniLog(wxWindow *w): 
  wxLogWindow(w,wxT("Log Messages"),false),topwin(w),visible(false)
{
  SetActiveTarget(this);
}

MuniLog::~MuniLog() {}


bool MuniLog::OnFrameClose(wxFrame *frame)
{
  wxCommandEvent e(wxEVT_COMMAND_MENU_SELECTED,ID_LOG);
  //  wxPostEvent(topwin,e);
  wxQueueEvent(topwin,e.Clone());
  return true;
}

bool MuniLog::Visibility()
{
  visible = ! visible;
  Show(visible);
  return visible;
}
