/*

  xmunipack - list threads

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

#include "event.h"
#include "dirscan.h"
#include <wx/wx.h>
#include <wx/dir.h>
#include <wx/thread.h>
#include <wx/filename.h>
#include <wx/filesys.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


// -- List file Scanner

DirScan::DirScan(wxEvtHandler *h, const wxString& p,
		 const wxString& f, size_t i):
  wxThread(wxTHREAD_DETACHED),handler(h),
  path(p == ""? "." : p), filter(f),id(i)
{
  wxASSERT(handler && ! path.IsEmpty() );
  wxLogDebug("DirScan::DirScan(): `"+path+"' "+f);
}


wxThread::ExitCode DirScan::Entry()
{
  wxDir dir(path);
  if( dir.IsOpened() ) {

    wxStopWatch sw;
    sw.Start();
    wxArrayString files;
    wxString filename;
    bool cont = dir.GetFirst(&filename,filter, wxDIR_FILES);
    while ( cont ) {
      if( TestDestroy() ) break;

      wxFileName file(path,filename);
      if( file.IsFileReadable() )
	files.Add(wxFileSystem::FileNameToURL(file));
      else
	wxLogDebug("DirScan: `"+filename+"' is unreadable.");

      //      wxLogDebug("%d",int(sw.Time()));
      if( sw.Time() > 100 ) {
	DirScanEvent ev(EVT_DIR_SCAN);
	ev.finish = false;
	ev.id = id;
	ev.files = files;
	wxQueueEvent(handler,ev.Clone());
	sw.Start();
	files.Empty();
      }

      cont = dir.GetNext(&filename);
      //      Sleep(100);
    }

    DirScanEvent ev(EVT_DIR_SCAN);
    ev.finish = true;
    ev.id = id;
    ev.files = files;
    wxQueueEvent(handler,ev.Clone());
  }

  //  wxLogDebug("DirScan::Entry()  FINISH");

  return (wxThread::ExitCode) 0;
}
