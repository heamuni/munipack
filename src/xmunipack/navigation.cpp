/*

  xmunipack - navigation for Browser

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

  Implementation notes:

  * The entry needs some kind of completion. I spend hours
    with activation AutoComplete(wxStringArray) which definitely
    does not work. AutoComplete(ComptereFactory*) has been left
    on future. The rest of code is commented. Generally standard
    command line habits interfere with TAB in GUI role.

  * Another attempt utilising a popup window on request is another
    attractive alternative to implement.


*/

#include "navigation.h"
#include <wx/wx.h>
#include <wx/dir.h>
//#include <wx/popupwin.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

#define NORMALIZE_FLAGS wxPATH_NORM_DOTS | wxPATH_NORM_TILDE  | wxPATH_NORM_ABSOLUTE

// MuniTextCtrl::MuniTextCtrl(wxWindow *w, const wxString& path, long style):
//   wxTextCtrl(w,wxID_ANY,path,wxDefaultPosition,wxDefaultSize,style)
// {
//   //  Bind(wxEVT_KEY_DOWN,&MuniTextCtrl::OnChar,this/*,cmd->GetId()*/);
//   //  Bind(wxEVT_TEXT,&MuniTextCtrl::OnText,this);

//   wxArrayString dirs = GetSubdirs(path);
//   bool x = AutoComplete(dirs);

// }

// void MuniTextCtrl::OnText(wxCommandEvent& event)
// {
//   if( event.GetString()  == "" ) return;
//   wxLogDebug("OnText: %d "+event.GetString(),IsModified());


//   wxFileName d(event.GetString());
//   if( d.DirExists() ) {
//     wxArrayString dirs = GetSubdirs(event.GetString());
//     bool x = AutoComplete(dirs);
//     wxLogDebug(".................. %d",x);
//   }
//   else {
//     wxArrayString dirs = GetSubdirs(event.GetString());
//     bool x = AutoComplete(dirs);
//     wxLogDebug(".................. %d",x);
//   }
//   wxLogDebug("****");
// }

// wxArrayString MuniTextCtrl::GetSubdirs(const wxString& path) const
// {
//   wxString a = path;
//   wxUniChar sep = wxFileName::GetPathSeparator();
//   int m = a.Find(sep,true);
//   //  wxLogDebug("%d %c",m,sep);
//   if( m != wxNOT_FOUND ) {
//     a = a.Left(m+1);
//     //      a = d.GetPath();
//   }
//   wxFileName d(path);
//   //  if( !d.DirExists() ) {
//     //    d.Normalize();
//     //    d.RemoveLastDir();
//   //  }
//   wxLogDebug(a+" path was:"+path);

//   wxArrayString dirs;
//   wxDir dir(a);
//   if( dir.IsOpened() ) {

//     wxString dirname;
//     bool cont = dir.GetFirst(&dirname,"",wxDIR_DIRS);
//     while ( cont ) {

//       wxFileName d(path);
//       d.AppendDir(dirname);

//       //      wxLogDebug(d.GetPath()+" "+path);
//       if( /*d.DirExists() &&*/ d.GetPath().StartsWith(path) ) {
// 	dirs.Add(d.GetPath());
// 	wxLogDebug("Navigation add: "+d.GetPath());
//       }

//       cont = dir.GetNext(&dirname);
//     }
//   }
//   return dirs;
// }

// void MuniTextCtrl::OnChar(wxKeyEvent& event)
// {
//   if(event.GetKeyCode() == WXK_DOWN /*&& event.ShiftDown()*/ ) {
//     wxLogDebug("OnChar: DOWN");

//     wxSize size = GetSize();
//     wxPoint pos = GetPosition();
//     int y = pos.y + size.GetHeight();
//     int h = GetCharHeight();
//     wxLogDebug("%d %d %d %d",pos.x,pos.y,size.GetWidth(),size.GetHeight());

//     /*
//     wxWindow *win = (wxWindow*) event.GetEventObject();
//     wxSize sz = win->GetSize();
//     */

//     wxPopupTransientWindow *w = new wxPopupTransientWindow(GetParent(),wxBORDER_DEFAULT);
//     wxPoint p = ClientToScreen( wxPoint(0,0) );
//     w->Position(pos,size);
//     //    w->Position(pos,sz);
//     w->Popup();
//     //    w->Show();

//     /*
//     wxMenu popup;
//     popup.Append(wxID_ANY,"x");
//     popup.Append(wxID_ANY,"y");

//     PopupMenu(&popup);
//     */
//     return;
//   }

//   event.Skip();
// }



MuniNavigation::MuniNavigation(wxWindow *w, MuniConfig *c):
  wxPanel(w,wxID_ANY),config(c),phase(0)
{
  cmd = new wxTextCtrl(this,wxID_ANY,"*.fits",wxDefaultPosition,wxDefaultSize,
		       wxTE_PROCESS_ENTER|wxWANTS_CHARS);

  int h = GetCharHeight();
  wxBitmap bmp(h,h);
  moon = new wxStaticBitmap(this,wxID_ANY,bmp);

  MuniArtIcons ico(wxART_BUTTON,wxSize(24,24));
  wxBitmap oicon = ico.Icon(wxART_FOLDER_OPEN);
  bopen = new wxBitmapButton(this,wxID_OPEN,oicon,wxDefaultPosition,wxDefaultSize,
			     wxBORDER_NONE);
  bopen->SetToolTip("Change directory");

  wxBoxSizer *topsizer = new wxBoxSizer(wxHORIZONTAL);
  topsizer->Add(bopen,wxSizerFlags().Border());
  topsizer->Add(cmd,wxSizerFlags(1).Border());
  topsizer->Add(moon,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).Border().ReserveSpaceEvenIfHidden());
  moon->Show(false);
  SetSizer(topsizer);

  Bind(wxEVT_TEXT_ENTER,&MuniNavigation::OnEnter,this);
}

void MuniNavigation::OnEnter(wxCommandEvent& event)
{
  //  wxLogDebug("OnEnter: "+event.GetString());

  GetPath(event.GetString(),pwd,mask);
  NavigationEvent e(EVT_NAVIGATION);
  e.pwd = pwd;
  e.mask = mask;
  wxQueueEvent(this,e.Clone());
}


void MuniNavigation::OnUpdateMoon(wxUpdateUIEvent& event)
{
  const int width = 48;
  const int height = 48;
  const int h = GetCharHeight();

  int day = 56*phase;
  wxASSERT(0 <= day && day < 56);
  wxImage moon_56frames(config->moon_56frames);
  wxRect rect(1+day*width,1,width,height);
  wxImage sub = moon_56frames.GetSubImage(rect);
  sub.Rescale(h,h);
  moon->SetBitmap(wxBitmap(sub));
}

void MuniNavigation::NewMoon()
{
  moon->Show(true);
  cmd->Enable(false);
  bopen->Enable(false);
  Bind(wxEVT_UPDATE_UI,&MuniNavigation::OnUpdateMoon,this,moon->GetId());
}

void MuniNavigation::HideMoon()
{
  moon->Show(false);
  cmd->Enable(true);
  bopen->Enable(true);
  Unbind(wxEVT_UPDATE_UI,&MuniNavigation::OnUpdateMoon,this,moon->GetId());
}

void MuniNavigation::PhaseMoon(double f)
{
  phase = f;
}

void MuniNavigation::ChangeDir(const wxString& p, const wxString& m)
{
  pwd = p;
  mask = m;
  //  wxLogDebug("+*+*+*+*+*+* "+pwd+"/"+mask+"    `"+cmd->GetLabel()+"'");

  wxFileName filename(pwd,mask);
  filename.Normalize(NORMALIZE_FLAGS);

  cmd->Clear();
  cmd->AppendText(filename.GetFullPath());
}


void MuniNavigation::GetPath(const wxString& fullpath,
			   wxString& path, wxString& mask) const
{
  int l = fullpath.Find('*',true);
  int k = fullpath.Find('?',true);

  if( l != wxNOT_FOUND || k != wxNOT_FOUND ) {
    // mask
    wxChar sep = wxFileName::GetPathSeparator();
    path = fullpath.BeforeLast(sep);
    mask = fullpath.AfterLast(sep);
  }
  else {
    wxFileName filename(fullpath);
    filename.Normalize(NORMALIZE_FLAGS);
    if( filename.DirExists() ) {
      path = filename.GetFullPath();
      mask = "";
    }
    else {
      path = filename.GetPath();
      mask = filename.GetFullName();
    }
  }
}
