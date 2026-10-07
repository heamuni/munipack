/*

  xmunipack - FITS extension side-list

  Copyright © 2019-22 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/dataview.h>
#include <wx/wx.h>
#include <wx/graphics.h>
#include <vector>

#define LABEL "Extension"

using namespace std;

MuniExtensionList::MuniExtensionList(wxWindow *w, const MuniConfig *c):
  wxDataViewListCtrl(w,wxID_ANY), config(c)
{
  // indexes
  idxtype[HDU_UNKNOWN] = 0;
  idxtype[HDU_HEAD] = 1;
  idxtype[HDU_IMAGE] = 2;
  idxtype[HDU_TABLE] = 3;

  // icon size
  int height = GetCharHeight();
  int size = int(height + 1);

  // icons
  icons.push_back(DrawIcon(size,HDU_HEAD)); // unknown
  icons.push_back(DrawIcon(size,HDU_HEAD));
  icons.push_back(DrawIcon(size,HDU_IMAGE));
  icons.push_back(DrawIcon(size,HDU_TABLE));

  // Moon
  phases = MoonPhases(height);

  Init();

  Bind(wxEVT_DATAVIEW_SELECTION_CHANGED,&MuniExtensionList::OnExtChanged,this);
}

void MuniExtensionList::Clear()
{
  DeleteAllItems();
  ClearColumns();
}

void MuniExtensionList::Init()
{
  wxASSERT(GetItemCount() == 0 && GetColumnCount() == 0);

  int chwidth = GetCharWidth();

  timer.Start();
  AppendIconTextColumn(LABEL,wxDATAVIEW_CELL_INERT,15*chwidth);
  AppendIconTextColumn(L"⇩",wxDATAVIEW_CELL_INERT,2*chwidth);
}

void MuniExtensionList::Finish()
{
  wxASSERT(GetItemCount() > 0 && GetColumnCount() == 2);

  wxDataViewColumn *col = GetColumn(1);
  DeleteColumn(col);
  SelectRow(0);
}


void MuniExtensionList::Append(const wxString& extname, int hdutype)
{
  wxASSERT(GetColumnCount() == 2);
  wxLogDebug("MuniExtensionList::Append: "+extname+": %d %d",hdutype,int(GetItemCount()));

  const size_t maxlen = 15;
  wxString label(extname);
  if( label.Len() > maxlen )
    label = label.Truncate(maxlen) + "...";
  wxVariant name;
  name << wxDataViewIconText(label, icons[hdutype]);

  wxVariant moonicon;
  moonicon << wxDataViewIconText("", phases[0]);

  wxVector<wxVariant> cols;
  cols.push_back(name);
  cols.push_back(moonicon);
  AppendItem(cols);
  timer.Start();

}

void MuniExtensionList::FinishProgress(int chdu)
{
  UpdateProgress(chdu,1.01);
}

void MuniExtensionList::UpdateProgress(int chdu, double f)
{
  int n = int(28*f);

  if( GetItemCount() == chdu ) {
    if( timer.Time() > 250 ) {
      wxVariant moonicon;
      moonicon << wxDataViewIconText("",phases[n]);
      SetValue(moonicon,chdu-1,1);
    }
  }

  /*
    There is a wx/GTK trouble.
    The warning:

(xmunipack:31084): Gtk-CRITICAL **: 12:24:47.105: gtk_widget_queue_draw_area: assertion 'width >= 0' failed

is issued shortly after start when something is probably uninitialised.
I found the only solution: Moon icon should be draw after some delay.

  */
}


void MuniExtensionList::ChangeSelection(int n)
{
  SelectRow(n);
}

void MuniExtensionList::OnExtChanged(wxDataViewEvent& event)
{
  if( ! IsOk() ) return;

  int n = GetSelectedRow();
  //  wxLogDebug("MuniExtensionList::OnExtChanged %d",n);

  if( n >= 0 ) {
    wxCommandEvent e(wxEVT_COMMAND_MENU_SELECTED,ID_EXTENSION_SELECT);
    e.SetInt(n);
    wxQueueEvent(GetParent(),e.Clone());
  }
}

bool MuniExtensionList::IsOk() const
{
  return GetItemCount() > 0 && GetColumnCount() == 1;
}


wxIcon MuniExtensionList::DrawIcon(int size, int hdutype) const
{
  wxBitmap bmp(size,size);
  wxMemoryDC mdc(bmp);

  wxGraphicsContext *gc = wxGraphicsContext::Create(mdc);
  if( gc ) {

    wxString ch;
    if( hdutype == HDU_IMAGE )
      ch = config->image_symbol;
    else if( hdutype == HDU_TABLE )
      ch = config->table_symbol;
    else if( hdutype == HDU_HEAD )
      ch = config->head_symbol;

    gc->SetBrush(wxBrush(wxSystemSettings::GetColour(wxSYS_COLOUR_LISTBOX)));
    gc->DrawRectangle(0,0,size,size);

    double w,h,u,v;
    gc->SetFont(*wxNORMAL_FONT,*wxBLACK);
    gc->GetTextExtent(ch,&w,&h,&u,&v);
    int x = wxMax((size - w) / 2,0);
    int y = wxMax((size - h) / 2,0);
    gc->DrawText(ch,x,y);

    delete gc;
  }
  mdc.SelectObject(wxNullBitmap);

  wxIcon icon;
  icon.CopyFromBitmap(bmp);
  return icon;
}


vector<wxIcon> MuniExtensionList::MoonPhases(int size) const
{
  const int width = 48;
  const int height = 48;

  wxImage moon_56frames(config->moon_56frames);

  vector<wxIcon> icons;
  for(int i = 0; i < 56; i++)
  {
     wxRect rect(1+i*width,1,width,height);
     wxImage sub = moon_56frames.GetSubImage(rect);
     wxBitmap bmp(sub.Rescale(size,size));
     wxIcon icon;
     icon.CopyFromBitmap(bmp);
     icons.push_back(icon);
  }
  return icons;

}
