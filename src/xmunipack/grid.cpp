/*

  xmunipack - table


  Copyright © 1997-2011 F.Hroch (hroch@physics.muni.cz)

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


  It will need a lot of work...

*/



#include "xmunipack.h"
#include <wx/wx.h>
#include <algorithm>

using namespace std;


MuniGrid::MuniGrid(wxWindow *w, MuniConfig *c): wxGrid(w,wxID_ANY) {}

bool MuniGrid::SetHdu(const FitsHdu& h)
{
  table = FitsTable(h);

  wxASSERT(table.IsOk());
  if( ! table.IsOk() ) return false;

  CreateGrid(table.Nrows(),table.Ncols());

  BeginBatch();

  for(int i = 0; i < table.Ncols(); i++) {
    wxString key;
    key.Printf("TTYPE%d",i+1);
    wxString label = table.GetKey(key);
    SetColLabelValue(i,label);
  }

  // fill table by blank values
  /*
  for(int i = 0; i < table.Ncols(); i++)
    for(int j = 0; j < table.Nrows(); j++)
      SetCellValue(j,i,wxEmptyString);
  */
  EndBatch();

  Bind(wxEVT_IDLE,&MuniGrid::OnIdle,this);

  rows_filled = 0;

  return true;
}


// implement via threading ???
void MuniGrid::OnIdle(wxIdleEvent& event)
{
  const int npart = 100;

  int imin = rows_filled;
  int imax = imin + min(npart,int(table.Nrows()-imin));

  BeginBatch();
  //  grid->AppendRows(imax - imin);
  for(int i = 0; i < table.Ncols(); i++)
    for(int j = imin; j < imax; j++)
      SetCellValue(j,i,table.Cell_str(j,i));
  EndBatch();

  rows_filled = imax;

  if( table.Nrows() > imax )
    event.RequestMore();
  else
    Unbind(wxEVT_IDLE,&MuniGrid::OnIdle,this);
}
