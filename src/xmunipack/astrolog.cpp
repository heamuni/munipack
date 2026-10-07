/*

  xmunipack - astrological symbols

  Copyright © 2015 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/wx.h>
#include <wx/datetime.h>

MuniAstrolog::MuniAstrolog()
{
  datetime.SetToCurrent();
}

MuniAstrolog::MuniAstrolog(const wxDateTime& dt)
{
  datetime.Set(dt.GetTm());
}

wxString MuniAstrolog::GetSign() const
{
  int y = datetime.GetYear();
  wxDateTime d(datetime);

  if( d.IsBetween(wxDateTime(20,wxDateTime::Mar,y),wxDateTime(20,wxDateTime::Apr,y)) )
    return L"♈";
  else if(d.IsBetween(wxDateTime(20,wxDateTime::Apr,y),wxDateTime(21,wxDateTime::May,y)))
    return L"♉";
  else if(d.IsBetween(wxDateTime(21,wxDateTime::May,y),wxDateTime(21,wxDateTime::Jun,y)))
    return L"♊";
  else if(d.IsBetween(wxDateTime(21,wxDateTime::Jun,y),wxDateTime(23,wxDateTime::Jul,y)))
    return L"♋";
  else if(d.IsBetween(wxDateTime(23,wxDateTime::Jul,y),wxDateTime(23,wxDateTime::Aug,y)))
    return L"♌";
  else if(d.IsBetween(wxDateTime(23,wxDateTime::Aug,y),wxDateTime(23,wxDateTime::Sep,y)))
    return L"♍";
  else if(d.IsBetween(wxDateTime(23,wxDateTime::Sep,y),wxDateTime(23,wxDateTime::Oct,y)))
    return L"♎";
  else if(d.IsBetween(wxDateTime(23,wxDateTime::Oct,y),wxDateTime(22,wxDateTime::Nov,y)))
    return L"♏";
  else if(d.IsBetween(wxDateTime(23,wxDateTime::Nov,y),wxDateTime(21,wxDateTime::Dec,y)))
    return L"♐";
  else if(d.IsBetween(wxDateTime(21,wxDateTime::Jan,y),wxDateTime(19,wxDateTime::Feb,y)))
    return L"♒";
  else if(d.IsBetween(wxDateTime(20,wxDateTime::Feb,y),wxDateTime(20,wxDateTime::Mar,y)))
    return L"♓";
  else  
    return L"♑";
}
