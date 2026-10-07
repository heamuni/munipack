
/*
   xmunipack - time conversions

   Copyright © 1997-2011, 2019-2022 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/tokenzr.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


FitsTime::FitsTime(const wxString& text)
{
  year = 0; month = 0; day = 0;
  hour = 0; minute = 0; second = 0; milisecond = 0;
  SetDateTime(text);
}

void FitsTime::SetDateTime(const wxString& text)
{
  if( text.IsEmpty() )
    return;

  wxString a;
  wxStringTokenizer tb(text,"T");
  if( tb.HasMoreTokens() )
    date = tb.GetNextToken().Trim();
  if( tb.HasMoreTokens() )
    time = tb.GetNextToken();

  SetDate(date);
  SetTime(time);
}

void FitsTime::SetDate(const wxString& xdate)
{
  wxString a;
  wxStringTokenizer tb(xdate,"-");
  if( tb.HasMoreTokens() ) {
    a = tb.GetNextToken().Trim();
    if( ! a.ToLong(&year) )
      year = 0;
  }
  if( tb.HasMoreTokens() ) {
    a = tb.GetNextToken().Trim();
    if( ! a.ToLong(&month) )
      month = 0;
  }
  if( tb.HasMoreTokens() ) {
    a = tb.GetNextToken().Trim();
    if( ! a.ToLong(&day) )
      day = 0;
  }
}

void FitsTime::SetTime(const wxString& xtime)
{
  double sec;
  wxString a;
  wxStringTokenizer tb(xtime,":");
  if( tb.HasMoreTokens() ) {
    a = tb.GetNextToken().Trim();
    if( ! a.ToLong(&hour) )
      hour = 0;
  }
  if( tb.HasMoreTokens() ) {
    a = tb.GetNextToken().Trim();
    if( ! a.ToLong(&minute) )
      minute = 0;
  }
  if( tb.HasMoreTokens() ) {
    a = tb.GetNextToken().Trim();
    if( ! a.ToDouble(&sec) )
      sec = 0.0;
  }

  second = long(sec);
  milisecond = long(1000.0*(sec - second));
}

double FitsTime::GetJd() const
{
  double y,m,d,j;

  //  wxLogDebug("%ld %ld %ld %ld %ld %ld.%ld",year,month,day,hour,minute,second,milisecond);

  d = day + (hour + (minute + (second + milisecond/1000.0)/60.0)/60.0)/24.0;

  y = year;
  if( y < 0 ) y = y + 1;

  if( month > 2 )
    m = month + 1;
  else {
    y = y - 1.0;
    m = month + 13;
  }

  j = int(365.25*y) + int(30.6001*m) + d + 1720994.5;
  if( d + 31*(m + 12*y) >=  15 + 31*(10 + 12*1582) ) {
    double a = int(y/100.0);
    j = j + 2.0 - a + int(a/4.0);
  }

  return j;
}
