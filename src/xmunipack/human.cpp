/*

  Prints numbers in human format

  Copyright © 2019 F.Hroch (hroch@physics.muni.cz)

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


  https://en.wikipedia.org/wiki/Metric_prefix

*/

#include "fits.h"
#include <wx/wx.h>
#include <cfloat>


wxString HumanFormat(double x)
{
  const wchar_t *prefin[] = { L"",L"m",L"μ",L"n",L"p",L"f",L"a",L"z",L"y" };
  const wchar_t *prefix[] = { L"",L"k",L"M",L"G",L"T",L"P",L"E",L"Z",L"Y" };
  const wchar_t *form[] = { L"%.3f", L"%.2f", L"%.1f" };

  wxString text;
  double w = fabs(x);

  if( w < DBL_EPSILON )
    return "0";

  double sgn = x / w;
  double d = log10(w);
  if( -24 <= d && d < 27 ) {
    int thorder = d / 3;
    if( d < 0 ) thorder = (d - 1) / 3;
    double q = w / pow(10.0,3*thorder);
    int order = log10(q);
    int o = d >= 0 ? 1 : -1;
    //    wxLogDebug(L"%e %e %d %f %d %f",x,w,thorder,q,order,d);
    wxASSERT(0 <= o*order && o*order < 3);
    wxString fmt;
    fmt.Printf("%s%s",form[o*order],"%s");
    if( d >= 0 ) {
      wxASSERT(0 <= thorder && thorder < 9);
      text.Printf(fmt,sgn*q,prefix[thorder]);
    }
    else {
      wxASSERT(0 <= -thorder && -thorder < 9);
      text.Printf(fmt,sgn*q,prefin[-thorder]);
    }
  }
  else
    text.Printf("%e",x);
  //    printf("%e %f %d\n",x,d,thorder);
  return text;
}


/*

 g++ `wx-config --cxxflags` human.cpp testhuman.cpp `wx-config --libs`

#include "fits.h"
#include <wx/wx.h>
#include <cfloat>

int main()
{
  wxLogDebug(HumanFormat(0.0));
  wxLogDebug(HumanFormat(1.0));
  wxLogDebug(HumanFormat(0.5));
  wxLogDebug(HumanFormat(5));

  wxLogDebug(HumanFormat(666.0));
  wxLogDebug(HumanFormat(2e5));
  wxLogDebug(HumanFormat(1.23456789e6));
  wxLogDebug(HumanFormat(1e12));

  wxLogDebug(HumanFormat(-1.0));
  wxLogDebug(HumanFormat(-666.0));
  wxLogDebug(HumanFormat(-2e5));
  wxLogDebug(HumanFormat(-1.23456789e6));
  wxLogDebug(HumanFormat(-1e12));

  wxLogDebug(HumanFormat(0.005));
  wxLogDebug(HumanFormat(2e-5));
  wxLogDebug(HumanFormat(1.23456789e-6));
  wxLogDebug(HumanFormat(1e-12));

  wxLogDebug(HumanFormat(-0.005));
  wxLogDebug(HumanFormat(-2e-5));
  wxLogDebug(HumanFormat(-1.23456789e-6));
  wxLogDebug(HumanFormat(-1e-12));

  wxLogDebug(HumanFormat(1.01e9));
  wxLogDebug(HumanFormat(0.99e9));

  wxLogDebug(HumanFormat(1.01e-9));
  wxLogDebug(HumanFormat(0.99e-9));


  wxLogDebug(HumanFormat(5000));
  wxLogDebug(HumanFormat(0.005));


}

*/
