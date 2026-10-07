/*

  Common functions

  Copyright © 2018-2026 F.Hroch (hroch@physics.muni.cz)

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

#include "common.h"
#include <wx/wx.h>
#include <wx/regex.h>
#include <wx/filename.h>

// ------------------------------------------------------------------
// Common functions

// it would be nice to create a fortran I/O interoperability functions
// (library)?


wxString GetString(const wxString& line)
{
  wxRegEx re(".+ = '(.*)'");
  // double quotes in string:  "'(''|[^'])*'" - perhaps unfunctional
  //  wxRegEx re(".+ = '(''|[^']*)'");
  wxASSERT(re.IsValid());
  if( re.Matches(line) )
    return re.GetMatch(line,1);
  else
    return "";
}

double GetDouble(const wxString& line)
{
  double x;
  wxString a = line.AfterFirst('=');

  if( a.ToCDouble(&x) )
    return x;

  wxLogFatalError("Failed to read the number: "+a);
  return 666; // formally
}

long GetLong(const wxString& line)
{
  long l;
  wxString a = line.AfterFirst('=');

  if( a.ToCLong(&l) )
    return l;

  wxLogFatalError("Failed to read the number: "+a);
  return 0; // formally
}

bool GetBool(const wxString& line)
{
  wxRegEx re(".+ = (.*)");
  wxASSERT(re.IsValid());
  if( re.Matches(line) ) {
    wxString a = re.GetMatch(line,1);
    a.Upper();
    if( a.StartsWith("T") )
      return true;
    else if( a.StartsWith("F") )
      return false;
  }
  wxLogFatalError("Failed to parse boolean (logical) expression: `"+line+"'.");
  return false; // formally
}

wxString GetFileType(const wxString& output)
{
  // determine file type of output by suffix
  wxString type;
  if( ! output.IsEmpty() ) {
    wxFileName fn(output);
    if( fn.IsOk() )
      type = fn.GetExt().Upper();
  }
  return type;
}
