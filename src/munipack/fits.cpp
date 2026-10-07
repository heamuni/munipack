/*

  FITS manipulations

  Copyright © 2011 - 2026 F.Hroch (hroch@physics.muni.cz)

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

  Note.

  Multiple use of AddOption is not supported. Therefore the

     fits -K NAXIS1 -K NAXIS2

  does not works and tokenizers(s) must be used on boths tails of pipes.

*/

#include "munipack.h"
#include <wx/wx.h>
#include <wx/tokenzr.h>


using namespace std;


MuniFits::MuniFits(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"fits")
{
  wxString mode = "STRUCTURE";
  wxString a;

  if( cmd.Found("K",&a) || cmd.Found("keys",&a) ) {

    mode = "HEADER";

    wxStringTokenizer tokenizer(a, ",");
    while ( tokenizer.HasMoreTokens() )
      Write("KEYWORD = '" + tokenizer.GetNextToken() + "'");
  }

  else if( cmd.Found("remove-keys",&a) ) {

    mode = "HEADER";

    wxStringTokenizer tokenizer(a, ",");
    while ( tokenizer.HasMoreTokens() )
      Write("REMKEY = '" + tokenizer.GetNextToken() + "'");
  }

  else if( cmd.Found("update") ) {

    mode = "HEADER";

    Write("UPDATE = T");

    if( cmd.Found("key",&a) )
      Write("KEY = '"+a+"'");

    if( cmd.Found("val",&a) ) {
      wxString b = a.Find("'") == wxNOT_FOUND ? "'" : "\"";
      Write("VALUE = "+b+a+b);
    }

    if( cmd.Found("com",&a) ) {
      wxString b = a.Find("'") == wxNOT_FOUND ? "'" : "\"";
      Write("COMMENT = "+b+a+b);
    }

    if( cmd.Found("templ",&a) )
      Write("TEMPL = '"+a+"'");
  }

  else if( cmd.Found("lh") || cmd.Found("header") ) {

    mode = "HEADER";

  }

  else if( cmd.Found("lt") || cmd.Found("table") ) {

    mode = "TABLE";

  }

  else if( cmd.Found("li") || cmd.Found("image") ) {

    mode = "IMAGE";

  }

  else if( cmd.Found("dump") ) {

    mode = "DUMP";

  }
  else if( cmd.Found("restore") ) {

    mode = "RESTORE";

  }
  else if( cmd.Found("cat") ) {

    mode = "CAT";

  }
  else if( cmd.Found("remove-extensions",&a) ) {

    mode = "EXTENSION";

    wxStringTokenizer tokenizer(a, ",");
    while ( tokenizer.HasMoreTokens() )
      Write("REMEXT = '" + tokenizer.GetNextToken() + "'");
  }

  Write("MODE = '" + mode + "'");

  if( cmd.Found("shell") )
    Write("KEYLIST = 'SHELL'");

  if( cmd.Found("value") )
    Write("KEYLIST = 'VALUE'");

  CommonOutputMultiple(cmd);

  if( ! cmd.Found("O") && (mode == "RESTORE" || mode == "DUMP") ) {
    // In this modes, .lst are by default

    SetAdvanced(true);
    if( mode == "RESTORE" ) {
      SetPattern("(.*).lst");
      SetMask("\\1.fits");
    }
    else if( mode == "DUMP" ) {
      SetPattern("(.*).fits");
      SetMask("\\1.lst");
    }
  }

  WriteFilesMultiple(cmd.GetFiles());
}
