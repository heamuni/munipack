/*

  VOTable convertor

  Copyright © 2010 - 2013, 2017-2026 F.Hroch (hroch@physics.muni.cz)

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


  Reference:

  http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html

*/

#include "common.h"
#include "votable.h"
#include <wx/wx.h>
#include <wx/app.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>
#include <cstdio>

using namespace std;

class VOTab: public wxAppConsole
{
public:
  bool OnInit();
  int OnRun();
private:
  wxString vooper();
};

IMPLEMENT_APP_CONSOLE(VOTab)

bool VOTab::OnInit()
{
  wxLog::DisableTimestamp();
  return true;
}

int VOTab::OnRun()
{
  wxString msg = vooper();

  // Looks strangle? See: cone.cpp.
  if( msg == "" )
    fprintf(stderr,"STOP 0\n");
  else
    fprintf(stderr,"STOP '%s'\n",static_cast<const char *>(msg.c_str()));
  return msg != "";
}


wxString VOTab::vooper()
{
  wxString output, type, file;
  wxString pt, mk, ak, dk;
  double pa, pd, ps, ml;
  long pw = -1, ph = -1;

  // limit logging
  wxLog::SetLogLevel(wxLOG_Message);
  wxLog::SetVerbose(); // important to activate previous setup

  wxFFileInputStream pstream(stdin);
  wxTextInputStream input(pstream);

  while(pstream.IsOk() && ! pstream.Eof()) {

    wxString line = input.ReadLine();

    if( line.StartsWith("VERBOSE") ) {
      if( GetBool(line) )
	wxLog::SetLogLevel(wxLOG_Debug);
      else
	wxLog::SetLogLevel(wxLOG_Message);
      wxLog::SetVerbose(); // important to activate previous setup
    }

    if( line.StartsWith("TYPE") )
      type = GetString(line);

    if( line.StartsWith("FILE") )
      file = GetString(line);

    if( line.StartsWith("PROJ TYPE") )
      pt = GetString(line);

    if( line.StartsWith("PROJ ALPHA") )
      pa = GetDouble(line);

    if( line.StartsWith("PROJ DELTA") )
      pd = GetDouble(line);

    if( line.StartsWith("PROJ SCALE") )
      ps = GetDouble(line);

    if( line.StartsWith("PROJ WIDTH") )
      pw = GetLong(line);

    if( line.StartsWith("PROJ HEIGHT") )
      ph = GetLong(line);

    if( line.StartsWith("MAG LIMIT") )
      ml = GetDouble(line);

    if( line.StartsWith("COL_RA") )
      ak = GetString(line);

    if( line.StartsWith("COL_DEC") )
      dk = GetString(line);

    if( line.StartsWith("COL_MAG") )
      mk = GetString(line);

    if( line.StartsWith("OUTPUT") )
      output = GetString(line);

  }


  if( type == "FITS" && ! output.IsEmpty() ) {
    FITStable vt(file);
    return vt.SaveFITS(output,true) ? "" : "FITS save failed.";
  }

  wxFFile ffile;
  if( output.IsEmpty() )
    ffile.Attach(stdout);
  else
    ffile.Open(output.c_str(),"w");

  if( ! ffile.IsOpened() )
    return "Failed to create or open `"+output+"'.";

  wxFFileOutputStream ostream(ffile);
  if( ! ostream.IsOk() )
    return "Failed to initialise an output stream.";

  if( type.IsEmpty() || type == "XML" ) {
    VOTable vt(file);
    return vt.SaveVO(ostream) ? "" : "XML save failed.";
  }

  else if( type == "CSV" ) {
    CSVtable vt(file);
    return vt.SaveCSV(ostream) ? "" : "CSV save failed.";
  }

  else if( type == "TXT" ) {
    TXTable vt(file);
    return vt.SaveTXT(ostream) ? "" : "Text save failed.";
  }
  else if( type == "SVG" ) {
    SVGcanvas vt(file);

    vt.SetProjection(pt);
    vt.SetProjectionCenter(pa,pd);
    vt.SetScale(ps);
    vt.SetMaglim(ml);
    vt.SetMagkey(mk);
    vt.SetAlphakey(ak);
    vt.SetDeltakey(dk);
    if( pw > 0 && ph > 0 )
      vt.SetCanvasSize(pw,ph);

    return vt.SaveSVG(ostream) ? "" : "SVG save failed.";
  }
  else
    return "Type `"+type+"' unrecognized.";

  return "Votable failed.";
}
