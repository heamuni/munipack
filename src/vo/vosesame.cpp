/*

  Virtual observatory Sesame client (to resolve names of astronomical objects)

  Copyright © 2019-20 F.Hroch (hroch@physics.muni.cz)

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
     http://vizier.u-strasbg.fr/vizier/doc/sesame.htx

*/

#include "../config.h"
#include "vosesame.h"
#include <wx/wx.h>
#include <wx/url.h>
#include <wx/tokenzr.h>
#include <wx/txtstrm.h>
#include <wx/protocol/http.h>

/*
#define STRASBOURG "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame?"
#define HARVARD    "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame?"
*/

VOSesame::VOSesame(const wxString& o):
  wxHTTP(), object(o), ra(-999.9), dec(-99.9), url(o), //url(STRASBOURG+object),
  istream(0)
{
  if( ! url.IsOk() ) {
    if( url.GetError() == wxURL_SNTXERR )
      errmsg = "Syntax error in the URL string.";
    else if( url.GetError() == wxURL_NOPROTO )
      errmsg = "Found no protocol which can get this URL.";
    else if( url.GetError() == wxURL_NOHOST )
      errmsg = "A host name is required for this protocol.";
    else if( url.GetError() == wxURL_NOPATH )
      errmsg = "A path is required for this protocol.";
  }
}

VOSesame::~VOSesame()
{
  delete istream;
}

const char* VOSesame::GetErrMsg() const
{
  return errmsg.c_str();
}


bool VOSesame::Connect(const wxString& server)
{
  SetFlags(wxSOCKET_WAITALL);
  SetHeader("User-Agent",wxString("Munipack/")+wxString(PACKAGE_VERSION));
  wxLogInfo("Connecting to VOSesame http://"+server+" ...");

  return wxHTTP::Connect(server);
}

bool VOSesame::Resolve()
{
  if( errmsg != "" ) return false;

  if( ! Connect(url.GetServer()) ) {
    errmsg = "Failed to connect server `http://"+url.GetServer()+"'.";
    return false;
  }

  if( ! Get(url.GetPath() + wxString("?") + url.GetQuery()) )
    return false;

  return true;
}


bool VOSesame::Get(const wxString& query)
{
  wxLogInfo("Sending query `"+query+"'...");
  istream = GetInputStream(query);
  wxASSERT(istream);

  wxLogInfo("HTTP response: %d",GetResponse());
  wxLogInfo("HTTP protocol state: %d",GetError());
  wxLogInfo("Content-Type (code):",GetContentType());
  wxLogInfo("Content-Type (head):",GetHeader("Content-Type"));
  wxLogInfo("Content-Length:",GetHeader("Content-Length"));

  if( GetResponse() != 200 ) {
    errmsg.Printf("HTTP response: %d (failed).",GetResponse());
    return false;
  }

  if( GetError() != wxPROTO_NOERR ) {
    errmsg.Printf("Network error %d",(int) GetError());
    return false;
  }

  wxLogInfo("Downloading and parsing response ...");

  wxTextInputStream text(*istream,"\r\n");
  while( ! istream->Eof() && istream->IsOk()) {
    wxString line = text.ReadLine();
    if( istream->LastRead() > 0 ) {
      wxLogInfo(line);

      if( line.StartsWith("%J ") ) {
	wxStringTokenizer tk(line);
	for(int n = 0; tk.HasMoreTokens(); n++) {
	  wxString item = tk.GetNextToken();
	  double q;
	  if( n == 1 && item.ToDouble(&q) )
	    ra = q;
	  if( n == 2 && item.ToDouble(&q) )
	    dec = q;
	}
      }

      if( line.StartsWith("#!") ) {
	if( line.Find("Nothing found") )
	  errmsg = "Required object not found.";
	else
	  errmsg = line.Mid(3);
	return false;
      }
    }
  }

  return ra > -999 && dec > -99;
}


bool VOSesame::Request()
{
  wxString server = url.GetServer();
  wxString query = url.GetPath() + wxString("?") + url.GetQuery();

  SetFlags(wxSOCKET_WAITALL);
  SetHeader("User-Agent",wxString("Munipack/")+wxString(PACKAGE_VERSION));
  wxLogInfo("Connecting to VOSesame http://"+server+" ...");

  if( wxHTTP::Connect(server) ) {
    wxLogInfo("Sending query `"+query+"'...");
    istream = GetInputStream(query);
    wxASSERT(istream);

    wxLogInfo("HTTP response: %d",GetResponse());
    wxLogInfo("HTTP protocol state: %d",GetError());
    wxLogInfo("Content-Type (code):",GetContentType());
    wxLogInfo("Content-Type (head):",GetHeader("Content-Type"));
    wxLogInfo("Content-Length:",GetHeader("Content-Length"));

    if( GetResponse() != 200 ) {
      errmsg.Printf("HTTP response: %d (failed).",GetResponse());
      return false;
    }

    if( GetError() != wxPROTO_NOERR ) {
      errmsg.Printf("Network error %d",(int) GetError());
      return false;
    }

    wxLogInfo("Downloading and parsing response ...");
    return true;
  }

  return false;
}
