/*

  Virtual Observatory Client

  Copyright © 2018-9, 2025 F.Hroch (hroch@physics.muni.cz)

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

/*

  This file contais unfinished implementation of http client
  via webrequiest interface

  https://docs.wxwidgets.org/3.2/classwx_web_request.html

  Unfortunatelly, is requires cURL which is unavailable
  for wxWidgets in my current system.

 */


#include "../config.h"
#include "voclient.h"
#include <wx/wx.h>
#include <wx/protocol/http.h>
#include <wx/wfstream.h>
#include <wx/ffile.h>
#include <wx/txtstrm.h>

VOclient::VOclient(): wxHTTP() {}

VOclient::VOclient(const wxString& server): wxHTTP()
{
  Connect(server);
}

bool VOclient::Connect(const wxString& server)
{
  SetFlags(wxSOCKET_WAITALL);

  SetHeader("User-Agent",wxString("Munipack/")+wxString(PACKAGE_VERSION));
  wxLogInfo("Connecting VO server http://"+server);

  return wxHTTP::Connect(server);
}

bool VOclient::Get(const wxString& path, const wxString& filename)
{
  wxLogInfo("Sending query "+path);
  wxInputStream *istream = GetInputStream(path);

  wxLogInfo("HTTP response: %d",GetResponse());
  wxLogInfo("HTTP protocol state: %d",GetError());
  wxLogInfo("Content-Type (code):",GetContentType());
  wxLogInfo("Content-Type (head):",GetHeader("Content-Type"));
  wxLogInfo("Content-Length:",GetHeader("Content-Length"));

  if( GetResponse() != 200 )
    wxLogError("HTTP response: %d (failed).",GetResponse());

  if( istream && GetResponse() == 200 && GetError() == wxPROTO_NOERR ) {

    wxLogInfo("Downloading data ("+filename+") ...");
    wxTextInputStream text(*istream,"\r\n");

    wxFFile file(filename,"w");
    wxFFileOutputStream ostream(file);
    wxTextOutputStream cout(ostream);

    while( ! istream->Eof() && istream->IsOk() ) {
      wxString line = text.ReadLine();
      if( istream->LastRead() > 0 )
	cout.WriteString(line+"\n");
    }
  }
  delete istream;

  if( GetError() != wxPROTO_NOERR )
    wxLogError("Network error %d",(int) GetError());

  return GetResponse() == 200 && GetError() == wxPROTO_NOERR;
}



/*
#include <wx/file.h>
#include <wx/webrequest.h>


bool VOclient::Query(wxEvtHandler *handler, const wxString& server,
		     const wxString& query, const wxString& filename): handler(0)
{
  wxLogInfo("Connecting "+server);
  wxLogInfo("Sending query "+path);

  wxString a = server + query;
  wxWebRequest request = wxWebSession::GetDefault().CreateRequest(this,a);
  request.SetStorate(Storage_None);

  wxWebResponse r = request.GetResponse();
  wxLogDebug(r.AsString());

  if ( !request.IsOk() ) {
    wxLogError(request.GetResponse().AsString());
    return false;
  }

  Bind(wxEVT_WEBREQUEST_STATE,&VOclient::OnWebRequestState,this);
  Bind(wxEVT_WEBREQUEST_DATA,&VOclient::OnWebRequestData,this);

  tmpfile = wxFile(filename,"w");
  if( ! tmpfile.IsOpened() ) {
    wxLogError("VOClient failed to create temporary file.");
    return false;
  }

  request.Start();
}

void VOclient::OnWebRequestState(wxWebRequestEvent& evt)
{
  switch (evt.GetState())
    {

    case wxWebRequest::State_Completed:
      {
	tmfile.Close();
	// finish this class
	break;
      }
    case wxWebRequest::State_Failed:
      wxLogError("Could not load logo: %s", evt.GetErrorDescription());
      break;
    }

}


void VOclient::OnWebRequestData(wxWebRequestEvent& evt)
{
  tmpfile.Write(evt.GetDataBufer(),evt.GetDataSize());
  wxLogDebug("VOclient::OnWebRequestData received: %d bytes",(int)evt.GetDataSize());
}
*/
