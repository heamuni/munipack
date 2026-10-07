/*

  Virtual observatory Sesame CLI client

  Copyright © 2019-2026 F.Hroch (hroch@physics.muni.cz)

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


#include "common.h"
#include "../config.h"
#include "vosesame.h"
#include <wx/wx.h>
#include <wx/txtstrm.h>
#include <wx/wfstream.h>
#include <wx/sstream.h>
#include <cstdio>

#if wxUSE_WEBREQUEST
#include <wx/webrequest.h>
#endif


#define UNISTRA    "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI?"
#define STRASBOURG "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame?"
#define USER_AGENT "Munipack/" PACKAGE_VERSION " (" wxVERSION_STRING  ")"

class Sesame: public wxAppConsole
{
  wxString object;

  int Print(wxInputStream *);
#if wxUSE_WEBREQUEST
  wxWebRequest request;
  void OnWebRequestState(wxWebRequestEvent& evt);
#endif

public:
  bool OnInit();
#if !wxUSE_WEBREQUEST
  int OnRun();
#endif
};
IMPLEMENT_APP_CONSOLE(Sesame)

bool Sesame::OnInit()
{
  bool verbose = false;

  // logging
  wxLog::DisableTimestamp();
  wxLog::SetLogLevel(wxLOG_Message);
  wxLog::SetVerbose(); // important to activate previous setup

  wxFFileInputStream istream(stdin);
  wxTextInputStream input(istream);

  while(istream.IsOk() && ! istream.Eof()) {

    wxString line = input.ReadLine();

    if( line.StartsWith("VERBOSE") ) {
      verbose = GetBool(line);
      if( verbose )
	wxLog::SetLogLevel(wxLOG_Debug);
      else
	wxLog::SetLogLevel(wxLOG_Message);
      wxLog::SetVerbose();
    }

    if( line.StartsWith("OBJECT") )
      object = GetString(line);
  }

  if( object.IsEmpty() ) {
    fprintf(stderr,"STOP 'An object name undefined.'\n");
    return false;
  }

#if wxUSE_WEBREQUEST
  wxString url(UNISTRA+object);

  request = wxWebSession::GetDefault().CreateRequest(this,url);
  request.SetHeader("User-Agent",USER_AGENT);
  request.SetStorage(wxWebRequest::Storage_Memory);
  if ( !request.IsOk() ) {
    wxLogError("Failed to create a Web request.");
    return false;
  }
  Bind(wxEVT_WEBREQUEST_STATE, &Sesame::OnWebRequestState, this);
  request.Start();
#endif

  return true;
}

#if wxUSE_WEBREQUEST

void Sesame::OnWebRequestState(wxWebRequestEvent& evt)
{
  wxWebResponse response = evt.GetResponse();
  int result;

  switch (evt.GetState()) {
  case wxWebRequest::State_Idle:
    wxLogDebug("WebRequest in idle.");
    break;
  case wxWebRequest::State_Active:
    wxLogDebug("WebRequest is active.");
    break;
  case wxWebRequest::State_Unauthorized:
    wxLogError("WebRequest is unauthorized.");
    ExitMainLoop();
    SetErrorExitCode(1);
    break;
  case wxWebRequest::State_Cancelled:
    wxLogDebug("WebRequest cancelled");
    fprintf(stderr,"STOP 'Cancelled'\n");
    ExitMainLoop();
    SetErrorExitCode(1);
    break;
  case wxWebRequest::State_Completed:
    if( response.IsOk() ) {
      wxLogDebug("Request: "+response.GetURL());
      wxLogDebug("HTTP response: %d ("+response.GetStatusText()+")",
		response.GetStatus());
      wxLogDebug("Content-Type: "+response.GetHeader("Content-Type"));
      wxLogDebug("Content-Length: "+response.GetHeader("Content-Length"));
      result = Print(response.GetStream());
    }
    else
      result = 1;
    fprintf(stderr,"STOP %d\n",result);
    ExitMainLoop();
    SetErrorExitCode(result);
    break;
  case wxWebRequest::State_Failed:
    if( response.IsOk() ) {
      wxLogDebug("Request: "+response.GetURL());
      if( response.GetStatus() > 0 )
	wxLogDebug("HTTP response: %d ("+response.GetStatusText()+")",
		   response.GetStatus());
    }
    wxLogInfo(evt.GetErrorDescription());
    fprintf(stderr,"STOP 'Failed'\n");
    ExitMainLoop();
    SetErrorExitCode(1);
    break;
  }
}

#else

/*
  This is based on wxHTTP class which does not support the https.
  I should not remove it as WebRequests (above) are only included
  by Debian Trixie 13 since '25, there are still too much older versions in use.
*/

int Sesame::OnRun()
{
  VOSesame sezam(STRASBOURG+object);
  if( sezam.Request() && Print(sezam.GetStream()) == 0 ) {
    fprintf(stderr,"STOP 0\n");
    return 0;
  }
  else {
    fprintf(stderr,"STOP '%s'\n",sezam.GetErrMsg());
    return 1;
  }
}

#endif

int Sesame::Print(wxInputStream *istream)
{
  double ra,dec;
  bool defined = false;

  wxTextInputStream text(*istream);
  while( ! istream->Eof() && istream->IsOk()) {
    wxString line = text.ReadLine();
    wxLogDebug("%s",line.c_str());
    // Looks stupid, but it does not thow "Assert failure"" failed in FromWChar():
    // trying to encode undefined Unicode character

    if( line.StartsWith("%J") && (line.Find("%J.E") == wxNOT_FOUND) ) {

      wxStringInputStream s(line.Mid(3));
      wxTextInputStream text(s);
      wxString a,b;
      text >> a >> b;
      defined = a.ToDouble(&ra) && b.ToDouble(&dec);
    }

    if( line.StartsWith("#!") ) {
      wxLogError(line.Mid(3));
      return 1;
    }
  }

  if( defined )
    fprintf(stdout,"%.7f %.7f \n",ra,dec);
  else
    fprintf(stderr,"Required object not found.\n");

  return defined ? 0 : 1;
}
