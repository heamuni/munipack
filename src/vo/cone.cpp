/*

  Virtual Observatory capable cone search

  Copyright © 2010-2026 F.Hroch (hroch@physics.muni.cz)

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

  Obsolete parts: voclient is obsolete and will be removed as OnRun().
  A connection can be only https based.

*/

#include "common.h"
#include "../config.h"
#include "votable.h"
#include "vocatconf.h"
#include "voclient.h"
#include <wx/wx.h>
#include <wx/app.h>
#include <wx/url.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>
#include <wx/filename.h>
#include <wx/filefn.h>
#include <fitsio.h>
#include <list>
#include <cstdio>
#include <cmath>

#if wxUSE_WEBREQUEST
#include <wx/webrequest.h>
#endif

#define USER_AGENT "Munipack/" PACKAGE_VERSION " (" wxVERSION_STRING  ")"

class Cone: public wxAppConsole
{
  bool PatchJohn(const wxString&) const;
  void PipePrint(const wxString&) const;
  void Stop(const wxString&) const;
  int Save() const;
  wxString url, output, catname, sort, type, tmpcone;
  bool patch, pipelog;

#if wxUSE_WEBREQUEST
  wxWebRequest request;
  void OnWebRequestState(wxWebRequestEvent& evt);
#endif

public:
  bool OnInit();
  int OnExit();

#if !wxUSE_WEBREQUEST
  int OnRun();
#endif


};
IMPLEMENT_APP_CONSOLE(Cone)


bool Cone::OnInit()
{
  const wxString amp("&");
  wxString ra, dec, sr, magmin, magmax;
  std::list<wxString> pars;
  bool verbose = false;

  patch = false;
  pipelog = false;
  catname = "UCAC4";
  tmpcone = wxFileName::CreateTempFileName("cone_");

  // limit logging
  wxLog::DisableTimestamp();
  wxLog::SetLogLevel(wxLOG_Message);
  wxLog::SetVerbose(); // important to activate previous setup

  // read input
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

    if( line.StartsWith("PIPELOG") )
      pipelog = GetBool(line);

    if( line.StartsWith("PATCH") )
      patch = GetBool(line);

    if( line.StartsWith("URL") )
      url = GetString(line);

    if( line.StartsWith("CATNAME") )
      catname = GetString(line);

    if( line.StartsWith("OUTPUT") )
      output = GetString(line);

    if( line.StartsWith("TYPE") )
      type = GetString(line);

    if( line.StartsWith("SORT") )
      sort = GetString(line);

    if( line.StartsWith("MAGMIN") )
      magmin = wxString::FromCDouble(GetDouble(line));

    if( line.StartsWith("MAGMAX") )
      magmax = wxString::FromCDouble(GetDouble(line));

    if( line.StartsWith("PAR") )
      pars.push_back(GetString(line));
  }

  if( output.IsEmpty() ) {
    Stop("Missing output filename.");
    return false;
  }

  if( type.IsEmpty() )
    type = GetFileType(output);

  wxASSERT(!url.IsEmpty());

  // form URL
  if( magmin != "" && magmax == "" )
    url += amp + sort + "=>" + magmin;
  else if ( magmax != "" && magmin == "" )
    url += amp + sort + "=<" + magmax;
  else if( magmin != "" && magmax != "" )
    url += amp + sort + "=" + magmin + ".." + magmax;

  for(std::list<wxString>::const_iterator a = pars.begin(); a != pars.end(); ++a)
    url += (a != pars.end() ? amp : "") + *a;


#if wxUSE_WEBREQUEST

  url.Replace("http://","https://");
    // this will be obsolete when voclient dependence will be removed.

  request = wxWebSession::GetDefault().CreateRequest(this,url);
  request.SetHeader("User-Agent",USER_AGENT);
  request.SetStorage(wxWebRequest::Storage_File);
  if ( !request.IsOk() ) {
    wxLogError("Failed to create a Web request.");
    return false;
  }
  Bind(wxEVT_WEBREQUEST_STATE, &Cone::OnWebRequestState, this);
  request.Start();
#endif


  return true;
}

int Cone::OnExit()
{
  if( wxFileExists(tmpcone) )
    wxRemoveFile(tmpcone);

  return wxAppConsole::OnExit();
}

#if wxUSE_WEBREQUEST

void Cone::OnWebRequestState(wxWebRequestEvent& evt)
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
      //wxLogDebug("Request: "+response.GetURL()); // fails for undefined Unicode character
      wxLogDebug("HTTP response: %d ("+response.GetStatusText()+")",
		response.GetStatus());
      wxLogDebug("Content-Type: "+response.GetHeader("Content-Type"));
      wxLogDebug("Content-Length: "+response.GetHeader("Content-Length"));
      if( wxCopyFile(response.GetDataFile(),tmpcone) )
	result = Save();
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
    wxLogError(evt.GetErrorDescription());
    fprintf(stderr,"STOP 'Failed'\n");
    ExitMainLoop();
    SetErrorExitCode(1);
    break;
  }
}

#else


int Cone::OnRun()
{
  wxURL u(url);
  if( ! u.IsOk() ) {
    if( u.GetError() == wxURL_SNTXERR )
      wxLogError("Syntax error in the URL string.");
    else if( u.GetError() == wxURL_NOPROTO )
      wxLogError("Found no protocol which can get this URL.");
    else if( u.GetError() == wxURL_NOHOST )
      wxLogError("A host name is required for this protocol.");
    else if( u.GetError() == wxURL_NOPATH )
      wxLogError("A path is required for this protocol.");
    Stop("Badly formatted URL string.");
    return 1;
  }

  wxString server = u.GetServer();
  wxString query = u.GetPath() + wxString("?") + u.GetQuery();

  if( pipelog ) PipePrint("Connecting "+server+"...");
  VOclient voc;
  //if( ! voc.Query(this,server,query,tmpcone) ) {
  if( ! voc.Connect(server) ) {
    Stop("Failed to connect VO server.");
    return 1;
  }

  if( pipelog ) PipePrint("Downloading data...");
  if( ! voc.Get(query,tmpcone) ) {
    Stop("Download failed.");
    return 1;
  }

  Save();

  Stop("");
  return 0;
}
#endif

int Cone::Save() const
{

  if( pipelog ) PipePrint("Parsing XML data...");

  VOTable vt(tmpcone);
  wxLogDebug("XML IsOK? %s", vt.IsOk() ? "T" : "F");

  if( ! vt.IsOk() ) {
    if( pipelog ) PipePrint("Parsing XML data failed.");
    Stop("Parsing XML data failed.");
    return 1;
  }

  if( vt.HasError() ) {
    Stop(vt.GetErrorMsg());
    return 1;
  }

  if( vt.IsEmpty() ) {
    wxLogWarning("No objects found by given constrains.");
    if( pipelog ) PipePrint("No objects found.");
    Stop("");
    return 0;
  }
  else
    wxLogDebug("Cone search: %d objects found.",vt.RecordCount());

  if( ! sort.IsEmpty() ) {
    if( pipelog ) PipePrint("Sorting XML data...");
    vt.Sort(sort);
  }

  if( pipelog ) PipePrint("Data are being saved ...");
  if( type == "FITS" ) {
    FITStable ft(vt);
    if( ft.SaveFITS(output,true) ) {
      if( patch ) {
	if( catname == "UCAC4" ) {
	  if( ! PatchJohn(output) ) {
	    Stop("Application of Johnson's patch failed.");
	    return 1;
	  }
	}
	else
	  wxLogWarning("Johnson's patch is applicable only on UCAC 4.");
      }
    }
    else {
      Stop("FITS save failed.");
      return 1;
    }
  }
  else {
    if( ! vt.Save(output) ) {
      Stop("Save to "+type+" failed.");
      return 1;
    }
  }

  if( pipelog ) wxPrintf("=CONE> %d objects found\n",vt.RecordCount());
  return 0;
}

bool Cone::PatchJohn(const wxString& fitsname) const
{
  /*
     Adds two new colums: Johnson RI filters to the table.
     Also removes the original Gunn ri columns.

     https://gaia.esac.esa.int/documentation/GDR1/Data_processing/chap_cu5phot/sec_phot_calibr.html

      http://www.sdss.org/dr4/algorithms/sdssUBVRITransform.html
      see Lupton (2005)

  */

  wxLogDebug("Transforming Gunn's ri to Johnson RI ..");

  fitsfile *fits = 0;
  int status = 0;

  if( fits_open_table(&fits,fitsname.char_str(),READWRITE,&status) ) {
    fits_report_error(stderr, status);
    return false;
  }

  long nrows;
  if( fits_get_num_rows(fits,&nrows,&status) ) {
    fits_report_error(stderr, status);
    return false;
  }

  float *rmag = new float[nrows];
  float *imag = new float[nrows];
  float *e_rmag = new float[nrows];
  float *e_imag = new float[nrows];
  char **f_rmag = new char*[nrows];
  char **f_imag = new char*[nrows];
  for(int i = 0; i < nrows; i++) {
    f_rmag[i] = new char[7];
    f_imag[i] = new char[7];
  }

  int ncol, anynul;
  fits_get_colnum(fits,CASEINSEN,(char *)"rmag",&ncol,&status);
  fits_read_col(fits,TFLOAT,ncol,1,1,nrows,NULL,rmag,&anynul,&status);
  fits_delete_col(fits,ncol,&status);
  fits_get_colnum(fits,CASEINSEN,(char *)"e_rmag",&ncol,&status);
  fits_read_col(fits,TFLOAT,ncol,1,1,nrows,NULL,e_rmag,&anynul,&status);
  fits_delete_col(fits,ncol,&status);
  fits_get_colnum(fits,CASEINSEN,(char *)"f_rmag",&ncol,&status);
  fits_read_col(fits,TSTRING,ncol,1,1,nrows,NULL,f_rmag,&anynul,&status);
  fits_delete_col(fits,ncol,&status);
  fits_get_colnum(fits,CASEINSEN,(char *)"imag",&ncol,&status);
  fits_read_col(fits,TFLOAT,ncol,1,1,nrows,NULL,imag,&anynul,&status);
  fits_delete_col(fits,ncol,&status);
  fits_get_colnum(fits,CASEINSEN,(char *)"e_imag",&ncol,&status);
  fits_read_col(fits,TFLOAT,ncol,1,1,nrows,NULL,e_imag,&anynul,&status);
  fits_delete_col(fits,ncol,&status);
  fits_get_colnum(fits,CASEINSEN,(char *)"f_imag",&ncol,&status);
  fits_read_col(fits,TSTRING,ncol,1,1,nrows,NULL,f_imag,&anynul,&status);
  fits_delete_col(fits,ncol,&status);

  if( status != 0 ) {
    fits_report_error(stderr, status);
    return false;
  }

  float *Rmag = new float[nrows];
  float *Imag = new float[nrows];
  float *e_Rmag = new float[nrows];
  float *e_Imag = new float[nrows];

  for(int i = 0; i < nrows; i++) {
    if( isnormal(rmag[i]) &&  isnormal(imag[i]) ) {
      float ri = rmag[i] - imag[i];
      Rmag[i] = rmag[i] - 0.2936*ri - 0.1439;
      Imag[i] = rmag[i] - 1.2444*ri - 0.3820;
      e_Rmag[i] = e_rmag[i];
      e_Imag[i] = e_imag[i];
    }
    else {
      Rmag[i] = 99.999;
      Imag[i] = 99.999;
      e_Rmag[i] = 9.999;
      e_Imag[i] = 9.999;
    }
    /*
    wxLogDebug("%d %f %f %f %f %d %d",i,Rmag[i],Imag[i],rmag[i],imag[i],
	      isnormal(rmag[i]),isnan(imag[i]));
    */
  }

  fits_insert_col(fits,ncol,(char *)"Rmag",(char *)"1E",&status);
  fits_write_col(fits,TFLOAT,ncol,1,1,nrows,Rmag,&status);
  ncol++;
  fits_insert_col(fits,ncol,(char *)"e_Rmag",(char *)"1E",&status);
  fits_write_col(fits,TFLOAT,ncol,1,1,nrows,e_Rmag,&status);
  ncol++;
  fits_insert_col(fits,ncol,(char *)"f_Rmag",(char *)"1A",&status);
  fits_write_col(fits,TSTRING,ncol,1,1,nrows,f_rmag,&status);
  ncol++;
  fits_insert_col(fits,ncol,(char *)"Imag",(char *)"1E",&status);
  fits_write_col(fits,TFLOAT,ncol,1,1,nrows,Imag,&status);
  ncol++;
  fits_insert_col(fits,ncol,(char *)"e_Imag",(char *)"1E",&status);
  fits_write_col(fits,TFLOAT,ncol,1,1,nrows,e_Imag,&status);
  ncol++;
  fits_insert_col(fits,ncol,(char *)"f_Imag",(char *)"1A",&status);
  fits_write_col(fits,TSTRING,ncol,1,1,nrows,f_imag,&status);

  fits_write_comment(fits,"",&status);
  fits_write_comment(fits," *** NOTICE ***",&status);
  fits_write_comment(fits,
      "Rmag, Imag in Johnson system, and related columns, has been derived",
      &status);
  fits_write_comment(fits,
      "from Gunn's ri magnitudes by the transformation, Lupton (2005): ",
      &status);
  fits_write_comment(fits,
      "http://www.sdss.org/dr4/algorithms/sdssUBVRITransform.html",&status);

  fits_close_file(fits, &status);
  fits_report_error(stderr,status);

  delete[] rmag; delete[] imag; delete[] e_rmag; delete[] e_imag;
  delete[] Rmag; delete[] Imag; delete[] e_Rmag; delete[] e_Imag;
  for(int i = 0; i < nrows; i++) {
    delete[] f_rmag[i];
    delete[] f_imag[i];
  }
  delete[] f_rmag; delete[] f_imag;

  return true;
}

void Cone::Stop(const wxString& msg) const
{
  // It is just strangle workaround.
  // The code simulates STOP commands which terminate of Fortran programs.
  // We have no faith to correctness of exit codes returned by child processes
  // because a correct terminate would leads to generate 'Child process
  // (PID XXX) still alive but pipe closed so generating a close notification'
  // under some GNU/Linux distributions (Fedora). Processing of STOP mark
  // eliminates use of process execution terminate codes.
  // One looks as a feature (bug) of WX.
  if( msg == "" )
    fprintf(stderr,"STOP 0\n");
  else
    fprintf(stderr,"STOP '%s'\n",static_cast<const char *>(msg.c_str()));
}

void Cone::PipePrint(const wxString& line) const
{
  wxPrintf("=CONE> "+line+"\n");
  fflush(stdout);
}
