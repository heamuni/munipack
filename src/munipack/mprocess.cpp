/*

  Launcher of an external process

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

*/

#include "mprocess.h"
#include <wx/wx.h>
#include <wx/app.h>
#include <wx/stream.h>
#include <wx/txtstrm.h>
#include <wx/wfstream.h>
#include <wx/log.h>
#include <wx/event.h>
#include <wx/process.h>
#include <wx/regex.h>
#include <wx/stopwatch.h>
#include <wx/datetime.h>
#include <cstring>

#define EXITCODE 13

MuniProcess::MuniProcess(wxEvtHandler *h, const wxString& c,
  			 const wxArrayString& args):
  wxProcess(h),command(c),argv(0),exitcode(EXITCODE),timer(this),tick(250),Index(0),
  echo(true)
{
  size_t argc = args.GetCount() + 2;
  argv = new wxChar*[argc];
  argv[0] = chardup(command.wc_str());
  for(size_t i = 0; i < args.GetCount(); i++)
    argv[i+1] = chardup(args[i].wc_str());
  argv[argc-1] = 0;

  SetRuntimeEnvironment();
  SetFitsKeys();

  Bind(wxEVT_TIMER,&MuniProcess::OnTimer,this);
  // Idle events are, perhaps, unimplemented for wxConsoleApps.
}

MuniProcess::~MuniProcess()
{
  Unbind(wxEVT_TIMER,&MuniProcess::OnTimer,this);
  deallocate_argv();
}

void MuniProcess::deallocate_argv()
{
  if( argv ) {
    for(size_t i = 0; argv[i] != 0; i++)
      delete[] argv[i];
    delete[] argv;
    argv = 0;
  }
}

void MuniProcess::SetCmdLine(int argc, wxChar **args)
{
  deallocate_argv();

  argv = new wxChar*[argc+1];
  for(int i = 0; i < argc; i++)
    argv[i] = chardup(args[i]);
  argv[argc] = 0;
}


// Echoes output of the childerns; Console utilities has enabled it
// be default, GUI has to change it switch-off.
void MuniProcess::SetEcho(bool e)
{
  echo = e;
}


void MuniProcess::OnStart()
{
  wxASSERT(command != "");

  wxLogDebug("Launching `" + command + "' ...");

  Redirect();
  long pid = wxExecute(argv,wxEXEC_ASYNC,this);

  if( pid <= 0 ) {
    wxLogError("Failed to launch the external command `" + command + "'.");
    return;
  }

  timer.Start(tick);
  stopwatch.Start();
  Send();

  wxASSERT(wxProcess::Exists(pid));
}

void MuniProcess::SetRuntimeEnvironment()
{
  // Switch-off buffering of gfortran's stdout and stderr.
  // We need this setting for on-the-fly parsing of outputs.
  wxSetEnv("GFORTRAN_UNBUFFERED_PRECONNECTED","Y");

  // set path for libexec, generally non-portable (!)
  wxString xpath;
  wxGetEnv("PATH",&xpath);

#ifdef MUNIPACK_LIBEXEC_DIR
  xpath = wxString(MUNIPACK_LIBEXEC_DIR ":") + xpath;
#endif

  wxString libexecpath;
  if( wxGetEnv("MUNIPACK_LIBEXEC_PATH",&libexecpath) )
    xpath = libexecpath + wxString(":") + xpath;

  wxSetEnv("PATH",xpath);

}

void MuniProcess::SetFitsKeys()
{
  // redefine FITS keywords by environment variables

  const char *keys[] = {
    "FITS_KEY_FILTER",
    "FITS_KEY_TEMPERATURE",
    "FITS_KEY_DATEOBS",
    "FITS_KEY_EXPTIME",
    "FITS_KEY_OBJECT",
    "FITS_KEY_SATURATE",
    "FITS_KEY_READNOISE",
    "FITS_KEY_GAIN",
    "FITS_KEY_AREA",
    "FITS_KEY_EPOCH",
    "FITS_KEY_LATITUDE",
    "FITS_KEY_LONGITUDE",
    "FITS_KEY_ALTITUDE",
    "FITS_KEY_AIRMASS",
    "FITS_KEY_TIME",
    0
  };

  for(size_t i = 0; keys[i] != 0; i++) {
    wxString var, key(keys[i]);
    if( wxGetEnv(key,&var) )
      Input.Add(key + " = '" + var + "'" );
  }
}


void MuniProcess::Write(const wxArrayString& i)
{
  wxASSERT(Index == 0);

  for(size_t l = 0; l < i.GetCount(); l++)
    Input.Add(i[l]);
}

void MuniProcess::Write(const char *line)
{
  wxASSERT(Index == 0);
  Input.Add(wxString(line));
}

void MuniProcess::Write(const wxString& fmt, ...)
{
  wxASSERT(Index == 0);

  wxString line;

  va_list par;
  va_start(par, fmt);
  line.PrintfV(fmt,par);
  va_end(par);

  Input.Add(line);
}

wxKillError MuniProcess::Kill(wxSignal sig, int flags)
{
  return wxProcess::Kill(GetPid(),sig,flags);
}

void MuniProcess::OnTimer(wxTimerEvent& event)
{
  //  wxLogDebug("MuniProcess::OnTimer");

  // Fill the input of a subprocess
  Send();

  /*
    To give a chance of GUI/CLI to be updated, we should periodically
    interrupt the data streams. The limit for each stream
    is given by time ticks and one would be short (about 10-50ms).
  */
  wxStopWatch sw;
  sw.Start();
  while( Receive() && sw.Time() < tick / 5)
    ;
}

void MuniProcess::OnTerminate(int pid, int status)
{
  timer.Stop();
  while( Receive() ) // recieve all remaining messages
    ;

  /*
    On terminate, we should test both the exit code and the error output.
    The testing of exit code only is inappropriate because the implementation
    of wxExecute returns -1 in case when waitpid gives 0 (no changes).

    The outline is little bit complicated by using of some
    Fortran utilities which prints STOP <NUMBER>[string] to indicate
    their return status.

    By the way, this code looks for 'STOP 0' string over final
    error output lines of every subprocces; the approach ensures us
    that the process finished correctly.

    There is also another exception. When user requested killing
    of a process (no STOP is emmited), we always returns non-zero.

    Shortly:

      Just when both 'status' and STOP indicates 0 (zero) the exit 0
      is returned. All others alternatives gives non-zero status.

  */

  if( status == 0 ) {

    for(int i = Error.GetCount() - 1; i >= 0; i--) {
      int code;
      if( StopLine(Error.Item(i),code) ) {
	exitcode = code;
	break;
      }
    }
  }
  else
    exitcode = status;

  wxTimeSpan ts(wxTimeSpan::Milliseconds(stopwatch.Time()));
  wxLogDebug("MuniProcess::OnTerminate: Exitcode: %d, "+
	     ts.Format("Elapsed time: %Hh %Mm %S.%ls"),exitcode);

  // if subprocess failed, dump logs (in debug mode)
  if( exitcode != 0 && echo ) {
    if( wxLog::GetLogLevel() == wxLOG_Debug )
      SaveErrorLog();
    else
      wxLogError("Execution failed (consider re-run with --verbose).");
  }

  // notify the parent process
  wxQueueEvent(this,new wxProcessEvent(wxEVT_END_PROCESS,pid,exitcode));

}

void MuniProcess::Send()
{
  // Sends the input stream to child

  // we're waiting in rest for child finish
  if( Index == Input.GetCount() ) return;

  // the input is feeded for a limited time only (see Recive())
  wxStopWatch sw;
  sw.Start();
  while( IsInputOpened() && Index < Input.GetCount() && sw.Time() < tick/5 ) {
    wxTextOutputStream out(*GetOutputStream());

    /*
      While documentation of wxTextOutputStream recommends notation:

      out << Input[Index] << endl;

      I found that some lines are missing, if the output buffer
      is temporary full: LastWrite() returns only 1 byte for
      lines with width > 1 as consequence of putting, and filling,
      the endl mark, that suppress to check of Input write status.

      The portability is holded due WriteString().
    */
    out.WriteString(Input[Index]+"\n");

    if( GetOutputStream()->GetLastError() == wxSTREAM_WRITE_ERROR )
      wxLogFatalError("mprocess send(): "
		      "generic write error on the last write call.");

    // if no byte is successfully written, then the output buffer is full
    // of a garbage; we're going to idle state to enable processing
    // of delivered data
    if( GetOutputStream()->LastWrite() == 0 ) break;
    // the assumption: a complete line at once, or nothing, is transfered

    wxLogDebug(Input[Index]);
    Index++;
  }

  // close the stream, no more input data are available
  if( Index == Input.GetCount() )
    CloseOutput();

}

bool MuniProcess::Receive()
{
  //  wxLogDebug("MuniProcess::Receive");
  bool debug = wxLog::GetLogLevel() == wxLOG_Debug;

  bool more = false;

  // Output stream
  if( IsInputAvailable() ) {
    wxTextInputStream out(*GetInputStream());
    wxString line = out.ReadLine();
    if( GetInputStream()->LastRead() > 0 ) {
      Output.Add(line);
      if( echo )
	fprintf(stdout,"%s\n",static_cast<const char *>(line.char_str()));
    }
    more = true;
  }

  // Error stream
  if( IsErrorAvailable() ) {
    wxTextInputStream err(*GetErrorStream());
    wxString line = err.ReadLine();
    if( GetErrorStream()->LastRead() > 0 ) {
      Error.Add(line);

      // filter lines with STOP 0 in non-debug mode
      int code;
      if( StopLine(line,code) && code == 0 && ! debug )
	;
      else
	fprintf(stderr,"%s\n",static_cast<const char *>(line.char_str()));
    }
    more = true;
  }

  return more;
}

bool MuniProcess::StopLine(const wxString& line, int& stopcode) const
{
  // Recognises 'STOP <NUMBER>' patterns, decode the number.

  stopcode = EXITCODE; // default code, applied for: STOP [string]

  // One filters lines with STOP mark
  wxRegEx res("^[[:space:]]*STOP",wxRE_DEFAULT|wxRE_ICASE);
  wxASSERT(res.IsValid());

  if( res.Matches(line) ) {
    // STOP is presented
    wxRegEx re("^[[:space:]]*STOP[[:space:]]+([[:digit:]])[[:space:]]*$",
	       wxRE_DEFAULT|wxRE_ICASE);
    wxASSERT(re.IsValid());

    if( re.Matches(line) ) {
      // STOP <NUMBER>?
      wxString a(re.GetMatch(line,1));
      long s;
      if( a.ToLong(&s) )
	stopcode = s;
    }
    return true;
  }

  return false;
}

void MuniProcess::SaveErrorLog() const
{
  wxString filename = command + ".log";

  wxLogError("Executing of the utility `" + command +
	     "' finished with a failure (see `"+filename+"'"
	     " for details).");

  wxFFileOutputStream logfile(filename);

  if( ! logfile.IsOk() ) {
    wxLogError("Failed to create `"+filename+"'");
    return;
  }

  wxTextOutputStream log(logfile);
  log.WriteString("=== "+command+" ===\n");

  log.WriteString("=> Command line: ");
  for(size_t i = 0; argv[i] != 0; i++)
    log.WriteString(wxString(argv[i])+" ");
  log.WriteString("\n");

  log.WriteString("=> Standard input:\n");
  for(size_t i = 0; i < Input.GetCount(); i++)
    log.WriteString(Input[i]+"\n");

  log.WriteString("=> Standard output:\n");
  for(size_t i = 0; i < Output.GetCount(); i++)
    log.WriteString(Output[i]+"\n");

  log.WriteString("=> Error output:\n");
  for(size_t i = 0; i < Error.GetCount(); i++)
    log.WriteString(Error[i]+"\n");

  /*
    There is also idea to save output into a structured document
    like XML to enable additional processing.
  */
}

wxChar *MuniProcess::chardup(const wxChar *src) const
{
  size_t size = (wxStrlen(src) + 1) * sizeof(wxChar);
  wxChar *dst = new wxChar[size];
  memcpy(dst,src,size);
  return dst;
}
