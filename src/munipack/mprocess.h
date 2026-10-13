/*

  Managing of external processes

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

#ifndef _MUNIPACK_CONFIG_H_
#define _MUNIPACK_CONFIG_H_

#include <wx/wx.h>
#include <wx/app.h>
#include <wx/process.h>
#include <wx/event.h>
#include <wx/stopwatch.h>
#include <queue>

class MuniProcess: public wxProcess
{
  const wxString command;
  wxChar **argv;
  int exitcode;
  wxTimer timer;
  int tick;
  wxStopWatch stopwatch;
  size_t Index;
  wxArrayString Input,Output,Error;
  bool echo;
  bool broken;

  wxChar *chardup(const wxChar *) const;
  void deallocate_argv();
  void Send();
  bool Receive();
  void OnTimer(wxTimerEvent&);
  void OnTerminate(int,int);
  bool StopLine(const wxString&, int&) const;
  void SaveErrorLog() const;

  void SetRuntimeEnvironment();
  void SetFitsKeys();

public:
  MuniProcess(wxEvtHandler *, const wxString&,
	      const wxArrayString& =wxArrayString());
  virtual ~MuniProcess();

  void SetCmdLine(int, wxChar **);

  virtual void OnPreProcess() {}
  virtual void OnStart();
  virtual void OnPostProcess() {}
  wxKillError Kill(wxSignal sig =wxSIGTERM, int flags =wxKILL_NOCHILDREN);

  void Write(const char *);
  void Write(const wxString&, ...);
  void Write(const wxArrayString&);
  int GetExitCode() const { return exitcode; }
  wxArrayString GetInput() const { return Input; }
  wxArrayString GetOutput() const { return Output; }
  wxArrayString GetErrors() const { return Error; }
  wxString GetCommand() const { return command; }
  void SetEcho(bool);

};

class MuniPipe: public wxEvtHandler
{
  wxEvtHandler *parent;
  std::queue<MuniProcess *> procs;
  wxArrayString OutputBuffer,ErrorBuffer;
  int exitcode;

  void OnFinish(wxProcessEvent&);

public:
  MuniPipe(wxEvtHandler * =0);
  virtual ~MuniPipe();

  void push(MuniProcess *);
  void Start();
  void Stop();
  wxArrayString GetOutput() const;
  wxArrayString GetErrors() const;

};

#endif
