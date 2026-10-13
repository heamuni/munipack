/*

  Munipack -- the command line interface

  Copyright © 2009 - 2026 F.Hroch (hroch@physics.muni.cz)

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
#include "vocatconf.h"
#include <wx/wx.h>
#include <wx/app.h>
#include <wx/cmdline.h>
#include <wx/process.h>


// Prints the error message and terminates munipack with the exit code 1.
// (wxLogFatalError aborts: its exit code, 134 or 139, looks like a crash.)
[[noreturn]] void MuniFatalError(const wxString&);

// The same with a format string (as wxString::Format).
template<typename T, typename... Args>
[[noreturn]] void MuniFatalError(const wxString& fmt, const T& a,
				 const Args&... args)
{
  MuniFatalError(wxString::Format(fmt,a,args...));
}


class MuniCmdLineParser: public wxCmdLineParser
{
  wxArrayString lstfiles;
  void InputByFile(const wxString&);

public:
  MuniCmdLineParser(int, wxChar **);
  MuniCmdLineParser(const wxString&, int, wxChar **);

  //  void Usage() const;
  wxArrayString GetFiles() const;

};


class Munipack: public wxAppConsole
{
  int exitcode;
  MuniProcess *proc;
  void OnFinish(wxProcessEvent&);

public:

  bool OnInit();
  int OnRun();
  int OnExit();

};


class MuniAction: public MuniProcess
{
  bool advanced;
  wxString pattern, mask, format, target;

protected:

  void SetCAF(const MuniCmdLineParser&);
  void CommonOptionsBitpix(const MuniCmdLineParser&);
  void CommonOptionsPhCorr(const MuniCmdLineParser&);
  void CommonOptionsEnableOverwrite(const MuniCmdLineParser&);
  void CommonOptionsPhq(const MuniCmdLineParser&);
  void CommonOutputMultiple(const MuniCmdLineParser&);
  void CommonOutputSingle(const MuniCmdLineParser&, const wxString& =wxEmptyString);
  void EnvironmentPhCorr();
  void WriteFilesMultiple(const wxArrayString&);
  void WriteFilesSingle(const wxArrayString&);

  wxString GetAdvanced(const wxString&) const;
  wxString GetOutput(const wxString&) const;
  wxString GetFormat(size_t) const;

  void SetAdvanced(bool b) { advanced = b; }
  void SetPattern(const wxString& a) { pattern = a; }
  void SetMask(const wxString& a) { mask = a; }
  void SetFormat(const wxString& a) { format = a; }
  bool SetTargetDir(const wxString&);
  bool FindDataDir(const wxString&, wxString*) const;

  void apstr(const wxString&, long&, wxString&) const;
  void dblstr(const wxString&, const wxString&, size_t&, wxString&) const;

public:
  MuniAction(wxEvtHandler *, const wxString&);

};

class MuniBias: public MuniAction
{
public:
  MuniBias(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniDark: public MuniAction
{
public:
  MuniDark(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniFlat: public MuniAction
{
public:
  MuniFlat(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniPhcorr: public MuniAction
{
public:
  MuniPhcorr(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniCone: public MuniAction
{
  void list_servers(const VOCatConf&) const;
  void list_catalogues(const VOCatConf&) const;

public:
  MuniCone(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniFits: public MuniAction
{
public:
  MuniFits(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniAstrometry: public MuniAction
{
public:
  MuniAstrometry(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniPhcal: public MuniAction
{
public:
  MuniPhcal(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniAphot: public MuniAction
{
public:
  MuniAphot(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniGphot: public MuniAction
{
public:
  MuniGphot(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniPsfphot: public MuniAction
{
public:
  MuniPsfphot(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniFind: public MuniAction
{
public:
  MuniFind(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniKombine: public MuniAction
{
public:
  MuniKombine(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniColouring: public MuniAction
{
public:
  MuniColouring(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniVotable: public MuniAction
{
public:
  MuniVotable(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniSesame: public MuniAction
{
public:
  MuniSesame(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniCross: public MuniAction
{
public:
  MuniCross(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniPhfotran: public MuniAction
{
public:
  MuniPhfotran(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniTimeseries: public MuniAction
{
public:
  MuniTimeseries(wxEvtHandler *, const MuniCmdLineParser&);
};

class MuniArtificial: public MuniAction
{
public:
  MuniArtificial(wxEvtHandler *, const MuniCmdLineParser&);
};
