/*

  Munipack - command line interface

  Copyright © 2010 - 2026 F.Hroch (hroch@physics.muni.cz)

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

#include "../config.h"
#include "munipack.h"
#include <wx/wx.h>
#include <wx/app.h>

wxIMPLEMENT_APP_CONSOLE(Munipack);

bool Munipack::OnInit()
{
  exitcode = 3;
  proc = 0;
  SetCLocale();

  // switch-off timestamps in logs
  wxLog::DisableTimestamp();

  // the default log level print errors and warnings only
  wxLog::SetLogLevel(wxLOG_Warning);

  // implemented actions
  const char *actions[] = {
    "bias", "dark", "flat", "phcorr", "colouring", "cone", "astrometry", "phcal",
    "votable", "aphot", "gphot", "find", "timeseries", "fits", "kombine",
    "phfotran", "cross", "artificial", "sesame", 0
  };

  // recognise the action
  wxString action;
  for(int i = 0; i < argc; i++) {
    for(int j = 0; actions[j] != 0; j++)
      if( argv[i] == actions[j] )
	action = actions[j];
  }

  if( action == "" ) {
    MuniCmdLineParser cmd(argc,argv);
    OnInitCmdLine(cmd);
    if( cmd.Parse(true) == 0 ) {
      if( cmd.Found("version") ) {
	wxPrintf("%s %s, %s\n",PACKAGE_NAME,PACKAGE_VERSION,PACKAGE_COPYLEFT);
	wxPrintf("%s\n\n",PACKAGE_DESCRIPTION);
	wxPrintf("This program comes with ABSOLUTELY NO WARRANTY;\n"
		 "for details, see the GNU General Public License,"
		 "version 3 or later.\n");
      }
      else
	cmd.Usage();
    }
    return false;
  }

  // process command-line
  MuniCmdLineParser cmd(action,argc,argv);
  OnInitCmdLine(cmd);

  if( cmd.Parse(true) == 0 ) {

    if( cmd.Found("verbose") )
      wxLog::SetLogLevel(wxLOG_Debug);

    if ( action == "bias" )
      proc = new MuniBias(this,cmd);
    else if ( action == "dark" )
      proc = new MuniDark(this,cmd);
    else if ( action == "flat" )
      proc = new MuniFlat(this,cmd);
    else if ( action == "phcorr" )
      proc = new MuniPhcorr(this,cmd);
    else if ( action == "cone" )
      proc = new MuniCone(this,cmd);
    else if ( action == "fits" )
      proc = new MuniFits(this,cmd);
    else if ( action == "astrometry" )
      proc = new MuniAstrometry(this,cmd);
    else if ( action == "phcal" )
      proc = new MuniPhcal(this,cmd);
    else if ( action == "aphot" )
      proc = new MuniAphot(this,cmd);
    else if ( action == "gphot" )
      proc = new MuniGphot(this,cmd);
    else if ( action == "find" )
      proc = new MuniFind(this,cmd);
    else if ( action == "kombine" )
      proc = new MuniKombine(this,cmd);
    else if ( action == "timeseries" )
      proc = new MuniTimeseries(this,cmd);
    else if ( action == "colouring" )
      proc = new MuniColouring(this,cmd);
    else if ( action == "votable" )
      proc = new MuniVotable(this,cmd);
    else if ( action == "sesame" )
      proc = new MuniSesame(this,cmd);
    else if ( action == "cross" )
      proc = new MuniCross(this,cmd);
    else if( action == "artificial" )
      proc = new MuniArtificial(this,cmd);

    Bind(wxEVT_END_PROCESS,&Munipack::OnFinish,this);

    return true;
  }

  return false;
}

int Munipack::OnRun()
{
  wxASSERT(proc);

  proc->OnStart();
  wxAppConsole::OnRun();
  return exitcode;
}

int Munipack::OnExit()
{
  wxASSERT(proc);

  Unbind(wxEVT_END_PROCESS,&Munipack::OnFinish,this);
  delete proc;
  return 0;
}

void Munipack::OnFinish(wxProcessEvent& event)
{
  ExitMainLoop();
  exitcode = event.GetExitCode();
}
