/*

  xmunipack - main

  Copyright © 2009-2015, 2019-22 F.Hroch (hroch@physics.muni.cz)

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

#include "xmunipack.h"
#include <wx/wx.h>
#include <wx/app.h>
#include <wx/cmdline.h>
#include <wx/utils.h>


wxIMPLEMENT_APP(XMunipack);

bool XMunipack::OnInit()
{
  config = 0;
  view = 0;
  browser = 0;

  // switch-off timestamps in log
  wxLog::DisableTimestamp();

  // default log level prints errors only
  wxLog::SetLogLevel(wxLOG_Error);

  // command line parameters
  wxCmdLineParser cmd(argc,argv);
  OnInitCmdLine(cmd);
  cmd.AddSwitch("","version","print version and license");
  cmd.AddParam("filename",wxCMD_LINE_VAL_STRING,wxCMD_LINE_PARAM_OPTIONAL);

  if( cmd.Parse() == 0 ) {

    if( cmd.Found("verbose") )
      wxLog::SetLogLevel(wxLOG_Debug);
    else
      wxSetAssertHandler(NULL);

    if( cmd.Found("version") ) {
      wxPrintf("%s %s, %s\n",PACKAGE_NAME,PACKAGE_VERSION,PACKAGE_COPYLEFT);
      wxPrintf("%s\n\n",PACKAGE_DESCRIPTION);
      wxPrintf("This program comes with ABSOLUTELY NO WARRANTY;\nfor details, see the GNU General Public License, version 3 or later.\n");
      return false;
    }

  }
  else
    return false;

  // look for passed file
  wxString file;
  if( cmd.GetParamCount() > 0 )
    file = cmd.GetParam();


  // initialisation
  wxHandleFatalExceptions();
  wxInitAllImageHandlers();

  /*
   https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html

"$XDG_CONFIG_HOME defines the base directory relative to which user-specific configuration files should be stored. If $XDG_CONFIG_HOME is either not set or empty, a default equal to $HOME/.config should be used.
"
  */

  wxString config_file;
  wxString xdg_config_home;
  if( wxGetEnv("XDG_CONFIG_HOME",&xdg_config_home) && xdg_config_home != "" ) {
    wxFileName filename(xdg_config_home,"xmunipack.conf");
    config_file = filename.GetFullPath();
  }
  else {
    wxString home;
    if( wxGetEnv("HOME",&home) && home != "" ) {
      wxFileName filename(home,"xmunipack.conf");
      filename.AppendDir(".config");
      config_file = filename.GetFullPath();
    }
  }
  config = new MuniConfig(config_file);

#ifdef __WXMAC__
  // all main windows?
  SetExitOnFrameDelete(true);
  wxMenuBar *menubar = new wxMenuBar;
  wxMenuBar::MacSetCommonMenuBar(menubar);
#endif


  if( ! file.IsEmpty() ) {

    view = new MuniView(NULL,config);
    view->Show();
    view->LoadFile(file);
    SetTopWindow(view);
  }
  else {

    browser = new MuniBrowser(NULL,config);
    browser->Show(true);

  }

  return true;
}


int XMunipack::OnExit()
{
  delete config;
  return 0;
}


void XMunipack::OnEventLoopEnter(wxEventLoopBase* WXUNUSED(loop))
{
  //wxLogDebug("XMunipack::OnEventLoopEnter");
  if( view )
    view->CreateFSWatch();
}
