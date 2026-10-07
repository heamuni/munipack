/*

  xmunipack - help and about

  Copyright © 2009-2013, 2019-26 F.Hroch (hroch@physics.muni.cz)

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
#include "help.h"
#include <wx/wx.h>
#include <wx/filename.h>
#include <wx/utils.h>
#include <wx/aboutdlg.h>

#define COPYLEFT PACKAGE_COPYLEFT "\n\n This program comes with ABSOLUTELY NO WARRANTY;\nfor details, see the GNU General Public License, version 3 or later."

MuniAbout::MuniAbout(const wxIcon& icon)
{
  wxAboutDialogInfo info;
  info.SetName(PACKAGE_NAME);
  info.SetVersion(PACKAGE_VERSION);
  info.SetDescription(PACKAGE_DESCRIPTION);
  info.SetCopyright(COPYLEFT);
  info.SetWebSite(PACKAGE_URL);
  info.SetIcon(icon);
  info.AddDeveloper("Filip Hroch");
  wxAboutBox(info);
}


bool MuniHelp(const wxString& munipack_html_dir, const wxString& subpage)
{
  wxString webpage(PACKAGE_URL);

  wxFileName path(munipack_html_dir,(subpage == "") ? "munipack.html" : subpage);
  if( path.FileExists() )
    webpage = "file://" + path.GetFullPath();
  else
    webpage = webpage + subpage;

  return wxLaunchDefaultBrowser(webpage,wxBROWSER_NOBUSYCURSOR);
}
