/*

  Virtual Observatory Client

  Copyright © 2018, 2025 F.Hroch (hroch@physics.muni.cz)

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


#include <wx/wx.h>
#include <wx/protocol/http.h>

class VOclient: public wxHTTP
{
 public:
  VOclient();
  VOclient(const wxString&);
  bool Connect(const wxString&);
  bool Get(const wxString&, const wxString&);

};


/*

The part below is prepared for implementation via webrequest.
It is required for support of https servers. Unfortunatelly,
they require cURL library which is left unsupported in my current
wxWidgets.


#include <wx/file.h>
#include <wx/webrequest.h>

class VOclient
{
  wxFile tmpfile;

  void OnWebRequestState(wxWebRequestEvent&);
  void OnWebRequestData(wxWebRequestEvent&);

public:
  bool Query(wxEvtHandler *, const wxString&, const wxString&, const wxString&);
};
*/
