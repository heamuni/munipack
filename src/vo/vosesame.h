/*

  Virtual observatory Sesame client (resolves astronomical names of objects)

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

*/

#include <wx/wx.h>
#include <wx/url.h>
#include <wx/stream.h>
#include <wx/protocol/http.h>

class VOSesame: public wxHTTP
{
  wxString object, errmsg;
  double ra, dec;
  wxURL url;
  wxInputStream *istream;

 public:
  VOSesame(const wxString&);
  virtual ~VOSesame();
  bool Connect(const wxString&);
  bool Resolve();
  bool Get(const wxString&);
  const char* GetErrMsg() const;
  double GetRa() const { return ra; }
  double GetDec() const { return dec; }
  bool Request();
  wxInputStream *GetStream() { return istream; }
};
