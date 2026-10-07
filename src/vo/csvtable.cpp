/*

  VOTable to CSV convertor

  Copyright © 2018-2026 F.Hroch (hroch@physics.muni.cz)

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


  Description of CSV format:

      http://en.wikipedia.org/wiki/Comma-separated_values

*/

#include "votable.h"
#include <wx/wx.h>
#include <wx/txtstrm.h>
#include <vector>

using namespace std;

CSVtable::CSVtable(const wxString& filename): VOTable(filename) {}

bool CSVtable::SaveCSV(wxOutputStream& output) const
{
  if( ! IsOk() ) return false;
  if( HasError() ) return false;

  if( resources.size() == 0 )
    return false;

  // We are loosers, only the first table is converted.
  const VOTableTable table(resources[0].table);

  //  wxFFileOutputStream output(file);
  wxTextOutputStream cout(output);

  vector<wxString> fs = table.GetFields();
  for(vector<wxString>::const_iterator i = fs.begin(); i != fs.end(); ++i) {
    if( i != fs.begin() ) cout << ",";
    wxString a(*i);
    a.Replace("\n",""); a.Replace("\r","");
    cout << a;
  }
  cout << endl;

  for(int n = 0; n < table.RecordCount(); n++) {

    vector<wxString> r = table.GetRecord(n);

    for(vector<wxString>::const_iterator i = r.begin(); i != r.end(); ++i) {
      if( i != r.begin() ) cout << ",";
      wxString a(*i);
      a.Replace("\n",""); a.Replace("\r","");
      cout << a;
    }

    cout << endl;
  }

  return true;
}
