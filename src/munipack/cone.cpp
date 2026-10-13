/*

  Cone search

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

/*

   http://www.usno.navy.mil/USNO/astrometry/optical-IR-prod/icas/vo_nofs

*/

#include "munipack.h"
#include "vocatconf.h"
#include <wx/wx.h>
#include <wx/filename.h>
#include <iostream>
#include <map>

using namespace std;


MuniCone::MuniCone(wxEvtHandler *h, const MuniCmdLineParser& cmd):
  MuniAction(h,"cone")
{
  wxString a, type;
  wxString vocat, sort;

  if( cmd.Found("vocat",&a) )
    vocat = a;
  else if( FindDataDir("VOcat_conf.xml",&a) )
    vocat = a;

  VOCatConf catalogs(vocat);

  if( ! catalogs.IsOk() ) {
    MuniFatalError("Failed to parse the definition file of catalogues: `" +
		   vocat+"'.");
    return;
  }

  if( cmd.Found("list-catalogues") ) {

    list_catalogues(catalogs);
    wxQueueEvent(this,new wxProcessEvent(wxEVT_END_PROCESS,GetPid(),0));
    return;
  }
  else if( cmd.Found("list-servers") ) {

    list_servers(catalogs);
    wxQueueEvent(this,new wxProcessEvent(wxEVT_END_PROCESS,GetPid(),0));
    return;
  }

  // the cone search
  CommonOutputSingle(cmd);

  double r,ra,dec,mag;

  if( cmd.Found("par",&a) )
    Write("PAR = '"+a+"'");

  if( cmd.Found("o",&a) || cmd.Found("output",&a) ) {
    wxFileName fn(a);
    if( fn.IsOk() )
      type = fn.GetExt();
  }

  if( cmd.Found("type",&a) )
    type = a;

  if( ! type.IsEmpty() )
    Write("TYPE = '"+type.Upper()+"'");

  map<wxString,wxString> pars;

  if( cmd.Found("r",&r) || cmd.Found("radius",&r) )
    pars["SR"] = wxString::FromCDouble(r);
  else
    pars["SR"] = "0.1";

  if ( cmd.GetParamCount() != 2 )
    MuniFatalError("Unspecified position? Use: RA DEC (both in degress, sexadecimal).");

  a = cmd.GetParam(0);
  if( a.ToDouble(&ra) )
    pars["RA"] = wxString::FromCDouble(ra);
  else
    MuniFatalError("Right Ascension unrecognized.");

  a = cmd.GetParam(1);
  if( a.ToDouble(&dec) )
    pars["DEC"] = wxString::FromCDouble(dec);
  else
    MuniFatalError("Declination unrecognized.");

  if( cmd.Found("server",&a) ) {
    bool b = catalogs.SetSite(a);
    if( ! b )
      MuniFatalError("Selected VO server unavailable (try --list-servers).");
  }

  if( cmd.Found("id",&a) ) {
    pars["ID"] = a;
    catalogs.UnSetCat();
  }

  if( cmd.Found("url",&a) ) {
    VOCatConf::ReplaceAll(a,pars);
    Write("URL = '"+a+"'");
  }
  else {

    if( cmd.Found("cat",&a) || cmd.Found("c",&a) ) {
      bool b = catalogs.SetCat(a);
      if( ! b )
	MuniFatalError("Required catalogue unavailable (try --list-catalogues).");
    }

    Write("URL = '"+catalogs.GetUrl(pars)+"'");
    sort = catalogs.GetSort();

    Write("CATNAME = '" + catalogs.GetName() + "'");
  }

  if( cmd.Found("s",&a) || cmd.Found("sort",&a) )
    sort = a;

  if( sort != "" )
    Write("SORT = '"+sort+"'");

  if( cmd.Found("magmin",&mag) )
    Write("MAGMIN = "+wxString::FromCDouble(mag));

  if( cmd.Found("magmax",&mag) )
    Write("MAGMAX = "+wxString::FromCDouble(mag));

  if( cmd.Found("Johnson-patch") )
    Write("PATCH = T");

  return;
}

void MuniCone::list_catalogues(const VOCatConf& catalogs) const
{
  cout << "Available catalogues:" << endl;
  vector<wxString> names(catalogs.GetNames());
  for(vector<wxString>::const_iterator n = names.begin(); n != names.end(); ++n ) {
    if( n - names.begin() == 1 ) {
      cout.width(2);
      cout << "  ";
    }
    cout.width(10);
    cout << *n;
  }
  cout << endl;

  vector<wxString> ucd;
  ucd.push_back("CAT");
  ucd.push_back("ID");
  ucd.push_back("POS_EQ_RA");
  ucd.push_back("POS_EQ_DEC");
  ucd.push_back("POS_EQ_PMRA");
  ucd.push_back("POS_EQ_PMDEC");
  ucd.push_back("PHOT_MAG");

  vector<VOCatResources>::const_iterator c;
  vector<VOCatResources> cats(catalogs.GetCatalogues());
  for(c = cats.begin(); c != cats.end(); ++c ){
    for(vector<wxString>::const_iterator u = ucd.begin(); u != ucd.end(); ++u) {
      if( u - ucd.begin() == 1 ) {
	cout.width(2);
	if( catalogs.GetName() == c->GetName() ) cout << "*"; else cout << " ";
      }
      cout.width(10);
      cout << c->GetLabel(*u);
    }
    cout << endl;
  }
}

void MuniCone::list_servers(const VOCatConf& catalogs) const
{
  cout << "Available Virtual observatory servers:" << endl;
  cout.width(20);
  cout << "Alias" << "   " << "URL" << endl;

  map<wxString,wxString>::const_iterator s;
  const map<wxString,wxString> sites(catalogs.GetSites());
  for(s = sites.begin(); s != sites.end(); ++s ) {
    wxString flag(catalogs.GetSite() == s->second ? "*" : " ");
    cout.width(20);
    cout << s->first << " " << flag << " " << s->second  << endl;
  }
  cout << "(Some servers works erratically.)" << endl;
}
