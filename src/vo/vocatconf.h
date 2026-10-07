/*

  VOTable parser specialized on default catalogues configuration


  Copyright © 2013 - 2026 F.Hroch (hroch@physics.muni.cz)

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


  Reference:

  http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html

*/

#ifndef _VOCATCONF_CONFIG_H_
#define _VOCATCONF_CONFIG_H_

#include <wx/wx.h>
#include <wx/xml/xml.h>
#include <map>
#include <vector>


class VOCatResources
{

  wxString name,query;
  std::map<wxString,wxString> labels;

public:

  VOCatResources();
  VOCatResources(const wxString&);
  VOCatResources(const wxString&, const wxString&, const std::map<wxString,wxString>&);
  bool IsOk() const;
  wxString GetLabel(const wxString&) const;
  wxString GetSort() const;
  wxString GetQuery() const { return query; }
  wxString GetName() const { return name; }

  std::map<wxString,wxString> GetLabels() const { return labels; }
};

class VOCatConf: public wxXmlDocument
{

  void ParseSites(const wxXmlNode *);
  void ParseCats(const wxXmlNode *);

  std::vector<VOCatResources> cats;
  std::map<wxString,wxString> sites;
  std::vector<wxString> name,ucd;
  wxString query;
  std::vector<VOCatResources>::const_iterator cat_current;
  std::map<wxString,wxString>::const_iterator site_current;

public:
  VOCatConf(const wxString& ="");
  bool IsOk() const;

  static int Replace(wxString& , const wxString, const wxString);
  static int ReplaceAll(wxString&, const std::map<wxString,wxString>&);
  wxString GetUrl(const std::map<wxString,wxString>&) const;
  VOCatResources GetCat() const;
  VOCatResources GetCat(const wxString&) const;
  wxString GetName() const;
  wxString GetSite() const;
  wxString GetSort() const;
  std::vector<wxString> GetNames() const { return name; }
  VOCatResources GetCatFile(const wxString&) const;
  VOCatResources GetCatFits(const wxString&) const;
  void UnSetCat();
  bool SetCat(const wxString&);
  bool SetSite(const wxString&);
  std::vector<VOCatResources> GetCatalogues() const { return cats; }
  std::map<wxString,wxString> GetSites() const { return sites; }

};

#endif
