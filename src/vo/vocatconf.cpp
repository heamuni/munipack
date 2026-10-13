/*

  VOTable parser specialized on default catalogues configuration

  Copyright © 2013-2026 F.Hroch (hroch@physics.muni.cz)

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


#include "vocatconf.h"
#include "votable.h"

#include <wx/wx.h>
#include <wx/xml/xml.h>
#include <wx/url.h>
#include <wx/regex.h>
#include <wx/filefn.h>
#include <map>
#include <vector>
#include <fitsio.h>


using namespace std;

// ---- VOCatResources

VOCatResources::VOCatResources() {}

VOCatResources::VOCatResources(const wxString& a): name(a) {}

VOCatResources::VOCatResources(const wxString& a, const wxString& q,
			       const std::map<wxString,wxString>& l):
  name(a),query(q),labels(l) {}

bool VOCatResources::IsOk() const
{
  return ! labels.empty() && ! query.IsEmpty();
}

wxString VOCatResources::GetLabel(const wxString& name) const
{
  map<wxString,wxString>::const_iterator l;
  for(l = labels.begin(); l != labels.end(); ++l) {
    if( l->first == name )
      return l->second;
  }

  return "";
}


wxString VOCatResources::GetSort() const
{
  return GetLabel("PHOT_MAG");
}




//  -----   VOCatConf


VOCatConf::VOCatConf(const wxString& path): wxXmlDocument()
{
  // an error leaves the configuration empty: IsOk() is false
  // (callers report it; wxLogFatalError aborted the program)
  if( ! wxFileExists(path) ) {
    wxLogError("VOCatConf::VOCatConf(): There is no file `"+path+"'.");
    return;
  }

  if( !( Load(path) && GetRoot()->GetName() == "VOTABLE" ) ) {
    wxLogError("File `"+path+"' not recognized as a VOTable.");
    return;
  }

  const wxXmlNode *node = GetRoot()->GetChildren();
  while (node) {

    if( node->GetName() == "RESOURCE" ) {

      if( node->GetAttribute("name") == "SITE" )
	ParseSites(node);

      if( node->GetAttribute("name") == "CATALOGUE" )
	ParseCats(node);
    }

    node = node->GetNext();
  }

  SetSite("CDS");
  SetCat("UCAC4");
}

bool VOCatConf::IsOk() const
{
  return ! cats.empty() && ! sites.empty() && wxXmlDocument::IsOk();
}


void VOCatConf::ParseSites(const wxXmlNode *rnode)
{
  const wxXmlNode *node = rnode->GetChildren();
  while (node) {

    if( node->GetName() == "TR" ) {

      vector<wxString> td;

      const wxXmlNode *tnode = node->GetChildren();
      while( tnode ) {
	td.push_back(tnode->GetNodeContent());
	tnode = tnode->GetNext();
      }

      wxASSERT(td.size() == 2);
      sites[td[0]] = td[1];
    }
    else

      ParseSites(node);

    node = node->GetNext();
  }

}

void VOCatConf::ParseCats(const wxXmlNode *rnode)
{
  const wxXmlNode *node = rnode->GetChildren();
  while (node) {

    if( node->GetName() == "FIELD" ) {

      name.push_back(node->GetAttribute("name"));
      ucd.push_back(node->GetAttribute("ucd"));

    }
    else if( node->GetName() == "PARAM" && node->GetAttribute("name") == "query" )

      query = node->GetAttribute("value");

    else if( node->GetName() == "TR" ) {


      vector<wxString> td;

      const wxXmlNode *tnode = node->GetChildren();
      while( tnode ) {
	td.push_back(tnode->GetNodeContent());
	tnode = tnode->GetNext();
      }
      wxASSERT(td.size() == ucd.size() && td.size() > 1);

      map<wxString,wxString> labels;
      vector<wxString>::const_iterator l;
      for(l = td.begin(); l != td.end(); ++l)
	labels[ucd[l-td.begin()]] = *l;

      wxString newq(query);
      Replace(newq,"ID",td[1]);

      VOCatResources r(td[0],newq,labels);
      if( ! r.IsOk() )
	wxLogError("Catalogue `"+td[0]+"' looks weird.");
      cats.push_back(r);

    }
    else

      ParseCats(node);

    node = node->GetNext();
  }

}

int VOCatConf::Replace(wxString& text, const wxString from, const wxString to)
{
  int ret = -1;

  wxRegEx regex("[$][{]"+from+"[}]");

  if( regex.IsValid() && regex.Matches(text) )
    ret = regex.Replace(&text,to);

  return ret;
}

int VOCatConf::ReplaceAll(wxString& text, const map<wxString,wxString>& replace)
{
  int cts = 0;

  map<wxString,wxString>::const_iterator r;
  for(r = replace.begin(); r != replace.end(); ++r)
    cts += Replace(text,r->first,r->second) == 1 ? 1 : 0 ;

  return cts; // returns number of succefull replacements
}


bool VOCatConf::SetCat(const wxString& name)
{
  vector<VOCatResources>::const_iterator c;
  for(c = cats.begin(); c != cats.end(); ++c) {
    wxString cat(c->GetName());
    if( ! name.IsEmpty() && cat.Find(name) == 0 ) {
      cat_current = c;
      return true;
    }
  }
  return false;
}

bool VOCatConf::SetSite(const wxString& name)
{
  map<wxString,wxString>::const_iterator s;
  for(s = sites.begin(); s != sites.end(); ++s) {
    wxString a(s->first);
    if( ! name.IsEmpty() && a.Find(name) == 0 ) {
      site_current = s;
      return true;
    }
  }
  return false;
}

VOCatResources VOCatConf::GetCat() const
{
  return *cat_current;
}

wxString VOCatConf::GetName() const
{
  return cat_current != cats.end() ? cat_current->GetName() : "";
}

wxString VOCatConf::GetSite() const
{
  return site_current->second;
}

wxString VOCatConf::GetSort() const
{
  if( cat_current != cats.end() )
    return cat_current->GetLabel("PHOT_MAG");
  else
    return "";
}

wxString VOCatConf::GetUrl(const map<wxString,wxString>& replace) const
{
  wxString url(cat_current != cats.end()? cat_current->GetQuery() : query);
  ReplaceAll(url,replace);
  wxURI uri(site_current->second+"/"+url);
  return uri.BuildUnescapedURI();
}

void VOCatConf::UnSetCat()
{
  cat_current = cats.end();
}

VOCatResources VOCatConf::GetCat(const wxString& name) const
{
  vector<VOCatResources>::const_iterator c;
  for(c = cats.begin(); c != cats.end(); ++c)
    if( c->GetName() == name )
      return *c;

  return VOCatResources();
}


VOCatResources VOCatConf::GetCatFile(const wxString& name) const
{
  VOTable ct(name);
  if( ct.IsOk() ) {

    wxString catname = ct.GetDescription();

    vector<VOCatResources>::const_iterator c;
    for(c = cats.begin(); c != cats.end(); ++c)
      if( catname.Find(c->GetName()) != wxNOT_FOUND  )
	return *c;
  }

  return VOCatResources();
}

VOCatResources VOCatConf::GetCatFits(const wxString& filename) const
{
  fitsfile *f;
  int status = 0;
  int nhdu = 0;
  int hdutype;
  char extname[FLEN_CARD],comment[FLEN_CARD];

  wxString catname;

  status = 0;
  if( fits_open_file(&f, filename.fn_str(), READONLY, &status) )
    fits_report_error(stderr, status);

  if( fits_get_num_hdus(f,&nhdu,&status) )
    fits_report_error(stderr, status);

  for(int k = 0; k < nhdu; k++) {

    if( fits_movabs_hdu(f,k+1,&hdutype,&status) )
      fits_report_error(stderr, status);

    fits_read_keyword(f,"EXTNAME",extname,comment,&status);
    if( status == 0 && (hdutype == ASCII_TBL || hdutype == BINARY_TBL) ) {
      catname = extname;
      break;
    }
    status = 0;
  }
  fits_close_file(f, &status);

  if( status == 0 && ! catname.IsEmpty() ) {
    vector<VOCatResources>::const_iterator c;
    for(c = cats.begin(); c != cats.end(); ++c) {
      if( catname.Find(c->GetLabel("ID")) != wxNOT_FOUND  )
	return *c;
    }
  }
  else
    fits_report_error(stderr, status);

  return VOCatResources();
}
