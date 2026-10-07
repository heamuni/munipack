/*

  VOTable parser

  Copyright © 2010-2015, 2017-2026 F.Hroch (hroch@physics.muni.cz)

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

  This is very simple VOTable parser. It is focused only
  on reading and interpretation of data from a VO server:

  * One TABLE in RESOURCE
  * Grouping anything is ignored.
  * Recursive use of RESOURCE is avoided.
  * Many parameters is empirically determinted (like error processig).

  I've doubts about determination of coordinate frame. Vizier offers
  a different structure for UCAC4 and UCAC5 catalogues. It means that
  determination of epoch of catalogue coordinates cann't be assumed
  as a simple parameter. Does it look as the maintainers proposed
  another way of processing?

*/

#include "common.h"
#include "votable.h"
#include <wx/wx.h>
#include <wx/xml/xml.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>
#include <wx/sstream.h>
#include <wx/file.h>
#include <wx/filename.h>
#include <wx/filefn.h>
#include <wx/tokenzr.h>
#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>

using namespace std;


wxString VOField::GetLabel() const
{
  wxASSERT(node);
  return node->GetAttribute("name");
}

wxString VOField::GetType() const
{
  wxASSERT(node);
  return node->GetAttribute("datatype");
}

wxString VOField::GetUnit() const
{
  wxASSERT(node);
  return node->GetAttribute("unit");
}

wxString VOField::GetArraySize() const
{
  wxASSERT(node);
  return node->GetAttribute("arraysize");
}


void VOTableTable::VOTableData(const wxXmlNode *root)
{
  wxASSERT(root);

  wxXmlNode *tr = root->GetChildren();
  while (tr) {

    if( tr->GetName() == "TR" ) {

      nrows++;
      rows.push_back(tr);

      wxXmlNode *td = tr->GetChildren();
      while (td) {
	elements.push_back(td->GetNodeContent());
	td = td->GetNext();
      }
    }
    tr = tr->GetNext();
  }
}


void VOTableTable::VOData(const wxXmlNode *root)
{
  wxASSERT(root);

  wxXmlNode *node = root->GetChildren();
  while (node) {

    if( node->GetName() == "TABLEDATA" ) {
      tablenode = node;
      VOTableData(node);
    }
    else
      wxLogDebug("Only TABLEDATA element is supported.");

    node = node->GetNext();
  }
}


VOTableTable::VOTableTable(const wxXmlNode *root): nrows(0)
{
  wxASSERT(root);

  const wxXmlNode *node = root->GetChildren();
  while (node) {

    //    wxLogDebug(node->GetName());

    if( node->GetName() == "DESCRIPTION" )
      description = node->GetNodeContent();

    else if( node->GetName() == "FIELD" )
      fields.push_back(VOField(node));

    else if( node->GetName() == "LINK" )
      links.push_back(node->GetNodeContent());

    else if( node->GetName() == "DATA" )
      VOData(node);

    node = node->GetNext();
  }
}

int VOTableTable::RecordCount() const
{
  return nrows;
}

vector<wxString> VOTableTable::GetRecord(int n) const
{
  vector<wxString> r;
  int ncol = fields.size();

  for(int i = 0; i < ncol; i++)
    r.push_back(elements[n*ncol+i]);

  return r;
}

vector<wxString> VOTableTable::GetColumn(int n) const
{
  int ncol = fields.size();

  //  wxLogDebug(_("%d ")+data.table.elements[1],ncol);

  vector<wxString> c;
  for(int i = 0; i < nrows; i++)
    c.push_back(elements[i*ncol+n]);

  return c;
}

vector<wxString> VOTableTable::GetFields() const
{
  vector<wxString> fl;
  for(size_t i = 0; i < fields.size(); i++)
    fl.push_back(fields[i].GetLabel());
  return fl;
}

vector<wxString> VOTableTable::GetTypes() const
{
  vector<wxString> fl;
  for(size_t i = 0; i < fields.size(); i++)
    fl.push_back(fields[i].GetType());
  return fl;
}

vector<wxString> VOTableTable::GetUnits() const
{
  vector<wxString> fl;
  for(size_t i = 0; i < fields.size(); i++)
    fl.push_back(fields[i].GetUnit());
  return fl;
}


VOResource::VOResource(const wxXmlNode *root): resource(root)
{
  wxASSERT(root);

  const wxXmlNode *node = root->GetChildren();
  while (node) {

    //    wxLogDebug(node->GetName());

    if( node->GetName() == "DESCRIPTION" )
      description = node->GetNodeContent();

    else if( node->GetName() == "TABLE" )
      table = VOTableTable(node);

    node = node->GetNext();
  }
}


VOTable::VOTable(): wxXmlDocument() {}

VOTable::VOTable(const wxString& filename): wxXmlDocument(filename)
{
  if( ! IsOk() )
    return;

  if( GetRoot()->GetName() != "VOTABLE" ) {
    wxLogError("The document does not looks as VOTable.");
    return;
  }

  ParseErrorInfo();
  ParseTable();
}


VOTable::VOTable(const wxXmlDocument& vt): wxXmlDocument(vt)
{
  if( ! IsOk() )
    return;

  if( GetRoot()->GetName() != "VOTABLE" ) {
    wxLogError("The document does not looks as VOTable.");
    return;
  }

  ParseErrorInfo();
  ParseTable();
}

bool VOTable::LoadVO(const wxString& filename)
{
  wxLogDebug("Loading data ("+filename+") ...");
  wxXmlDocument::Load(filename,"UTF-8", wxXMLDOC_KEEP_WHITESPACE_NODES);

  if( ! IsOk() )
    return false;

  if( GetRoot()->GetName() != "VOTABLE" ) {
    wxLogError("The document does not looks as VOTable.");
    return false;
  }

  ParseErrorInfo();
  ParseTable();
  return IsOk();
}

void VOTable::ParseTable()
{
  if( ! IsOk() ) return;
  if( HasError() ) return;

  wxLogDebug("Parsing data ...");

  resources.clear();

  wxASSERT(IsOk() && GetRoot()->GetName() == "VOTABLE");

  const wxXmlNode *node = GetRoot()->GetChildren();
  while (node) {

    if( node->GetName() == "RESOURCE" )
      resources.push_back(VOResource(node));

    node = node->GetNext();
  }
}


bool VOTable::IsEmpty() const
{
  return resources.size() == 0 ||
    ( resources.size() > 0 && resources[0].table.RecordCount() == 0);
}


int VOTable::RecordCount() const
{
  if( resources.size() > 0 )
    return resources[0].table.RecordCount();
  else
    return 0;
}

bool VOTable::SaveVO(const wxString& filename) const
{
  if( filename.IsEmpty() ) {
    wxFFile ffile;
    ffile.Attach(stdout);
    wxFFileOutputStream ostream(ffile);
    return Save(ostream);
  }

  return wxXmlDocument::Save(filename, wxXML_NO_INDENTATION);
}

bool VOTable::SaveVO(wxOutputStream& ostream) const
{
  return wxXmlDocument::Save(ostream, wxXML_NO_INDENTATION);
}


void VOTable::ShowStructure() const
{
  // Diagnostic routine printing structure of VOTable

  wxASSERT(IsOk());

  wxLogDebug(GetRoot()->GetName());
  ShowNode(GetRoot()->GetChildren());
}

void VOTable::ShowNode(const wxXmlNode *node) const
{
  wxASSERT(IsOk());

  if( ! node ) return;
  while (node) {
    if( ! (node->GetName() == "TR" || node->GetName() == "TD" ||
	   node->GetType() == wxXML_TEXT_NODE ||
	   node->GetType() == wxXML_COMMENT_NODE  ) ) {
      wxString vata('\t',node->GetDepth(GetRoot())+1);
      if( node->GetName() == "DESCRIPTION" )
	wxLogDebug(vata+node->GetName()+"(%d) "+node->GetNodeContent(),
		  node->GetType());
      else
	wxLogDebug(vata+node->GetName()+"(%d)",node->GetType());
    }
    ShowNode(node->GetChildren());
    node = node->GetNext();
  }

}

wxString VOTable::GetDescription() const
{
  wxASSERT(IsOk());

  wxString description;

  const wxXmlNode *node = GetRoot()->GetChildren();
  while( node ) {
    if( node->GetName() == "DESCRIPTION" ) {
      description = node->GetNodeContent();
      return description;
    }
    node = node->GetNext();
  }
  return description;
}



void VOTable::ParseErrorInfo()
{
  wxASSERT(IsOk());

  infomsg = "";
  const wxXmlNode *node = GetRoot()->GetChildren();
  while( node ) {

    if( node->GetName() == "INFO" ) {
      if( node->GetAttribute("name") == "Error" ) {
	infomsg = node->GetAttribute("value");
	return;
      }
    }
    node = node->GetNext();
  }
}

vector<wxString> VOTable::GetQueryPar() const
{
  wxASSERT(IsOk());

  vector<wxString> pars;

  const wxXmlNode *node = GetRoot()->GetChildren();
  while( node ) {

    if( node->GetName() == "INFO" ) {
      if( node->GetAttribute("name") == "queryParameters" ) {
	wxStringTokenizer tokenizer(node->GetNodeContent(),"\r\n");
	while ( tokenizer.HasMoreTokens() )
	  pars.push_back(tokenizer.GetNextToken());
	return pars;
      }
    }
    node = node->GetNext();
  }
  return pars;
}




// needs to add other types
bool cmp(pair<int,double> a, pair<int,double> b)
{
  //  wxLogDebug("%f %f",b.second,a.second);
  return b.second > a.second;
}

bool VOTable::Sort(const wxString& key, const size_t index)
{
  wxLogDebug("Sorting by "+key+" ... ");

  if( ! ( resources.size() > index ) ) {
    wxLogDebug("Sorting skipped (no data).");
    return true;
  }

  const VOTableTable table(resources[index].table);

  if( ! (table.RecordCount() > 0 ) ) {
    wxLogDebug("Sorting skipped (empty table).");
    return true;
  }

  int n = -1;
  vector<wxString> fs = table.GetFields();
  for(vector<wxString>::const_iterator i = fs.begin(); i != fs.end(); ++i) {
    if( *i == key ) {
      n = i - fs.begin();
      break;
    }
  }

  if( n == -1 )
    return false;

  vector< pair<int,double> > d;

  vector<wxString> col = table.GetColumn(n);

  for(vector<wxString>::const_iterator i = col.begin(); i != col.end(); ++i) {
    wxString a(*i);
    double x;
    if( ! a.IsEmpty() && a.ToDouble(&x) )
      d.push_back(make_pair(i-col.begin(),x));
    else
      // undefined values, works correctly only for magnitudes
      d.push_back(make_pair(i-col.begin(),99.999));;
  }

  sort(d.begin(),d.end(),cmp);

  /*
  for(vector< pair<int,double> >::const_iterator i = d.begin(); i!=d.end();++i)
    wxLogDebug("%d %f",i->first,i->second);
  */

  wxXmlNode *t = table.tablenode;
  wxXmlNode *nt = new wxXmlNode(t->GetParent(),t->GetType(),
				t->GetName(),t->GetContent(),
				t->GetAttributes(),t->GetNext(),
				t->GetLineNumber());

  vector<wxXmlNode *> rows(table.rows);
  for(vector< pair<int,double> >::const_iterator i=d.begin(); i!=d.end(); ++i){
    nt->AddChild(rows[i->first]);
  }

  wxXmlNode *parent = t->GetParent();
  parent->RemoveChild(t);

  ParseTable();
  return true;
}
