/*

  VOTable parser


  Copyright © 2010-2026 F.Hroch (hroch@physics.muni.cz)

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

#include <wx/wx.h>
#include <wx/xml/xml.h>
#include <vector>

class VOField
{
  const wxXmlNode *node;

 public:

 VOField(): node(0) {}
 VOField(const wxXmlNode *nod): node(nod) {}

  wxString GetLabel() const;
  wxString GetType() const;
  wxString GetUnit() const;
  wxString GetArraySize() const;
};

class VOTableTable {
public:
 VOTableTable(): nrows(0), tablenode(0) {}
  VOTableTable(const wxXmlNode *node);
  int RecordCount() const;
  std::vector<wxString> GetRecord(int) const;
  std::vector<wxString> GetColumn(int) const;

  std::vector<wxString> GetFields() const;
  std::vector<wxString> GetUnits() const;
  std::vector<wxString> GetTypes() const;

  //private:

  wxString description;
  std::vector<VOField> fields;
  std::vector<wxString> links;

  int nrows;
  wxXmlNode *tablenode;
  std::vector<wxString> elements;
  std::vector<wxXmlNode *> rows;

  void VOData(const wxXmlNode *);
  void VOTableData(const wxXmlNode *);

};


class VOResource
{
  const wxXmlNode *resource;

public:
  VOResource(): resource(0) {}
  VOResource(const wxXmlNode *);

  wxString description;
  VOTableTable table;

};


class VOTable: public wxXmlDocument
{
public:
  VOTable();
  VOTable(const wxString&);
  VOTable(const wxXmlDocument&);
  bool LoadVO(const wxString&);


  wxString GetDescription() const;
  std::vector<VOResource> GetResources() const { return resources; }

  int RecordCount() const;
  bool IsEmpty() const;

  bool HasError() const { return ! infomsg.IsEmpty(); }
  wxString GetErrorMsg() const { return infomsg; }
  std::vector<wxString> GetQueryPar() const;

  bool Sort(const wxString&, const size_t index=0);

  bool SaveVO(const wxString&) const;
  bool SaveVO(wxOutputStream&) const;

  void ShowStructure() const;

 protected:

  wxString infomsg;
  std::vector<VOResource> resources;

 private:

  void ParseTable();
  void ParseErrorInfo();
  void ShowNode(const wxXmlNode *node) const;

};


class SVGcanvas: public VOTable
{
  wxString proj_type, mag_key, alpha_key, delta_key;

  double proj_alpha, proj_delta, proj_scale, mag_limit;
  long canvas_width, canvas_height;

  double ToDouble(const wxString& a) const;

 public:
  SVGcanvas(const wxString&);
  bool SaveSVG(wxOutputStream&) const;

  void SetProjection(const wxString&);
  void SetProjectionCenter(double,double);
  void SetCanvasSize(long,long);
  void SetScale(double);
  void SetMaglim(double);
  void SetMagkey(const wxString&);
  void SetAlphakey(const wxString&);
  void SetDeltakey(const wxString&);

};

class FITStable: public VOTable
{

  std::vector<wxString> FitsTypes(const std::vector<VOField>& fields) const;
  char **GetArray(const std::vector<wxString>&) const;
  wxString GetExtname() const;
  wxString GetQueryRequest() const;
  bool GetQueryInfo(wxString& , wxString&) const;
  wxString GetInfo(const wxXmlNode *, const wxString& ) const;
  bool GetEpoch(double&, wxString&) const;
  const wxXmlNode *GetResource() const;

 public:
  FITStable(const wxString&);
  FITStable(const VOTable&);
  bool SaveFITS(const wxString&, bool) const;

};

class TXTable: public VOTable
{

 public:
  TXTable(const wxString&);
  bool SaveTXT(wxOutputStream& output) const;

};

class CSVtable: public VOTable
{

 public:
  CSVtable(const wxString&);
  bool SaveCSV(wxOutputStream& output) const;

};
