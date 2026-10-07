/*

  VOTable to FITS convertor

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

*/

#include "votable.h"
#include "mfitsio.h"
#include <wx/wx.h>
#include <wx/txtstrm.h>
#include <wx/sstream.h>
#include <wx/tokenzr.h>
#include <iostream>
#include <fitsio.h>
#include <vector>

using namespace std;


FITStable::FITStable(const wxString& filename): VOTable(filename) {}

FITStable::FITStable(const VOTable& name): VOTable(name) {}

// In fact, to save of cone search results
bool FITStable::SaveFITS(const wxString& fitsname, bool clobber) const
{
  if( ! IsOk() ) return false;
  if( HasError() ) return false;

  fitsfile *file = 0;
  int status = 0;

  fits_create_file(&file,fitsname.char_str(),&status);
  if( status == FILE_NOT_CREATED && clobber ) {
    mfitsio_unlink(string(fitsname));
    status = 0;
    fits_create_file(&file,fitsname.char_str(),&status);
  }
  if( status ) {
    fits_report_error(stderr,status);
    return false;
  }

  // the dumb table
  fits_create_img(file,8,0,0,&status);

  wxString qvalue, constrain;
  if( GetQueryInfo(qvalue,constrain) ) {
    wxString com = "Cone search parameters: "+qvalue;
    const char *s = com.fn_str();
    fits_write_comment(file,(char *) s,&status);
    com = "Cone search constrains: "+constrain;
    const char *ss = com.fn_str();
    fits_write_comment(file,(char *) ss,&status);
  }

  fits_write_comment(file,"VO query:",&status);
  const wxString q(GetQueryRequest());
  fits_write_comment(file,q.fn_str(),&status);

  fits_write_comment(file,"VO server description:",&status);
  wxStringTokenizer tokenizer(GetDescription(),"\r\n");
  while ( tokenizer.HasMoreTokens() ) {
    wxString description = tokenizer.GetNextToken();
    fits_write_comment(file,description.fn_str(),&status);
  }

  for(size_t ir = 0; ir < resources.size(); ir++) {

    // tables
    const VOTableTable table(resources[ir].table);

    vector<wxString> fs = table.GetFields();
    vector<wxString> ft = table.GetTypes();
    vector<wxString> fu = table.GetUnits();

    vector<wxString> tf = FitsTypes(table.fields);
    vector<wxString> tu = fu;

    // VizieR catalogues reports magnitude errors (e_Vmag) in an obscure
    // type as characters (or integer shorts in past) like '02' for 0.02 (!?).
    // These values are converted to more convenient form.
    for(size_t i = 0; i < tf.size(); i++) {
      if( fu[i] == "cmag" ) {
	ft[i] = "float";
	tu[i] == "mag";
	tf[i] = "1E";
      }
    }

    char **ttype = GetArray(resources[ir].table.GetFields());
    char **tform = GetArray(tf);
    char **tunit = GetArray(tu);

    const wxString extname(GetExtname());
    fits_create_tbl(file,BINARY_TBL,0,fs.size(),ttype,tform,tunit,
		    extname.fn_str(),&status);

    for(size_t j = 0; j < fs.size(); j++) {
      delete[] ttype[j];
      delete[] tform[j];
      delete[] tunit[j];
    }
    delete[] ttype;
    delete[] tform;
    delete[] tunit;

    string rurl(GetInfo(GetResource()->GetChildren(),"reference_url"));
    fits_write_key(file,TSTRING,"REFURL",(char *) rurl.c_str(),
		   "Reference URL",&status);

    double epoch = 2000.0;
    wxString radesys;
    if( GetEpoch(epoch,radesys) ) {
      const char *s = radesys.fn_str();
      fits_write_key(file,TSTRING,"RADESYS",(char*)s,
		     "reference coordinate system",&status);
      fits_write_key(file,TDOUBLE,"JEPOCH",&epoch,"reference time (Julian epoch) ",
		     &status);
      // keywords must match ones in lib/titsio.f08
    }

    fits_write_comment(file,resources[ir].description.fn_str(),&status);

    for(size_t i = 0; i < fs.size() && status == 0 ; i++) {

      vector<wxString> c = table.GetColumn(i);
      size_t n = c.size();

      if( ft[i] == "boolean" ) {
	char *x = new char[n];
	for(size_t j = 0; j < n; j++) {
	  if( c[j]== "T" || c[j]== "t" || c[j]== "1" )
	    x[j] = 'T';
	  else if( c[j]== "F" || c[j]== "f" || c[j]== "0" )
	    x[j] = 'F';
	  else {
	    x[j] = ' ';
	    if( c[j] != "" ) {
	      wxLogError("Failed convert element `"+c[j]+"' at row=%d col="+
			 fs[i]+" to boolean.",(int)j);
	    }
	  }
	}
	fits_write_col(file,TLOGICAL,i+1,1,1,n,x,&status);
	delete[] x;
      }

      else if( ft[i] == "bit" ) {
	char *x = new char[n];
	for(size_t j = 0; j < n; j++) {
	  long t;
	  if( c[j].ToLong(&t) && (t == 0 || t == 1) )
	    x[j] = t;
	  else {
	    x[j] = 0;
	    if( c[j] != "" ) {
	      wxLogError("Failed convert element `"+c[j]+"' at row=%d col="+
			 fs[i]+" to bit(s).",(int)j);
	    }
	  }
	}
	fits_write_col(file,TBIT,i+1,1,1,n,x,&status);
	delete[] x;
      }

      else if( ft[i] == "unsignedByte" ) {
	unsigned char *x = new unsigned char[n];
	for(size_t j = 0; j < n; j++) {
	  long t;
	  if(  c[j].ToLong(&t) && (0 <= t && t < 256) )
	    x[j] = t;
	  else {
	    x[j] = 0;
	    if( c[j] != "" ) {
	      wxLogError("Failed convert element `"+c[j]+"' at row=%d col="+
			 fs[i]+" to uchar.",(int)j);
	    }
	  }
	}
	fits_write_col(file,TBYTE,i+1,1,1,n,x,&status);
	delete[] x;
      }

      else if( ft[i] == "short" ) {
	short *x = new short[n];
	for(size_t j = 0; j < n; j++) {
	  long t;
	  if( c[j].ToLong(&t) )
	    x[j] = t;
	  else {
	    x[j] = 0;
	    if( c[j] != "" ) {
	      wxLogError("Failed convert element `"+c[j]+"' at row=%d col="+
			 fs[i]+" to short.",(int)j);
	    }
	  }
	}
	fits_write_col(file,TSHORT,i+1,1,1,n,x,&status);
	delete[] x;
      }

      else if( ft[i] == "int" ) {
	int *x = new int[n];
	for(size_t j = 0; j < n; j++) {
	  long t;
	  if( c[j].ToLong(&t) )
	    x[j] = t;
	  else {
	    x[j] = 0;
	    if( c[j] != "" ) {
	      wxLogError("Failed convert element `"+c[j]+"' at row=%d col="+
			 fs[i]+" to int.",(int)j);
	    }
	  }
	}
	fits_write_col(file,TINT,i+1,1,1,n,x,&status);
	delete[] x;
      }

      else if( ft[i] == "long" ) {
	long *x = new long[n];
	for(size_t j = 0; j < n; j++) {
	  long t;
	  if( c[j].ToLong(&t) )
	    x[j] = t;
	  else {
	    x[j] = 0;
	    if( c[j] != "" ) {
	      wxLogError("Failed convert element `"+c[j]+"' at row=%d col="+
			 fs[i]+" to long.",(int)j);
	    }
	  }
	}
	fits_write_col(file,TLONG,i+1,1,1,n,x,&status);
	delete[] x;
      }

      else if( ft[i] == "char" || ft[i] == "unicodeChar" ) {
	char **x = GetArray(c);
	fits_write_col(file,TSTRING,i+1,1,1,n,x,&status);
	for(size_t j = 0; j < n; j++)
	  delete[] x[j];
	delete[] x;
      }

      else if( ft[i] == "float" ) {
	float scale = fu[i] == "cmag" ? 0.01 : 1.0;
	float *x = new float[n];
	for(size_t j = 0; j < n; j++) {
	  double t;
	  if( c[j].ToDouble(&t) )
	    x[j] = scale*t;
	  else {
	    x[j] = nan("");
	    if( c[j] != "" ) {
	      wxLogError("Failed convert element `"+c[j]+"' at row=%d col="+
			 fs[i]+" to float.",(int)j);
	    }
	  }
	}
	fits_write_col(file,TFLOAT,i+1,1,1,n,x,&status);
	delete[] x;
      }

      else if( ft[i] == "double" ) {
	double scale = fu[i] == "cmag" ? 0.01 : 1.0;
	double *x = new double[n];
	for(size_t j = 0; j < n; j++) {
	  double t;
	  if( c[j].ToDouble(&t) )
	    x[j] = scale*t;
	  else {
	    x[j] = nan("");
	    if( c[j] != "" ) {
	      wxLogError("Failed convert element `"+c[j]+"' at row=%d col="+
			 fs[i]+" to double.",(int)j);
	    }
	  }
	}
	fits_write_col(file,TDOUBLE,i+1,1,1,n,x,&status);

	delete[] x;
      }

      else if( ft[i] == "floatComplex" ) {
	float *x = new float[2*n];
	for(size_t j = 0; j < n; j++) {
	  wxStringInputStream is(c[j]);
	  wxTextInputStream ts(is);
	  ts >> x[2*j] >> x[2*j+1];
	}
	fits_write_col(file,TCOMPLEX,i+1,1,1,2*n,x,&status);
	delete[] x;
      }
      else if( ft[i] == "doubleComplex" ) {
	double *x = new double[2*n];
	for(size_t j = 0; j < n; j++) {
	  wxStringInputStream is(c[j]);
	  wxTextInputStream ts(is);
	  ts >> x[2*j] >> x[2*j+1];
	}
	fits_write_col(file,TDBLCOMPLEX,i+1,1,1,2*n,x,&status);
	delete[] x;
      }

    }
  }

  fits_close_file(file, &status);
  fits_report_error(stderr,status);

  return status == 0;
}

vector<wxString> FITStable::FitsTypes(const vector<VOField>& fields) const
{
  vector<wxString> tform;
  for(vector<VOField>::const_iterator f = fields.begin(); f !=fields.end();++f){

    wxString datatype = f->GetType();
    wxString arraysize = f->GetArraySize();

    if( datatype == "char" ) {
      if( ! arraysize.empty() ) {
	arraysize.Replace("*","");
	long n;
	if( arraysize.ToLong(&n) ) {
	  wxString a;
	  a.Printf("%d",int(n));
	  tform.push_back(a+"A");
	}
      }
      else
	tform.push_back("1A");
    }
    else if( datatype == "double" )
      tform.push_back("1D");
    else if( datatype == "boolean" )
      tform.push_back("1L");
    else if( datatype == "bit" )
      tform.push_back("1X");
    else if( datatype == "unsignedByte" )
      tform.push_back("1B");
    else if( datatype == "short" )
      tform.push_back("1I");
    else if( datatype == "int" || datatype == "long" )
      tform.push_back("1J");
    else if( datatype == "longlong" )
      tform.push_back("1K");
    else if( datatype == "float" )
      tform.push_back("1E");
    else if( datatype == "floatComplex" )
      tform.push_back("1C");
    else if( datatype == "doubleComplex" )
      tform.push_back("1M");
    else if( datatype == "unicodeChar" )
      tform.push_back("");
    else
      wxLogFatalError("Unsupported data type");

  }
  return tform;
}

char **FITStable::GetArray(const vector<wxString>& vec) const
{
  int n = vec.size();
  char **x = new char*[n];

  for(int j = 0; j < n; j++) {
    size_t l = vec[j].Len();
    char *a = new char[l + 1];
    if( l > 0 ) {
      for(size_t i = 0; i < l; i++)
	a[i] = vec[j].GetChar(i);
      a[l] = '\0';
    }
    else
      a[0] = '\0';
    x[j] = a;
  }
  return x;
}

wxString FITStable::GetQueryRequest() const
{
  const wxXmlNode *node = GetRoot()->GetChildren();
  while (node) {
    if( node->GetName() == "INFO" ) {
      wxString a(node->GetAttribute("name"));
      if( a == "request" )
	return node->GetAttribute("value");
    }
    node = node->GetNext();
  }
  return "";
}

const wxXmlNode *FITStable::GetResource() const
{
  const wxXmlNode *node = GetRoot()->GetChildren();
  while (node) {
    if( node->GetName() == "RESOURCE" )
      return node;
    node = node->GetNext();
  }
  return 0;
}

wxString FITStable::GetExtname() const
{
  const wxXmlNode *resource = GetResource();
  return resource->GetAttribute("name");
}


wxString FITStable::GetInfo(const wxXmlNode *node, const wxString& name) const
{
  while (node) {
    if( node->GetName() == "INFO" ) {
      if( node->GetAttribute("name") == name )
	return node->GetAttribute("value");
    }
    node = node->GetNext();
  }
  return "";
}


bool FITStable::GetEpoch(double& epoch, wxString& system) const
{
  const wxXmlNode *node = GetResource()->GetChildren();
  while (node) {
    if( node->GetName() == "COOSYS" ) {
      wxString a(node->GetAttribute("epoch"));
      if( a != "" ) {
	system = node->GetAttribute("system");
	if( a.ToDouble(&epoch) )
	  return true;
      }
    }
    node = node->GetNext();
  }
  return false;
}

bool FITStable::GetQueryInfo(wxString& value, wxString& constrain) const
{
  const wxXmlNode *node = GetRoot()->GetChildren();
  while (node) {
    if( node->GetName() == "INFO" && node->GetAttribute("ID") == "Target" ) {
      value = node->GetAttribute("value");
      constrain = node->GetNodeContent();
      return true;
    }
    node = node->GetNext();
  }
  return false;
}
