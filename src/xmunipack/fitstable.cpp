/*

  xmunipack - fits tables

  Copyright © 2021-22 F.Hroch (hroch@physics.muni.cz)

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


#include "fits.h"
#include <cfloat>
#include <wx/wx.h>
#include <wx/regex.h>
#include <algorithm>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#include <wx/txtstrm.h>
#include <wx/wfstream.h>
#endif


using namespace std;

// ?? move to fits.cpp ???
class FitsTableColumnData: public wxObjectRefData
{
public:
  FitsTableColumnData();
  FitsTableColumnData(int, long, float *);
  FitsTableColumnData(int, long, int *);
  FitsTableColumnData(int, long, char **);
  FitsTableColumnData(int, long, char *);
  FitsTableColumnData(int, long, double *);
  FitsTableColumnData(int, long, short *);
  FitsTableColumnData(int, long, long *);
  FitsTableColumnData(int, long, bool *);
  FitsTableColumnData(const FitsTableColumnData&);
  FitsTableColumnData& operator = (const FitsTableColumnData&);
  virtual ~FitsTableColumnData();

  bool operator == (const FitsTableColumnData&) const;

  int typecode;
  long nrows;
  bool *otable;
  char *btable;
  short *stable;
  int *itable;
  long *ltable;
  float *ftable;
  double *dtable;
  char **ctable;

};

class FitsTableData: public wxObjectRefData
{
public:
  FitsTableData();
  FitsTableData(const std::vector<FitsTableColumn>&);
  bool operator == (const FitsTableData&) const;

  std::vector<FitsTableColumn> columns;

};




FitsTableColumnData::FitsTableColumnData():
  typecode(0),nrows(0),otable(0),btable(0),stable(0),itable(0),
  ltable(0),ftable(0),dtable(0),ctable(0) {}

FitsTableColumnData::FitsTableColumnData(int t, long n, float *d):
  typecode(t),nrows(n),otable(0),btable(0),stable(0),itable(0),
  ltable(0),ftable(d),dtable(0),ctable(0) {}

FitsTableColumnData::FitsTableColumnData(int t, long n, int *d):
  typecode(t),nrows(n),otable(0),btable(0),stable(0),itable(d),
  ltable(0),ftable(0),dtable(0),ctable(0) {}

FitsTableColumnData::FitsTableColumnData(int t, long n, char **d):
  typecode(t),nrows(n),otable(0),btable(0),stable(0),itable(0),
  ltable(0),ftable(0),dtable(0),ctable(d) {}

FitsTableColumnData::FitsTableColumnData(int t, long n, char *d):
  typecode(t),nrows(n),otable(0),btable(d),stable(0),itable(0),
  ltable(0),ftable(0),dtable(0),ctable(0) {}

FitsTableColumnData::FitsTableColumnData(int t, long n, bool *d):
  typecode(t),nrows(n),otable(d),btable(0),stable(0),itable(0),
  ltable(0),ftable(0),dtable(0),ctable(0) {}

FitsTableColumnData::FitsTableColumnData(int t, long n, short *d):
  typecode(t),nrows(n),otable(0),btable(0),stable(d),itable(0),
  ltable(0),ftable(0),dtable(0),ctable(0) {}

FitsTableColumnData::FitsTableColumnData(int t, long n, long *d):
  typecode(t),nrows(n),otable(0),btable(0),stable(0),itable(0),
  ltable(d),ftable(0),dtable(0),ctable(0) {}

FitsTableColumnData::FitsTableColumnData(int t, long n, double *d):
  typecode(t),nrows(n),otable(0),btable(0),stable(0),itable(0),
  ltable(0),ftable(0),dtable(d),ctable(0) {}

FitsTableColumnData::~FitsTableColumnData()
{
  delete[] ftable;
  delete[] itable;
  delete[] dtable;
  delete[] ltable;
  delete[] stable;
  delete[] btable;
  delete[] otable;

  if( ctable ) {
    for(int i = 0; i < nrows; i++)
      delete[] ctable[i];
    delete[] ctable;
  }
}


FitsTableColumnData::FitsTableColumnData(const FitsTableColumnData& copy)
{
  wxFAIL_MSG("FitsTableColumn ---- WE ARE REALY NEED COPY CONSTRUCTOR ----");

}

FitsTableColumnData& FitsTableColumnData::operator = (const FitsTableColumnData& other)
{
 wxFAIL_MSG("** FitsTableColumn: WE ARE REALLY NEED ASSIGNMENT CONSTRUCTOR");
 return *this;
}


FitsTableColumn::FitsTableColumn() {}

FitsTableColumn::~FitsTableColumn() {}

FitsTableColumn::FitsTableColumn(int t, long nr, float *d)
{
  SetRefData(new FitsTableColumnData(t,nr,d));
}

FitsTableColumn::FitsTableColumn(int t, long nr, double *d)
{
  SetRefData(new FitsTableColumnData(t,nr,d));
}

FitsTableColumn::FitsTableColumn(int t, long nr, long *d)
{
  SetRefData(new FitsTableColumnData(t,nr,d));
}

FitsTableColumn::FitsTableColumn(int t, long nr, int *d)
{
  SetRefData(new FitsTableColumnData(t,nr,d));
}

FitsTableColumn::FitsTableColumn(int t, long nr, short *d)
{
  SetRefData(new FitsTableColumnData(t,nr,d));
}

FitsTableColumn::FitsTableColumn(int t, long nr, char *d)
{
  SetRefData(new FitsTableColumnData(t,nr,d));
}

FitsTableColumn::FitsTableColumn(int t, long nr, bool *d)
{
  SetRefData(new FitsTableColumnData(t,nr,d));
}

FitsTableColumn::FitsTableColumn(int t, long nr, char **d)
{
  SetRefData(new FitsTableColumnData(t,nr,d));
}

FitsTableColumn::FitsTableColumn(int typecode, long repeat, long width, long nrows)
{
  if( typecode == TSTRING ) {
    char **a = new char*[nrows];
    for(int i = 0; i < nrows; i++)
      a[i] = new char[width];
    SetRefData(new FitsTableColumnData(typecode,nrows,a));
  }
  else if( typecode == TLOGICAL ) {
    bool *b = new bool[nrows];
    SetRefData(new FitsTableColumnData(typecode,nrows,b));
  }
  else if( typecode == TBYTE || typecode == TBIT ) {
    char *a = new char[nrows];
    SetRefData(new FitsTableColumnData(typecode,nrows,a));
  }
  else if( typecode == TSHORT ) {
    short *d = new short[nrows];
    SetRefData(new FitsTableColumnData(typecode,nrows,d));
  }
  else if( typecode == TLONG || typecode == TINT32BIT ) {
    long *d = new long[nrows];
    SetRefData(new FitsTableColumnData(typecode,nrows,d));
  }
  else if( typecode == TFLOAT ) {
    float *d = new float[nrows];
    SetRefData(new FitsTableColumnData(typecode,nrows,d));
  }
  else if( typecode == TDOUBLE ) {
    double *d = new double[nrows];
    SetRefData(new FitsTableColumnData(typecode,nrows,d));
  }
}

bool FitsTableColumnData::operator == (const FitsTableColumnData& other) const
{
  return true; // !!!!!!!!!!!!!!!!! IMPLEMENT !!!!!!!!!!!!!!!1111
}

/*
FitsTableColumn::FitsTableColumn(const FitsTableColumn& copy)
{
}

FitsTableColumn& FitsTableColumn::operator = (const FitsTableColumn& other)
{
  wxFAIL_MSG("** FitsTableData: WE ARE REALLY NEED ASSIGNMENT CONSTRUCTOR");
  return *this;
}
*/

int FitsTableColumn::GetColType() const
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->typecode;
};

long FitsTableColumn::Nrows() const
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->nrows;
};


const float *FitsTableColumn::GetCol_float() const
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->ftable;
}

const double *FitsTableColumn::GetCol_double() const
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->dtable;
}

const long *FitsTableColumn::GetCol_long() const
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->ltable;
}

const char **FitsTableColumn::GetCol_char() const
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return (const char **) data->ctable;
}


wxObjectRefData *FitsTableColumn::CreateRefData() const
{
  return new FitsTableColumnData;
}


wxObjectRefData *FitsTableColumn::CloneRefData(const wxObjectRefData *that) const
{
  const FitsTableColumnData *olddata = static_cast<const FitsTableColumnData *>(that);
  wxASSERT(olddata);

  /*
  FitsTableData *newdata = new FitsTableData;
  newdata->ncols = olddata->ncols;
  newdata->nrows = olddata->nrows;
  long nelem = newdata->ncols*newdata->nrows;
  newdata->table = new float[nelem];
  memcpy(newdata->table,olddata->table,nelem*sizeof(float));
  return newdata;
  */

  FitsTableColumnData *newdata = new FitsTableColumnData;
  newdata->typecode = olddata->typecode;
  newdata->nrows = olddata->nrows;
  long nr = olddata->nrows;

  if( olddata->itable ) {
    newdata->itable = new int[nr];
    memcpy(newdata->itable,olddata->itable,nr*sizeof(int));
  }
  else if( olddata->ltable ) {
    newdata->ltable = new long[nr];
    memcpy(newdata->ltable,olddata->ltable,nr*sizeof(long));
  }
  else if( olddata->otable ) {
    newdata->otable = new bool[nr];
    memcpy(newdata->otable,olddata->otable,nr*sizeof(bool));
  }
  else if( olddata->btable ) {
    newdata->btable = new char[nr];
    memcpy(newdata->btable,olddata->btable,nr*sizeof(char));
  }
  else if( olddata->stable ) {
    newdata->stable = new short[nr];
    memcpy(newdata->stable,olddata->stable,nr*sizeof(short));
  }
  else if( olddata->ftable ) {
    newdata->ftable = new float[nr];
    memcpy(newdata->ftable,olddata->ftable,nr*sizeof(float));
  }
  else if( olddata->dtable ) {
    newdata->dtable = new double[nr];
    memcpy(newdata->dtable,olddata->dtable,nr*sizeof(double));
  }
  else if( olddata->ctable ) {
    newdata->ctable = new char*[nr];
    for(int i = 0; i < olddata->nrows; i++)
      newdata->ctable[i] = wxStrdup(olddata->ctable[i]);
  }
  else
    wxFAIL_MSG("FitsTableColumn::CloneRefData: Unknown data type.");

  return newdata;
}

FitsTableColumn FitsTableColumn::Copy() const
{
  FitsTableColumn new_col;
  new_col.m_refData = CloneRefData(m_refData);
  return new_col;
}


char **FitsTableColumn::GetData_string()
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->ctable;
}

char *FitsTableColumn::GetData_char()
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->btable;
}

bool *FitsTableColumn::GetData_bool()
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->otable;
}

short *FitsTableColumn::GetData_short()
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->stable;
}

int *FitsTableColumn::GetData_int()
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->itable;
}

long *FitsTableColumn::GetData_long()
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->ltable;
}

float *FitsTableColumn::GetData_float()
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->ftable;
}

double *FitsTableColumn::GetData_double()
{
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data);
  return data->dtable;
}




FitsTableData::FitsTableData() {}
FitsTableData::FitsTableData(const std::vector<FitsTableColumn>& cols): columns(cols) {}

bool FitsTableData::operator == (const FitsTableData& other) const
{
  if( columns.size() != other.columns.size() )
    return false;

  /*
  for(size_t i = 0; i < columns.size(); i++) {
    if( ! (columns[i] == other.columns[i]) )
      return false;
  }
  */

  return true;
}


// FitsTableData::FitsTableData(): ncols(0),columns(0) {}

// /*
// FitsTableData::FitsTableData(long nr, long nc, float *t):
//   nrows(nr),ncols(nc),table(t) {}
// */

// FitsTableData::FitsTableData(long nc): ncols(nc)
// {
//   columns = new FitsTableColumn[nc];
//   for(int i = 0; i < nc; i++)
//     columns[i] = 0;
// }

// FitsTableData::FitsTableData(const FitsTableData& copy)
// {
//   wxFAIL_MSG("FitsTableData ---- WE ARE REALY NEED COPY CONSTRUCTOR ----");
// }

// FitsTableData& FitsTableData::operator = (const FitsTableData& other)
// {
//  wxFAIL_MSG("** FitsTableData: WE ARE REALLY NEED ASSIGNMENT CONSTRUCTOR");
//  return *this;
// }


// FitsTableData::~FitsTableData()
// {
//   delete[] columns;

//   /*
//   for(int i = 0; i < ncols; i++) {
//     if( typecode[i] == TFLOAT ) {
//       float *d = (float *) table[i];
//       delete[] d;
//     }
//     else if( typecode[i] == TLONG ) {
//       long *d = (long *) table[i];
//       delete[] d;
//     }
//     else if( typecode[i] == TSTRING ) {
//       char *a = (char *) table[i];
//       delete[] a;
//     }
//   */

//     //    delete[] table[i];
//   //  }
//   //  delete[] typecode;
// }

/*
void FitsTableData::InsertColumn(long k, long nr, float *d)
{
  wxASSERT(0 <= k && k < ncols);
  columns[k] = FitsDataColumn(nr,d);
}

void FitsTableData::InsertColumn(long k, long nr, int *d)
{
  wxASSERT(0 <= k && k < ncols);
  columns[k] = FitsDataColumn(nr,d);
}

void FitsTableData::InsertColumn(long k, long nr, char **d)
{
  wxASSERT(0 <= k && k < ncols);
  columns[k] = FitsDataColumn(nr,d);
}
*/


// ----------    FitsTable


FitsTable::FitsTable(const FitsHeader& h, int ht,
		     const vector<FitsTableColumn>& cols):
  FitsHdu(h,ht)
{
  SetRefData(new FitsTableData(cols));
  //  wxLogDebug("%d %d",(int)Nrows(),(int)Ncols());
}


FitsTable::FitsTable(const FitsHdu& h):
  FitsHdu(h.GetHeader(),h.GetHduType())
{
  wxASSERT(h.Type() == HDU_TABLE && h.IsOk());
  Ref(h);

  //  wxLogDebug("FitsTable::FitsTable(const FitsHdu& h): ");

  //  wxLogDebug("%d %d",(int)Nrows(),(int)Ncols());

  //  fits_type = (const FitsTable) h.fits_type;

  /*
  fits_type = static_cast<FitsTable>(h).fits_type;
  columns = static_cast<FitsTable>(h).columns;
  */

  /*
  const FitsTable *t = static_cast<const FitsTable *>(&h);
  wxASSERT(t && t->IsOk());
  columns = t->columns;
  */

  /*
  FitsTableData *data = static_cast<FitsTableData *>(h.GetRefData());
  if( data )
    SetRefData(data);
  */
}

FitsTable::FitsTable(const FitsTable& t): FitsHdu(t)
{
  wxASSERT(t.IsOk());
  Ref(t);
}


FitsTable& FitsTable::operator = (const FitsTable& other)
{
  //  wxFAIL_MSG("*** FitsTable: WE ARE REALLY NEED ASSIGNMENT CONSTRUCTOR");
  if( this != & other ) {

    header = other.header;
    htype = other.htype;
    type = other.type;
    flavour = other.flavour;
    modified = other.modified;
    Ref(other);
  }
  return *this;
}

bool FitsTable::operator == (const FitsTable& other) const
{
  if (m_refData == other.m_refData)
    return true;

  if (!m_refData || !other.m_refData)
    return false;

  return ( *(FitsTableData*)m_refData == *(FitsTableData*)other.m_refData );
}

bool FitsTable::operator != (const FitsTable& other) const
{
  return  !(*this == other);
}

wxObjectRefData *FitsTable::CreateRefData() const
{
  return new FitsTableData;
}


wxObjectRefData *FitsTable::CloneRefData(const wxObjectRefData *that) const
{
  const FitsTableData *olddata = static_cast<const FitsTableData *>(that);
  wxASSERT(olddata);

  FitsTableData *newdata = new FitsTableData;
  newdata->columns = olddata->columns;

  return newdata;
}

//   /*
//   FitsTableData *newdata = new FitsTableData;
//   newdata->ncols = olddata->ncols;
//   newdata->nrows = olddata->nrows;
//   long nelem = newdata->ncols*newdata->nrows;
//   newdata->table = new float[nelem];
//   memcpy(newdata->table,olddata->table,nelem*sizeof(float));
//   return newdata;
//   */

//   FitsTableData *newdata = new FitsTableData;
//   newdata->ncols = olddata->ncols;
//   newdata->columns = new FitsDataColumn[newdata->ncols];
//   for(int i = 0; i < olddata->ncols; i++) {
//     newdata->columns[i] = olddata->columns[i];



//   /*
//   int len = 1;
//   for(int i = 0; i < olddata->ncols; i++) {
//     if( olddata->typecode[i] == TFLOAT )
//       len = len*sizeof(float);
//     else if( olddata->typecode[i] == TLONG )
//       len = len*sizeof(long);
//     else if( olddata->typecode[i] == TINT )
//       len = len*sizeof(int);
//     else if( olddata->typecode[i] == TLOGICAL )
//       len = len*sizeof(bool);
//     else
//       wxLogDebug("FitsTable::CloneRefData: Type `%d' not recognized.",olddata->typecode[i]);
//   }

//   long nelem = len*olddata->nrows;

//   newdata->typecode = new int[olddata->ncols];
//   memcpy(newdata->typecode,olddata->typecode,olddata->ncols*sizeof(int));

//   void *table = new void*[olddata->ncols];
//   for(int i = 0; i < olddata->ncols; i++) {

//     if( newdata->typecode[i] == TFLOAT ) {
//       table[i] = (void *) new float[olddata->nrows];
//       memcpy(table,olddata->table[i],olddata->nrows*sizeof(float));
//     }



//   }


//   newdata->table = new void**[];

//   memcpy(newdata->table,olddata->table,nelem);
//   return newdata;
//   */
// }


FitsTable FitsTable::Copy() const
{
  FitsTable new_table;
  new_table.m_refData = CloneRefData(m_refData);
  return new_table;
}



long FitsTable::Nrows() const
{
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  if( ! data )
    return 0;

  const vector<FitsTableColumn> columns(data->columns);
  long nrows = 0;

  wxASSERT(columns.size() > 0 );

  /*
  for(size_t i = 0; i < columns.size(); i++)
    nrows = columns[i].Nrows();
  */

  if( columns.size() > 0 )
    nrows = columns[0].Nrows();

  //    wxLogDebug("nrows = %d",(int)nrows);

  return nrows;

  /*
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  if( data )
    return data->nrows;
  else
    return 0;
  */
}


int FitsTable::Ncols() const
{
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  const vector<FitsTableColumn> columns(data->columns);
  return columns.size();
  /*
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  if( data )
    return data->ncols;
  else
    return 0;
  */
}


long FitsTable::Naxes(int n) const
{
  //  wxLogDebug("Don't use Naxes for tables!");
  /*
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  wxASSERT(data);
  */
  switch(n) {
  case 0:  return Ncols(); break;
  case 1:  return Nrows(); break;
  default: return 0;           break;
  }
}


bool FitsTable::IsOk() const
{
  /*
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  if( ! data ) return false;
  return m_refData && data->nrows > 0 && data->ncols > 0;
  */
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  bool b = data != 0;
  if( data ) {
    const vector<FitsTableColumn> columns(data->columns);
    b = ! columns.empty();
  }
  return b;
}


int FitsTable::GetColType(int c) const
{
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  const vector<FitsTableColumn> columns(data->columns);
  /*
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  wxASSERT(data && 0 <= x && x < data->ncols);
  return data->typecode + c;
  */
  wxASSERT(0 <= c && c < (int) columns.size());
  return columns[c].GetColType();
}


// const float *FitsTable::GetCol_float(int c) const
// {
//   /*
//   FitsTableData *data = static_cast<FitsTableData *>(m_refData);
//   wxASSERT(data && 0 <= x && x < data->ncols && data->typecode[c] == TFLOAT);
//   return (float *) data->table[c];
//   */
//   /*
//   wxASSERT(0 <= c && c < (int) columns.size());
//   return columns[c].GetCol_float();
//   */
//   //  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
//   FitsTableData *data = static_cast<FitsTableData *>(m_refData);
//   const vector<FitsTableColumn> columns(data->columns);
//   return data->ftable;
// }

// const int *FitsTable::GetCol_int(int c) const
// {
//   /*
//   FitsTableData *data = static_cast<FitsTableData *>(m_refData);
//   wxASSERT(data && 0 <= x && x < data->ncols && data->typecode[c] == TLONG);
//   return (int *) data->table[c];
//   */
//   wxASSERT(0 <= c && c < (int) columns.size());
//   //  return columns[c].GetCol_int();
//   FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
//   return data->itable;
// }

// const char **FitsTable::GetCol_char(int c) const
// {
//   /*
//   FitsTableData *data = static_cast<FitsTableData *>(m_refData);
//   wxASSERT(data && 0 <= x && x < data->ncols && data->typecode[c] == TSTRING);
//   return (char **) data->table[c];
//   */
//   FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
//   return (const char **) data->ctable;
//   /*
//   wxASSERT(0 <= c && c < (int) columns.size());
//   return columns[c].GetCol_char();
//   */
// }

FitsTableColumn FitsTable::GetColumn(int k) const
{
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  const vector<FitsTableColumn> columns(data->columns);
  wxASSERT(0 <= k && k < (int) columns.size());
  return columns[k];
}


FitsTableColumn FitsTable::GetColumn(const wxString& label) const
{
  return GetColumn(GetColIndex(label));
}

int FitsTable::GetColIndex(const wxString& label) const
{
  for(size_t i = 0; i < GetCount(); i++) {
    wxString key,val,com;

    FitsHeader::ParseRecord(Item(i),key,val,com);

    if( label.StartsWith(val) ) {

      wxRegEx re("TTYPE([0-9]+)");
      wxASSERT(re.IsValid());

      if( re.Matches(key) ) {
	wxString a(re.GetMatch(key,1));
	long i;
	if( a.ToLong(&i) )
	  return i - 1;
      }
    }
  }
  return -1;
}


wxArrayString FitsTable::GetColLabels() const
{
  wxArrayString cols;

  for(size_t i = 0; i < GetCount(); i++) {
    wxString key,val,com;

    FitsHeader::ParseRecord(Item(i),key,val,com);

    if( key.StartsWith("TTYPE") )
      cols.Add(val);
  }
  return cols;
}


/*
const float *FitsTable::GetColumn(int x) const
{
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  wxASSERT(data && 0 <= x && x < data->ncols);
  return data->table+x*data->nrows;
}
*/


void FitsTable::GetStarChart(wxOutputStream& output)
{
  wxTextOutputStream cout(output);

  cout << "<svg xmlns=\"http://www.w3.org/2000/svg\">" << endl;

  for(int i = 0; i < Naxes(1); i++) {

    float x = Pixel(0,i);
    float y = Pixel(1,i);

    float f = Pixel(11,i);

    if( f > 0 ) {
      float r = wxMax(f/1e4,10.0);
      cout << "<circle cx=\"" << x << "\" cy=\"" << y
	   << "\" r=\"" << r << "\"/>" << endl;
    }
  }

  cout << "</svg>" << endl;
}


// inline float FitsTable::Cell(int x, int y) const
// {
//   FitsTableData *data = static_cast<FitsTableData *>(m_refData);
//   wxASSERT(data);
//   const std::vector<FitsTableColumn> columns(data->columns);

//   wxASSERT(0 <= y && y < (int) columns.size());
//   return columns[y].Cell(x);
// 	   //  FitsTableColumnData data(columns[y].GetData());
//   // FitsTableColumn col(columns[y]);
//   // FitsTableColumnData *data = static_cast<FitsTableColumnData *>(col.m_refData);
//   // 	   //  wxASSERT(data && 0 <= x && x < data->nrows);
//   // 	   //  wxASSERT(0 <= x && x < data.Nrows());

//   // if( data->ftable )
//   //   return data->ftable[x];
//   // else
//   //   return 0.0;

//   /*
//   FitsTableData *data = static_cast<FitsTableData *>(m_refData);
//   wxASSERT(data && 0 <= x && x < data->ncols && 0 <= y && y < data->nrows);
//   float *p = (float *) data->table[x*data->nrows + y];
//   return (float) *p;
//   */
//   //return data->Cell_float(x,y);
// }


float FitsTable::Cell(int row, int col) const
{
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  wxASSERT(data);
  const std::vector<FitsTableColumn> columns(data->columns);

  wxASSERT(0 <= col && col < (int) columns.size());
  return columns[col].Cell(row);
}

wxString FitsTable::Cell_str(int x, int y) const
{
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  wxASSERT(data);
  const std::vector<FitsTableColumn> columns(data->columns);

  wxASSERT(0 <= y && y < (int) columns.size());
  return columns[y].Cell_str(x);

  /*
  wxASSERT(0 <= y && y < (int) columns.size());
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(columns[y].m_refData);
  wxASSERT(data && 0 <= x && x < data->nrows);

  wxString a;

  if( data->ftable )
    a.Printf("%g",data->ftable[x]);
  else if( data->itable )
    a.Printf("%d",data->itable[x]);
  else if( data->ctable )
    a = data->ctable[x];
  */

  /*
  switch(data->typecode[y]) {
  case TBYTE:
  case TINT:
  case TLONG:
    a.Printf("%d",(int **) data->table[x*data->nrows + y]); break;

  case TFLOAT:
     break;

  case TDOUBLE:
    a.Printf("%lf",(double **) data->table[x*data->nrows + y]); break;

  case TSTRING:
    a = wxString((char *) data->table[x*data->nrows + y]); break;

  }
  */

  //  wxString a;
  //  a.Printf("%g",Cell(x,y));
  //  return a;
}

float FitsTableColumn::Cell(int x) const
{
  //  wxASSERT(0 <= y && y < (int) columns.size());
	   //  FitsTableColumnData data(columns[y].GetData());
	   //  FitsTableColumn col(columns[y]);
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
	   //  wxASSERT(data && 0 <= x && x < data->nrows);
  wxASSERT(data && 0 <= x && x < data->nrows);

	   //if( data->ftable )
  //  return data->ftable[x];
	   //  else
	   //    return 0.0;


  if( data->ftable )
    return data->ftable[x];
  else if( data->dtable )
    return data->dtable[x];
  else if( data->stable )
    return data->stable[x];
  else if( data->ltable )
    return data->ltable[x];
  else if( data->itable )
    return data->itable[x];
  else if( data->btable )
    return data->btable[x];
  else if( data->otable )
    return data->otable[x] ? 1 : 0;
  else if( data->ctable ) {
    wxFAIL_MSG("FitsTableColumn::Cell: char to float unimplemented.");
  }

  return 0.0;

  /*
  FitsTableData *data = static_cast<FitsTableData *>(m_refData);
  wxASSERT(data && 0 <= x && x < data->ncols && 0 <= y && y < data->nrows);
  float *p = (float *) data->table[x*data->nrows + y];
  return (float) *p;
  */
  //return data->Cell_float(x,y);
}

wxString FitsTableColumn::Cell_str(int x) const
{
  /*
  wxASSERT(0 <= y && y < (int) columns.size());
  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data && 0 <= x && x < data->nrows);
  */

  FitsTableColumnData *data = static_cast<FitsTableColumnData *>(m_refData);
  wxASSERT(data && 0 <= x && x < data->nrows);

  wxString a;

  if( data->ftable )
    a.Printf("%g",data->ftable[x]);
  else if( data->dtable )
    a.Printf("%lg",data->dtable[x]);
  else if( data->stable )
    a.Printf("%d",data->stable[x]);
  else if( data->ltable )
    a.Printf("%ld",data->ltable[x]);
  else if( data->itable )
    a.Printf("%d",data->itable[x]);
  else if( data->btable )
    a.Printf("%d",data->btable[x]);
  else if( data->otable )
    a = data->otable[x] ? "T" : "F";
  else if( data->ctable )
    a = data->ctable[x];

//wxLogDebug("cell_str: "+a);
  return a;
}

// class FitsTableData: public wxObjectRefData
// {
// public:
//   FitsTableData();
//   //  FitsTableData(long, long, float *);
//   FitsTableData(long);
//   FitsTableData(const FitsTableData&);
//   FitsTableData& operator = (const FitsTableData&);
//   virtual ~FitsTableData();

//   void InsertColumn(long, long, float *);
//   void InsertColumn(long, long, int *);
//   void InsertColumn(long, long, char **);

//   inline float Cell_float(int,int) const;
//   //  inline int Cell_int(int,int) const;
//   //  inline char *Cell_char(int,int) const;
//   //  inline bool Cell_bool(int,int) const;

//   //  long ncols;
//   //  float *table;
//   //  int *typecode;
//   //  void *table;
//   //  FitsTableColumn *columns;
//   std::vector<FitsTableColumn> columns;
// };

// inline float FitsTableData::Cell_float(int x, int y) const
// {
//   wxASSERT(table && 0 <= x && x < ncols && 0 <= y && y < nrows && table[y] && typecode[y] == TFLOAT);
//   float *p = (float *) table[x*nrows + y];
//   return *p;
// 	   //   return (float) table[x*nrows + y];
// }
