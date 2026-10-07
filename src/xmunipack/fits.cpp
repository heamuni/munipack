/*

  xmunipack - fits implementation

  Copyright © 2011 - 2026 F.Hroch (hroch@physics.muni.cz)

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
#include <string>
#include <cstdio>
#include <fitsio.h>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <wx/wx.h>
#include <wx/tokenzr.h>
#include <wx/regex.h>
#include <wx/uri.h>
#include <wx/filesys.h>
#include <wx/txtstrm.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


using namespace std;

// we know the types of images:
// * 1d (spectrum)
// * 2d grayscale
// * 3d colour (3d, unique sizes, header with CSPACE keyword)
// * 3d (3d spectrum) ?
// * multilayer (more layers, arbitrary sizes)
//




// ---  FitsHeader (HDU)

int FitsHeader::GetNaxis() const
{
  wxString l = GetKey("NAXIS");
  long n;
  if( l.ToLong(&n) )
    return n;
  else
    return 0;
}



int FitsHeader::Bitpix() const
{
  wxString l = GetKey("BITPIX");
  long b;
  if( l.ToLong(&b) )
    return b;
  else
    return 0;
}


wxString FitsHeader::Bitpix_str() const
{
  int bitpix = Bitpix();

  switch (bitpix) {
  case 8: return   "0 .. 255";
  case 16: return  "0 .. 65535";
  case 32: return  "0 .. 4294967296";
  case -32:return L"0 .. ± 10^38, 6 digits";
  default:
    wxString a;
    a.Printf("Bitpix = %d",bitpix);
    return a;
  }
}

bool FitsHeader::ParseRecord(const wxString& record, wxString& key,
			     wxString& value, wxString& comment)
{
  int status = 0;
  const wxScopedCharBuffer name(record.ToAscii());
  char *card = (char *) name.data();
  char keyname[FLEN_CARD],cval[FLEN_CARD],com[FLEN_CARD];
  int keylen;

  fits_get_keyname((char*)card,keyname,&keylen,&status);
  fits_parse_value((char*)card,cval,com,&status);

  if( status == 0 ) {
    key = keyname;
    comment = com;
    value = cval;

    char *i1 = index(cval,'\'');
    char *i2 = rindex(cval,'\'');
    if( i1 != NULL && i2 != NULL && i1 != i2 ) {
      int l = i2 - i1 - 1;
      for(int i = 0; i < l; i++)
	cval[i] = *(i1 + i + 1);
      cval[l] = '\0';
    }

    value = cval;
    value.Trim();


    return true;
  }
  else
    return false;
}

bool FitsHeader::FindKey(const wxString& keyword, wxString& value, wxString& comment) const
{
  for(size_t i = 0; i < GetCount(); i++) {
    wxString key,val,com;
    if( ParseRecord(Item(i),key,val,com) ) {
      wxStringTokenizer tb = wxStringTokenizer(keyword,",");
      while( tb.HasMoreTokens() ) {
	wxString k = tb.GetNextToken();
	if( key == k ) {
	  value = val; comment = com;
	  return true;
	}
      }
    }
  }
  return false;
}

wxString FitsHeader::GetKey(const wxString& key) const
{
  wxString l,c;
  if( FindKey(key,l,c) )
    return l;
  return wxEmptyString;
}

wxString FitsHeader::GetUnit(const wxString& key) const
{
  wxString l,c;
  if( FindKey(key,l,c) ) {
    wxRegEx re("^\\s*\\[.*\\].*");
    wxASSERT(re.IsValid());
    if( re.Matches(c) ) {
      size_t start,len;
      re.GetMatch(&start,&len);
      wxString s = c.SubString(start,len);
      int i1 = s.Index('[') + 1;
      int i2 = s.Index(']') - 1;
      return c.SubString(i1,i2);
    }
  }
  return wxEmptyString;
}

wxString FitsHeader::Exposure_str(const wxString& key) const
{
  wxString exp = GetKey(key);
  double e;
  if( ! exp.IsEmpty() && exp.ToDouble(&e) ) {
    wxString line;
    line.Printf("%g",e);
    return line;
  }
  else
    return exp;
}

hdu_type FitsHeader::GetType(int htype) const
{
  hdu_type hdutype = HDU_UNKNOWN;

  if( htype == IMAGE_HDU ) {
    wxString axes = GetKey("NAXIS");
    long naxes;
    if( axes != "" && axes.ToLong(&naxes) )
      hdutype = naxes > 0 ? HDU_IMAGE : HDU_HEAD;
    else
      wxLogError("FitsHeader::GetType(): NAXIS undefined.");
  }
  else if ( htype == ASCII_TBL || htype == BINARY_TBL )
    hdutype = HDU_TABLE;

  return hdutype;
}

hdu_flavour FitsHeader::GetFlavour(int htype, hdu_type hdutype) const
{
  hdu_flavour subtype = HDU_DUMMY;

  if( hdutype == HDU_HEAD )
    subtype = HDU_DUMMY;

  else if( hdutype == HDU_IMAGE ) {

    wxString axis = GetKey("NAXIS");
    long naxis;
    if( axis.ToLong(&naxis) ) {
      if( naxis == 1 ) {
	wxString axis1 = GetKey("NAXIS1");
	long naxis1;
	if( axis1.ToLong(&naxis1) && naxis1 > 0 )
	  subtype = HDU_IMAGE_LINE;
	else
	  subtype = HDU_IMAGE_UNKNOWN;
      }
      else if( naxis == 2 ) {
	wxString axis1 = GetKey("NAXIS1");
	wxString axis2 = GetKey("NAXIS2");
	long naxis1, naxis2;
	if( axis1.ToLong(&naxis1) && naxis1 > 0 &&
	    axis2.ToLong(&naxis2) && naxis2 > 0)
	  subtype = HDU_IMAGE_FRAME;
	else
	  subtype = HDU_IMAGE_UNKNOWN;
      }
      else if ( naxis == 3 ) {
	wxString cspace = GetKey("CSPACE");
	wxString axis3 = GetKey("NAXIS3");
	long naxis3;
	if( axis3.ToLong(&naxis3) && naxis3 == 3 && cspace.Find("XYZ")!=wxNOT_FOUND)
	  subtype = HDU_IMAGE_COLOUR;
	else
	  subtype = HDU_IMAGE_CUBE;
      }
      else
	subtype = HDU_IMAGE_UNKNOWN;
    }
    else
      subtype = HDU_IMAGE_UNKNOWN;
  }
  else if( hdutype == HDU_TABLE ) {
    switch(htype) {
    case ASCII_TBL: subtype = HDU_TABLE_BIN; break;
    case BINARY_TBL:subtype = HDU_TABLE_ASCII; break;
    default: subtype = HDU_TABLE_UNKNOWN;
    }
  }

  return subtype;
}

wxString FitsHeader::GetType_str(hdu_type hdutype)
{
  switch (hdutype) {
  case HDU_HEAD:   return "Head";
  case HDU_IMAGE:  return "Image";
  case HDU_TABLE:  return "Table";
  case HDU_UNKNOWN:return "HDU unknown";
  }
  return "HDU unknown";
}

wxString FitsHeader::GetFlavour_str(hdu_flavour flavour)
{
  switch (flavour) {
  case HDU_IMAGE_LINE: return "Line";
  case HDU_IMAGE_FRAME: return "Grey picture";
  case HDU_IMAGE_CUBE: return "Cube";
  case HDU_IMAGE_COLOUR: return "Colour picture";
  case HDU_IMAGE_UNKNOWN: return "Unknown picture";
  case HDU_TABLE_BIN: return "Binary table";
  case HDU_TABLE_ASCII: return "ASCII table";
  case HDU_TABLE_UNKNOWN: return "Unknown table";
  case HDU_DUMMY: return "Dummy";
  }
  return "Flavour unknown";
}


// ------ FitsHDU

FitsHdu::FitsHdu():
  htype(IMAGE_HDU),type(HDU_UNKNOWN), flavour(HDU_DUMMY), modified(false) {}

FitsHdu::FitsHdu(const FitsHeader& h, int t):
  header(h),htype(t),type(header.GetType(htype)),
  flavour(header.GetFlavour(htype,type)), modified(false) {}

FitsHdu::FitsHdu(const FitsHeader& h, int ht, hdu_type t, hdu_flavour f):
  header(h),htype(ht),type(t),flavour(f),modified(false) {}

size_t FitsHdu::GetCount() const
{
  return header.GetCount();
}

wxString FitsHdu::Item(size_t i) const
{
  return header.Item(i);
}

wxString FitsHdu::GetKey(const wxString& a) const
{
  return header.GetKey(a);
}

long FitsHdu::GetKeyLong(const wxString& key) const
{
  wxString l = GetKey(key);
  long i;
  return l.ToLong(&i) ? i : 0;
}

double FitsHdu::GetKeyDouble(const wxString& key) const
{
  wxString l = GetKey(key);
  double x;
  return l.ToDouble(&x) ? x : 0.0;
}

wxString FitsHdu::GetUnit(const wxString& a) const
{
  return header.GetUnit(a);
}

int FitsHdu::Bitpix() const
{
  return header.Bitpix();
}

wxString FitsHdu::Bitpix_str() const
{
  return header.Bitpix_str();
}

wxString FitsHdu::Exposure_str(const wxString& a) const
{
  return header.Exposure_str(a);
}

int FitsHdu::GetHduType() const
{
  return htype;
}

hdu_type FitsHdu::Type() const
{
  return type;
}

wxString FitsHdu::Type_str() const
{
  return header.GetType_str(type);
}

hdu_flavour FitsHdu::Flavour() const
{
  return flavour;
}

wxString FitsHdu::Flavour_str() const
{
  return header.GetFlavour_str(flavour);
}

int FitsHdu::Naxis() const
{
  return 0;
}

long FitsHdu::Naxes(int n) const
{
  return 0;
}

long FitsHdu::Width() const
{
  return Naxes(0);
}

long FitsHdu::Height() const
{
  return Naxes(1);
}

long FitsHdu::GetWidth() const
{
  return Naxes(0);
}

long FitsHdu::GetHeight() const
{
  return Naxes(1);
}

long FitsHdu::GetDepth() const
{
  switch (Naxis()) {
  case 2:  return 1;
  case 3:  return Naxes(2);
  default: return 0;
  }
}


bool FitsHdu::IsOk() const
{
  return ! header.IsEmpty();
}

bool FitsHdu::IsColour() const
{
  return type == HDU_IMAGE && flavour == HDU_IMAGE_COLOUR;
}

bool FitsHdu::IsDisplayImplemented() const
{
  return type == HDU_IMAGE &&
    (flavour == HDU_IMAGE_FRAME || flavour == HDU_IMAGE_COLOUR);
}

wxString FitsHdu::GetExtname() const
{
  return GetKey("EXTNAME");
}

wxString FitsHdu::GetControlLabel() const
{
  wxString label = GetKey("EXTNAME");
  if( label.IsEmpty() )
    label = Type_str();
  return label;
}

FitsHeader FitsHdu::GetHeader() const
{
  return header;
}


bool FitsHdu::GetWCS(double& xpix, double& ypix, double& alpha, double& delta, double& scale, double& angle, double& reflex) const
{
  // projection is silently ignored
  wxString type(GetKey("CTYPE1"));

  if( type.IsEmpty() ) return false;

  xpix = GetKeyDouble("CRPIX1");
  ypix = GetKeyDouble("CRPIX2");

  alpha = GetKeyDouble("CRVAL1");
  delta = GetKeyDouble("CRVAL2");

  double cd11,cd12,cd21,cd22;
  cd11 = GetKeyDouble("CD1_1");
  cd12 = GetKeyDouble("CD1_2");
  cd21 = GetKeyDouble("CD2_1");
  cd22 = GetKeyDouble("CD2_2");

  reflex = -cd11*cd22 >= 0 ? 1.0 : -1.0;
  scale = 1.0/sqrt(cd11*cd11 + cd12*cd12);
  angle = (45.0/atan(1.0)) * atan2(cd21,cd22);

  return true;
}


int FitsHdu::hdr2str(int nocomments, char **exclist, int nexc,
		     char **blockhead, int *nkeys, int *status) const
{
  if( *status != 0 ) return *status;

  *nkeys = GetCount();
  *blockhead = (char *) malloc(*nkeys*80+1);

  if( *blockhead == NULL ) return 1;

  char *h = *blockhead;
  for(size_t i = 0; i < header.GetCount(); i++) {
    wxString a = header[i];
    const char *b = (const char *) a.fn_str();
    size_t blen = strlen(b);
    strncpy(h,b,80);
    memset(h + blen,' ',80 - blen + 1);
    h = h + 80;
  }
  /*
   *h = '\0';
   fprintf(stderr, "head: `%s' \n", *blockhead);
  */

  return 0;
}


// ----------    FitsFile

FitsFile::FitsFile(): status(false),type(FITS_UNKNOWN) {}

FitsFile::FitsFile(const wxString& name):
  filename(name),status(false),type(FITS_UNKNOWN)
{
  fitsfile *f;
  int stat = 0;
  int dummy, htype, bitpix, naxis;

  int nhdu = 0;

  // open file
  stat = 0;
  fits_open_file(&f, name.fn_str(), READONLY, &stat);
  if( stat ) goto crash;

  fits_get_num_hdus(f,&nhdu,&stat);
  if( stat ) goto crash;

  for(int k = 0; k < nhdu; k++) {

    fits_movabs_hdu(f,k+1,&htype,&stat);
    if( stat ) goto crash;

    // load header
    int nhead;
    char h[FLEN_CARD];
    FitsHeader head;
    fits_get_hdrspace(f,&nhead,&dummy,&stat);
    for(int n = 0; stat == 0 && n < nhead; n++) {
      if( fits_read_record(f,n+1,h,&stat) == 0 )
	head.Add(wxString(h,wxConvUTF8));
    }
    if( stat ) goto crash;

    // load data
    if( htype == IMAGE_HDU ) {

      fits_get_img_type(f,&bitpix,&stat);
      fits_get_img_dim(f,&naxis,&stat);
      if( stat ) goto crash;

      if( naxis > 0 ) {

	long *naxes = new long[naxis];
	fits_get_img_size(f,naxis,naxes,&stat);
	if( stat ) { delete[] naxes; goto crash; }

	long ndata = 1; for(int i = 0; i < naxis; i++ ) ndata = ndata*naxes[i];
	long firstelem = 1;
	float nullval = 0.0;

	float *image = new float[ndata];
	wxASSERT(image);
	fits_read_img(f,TFLOAT,firstelem,ndata,&nullval,image,&dummy,&stat);
	if( stat ) { delete[] naxes; delete[] image; goto crash; }

	hdu.push_back(FitsArray(head,htype,naxis,naxes,image));
      }
      else {
	hdu.push_back(FitsHdu(head,htype));
      }

    }
    else if( htype == ASCII_TBL || htype == BINARY_TBL ) {

      long nrows, ncols;
      int nc;
      fits_get_num_rows(f,&nrows,&stat);
      fits_get_num_cols(f,&nc,&stat);
      if( stat ) goto crash;

      ncols = nc;
      //      float *table = new float[nrows*ncols];
      //      void **table = new void*[ncols];

      //      FitsTable table(head,htype,nrows,ncols);
      vector<FitsTableColumn> cols;

      long frow = 1, felem = 1;
      for(int k = 0; k < ncols; k++ ) {

	int colnum = k + 1;
	int typecode;
	long repeat, width;

	fits_get_coltype(f,colnum,&typecode,&repeat,&width,&stat);
	if( stat ) goto crash;

	if( typecode == TSTRING ) {
	  int width;
	  fits_get_col_display_width(f,colnum,&width,&stat);
	  char **a = new char*[nrows];
	  for(int i = 0; i < nrows; i++)
	    a[i] = new char[width];
	  char nullval[width];
	  strcpy(nullval," ");
	  fits_read_col(f,TSTRING,colnum, frow, felem, nrows, &nullval,
			a,&dummy,&stat);
	  cols.push_back(FitsTableColumn(typecode,nrows,a));
	}
	else if( typecode == TLOGICAL ) {
	  bool *b = new bool[nrows];
	  bool nullval = false;
	  fits_read_col(f,TLOGICAL,colnum, frow, felem, nrows, &nullval,
			b,&dummy,&stat);
	  cols.push_back(FitsTableColumn(typecode,nrows,b));
	}
	else if( typecode == TBYTE || typecode == TBIT ) {
	  char *b = new char[nrows];
	  char nullval = 0;
	  fits_read_col(f,TBYTE,colnum, frow, felem, nrows, &nullval,
			b,&dummy,&stat);
	  cols.push_back(FitsTableColumn(typecode,nrows,b));
	}
	else if( typecode == TSHORT ) {
	  short *d = new short[nrows];
	  short nullval = 0;
	  fits_read_col(f,TSHORT,colnum, frow, felem, nrows, &nullval,
			d,&dummy,&stat);
	  cols.push_back(FitsTableColumn(typecode,nrows,d));
	}
	else if( typecode == TLONG || typecode == TINT32BIT ) {
	  long *d = new long[nrows];
	  long nullval = 0;
	  fits_read_col(f,TLONG,colnum, frow, felem, nrows, &nullval,
			d,&dummy,&stat);
	  cols.push_back(FitsTableColumn(typecode,nrows,d));
	}
	else if( typecode == TFLOAT ) {
	  float *d = new float[nrows];
	  float nullval = 0.0;
	  fits_read_col(f,TFLOAT,colnum, frow, felem, nrows, &nullval,
			d,&dummy,&stat);
	  cols.push_back(FitsTableColumn(typecode,nrows,d));
	}
	else if( typecode == TDOUBLE ) {
	  double *d = new double[nrows];
	  double nullval = 0.0;
	  fits_read_col(f,TDOUBLE,colnum, frow, felem, nrows, &nullval,
			d,&dummy,&stat);
	  cols.push_back(FitsTableColumn(typecode,nrows,d));
	}
	else
	  wxLogDebug("FitsFile::FitsFile: The type code `%d' not implemented yet.",typecode);

	/*
	fits_read_col(f, TFLOAT,colnum, frow, felem, nrows, &nullval,
		      table+k*nrows,&dummy,&stat);
	*/
#ifdef __WXDEBUG__
	if( stat ) { wxLogDebug("FITS unknown fail: %d (check!)",stat); }
#endif
	//	stat = 0;
      }
      //      FitsTable t(head,htype, nrows, ncols, table);
      //      hdu.push_back(t);
      //      FitsTable table(head,htype,cols);
      //      hdu.push_back(table);
      wxASSERT(stat == 0);
      hdu.push_back(FitsTable(head,htype,cols));
    }
    else {
      hdu.push_back(FitsHdu(head,htype));
    }

    // don't delete image and naxes (both will deleted by FitsHdu) !!!!!
  }

  // close file
  fits_close_file(f, &stat);
  if( stat ) goto crash;

  // classify image
  Recognise();
  status = true;
  return;

 crash:

  // save error description
  char emsg[FLEN_ERRMSG];
  while( fits_read_errmsg(emsg) )
    errmsg.Add(wxString(emsg,wxConvUTF8));

  char msg[FLEN_STATUS];
  fits_get_errstatus(stat,msg);
  smsg = wxString(msg,wxConvUTF8);
  wxLogDebug("FitsFile: " + smsg);

  filename = "";
  status = false;
  fits_close_file(f, &stat);
}

FitsFile::FitsFile(const wxString& name, const vector<FitsHdu>& h):
  filename(name),status(true),type(FITS_UNKNOWN),hdu(h)
{
  Recognise();
}


FitsFile::~FitsFile()
{
  hdu.clear();
}

bool FitsFile::Status() const
{
  return status;
}

size_t FitsFile::HduCount() const
{
  return hdu.size();
}

int FitsFile::size() const
{
  return hdu.size();
}

FitsHdu FitsFile::Hdu(size_t n) const
{
  wxASSERT(0 <= n && n < hdu.size());
  return hdu[n];
}

FitsHdu FitsFile::FindHdu(const wxString& extname) const
{
  for(size_t n = 0; n < hdu.size(); n++) {
    if( hdu[n].GetExtname() == extname )
      return hdu[n];
  }
  return FitsHdu();
}


void FitsFile::Recognise()
{
  type = FITS_UNKNOWN;
  int nhdu = hdu.size();

  if( nhdu == 1 && hdu[0].Type() == HDU_IMAGE ) {

    if( hdu[0].Flavour() == HDU_IMAGE_FRAME )
      // simple, grayscale image
      type = FITS_GREY;

    else if( hdu[0].Flavour() == HDU_IMAGE_COLOUR )
      // colour image
      type = FITS_COLOUR;
  }
  else if( nhdu > 1 )
    // more hdus
    type = FITS_MULTI;

}

int FitsFile::Type() const
{
  return type;
}

wxString FitsFile::Type_str() const
{
  switch (type) {
  case FITS_GREY:  return "Grey image";
  case FITS_COLOUR:return "Colour image";
  case FITS_MULTI: return "Multi-extension";
  default: return "";
  }
}


bool FitsFile::HasImage() const
{
  // locate image hdu
  for(size_t k = 0; k < HduCount(); k++)
    if( hdu[k].Type() == HDU_IMAGE )
      return true;
  return false;
}

bool FitsFile::HasFind() const
{
  // locate FIND extension
  for(size_t k = 0; k < HduCount(); k++)
    if( hdu[k].GetExtname() == FINDEXTNAME )
      return true;
  return false;
}

bool FitsFile::HasPhotometry() const
{
  // locate APHOT extension
  for(size_t k = 0; k < HduCount(); k++)
    if( hdu[k].GetExtname() == APEREXTNAME )
      return true;
  return false;
}

bool FitsFile::HasPhcal() const
{
  // locate photometry extension
  for(size_t k = 0; k < HduCount(); k++)
    if( hdu[k].GetExtname() == PHOTOEXTNAME )
      return true;
  return false;
}

wxString FitsFile::GetURL() const
{
  return wxFileSystem::FileNameToURL(filename);
}

wxString FitsFile::GetName() const
{
  wxFileName name(filename);
  return name.GetName();
}

wxString FitsFile::GetFullName() const
{
  wxFileName name(filename);
  return name.GetFullName();
}

wxString FitsFile::GetPath() const
{
  wxFileName name(filename);
  if( name.GetPath().IsEmpty() )
    name.AssignCwd();
  return name.GetPath();
}

wxString FitsFile::GetFullPath() const
{
  wxFileName name(filename);
  if( name.GetPath().IsEmpty() ) {
    name.AssignCwd();
    name.SetFullName(filename);
  }
  return name.GetFullPath();
}

bool FitsFile::IsOk() const
{
  return ! filename.IsEmpty();
}

bool FitsFile::IsModified() const
{
  for(size_t k = 0; k < HduCount(); k++) {
    if( Hdu(k).IsModified() )
      return true;
  }
  return false;
}


wxArrayString FitsFile::GetErrorMessage() const
{
  return errmsg;
}

wxString FitsFile::GetErrorDescription() const
{
  return smsg;
}

void FitsFile::Clear()
{
  filename.Clear();
  status = false;
  type = FITS_UNKNOWN;
  hdu.clear();
  errmsg.Clear();
  smsg.Clear();
}


FitsFile::FitsFile(const FitsHdu& h):
  status(true),type(h.Type())
{
  hdu.push_back(h);
}


bool FitsFile::Save(const wxString& name)
{
  int status = 0;
  fitsfile *f;

  fits_create_file(&f,name.fn_str(),&status);
  if( status != 0 ) return false;

  for(size_t k = 0; k < HduCount(); k++) {
    int type = Hdu(k).Type();

    if( type == HDU_IMAGE ) {

      const FitsArray image(Hdu(k));

      long naxis = image.Naxis();
      long *naxes = new long[naxis];
      long fpixel = 1;
      long nelements = 1;
      for(int i = 0; i < naxis; i++) {
	naxes[i] = image.Naxes(i);
	nelements = nelements * naxes[i];
      }

      fits_create_img(f, FLOAT_IMG, naxis, naxes, &status);
      merge_head(f,Hdu(k),&status);
      fits_write_img(f, TFLOAT, fpixel, nelements, (float *) image.PixelData(), &status);


    }
    else if( type == HDU_TABLE ) {

      const FitsTable table(Hdu(k));

      int tfields = table.Ncols();
      char **ttype = new char*[tfields];
      char **tform = new char*[tfields];
      char **tunit = new char*[tfields];

      for(int i = 0; i < table.Ncols(); i++) {
	wxString key;
	key.Printf("TTYPE%d",i+1);
	ttype[i] = wxStrdup(table.GetKey(key));
	key.Printf("TFORM%d",i+1);
	tform[i] = wxStrdup(table.GetKey(key));
	key.Printf("TUNIT%d",i+1);
	tunit[i] = wxStrdup(table.GetKey(key));
      }

      wxString extname(table.GetExtname());

      fits_create_tbl(f,BINARY_TBL,table.Nrows(),table.Ncols(), ttype, tform,
		      tunit, extname.fn_str(), &status);

      delete[] ttype;
      delete[] tform;
      delete[] tunit;

      merge_head(f,table,&status);

      /*
      for(size_t i = 0; i < table.GetCount(); i++) {
	fits_write_record(f,table.Item(i).fn_str(),&status);
      }
      */

      //      wxLogDebug("%d",(int) status);

      long firstrow  = 1;
      long firstelem = 1;
      long nelements = table.Nrows();

      for(int colnum = 0; colnum < table.Ncols(); colnum++) {

	int typecode = table.GetColType(colnum);
	//	wxString form(tform[colnum]);

	//	wxLogDebug(form+" %d",(int)status);
	const FitsTableColumn col(table.GetColumn(colnum));

	if( typecode == TSTRING ) {
	  /*
	  char **a = new char*[nelements];
	  for(int i = 0; i < nelements; i++)
	    a[i] = strdup("");
	  */
	  char **a = (char **) col.GetCol_char();
	  fits_write_col(f, TSTRING, colnum+1, firstrow, firstelem,
			 nelements, a, &status);
//  delete[] a;
	}
	//	else if( form.Find("D") != wxNOT_FOUND ) {
	else if( typecode == TDOUBLE ) {
	  /*
	  double *u = new double[nelements];
	  const float *col = table.GetCol(colnum);
	  for(int i = 0; i < nelements; i++)
	    u[i] = col[i];
	  */
	  //	  wxLogDebug("%f %d %d",col[0],(int)nelements, (int)status);
	  double *d = (double *) col.GetCol_double();
	  fits_write_col(f, TDOUBLE, colnum+1, firstrow, firstelem,
			 nelements, d, &status);
	  //	  delete[] u;
	}
	//	else if( form.Find("J") != wxNOT_FOUND ) {
	else if( typecode == TLONG ) {
	  /*
	  long *d = new long[nelements];
	  const float *col = table.GetCol(colnum);
	  for(int i = 0; i < nelements; i++)
	    d[i] = (int) col[i];
	  */
	  long *d = (long *) col.GetCol_long();
	  fits_write_col(f, TLONG, colnum+1, firstrow, firstelem,
			 nelements, d, &status);
	  //	  delete[] d;
	}
      }
    }
  }

  fits_close_file(f,&status);


  bool s = status == 0;

  if( s )
    filename = name;

  return s;
}

int FitsFile::merge_head(fitsfile *f, const FitsHdu& hdu, int *status) const
{
  int nhead,dummy;
  char h[FLEN_CARD];
  wxArrayString head;

  fits_get_hdrspace(f,&nhead,&dummy,status);

  for(int n = 0; *status == 0 && n < nhead; n++) {
    if( fits_read_record(f,n+1,h,status) == 0 )
      head.Add(wxString(h,wxConvUTF8));
  }

  for(size_t i = 0; i < hdu.GetCount(); i++) {

    wxString record(hdu.Item(i));
    bool presented = false;
    wxString xkey,ykey,value,com;

    FitsHeader::ParseRecord(record,xkey,value,com);
    for(size_t j = 0; j < head.GetCount(); j++) {
      FitsHeader::ParseRecord(head[j],ykey,value,com);
      if( xkey == ykey ) {
	presented = true;
	break;
      }
    }

    if( ! presented )
      fits_write_record(f,record.fn_str(),status);
  }

  return *status;
}



// auxiliary  functions

bool FitsCopyFile(const wxString& in, const wxString& out)
{
  fitsfile *fin,*fout;
  int stat = 0;
  bool result = false;

  fits_open_file(&fin, in.fn_str(), READONLY, &stat);
  fits_create_file(&fout, out.fn_str(), &stat);
  result = fits_copy_file(fin,fout,1,1,1,&stat) == 0;
  fits_close_file(fout,&stat);
  fits_close_file(fin,&stat);

  fits_report_error(stderr,stat);

  return result && stat == 0;
}

bool FitsCopyHdu(const wxString& in, const wxString& out, const wxString& hdu)
{
  fitsfile *fin,*fout;
  int stat = 0;
  int hdutype;
  bool result = false;
  const char *extname = hdu.ToAscii();

  fits_open_file(&fin, in.fn_str(), READONLY, &stat);
  fits_open_file(&fout, out.fn_str(), READWRITE, &stat);
  fits_movnam_hdu(fin,BINARY_TBL,(char*)extname,0,&stat);
  if( fits_movnam_hdu(fout,BINARY_TBL,(char*)extname,0,&stat) == 0 )
    fits_delete_hdu(fout,&hdutype,&stat);
  if( stat == BAD_HDU_NUM )
    stat = 0;
  result = fits_copy_hdu(fin,fout,0,&stat) == 0;
  fits_close_file(fout,&stat);
  fits_close_file(fin,&stat);

  fits_report_error(stderr,stat);

  return result && stat == 0;
}


wxArrayString FitsColumns(const wxString& name)
{
  fitsfile *f;
  int stat = 0;
  int keysexist, morekeys;
  char keyname[FLEN_KEYWORD],value[FLEN_VALUE],com[FLEN_COMMENT];
  wxArrayString cols;

  fits_open_table(&f, name.fn_str(), READONLY, &stat);
  fits_get_hdrspace(f,&keysexist,&morekeys,&stat);
  for(int i = 1; i <= keysexist; i++) {
    fits_read_keyn(f,i,keyname,value,com,&stat);
    if( stat == 0 && strncmp(keyname,"TTYPE",5) == 0 ) {

      char *i1 = index(value,'\'');
      char *i2 = rindex(value,'\'');
      if( i1 != NULL && i2 != NULL && i1 != i2 ) {
	int l = i2 - i1 - 1;
	for(int i = 0; i < l; i++)
	  value[i] = *(i1 + i + 1);
	value[l] = '\0';
      }
      cols.Add(value);
    }
  }
  fits_close_file(f,&stat);
  fits_report_error(stderr,stat);

  return cols;
}
