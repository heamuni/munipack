/*
   xmunipack - value units (intensity, flux, etc..) conversions

   Copyright © 1997-2013, 2019-22 F.Hroch (hroch@physics.muni.cz)

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

   This routine represents physical value intensity and provides
   its represenattion in various physical units.

   Currently only HST magnitude - photometry convention is used.

*/


#include "fits.h"
#include <wx/wx.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


using namespace std;

FitsValue::FitsValue(): type(UNIT_COUNT),hascal(false) {}

FitsValue::FitsValue(const FitsArray& a, const wxString& phsystemfile,
		     const wxString& f1, const wxString& f2,
		     const wxString& f3):
  type(UNIT_COUNT),hascal(false),array(a),fits_key_area(f1),
  fits_key_exptime(f2), fits_key_filter(f3)
{
  hascal = Init();
  if( hascal )
    Init_phsystem(phsystemfile);

  hascal = ! phconv.empty(); // photosystem is missing or an unknown
}

void FitsValue::Init_phsystem(const wxString& phsystemfile)
{
  phsystems = FitsPhotosystems(phsystemfile);
  if( phsystems.IsOk() ) {
    for(size_t i = 0; i < filter.size(); i++) {

      PhotoFilter phfilter = phsystems.GetFilter(photsys,filter);
      if( phfilter.IsOk() )
	phconv.push_back(PhotoConv(phfilter,area,exptime,scale));
    }
  }
}


bool FitsValue::Init()
{
  double sc11, sc12;

  area = GetKeyDouble(array,fits_key_area);
  exptime = GetKeyDouble(array,fits_key_exptime);
  sc11 = GetKeyDouble(array,"CD1_1");
  sc12 = GetKeyDouble(array,"CD1_2");
  scale = 3600.0*sqrt(sc11*sc11 + sc12*sc12);

  filter = array.GetKey(fits_key_filter);
  photsys = array.GetKey("PHOTSYS");

  return ! photsys.IsEmpty();
}

double FitsValue::GetKeyDouble(const FitsHdu& a, const wxString& key) const
{
  if( ! key.IsEmpty() ) {
    wxString s = a.GetKey(key);

    double x;
    if( !s.IsEmpty() &&  s.ToDouble(&x) )
      return x;
  }
  return 1.0;
}


void FitsValue::SetType(int t)
{
  if( hascal ) {
    if( UNIT_FIRST + 1 < t && t < UNIT_LAST )
      type = static_cast<units_type>(t);
    else
      type = UNIT_PHOTON;
  }
  else {
    if( UNIT_FIRST < t && t < UNIT_PHOTON )
      type = static_cast<units_type>(t);
    else
      type = UNIT_COUNT;
  }
}

void FitsValue::SetType(const wxString& a)
{
  for(int i = UNIT_FIRST+1; i < UNIT_LAST; i++)
    if( a == Label_str(i) ) {
      SetType(i);
      return;
    }
}


units_type FitsValue::GetType() const
{
  return type;
}

wxString FitsValue::GetName() const
{
  return Label_str(type);
}


wxString FitsValue::GetUnit() const
{
  return Units_str(type);
}

wxString FitsValue::ToString(double d, int k) const
{
  wxString a;

  if( type == UNIT_PHOTON || type == UNIT_COUNT ) {
    if( array.Bitpix() > 0 )
      a.Printf("%d",int(d+0.5));
    else
      a = HumanFormat(d);
  }

  if( hascal ) {

    if( type == UNIT_MAG )

      a.Printf("%.5g",phconv[k].GetMag(d));

    else if( type == UNIT_INTENSITY )

      a = HumanFormat(phconv[k].GetIntensity(d));

  }

  return a;
}

wxString FitsValue::Get_str(int i, int j) const
{
  wxASSERT(array.IsOk() && array.Naxis() == 2);
  if( 0 <= i && i < array.Width() && 0 <= j && j < array.Height() ) {
    double c = array.Pixel(i,j);
    return ToString(c,0);
  }
  return "";
}

vector<wxString> FitsValue::Get_str(int i, int j, const vector<int>& z) const
{
  wxASSERT(array.IsOk() && array.Naxis() == 3);

  vector<wxString> labels;
  if( 0 <= i && i < array.Width() && 0 <= j && j < array.Height() &&
      z.size() == 3 && array.IsColour() ) {

    float X = array.Pixel(i,j,2);
    float Y = array.Pixel(i,j,1);
    float Z = array.Pixel(i,j,0);

    float *d[3] = { &X, &Y, &Z };

    for(size_t i = 0; i < 3; i++) {
      double q = *d[i];
      wxString a;
      a.Printf(fabs(q) < 1e5 ? "% 7.1f" : "%.2e",q);
      labels.push_back(a);
    }
    return labels;
  }
  else {
    vector<wxString> labels = { "", "", ""};
    return labels;
  }
}

wxString FitsValue::Label_str(int n)
{
  switch(n){
  case PHQUANTITY_COUNT:     return "Count";
  case PHQUANTITY_PHOTON:    return "Photon";
  case PHQUANTITY_MAG:       return "Magnitude";
  case PHQUANTITY_INTENSITY: return "Intensity";
  default:                   return  "???";
  }
}

wxArrayString FitsValue::Label_str()
{
  wxArrayString a;

  for(int i = PHQUANTITY_FIRST+1; i < PHQUANTITY_LAST; i++)
    a.Add(Label_str(i));

  return a;
}

wxString FitsValue::Units_str(int n)
{
  switch(n){
  case UNIT_COUNT:     return  "";
  case UNIT_PHOTON:    return  "";
  case UNIT_MAG:       return L"mag/arcsec²";
  case UNIT_INTENSITY: return L"eV/s/m²/arcsec²";
  default:             return  "";
  }
}

wxArrayString FitsValue::Units_str()
{
  wxArrayString a;

  for(int i = UNIT_FIRST+1; i < UNIT_LAST; i++)
    a.Add(Units_str(i));

  return a;
}


// ----- Photometric Systems

const double PhotoConv::cspeed = 299792458.0;
const double PhotoConv::hplanck = 6.62606957e-34;
const double PhotoConv::evolt = 1.6021766208e-19;  // = 1eV
const double PhotoConv::sqrtpi2 = 1.2533141373155001; // = sqrt(pi/2)

PhotoConv::PhotoConv(const PhotoFilter& phfilter, double a, double e, double s):
  area(a), exptime(e), scale(s)
{
  flam = phfilter.flam;
  leff = phfilter.leff;
  lwidth = phfilter.lwidth;
}

double PhotoConv::intensity(double photon) const
{
  double phi = photon / (area * exptime * scale * scale);
  return phi * hplanck * cspeed / leff;
}

double PhotoConv::mag(double f, double f0) const
{
  if( f0 <= 0.0 || f <= 0.0) return 99.999;
  double x = f / f0;
  return -2.5*log10(x);
}

double PhotoConv::GetIntensity(double photon) const
{
  return intensity(photon) / evolt;
}

double PhotoConv::GetMag(double photon) const
{
  return mag(intensity(photon), sqrtpi2 * flam * (1e9*lwidth));
}


FitsPhotosystems::FitsPhotosystems(const wxString& name)
{
  fitsfile *f;
  int status = 0;
  float nullval = 0.0;
  int dummy,htype;
  char extname[FLEN_CARD];//,comment[FLEN_CARD];
  const char *cols[] = {"FILTER","LAM_EFF","LAM_FWHM","FLAM_REF"};
  int col[4];

  int nhdu = 0;

  // open file
  status = 0;
  fits_open_file(&f, name.fn_str(), READONLY, &status);
  fits_get_num_hdus(f,&nhdu,&status);
  if( status ) goto crash;

  for(int k = 1; k < nhdu; k++) { // the first HDU is skipped

    fits_movabs_hdu(f,k+1,&htype,&status);
    if( status ) goto crash;

    fits_read_key(f,TSTRING,"EXTNAME",extname,NULL,&status);

    //    fits_read_keyword(f,"EXTNAME",extname,comment,&status);
    // removing aphostrophes
    /*
    wxString e(extname);
    wxString e1 = e.SubString(1,e.Len()-2);
    e = e1.Trim();
    strcpy(extname,e.ToAscii());
    */

    if( htype == BINARY_TBL ) {

      long nrows;
      int ncols;
      fits_get_num_rows(f,&nrows,&status);
      fits_get_num_cols(f,&ncols,&status);
      if( status ) goto crash;

      for(int i = 0; i < 4; i++)
	fits_get_colnum(f,CASESEN,(char*)cols[i],&col[i],&status);

      double *leff = new double[nrows];
      double *lwidth = new double[nrows];
      double *flam = new double[nrows];

      long frow = 1, felem = 1;

      int width;
      fits_get_col_display_width(f,col[0],&width,&status);
      char **fs = new char*[nrows];
      for(int i = 0; i < nrows; i++)
	fs[i] = new char[width];
      fits_read_col(f,TSTRING,col[0], frow, felem, nrows, &nullval,fs,
		    &dummy, &status);
      fits_read_col(f,TDOUBLE,col[1], frow, felem, nrows, &nullval,leff,
		    &dummy, &status);
      fits_read_col(f,TDOUBLE,col[2], frow, felem, nrows, &nullval,lwidth,
		    &dummy,&status);
      fits_read_col(f,TDOUBLE,col[3], frow, felem, nrows, &nullval,flam,
		    &dummy,&status);
      if( status ) goto crash;

      vector<PhotoFilter> phfs;
      for(int i = 0; i < nrows; i++)
	phfs.push_back(PhotoFilter(fs[i],leff[i],lwidth[i],flam[i]));

      phsystems.push_back(Photosys(extname,phfs));

      delete[] leff;
      delete[] lwidth;
      delete[] flam;
      for(int i = 0; i < nrows; i++)
	delete[] fs[i];
      delete[] fs;

    }
  }

  fits_close_file(f, &status);

 crash:

  fits_report_error(stderr,status);
}

bool FitsPhotosystems::IsOk() const
{
  return ! phsystems.empty();
}

PhotoFilter FitsPhotosystems::GetFilter(const wxString& name,
					const wxString& filter) const
{
  for(size_t i = 0; i < phsystems.size(); i++) {
    if( phsystems[i].GetName() == name ) {
      return phsystems[i].GetFilter(filter);
    }
  }
  return PhotoFilter();
}


PhotoFilter Photosys::GetFilter(const wxString& filter) const
{
  for(size_t i = 0; i < filters.size(); i++) {
    if( filters[i].name == filter ) {
      return filters[i];
    }
  }
  return PhotoFilter();
}
