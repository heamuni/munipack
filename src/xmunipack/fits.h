/*

  XMunipack - FITS headers

  Copyright © 2009-2026 F.Hroch (hroch@physics.muni.cz)

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

#ifndef _FITS_H_
#define _FITS_H_

#include <wcslib.h>
#include <fitsio.h>
#include <wx/wx.h>
#include <vector>
#include <algorithm>
#include <cmath>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


#define FINDEXTNAME "FIND"
#define APEREXTNAME "APERPHOT"
#define PHOTOEXTNAME "PHOTOMETRY"
#define FITS_LABEL_X "X"
#define FITS_LABEL_Y "Y"
#define FITS_LABEL_APCOUNT "APCOUNT"
#define FITS_COL_PHOTONS "PHOTONS"
#define FITS_KEY_CSPACE "CSPACE"

enum fits_type {
  FITS_GREY, FITS_COLOUR, FITS_MULTI, FITS_3D,
  FITS_SCI, FITS_FLAT, FITS_DARK, FITS_BIAS, FITS_UNKNOWN
};

enum hdu_type {
  HDU_UNKNOWN, HDU_HEAD, HDU_IMAGE, HDU_TABLE
};

enum hdu_flavour {
  HDU_IMAGE_LINE, HDU_IMAGE_FRAME, HDU_IMAGE_CUBE, HDU_IMAGE_COLOUR,
  HDU_IMAGE_UNKNOWN,
  HDU_TABLE_BIN, HDU_TABLE_ASCII, HDU_TABLE_UNKNOWN,
  HDU_DUMMY
};

enum zoom_type {
  ZOOM_FIRST,
  ZOOM_SHRINK, ZOOM_ZOOM,
  ZOOM_LAST
};

enum itt_type {
  ITT_FIRST,
  ITT_LINE, ITT_ASINH, ITT_TANH, ITT_SNLIKE, ITT_SQR, ITT_PHOTO,
  ITT_LAST
};

enum tone_kind {
  TONE_KIND_FIRST,
  TONE_KIND_ABS, TONE_KIND_REL,
  TONE_KIND_LAST
};

enum colour_palette {
  PAL_FIRST,
  PAL_GREY, PAL_SEPIA, PAL_ROYAL,
  PAL_COLOUR, PAL_HIGHLIGHT, PAL_RAINBOW,
  PAL_MADNESS, PAL_VGA, PAL_AIPS0,
  PAL_LAST
};

enum colour_space {
  COLOURSPACE_FIRST,
  COLOURSPACE_XYZ,
  COLOURSPACE_LAST
};

enum colour_temperature {
  COLOURTEMP_D65
};

enum tune_operation {
  OP_TUNE_ZOOM=2, OP_TUNE_SCALE=4, OP_TUNE_ITT=8, OP_TUNE_PAL=16,
  OP_TUNE_COLOUR=32, OP_TUNE_RGB=64
};

enum units_type {
  UNIT_FIRST,
  UNIT_COUNT, UNIT_PHOTON, UNIT_MAG, UNIT_INTENSITY,
  UNIT_LAST
};

enum phquantity_type {
  PHQUANTITY_FIRST,
  PHQUANTITY_COUNT, PHQUANTITY_PHOTON, PHQUANTITY_MAG, PHQUANTITY_INTENSITY,
  PHQUANTITY_LAST
};

enum coords_type {
  COO_FIRST,
  COO_PIXEL, COO_EQDEG, COO_EQSIX,
  COO_LAST
};



// -- content of HDUs


class FitsHeader: public wxArrayString
{
public:
  static bool ParseRecord(const wxString&,wxString&,wxString&,wxString&);
  bool FindKey(const wxString&,wxString&, wxString&) const;
  bool FindVal(const wxString&,wxString&, wxString&) const;
  wxString GetKey(const wxString&) const;
  wxString GetUnit(const wxString&) const;
  int Bitpix() const;
  wxString Bitpix_str() const;
  wxString Exposure_str(const wxString&) const;
  virtual int GetNaxis() const;
  virtual hdu_type GetType(int) const;
  virtual hdu_flavour GetFlavour(int,hdu_type) const;
  static wxString GetType_str(hdu_type);
  static wxString GetFlavour_str(hdu_flavour);
  bool IsOk() const;
};


class FitsHdu: public wxObject
{
protected:
  FitsHeader header;
  int htype;
  hdu_type type;
  hdu_flavour flavour;
  bool modified;

public:
  FitsHdu();
  FitsHdu(const FitsHeader& h, int);
  FitsHdu(const FitsHeader& h, int, hdu_type, hdu_flavour);

  size_t GetCount() const;
  wxString Item(size_t i) const;

  virtual wxString GetKey(const wxString& a) const;
  virtual long GetKeyLong(const wxString& a) const;
  virtual double GetKeyDouble(const wxString& a) const;
  virtual wxString GetUnit(const wxString& a) const;
  virtual int Bitpix() const;
  virtual wxString Bitpix_str() const;
  virtual wxString Exposure_str(const wxString& a) const;
  virtual wxString GetExtname() const;
  virtual wxString GetControlLabel() const;
  virtual FitsHeader GetHeader() const;

  virtual int GetHduType() const;
  virtual hdu_type Type() const;
  virtual wxString Type_str() const;
  virtual hdu_flavour Flavour() const;
  virtual wxString Flavour_str() const;

  // templates and helpers for derived classes (Array + Table)
  virtual int Naxis() const;
  virtual long Naxes(int n) const;
  virtual long Width() const; // obsolete
  virtual long Height() const; // obsolete
  virtual long GetWidth() const;
  virtual long GetHeight() const;
  virtual long GetDepth() const;
  virtual bool IsOk() const;
  virtual bool IsColour() const;
  virtual bool IsModified() const { return modified; }
  virtual bool IsDisplayImplemented() const;

  virtual inline float Pixel(int) const { return 0; }
  virtual inline float Pixel(int, int) const { return 0; }
  virtual inline wxString Pixel_str(int, int) const;

  virtual bool GetWCS(double&,double&,double&,double&,double&,double&,
		      double&) const;
  virtual int hdr2str(int, char**, int, char**, int *, int *) const;
};


inline wxString FitsHdu::Pixel_str(int x, int y) const
{
  wxString a;
  a.Printf("%g",Pixel(x,y));
  return a;
}




class FitsArray: public FitsHdu
{
public:
  FitsArray() {}
  FitsArray(const FitsHeader&, long, long);
  FitsArray(const FitsHeader&, long, long, long);
  FitsArray(const FitsHeader&, int, int, long *, float *);
  FitsArray(const FitsHdu&);
  FitsArray(const FitsArray&);
  virtual ~FitsArray();

  FitsArray& operator = (const FitsArray&);
  bool operator == (const FitsArray&) const;
  bool operator != (const FitsArray&) const;

  int Naxis() const;
  long Naxes(const int n) const;
  long* Naxes() const;
  float Pixel(int) const;
  float Pixel(int, int) const;
  float Pixel(int, int, int) const;
  const float *PixelData() const;
  float *GetData() const;
  FitsArray GetGrid(int, int) const;
  FitsArray Plane(int) const;
  long Npixels() const;

  bool IsOk() const;

  void GetSubImage(int,int,FitsArray&) const;
  void SetSubImage(int,int,const FitsArray&);

  // scheduled for the accessor FitsImage
  bool HasWCS() const;

protected:

  wxObjectRefData* CreateRefData() const;
  wxObjectRefData *CloneRefData(const wxObjectRefData *) const;

};


class FitsArrayStat: public FitsArray
{
  bool minmax_initialised, medmad_initialised;
  int skip;
  float med, mad, xmin, xmax;

  void Setup_minimax();
  void Setup_medmad();

public:
  FitsArrayStat(const FitsArray&, int=0);

  void InitStat();

  float GetMed();
  float GetMad();
  float GetMin();
  float GetMax();
  static float QMed(long, float *, int);

};


class FitsTableColumn: public wxObject
{

protected:

  wxObjectRefData* CreateRefData() const;
  wxObjectRefData *CloneRefData(const wxObjectRefData *) const;

public:
  FitsTableColumn();
  FitsTableColumn(int, long, double *);
  FitsTableColumn(int, long, float *);
  FitsTableColumn(int, long, int *);
  FitsTableColumn(int, long, char **);
  FitsTableColumn(int, long, char *);
  FitsTableColumn(int, long, bool *);
  FitsTableColumn(int, long, short *);
  FitsTableColumn(int, long, long *);
  FitsTableColumn(int, long, long, long);

  //  FitsTableColumn(const FitsTableColumn&);
  //  FitsTableColumn& operator = (const FitsTableColumn&);
  FitsTableColumn Copy() const;
  virtual ~FitsTableColumn();

  int GetColType() const;
  long Nrows() const;
  float Cell(int x) const;
  wxString Cell_str(int x) const;

  const float *GetCol_float() const;
  const double *GetCol_double() const;
  const char **GetCol_char() const;
  const long *GetCol_long() const;

  char **GetData_string();
  char *GetData_char();
  bool *GetData_bool();
  int *GetData_int();
  short *GetData_short();
  long *GetData_long();
  float *GetData_float();
  double *GetData_double();

};


class FitsTable: public FitsHdu
{

protected:

  wxObjectRefData* CreateRefData() const;
  wxObjectRefData *CloneRefData(const wxObjectRefData *) const;

public:
  FitsTable() {}
  FitsTable(const FitsHeader&,int,const std::vector<FitsTableColumn>&);
  FitsTable(const FitsHdu&);
  FitsTable(const FitsTable&);
  FitsTable Copy() const;

  FitsTable& operator = (const FitsTable&);
  bool operator == (const FitsTable&) const;
  bool operator != (const FitsTable&) const;

  float Cell(int,int) const;
  wxString Cell_str(int, int) const;
  //  const float *GetCol_float(int) const;
  //  const char **GetCol_char(int) const;
  //  const int *GetCol_int(int) const;
  int GetColType(int) const;
  bool IsOk() const;
  FitsTableColumn GetColumn(int) const;
  FitsTableColumn GetColumn(const wxString&) const;
  int GetColIndex(const wxString&) const;
  wxArrayString GetColLabels() const;

  long Nrows() const;
  int Ncols() const;
  int Naxis() const { return 2; }
  long Naxes(int n) const;
  //  inline float Pixel(int, int) const;

  void GetStarChart(wxOutputStream&);

};



// ---  FITS file

class FitsFile
{
public:

  FitsFile();
  FitsFile(const wxString&);
  FitsFile(const FitsHdu&);
  FitsFile(const wxString&, const std::vector<FitsHdu>&);
  virtual ~FitsFile();
  void Clear();

  bool Status() const;
  size_t HduCount() const;
  int size() const;
  FitsHdu Hdu(size_t n) const;
  FitsHdu FindHdu(const wxString&) const;
  //  FitsArray GetArray(size_t n) const;
  //  FitsTable GetTable(size_t n) const;

  int Type() const;
  wxString Type_str() const;
  bool HasImage() const;
  bool HasFind() const;
  bool HasPhotometry() const;
  bool HasPhcal() const;
  wxString GetURL() const;

  bool IsOk() const;
  bool IsModified() const;
  wxString GetName() const;
  wxString GetFullName() const;
  wxString GetPath() const;
  wxString GetFullPath() const;

  wxArrayString GetErrorMessage() const;
  wxString GetErrorDescription() const;

  bool Save(const wxString&);

private:

  wxString filename;
  bool status;
  int type;
  std::vector<FitsHdu> hdu;
  wxArrayString errmsg;
  wxString smsg;

  void Recognise();
  int merge_head(fitsfile *, const FitsHdu&, int *) const;

};


// --- meta classes

class FitsMetaHdu: public FitsHeader
{
  int htype;
  hdu_type type;
  hdu_flavour subtype;
  long nrows;
  int ncols;
  std::vector<long> naxes;
  wxImage icon;

public:
  FitsMetaHdu(const FitsHeader&, int);
  FitsMetaHdu(const FitsHeader&, int, int, const long *);
  FitsMetaHdu(const FitsHeader&, int, long, int);
  FitsMetaHdu(const FitsHdu&, const wxImage&);

  size_t Naxis() const;
  long Naxes(size_t n) const;
  long Width() const;
  long Height() const;
  long Nrows() const;
  int Ncols() const;
  hdu_type Type() const;
  hdu_flavour SubType() const;
  wxImage GetIcon() const;
  void SetIcon(const wxImage&);

  wxString Type_str() const;
  wxString SubType_str() const;

  std::vector<long> GetNaxes() const;
  wxString GetControlLabel() const;

};

class FitsMeta
{
  wxString url;
  int status;
  wxString type_str;
  int type;
  std::vector<FitsMetaHdu> hdu;
  wxImage icon;
  wxULongLong size;

public:
  FitsMeta();
  FitsMeta(const wxString&);
  FitsMeta(const FitsFile&, const wxImage&, const std::vector<wxImage>&);
  FitsMeta(const wxString&, const std::vector<FitsMetaHdu>&, wxULongLong);
  void Clear();

  FitsMetaHdu Hdu(size_t) const;
  size_t HduCount() const;
  FitsMetaHdu *GetHdu(size_t);

  wxString Mtime() const;

  wxString GetURL() const;
  wxImage GetIcon() const;
  wxString GetKeys(const wxString &) const;
  wxString GetDateobs(const wxString &) const;
  wxString GetDate(const wxString &) const;
  wxString GetTime(const wxString &) const;
  wxString GetExposure(const wxString &) const;
  wxString GetFilter(const wxString &) const;

  int Type() const;
  wxString Type_str() const;

  bool IsOk() const;
  wxString GetName() const;
  wxString GetPath() const;
  wxULongLong GetSize() const;
  wxString GetFullName() const;
  wxString GetFullPath() const;
  wxString GetHumanReadableSize() const;

  void SetURL(const wxString&);
  void SetIcon(const wxImage&);
  void SetIconList(const std::vector<wxImage>&);

};

/*
class MuniFits: public FitsFile
{
  const FitsFile fitsfile;
  std::vector<FitsTone> tones;

public:
  MuniFits(const FitsFile&);

  const GetHeader(size_t) const;
  const GetArray(size_t) const;
  const GetTable(size_t) const;

};
*/

// image operations

class EmpiricalCDF
{
  long ncdf;
  float *xcdf, *ycdf;

 public:
  EmpiricalCDF(): ncdf(0),xcdf(0),ycdf(0) {}
  EmpiricalCDF(long, const float*);
  EmpiricalCDF(const EmpiricalCDF&);
  EmpiricalCDF& operator=(const EmpiricalCDF&);
  virtual ~EmpiricalCDF();

  bool IsOk() const;
  float GetQuantile(float) const;
  float GetInverse(float) const;
  void GetPositiveQuantile(float&, float&) const;
};


class FitsHisto: public wxObject
{
  double wbin,xmin,xmax;
  long cmax,cmin,ntotal;

  void Create(int, int, const float*);

protected:
  wxObjectRefData* CreateRefData() const;
  wxObjectRefData *CloneRefData(const wxObjectRefData *) const;

public:
  FitsHisto();
  FitsHisto(long, const float *d, double=0.999,int=0);
  FitsHisto(const FitsHisto&);
  FitsHisto& operator = (const FitsHisto&);

  int NBins() const;
  int Count(int n) const;
  float Freq(int n) const;
  float Bin(int n) const;
  float BinWidth() const;
  float BinMin() const;
  float BinMax() const;
  int MaxCount() const;
  int MinCount() const;
  float MaxFreq() const;
  float MinFreq() const;
  float MinVal() const;
  float MaxVal() const;
  long Ntotal() const;

  bool IsOk() const;

};


class FitsGeometry: public FitsArray
{
  float MeanLine(int, int) const;
  float MeanRect(int, int, int, int) const;
  float MeanRect_debug(int, int, int, int) const;
  float MeanRect(const float *, long, long, int, int, int, int) const;
  float MeanSquare(const float *, int, int, int, int, int) const;

 public:
  FitsGeometry(const FitsArray&);

  void ShrinkSubwin(int, int,int,int,int,int,int,float*) const;
  FitsArray GetSubArray(int,int,int,int);
  void SetSubArray(int,int,const FitsArray&);

};


class FitsTone
{
  bool initialised;
  float black, qblack, sense, rsense, refsense, refblack, refqblack;
  EmpiricalCDF cdf_back;

 public:

  FitsTone();
  FitsTone(const FitsArray&);
  FitsTone(const FitsTone&);
  FitsTone& operator = (const FitsTone&);
  FitsTone(long, const float *);
  bool IsOk() const;
  void Setup(long, const float *,bool);
  inline float Scale(float) const;
  float *Scale(long, const float *) const;
  void SetBlack(float);
  void SetQblack(float);
  void SetSense(float);
  void SetRsense(float);
  float GetSense() const { return sense; }
  float GetBlack() const { return black; }
  float GetQblack() const { return qblack; }
  float GetRsense() const { return rsense; }
  void Reset();

};

inline float FitsTone::Scale(float x) const
{
  return (x - black)/sense;
}


// -- Intensity transfer table (ITT)
class FitsItt
{
  int itt;
  float (FitsItt::*F)(float) const;
  float itt_line(float) const;
  float itt_snlike(float) const;
  float itt_tanh(float) const;
  float itt_asinh(float) const;
  float itt_photo(float) const;
  float itt_sqr(float) const;

public:
  FitsItt();
  FitsItt(int);

  void SetItt(int);
  void SetItt(const wxString&);
  int GetItt() const;
  wxString GetItt_str() const;
  inline float Scale(float) const;
  float *Scale(long, const float *) const;

  static wxString Type_str(int);
  static wxArrayString Type_str();

  bool IsLinear() const;
  float InvScale(float) const;

};

inline float FitsItt::Scale(float a) const
{
  return (this->*F)(a);
}

class FitsColour
{
  wxString cspace;
  int ispace;
  float saturation;
  bool nitevision;
  float meso_level, meso_width;
  int ncolours, nbands;
  float *trafo,*level,*weight;
  float Sn,uw,vw;

public:
  FitsColour();
  FitsColour(const wxString&,const FitsArray&);
  FitsColour(const FitsColour&);
  FitsColour& operator = (const FitsColour&);
  virtual ~FitsColour();

  void Init(const wxString&,const FitsArray&);

  void SetWhitePoint(float,float);
  void SetSaturation(float);
  void SetNiteVision(bool);
  void SetMesoLevel(float);
  void SetMesoWidth(float);
  void SetTrans(const wxString&);
  void SetTrans(const wxString&,const wxString&);
  void SetTrans(int,int);
  void SetTrans(int,int,float);
  void SetLevel(int,float);
  void SetWeight(int,float);
  float GetWeight(int) const;
  float GetLevel(int) const;
  float GetTrans(int,int) const;
  int GetColours() const;
  int GetBands() const;
  inline void Luv_XYZ(float,float,float,float&, float&, float&) const;
  inline void XYZ_Luv(float,float,float,float&, float&, float&) const;

  void Instr_XYZ(long,size_t,const float **,float*,float*,float*);
  void TuneColours(long,float*,float*,float*) const;
  void TuneColours(float&,float&,float&) const;

  float Scotopic(float,float,float) const;

  float GetWhitePointX() const;
  float GetWhitePointY() const;
  float GetSaturation() const { return saturation; }
  bool GetNiteVision() const { return nitevision; }
  float GetMesoLevel() const { return meso_level; }
  float GetMesoWidth() const { return meso_width; }
  void NiteVision(float,float,float&,float&,float&) const;
  wxString GetColourspace() const;

  static wxString Type_str(int);
  static wxArrayString Type_str();

};

inline void FitsColour::XYZ_Luv(float X, float Y, float Z,
			 float& L, float& u, float& v) const
{
  float s,u1,v1,y,t;

  s = X + 15.0f*Y + 3.0f*Z;
  if( s > 0.0f ) {
    u1 = 4.0f*X/s;
    v1 = 9.0f*Y/s;
  }
  else {
    u1 = uw;
    v1 = vw;
  }

  y = Y / 100.0f;
  if( y > 0.0088565f /* = powf(6.0/29.0,3)*/ )
    L = 116.0f*cbrtf(y) - 16.0f;
  else
    L = 903.30f*y;   /* = pow(26.0/3.0,3) */

  t = 13.0f*L;
  u = t*(u1 - uw);
  v = t*(v1 - vw);
}

inline void FitsColour::Luv_XYZ(float L, float u, float v,
			 float& X, float& Y, float& Z) const
{
  float u1,v1,s,t,w;

  s = 13.0f*L;
  if( s > 0.0f ) {
    u1 = u/s + uw;
    v1 = v/s + vw;
  }
  else {
    u1 = uw;
    v1 = vw;
  }

  if( L > 8.0f ) {
    w = (L + 16.0f)/116.0f;
    Y = 100.0f*w*w*w;
  }
  else
    Y = L*0.11071f; /* = 100*powf(3.0/29.0,3) */

  t = Y / (4.0f*v1);
  X = t*(9.0f*u1);
  Z = t*(12.0f - 3.0f*u1 - 20.0f*v1);
}


// ---  FitsPalette

class FitsPalette
{
  int pal;
  bool inverse;
  int npal, npal1;
  unsigned char *rpal,*gpal,*bpal, *pals;

  void Create_Grey();
  void Create_Sepia();
  void Create_Royal();
  void Create_Colour();
  void Create_Highlight();
  void Create_Rainbow();
  void Create_Madness();
  void Create_VGA();
  void Create_AIPS0();
  void CreatePalette();
  void HSL_RGB(float, float, float, unsigned char&, unsigned char&, unsigned char&);
  float f(int,float,float,float) const;

public:
  FitsPalette(int =PAL_GREY, bool =false);
  FitsPalette(const FitsPalette&);
  FitsPalette& operator = (const FitsPalette&);
  virtual ~FitsPalette();

  void SetPalette(int);
  void SetPalette(const wxString&);
  void SetInverse(bool=false);
  int GetColours() const;
  int GetPalette() const;
  wxString GetPalette_str() const;
  bool GetInverse() const;

  unsigned char R(int i) const;
  unsigned char G(int i) const;
  unsigned char B(int i) const;
  inline void RGB(float, unsigned char&, unsigned char&, unsigned char&) const;
  inline void RGB(float, unsigned char *) const;
  unsigned char *RGB(long,float *);

  static wxString Type_str(int);
  static wxArrayString Type_str();

};


inline void FitsPalette::RGB(float f, unsigned char& r, unsigned char& g,
			     unsigned char& b) const
{
  int l = std::max(std::min(int(255.0*f + 0.5f),npal-1),0);
 // wxASSERT(0 <= l && l < npal);
  r = rpal[l];
  g = gpal[l];
  b = bpal[l];
}

inline void FitsPalette::RGB(float g, unsigned char *rgb) const
{
  int l = npal1 * g;
  l = l < 0 ? 0 : l;
  l = l > npal1 ? npal1 : l;
  memcpy(rgb,pals+3*l,3);
}



// projections

class FitsProjection
{
  wxString type;
  double acen, dcen, xcen, ycen, scale, angle, reflex;

  void gnomon(double,double,double&,double&) const;
  void ignomon(double,double,double&,double&) const;

public:
  FitsProjection();
  FitsProjection(const wxString& t, double a,double d,double x,double y,
		 double s,double r, double z);
  FitsProjection(const wxString& t, double a,double d, double x, double y,
		 double cd11, double cd12, double cd21, double cd22);

  void ad2xy(double,double,double&,double&) const;
  void xy2ad(double,double,double&,double&) const;
  double xy2dist(double,double,double,double) const;
  double GetScale() const;
  bool IsOk() const;
};


// -- FitsCoo

class FitsCoo: public FitsArray
{
  struct wcsprm *wcs;
  int nwcs;
  coords_type type;
  bool haswcs;
  int digits;

public:
  FitsCoo();
  FitsCoo(const FitsArray&);
  FitsCoo(const FitsCoo&);
  virtual ~FitsCoo();

  FitsCoo& operator = (const FitsCoo&);

  void SetType(int);
  void SetType(const wxString&);
  coords_type GetType() const;

  void GetEq(int,int,double&,double&) const;
  void RaSix(double,int&,int&,double&) const;
  void DecSix(double,char&,int&,int&,double&) const;

  void GetPix(int,int,wxString&,wxString&) const;
  void GetCoo(int,int,wxString&,wxString&) const;

  bool HasWCS() const;

  static wxString Label_str(int);
  static wxArrayString Label_str();
};


// Photometric Systems

class PhotoFilter
{
public:
 PhotoFilter(): leff(0.0), lwidth(0.0), flam(-1.0) {}
 PhotoFilter(const wxString& n, double le, double lw, double f0):
  name(n),leff(le),lwidth(lw),flam(f0) {}

  bool IsOk() const { return flam > 0.0; }

  wxString name;
  double leff, lwidth, flam;

};

class Photosys
{
  wxString name;
  std::vector<PhotoFilter> filters;

public:
  Photosys(const wxString& n, const std::vector<PhotoFilter>& f):
    name(n), filters(f) {}

  PhotoFilter GetFilter(const wxString&) const;
  wxString GetName() const { return name; }
};

class FitsPhotosystems
{
  std::vector<Photosys> phsystems;

public:
  FitsPhotosystems() {}
  FitsPhotosystems(const wxString&);

  bool IsOk() const;
  PhotoFilter GetFilter(const wxString&, const wxString&) const;

};

class PhotoConv
{
  double flam, leff, lwidth, area, exptime, scale;
  static const double cspeed, hplanck, evolt, sqrtpi2, STspflux;

  double mag(double,double) const;
  double intensity(double) const;

 public:
  PhotoConv() {}
  PhotoConv(const PhotoFilter&, double, double, double);

  double GetIntensity(double) const;
  double GetMag(double) const;

};


class FitsValue
{
public:
  FitsValue();
  FitsValue(const FitsArray&, const wxString&, const wxString&,
	    const wxString&, const wxString&);

  void SetType(int);
  void SetType(const wxString&);
  units_type GetType() const;
  wxString GetName() const;
  wxString GetUnit() const;

  wxString Get_str(int,int) const;
  std::vector<wxString> Get_str(int,int,const std::vector<int>&) const;

  static wxString Label_str(int);
  static wxArrayString Label_str();
  static wxString Units_str(int);
  static wxArrayString Units_str();

  bool HasPhcal() const { return hascal; }

private:

  units_type type;
  bool hascal;
  FitsArray array;
  double area, exptime, scale;
  wxString filter;
  wxString photsys, fits_key_area, fits_key_exptime, fits_key_filter;
  FitsPhotosystems phsystems;
  std::vector<PhotoConv> phconv;

  bool Init();
  void Init_phsystem(const wxString&);
  wxString ToString(double,int) const;
  double GetKeyDouble(const FitsHdu&, const wxString&) const;

};


class FitsTime
{
  long year, month, day, hour, minute, second, milisecond;
  wxString date,time;

public:
  FitsTime(const wxString&);

  void SetDateTime(const wxString& );
  void SetDate(const wxString& );
  void SetTime(const wxString& );

  double GetJd() const;
  wxString GetDate() const { return date; }
  wxString GetTime() const { return time; }

};



// auxiliary  functions

bool FitsCopyFile(const wxString&,const wxString&);
bool FitsCopyHdu(const wxString&,const wxString&,const wxString&);
//wxArrayString FitsColumns(const wxString&);

wxString HumanFormat(double);


#endif
