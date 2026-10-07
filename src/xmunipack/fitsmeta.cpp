/*

  xmunipack - fits meta-info


  Copyright © 2009-2011, 2019-2022 F.Hroch (hroch@physics.muni.cz)

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
#include <fitsio.h>
#include <math.h>
#include <vector>
#include <wx/wx.h>
#include <wx/uri.h>
#include <wx/filesys.h>
#include <wx/tokenzr.h>
#include <cstdio>
#include <cstring>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif


// --------  FitsMetaHdu

// by single HDU
FitsMetaHdu::FitsMetaHdu(const FitsHeader& head, int t):
  FitsHeader(head),htype(t),type(head.GetType(htype)),
  subtype(head.GetFlavour(htype,type)),nrows(-1),ncols(-1) {}

FitsMetaHdu::FitsMetaHdu(const FitsHeader& head, int t, int maxis,const long *maxes):
  FitsHeader(head),htype(t),type(head.GetType(htype)),
  subtype(head.GetFlavour(htype,type)), nrows(-1),ncols(-1)
{
  for(int i = 0; i < maxis; i++)
    naxes.push_back(maxes[i]);
}

FitsMetaHdu::FitsMetaHdu(const FitsHeader& head, int ht, long nr, int nc):
  FitsHeader(head),htype(ht), type(head.GetType(htype)),
  subtype(head.GetFlavour(htype,type)), nrows(nr),ncols(nc) {}

// by full FITS file
FitsMetaHdu::FitsMetaHdu(const FitsHdu& hdu, const wxImage& ico):
  FitsHeader(hdu.GetHeader()), htype(hdu.GetHduType()),
  type(hdu.Type()),subtype(hdu.Flavour()), nrows(0),ncols(0),icon(ico)
{
  for(size_t i = 0; i < hdu.GetCount(); i++)
    Add(hdu.Item(i));

  if( type == HDU_IMAGE ) {
    const FitsArray array(hdu);
    for(int i = 0; i < array.Naxis(); i++)
      naxes.push_back(array.Naxes(i));
  }
  else if( type == HDU_TABLE ) {
    const FitsTable table(hdu);
    ncols = table.Width();
    nrows = table.Height();
  }
}

size_t FitsMetaHdu::Naxis() const
{
  return naxes.size();
}

long FitsMetaHdu::Naxes(size_t n) const
{
  if( 0 <= n && n < naxes.size() )
    return naxes[n];
  else
    return 0;
}

std::vector<long> FitsMetaHdu::GetNaxes() const
{
  return naxes;
}

long FitsMetaHdu::Width() const
{
  return Naxes(0);
}

long FitsMetaHdu::Height() const
{
  return Naxes(1);
}

long FitsMetaHdu::Nrows() const
{
  return nrows;
}

int FitsMetaHdu::Ncols() const
{
  return ncols;
}

hdu_type FitsMetaHdu::Type() const
{
  return type;
}

hdu_flavour FitsMetaHdu::SubType() const
{
  return subtype;
}


wxString FitsMetaHdu::Type_str() const
{
  return GetType_str(type);
}

wxString FitsMetaHdu::SubType_str() const
{
  return GetFlavour_str(subtype);
}

wxImage FitsMetaHdu::GetIcon() const
{
  wxASSERT(icon.IsOk());
  return icon;
}

void FitsMetaHdu::SetIcon(const wxImage& i)
{
  wxASSERT(i.IsOk());
  icon = i;
}

wxString FitsMetaHdu::GetControlLabel() const
{
  wxString label = GetKey("EXTNAME");
  if( label.IsEmpty() )
    label = Type_str();
  return label;
}


// ------- FitsMeta


FitsMeta::FitsMeta():status(1),type(FITS_UNKNOWN),size(wxInvalidSize) {}

FitsMeta::FitsMeta(const wxString& name):
  url(name),status(1),type(FITS_UNKNOWN),size(wxInvalidSize)
{
  wxFileName fn(wxFileSystem::URLToFileName(url));
  size = fn.GetSize();
}

FitsMeta::FitsMeta(const FitsFile& fits, const wxImage& ico,
		   const std::vector<wxImage>& list):
  url(fits.GetURL()),status(0),type_str(fits.Type_str()),type(fits.Type()),
  icon(ico), size(wxInvalidSize)
{
  wxASSERT(fits.Status() && list.size() == fits.HduCount());

  for(size_t k = 0; k < fits.HduCount(); k++)
    hdu.push_back(FitsMetaHdu(fits.Hdu(k),list[k]));

  wxFileName fn(wxFileSystem::URLToFileName(url));
  size = fn.GetSize();
}


FitsMeta::FitsMeta(const wxString& u, const std::vector<FitsMetaHdu>& h, wxULongLong s):
  url(u),status(0),type(FITS_UNKNOWN),hdu(h), size(s)
{
  if( hdu.size() == 1 && hdu[0].Type() == HDU_IMAGE ) {

    if( hdu[0].SubType() == HDU_IMAGE_FRAME ) {
      // simple, grayscale image
      type = FITS_GREY;
      type_str = "Grey image";
    }
    else  if( hdu[0].SubType() == HDU_IMAGE_COLOUR ) {
      // colour image
      type = FITS_COLOUR;
      type_str = "Colour image";
    }
  }
  else {
    // more hdus
    type = FITS_MULTI;
    type_str = "Multi-extension";
  }
}


wxImage FitsMeta::GetIcon() const
{
  return icon;
}

void FitsMeta::SetIcon(const wxImage& i)
{
  wxASSERT(i.IsOk());
  icon = i;
}

void FitsMeta::SetIconList(const std::vector<wxImage>& ilist)
{
  wxASSERT(ilist.size() == hdu.size());
  for(size_t k = 0; k < hdu.size(); k++)
    hdu[k].SetIcon(ilist[k]);
}

size_t FitsMeta::HduCount() const
{
  return hdu.size();
}

FitsMetaHdu FitsMeta::Hdu(size_t n) const
{
  wxASSERT(0 <= n && n < hdu.size() );
  return hdu[n];
}

FitsMetaHdu *FitsMeta::GetHdu(size_t n)
{
  wxASSERT(0 <= n && n < hdu.size() );
  return &hdu[n];
}

int FitsMeta::Type() const
{
  return type;
}

wxString FitsMeta::Type_str() const
{
  return type_str;
}


wxString FitsMeta::Mtime() const
{
  wxURI uri(url);
  wxFileName name(wxURI::Unescape(uri.GetPath()));

  wxDateTime t = name.GetModificationTime();
  if( t.IsValid() )
    return t.FormatDate()+" "+t.FormatTime();
  else
    return wxEmptyString;
}


wxString FitsMeta::GetKeys(const wxString& key) const
{
  wxString a;
  for(size_t k = 0; k < hdu.size(); k++) {
    wxString l = hdu[k].GetKey(key);
    if( a.IsEmpty() )
      a = l;
    else
      a += ";" + l;
  }
  return a;
}

wxString FitsMeta::GetDateobs(const wxString& key) const
{
  return GetKeys(key);
}

wxString FitsMeta::GetDate(const wxString& key) const
{
  FitsTime dateobs(GetKeys(key));
  return dateobs.GetDate();
}

wxString FitsMeta::GetTime(const wxString& key) const
{
  FitsTime dateobs(GetKeys(key));
  return dateobs.GetTime();
}

wxString FitsMeta::GetFilter(const wxString& key) const
{
  return type == FITS_COLOUR ? "XYZ" : GetKeys(key);
}

wxString FitsMeta::GetExposure(const wxString& key) const
{
  const wchar_t *prefin[] = { L"",L"m",L"μ",L"n",L"p",L"f",L"a",L"z",L"y" };
  const wchar_t *prefix[] = { L"",L"k",L"M",L"G",L"T",L"P",L"E",L"Z",L"Y" };
  const wxString fmt = "%g%ss";

  wxString line;
  wxStringTokenizer tokenizer(GetKeys(key),";");
  while ( tokenizer.HasMoreTokens() ) {

    if( line != "" ) line += ";";
    wxString l = tokenizer.GetNextToken();
    double e;
    if( ! l.IsEmpty() && l.ToDouble(&e) ) {
      if( e > 0 ) {
	wxString a;
	double d = log10(e);
	int thorder = d >= 0 ?  d / 3 : (d - 1) / 3;
	double q = e / pow(10.0,3*thorder);
	wxString text;
	if( d >= 0 ) {
	  wxASSERT(0 <= thorder && thorder < 9);
	  text.Printf(fmt,q,prefix[thorder]);
	}
	else {
	  wxASSERT(0 <= -thorder && -thorder < 9);
	  text.Printf(fmt,q,prefin[-thorder]);
	}
	line += text;
      }
      else
	line += "0s";
    }
    else
      line += l;
  }
  return line;
}

wxString FitsMeta::GetURL() const
{
  return url;
}

void FitsMeta::SetURL(const wxString& a)
{
  url = a;
}

bool FitsMeta::IsOk() const
{
  return status == 0;
}


wxString FitsMeta::GetName() const
{
  wxURI uri(url);
  wxFileName name(uri.GetPath());
  return name.GetName();
}

wxString FitsMeta::GetPath() const
{
  wxURI uri(url);
  wxFileName name(uri.GetPath());
  return name.GetPath();
}

wxString FitsMeta::GetFullPath() const
{
  wxURI uri(url);
  wxFileName name(uri.GetPath());
  return name.GetFullPath();
}

wxString FitsMeta::GetFullName() const
{
  wxURI uri(url);
  wxFileName name(uri.GetPath());
  return name.GetFullName();
}

wxString FitsMeta::GetHumanReadableSize() const
{
  return wxFileName::GetHumanReadableSize(size);
}


wxULongLong FitsMeta::GetSize() const
{
  return size;
}

void FitsMeta::Clear()
{
  url.Clear();
  type_str.Clear();
  type = FITS_UNKNOWN;
  hdu.clear();
  icon = wxImage();
  size = wxInvalidSize;
}
