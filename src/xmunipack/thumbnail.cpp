/*

  xmunipack - save/load thumbnails

  Copyright © 2009-2011 F.Hroch (hroch@physics.muni.cz)

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

  Note.

  Thumbnails can be stored as files in a filesystem and also serialized
  to a stream.


  In filesystem storage, the files are: XML description with *.fm suffix,
  and a set if png files cointainig icons. Pointers to icon names are
  in XML file. Example:

   thumb_23.fm    .. the xml part (description)
   thumb_23.png   .. the global icon
   thumb_23.0.png .. the first hdu icon


  In stream version, the stream for a single FitsMeta has the structure:

  15781    .. the length of XML part in bytes
  <XML>
    ...    .. the xml part
  </XML>
  20500    .. the length of global icon in bytes
  [PNG]    .. the global icon
  20500    .. the length of first hdu icon in bytes
  [PNG]    .. the icon
   ...

  The length of parts are stored by wxWidget in an architecture independed
  way as is described in documentation.

*/


#include "xmunipack.h"
#include <wx/wx.h>
#include <wx/filename.h>
#include <wx/utils.h>
#include <wx/xml/xml.h>
#include <wx/filesys.h>
#include <wx/tokenzr.h>
#include <wx/txtstrm.h>
#include <wx/sstream.h>
#include <wx/datstrm.h>
#include <wx/mstream.h>
#include <vector>

using namespace std;


// --- xMemoryInputStream

class xMemoryInputStream
{
public:
  xMemoryInputStream(wxInputStream& s): size(0),data(0),stream(s),mstream(0)
  {
    wxDataInputStream ds(stream);
    ds >> size;
    if( stream.Eof() ) return;
    data = new char[size];
    stream.Read(data,size);
    mstream = new wxMemoryInputStream(data,size);
  }
  virtual ~xMemoryInputStream() { delete[] data; delete mstream; }
  wxMemoryInputStream& Cutout() { return *mstream; }
  bool IsOk() const { return mstream; }

private:
  unsigned int size;
  char *data;
  wxInputStream& stream;
  wxMemoryInputStream *mstream;
};


// ----- MuniThumbnail

MuniThumbnail::MuniThumbnail(const FitsMeta& fm)
{
  SetMeta(fm);
}

MuniThumbnail::MuniThumbnail(const wxString& filename)
{
  Load(filename);
}

MuniThumbnail::MuniThumbnail(wxInputStream& stream)
{
  Load(stream);
}



void MuniThumbnail::Load(wxInputStream& stream)
{
  xMemoryInputStream ms(stream);
  if( ! ms.IsOk() ) return;

  /*
  wxDataInputStream ds(stream);
  unsigned int xs;
  ds >> xs;
  if( stream.Eof() ) return;
  char *b = new char[xs];
  stream.Read(b,xs);
  wxMemoryInputStream ms(b,xs);
  */

  wxXmlDocument xml;
  xml.Load(ms.Cutout()/*stream*/);

  if( xml.IsOk() ) {
    meta = ParseXML(xml);
    wxASSERT(meta.IsOk());
    //    wxLogDebug(url);
    /*
    ds >> xs;
    char *bi = new char[xs];
    stream.Read(bi,xs);
    wxMemoryInputStream mi(bi,xs);
    */
    xMemoryInputStream mi(stream);

    //    /*wxImage*/ icon.LoadFile(/*stream*/mi,wxBITMAP_TYPE_PNG);
    wxImage i(mi.Cutout(),wxBITMAP_TYPE_PNG);
    wxASSERT(i.IsOk());
    meta.SetIcon(i);

    //for(vector<MuniThumbnailHdu>::iterator h =hdu.begin(); h !=hdu.end();++h){
    //    vector<FitsMetaHdu> hdu = meta.GetHdu();
    //    vector<FitsMetaHdu> *hdu = meta.GetHduX();
    //    for(vector<FitsMetaHdu>::iterator h =hdu->begin();h !=hdu->end(); ++h){
    for(size_t k = 0; k < meta.HduCount(); k++) {
      FitsMetaHdu *h = meta.GetHdu(k);
      /*
      ds >> xs;
      char *bh = new char[xs];
      stream.Read(bh,xs);
      wxMemoryInputStream mh(bh,xs);
      */
      xMemoryInputStream mh(stream);

      wxImage i(mh.Cutout()/*stream*/,wxBITMAP_TYPE_PNG);
      wxASSERT(i.IsOk());
      h->SetIcon(i);
      //      delete[] bh;
    }
    //    delete[] bi;
  }

  //  delete[] b;

  wxASSERT(meta.IsOk());

  /*
  wxXmlDocument thumb;
  if( thumb.Load(stream) )
      LoadXML(thumb);
  */
}


void MuniThumbnail::Load(const wxString& filename)
{
  if( ! wxFileName::FileExists(filename) ) return;

  wxXmlDocument xml;
  xml.Load(filename);
  if( xml.IsOk() ) {
    meta = ParseXML(xml);

    wxFileName f(wxFileSystem::URLToFileName(icon/*url*/));
    wxImage i(f.GetFullPath());
    meta.SetIcon(i);
    //    wxLogDebug(f.GetFullPath());

    //for(vector<MuniThumbnailHdu>::iterator h =hdu.begin();h !=hdu.end(); ++h){
    //    vector<FitsMetaHdu> *hdu = meta.GetHduX();
    //    for(vector<FitsMetaHdu>::iterator h =hdu->begin();h !=hdu->end(); ++h){
      //      wxFileName f(wxFileSystem::URLToFileName(h->ImageURL()));
    for(size_t k = 0; k < meta.HduCount(); k++) {
      FitsMetaHdu *h = meta.GetHdu(k);
      wxFileName f(wxFileSystem::URLToFileName(icons[k/*h-hdu->begin()*/]));
      wxImage i(f.GetFullPath());
      h->SetIcon(i);
    }
  }

  //  wxASSERT(url.IsEmpty());
  //  wxLogDebug(_("url: ")+url);
  /*
  if( wxFileName::FileExists(filename) ) {
    wxXmlDocument thumb(filename);
    if( thumb.IsOk() )
      LoadXML(thumb);
  }
  */
}

void MuniThumbnail::Save(wxOutputStream& stream)
{
  //  wxString buf(' ',80);

  wxCountingOutputStream b;
  //  wxXmlDocument thumb = CreateXML(wxEmptyString);
  //  thumb.Save(b);
  wxXmlDocument xml = CreateXML(meta,wxEmptyString);
  xml.Save(b);
  wxLogDebug(_("MuniThumbnail::Save %d"),(int) b.GetSize());
  wxDataOutputStream x(stream);
  x << (unsigned int) b.GetSize();
  //  stream << "XML " << b.GetSize() << "\n";
  xml.Save(stream);

  wxCountingOutputStream bi;
  wxImage i = meta.GetIcon();
  wxASSERT(i.IsOk());

  i.SaveFile(bi,wxBITMAP_TYPE_PNG);
  //  stream << "PNG " << bi.GetSize() << "\n";
  x << (unsigned int) bi.GetSize();
  i.SaveFile(stream,wxBITMAP_TYPE_PNG);
  //  icon.SaveFile(stream,wxBITMAP_TYPE_PNG);

  for(size_t k = 0; k < meta.HduCount(); k++) {
    FitsMetaHdu hdu = meta.Hdu(k);
    wxImage i = hdu.GetIcon();
    wxASSERT(i.IsOk());

    wxCountingOutputStream bb;
    i.SaveFile(bb,wxBITMAP_TYPE_PNG);
    x << (unsigned int) bb.GetSize();
    i.SaveFile(stream,wxBITMAP_TYPE_PNG);
  }


//   for(vector<MuniThumbnailHdu>::const_iterator h=hdu.begin();h!=hdu.end();++h){
//     wxImage i(h->GetImage());
//     wxCountingOutputStream bb;
//     i.SaveFile(bb,wxBITMAP_TYPE_PNG);
//     //    stream << "PNG " << bb.GetSize() << "\n";
//     x << (unsigned int) bb.GetSize();
//     i.SaveFile(stream,wxBITMAP_TYPE_PNG);
//   }
  //  return b.GetSize();
}

bool MuniThumbnail::IsOk() const
{
  //  wxLogDebug(_("%d %d %d"),icon.IsOk(),url.IsEmpty(),hdu.size());
  //  return icon.IsOk() && ! url.IsEmpty() && hdu.size() > 0;
  return meta.GetIcon().IsOk() && ! meta.GetURL().IsEmpty() &&
    meta.HduCount() > 0 && meta.IsOk();
}

wxString MuniThumbnail::GetURL() const
{
  return meta.GetURL();
  //  return url;
}

FitsMeta MuniThumbnail::ParseXML(const wxXmlDocument& thumb)
{
  if( thumb.IsOk() && thumb.GetRoot()->GetName() == "munipack" ){

    std::vector<FitsMetaHdu> hdu;
    wxString url, type_str;
    wxULongLong size = wxInvalidSize;

    // clean shared
    icon.Clear();
    icons.Clear();

    wxXmlNode *child = thumb.GetRoot()->GetChildren();
    while (child) {

      if (child->GetName() == "fitsmeta") {

	//	wxString xiconame;

	wxXmlAttribute *prop = child->GetAttributes();
	while(prop) {

	  if( prop->GetName() == "url" )
	    url = prop->GetValue();
	  else if( prop->GetName() == "type" ) {
	    type_str = prop->GetValue();
	    /*
	    if( type_str == wxT("Gray image") )
	      type = FITS_GRAY;
	    else if( type_str == wxT("Color image") )
	      type = FITS_COLOR;
	    else if( type_str == wxT("Multi layer") )
	      type = FITS_MULTI;
	    */
	  }
	  else if( prop->GetName() == "icon" ) {
	    icon = prop->GetValue();
	    //	    iconurl = prop->GetValue();
	    /*
	    wxFileName f(wxFileSystem::URLToFileName(prop->GetValue()));
	    xiconame = f.GetFullPath();
	    */
	  }
	  else if( prop->GetName() == "size" ) {
	    wxString x = prop->GetValue();
	    // this conversions doesn't works for ULongLong under 64-bit
	    // compiler and system
	    //	    wxULongLong_t n;
	    //	    if( x.ToULongLong(&n) )
	    unsigned long n;
	    if( x.ToULong(&n) )
	      size = n;
	  }
	  prop = prop->GetNext();
	}

	// wxLogDebug(_("MuniThumbnail::LoadFile iconname ")+xiconame);
	/*
	wxImage i;
	if( ! xiconame.IsEmpty() && i.LoadFile(xiconame) && i.IsOk() )
	  icon = i;
	*/

	wxXmlNode *dchild = child->GetChildren();
	while (dchild) {

	  if( dchild->GetName() == "hdu") {

	    wxArrayString head;
	    //	    wxString imageurl;
	    //	    wxImage image;
	    wxString htype,hstype;
	    int ncols = 0;
	    long nrows = 0;
	    vector<long> naxes;


	    //	    wxXmlProperty *prop = dchild->GetProperties();
	    wxXmlAttribute *prop = dchild->GetAttributes();
	    while(prop) {

	      if( prop->GetName() == "type" ) {
		htype = prop->GetValue();
	      }
	      else if( prop->GetName() == "subtype" ) {
		hstype = prop->GetValue();
	      }
	      else if( prop->GetName() == "ncols" ) {
		wxString x = prop->GetValue();
		long n;
		if( x.ToLong(&n) )
		  ncols = n;
	      }
	      else if( prop->GetName() == "nrows" ) {
		wxString x = prop->GetValue();
		long n;
		if( x.ToLong(&n) )
		  nrows = n;
	      }
	      else if( prop->GetName() == "icon" ) {
		/*
		wxFileName f(wxFileSystem::URLToFileName(prop->GetValue()));
		wxString a(f.GetFullPath());
		wxImage i;
		*/
		//		  wxLogDebug(_("MuniThumbnail::LoadFile ")+a);
		/*
		if( ! a.IsEmpty() && i.LoadFile(a) && i.IsOk() )
		  image = i;
		*/
		//		imageurl = prop->GetValue();
		icons.Add(prop->GetValue());
	      }
	      prop = prop->GetNext();
	    }

	    wxXmlNode *hchild = dchild->GetChildren();
	    while (hchild) {

	      if( hchild->GetName() == "naxes" ) {

		wxString b = hchild->GetNodeContent();

		wxStringInputStream buf(b);
		wxTextInputStream xbuf(buf);
		while(! buf.Eof()) {
		  int n;
		  xbuf >> n;
		  naxes.push_back(long(n));
		}
	      }
	      else if( hchild->GetName() == "fitshead") {

		wxXmlNode *hchildx = hchild->GetChildren();
		if( hchildx->GetName() == "cdata") {
		  wxStringTokenizer t(hchildx->GetContent(),"\r\n");
		  while(t.HasMoreTokens() )
		    head.Add(t.GetNextToken());
		}
	      }

	      hchild = hchild->GetNext();
	    }

	    FitsMetaHdu h(head,wxImage(),ncols,nrows,naxes,htype,hstype);
	    //MuniThumbnailHdu h(head,imageurl,ncols,nrows,naxes,htype,hstype);
	    hdu.push_back(h);

	  }
	  dchild = dchild->GetNext();
	}
      }
      child = child->GetNext();
    }
    return FitsMeta(url,type_str,hdu,wxImage(),size);
  }





//   if( wxFileName::FileExists(filename) ) {
//     if( thumb.Load(filename) && thumb.GetRoot()->GetName() == wxT("munipack") ){

//       wxXmlNode *child = thumb.GetRoot()->GetChildren();
//       while (child) {

// 	if (child->GetName() == wxT("fitsmeta")) {

// 	  wxString xiconame;

// 	  wxXmlProperty *prop = child->GetProperties();
// 	  while(prop) {

// 	    if( prop->GetName() == wxT("url") )
// 	      url = prop->GetValue();
// 	    else if( prop->GetName() == wxT("type") ) {
// 	      type_str = prop->GetValue();
// 	      if( type_str == wxT("Gray image") )
// 		type = FITS_GRAY;
// 	      else if( type_str == wxT("Color image") )
// 		type = FITS_COLOR;
// 	      else if( type_str == wxT("Multi layer") )
// 		type = FITS_MULTI;
// 	    }
// 	    else if( prop->GetName() == wxT("icon") ) {
// 	      wxFileName f(wxFileSystem::URLToFileName(prop->GetValue()));
// 	      xiconame = f.GetFullPath();
// 	    }

// 	    prop = prop->GetNext();
// 	  }

// 	  wxImage i;
// 	  // wxLogDebug(_("MuniThumbnail::LoadFile iconname ")+xiconame);
// 	  if( ! xiconame.IsEmpty() && i.LoadFile(xiconame) && i.IsOk() )
// 	    icon = i;

// 	  wxXmlNode *dchild = child->GetChildren();
// 	  while (dchild) {

// 	    if( dchild->GetName() == wxT("hdu")) {

// 	      wxArrayString head;
// 	      wxImage image;
// 	      wxString htype,hstype;
// 	      int ncols = 0;
// 	      long nrows = 0;
// 	      vector<long> naxes;

// 	      wxXmlProperty *prop = dchild->GetProperties();
// 	      while(prop) {

// 		if( prop->GetName() == wxT("type") ) {
// 		  htype = prop->GetValue();
// 		}
// 		else if( prop->GetName() == wxT("subtype") ) {
// 		  hstype = prop->GetValue();
// 		}
// 		else if( prop->GetName() == wxT("ncols") ) {
// 		  wxString x = prop->GetValue();
// 		  long n;
// 		  if( x.ToLong(&n) )
// 		    ncols = n;
// 		}
// 		else if( prop->GetName() == wxT("nrows") ) {
// 		  wxString x = prop->GetValue();
// 		  long n;
// 		  if( x.ToLong(&n) )
// 		    nrows = n;
// 		}
// 		else if( prop->GetName() == wxT("icon") ) {
// 		  wxFileName f(wxFileSystem::URLToFileName(prop->GetValue()));
// 		  wxString a(f.GetFullPath());
// 		  wxImage i;
// 		  //		  wxLogDebug(_("MuniThumbnail::LoadFile ")+a);
// 		  if( ! a.IsEmpty() && i.LoadFile(a) && i.IsOk() )
// 		    image = i;
// 		}
// 		prop = prop->GetNext();
// 	      }

// 	      wxXmlNode *hchild = dchild->GetChildren();
// 	      while (hchild) {

// 		if( hchild->GetName() == wxT("naxes") ) {

// 		  wxString b = hchild->GetNodeContent();

// 		  wxStringInputStream buf(b);
// 		  wxTextInputStream xbuf(buf);
// 		  while(! buf.Eof()) {
// 		    int n;
// 		    xbuf >> n;
// 		    naxes.push_back(long(n));
// 		  }
// 		}
// 		else if( hchild->GetName() == wxT("fitshead")) {

// 		  wxXmlNode *hchildx = dchild->GetChildren();
// 		  if( hchildx->GetName() == wxT("cdata")) {
// 		    wxStringTokenizer t(hchildx->GetContent(),wxT("\r\n"));
// 		    while(t.HasMoreTokens() )
// 		      head.Add(t.GetNextToken());
// 		  }
// 		}

// 		hchild = hchild->GetNext();
// 	      }

// 	      FitsMetaHdu h(head,image,ncols,nrows,naxes,htype,hstype);
// 	      hdu.push_back(h);

// 	    }
// 	    dchild = dchild->GetNext();
// 	  }
// 	}
// 	child = child->GetNext();
//       }
//       return;
//     }
//   }


  return FitsMeta();
 //wxLogDebug(wxT("MuniArchive::LoadMeta !!!!!!!!!!! Missing file: ")+filename);
}



wxString MuniThumbnail::CreateIconame(const wxString& filename,
				      const wxString& suf)
{
  wxFileName xicon(filename);
  if( suf.IsEmpty() )
    xicon.SetExt("png");
  else
    xicon.SetExt(suf+".png");
  return xicon.GetFullPath();
}


wxXmlDocument MuniThumbnail::CreateXML(const FitsMeta& fitsmeta,
				       const wxString& filename)
{
  wxXmlDocument thumb;

  wxXmlNode *root = new wxXmlNode(0,wxXML_ELEMENT_NODE,"munipack");
  root->AddAttribute("version",VERSION);
  thumb.SetRoot(root);

  wxXmlNode *meta = new wxXmlNode(0,wxXML_ELEMENT_NODE,"fitsmeta");
  wxString x;
  meta->AddAttribute("url",/*url*//*GetURL()*/fitsmeta.GetURL());
  meta->AddAttribute("type",/*type_str*//*Type_str()*/fitsmeta.Type_str());
  /*
  wxFileName xicon(filename);
  xicon.SetExt(wxT("png"));
  wxString xiconame = xicon.GetFullPath();
  wxImage icon(GetIcon());
  icon.SaveFile(xiconame);
  */
  if( !filename.IsEmpty() ) {
    wxASSERT(!icon.IsEmpty());
    //    wxString xiconame = CreateIconame(filename);
    //    meta->AddProperty(wxT("icon"),wxFileSystem::FileNameToURL(xiconame));
    meta->AddAttribute("icon",wxFileSystem::FileNameToURL(icon));
  }
  //  x.Printf(wxT("%d"),(long long unsigned int)fitsmeta.GetSize());
  x << fitsmeta.GetSize();
  meta->AddAttribute("size",x);
  root->AddChild(meta);

  //  vector<FitsMetaHdu> hdu = fitsmeta.GetHdu();
  //  for(vector<FitsMetaHdu>::iterator h =hdu.begin();h !=hdu.end(); ++h){
    //for(size_t k = 0; k < hdu.size(); k++) {
  for(size_t k = 0; k < fitsmeta.HduCount(); k++) {
    //size_t k = h-hdu.begin();
    FitsMetaHdu h = fitsmeta.Hdu(k);

    wxXmlNode *xhdu = new wxXmlNode(0,wxXML_ELEMENT_NODE,"hdu");
    wxString num;
    num.Printf("%d",(int) k);
    xhdu->AddAttribute("index",num);
    xhdu->AddAttribute("type",/*hdu[k].*/h.Type_str());
    xhdu->AddAttribute("subtype",/*hdu[k].*/h.SubType_str());
    x.Printf("%d",/*hdu[k].*/(int) h.Ncols());
    xhdu->AddAttribute("ncols",x);
    x.Printf("%d",/*hdu[k].*/(int) h.Nrows());
    xhdu->AddAttribute("nrows",x);
    /*
    wxFileName xicon(filename);
    xicon.SetExt(num+wxT(".png"));
    wxString xiconame = xicon.GetFullPath();
    wxImage icon = hdu[k].GetIcon();
    icon.SaveFile(xiconame);
    */
    if( !filename.IsEmpty() ) {
      wxASSERT(!icons[k].IsEmpty());
      //      wxString xiconame = CreateIconame(filename,num);
      //      xhdu->AddProperty(wxT("icon"),wxFileSystem::FileNameToURL(xiconame));
      xhdu->AddAttribute("icon",wxFileSystem::FileNameToURL(icons[k]));
    }

    wxStringOutputStream buf;
    wxTextOutputStream xbuf(buf);
    vector<long> naxes = /*hdu[k].*/h.GetNaxes();
    for(size_t i = 0; i < naxes.size(); i++) {
      xbuf << int(naxes[i]);
      if( long(i + 1) < naxes[i] )
	xbuf << " ";

      /*
      for(size_t i = 0; i < hdu[k].Naxis(); i++) {
      xbuf << int(hdu[k].Naxes(i));
      if( i + 1 < hdu[k].Naxis() )
	xbuf << wxT(" ");
      */
    }

    wxXmlNode *xh = new wxXmlNode(0,wxXML_ELEMENT_NODE,"naxes");
    wxXmlNode *tn = new wxXmlNode(0,wxXML_TEXT_NODE,wxEmptyString,
				  buf.GetString());
    xh->AddChild(tn);
    xhdu->AddChild(xh);

    wxString hh;
    //    MuniThumbnailHdu hh(hdu[k]);
    for(size_t i = 0; i < /*hh.*/h.GetCount(); i++ )
      hh.Append(/*hh[i]*/h.Item(i)+"\n");

    // \n is the most ideal character to separate of strings. The
    // output xml looks as we are expecting and one is not allowed
    // in FITS header.

    wxXmlNode *fh = new wxXmlNode(0,wxXML_ELEMENT_NODE,"fitshead");
    wxXmlNode *t = new wxXmlNode(0,wxXML_CDATA_SECTION_NODE,wxEmptyString,hh);
    fh->AddChild(t);
    xhdu->AddChild(fh);

    meta->AddChild(xhdu);
  }

  return thumb;
}

void MuniThumbnail::Save(const wxString& filename)
{

//   wxXmlDocument thumb;

//   wxXmlNode *root = new wxXmlNode(0,wxXML_ELEMENT_NODE,wxT("munipack"));
//   root->AddProperty(wxT("version"),wxT(VERSION));
//   thumb.SetRoot(root);

//   wxXmlNode *meta = new wxXmlNode(0,wxXML_ELEMENT_NODE,wxT("fitsmeta"));
//   wxString x;
//   meta->AddProperty(wxT("url"),GetURL());
//   meta->AddProperty(wxT("type"),Type_str());
//   wxFileName xicon(filename);
//   xicon.SetExt(wxT("png"));
//   wxString xiconame = xicon.GetFullPath();
//   wxImage icon(GetIcon());
//   icon.SaveFile(xiconame);
//   meta->AddProperty(wxT("icon"),wxFileSystem::FileNameToURL(xiconame));
//   root->AddChild(meta);

//   for(size_t k = 0; k < hdu.size(); k++) {

//     wxXmlNode *xhdu = new wxXmlNode(0,wxXML_ELEMENT_NODE,wxT("hdu"));
//     wxString num;
//     num.Printf(wxT("%d"),k);
//     xhdu->AddProperty(wxT("index"),num);
//     xhdu->AddProperty(wxT("type"),hdu[k].Type_str());
//     xhdu->AddProperty(wxT("subtype"),hdu[k].SubType_str());
//     x.Printf(wxT("%d"),hdu[k].Ncols());
//     xhdu->AddProperty(wxT("ncols"),x);
//     x.Printf(wxT("%d"),hdu[k].Nrows());
//     xhdu->AddProperty(wxT("nrows"),x);
//     wxFileName xicon(filename);
//     xicon.SetExt(num+wxT(".png"));
//     wxString xiconame = xicon.GetFullPath();
//     wxImage icon = hdu[k].GetIcon();
//     icon.SaveFile(xiconame);
//     xhdu->AddProperty(wxT("icon"),wxFileSystem::FileNameToURL(xiconame));

//     wxStringOutputStream buf;
//     wxTextOutputStream xbuf(buf);
//     for(size_t i = 0; i < hdu[k].Naxis(); i++) {
//       xbuf << int(hdu[k].Naxes(i));
//       if( i + 1 < hdu[k].Naxis() )
// 	xbuf << wxT(" ");
//     }

//     wxXmlNode *xh = new wxXmlNode(0,wxXML_ELEMENT_NODE,wxT("naxes"));
//     wxXmlNode *tn = new wxXmlNode(0,wxXML_TEXT_NODE,wxEmptyString,
// 				  buf.GetString());
//     xh->AddChild(tn);
//     xhdu->AddChild(xh);

//     wxString h;
//     FitsHeader hh(hdu[k]);
//     for(size_t i = 0; i < hh.GetCount(); i++ )
//       h.Append(hh[i]+wxT("\n"));

//     // \n is the most ideal character to separate of strings. The
//     // output xml looks as we are expecting and one is not allowed
//     // in FITS header.

//     wxXmlNode *fh = new wxXmlNode(0,wxXML_ELEMENT_NODE,wxT("fitshead"));
//     wxXmlNode *t = new wxXmlNode(0,wxXML_CDATA_SECTION_NODE,wxEmptyString,h);
//     fh->AddChild(t);
//     xhdu->AddChild(fh);

//     meta->AddChild(xhdu);
//   }

  wxASSERT(meta.IsOk());

  icon.Clear();
  icon = CreateIconame(filename);

  icons.Clear();
  for(size_t k = 0; k < meta.HduCount(); k++) {
    wxString num;
    num.Printf("%d",(int) k);
    icons.Add(CreateIconame(filename,num));
  }

  wxXmlDocument xml = CreateXML(meta,filename);
  //  xml = CreateXML(filename);
  //  thumb.Save(filename);
  xml.Save(filename);

  //  wxImage icon(GetIcon());
  //  wxString xiconame = CreateIconame(filename);
  //  icon.SaveFile(xiconame);
  wxImage i = meta.GetIcon();
  i.SaveFile(icon);

  for(size_t k = 0; k < meta.HduCount(); k++) {
    FitsMetaHdu hdu = meta.Hdu(k);
    wxImage i = hdu.GetIcon();
    i.SaveFile(icons[k]);
  }

  /*
  for(vector<MuniThumbnailHdu>::const_iterator h=hdu.begin();h!=hdu.end();++h){
     wxString num;
     num.Printf(wxT("%d"),h-hdu.begin());
     wxString xiconame = CreateIconame(filename,num);
     wxImage i(h->GetImage());
     i.SaveFile(xiconame);
  }
  */

  /*
  for(size_t k = 0; k < hdu.size(); k++) {
     wxImage icon = hdu[k].GetIcon();
     icon.SaveFile(xiconame);
  }
  */
}


FitsMeta MuniThumbnail::GetMeta() const
{
  /*
  vector<FitsMetaHdu> hdus;

  for(vector<MuniThumbnailHdu>::const_iterator h=hdu.begin();h!=hdu.end();++h){
    FitsMetaHdu hh(*h,h->GetImage(),h->Ncols(),h->Nrows(),
		  h->Naxes(),h->Type_str(),h->SubType_str());
    hdus.push_back(hh);
  }

  return FitsMeta(url,type_str,hdus,icon);
  */
  return meta;
}


void MuniThumbnail::SetMeta(const FitsMeta& fm)
{
  meta = fm;
  /*
  url = fm.GetURL();
  type_str = fm.Type_str();
  icon = fm.GetIcon();
  vector<FitsMetaHdu> hdus(fm.GetHdu());

  hdu.clear();

  for(vector<FitsMetaHdu>::const_iterator h =hdus.begin(); h !=hdus.end(); ++h){
    MuniThumbnailHdu hh(*h,wxEmptyString,h->Ncols(),h->Nrows(),
		  h->GetNaxes(),h->Type_str(),h->SubType_str());
    hh.SetImage(h->GetIcon());
    hdu.push_back(hh);
  }
  */
}
