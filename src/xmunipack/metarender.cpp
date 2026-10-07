/*

  xmunipack - list FitsMeta render

  Copyright © 2022 F.Hroch (hroch@physics.muni.cz)

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

#include "list.h"
#include "fitsdisplay.h"
#include <wx/wx.h>
#include <wx/thread.h>
#include <wx/filesys.h>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

// ----  MetaRender

MetaRender::MetaRender(wxEvtHandler *eh, const size_t i, const wxArrayString& f,
		       const wxSize& s):
  wxThread(wxTHREAD_DETACHED),handler(eh),id(i),files(f),size(s)
{
  wxASSERT(handler);
  wxLogDebug("MetaRender::MetaRender() %d",int(files.GetCount()));
}

void *MetaRender::Entry()
{
  bool imake = size.GetWidth() > 0;

  for(size_t i = 0; i < files.GetCount(); i++) {
    if( TestDestroy() ) break;

    array = FitsArray();
    hdu.clear();
    wxImage icon;
    FitsMeta meta;
    wxFileName fn(wxFileSystem::URLToFileName(files[i]));
    wxULongLong filesize = fn.GetSize();
    wxString filename = fn.GetFullPath();

    if( MagicFile(filename) ) {
      if( FitsOpen(filename,imake) ) {
	if( imake ) {
	  if( array.IsOk() ) {
	    if( array.IsDisplayImplemented() ) {
	      FitsTone tone(array);
	      FitsDisplay display(array.IsColour());
	      display.SetTone(tone);
	      icon = display.MakeIcon(array,size.GetWidth(),size.GetHeight());
	    }
	    else
	      wxLogDebug("MetaRender::Entry(): `"+files[i]+
			 "' is corrupted or have unimplemeted shape of the array.");
	  }
	}
	meta = FitsMeta(files[i],hdu,filesize);
      }
      else
	wxLogDebug("MetaRender::Entry(): `"+files[i]+"' failed to open as FITS.");
    }
    else
      wxLogDebug("MetaRender::Entry() on "+files[i]+"'; the magic say: no FITS.");


    MetaOpenEvent ev(EVT_META_OPEN,ID_METALOAD_PHASE);
    ev.id = id;
    ev.index = i;
    ev.meta = meta;
    ev.icon = icon;
    wxQueueEvent(handler,ev.Clone());
  }

  MetaOpenEvent ev(EVT_META_OPEN,ID_METALOAD_FINISH);
  ev.id = id;
  wxQueueEvent(handler,ev.Clone());

  return (wxThread::ExitCode) 0;
}


bool MetaRender::FitsOpen(const wxString& filename, bool pick_array)
{
  fitsfile *f;
  int status = 0;
  int dummy, htype, bitpix, naxis;

  hdu.clear();
  array = FitsArray();
  int nhdu = 0;

  // open file
  status = 0;
  fits_open_file(&f, filename.fn_str(), READONLY, &status);
  if( status ) goto crash;

  fits_get_num_hdus(f,&nhdu,&status);
  if( status ) goto crash;

  for(int k = 0; k < nhdu; k++) {

    fits_movabs_hdu(f,k+1,&htype,&status);
    if( status ) goto crash;

    // load header
    int nhead;
    char h[FLEN_CARD];
    FitsHeader head;
    fits_get_hdrspace(f,&nhead,&dummy,&status);
    for(int n = 0; status == 0 && n < nhead; n++) {
      if( fits_read_record(f,n+1,h,&status) == 0 )
	head.Add(wxString(h,wxConvUTF8));
    }
    if( status ) goto crash;

    if( htype == IMAGE_HDU ) {

      fits_get_img_type(f,&bitpix,&status);
      fits_get_img_dim(f,&naxis,&status);
      if( status ) goto crash;

      int hdutype = naxis > 0 ? HDU_IMAGE : HDU_HEAD;
      if( hdutype == HDU_HEAD )
	hdu.push_back(FitsMetaHdu(head,htype));
      else {
	long *naxes = new long[naxis];
	fits_get_img_size(f,naxis,naxes,&status);
	if( status ) { delete[] naxes; goto crash; }
	hdu.push_back(FitsMetaHdu(head,htype,naxis,naxes));

	if( pick_array && !array.IsOk() ) {
	  // this image becomes icon

	  long *naxes = new long[naxis];
	  fits_get_img_size(f,naxis,naxes,&status);
	  if( status ) { delete[] naxes; goto crash; }

	  long ndata = 1; for(int i = 0; i < naxis; i++ ) ndata = ndata*naxes[i];
	  long firstelem = 1;
	  float nullval = 0.0;

	  float *image = new float[ndata];
	  wxASSERT(image);
	  fits_read_img(f,TFLOAT,firstelem,ndata,&nullval,image,&dummy,&status);
	  if( status ) { delete[] naxes; delete[] image; goto crash; }

	  array = FitsArray(head,htype,naxis,naxes,image);
	}
      }
    }
    else if( htype == ASCII_TBL || htype == BINARY_TBL ) {
      long nrows;
      int ncols;
      fits_get_num_rows(f,&nrows,&status);
      fits_get_num_cols(f,&ncols,&status);
      hdu.push_back(FitsMetaHdu(head,htype,nrows,ncols));
    }
    else
      wxLogError("FitsStream: `"+filename+"' unrecognized htype for HDU #%d",k);
  } // for over HDUs

 crash:

  fits_close_file(f, &status);
  if( status ) {
    char emsg[FLEN_ERRMSG];
    if( status ) wxLogDebug("FitsMeta():");
    while( fits_read_errmsg(emsg) )
      wxLogDebug(wxString(emsg,wxConvUTF8));
  }
  return status == 0;
}

bool MetaRender::MagicFile(const wxString& filename) const
{
  const char *record = "SIMPLE  =                    T";
  const size_t size = strlen(record);
  bool result = false;

  FILE *f = fopen(filename.fn_str(),"r");
  if( f != NULL ) {
    char buf[size+1];
    if( fread(buf,1,size,f) == size ) {
      buf[size] = '\0';
      result = strcmp(buf,record) == 0;
    }
    fclose(f);
  }
  return result;
}
