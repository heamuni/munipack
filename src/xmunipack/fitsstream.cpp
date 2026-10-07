/*

  xmunipack - FITS stream for input

  Copyright © 2019-22 F.Hroch (hroch@physics.muni.cz)

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

#include "xmunipack.h"
#include "fits.h"
#include "fitsdisplay.h"
#include <fitsio.h>
#include <wx/wx.h>
#include <vector>
#include <algorithm>

using namespace std;

// --- FitsStream -----------------
// A serial loader of FITS files -- a file load process is breaked on chunks,
// and an event is issued when the chunk is successfully loaded.

/*
  FITSes can be large files. Ones can be downloaded over the network;
  ones should be parsed, and image values should be computed.

  The stream loader separate of the loading procedure:

  * Headers are red, and passed, as whole, and the calling procedure
    is notified.

  * Images are read by blocks, every block is reported and passed up.

  * Tables are read as blocks, current row order is reported. Table data
    are unavailable during load.

Consider:

  * The preview handling leads to a complex code (also in View);
    I can consider to replace it by a plain progress with percents
    only (no frame display).

  * Utilising of stream driver of cfitsio can be important for
    network downloads.

*/


FitsStream::FitsStream(wxEvtHandler *h, const wxString& fn, int is):
  wxThread(wxTHREAD_JOINABLE),handler(h),filename(fn),
  shrink(0),icon_size(is),stop(false),stop_wait(false)
{
  wxASSERT(handler);
}

void FitsStream::Stop()
{
  stop = true;
}

void FitsStream::StopAndWait()
{
  stop = true;
  stop_wait = true;
}


wxThread::ExitCode FitsStream::Entry()
{
  fitsfile *fits;
  int status, htype;
  vector<FitsHdu> hdu;

  // open file
  status = 0;
  fits_open_file(&fits,filename.fn_str(), READONLY, &status);

  if( status == 0 ) {
    FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_OPEN);
    ev.filename = filename;
    wxQueueEvent(handler,ev.Clone());
  }

  //wxMilliSleep(20000);

  // iterate over HDUs
  int nhdu = 0;
  int chdu;
  while ( status == 0 && stop == false ) {

    // We're terminating, if there're no more HDUs to read
    int mstatus = 0;
    fits_movabs_hdu(fits,nhdu+1,&htype,&mstatus);
    if( mstatus ) break;
    fits_get_hdu_num(fits,&chdu);
    nhdu = nhdu + 1;

    // load header
    int nhead;
    char h[FLEN_CARD];
    FitsHeader head;
    fits_get_hdrspace(fits,&nhead,NULL,&status);
    for(int n = 0; n < nhead && status == 0; n++) {
      if( fits_read_record(fits,n+1,h,&status) == 0 )
	head.Add(wxString::FromAscii(h));
    }
    if( status || stop ) break;

    // determine HDU type
    hdu_type hdutype = head.GetType(htype);
    hdu_flavour flavour = head.GetFlavour(htype,hdutype);

    FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_HEAD);
    ev.head = head;
    ev.extname = head.GetKey("EXTNAME");
    ev.hdutype = hdutype;
    ev.flavour = flavour;
    ev.chdu = chdu;
    wxQueueEvent(handler,ev.Clone());

    // decode current HDU
    if( htype == IMAGE_HDU ) {

      int bitpix, naxis;

      fits_get_img_type(fits,&bitpix,&status);
      fits_get_img_dim(fits,&naxis,&status);
      if( status ) break;

      // header only
      if( hdutype == HDU_HEAD ) {

	hdu.push_back(FitsHdu(head,htype));

	FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_HEADONLY);
	ev.chdu = chdu;
	wxQueueEvent(handler,ev.Clone());

      }
      // full image
      else {

	float *nullval = 0;
	int anynul;
	long *naxes = new long[naxis];
	fits_get_img_size(fits,naxis,naxes,&status);
	if( status ) { delete[] naxes; break; }

	FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_IMGHEAD);
	ev.chdu = chdu;
	ev.bitpix = bitpix;
	ev.naxis = naxis;
	vector<long> maxes(naxes,naxes+naxis);
	ev.naxes = maxes;
	ev.progress = 0;
	wxQueueEvent(handler,ev.Clone());

       	long ndata = 1;
	for(int i = 0; i < naxis; i++ ) ndata = ndata*naxes[i];
	float *image = new float[ndata];
	if( image == 0 ) {
	  status = -1;
	  wxLogError("Failed to allocate memory for an image.");
	}

	if( hdutype == HDU_IMAGE &&
	    (flavour == HDU_IMAGE_FRAME || flavour == HDU_IMAGE_COLOUR) ) {
	  // regular images
	  wxASSERT((naxis == 2 || naxis == 3) && naxes[0] > 0);
	  long coffset = naxes[0] * naxes[1];
	  int maxk = naxis == 3 ? 3 : 1;
	  for(int k = 0; k < maxk; k++) {
	    long maxrows = wxMax((NIOBUF*IOBUFLEN) / (sizeof(float)*naxes[0]),1L);
	    wxASSERT(maxrows > 0);
	    long jrow = 0;
	    for(long j = 0; j < naxes[1]; j += maxrows) {
	      long nrows = std::min(maxrows,naxes[1]-j);
	      long nelements = nrows*naxes[0];
	      wxASSERT(nelements > 0);

	      long offset = naxes[0]*j + k*coffset;
	      long fpixel = offset + 1;
	      float *ptr = image + offset;

	      fits_read_img(fits,TFLOAT,fpixel,nelements,nullval,
			    ptr,&anynul,&status);
	      if( status ) goto crash;

	      FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_IMAGE);
	      ev.chdu = chdu;
	      ev.cube = k+1;
	      ev.crow = j+1;
	      ev.progress = double(j+nrows)/double(maxes[0]);
	      wxQueueEvent(handler,ev.Clone());

	      if( chdu == 1 ) {

		if( j == 0 && tone.IsOk() == false ) {
		  tone = FitsTone(nelements,ptr);
		  tone.SetQblack(0.0);
		  tone.SetRsense(0.05);
		}

		if( shrink > 0 && j > jrow + shrink ) {
		  FitsStreamEvent pev(EVT_FITS_STREAM,ID_FITSLOAD_IMGPREVIEW);
		  pev.chdu = chdu;
		  pev.cube = k+1;
		  pev.crow = j+1;
		  long off = naxes[0]*jrow + k*coffset;
		  pev.preview = RenderPreview(naxes[0],j-jrow+1,image+off);
		  pev.naxes = maxes;
		  pev.progress = double(j+nrows)/double(maxes[0]);
		  wxQueueEvent(handler,pev.Clone());
		  jrow = j;
		}
	      }

	      if( stop ) goto crash;

	      // wxMilliSleep(250);
	    }
	  }
	}
	else {
	  // line, cube, no preview, needs testing

	  long bufsize = NIOBUF*IOBUFLEN;
	  long blocksize= bufsize / sizeof(float);
	  long nblocks = ndata / blocksize + 1;
	  //	  wxLogDebug("%d %d %dx%d",int(bufsize),int(nblocks),int(naxes[0]),int(naxes[1]));
	  for(long j = 0; j < nblocks; j++) {
	    long nelements = std::min(blocksize,ndata-j*blocksize);
	    long offset = j*blocksize;
	    long fpixel = offset + 1;
	    float *ptr = image + offset;
	    //	    wxLogDebug("%d %d",int(j),int(nelements));

	    fits_read_img(fits,TFLOAT,fpixel,nelements,nullval,
			  ptr,&anynul,&status);

	    if( status ) goto crash;
	    if( stop ) goto crash;
	    //     wxMilliSleep(25);
	  }
	}
      crash:

	if( status == 0 && stop == false  ) {
	  hdu.push_back(FitsArray(head,htype,naxis,naxes,image));
	}
	else {
	  delete[] naxes;
	  delete[] image;
	  break;
	}
      }
    }
    else if( htype == ASCII_TBL || htype == BINARY_TBL ) {

      long nrows, ncols, srows;
      int nc;
      fits_get_num_rows(fits,&nrows,&status);
      fits_get_num_cols(fits,&nc,&status);
      if( status ) break;
      ncols = nc;
      vector<FitsTableColumn> cols;

      // initialise full table
      for(int k = 0; k < ncols && status == 0; k++ ) {

	int colnum = k + 1;
	int typecode;
	long repeat, width;

	fits_get_coltype(fits,colnum,&typecode,&repeat,&width,&status);
	if( status == 0 )
	  cols.push_back(FitsTableColumn(typecode,repeat,width,nrows));
      }
      if( status ) break;

      FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_TBLHEAD);
      ev.nrows = nrows;
      ev.ncols = ncols;
      ev.chdu = chdu;
      wxQueueEvent(handler,ev.Clone());

      fits_get_rowsize(fits,&srows,&status);
      for(long row = 0; row < nrows && stop == false; row = row + srows) {
	long felem = 1;
	long frow = row + 1;
	long nelem = min(srows,nrows-row);

	for(int col = 0; col < ncols; col++ ) {

	  int colnum = col + 1;
	  int typecode, anynul;
	  long repeat, width;

	  fits_get_coltype(fits,colnum,&typecode,&repeat,&width,&status);

	  if( typecode == TSTRING ) {
	    char **a = cols[col].GetData_string();
	    char *nullval = 0;
	    fits_read_col(fits,TSTRING,colnum, frow, felem, nelem, &nullval,
			  a+row,&anynul,&status);
	  }
	  else if( typecode == TLOGICAL ) {
	    bool *b = cols[col].GetData_bool();
	    bool nullval = false;
	    fits_read_col(fits,TLOGICAL,colnum, frow, felem, nelem, &nullval,
			b+row,&anynul,&status);
	  }
	  else if( typecode == TBYTE || typecode == TBIT ) {
	    char *b = cols[col].GetData_char();
	    char nullval = ' ';
	    fits_read_col(fits,TBYTE,colnum, frow, felem, nelem, &nullval,
			  b+row,&anynul,&status);
	  }
	  else if( typecode == TSHORT ) {
	    short *d = cols[col].GetData_short();
	    short nullval = 0;
	    fits_read_col(fits,TSHORT,colnum, frow, felem, nelem, &nullval,
			   d+row,&anynul,&status);
	    }
	  else if( typecode == TLONG || typecode == TINT32BIT ) {
	    long *d = cols[col].GetData_long();
	    long nullval = 0;
	    fits_read_col(fits,TLONG,colnum, frow, felem, nelem, &nullval,
			  d+row,&anynul,&status);
	  }
	  else if( typecode == TFLOAT ) {
	    float *d = cols[col].GetData_float();
	    float nullval = 0.0;
	    fits_read_col(fits,TFLOAT,colnum, frow, felem, nelem, &nullval,
			  d+row,&anynul,&status);
	  }
	  else if( typecode == TDOUBLE ) {
	    double *d = cols[col].GetData_double();
	    double nullval = 0.0;
	    fits_read_col(fits,TDOUBLE,colnum, frow, felem, nelem, &nullval,
			  d+row,&anynul,&status);
	  }
	  else
	    wxLogWarning("FitsStream: The type code `%d' not implemented yet.",
			 typecode);

	  if( status || stop ) goto clash;

	  /*
             Tables are read as blocks of rows according to cFISTIO recommendations.
	     We are notifying calling thread by current row. No data are passed --
             tables shouldn't been the first HDU. The only first HDU will be
             displayed during FITS opening.
	  */

	}

	FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_TABLE);
	ev.chdu = chdu;
	ev.crow = row + 1;
	ev.nrows = nrows;
	ev.ncols = ncols;
	ev.progress = double(row+nelem)/double(nrows);
	wxQueueEvent(handler,ev.Clone());

	//	wxMilliSleep(1000);
      }
    clash:

      hdu.push_back(FitsTable(head,htype,cols));
    }
    else {
      wxLogWarning("FitsStream: Unrecognized htype code `%d'.",htype);
    }

    FitsStreamEvent fv(EVT_FITS_STREAM,ID_FITSLOAD_HDUFIN);
    ev.chdu = chdu;
    wxQueueEvent(handler,fv.Clone());

  }

  // close file
  fits_close_file(fits, &status);

  //  wxLogDebug("********** %d %d %d",status,stop,stop_wait);

  if( status == 0 ) {

    if( stop_wait )
      ;

    else if( stop ) {

      FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_INTERRUPT);
      ev.filename = filename;
      ev.errmsg = "An interrupt on request.";
      wxQueueEvent(handler,ev.Clone());

    } else {

      FitsFile fitsfile(filename,hdu);

      vector<FitsTone> tones;
      for(size_t k = 0; k < hdu.size(); k++)
	if( hdu.at(k).IsDisplayImplemented() ) {
	  FitsTone tone(FitsArray(hdu.at(k)));
	  tones.push_back(tone);
	}
	else {
	  FitsTone tone;
	  tones.push_back(tone);
	}

      vector<wxImage> icons;
      make_icons(fitsfile,tones,icons);

      FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_FINISH);
      ev.filename = filename;
      ev.fitsfile = fitsfile;
      ev.icons = icons;
      ev.tones = tones;
      wxQueueEvent(handler,ev.Clone());
    }

  }
  else {

    // save error description
    char emsg[FLEN_ERRMSG];
    char msg[FLEN_STATUS];
    wxArrayString errlog;

    while( fits_read_errmsg(emsg) )
      errlog.Add(wxString::FromAscii(emsg));
    fits_get_errstatus(status,msg);

    FitsStreamEvent ev(EVT_FITS_STREAM,ID_FITSLOAD_CRASH);
    ev.filename = filename;
    ev.errmsg = wxString::FromAscii(msg);
    ev.errlog = errlog;
    wxQueueEvent(handler,ev.Clone());

  }

  return (wxThread::ExitCode) 0;
}

wxImage FitsStream::RenderPreview(int width, int height, const float *ptr)
{
  wxASSERT(shrink > 0);

  double zoom = 1.0 / shrink;

  int canvas_width = int(width * zoom + 1);
  int canvas_height = int(height * zoom + 1);
  wxASSERT(canvas_width > 0 && canvas_height > 0);

  int npix = canvas_width*canvas_height;
  int nbytes = 3*npix*sizeof(unsigned char);
  unsigned char *d = (unsigned char *) malloc(nbytes);
  memset(d,0,nbytes);

  for(int j = 0; j < canvas_height; j++) {
    int jj = shrink*j;
    unsigned char *line = d + 3*j*canvas_width;
    const float *pline = ptr + wxMax(height-1-jj,0)*width;
    for(long i = 0; i < canvas_width; i++) {
      int ii = wxMin(shrink*i,width);
      //wxLogDebug("%d %d %d %d",int(i),int(j),int(ii),int(jj));
      float f = tone.Scale(pline[ii]);
      unsigned char rgb = std::min(std::max(int(255*f),0),255);
      memset(line+3*i,rgb,3);
    }
  }
  return wxImage(canvas_width,canvas_height,d);
}


void FitsStream::SetZoom(float z)
{
  if( shrink > 0 ) return;

  if( z < 1.01 )
    shrink = int(1.0 / z + 0.5);
  else
    shrink = 1;
}


void FitsStream::make_icons(const FitsFile& fits, const vector<FitsTone>& tones,
			    vector<wxImage>& icons)
{
  wxLogDebug("FitsStream:: making icons");
  wxASSERT(fits.IsOk());

  // for images, estimate stat. parameters of all HDUs, and make icons
  for(size_t k = 0; k < fits.HduCount(); k++) {
    if( fits.Hdu(k).IsDisplayImplemented() ) {
      FitsArray array(fits.Hdu(k));
      FitsDisplay display(array.IsColour());
      display.SetTone(tones[k]);
      wxImage icon = display.MakeIcon(array,icon_size,icon_size);
      icons.push_back(icon);
    }
    else {
      icons.push_back(wxImage());
    }
  }
}
