/*

  Growth-Curve Photometry Report

  Copyright © 2016 F.Hroch (hroch@physics.muni.cz)

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

#include <fitsio.h>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>


//#define DEBUG 1

#define C(a) const_cast<char*>((a))

using namespace std;


float farray2d(float *pic, int w, int h, int x, int y, int d)
{
  int xdim,ydim,ndim;
  int x1 = x, y1 = y, x2 = x + d, y2 = y + d;
  int i, j, l, n;
  double s, d2;

  assert(pic != NULL);
  assert(d >= 1);

  xdim = w;
  ydim = h;

  if( d == 1 ) {

    assert((0 <= x && x < xdim) && (0 <= y && y < ydim));
    return(pic[y*xdim + x]);

  }
  else if( d > 1 ) {

    d2 = d*d;

    if( (0 <= x1 && x1 < xdim) && (0 <= y1 && y1 < ydim) &&
        (0 <= x2 && x2 < xdim) && (0 <= y2 && y2 < ydim)) {

      s = 0.0;
      for(j = y1; j < y2; j++) {
        float *nrow = pic + j*xdim;
        for(i = x1; i < x2; i++)
          s = s + *(nrow + i);
      }
      return(s/d2);
    }
    else {
      s = 0.0;
      ndim = w*h;
      n = 0;
      for(j = y1; j < y2; j++) {
        for(i = x1; i < x2; i++) {
	  l = i+j*xdim;
	  if( 0 <= l && l < ndim ) {
	    s = s + pic[l];
	    n++;
	  }
	}
      }
      if( n > 0 )
	return(s/n);
      else
	return(s);
    }
  }
  return(0.0);
}


float qmed(int n, float *a, int k)
{
  float w,x;
  int l,r,i,j;

  x = 0.0;

  l = 0;
  r = n - 1;

  while( l < r ) {
    x = a[k];
    i = l;
    j = r;

    do {

      while( a[i] < x )
	i++;

      while( x < a[j] )
	j--;

      if( i <= j ) {
	w = a[i];
	a[i] = a[j];
	a[j] = w;
	i++;
	j--;
      }

    } while ( i <= j );

    if( j < k ) l = i;
    if( k < i ) r = j;
  }

  return(x);
}


float extinvdist(float mean, float sig)
{
  /*
     returns random variable from exponential distribution
     https://en.wikipedia.org/wiki/Laplace_distribution
  */

  float x,sign;
  const float rmax = RAND_MAX;

  x = rand() / rmax;
  x = x - 0.5;
  sign = x >= 0 ? 1 : -1;
  x = mean - sig*sign*log10(1-2*fabs(x));
  return(x);
}


int fitsflush(string fitsname, int fd, int gd, int gs, int dd, int gc, int pd,
	      int *nx, int *ny, int *rmax, bool *psf)
{
  const long firstelem = 1;
  float nullval = 0.0;
  char nullvac = 0;
  const long frow = 1, felem = 1;

  FILE *fc,*ff,*fgd, *fdd, *fpsf;
  fitsfile *f;
  int dummy, bitpix, naxis,i,j,n,ii,jj,zoom, nstep, col, naper, stat, nc,
    col_x, col_y, col_f;
  long nrows,ndata;
  long *naxes;
  float *data,*xdata, *grow, *res, *raper, *xcoo, *ycoo, x,mad,med,s,t,anu1,rf90;
  char *flag;
  unsigned char obr[4],g;
  char com[80],key[80];

  stat = 0;
  fits_open_file(&f,fitsname.c_str(), READONLY, &stat);

  assert(stat == 0);

  // image
  fits_get_img_type(f,&bitpix,&stat);
  fits_get_img_dim(f,&naxis,&stat);
  naxes = new long[naxis];
  fits_get_img_size(f,naxis,naxes,&stat);
  ndata = 1; for(i = 0; i < naxis; i++ ) ndata = ndata*naxes[i];
  data = new float[ndata];
  fits_read_img(f,TFLOAT,firstelem,ndata,&nullval,data,&dummy,&stat);

  // estimate zoom, final height is limited by 800
  zoom = 1;
  n = naxes[1] / 800;
  if( n > 0 )
    zoom = n + 1;

  // size of plotting area
  *nx = naxes[0] / zoom;
  *ny = naxes[1] / zoom;

  // get size of annulus
  fits_movnam_hdu(f,BINARY_TBL,C("APERPHOT"),0,&stat);
  fits_read_key(f,TFLOAT,"ANNULUS1",&anu1,com,&stat);

  // grow curve
  fits_movnam_hdu(f,BINARY_TBL,C("GROWFUNC"),0,&stat);
  fits_get_num_rows(f,&nrows,&stat);
  raper = new float[nrows];
  grow = new float[nrows];
  fits_get_colnum(f,CASEINSEN,C("R"),&col,&stat);
  fits_read_col(f,TFLOAT,col, frow, felem, nrows, &nullval, raper,&dummy,&stat);
  fits_get_colnum(f,CASEINSEN,C("GROWCURVE"),&col,&stat);
  fits_read_col(f,TFLOAT,col, frow, felem, nrows, &nullval, grow,&dummy,&stat);
  fc = fdopen(gc,"w");
  for(int i = 0; i < nrows; i++)
    fprintf(fc,"%.3f %.3f\n",raper[i],grow[i]);
  fclose(fc);
  *rmax = round(raper[nrows-1]) + 1;
  delete[] raper;
  delete[] grow;

  // stars used for grow-curve, radius of flux of 90%
  fits_movnam_hdu(f,BINARY_TBL,C("GROWPHOT"),0,&stat);
  fits_read_key(f,TFLOAT,"RADFLX90",&rf90,com,&stat);
  fits_get_num_rows(f,&nrows,&stat);
  fits_get_colnum(f,CASEINSEN,C("X"),&col_x,&stat);
  fits_get_colnum(f,CASEINSEN,C("Y"),&col_y,&stat);
  fits_get_colnum(f,CASEINSEN,C("GROWFLAG"),&col_f,&stat);
  xcoo = new float[nrows];
  ycoo = new float[nrows];
  flag = new char[nrows];
  fits_read_col(f,TFLOAT,col_x, frow, felem, nrows, &nullval, xcoo,&dummy,&stat);
  fits_read_col(f,TFLOAT,col_y, frow, felem, nrows, &nullval, ycoo,&dummy,&stat);
  fits_read_col(f,TBYTE,col_f, frow, felem, nrows, &nullvac, flag,&dummy,&stat);

  ff = fdopen(gs,"w");
  for(int i = 0; i < nrows; i++) {
    if( flag[i] == 1 )
      fprintf(ff,"%.1f %.1f %.2f\n",xcoo[i]/zoom,ycoo[i]/zoom,anu1/zoom);
  }
  fclose(ff);

  // empirical grow-curves
  fits_movnam_hdu(f,BINARY_TBL,C("GROWDATA"),0,&stat);

  fits_read_key(f,TINT,"NAPER",&naper,com,&stat);
  raper = new float[naper];
  for(int i = 1; i <= naper; i++) {
    sprintf(key,"%s%d","APER",i);
    fits_read_key(f,TFLOAT,key,&raper[i-1],com,&stat);
  }

  fgd = fdopen(gd,"w");
  fdd = fdopen(dd,"w");

  fits_get_num_rows(f,&nrows,&stat);
  grow = new float[nrows];
  res = new float[nrows];
  for(int i = 0; i < naper; i++) {
    n = i + 1;
    sprintf(key,"%s%d","GROWCURVE",n);
    fits_get_colnum(f,CASEINSEN,key,&col,&stat);
    fits_read_col(f,TFLOAT,col, frow, felem, nrows, &nullval, grow,&dummy,&stat);

    sprintf(key,"%s%d","RESGROW",n);
    fits_get_colnum(f,CASEINSEN,key,&col,&stat);
    fits_read_col(f,TFLOAT,col, frow, felem, nrows, &nullval, res,&dummy,&stat);

    for(int j = 0; j < nrows; j++) {
      float r = raper[i] + extinvdist(0.0,0.03*raper[i]);
      /*
      double dx = xcoo[j] - int(xcoo[j]);
      double dy = ycoo[j] - int(ycoo[j]);
      double r = raper[i] + sqrt(dx*dx + dy*dy)*0;
      fprintf(stderr,"r: %f %f %f %f %f %f\n",xcoo[j],ycoo[j],dx,dy,r,sqrt(dx*dx + dy*dy));
      */
      if( flag[j] == 1 )
	fprintf(fgd,"%5.2f %7.4f %7.4f\n",r,grow[j],res[j]);
      else
	fprintf(fdd,"%5.2f %7.4f %7.4f\n",r,grow[j],res[j]);
    }
  }

  // PSF (artificial data only)
  fpsf = fdopen(pd,"w");

  fits_movnam_hdu(f,IMAGE_HDU,C("PSF"),0,&stat);
  if( stat == 0 ) {
    *psf = true;

    // PSF is presented

    int bitpix, naxis, zoom;

    fits_read_key(f,TINT,"ZOOM",&zoom,com,&stat);

    fits_get_img_type(f,&bitpix,&stat);
    fits_get_img_dim(f,&naxis,&stat);
    long *naxes = new long[naxis];
    fits_get_img_size(f,naxis,naxes,&stat);
    long ndata = 1; for(i = 0; i < naxis; i++ ) ndata = ndata*naxes[i];
    float *psf = new float[ndata];
    fits_read_img(f,TFLOAT,firstelem,ndata,&nullval,psf,&dummy,&stat);

    long *cens = new long[naxis];
    for(int i = 0; i < naxis; i++) cens[i] = naxes[i] / 2 + 1;

    fprintf(fpsf,"0 0\n");

    for(int n = 1; n < zoom* *rmax; n++) {
      float f = 0;
      for(int i = -n; i <= n; i++) {
	int l = i + cens[0];
	for(int j = -n; j <= n; j++) {
	  int k = j + cens[1];
	  if( (i*i + j*j <= n*n) && (0 < l && l < naxes[0]) && (0 < k && k < naxes[1]) )
	    f = f + psf[l+k*naxes[0]];
	}
      }
      fprintf(fpsf,"%.1f %.6f\n",float(n)/zoom,f/(zoom*zoom));
    }

    delete[] naxes;
    delete[] cens;
    delete[] psf;
  }
  else {
    stat = 0;
    fprintf(fpsf,"# This file is intentionally left empty.\n");
  }
  fclose(fpsf);

  delete[] grow;
  delete[] res;
  delete[] flag;
  delete[] raper;
  delete[] xcoo;
  delete[] ycoo;
  fclose(fgd);
  fclose(fdd);

  fits_close_file(f, &stat);

  // estimate levels
  nstep = ndata / 10000;
  nc = ndata / nstep;
  xdata = new float[nc];

  n = 0;
  for(int i = 0; n < nc && i < ndata-nstep-1; i = i + nstep)
    xdata[n++] = data[i];
  med = qmed(n,xdata,n/2+1);

  n = 0;
  for(int i = 0; n < nc && i < ndata-nstep-1; i = i + nstep) {
    x = data[i] - med;
    if( x > 0 )
      xdata[n++] = x;
  }
  mad = qmed(n,xdata,n/2+1);
  delete[] xdata;

  // scaling parameters
  t = med + 3.0 * mad;
  s = 10*mad;
  if( mad < 0.2 ) s = 1;

#ifdef DEBUG
  fprintf(stderr,"%ld %ld %d %d %.1f %.1f\n",naxes[0],naxes[1],stat,zoom,mad,med);
#endif

  // flushing image to AVS frame ...
  for(i = 0; i < 2; i++) {
    n = naxes[i] / zoom;
    ssize_t bytes = write(fd,&n,4);
    if( bytes != 4 ) {
#ifdef DEBUG
      fprintf(stderr,"%s: Flushing to AVS failed.\n",__FILE__);
#endif
    }
  }

  /*
    Gnuplot documenation says "2 longs (xwidth, ywidth) followed by a stream of pixels,
    each with four bytes of information alpha/red/green/blue." which does not works.
    Lenghts of xwidth and ywidth are 'int' (4 bytes). By my opinion, the inconsistency
    has origin in 16/32/64-bit coding of integers on various platforms.

    `help avs` or  http://gnuplot.sourceforge.net/docs_4.2/node109.html
  */

  for(j = naxes[1]/zoom - 1; j >=0; j--) {
    jj = j * zoom;

    for(i=0; i < naxes[0]/zoom; i++) {
      ii = i * zoom;

      // zoom
      x = farray2d(data,naxes[0],naxes[1],ii,jj,zoom);

      // intensity scaling
      x = (x - t) / s;
      x = (1 + erf(2*x - 1)) / 2.0;

      // gamma function by sRGB
      if( x <= 3.1308e-3 )
	x = 12.92*x;
      else if( x <= 1 )
	x = 1.055*pow(x,1.0/2.4) - 0.055;

      // convert to interval 0 - 255,
      // because results of erf are limited, cut-off is not useless
      g = 255*x;

      // set up RGB, invert intensity
      obr[0] = 255; // alpha channel
      for(n = 1; n <= 3; n++)
	obr[n] = 255 - g;
      ssize_t bytes = write(fd,obr,4);
      if( bytes != 4 ) {
#ifdef DEBUG
	fprintf(stderr,"%s: Flushing to gnuplot failed.\n",__FILE__);
#endif
    }
    }
  }

  delete[] naxes;
  delete[] data;

  return 0;
}

void report(string fitsname)
{
  string pngname = fitsname;
  string suf = ".fits";
  size_t n = fitsname.rfind(suf);
  if( n != string::npos )
    pngname.replace(n,suf.length(),".png");
  else
    pngname = fitsname + ".png";

#ifdef DEBUG
  char avsname[10] = "/tmp/x";
  char dname[10] = "/tmp/d";
  char ddname[10] = "/tmp/b";
  char gsname[13] = "/tmp/g";
  char cname[13] = "/tmp/f";
  char growname[10] = "/tmp/grow";
  char gplname[10] = "/tmp/gpl";
  FILE *df = fopen(avsname,"w");
  FILE *dp = fopen(dname,"w");
  FILE *ddp = fopen(ddname,"w");
  FILE *dg = fopen(gsname,"w");
  FILE *dc = fopen(cname,"w");
  FILE *ds = fopen(growname,"w");
  int fd = fileno(df);
  int gd = fileno(dp);
  int dd = fileno(ddp);
  int gs = fileno(dg);
  int gc = fileno(dc);
  int gf = fileno(ds);
#else
  char avsname[10] = "avsXXXXXX";
  char dname[11] = "dataXXXXXX";
  char ddname[10] = "datXXXXXX";
  char gsname[13] = "gstarsXXXXXX";
  char cname[13] = "funXXXXXX";
  char growname[11] = "growXXXXXX";
  char gplname[10] = "gplXXXXXX";
  int fd = mkstemp(avsname);
  int gd = mkstemp(dname);
  int dd = mkstemp(ddname);
  int gs = mkstemp(gsname);
  int gc = mkstemp(cname);
  int gf = mkstemp(growname);
#endif


  int nx,ny,rmax;
  bool psf = false;

  fitsflush(fitsname,fd,gd,gs,dd,gc,gf,&nx,&ny,&rmax,&psf);

  ofstream gpl(gplname);
  gpl << "set output '" << pngname << "'" << endl;
  gpl << "set term 'png' truecolor size " << 2*nx << "," << ny  << endl;
  gpl << "set multiplot title 'Growth-curve (" << fitsname << ")'" << endl;
  gpl << "set style line 99 lt 2 lc rgb 'gold' lw 2" << endl;
  gpl << "set style line 98 lt 2 lc rgb 'gray' lw 4" << endl;
  gpl << "unset key" << endl;
  gpl << "set lmargin 5" << endl;
  gpl << "set rmargin 2" << endl;
  gpl << "set origin 0.5,0.3" << endl;
  gpl << "set size 0.5,0.7" << endl;
  gpl << "set key right bottom" << endl;
  gpl << "set tmargin 2" << endl;
  gpl << "set bmargin 0" << endl;
  gpl << "unset xtics" << endl;
  gpl << "set xrange[0:" << rmax << "]" << endl;
  gpl << "set yrange[0:1.618]" << endl;
  gpl << "set ytics 0,0.5,1.5 format '%0.1f'" << endl;
  gpl << "set grid" << endl;
  gpl << "plot ";
  if( psf ) gpl << "'" << growname << "' t 'the genuine' with lines ls 98,";
  gpl << "'" << cname << "' t 'growth-curve' with lines ls 99, '" << ddname
      << "' t 'all stars' with points pt 6 lc rgb '#87CEFA', '" << dname
      << "' t 'designed by' with points pt 13 lc rgb 'blue'" << endl;
  gpl << "set origin 0.5,0" << endl;
  gpl << "set size 0.5,0.3" << endl;
  gpl << "unset key" << endl;
  gpl << "set grid" << endl;
  gpl << "set xlabel 'distance from a star center in pixels'" << endl;
  gpl << "set bmargin 3" << endl;
  gpl << "set tmargin 0" << endl;
  gpl << "set xtics nomirror" << endl;
  gpl << "set xrange[0:" << rmax << "]" << endl;

  // The second graph shows detail of grow-curve when PSF is available
  // (artificial data) and residuals for a real data.
  if( psf ) {
    gpl << "set yrange[0.95:1.01]" << endl;
    gpl << "set ytics 0.95,0.01,1.0 format '%0.2f'" << endl;
    gpl << "plot '" << growname << "' with lines ls 98,"
	<< "'" << dname  << "' with points pt 13 lc rgb 'blue'" << endl;
  }
  else {
    gpl << "set yrange[-0.21:0.21]" << endl;
    gpl << "set ytics -0.1,0.1,0.1 format '%0.1f'" << endl;
    gpl << "plot '" << ddname << "' u 1:3 with points pt 6 lc rgb '#87CEFA', '"
	<< dname << "' u 1:3 with points pt 13 lc rgb 'blue'" << endl;
  }
  gpl << "set origin 0,0" << endl;
  gpl << "set size 0.5,1" << endl;
  gpl << "set tmargin 2" << endl;
  gpl << "unset border" << endl;
  gpl << "unset xtics" << endl;
  gpl << "unset ytics" << endl;
  gpl << "set xrange[0:" << nx << "]" << endl;
  gpl << "set yrange[0:" << ny << "]" << endl;
  gpl << "unset xlabel" << endl;
  gpl << "plot '" << avsname << "' binary filetype=avs with rgbimage,'"
      << gsname << "' with circles lc rgb 'blue' fs transparent solid 0.2 noborder"
      << endl;

  string com("gnuplot");
  com = com + " " + gplname;
  int e = system(com.c_str());
  if( e != 0  )
    fprintf(stderr,"%s: gnuplot invoke failed.\n",__FILE__);

#ifndef DEBUG
  unlink(avsname); unlink(dname); unlink(ddname);
  unlink(gsname); unlink(cname); unlink(growname); unlink(gplname);
#endif

  close(fd); close(gd); close(dd); close(gs); close(gc); close(gf);
}


int main()
{
  string line;
  while( getline(cin,line) ) {
#ifdef DEBUG
    cout << line << endl;
#endif
    if( line.find("FILE") != string::npos && line.find("=") != string::npos ) {
      size_t l = line.find("'") + 1;
      size_t n = line.find("'",l) - l;
      report(line.substr(l,n));
    }
  }

  fprintf(stderr,"STOP 0\n");
  return 0;
}
