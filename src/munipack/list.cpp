/*

  Listing

  Copyright © 2011-2, 2018 F.Hroch (hroch@physics.muni.cz)

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

// listing of all usefull info from FITS

// light curve, absolute, relative  magnitudes, intensity, sky, jd,
//              airmass, helcor ...?
// all calibrated !


#include "munipack.h"
#include <wx/wx.h>
#include <wx/tokenzr.h>
#include <wx/intl.h>
#include <vector>

using namespace std;

bool Munipack::list(MuniCmdLineParser& cmd)
{
  MuniProcess *action = new MuniProcess(&pipe,"list");
  pipe.push(action);

  wxString cooformat("SPHERICAL");
  wxString a;
  double x;
  wxArrayString cols,keys;
  vector<size_t> icoo;

  // detrmine character used as decimal point
  wxString decimal = wxLocale::GetInfo(wxLOCALE_DECIMAL_POINT,wxLOCALE_CAT_NUMBER);
  wxString sep = decimal == "." ? ",;" : ";";

  CommonOutputSingle(action,cmd);

  if( cmd.Found("c",&a) || cmd.Found("create",&a) )
    action->Write("MODE = '" + a + "'");

  /*
  if( cmd.Found("coo",&a) )
    action << ("COOTYPE = '" + a + "'");
  */

  if( cmd.Found("col",&a) ) {
    wxStringTokenizer tokenizer(a, ",");
    while ( tokenizer.HasMoreTokens() )
      cols.Add(tokenizer.GetNextToken());
  }

  if( cmd.Found("key",&a) ) {
    wxStringTokenizer tokenizer(a, ",");
    while ( tokenizer.HasMoreTokens() )
      keys.Add(tokenizer.GetNextToken());
  }

  //    action << ("COLUMN = '" + a + "'");

  /*
  if( cmd.Found("mag") )
    action << wxString("MAG = 1");

  if( cmd.Found("full") )
    action << wxString("FULL_LISTING = 1");

  if( cmd.Found("pos") )
    action << wxString("POSITION_LISTING = 1");
  */

  /*
  if( cmd.Found("relative",&a) )
    action << wxString("RELATIVE = 1");
  */

  /*
  if( cmd.Found("geo",&a) )
    action << ("GEOGRAPHICAL = " + a);
  */

  /*
  if( cmd.Found("longitude",&x) )
    action->Write("LONGITUDE = %f",x);

  if( cmd.Found("latitude",&x) )
    action->Write("LATITUDE = %f",x);

  if( cmd.Found("ra",&x) )
    action->Write("RA = %f",x);

  if( cmd.Found("dec",&x) )
    action->Write("DEC = %f",x);
  */


  /*
  if( cmd.Found("longitude",&x) )
    action->Write("LONGITUDE = %f",x);

  if( cmd.Found("latitude",&x) )
    action->Write("LATITUDE = %f",x);

  if( cmd.Found("alpha",&x) )
    action->Write("ALPHA = %f",x);

  if( cmd.Found("delta",&x) )
    action->Write("DELTA = %f",x);

  if( cmd.Found("epoch",&x) )
    action->Write("EPOCH = %f",x);

  if( cmd.Found("period",&x) )
    action->Write("PERIOD = %f",x);
  */

  /*
  if( cmd.Found("equ",&a) )
    action << ("EQUATORIAL = " + a);
  */

  /*
  if( cmd.Found("ele",&a) )
    action << ("LC_ELEMENTS = " + a);
  */

  if( cmd.Found("tol",&x) )
    action->Write("TOL = %f",x);

  if( cmd.Found("print-filename") )
    action->Write("PRINT_FILENAME = T");

  if( cmd.Found("mag") )
    action->Write("MAGNITUDES = T");

  if( cmd.Found("flux") )
    action->Write("FLUXES = T");

  if( cmd.Found("calibr") )
    action->Write("CALIBRATED = T");

  if( cmd.Found("diffmag") )
    action->Write("DIFFMAG = T");

  if( cmd.Found("zeromag",&x) )
    action->Write("ZEROMAG = %f",x);

  if( cmd.Found("photflux",&x) )
    action->Write("PHOTFLUX = %f",x);

  if( cmd.Found("epoch",&x) )
    action->Write("EPOCH = %f",x);

  if( cmd.Found("aperture",&x) )
    action->Write("APERTURE = %f",x);

  if( cmd.Found("s") || cmd.Found("spherical") )
    cooformat = "SPHERICAL";

  if( cmd.Found("p") || cmd.Found("pixels") )
    cooformat = "CARTESIAN";

  /*
  if( cmd.Found("zeromag",&x) )
    action->Write("ZEROMAG = %f",x);

  if( cmd.Found("photflux",&x) )
    action->Write("PHOTFLUX = %f",x);
  */


  /*
  if( cmd.Found("index") )
    cooformat = "INDEX";
  */

  /*
  if( cmd.Found("k",&a) ) {

    action << wxString("HEADER = 1");

    wxStringTokenizer tokenizer(a, ",");
    while ( tokenizer.HasMoreTokens() )
      action << ("KEYWORD = '" + tokenizer.GetNextToken() + "'");
  }
  */

  if( cols.GetCount() > 0 ) {
    action->Write("NCOLUMNS = %ld",long(cols.GetCount()));
    for(size_t i = 0; i < cols.GetCount(); i++)
      action->Write("COLUMNS = '" + cols[i] + "'");
  }

  if( keys.GetCount() > 0 ) {
    action->Write("NKEYWORDS = %ld",long(keys.GetCount()));
    for(size_t i = 0; i < keys.GetCount(); i++)
      action->Write("KEYWORDS = '" + keys[i] + "'");
  }

  wxASSERT(!cooformat.IsEmpty());

  if( cooformat == "INDEX" ) {

    vector<long> idx;

    for(size_t i = 1; i < cmd.GetParamCount(); i++) {
      wxString l(cmd.GetParam(i));
      long d;
      if( l.ToLong(&d) ) {
	idx.push_back(d);
	icoo.push_back(i);
      }
    }

    if( idx.size() > 0 ) {
      action->Write("NINDEX = %ld",long(idx.size()));
      for(size_t i = 0; i < idx.size(); i++)
	action->Write("INDEX = %ld",idx[i]);
    }
  }

  if( cooformat == "CARTESIAN" || cooformat == "SPHERICAL" ) {

    vector<double> x,y,vx,vy;
    vector<double> w;

    for(size_t i = 1; i < cmd.GetParamCount(); i++) {
      wxString l(cmd.GetParam(i));

      w.clear();
      wxStringTokenizer tokenizer(l,sep);
      while( tokenizer.HasMoreTokens() ){
        wxString token = tokenizer.GetNextToken();
	double q;
	if( token.ToDouble(&q) )
	  w.push_back(q);
      }

      if( w.size() == 2 ) {
	x.push_back(w[0]);
	y.push_back(w[1]);
	vx.push_back(0.0);
	vy.push_back(0.0);
	icoo.push_back(i);
      }
      else if( w.size() == 4 ) {
	x.push_back(w[0]);
	y.push_back(w[1]);
	vx.push_back(w[2]);
	vy.push_back(w[3]);
	icoo.push_back(i);
      }

    }

    if( x.size() > 0 ) {
      action->Write(wxString("N"+cooformat+" = %ld"),long(x.size()));
      for(size_t i = 0; i < x.size(); i++) {
	/*
	wxString a;
	a.Printf(wxString(),
		 x[i],y[i],vx[i],vy[i]);
	*/
	action->Write(cooformat+" = %.15f %.15f %.15e %.15e",x[i],y[i],vx[i],vy[i]);
      }
    }
  }

  vector<wxString> files;
  for(size_t i = 1; i < cmd.GetParamCount(); i++) {

    bool found = false;
    for(size_t j = 0; j < icoo.size(); j++)
      if( icoo[j] == i )
	found = true;
    if( ! found )
      files.push_back(cmd.GetParam(i));


    // wxString l = cmd.GetParam(i);
    // long d;

    // if( l.Find(',') != wxNOT_FOUND )
    //   action << ("RECT = " + l + ",0,0");

    // else if( l.ToLong(&d) )
    //   action << ("INDEX = " + l);

    // else
    //   action->WriteFiles(l);
    //   //      action << ("FILE = '" + l + "'");
  }

  if( ! files.empty() ) {
    //action->Write("NFILES = %ld",long(files.size()));
    action->Write("NFILES = %d",static_cast<int>(files.size()));
    for(size_t i = 0; i < files.size(); i++)
      action->WriteFile(files[i]);
    //      WriteFile(action,files[i]);
  }

  //  pipe.Start();

  return true;
}
