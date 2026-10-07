/*

  Astrometrical calibration of FITS image

  Copyright © 2010-3 F.Hroch (hroch@physics.muni.cz)

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

  1. create list to match:
     * projection from catalogue
     * photometry of ref
  2. for every image:
     * photometry
     * match
     * astrometry

*/

#include "munipack.h"
#include <wx/wx.h>
#include <wx/filename.h>
#include <iostream>
#include <map>

using namespace std;

//#define TESTING 1

#ifdef TESTING

void Munipack::calibrate(MuniCmdLineParser& cmd)
{
  MuniProcess *p = new MuniProcess(&pipe,"uname -a");
  pipe.push(p);

  MuniProcess *f = new MuniProcess(&pipe,"date");
  pipe.push(f);

  MuniProcess *g = new MuniProcess(&pipe,"sleep 3");
  pipe.push(g);

  MuniProcess *c = new MuniProcess(&pipe,"cat");
  c->Write("ssss");
  pipe.push(c);

  pipe.Start();

  return true;
}

#else

void Munipack::calibrate(MuniCmdLineParser& cmd)
{
  size_t cidx = 0;
  wxString tmpname = wxFileName::CreateTempFileName("munipack-calibrate-");
  if( tmpname.IsEmpty() )
    wxLogFatalError("Failed to create temporary files.");
    
  MuniProcess *action = new MuniProcess(&pipe,"cone");
  pipe.push(action);

  action->Write("OUTPUT = '!" + tmpname + "'");
  action->Write("TYPE = 'FITS'");
  action->Write("SORT = 'R'");
  action->Write("URL = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=UCAC-3&'");

  wxString a;
  if( cmd.Found("r",&a) )
    action->Write("SR = '" + a + "'");
  else
    action->Write("SR = '0.15'");

  for(size_t i = 1; i < cmd.GetParamCount(); i++) {
    wxString p(cmd.GetParam(i));
    int pos = p.Find(',');
    if( pos != wxNOT_FOUND ) {
      cidx = i;
      action->Write("RA = '" + p.substr(0,pos) + "'");
      action->Write("DEC = '" + p.substr(pos+1) + "'");
    }
  }

  MuniProcess *afit = new MuniProcess(&pipe,"astrofit");
  pipe.push(afit);

  afit->Write("CAT = '"+tmpname+"'");

  for(size_t i = 1; i < cmd.GetParamCount(); i++)
    if( i != cidx )
      afit->Write("FILE = '" + cmd.GetParam(i) + "'");

  MuniProcess *rm = new MuniProcess(&pipe,"rm -f "+tmpname);
  pipe.push(rm);  

  pipe.Start();

  return true;
}

#endif

// int main()
// {

//   if( catalogue ) {

//     // cone search
//     MuniOperation cone(ra,dec,..);
//     if( ! catalog.IsOk() )
//       return 1;
    
//     // convert VO table to txt
//     MuniOperation votxt();
    

//     // projection of sperical to rectangular coordinates
//     MuniOperation proj();

//   }
//   else if ( ref ) {
    
//     // photometry of reference image
//     MuniOperation ref();

//   }
  

//   //  forall

//   // looking for stars
//   MuniOperation photo();

//   // merge lists

//   // match to list
//   MuniOperation match();

//   // preliminary astrometry
//   MuniOperation astrometry();

//   // find all possible stars

//   // final astrometry
//   MuniOperation astrometry();

//   // calibrate images
//   MuniOperation fits();


//   // finish

