/*

  Implements FITS file name utilities (taylored by wxFileName):

  * GNU backup-file conventions

  * file-existing


  Copyright © 2013, 2018 F.Hroch (hroch@physics.muni.cz)

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


#include "fitsfn.h"
#include <wx/wx.h>
#include <wx/regex.h>
#include <wx/filename.h>
#include <wx/dir.h>
#include <fitsio.h>

FITSFileName::FITSFileName(const wxString& name): filename(name)
{
  char *filter = 0, *binspec = 0, *colspec = 0, *pixspec = 0;
  char ftype[FLEN_FILENAME], infile[FLEN_FILENAME], outfile[FLEN_FILENAME],
    extspec[FLEN_FILENAME];
  const char *fn = filename.fn_str();

  status = 0;
  fits_parse_input_filename((char *)fn,ftype,infile,outfile,extspec,
			    filter,binspec,colspec,pixspec,&status);

  //  wxLogFatalError("%d %s %s %s %s "+filename,(int)status,filetype,
  //		  infile,outfile,extspec);

  if( status == 0 ) {
    filetype = ftype;
    fullpath = infile;
    tmpname = outfile;
    suffix = extspec;
  }
}

wxString FITSFileName::GetFullName() const
{
  wxFileName fname(fullpath);
  wxString output = fname.GetFullName();

  return output;
}
