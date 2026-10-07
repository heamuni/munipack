/*

  FITS related utility


  Copyright © 2011-8 F.Hroch (hroch@physics.muni.cz)

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


  Plan:

  * implement adding and removing COMMENT,HISTORY records


*/

#include "fits.h"
#include "mfitsio.h"
#include "fortranio.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

int main()
{
  int mode = ID_NULL;
  int hmode = ID_MODE_NULL;
  string key, val, com, templ;
  vector<string> keywords, extnames;
  int keylist = ID_FULL;
  bool exitus = true;

  while( cin.good() ) {

    string line;
    getline(cin,line);

    if( line.length() == 0 )
      continue;

    size_t eq = line.find('=');
    if( eq == string::npos ) {
      cerr  << line << endl;
      cerr << "STOP 'Malformed input record.'" << endl;
      return 1;
    }

    string value = line.substr(eq+1);

    if( line.find("MODE") != string::npos ) {

      if( value.find("STRUCTURE") != string::npos )
	mode = ID_STRUCTURE;
      else if ( value.find("HEADER") != string::npos )
	mode = ID_HEADER;
      else if ( value.find("TABLE") != string::npos )
	mode = ID_TABLE;
      else if ( value.find("IMAGE") != string::npos )
	mode = ID_IMAGE;
      else if ( value.find("DUMP") != string::npos )
	mode = ID_DUMP;
      else if ( value.find("RESTORE") != string::npos )
	mode = ID_RESTORE;
      else if ( value.find("CAT") != string::npos )
	mode = ID_CAT;
      else if ( value.find("EXTENSION") != string::npos )
	mode = ID_EXTENSION;

    }

    if( line.find("KEYWORD") != string::npos ) {

      hmode = ID_MODE_PRINT;
      size_t i = value.find('\'');
      size_t j = value.rfind('\'');
      size_t n = j - i - 1;
      keywords.push_back(value.substr(i+1,n));

    }

    if( line.find("REMKEY") != string::npos ) {

      hmode = ID_MODE_REMOVE;
      size_t i = value.find('\'');
      size_t j = value.rfind('\'');
      size_t n = j - i - 1;
      keywords.push_back(value.substr(i+1,n));

    }

    if( line.find("UPDATE") != string::npos ) {

      hmode = ID_MODE_UPDATE;

    }

    if( line.find("KEY") != string::npos ) {

      key = forstr(value);

    }

    if( line.find("VALUE") != string::npos ) {

      val = forstr(value);

    }

    if( line.find("COMMENT") != string::npos ) {

      com = forstr(value);

    }

    if( line.find("TEMPL") != string::npos ) {

      templ = forstr(value);

    }

    if( line.find("REMEXT") != string::npos ) {

      hmode = ID_MODE_REMOVE;
      size_t i = value.find('\'');
      size_t j = value.rfind('\'');
      size_t n = j - i - 1;
      extnames.push_back(value.substr(i+1,n));

    }

    if( line.find("KEYLIST") != string::npos ) {

      if( value.find("SHELL") != string::npos )
	keylist = ID_SHELL;

      else if( value.find("VALUE") != string::npos )
	keylist = ID_VALUE;

    }

    if( line.find("FILE") != string::npos &&
	line.find("NFILES") == string::npos ) {

      vector<string> items(strftok(value));
      assert(items.size() == 2);

      string filename(items[0]);
      string output(items[1]);

      // begin of this file processing
      int ret = -1;

      if( mode == ID_STRUCTURE )
	ret = structure(filename);

      else if( mode == ID_HEADER ) {
	if( hmode == ID_MODE_PRINT || hmode == ID_MODE_NULL )
	  ret = header_print(filename,keywords,keylist);
	else if( hmode == ID_MODE_REMOVE || hmode == ID_MODE_UPDATE ) {

	  string name;
	  if( filename == output )
	    name = filename;
	  else {
	    mfitsio_unlink(output);
	    mfitsio_copy(filename,output);
	    name = output;
	  }

	  if( hmode == ID_MODE_REMOVE )
	    ret = header_remove(name,keywords);
	  else if( hmode == ID_MODE_UPDATE ) {
	    if( templ == "" )
	      ret = header_update(name,key,val,com);
	    else
	      ret = header_template(name,templ);
	  }
	}
      }

      else if( mode == ID_TABLE )
	ret = table(filename);

      else if( mode == ID_IMAGE )
	ret =  image(filename);

      else if( mode == ID_DUMP )
	ret =  dump(filename,output);

      else if( mode == ID_RESTORE )
	ret =  restore(filename,output);

      else if( mode == ID_EXTENSION )
	ret =  ext_remove(filename,extnames,output);

      else if( mode == ID_CAT )
	ret =  cat(filename,output);

      if( ret != 0 ) {
	cerr << "Processing of `" << filename <<
	  "' failed with FITS return code: " << ret << endl;
	exitus = false;
      }
      // end of process file
    }
  }

  if( exitus ) {
    cerr << "STOP 0" << endl;
    return 0;
  }
  else {
    cerr << "STOP 'An error in fits occurred.'" << endl;
    return 1;
  }
}
