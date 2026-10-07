/*

  FITS related utility


  Copyright © 2011-3, 2017-8 F.Hroch (hroch@physics.muni.cz)

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
#include "fortranio.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <fitsio.h>

#include <cassert>

using namespace std;


// removes trailing spaces, "IMAGETYP='LIGHT T '  " -> "IMAGETYP='LIGHT T'"
string rtrim(const string& str)
{
  char *a = strdup(str.c_str());

  if( strlen(a) > 0 ) {

    // removes trailing spaces
    for(char *l = a + strlen(a)-1; l != a; l--) {
      if( *l == ' ' )
	*l = '\0';
      else
	break;
    }

    // removes trailing spaces inside appostrophes
    char *l = a + strlen(a);
    while( l-- != a ) {
      if( *l == ' ' && *(l+1) == '\'' ) {
	*l = '\''; *(l+1)='\0';
      }
    }
  }

  string out(a);
  free(a);

  return out;
}


int type_recognize(const string& value)
{
  double x;

  if( value == "T" || value == "F" )
    return TYPE_LOGICAL;

  if( sscanf(value.c_str(),"%lf",&x) == 1 ) {
    for(string::const_iterator l = value.begin(); l != value.end(); l++)
      if( *l == '.' || *l == 'e' || *l == 'E' )
	return TYPE_REAL;

    return TYPE_INT;
  }

  return TYPE_STRING;
}



int header_print(const string& filename, const vector<string>& keywords, int keylist)
{
  fitsfile *f;
  int status = 0;

  status = 0;
  if( fits_open_file(&f, filename.c_str(), READONLY, &status) )
    fits_report_error(stderr, status);

  if( ! keywords.empty() ) {
    /* print selected keys */

    int ninc = keywords.size();
    char **inclist = (char **) malloc(sizeof(char *)*ninc);

    for(int i = 0; i < ninc; i++)
      inclist[i] = strdup(keywords[i].c_str());

    char record[FLEN_CARD];
    fits_read_record(f,0,record,&status); // useless, but resets init position

    if( keylist == ID_FULL ) {

      while( fits_find_nextkey(f,inclist,ninc,NULL,0,record,&status) == 0 )
	cout << record << endl;

    }
    else if( keylist == ID_SHELL || keylist == ID_VALUE ) {

      int keylen;
      char keyname[FLEN_KEYWORD],value[FLEN_VALUE],comment[FLEN_COMMENT];

      while( fits_find_nextkey(f,inclist,ninc,NULL,0,record,&status) == 0 ) {
	fits_get_keyname(record,keyname,&keylen,&status);
	fits_parse_value(record,value,comment,&status);
	if( keylist == ID_SHELL )
	  cout << keyname << "=" << rtrim(value) << endl;
	else if( keylist == ID_VALUE )
	  cout << rtrim(value) << endl;
      }
    }

    for(int i = 0; i < ninc; i++)
      free(inclist[i]);
    free(inclist);

    if( status == KEY_NO_EXIST )
      status = 0;

  }
  else {

    /* print full header */
    int nhead;
    char record[FLEN_CARD];

    fits_get_hdrspace(f,&nhead,NULL,&status);
    for(int i = 0; status == 0 && i < nhead; i++) {
      if( fits_read_record(f,i+1,record,&status) == 0 )
	cout << record << endl;
    }

  }

  fits_close_file(f, &status);

  if( status != 0 )
    fits_report_error(stderr, status);

  return status;
}

int header_remove(const string& filename, const vector<string>& keywords)
{
  fitsfile *f;
  int status = 0;

  status = 0;
  if( fits_open_file(&f, filename.c_str(), READWRITE, &status) )
    fits_report_error(stderr, status);

  vector<string>::const_iterator i;
  for( i = keywords.begin(); i < keywords.end(); i++ ) {
    if( fits_delete_key(f,i->c_str(),&status) != 0 ) {
      fits_report_error(stderr, status);
      status = 0;
    }
  }

  fits_close_file(f, &status);

  if( status != 0 )
    fits_report_error(stderr, status);

  return status;
}

int header_update(const string& filename, const string& keyword,
		  const string& value, const string& comment)
{
  fitsfile *f;
  int status = 0;

  int type = type_recognize(value);

  status = 0;
  if( fits_open_file(&f, filename.c_str(), READWRITE, &status) )
    fits_report_error(stderr, status);

  if( type == TYPE_INT ) {
    stringstream s(value);
    int i;
    s >> i;
    fits_update_key(f,TINT,const_cast<char*>(keyword.c_str()),&i,
		    const_cast<char*>(comment.c_str()),&status);
  }
  else if( type == TYPE_REAL ) {
    stringstream s(value);
    double x;
    s >> x;
    fits_update_key(f,TDOUBLE,const_cast<char*>(keyword.c_str()),&x,
		    const_cast<char*>(comment.c_str()),&status);
  }
  else if( type == TYPE_LOGICAL ) {
    stringstream s(value);
    int i = value == "T" ? 1 : 0;
    fits_update_key(f,TLOGICAL,const_cast<char*>(keyword.c_str()),&i,
		    const_cast<char*>(comment.c_str()),&status);
  }
  else if( type == TYPE_STRING ) {

    string val = forstr(value);
    fits_update_key(f,TSTRING,const_cast<char*>(keyword.c_str()),
		    const_cast<char*>(val.c_str()),
		    const_cast<char*>(comment.c_str()),&status);
  }

  fits_close_file(f, &status);

  if( status != 0 )
    fits_report_error(stderr, status);


  return status;
}

int header_template(const string& filename, const string& templ)
{
  fitsfile *f;
  int status = 0;
  string line;
  int keytype, keylen;
  char card[FLEN_CARD], key[FLEN_KEYWORD];

  status = 0;
  if( fits_open_file(&f, filename.c_str(), READWRITE, &status) )
    fits_report_error(stderr, status);

  ifstream fin(templ.c_str(),ifstream::in);
  if( ! fin.good() ) {
    cerr << "Failed to open `" << filename << "'." << endl;
    return 1;
  }

  while( fin.good() && status == 0 ) {

    getline(fin,line);
    if( fin.eof() ) break;

    fits_parse_template((char *)line.c_str(),&card[0],&keytype,&status);
    if( status != 0 )
      cerr << "Failed to parse `" << line << "'." << endl;
    fits_get_keyname(card,&key[0],&keylen,&status);
    fits_update_card(f,key,card,&status);

  }

  fits_close_file(f, &status);

  if( status != 0 )
    fits_report_error(stderr, status);

  return status;
}
