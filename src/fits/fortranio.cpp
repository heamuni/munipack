/*

  Fortran text parsers


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

#include "fortranio.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>

using namespace std;


// string tokenizer recognizing Fortran strings
//  * delimiters are by norm, only ' ,;\t' acceptable
//  * double apostrophes are treated as the single character
//  * strings starts and finish with a single apostrophe character
vector<string> strftok(const string& fstring, const string& delimiters)
{
  const char *a = fstring.c_str();
  const char *delim = delimiters.size() > 0 ? delimiters.c_str() : " ,;\t\n\0";
  vector<string> tokens;

  char *b = new char[strlen(a)+1];

  bool open = false;
  bool str = false;
  bool fin = false;
  char apho = ' ';

  int l = 0;
  int i = 0;
  for(;;) {

    if( strchr(delim,a[i]) != NULL ) {

      if( str )
	// the separator inside string appostrophes
	b[l++] = a[i];
      else {
	if( open )
	  // non-string value finished
	  fin = true;
      }
    }
    else {


      if( ((a[i] == '\'' || a[i] == '"') && apho == ' ') ||
	  (apho != ' ' && a[i] == apho)  ) {
	apho = a[i];

	// appostroph encountered, what's now?
	if( ! str ) {
	  str = true;
	  open = true;
	}
	else {
	  // inside string

	  if( a[i+1] != '\0' && ( a[i+1] == '\'' || a[i+1] == '"') ) {
	    //	  if( a[i+1] != '\0' &&  apho != ' ' && a[i+1] == apho ) {
	    // inside string, double appostrophes are treated as a single char
	    i++;
	    b[l++] = a[i];
	  }
	  else if( a[i+1] != '\0' && strchr(delim,a[i+1]) != NULL ) {
	    // finishing string
	    str = false;
	    open = false;
	    fin = true;
	  }
	  else if( a[i+1] == '\0' ) {
	    // finishing string
	    str = false;
	    open = false;
	    fin = true;
	  }
	}
      }
      else {
	// non-empty, separator character

	open = true;
	b[l++] = a[i];
      }

    }

    if( fin ) {
      b[l++] = '\0';
	//	cerr << "number:" << b << endl;
      tokens.push_back(b);
      l = 0;
      open = false;
      str = false;
      fin = false;
      apho = ' ';
    }

    if( a[i] == '\0' ) break;
    i++;
  }

  delete[] b;

  return tokens;
}

string forstr(const string& a)
{
  vector<string> s = strftok(a," \t\n\0");
  if( s.size() == 1)
    return s[0];
  else
    return "";
}
