/*

  FITS related utility


  Copyright © 2011-3, 2016-8 F.Hroch (hroch@physics.muni.cz)

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

#include <string>
#include <vector>

enum { ID_NULL, ID_STRUCTURE, ID_HEADER, ID_TABLE, ID_IMAGE, ID_DUMP,
       ID_RESTORE, ID_EXTENSION, ID_CAT };
enum { ID_MODE_NULL, ID_MODE_PRINT, ID_MODE_REMOVE, ID_MODE_UPDATE };
enum { ID_FULL, ID_SHELL, ID_VALUE };
enum { TYPE_STRING, TYPE_LOGICAL, TYPE_INT, TYPE_REAL };

int main();
int structure(const std::string&);
int header_print(const std::string&, const std::vector<std::string>&,int);
int header_remove(const std::string&, const std::vector<std::string>&);
int header_update(const std::string&, const std::string&, const std::string&,
		  const std::string&);
int header_template(const std::string&, const std::string&);
int table(const std::string&);
int image(const std::string&);
int dump(const std::string&,const std::string&);
int restore(const std::string&,const std::string&);
int cat(const std::string&,const std::string&);
std::string rtrim(const std::string&);
int type_recognize(const std::string&);
int ext_remove(const std::string&, const std::vector<std::string>&,
	       const std::string&);
