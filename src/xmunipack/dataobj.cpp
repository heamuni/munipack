/*

  xmunipack - implementation of data objects for Dnd and clipboard

  Copyright © 2010-2011 F.Hroch (hroch@physics.muni.cz)

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
#include <vector>
#include <wx/wx.h>
#include <wx/dataobj.h>
#include <wx/mstream.h>

using namespace std;

class MuniDataFormat: public wxDataFormat
{
public:
  MuniDataFormat() { SetId("MUNIPACK_METAFITS"); }
};


// ---- MuniDataObjectMeta

MuniDataObjectMeta::MuniDataObjectMeta(): 
  len(0),data(0)
{
  SetFormat(MuniDataFormat());
}

MuniDataObjectMeta::MuniDataObjectMeta(const vector<FitsMeta>& mlist): 
  len(0),data(0)
{
  SetFormat(MuniDataFormat());

  wxMemoryOutputStream b;
  for(vector<FitsMeta>::const_iterator i = mlist.begin();i != mlist.end(); ++i){
    MuniThumbnail th(*i);
    th.Save(b);
  }

  len = b.GetSize();
  wxASSERT(len > 0);
  data = new char[len];
  b.CopyTo(data,len);
}

MuniDataObjectMeta::~MuniDataObjectMeta()
{
  delete[] data;
}


// copy & assignement consttructor?

size_t MuniDataObjectMeta::GetDataSize() const
{
  return len;
}

bool MuniDataObjectMeta::GetDataHere(void *d) const
{
  wxASSERT(data && len > 0);
  memcpy(d,data,len);
  return true;
}

bool MuniDataObjectMeta::SetData(size_t size, const void *d)
{
  delete[] data;
  data = new char[size];
  len = size;
  memcpy(data,d,size);

  return true;
}


vector<FitsMeta> MuniDataObjectMeta::GetMetafitses() const
{
  wxMemoryInputStream b(data,len);

  vector<FitsMeta> mlist;
  while(true) {
    MuniThumbnail th(b);
    if( ! th.IsOk() ) break;
    mlist.push_back(th.GetMeta());
  }

  return mlist;
}

