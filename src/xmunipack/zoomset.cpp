/*

  xmunipack - zoom set

  Copyright © 2021-22 F.Hroch (hroch@physics.muni.cz)

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

  A new numerial type MuniZoomSet is introduced.
  The purpose of the complex structure is handling of
  specific set of values for zoom 4:1, 2:1 ... 1:8,
  and to implement operations on the set.

*/


#include "zoomer.h"


MuniZoomSet::MuniZoomSet()
{
  Init();
}

void MuniZoomSet::Init()
{
  index = 5;
  nelements = 11;
  table = new double[nelements];
  table[0] = 1.0 /16.0;
  table[1] = 1.0 / 8.0;
  table[2] = 1.0 / 4.0;
  table[3] = 1.0 / 3.0;
  table[4] = 1.0 / 2.0;
  table[5] = 1.0;
  table[6] = 2.0;
  table[7] = 3.0;
  table[8] = 4.0;
  table[9] = 8.0;
  table[10]=16.0;
}

MuniZoomSet::MuniZoomSet(double zoom)
{
  Init();

  if( SetZoom(zoom) )
    ;
  else
    wxFAIL_MSG("MuniZoomSet(double) failed");
}

MuniZoomSet::MuniZoomSet(int i)
{
  Init();

  if( SetZoom(i) )
    ;
  else
    wxFAIL_MSG("MuniZoomSet(int) failed");
}


MuniZoomSet::~MuniZoomSet()
{
  delete[] table;
}

MuniZoomSet::MuniZoomSet(const MuniZoomSet& other)
{
  nelements = other.nelements;
  table = new double[nelements];
  std::copy(other.table,other.table+nelements,table);
}

MuniZoomSet& MuniZoomSet::operator = (const MuniZoomSet& other)
{
  if( this != &other ) {
    int n = other.nelements;
    double *t = new double[nelements];
    std::copy(other.table,other.table+nelements,t);

    if( nelements > 0 )
      delete[] table;

    nelements = n;
    table = t;
  }
  return *this;
}

int MuniZoomSet::GetIndex() const
{
  return index;
}

double MuniZoomSet::GetZoom()
{
  wxASSERT(0 <= index && index < nelements && table);
  return table[index];
}

bool MuniZoomSet::SetZoom(int i)
{
  if( 0 <= i && i < nelements ) {
    index = i;
    return true;
  }
  else
    return false;
}

bool MuniZoomSet::SetZoom(double z)
{
  int n = -1;
  for(int i = 0; i < nelements; i++) {
    if( fabs(z - table[i]) < 0.001 ) {
      n = i;
      break;
    }
  }

  bool found = n != -1;
  if( found ) index = n;

  return found;
}


double MuniZoomSet::SetBestFitZoom(const wxSize& window, const wxSize& image)
{
  index = 0;
  for(int i = 5; i >= 0; i--) {
    bool horizontal = window.GetWidth() >= image.GetWidth()*table[i];
    bool vertical = window.GetHeight() >= image.GetHeight()*table[i];
    if( horizontal && vertical ) {
      index = i;
      break;
    }
  }

  /*
    This part check whatever images with the selected zoom have both sizes
    non-zero. Otherwise, the zoom is reduced.
  */
  while( (int(image.GetWidth()*table[index]) < 1 ||
	  int(image.GetHeight()*table[index]) < 1 ) && table[index] < 1.001 )
    index++;

  wxLogDebug("SetBestFitZoom %d %d %d %d",index,index,
	     int(image.GetWidth()*table[index]),int(image.GetHeight()*table[index]));

  return table[index];
}

double MuniZoomSet::SetNormalZoom()
{
  for(int i = 0; i < nelements; i++)
    if( fabs(table[i] - 1) < 1e-3 ) {
      index = i;
      break;
    }

  return table[index];
}

double MuniZoomSet::Increase()
{
  if( index < nelements-2 )
    index++;

  return table[index];
}

double MuniZoomSet::Decrease()
{
  if( index > 0)
    index--;

  return table[index];
}
