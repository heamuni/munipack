/*

  xmunipack - tone profiles

  Copyright © 2009-2011, 2018-21 F.Hroch (hroch@physics.muni.cz)

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
#include <cmath>

FitsItt::FitsItt(): itt(ITT_LINE),F(&FitsItt::itt_line) {}

FitsItt::FitsItt(int i) { SetItt(i); }

void FitsItt::SetItt(int i)
{
  wxASSERT(ITT_FIRST < i && i < ITT_LAST);
  itt = i;

  switch(itt) {
  case ITT_LINE:  F = &FitsItt::itt_line; break;
  case ITT_SNLIKE:F = &FitsItt::itt_snlike; break;
  case ITT_ASINH: F = &FitsItt::itt_asinh; break;
  case ITT_TANH:  F = &FitsItt::itt_tanh; break;
  case ITT_PHOTO: F = &FitsItt::itt_photo;break;
  case ITT_SQR:   F = &FitsItt::itt_sqr;  break;
  }
}

void FitsItt::SetItt(const wxString& a)
{
  for(int i = ITT_FIRST+1; i < ITT_LAST; i++)
    if( a == Type_str(i) ) {
      SetItt(i);
      return;
    }
  wxFAIL_MSG("*** FitsItt::SetItt(const wxString&): unreachable line reached.");
}


int FitsItt::GetItt() const
{
  return itt;
}

wxString FitsItt::GetItt_str() const
{
  return Type_str(itt);
}

bool FitsItt::IsLinear() const
{
  return itt == ITT_LINE;
}

wxString FitsItt::Type_str(int n)
{
  switch(n){
  case ITT_LINE:  return "Linear";
  case ITT_SNLIKE:return "S/N like";
  case ITT_SQR:   return "Square";
  case ITT_ASINH: return "asinh";
  case ITT_TANH:  return "tanh";
  case ITT_PHOTO: return "Photo";
  default:        return "";
  }
}

wxArrayString FitsItt::Type_str()
{
  wxArrayString a;

  for(int i = ITT_FIRST+1; i < ITT_LAST; i++)
    a.Add(Type_str(i));

  return a;
}

float FitsItt::itt_line(float x) const
{
  return x;
}

float FitsItt::itt_snlike(float x) const
{
  float t = 6.66f*x;
  return t > 0.0f ? t / asinh(t) - 1.0f : 0.0f;
}

float FitsItt::itt_sqr(float x) const
{
  return x*x;
}

float FitsItt::itt_tanh(float x) const
{
  return tanh(x);
}

float FitsItt::itt_asinh(float x) const
{
  return asinh(x);
}

float FitsItt::itt_photo(float x) const
{
  return (1.0f + tanh(2.0f*(x - 1.0f)))/2.0f;
}


float *FitsItt::Scale(long n, const float *a) const
{
  float *f = new float[n];
  for(long i = 0; i < n; i++)
    f[i] = Scale(a[i]);
  return f;
}


float FitsItt::InvScale(float f) const
{
  wxASSERT(ITT_FIRST < itt && itt < ITT_LAST);

  if( itt == ITT_LINE )
    return f;
  else if( itt == ITT_SNLIKE )
    return (f - 1.0)/f;
  else if( itt == ITT_SQR )
    return f > 0.0 ? sqrt(f) : 0.0;
  else if( itt == ITT_TANH  )
    return fabs(f) < 1.0 ? atanh(f) : copysign(1.0,f);
  else if( itt == ITT_ASINH )
    return sinh(f);
  else if( itt == ITT_PHOTO  ) {
    if( f < 0.0 )
      return 0.0;
    else if( f > 1.0 )
      return 1.0;
    else
      return 1.0 + atanh(2.0*f-1.0) / 2.0;
  }

  return 1.0;
}
