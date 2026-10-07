/*

  xmunipack - dummy display for non-implemented frames

  Copyright © 2022 F.Hroch (hroch@physics.muni.cz)

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

#include "view.h"
#include <wx/wx.h>


MuniDummy::MuniDummy(wxWindow *w, const MuniConfig *c):
  wxPanel(w),config(c)
{
  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->AddStretchSpacer(3);

  wxStaticText *ch = new wxStaticText(this,wxID_ANY,"");
  ch->SetLabelMarkup(L"<span size=\"xx-large\">🚫</span>");
  topsizer->Add(ch,wxSizerFlags().Center().DoubleBorder());

  wxStaticText *label = new wxStaticText(this,wxID_ANY,"");
  label->SetLabelMarkup("<span fgcolor=\"black\" font-weight=\"bold\">An unimplemented feature</span>");
  topsizer->Add(label,wxSizerFlags().Centre().DoubleBorder());

  wxStaticText *description = new wxStaticText(this,wxID_ANY,"This kind of HDU can not be displayed as a regular picture.");
  description->SetForegroundColour(wxColour(128,128,128));
  topsizer->Add(description,wxSizerFlags().Center());
  topsizer->AddStretchSpacer(3);
  SetSizer(topsizer);
}
