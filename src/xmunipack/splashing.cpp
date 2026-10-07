/*

  xmunipack - splasing


  Copyright © 2018, 2021 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/wx.h>
#include <wx/animate.h>


MuniSplashing::MuniSplashing(wxWindow *w, const MuniConfig *config):
  wxWindow(w,wxID_ANY),logo(wxBitmap(config->munipack_icon)),anim(0)
{
  //    wxLogDebug("Splashing...");

  wxFont bf(*wxNORMAL_FONT);
  bf.SetWeight(wxFONTWEIGHT_BOLD);
  bf.Scale(2.0);

  anim = new wxAnimationCtrl(this,wxID_ANY,config->throbber);
  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->AddStretchSpacer();
  topsizer->Add(new wxStaticBitmap(this,wxID_ANY,logo),
		wxSizerFlags().Center().Border());
  wxStaticText *label = new wxStaticText(this,wxID_ANY,"Munipack");
  label->SetFont(bf);
  label->SetForegroundColour(wxColour(128,128,128));
  topsizer->Add(label,wxSizerFlags().Center().DoubleBorder());
  topsizer->Add(anim,wxSizerFlags().Center().Border());
  topsizer->AddStretchSpacer();
  SetSizer(topsizer);
  anim->Show(false);
}

void MuniSplashing::Play()
{
  wxASSERT(anim);
  anim->Play();
  anim->Show(true);
  Layout();
}

void MuniSplashing::Stop()
{
  wxASSERT(anim);
  anim->Stop();
  anim->Show(false);
  Layout();
}
