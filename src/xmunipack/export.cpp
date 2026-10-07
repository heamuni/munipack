/*

  XMunipack -- export FITS images dialog

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

#include "fits.h"
#include "fitsdisplay.h"
#include "view.h"
#include "export.h"
#include <wx/wx.h>

// -- FileShrinkControls
// helper class with the shrink option

class FileShrinkControls: public wxPanel
{
  double shrink;

  void OnChoice(wxCommandEvent&);
public:
  FileShrinkControls(wxWindow *);
  double GetShrink() const { return shrink; }
};

FileShrinkControls::FileShrinkControls(wxWindow *w): wxPanel(w), shrink(1)
{
  wxSizerFlags flags;
  flags.Align(wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT);

  wxBoxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);
  /*  sizer->Add(new wxStaticText(this,wxID_ANY,"Size:"),flags);*/

  wxString *shrinks = new wxString[6];
  shrinks[0] = "Original size";
  shrinks[1] = L"Shrink 2𐄂";
  shrinks[2] = L"Shrink 3𐄂";
  shrinks[3] = L"Shrink 4𐄂";
  shrinks[4] = L"Shrink 8𐄂";
  shrinks[5] = L"Shrink 16𐄂";

  wxChoice *ch =
    new wxChoice(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,6,shrinks);
  ch->SetSelection(0);
  sizer->Add(ch,flags);
  SetSizerAndFit(sizer);

  Bind(wxEVT_CHOICE,&FileShrinkControls::OnChoice,this);
}

void FileShrinkControls::OnChoice(wxCommandEvent& event)
{
  switch(event.GetSelection()) {
  case 0: shrink = 1; break;
  case 1: shrink = 2; break;
  case 2: shrink = 3; break;
  case 3: shrink = 4; break;
  case 4: shrink = 8; break;
  case 5: shrink = 16; break;
  }
}

static wxWindow *CreateShrinkChoice(wxWindow *w)
{
  return new FileShrinkControls(w);
}


// ----  MuniExportDialog --------------------------------------

MuniExportDialog::MuniExportDialog(wxWindow *w, const wxString& t, const wxString& d,
				   const wxString& fn, const wxString& ms, long s):
  wxFileDialog(w,t,d,fn,ms,s)
{
  SetExtraControlCreator(&CreateShrinkChoice);
}

bool MuniExportDialog::Save(const FitsArray& array, const MuniDisplay *mdisplay)
{
  int shrink = static_cast<FileShrinkControls*>(GetExtraControl())->GetShrink();

  FitsDisplay display(array.IsColour());
  display.SetShrink(shrink);
  display.SetTone(mdisplay->GetTone());
  display.SetItt(mdisplay->GetItt());
  display.SetPalette(mdisplay->GetPalette());
  display.SetColour(mdisplay->GetColour());
  wxImage image = display.ConvertTowxImage(array);
  return image.SaveFile(GetFilename());
}
