/*

  xmunipack - head


  Copyright © 2009-2011, 2020-22 F.Hroch (hroch@physics.muni.cz)

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
#include "help.h"
#include <wx/wx.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>

MuniHead::MuniHead(wxWindow *w, const MuniConfig *c):
  wxTextCtrl(w,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
	     wxTE_MULTILINE|wxTE_READONLY|wxTE_DONTWRAP),config(c)
{
  const wxFont *ft = wxNORMAL_FONT;
  wxFont tt = wxFont(ft->GetPointSize(),wxFONTFAMILY_MODERN,ft->GetStyle(),
	      ft->GetWeight());
  wxTextAttr ta;
  ta.SetFont(tt);
  SetDefaultStyle(ta);
}


bool MuniHead::SetHdu(const FitsHdu& h)
{
  head = h;

  wxASSERT(head.IsOk());
  if( ! head.IsOk() ) return false;

  Clear();
  for(size_t i = 0; i < head.GetCount(); i++) {

      wxString line(head.Item(i));

      SetDefaultStyle(wxTextAttr(*wxBLUE));
      AppendText(line.SubString(0,7));

      SetDefaultStyle(wxTextAttr(*wxBLACK));
      AppendText(line.SubString(8,line.Len()));

      AppendText("\n");
  }

  ShowPosition(0);
  Layout();

  return true;
}


// -----   MuniHeader ---------------------------------

MuniHeader::MuniHeader(wxWindow *w, MuniConfig *c):
  wxFrame(w,wxID_ANY,"Header",wxDefaultPosition,c->header_size),config(c)
{
  SetIcon(config->munipack_icon);

  menuFile = new wxMenu;
  menuFile->Append(ID_EXPORT,"Export As...");
  menuFile->AppendSeparator();
  menuFile->Append(wxID_CLOSE);

  /*
  menuEdit =  new wxMenu;
  menuEdit->Append(wxID_CUT);
  menuEdit->Append(wxID_COPY);
  menuEdit->Append(wxID_PASTE);
  */

  wxMenu *menuHelp = new wxMenu;
  menuHelp->Append(wxID_ABOUT);

  wxMenuBar *menuBar = new wxMenuBar;
  menuBar->Append(menuFile,"&File");
  //  menuBar->Append(menuEdit,"&Edit");
  menuBar->Append(menuHelp,"&Help");
  SetMenuBar(menuBar);

  header = new MuniHead(this,config);

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(header,wxSizerFlags(1).Expand());
  SetSizer(topsizer);
  SetSizeHints(wxSize(30*GetCharWidth(),7*GetCharHeight()));

  // init state
  /*
  menuFile->Enable(ID_EXPORT,false);
  menuEdit->Enable(wxID_CUT,false);
  menuEdit->Enable(wxID_COPY,false);
  menuEdit->Enable(wxID_PASTE,false);
  */

  Bind(wxEVT_CLOSE_WINDOW,&MuniHeader::OnClose,this);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniHeader::FileClose,this,wxID_CLOSE);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniHeader::FileExport,this,ID_EXPORT);
  /*
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniHeader::EditCut,this,wxID_CUT);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniHeader::EditCopy,this,wxID_COPY);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniHeader::EditPaste,this,wxID_PASTE);
  */
}

void MuniHeader::OnClose(wxCloseEvent& event)
{
  config->header_size = GetSize();
  Destroy();
}

void MuniHeader::FileClose(wxCommandEvent& event)
{
  Close();
}

void MuniHeader::SetHdu(const FitsHdu& h)
{
  wxASSERT(header);
  head = h;
  header->SetHdu(head);

  menuFile->Enable(ID_EXPORT,true);
  /*
  menuEdit->Enable(wxID_CUT,true);
  menuEdit->Enable(wxID_COPY,true);
  menuEdit->Enable(wxID_PASTE,true);
  */
}

void MuniHeader::FileExport(wxCommandEvent& event)
{
  wxFileDialog select(this,"Export As",wxEmptyString,wxEmptyString,
		      "TEXT files (*.txt)|*.txt",
		      wxFD_SAVE|wxFD_OVERWRITE_PROMPT|wxFD_CHANGE_DIR);

  if (select.ShowModal() == wxID_OK ) {

    wxBusyCursor wait;

    wxFileOutputStream output(select.GetPath());
    wxTextOutputStream cout(output);

    for(size_t i = 0; i < head.GetCount(); i++)
      cout << head.Item(i) << endl;

    // Create FITS export module !!!!!
  }
}

void MuniHeader::EditCut(wxCommandEvent& event)
{
  wxLogDebug("MuniHeader::EditCut not implemented");
}

void MuniHeader::EditCopy(wxCommandEvent& event)
{
  wxLogDebug("MuniHeader::EditCopy not implemented");
}

void MuniHeader::EditPaste(wxCommandEvent& event)
{
  wxLogDebug("MuniHeader::EditPaste not implemented");
}
