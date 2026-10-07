/*

  xmunipack - browser search entry


  Copyright © 1997-2011 F.Hroch (hroch@physics.muni.cz)

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
#include <algorithm>
#include <vector>
#include <wx/wx.h>
#include <wx/srchctrl.h>
#include <wx/regex.h>
#include <wx/tokenzr.h>


using namespace std;


MuniBrowserSearch::MuniBrowserSearch(wxWindow *w, wxWindowID id, 
				     const wxString& value, 
				     const wxPoint& pos, const wxSize& size, 
				     long style):
  wxSearchCtrl(w,id,value,pos,size,style)
{
  type = ID_SELECT_BY_KEY;

  wxMenu *smenu = new wxMenu;
  smenu->AppendRadioItem(ID_SELECT_BY_KEY,"Search in Metadata");
  smenu->AppendRadioItem(ID_SELECT_BY_NAME,"Search by Filename");
  smenu->AppendRadioItem(ID_SELECT_BY_ADVANCE,"Advanced Search");
  SetMenu(smenu);

  Bind(wxEVT_COMMAND_TEXT_ENTER,&MuniBrowserSearch::OnSearchEnter,this);
  Bind(wxEVT_COMMAND_SEARCHCTRL_SEARCH_BTN,&MuniBrowserSearch::OnSearchButton,this);
  Bind(wxEVT_COMMAND_SEARCHCTRL_CANCEL_BTN,&MuniBrowserSearch::OnSearchFinish,this);
  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniBrowserSearch::OnSearchUpdate,this);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniBrowserSearch::OnSearchMenu,this);
  Bind(wxEVT_UPDATE_UI,&MuniBrowserSearch::OnUpdateUI,this);
}

void MuniBrowserSearch::OnSearchUpdate(wxCommandEvent& event)
{
  ShowCancelButton(true);
}

void MuniBrowserSearch::OnUpdateUI(wxUpdateUIEvent& event)
{
  switch(type) {
  case ID_SELECT_BY_KEY:
    SetToolTip("Type a string to search for.");
    break;
  case ID_SELECT_BY_NAME:
    SetToolTip("Type a file name to search for.");
    break;
  case ID_SELECT_BY_ADVANCE:
    SetToolTip("Type a regular expression to search for. Use key: to specify a keyword. Example: M.*31 key:object"); 
    break;    
  }
}

void MuniBrowserSearch::OnSearchMenu(wxCommandEvent& event)
{
  type = event.GetId();
}

void MuniBrowserSearch::OnSearchFinish(wxCommandEvent& event)
{
  SetValue(wxEmptyString);
  ShowCancelButton(false);
  wxQueueEvent(GetParent(),event.Clone());
}

void MuniBrowserSearch::OnSearchButton(wxCommandEvent& event)
{
  // this empty code blocks propagation to top win (browser)
}

void MuniBrowserSearch::OnSearchEnter(wxCommandEvent& event)
{
  muster = event.GetString();
  wxLogDebug(muster);

  wxCommandEvent e(wxEVT_COMMAND_SEARCHCTRL_SEARCH_BTN,wxID_FIND);
  wxQueueEvent(GetParent(),e.Clone());
}

vector<long> MuniBrowserSearch::Find(const std::vector<FitsMeta>& list) const
{
  if( type == ID_SELECT_BY_NAME )
    return FindByName(list);
  else if( type == ID_SELECT_BY_KEY )
    return FindByKey(list);
  else if( type == ID_SELECT_BY_ADVANCE )
    return FindByAdv(list);
  else {
    wxFAIL_MSG("REACHED !!! ");
    vector<long> index;
    return index;
  }
}

vector<long> MuniBrowserSearch::FindByName(const vector<FitsMeta>& list) const
{
  vector<long> index;
  for(vector<FitsMeta>::const_iterator i=list.begin();i!=list.end();++i) {
    wxString name = i->GetName();
    if( name.Find(muster) != wxNOT_FOUND )
      index.push_back(i - list.begin());
  }
  return index;
}

vector<long> MuniBrowserSearch::FindByKey(const vector<FitsMeta>& list) const
{
  vector<long> index;
  for(vector<FitsMeta>::const_iterator m = list.begin(); m != list.end(); ++m ){

    for(size_t i = 0; i < m->HduCount(); i++) {
      FitsHeader h(m->Hdu(i));
      for(size_t l = 0; l < h.GetCount(); l++) {
	if( h[l].Find(muster) != wxNOT_FOUND ) {
	  index.push_back(m - list.begin());
	  goto xxx;
	}
      }
    }
  xxx: ;
  }
  return index;
}


vector<long> MuniBrowserSearch::FindByAdv(const vector<FitsMeta>& list) const
{
  //  parse mask
  wxString key,val = muster;
  wxStringTokenizer tkz(muster,", ");
  while ( tkz.HasMoreTokens() ) {
    wxString token = tkz.GetNextToken();

    if( token.Find("key:") ) {
      wxStringTokenizer tkz(token,":");
      if( tkz.HasMoreTokens() ) {
	tkz.GetNextToken();
	if( tkz.HasMoreTokens() )
	  key = tkz.GetNextToken();
      }
    }
    else 
      val = token;

  }
  wxLogDebug(key+val);

  vector<long> index;
  for(vector<FitsMeta>::const_iterator m = list.begin(); m != list.end(); ++m ){

    for(size_t i = 0; i < m->HduCount(); i++) {

      wxRegEx re(val,wxRE_EXTENDED|wxRE_ICASE);

      if( ! key.IsEmpty() ) {
	FitsHeader h(m->Hdu(i));
	wxString c = h.GetKey(key);

	if( re.Matches(c) ) {
	  index.push_back(m - list.begin());
	  goto xxx;
	}
      }
      else {

	for(size_t i = 0; i < m->HduCount(); i++) {
	  FitsHeader h(m->Hdu(i));
	  for(size_t l = 0; l < h.GetCount(); l++) {
	    if( re.Matches(h[l]) ) {
	      index.push_back(m - list.begin());
	      goto xxx;
	    }
	  }
	}
      }
    }
  xxx: ;
  }
  return index;
}


