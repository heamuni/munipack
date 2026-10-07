/*

  xmunipack - average of darks, biases and flats


  Copyright © 2010-2012 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/valgen.h>
#include <vector>

using namespace std;


// --- Options

class MuniAverageOptions: public wxDialog
{
public:
  MuniAverageOptions(wxWindow *, MuniConfig *);
  wxString GetBitpix() const;
  bool GetRobust() const;
  wxString GetLevel() const;

private:

  MuniConfig *config;
  wxRadioButton *bitpix0, *bitpix1;
  wxCheckBox *rbutt;
  wxTextCtrl *xlevel;
  bool robust,bitpix_16bit,bitpix_float;
  wxString level;

  void Init();
  void CreateControls();
  void OnUpdateUI(wxUpdateUIEvent&);

};

MuniAverageOptions::MuniAverageOptions(wxWindow *w, MuniConfig *c):
  wxDialog(w,wxID_ANY,"Averaging Options"),config(c)
{
  SetIcon(config->munipack_icon);
  EnableCloseButton(false);

  Init();
  CreateControls();
}

void MuniAverageOptions::Init()
{
  bitpix_16bit = false;
  bitpix_float = true;
  robust = true;
  level = "1";
}


void MuniAverageOptions::CreateControls()
{
  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);

  wxFlexGridSizer *gsizer = new wxFlexGridSizer(2);
  gsizer->AddGrowableCol(1);

  // output file
  gsizer->Add(new wxStaticText(this,wxID_ANY,"Output FITS:"),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));

  wxBoxSizer *psizer = new wxBoxSizer(wxVERTICAL);

  bitpix0 = new wxRadioButton(this,wxID_ANY,"16-bit",wxDefaultPosition,
			    wxDefaultSize,wxRB_GROUP);
  bitpix1 = new wxRadioButton(this,wxID_ANY,"float");
  bitpix0->SetToolTip("Select representation of numbers in output image. The float numbers are intended for general usage. 16-bit numbers are less-precise but saves some space.");
  bitpix1->SetToolTip("Select representation of numbers in output image. The float numbers are intended for general usage. 16-bit numbers are less-precise but saves some space.");
  psizer->Add(bitpix0,wxSizerFlags());
  psizer->Add(bitpix1,wxSizerFlags());
  gsizer->Add(psizer,wxSizerFlags().Border());

  xlevel = new wxTextCtrl(this,wxID_ANY,"1");
  xlevel->SetToolTip("Sets the mean level of the output frame (valid only for flat-field).");
  gsizer->Add(new wxStaticText(this,wxID_ANY,"Level:"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  gsizer->Add(xlevel,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).Expand().Border());

  rbutt = new wxCheckBox(this,wxID_ANY,"Be Robust");
  rbutt->SetToolTip("Selects algorithm to estimate of the mean. Check for the robust mean (recommended). Unchecked for the arithmetical mean (faster). One is valid only for darks and biases.");
  gsizer->Add(new wxStaticText(this,wxID_ANY,"Algorithm:"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  gsizer->Add(rbutt,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL).Border());
  topsizer->Add(gsizer,wxSizerFlags().Border().Expand());

  topsizer->Add(CreateButtonSizer(wxOK|wxCANCEL),wxSizerFlags().Expand().Border());
  SetSizerAndFit(topsizer);

  // data
  bitpix0->SetValidator(wxGenericValidator(&bitpix_16bit));
  bitpix1->SetValidator(wxGenericValidator(&bitpix_float));
  rbutt->SetValidator(wxGenericValidator(&robust));
  xlevel->SetValidator(wxTextValidator(wxFILTER_NONE,&level));

  Connect(xlevel->GetId(),wxEVT_UPDATE_UI,
	  wxUpdateUIEventHandler(MuniAverageOptions::OnUpdateUI));
}

void MuniAverageOptions::OnUpdateUI(wxUpdateUIEvent& event)
{
  if( bitpix0->GetValue() )
    event.SetText("10000");

  else if( bitpix1->GetValue() )
    event.SetText("1");
}


wxString MuniAverageOptions::GetBitpix() const
{
  if( bitpix_16bit )
    return "16";

  if( bitpix_float )
    return "-32";

  return wxEmptyString;
}

bool MuniAverageOptions::GetRobust() const
{
  return robust;
}

wxString MuniAverageOptions::GetLevel() const
{
  return level;
}


// --  Average


MuniAverage::MuniAverage(wxWindow *w, wxWindowID id, long t, MuniConfig *c):
  MuniListWindow(w,id,t,c), config(c), pipe(this), timer(this)
{
  wxWindow *window = new wxPanel(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
				 wxTAB_TRAVERSAL|wxBORDER_THEME);

  MuniArtIcons icons(wxART_CMN_DIALOG,wxSize(48,48));

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);

  wxBoxSizer *isizer = new wxBoxSizer(wxHORIZONTAL);


  wxFlexGridSizer *gsizer = new wxFlexGridSizer(2);
  gsizer->AddGrowableCol(1);

  gsizer->Add(new wxStaticText(window,wxID_ANY,"Average:"),
	      wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  btype[0] = new wxRadioButton(window,wxID_ANY,"Bias",
					    wxDefaultPosition,wxDefaultSize,
					    wxRB_GROUP);
  btype[1] = new wxRadioButton(window,wxID_ANY,"Dark");
  btype[2] = new wxRadioButton(window,wxID_ANY,"Flat");

  wxBoxSizer *btypes = new wxBoxSizer(wxHORIZONTAL);
  for(int i = 0; i < 3; i++)
    btypes->Add(btype[i]);

  gsizer->Add(btypes,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL));

  gsizer->Add(new wxStaticText(window,wxID_ANY,"Filename:"),
	      wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  flabel = new wxTextCtrl(window,wxID_ANY);
  gsizer->Add(flabel,wxSizerFlags(1).Expand().Border());


  gsizer->Add(new wxStaticText(window,wxID_ANY,"Destination:"),
	      wxSizerFlags().Border(wxLEFT).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));

  dirpic = new wxDirPickerCtrl(window,wxID_ANY);
  gsizer->Add(dirpic,wxSizerFlags(1).Expand().Border(wxLEFT|wxRIGHT));

  isizer->Add(gsizer,wxSizerFlags(1).Expand());

  MuniThumbCanvas *th = new MuniThumbCanvas(window,icons.Icon(wxART_EXECUTABLE_FILE));
  isizer->Add(th,wxSizerFlags().Center().TripleBorder());

  topsizer->Add(isizer,wxSizerFlags().Expand().Border());

  wxButton *bcre = new wxButton(window,ID_MDARK_CREATE,"Create");
  wxButton *bdet = new wxButton(window,wxID_ANY,"Options...");
  gauge = new wxGauge(window,wxID_ANY,1);
  label = new wxStaticText(window,wxID_ANY,"Label");

  wxSizer *butt = new wxBoxSizer(wxHORIZONTAL);
  butt->Add(bcre,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL));
  butt->Add(bdet,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL));
  butt->Add(gauge,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL));
  butt->Add(label,wxSizerFlags(1).Border().Align(wxALIGN_CENTER_VERTICAL));
  topsizer->Add(butt,wxSizerFlags().Border().Expand());

  window->SetSizer(topsizer);

  AddWindow(window);


  // data
  filename[0] = "bias.fits";
  filename[1] = "dark.fits";
  filename[2] = "flat.fits";
  xtype = 1;
  robust = true;
  level = "1";
  btype[1]->SetValue(true);
  bcre->Enable(false);
  flabel->SetValue(filename[xtype]);
  gauge->Show(false);
  label->Show(false);

  Bind(wxEVT_COMMAND_RADIOBUTTON_SELECTED,&MuniAverage::OnBtype,this,btype[0]->GetId());
  Bind(wxEVT_COMMAND_RADIOBUTTON_SELECTED,&MuniAverage::OnBtype,this,btype[1]->GetId());
  Bind(wxEVT_COMMAND_RADIOBUTTON_SELECTED,&MuniAverage::OnBtype,this,btype[2]->GetId());
  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniAverage::OnFlabel,this,flabel->GetId());
  Bind(wxEVT_UPDATE_UI,&MuniAverage::OnUpdateButt,this,bcre->GetId());
  Bind(wxEVT_COMMAND_DIRPICKER_CHANGED,&MuniAverage::OnDirname,this,dirpic->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAverage::OnOptions,this,bdet->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAverage::OnCreate,this,bcre->GetId());
}


void MuniAverage::OnBtype(wxCommandEvent& event)
{
  for(int i = 0; i < 3; i++)
    if( event.GetId() == btype[i]->GetId() )
      xtype = i;

  flabel->SetValue(filename[xtype]);
}

void MuniAverage::OnFlabel(wxCommandEvent& event)
{
  filename[xtype] = event.GetString();
}

void MuniAverage::OnUpdateButt(wxUpdateUIEvent& event)
{
  event.Enable(GetItemCount() > 0);
}

void MuniAverage::OnDirname(wxFileDirPickerEvent& event)
{
  dirname = event.GetPath();
  wxLogDebug(dirname);
}

void MuniAverage::OnOptions(wxCommandEvent& event)
{
  MuniAverageOptions opt(this,config);
  if( opt.ShowModal() == wxID_OK ) {
    bitpix = opt.GetBitpix();
    robust = opt.GetRobust();
    level = opt.GetLevel();
  }
}


void MuniAverage::OnCreate(wxCommandEvent& event)
{
  wxString command;

  switch(xtype) {
  case 0: command = "adark"; break;
  case 1: command = "adark"; break;
  case 2: command = "aflat"; break;
  }

  MuniProcess *com = new MuniProcess(&pipe,command);

  wxFileName f(dirname,filename[xtype]);
  wxString rob = robust ? "'R'" : "'A'";

  com->Write("OUTPUT = '"+f.GetFullPath()+"'");

  if( xtype == 2 ) {
    if( ! level.IsEmpty() )
      com->Write("Output level = " + level);
    com->Write("Grid step = 1");
  }
  else
    com->Write("MEAN = " + rob);

  if( ! bitpix.IsEmpty() ) {
    if( bitpix == "float" )
      com->Write("BITPIX = -32");
    else if( bitpix == "16-bit" )
      com->Write("BITPIX = 16");
    else
      com->Write("BITPIX = -32");
  }
  else
    com->Write("BITPIX = -32");

  vector<FitsMeta>::const_iterator m;
  vector<FitsMeta> metalist = GetAllMeta();
  for(m = metalist.begin(); m != metalist.end(); ++m){
    wxASSERT(m->IsOk());
    //    com->WriteFiles(m->GetURL());
    com->Write("FILE = '"+m->GetFullPath()+"'");
    wxLogDebug(m->GetURL());
    //    com->Write("'"+m->GetURL()+"'");
  }

  pipe.push(com);

  timer.Start(200);
  gauge->SetRange(metalist.size());
  gauge->Show(true);
  label->Show(true);
  label->SetLabel("Creating average ...");

  Bind(wxEVT_END_PROCESS,&MuniAverage::OnFinish,this);
  Bind(wxEVT_TIMER,&MuniAverage::OnUpdate,this);

  pipe.Start();

  Layout();
}

void MuniAverage::OnFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniAverage::OnFinish");

  timer.Stop();
  gauge->Show(false);
  label->Show(false);

  Unbind(wxEVT_END_PROCESS,&MuniAverage::OnFinish,this);
  Unbind(wxEVT_TIMER,&MuniAverage::OnUpdate,this);

  if( event.GetExitCode() == 0 ) {

    wxFileName f(dirname,filename[xtype]);
    FitsFile fits(f.GetFullPath());
    if( fits.Status() ) {
      /*
      MuniIcon micon(fits,config->display_colourspace,config->cdatafile,
		     config->default_icon,config->table_icon,config->head_icon,
		     config->icon_size);
      FitsMeta meta(fits,micon.GetIcon(),micon.GetList());
      if( meta.IsOk() )
	AddMeta(meta);
      */
    }
  }

  Layout();
}

void MuniAverage::OnUpdate(wxTimerEvent& event)
{
  wxLogDebug("MuniAverage::OnUpdate ");
  wxArrayString out(pipe.GetOutput());

  int n = out.GetCount();
  if( n <= 0 || n > gauge->GetRange() )
    gauge->Pulse();
  else
    gauge->SetValue(n);

  // if( out.GetCount() > 0 ) {
  //   wxLogDebug("MuniAverage::OnUpdate "+out.Last());
  //   if( out.GetCount() < gauge->GetRange() )
  //     gauge->SetValue(out.GetCount());
  //   else {
  //     // gauge->SetValue(gauge->GetRange());
  //     gauge->Pulse();
  //   }
  // }

}
