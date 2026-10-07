/*

  rawimport - import of raw files acquired by digital photographic cameras


  Copyright © 2010-2011, 2022-26 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/collpane.h>

using namespace std;


// ----- MuniProcessRaw ----

class MuniProcessRaw: public MuniProcess
{
private:

  wxWindow *win;
  wxString output;
  int idx;

  void OnPostProcess();

public:

  MuniProcessRaw(wxEvtHandler *h, wxWindow *w, const wxString& c, const wxArrayString& args):
    MuniProcess(h,c,args), win(w),idx(0) {}

  void SetOutput(const wxString& o) { output = o; }
  void SetIndex(int i) { idx = i; }

};

void MuniProcessRaw::OnPostProcess()
{
  wxASSERT(!output.IsEmpty());

  // check correctness of the file?

  wxASSERT(win);
  wxASSERT(dynamic_cast<MuniImportRaw *>(win));
  dynamic_cast<MuniImportRaw *>(win)->UpdateGauge(idx);
  dynamic_cast<MuniImportRaw *>(win)->LoadFile(output);
}



// ---- MuniImportRaw -----

MuniImportRaw::MuniImportRaw(wxWindow *w, const MuniImportRawOptions& opt,
			     const wxArrayString& files):
  wxDialog(w,wxID_ANY,"Raw Photos Import",wxDefaultPosition,wxDefaultSize,wxCAPTION),
  fcount(files.GetCount()),pipe(this)
{
  CreateControls();
  CreatePipe(opt,files);

  // set layout !! ???

  Bind(wxEVT_CLOSE_WINDOW,&MuniImportRaw::OnClose,this);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniImportRaw::OnCancel,this,wxID_CANCEL);
  Bind(wxEVT_END_PROCESS,&MuniImportRaw::OnFinish,this);

  pipe.Start();
  UpdateGauge(1);
  Layout();
}

MuniImportRaw::~MuniImportRaw()
{
  wxLogDebug("MuniImportRaw::~MuniImportRaw()");
}

void MuniImportRaw::CreateControls()
{
  // Controls
  label = new wxStaticText(this,wxID_ANY,"          Importing raw photos...             ");
  gauge = new wxGauge(this,wxID_ANY,fcount);
  gauge->SetValue(0);

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(label,wxSizerFlags().Align(wxALIGN_CENTER_HORIZONTAL).DoubleBorder());
  topsizer->Add(gauge,wxSizerFlags().Expand().Border());
  topsizer->Add(CreateButtonSizer(wxCANCEL),wxSizerFlags().Expand().Border());
  SetSizerAndFit(topsizer);
}

void MuniImportRaw::CreatePipe(const MuniImportRawOptions& opt, const wxArrayString& files)
{
  // Actions
  wxString type = opt.GetType();
  wxString band = opt.GetBand();
  wxString bitpix = opt.GetBitpix();
  wxString overw = opt.GetOverWrite();
  wxString dark = opt.GetDarkframe();
  wxString dcopt = opt.GetDcoptions();
  wxString interpol = opt.GetInterpol();
  wxString dir = opt.GetDirs();


  if( ! dark.IsEmpty() ) {
    dcdark = wxFileName::CreateTempFileName("dcdark");
    wxASSERT(! dcdark.IsEmpty());
    wxCopyFile(dark,dcdark);

    wxArrayString args;
    args.Add("-D"); args.Add("-4"); args.Add("-j"); args.Add("-t"); args.Add("0"); args.Add(dcdark);

    MuniProcess *com = new MuniProcess(&pipe,"dcraw",args);
    pipe.push(com);
  }

  for(size_t i = 0; i < files.GetCount(); i++) {


    wxFileName rname(files[i]);
    wxFileName fname(dir,rname.GetName(),"fits");

    wxArrayString args;

    args.Add("-t");
    args.Add(type);

    if( ! band.IsEmpty() ) {
      args.Add("-c");
      args.Add(band);
    }

    if( ! bitpix.IsEmpty() ) {
      args.Add("-b");
      args.Add(bitpix);
    }

    wxString opt = dcopt;

    if( ! dcdark.IsEmpty() )
      opt += " -K " + dcdark + ".pgm";

    if( ! interpol.IsEmpty() )
      opt += " " + interpol;

    if( ! opt.IsEmpty() ) {
      args.Add("-X");
      args.Add(opt);
    }

    args.Add("-o");
    args.Add(overw + fname.GetFullPath());
    args.Add(files[i]);

    MuniProcessRaw *com = new MuniProcessRaw(&pipe,this,"rawtran",args);
    com->SetOutput(fname.GetFullPath());
    com->SetIndex(i+1);
    pipe.push(com);
  }
}

void MuniImportRaw::OnFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniImportRaw::OnFinish");

  // remove dcdark
  if( ! dcdark.IsEmpty() ) {
    wxRemoveFile(dcdark);
    wxRemoveFile(dcdark+".pgm");
  }

  Destroy();
}

void MuniImportRaw::OnClose(wxCloseEvent& event)
{
  // One doesn't receive the proper signal ?
  wxLogDebug("MuniImportRaw::OnClose");
  pipe.Stop();
}

void MuniImportRaw::OnCancel(wxCommandEvent& event)
{
  wxLogDebug("MuniImportRaw::OnCancel");
  pipe.Stop();
}

void MuniImportRaw::UpdateGauge(int value)
{
  const wxString blank("          ");
  wxString a;
  a.Printf(blank+"Importing %d of %d files."+blank,value,fcount);

  label->SetLabel(a);
  gauge->SetValue(value);
}

void MuniImportRaw::LoadFile(const wxString& file)
{
  wxArrayString a;
  a.Add(file);
  dynamic_cast<MuniBrowser *>(GetParent())->FileLoad(a);
}


// --- MuniImportRawOptions ----

MuniImportRawOptions::MuniImportRawOptions(wxWindow *w, MuniConfig *c):
  wxDialog(w,wxID_ANY,"Raw Pictures Import"),config(c)
{
  SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);
  SetIcon(config->munipack_icon);
  EnableCloseButton(false);

  Init();
  CreateControls();
}

void MuniImportRawOptions::Init()
{
  fchoices.Add("Red Cone Eye");          fopt.Add("X");
  fchoices.Add("Green Cone Eye");        fopt.Add("Y");
  fchoices.Add("Blue Cone Eye");         fopt.Add("Z");
  fchoices.Add("Rod Eye");               fopt.Add("s");
  fchoices.Add("Unfiltered");            fopt.Add("u");
  fchoices.Add("Landolt R");             fopt.Add("R");
  fchoices.Add("Landolt V");             fopt.Add("V");
  fchoices.Add("Landolt B");             fopt.Add("B");
  fchoices.Add("Instrumental R");        fopt.Add("Ri");
  fchoices.Add("Instrumental G");        fopt.Add("Gi");
  fchoices.Add("Instrumental B");        fopt.Add("Bi");
  fchoices.Add("Instrumental G1");       fopt.Add("Gi1");
  fchoices.Add("Instrumental G2");       fopt.Add("Gi2");

  ichoices.Add("Adaptive Homogeneity-Directed (AHD)"); iopt.Add("-q 3");
  ichoices.Add("Patterned Pixel Grouping (PPG)");      iopt.Add("-q 2");
  ichoices.Add("Variable Number of Gradients (VNG)");  iopt.Add("-q 1");
  ichoices.Add("Bilinear");                            iopt.Add("-q 0");
  ichoices.Add("Interpolate RGB as four colors");      iopt.Add("-f");
  ichoices.Add("Output a half-size color image");      iopt.Add("-h");

  type_colour = true;
  type_grey = false;
  filter = 1;
  bitpix_16bit = true;
  bitpix_float = false;
  over = true;
  dark = wxEmptyString;
  dcopt = wxEmptyString;
  interpol = 0;
}

void MuniImportRawOptions::CreateControls()
{
  wxFont bf(*wxNORMAL_FONT);
  bf.SetWeight(wxFONTWEIGHT_BOLD);

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);

  wxBoxSizer *ysizer = new wxBoxSizer(wxVERTICAL);

  wxBoxSizer *osizer = new wxBoxSizer(wxHORIZONTAL);

  wxImage import_raw_icon = MuniConfig::LoadImage("gphoto.png");
  MuniThumbCanvas *th = new MuniThumbCanvas(this,import_raw_icon);
  osizer->Add(th,wxSizerFlags().Center().DoubleBorder());

  wxStaticText *title = new wxStaticText(this,wxID_ANY,"Digital Camera RAW Pictures");
  title->SetFont(bf);
  osizer->Add(title,wxSizerFlags().DoubleBorder().Align(wxALIGN_CENTER_VERTICAL));
  ysizer->Add(osizer,wxSizerFlags().Center());

  wxFlexGridSizer *csizer = new wxFlexGridSizer(2);
  csizer->AddGrowableCol(1);
  csizer->Add(new wxStaticText(this,wxID_ANY,"Conversion:"),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));

  wxBoxSizer *tsizer = new wxBoxSizer(wxVERTICAL);

  type0 = new wxRadioButton(this,wxID_ANY,"Colour",wxDefaultPosition,
			    wxDefaultSize,wxRB_GROUP);
  type1 = new wxRadioButton(this,wxID_ANY,"Grey");
  tsizer->Add(type0,wxSizerFlags());
  tsizer->Add(type1,wxSizerFlags());

  filters = new wxChoice(this,wxID_ANY,wxDefaultPosition,wxDefaultSize,
			 fchoices);
  filters->SetToolTip("Select the output spectral band.");
  wxBoxSizer *fbox = new wxBoxSizer(wxHORIZONTAL);
  fbox->Add(new wxStaticText(this,wxID_ANY,"Filters:"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).Border(wxLEFT));
  fbox->Add(filters,wxSizerFlags());
  tsizer->Add(fbox,wxSizerFlags().Center());
  csizer->Add(tsizer,wxSizerFlags().Border());

  csizer->Add(new wxStaticText(this,wxID_ANY,"Output:"),
	      wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  dirs = new wxDirPickerCtrl(this,wxID_ANY);
  dirs->SetToolTip("Specify a directory where the output files will be stored.");
  csizer->Add(dirs,wxSizerFlags().Expand().Border());

  ysizer->Add(csizer,wxSizerFlags().Right());

  topsizer->Add(ysizer,wxSizerFlags().Border());

  wxCollapsiblePane *collPane = new wxCollapsiblePane(this,wxID_ANY,"Details");
  wxWindow *win = collPane->GetPane();

  // detail panel
  wxFlexGridSizer *gsizer = new wxFlexGridSizer(2);
  gsizer->AddGrowableCol(1);

  // output file
  gsizer->Add(new wxStaticText(win,wxID_ANY,"Output FITS:"),
	      wxSizerFlags().Border().Align(wxALIGN_RIGHT));

  wxBoxSizer *psizer = new wxBoxSizer(wxVERTICAL);

  bitpix0 = new wxRadioButton(win,wxID_ANY,"16-bit",wxDefaultPosition,
			    wxDefaultSize,wxRB_GROUP);
  bitpix1 = new wxRadioButton(win,wxID_ANY,"float");
  wxString tooltip("Select representation of numbers in output image. The float numbers are intended for general usage and additional processing. 16-bit numbers are less-precise but saves some space.");
  bitpix0->SetToolTip(tooltip);
  bitpix1->SetToolTip(tooltip);
  psizer->Add(bitpix0,wxSizerFlags());
  psizer->Add(bitpix1,wxSizerFlags());

  overs = new wxCheckBox(win,wxID_ANY,"Overwrite file(s)");
  overs->SetToolTip("Check it, if you want overwrite output files.");
  psizer->Add(overs,wxSizerFlags().Right());
  gsizer->Add(psizer);

  // interpolation
  gsizer->Add(new wxStaticText(win,wxID_ANY,"Interpolation:"),
	      wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  interpols = new wxChoice(win,wxID_ANY,wxDefaultPosition,wxDefaultSize,
			 ichoices);
  interpols->SetToolTip("Select method for color interpolation.");
  gsizer->Add(interpols,wxSizerFlags().Border());

  // dark frame
  gsizer->Add(new wxStaticText(win,wxID_ANY,"Dark Frame:"),
	      wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  darks = new wxFilePickerCtrl(win,wxID_ANY,wxEmptyString,
			       "Select a Dark Frame",
			       "RAW files ("+config->rawmask+")|"+
			       config->rawmask+"|All files (*)|*");
  darks->SetToolTip("Specify file to be subtracted from exposures. The file is usually grabbed as a dark frame exposure (shot with no light) at same ambient temperature.");
  gsizer->Add(darks,wxSizerFlags().Expand().Border());

  // options
  gsizer->Add(new wxStaticText(win,wxID_ANY,"Options for dcraw:"),
	      wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  dcopts = new wxTextCtrl(win,wxID_ANY);
  dcopts->SetToolTip("Write any additional parameters. It'll passed directly to dcraw.");
  gsizer->Add(dcopts,wxSizerFlags().Expand().Border());

  win->SetSizer(gsizer);
  gsizer->SetSizeHints(win);
  topsizer->Add(collPane,wxSizerFlags().Border(wxLEFT|wxRIGHT|wxBOTTOM));

  topsizer->Add(CreateButtonSizer(wxOK|wxCANCEL),
		wxSizerFlags().Expand().Border());
  SetSizerAndFit(topsizer);


  // data
  type0->SetValidator(wxGenericValidator(&type_colour));
  type1->SetValidator(wxGenericValidator(&type_grey));
  bitpix0->SetValidator(wxGenericValidator(&bitpix_16bit));
  bitpix1->SetValidator(wxGenericValidator(&bitpix_float));
  filters->SetValidator(wxGenericValidator(&filter));
  overs->SetValidator(wxGenericValidator(&over));
  interpols->SetValidator(wxGenericValidator(&interpol));
  dcopts->SetValidator(wxTextValidator(wxFILTER_NONE,&dcopt));
  //  darks->SetValidator(wxGenericValidator(&dark));
  // NOTE: Validator for FilePicker DOES NOT work !


  // Update UI
  Bind(wxEVT_UPDATE_UI,&MuniImportRawOptions::OnUpdateUI,this,filters->GetId());
  Bind(wxEVT_COMMAND_FILEPICKER_CHANGED,&MuniImportRawOptions::OnDarks,this,darks->GetId());
  Bind(wxEVT_COMMAND_DIRPICKER_CHANGED,&MuniImportRawOptions::OnDirs,this,dirs->GetId());
}

void MuniImportRawOptions::OnDarks(wxFileDirPickerEvent& event)
{
  dark = event.GetPath();
}

void MuniImportRawOptions::OnDirs(wxFileDirPickerEvent& event)
{
  dir = event.GetPath();
}

void MuniImportRawOptions::OnUpdateUI(wxUpdateUIEvent& event)
{
  event.Enable(type1->GetValue());
}

wxString MuniImportRawOptions::GetType() const
{
  if( type_colour )
    return "0";

  else {

    if( fchoices[filter].StartsWith("Instrumental") )
      return "1";
    else
      return "0";

  }
}

wxString MuniImportRawOptions::GetBand() const
{
  if( type_grey )
    return fopt[filter];
  else
    return wxEmptyString;
}

wxString MuniImportRawOptions::GetBitpix() const
{
  if( bitpix_16bit )
    return "16";

  if( bitpix_float )
    return "-32";

  return wxEmptyString;
}

wxString MuniImportRawOptions::GetOverWrite() const
{
  if( over )
    return "!";
  else
    return wxEmptyString;
}

wxString MuniImportRawOptions::GetInterpol() const
{
  return iopt[interpol];
}

wxString MuniImportRawOptions::GetDarkframe() const
{
  return dark;
}

wxString MuniImportRawOptions::GetDcoptions() const
{
  return dcopt;
}

wxString MuniImportRawOptions::GetDirs() const
{
  return dir;
}
