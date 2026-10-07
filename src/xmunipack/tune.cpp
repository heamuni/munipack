/*

  xmunipack - tune panel

  Copyright © 2021-2022 F.Hroch (hroch@physics.muni.cz)

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

#include "tune.h"
#include <wx/wx.h>
#include <vector>
#include <limits>
#include <cfloat>

#define TUNE_FRAME_STYLE wxCAPTION | wxCLOSE_BOX | wxFRAME_TOOL_WINDOW

#define TOOL_TITLE  "Tune display"
#define HIST_WIDTH  600
#define HIST_HEIGHT 370
#define SCALE_FMT "%.3g"

// ---- MuniTune

MuniTune::MuniTune(wxWindow *w, wxWindowID id, const wxPoint& pos,
		   const wxSize& size, int icon_size,
		   const FitsArray& a, const FitsTone& t,
		   const FitsItt& i, const FitsPalette& p):
  wxFrame(w,id,TOOL_TITLE,pos,size,TUNE_FRAME_STYLE),
  array(a),black(t.GetBlack()),sense(t.GetSense()),
  itt(i.GetItt()),palette(p.GetPalette()),inverse(p.GetInverse()),qrange(0.999),
  undo(false), nitevision(false),pltpal(0),pltsat(0),pltnite(0),
  thread_histo(NULL),shutdown(false)
{
  wxNotebook *book = new wxNotebook(this,wxID_ANY);

  // prescale
  wxPanel *spanel = new wxPanel(book);
  FeedScaleTab(spanel);
  book->AddPage(spanel,"Scale");

  // ITT
  wxPanel *ipanel = new wxPanel(book);
  FeedIttTab(ipanel);
  book->AddPage(ipanel,"Tone");

  // Palette
  wxPanel *ppanel = new wxPanel(book);
  FeedPaletteTab(ppanel);
  book->AddPage(ppanel,"Palette");

  // arrange all
  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(book,wxSizerFlags().Center());
  SetSizerAndFit(topsizer);

  // events
  Bind(wxEVT_CLOSE_WINDOW,&MuniTune::OnClose,this);
  Bind(wxEVT_CHAR_HOOK,&MuniTune::OnEscape,this);

  RenderHisto();
}

MuniTune::MuniTune(wxWindow *w, wxWindowID id, const wxPoint& pos,
		   const wxSize& size, int icon_size, const FitsArray& a,
		   const FitsTone& t, const FitsItt& i,
		   const FitsColour& c):
  wxFrame(w,id,TOOL_TITLE,pos,size,TUNE_FRAME_STYLE),
  array(a), black(t.GetBlack()),sense(t.GetSense()),
  itt(i.GetItt()),palette(0),saturation(c.GetSaturation()),
  x_white(c.GetWhitePointX()),y_white(c.GetWhitePointY()),
  meso_level(c.GetMesoLevel()),meso_width(c.GetMesoWidth()),qrange(0.999),
  undo(false),nitevision(c.GetNiteVision()),pltpal(0),pltsat(0),pltnite(0),
  thread_histo(NULL),shutdown(false)
{
  wxNotebook *book = new wxNotebook(this,wxID_ANY);

  // prescale
  wxPanel *spanel = new wxPanel(book);
  FeedScaleTab(spanel);
  book->AddPage(spanel,"Scale");

  // ITT
  wxPanel *ipanel = new wxPanel(book);
  FeedIttTab(ipanel);
  book->AddPage(ipanel,"Tone");

  // Colour
  wxPanel *cpanel = new wxPanel(book);
  FeedColourTab(cpanel);
  book->AddPage(cpanel,"Colours");

  // Night vision
  wxPanel *npanel = new wxPanel(book);
  FeedNiteTab(npanel);
  book->AddPage(npanel,"Nite");

  // arrange all
  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(book,wxSizerFlags().Center());
  SetSizerAndFit(topsizer);

  Bind(wxEVT_CLOSE_WINDOW,&MuniTune::OnClose,this);
  Bind(wxEVT_CHAR_HOOK,&MuniTune::OnEscape,this);

  RenderHisto();
}

MuniTune::~MuniTune()
{
  //  wxLogDebug("MuniTune::~MuniTune()");
  if( thread_histo ) {
    thread_histo->Delete();
  }
}


void MuniTune::OnClose(wxCloseEvent& event)
{
  if( thread_histo ) {
    thread_histo->Delete();
    event.Veto();
    shutdown = true;
    return;
  }

  //  wxLogDebug("MuniTune::OnClose()");
  wxQueueEvent(GetParent(),event.Clone());
}

void MuniTune::OnEscape(wxKeyEvent& event)
{
  if( event.GetKeyCode() == WXK_ESCAPE )
    Close();
  else
    event.Skip();
}

void MuniTune::FeedScaleTab(wxPanel *tpanel)
{
  wxSizerFlags vsizer(1);
  vsizer.Align(wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT).Expand();

  int w = 12*GetCharWidth();

  entry_sense = new wxTextCtrl(tpanel,ID_TONE_SENSE,ToString(SCALE_FMT,sense),
			       wxDefaultPosition,wxSize(w,-1),wxTE_PROCESS_ENTER);
  entry_black = new wxTextCtrl(tpanel,ID_TONE_BLACK,ToString(SCALE_FMT,black),
			       wxDefaultPosition,wxSize(w,-1),wxTE_PROCESS_ENTER);
  entry_sense->SetToolTip("Sensitivity adjusts the displayed range of intensities.");
  entry_black->SetToolTip("Black point is an intensity corresponding "
			  "to black colour on the display.");

  wxSizerFlags label_flags;
  label_flags.Align(wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT);
  wxSizerFlags icon_flags;
  icon_flags.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER);
  wxSizerFlags entry_flags;
  entry_flags.Align(wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT|wxTOP);
  wxSizerFlags button_flags;
  button_flags.Align(wxALIGN_CENTER_VERTICAL);

  wxFlexGridSizer *grid = new wxFlexGridSizer(4);
  grid->Add(new wxStaticText(tpanel,wxID_ANY,"Sensitivity:"),label_flags);
  grid->Add(new wxStaticText(tpanel,wxID_ANY,L"☾"),icon_flags);
  grid->Add(entry_sense,entry_flags);
  grid->Add(new wxStaticText(tpanel,wxID_ANY,L"☼"),icon_flags);
  grid->Add(new wxStaticText(tpanel,wxID_ANY,"Black:"),label_flags);
  grid->Add(new wxStaticText(tpanel,wxID_ANY,L"◑"),icon_flags);
  grid->Add(entry_black,entry_flags);
  grid->Add(new wxStaticText(tpanel,wxID_ANY,L"○"),icon_flags);

  wxButton *apply = new wxButton(tpanel,wxID_APPLY);
  wxButton *bundo = new wxButton(tpanel,wxID_UNDO,L"⎌",wxDefaultPosition,
				 wxDefaultSize,wxBU_EXACTFIT);
  bundo->SetToolTip("Undo to the initial setup");

  wxBoxSizer *ugrid = new wxBoxSizer(wxHORIZONTAL);
  ugrid->AddStretchSpacer(1);
  ugrid->Add(grid);
  ugrid->Add(apply,button_flags.DoubleBorder());
  ugrid->Add(bundo,button_flags.Border());

  // histogram
  plthist = new MuniPlotHisto(tpanel,wxSize(HIST_WIDTH,HIST_HEIGHT),black,sense);

  // pixel ranges
  rb99 = new wxRadioButton(tpanel,wxID_ANY,"99%",wxDefaultPosition,wxDefaultSize,
			  wxRB_GROUP);
  rb999 = new wxRadioButton(tpanel,wxID_ANY,"99.9%");
  rbfull = new wxRadioButton(tpanel,wxID_ANY,"Full");
  rb999->SetValue(true);

  wxBoxSizer *rbsizer = new wxBoxSizer(wxHORIZONTAL);
  rbsizer->Add(new wxStaticText(tpanel,wxID_ANY,"Quantile range:"));
  rbsizer->Add(rb99,wxSizerFlags().Border(wxLEFT));
  rbsizer->Add(rb999,wxSizerFlags().Border(wxLEFT));
  rbsizer->Add(rbfull,wxSizerFlags().Border(wxLEFT));

  wxBoxSizer *isizer = new wxBoxSizer(wxVERTICAL);
  isizer->Add(ugrid,wxSizerFlags().Border().Expand());
  isizer->Add(plthist,wxSizerFlags(1));
  isizer->Add(rbsizer,wxSizerFlags().Border().Centre());
  tpanel->SetSizer(isizer);

  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniTune::OnApply,this,wxID_APPLY);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniTune::OnUndo,this,wxID_UNDO);
  Bind(wxEVT_UPDATE_UI,&MuniTune::OnUpdateApply,this,wxID_APPLY);
  Bind(wxEVT_UPDATE_UI,&MuniTune::OnUpdateUndo,this,wxID_UNDO);
  Bind(wxEVT_TEXT_ENTER,&MuniTune::OnEnterBlack,this,entry_black->GetId());
  Bind(wxEVT_TEXT_ENTER,&MuniTune::OnEnterSense,this,entry_sense->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnQrange,this,rb99->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnQrange,this,rb999->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnQrange,this,rbfull->GetId());
  Bind(EVT_HISTO,&MuniTune::OnHistoFinish,this,ID_HISTO_FINISH);
}

void MuniTune::FeedIttTab(wxPanel *tpanel)
{
  itt_line = new wxRadioButton(tpanel,wxID_ANY,FitsItt::Type_str(ITT_LINE),
			       wxDefaultPosition,wxDefaultSize,wxRB_GROUP);
  itt_tanh = new wxRadioButton(tpanel,wxID_ANY,FitsItt::Type_str(ITT_TANH));
  itt_asinh = new wxRadioButton(tpanel,wxID_ANY,FitsItt::Type_str(ITT_ASINH));
  itt_snlike = new wxRadioButton(tpanel,wxID_ANY,FitsItt::Type_str(ITT_SNLIKE));
  itt_sqr = new wxRadioButton(tpanel,wxID_ANY,FitsItt::Type_str(ITT_SQR));
  itt_photo = new wxRadioButton(tpanel,wxID_ANY,FitsItt::Type_str(ITT_PHOTO));
  SetIttRadio();
  itt_line->SetToolTip("Displayed intensities are linearly proportional "
		       "to the input values");
  itt_tanh->SetToolTip("This sigmoid intensity curve is equivalent to scaled "
		       "logistics function simulating some perception with "
		       "saturation.");
  itt_asinh->SetToolTip("A wide range asinh() function keeps faint details "
			"white the bright objects are still visible.");
  itt_snlike->SetToolTip("This function roughly approximates Signal-to-Noise "
			 "ratio. One offers a nice contrast.");
  itt_sqr->SetToolTip("Square of intensities does emphasise of faint objects. "
		      " Good for varying background.");
  itt_photo->SetToolTip("It is simulation of the classic photo-emulsion having "
			"a low dynamical range. Great for vintage photos.");

  pltitt = new MuniPlotItt(tpanel,wxSize(HIST_WIDTH,HIST_HEIGHT),FitsItt(itt));

  wxSizerFlags radio_flags;
  radio_flags.Border().DoubleBorder(wxLEFT);

  wxFlexGridSizer *grid = new wxFlexGridSizer(2);
  grid->Add(itt_line,radio_flags);
  grid->Add(itt_snlike,radio_flags);
  grid->Add(itt_asinh,radio_flags);
  grid->Add(itt_tanh,radio_flags);
  grid->Add(itt_sqr,radio_flags);
  grid->Add(itt_photo,radio_flags);

  wxSizerFlags label_flags;
  label_flags.Align(wxALIGN_CENTER_VERTICAL).Border();

  wxBoxSizer *fsizer = new wxBoxSizer(wxHORIZONTAL);
  fsizer->Add(new wxStaticText(tpanel,wxID_ANY,"Function:"),label_flags);
  fsizer->Add(grid,wxSizerFlags().Align(wxALIGN_CENTER|wxALIGN_CENTER_VERTICAL));

  wxBoxSizer *isizer = new wxBoxSizer(wxVERTICAL);
  isizer->Add(fsizer,wxSizerFlags(1).Centre());
  isizer->Add(pltitt);
  tpanel->SetSizer(isizer);

  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnItt,this,itt_line->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnItt,this,itt_asinh->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnItt,this,itt_tanh->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnItt,this,itt_snlike->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnItt,this,itt_sqr->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnItt,this,itt_photo->GetId());

}

void MuniTune::FeedColourTab(wxPanel *panel)
{
  wp_d65 = new wxRadioButton(panel,wxID_ANY,"D65",wxDefaultPosition,wxDefaultSize,
			     wxRB_GROUP);
  wp_d50 = new wxRadioButton(panel,wxID_ANY,"D50");
  wp_tungsten = new wxRadioButton(panel,wxID_ANY,"Tungsten");
  wp_flour = new wxRadioButton(panel,wxID_ANY,"Fluorescent");
  wp_led = new wxRadioButton(panel,wxID_ANY,"LED");
  wp_equal = new wxRadioButton(panel,wxID_ANY,"E");
  SetWhitePointRadio();
  wp_d65->SetToolTip("The most common D65 illuminant represents average "
		     "spectral distribution during bright days of our lives");
  wp_d50->SetToolTip("Lightly reddered D65; suitable for printing");
  wp_tungsten->SetToolTip("The old school light bulb with a tungsten (wolfram) "
			  "filament");
  wp_flour->SetToolTip("A common floursescent light tube");
  wp_led->SetToolTip("LED based light sources");
  wp_equal->SetToolTip("An artificial light source with equal energy spectral "
		       "distribution");

  wxSizerFlags radio_flags;
  radio_flags.Border().DoubleBorder(wxLEFT);

  wxFlexGridSizer *grid = new wxFlexGridSizer(3);
  grid->Add(wp_d65,radio_flags);
  grid->Add(wp_d50,radio_flags);
  grid->Add(wp_tungsten,radio_flags);
  grid->Add(wp_flour,radio_flags);
  grid->Add(wp_led,radio_flags);
  grid->Add(wp_equal,radio_flags);

  int w = 22*GetCharWidth();
  int s = int(10*saturation + 0.5);
  slider_sat = new wxSlider(panel,ID_COLOUR_SATUR,s,0,20,wxDefaultPosition,
			    wxSize(w,-1),wxSL_BOTTOM);
  slider_sat->SetToolTip("Colour saturation sets amount of colours in the image.");

  pltsat = new MuniPlotSaturation(panel,wxSize(HIST_WIDTH,HIST_HEIGHT),
				  saturation,x_white,y_white);

  wxSizerFlags label_flags;
  label_flags.Align(wxALIGN_CENTER_VERTICAL).Border();
  wxSizerFlags entry_flags;
  entry_flags.Align(wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT|wxBOTTOM);

  wxBoxSizer *wsizer = new wxBoxSizer(wxHORIZONTAL);
  wsizer->Add(new wxStaticText(panel,wxID_ANY,"Illuminant:"),label_flags);
  wsizer->Add(grid,wxSizerFlags().Align(wxALIGN_CENTER|wxALIGN_CENTER_VERTICAL));

  wxBoxSizer *asizer = new wxBoxSizer(wxHORIZONTAL);
  asizer->Add(new wxStaticText(panel,wxID_ANY,"Saturation:"),label_flags);
  asizer->Add(slider_sat,entry_flags);

  wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);
  sizer->Add(wsizer,wxSizerFlags().Center().Border().DoubleBorder(wxTOP|wxBOTTOM));
  sizer->Add(asizer,wxSizerFlags().Center().Border());
  sizer->Add(pltsat,wxSizerFlags(1));
  panel->SetSizer(sizer);

  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnWhitePoint,this,wp_d65->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnWhitePoint,this,wp_d50->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnWhitePoint,this,wp_tungsten->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnWhitePoint,this,wp_flour->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnWhitePoint,this,wp_led->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnWhitePoint,this,wp_equal->GetId());
  Bind(wxEVT_SLIDER,&MuniTune::OnSat,this,slider_sat->GetId());
}

// Night vision
void MuniTune::FeedNiteTab(wxPanel *panel)
{
  nitecheck = new wxCheckBox(panel,ID_CHECK_NITE,"Night vision");
  nitecheck->SetValue(nitevision);
  nitecheck->SetToolTip("Simulates human vision under low-light conditions.");

  int w = 12*GetCharWidth();
  entry_meso_level = new wxTextCtrl(panel,ID_COLOUR_MESO_LEVEL,
				    ToString(SCALE_FMT,meso_level),
				    wxDefaultPosition,wxSize(w,-1),
				    wxTE_PROCESS_ENTER);
  entry_meso_level->SetToolTip("Mesopic level is an intensity in which "
			       "human photopic (daylight) vision does down-grade"
			       " to scotopic (nightlight) vision.");
  entry_meso_width = new wxTextCtrl(panel,ID_COLOUR_MESO_WIDTH,
				    ToString(SCALE_FMT,meso_width),
				    wxDefaultPosition,wxSize(w,-1),
				    wxTE_PROCESS_ENTER);
  entry_meso_width->SetToolTip("The width of mesotopic layer "
                               "(a transition regime connecting "
			       "the scotopic and photopic vision).");

  wxButton *apply = new wxButton(panel,wxID_ANY,"Apply");

  wxSizerFlags label_flags;
  label_flags.Align(wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT);
  wxSizerFlags entry_flags;
  entry_flags.Align(wxALIGN_CENTER_VERTICAL).Border(wxLEFT|wxRIGHT|wxBOTTOM);
  wxSizerFlags button_flags;
  button_flags.Align(wxALIGN_CENTER_VERTICAL);

  wxFlexGridSizer *grid = new wxFlexGridSizer(2);
  grid->Add(new wxStaticText(panel,wxID_ANY,"Mesopic level:"),label_flags);
  grid->Add(entry_meso_level,entry_flags);
  grid->Add(new wxStaticText(panel,wxID_ANY,"Mesopic width:"),label_flags);
  grid->Add(entry_meso_width,entry_flags);

  wxBoxSizer *ugrid = new wxBoxSizer(wxHORIZONTAL);
  ugrid->AddStretchSpacer(1);
  ugrid->Add(grid);
  ugrid->Add(apply,button_flags.DoubleBorder());

  pltnite = new MuniPlotNite(panel,wxSize(HIST_WIDTH,HIST_HEIGHT),
			     meso_level,meso_width,nitevision);

  wxBoxSizer *tsizer = new wxBoxSizer(wxVERTICAL);
  tsizer->Add(nitecheck,wxSizerFlags().Center().DoubleBorder(wxTOP));
  tsizer->Add(ugrid,wxSizerFlags().Align(wxALIGN_CENTER).Border());
  tsizer->Add(pltnite,wxSizerFlags(1));

  panel->SetSizer(tsizer);

  Bind(wxEVT_COMMAND_CHECKBOX_CLICKED,&MuniTune::OnCheckNite,this,ID_CHECK_NITE);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniTune::OnApplyNite,this,apply->GetId());
  Bind(wxEVT_UPDATE_UI,&MuniTune::OnUpdateNite,this,ID_COLOUR_MESO_WIDTH);
  Bind(wxEVT_UPDATE_UI,&MuniTune::OnUpdateNite,this,ID_COLOUR_MESO_LEVEL);
  Bind(wxEVT_UPDATE_UI,&MuniTune::OnUpdateNiteApply,this,apply->GetId());
  Bind(wxEVT_TEXT_ENTER,&MuniTune::OnEnterMesoLevel,this,ID_COLOUR_MESO_LEVEL);
  Bind(wxEVT_TEXT_ENTER,&MuniTune::OnEnterMesoWidth,this,ID_COLOUR_MESO_WIDTH);
}

void MuniTune::FeedPaletteTab(wxPanel *ppanel)
{
  pal_grey = new wxRadioButton(ppanel,wxID_ANY,FitsPalette::Type_str(PAL_GREY),
			       wxDefaultPosition,wxDefaultSize,wxRB_GROUP);
  pal_sepia = new wxRadioButton(ppanel,wxID_ANY,FitsPalette::Type_str(PAL_SEPIA));
  pal_royal = new wxRadioButton(ppanel,wxID_ANY,FitsPalette::Type_str(PAL_ROYAL));
  pal_colour = new wxRadioButton(ppanel,wxID_ANY,FitsPalette::Type_str(PAL_COLOUR));
  pal_highlight = new wxRadioButton(ppanel,wxID_ANY,
				    FitsPalette::Type_str(PAL_HIGHLIGHT));
  pal_rainbow = new wxRadioButton(ppanel,wxID_ANY,
				     FitsPalette::Type_str(PAL_RAINBOW));
  pal_madness = new wxRadioButton(ppanel,wxID_ANY,
				  FitsPalette::Type_str(PAL_MADNESS));
  pal_vga = new wxRadioButton(ppanel,wxID_ANY,FitsPalette::Type_str(PAL_VGA));
  pal_aips0 = new wxRadioButton(ppanel,wxID_ANY,FitsPalette::Type_str(PAL_AIPS0));

  switch(palette){
  case PAL_GREY:      pal_grey->SetValue(true); break;
  case PAL_SEPIA:     pal_sepia->SetValue(true); break;
  case PAL_ROYAL:     pal_royal->SetValue(true); break;
  case PAL_COLOUR:    pal_colour->SetValue(true); break;
  case PAL_HIGHLIGHT: pal_highlight->SetValue(true); break;
  case PAL_RAINBOW:   pal_rainbow->SetValue(true); break;
  case PAL_MADNESS:   pal_madness->SetValue(true); break;
  case PAL_VGA:       pal_vga->SetValue(true); break;
  case PAL_AIPS0:     pal_aips0->SetValue(true); break;
  }

  invcheck = new wxCheckBox(ppanel,ID_CHECK_INVERSE,"Invert");
  invcheck->SetValue(inverse);

  wxSizerFlags radio_flags;
  radio_flags.Border().DoubleBorder(wxLEFT);

  wxFlexGridSizer *grid = new wxFlexGridSizer(3);
  grid->Add(pal_grey,radio_flags);
  grid->Add(pal_sepia,radio_flags);
  grid->Add(pal_royal,radio_flags);
  grid->Add(pal_colour,radio_flags);
  grid->Add(pal_highlight,radio_flags);
  grid->Add(pal_rainbow,radio_flags);
  grid->Add(pal_madness,radio_flags);
  grid->Add(pal_vga,radio_flags);
  grid->Add(pal_aips0,radio_flags);

  wxSizerFlags block_flags;
  block_flags.Align(wxALIGN_CENTER_VERTICAL).TripleBorder();

  wxBoxSizer *fsizer = new wxBoxSizer(wxHORIZONTAL);
  fsizer->Add(grid,wxSizerFlags().Align(wxALIGN_CENTER|wxALIGN_CENTER_VERTICAL));
  fsizer->Add(invcheck,block_flags);

  pltpal = new MuniPlotPalette(ppanel,wxSize(HIST_WIDTH,HIST_HEIGHT),
			       palette,inverse,itt,black,sense);

  wxBoxSizer *tsizer = new wxBoxSizer(wxVERTICAL);
  tsizer->Add(fsizer,wxSizerFlags(1).Centre());
  tsizer->Add(pltpal);

  ppanel->SetSizer(tsizer);

  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnPal,this,pal_grey->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnPal,this,pal_sepia->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnPal,this,pal_royal->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnPal,this,pal_colour->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnPal,this,pal_highlight->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnPal,this,pal_rainbow->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnPal,this,pal_madness->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnPal,this,pal_vga->GetId());
  Bind(wxEVT_RADIOBUTTON,&MuniTune::OnPal,this,pal_aips0->GetId());
  Bind(wxEVT_COMMAND_CHECKBOX_CLICKED,&MuniTune::OnCheckInverse,this,
       ID_CHECK_INVERSE);
}

void MuniTune::OnPal(wxCommandEvent& event)
{
  if( event.GetId() == pal_grey->GetId() )
    palette = PAL_GREY;
  else if( event.GetId() == pal_sepia->GetId() )
    palette = PAL_SEPIA;
  else if( event.GetId() == pal_royal->GetId() )
    palette = PAL_ROYAL;
  else if( event.GetId() == pal_colour->GetId() )
    palette = PAL_COLOUR;
  else if( event.GetId() == pal_highlight->GetId() )
    palette = PAL_HIGHLIGHT;
  else if( event.GetId() == pal_rainbow->GetId() )
    palette = PAL_RAINBOW;
  else if( event.GetId() == pal_madness->GetId() )
    palette = PAL_MADNESS;
  else if( event.GetId() == pal_vga->GetId() )
    palette = PAL_VGA;
  else if( event.GetId() == pal_aips0->GetId() )
    palette = PAL_AIPS0;
  else
    wxFAIL_MSG("*** MuniTune::OnPal(): unreachable line");

  pltpal->SetPalette(palette);

  MuniTunePalEvent e(EVT_TUNE_PAL);
  e.palette = palette;
  e.inverse = inverse;
  wxQueueEvent(GetParent(),e.Clone());
}

void MuniTune::OnCheckInverse(wxCommandEvent& event)
{
  inverse = event.IsChecked();

  MuniTunePalEvent ev(EVT_TUNE_PAL,ID_PALETTE_INVERSE);
  ev.palette = palette;
  ev.inverse = inverse;
  wxQueueEvent(GetParent(),ev.Clone());

  pltpal->SetInverse(inverse);
}

void MuniTune::OnCheckNite(wxCommandEvent& event)
{
  nitevision = event.IsChecked();
  pltnite->SetEnabled(nitevision);

  MuniTuneNiteEvent ev(EVT_TUNE_NITE,ID_COLOUR_NITEVISION);
  ev.nite = nitevision;
  ev.level = meso_level;
  ev.width = meso_width;
  wxQueueEvent(GetParent(),ev.Clone());
}


void MuniTune::OnApplyNite(wxCommandEvent& event)
{
  NiteUpdated();
}

void MuniTune::OnEnterMesoLevel(wxCommandEvent& event)
{
  NiteUpdated();
}

void MuniTune::OnEnterMesoWidth(wxCommandEvent& event)
{
  NiteUpdated();
}


void MuniTune::NiteUpdated()
{
  wxString uentry = entry_meso_level->GetValue();
  wxString ventry = entry_meso_width->GetValue();
  double u,v;

  if( uentry.ToDouble(&u) && ventry.ToDouble(&v) && v > 0.0 &&
      fabs(u) < FLT_MAX && fabs(v) < FLT_MAX ) {
    meso_level = u;
    meso_width = v;
    pltnite->SetMeso(meso_level,meso_width);

    MuniTuneNiteEvent ev(EVT_TUNE_NITE);
    ev.nite = nitevision;
    ev.level = meso_level;
    ev.width = meso_width;
    wxQueueEvent(GetParent(),ev.Clone());
  }
}

void MuniTune::OnSat(wxCommandEvent& event)
{
  wxASSERT(event.GetId() == ID_COLOUR_SATUR);

  saturation = 0.1*event.GetInt();

  pltsat->SetSaturation(saturation);

  MuniTuneColourEvent ev(EVT_TUNE_COLOUR);
  ev.saturation = saturation;
  ev.x_white = x_white;
  ev.y_white = y_white;
  wxQueueEvent(GetParent(),ev.Clone());
}

void MuniTune::OnApply(wxCommandEvent& event)
{
  double u,v;
  wxString bentry = entry_black->GetValue();
  wxString sentry = entry_sense->GetValue();
  if( bentry.ToDouble(&u) && sentry.ToDouble(&v) && v > 0.0 &&
      fabs(u) < FLT_MAX && fabs(v) < FLT_MAX ) {
    black = u;
    sense = v;
    ScaleUpdated();
    undo = true;
  }
  else
    wxFAIL_MSG("*** MuniTune::OnApplyScale(): Values check fail");
}

void MuniTune::OnUndo(wxCommandEvent& event)
{
  MuniTuneScaleEvent e(EVT_TUNE_SCALE,event.GetId());
  e.reset = true;
  wxQueueEvent(GetParent(),e.Clone());
  undo = false;
}

void MuniTune::OnQrange(wxCommandEvent& event)
{
  if( event.GetId() == rb99->GetId() )
    qrange = 0.99;
  else if( event.GetId() == rb999->GetId() )
    qrange = 0.999;
  else if( event.GetId() == rbfull->GetId() )
    qrange = 1.0;

  RenderHisto();
}

void MuniTune::OnItt(wxCommandEvent& event)
{
  if( event.GetId() == itt_line->GetId() )
    itt = ITT_LINE;
  else if( event.GetId() == itt_asinh->GetId() )
    itt = ITT_ASINH;
  else if( event.GetId() == itt_tanh->GetId() )
    itt = ITT_TANH;
  else if( event.GetId() == itt_snlike->GetId() )
    itt = ITT_SNLIKE;
  else if( event.GetId() == itt_sqr->GetId() )
    itt = ITT_SQR;
  else if( event.GetId() == itt_photo->GetId() )
    itt = ITT_PHOTO;
  else
    wxFAIL_MSG("*** MuniTune::OnItt(): unreachable line");

  pltitt->SetItt(itt);
  if( pltpal )
    pltpal->SetItt(itt);

  MuniTuneIttEvent e(EVT_TUNE_ITT);
  e.itt = itt;
  wxQueueEvent(GetParent(),e.Clone());
}


void MuniTune::OnWhitePoint(wxCommandEvent& event)
{
  // https://en.wikipedia.org/wiki/Standard_illuminant

  if( event.GetId() == wp_d65->GetId() ) {
    x_white = 0.31271; y_white = 0.32902;
  }
  else if( event.GetId() == wp_d50->GetId() ) {
    x_white = 0.34567; y_white = 0.35850;
  }
  else if( event.GetId() == wp_tungsten->GetId() ) {
    x_white = 0.44757; y_white = 0.40745;
  }
  else if( event.GetId() == wp_flour->GetId() ) {
    x_white = 0.40910; y_white = 0.39430;  // F3
  }
  else if( event.GetId() == wp_led->GetId() ) {
    x_white = 0.4557; y_white = 0.4211;    // LED-RGB1
  }
  else if( event.GetId() == wp_equal->GetId() ) {
    x_white = 0.33333; y_white = 0.33333;
  }

  pltsat->SetWhitePoint(x_white,y_white);

  MuniTuneColourEvent e(EVT_TUNE_COLOUR);
  e.saturation = saturation;
  e.x_white = x_white;
  e.y_white = y_white;
  wxQueueEvent(GetParent(),e.Clone());
}

void MuniTune::OnUpdateNite(wxUpdateUIEvent& event)
{
  event.Enable(nitevision);
}

void MuniTune::OnUpdateNiteApply(wxUpdateUIEvent& event)
{
  double u,v;
  wxString b = entry_meso_level->GetLineText(0);
  wxString s = entry_meso_width->GetLineText(0);
  event.Enable(b.ToDouble(&u) && s.ToDouble(&v) && v > 0.0 && nitevision &&
	       fabs(u) < FLT_MAX && fabs(v) < FLT_MAX);
}

void MuniTune::OnUpdateUndo(wxUpdateUIEvent& event)
{
  event.Enable(undo);
}

void MuniTune::OnUpdateApply(wxUpdateUIEvent& event)
{
  double u,v;
  wxString b = entry_black->GetLineText(0);
  wxString s = entry_sense->GetLineText(0);
  event.Enable(b.ToDouble(&u) && s.ToDouble(&v) && v > 0.0 &&
	       fabs(u) < FLT_MAX && fabs(v) < FLT_MAX);
}


void MuniTune::OnEnterBlack(wxCommandEvent& event)
{
  wxASSERT(event.GetId() == ID_TONE_BLACK);
  ScaleUpdated();
  undo = true;
}

void MuniTune::OnEnterSense(wxCommandEvent& event)
{
  wxASSERT(event.GetId() == ID_TONE_SENSE);
  ScaleUpdated();
  undo = true;
}

void MuniTune::ScaleUpdated()
{

  double u,v;
  wxString b = entry_black->GetLineText(0);
  wxString s = entry_sense->GetLineText(0);

  if( b.ToDouble(&u) && s.ToDouble(&v) && v > 0.0 &&
      fabs(u) < FLT_MAX && fabs(v) < FLT_MAX ) {
    black = u;
    sense = v;

    plthist->SetScale(black,sense);
    if( pltpal )
      pltpal->SetScale(black,sense);

    MuniTuneScaleEvent e(EVT_TUNE_SCALE);
    e.black = black;
    e.sense = sense;
    wxQueueEvent(GetParent(),e.Clone());
  }
}

wxString MuniTune::ToString(const char *fmt, double x) const
{
  wxString a;
  a.Printf(fmt,x);
  return a;
}

void MuniTune::SetIttRadio()
{
  if( itt == ITT_LINE )
    itt_line->SetValue(true);
  else if( itt == ITT_ASINH )
    itt_asinh->SetValue(true);
  else if( itt == ITT_TANH )
    itt_tanh->SetValue(true);
  else if( itt == ITT_SNLIKE )
    itt_snlike->SetValue(true);
  else if( itt == ITT_SQR )
    itt_sqr->SetValue(true);
  else if( itt == ITT_PHOTO )
    itt_photo->SetValue(true);
}

void MuniTune::SetWhitePointRadio()
{
  const float e = 1e-3;

  if( fabs(x_white - 0.31271) < e && fabs(y_white - 0.32902) < e )
    wp_d65->SetValue(true);
  else if( fabs(x_white - 0.34567) < e && fabs(y_white - 0.35850) < e )
    wp_d50->SetValue(true);
  else if( fabs(x_white - 0.44757) < e && fabs(y_white - 0.40745) < e )
    wp_tungsten->SetValue(true);
  else if( fabs(x_white - 0.40910) < e && fabs(y_white - 0.39430) < e )
    wp_flour->SetValue(true);
  else if( fabs(x_white - 0.4557) < e && fabs(y_white - 0.4211) < e )
    wp_led->SetValue(true);
  else if( fabs(x_white - 0.3333) < e && fabs(y_white - 0.3333) < e )
    wp_equal->SetValue(true);

}


void MuniTune::SetScale(double b, double s)
{
  black = b;
  sense = s;
  entry_black->SetValue(ToString(SCALE_FMT,black));
  entry_sense->SetValue(ToString(SCALE_FMT,sense));
  plthist->SetScale(black,sense);
  if( pltpal )
    pltpal->SetScale(black,sense);
}


void MuniTune::RenderHisto()
{
  if( thread_histo )
    thread_histo->Delete();

  thread_histo = new MuniHistoThread(this,array,qrange);
  if( thread_histo->Run() != wxTHREAD_NO_ERROR ) {
    wxLogError("Can't create the thread!");
    delete thread_histo;
    thread_histo = NULL;
  }
}


void MuniTune::OnHistoFinish(MuniHistoEvent& e)
{
  //  wxLogDebug("MuniTune::OnHistoFinish()");

  thread_histo = NULL;

  if( shutdown ) {
    Close();
    return;
  }

  wxASSERT(plthist);
  plthist->SetHisto(e.hist);

  if( pltnite )
    pltnite->SetHisto(e.hist);
}


// ---------------------------------------------------------------------

MuniHistoThread::MuniHistoThread(wxEvtHandler *h, const FitsArray& a, double q):
  wxThread(wxTHREAD_DETACHED),handler(h),array(a),qrange(q)
{
  wxASSERT(handler);
}

wxThread::ExitCode MuniHistoThread::Entry()
{
  size_t npix = array.GetWidth()*array.GetHeight();
  const int maxdat = 262144; // =2**18
  const int stride = std::max(int(log(float(npix)/float(maxdat)) / log(2.0)),1);

  //wxLogDebug("stride=%d %d %d",stride,int(array.GetWidth()),int(array.GetHeight()));

  int idx = array.IsColour() ? 1 : 0;
  FitsArray a(array.GetGrid(idx,stride));
  FitsHisto hist(a.Npixels(),a.PixelData(),qrange);

  if( TestDestroy() )
    return (wxThread::ExitCode) 1;

  wxLogDebug("MuniHistoThread()");

  // Only the intensity component histogram is kept.
  MuniHistoEvent e(EVT_HISTO,ID_HISTO_FINISH);
  e.hist = hist;
  wxQueueEvent(handler,e.Clone());

  //  wxLogDebug("MuniHistoThread::Entry() FINISH");
  return (wxThread::ExitCode) 0;
}
