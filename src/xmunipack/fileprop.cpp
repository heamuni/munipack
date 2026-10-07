/*

  xmunipack - file properties dialog

  Copyright © 2009-2013, 2019-2022 F.Hroch (hroch@physics.muni.cz)

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

#include "fileprop.h"
#include "fits.h"
#include <wx/wx.h>
#include <wx/propdlg.h>
#include <wx/bookctrl.h>
#include <wx/filename.h>
#include <wx/uri.h>
#include <vector>

#ifdef __WXDEBUG__
#include <wx/debug.h>
#endif

using namespace std;


MuniFileProperties::MuniFileProperties(wxWindow *w, const MuniConfig *config,
				       const FitsMeta& meta, const wxImage& cover,
				       const std::vector<wxImage>& ilist):
  wxDialog(w,wxID_ANY,"Properties of "+wxURI::Unescape(meta.GetName()))
{
  SetIcon(config->munipack_icon);

  wxSizerFlags label_flags, value_flags;
  label_flags.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).DoubleBorder(wxRIGHT);
  value_flags.Align(wxALIGN_CENTER_VERTICAL).Expand();

  wxString a;

  wxArrayString flabels;
  flabels.Add("Name");
  flabels.Add("Format");
  flabels.Add("Type");
  flabels.Add("Hdu(s)");
  flabels.Add("Size");
  flabels.Add("Access");
  flabels.Add("Modified");
  flabels.Add("Location");

  wxArrayString fvals;
  // name
  fvals.Add(wxURI::Unescape(meta.GetFullName()));

  // format
  fvals.Add("FITS");

  // type
  fvals.Add(meta.Type_str());

  // HDUs
  a.Printf("%d",(int) meta.HduCount());
  fvals.Add(a);

  // bytes
  fvals.Add(meta.GetHumanReadableSize());

  // access
  wxFileName fn(wxURI::Unescape(meta.GetFullPath()));
  if( fn.IsFileReadable() && fn.IsFileWritable() )
    fvals.Add("Read-Write");
  else if( fn.IsFileReadable() && ! fn.IsFileWritable() )
    fvals.Add("Read only");
  else if( ! fn.IsFileReadable() &&  fn.IsFileWritable() )
    fvals.Add("Write only");
  else
    fvals.Add(wxEmptyString);

  // mtime
  fvals.Add(meta.Mtime());

  // location
  fvals.Add(wxURI::Unescape(meta.GetPath()));

  wxNotebook *book = new wxNotebook(this,wxID_ANY);

  wxPanel* gpanel = new wxPanel(book);

  wxBoxSizer *gtop = new wxBoxSizer(wxHORIZONTAL);
  wxBoxSizer *gleft = new wxBoxSizer(wxVERTICAL);
  MuniThumbCanvas *th = new MuniThumbCanvas(gpanel,cover);
  gleft->Add(th,wxSizerFlags().Border().Center());
  gtop->Add(gleft,wxSizerFlags().Border(wxALL));

  wxSizerFlags flags; flags.Align(wxALIGN_BOTTOM);
  wxSizerFlags lflags; flags.Align(wxALIGN_BOTTOM).Border(wxRIGHT);

  wxFlexGridSizer *gright = new wxFlexGridSizer(2);

  // create the table
  for(size_t i = 0; i < flabels.GetCount(); i++) {

    wxString label("<span fgcolor=\"gray\">"+flabels.Item(i)+"</span>");
    wxStaticText *l = new wxStaticText(gpanel, wxID_ANY, "");
    l->SetLabelMarkup(label);
    gright->Add(l,label_flags);
    wxStaticText *x = new wxStaticText(gpanel, wxID_ANY, fvals.Item(i));
    gright->Add(x,value_flags);

  }

  gtop->Add(gright,wxSizerFlags().Border(wxALL));

  gpanel->SetSizer(gtop);

  book->AddPage(gpanel,"File");

  for(size_t k = 0; k < meta.HduCount(); k++ ) {

    FitsMetaHdu hdu = meta.Hdu(k);

    wxArrayString hlabels, hvals;
    wxString title("Component");

    // extname
    a = hdu.GetKey("EXTNAME");
    if( ! a.IsEmpty() ) {
      hlabels.Add("Extname");
      hvals.Add(a);
    }

    // obj
    hlabels.Add("Object");
    hvals.Add(hdu.GetKey(config->fits_key_object));

    // obs
    hlabels.Add("Observer");
    hvals.Add(hdu.GetKey(config->fits_key_observer));

    if( hdu.Type() == HDU_IMAGE ) {

      title = "Image";

      // axes
      hlabels.Add("Axes");
      a.Printf("%d",(int) hdu.Naxis());
      hvals.Add(a);

      switch(hdu.SubType()) {

      case HDU_IMAGE_LINE:
	hlabels.Add("Length");
	a.Printf("%d pix",(int) hdu.Width());
	hvals.Add(a);
	break;

      case HDU_IMAGE_FRAME:
	hlabels.Add("Width");
	hlabels.Add("Height");
	a.Printf("%d pix",(int) hdu.Width());
	hvals.Add(a);
	a.Printf("%d pix",(int) hdu.Height());
	hvals.Add(a);
	break;

      case HDU_IMAGE_COLOUR:
	hlabels.Add("Bands");
	hlabels.Add("Width");
	hlabels.Add("Height");
	a.Printf("%d",(int) hdu.Naxes(2));
	hvals.Add(a);
	a.Printf("%d pix",(int) hdu.Width());
	hvals.Add(a);
	a.Printf("%d pix",(int) hdu.Height());
	hvals.Add(a);
	break;

      case HDU_IMAGE_CUBE:
	hlabels.Add("Width");
	hlabels.Add("Height");
	hlabels.Add("Depth");
	a.Printf("%d pix",(int) hdu.Width());
	hvals.Add(a);
	a.Printf("%d pix",(int) hdu.Height());
	hvals.Add(a);
	a.Printf("%d pix",(int) hdu.Naxes(2));
	hvals.Add(a);
	break;

      default: break;
      }
    }
    else if( hdu.Type() == HDU_TABLE ) {

      title = "Table";

      hlabels.Add("Rows");
      hlabels.Add("Columns");
      a.Printf("%d",(int) hdu.Nrows());
      hvals.Add(a);
      a.Printf("%d",(int) hdu.Ncols());
      hvals.Add(a);
    }
    else if( hdu.Type() == HDU_HEAD ) {

      title = "Head";

    }

    // date & time
    hlabels.Add("Date");
    hlabels.Add("Time");
    FitsTime t(hdu.GetKey(config->fits_key_dateobs));
    hvals.Add(t.GetDate());
    hvals.Add(t.GetTime());

    // filter
    hlabels.Add("Filter");
    hvals.Add(hdu.GetKey(config->fits_key_filter));

    // exp
    hlabels.Add("Exposure");
    a = hdu.Exposure_str(config->fits_key_exptime);
    if( ! a.IsEmpty() )
       a = a + " sec";
    hvals.Add(a);

    // type
    hlabels.Add("Type");
    hvals.Add(hdu.Type_str() + " (" + hdu.SubType_str() + ")");

    // range
    hlabels.Add("Range");
    hvals.Add(hdu.Bitpix_str());


    wxPanel *panel = new wxPanel(book);

    wxBoxSizer *top = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *left = new wxBoxSizer(wxVERTICAL);

    MuniThumbCanvas *tb = new MuniThumbCanvas(panel,ilist[k]);
    left->Add(tb,wxSizerFlags().Border().Center());
    top->Add(left,wxSizerFlags().Border(wxALL));

    wxFlexGridSizer *right = new wxFlexGridSizer(2);

    // create the table
    for(size_t i = 0; i < hlabels.GetCount(); i++) {

      wxString label("<span fgcolor=\"gray\">"+hlabels.Item(i)+"</span>");
      wxStaticText *l = new wxStaticText(panel, wxID_ANY,"");
      l->SetLabelMarkup(label);
      right->Add(l,label_flags);
      wxStaticText *x = new wxStaticText(panel, wxID_ANY, hvals.Item(i));
      right->Add(x,value_flags);

    }

    top->Add(right,wxSizerFlags().Border(wxALL));

    panel->SetSizer(top);

    if( meta.HduCount() == 1 )
      a = "Image";
    else {
      a = hdu.GetKey("EXTNAME");
      if( a.IsEmpty() )
	a.Printf(title+" %d",(int) k);
    }
    book->AddPage(panel,a);

  }

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(book,wxSizerFlags(1).Expand().Border());
  wxSizer *butt = CreateButtonSizer(wxCLOSE);
  if( butt ) {
    topsizer->Add(butt,wxSizerFlags().Right().Border(wxBOTTOM));
    SetEscapeId(wxID_CLOSE);
  }
  SetSizerAndFit(topsizer);

}


// ----  MuniDirProp


MuniDirProperties::MuniDirProperties(wxWindow *w, const MuniConfig *config,
				     const vector<FitsMeta>& flist):
  wxDialog(w,wxID_ANY,"Properties")
{
  SetIcon(config->munipack_icon);
  EnableCloseButton();

  wxSizerFlags label_flags, value_flags;
  label_flags.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT).DoubleBorder(wxRIGHT);
  value_flags.Align(wxALIGN_CENTER_VERTICAL).Expand();


  wxNotebook *book = new wxNotebook(this,wxID_ANY);

  wxPanel* gpanel = new wxPanel(book);

  wxBoxSizer *gtop = new wxBoxSizer(wxHORIZONTAL);
  wxBoxSizer *gleft = new wxBoxSizer(wxVERTICAL);

  MuniArtIcons icons(wxART_MESSAGE_BOX,wxSize(48,48));

  MuniThumbCanvas *th = new MuniThumbCanvas(gpanel,icons.Icon(wxART_FOLDER));
  gleft->Add(th,wxSizerFlags().Center());
  gtop->Add(gleft,wxSizerFlags().Border(wxALL));

  wxString a;
  a.Printf("%d",(int) flist.size());

  wxFlexGridSizer *gright = new wxFlexGridSizer(2);
  gright->AddGrowableCol(1);
  wxStaticText *l1 = new wxStaticText(gpanel, wxID_ANY, "");
  l1->SetLabelMarkup("<span fgcolor=\"gray\">Total images</span>");
  gright->Add(l1,label_flags);
  gright->Add(new wxStaticText(gpanel, wxID_ANY,a),value_flags);
  wxStaticText *l2 = new wxStaticText(gpanel, wxID_ANY, "");
  l2->SetLabelMarkup("<span fgcolor=\"gray\">Total size</span>");
  gright->Add(l2,label_flags);

  wxULongLong size = 0;
  for(vector<FitsMeta>::const_iterator i = flist.begin(); i != flist.end();++i){
    wxULongLong s = i->GetSize();
    if( s == wxInvalidSize ) {
      size = s;
      break;
    }
    size += s;
  }

  a = wxFileName::GetHumanReadableSize(size);
  wxStaticText *ss = new wxStaticText(gpanel, wxID_ANY, a);
  gright->Add(ss,value_flags);
  gtop->Add(gright,wxSizerFlags(1).Center().Border(wxALL));

  gpanel->SetSizer(gtop);

  book->AddPage(gpanel,"General");

  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);
  topsizer->Add(book,wxSizerFlags(1).Expand().Border());
  wxSizer *butt = CreateButtonSizer(wxCLOSE);
  if( butt ) {
    topsizer->Add(butt,wxSizerFlags().Right().Border());
    SetEscapeId(wxID_CLOSE);
  }
  SetSizerAndFit(topsizer);
}
