/*

  xmunipack - astrometry

  Copyright © 2011-5,2017 F.Hroch (hroch@physics.muni.cz)

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
#include <wx/collpane.h>
#include <wx/statline.h>
#include <wx/regex.h>
#include <wx/utils.h>
#include <wx/sstream.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>
#include <wx/filefn.h>
#include <wx/richtooltip.h>
#include <list>
#include <cfloat>


#define ARCUAS L"μas"
#define ARCMAS  "mas"
#define ARCSEC  "arcsec"
#define ARCMIN  "arcmin"
#define ARCDEG  "deg"

#define ARC_UAS ARCUAS "/pix"
#define ARC_MAS ARCMAS "/pix"
#define ARC_SEC ARCSEC "/pix"
#define ARC_MIN ARCMIN "/pix"
#define ARC_DEG ARCDEG "/pix"

#define MAXGRANGE 1000

using namespace std;

// --- MuniAstrometryLayer


class MuniAstrometryLayer
{
public:
  MuniAstrometryLayer(int w, int h):
    ltic(10),gold(wxColour(255,215,0,240)),DarkOrange2(wxColour(238,118,0)),
    width(w),height(h),alpha(0.0),delta(0.0),xoff(0.0),yoff(0.0),scale(1.0),
    angle(0.0),reflex(1.0),maglim(18.0) {}
  MuniLayer GetLayer() const;
  bool IsEmpty() const { return objects.empty(); }

  void SetProj(const wxString&,double,double,double,double,double,double,
	       double);
  void DrawStars(const FitsTable& table, const wxString&,const wxString&,
		 const wxString&,const wxString&,const wxString&);
  void DrawCompas();
  void DrawMatch(double,const std::vector<double>&,const std::vector<double>&);
  void Draw3Space(int, const std::vector<int>&);
  void Draw3Tri(const std::vector<double>&,const std::vector<double>&);

private:

  const int ltic;
  const wxColour gold,DarkOrange2;
  std::vector<wxObject *> objects;
  int width, height;
  wxString proj;
  double alpha,delta,xoff,yoff,scale,angle,reflex,maglim;
  wxBitmap tribit;


};

MuniLayer MuniAstrometryLayer::GetLayer() const
{
  MuniLayer layer(ID_ASTROMETRY,objects);
  return layer;
}

void MuniAstrometryLayer::SetProj(const wxString& t,double a,double d,double x,double y,double s,double f, double z)
{
  proj = t;
  alpha = a;
  delta = d;
  xoff = x;
  yoff = y;
  scale = s;
  angle = f;
  reflex = z;
  //wxLogDebug(proj+" %f %f %f %f %f %f",a,d,x,y,s,f);
}

void MuniAstrometryLayer::DrawStars(const FitsTable& table,
				    const wxString& label_ra,
				    const wxString& label_dec,
				    const wxString& label_pmra,
				    const wxString& label_pmdec,
				    const wxString& label_mag)
{
  objects.push_back(new MuniDrawPen(wxPen(gold,1.0)));
  objects.push_back(new MuniDrawFont(*wxNORMAL_FONT,gold));
  objects.push_back(new MuniDrawBrush(wxColour(30,30,205,190)));

  // stars
  if( table.IsOk() ) {

    int n = table.Nrows();

    int ira = table.GetColIndex(label_ra);
    int idec = table.GetColIndex(label_dec);
    int imag = table.GetColIndex(label_mag);
    int ipmra = table.GetColIndex(label_pmra);
    int ipmdec = table.GetColIndex(label_pmdec);

    if( ira >= 0 && idec >= 0 ) {

      const double *a = table.GetColumn(label_ra).GetCol_double();
      const double *d = table.GetColumn(label_dec).GetCol_double();
      const float *mag = imag >= 0 ? table.GetColumn(label_mag).GetCol_float() : 0;
      double *ppma = ipmra  >= 0 ?
	(double *)table.GetColumn(label_pmra).GetCol_double() : 0;
      double *ppmd = ipmdec >= 0 ?
	(double *)table.GetColumn(label_pmdec).GetCol_double() : 0;

      objects.push_back(new MuniDrawPen(wxPen(gold,1.8)));
      objects.push_back(new MuniDrawBrush(wxColour(90,90,255,190)));

      FitsProjection pr(proj,alpha,delta,xoff,yoff,scale,angle,reflex);

      double dt = 2012.8 - 2000.0;
      //      dt = 0.0;

      for(int i = 0; i < n; i++) {

	if( wxFinite(a[i]) && wxFinite(d[i]) && wxFinite(ppma[i]) &&
	    wxFinite(ppmd[i]) && wxFinite(mag[i]) ) {

	  double x,y;
	  double pma = ppma ? dt*ppma[i]/3600000.0 : 0.0;
	  double pmd = ppmd ? dt*ppmd[i]/3600000.0 : 0.0;

	  //wxPrintf("pm: %f %f %f %f %f %f\n",a[i],d[i],pma,pmd,ppma[i],ppmd[i]);

	  pr.ad2xy(a[i]+pma,d[i]+pmd,x,y);

	  //wxLogDebug("a,d: %f %f, x,y: %f %f",a[i]+pma,d[i]+pmd,x,y);

	  if( 0 < x && x < width && 0 < y && y < width ) {
	    if( mag ) {
	      double r = 1.5*pow(10.0,0.11*(maglim - mag[i]));
	      //	  double r = 3.0*pow(maglim - mag[i],0.333);
	      //	  if( mag[i] > 1.0 )
	      //	  wxLogDebug("%f %f %f",x,y,r);
	      objects.push_back(new MuniDrawCircle(x,y,r));
	    }
	    else
	      objects.push_back(new MuniDrawCircle(x,y,1.0));
	  }
	}
      }
    }
  }


  // scale
  /*
  const int xscale = 50;
  const int yscale = 30;
  double s = 60.0/Scale(scale);
  objects.push_back(new MuniDrawLine(xscale,yscale,xscale+s,yscale));
  objects.push_back(new MuniDrawText(xscale+s/2,yscale,"1'"));
  */
  //  wxLogDebug("%f %f",60.0/scale,s);

}


// central cross
void MuniAstrometryLayer::DrawCompas()
{
  objects.push_back(new MuniDrawLine(xoff-ltic,yoff,xoff+ltic,yoff));
  objects.push_back(new MuniDrawLine(xoff,yoff-ltic,xoff,yoff+ltic));

  const int rrose = 137;
  const double rad = 57.3;

  // compass rose
  for(int i = 0; i < 36; i++) {
    double f = 10.0*i;
    double r1 = rrose;
    double r2 = rrose + ltic;
    double x1 = xoff + r1*cos(f/rad);
    double y1 = yoff + r1*sin(f/rad);
    double x2 = xoff + r2*cos(f/rad);
    double y2 = yoff + r2*sin(f/rad);
    objects.push_back(new MuniDrawLine(x1,y1,x2,y2));
  }

  // big cross and cardinal points
  const char *cp[4] = { "W", "N", "E", "S" };
  for(int i = 0; i < 4; i++) {
    double f = 90.0*i - angle;
    double r1 = 3.0*ltic;
    double r2 = rrose - 4.0*ltic;
    double x1 = xoff + r1*cos(f/rad)*reflex;
    double y1 = yoff + r1*sin(f/rad);
    double x2 = xoff + r2*cos(f/rad)*reflex;
    double y2 = yoff + r2*sin(f/rad);
    objects.push_back(new MuniDrawLine(x1,y1,x2,y2));
    objects.push_back(new MuniDrawText(x2,y2,cp[i]));
  }
}

void MuniAstrometryLayer::DrawMatch(double sig, const std::vector<double>& x,const std::vector<double>& y)
{
  // print stars with position uncertainity by sig
  objects.push_back(new MuniDrawPen(wxPen(DarkOrange2,1.0)));
  objects.push_back(new MuniDrawBrush(*wxTRANSPARENT_BRUSH));
  for(size_t i = 0; i < x.size(); i++)
    objects.push_back(new MuniDrawCircle(x[i],y[i],3*sig));

  wxColour orange(DarkOrange2);
  wxColour colour(orange.Red(),orange.Green(),orange.Blue(),137);
  wxPen pen(colour);
  objects.push_back(new MuniDrawPen(pen));

  for(size_t i = 1; i < x.size(); i++)
    objects.push_back(new MuniDrawLine(x[i-1],y[i-1],x[i],y[i]));

  objects.push_back(new MuniDrawFont(*wxSMALL_FONT,colour));
  for(size_t i = 0; i < x.size(); i++) {
    wxString t;
    t.Printf("%d",(int)i+1);
    objects.push_back(new MuniDrawText(x[i]-sig-3,y[i]-sig-5,t));
  }
}

void MuniAstrometryLayer::Draw3Space(int nhist, const std::vector<int>& hist)
{
  wxASSERT(nhist*nhist == (int)hist.size());

  // triangle space
  int hmax = 0;
  for(size_t i = 0; i < hist.size(); i++)
    if( hist[i] > hmax )
      hmax = hist[i];

  // pix can be adjusted to get zoomed triangle space
  int pix = 1;
  int uvdim = pix*nhist;
  const int blevel = 64;
  wxImage img(uvdim,uvdim);
  for(int i = 0; i < nhist; i++) {
    for(int j = 0; j < nhist; j++) {
      wxRect rect(i*pix,uvdim-(j+1)*pix,pix,pix);
      double x = double(hist[j*nhist+i])/double(hmax);
      // Logistic function shifted and scaled
      x = 2.0/(1.0 + exp(-x/0.1618)) - 1.0;
      // color table: dark is blue, light is white
      int r = int(255*x);
      int g = r;
      int b = blevel + int((255-blevel)*x);
      img.SetRGB(rect,r,g,b);
    }
  }

  tribit = wxBitmap(img);
  wxMemoryDC mdc(tribit);
  if( mdc.IsOk() ) {
    mdc.SetBrush(*wxTRANSPARENT_BRUSH);
    mdc.SetPen(wxPen(DarkOrange2));
    mdc.DrawRectangle(0,0,uvdim,uvdim);
  }
  mdc.SelectObjectAsSource(wxNullBitmap);

  double x = width - 10;
  double y = height - 10;

  objects.push_back(new MuniDrawBitmap(tribit,x-uvdim,y-uvdim,uvdim,uvdim));
}

void MuniAstrometryLayer::Draw3Tri(const std::vector<double>& u, const std::vector<double>& v)
{
  wxASSERT(tribit.IsOk() && tribit.GetWidth() == tribit.GetHeight());

  int uvdim = tribit.GetWidth();
  const int q = 3;

  wxMemoryDC mdc(tribit);
  if( mdc.IsOk() ) {
    mdc.SetBrush(*wxTRANSPARENT_BRUSH);
    mdc.SetPen(wxPen(DarkOrange2));

    for(size_t i = 0; i < u.size(); i++) {
      int uu = int(uvdim*u[i]);
      int vv = int(uvdim*(1.0 - v[i]));

      mdc.DrawLine(uu-q,vv,uu+q,vv);
      mdc.DrawLine(uu,vv-q,uu,vv+q);
    }

    mdc.SetFont(*wxSMALL_FONT);
    mdc.SetTextForeground(DarkOrange2);
    wxString uv("u-v space");
    wxSize ps = mdc.GetTextExtent(uv);
    mdc.DrawText(uv,(uvdim - ps.GetWidth())/2,uvdim - ps.GetHeight());
  }

  mdc.SelectObjectAsSource(wxNullBitmap);

  double x = width - 10;
  double y = height - 10;

  objects.push_back(new MuniDrawBitmap(tribit,x-uvdim,y-uvdim,uvdim,uvdim));
}


// --- MuniAstrometry

MuniAstrometry::MuniAstrometry(wxWindow *w, MuniConfig *c):
  wxDialog(w,wxID_ANY,"Astrometry",wxDefaultPosition,wxDefaultSize,
	   wxDEFAULT_DIALOG_STYLE|wxRESIZE_BORDER),
  config(c),
  timer(this),pipe(this),
  nhist(0),
  maglim(18.0),
  proj(config->astrometry_proj),
  init_par(false),init_file(false),init_ref(false),
  running(false),unsaved(false),
  draw_overlay(true),edited(false),
  relative(false), tmpcat(false), readonly(false), showtooltip(true),page(0)
{
  SetIcon(config->munipack_icon);
  EnableCloseButton(false);

  CreateControls();

  calbutt->Enable(false);
}


MuniAstrometry::~MuniAstrometry()
{
  // delete temps

  config->astrometry_proj = proj;

  RemoveWorkingFile();

  EraseTemp();

  wxLogDebug("MuniAstrometry::~MuniAstrometry()");
}

void MuniAstrometry::EraseTemp()
{
  if( tmpcat && !catfile.IsEmpty() && wxFileExists(catfile) )
    wxRemoveFile(catfile);
}

bool MuniAstrometry::FinishClean(wxCommandEvent& event)
{
  wxLogDebug("MuniAstrometry::FinishClean %c %d",running ? 'T' : 'F',event.GetId());

  if( running ) {
    pipe.Stop();
    wxQueueEvent(this,event.Clone());

    return true;
  }
  else {

    // to clear of map
    MuniAstrometryEvent ev(EVT_ASTROMETRY);
    ev.astrometry = false;
    wxQueueEvent(GetParent(),ev.Clone());

    // This is important for modeless dialogs which must notify parents
    // (else the close event is stopped here). And potentialy pass outputs.
    wxCommandEvent e(EVT_FINISH_DIALOG,this->GetId());
    e.SetString(file);
    e.SetInt(event.GetId());
    wxQueueEvent(GetParent(),e.Clone());

    return false;
  }
}


void MuniAstrometry::OnStdButton(wxCommandEvent& event)
{
  /*
  if( edited && event.GetId() == wxID_APPLY ) {
    wxQueueEvent(savebutt,new wxCommandEvent());
    wxQueueEvent(this,event.Clone());
    return;
  }
  */

  if( unsaved && event.GetId() == wxID_APPLY ) {
    unsaved = false;
    backupfile = file+config->backup_suffix;
    wxRenameFile(file,backupfile);
    wxRenameFile(workingfile,file);
    workingfile.Clear();
  }

  if( FinishClean(event) )
    return;

  SetReturnCode(event.GetId());
  event.Skip();
}

void MuniAstrometry::CreateControls()
{
  MuniAstrolog astrolog;

  wxSizerFlags sf, lf, cf, sl, rl;
  sf.Align(wxALIGN_CENTER_VERTICAL);
  lf.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT);
  cf.Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT);
  rl.Align(wxALIGN_CENTER_VERTICAL).Border(wxRIGHT);
  sl.Border().Expand();

  wxStaticText *label;

  wxFont bf(*wxNORMAL_FONT);
  bf.SetWeight(wxFONTWEIGHT_BOLD);

  // title
  wxBoxSizer *topsizer = new wxBoxSizer(wxVERTICAL);


  wxNotebook *book = new wxNotebook(this,wxID_ANY);

  wxPanel* manpanel = new wxPanel(book);
  wxBoxSizer *mansizer = new wxBoxSizer(wxVERTICAL);

  wxStaticBoxSizer *psizer = new wxStaticBoxSizer(wxVERTICAL,manpanel,"Projection");

  wxArrayString projections;
  projections.Add("Identity");
  projections.Add("Gnomonic");
  wxChoice *chproj = new wxChoice(manpanel,wxID_ANY,wxDefaultPosition,wxDefaultSize,projections);
  chproj->SetSelection(chproj->FindString(proj));
  psizer->Add(chproj,cf.Center());
  mansizer->Add(psizer,wxSizerFlags().Expand().Border());

  wxStaticBoxSizer *tsizer = new wxStaticBoxSizer(wxVERTICAL,manpanel,"Transformation");


  wxFlexGridSizer *grid = new wxFlexGridSizer(2);
  grid->AddGrowableCol(1);

  label = new wxStaticText(manpanel,wxID_ANY,"Scale:");
  grid->Add(label,lf);

  wscale = new wxSpinCtrlDouble(manpanel,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
				wxSP_ARROW_KEYS|wxTE_PROCESS_ENTER,0.0,100.0,1.0,0.01);
  sunit = new wxButton(manpanel,wxID_ANY,ARC_SEC,wxDefaultPosition,wxDefaultSize,wxBU_EXACTFIT);

  wxBoxSizer *ss = new wxBoxSizer(wxHORIZONTAL);
  ss->Add(wscale);
  ss->Add(sunit);
  grid->Add(ss,cf.Left());

  label = new wxStaticText(manpanel,wxID_ANY,"Rotation:");
  grid->Add(label,lf);
  ss = new wxBoxSizer(wxHORIZONTAL);
  wangle = new wxSpinCtrlDouble(manpanel,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
				wxSP_ARROW_KEYS|wxSP_WRAP|wxTE_PROCESS_ENTER,-360.0,360.0,0.0,0.1);
  ss->Add(wangle);
  ss->Add(new wxStaticText(manpanel,wxID_ANY,L"[°]"),cf);
  grid->Add(ss,cf);

  label = new wxStaticText(manpanel,wxID_ANY,L"Center in α:");
  grid->Add(label,lf);
  ss = new wxBoxSizer(wxHORIZONTAL);
  acenter = new wxSpinCtrlDouble(manpanel,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
				wxSP_ARROW_KEYS|wxTE_PROCESS_ENTER,0.0,360.0,0.0,0.001);
  ss->Add(acenter);
  ss->Add(new wxStaticText(manpanel,wxID_ANY,L"[°]"),cf);
  grid->Add(ss,cf);

  label = new wxStaticText(manpanel,wxID_ANY,L"Center in δ:");
  grid->Add(label,lf);
  ss = new wxBoxSizer(wxHORIZONTAL);
  dcenter = new wxSpinCtrlDouble(manpanel,wxID_ANY,"",wxDefaultPosition,wxDefaultSize,
				wxSP_ARROW_KEYS|wxTE_PROCESS_ENTER,-90.0,90.0,0.0,0.001);
  ss->Add(dcenter);
  ss->Add(new wxStaticText(manpanel,wxID_ANY,L"[°]"),cf);
  grid->Add(ss,cf);

  ss = new wxBoxSizer(wxHORIZONTAL);
  reflex_checkbox = new wxCheckBox(manpanel,wxID_ANY,"Reflex");
  ss->Add(reflex_checkbox);
  grid->Add(ss,cf);

  savebutt = new wxButton(manpanel,wxID_ANY,"Save to header");
  savebutt->Enable(false);

  wxBoxSizer *st = new wxBoxSizer(wxHORIZONTAL);
  st->Add(grid,wxSizerFlags().Border().Center());
  st->Add(savebutt,wxSizerFlags().Border().Align(wxALIGN_RIGHT|wxALIGN_BOTTOM));

  tsizer->Add(st,wxSizerFlags().Border().Right());

  mansizer->Add(tsizer,wxSizerFlags().Expand().Border());

  // accurate

  overlayid = new wxStaticText(manpanel,wxID_ANY,"Overlay: "+astrolog.GetSign());
  overlay_check = new wxCheckBox(manpanel,wxID_ANY,"Draw");
  overlay_check->SetValue(draw_overlay);
  overlay_check->Enable(false);
  //  init_params->SetValue(init_ref);
  //  init_params->SetToolTip("Parameters are initialised from Reference (averaged coordinates as its center, scale as the second moment).");
  wxButton *catbutt = new wxButton(manpanel,wxID_ANY,L"Specify…");

  wxBoxSizer *oversizer = new wxBoxSizer(wxHORIZONTAL);
  oversizer->Add(overlayid,cf.Border());
  oversizer->AddStretchSpacer();
  oversizer->Add(overlay_check,lf);
  oversizer->Add(catbutt,lf.Border());
  mansizer->Add(oversizer,wxSizerFlags().Expand());

  manpanel->SetSizer(mansizer);

  book->AddPage(manpanel,"Edit");


  wxPanel* autopanel = new wxPanel(book);
  autosizer = new wxBoxSizer(wxVERTICAL);

  wxBoxSizer *rhsizer = new wxBoxSizer(wxHORIZONTAL);
  wxStaticText *rt = new wxStaticText(autopanel,wxID_ANY,"Reference:");
  rt->SetLabelMarkup("<small>Reference:</small>");
  rhsizer->Add(rt,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT));

  refcatid = new wxStaticText(autopanel,wxID_ANY,astrolog.GetSign());
  rhsizer->Add(refcatid,wxSizerFlags(1).Border().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT));

  wxButton *refbutt = new wxButton(autopanel,wxID_ANY,L"Specify…");
  rhsizer->Add(refbutt,wxSizerFlags().Border().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  autosizer->Add(rhsizer,wxSizerFlags().Expand());

  wxBoxSizer *msizer = new wxBoxSizer(wxHORIZONTAL);
  infolabel = new wxStaticText(autopanel,wxID_ANY,"");
  info = new wxStaticText(autopanel,wxID_ANY,"");
  msizer->Add(infolabel,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT).Border(wxRIGHT));
  msizer->Add(info,wxSizerFlags(1).Align(wxALIGN_LEFT));
  autosizer->Add(msizer,wxSizerFlags().Expand().Border());


  wxBoxSizer *fsizer = new wxBoxSizer(wxHORIZONTAL);
  proglabel = new wxStaticText(autopanel,wxID_ANY,"");
  fsizer->Add(proglabel,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT).Border(wxRIGHT));
  proginfo = new wxStaticText(autopanel,wxID_ANY,"");
  fsizer->Add(proginfo,wxSizerFlags(1).Align(wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT));
  autosizer->Add(fsizer,wxSizerFlags().Expand().Border());

  wxBoxSizer *prsizer = new wxBoxSizer(wxHORIZONTAL);

  progress = new wxGauge(autopanel,wxID_ANY,MAXGRANGE);
  prsizer->Add(progress,wxSizerFlags(1).Expand().Border());

  calbutt = new wxButton(autopanel,wxID_EXECUTE);
  stopbutt = new wxButton(autopanel,wxID_STOP);
  stopbutt->SetToolTip("Please be patient.");
  stopbutt->Hide();
  prsizer->Add(calbutt,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  prsizer->Add(stopbutt,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT));
  autosizer->Add(prsizer,wxSizerFlags().Expand().Border());


  autosizer->Add(new wxStaticLine(autopanel,wxID_ANY),sl);

  astropt = new MuniAstrometryOptions(autopanel,config);
  autosizer->Add(astropt,wxSizerFlags().Center().Border());

  autopanel->SetSizer(autosizer);

  book->AddPage(autopanel,"Calibrate");


  wxPanel* rempanel = new wxPanel(book);
  wxBoxSizer *remsizer = new wxBoxSizer(wxVERTICAL);

  remsizer->AddStretchSpacer();

  remsizer->Add(new wxStaticText(rempanel,wxID_ANY,"Do you want to remove the astrometry calibration from this frame ?"),wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER).DoubleBorder());

  rembutt = new wxButton(rempanel,wxID_ANY,"Remove");
  remsizer->Add(rembutt,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER));
  rembutt->Enable(false);

  remsizer->AddStretchSpacer();

  rempanel->SetSizer(remsizer);

  book->AddPage(rempanel,"Remove");

  topsizer->Add(book,wxSizerFlags(1).Expand().DoubleBorder());

  wxStdDialogButtonSizer *buttsize =
    static_cast<wxStdDialogButtonSizer *>(CreateButtonSizer(wxAPPLY|wxCANCEL));

  if( buttsize ) {
    topsizer->Add(buttsize,wxSizerFlags().Right().Border(wxBOTTOM));
    SetAffirmativeId(wxID_APPLY);
  }

  SetSizerAndFit(topsizer);

  Bind(wxEVT_UPDATE_UI,&MuniAstrometry::OnUpdateUI,this);
  Bind(wxEVT_IDLE,&MuniAstrometry::OnIdle,this);
  Bind(wxEVT_COMMAND_NOTEBOOK_PAGE_CHANGED,&MuniAstrometry::OnBookChange,this,book->GetId());

  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometry::OnStdButton,this,wxID_APPLY);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometry::OnStdButton,this,wxID_CANCEL);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometry::OnOverlay,this,catbutt->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometry::OnReference,this,refbutt->GetId());
  Bind(wxEVT_COMMAND_CHECKBOX_CLICKED,&MuniAstrometry::OnDrawOverlay,this,overlay_check->GetId());

  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometry::OnCalibrate,this,wxID_EXECUTE);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometry::OnCalibrateStop,this,wxID_STOP);
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometry::OnRemove,this,rembutt->GetId());
  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometry::OnSaveWCS,this,savebutt->GetId());

  Bind(wxEVT_COMMAND_BUTTON_CLICKED,&MuniAstrometry::OnPopScaleUnit,this,sunit->GetId());
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniAstrometry::OnScaleUnit,this,ID_SCALE_ARCUAS);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniAstrometry::OnScaleUnit,this,ID_SCALE_ARCMAS);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniAstrometry::OnScaleUnit,this,ID_SCALE_ARCSEC);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniAstrometry::OnScaleUnit,this,ID_SCALE_ARCMIN);
  Bind(wxEVT_COMMAND_MENU_SELECTED,&MuniAstrometry::OnScaleUnit,this,ID_SCALE_ARCDEG);

  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniAstrometry::OnTextDouble,this,wangle->GetId());
  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniAstrometry::OnTextDouble,this,wscale->GetId());
  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniAstrometry::OnTextDouble,this,acenter->GetId());
  Bind(wxEVT_COMMAND_TEXT_UPDATED,&MuniAstrometry::OnTextDouble,this,dcenter->GetId());
  Bind(wxEVT_COMMAND_CHECKBOX_CLICKED,&MuniAstrometry::OnReflex,this,
       reflex_checkbox->GetId());

  Bind(wxEVT_COMMAND_TEXT_ENTER,&MuniAstrometry::OnTextDouble,this,wangle->GetId());
  Bind(wxEVT_COMMAND_TEXT_ENTER,&MuniAstrometry::OnTextDouble,this,wscale->GetId());
  Bind(wxEVT_COMMAND_TEXT_ENTER,&MuniAstrometry::OnTextDouble,this,acenter->GetId());
  Bind(wxEVT_COMMAND_TEXT_ENTER,&MuniAstrometry::OnTextDouble,this,dcenter->GetId());

  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniAstrometry::OnSpinDouble,this,wscale->GetId());
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniAstrometry::OnSpinDouble,this,wangle->GetId());
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniAstrometry::OnSpinDouble,this,acenter->GetId());
  Bind(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED,&MuniAstrometry::OnSpinDouble,this,dcenter->GetId());
  Bind(wxEVT_COMMAND_CHOICE_SELECTED,&MuniAstrometry::OnChoiceProj,this,chproj->GetId());

  ids.push_back(chproj->GetId());
  ids.push_back(wscale->GetId());
  ids.push_back(sunit->GetId());
  ids.push_back(wangle->GetId());
  ids.push_back(acenter->GetId());
  ids.push_back(dcenter->GetId());
  ids.push_back(savebutt->GetId());
  ids.push_back(overlay_check->GetId());
  ids.push_back(catbutt->GetId());
  ids.push_back(refbutt->GetId());
  ids.push_back(rembutt->GetId());
  ids.push_back(reflex_checkbox->GetId());
}

void MuniAstrometry::OnUpdateUI(wxUpdateUIEvent& e)
{
  //  wxLogDebug("MuniAstrometry::OnUpdateUI");

  for(list<int>::const_iterator i = ids.begin(); i != ids.end(); ++i)
    FindWindow(*i)->Enable(!running);

  wxASSERT(astropt);
  astropt->Enable(!running);

  savebutt->Enable(edited);

  wxASSERT(FindWindow(wxID_APPLY));
  FindWindow(wxID_APPLY)->Enable(unsaved && !readonly);

  wxASSERT(FindWindow(wxID_CANCEL));
  FindWindow(wxID_CANCEL)->Enable(!running);

  if( init_file ) {
    init_file = false;

    Reset();

    if( array.GetWCS(xoff,yoff,alpha,delta,scale,angle,reflex) ) {

      wscale->SetValue(Scale(scale));
      //      wangle->SetValue(Period(angle));
      wangle->SetValue(angle);
      acenter->SetValue(alpha);
      dcenter->SetValue(delta);
      reflex_checkbox->SetValue(reflex < 0);

      init_par = true;

      rembutt->Enable();

      DrawOverlay(coverlay);
    }
  }

  if( stars.IsOk() && (! catfile.IsEmpty() || ! reffile.IsEmpty()) ) {
    calbutt->Enable();
  }


}

void MuniAstrometry::OnIdle(wxIdleEvent&)
{
  if( showtooltip ) {
    showtooltip = false;

    wxASSERT(FindWindow(wxID_APPLY));
    if( readonly ) {
      wxRichToolTip tip("Read Only Data",
			"This file is not writable.\n"
			"Any results will not saved.");
      tip.SetIcon(wxICON_WARNING);
      tip.ShowFor(FindWindow(wxID_APPLY));
    }
  }

}

void MuniAstrometry::OnChoiceProj(wxCommandEvent& event)
{
  proj = event.GetString();
}

void MuniAstrometry::OnPopScaleUnit(wxCommandEvent& event)
{
  wxMenu popup;
  popup.AppendRadioItem(ID_SCALE_ARCMAS,ARC_UAS);
  popup.AppendRadioItem(ID_SCALE_ARCMAS,ARC_MAS);
  popup.AppendRadioItem(ID_SCALE_ARCSEC,ARC_SEC);
  popup.AppendRadioItem(ID_SCALE_ARCMIN,ARC_MIN);
  popup.AppendRadioItem(ID_SCALE_ARCDEG,ARC_DEG);

  wxString l(sunit->GetLabel());
  if( l == ARC_UAS )
    popup.Check(ID_SCALE_ARCUAS,true);
  else if( l == ARC_MAS )
    popup.Check(ID_SCALE_ARCMAS,true);
  else if( l == ARC_SEC )
    popup.Check(ID_SCALE_ARCSEC,true);
  else if( l == ARC_MIN )
    popup.Check(ID_SCALE_ARCMIN,true);
  else if( l == ARC_DEG )
    popup.Check(ID_SCALE_ARCDEG,true);

  PopupMenu(&popup);
}

void MuniAstrometry::OnScaleUnit(wxCommandEvent& event)
{
  switch(event.GetId()){
  case ID_SCALE_ARCUAS: sunit->SetLabel(ARC_UAS); wscale->SetValue(3600.e6*scale); return;
  case ID_SCALE_ARCMAS: sunit->SetLabel(ARC_MAS); wscale->SetValue(3600.e3*scale); return;
  case ID_SCALE_ARCSEC: sunit->SetLabel(ARC_SEC); wscale->SetValue(3600.*scale); return;
  case ID_SCALE_ARCMIN: sunit->SetLabel(ARC_MIN); wscale->SetValue(60.*scale);   return;
  case ID_SCALE_ARCDEG: sunit->SetLabel(ARC_DEG); wscale->SetValue(scale);       return;
  }
}

void MuniAstrometry::OnInitRef(wxCommandEvent& event)
{
  init_ref = event.IsChecked();
}

void MuniAstrometry::SetFile(const wxString& f, const FitsArray& a)
{
  file = f;
  array = a;

  init_file = true;

  readonly = ! wxFileName::IsFileWritable(file);

  workingfile= wxFileName::CreateTempFileName("xmunipack-astrometry_");
  wxRemoveFile(workingfile);
  FitsCopyFile(file,workingfile);
}

void MuniAstrometry::SetDetectedSources(const FitsTable& t)
{
  stars = t;
}


wxString MuniAstrometry::GetBackup() const
{
  return backupfile;
}

void MuniAstrometry::Reset()
{
  wscale->SetValue("1.0");
  wangle->SetValue("0.0");
  acenter->SetValue("0.0");
  dcenter->SetValue("0.0");
  reflex_checkbox->SetValue(false);

  if( array.IsOk() ) {
    wxASSERT(array.Naxis() == 2);
    xoff = array.Naxes(0)/2.0;
    yoff = array.Naxes(1)/2.0;
  }
}

FitsTable MuniAstrometry::LoadCatalogue(const wxString& filename)
{

  FitsTable table;

  FitsFile fits(filename);

  if( fits.IsOk() ) {
    for(size_t i = 0; i < fits.HduCount(); i++) {

      if( ! fits.IsOk() || fits.Hdu(i).Type() != HDU_TABLE ) continue;

      table = fits.Hdu(i);
      wxLogDebug("%d stars in list.",(int)table.Nrows());
    }
  }

  return table;
}


void MuniAstrometry::InitByCatalogue(const FitsTable& catalogue)
{
  // finds coordinate keys
  int ira = catalogue.GetColIndex(label_ra);
  int idec = catalogue.GetColIndex(label_dec);

  if( ! init_par && (ira >= 0 && idec >= 0) ) {

    double xx,yy,sx,sy;
    xx = 0; yy = 0; sx = 0; sy = 0;

    int n = catalogue.Nrows();
    const double *a = catalogue.GetColumn(label_ra).GetCol_double();
    const double *d = catalogue.GetColumn(label_dec).GetCol_double();


    for(int i = 0; i < n; i++) {

      double x = a[i];
      double y = d[i];
      xx = xx + x;
      yy = yy + y;
      sx = sx + x*x;
      sy = sy + y*y;

    }

    if( n > 0 ) {

      alpha = xx/n;
      delta = yy/n;

      acenter->SetValue(alpha);
      dcenter->SetValue(delta);

      sx = sqrt(sx/n);
      //	sy = sqrt(sy/n);
      //	wxLogDebug("%f %d %f",sy,n,delta);
      sy = sqrt(sy/n - delta*delta);

      double r = sy; // in declination, there is no scale contraction

      if( yoff ) {
	scale = yoff/r;
	//	  wxLogDebug("%f %f %f",(float)alpha,(float)delta,(float)scale);
	init_par = true;
      }
      return;
      }

  }

  if( init_par ) {

    wscale->SetValue(Scale(scale));
  }
}

void MuniAstrometry::OnSpinDouble(wxSpinDoubleEvent& event)
{

  if( event.GetId() == wangle->GetId() )
    angle = event.GetValue();
  else if( event.GetId() == acenter->GetId() )
    alpha = event.GetValue();
  else if( event.GetId() == dcenter->GetId() )
    delta = event.GetValue();
  else if( event.GetId() == wscale->GetId() ) {
    scale = Scale(event.GetValue());
  }

  DrawOverlay(coverlay);
  edited = true;
}

void MuniAstrometry::DrawOverlay(const FitsTable& t)
{
  if( draw_overlay && t.IsOk() ) {

    MuniAstrometryEvent ev(EVT_ASTROMETRY);
    ev.astrometry = true;
    ev.proj = proj;
    ev.scale = scale;
    ev.angle = angle;
    ev.reflex = reflex;
    ev.xcen = xoff;
    ev.ycen = yoff;
    ev.acen = alpha;
    ev.dcen = delta;
    MuniAstrometryLayer l(array.Naxes(0),array.Naxes(1));
    l.SetProj(proj,alpha,delta,xoff,yoff,scale,angle,reflex);
    l.DrawStars(t,label_ra,label_dec,label_pmra,label_pmdec,label_mag);
    if( page == 0 )
     l.DrawCompas();
    ev.layer = l.GetLayer();
    wxQueueEvent(GetParent(),ev.Clone());
  }
  else {
    MuniAstrometryEvent ev(EVT_ASTROMETRY);
    ev.astrometry = false;
    wxQueueEvent(GetParent(),ev.Clone());
  }

}

void MuniAstrometry::OnTextDouble(wxCommandEvent& event)
{
  double x;

  if( ! event.GetString().ToDouble(&x) ) return;

  if( event.GetId() == wangle->GetId() )
    angle = x;
  else if( event.GetId() == acenter->GetId() )
    alpha = x;
  else if( event.GetId() == dcenter->GetId() )
    delta = x;
  else if( event.GetId() == wscale->GetId() ) {
    scale = Scale(x);
  }

  DrawOverlay(coverlay);
  edited = true;
}

void MuniAstrometry::OnReflex(wxCommandEvent& event)
{
  reflex = event.IsChecked() ? -1 : 1;
  DrawOverlay(coverlay);
  edited = true;
}

void MuniAstrometry::RunProcessing()
{
  wxASSERT(astropt);

  MuniProcess *astrometry = new MuniProcess(&pipe,"astrometry");
  pipe.push(astrometry);

  astrometry->Write("PIPELOG = T");
  astrometry->Write("PROJECTION = '"+proj.Upper()+"'");
  astrometry->Write("WCSSAVE = T");
  astrometry->Write("AUNITS = '"+astropt->GetOutputUnits()+"'");
  astrometry->Write("SIG = %e",astropt->GetSig());
  astrometry->Write("FSIG = %e",astropt->GetFSig());

  if( astropt->GetMatchType() == ID_ASTRO_NEAR ) {
    wxString a;
    astrometry->Write("MATCH = 'NEARLY'");
    astrometry->Write("INITPAR = F");
    astrometry->Write("SCALE = %25.15e",scale);
    astrometry->Write("ANGLE = %20.15f",angle);
    a.Printf("CRVAL = %20.15f %20.15f",alpha,delta);
    astrometry->Write(a);
  }
  else if( astropt->GetMatchType() == ID_ASTRO_MATCH ) {
    astrometry->Write("MATCH = 'BACKTRACKING'");
    astrometry->Write("MINMATCH = %ld",long(astropt->GetMinMatch()));
    astrometry->Write("MAXMATCH = %ld",long(astropt->GetMaxMatch()));
    if( astropt->GetFullMatch() ) astrometry->Write("FULLMATCH = T");
  }

  if( label_ra != "" && label_dec != "" ) {
    astrometry->Write("COL_RA = '" + label_ra + "'");
    astrometry->Write("COL_DEC = '" + label_dec + "'");
    if( label_pmra != "" )  astrometry->Write("COL_PMRA = '" + label_pmra + "'");
    if( label_pmdec != "" ) astrometry->Write("COL_PMDEC = '" + label_pmdec + "'");
    if( label_mag != "" ) astrometry->Write("COL_MAG = '" + label_mag + "'");
    astrometry->Write("CAT = '"+catfile+"'");
  }
  else {
    if( relative )
      astrometry->Write("REL = '"+reffile+"'");
    else
      astrometry->Write("REF = '"+reffile+"'");
  }

  astrometry->Write("FILE = '"+workingfile+"' ''");
}

void MuniAstrometry::OnCalibrate(wxCommandEvent& event)
{
  Bind(wxEVT_END_PROCESS,&MuniAstrometry::OnCalibrateFinish,this);
  Bind(wxEVT_TIMER,&MuniAstrometry::OnTimer,this);

  wxBeginBusyCursor();

  calbutt->Hide();
  stopbutt->Show();
  info->SetLabel("");
  infolabel->SetLabel("");
  proglabel->SetLabelMarkup("<b>Status:</b>");
  proginfo->SetLabel("Astrometry is being started ...");
  autosizer->Layout();

  RunProcessing();

  running = true;
  parsing = false;
  output_index = 0;
  nhist = 0;
  timer.Start(250);
  start.SetToCurrent();
  pipe.Start();
}

void MuniAstrometry::OnCalibrateStop(wxCommandEvent& event)
{
  pipe.Stop();
}

void MuniAstrometry::OnCalibrateFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniAstrometry::OnFinish");

  ParseProcessing(pipe.GetOutput());

  timer.Stop();
  stopbutt->Hide();
  calbutt->Show();
  running = false;
  wxEndBusyCursor();

  Unbind(wxEVT_END_PROCESS,&MuniAstrometry::OnCalibrateFinish,this);
  Unbind(wxEVT_TIMER,&MuniAstrometry::OnTimer,this);

  // check status
  if( event.GetExitCode() == 0 )
    unsaved = true;
  else {
    infolabel->SetLabel("");
    info->SetLabel("");
    proglabel->SetLabel("Status:");
    proginfo->SetLabelMarkup("<i>Interrupted on request.</i>");
    progress->SetValue(0);
  }
  autosizer->Layout();
}


// parse output of astrometry fit and fill the data table
void MuniAstrometry::OnTimer(wxTimerEvent& event)
{
  ParseProcessing(pipe.GetOutput());
}

void MuniAstrometry::ParseProcessing(const wxArrayString& out)
{
  //  wxLogDebug("MuniAstrometry::ParseProcessing");
  wxASSERT(astropt);

  if( parsing ) return;
  parsing = true;

  wxRegEx re("^=(.*)> (.+)");
  wxASSERT(re.IsValid());

  bool afit = false;
  bool match = false;

  for(size_t idx = output_index; idx < out.GetCount(); idx++) {

    if( re.Matches(out[idx]) ) {

      wxString key(re.GetMatch(out[idx],1));
      wxString value(re.GetMatch(out[idx],2));

      if( key == "ASTROMETRY" ) {

	if( value.Find("Start") != wxNOT_FOUND ) {
	  progress->SetValue(0);
	  proginfo->SetLabel("Starting...");

	  nhist = 0;
	  hist.clear();
	  u.clear(); v.clear(), x.clear(); y.clear();
	}
	else if( value.Find("Finish") != wxNOT_FOUND ) {
	  progress->SetValue(MAXGRANGE);
	}

      }
      else if( key == "ASTROFIT" ) {

	if( value.Find("Start") != wxNOT_FOUND ) {

	}
	else if( value.Find("Final") != wxNOT_FOUND ) {

	  wxStringInputStream ss(value.Mid(value.Find("Final")+5));
	  wxTextInputStream t(ss);

	  t >> rms;
	}
	else if( value.Find("Finish") != wxNOT_FOUND ) {

	  proglabel->SetLabelMarkup("<small>Astrometry:</small>");

	  if( value.Find("Success") != wxNOT_FOUND ) {

	    // update parameters in Edit
	    wscale->SetValue(Scale(scale));
	    wangle->SetValue(angle);
	    acenter->SetValue(alpha);
	    dcenter->SetValue(delta);
	    //	    reflex_checkbox->SetValue();

	    wxString a;
	    a.Printf(" Success (RMS = %.2f "+astropt->GetOutputUnits()+")",
		     Scale(1.0/rms));
	    proginfo->SetLabel(a);

	  }
	  else if( value.Find("Fail") != wxNOT_FOUND ) {
	    proginfo->SetLabelMarkup("<i>Failed</i>");
	  }
	}
      }

      else if( key == "ASTROMATCH" ) {

	if( value.Find("Start") != wxNOT_FOUND ) {
	  proglabel->SetLabelMarkup("<b>Search:</b>");
	  infolabel->SetLabelMarkup("<small>Match:</small>");
	}
	else if( value.Find("Finish") != wxNOT_FOUND ) {

	  proglabel->SetLabel("");
	  proginfo->SetLabel("");

	  if( value.Find("Success") != wxNOT_FOUND ) {
	  }
	  else if( value.Find("Fail") != wxNOT_FOUND ) {
	    info->SetLabelMarkup("<i>Not found</i>");
	  }

	}
      }
      else if( key == "WCSAVE" ) {

	if( value.Find("Start") != wxNOT_FOUND ) {
	  proglabel->SetLabel("WCSSave:");
	  proginfo->SetLabel("Saving ...");
	}
	else if( value.Find("Finish") != wxNOT_FOUND ) {

	  infolabel->SetLabel("WCSSave:");
	  if( value.Find("Success") != wxNOT_FOUND ) {
	    proginfo->SetLabel("Success");
	  }
	  else if( value.Find("Fail") != wxNOT_FOUND ) {
	    proginfo->SetLabel("Failed.");
	  }

	}

      }
      else if( key == "AFIT" || key == "RFIT" || key == "MFIT" ) {

	afit = true;

	wxStringInputStream ss(value);
	wxTextInputStream t(ss);

	double xcen,ycen;
	t >> s0 >> scale >> angle >> alpha >> delta >> xcen >> ycen;
      }


      else if( key == "MPROGRESO" ) {

	wxStringInputStream ss(value);
	wxTextInputStream t(ss);

	int ns, ntot;
	t >> ns >> ntot;

	double r = double(ns)/double(ntot);

	wxDateTime now;
	now.SetToCurrent();
	wxTimeSpan d = now.Subtract(start);
	wxLongLong sec = d.GetSeconds();

	wxString ts;

	//	wxLogDebug("%f %f",r,sec.ToDouble());
	// print in human time units
	if( r > 0.0 && sec > 1 ) {

	  double vel = r/sec.ToDouble(); // velocity in fractions per second
	  wxLongLong est;
	  est.Assign((1 - r)/vel);

	  wxTimeSpan dt = wxTimeSpan::Seconds(est);
	  int d = dt.GetDays();
	  int h = dt.GetHours();
	  int m = dt.GetMinutes();
	  wxLongLong s = dt.GetSeconds();

	  //	  wxLogDebug("%d %d %d "+s.ToString(),d,h,m);

	  wxString ss;
	  if( d > 0 ) {
	    if( d <= 2 ) {
	      ss.Printf("%d day %d hour",d,h-24*d);
	    }
	    else
	      ss.Printf("%d days",d);
	  }
	  else if( h > 0 ) {
	    if( h <= 2 ) {
	      int mm = m - 60*h;
	      if( mm > 1 ) mm = 10 * (mm / 10);
	      ss.Printf("%d hour %d min",h,mm);
	    }
	    else
	      ss.Printf("%d hours",h);
	  }
	  else if( m > 0 ) {
	    if( m <= 2 ) {
	      int sec = int(s.ToLong())-60*m;
	      if( sec > 15 ) sec = 10 * (sec / 10);
	      ss.Printf("%d min %d sec",m,sec);
	    }
	    else
	      ss.Printf("%d min",m);
	  }
	  else {
	    long sec = s.ToLong();
	    if( sec > 5 ) sec = 5 *(sec / 5);
	    ss.Printf("%d sec",int(sec));
	  }

	  ts = "(finish up to " + ss + ")";
	}
	else
	  ts = "";


	wxString a;
	a.Printf(" %d%% completed "+ts,int(100.0*r));

	proginfo->SetLabel(a);
	progress->SetValue(wxMin(int(MAXGRANGE*r),MAXGRANGE));

      }

      else if( key == "MPROGRES2" ) {

	match = true;

	wxStringInputStream ss(value);
	wxTextInputStream t(ss);

	int n;
	t >> amin >> dmin >> fmin >> n;

	x.clear(); y.clear();
	if( n > 0 ) {
	  x.resize(n); y.resize(n);
	  for(int i = 0; i < n; i++)
	    t >> x[i] >> y[i];
	}

	int m;
	t >> m;
	u.clear(); v.clear();
	if( m > 0 ) {
	  u.resize(m); v.resize(m);
	  for(int i = 0; i < m; i++)
	    t >> u[i] >> v[i];
	}

	wxString b;
	b.Printf(L" Sequence %d stars (Δd = %.2f‰, Δf = %.1f, Δφ = %.2f°)",
		 n,1000.0*dmin,fmin,57.3*amin);
	/*
	if( amin < probability )
	  b.Printf(L"p=%.3f σ=%.5f (Excellent)",pmin,dmin);
	else {
	  wxString c(L"∞");
	  if( dmin > DBL_EPSILON ) c.Printf("%.5g",dmin);
	  b.Printf(L"p=%.3f σ="+c+" (Poor)",pmin);
	}
	*/
	info->SetLabel(b);
      }

      else if( key == "MTRI" ) {

	wxStringInputStream ss(value);
	wxTextInputStream t(ss);

	t >> nhist;
	hist.clear();
	hist.resize(nhist*nhist);
	for(int i = 0; i < nhist*nhist && ss.IsOk(); i++)
	  t >> hist[i];
      }

    }
  }
  if( out.GetCount() > 0 )
    output_index = out.GetCount() - 1;

  if( nhist > 0 && (match || afit) ) {

    MuniAstrometryLayer layer(array.Naxes(0),array.Naxes(1));
    MuniAstrometryEvent ev(EVT_ASTROMETRY);
    ev.astrometry = true;

    layer.Draw3Space(nhist,hist);
    if( afit ) {
      ev.proj = proj;
      ev.scale = scale;
      ev.angle = angle;
      ev.xcen = xoff;
      ev.ycen = yoff;
      ev.acen = alpha;
      ev.dcen = delta;
    }
    layer.SetProj(proj,alpha,delta,xoff,yoff,scale,angle,reflex);
    //    layer.DrawStars(catalogue,label_ra,label_dec,label_pmra,label_pmdec,label_mag);

    layer.Draw3Tri(u,v);

    if( x.size() > 0 ) {
      layer.DrawMatch(astropt->GetSig(),x,y);
    }
    if( ! layer.IsEmpty() ) {
      ev.layer = layer.GetLayer();
      wxQueueEvent(GetParent(),ev.Clone());
    }
  }

  autosizer->Layout();
  parsing = false;
}


void MuniAstrometry::OnSaveWCS(wxCommandEvent& event)
{
  if( running ) { wxFAIL_MSG("Process in progress?"); }


  Bind(wxEVT_END_PROCESS,&MuniAstrometry::OnSaveWCSFinish,this);

  wxBeginBusyCursor();
  running = true;

  MuniProcess *action = new MuniProcess(&pipe,"wcsupdate");
  pipe.push(action);

  action->Write("WCSSAVE = T");
  action->Write("PROJECTION = '"+proj.Upper()+"'");

  wxString a;
  a.Printf("CRVAL = %20.15f %20.15f",alpha,delta);
  action->Write(a);
  action->Write("SCALE = %25.15e",scale);
  action->Write("ANGLE = %20.15f",angle);
  action->Write("FILE = '"+workingfile+"' ''");

  pipe.Start();
}

void MuniAstrometry::OnSaveWCSFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniAstrometry::OnSaveFinish");

  Unbind(wxEVT_END_PROCESS,&MuniAstrometry::OnSaveWCSFinish,this);
  wxEndBusyCursor();
  running = false;

  if( event.GetExitCode() == 0 ) {
    unsaved = true;
  }
  edited = false;
}


void MuniAstrometry::OnRemove(wxCommandEvent& event)
{
  if( running ) { wxFAIL_MSG("Process in progress?"); }

  Bind(wxEVT_END_PROCESS,&MuniAstrometry::OnRemoveFinish,this);

  wxBeginBusyCursor();
  running = true;

  MuniProcess *action = new MuniProcess(&pipe,"wcsremove");
  pipe.push(action);

  action->Write("WCSSAVE = T");
  action->Write("FILE= '"+workingfile+"' ''");

  pipe.Start();
}

void MuniAstrometry::OnRemoveFinish(wxProcessEvent& event)
{
  wxLogDebug("MuniAstrometry::OnRemoveFinish");

  Unbind(wxEVT_END_PROCESS,&MuniAstrometry::OnRemoveFinish,this);
  wxEndBusyCursor();
  running = false;
  unsaved = true;
}

void MuniAstrometry::RemoveWorkingFile()
{
  if( ! workingfile.IsEmpty() && wxFileExists(workingfile) ) {
    wxRemoveFile(workingfile);
    workingfile.Clear();
  }
}


double MuniAstrometry::Scale(double s) const
{
  if( sunit->GetLabel() == ARC_UAS )
    return 3600.0e6/s;
  else if( sunit->GetLabel() == ARC_MAS )
    return 3600.0e3/s;
  else if( sunit->GetLabel() == ARC_SEC )
    return 3600.0/s;
  else if( sunit->GetLabel() == ARC_MIN )
    return 60.0/s;
  else
    return s;
}

double MuniAstrometry::Period(double s) const
{
  if( s > 360.0 ) {
    int n = int(s/360.0);
    s = s - n*360.0;
  }
  else if ( s < 0.0 ) {
    int n = int(s/360.0);
    s = s + n*360.0;
  }
  return s;
}

void MuniAstrometry::OnOverlay(wxCommandEvent& event)
{
  MuniSelectSource ms(this,config);
  if( ms.ShowModal() == wxID_OK ) {
    wxString coverfile = ms.GetPath();
    label_ra = ms.GetLabelRA();
    label_dec = ms.GetLabelDec();
    label_pmra = ms.GetLabelPMRA();
    label_pmdec = ms.GetLabelPMDec();
    label_mag = ms.GetLabelMag();
    overlayid->SetLabel(ms.GetId());
    coverlay = LoadCatalogue(coverfile);
    DrawOverlay(coverlay);
    overlay_check->Enable();
    Layout();
  }
}

void MuniAstrometry::OnReference(wxCommandEvent& event)
{
  MuniSelectSource ms(this,config,true);
  if( ms.ShowModal() == wxID_OK ) {

    EraseTemp();

    if( ms.GetType() == ID_ASTRO_CAT ) {

      tmpcat = ms.IsTemporary();
      if( tmpcat ) {
	catfile = wxFileName::CreateTempFileName("xmunipack-astrometry_");
	wxCopyFile(ms.GetPath(),catfile);
      }
      else
	catfile = ms.GetPath();
      label_ra = ms.GetLabelRA();
      label_dec = ms.GetLabelDec();
      label_pmra = ms.GetLabelPMRA();
      label_pmdec = ms.GetLabelPMDec();
      label_mag = ms.GetLabelMag();
      refcatid->SetLabel(ms.GetId());
      catalogue = LoadCatalogue(catfile);
      if( init_par )
	InitByCatalogue(catalogue);
      if( ! coverlay.IsOk() ) {
	overlayid->SetLabel(ms.GetId());
	coverlay = LoadCatalogue(catfile);
	DrawOverlay(coverlay);
	overlay_check->Enable();
      }
      Layout();
    }
    else if( ms.GetType() == ID_ASTRO_REF ) {
      reffile = ms.GetPath();
      relative = ms.GetRelative();
      label_ra.Clear();
      label_dec.Clear();
      label_pmra.Clear();
      label_pmdec.Clear();
      label_mag.Clear();
      refcatid->SetLabel(ms.GetId());
      Layout();
    }
  }
}

void MuniAstrometry::OnDrawOverlay(wxCommandEvent& event)
{
  draw_overlay = event.IsChecked();
  DrawOverlay(coverlay);
}

void MuniAstrometry::OnBookChange(wxBookCtrlEvent& event)
{
  page = event.GetSelection();
  DrawOverlay(coverlay);

  if( page == 1 && ! stars.IsOk() ) {
    wxRichToolTip tip("Photometry unavailable",
		      "Already detected stars are required.");
    tip.SetIcon(wxICON_WARNING);
    tip.ShowFor(calbutt);
  }
}
