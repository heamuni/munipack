/*

  Munipack - the action

  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)

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

#include "mprocess.h"
#include "munipack.h"
#include "fitsfn.h"
#include <wx/wx.h>
#include <wx/cmdline.h>
#include <wx/regex.h>
#include <wx/tokenzr.h>
#include <wx/filename.h>
#include <thread>

// by default, switch off the advance processing
#define ADVANCE false

// the default regexp pattern
#define PATTERN "(.+)\\.(.+)"

MuniAction::MuniAction(wxEvtHandler *h, const wxString& c):
  MuniProcess(h,c), advanced(ADVANCE), pattern(PATTERN) {}

void MuniAction::SetCAF(const MuniCmdLineParser& cmd)
{
  /* coarray wrapper */
  wxArrayString args;
  wxString caf_run, caf_nproc, xproc, a;
  long n;
  int mproc = 3; // by default, 3 processes
  bool caf_in_action = false;

  // * run-time check of CAF_RUN
  // * check for available processors
  // * default processes by MUNIPACK_CAF_MPROC
  // * --Xcaf "args...", OnCmd("args ..")

#ifdef MUNIPACK_CAF_RUN
  caf_run = wxString(MUNIPACK_CAF_RUN);
#endif

  if( wxGetEnv("MUNIPACK_CAF_RUN",&a) ) caf_run = a;
  caf_in_action = caf_run != "";

  if( wxGetEnv("MUNIPACK_CAF_MPROC",&a) && a.ToLong(&n) )
    mproc = n;

  if( cmd.Found("mproc",&n) )
    mproc = n;

  if( wxGetEnv("MUNIPACK_CAF_XPROC",&a) )
    xproc = a;

  if( cmd.Found("Xproc",&a) )
    xproc = a;

  int maxproc = std::thread::hardware_concurrency();
  int nproc = std::max(1,std::min(std::max(1,maxproc),mproc));

  // One process runs the utility directly (built for a single process);
  // more processes run its parallel variant (MUNIPACK_CAF_SUFFIX, built
  // with the MPI coarray library) by the launcher.
  if( caf_in_action && nproc > 1 ) {
    caf_nproc << nproc;

    args.push_back(caf_run);
    args.push_back("-n");
    args.push_back(caf_nproc);

    if( xproc != "" ) {
      wxArrayString cproc =
	wxCmdLineParser::ConvertStringToArgs(xproc,wxCMD_LINE_SPLIT_UNIX);
      for(size_t i = 0; i < cproc.GetCount(); i++)
	args.push_back(cproc[i]);
    }

#ifdef MUNIPACK_CAF_SUFFIX
    args.push_back(GetCommand() + MUNIPACK_CAF_SUFFIX);
#else
    args.push_back(GetCommand());
#endif
  }
  else
    args.push_back(GetCommand());

  size_t argc = args.GetCount();
  wxChar **argv = new wxChar*[argc];
  for(size_t i = 0; i < argc; i++)
    argv[i] = const_cast<wxChar *>(args[i].wx_str());

  SetCmdLine(argc,argv);

  delete[] argv;
}

void MuniAction::CommonOutputSingle(const MuniCmdLineParser& cmd, const wxString& name)
{
  wxString a;

  if( cmd.Found("verbose") )
    Write("VERBOSE = T");

  wxString output;
  if( cmd.Found("output",&a) || cmd.Found("o",&a) )
    output = a;
  else if( name != "" )
    output = name + ".fits";
  else
    output = GetCommand() + ".fits";

  //  WriteOutput(output);
  Write("OUTPUT = '" + output + "' ");
}

void MuniAction::CommonOutputMultiple(const MuniCmdLineParser& cmd)
{
  wxString a;

  if( cmd.Found("verbose") )
    Write("VERBOSE = T");

  if( cmd.Found("target-directory",&a) || cmd.Found("t",&a) ) {
    if( ! SetTargetDir(a) )
      MuniFatalError("Directory `"+target+"' does not exists.");
  }

  if( cmd.Found("O") || cmd.Found("mask",&a) ) {

    SetAdvanced(true);

    if( cmd.Found("pattern",&a) )
      SetPattern(a);

    if( cmd.Found("mask",&a) )
      SetMask(a);

    if( cmd.Found("format",&a) )
      SetFormat(a);
  }
}

void MuniAction::CommonOptionsBitpix(const MuniCmdLineParser& cmd)
{
  long n;
  if( cmd.Found("bitpix",&n) || cmd.Found("B",&n) )
    Write("BITPIX = %ld",n);
}

void MuniAction::CommonOptionsEnableOverwrite(const MuniCmdLineParser& cmd)
{
  if( cmd.Found("enable-overwrite") )
    Write("OVERWRITE = T");
}

void MuniAction::CommonOptionsPhCorr(const MuniCmdLineParser& cmd)
{
  wxString a;
  double x;

  // command line options
  if( cmd.Found("bitmask",&a) )
    Write("MASK = '"+a+"'");

  if( cmd.Found("st",&x) || cmd.Found("saturate",&x) )
    Write("SATURATE = %lf",x);

  if( cmd.Found("th",&x) || cmd.Found("threshold",&x) )
    Write("THRESHOLD = %lf",x);
}

void MuniAction::EnvironmentPhCorr()
{
  // environment variables
  const char *envvar[] = {
    "MUNIPACK_TEMPERATURE_TOLERANCE",
    "MUNIPACK_EXPTIME_TOLERANCE",
    0
  };

  for(size_t i = 0; envvar[i] != 0; i++) {
    wxString val, var(envvar[i]);
    double x;
    if( wxGetEnv(var,&val) && val.ToDouble(&x) )   // for numbers
      Write(var + " = %.1le",x);
  }
}

void MuniAction::CommonOptionsPhq(const MuniCmdLineParser& cmd)
{
  wxString a;

  if( cmd.Found("q",&a) || cmd.Found("quantity",&a) ) {

    wxString l;
    long n;

    apstr(a,n,l);
    Write("NQUANTITIES = %d",n);
    Write("QUANTITIES = " + l);
  }
}


void MuniAction::WriteFilesMultiple(const wxArrayString& files)
{
  Write("NFILES = %d",static_cast<int>(files.GetCount()));

  for(size_t n = 0; n < files.GetCount(); n++) {

    // splits input 'line' by filenames: source [target]
    // filenames with spaces are taken on more pieces than expected.
    wxArrayString filenames;
    wxStringTokenizer tokenizer(files[n]);
    while ( tokenizer.HasMoreTokens() ) {
      wxString name = tokenizer.GetNextToken();
      filenames.push_back(name);
    }

    wxString filename, output;

    if( filenames.size() == 1 ) {
      filename = filenames[0];
      if( advanced ) {
	if( format == "" )
	  output = GetAdvanced(filename);
	else
	  output = GetFormat(n);
      }
      else
	output = GetOutput(filename);
    }
    else if( filenames.size() == 2 ) {
      filename = filenames[0];
      output = filenames[1];
    }
    else
      MuniFatalError("Required the format of input file: source-file [target-file]");

    Write("FILE = '" + filename + "' '" + output + "'");
  }
}


void MuniAction::WriteFilesSingle(const wxArrayString& files)
{
  Write("NFILES = %d",static_cast<int>(files.size()));
  for(size_t n = 0; n < files.size(); n++)
    Write("FILE = '" + files[n] + "'");
}

void MuniAction::apstr(const wxString& a, long& n, wxString& l) const
{
  n = 0;
  l = "";

  wxStringTokenizer tokenizer(a,",");
  while( tokenizer.HasMoreTokens() ) {
    l += " '"+tokenizer.GetNextToken()+"'";
    n++;
  }
}

void MuniAction::dblstr(const wxString& a, const wxString &sep,
		      size_t& n, wxString& l) const
{
  n = 0;
  l = "";

  wxStringTokenizer objtokenizer(a,sep);
  while( objtokenizer.HasMoreTokens() ) {
    wxString c(objtokenizer.GetNextToken());
    double b;
    if( c.ToDouble(&b) ) {
      l += " " + c;
      n++;
    }
    else {
      MuniFatalError("Failed to interpret `"+c+"' as a numerical value.");
      return;
    }
  }
}



// -----

wxString MuniAction::GetAdvanced(const wxString& filename) const
{
  if( ! advanced ) return "";

  FITSFileName fitsfile(filename);
  wxString fullname = fitsfile.GetFullName();
  wxString output;

  if( mask == "" )
    MuniFatalError("Regular expression mask is empty.");

  wxRegEx rf(pattern);
  if( ! rf.IsValid() )
    MuniFatalError("Regular expression pattern `"+pattern+"' is invalid.");

  if( rf.Matches(filename) ) {

    output = fullname;
    int e = rf.Replace(&output,mask);
    if( e == -1 )
      MuniFatalError("Regular expression `"+mask+"'failed to replace.");
    if( e == 0 )
      wxLogWarning("No replace for regular expression mask: `"+mask+"'.");
  }
  else
    MuniFatalError("Regular expression has no match on the filename `"
		   + filename+"'.");

  return output;
}

wxString MuniAction::GetOutput(const wxString& name) const
{
  FITSFileName ftname(name);
  wxString filename(ftname.GetFullName());

  if( target == "" )
    return filename;
  else {
    wxFileName src(filename);
    wxFileName dest(target,src.GetFullName());
    return dest.GetFullPath();
  }
}

wxString MuniAction::GetFormat(size_t n) const
{
  wxASSERT(format != "");

  wxString output;
  output.Printf(format,static_cast<int>(n));
  return output;
}


bool MuniAction::SetTargetDir(const wxString& t)
{
  target = t;

  // check existence of the target directory
  wxFileName d;
  d.AssignDir(target);
  return d.DirExists();
}

bool MuniAction::FindDataDir(const wxString& file, wxString *filepath) const
{
  wxString dir, datadir;

#ifdef MUNIPACK_DATA_DIR
  dir = wxString(MUNIPACK_DATA_DIR);
#endif

  if( wxGetEnv("MUNIPACK_DATA_DIR",&datadir) )
    dir = datadir;

  wxFileName path(dir,file);
  if( path.FileExists() ) {
    *filepath = path.GetFullPath();
    return true;
  }

  wxLogError("Munipack::FindDataDir(): an improper setup of MUNIPACK_DATA_DIR=`"
	     +dir+"' for file: `"+file+"'");
  return false;
}
