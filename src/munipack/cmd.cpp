/*

  Munipack - command line parser

  Copyright © 2010 - 2026 F.Hroch (hroch@physics.muni.cz)

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

#include "munipack.h"
#include <wx/wx.h>
#include <wx/app.h>
#include <wx/cmdline.h>
#include <wx/regex.h>
#include <wx/tokenzr.h>
#include <wx/filename.h>
#include <wx/file.h>
#include <wx/wfstream.h>
#include <wx/txtstrm.h>

#define COMMON_OPTIONS_SINGLE \
  AddOption("o","output","output file");

#define COMMON_OPTIONS_ENABLE_OVERWRITE	\
  AddSwitch("","enable-overwrite","overwrite existing files");

#define COMMON_OPTIONS_MULTIPLE \
  AddOption("t","target-directory","set target directory");		\
  AddSwitch("O","","switch-on advanced output filenames processing");\
  AddOption("","pattern","matching pattern (for -O)");	\
  AddOption("","mask","mask for output file(s) (for -O)");	\
  AddOption("","format","format for output files (for -O)");

#define COMMON_OPTIONS_BITPIX \
  AddOption("B","bitpix","bitpix (-32,32,16,8)",wxCMD_LINE_VAL_NUMBER);

#define COMMON_OPTIONS_PHCORR \
  AddOption("bitmask",wxEmptyString,"mask frame"); \
  AddOption("st","saturate","saturation in counts",wxCMD_LINE_VAL_DOUBLE); \
  AddOption("th","threshold","threshold in counts",wxCMD_LINE_VAL_DOUBLE);

#define COMMON_OPTIONS_PHQ \
    AddOption("q","quantity","select: PHOTRATE,FNU,FLAM,ABMAG,STMAG");
//    AddOption("q","quantity","instrumental: COUNT,RATE,MAG\n\tcalibrated: PHOTON,PHOTRATE,PHOTNU,PHOTLAM,FLUX,FNU,FLAM,MAG,ABMAG,STMAG");

#ifdef HAVE_CAF
#define COMMON_OPTIONS_PARALLEL \
  AddOption("mproc","mproc","number of parallel processes",wxCMD_LINE_VAL_NUMBER); \
  AddOption("","Xproc","additional options for parallel processing");
#endif


MuniCmdLineParser::MuniCmdLineParser(int xargc, wxChar *xargv[]):
  wxCmdLineParser()
{
    //    SetLogo(INFOTEXT " " XVERSION ", " COPYLEFT "\n" INFOTEXTFULL "\n");
  AddUsageText(
	       /*
                 "\nGeneral Actions:\n"
		 " calibrate\tcalibrate\n"
		 */
		 "\nActions:\n"
		 " bias\t\tAveraged bias frame\n"
		 " dark\t\tAveraged dark frame\n"
		 " flat\t\tAveraged flat-field frame\n"
		 " phcorr\t\tApply photometric corrections\n"
		 " colouring\tCompose colour frames\n"
		 " kombine\tCombine frames (deep exposures, mosaics)\n"
		 " cone\t\tCone search in a Virtual Observatory catalogue\n"
		 " astrometry\tAstrometry calibration\n"
		 " find\t\tSearch for stars on frames\n"
		 " aphot\t\tAperture photometry\n"
		 " gphot\t\tGrowth-Curve photometry\n"
		 " phcal\t\tPhotometry calibration\n"
		 " artificial\tArtificial frames\n"
		 " phfotran\tDetermine photometry system transformation\n"
		 " timeseries\tList arbitrary quantity to time-series\n"
		 " sesame\t\tAstronomical name resolver\n"
		 "\nLow Level Actions:\n"
		 " votable\tConversions of VOTables\n"
		 " fits\t\tOperations on FITS files\n"
		 " cross\t\tCross-match on FITS tables\n"
	       );

  AddSwitch("","version","print version and license");
  AddParam("action",wxCMD_LINE_VAL_STRING,wxCMD_LINE_PARAM_OPTIONAL);

  SetCmdLine(xargc,xargv);
}


MuniCmdLineParser::MuniCmdLineParser(const wxString& action, int xargc,
				     wxChar *xargv[]):
  wxCmdLineParser()
{
  wxASSERT(action != "" && xargc > 0); // Requires the action, at least

  // filters command-line options:
  // * actions are not passed
  // * remove stdin "-" switch
  // * remove @in item
  // * read input files and defines files

  wxChar **argv = new wxChar*[xargc-1];
  int argc = 0;
  for(int i = 0; i < xargc; i++) {
    wxString arg(xargv[i]);
    if( arg == action )
      // the action itself is omitted
      ;
    else if( arg == "-" || arg.substr(0,1) == "@" )
      // Look for - or @in option and eventually read the input files
      InputByFile(arg.substr(1));
    else
      argv[argc++] = xargv[i];
  }

  SetCmdLine(argc,argv);

  delete[] argv;


  /* There two alternatives how to handle the action as parameters:
   a) Remove it from options: It does't appear as a parameter in help prints.
   b) Include it in
       AddParam("artificial",wxCMD_LINE_VAL_STRING,wxCMD_LINE_OPTION_MANDATORY);
      as it appear in the end of print list. They should be ignored later
      during options proccessing.
  */

  // The text represents a brace for remove original Usage text
  //  AddUsageText("#");

  if( action == "bias" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)");
    AddUsageText("* Averaged bias frame.");
    COMMON_OPTIONS_PHCORR
    COMMON_OPTIONS_SINGLE
    COMMON_OPTIONS_BITPIX
    AddParam("file(s)",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }
  else if( action == "dark" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)");
    AddUsageText("* Averaged dark frame.");
    AddOption("bias","","bias frame");
    COMMON_OPTIONS_PHCORR
    COMMON_OPTIONS_SINGLE
    COMMON_OPTIONS_BITPIX
    AddParam("file(s)",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "flat" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)");
    AddUsageText("* Averaged flat-field frame.");

    AddOption("gain","","provide value of gain (rather by FITS header)",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("bias",wxEmptyString,"bias frame");
    AddOption("dark",wxEmptyString,"dark frame");
    AddOption("xdark",wxEmptyString,"dark frame exposure time factor",
	      wxCMD_LINE_VAL_DOUBLE);
    COMMON_OPTIONS_PHCORR
    AddOption("","approximation","accuracy: basic, standard (default)");
    COMMON_OPTIONS_SINGLE
    COMMON_OPTIONS_BITPIX
    AddParam("file(s)",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "phcorr" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)");
    AddUsageText("* Apply photometric corrections (gain, bias, dark, flat-field).");

    AddOption("gain","","provide value of gain (rather than FITS keyword)",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("flat","","flat-field frame");
    AddOption("bias","","bias frame");
    AddOption("dark","","dark frame");
    AddOption("xdark","","dark frame exposure time factor",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("bitmask","","mask frame");
    AddOption("xbitmask","","mask frame as ZERO, MEAN(default), MEDIAN");
    AddOption("box","","box size for MEDIAN bitmask",wxCMD_LINE_VAL_NUMBER);
    AddSwitch("","normalise-flat","normalise flat-field prior to use");
    COMMON_OPTIONS_ENABLE_OVERWRITE
    COMMON_OPTIONS_MULTIPLE
    COMMON_OPTIONS_BITPIX
    AddParam("file(s)",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "colouring" ) {

    AddUsageText("Usage: munipack "+action+" [options] files");
    AddUsageText("* Compose colour frames");

    AddOption("c","cspace-input","input colour-space, default: ask oracle");
    AddOption("","cspace-output","output colour-space, default:'CIE 1931 XYZ'");
    AddOption("w","weights","w1,w2,.. weights of images");
    AddOption("q","ctphs","ctph1,ctph2,.. ctphs of images");
    AddOption("b","backs","b1,b2,.. backgrounds of images");
    AddSwitch("","disable-back","don't estimate background of images");
    AddOption("","white-spot","x,y coordinates of centre of white spot");
    AddOption("","white-star","x,y coordinates of centre of white star");
    AddOption("","white-radius","for the spot or the star in pixels, default: 7",
	      wxCMD_LINE_VAL_DOUBLE);
    AddSwitch("","list","available colour-space transformations");
    AddOption("","ctable","colour table transformation definition");
    AddOption("","phsystab","photometry systems definition table");
    COMMON_OPTIONS_SINGLE
    COMMON_OPTIONS_BITPIX
    AddParam("files",wxCMD_LINE_VAL_STRING,
  	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "cone" ) {

    AddUsageText("Usage: munipack "+action+" [options] [--] RA DEC");
    AddUsageText("* Cone search in catalogues of Virtual Observatory.");

    AddOption("c","cat","catalogue (see --list-catalogues)");
    AddOption("r","radius","radius for cone search [deg] (default: 6 arcmin)",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","magmin","bright object magnitude limit",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","magmax","faint object magnitude limit",wxCMD_LINE_VAL_DOUBLE);
    AddSwitch("","Johnson-patch",
	      "convert Gunn's ri magnitudes to Johnson RI (UCAC4 only)");
    AddOption("s","sort","sort by the column");
    AddOption("","par","add optional parameters");
    AddOption("","url","service URL");
    AddOption("","id","catalogue identifier");
    AddOption("","server","Virtual observatory server (see --list-servers)");
    AddSwitch("","list-catalogues","show predefined catalogues");
    AddSwitch("","list-servers","show available servers");
    AddOption("","vocat","VO catalogue config file (column IDs)");
    AddOption("","type","type of output file: fits,votable");
    COMMON_OPTIONS_SINGLE
    AddParam("ra",wxCMD_LINE_VAL_DOUBLE,wxCMD_LINE_PARAM_OPTIONAL);
    AddParam("dec",wxCMD_LINE_VAL_DOUBLE,wxCMD_LINE_PARAM_OPTIONAL);
  }

  else if( action == "votable" ) {

    AddUsageText("Usage: munipack "+action+" [options] file");
    AddUsageText("* Conversions of VOTables.");

    AddOption("","type","type of output file: fits,csv,txt,svg,xml");
    AddOption("pt","","projection: TAN (gnomonic, default), ..");
    AddOption("pa","","Right Ascension of projection centre",wxCMD_LINE_VAL_DOUBLE);
    AddOption("pd","","Declination of projection centre",wxCMD_LINE_VAL_DOUBLE);
    AddOption("ps","","scale of projection",wxCMD_LINE_VAL_DOUBLE);
    AddOption("pw","","canvas width",wxCMD_LINE_VAL_NUMBER);
    AddOption("ph","","canvas height",wxCMD_LINE_VAL_NUMBER);
    AddOption("ml","","magnitude limit",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","col-ra","Right Ascension column");
    AddOption("","col-dec","Declination column");
    AddOption("","col-mag","magnitude column");
    COMMON_OPTIONS_SINGLE
    AddParam("file",wxCMD_LINE_VAL_STRING,wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "sesame" ) {

    AddUsageText("Usage: munipack "+action+" [options] object");
    AddUsageText("* Astronomical objects name resolver by Virtual Observatory.");
    AddParam("object",wxCMD_LINE_VAL_STRING);
  }

  else if( action == "astrometry" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)[,results(s)]");
    AddUsageText("* Astrometry calibration.");

    AddOption("m","mode","match, cross, sequence or manual");
    AddOption("c","cat","catalogue");
    AddOption("r","ref","reference frame");
    AddOption("R","rel","reference frame (relative)");
    AddOption("p","projection","projection: TAN (gnomonic, default), ..");
    AddOption("","col-ra","Right Ascension column");
    AddOption("","col-dec","Declination column");
    AddOption("","col-pm-ra","Proper motion in Right Ascension column");
    AddOption("","col-pm-dec","Proper motion in Declination column");
    AddOption("","col-epoch-ra","Reference epoch in Right Ascension column");
    AddOption("","col-epoch-dec","Reference epoch in Declination column");
    AddOption("","col-mag","Magnitude-like column");
    AddOption("","xcen","centre of frame [pix] (default: width/2)",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","ycen","centre of frame [pix] (default: height/2)",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","rcen","centre of FOV in Right Ascension [deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","dcen","center of FOV in Declination [deg]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","scale","scale [pix/deg]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","angle","rotation angle [deg], clockwise positive",
	      wxCMD_LINE_VAL_DOUBLE);
    AddSwitch("","reflex","set reflection");
    AddOption("","radesys","reference frame (default: IRCS)");
    AddOption("","epoch","reference epoch",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","fit","robust (default) or squares");
    AddOption("","units","output units: deg,arcmin,arcsec,mas,uas,pix");
    AddOption("","sig","errors of centroids on frame [pix]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","sigcat","coordinate errors of cat/ref [deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","flxmin","minimum of relative fluxes, default=0.8",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","flxmax","maximum of relative fluxes, default=1.3",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","minmatch","minimal length of match sequence",
	      wxCMD_LINE_VAL_NUMBER);
    AddOption("","maxmatch","maximal length of match sequence",
	      wxCMD_LINE_VAL_NUMBER);
    AddOption("","maxstars","total number of astrometry stars",
	      wxCMD_LINE_VAL_NUMBER);
    AddOption("","luckymatch","minimal reliable sequence, default=3",
	      wxCMD_LINE_VAL_NUMBER);
    AddOption("","seq1","sequence for reference");
    AddOption("","seq2","sequence for calibrated");
    AddSwitch("","disable-lucky-match","finish at first success match");
    AddSwitch("","disable-flux-check","fluxes are not used during matching");
    AddSwitch("","disable-rms-check","disable RMS check");
    AddSwitch("","disable-save","don't save calibration to header");
    AddSwitch("","remove","remove complete astrometry calibration");
    AddSwitch("","report","add FITS extension with a report");
    AddSwitch("","show-defaults","show defaults");
    AddOption("","vocat","VO catalogue config file (column IDs)");
    COMMON_OPTIONS_MULTIPLE
#ifdef HAVE_CAF
    COMMON_OPTIONS_PARALLEL
#endif
    AddParam("file(s)[,results(s)]",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "find" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)");
    AddUsageText("* Search for stars on frames.");

    AddOption("f","fwhm","FWHM in pixels",wxCMD_LINE_VAL_DOUBLE);
    AddOption("th","threshold","threshold in sigmas above background",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","saturate","saturation in ADU",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","lothresh","lower for threshold in sigmas",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","rndlo","lower for round",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","rndhi","higher for round",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","shrplo","lower for sharp",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","shrphi","higher for sharp",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","read-noise","read noise in ADU (default: frame header)",
	      wxCMD_LINE_VAL_DOUBLE);
    COMMON_OPTIONS_MULTIPLE
    AddParam("file(s)",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "aphot" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)");
    AddUsageText("* Aperture photometry");

    AddOption("","apertures","a1,a2,... define aperture radii in pixels");
    AddOption("","ring","ri,ro inner and outer sky ring radius in pixels");
    AddSwitch("","enable-ellipticity","enable elliptic apertures");
    AddOption("","eccentricity","of apertures (by default: circle)",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","inclination","inclination of major semiaxis, degrees",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","coo","coordinates of objects");
    COMMON_OPTIONS_MULTIPLE
    AddParam("file(s)",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "gphot" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)[,results(s)]");
    AddUsageText("* Growth-curve photometry");

    AddOption("th","threshold","threshold in sigmas above background",wxCMD_LINE_VAL_DOUBLE);
    COMMON_OPTIONS_MULTIPLE
    AddParam("file(s)[,results(s)]",wxCMD_LINE_VAL_STRING,wxCMD_LINE_PARAM_MULTIPLE);
  }

    /*
  case ID_CALIBRATE:
    AddParam("calibrate",wxCMD_LINE_VAL_STRING);
    AddOption("r",wxEmptyString,"search radius in degrees");
    AddParam("ra,dec file(s)",wxCMD_LINE_VAL_STRING,
  	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
    break;
    */

  else if( action == "timeseries" ) {

    AddUsageText("Usage: munipack "+action+" [options] ... [file(s)]");
    AddUsageText("* List keywords or objects to time-series.");

    AddOption("c","coo","coordinates of objects");
    AddOption("","cat","coordinates by FITS catalogue");
    AddOption("l","col","column label(s) to list");
    AddOption("K","keys","key[,key] header key(s) to list");
    AddOption("T","time-type","time types:  JD (default), MJD, HJD, PHASE");
    AddOption("","time-stamp","MID (default), BEGIN, END");
    AddOption("","coo-type","DEG (default), RECT");
    AddOption("","coo-col","column label(s) of coordinates (object search)");
    AddOption("","tol","crosmatch coordinate search limit [deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","lc-epoch","reference time for light curve elements [JD]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","lc-period","period for light curve elements [days]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","extname","specify FITS extension name");
    AddSwitch("","stdout","list results to standard output");
    AddSwitch("","enable-filename","print filename");
    AddSwitch("","disable-timetype","disable print of time-type quantity");
    AddSwitch("","enable-horizon","print horizontal coordinates");
    AddSwitch("","enable-airmass","print air-mass");
    COMMON_OPTIONS_SINGLE
    AddParam("[file(s)]",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL | wxCMD_LINE_PARAM_MULTIPLE);
  }

  // case ID_LIST:
  //   AddParam("list",wxCMD_LINE_VAL_STRING);
  //   //    AddOption("c","create","cat (star catalogue), lc (light curve, default)");
  //   AddOption("T","time-type","JD-Julian day (default), MJD-modified JD, HJD-heliocentric JD (+), PHASE(*)");
  //   AddOption(wxEmptyString,"time-stamp","mid (default), begin, end");
  //   AddOption("C","coordinate-type","S-spherical (default), P-cartesian, pixels");
  //   //    AddSwitch("s","spherical","spherical coordinates (default)");
  //   //    AddSwitch("p","pixels","cartesian coordinates");
  //   AddSwitch(wxEmptyString,"index","specify index");
  //   AddOption(wxEmptyString,"col","column(s) to list");
  //   AddOption(wxEmptyString,"key","header keyword(s)");
  //   //    AddSwitch(wxEmptyString,"mag","output in instrumental magnitudes");
  //   //    AddSwitch(wxEmptyString,"full","full list");
  //   //    AddSwitch(wxEmptyString,"pos","position list");
  //   AddOption(wxEmptyString,"longitude","geographical longitude of station (+east)[deg]",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"latitude","geographical latitude of station (+north)[deg]",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"ra","Right Ascension [deg]",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"dec","Declination [deg]",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"minim","basic point of light curve elements in JD (*)",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"period","period light curve elements in days (*)",wxCMD_LINE_VAL_DOUBLE);
  //   //    AddOption(wxEmptyString,"coo","coordinates: horizontal, pixels, equatorial");
  //   AddSwitch(wxEmptyString,"mag","output in magnitudes instead counts");
  //   AddSwitch(wxEmptyString,"flux","output in fluxes instead counts");
  //   AddSwitch(wxEmptyString,"calibr","calibrated quantities (for --mag or --flux)");
  //   AddSwitch(wxEmptyString,"diffmag","differential magnitudes (for --mag)");
  //   /*
  //   AddOption(wxEmptyString,"zeromag","zero magnitude (along with --mag, --calibr)",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"photflux","flux of 1 ADU (count) source (--flux,--calibr)",wxCMD_LINE_VAL_DOUBLE);
  //   */
  //   AddOption(wxEmptyString,"zeromag","zero magnitude (along with -Q mag, --calibr)",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"photflux","flux of 1 ADU (count) source (-Q flux,--calibr)",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"epoch","reference Julian date (proper motion)",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"naperture","aperture number",wxCMD_LINE_VAL_NUMBER);
  //   AddOption(wxEmptyString,"aperture","radius of aperture [pix]",wxCMD_LINE_VAL_DOUBLE);
  //   AddOption(wxEmptyString,"tol","search radius [deg]",wxCMD_LINE_VAL_DOUBLE);
  //   AddSwitch(wxEmptyString,"stdout","results print also to standard output");
  //   AddSwitch(wxEmptyString,"print-filename","print filenames");
  //   AddOption(wxEmptyString,"file","read coordinates from the file, (-) for stdin");
  //   COMMON_OPTIONS_SINGLE
  //   AddParam("[a,d[,mua,mud]]|[x,y[,mux,muy]]... file(s)",wxCMD_LINE_VAL_STRING,
  // 	     wxCMD_LINE_PARAM_MULTIPLE);
  //   break;

  else if( action == "fits" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)");
    AddUsageText("* Operations on FITS files.");

    AddSwitch("lh","header","list header");
    AddSwitch("lt","table","list table");
    AddSwitch("li","image","list image");
    AddOption("K","keys","key[,key,..] print header values by keyword(s)");
    AddSwitch("","value","print only values during keywords print");
    AddSwitch("","shell","shell-friendly format of keywords print");
    AddSwitch("","update","add or update header record");
    AddOption("","key","specify keyword of updated record");
    AddOption("","val","specify value of updated record");
    AddOption("","com","specify comment of updated record");
    AddOption("","templ","update records by this file");
    AddOption("","remove-keys","key[,key,..] remove keyword(s)");
    AddOption("","remove-extensions","extname[,extname,..] remove extension(s)");
    AddSwitch("","dump","dump FITS to plain text");
    AddSwitch("","restore","restore FITS from plain text");
    AddSwitch("","cat","cat (copy) input on output");
    COMMON_OPTIONS_MULTIPLE
    AddParam("[file(s)]",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "cross" ) {

    AddUsageText("Usage: munipack "+action+" [options] file1 file2");
    AddUsageText("* Cross-match of FITS tables.");
    AddOption("","tol","search radius [deg]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","ftol","relative flux tolerance",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","col-ra","Right Ascension columns (.,.)");
    AddOption("","col-dec","Declination columns (.,.)");
    AddOption("","col-pm-ra","Proper motion in Right Ascension columns (.,.)");
    AddOption("","col-pm-dec","Proper motion in Declination columns (.,.)");
    AddOption("","col-epoch-ra","Reference epoch in Right Ascension column");
    AddOption("","col-epoch-dec","Reference epoch in Declination column");
    AddOption("","col-mag","magnitude columns (.,.)");
    AddOption("","vocat","VO catalogue config file (column IDs)");
    COMMON_OPTIONS_SINGLE
    AddParam("file1 file2",wxCMD_LINE_VAL_STRING,wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "kombine" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)");
    AddUsageText("* Compose frames (deep exposures, mosaics).");

    AddOption("p","projection","projection: TAN (gnomonic, default), ..");
    AddOption("","xcen","centre of frame [pix] (default: width/2)",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","ycen","centre of frame [pix] (default: height/2)",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","rcen","centre of FOV in Right Ascension [deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","dcen","centre of FOV in Declination [deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","pm-ra","proper motion in Right Ascension [deg/day]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","pm-dec","proper motion in Declination [deg/day]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","pm-jdref","proper motion reference JD",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","width","width of output [pix]",wxCMD_LINE_VAL_NUMBER);
    AddOption("","height","height of output [pix]",wxCMD_LINE_VAL_NUMBER);
    AddOption("","scale","scale [pix/deg]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","angle","position angle [deg]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","reflex","[yes|no] set/unset reflection");
    AddOption("i","interpol","[near,bilinear,bicubic,bi3conv] interpol");
    AddSwitch("","disable-back","does not subtract background");
    AddSwitch("","arith","compute average by arithmetical mean");
    COMMON_OPTIONS_SINGLE
    COMMON_OPTIONS_BITPIX
    AddParam("file(s)",wxCMD_LINE_VAL_STRING,
	     wxCMD_LINE_PARAM_OPTIONAL|wxCMD_LINE_PARAM_MULTIPLE);
  }

  else if( action == "phcal" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)[,results(s)]");
    AddUsageText("* Photometry calibration.");

    AddOption("C","cal","specify the calibration ratio(s) by hand");
    AddOption("c","cat","reference photometry catalogue");
    AddOption("r","ref","reference frame");
    AddOption("f","filters","list of filters");
    AddOption("th","threshold","select stars with the threshold above sky",wxCMD_LINE_VAL_DOUBLE);
    AddOption("e","maxerr","select stars with the maximal relative error",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","apcorr","aperture correction",wxCMD_LINE_VAL_DOUBLE);
    AddOption(wxEmptyString,"col-ra","Right Ascension column");
    AddOption(wxEmptyString,"col-dec","Declination column");
    AddOption(wxEmptyString,"col-mag","magnitude column(s)");
    AddOption(wxEmptyString,"col-magerr","magnitude std. error column(s)");
    AddOption("","photsys-ref","reference photometry system (catalogue)");
    AddOption("","photsys-instr","instrumental photometry system (frames)");
    AddOption(wxEmptyString,"area","area of input aperture [m2]",wxCMD_LINE_VAL_DOUBLE);
    COMMON_OPTIONS_PHQ
    AddOption(wxEmptyString,"tratab","instrumental to reference photo-system table");
    AddOption(wxEmptyString,"phsystab","photometry systems definition table");
    AddSwitch(wxEmptyString,"advanced","advanced format (additional extensions included)");
    AddSwitch(wxEmptyString,"list","show available photometry systems");
    AddOption(wxEmptyString,"tol","search radius [deg]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","saper","select aperture number",wxCMD_LINE_VAL_NUMBER);
    AddOption("","vocat","VO catalogue config file (column IDs)");
    COMMON_OPTIONS_MULTIPLE
    AddParam("file(s)[,results(s)]",wxCMD_LINE_VAL_STRING,
  	     wxCMD_LINE_PARAM_MULTIPLE|wxCMD_LINE_PARAM_OPTIONAL);
  }

  else if( action == "phfotran" ) {

    AddUsageText("Usage: munipack "+action+" [options] file(s)");
    AddUsageText("* Determine photometry system colour transformation.");

    AddOption("c","cat","reference photometry catalogue");
    AddOption("f","filters","list of filters");
    AddOption("C","cal","calibration ratios in the filters");
    AddOption("E","extin","extinction coefficients k (exp(k*X) in filters");
    AddOption(wxEmptyString,"col-ra","Right Ascension column");
    AddOption(wxEmptyString,"col-dec","Declination column");
    AddOption(wxEmptyString,"col-mag","magnitude column(s)");
    AddOption(wxEmptyString,"col-magerr","magnitude std. error column(s)");
    AddOption("","photsys-ref","reference photometry system (catalogue)");
    AddOption("","photsys-instr","instrumental photometry system (frames)");
    AddOption(wxEmptyString,"area","area of input aperture [m2]",wxCMD_LINE_VAL_DOUBLE);
    AddOption(wxEmptyString,"phsystab","photometry systems definition table");
    AddSwitch(wxEmptyString,"list","show available photometry systems");
    AddOption(wxEmptyString,"tol","search radius [deg]",wxCMD_LINE_VAL_DOUBLE);
    COMMON_OPTIONS_SINGLE
    AddParam("file(s)",wxCMD_LINE_VAL_STRING,
  	     wxCMD_LINE_PARAM_MULTIPLE|wxCMD_LINE_PARAM_OPTIONAL);
  }

  else if( action == "artificial" ) {

    AddUsageText("Usage: munipack "+action+" [options]");
    AddUsageText("* Photometrically realistic artificial frames.");

    AddOption("","psf","PSF function: SEEING, MOFFAT, GAUSS");
    AddOption("","psf-file","PSF by file");
    AddOption("","spread","seeing spread method: FFT, RANDOM, AUTO");
    AddOption("","gain","gain in e-/ADU",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","saturate","saturation in ADU",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","cosmics-freq","frequency of cosmic rays [1/s]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","cosmics-cmax","electrons in the maximum of cosmics rays",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","cosmics-len","cosmics rays lengths",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","hwhm","PSF width [pix]: sigma for GAUSS and SEEING, "
	      "alpha for MOFFAT",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","airy","radius of Airy spot [pix]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","beta","Moffat exponent",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","eccentricity","of ellipsis (by default: circle)",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","inclination","inclination of major semiaxis, degrees",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","maglim","magnitude limit",wxCMD_LINE_VAL_NUMBER);
    AddOption("","lc-model","light curve model: Cepheid,Algol,wave(default)");
    AddOption("","lc-table","light curve model in the table");
    AddOption("","lc-fourier","light curve model by Fourier coefficients");
    AddOption("","lc-mag","mean magnitude of the variable object",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","lc-amp","amplitude of the curve",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","lc-jd0","reference Julian date for the light curve",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","lc-per","period for the light curve [days]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","lc-ra","variable object Right Ascension [deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","lc-dec","variable object Declination [deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","sky-mag","sky brightness [mag/arcsec2]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","sky-grad-x","sky brightness gradient in x [mag/arcsec2/pix],"
	      " positive is fainter",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","sky-grad-y","sky brightness gradient in y [mag/arcsec2/pix],"
	      " positive is fainter",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","area","area of input aperture [m2]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","diameter","diameter of input aperture [m]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","exptime","exposure time [s]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","qeff","quantum efficiency",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","dark","dark current [e-/s/pix]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","readnoise","read-out noise [e- rms]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","gnoise","(deprecated) additional background rate "
	      "[ph/s/m2/arcsec2], merged in quadrature with the sky",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","seed","random generator seed (reproducible frames)",
	      wxCMD_LINE_VAL_NUMBER);
    AddSwitch("","atmosphere","apply atmosphere modelling: extinction+seeing");
    AddOption("","extk","extinction coefficient",wxCMD_LINE_VAL_DOUBLE);
    AddOption("","long","geographic longitude of station (+east)[deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","lat","geographic latitude of station (+north)[deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","date","initial date as YYYY-MM-DD");
    AddOption("","time","initial time as HH:MM:SS");
    AddOption("","count","total number of generated files",
	      wxCMD_LINE_VAL_NUMBER);
    AddOption("","timestep","timestep in seconds",wxCMD_LINE_VAL_DOUBLE);
    AddOption("f","filter","photometry filter");

    AddOption("c","cat","reference photometry catalogue");
    AddOption("","galcat","galaxy catalogue");

    AddOption("","fov","field of view [deg]",wxCMD_LINE_VAL_DOUBLE);
    AddOption("p","projection","projection: TAN (gnomonic, default), ..");
    AddOption(wxEmptyString,"rcen","centre of FOV in Right Ascension [deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption(wxEmptyString,"dcen","centre of FOV in Declination [deg]",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption(wxEmptyString,"scale","scale [pix/deg]",wxCMD_LINE_VAL_DOUBLE);
    AddOption(wxEmptyString,"angle","rotation angle [deg], clockwise positive",
	      wxCMD_LINE_VAL_DOUBLE);
    AddOption("","width","width of output [pix]",wxCMD_LINE_VAL_NUMBER);
    AddOption("","height","height of output [pix]",wxCMD_LINE_VAL_NUMBER);
    AddOption("","col-ra","Right Ascension column");
    AddOption("","col-dec","Declination column");
    AddOption("","col-pm-ra","Proper motion in Right Ascension column");
    AddOption("","col-pm-dec","Proper motion in Declination column");
    AddOption("","col-epoch-ra","Reference epoch in Right Ascension column");
    AddOption("","col-epoch-dec","Reference epoch in Declination column");
    AddOption("","col-mag","Magnitude-like column");
    AddSwitch("","disable-noise","noise free frames (no photon and read noise)");
    AddOption("","photsys","photometry system");
    AddOption("","phsystab","photometry systems definition table");
    AddOption("","fits-keys-file","add FITS keys into a header");
    AddOption("","flat-model-file","a FITS file with flat model");
    AddOption("","vocat","VO catalogue config file (column IDs)");
    AddOption("","mask","mask for output files (artificial_???.fits)");
    COMMON_OPTIONS_BITPIX
#ifdef HAVE_CAF
    COMMON_OPTIONS_PARALLEL
#endif
  }

}

/*
void MuniCmdLineParser::Usage() const
{
  wxString u(GetUsageString());

  // The default Usage prints no actions, we are remove one and replace
  // with more suitable form. Unfortunatelly, options list is omitted.
  // I don't know which way is suitable.
  wxRegEx re("(Usage:.*[#]\n)(.*)");
  wxASSERT(re.IsValid());
  if( re.Matches(u) ) {
    wxString t = re.GetMatch(u);
    re.Replace(&t,"\\2");
    wxPrintf(t);
  }
  else
    wxPrintf(GetUsageString());

}
*/

wxArrayString MuniCmdLineParser::GetFiles() const
{
  wxArrayString files;

  for(size_t i = 0; i < GetParamCount(); i++)
    files.push_back(GetParam(i));

  for(size_t i = 0; i < lstfiles.GetCount(); i++)
    files.push_back(lstfiles[i]);

  return files;
}


void MuniCmdLineParser::InputByFile(const wxString& lstname)
{
  wxFile lstfile;
  bool is_stdin = lstname == "";

  if( is_stdin )
    lstfile.Attach(wxFile::fd_stdin);
  else
    lstfile.Open(lstname);

  if( lstfile.IsOpened() ) {
    wxFileInputStream input(lstfile);
    wxTextInputStream text(input);
    while( input.IsOk() && input.CanRead() ) {
      wxString line(text.ReadLine());
      wxString filename = line.BeforeFirst('#');
      filename.Trim();
      if( ! filename.IsEmpty() )
	lstfiles.push_back(filename);
    }

    if( is_stdin )
      // Stdin must be pulled back to the parent process,
      // else sockets of calling process are unavailable.
      lstfile.Detach();

    lstfile.Close();
  }

  if( lstfile.GetLastError() != 0 ) {
    if( is_stdin )
      wxLogFatalError("Failed to read the standard input.");
    else
      wxLogFatalError("Failed to read the input file `"+lstname+"'.");
  }
}
