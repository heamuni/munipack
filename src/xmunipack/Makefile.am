#
#  Copyright © 2009-2026 F.Hroch (hroch@physics.muni.cz)
#
#  This file is part of Munipack.
#
#  Munipack is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Munipack is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
#

bin_PROGRAMS = xmunipack

noinst_HEADERS = aphot.h dismag.h export.h icon.h mtool.h thumbcanvas.h  \
                 xmunipack.h browser.h display.h fileprop.h legend.h \
                 navigation.h tune.h zoomer.h caption.h dispreview.h \
                 fitsdisplay.h list.h plot.h tuner.h console.h enum.h \
                 fits.h Luv_Gamut.h preferences.h types.h dirscan.h \
                 event.h help.h mconfig.h render.h view.h

xmunipack_SOURCES = xmunipack.cpp fits.cpp itt.cpp fitscolour.cpp \
       palette.cpp fitsmeta.cpp browser.cpp list.cpp navigation.cpp \
       dirscan.cpp metarender.cpp fileprop.cpp render.cpp \
       splashing.cpp extlist.cpp view.cpp config.cpp dummy.cpp \
       display.cpp fitsstream.cpp \
       time.cpp value.cpp grid.cpp \
       tune.cpp tuner.cpp export.cpp \
       fitsdisplay.cpp legend.cpp \
       fitsgeometry.cpp fitshisto.cpp fitsstat.cpp fitstone.cpp ecdf.cpp \
       fitsarray.cpp fitstable.cpp event.cpp head.cpp icon.cpp help.cpp \
       fitscoo.cpp \
       plplot.cpp rawimport.cpp cone.cpp \
       preferences.cpp calibrate.cpp selectsrc.cpp \
       find.cpp findplot.cpp astrometry.cpp astropt.cpp astrolog.cpp \
       astrometer.cpp thumbcanvas.cpp mtool.cpp \
       human.cpp aphotaper.cpp aphotdisplay.cpp aphot.cpp \
       starlayer.cpp caption.cpp dismag.cpp dispreview.cpp zoomer.cpp zoomset.cpp \
       articons.cpp

#log.cpp photometry.cpp

iconsdir = $(datarootdir)/pixmaps/munipack
xdgicondir = $(datarootdir)/icons/hicolor/128x128/apps
desktopdir = $(datarootdir)/applications
metadir = $(datarootdir)/metainfo
XFLAGS = -DMUNIPACK_LIBEXEC_DIR=\"$(pkglibexecdir)\"            \
         -DMUNIPACK_DATA_DIR=\"$(pkgdatadir)\"                  \
         -DMUNIPACK_ICON_DIR=\"$(iconsdir)\"                    \
         -DMUNIPACK_HTML_DIR=\"$(htmldir)\"

xmunipack_CPPFLAGS = -I$(top_srcdir)/vo -I$(top_srcdir)/wxMathPlot \
                     -I$(top_srcdir)/munipack -I$(top_srcdir)/lib \
                     $(WX_CXXFLAGS_GUI) $(PLplot_CFLAGS) $(WCSLIB_CFLAGS) \
                     $(XFLAGS)
xmunipack_LDADD = ../wxMathPlot/libwxMathPlot.a \
                  ../munipack/libmprocess.a \
                  ../vo/libvocatconf.a \
                  ../vo/libvotable.a \
                  ../lib/libmfitsio.a \
                  $(PLplot_LIBS) $(WX_LIBS_GUI) $(WCSLIB_LIBS) $(LIBCFITSIO) $(LIBM)

dist_man_MANS = xmunipack.1

include image_list.mk
icons_DATA = $(image_list)
desktop_DATA = xmunipack.desktop
xdgicon_DATA = icons/munipack_icon.png
meta_DATA = cz.muni.physics.munipack.xmunipack.metainfo.xml
EXTRA_DIST = $(image_list) $(desktop_DATA) $(xdgicon_DATA) $(meta_DATA)


#Xmunipack.app: Info.plist xmunipack
#	-mkdir -p Xmunipack.app
#	-mkdir -p Xmunipack.app/Contents
#	-mkdir -p Xmunipack.app/Contents/MacOS
#	-mkdir -p Xmunipack.app/Contents/Resources
#	cp Info.plist Xmunipack.app/Contents/
#	echo -n 'APPL????' > Xmunipack.app/Contents/PkgInfo
#	cp xmunipack Xmunipack.app/Contents/MacOS/
#	cp icons/munipack_icon.png Xmunipack.app/Contents/Resources/Xmunipack.png
