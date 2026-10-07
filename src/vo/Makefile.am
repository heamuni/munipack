#
#  Copyright © 2010-2026 F.Hroch (hroch@physics.muni.cz)
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
# Obsolete!
# Parts based on wxHTTP of sesame and cone are obsolete, they are kept
# due compatibility and will be removed: vosesame.cpp, voclient.cpp.
#

XFLAGS= -DMUNIPACK_DATA_DIR=\"$(pkgdatadir)\"

AM_CPPFLAGS = -DwxUSE_GUI=0 -I$(top_srcdir)/lib $(WX_CXXFLAGS_VO) $(XFLAGS)
LDADD = $(WX_LIBS_VO) $(LIBCFITSIO) $(LIBM)
VOTABLIB = libvotable.a ../lib/libmfitsio.a

pkglibexec_PROGRAMS = cone votable sesame

noinst_HEADERS = vocatconf.h votable.h voclient.h vosesame.h common.h
noinst_LIBRARIES = libvotable.a libvocatconf.a

cone_SOURCES = common.cpp voclient.cpp cone.cpp
votable_SOURCES = common.cpp votab.cpp
sesame_SOURCES = vosesame.cpp common.cpp sesame.cpp
libvocatconf_a_SOURCES = vocatconf.cpp
libvotable_a_SOURCES = fitstable.cpp txtable.cpp csvtable.cpp \
     svgcanvas.cpp common.cpp votable.cpp
votable_LDADD = $(VOTABLIB) $(LDADD)
cone_LDADD = $(VOTABLIB) $(LDADD)

dist_pkgdata_DATA = VOcat_conf.xml
