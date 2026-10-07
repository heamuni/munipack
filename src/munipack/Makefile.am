#
#  Copyright © 2010 - 2026 F.Hroch (hroch@physics.muni.cz)
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

AM_CPPFLAGS = -I$(top_srcdir)/lib -I$(top_srcdir)/vo \
              -DwxUSE_GUI=0 $(WX_CXXFLAGS_CLI) \
              -DMUNIPACK_LIBEXEC_DIR=\"$(pkglibexecdir)\" \
              -DMUNIPACK_DATA_DIR=\"$(pkgdatadir)\"

if CAF
CAF_CPPFLAGS = -DHAVE_CAF=1 -DMUNIPACK_CAF_RUN=\"$(cafrun)\"
endif

bin_PROGRAMS = munipack

pkglibexec_PROGRAMS = growreport

noinst_HEADERS = fitsfn.h mprocess.h munipack.h
noinst_LIBRARIES = libmprocess.a

libmprocess_a_SOURCES = mpipe.cpp mprocess.cpp

munipack_CPPFLAGS = $(AM_CPPFLAGS) $(CAF_CPPFLAGS)
munipack_SOURCES = munipack.cpp cmd.cpp action.cpp fitsfn.cpp \
                   bias.cpp dark.cpp flat.cpp phcorr.cpp \
                   colouring.cpp cone.cpp timeseries.cpp fits.cpp \
                   astrometry.cpp phcal.cpp votable.cpp find.cpp \
		   kombine.cpp aphot.cpp gphot.cpp  \
                   phfotran.cpp cross.cpp artificial.cpp sesame.cpp
# calibrate.cpp list.cpp grow_report.cpp
munipack_LDADD = libmprocess.a ../vo/libvocatconf.a \
                 ../vo/libvotable.a ../lib/libmfitsio.a \
                 $(WX_LIBS_VO) $(WX_LIBS_CLI) $(LIBCFITSIO) $(LIBM)

growreport_SOURCES = growreport.cpp
growreport_LDADD = $(LIBCFITSIO) $(LIBM)

dist_man_MANS = munipack.1
