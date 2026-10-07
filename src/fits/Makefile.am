#
#  Copyright © 2011-2025 F.Hroch (hroch@physics.muni.cz)
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

AM_CPPFLAGS = -I$(top_srcdir)/lib

#libexec_PROGRAMS = fitshead

#fitshead_SOURCES = fitshead.cpp
#fitshead_LDFLAGS = -lcfitsio -lm

noinst_HEADERS = fits.h fortranio.h

pkglibexec_PROGRAMS = fits
fits_SOURCES = fits.cpp structure.cpp header.cpp table.cpp image.cpp \
                 dump.cpp restore.cpp fortranio.cpp extremove.cpp cat.cpp
fits_LDADD   = ../lib/libmfitsio.a $(LIBCFITSIO) $(LIBM)
