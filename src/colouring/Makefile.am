#
#  Copyright © 1997-2026 F.Hroch (hroch@physics.muni.cz)
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

AM_FCFLAGS = $(FC_MODINC). $(FC_MODINC)../lib $(AM_FC_MODINC)
ldir = ../lib
LDADD =  $(ldir)/libphotometry.a $(ldir)/libtitsio.a \
        $(LIBOAKLEAF) $(LIBCFITSIO) $(LIBM)
#  $(ldir)/libminpacks.a      $(LIBOAKLEAF) $(LIBCFITSIO) $(LIBMINPACK) $(LIBM)

pkglibexec_PROGRAMS = colouring

.NOTPARALLEL:

colouring_SOURCES = colourfits.f08 ctables.f08 colour.f08 colouring.f08

dist_pkgdata_DATA = ctable.dat
