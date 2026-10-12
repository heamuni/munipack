#
#  Copyright © 1997-2025 F.Hroch (hroch@physics.muni.cz)
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
LDADD = ../lib/libtitsio.a $(LIBOAKLEAF) $(LIBCFITSIO)

pkglibexec_PROGRAMS = dark flat phcorr

.NOTPARALLEL:

dark_SOURCES = fitscorr.f08 dark.f08
dark_LDADD = ../lib/libclipmean.a $(LDADD)
flat_SOURCES = fitscorr.f08 flat.f08
phcorr_SOURCES = fitscorr.f08 phcorr.f08
