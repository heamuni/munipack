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
FITSLIBS = $(ldir)/libtitsio.a $(LIBCFITSIO)

pkglibexec_PROGRAMS = astrometry cross

.NOTPARALLEL:

astrometry_FCFLAGS = $(AM_FCFLAGS) $(CAF_FCFLAGS)
astrometry_SOURCES = match.f08 nearmatch.f08 estimator.f08 robfit.f08 deform.f08 \
                     astrofit.f08 astromatch.f08 astrofits.f08 astrometry.f08
# lsfit.f08 absfit.f08
astrometry_LDADD = $(ldir)/libastrosphere.a $(ldir)/libastrotrafo.a \
                   $(LIBOAKLEAF) $(WCSLIB_LIBS) $(FITSLIBS) $(CAF_LIBS)

cross_SOURCES = crossmatch.f95 cross.f08
cross_LDADD = $(ldir)/libastrosphere.a $(FITSLIBS)
