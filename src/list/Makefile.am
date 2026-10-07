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

AM_FCFLAGS = $(FC_MODINC). $(FC_MODINC)../lib
ldir = ../lib

pkglibexec_PROGRAMS = timeseries

.NOTPARALLEL:

timeseries_SOURCES = fitscat.f08 fitsts.f08 savelc.f08 timeseries.f08
timeseries_LDADD   = $(ldir)/libastrosphere.a $(ldir)/libtitsio.a \
                     $(LIBCFITSIO) $(LIBM)
