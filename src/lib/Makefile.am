#
#  munipack libraries
#
#  Copyright © 1997 - 2026 F.Hroch (hroch@physics.muni.cz)
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
#

noinst_LIBRARIES = libminpacks.a libtitsio.a libsort.a \
                   liblmin.a libstat.a libfmm.a libmfitsio.a \
                   libastrotrafo.a libastrosphere.a libphotometry.a \
                   libpppack.a libnoise.a libpsfmodel.a
#                   libphio.a libpppack.a libnoise.a

.NOTPARALLEL:

libminpacks_a_SOURCES = minpack.f90 minpacks.f95
libtitsio_a_SOURCES = cfitsio.f08 fitsio.f08 fitsio_mmviii.f08 trajd.f08 titsio.f08
libfmm_a_SOURCES = fmin.f95 zeroin.f95 fmm.f95
libsort_a_SOURCES = qsort.f95 ssort.f95
liblmin_a_SOURCES = nelmin.f95
libstat_a_SOURCES = mean.f95
libmfitsio_a_SOURCES = mfitsio.cpp
#libastrotrafo_a_SOURCES = projections.f08 astrotrafo.f08
libastrotrafo_a_SOURCES = wcslib.f08 astrowcs.f08
libastrosphere_a_SOURCES = trajd.f08 astrosphere.f08
libphotometry_a_SOURCES = photoconv.f08 phsysfits.f08
#libphio_a_SOURCES = phio.f08
libpppack_a_SOURCES = pppack.f90
libnoise_a_SOURCES = noise.f08
libpsfmodel_a_SOURCES = psfmodel.f08

noinst_HEADERS = mfitsio.h
EXTRA_DIST = fversion.inc
