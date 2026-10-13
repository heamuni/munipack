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

AM_FCFLAGS = $(FC_MODINC). $(FC_MODINC)../lib $(CAF_FCFLAGS)
ldir = ../lib

pkglibexec_PROGRAMS = artificial
if CAFMPI
pkglibexec_PROGRAMS += artificial-caf
endif

.NOTPARALLEL:

artificial_SOURCES = ftransform.f08 spray.f08 catio.f08 lcio.f08 generator.f08 \
                     artpicture.f95 lcparams.f08 lcurves.f08 marker.f08 \
                     galaxy.f08 cosmics.f08 flatmodule.f08 fitskeys.f08 \
                     artificial.f08
artificial_LDADD   = $(ldir)/libpsfmodel.a $(ldir)/libnoise.a $(ldir)/libastrotrafo.a \
                     $(ldir)/libphotometry.a $(ldir)/libastrosphere.a \
                     $(ldir)/libpppack.a $(ldir)/libtitsio.a \
                     $(CAF_LIBS) $(LIBCFITSIO) $(LIBM)

# The parallel variant (--mproc > 1, launched by mpirun): the main program,
# the only source with coarrays, compiled for the MPI coarray library and
# linked with the modules of artificial (built for one process as before).
artificial_caf_SOURCES = artificial.f08
artificial_caf_FCFLAGS = $(FC_MODINC). $(FC_MODINC)../lib $(CAFMPI_FCFLAGS)
artificial_caf_LDADD = ftransform.$(OBJEXT) spray.$(OBJEXT) catio.$(OBJEXT) \
                       lcio.$(OBJEXT) generator.$(OBJEXT) artpicture.$(OBJEXT) \
                       lcparams.$(OBJEXT) lcurves.$(OBJEXT) marker.$(OBJEXT) \
                       galaxy.$(OBJEXT) cosmics.$(OBJEXT) flatmodule.$(OBJEXT) \
                       fitskeys.$(OBJEXT) \
                       $(ldir)/libpsfmodel.a $(ldir)/libnoise.a $(ldir)/libastrotrafo.a \
                       $(ldir)/libphotometry.a $(ldir)/libastrosphere.a \
                       $(ldir)/libpppack.a $(ldir)/libtitsio.a \
                       $(CAFMPI_LIBS) $(LIBCFITSIO) $(LIBM)
# the modules are compiled for artificial
artificial_caf-artificial.$(OBJEXT): artificial.$(OBJEXT)
