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
FITSADD = $(ldir)/libtitsio.a $(LIBCFITSIO)

pkglibexec_PROGRAMS = find aphot gphot phcal phfotran psfphot
noinst_PROGRAMS = aphottester

.NOTPARALLEL:

find_SOURCES = daosky.f95 daofind.f95 fitsfind.f08 find.f08
find_LDADD   = $(ldir)/libsort.a $(LIBOAKLEAF) $(LIBMINPACK) $(FITSADD)

aphot_SOURCES = daosky.f95 daofotometr.f08 fitsaphot.f08 exscatter.f08 \
                apee.f08 aphot.f08
aphot_LDADD   = $(ldir)/libpsfmodel.a $(LIBOAKLEAF) $(LIBMINPACK) $(FITSADD)

aphottester_SOURCES = daofotometr.f08 fitsaphot.f08 aphottester.f08
aphottester_LDADD   = $(LIBOAKLEAF) $(LIBMINPACK) $(FITSADD)

gphot_SOURCES = growinit.f95 growmodel.f95 growfit.f95\
                growreport.f95 growcurve.f95 fitsgphot.f08 gphot.f08
gphot_LDADD   = $(ldir)/liblmin.a $(ldir)/libpppack.a $(ldir)/libminpacks.a \
                $(ldir)/libfmm.a $(LIBOAKLEAF) $(LIBMINPACK) $(FITSADD)

phcal_SOURCES =  fitsfotran.f95 rational.f95 rplane.f95 \
                 fotran.f95 mfits.f08 calibre.f95 jam.f95 colorex.f95 \
                 sfits.f08 phcal.f08
phcal_LDADD = $(ldir)/liblmin.a $(ldir)/libsort.a $(ldir)/libastrotrafo.a \
              $(ldir)/libphotometry.a  $(WCSLIB_LIBS) \
              $(ldir)/libastrosphere.a $(LIBOAKLEAF) $(LIBMINPACK) $(FITSADD)

phfotran_SOURCES = mfits.f08 jam.f95 fitsfotran.f95 rplane.f95 \
                   rational.f95 fotran.f95 phfotran.f95
phfotran_LDADD   = $(ldir)/libastrotrafo.a $(ldir)/libphotometry.a \
                   $(ldir)/libminpacks.a $(ldir)/liblmin.a \
                   $(ldir)/libastrosphere.a \
                   $(LIBOAKLEAF) $(LIBMINPACK) $(FITSADD)
# $(ldir)/libphio.a

psfphot_SOURCES = fitspsfphot.f08 psffit.f08 psfempir.f08 psfphot.f08
psfphot_LDADD   = $(ldir)/libpsfmodel.a $(FITSADD)
# the module dependencies (for parallel builds)
psffit.$(OBJEXT): fitspsfphot.$(OBJEXT)
psfempir.$(OBJEXT): psffit.$(OBJEXT)
psfphot.$(OBJEXT): psfempir.$(OBJEXT)

dist_pkgdata_DATA = photosystems.lst
dist_html_DATA = grow_report.sh grow_report.py

LSD = "MODE = 'RESTORE'\nFILE = 'photosystems.lst' 'photosystems.fits'\n"
install-data-hook: photosystems.lst
	cd $(DESTDIR)$(pkgdatadir) && \
           echo -e $(LSD) | $(DESTDIR)$(pkglibexecdir)/fits

uninstall-local:
	rm -f $(DESTDIR)$(pkgdatadir)/photosystems.fits
