#
#  Copyright © 2018-2026 F.Hroch (hroch@physics.muni.cz)
#
#  This file is part of Oakleaf.
#
#  Oakleaf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Oakleaf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
#

SUBDIRS = fmm minpack simplex
AM_FCFLAGS = $(FC_MODINC)fmm $(FC_MODINC)minpack $(FC_MODINC)simplex

#include_HEADERS = oakleaf.${FC_MODEXT}
install-exec-hook: oakleaf.${FC_MODEXT}
	$(MKDIR_P) $(DESTDIR)$(includedir) \
           && $(install_sh) oakleaf.${FC_MODEXT} $(DESTDIR)$(includedir)
uninstall-hook:
	rm -f $(DESTDIR)$(includedir)/oakleaf.${FC_MODEXT}


lib_LTLIBRARIES = liboakleaf.la

.NOTPARALLEL:

rfun_hampel_BUILD = rfun_hampel_real32.f08 rfun_hampel_real64.f08
rfun_tukey_BUILD = rfun_tukey_real32.f08 rfun_tukey_real64.f08
rfun_huber_BUILD = rfun_huber_real32.f08 rfun_huber_real64.f08
rnfun_square_BUILD = rnfun_square_real32.f08 rnfun_square_real64.f08
rlike_BUILD = rlike_real32.f08 rlike_real64.f08
q2eig_BUILD = q2eig_real32.f08 q2eig_real64.f08
qsort_BUILD = qsort_real32.f08 qsort_real64.f08
qmed_BUILD = qmed_real32.f08 qmed_real64.f08
qmean_BUILD = qmean_real32.f08 qmean_real64.f08
fmin_BUILD = fmin_real32.f08 fmin_real64.f08
qtrstep_BUILD = qtrstep_real32.f08 qtrstep_real64.f08
sqpfmean_BUILD = sqpfmean_real32.f08 sqpfmean_real64.f08
rmean_BUILD = rmean_real32.f08 rmean_real64.f08
median_BUILD = median_real32.f08 median_real64.f08

rfun_hampel_real32.f08: rfun_hampel.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
rfun_tukey_real32.f08: rfun_tukey.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
rfun_huber_real32.f08: rfun_huber.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
rnfun_square_real32.f08: rnfun_square.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
rlike_real32.f08: rlike.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
qsort_real32.f08: qsort.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
qmed_real32.f08: qmed.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
qmean_real32.f08: qmean.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
fmin_real32.f08: fmin.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
q2eig_real32.f08: q2eig.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
qtrstep_real32.f08: qtrstep.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
sqpfmean_real32.f08: sqpfmean.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
rmean_real32.f08: rmean.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
median_real32.f08: median.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@

rfun_hampel_real64.f08: rfun_hampel.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
rfun_tukey_real64.f08: rfun_tukey.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
rfun_huber_real64.f08: rfun_huber.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
rnfun_square_real64.f08: rnfun_square.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
rlike_real64.f08: rlike.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
qsort_real64.f08: qsort.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
qmed_real64.f08: qmed.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
qmean_real64.f08: qmean.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
fmin_real64.f08: fmin.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
q2eig_real64.f08: q2eig.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
qtrstep_real64.f08: qtrstep.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
sqpfmean_real64.f08: sqpfmean.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
rmean_real64.f08: rmean.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@
median_real64.f08: median.f08
	$(SED) 's/REAL64/REAL64/g' < $< >$@

if REAL128

AM_FCFLAGS += -DHAVE_REAL128
rfun_hampel_BUILD += rfun_hampel_real128.f08
rfun_tukey_BUILD += rfun_tukey_real128.f08
rfun_huber_BUILD += rfun_huber_real128.f08
rnfun_square_BUILD += rnfun_square_real128.f08
rlike_BUILD += rlike_real128.f08
q2eig_BUILD += q2eig_real128.f08
qsort_BUILD += qsort_real128.f08
qmed_BUILD += qmed_real128.f08
qmean_BUILD += qmean_real128.f08
fmin_BUILD += fmin_real128.f08
qtrstep_BUILD += qtrstep_real128.f08
sqpfmean_BUILD += sqpfmean_real128.f08
rmean_BUILD += rmean_real128.f08
median_BUILD += median_real128.f08

rfun_hampel_real128.f08: rfun_hampel.f08
	$(SED) 's/REAL64/REAL128/g'  < $< >$@
rfun_tukey_real128.f08: rfun_tukey.f08
	$(SED) 's/REAL64/REAL128/g'  < $< >$@
rfun_huber_real128.f08: rfun_huber.f08
	$(SED) 's/REAL64/REAL128/g'  < $< >$@
rnfun_square_real128.f08: rnfun_square.f08
	$(SED) 's/REAL64/REAL128/g'  < $< >$@
rlike_real128.f08: rlike.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
qsort_real128.f08: qsort.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
qmed_real128.f08: qmed.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
qmean_real128.f08: qmean.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
fmin_real128.f08: fmin.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
q2eig_real128.f08: q2eig.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
qtrstep_real128.f08: qtrstep.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
sqpfmean_real128.f08: sqpfmean.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
rmean_real128.f08: rmean.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
median_real128.f08: median.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@

endif


BUILT_SOURCES = $(qsort_BUILD) $(qmed_BUILD) $(qmean_BUILD) $(fmin_BUILD) \
                $(rfun_hampel_BUILD) $(rfun_BUILD) $(rlike_BUILD) $(q2eig_BUILD) \
                $(qtrstep_BUILD) $(sqpfmean_BUILD) $(rmean_BUILD) $(median_BUILD)

liboakleaf_la_SOURCES = $(rfun_hampel_BUILD) $(rfun_tukey_BUILD) \
                        $(rfun_huber_BUILD) $(rnfun_square_BUILD) rfun.f08 \
			$(fmin_BUILD)  \
                        $(rlike_BUILD) \
                        $(q2eig_BUILD)  \
                        $(qtrstep_BUILD) \
                        $(sqpfmean_BUILD)  \
                        $(qsort_BUILD) quicksort.f08 \
                        $(qmed_BUILD)  quickmedian.f08 \
                        $(qmean_BUILD) quantilmean.f08  \
                        $(median_BUILD) medians.f08 \
                        $(rmean_BUILD)  \
                        iscale.f08 fmean.f08 \
                        cholesky.f08 qrhouse.f08 \
                        rline.f08 \
                        bfgs.f08 huberfit.f08 \
                        oakleaf.f08
liboakleaf_la_LIBADD = fmm/libfmm.la simplex/libsimplex.la \
                       minpack/libminpack90.la $(LIBS)

EXTRA_DIST = rfun_hampel.f08 rfun_tukey.f08 rfun_huber.f08 \
        rnfun_square.f08 qsort.f08 qmed.f08 qmean.f08 median.f08 \
	rlike.f08 fmin.f08 q2eig.f08 qtrstep.f08 sqpfmean.f08 rmean.f08

# https://autotools.info/libtool/version.html
liboakleaf_la_LDFLAGS = -version-info 2:0:0

#  trust_step.f08 sr1.f08 rmean_rfun.f08 trmean.f08 rmean2.f08
#   rmean.f08 rwmean.f08 zeroin.f08 madmean.f08
# sqps.f08 phcommon.f08 phnormal.f08 phinit.f08 phkmean.f08
# phopt.f08 phiopt.f08 phline.f08 phrate.f08

CLEANFILES = $(qsort_BUILD) $(qmed_BUILD) $(qmean_BUILD) $(fmin_BUILD) \
             $(rfun_hampel_BUILD) $(rfun_BUILD) $(rlike_BUILD) $(q2eig_BUILD) \
             $(qtrstep_BUILD) $(sqpfmean_BUILD) $(rmean_BUILD) $(median_BUILD)

clean-local:
	rm -f *.$(FC_MODEXT)
