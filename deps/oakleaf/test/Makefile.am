#
#  Copyright © 2018-2025 F.Hroch (hroch@physics.muni.cz)
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

AM_FCFLAGS = $(FC_MODINC)../src $(FC_MODINC)../src/fmm
LDADD = -L../src -loakleaf

#SUBDIRS = qr

# testphrate.sh

.NOTPARALLEL:

testrmean_BUILD = testmean_real32.f08
noise_BUILD = noise_real32.f08
suite_BUILD = suite_real32.f08
# testphrate

testmean_real32.f08: testmean.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
noise_real32.f08: noise.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@
suite_real32.f08: suite.f08
	$(SED) 's/REAL64/REAL32/g' < $< >$@

#testphrate_SOURCES = noise.f08 toolbox.f08 suite.f08 romberg.f08 testphrate.f08
testmean_SOURCES = noise.f08 toolbox.f08 suite.f08 testmean.f08
testmean_real32_SOURCES = noise_real32.f08 suite_real32.f08 testmean_real32.f08
testfmean_SOURCES = noise.f08 toolbox.f08 suite.f08 testfmean.f08
testrline_SOURCES = noise.f08 toolbox.f08 suite.f08 testrline.f08
#testsqps_SOURCES = testsqps.f08
#testqrhouse_SOURCES = testqrhouse.f08

TESTS = testmean.sh testmean_real32.sh testfmean.sh testrline.sh
BUILT_SOURCES = $(testrmean_BUILD) $(noise_BUILD) $(suite_BUILD)
check_PROGRAMS = testmean testmean_real32 testfmean testrline

if REAL128

testrmean_BUILD += testmean_real128.f08
noise_BUILD += noise_real128.f08
suite_BUILD += suite_real128.f08
testmean_real128_SOURCES = noise_real128.f08 suite_real128.f08 testmean_real128.f08

testmean_real128.f08: testmean.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
noise_real128.f08: noise.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@
suite_real128.f08: suite.f08
	$(SED) 's/REAL64/REAL128/g' < $< >$@

TESTS += testmean_real128.sh
check_PROGRAMS += testmean_real128
endif

EXTRA_DIST = $(TESTS)
CLEANFILES = $(testrmean_BUILD) $(noise_BUILD) $(suite_BUILD)
clean-local:
	rm -f *.out
