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

SUFFIXES: .shtml .html

PAGES = oakleaf.html fortran.html rmean.html order.html build.html \
        status.html changes.html bfgs.html huberfit.html

EXTRA_DIST = oakleaf.shtml fortran.shtml rmean.shtml order.shtml build.shtml \
        status.shtml changes.shtml bfgs.shtml huberfit.shtml \
        meta.shtml foot.shtml head.shtml render.awk

SUPPLEMENT = oakleaf.css favicon.ico
#example_rmean.f08 example_fmean.f08

dist_html_DATA = $(PAGES) $(SUPPLEMENT)

.shtml.html:
	$(AWK) -f render.awk $< > $@

CLEANFILES = $(PAGES)
