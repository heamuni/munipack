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

# ##
# process #include directores in .shtml files
# in fashion of Apache's SSI includes
#
# the tags: <!-- #include virtual="*" --> are replaced by content of a file
#

SUFFIXES: .shtml .html

TARGETS = munipack.html features.html download.html docs.html notfound.html \
         basic.html external.html classics.html colourspace.html \
         colourfits.html itone.html guide.html faq.html vobs.html history.html \
         modules.html credits.html pererr.html lctut.html kombitut.html \
         phcorrtut.html colourimage.html colorcal.html colcal.html cmd.html \
         chart.html \
         man_intro.html man_help.html man_colouring.html \
         man_bias.html man_dark.html man_flat.html man_phcorr.html \
         man_cone.html man_astrometry.html man_phcal.html man_aphot.html \
         man_gphot.html man_psfphot.html \
         man_list.html man_fits.html man_votable.html man_com.html  \
         man_env.html man_kombine.html man_timeseries.html man_setup.html \
         man_backup.html man_i18n.html \
         dataform_tmseries.html dataform_photometry.html \
         dataform_astrometry.html dataform_artgalaxy.html dataform_artpsf.html \
         dataform_photosys.html dataform_phfotran.html dataform_proc.html \
         dataform_biasdark.html dataform_flat.html mosaics.html \
         phcalibration.html phoverview.html astoverview.html deka154.html \
         man_find.html man_phfotran.html man_cross.html fwhm.html changes.html \
         appimage.html install.html version.html SourceInstallation.html \
         grow.html dataform_growrep.html man_artificial.html artific.html \
         shcripts.html flatfielding.html man_sesame.html

PADDING = munipack.css favicon.ico henden_m67.sh m67.sh m27.sh \
          robots.txt 0716+71.lst timeserie_cat.lst lctut.sh ucac_jmuc.py \
          artificial.sh kombitut.sh mosaic.sh Kepler-21b.lst fourier.lst \
          artgalaxy.lst artpsf.lst artpsf.f08

EXTRA_DIST = munipack.shtml features.shtml download.shtml docs.shtml notfound.shtml \
         basic.shtml external.shtml classics.shtml colourspace.shtml \
         colourfits.shtml itone.shtml guide.shtml faq.shtml vobs.shtml \
         history.shtml \
         modules.shtml credits.shtml pererr.shtml lctut.shtml kombitut.shtml \
         phcorrtut.shtml colourimage.shtml colorcal.shtml colcal.shtml cmd.shtml \
         chart.shtml \
         man_intro.shtml man_help.shtml man_colouring.shtml \
         man_bias.shtml man_dark.shtml man_flat.shtml man_phcorr.shtml \
         man_cone.shtml man_astrometry.shtml man_phcal.shtml man_aphot.shtml \
         man_gphot.shtml man_psfphot.shtml \
         man_list.shtml man_fits.shtml man_votable.shtml man_com.shtml  \
         man_env.shtml man_kombine.shtml man_timeseries.shtml man_setup.shtml \
         man_backup.shtml man_i18n.shtml \
         dataform_tmseries.shtml dataform_photometry.shtml \
         dataform_astrometry.shtml dataform_artgalaxy.shtml dataform_artpsf.shtml \
         dataform_photosys.shtml dataform_phfotran.shtml dataform_proc.shtml \
         dataform_biasdark.shtml dataform_flat.shtml mosaics.shtml \
         phcalibration.shtml phoverview.shtml astoverview.shtml deka154.shtml \
         man_find.shtml man_phfotran.shtml man_cross.shtml fwhm.shtml changes.shtml \
         appimage.shtml install.shtml version.shtml SourceInstallation.shtml \
         grow.shtml dataform_growrep.shtml man_artificial.shtml artific.shtml \
         shcripts.shtml flatfielding.shtml man_sesame.shtml \
         meta.shtml foot.shtml head.shtml include.awk

include image_list.mk

dist_html_DATA = $(TARGETS) $(image_list) $(PADDING)

.shtml.html:
	$(AWK) -f include.awk $< > $@

#install-exec-hook:
#	cd $(docdir) && $(LN_S) munipack.html index.html
# We got a problem when the link already exists.

CLEANFILES = $(TARGETS)
