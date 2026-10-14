
#
#  Copyright © 2025 F.Hroch (hroch@physics.muni.cz)
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

# switch-on for a debuging
set -x

# I lost 6h of investigation how to provide the test data
# for tests of `make distcheck', as I overloooked a small notice
# in the end of the page `Scripts-based Testsuites' of the automake manual.
# The `${srcdir}' is the helpfull variable.
#
# This check script lye inside of the automake machinery. All utilities
# are uninstalled yet, their using requires update of the libexec, so
# it will works only as the test driver.

echo "srcdir:" ${srcdir}
PATH="$PATH:../munipack"; export PATH
MUNIPACK_DATA_DIR="${srcdir}/../vo"; export MUNIPACK_DATA_DIR

FILE=barnard-calibrated.fits
CONE=cone_barnard.fits

test -f ${FILE} && test -f ${CONE}
if test $? -eq 1; then
    FILES_AVAILABLE="no"
    cp ${srcdir}/${FILE} ${srcdir}/${CONE} .
fi

# exit imediatelly, when any command fails
set -e

MUNIPACK_LIBEXEC_PATH="../photometry"; export MUNIPACK_LIBEXEC_PATH
munipack find -f 5 $FILE
munipack aphot $FILE

# munipack cone -r 0.1 --par 'Vmag=<15' -o cone_barnard.fits -- 269.44 4.69
MUNIPACK_LIBEXEC_PATH="../astrometry"; export MUNIPACK_LIBEXEC_PATH
munipack astrometry -c $CONE $FILE

cp ${srcdir}/../photometry/photosystems.lst .
echo -e "MODE = 'RESTORE'\nFILE = 'photosystems.lst' 'photosystems.fits'\n" |\
    ../fits/fits

MUNIPACK_LIBEXEC_PATH="../photometry"; export MUNIPACK_LIBEXEC_PATH
munipack phcal --photsys-ref Johnson --area 0.3 \
         --col-ra RAJ2000 --col-dec DEJ2000 \
         -f V --col-mag Vmag --col-magerr e_Vmag -c $CONE \
         -O --mask phcal.fits \
	 --phsystab photosystems.fits $FILE

MUNIPACK_LIBEXEC_PATH="../fits"; export MUNIPACK_LIBEXEC_PATH
# CTPH within twice its error of 11.5 (the error is ~1.3 on this frame:
# a window of 11 .. 12 failed on legitimate changes of the star selection)
CTPH=`munipack fits -K CTPH --value phcal.fits\[PHOTOMETRY\]`
CTPHERR=`munipack fits -K CTPHERR --value phcal.fits\[PHOTOMETRY\]`
echo "CTPH = $CTPH +- $CTPHERR"
STATUS=`echo $CTPH $CTPHERR | awk '{ d = $1 - 11.5; if(d < 0) d = -d; if(0 < $2 && $2 < 3 && d < 2*$2) {print "OK";}}'`

CLEAN_FILES="phcal.fits gmon.out photosystems.lst photosystems.fits"
if test x$FILES_AVAILABLE = "xno"; then
    CLEAN_FILES="$CLEAN_FILES $FILE $CONE"
fi
rm -f $CLEAN_FILES

test x$STATUS = xOK && exit 0
