#!/bin/sh

# This script creates light curves and other related
# products by using Munipack data for 0716+71
#
# See detailed description in lctut.html

# data
wget -nv ftp://munipack.physics.muni.cz/pub/munipack/munipack-data-0716.tar.gz
tar zxf munipack-data-0716.tar.gz
DATADIR=$(pwd)/munipack-data-0716/

mkdir workdir-0716/
cd workdir-0716/

# preprocessing
munipack dark -o d30.fits $DATADIR/d30_*.fits
munipack dark -o d120.fits $DATADIR/d120_*.fits
munipack flat -o f_V.fits -dark d30.fits $DATADIR/f30_*V.fits
munipack flat -o f_R.fits -dark d30.fits $DATADIR/f30_*R.fits

OPT="--enable-overwrite"
munipack phcorr -dark d120.fits -flat f_V.fits $OPT -t . $DATADIR/0716_*V.fits
munipack phcorr -dark d120.fits -flat f_R.fits $OPT -t . $DATADIR/0716_*R.fits

# aperture photometry and star detection
munipack find -f 3 0716_*.fits
munipack aphot 0716_*.fits

# astrometry calibration
munipack cone --Johnson-patch --magmin 11 --magmax 15 -r 0.2 -- 110.5 71.3
munipack astrometry -c cone.fits 0716_*.fits

# photometry calibration
wget -nv http://munipack.physics.muni.cz/0716+71.lst
munipack fits --restore 0716+71.lst

# SAPER=5 minimises dispersion of A comparison star
# SAPER=3 minimises dispersion of the blasar: PSF is actually need.
SAPER="--saper 5"
#THRESH="-th 0.1 -e 0.01"
for FILTER in V R; do
    munipack kombine --disable-back -o 0716.${FILTER}.fits \
	     0716_04[0-7]${FILTER}.fits
  #cp 0716_042${FILTER}.fits 0716.${FILTER}.fits
  munipack find -f 3 0716.${FILTER}.fits
  munipack aphot 0716.${FILTER}.fits
  munipack astrometry -c cone.fits 0716.${FILTER}.fits
  munipack phcal --photsys-ref Johnson -f ${FILTER} --area 0.3 \
	   -c 0716+71.fits --col-mag ${FILTER}mag --col-magerr e_${FILTER}mag \
	   $THRESH $SAPER 0716.${FILTER}.fits
done


for FILTER in V R; do

  # catalogue
  for A in 0716_*${FILTER}.fits; do
      munipack phcal \
        --photsys-ref Johnson --area 0.3 \
        -f ${FILTER} --col-mag ${FILTER}mag --col-magerr e_${FILTER}mag  \
        -c 0716+71.fits \
        $THRESH $SAPER  \
	-O --mask '\!\1_catcal.\2' $A;
  done

  # reference frame
  for A in 0716_*${FILTER}.fits; do
      munipack phcal \
        --photsys-ref Johnson -f ${FILTER} --area 0.3 \
	$THRESH $SAPER -r 0716.${FILTER}_cal.fits \
        -O --mask '\!\1_refcal.\2' $A;
  done

  # ucac
  for A in 0716_*${FILTER}.fits; do
      munipack phcal \
        --photsys-ref Johnson --area 0.3 --col-ra RAJ2000 --col-dec DEJ2000 \
        -f ${FILTER} --col-mag ${FILTER}mag --col-magerr e_${FILTER}mag  \
	-c cone.fits $THRESH $SAPER \
        -O --mask '\!\1_ucacal.\2' $A;
  done

  # manual
  for A in 0716_*${FILTER}.fits; do
    munipack phcal -C 1 --photsys-ref Johnson -f ${FILTER} \
        $SAPER -O --mask '\!\1_mancal.\2' $A;
  done

done

# light curves
for FILTER in V R; do
   for TYPE in man cat ref uca; do
      munipack timeseries -c "110.473,71.343 110.389,71.322 110.468,71.305" \
               -l MAG,MAGERR --stdout  0716_*${FILTER}_${TYPE}cal.fits | \
	  sort -n > ${TYPE}cal_${FILTER}
   done
done


# plots, V filter
gnuplot <<EOF
set key bottom left
set title "0716+71"
set ylabel "magnitude in V"
set xlabel "JD - 2453409"
set yrange[14.8:13.5]
set term svg dynamic
set output 'lc0716_V.svg'
plot 'mancal_V' u (\$1-2453409):(\$2-\$4+11.8) t "diff(0716+71 - A) + 11.8", \
     'catcal_V' u (\$1-2453409):(\$2-0.2) t "std.field - 0.2" ls 1 lc 2, \
     'refcal_V' u (\$1-2453409):2 t "frame" ls 1 lc 3, \
     'ucacal_V' u (\$1-2453409):(\$2+0.2) t "UCAC4 + 0.2" ls 1 lc 4
EOF


gnuplot <<EOF
set key top left
set title "Calibration star A"
set ylabel "magnitude in V"
set xlabel "JD - 2453409"
set yrange[12:11]
set ytics 0.1
set term svg dynamic
set output 'comp0716_V.svg'
plot 'mancal_V' u (\$1-2453409):(\$4-\$6+12.65) t "diff(A-B)+12.65", \
     'catcal_V' u (\$1-2453409):(\$4-0.1) t "std.field - 0.1" ls 1 lc 2, \
     'refcal_V' u (\$1-2453409):4 t "frame" ls 1 lc 3, \
     'ucacal_V' u (\$1-2453409):(\$4+0.2) t "UCAC4 + 0.2" ls 1 lc 4
EOF

# plots, R filter
gnuplot <<EOF
set key bottom left
set title "0716+71"
set ylabel "magnitude in R"
set xlabel "JD - 2453409"
set yrange[14.5:13.0]
set term svg dynamic
set output 'lc0716_R.svg'
plot 'mancal_R' u (\$1-2453409):(\$2-\$4+11.7) t "diff(0716+71 - A) + 11.7", \
     'catcal_R' u (\$1-2453409):(\$2-0.2) t "std.field - 0.2" ls 1 lc 2, \
     'refcal_R' u (\$1-2453409):2 t "frame" ls 1 lc 3, \
     'ucacal_R' u (\$1-2453409):(\$2+0.1) t "UCAC4 + 0.1" ls 1 lc 4
EOF


gnuplot <<EOF
set key top left
set title "Calibration star A"
set ylabel "magnitude in R"
set xlabel "JD - 2453409"
set yrange[11.7:10.7]
set ytics 0.1
set term svg dynamic
set output 'comp0716_R.svg'
plot 'mancal_R' u (\$1-2453409):(\$4-\$6+12.5) t "diff(A-B+12.5)", \
     'catcal_R' u (\$1-2453409):(\$4-0.1) t "std.field - 0.1" ls 1 lc 2, \
     'refcal_R' u (\$1-2453409):4 t "frame" ls 1 lc 3, \
     'ucacal_R' u (\$1-2453409):(\$4+0.1) t "UCAC4 + 0.1" ls 1 lc 4
EOF
