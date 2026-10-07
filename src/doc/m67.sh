

#wget https://integral.physics.muni.cz/ftp/munipack/munipack-data-m67.tar.gz
tar zxf munipack-data-m67.tar.gz

DATADIR=$(pwd)/munipack-data-m67/

mkdir workdir-m67/
cd workdir-m67/

# darks
munipack dark -o d7.fits $DATADIR/d7_*.fits &
munipack dark -o d30.fits $DATADIR/d30_*.fits &
wait

# flats
munipack flat -o fB.fits -dark d7.fits $DATADIR/flat_*B.fits &
munipack flat -o fV.fits -dark d7.fits $DATADIR/flat_*V.fits &
munipack flat -o fR.fits -dark d7.fits $DATADIR/flat_*R.fits &
wait

# precorr
OPT="--enable-overwrite"
munipack phcorr -dark d30.fits -flat fB.fits $OPT -t . $DATADIR/m67_*B.fits &
munipack phcorr -dark d30.fits -flat fV.fits $OPT -t . $DATADIR/m67_*V.fits &
munipack phcorr -dark d30.fits -flat fR.fits $OPT -t . $DATADIR/m67_*R.fits &
wait

# find, aperture photometry
munipack find -th 10 -f 6 m67_*.fits
munipack aphot m67_*.fits

# astrometry
RA=132.8304
DEC=11.7771
munipack cone -r 0.15 --magmax 15 --Johnson-patch -- $RA $DEC

munipack astrometry -c cone.fits m67_??B.fits &
munipack astrometry -c cone.fits m67_??V.fits &
munipack astrometry -c cone.fits m67_??R.fits &
wait

# collect
RADEC="--rcen $RA --dcen $DEC"
munipack kombine $RADEC -o m67_B.fits m67_??B.fits &
munipack kombine $RADEC -o m67_V.fits m67_??V.fits &
munipack kombine $RADEC -o m67_R.fits m67_??R.fits &
wait

# find, aperture photometry
munipack find -th 10 -f 6 m67_?.fits
munipack aphot m67_?.fits

# photometry calibration by UCAC4
munipack phcal --photsys-ref Johnson -f B --area 0.3 --col-mag Bmag \
	 --col-magerr e_Bmag -c cone.fits m67_B.fits
munipack phcal --photsys-ref Johnson -f V --area 0.3 --col-mag Vmag \
	 --col-magerr e_Vmag -c cone.fits m67_V.fits
munipack phcal --photsys-ref Johnson -f R --area 0.3 --col-mag Rmag \
	 --col-magerr e_Rmag -c cone.fits m67_R.fits
