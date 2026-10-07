

echo "Processing data of M27 to get colour frame..."

wget ftp://integral.physics.muni.cz/pub/munipack/munipack-data-m27.tar.gz
tar zxf munipack-data-m27.tar.gz

DATADIR=$(pwd)/munipack-data-m27/

mkdir workdir-m27/
cd workdir-m27/

munipack dark -o d60.fits $DATADIR/d60_*.fits
munipack dark -o d40.fits $DATADIR/d40_*.fits

OPT="--enable-overwrite"
munipack phcorr -dark d60.fits $OPT -t . $DATADIR/m27_*B.fits
munipack phcorr -dark d40.fits $OPT -t . $DATADIR/m27_*[VR].fits

munipack find -f 5 m27_*.fits
munipack aphot m27_*.fits

munipack cone --Johnson-patch -r 0.1 --magmin 10 --magmax 13 -- 299.87 22.71

munipack astrometry m27_*B.fits &
munipack astrometry m27_*V.fits &
munipack astrometry m27_*R.fits &
wait

munipack kombine -o m27_B.fits --rcen 299.9 --dcen 22.72 m27_*B.fits &
munipack kombine -o m27_V.fits --rcen 299.9 --dcen 22.72 m27_*V.fits &
munipack kombine -o m27_R.fits --rcen 299.9 --dcen 22.72 m27_*R.fits &
wait

munipack find -f 5 m27_?.fits
munipack aphot m27_?.fits

for F in B V R; do
    munipack phcal --photsys-ref Johnson --area 0.2 -c cone.fits \
	     -f $F --col-mag ${F}mag --col-magerr e_${F}mag m27_${F}.fits
done

munipack colouring -o m27.fits m27_B_cal.fits m27_V_cal.fits m27_R_cal.fits
