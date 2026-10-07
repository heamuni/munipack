
# prepare data
wget ftp://integral.physics.muni.cz/pub/munipack/munipack-data-m51.tar.gz
tar zxf munipack-data-m51.tar.gz

DATADIR=$(pwd)/munipack-data-m51/

mkdir workdir-m51/
cd workdir-m51/

# photometric corrections
export FITS_KEY_TEMPERATURE="CCD-TEMP"
munipack dark -o dark120.fits $DATADIR/dark_005?.fits $DATADIR/dark_008?.fits
munipack dark -o dark10.fits $DATADIR/dark_003?.fits $DATADIR/dark_004?.fits
munipack flat -o flat_Green.fits -dark dark10.fits $DATADIR/flat_Green_*.fits
munipack phcorr -dark dark120.fits -flat flat_Green.fits \
	 --enable-overwrite  -t . $DATADIR/m51_Green_*.fits

# find stars, aperture photometry
munipack find -f 7 -th 7 m51_Green_*.fits
munipack aphot m51_Green_*.fits

# astrometry
RA="202.47"
DEC="+47.2"
munipack cone -r 0.2 --magmax 18 -- $RA $DEC
munipack astrometry m51_Green_00*.fits

# kombine
munipack kombine --rcen $RA --dcen $DEC --width 1000 --height 1000 \
	 -o m51.fits m51_Green_00??.fits

# update filters, photometric calibration
for A in m51_Green_00??.fits; do
    munipack fits --update --key FILTER --val "'V'" $A
    munipack phcal --photsys-ref Johnson --area 0.3 -c cone.fits \
             --col-ra RAJ2000 --col-dec DEJ2000 \
	     -f V --col-mag Vmag --col-magerr e_Vmag $A
done

# kombine of calibrated frames
munipack kombine --rcen $RA --dcen $DEC --width 1000 --height 1000 \
	 -o m51_cal.fits m51_Green_00??_cal.fits
