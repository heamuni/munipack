
# prepare data
wget ftp://munipack.physics.muni.cz/pub/munipack/munipack-data-bubble.tar.gz
tar zxf munipack-data-bubble.tar.gz

DATADIR=$(pwd)/munipack-data-bubble/

mkdir workdir-bubble/
cd workdir-bubble/

# precorrections
munipack dark $DATADIR/d120_*.fits
munipack phcorr -gain 2.3 -dark dark.fits --enable-overwrite -t . $DATADIR/ngc7635_*.fits

# find stars, aperture photometry
munipack find ngc7635_*.fits
munipack aphot ngc7635_*.fits

# astrometry
munipack cone -r 0.1 350.20 61.20
munipack astrometry ngc7635_*.fits

# kombine
munipack kombine -o bubble.fits ngc7635_*.fits
