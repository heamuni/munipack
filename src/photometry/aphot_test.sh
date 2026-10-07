
if [ -f aphottester.f08 ]; then
    echo "Use another working place."
    exit 0
fi

PATH=$PATH:${0%aphot_test.sh}

if [ ! -f cone.fits ]; then
    munipack cone -r 0.1 -- 110.473 71.343
fi

HWHM=2
ECC=0.9
INCL=60

munipack artificial --psf gauss --hwhm $HWHM --disable-noise \
	 --eccentricity $ECC --inclination $INCL  --sky-mag 18 \
	 -c cone.fits --rcen 110.473 --dcen 71.343 --fov 0.1 \
	 --width 512 --height 384 --col-mag Vmag --mask art.fits

munipack find art.fits
munipack aphot --verbose --eccentricity $ECC --inclination $INCL  art.fits

aphottester $HWHM $ECC
