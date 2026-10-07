
# Examples for both artificial tutorial and manual page

# Cover Image -- Blazar 0716+71
munipack cone -r 0.3 -- 110.5 71.35
munipack artificial --qeff 0.02 -c cone.fits --exptime 120 \
	 --hwhm 0.6 --diameter 0.6 --sky-grad-x 0.003 \
         --rcen 110.476 --dcen 71.346 --scale 5.658e-4 --width 510 --height 360 \
	 --verbose
fitspng -o 0716+71.png artificial.fits

# merged by hand with 0716_035R.fits
