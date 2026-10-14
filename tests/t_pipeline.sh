# The main route on a simulated field with a known truth: artificial ->
# find -> aphot -> astrometry (from scratch) -> phcal, and psfphot -> phcal.
# 150 stars V = 10-16 at Dec +30, Gaussian PSF (sigma 1.5 px), a detector
# with gain, read noise, bias and dark.
# Covers: HWHM is a true HWHM (BUGS P8), RA residuals with cos(dec) (S7),
# apertures not invalidated by isolated neighbours (P29), errors by the
# gain (P34), calibrated MAG against the catalogue, psfphot (P18, P35, P42).

mkcat cat.fits 150 11 180 30 0.08 10 16 0.004
munipack artificial -c cat.fits $CATCOLS -f V --rcen 180 --dcen 30 \
    --scale 3600 --width 640 --height 640 --exptime 60 --psf GAUSS --hwhm 1.5 \
    --gain 2 --readnoise 5 --bias 100 --dark 0.01 --seed 2 --no-watermark \
    --mask frame.fits > /dev/null 2>&1

munipack find -f 3.5 frame.fits > /dev/null 2>&1
check "stars found" "$(nrows 'frame.fits[FIND]')" 145 155

munipack aphot frame.fits > /dev/null 2>&1
# the true HWHM is 1.1774 sigma = 1.766 px (sigma 1.5 before P8), pixelised
check "APERPHOT HWHM [px]" "$(key 'frame.fits[APERPHOT]' HWHM)" 1.65 1.95

# Astrometry from scratch: the WCS of artificial removed.
cp frame.fits psf.fits
munipack astrometry --remove frame.fits > /dev/null 2>&1
[ -z "$(key frame.fits CRVAL1)" ] || fail "astrometry --remove left CRVAL1"
munipack astrometry -c cat.fits $CATCOLS frame.fits > /dev/null 2>&1
read -r cr1 cr2 px1 px2 c11 c12 c21 c22 <<< "$(munipack fits -K \
    CRVAL1,CRVAL2,CRPIX1,CRPIX2,CD1_1,CD1_2,CD2_1,CD2_2 --value frame.fits | tr '\n' ' ')"
# at CRPIX = 320, 320 the truth is 180, +30; 0.05 px = 1.4e-5 deg
check "CRPIX1" "$px1" 320 320
check "CRVAL1 offset [px]" "$(awk -v v="$cr1" 'BEGIN { print (v - 180) * 3600 * cos(30 * atan2(0, -1) / 180) }')" -0.05 0.05
check "CRVAL2 offset [px]" "$(awk -v v="$cr2" 'BEGIN { print (v - 30) * 3600 }')" -0.05 0.05
check "scale [pix/deg]" "$(awk -v a="$c11" -v b="$c12" -v c="$c21" -v d="$c22" \
    'BEGIN { x = a * d - b * c; print 1 / sqrt(x < 0 ? -x : x) }')" 3599.5 3600.5
check "rotation CD1_2/CD2_2" "$(awk -v b="$c12" -v d="$c22" 'BEGIN { print b / d }')" -3e-4 3e-4
check "CD1_1 < 0 (east left)" "$(awk -v a="$c11" 'BEGIN { print (a < 0) }')" 1 1

# Calibration against the catalogue (the frame is rendered from V exactly).
pcal() {
    munipack phcal --photsys-ref Johnson -f V --col-ra RAJ2000 --col-dec DEJ2000 \
        --col-mag Vmag --col-magerr e_Vmag -c cat.fits "$@" > /dev/null 2>&1
}
calcheck() {   # FILE LABEL
    match "cat.fits[1]" "RAJ2000 DEJ2000 Vmag" "$1[PHOTOMETRY]" \
        "RAJ2000 DEJ2000 MAG MAGERR" 3e-4 > m.txt
    check "$2: stars matched" "$(wc -l < m.txt)" 145 150
    # all valid but the faintest few (a 16 mag star at S/N ~ 50 is valid)
    check "$2: valid MAG" "$(awk '$2 < 99 && $3 > 0' m.txt | wc -l)" 145 150
    awk '$2 < 99 && $3 > 0' m.txt > v.txt
    read -r med sd n <<< "$(awk '{ print $2 - $1 }' v.txt | stats)"
    check "$2: median MAG - V" "$med" -0.003 0.003
    check "$2: robust sd MAG - V" "$sd" 0 0.005
    # the reported errors describe the scatter
    read -r med sd n <<< "$(awk '{ print ($2 - $1) / $3 }' v.txt | stats)"
    check "$2: robust sd of (MAG - V)/MAGERR" "$sd" 0.7 1.4
}
pcal frame.fits
calcheck frame_cal.fits "aphot"

# PSF photometry on the same frame (the WCS of artificial kept).
munipack psfphot psf.fits > /dev/null 2>&1
check "psfphot stars" "$(nrows 'psf.fits[PSFPHOT]')" 145 155
pcal --photometry psf psf.fits
check "PHOTOTYP" "$(key 'psf_cal.fits[PHOTOMETRY]' PHOTOTYP | grep -c PSF)" 1 1
calcheck psf_cal.fits "psfphot"
