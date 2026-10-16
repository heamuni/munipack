# The empirical PSF (psfphot --psf EMPIRICAL) on a simulated field with a
# known truth: 150 stars V = 10-16, a Moffat PSF (beta 2.5) with wide
# wings, the core fitted, the radial variation, then phcal.
# Covers: the table radius by the growth curve and the normalisation to
# an aperture (BUGS P50), the radial terms (P51), PSFMODERR in the errors
# (P52), the cross-validated regularisation (P53).

mkcat cat.fits 150 13 180 30 0.08 10 16 0.006
munipack artificial -c cat.fits $CATCOLS -f V --rcen 180 --dcen 30 \
    --scale 3600 --width 640 --height 640 --exptime 60 --psf MOFFAT \
    --hwhm 2.5 --beta 2.5 --gain 2 --readnoise 5 --bias 100 --dark 0.01 \
    --seed 5 --no-watermark --mask frame.fits > /dev/null 2>&1
munipack find -f 3.5 frame.fits > /dev/null 2>&1
munipack aphot frame.fits > /dev/null 2>&1

munipack psfphot --psf EMPIRICAL --psf-var-radial 1 frame.fits > /dev/null 2>&1
check "psfphot stars" "$(nrows 'frame.fits[PSFPHOT]')" 145 155
check "PSFVARR" "$(key 'frame.fits[PSFPHOT]' PSFVARR)" 1 1
# the radius of 99.5% of the flux of beta 2.5 is 5.1 FWHM (4 FWHM before)
check "PSFRAD / FWHM" "$(awk -v r="$(key 'frame.fits[PSFPHOT]' PSFRAD)" \
    -v f="$(key 'frame.fits[PSFPHOT]' FWHM)" 'BEGIN { print r / f }')" 4.5 5.8
check "PSFNORMN (stars of the normalisation)" "$(key 'frame.fits[PSFPHOT]' PSFNORMN)" 10 100
check "PSFSMOOT" "$(key 'frame.fits[PSFPHOT]' PSFSMOOT)" 0.001 10
[ -n "$(colnum 'frame.fits[PSFPHOT]' PSFMODERR)" ] || fail "no PSFMODERR column"
ok "PSFMODERR column"
# PSFMODERR is a part of PSFCOUNTERR
check "PSFMODERR <= PSFCOUNTERR" "$(column 'frame.fits[PSFPHOT]' PSFCOUNTERR PSFMODERR | \
    awk '$1 > 0 && $2 > $1 * 1.000001 { n++ } END { print n + 0 }')" 0 0

munipack phcal --photometry psf --photsys-ref Johnson -f V --col-ra RAJ2000 \
    --col-dec DEJ2000 --col-mag Vmag --col-magerr e_Vmag -c cat.fits \
    frame.fits > /dev/null 2>&1
match "cat.fits[1]" "RAJ2000 DEJ2000 Vmag" "frame_cal.fits[PHOTOMETRY]" \
    "RAJ2000 DEJ2000 MAG MAGERR" 3e-4 > m.txt
check "stars matched" "$(wc -l < m.txt)" 140 150
awk '$2 < 99 && $3 > 0' m.txt > v.txt
read -r med sd n <<< "$(awk '{ print $2 - $1 }' v.txt | stats)"
check "median MAG - V" "$med" -0.003 0.003
check "robust sd MAG - V" "$sd" 0 0.006
read -r med sd n <<< "$(awk '{ print ($2 - $1) / $3 }' v.txt | stats)"
check "robust sd of (MAG - V)/MAGERR" "$sd" 0.7 1.4
