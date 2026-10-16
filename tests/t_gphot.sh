# Growth-curve photometry: gphot on a simulated field with a known truth
# (the CATALOGUE extension of artificial: COUNT/GAIN is the total flux in
# ADU), then phcal --photometry grow.
# 200 stars V = 10-16.5 at Dec +30, Gaussian PSF (sigma 1.5 px), a
# ground-like detector (gain 2, read noise 5 e-).
# Covers: GCOUNT of every star by the enclosed energy of the frame's growth
# curve (BUGS P24: no valid GCOUNT with an extended PSF), GCOUNTERR with
# the error of the growth curve (bright stars' errors were too small),
# the corrected sky, phcal of GCOUNT, --legacy, too few stars.

mkcat cat.fits 200 21 180 30 0.1 10 16.5 0.004
munipack artificial -c cat.fits $CATCOLS -f V --rcen 180 --dcen 30 \
    --scale 3600 --width 800 --height 800 --exptime 60 --psf GAUSS --hwhm 1.5 \
    --gain 2 --readnoise 5 --bias 100 --dark 0.01 --sky-mag 21 --seed 5 \
    --no-watermark --mask frame.fits > /dev/null 2>&1
munipack find -f 3.5 frame.fits > /dev/null 2>&1
munipack aphot --apertures 1,1.5,2,2.6,3.4,4.5,6,8,10.5,14,18 --ring 20,30 \
    frame.fits > /dev/null 2>&1

cp frame.fits g.fits
munipack gphot g.fits > /dev/null 2>&1
check "GROWMETH is EE" "$(key 'g.fits[GROWPHOT]' GROWMETH | grep -c EE)" 1 1
check "stars of the growth curve" "$(key 'g.fits[GROWPHOT]' EENSTAR)" 30 200
# the true HWHM is 1.1774 sigma = 1.766 px (pixelised)
check "GROWPHOT HWHM [px]" "$(key 'g.fits[GROWPHOT]' HWHM)" 1.6 2.1
n=$(nrows 'g.fits[GROWPHOT]')
column 'g.fits[GROWPHOT]' X Y GCOUNT GCOUNTERR SKY SKYERR > g.txt
check "valid GCOUNT [%]" "$(awk -v n="$n" '$3 > 0 && $4 > 0 { k++ } END { print 100 * k / n }' g.txt)" 97 100

# join with the truth by the positions (1 px), GCOUNT / truth and the pull
gain=$(key frame.fits GAIN)
column 'frame.fits[CATALOGUE]' X Y COUNT | sed 's/^/C /' > c.txt
awk -v g="$gain" '
    $1 == "C" { n++; cx[n] = $2; cy[n] = $3; ct[n] = $4 / g; next }
    $3 > 0 && $4 > 0 { best = 0; bd = 1
        for (i = 1; i <= n; i++) { d = ($1 - cx[i])^2 + ($2 - cy[i])^2; if (d < bd) { bd = d; best = i } }
        if (best) print ct[best], $3, $4 }' c.txt g.txt > m.txt
check "stars joined with the truth" "$(wc -l < m.txt)" 180 200
read -r med sd nb <<< "$(awk '$2 / $3 > 50 { print $2 / $1 - 1 }' m.txt | stats)"
check "median GCOUNT/truth - 1 at S/N > 50" "$med" -0.003 0.003
read -r med sd n <<< "$(awk '{ print ($2 - $1) / $3 }' m.txt | stats)"
check "median pull (GCOUNT - truth)/GCOUNTERR" "$med" -0.3 0.3
check "robust sd of the pull" "$sd" 0.8 1.25
read -r med sd n <<< "$(awk '$2 / $3 > 300 { print ($2 - $1) / $3 }' m.txt | stats)"
check "bright stars (S/N > 300)" "$n" 15 200
check "robust sd of the pull of bright stars" "$sd" 0.7 1.3

# the sky of the star by the ring and the apertures: SKY + b of the
# growth curve, b within the error of the sky of the ring, the error of
# the sum smaller than that of the ring (P6, P38)
column 'g.fits[APERPHOT]' SKY SKYERR > a.txt
paste g.txt a.txt | awk '$6 > 0 { print ($5 - $7) / $8 }' > s.txt
read -r med sd n <<< "$(stats < s.txt)"
check "median (SKY - SKY of APERPHOT)/SKYERR of APERPHOT" "$med" -0.5 0.5
check "SKYERR > SKYERR of APERPHOT" "$(paste g.txt a.txt | \
    awk '$6 > 0 && $6 > $8 * 1.0001 { k++ } END { print k + 0 }')" 0 0

# calibration by GCOUNT
munipack phcal --photometry grow --photsys-ref Johnson -f V --col-ra RAJ2000 \
    --col-dec DEJ2000 --col-mag Vmag --col-magerr e_Vmag -c cat.fits g.fits \
    > /dev/null 2>&1
check "PHOTOTYP" "$(key 'g_cal.fits[PHOTOMETRY]' PHOTOTYP | grep -c GPHOT)" 1 1
match "cat.fits[1]" "RAJ2000 DEJ2000 Vmag" "g_cal.fits[PHOTOMETRY]" \
    "RAJ2000 DEJ2000 MAG MAGERR" 3e-4 | awk '$2 < 99 && $3 > 0' > v.txt
check "calibrated stars" "$(wc -l < v.txt)" 180 200
read -r med sd n <<< "$(awk '$3 < 0.01 { print $2 - $1 }' v.txt | stats)"
check "median MAG - V (MAGERR < 0.01)" "$med" -0.003 0.003
check "robust sd MAG - V (MAGERR < 0.01)" "$sd" 0 0.006
read -r med sd n <<< "$(awk '{ print ($2 - $1) / $3 }' v.txt | stats)"
check "robust sd of (MAG - V)/MAGERR" "$sd" 0.7 1.3

# the growth-curve fit of the versions before
cp frame.fits l.fits
munipack gphot --legacy l.fits > /dev/null 2>&1
check "--legacy: GROWMETH" "$(key 'l.fits[GROWPHOT]' GROWMETH | grep -c LEGACY)" 1 1
expect_fail "--threshold without --legacy" munipack gphot -th 0.2 frame.fits

# too few stars for a growth curve: an error, not a crash
mkcat few.fits 2 3 180 30 0.05 12 13 0.004
munipack artificial -c few.fits $CATCOLS -f V --rcen 180 --dcen 30 \
    --scale 3600 --width 400 --height 400 --exptime 60 --psf GAUSS --hwhm 1.5 \
    --gain 2 --readnoise 5 --bias 100 --seed 6 --no-watermark --mask few2.fits \
    > /dev/null 2>&1
munipack find -f 3.5 few2.fits > /dev/null 2>&1
munipack aphot few2.fits > /dev/null 2>&1
expect_fail "gphot with two stars" munipack gphot few2.fits
