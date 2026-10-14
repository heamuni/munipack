# astrometry from scratch on rotated frames at Dec +60 (cos(dec) = 0.5):
# rotation and the CD matrix (BUGS S1-S3), RA residuals with cos(dec) (S7),
# the scale in pix/deg (A1, A18); and the same solution with --mproc 2.

mkcat cat.fits 80 21 30 60 0.08 10 15 0.004
for ang in 30 -120; do
    f=rot$ang.fits
    munipack artificial -c cat.fits $CATCOLS -f V --rcen 30 --dcen 60 \
        --scale 3600 --angle $ang --width 600 --height 600 --exptime 30 \
        --psf GAUSS --hwhm 1.2 --seed 4 --no-watermark --mask $f > /dev/null 2>&1
    true_cd=$(munipack fits -K CD1_1,CD1_2,CD2_1,CD2_2 --value $f | tr '\n' ' ')
    munipack find -f 3 $f > /dev/null 2>&1
    munipack aphot $f > /dev/null 2>&1
    munipack astrometry --remove $f > /dev/null 2>&1
    cp $f mp$ang.fits
    munipack astrometry -c cat.fits $CATCOLS $f > /dev/null 2>&1
    read -r cr1 cr2 c11 c12 c21 c22 <<< "$(munipack fits -K \
        CRVAL1,CRVAL2,CD1_1,CD1_2,CD2_1,CD2_2 --value $f | tr '\n' ' ')"
    check "angle $ang: CRVAL1 offset [px]" "$(awk -v v="$cr1" 'BEGIN { print (v - 30) * 3600 * 0.5 }')" -0.05 0.05
    check "angle $ang: CRVAL2 offset [px]" "$(awk -v v="$cr2" 'BEGIN { print (v - 60) * 3600 }')" -0.05 0.05
    # the largest difference of the CD elements, in px at the frame edge (300 px)
    check "angle $ang: CD matrix error [px at 300 px]" "$(echo "$true_cd $c11 $c12 $c21 $c22" | awk '{
        m = 0; for (i = 1; i <= 4; i++) { d = $i - $(i + 4); d = d < 0 ? -d : d; if (d > m) m = d }
        print m * 3600 * 300 }')" 0 0.05

    # parallel astrometry gives the same solution
    munipack astrometry --mproc 2 -c cat.fits $CATCOLS mp$ang.fits > /dev/null 2>&1
    [ "$(munipack fits -K CRVAL1,CRVAL2,CD1_1,CD1_2,CD2_1,CD2_2 --value mp$ang.fits)" = \
      "$(munipack fits -K CRVAL1,CRVAL2,CD1_1,CD1_2,CD2_1,CD2_2 --value $f)" ] \
        || fail "angle $ang: --mproc 2 gives another solution"
    ok "angle $ang: --mproc 2 identical"
done
