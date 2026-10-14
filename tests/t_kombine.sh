# kombine: four frames shifted in RA by 5-16 px stacked onto the grid of a
# reference frame. The stack is the sum (BUGS S20), the stars keep their
# flux (S14) and their positions by the WCS (S9-S13), the exposure time
# is the total.

# isolated stars (30 px apart): the test is about kombine, not crowding
mkcat cat.fits 40 31 100 10 0.06 10 14 0.008
art() {   # RA SEED FILE
    munipack artificial -c cat.fits $CATCOLS -f V --rcen $1 --dcen 10 \
        --scale 3600 --width 500 --height 500 --exptime 30 --psf GAUSS \
        --hwhm 1.3 --gain 1 --readnoise 3 --seed $2 --no-watermark \
        --mask $3 > /dev/null 2>&1
}
art 100 9 ref.fits
art 99.9955 1 k1.fits
art 99.9985 2 k2.fits
art 100.0015 3 k3.fits
art 100.0045 4 k4.fits

# the grid of ref.fits (artificial: CRPIX = NAXIS/2)
munipack kombine --rcen 100 --dcen 10 --xcen 250 --ycen 250 --width 500 \
    --height 500 --scale 3600 -i bilinear -o sum.fits \
    k1.fits k2.fits k3.fits k4.fits > /dev/null 2>&1
check "EXPTIME of the stack" "$(key sum.fits EXPTIME)" 120 120

for f in ref sum; do
    munipack find -f 3 $f.fits > /dev/null 2>&1
    munipack aphot $f.fits > /dev/null 2>&1
done
{ column "ref.fits[APERPHOT]" X Y APCOUNT8 | sed 's/^/R /'
  column "sum.fits[APERPHOT]" X Y APCOUNT8 | sed 's/^/S /'; } | awk '
    $1 == "R" { n++; x[n] = $2; y[n] = $3; f[n] = $4; next }
    { for (i = 1; i <= n; i++)
          if (($2 - x[i])^2 + ($3 - y[i])^2 < 1 && f[i] > 0 && $4 > 0)
              print $4 / f[i], $2 - x[i], $3 - y[i] }' > pairs.txt
check "stars in both" "$(wc -l < pairs.txt)" 36 40
check "median flux ratio stack/single" "$(cut -d' ' -f1 pairs.txt | stats | cut -d' ' -f1)" 3.96 4.04
check "median x offset [px]" "$(cut -d' ' -f2 pairs.txt | stats | cut -d' ' -f1)" -0.03 0.03
check "median y offset [px]" "$(cut -d' ' -f3 pairs.txt | stats | cut -d' ' -f1)" -0.03 0.03
