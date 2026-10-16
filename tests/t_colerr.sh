# A colour term against a catalogue whose colours are noisy: the true
# B - V of 120 stars is uniform in 0 .. 1.5 (sd 0.43), the catalogue V and
# B are the truth plus their errors 0.01 and 0.2 mag (quoted in e_Vmag,
# e_Bmag), so the catalogue B - V scatters by 0.2 about the true one
# (the reliability ratio var(true)/var(observed) 0.82). The V frame is
# rendered from X = V + 0.3 (B - V) of the true magnitudes, so k1 = -0.3
# exactly. The fit on the noisy colours is biased towards zero by
# the errors (the regression dilution, -0.24, six errors of k1); with
# the error columns phcal corrects it (BUGS P57), --color-err-correct no
# does not (the previous version: -0.24).

mkcat truth.fits 120 61 200 -20 0.07 10 15 0.006
column "truth.fits[1]" RAJ2000 DEJ2000 Vmag Bmag \
    | awk '{ print $1, $2, $3 + 0.3 * ($4 - $3), $4 }' | catfits inst.fits
# the catalogue measured with its errors (Gaussian, Park-Miller and
# Box-Muller: reproducible)
column "truth.fits[1]" RAJ2000 DEJ2000 Vmag Bmag | awk '
    function rnd() { s = (s * 16807) % 2147483647; return s / 2147483647 }
    function gauss() { return sqrt(-2 * log(rnd())) * cos(2 * pi * rnd()) }
    BEGIN { s = 4242; pi = atan2(0, -1) }
    { printf "%s %s 0 0 %.4f 0.01 %.4f 0.2\n", $1, $2, $3 + 0.01 * gauss(), $4 + 0.2 * gauss() }' \
    | tabfits cat.fits "RAJ2000 D deg DEJ2000 D deg pmRA D mas/yr pmDE D mas/yr \
    Vmag E - e_Vmag E - Bmag E - e_Bmag E -"
check "catalogue stars" "$(nrows "cat.fits[1]")" 120 120

munipack artificial -c inst.fits $CATCOLS -f V --rcen 200 --dcen -20 \
    --scale 3600 --width 540 --height 540 --exptime 60 --psf GAUSS --hwhm 1.4 \
    --gain 2 --readnoise 4 --bias 100 --seed 5 --no-watermark --mask V.fits \
    > /dev/null 2>&1
munipack find -f 3.3 V.fits > /dev/null 2>&1
munipack aphot V.fits > /dev/null 2>&1

C="--photsys-ref Johnson -f V --col-ra RAJ2000 --col-dec DEJ2000 -c cat.fits \
    --col-mag Vmag --col-magerr e_Vmag"
pull() {   # pull DIR: (COLTK1 + 0.3) / COLTK1E
    awk -v k="$(key "$1/V.fits[PHOTOMETRY]" COLTK1)" \
        -v e="$(key "$1/V.fits[PHOTOMETRY]" COLTK1E)" 'BEGIN { print (k + 0.3) / e }'
}

# corrected for the errors of the colours (the default)
mkdir ce cn nx
munipack phcal $C --color-term Bmag,Vmag,e_Bmag,e_Vmag -t ce V.fits > /dev/null 2>&1
h="ce/V.fits[PHOTOMETRY]"
check "COLTNSTA" "$(key "$h" COLTNSTA)" 110 120
check "COLTK1 (-0.3)" "$(key "$h" COLTK1)" -0.34 -0.26
check "(COLTK1 + 0.3) / COLTK1E" "$(pull ce)" -2.5 2.5
check "COLTDILU" "$(key "$h" COLTDILU | grep -c T)" 1 1
check "COLTLAM (0.82)" "$(key "$h" COLTLAM)" 0.7 0.92

# not corrected: k1 diluted towards zero by the errors of the colours
munipack phcal $C --color-term Bmag,Vmag,e_Bmag,e_Vmag --color-err-correct no \
    -t cn V.fits > /dev/null 2>&1
h="cn/V.fits[PHOTOMETRY]"
check "not corrected: COLTDILU" "$(key "$h" COLTDILU | grep -c F)" 1 1
check "not corrected: COLTLAM as corrected" "$(key "$h" COLTLAM)" \
    "$(key "ce/V.fits[PHOTOMETRY]" COLTLAM)" "$(key "ce/V.fits[PHOTOMETRY]" COLTLAM)"
check "not corrected: (COLTK1 + 0.3) / COLTK1E" "$(pull cn)" 2.5 10

# without the error columns: no correction, no keys of it
munipack phcal $C --color-term Bmag,Vmag -t nx V.fits > /dev/null 2>&1
[ -z "$(key "nx/V.fits[PHOTOMETRY]" COLTDILU)" ] || fail "COLTDILU without the error columns"
ok "no COLTDILU without the error columns"
check "no error columns: (COLTK1 + 0.3) / COLTK1E" "$(pull nx)" 2.5 10

expect_fail "--color-err-correct maybe" munipack phcal $C \
    --color-term Bmag,Vmag,e_Bmag,e_Vmag --color-err-correct maybe V.fits
