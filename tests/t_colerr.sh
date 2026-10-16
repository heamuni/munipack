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
# (COLTLAM by the weights of the fit, which depend on its slope a little)
lc=$(key "ce/V.fits[PHOTOMETRY]" COLTLAM)
check "not corrected: COLTLAM near the corrected one" "$(key "$h" COLTLAM)" \
    "$(awk -v l="$lc" 'BEGIN { print l - 0.01 }')" "$(awk -v l="$lc" 'BEGIN { print l + 0.01 }')"
check "not corrected: (COLTK1 + 0.3) / COLTK1E" "$(pull cn)" 2.5 10

# without the error columns: no correction, no keys of it
munipack phcal $C --color-term Bmag,Vmag -t nx V.fits > /dev/null 2>&1
[ -z "$(key "nx/V.fits[PHOTOMETRY]" COLTDILU)" ] || fail "COLTDILU without the error columns"
ok "no COLTDILU without the error columns"
check "no error columns: (COLTK1 + 0.3) / COLTK1E" "$(pull nx)" 2.5 10

# the band of the calibration in the colour (V of B - V; QUVIK: u of
# u - g against u) and noisy: V of the catalogue with errors of 0.15 mag,
# B 0.01, the error of V in both y and the colour, cov = -0.15**2, which
# biases k1 further (-0.375 uncorrected); the correction subtracts it.
# The covariance is by the error of the band in the colour: a magnitude
# error column of 0.01 for the same V (e_Vlow) gives the same k1 and
# a warning (by that error the covariance was left out: about -0.42)
column "truth.fits[1]" RAJ2000 DEJ2000 Vmag Bmag | awk '
    function rnd() { s = (s * 16807) % 2147483647; return s / 2147483647 }
    function gauss() { return sqrt(-2 * log(rnd())) * cos(2 * pi * rnd()) }
    BEGIN { s = 777; pi = atan2(0, -1) }
    { printf "%s %s 0 0 %.4f 0.15 %.4f 0.01 0.01\n", $1, $2, $3 + 0.15 * gauss(), $4 + 0.01 * gauss() }' \
    | tabfits catv.fits "RAJ2000 D deg DEJ2000 D deg pmRA D mas/yr pmDE D mas/yr \
    Vmag E - e_Vmag E - Bmag E - e_Bmag E - e_Vlow E -"
CV="--photsys-ref Johnson -f V --col-ra RAJ2000 --col-dec DEJ2000 -c catv.fits --col-mag Vmag"
mkdir sv svn svl
munipack phcal $CV --col-magerr e_Vmag --color-term Bmag,Vmag,e_Bmag,e_Vmag -t sv V.fits > /dev/null 2>&1
check "noisy V in y and B - V: (COLTK1 + 0.3) / COLTK1E" "$(pull sv)" -2.5 2.5
munipack phcal $CV --col-magerr e_Vmag --color-term Bmag,Vmag,e_Bmag,e_Vmag \
    --color-err-correct no -t svn V.fits > /dev/null 2>&1
check "noisy V, not corrected: COLTK1 (-0.375)" "$(key "svn/V.fits[PHOTOMETRY]" COLTK1)" -0.45 -0.34
munipack phcal $CV --col-magerr e_Vlow --color-term Bmag,Vmag,e_Bmag,e_Vmag -t svl V.fits \
    > /dev/null 2> svl/err.txt
check "noisy V, --col-magerr of 0.01: (COLTK1 + 0.3) / COLTK1E" "$(pull svl)" -2.5 2.5
check "noisy V, --col-magerr of 0.01: the warning" "$(grep -c 'differ' svl/err.txt)" 1 9

# the zero point: CTPH of the colour 0.75 against the one of the exact
# catalogue (the truth, no error columns), orders 1 and 2; the order 2 of
# these colours is capped (the square of the colour hardly wider than its
# errors), which is reported
mkdir zt zc zt2 zc2
CT="--photsys-ref Johnson -f V --col-ra RAJ2000 --col-dec DEJ2000 -c truth.fits --col-mag Vmag"
dz() {   # dz A B: (CTPH(A) - CTPH(B)) / the errors
    awk -v a="$(key "$1/V.fits[PHOTOMETRY]" CTPH)" -v ea="$(key "$1/V.fits[PHOTOMETRY]" CTPHERR)" \
        -v b="$(key "$2/V.fits[PHOTOMETRY]" CTPH)" -v eb="$(key "$2/V.fits[PHOTOMETRY]" CTPHERR)" \
        'BEGIN { print (a - b) / sqrt(ea * ea + eb * eb) }'
}
munipack phcal $CT --color-term Bmag,Vmag --color-ref 0.75 -t zt V.fits > /dev/null 2>&1
munipack phcal $C --color-term Bmag,Vmag,e_Bmag,e_Vmag --color-ref 0.75 -t zc V.fits \
    > /dev/null 2>&1
check "CTPH of the colour 0.75 against the truth (in errors)" "$(dz zc zt)" -3 3
munipack phcal $CT --color-term Bmag,Vmag --color-ref 0.75 --color-order 2 -t zt2 V.fits \
    > /dev/null 2>&1
munipack phcal $C --color-term Bmag,Vmag,e_Bmag,e_Vmag --color-ref 0.75 --color-order 2 \
    -t zc2 V.fits > /dev/null 2> zc2/err.txt
h="zc2/V.fits[PHOTOMETRY]"
check "order 2: the cap reported" "$(grep -c capped zc2/err.txt)" 1 9
check "order 2: COLTDILU" "$(key "$h" COLTDILU | grep -c T)" 1 1
check "order 2: COLTCAP" "$(key "$h" COLTCAP)" 0.01 0.99
check "order 2: (COLTK1 + 0.3) / COLTK1E" "$(pull zc2)" -3 3
check "order 2: COLTK2 (0) / COLTK2E" "$(awk -v k="$(key "$h" COLTK2)" \
    -v e="$(key "$h" COLTK2E)" 'BEGIN { print k / e }')" -3 3
check "order 2: CTPH against the truth (in errors)" "$(dz zc2 zt2)" -3 3

# two colours sharing the noisy band of the calibration: V (0.1 mag) in
# y, B - V and V - R (cov -0.1**2), R - V independent of B - V; the frame
# of X = V + 0.3 (B - V) - 0.2 (V - R), so k = -0.3, +0.2
column "truth.fits[1]" RAJ2000 DEJ2000 Vmag Bmag | awk '
    function rnd() { s = (s * 16807) % 2147483647; return s / 2147483647 }
    BEGIN { s = 99 }
    { print $1, $2, $3, $4, $3 - 0.3 - 0.8 * rnd() }' > vbr.txt
awk '{ print $1, $2, $3 + 0.3 * ($4 - $3) - 0.2 * ($3 - $5), $4 }' vbr.txt | catfits inst3.fits
munipack artificial -c inst3.fits $CATCOLS -f V --rcen 200 --dcen -20 \
    --scale 3600 --width 540 --height 540 --exptime 60 --psf GAUSS --hwhm 1.4 \
    --gain 2 --readnoise 4 --bias 100 --seed 5 --no-watermark --mask V3.fits \
    > /dev/null 2>&1
munipack find -f 3.3 V3.fits > /dev/null 2>&1
munipack aphot V3.fits > /dev/null 2>&1
awk '
    function rnd() { s = (s * 16807) % 2147483647; return s / 2147483647 }
    function gauss() { return sqrt(-2 * log(rnd())) * cos(2 * pi * rnd()) }
    BEGIN { s = 555; pi = atan2(0, -1) }
    { printf "%s %s 0 0 %.4f 0.1 %.4f 0.01 %.4f 0.01\n", $1, $2, $3 + 0.1 * gauss(), \
        $4 + 0.01 * gauss(), $5 + 0.01 * gauss() }' vbr.txt \
    | tabfits cat3.fits "RAJ2000 D deg DEJ2000 D deg pmRA D mas/yr pmDE D mas/yr \
    Vmag E - e_Vmag E - Bmag E - e_Bmag E - Rmag E - e_Rmag E -"
C3="--photsys-ref Johnson -f V --col-ra RAJ2000 --col-dec DEJ2000 -c cat3.fits \
    --col-mag Vmag --col-magerr e_Vmag --color-term Bmag,Vmag,e_Bmag,e_Vmag \
    --color-term Vmag,Rmag,e_Vmag,e_Rmag"
mkdir t3 t3n
munipack phcal $C3 -t t3 V3.fits > /dev/null 2>&1
munipack phcal $C3 --color-err-correct no -t t3n V3.fits > /dev/null 2>&1
p3() {   # p3 DIR KEY TRUTH: the pull of a coefficient
    awk -v k="$(key "$1/V3.fits[PHOTOMETRY]" $2)" -v e="$(key "$1/V3.fits[PHOTOMETRY]" ${2}E)" \
        -v t="$3" 'BEGIN { print (k - t) / e }'
}
check "two colours: (COLTK1 + 0.3) / COLTK1E" "$(p3 t3 COLTK1 -0.3)" -3 3
check "two colours: (COLT2K1 - 0.2) / COLT2K1E" "$(p3 t3 COLT2K1 0.2)" -3 3
check "two colours: COLTLAM2" "$(key "t3/V3.fits[PHOTOMETRY]" COLTLAM2)" 0.5 0.99
check "two colours, not corrected: COLTDILU" \
    "$(key "t3n/V3.fits[PHOTOMETRY]" COLTDILU | grep -c F)" 1 1
check "two colours, not corrected: |pull of COLT2K1|" \
    "$(p3 t3n COLT2K1 0.2 | awk '{ print ($1 < 0 ? -$1 : $1) }')" 3 100

expect_fail "--color-err-correct maybe" munipack phcal $C \
    --color-term Bmag,Vmag,e_Bmag,e_Vmag --color-err-correct maybe V.fits
