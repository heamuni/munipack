# phcal --syserr auto with few calibration stars (BUGS P54): the frame is
# rendered from the catalogue, the calibration catalogue scatters around
# it by 0.05 mag (6 stars) or 0.03 mag (30 stars). The header gets
# the upper error PHSYSERE and the number of stars PHSYSNST, fewer than 10
# stars a warning (stderr, COMMENT); a fixed --syserr writes neither.

mkcat cat.fits 30 17 180 0 0.035 11 14 0.004
munipack artificial -c cat.fits $CATCOLS -f V --rcen 180 --dcen 0 \
    --scale 3600 --width 300 --height 300 --exptime 60 --psf GAUSS --hwhm 1.5 \
    --gain 2 --readnoise 5 --bias 100 --seed 3 --no-watermark \
    --mask frame.fits > /dev/null 2>&1
munipack find -f 3.5 frame.fits > /dev/null 2>&1
munipack aphot frame.fits > /dev/null 2>&1

# the calibration catalogues: 0.05 mag scatter of the first 6 stars,
# 0.03 mag of all 30
column "cat.fits[1]" RAJ2000 DEJ2000 Vmag Bmag | awk 'NR <= 6 {
    split("0.05 -0.04 0.06 -0.07 0.02 -0.03", d, " "); print $1, $2, $3 + d[NR], $4 }' \
    | catfits few.fits
column "cat.fits[1]" RAJ2000 DEJ2000 Vmag Bmag | awk '{
    d = 0.03 * (NR % 2 ? 1 : -1) * (0.5 + (NR % 5) / 4); print $1, $2, $3 + d, $4 }' \
    | catfits many.fits

pcal() {   # CATALOGUE DIR [OPTIONS]
    local c=$1 d=$2
    shift 2
    mkdir -p $d
    munipack phcal --photsys-ref Johnson -f V --col-ra RAJ2000 --col-dec DEJ2000 \
        --col-mag Vmag --col-magerr e_Vmag -c $c -t $d "$@" frame.fits \
        > $d/out.txt 2> $d/err.txt
}
comments() { munipack fits --header "$1[PHOTOMETRY]" | grep -c "$2"; }

pcal few.fits few
check "6 stars: PHSYSNST" "$(key 'few/frame.fits[PHOTOMETRY]' PHSYSNST)" 6 6
check "6 stars: PHSYSERR [mag]" "$(key 'few/frame.fits[PHOTOMETRY]' PHSYSERR)" 0.02 0.10
check "6 stars: PHSYSERE [mag]" "$(key 'few/frame.fits[PHOTOMETRY]' PHSYSERE)" 0.01 0.10
check "6 stars: warning" "$(grep -c 'estimated from 6 stars only' few/err.txt)" 1 1
check "6 stars: COMMENT interval" "$(comments few/frame.fits '68% interval')" 1 1
check "6 stars: COMMENT warning" "$(comments few/frame.fits 'Warning: the systematic')" 1 1
# CTPHERR with the uncertainty of the estimated PHSYSERR (P55): at least
# t(5) = 1.11 times the error of the mean of 6 stars of the scatter PHSYSERR
# (the errors of the stars, a few mmag, are negligible; it was 1.0)
check "6 stars: CTPHERR / (PHSYSERR/sqrt 6)" "$(awk -v t="$(key 'few/frame.fits[PHOTOMETRY]' CTPH)" \
    -v dt="$(key 'few/frame.fits[PHOTOMETRY]' CTPHERR)" \
    -v s="$(key 'few/frame.fits[PHOTOMETRY]' PHSYSERR)" \
    'BEGIN { print 1.0857 * dt / t * sqrt(6) / s }')" 1.06 3

pcal few.fits fixed --syserr 0.05
check "fixed: PHSYSERR [mag]" "$(key 'fixed/frame.fits[PHOTOMETRY]' PHSYSERR)" 0.05 0.05
[ -z "$(key 'fixed/frame.fits[PHOTOMETRY]' PHSYSERE)" ] || fail "fixed --syserr wrote PHSYSERE"
[ -z "$(key 'fixed/frame.fits[PHOTOMETRY]' PHSYSNST)" ] || fail "fixed --syserr wrote PHSYSNST"
check "fixed: no warning" "$(grep -c 'Warning: the systematic' fixed/err.txt)" 0 0

pcal many.fits many
check "30 stars: PHSYSNST" "$(key 'many/frame.fits[PHOTOMETRY]' PHSYSNST)" 28 30
check "30 stars: PHSYSERR [mag]" "$(key 'many/frame.fits[PHOTOMETRY]' PHSYSERR)" 0.02 0.04
check "30 stars: PHSYSERE [mag]" "$(key 'many/frame.fits[PHOTOMETRY]' PHSYSERE)" 0.002 0.01
check "30 stars: no warning" "$(grep -c 'Warning: the systematic' many/err.txt)" 0 0
