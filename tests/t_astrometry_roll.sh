# astrometry by the pointing and the scale with the roll unknown (BUGS S29):
# a header of a wrong rotation, the pointing by --rcen, --dcen, --scale
# without a WCS, --roll-known, --reflex-unknown, the units of --wcs-offset,
# a wrong sky with the large search; the blind match by OpenMP threads
# gives the same solution for any number of them (S30); the catalogue
# of a cone larger than the frame by its pointing (S31).

# the solution against the truth: "dCRVAL1 dCRVAL2 dCD" in px (CD at 300 px)
truth() { munipack fits -K CRVAL1,CRVAL2,CD1_1,CD1_2,CD2_1,CD2_2 --value "$1" | tr '\n' ' '; }
err() {
    echo "$1 $(truth "$2")" | awk '{ pi = atan2(0, -1)
        d1 = ($7 - $1) * 3600 * cos($2 * pi / 180); d2 = ($8 - $2) * 3600
        m = 0; for (i = 3; i <= 6; i++) { d = $i - $(i + 6); d = d < 0 ? -d : d; if (d > m) m = d }
        printf "%.4f %.4f %.4f\n", d1, d2, m * 3600 * 300 }'
}
checkerr() {
    read -r e1 e2 ecd <<< "$(err "$2" "$3")"
    check "$1: CRVAL1 [px]" "$e1" -0.05 0.05
    check "$1: CRVAL2 [px]" "$e2" -0.05 0.05
    check "$1: CD [px at 300 px]" "$ecd" 0 0.05
}
meth() { key "$1" ASTRMETH | tr -d "' "; }
art() {   # art CAT OUT SIZE [options]: a frame of 1"/px at RA 30, Dec 60, rotated by 20 deg
    local c=$1 o=$2 s=$3; shift 3
    munipack artificial -c "$c" $CATCOLS -f V --rcen 30 --dcen 60 --scale 3600 --angle 20 \
        --width "$s" --height "$s" --exptime 30 --psf GAUSS --hwhm 1.2 --seed 4 \
        --no-watermark --mask "$o" "$@" > /dev/null 2>&1
    munipack find -f 3 "$o" > /dev/null 2>&1
    munipack aphot "$o" > /dev/null 2>&1
}
# perturb IN OUT DX DY ROT SCALE FLIP: the WCS of the header moved by DX, DY
# [px] (CRVAL by the WCS), rotated by ROT [deg], scaled, mirrored (FLIP -1)
perturb() {
    cp "$1" "$2"
    read -r r1 r2 c11 c12 c21 c22 <<< "$(munipack fits -K CRVAL1,CRVAL2,CD1_1,CD1_2,CD2_1,CD2_2 --value "$2" | tr '\n' ' ')"
    read -r n1 n2 n11 n12 n21 n22 <<< "$(echo "$r1 $r2 $c11 $c12 $c21 $c22 $3 $4 $5 $6 $7" | awk '{ pi = atan2(0, -1)
        a = $9 * pi / 180; s = $10; c = cos(a); n = sin(a); f = $11
        # the pointing moved by the old CD, the new CD = s CD R(a) diag(f, 1)
        da = ($3 * $7 + $4 * $8) / cos($2 * pi / 180); dd = $5 * $7 + $6 * $8
        printf "%.10f %.10f %.12e %.12e %.12e %.12e\n", $1 + da, $2 + dd,
            f * s * ($3 * c + $4 * n), s * (-$3 * n + $4 * c),
            f * s * ($5 * c + $6 * n), s * (-$5 * n + $6 * c) }')"
    for kv in CRVAL1=$n1 CRVAL2=$n2 CD1_1=$n11 CD1_2=$n12 CD2_1=$n21 CD2_2=$n22; do
        munipack fits --update --key "${kv%%=*}" --val "${kv#*=}" "$2" > /dev/null
    done
}

mkcat cat.fits 150 31 30 60 0.12 10 15.5 0.004
art cat.fits f.fits 600
t=$(truth f.fits)

# 1. the header rotated by 137 deg (beyond --wcs-rot), moved by 50 px, scaled
# by 1.5%: the guided match of any rotation (ASTRMETH = POINTING)
perturb f.fits r.fits 30 -40 137 1.015 1
cp r.fits k.fits
munipack astrometry -c cat.fits $CATCOLS r.fits > /dev/null 2>&1 || fail "the rotated header: not solved"
[ "$(meth r.fits)" = "POINTING" ] || fail "the rotated header: ASTRMETH = '$(meth r.fits)', expected POINTING"
ok "the header rotated by 137 deg: solved by the guided match of any rotation (POINTING)"
checkerr "any rotation" "$t" r.fits

# 2. --roll-known: the rotation of the header is trusted, the blind match follows
munipack astrometry --roll-known -c cat.fits $CATCOLS k.fits > /dev/null 2>&1 || fail "--roll-known: not solved"
[ "$(meth k.fits)" = "BLIND" ] || fail "--roll-known: ASTRMETH = '$(meth k.fits)', expected BLIND"
ok "--roll-known: the blind match after the failed rotation of the header"
checkerr "--roll-known (blind)" "$t" k.fits

# 3. no WCS in the header: the pointing by --rcen, --dcen (30" off) and --scale
# (1% off), any rotation
cp f.fits p.fits
munipack astrometry --remove p.fits > /dev/null 2>&1
read -r rc dc <<< "$(munipack fits -K CRVAL1,CRVAL2 --value f.fits | tr '\n' ' ' | awk '{ printf "%.8f %.8f", $1 + 0.0075, $2 - 0.0080 }')"
munipack astrometry --rcen "$rc" --dcen "$dc" --scale 3564 -c cat.fits $CATCOLS p.fits > /dev/null 2>&1 \
    || fail "the pointing by --rcen, --dcen, --scale: not solved"
[ "$(meth p.fits)" = "POINTING" ] || fail "the pointing: ASTRMETH = '$(meth p.fits)', expected POINTING"
ok "the pointing by --rcen, --dcen, --scale without a WCS: POINTING"
checkerr "the pointing" "$t" p.fits

# 4. the header mirrored and rotated: --reflex-unknown tries both reflections
perturb f.fits m.fits -20 25 75 0.99 -1
munipack astrometry --reflex-unknown --roll-unknown -c cat.fits $CATCOLS m.fits > /dev/null 2>&1 \
    || fail "--reflex-unknown: not solved"
[ "$(meth m.fits)" = "POINTING" ] || fail "--reflex-unknown: ASTRMETH = '$(meth m.fits)', expected POINTING"
ok "the header mirrored: --reflex-unknown solves by the guided match"
checkerr "--reflex-unknown" "$t" m.fits

# 5. the units of --wcs-offset: arcmin, %; unknown units are an error
perturb f.fits u.fits 200 150 -100 1.01 1
munipack astrometry --wcs-offset 5arcmin -c cat.fits $CATCOLS u.fits > /dev/null 2>&1 \
    || fail "--wcs-offset 5arcmin: not solved"
[ "$(meth u.fits)" = "POINTING" ] || fail "--wcs-offset 5arcmin: ASTRMETH = '$(meth u.fits)', expected POINTING"
ok "the header moved by 250 px, --wcs-offset 5arcmin: POINTING"
checkerr "--wcs-offset 5arcmin" "$t" u.fits
cp r.fits v.fits
expect_fail "--wcs-offset 3parsec" munipack astrometry --wcs-offset 3parsec -c cat.fits $CATCOLS v.fits
expect_fail "--roll-known --roll-unknown" munipack astrometry --roll-known --roll-unknown -c cat.fits $CATCOLS v.fits

# 6. a wrong sky (other stars at the same place) with the large search
# (any rotation, both reflections, the offset of half of the frame): no solution
mkcat wrong.fits 150 77 30 60 0.12 10 15.5 0.004
perturb f.fits w.fits 30 -40 137 1.015 1
expect_fail "wrong sky, the large search" munipack astrometry --roll-unknown --reflex-unknown \
    --wcs-offset 50% -c wrong.fits $CATCOLS w.fits
[ "$(meth w.fits)" = "" ] || fail "wrong sky: a WCS saved (ASTRMETH = '$(meth w.fits)')"

# 7. a cone of 0.3 deg (3.4 x the frame, 900 stars) by the pointing: the stars
# near the frame are matched; without the pointing, the brightest of the cone
mkcat cone.fits 900 5 30 60 0.3 10 15.5 0.003
art cone.fits c.fits 600
tc=$(truth c.fits)
munipack astrometry --remove c.fits > /dev/null 2>&1
cp c.fits c2.fits
read -r rc dc <<< "$(munipack fits -K CRVAL1,CRVAL2 --value f.fits | tr '\n' ' ' | awk '{ printf "%.8f %.8f", $1 - 0.006, $2 + 0.005 }')"
munipack astrometry --rcen "$rc" --dcen "$dc" --scale 3600 -c cone.fits $CATCOLS c.fits > /dev/null 2>&1 \
    || fail "the cone by the pointing: not solved"
check "the cone by the pointing: stars used" "$(munipack fits -lh c.fits | awk -F'= ' '/Objects used/ { print $2 + 0 }')" 60 250
checkerr "the cone by the pointing" "$tc" c.fits

# 8. the blind match by 1 and 4 OpenMP threads (the frame of the cone without
# the pointing: several chunk pairs): the same solution
cp c2.fits c4.fits
OMP_NUM_THREADS=1 munipack astrometry --blind --scale 3600 -c cone.fits $CATCOLS c2.fits > /dev/null 2>&1 || true
OMP_NUM_THREADS=4 munipack astrometry --blind --scale 3600 -c cone.fits $CATCOLS c4.fits > /dev/null 2>&1 || true
[ "$(truth c2.fits)" = "$(truth c4.fits)" ] || fail "1 and 4 threads: different solutions"
[ "$(key c2.fits ASTRMETH)" = "$(key c4.fits ASTRMETH)" ] || fail "1 and 4 threads: different outcomes"
ok "the blind match by 1 and 4 threads: the same outcome ($(meth c2.fits)) and solution"
