# astrometry by the WCS in the header (BUGS S21), the scale prior (S22),
# the significance of a match against a wrong sky (S23, S16), few-star
# frames and the runs of several frames (S24), the proper motions with
# cos(dec) and without the epoch of the catalogue (S25), the catalogue
# in any order (S26).

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

mkcat cat.fits 150 31 30 60 0.12 10 15.5 0.004
art cat.fits f.fits 600
t=$(truth f.fits)

# 1. the header moved by 75 px, rotated by 2.5 deg, scaled by 2.5%: the WCS-guided match
cp f.fits p.fits
read -r p1 p2 c11 c12 c21 c22 <<< "$(munipack fits -K CRPIX1,CRPIX2,CD1_1,CD1_2,CD2_1,CD2_2 --value p.fits | tr '\n' ' ')"
read -r n1 n2 n11 n12 n21 n22 <<< "$(echo "$p1 $p2 $c11 $c12 $c21 $c22" | awk '{ pi = atan2(0, -1)
    a = 2.5 * pi / 180; s = 1.025; c = cos(a); n = sin(a)
    printf "%.3f %.3f %.12e %.12e %.12e %.12e\n", $1 + 45, $2 - 60,
        s * ($3 * c + $4 * n), s * (-$3 * n + $4 * c), s * ($5 * c + $6 * n), s * (-$5 * n + $6 * c) }')"
for kv in CRPIX1=$n1 CRPIX2=$n2 CD1_1=$n11 CD1_2=$n12 CD2_1=$n21 CD2_2=$n22; do
    munipack fits --update --key "${kv%%=*}" --val "${kv#*=}" p.fits > /dev/null
done
cp p.fits b.fits
munipack astrometry -c cat.fits $CATCOLS p.fits > /dev/null 2>&1 || fail "the perturbed header: not solved"
[ "$(meth p.fits)" = "WCS" ] || fail "ASTRMETH = '$(meth p.fits)', expected WCS"
ok "the perturbed header: solved by the WCS-guided match (ASTRMETH = WCS)"
checkerr "WCS-guided" "$t" p.fits

# 2. the same frame by the blind match (--blind ignores the header)
munipack astrometry --blind -c cat.fits $CATCOLS b.fits > /dev/null 2>&1 || fail "--blind: not solved"
[ "$(meth b.fits)" = "BLIND" ] || fail "--blind: ASTRMETH = '$(meth b.fits)'"
ok "--blind: ASTRMETH = BLIND"
checkerr "blind" "$t" b.fits

# 3. a wrong sky (a catalogue of other stars at the same place): no solution,
# by the header (and the blind match after it) and by the blind match of any scale
mkcat wrong.fits 150 77 30 60 0.12 10 15.5 0.004
cp f.fits w1.fits; cp f.fits w2.fits
expect_fail "wrong sky, by the header" munipack astrometry -c wrong.fits $CATCOLS w1.fits
expect_fail "wrong sky, blind of any scale" munipack astrometry --blind --scale-free -c wrong.fits $CATCOLS w2.fits
[ -z "$(meth w1.fits)" ] && [ -z "$(meth w2.fits)" ] || fail "wrong sky: a WCS saved"

# 4. a small frame (200 x 200 px, about 15 stars) with the catalogue of its
# footprint and a margin of 20%, the WCS removed: blind; and with --scale
mkcat small.fits 30 9 30 60 0.04 10 15.5 0.004
art small.fits s.fits 200
ts=$(truth s.fits)
munipack astrometry --remove s.fits > /dev/null 2>&1
cp s.fits s2.fits
munipack astrometry -c small.fits $CATCOLS s.fits > /dev/null 2>&1 || fail "the small frame: not solved"
check "the small frame: stars used" "$(munipack fits -lh s.fits | awk -F'= ' '/Objects used/ { print $2 + 0 }')" 5 40
checkerr "small frame" "$ts" s.fits
munipack astrometry --scale 3600 -c small.fits $CATCOLS s2.fits > /dev/null 2>&1 || fail "--scale: not solved"
checkerr "small frame, --scale" "$ts" s2.fits

# 5. five stars at least: --minmatch 3 is an error; a frame of four stars fails
# without stopping the run: the frames after it are calibrated
column "cat.fits[1]" RAJ2000 DEJ2000 Vmag Bmag | sort -g -k3 | head -4 | catfits four.fits
art four.fits r4.fits 600
cp f.fits r1.fits; cp f.fits r2.fits
expect_fail "--minmatch 3" munipack astrometry --minmatch 3 -c cat.fits $CATCOLS r1.fits
log=$(munipack astrometry -c cat.fits $CATCOLS r1.fits r4.fits r2.fits 2>&1) && fail "three frames: exit 0"
grep -q "too few stars" <<< "$log" || fail "four stars: no error 'too few stars'"
for f in r1 r2; do
    [ "$(meth $f.fits)" = "WCS" ] || fail "three frames: $f not calibrated"
done
ok "three frames, the middle one of four stars: it fails (too few stars), the others calibrated"

# 6. proper motions (+-400 mas/yr) at Dec 60 (cos = 0.5): a frame of 2030
# against the catalogue of 2000 (JEPOCH = 2000, pmRA = mu_alpha cos(dec))
column "cat.fits[1]" RAJ2000 DEJ2000 Vmag | awk '
    function rnd() { s = (s * 16807) % 2147483647; return s / 2147483647 }
    BEGIN { s = 12345; pi = atan2(0, -1) }
    { pa = 800 * rnd() - 400; pd = 800 * rnd() - 400
      print $1, $2, pa, pd, $3 > "pm.txt"
      print $1 + 30 * pa / 3.6e6 / cos($2 * pi / 180), $2 + 30 * pd / 3.6e6, $3, $3 + 0.5 }' \
    | catfits cat2030.fits
tabfits pm0.fits "RAJ2000 D deg DEJ2000 D deg pmRA D mas/yr pmDE D mas/yr Vmag E -" < pm.txt
# the epoch of the catalogue in the table header
awk '/^EXTNAME/ { print "JEPOCH  =               2000.0" } { print }' pm0.lst > pm.lst
munipack fits --restore pm.lst > /dev/null
[ "$(key "pm.fits[1]" JEPOCH)" != "" ] || fail "no JEPOCH in pm.fits"
art cat2030.fits m.fits 600 --date 2030-01-01
tm=$(truth m.fits)
munipack astrometry --remove m.fits > /dev/null 2>&1
munipack astrometry -c pm.fits --col-ra RAJ2000 --col-dec DEJ2000 --col-pm-ra pmRA \
    --col-pm-dec pmDE --col-mag Vmag m.fits > /dev/null 2>&1 || fail "proper motions: not solved"
checkerr "proper motions to 2030" "$tm" m.fits

# without the epoch of the catalogue: a warning (stderr and the header)
awk '{ print $1, $2, $3 / 20, $4 / 20, $5 }' pm.txt \
    | tabfits pmne.fits "RAJ2000 D deg DEJ2000 D deg pmRA D mas/yr pmDE D mas/yr Vmag E -"
cp f.fits n.fits
log=$(munipack astrometry -c pmne.fits --col-ra RAJ2000 --col-dec DEJ2000 --col-pm-ra pmRA \
    --col-pm-dec pmDE --col-mag Vmag n.fits 2>&1) || fail "no epoch: not solved"
grep -q "proper motions of .* ignored" <<< "$log" || fail "no epoch: no warning"
hdr=$(munipack fits -lh n.fits)
grep -q "proper motions ignored" <<< "$hdr" || fail "no epoch: no COMMENT"
ok "no epoch of the catalogue: the warning on stderr and in the header"

# 7. the catalogue in any order: reversed, the same solution
column "cat.fits[1]" RAJ2000 DEJ2000 Vmag Bmag | sort -g -r -k3 | catfits rev.fits
cp f.fits v.fits
munipack astrometry -c rev.fits $CATCOLS v.fits > /dev/null 2>&1 || fail "reversed catalogue: not solved"
checkerr "reversed catalogue" "$t" v.fits
