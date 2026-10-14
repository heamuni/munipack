# Front end: help and version, exit status 1 on errors (not 139 or 13,
# BUGS A24), rejected invalid options (A14, A16, A17), phfotran dispatched
# (P45).

munipack --version > /dev/null || fail "--version"
munipack --help > /dev/null 2>&1 || fail "--help"
for a in bias dark flat phcorr kombine astrometry find aphot gphot psfphot \
         phcal artificial phfotran timeseries fits cross; do
    munipack $a --help > /dev/null 2>&1 || fail "$a --help"
done
ok "help of every action"

expect_fail "unknown action" munipack nonsense
expect_fail "missing file" munipack find no-such-file.fits
expect_fail "artificial --psf FOO (A14)" munipack artificial --psf FOO --mask x.fits
expect_fail "artificial --spread FOO (A14)" munipack artificial --spread FOO --mask x.fits
expect_fail "artificial --mproc 0 (A17)" munipack artificial --mproc 0 --mask x.fits
expect_fail "artificial --rcen without --dcen (A16)" munipack artificial --rcen 10 --mask x.fits
[ ! -e x.fits ] || fail "a rejected artificial wrote x.fits"

# phfotran is dispatched (P45): an error about its input, not an unknown action
out=$(munipack phfotran 2>&1) && fail "phfotran without input succeeded"
echo "$out" | grep -qi "unknown\|unrecognized" && fail "phfotran not dispatched: $out"
ok "phfotran dispatched"
