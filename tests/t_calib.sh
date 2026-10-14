# Calibration frames on a QUVIK-like detector (gain 1, read noise 0.4 e-,
# pedestal 100 ADU): master bias; master dark of 0.13 e-/px, where the
# old robust mean was ~16% low (BUGS P31); phcorr subtracts bias and dark
# and writes the noise model in electrons (P39, P41).

mkcat empty.fits 1 1 180 0 0.001 30 30
art() {
    munipack artificial -c empty.fits $CATCOLS -f V --rcen 180 --dcen 0 \
        --width 150 --height 150 --no-watermark --gain 1 --readnoise 0.4 \
        --bias 100 "$@" > /dev/null 2>&1
}
for i in $(seq 1 10); do
    art --sky-mag 99 --exptime 0.001 --seed $i --mask b$i.fits
    art --sky-mag 99 --exptime 100 --dark 0.0013 --seed $((20 + i)) --mask d$i.fits
done

munipack bias -o bias.fits b*.fits > /dev/null 2>&1
check "master bias [ADU]" "$(imstats bias.fits | cut -d' ' -f1)" 99.995 100.005

munipack dark -bias bias.fits -o dark.fits d*.fits > /dev/null 2>&1
check "DARKMEAN" "$(key dark.fits DARKMEAN | grep -c CLIP)" 1 1
check "master dark [e-/px] (0.13)" "$(imstats dark.fits | cut -d' ' -f1)" 0.127 0.133
check "EXPTIME of the dark" "$(key dark.fits EXPTIME)" 100 100

# a frame with sky (3638.5 e-/px above the pedestal, noiseless) and dark
art --sky-mag 21 --exptime 100 --seed 99 --disable-noise --mask skyonly.fits
sky=$(imstats skyonly.fits | awk '{ print $1 - 100 }')
art --sky-mag 21 --exptime 100 --dark 0.0013 --seed 99 --mask sci.fits
mkdir out
munipack phcorr -bias bias.fits -dark dark.fits -t out sci.fits > /dev/null 2>&1
check "phcorr: level - sky [e-/px]" "$(imstats out/sci.fits | awk -v s="$sky" '{ print $1 - s }')" -1 1
check "phcorr: GAIN" "$(key out/sci.fits GAIN)" 1 1
[ -z "$(key out/sci.fits BIASLEV)" ] || fail "phcorr: BIASLEV kept on a frame in electrons"
ok "phcorr: no BIASLEV"
