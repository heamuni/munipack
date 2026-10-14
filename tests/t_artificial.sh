# artificial: reproducible frames by --seed (BUGS A2), the detector noise
# model (A13, A19, §2.1), no sky with --sky-mag 99 (A33), integer frames
# rounded (A34), the number of random stars by the frame area (A23),
# --mproc (A25).

# A field without stars (one star of 30 mag): only the detector.
mkcat empty.fits 1 1 180 0 0.001 30 30
art() { munipack artificial -c empty.fits $CATCOLS -f V --rcen 180 --dcen 0 \
          --width 200 --height 200 --no-watermark "$@" > /dev/null 2>&1; }

# Reproducibility: the same seed gives the same frame, another seed another one.
art --seed 7 --mask s1.fits
art --seed 7 --mask s2.fits
art --seed 8 --mask s3.fits
[ "$(imsum s1.fits)" = "$(imsum s2.fits)" ] || fail "--seed 7 twice: frames differ"
[ "$(imsum s1.fits)" != "$(imsum s3.fits)" ] || fail "--seed 7 and 8: frames equal"
ok "--seed reproducible"

# Detector only: dark 0.5 e-/s/px over 100 s, read noise 5 e-, gain 2 e-/ADU,
# bias 100 ADU: mean 100 + 50/2 = 125 ADU, sd sqrt(50 + 25)/2 = 4.330 ADU.
art --seed 1 --sky-mag 99 --exptime 100 --dark 0.5 --readnoise 5 --gain 2 \
    --bias 100 --mask det.fits
read -r mean sd n <<< "$(imstats det.fits)"
check "detector mean [ADU]" "$mean" 124.9 125.1
check "detector sd [ADU]" "$sd" 4.24 4.42
check "GAIN" "$(key det.fits GAIN)" 2 2
check "RDNOISE" "$(key det.fits RDNOISE)" 5 5

# Without noise: the bias plus the dark (and the light of the 30 mag star).
art --seed 1 --sky-mag 99 --exptime 100 --dark 0.5 --gain 2 --bias 100 \
    --disable-noise --mask flat.fits
read -r mean sd n <<< "$(imstats flat.fits)"
check "noiseless mean [ADU]" "$mean" 124.999 125.001
check "noiseless sd [ADU]" "$sd" 0 0.01

# --sky-mag 99 alone: no sky, no dark: zeros (A33 added ~1 e-/s/px).
art --seed 1 --sky-mag 99 --exptime 100 --mask nosky.fits
read -r mean sd n <<< "$(imstats nosky.fits)"
check "--sky-mag 99 mean" "$mean" 0 0.001

# 16-bit frame: rounded, not truncated (A34): the mean stays 125 (truncation
# would give 124.5).
art --seed 1 --sky-mag 99 --exptime 100 --dark 0.5 --readnoise 5 --gain 2 \
    --bias 100 -B 16 --mask det16.fits
read -r mean sd n <<< "$(imstats det16.fits)"
check "16-bit mean [ADU]" "$mean" 124.9 125.1
check "BITPIX" "$(key det16.fits BITPIX)" 16 16

# Random star field (no catalogue): a density per square degree, so a frame
# of 4x the area has ~4x the stars (A23: 239 on any frame).
munipack artificial --seed 3 --mask r1.fits > /dev/null 2>&1
munipack artificial --seed 3 --width 1332 --height 1332 --mask r2.fits > /dev/null 2>&1
check "random stars on 666^2" "$(nrows 'r1.fits[CATALOGUE]')" 1 20
check "random stars on 1332^2" "$(nrows 'r2.fits[CATALOGUE]')" 15 60

# Parallel processes give bit-identical frames (A25).
art --seed 4 --count 2 --mask 'p1_?.fits'
art --seed 4 --count 2 --mproc 2 --mask 'p2_?.fits'
for i in 1 2; do
    [ -f p2_$i.fits ] || fail "--mproc 2: no p2_$i.fits"
    [ "$(imsum p1_$i.fits)" = "$(imsum p2_$i.fits)" ] || fail "--mproc 2: frame $i differs"
done
ok "--mproc 2 bit-identical"
