# aphot beside brighter neighbours (BUGS P29): apertures that hold a
# neighbour's light are invalidated, apertures clear of it are kept.
# Designed pairs on a 400^2 frame (offsets from the centre in px, V), with
# isolated stars for the mean profile of the frame. Gaussian PSF
# (sigma 1.5 px): a star 12-20 px from a brighter one keeps its aperture
# of 5.1 px with the flux of an isolated star (all its apertures were
# invalidated by the neighbour's growth curve, which held the star itself);
# a 14 mag star 7 px from a 9 mag one is invalidated (it was 4 mag too
# bright before P29). Moffat PSF (beta 1.5, heavy wings): the same.

cat > stars.txt << EOF
-120 -120 11.0
-120 -108 11.5
0 -120 10.0
0 -104 12.5
120 -120 11.0
120 -100 11.0
-120 60 9.0
-120 66 15.0
60 60 9.0
60 67 14.0
-150 160 11.0
-40 160 12.0
70 160 11.5
150 0 11.0
-150 0 12.0
EOF
awk '{ printf "%.10f %.10f %s %s\n", 150 - $1 / 3600, $2 / 3600, $3, $3 + 0.5 }' \
    stars.txt | catfits cat.fits

# ap FILE DX DY COL: the column of the star at the offset (nan when not found)
ap() {
    column "$1[APERPHOT]" X Y "$4" | awk -v x=$((200 + $2)) -v y=$((200 + $3)) '
        ($1 - x)^2 + ($2 - y)^2 < 2 { v = $3 } END { print (v == "" ? "nan" : v) }'
}

for psf in GAUSS MOFFAT; do
    if [ $psf = GAUSS ]; then opt="--hwhm 1.5"; else opt="--hwhm 2.7 --beta 1.5"; fi
    munipack artificial -c cat.fits $CATCOLS -f V --rcen 150 --dcen 0 \
        --scale 3600 --width 400 --height 400 --exptime 60 --psf $psf $opt \
        --gain 2 --readnoise 5 --bias 100 --seed 3 --no-watermark \
        --mask $psf.fits > /dev/null 2>&1
    munipack find -f 3.5 $psf.fits > /dev/null 2>&1
    munipack aphot $psf.fits > /dev/null 2>&1

    # the faint star beside a bright one: invalid (aperture 7 = 5.1 px),
    # or not detected at all (find rejects detections in the wings, P25)
    v=$(ap $psf.fits 60 67 APCOUNT7)
    [ "$v" = nan ] && v=-1
    check "$psf: 14 mag at 7 px from 9 mag, APCOUNT7 (invalid or not found)" "$v" -1 -1
    # isolated stars and the bright member of the pair: valid
    for s in "150 0" "-150 0" "-150 160" "60 60"; do
        check "$psf: star at $s, APCOUNT7 > 0" "$(ap $psf.fits $s APCOUNT7)" 1 1e12
    done
done

# Gaussian: the fainter member of each pair has the flux of an isolated star
# of its magnitude (the 11 mag star at 150,0) within 0.5%.
ref=$(ap GAUSS.fits 150 0 APCOUNT7)
for s in "-120 -108 11.5" "0 -104 12.5" "120 -100 11.0"; do
    set -- $s
    check "GAUSS: V=$3 at $1,$2 (12-20 px from a brighter star), APCOUNT7/expected" \
        "$(awk -v c="$(ap GAUSS.fits $1 $2 APCOUNT7)" -v r="$ref" -v v=$3 \
            'BEGIN { print c / (r * 10^(-0.4 * (v - 11))) }')" 0.995 1.005
done
