# Colour terms: an instrument whose B, V differ from the catalogue ones,
# B_i = B - 0.10 (B-V), V_i = V + 0.04 (B-V). phcal --color-term recovers
# the colour term of one filter, phfotran --model mag both, and
# phcal --tratab with the magnitude table leaves no colour slope.

mkcat cat.fits 80 43 200 -20 0.07 10 15 0.008
column "cat.fits[1]" RAJ2000 DEJ2000 Vmag Bmag \
    | awk '{ c = $4 - $3; print $1, $2, $3 + 0.04 * c, $4 - 0.10 * c }' \
    | catfits inst.fits
for flt in B V; do
    munipack artificial -c inst.fits --col-ra RAJ2000 --col-dec DEJ2000 \
        --col-pm-ra pmRA --col-pm-dec pmDE --col-mag ${flt}mag -f $flt \
        --rcen 200 --dcen -20 --scale 3600 --width 540 --height 540 \
        --exptime 60 --psf GAUSS --hwhm 1.4 --gain 2 --readnoise 4 --bias 100 \
        --seed 7 --no-watermark --mask $flt.fits > /dev/null 2>&1
    munipack find -f 3.3 $flt.fits > /dev/null 2>&1
    munipack aphot $flt.fits > /dev/null 2>&1
done

C="--photsys-ref Johnson --col-ra RAJ2000 --col-dec DEJ2000 -c cat.fits"

# one filter: B against the catalogue B with the colour B - V
mkdir ct
munipack phcal $C -f B --col-mag Bmag --col-magerr e_Bmag \
    --color-term Bmag,Vmag,e_Bmag,e_Vmag -t ct B.fits > /dev/null 2>&1
check "B: COLTK1 (0.10)" "$(key "ct/B.fits[PHOTOMETRY]" COLTK1)" 0.095 0.105
check "B: COLORSRC column" "$(colnum "ct/B.fits[PHOTOMETRY]" COLORSRC)" 1 99
match "cat.fits[1]" "RAJ2000 DEJ2000 Bmag Vmag" \
    "ct/B.fits[PHOTOMETRY]" "RAJ2000 DEJ2000 MAG" 3e-4 \
    | awk '$3 < 99 { print $3 - $1, $1 - $2 }' > res.txt
check "B: stars" "$(wc -l < res.txt)" 76 80
check "B: colour slope [mag/mag]" "$(awk '{ x = $2; y = $1; n++; sx += x; sy += y
    sxx += x * x; sxy += x * y } END { print (n * sxy - sx * sy) / (n * sxx - sx * sx) }' res.txt)" \
    -0.003 0.003

# without the colour term the slope is that of the instrument
mkdir nc
munipack phcal $C -f B --col-mag Bmag --col-magerr e_Bmag -t nc B.fits > /dev/null 2>&1
[ -z "$(key "nc/B.fits[PHOTOMETRY]" COLTK1)" ] || fail "COLTK1 without --color-term"
ok "no colour keywords without --color-term"

expect_fail "--color-term with -r" munipack phcal --photsys-ref Johnson \
    --color-term Bmag,Vmag -r ct/B.fits B.fits

# phfotran --model mag on both filters, phcal --tratab
munipack phfotran $C -f B,V --col-mag Bmag,Vmag --col-magerr e_Bmag,e_Vmag \
    --model mag -o tra.fits B.fits V.fits > /dev/null 2>&1
check "TRAMODEL" "$(key tra.fits TRAMODEL | grep -c MAG)" 1 1
k=$(munipack fits -lt "tra.fits[FOTRAN_MAG]" | awk '{ print $7 }' | tr '\n' ' ')
check "phfotran K1 of B (0.10)" "$(echo $k | cut -d' ' -f1)" 0.095 0.105
check "phfotran K1 of V (-0.04)" "$(echo $k | cut -d' ' -f2)" -0.045 -0.035
mkdir tr
munipack phcal $C -f B,V --col-mag Bmag,Vmag --col-magerr e_Bmag,e_Vmag \
    --tratab tra.fits -t tr B.fits V.fits > /dev/null 2>&1
for flt in B V; do
    match "cat.fits[1]" "RAJ2000 DEJ2000 ${flt}mag Bmag Vmag" \
        "tr/$flt.fits[PHOTOMETRY]" "RAJ2000 DEJ2000 MAG" 3e-4 \
        | awk '$4 < 99 { print $4 - $1, $2 - $3 }' > res.txt
    read -r med sd n <<< "$(cut -d' ' -f1 res.txt | stats)"
    check "$flt --tratab: robust sd" "$sd" 0 0.004
    check "$flt --tratab: colour slope" "$(awk '{ x = $2; y = $1; n++; sx += x; sy += y
        sxx += x * x; sxy += x * y } END { print (n * sxy - sx * sy) / (n * sxx - sx * sx) }' res.txt)" \
        -0.003 0.003
done
