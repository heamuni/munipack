# phcal with two filters: B and V frames of a field with B - V = 0 .. 1.5.
# The order of -f (and of the columns and files) does not matter (BUGS
# P14), each band matches its catalogue magnitudes without a colour term
# (P43).

mkcat cat.fits 80 41 200 -20 0.07 10 15 0.008
for flt in B V; do
    munipack artificial -c cat.fits --col-ra RAJ2000 --col-dec DEJ2000 \
        --col-pm-ra pmRA --col-pm-dec pmDE --col-mag ${flt}mag -f $flt \
        --rcen 200 --dcen -20 --scale 3600 --width 540 --height 540 \
        --exptime 60 --psf GAUSS --hwhm 1.4 --gain 2 --readnoise 4 --bias 100 \
        --seed 7 --no-watermark --mask $flt.fits > /dev/null 2>&1
    munipack find -f 3.3 $flt.fits > /dev/null 2>&1
    munipack aphot $flt.fits > /dev/null 2>&1
done

mkdir bv vb
munipack phcal --photsys-ref Johnson -f B,V --col-ra RAJ2000 --col-dec DEJ2000 \
    --col-mag Bmag,Vmag --col-magerr e_Bmag,e_Vmag -c cat.fits -t bv \
    B.fits V.fits > /dev/null 2>&1
munipack phcal --photsys-ref Johnson -f V,B --col-ra RAJ2000 --col-dec DEJ2000 \
    --col-mag Vmag,Bmag --col-magerr e_Vmag,e_Bmag -c cat.fits -t vb \
    V.fits B.fits > /dev/null 2>&1

for flt in B V; do
    check "$flt: FILTER" "$(key "bv/$flt.fits[PHOTOMETRY]" FILTER | grep -c "'$flt")" 1 1
    [ "$(munipack fits -lt "bv/$flt.fits[PHOTOMETRY]" | md5sum)" = \
      "$(munipack fits -lt "vb/$flt.fits[PHOTOMETRY]" | md5sum)" ] \
        || fail "$flt: -f B,V and -f V,B differ"
    ok "$flt: -f B,V and -f V,B identical"

    # residual MAG - catalogue and its slope against B - V
    match "cat.fits[1]" "RAJ2000 DEJ2000 ${flt}mag Bmag Vmag" \
        "bv/$flt.fits[PHOTOMETRY]" "RAJ2000 DEJ2000 MAG" 3e-4 \
        | awk '$4 < 99 { print $4 - $1, $2 - $3 }' > res.txt
    check "$flt: stars" "$(wc -l < res.txt)" 76 80
    read -r med sd n <<< "$(cut -d' ' -f1 res.txt | stats)"
    check "$flt: median MAG - $flt" "$med" -0.003 0.003
    check "$flt: robust sd" "$sd" 0 0.004
    check "$flt: colour term [mag/mag]" "$(awk '{ x = $2; y = $1; n++; sx += x; sy += y
        sxx += x * x; sxy += x * y } END { print (n * sxy - sx * sy) / (n * sxx - sx * sx) }' res.txt)" \
        -0.003 0.003
done
