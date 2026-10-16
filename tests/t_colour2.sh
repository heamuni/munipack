# Two colour terms with real stellar colours: 61 stars of SDSS DR16
# (VizieR V/154/sdss16, primary, clean, class star) in 0.27 x 0.27 deg
# about RA 180, Dec 0, 14.5 < g < 19.5, e_u < 0.08 (u - g 0.2-2.8,
# g - r 0.0-1.5, correlated by 0.96 as for any field). The frame in g
# (Gunn) is rendered from X = g + 0.3 (u - g) - 0.2 (g - r) of
# the catalogue, a band of an exactly known two-colour relation, and
# calibrated against the catalogue g with the colours u - g and g - r:
# m_cat - m_inst = ZP + k1 (u - g) + k2 (g - r), k1 = -0.3, k2 = +0.2
# exactly (one colour cannot follow it: 0.02-0.03 mag left about
# the best function of u - g). Three stars have no u in the catalogue
# used for the calibration: their u - g is predicted from g - r.

# RA, Dec [deg], u, e_u, g, e_g, r, e_r (SDSS DR16, Ahumada et al. 2020)
cat > rows.txt <<'END'
180.062985 -0.026339 15.903 0.031 14.687 0.018 14.291 0.017
180.132989 -0.131980 15.964 0.044 14.768 0.019 14.357 0.017
180.009773 -0.037720 17.073 0.032 15.300 0.018 14.641 0.017
180.124721 0.112205 17.752 0.020 15.452 0.021 14.423 0.014
179.966993 0.130989 16.004 0.013 15.788 0.019 15.800 0.013
180.073485 -0.037063 18.037 0.033 15.842 0.018 15.061 0.017
179.899777 -0.118683 17.333 0.019 16.180 0.016 15.811 0.017
180.039175 0.008801 17.340 0.019 16.216 0.021 15.889 0.014
179.917076 -0.115487 17.914 0.020 16.223 0.016 15.641 0.017
179.936477 -0.089697 19.030 0.024 16.334 0.016 15.159 0.017
180.025402 0.058744 17.790 0.020 16.347 0.021 15.793 0.014
179.930830 -0.027583 17.930 0.020 16.349 0.016 15.785 0.017
180.064752 -0.112689 18.183 0.033 16.364 0.018 15.665 0.017
179.946334 0.079222 17.947 0.018 16.662 0.019 16.111 0.013
180.069173 -0.071447 19.352 0.037 16.727 0.018 15.536 0.017
179.936652 -0.079587 18.181 0.021 16.836 0.016 16.373 0.017
180.117809 -0.103409 19.337 0.037 16.889 0.018 15.391 0.017
180.044037 -0.105999 18.164 0.033 16.937 0.018 16.463 0.017
179.984557 -0.001032 18.269 0.033 17.058 0.018 16.766 0.017
180.116588 -0.085720 19.304 0.037 17.097 0.018 16.199 0.017
179.902446 -0.127919 19.910 0.032 17.126 0.016 15.719 0.017
179.980871 -0.003932 18.178 0.021 17.156 0.016 16.746 0.017
180.092183 0.107178 18.636 0.025 17.349 0.021 16.766 0.014
179.984593 -0.128464 18.385 0.033 17.415 0.018 17.100 0.017
179.943092 -0.003132 18.500 0.022 17.498 0.016 17.185 0.017
180.052720 -0.014468 18.824 0.034 17.557 0.018 17.111 0.017
179.987833 0.096221 19.931 0.042 17.638 0.020 16.638 0.013
179.942496 -0.106783 20.306 0.039 17.677 0.016 16.540 0.017
180.040190 -0.045381 20.369 0.048 17.710 0.018 16.324 0.017
179.878392 0.096189 20.441 0.057 17.741 0.020 16.284 0.013
179.921126 0.002214 19.089 0.035 17.804 0.020 17.261 0.013
180.096719 -0.047520 18.902 0.035 17.842 0.018 17.499 0.017
180.082901 -0.096513 20.176 0.045 17.926 0.018 17.095 0.017
179.871649 0.063029 19.560 0.034 18.008 0.020 17.326 0.013
179.953318 -0.051059 20.796 0.052 18.033 0.016 16.630 0.017
179.943479 0.109439 19.300 0.030 18.056 0.020 17.596 0.013
180.121003 -0.057730 19.456 0.037 18.092 0.018 17.602 0.017
179.907509 0.121133 20.831 0.076 18.108 0.020 16.810 0.013
179.921150 -0.055557 18.696 0.022 18.210 0.016 18.118 0.017
180.119706 -0.011615 19.551 0.038 18.248 0.019 17.814 0.018
180.074414 -0.071685 20.034 0.043 18.331 0.019 17.745 0.018
180.028290 0.016018 19.365 0.033 18.344 0.021 18.056 0.015
179.931392 0.101432 20.792 0.073 18.344 0.020 17.278 0.013
180.132060 -0.054354 19.329 0.048 18.375 0.019 18.086 0.018
180.006280 -0.035553 19.470 0.037 18.377 0.019 17.934 0.018
179.940788 -0.077720 20.937 0.057 18.398 0.016 17.373 0.017
180.124248 -0.052382 19.893 0.041 18.430 0.019 17.770 0.018
179.948948 -0.103959 20.935 0.057 18.480 0.016 17.477 0.017
179.942781 -0.115114 20.879 0.056 18.701 0.017 17.811 0.017
179.981544 -0.065624 20.319 0.047 18.774 0.019 18.139 0.018
179.907568 0.036407 19.825 0.040 18.885 0.021 18.557 0.015
180.079117 -0.127770 20.812 0.059 18.939 0.019 18.221 0.018
179.869545 -0.045362 19.947 0.033 19.020 0.017 18.581 0.018
179.975284 0.083483 20.834 0.076 19.074 0.021 18.343 0.014
179.999710 -0.119677 20.235 0.046 19.138 0.020 18.682 0.018
179.955711 -0.106917 20.313 0.040 19.226 0.018 18.852 0.019
180.116628 0.112413 20.238 0.054 19.307 0.023 18.945 0.017
179.895045 -0.012717 20.252 0.039 19.311 0.018 18.982 0.019
180.038981 -0.082160 20.384 0.048 19.375 0.020 19.149 0.019
179.966404 0.004967 20.266 0.053 19.403 0.022 19.089 0.016
180.023111 0.066639 20.621 0.068 19.412 0.023 18.954 0.017
END
awk '{ x = $5 + 0.3 * ($3 - $5) - 0.2 * ($5 - $7)
       print $1, $2, 0, 0, $3, $4, $5, $6, $7, $8, x }' rows.txt > rowsx.txt
tabfits truth.fits "RAJ2000 D deg DEJ2000 D deg pmRA D mas/yr pmDE D mas/yr \
    upmag E - e_upmag E - gpmag E - e_gpmag E - rpmag E - e_rpmag E - Xmag E -" < rowsx.txt
# the catalogue of the calibration: u missing for the stars 5, 25, 45
awk 'NR == 5 || NR == 25 || NR == 45 { $5 = 99.999; $6 = 9.999 } { print }' rowsx.txt \
    | tabfits cat.fits "RAJ2000 D deg DEJ2000 D deg pmRA D mas/yr pmDE D mas/yr \
    upmag E - e_upmag E - gpmag E - e_gpmag E - rpmag E - e_rpmag E - Xmag E -"
check "catalogue stars" "$(nrows "cat.fits[1]")" 61 61

munipack artificial -c truth.fits --col-ra RAJ2000 --col-dec DEJ2000 \
    --col-pm-ra pmRA --col-pm-dec pmDE --col-mag Xmag -f g --photsys Gunn \
    --rcen 180 --dcen 0 --scale 3600 --width 1024 --height 1024 \
    --exptime 60 --psf GAUSS --hwhm 1.4 --gain 2 --readnoise 4 --bias 100 \
    --seed 9 --no-watermark --mask x.fits > /dev/null 2>&1
munipack find -f 3.3 x.fits > /dev/null 2>&1
munipack aphot x.fits > /dev/null 2>&1

C="--photsys-ref Gunn -f g --col-ra RAJ2000 --col-dec DEJ2000 -c cat.fits \
    --col-mag gpmag --col-magerr e_gpmag"
# resid DIR: MAG - g, u - g, g - r, MAGERR, COLORSRC of the stars with
# a catalogue u (the true u for the others)
resid() {
    match "truth.fits[1]" "RAJ2000 DEJ2000 upmag gpmag rpmag" \
        "$1/x.fits[PHOTOMETRY]" "RAJ2000 DEJ2000 MAG MAGERR COLORSRC" 3e-4 \
        | awk '$4 < 99 { print $4 - $2, $1 - $2, $2 - $3, $5, $6 }'
}
slope() {   # slope COLUMN < resid: the slope of MAG - g against a colour
    awk -v c="$1" '{ x = $c; y = $1; n++; sx += x; sy += y; sxx += x * x; sxy += x * y }
        END { print (n * sxy - sx * sy) / (n * sxx - sx * sx) }'
}

# the colours without errors (the catalogue is exact here): the plain fit
mkdir t2 t2e t1
munipack phcal $C --color-term upmag,gpmag --color-term gpmag,rpmag \
    -t t2 x.fits > /dev/null 2>&1
h="t2/x.fits[PHOTOMETRY]"
k1=$(key "$h" COLTK1); e1=$(key "$h" COLTK1E)
k2=$(key "$h" COLT2K1); e2=$(key "$h" COLT2K1E)
check "COLTNTRM" "$(key "$h" COLTNTRM)" 2 2
check "COLTERM2 is g - r" "$(key "$h" COLTERM2 | grep -c 'gpmag-rpmag')" 1 1
check "COLTK1 (-0.3)" "$k1" -0.33 -0.27
check "COLT2K1 (+0.2)" "$k2" 0.15 0.25
check "(COLTK1 + 0.3) / COLTK1E" "$(awk -v k="$k1" -v e="$e1" 'BEGIN { print (k + 0.3) / e }')" -3 3
check "(COLT2K1 - 0.2) / COLT2K1E" "$(awk -v k="$k2" -v e="$e2" 'BEGIN { print (k - 0.2) / e }')" -3 3
check "COLTC12 (u - g and g - r correlated)" "$(key "$h" COLTC12)" -1 -0.5
check "COLORSRC, COLOR2 columns" "$(colnum "$h" COLOR2ERR)" 1 99
resid t2 > res.txt
check "stars" "$(wc -l < res.txt)" 55 61
read -r med sd2 n <<< "$(awk '$5 == 3 { print $1 }' res.txt | stats)"
check "MAG - g: robust sd" "$sd2" 0 0.012
check "slope against u - g" "$(awk '$5 == 3' res.txt | slope 2)" -0.01 0.01
check "slope against g - r" "$(awk '$5 == 3' res.txt | slope 3)" -0.02 0.02
# the stars without u: COLORSRC 2 (g - r only), u - g predicted from g - r
check "stars without u (COLORSRC 2)" "$(awk '$5 == 2' res.txt | wc -l)" 2 3
check "their |MAG - g - median| / MAGERR (max)" "$(awk -v m="$med" '$5 == 2 {
    z = ($1 - m) / $4; if (z < 0) z = -z; if (z > zm) zm = z } END { print zm + 0 }' res.txt)" 0 3

# with the errors of the colours: their correlation (g shared) in the errors
# (the catalogue is exact here, so the correction for the errors of
# the colours, BUGS P57, over-corrects a little: lambda 0.99, within
# the checks; a real catalogue measured with these errors is the case of
# t_colerr)
munipack phcal $C --color-term upmag,gpmag,e_upmag,e_gpmag \
    --color-term gpmag,rpmag,e_gpmag,e_rpmag -t t2e x.fits > /dev/null 2>&1
h="t2e/x.fits[PHOTOMETRY]"
check "errors: (COLTK1 + 0.3) / COLTK1E" "$(awk -v k="$(key "$h" COLTK1)" \
    -v e="$(key "$h" COLTK1E)" 'BEGIN { print (k + 0.3) / e }')" -3 3
check "errors: (COLT2K1 - 0.2) / COLT2K1E" "$(awk -v k="$(key "$h" COLT2K1)" \
    -v e="$(key "$h" COLT2K1E)" 'BEGIN { print (k - 0.2) / e }')" -3 3

# one colour (u - g) leaves the g - r part of the relation
munipack phcal $C --color-term upmag,gpmag -t t1 x.fits > /dev/null 2>&1
[ -z "$(key "t1/x.fits[PHOTOMETRY]" COLTNTRM)" ] || fail "COLTNTRM with one colour"
ok "no keys of a second colour with one --color-term"
read -r med sd1 n <<< "$(resid t1 | awk '$5 == 1 { print $1 }' | stats)"
check "one colour: robust sd" "$sd1" "$(awk -v s="$sd2" 'BEGIN { print 2 * s }')" 0.1

expect_fail "dependent colours" munipack phcal $C --color-term upmag,gpmag \
    --color-term gpmag,upmag x.fits
expect_fail "four colour terms" munipack phcal $C --color-term upmag,gpmag \
    --color-term gpmag,rpmag --color-term rpmag,Xmag --color-term upmag,Xmag x.fits
expect_fail "--color-order of three values for two terms" munipack phcal $C \
    --color-term upmag,gpmag --color-term gpmag,rpmag --color-order 1,2,1 x.fits
