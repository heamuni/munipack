
# A variable star
RA=150
DEC=50
JD0=2457948.5
LCID="--lc-fourier fourier.fits"
LCID="--lc-table Kepler-21b.fits"
LCID="--lc-model Cepheid"
munipack artificial --verbose --mask 'art_??.fits' --fov 0.3 --count 48 \
	 --rcen $RA --dcen $DEC \
	 --lc-ra $RA --lc-dec $DEC \
 	 $LCID --lc-mag 12 --lc-amp 1 --lc-jd0 $JD0 --lc-per 0.1 \
	 --date 2017-07-14 --time 00:00:00  --timestep 180 --exptime 90
LS="art_??.fits"
munipack find $LS
munipack aphot $LS
cp art_01.fits art_01_grow.fits
munipack gphot art_01_grow.fits

# aperture correction
R90=$(munipack fits --shell -K RADFLX90 art_01_grow.fits\[GROWPHOT\] | \
	     awk -F '=' '{print $2;}')
BEGIN='BEGIN { q = -1; }'
BODY='{if($0~/^APER/ && $2 >= r && q < 0) {q=$2;gsub("APER","",$1);print $1;}}'
SAPER=$(munipack fits -lh art_01_grow.fits\[APERPHOT\] | \
	       awk -F '[=/]' -v r=$R90 "$BEGIN $BODY")
APCORR=$(munipack fits -lt art_01_grow.fits\[GROWCURVE\] | \
	       awk -v saper=$SAPER  '{if( NR == saper ) print 1/$2;}')

munipack phcal -C 1 --saper $SAPER --photsys-ref Johnson -f V \
	 art_01_grow.fits,man_01.fits
for A in $LS; do
    munipack phcal --photsys-ref Johnson -f V -r man_01.fits \
	     --apcorr $APCORR --saper $SAPER $A;
done
echo "# JD MAG " > lc
munipack timeseries -c MAG,MAGERR --stdout ${RA},${DEC} art_*_cal.fits >> lc

echo "# JD MAG " > lcmodel
BEGIN='BEGIN { JD=0; MAG=99.999; }'
BODY='{ if($0~/^JD/) JD=$2; if($0~/current magnitude/) MAG=$2; }'
END='END {print JD MAG;}'
for A in $LS; do
    munipack fits -lh $A | awk -F '[=/]' "$BEGIN $BODY $END" >> lcmodel
done

paste lc lcmodel > lcd

gnuplot <<EOF
set key right bottom
set xlabel "Julian Day - $JD0"
set ylabel "magnitude"
set title "Model $LCID"
#set yrange[13.7:11.7]
set size ratio 0.618
set term svg dynamic
set output 'lc.svg'
plot 'lc' u (\$1-$JD0):2 t "data art" ls 3 lc rgb '#8AB8E6', \
     'lcmodel' u (\$1-$JD0):2 t "model" w l lc rgb '#000000'
EOF

gnuplot <<EOF
set key right bottom
set xlabel "Julian Day - $JD0"
set ylabel "difference magnitude"
set title "Verification $LCID"
#set yrange[13.7:11.7]
set size ratio 0.618
set term svg dynamic
set output 'lcd.svg'
plot 'lcd' u (\$1-$JD0):(\$2-\$5) t "diff." ls 3 lc rgb '#8AB8E6'
EOF

echo "Result kept in 'lc.svg' and 'lcd.svg'."

exit 0


# BL Lac
#RA=330.68
#DEC=42.27
#munipack cone -r 0.2 -- $RA $DEC
#munipack artificial --verbose --fov 0.3 -c cone.fits --rcen $RA --dcen $DEC \
#	 --mask 'art_??.fits' --lat 50 --long -15 --date 2016-08-29 \
#	 --time 00:27:03.123 --count 50 --extk 0.1 --timestep 600 --exptime 10 \
#	 --atmosphere

# sparse in UMa
HWHM=0.5
RA=140
DEC=40
COUNT=55
TIMESTEP=600
#COUNT=5
#TIMESTEP=6600
munipack cone -r 0.2 -- $RA $DEC
munipack artificial --verbose --fov 0.3 -c cone.fits --rcen $RA --dcen $DEC \
         --mask 'art_??.fits' --lat 50 --long -15 --date 2016-11-29 \
         --time 04:47:03.123 --count $COUNT --timestep $TIMESTEP --exptime 60 \
         --hwhm $HWHM --extk 0.1 --atmosphere

SCALE=$(echo 0.3*3600/666 | bc -l) # update for FOV, size changes
LS="art_??.fits"

munipack find -f 2 $LS
munipack aphot $LS

# construction of growth-curve as product of calibration in various apertures
for A in $LS; do
    B=${A%fits}dat
    CAL=${A%.fits}_catcal.fits
    rm -f $B
    for I in `seq 1 12`; do
	munipack phcal --photsys-ref Johnson -f V --col-mag f.mag \
		 --col-magerr e_Vmag -c cone.fits -O --mask '\!\1_catcal.\2' \
		 --verbose --saper $I $A;
	M=$(munipack fits -K AIRMASS --shell ${CAL} | awk -F= '{print $2;}')
	E=$(munipack fits -lh ${CAL} | \
		   awk -v a=$M '{if( /Extinction/ ) print exp($5*(a-1));}')
	C=$(munipack fits -K CTPH --shell ${CAL}\[PHOTOMETRY\] | \
		   awk -F= -v e=$E '{print 1/$2*e;}')
	R=$(munipack fits -K APER --shell ${CAL}\[PHOTOMETRY\] | \
		   awk -F= -v s=$SCALE '{print $2*3600/s;}')
	echo $R $C >> $B
    done
done

# aperture photometry
for A in $LS; do
    munipack phcal --photsys-ref Johnson -f V --col-mag f.mag --col-magerr e_Vmag \
	     -c cone.fits -O --mask '\!\1_catcal.\2' --verbose $A;
done
rm -f ext_A
for A in art_*_catcal.fits; do
    X=$(munipack fits -K AIRMASS --shell $A | awk -F= '{print $2;}')
    C=$(munipack fits -K CTPH --shell ${A}\[PHOTOMETRY\] | awk -F= '{print $2;}')
    echo $X $C >> ext_A
done

# growth-curve photometry
munipack gphot --verbose $LS
for A in $LS; do
    munipack phcal --photsys-ref Johnson -f V --col-mag f.mag --col-magerr e_Vmag \
	     -c cone.fits -O --mask '\!\1_catcal.\2' --verbose $A;
done
rm -f ext_G
for A in art_*_catcal.fits; do
    D=${A%fits}dat
    E=${A%_catcal.fits}.dat
    X=$(munipack fits -K AIRMASS --shell $A | awk -F= '{print $2;}')
    C=$(munipack fits -K CTPH --shell ${A}\[PHOTOMETRY\] | \
	       awk -F= '{print $2;}')
    echo $X $C >> ext_G
    munipack fits -lt ${A%_catcal.fits}.fits\[GROWCURVE\] > $D

    paste $D $E > F
    PNG=${A%fits}png
    gnuplot <<EOF
    set output '$PNG'
    set term png
    set multiplot title 'Grows'
    set lmargin 6
    set rmargin 2
    set origin 0,0.3
    set size 1,0.7
    set tmargin 2
    set bmargin 0
    unset xtics
    set key bottom
    set yrange[0:1.2]
    set xrange[0:21]
    #set ylabel "attenuation"
    plot '$D' u 1:2 t "curve", '$E' t "by cal"
    set origin 0,0
    set size 1,0.3
    set tmargin 0
    set bmargin 3
    unset key
    set grid
    set xlabel "aperture [pix]"
    set xtics 0,5,20 format '%0.0f' nomirror
    set ytics -0.005,0.005,0.005 format '%0.3f'
    set yrange[-0.008:0.008]
    plot 'F' u 1:(\$2-\$6)
EOF
done


gnuplot <<EOF
set output 'att.png'
set term png
set multiplot title 'Extinction'
set lmargin 6
set rmargin 2
set origin 0,0.3
set size 1,0.7
set tmargin 2
set bmargin 0
unset xtics
set key bottom
set ytics 1,0.1,1.5 format '%0.1f'
set yrange[0.9:1.5]
set xrange[0.9:4.5]
set ylabel "attenuation"
plot 'ext_A' t "aperture", 'ext_G' t "grow", exp(0.1*(x-1))
set origin 0,0
set size 1,0.3
set tmargin 0
set bmargin 3
unset key
set xtics 1,1,5 format '%0.0f' nomirror
set grid
set ytics -0.01,0.01,0.1 format '%0.2f'
set yrange[-0.011:0.011]
set xlabel "airmass"
plot 'ext_A' u 1:(exp(0.1*(\$1-1))-\$2), 'ext_G' u 1:(exp(0.1*(\$1-1))-\$2)
EOF
