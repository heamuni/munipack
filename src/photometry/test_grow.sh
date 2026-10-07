
LS="art_??.fits"
LS="0716_*V.fits"


# construction of growth-curve as the product of calibration in various apertures
for A in $LS; do
    B=${A%fits}dat
    CAL=${A%.fits}_catcal.fits
    rm -f $B
    for I in `seq 1 12`; do
	munipack phcal --photsys-ref Johnson -f V --col-mag f.mag --col-magerr e_Vmag \
		 -c cone.fits -O --mask '\!\1_catcal.\2' --verbose --saper $I $A;
#	M=$(munipack fits -K AIRMASS --shell ${CAL} | awk -F= '{print $2;}')
#	E=$(munipack fits -lh ${CAL} | \
#		   awk -v a=$M '{if( /Extinction/ ) print exp($5*(a-1));}')
	C=$(munipack fits -K CTPH --shell ${CAL}\[PHOTOMETRY\] | \
		   awk -F= -v e=$E '{print 1/$2;}')
#	R=$(munipack fits -K APER --shell ${CAL}\[PHOTOMETRY\] | \
#		   awk -F= -v s=$SCALE '{print $2*3600/s;}')
	R=$(munipack fits -lh ${A}\[APERPHOT\] | \
		   awk -v i=$I '{if($1 ~ "^APER" i "$") print $3;}')
	echo $R $C >> $B
    done
    X=$(awk '{} END {print $2;}' < $B)
    awk -v X=$X '{print $1,$2/X;}' < $B > /tmp/x
    mv /tmp/x $B
done


# growth-curve photometry
munipack gphot --verbose $LS
for A in $LS; do
    munipack phcal --photsys-ref Johnson -f V --col-mag f.mag --col-magerr e_Vmag \
	     -c cone.fits -O --mask '\!\1_catcal.\2' --verbose $A;
done

rm -f ext_G
for L in $LS; do
    A=${L%.fits}_catcal.fits
    D=${A%fits}dat
    E=${A%_catcal.fits}.dat
    X=$(munipack fits -K AIRMASS --shell $A | awk -F= '{print $2;}')
    C=$(munipack fits -K CTPH --shell ${A}\[PHOTOMETRY\] | awk -F= '{print $2;}')
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
    set ytics -0.05,0.05,0.05 format '%0.2f'
#    set yrange[-0.008:0.008]
    set yrange[-0.06:0.06]
    plot 'F' u 1:(\$2-\$6)
EOF
done
