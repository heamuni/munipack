
echo "Laplace distribution:"

NAME='/tmp/laplace'
HISTNAME=${NAME}_hist.dat
PNGNAME=${NAME}_hist.png
CDFDAT=${NAME}_cdf.dat
CDFPNG=${NAME}_cdf.png

laplace $HISTNAME $CDFDAT

# histogram
gnuplot <<EOF
set term png
set output '${PNGNAME}'
set style line 1 lt rgb "salmon3" lw 1
set style line 2 lt rgb "black" lw 1
set style fill solid 0.25 border
set key
set xrange[-5.5:5.5]
#set yrange[0:0.405]
#set xtics nomirror out ("\$-3\$" -3, "\$-2\$" -2, "\$-1\$" -1, "\$0\$" 0,"\$1\$" 1, "\$2\$" 2, "\$3\$" 3)
#set x2tics -3,1,3 format ""
#set ytics ("\$0\$" 0, "\$\\\\sfrac{1}{10}\$" 0.1, "\$\\\\sfrac{2}{10}\$" 0.2, "\$\\\\sfrac{3}{10}\$" 0.3, "\$\\\\sfrac{4}{10}\$" 0.4 )
set ylabel "Frequency"
#plot '${HISTNAME}' with boxes ls 1 t '\$\\epsilon = \\sfrac{1}{10}\$', exp(-x**2/2)/sqrt(2*pi) ls 2 t '\$\\mathcal{L}(0,1)\$
f(x,t,s) = exp(-abs(x-t)/s)/(2*s)
plot '${HISTNAME}' with boxes ls 1, f(x,0,1) ls 2 t 'L(0,1)'
EOF


echo " * histogram saved to $PNGNAME"

# CDF
gnuplot <<EOF
#set term epslatex size 4.5in,2.7in solid color
set term png
set output '$CDFPNG'
set style line 1 lt rgb "black" lw 1
set style line 2 lt rgb "red" lw 2
set key bottom
set xrange[-5.5:5.5]
#set xtics format "\$%.0f\$"
#set xtics -3,1,3
#set ytics ("\$0\$" 0, "\$\\\\sfrac{1}{4}\$" 0.25, "\$\\\\sfrac{1}{2}\$" 0.5, "\$\\\\sfrac{3}{4}\$" 0.75, "\$1\$" 1)
#set xlabel "\$x\$"
set ylabel "Probability"
p(x,t,s) = (1 + sgn(x-t)*(1 - exp(-abs(x-t)/s)))/2
plot '${CDFDAT}' w histeps lc "#52b8ff" lw 2, p(x,0,1) ls 1
EOF

echo " * CDF saved to $CDFPNG"
