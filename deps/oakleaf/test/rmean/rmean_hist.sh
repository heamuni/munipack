
set -xe

SRC=../../src/
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
gfortran $FLAGS $SRC/qsort.f95 $SRC/qmean.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 ../noise.f08 ../toolbox.f08 robnewton.f08 rmean_hist.f08 -o rmean_hist

./rmean_hist


NDAT=$(cat /tmp/rmean_hist.par)

echo NDAT= $NDAT

gnuplot <<EOF
#set output 'rmean_hist_t.png'
#set term pngcairo size 800,400
set term epslatex size 4in,1in solid colour #monochrome
set output 'rmean_hist_t.tex'
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "blue" lw 1
set xtics -3,1 format "" #format "\$%.0f\$"
#set ytics 0,1,4  format "\$%.0f\$"
set ytics format ""
#set tmargin 0.1
#set bmargin 0
#set rmargin 0.1
#set lmargin 2
#set xrange[-1.5:2.25]
set xrange[-2.1:3.1]
sig = 1/sqrt($NDAT)
plot '/tmp/rmean_t' ls 1 with boxes t "rob", exp(-(x/sig)**2/2)/sqrt(2*pi)/sig t "N(0,sig)" ls 2
EOF

gnuplot <<EOF
set output 'rmean_hist_sig.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
#set xtics 1.9,0.1,2.2
#set ytics format ""
#set xrange[1.85:2.25]
sig = 1/sqrt(2*$NDAT)
plot '/tmp/rmean_stdsig' ls 1 with boxes t "rob", exp(-((x-1)/sig)**2/2)/sqrt(2*pi)/sig t "N(1,2*sig)"
EOF

gnuplot <<EOF
set output 'rmean_hist_s.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
#set xtics 1.9,0.1,2.2
#set ytics format ""
#set xrange[1.85:2.25]
plot '/tmp/rmean_stderr' ls 1 with boxes t "rob"
EOF

#for A in rmean_hist_t; do
#    E=${A}.eps
#    if [ -f $E ]; then
#        epstopdf $E
#        rm -f $E
#	mv ${A}.tex ${A}.pdf ~/rontrib/
#    fi
#done
