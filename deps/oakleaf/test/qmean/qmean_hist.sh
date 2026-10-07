
set -xe

SRC=../../src
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
gfortran $FLAGS $SRC/qsort.f08 $SRC/qmed.f08 $SRC/qmean.f08 \
	 ../noise.f08 ../toolbox.f08 \
	 qmean_hist.f08 -o qmean_hist

./qmean_hist


NDAT=$(cat /tmp/qmean_hist.par)

echo NDAT=$NDAT

gnuplot <<EOF
set output 'qmean_hist_median.png'
set term pngcairo size 800,400
#set term epslatex size 4in,1in solid colour
#set output 'qmean_hist_median.tex'
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "blue" lw 1
#set xtics -2,1 format ""
#set ytics format ""
#set xrange[-2.1:3.1]
sig = 1/sqrt($NDAT)
plot '/tmp/qmean_median' ls 1 with boxes t "rob", exp(-(x/sig)**2/2)/sqrt(2*pi)/sig t "N(0,sig)" ls 2
EOF

gnuplot <<EOF
set output 'qmean_hist_stdsig.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
#set xtics 1.9,0.1,2.2
#set ytics format ""
#set xrange[1.85:2.25]
sig = 1/sqrt(2*$NDAT)
plot '/tmp/qmean_stdsig' ls 1 with boxes t "rob", exp(-((x-1)/sig)**2/2)/sqrt(2*pi)/sig t "N(1,2*sig)"
EOF

gnuplot <<EOF
set output 'qmean_hist_stderr.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
#set xtics 1.9,0.1,2.2
#set ytics format ""
#set xrange[1.85:2.25]
plot '/tmp/qmean_stderr' ls 1 with boxes t "rob"
EOF
