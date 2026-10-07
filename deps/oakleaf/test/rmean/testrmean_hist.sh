
set -xe

SRC=../../src/
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
gfortran $FLAGS $SRC/qsort.f95 $SRC/qmean.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 ../noise.f08 ../toolbox.f08 robnewton.f08 rmean_hist.f08 -o rmean_hist

./rmean_hist


NDAT=$(cat /tmp/rmean_hist.par)

echo NDAT= $NDAT

gnuplot <<EOF
set output 'rmean_hist_t.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
#set xtics 1.9,0.1,2.2
#set ytics format ""
#set xrange[1.85:2.25]
sig = 1/sqrt($NDAT)
plot '/tmp/rmean_t' ls 1 with boxes t "rob", exp(-(x/sig)**2/2)/sqrt(2*pi)/sig t "N(0,sig)"
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
