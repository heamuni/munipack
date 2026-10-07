
set -xe

SRC=../../src/
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
gfortran $FLAGS $SRC/qsort.f95 $SRC/qmean.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 ../noise.f08 robnewton.f08 rmeanfun.f08 -o rmeanfun

./rmeanfun

PSHAMPEL=/tmp/rmeanfun
TUKEY=/tmp/rmeanfun_tukey
QUANT=/tmp/rmeanfun_quant
ECDF=/tmp/rmean_ecdf
ECDFmad=/tmp/rmean_ecdf_mad

gnuplot <<EOF
set term pngcairo size 800,400
set output 'rmeanfun.png'
set key right bottom
plot '$PSHAMPEL' u 1:3 w l, '$TUKEY' u 1:3 w l, \
     '$QUANT' w l, '$QUANT' u 3:4 w l, '$QUANT' u 5:6 w l
EOF

gnuplot <<EOF
set term pngcairo size 800,1200
set output 'rmeanfun_multi.png'
set key right bottom
set xrange[-10:10]

set multiplot
set size 1,0.3333
set origin 0,0.6666
set xtics format ""
set ytics 0,2,9 format "% 0.0f"
set tmargin 1
set bmargin 0
plot '$PSHAMPEL' u 1:2 w l t "Hampel", '$TUKEY' u 1:2 w l t "Tukey"

set origin 0,0.3333
set tmargin 0
#set xtics format "%0.0f"
#set ytics -2,1,2 format "% 0.0f"
unset key
plot '$PSHAMPEL' u 1:3 w l t "Hampel", '$TUKEY' u 1:3 w l t "Tukey", \
     '$QUANT' w l, '$QUANT' u 3:4 w l, '$QUANT' u 5:6 w l

set origin 0,0
set tmargin 0
set bmargin 2
set xtics format "%0.0f"
#set ytics -2,1,2 format "% 0.0f"
unset key
plot '$PSHAMPEL' u 1:4 w l t "Hampel", '$TUKEY' u 1:4 w l t "Tukey"
EOF

gnuplot <<EOF
set term pngcairo size 800,400
set output 'rmean_cdf.png'
set key right bottom
set xrange[-3.1:3.1]
plot '$ECDF' w steps, '$ECDFmad' w steps
EOF
