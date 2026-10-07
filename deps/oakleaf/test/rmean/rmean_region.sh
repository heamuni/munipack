
set -xe

SRC=../../src/
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
gfortran $FLAGS $SRC/qsort.f95 $SRC/qmean.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 ../noise.f08 robnewton.f08 rmean_region.f08 -o rmean_region

./rmean_region

PSFUN=/tmp/rmean_region
QUANT=/tmp/rmean_region_quant
RECT=$(cat /tmp/rmean_region_conv)
LABELS=/tmp/rmean_region_labels
LABELS0=/tmp/rmean_region_labels0
LABELS1=/tmp/rmean_region_labels1
LABELS2=/tmp/rmean_region_labels2

cat >$LABELS <<EOF
"\\\\footnotesize Trust-region \\\\normalsize" -4 2.5 -90
"\\\\footnotesize Newton \\\\normalsize" 0 3.5 -90
EOF

cat >$LABELS0 <<EOF
\$\\\\varrho(x)\$ -8.7 4.2
EOF

cat >$LABELS1 <<EOF
\$\\\\psi(x)\$ -8.7 2
EOF

cat >$LABELS2 <<EOF
\$\\\\psi'(x)\$ -8.7 2.5
EOF

gnuplot <<EOF
#set term pngcairo size 800,1200
#set output 'rmean_region_tukey.png'
set term epslatex size 4in,3in solid #monochrome
set output 'rmean_region_tukey.tex'
set key left bottom
set xrange[-10:10]
set style line 1 lt rgb "blue" lw 4
set style line 2 lt rgb "black" lw 6
set style fill transparent solid 0.50 border

set multiplot
set size 1,0.3333
set origin 0,0.6666
#set title "Convergence regions for Hampel's redescending estimator"
set tmargin 1
set yrange[-0.1:]
set xtics format ""
set ytics 0,2,9 format "% 0.0f"
set bmargin 0
#set key left bottom
unset key
$RECT

plot '$PSFUN' u 1:2 w l t "\$\\\\varrho(x)\$" ls 1, \
     '$LABELS' u 2:3:1 with labels rotate t "", \
     '$LABELS0' u 2:3:1 with labels t ""

set origin 0,0.3333
set tmargin 0
set autoscale
set xrange[-10:10]
#set xtics format "%0.0f"
set ytics -2,2,2 format "% 0.0f"
#unset key
set key right bottom
unset title
plot '$PSFUN' u 1:3 w l t "" ls 1, '$LABELS1' u 2:3:1 with labels t "" , \
     '$QUANT' u 1:2:3 w xerr t "quantiles" ls 2

set origin 0,0
set tmargin 0
set bmargin 2
set xtics format "%0.0f"
set ytics -2,2,2 format "% 0.0f"
unset key
plot '$PSFUN' u 1:4 w l t "\$\\\\psi'(x)\$" ls 1, \
     '$LABELS2' u 2:3:1 with labels t ""
EOF

for A in *.tex; do
    E=${A%.tex}.eps
    if [ -f $E ]; then
	epstopdf $E
	rm -f $E
    fi
#    sed -i s/".eps"/".pdf"/g ${A}.tex
done
