
set -xe

SRC=../../src/
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
gfortran $FLAGS $SRC/qsort.f95 $SRC/qmean.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 $SRC/fmm/fmm.f08 $SRC/fmm/zeroin.f95  ../noise.f08 ../toolbox.f08 robnewton.f08 rmean_newton.f08 -o rmean_newton

./rmean_newton

ITER00_5='/tmp/rmean_newton_iter_00_5'
ITER00_10='/tmp/rmean_newton_iter_00_10'
ITER00_20='/tmp/rmean_newton_iter_00_20'
ITER00_100='/tmp/rmean_newton_iter_00_100'
ITER00_10000='/tmp/rmean_newton_iter_00_10000'

ITER20_5='/tmp/rmean_newton_iter_20_5'
ITER20_10='/tmp/rmean_newton_iter_20_10'
ITER20_20='/tmp/rmean_newton_iter_20_20'
ITER20_100='/tmp/rmean_newton_iter_20_100'
ITER20_10000='/tmp/rmean_newton_iter_20_10000'


gnuplot <<EOF
#set term pngcairo size 800,400
#set output 'rmean_newton_iter_00.png'
set term epslatex size 2.8in,1.4in solid monochrome
set output 'rmean_newton_iter_00.tex'

#unset key
#unset colorb
#set colorb horizontal user origin .05, .05 size .9, .05
#set palette cubehelix
set grid y
set boxwidth 0.25 abs
set style fill solid border lt -1

set xrange[0.5:8]
set ytics 0,50,50
set yrange[0:100]
unset key
#set xtics ""
#set xtics scale 0 ()

#riter=15
#category = "5 10 20 100 10000"
#ClusterSize = 1
#set for [i=1:riter] xtics add (i 5+(i-1)*ClusterSize)
#xcoord(i) =  i*ClusterSize + column(1)
#color(i)  = i == 0 ? 0 : i == 1 ? 0.5 : 1.0

set multiplot
#set xtics format ""
set size 1,0.333
set tmargin 0
set bmargin 0

set origin 0,0
#set xtics
set title offset -30,-2.5 "N=5"
plot '$ITER00_5' u (\$1-0.25):2 w boxes t "tr", \
     '$ITER00_5' u (\$1-0.0):3 w boxes t "r",\
     '$ITER00_5' u (\$1+0.25):4 w boxes t "z",

set origin 0,0.333
set bmargin 0
set xtics format ""
set title offset -30,-2.5 "N=10"
plot '$ITER00_10' u (\$1-0.25):2 w boxes t "tr", \
     '$ITER00_10' u (\$1-0.0):3 w boxes t "r",\
     '$ITER00_10' u (\$1+0.25):4 w boxes t "z"

#set origin 0,0.666
#plot '$ITER00_20' u (\$1-0.25):2 w boxes t "tr", \
#     '$ITER00_20' u (\$1-0.0):3 w boxes t "r",\
#     '$ITER00_20' u (\$1+0.25):4 w boxes t "z"

set origin 0,0.666
set tmargin 0.5
set key
set title offset -30,-2.5 "N=100"
plot '$ITER00_100' u (\$1-0.25):2 w boxes t "tr", \
     '$ITER00_100' u (\$1-0.0):3 w boxes t "r",\
     '$ITER00_100' u (\$1+0.25):4 w boxes t "z"

#set origin 0,0.8
#set tmargin 0.5
#set key
#plot '$ITER00_10000' u (\$1-0.25):2 w boxes t "tr", \
#     '$ITER00_10000' u (\$1-0.0):3 w boxes t "r",\
#     '$ITER00_10000' u (\$1+0.25):4 w boxes t "z"

#plot for [i=0:7] '$ITER00' \
#     using (xcoord(i)):(column(i+2)):(column(1)/3) with boxes lc palette z

EOF


gnuplot <<EOF
#set term pngcairo size 800,400
#set output 'rmean_newton_iter_20.png'
set term epslatex size 2.8in,1.4in solid monochrome
set output 'rmean_newton_iter_20.tex'

set grid y
set boxwidth 0.25 abs
set style fill solid border lt -1

set xrange[0.5:8]
set ytics 0,50,50
set yrange[0:100]
unset key

set multiplot
set size 1,0.333
set tmargin 0
set bmargin 0

set origin 0,0
set title offset -30,-2.5 "\$N=5\$"
plot '$ITER20_5' u (\$1-0.25):2 w boxes t "tr", \
     '$ITER20_5' u (\$1-0.0):3 w boxes t "r",\
     '$ITER20_5' u (\$1+0.25):4 w boxes t "z",

set origin 0,0.333
set bmargin 0
set xtics format ""
set title offset -30,-2.5 "\$N=10\$"
plot '$ITER20_10' u (\$1-0.25):2 w boxes t "tr", \
     '$ITER20_10' u (\$1-0.0):3 w boxes t "r",\
     '$ITER20_10' u (\$1+0.25):4 w boxes t "z"

set origin 0,0.666
set tmargin 0.5
set key
set title offset -30,-2.5 "\$N=100\$"
plot '$ITER20_100' u (\$1-0.25):2 w boxes t "tr", \
     '$ITER20_100' u (\$1-0.0):3 w boxes t "r",\
     '$ITER20_100' u (\$1+0.25):4 w boxes t "z"

EOF

for A in rmean_newton_iter_00 rmean_newton_iter_20; do
    E=${A}.eps
    if [ -f $E ]; then
        epstopdf $E
        rm -f $E
    fi
    mv ${A}.tex ${A}.pdf ~/rontrib/
done
