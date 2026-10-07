
set -xe

SRC=../../src/
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
gfortran $FLAGS $SRC/qsort.f95 $SRC/qmean.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 ../noise.f08 rfun_multimin.f08 -o rfun_multimin

gfortran $FLAGS $SRC/qsort.f95 $SRC/qmean.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 ../noise.f08 rfun_secmin.f08 -o rfun_secmin

gfortran $FLAGS $SRC/qsort.f95 $SRC/qmean.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 ../noise.f08 robnewton.f08 rfun_multimed.f08 -o rfun_multimed

gfortran $FLAGS $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 ../noise.f08 rfun_multi3d.f08 -o rfun_multi3d


./rfun_multimin

./rfun_secmin

#./rfun_multimed

#./rfun_multi3d

MFUN='/tmp/rfun_multifun'

gnuplot <<EOF
#set output 'rfun_multifun.png'
#set term pngcairo size 400,400
set term epslatex size 2.8in,1.5in solid monochrome
set output 'rfun_multifun.tex'
#set key bottom right
unset key
set tmargin 0
set rmargin 0.1
set bmargin 0
set lmargin 2
set style line 1 lt rgb "black" lw 3
set ytics 0,5,10 format "\$%.0f\$"
set xtics -10,10,10 format "\$%.0f\$"
set xrange[-20.2:20.2]
set yrange[-0.1:11]
plot '$MFUN' w l ls 1 t 'rfun', for [i=3:7] '$MFUN' u 1:i w l t ""
EOF


gnuplot <<EOF
#set output 'rfun_convolve.png'
#set term pngcairo size 400,400
set term epslatex size 2.8in,1.5in solid monochrome
set output 'rfun_convolve.tex'
unset key
set tmargin 0
set rmargin 0.1
set bmargin 0
set lmargin 2
set style line 1 lt rgb "black" lw 3
set ytics 0,5,10  format "\$%.0f\$"
set xtics -10,10,10 format "\$%.0f\$"
set xrange[-20.2:20.2]
set yrange[-0.1:11]
plot '/tmp/rfun_secmin_convol' w l ls 1
EOF

gnuplot <<EOF
#set output 'multimin3d.png'
#set term pngcairo size 640,400
set term epslatex size 2.8in,2in solid monochrome
set output 'multimin3d.tex'
set style fill transparent solid 0.25 border
#set style fill solid
set boxdepth 0.3
set pm3d depthorder base
set walls z0
set xyplane at 0
set view 60, 40
set logscale x
set xrange[1:2e3]
set tmargin 0
set rmargin 0.1
set bmargin -2
set lmargin 2
set xtics format "\$%.0f\$"
set ytics 0,10,30 format "\$%.0f\$"
set ztics 1,1,2 format "\$%.0f\$"
set xlabel "\$N\$"
set ylabel "\$t_c\$"
unset key
splot '/tmp/rfun_multi3d' w boxes
EOF

for A in multimin3d rfun_multifun rfun_convolve; do
    E=${A}.eps
    if [ -f $E ]; then
        epstopdf $E
        rm -f $E
    fi
    mv ${A}.tex ${A}.pdf $HOME/rontrib/
done
