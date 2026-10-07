
set -xe

FFLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
SRC=../../src/
gfortran $FFLAGS $SRC/cholesky.f08 $SRC/trust_step.f08 $SRC/sr1.f08  $SRC/rfun.f08 ../noise.f08  -o testsr1line  testsr1line.f08

./testsr1line

FCN=/tmp/sr1line
STEP=/tmp/sr1line_step
ITER=/tmp/sr1line_approx_2


gnuplot <<EOF

unset key
set view map
set pm3d at b
#set samples 100; set isosamples 100
set cntrparam levels incr 0,1,30
#set cntrlabel onecolor
set palette rgb 33,13,10
set term pngcairo size 800,800
set output 'testsr1line.png'
set title "Test SR1 Tukey line"

set yrange[-5:5]
set xrange[-5:5]
#set zrange[:20]
set multiplot

set lmargin at screen 0.1
set rmargin at screen 0.85
set tmargin at screen 0.1
set bmargin at screen 0.9

set contour base
#set origin 0,0
#set size 1,1
#set linetype 11 lc rgb "yellow" lw 1
splot '$FCN' u 1:2:3 w pm3d lc "black"

#set origin 0,0
#set size 1,1
unset contour
unset colorbox
unset view
unset pm3d
unset xtics
unset ytics
set style fill transparent solid 0.0 border
set style arrow 1 head filled size screen 0.015,15,45 ls 1 lc "white"
#splot '$STEP' u 2:3:4 w points pt 7 lc rgb "red", '$STEP' u 2:3:4 w l,\
#      '$STEP' u 2:3:4:5 w circles lc "blue"
plot '$STEP' u 2:3 w points pt 7 lc rgb "red", \
     '$STEP' u 2:3:6 w circles lc "gray", '$STEP' u 2:3:7:8 w  vectors as 1
EOF


gnuplot <<EOF

unset key
set view map
set pm3d at b
#set samples 100; set isosamples 100
set cntrparam levels incr 0,1,30
#set cntrlabel onecolor
set palette rgb 33,13,10
set term pngcairo size 800,800
set output 'testsr1line_approx.png'
set title "SR1 Tukey line"

set yrange[-5:5]
set xrange[-5:5]
#set zrange[:20]

#set contour base
#set linetype 11 lc rgb "yellow" lw 1
splot '$FCN' w pm3d lc "black"
set obj 4 rect at -4,0 size 3,1
EOF
