
RFUN="hampel"

set -xe

SRC=../../src/
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"

gfortran $FLAGS $SRC/cholesky.f08 $SRC/fmm/zeroin.f95 $SRC/fmm/fmin.f95 $SRC/fmm/fmm.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/trust_step.f08 $SRC/sr1.f08 $SRC/sr1meanfun.f08 $SRC/sr1mean.f08 \
	 ../noise.f08 \
	 -o sr1_region sr1_region.f08

./sr1_region $RFUN

FCN=/tmp/sr1region.dat
REG=/tmp/sr1region_region.dat

gnuplot <<EOF

unset key
set view map
set pm3d at b
unset surface
#set samples 100; set isosamples 100
set cntrparam levels incr 0,1,30
#set cntrlabel onecolor
set palette rgb 33,13,10
set term pngcairo size 800,400
set output 'sr1region.png'
#set term epslatex size 4in,3in solid colour
#set output 'testsr1mean.tex'
set title "Convergence region for $RFUN"

set contour base
set cntrparam levels incr 0,0.5,6
set table "/tmp/contour.dat"
splot [:][:][:6] '$FCN' w l lc rgb "black"
unset table

set contour base
set cntrparam levels incr 0,0.5,6
set table "/tmp/contour_region.dat"
splot [:][:][:] '$FCN' w l lc rgb "black"
unset table

set yrange[0:5]
set xrange[-5:5]
#set zrange[0:10]

set multiplot

set lmargin at screen 0.1
set rmargin at screen 0.85
set tmargin at screen 0.1
set bmargin at screen 0.9

#set contour base
#set origin 0,0
#set size 1,1
#set linetype 11 lc rgb "yellow" lw 1
set cbtics 0,1
splot [][][:3] '$FCN' w pm3d lc "black"

unset xtics
unset ytics
unset colorbox

plot '/tmp/contour_region.dat' w l lc rgb "black"

#     '/tmp/contour_reg.dat' w l lc rgb "black"
#plot '/tmp/contour_reg.dat' w l lc rgb "black"


#set origin 0,0
#set size 1,1
unset contour
unset colorbox
unset view
unset pm3d
set style fill transparent solid 0.0 border
set style arrow 1 head filled size screen 0.015,15,45 ls 1 lc "white"
plot '$REG' u 1:2:3 with filledcu fs transparent pattern 7, \
     '$REG' u 1:2 w l, '$REG' u 1:3 w l
#plot '$REG' u 1:2:3 lw 1 lc rgb "gray40" w filledcu fs transparent pattern 3,\
#     '$REG' u 1:2 w l lw 1 lc rgb "gray40"

EOF
