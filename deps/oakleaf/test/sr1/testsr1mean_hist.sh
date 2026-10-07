
TERM="set term pngcairo size 1600,800"
GOUT="set output 'testsr1mean_hist_land.png'"


set -xe

SRC=../../src/
FLAGS="-Wall -O2 -g -p -cpp -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
gfortran $FLAGS $SRC/qsort_real64.f08 $SRC/qmed_real64.f08 $SRC/qmean_real64.f08 $SRC/q2eig.f08 $SRC/qrhouse.f08 $SRC/fmin.f08 $SRC/rfun_hampel_real64.f08 $SRC/rfun_tukey_real64.f08 $SRC/rfun_huber_real64.f08 $SRC/rnfun_square.f08 $SRC/rfun_real64.f08 $SRC/qtrstep.f08 $SRC/rlike.f08 $SRC/sqpfmean.f08 $SRC/rmean.f08  ../noise.f08 ../toolbox.f08 -o testsr1mean_hist testsr1mean_hist.f08

#gfortran $FLAGS $SRC/fmm/zeroin.f95 $SRC/fmm/fmin.f95 $SRC/fmm/fmm.f08 $SRC/rfun.f08 $SRC/rfun_hampel.f08 $SRC/rlike.f08 $SRC/q2eig.f08 $SRC/qtrstep.f08 $SRC/zeroin.f08 $SRC/qsort.f08 $SRC/qmean.f08 $SRC/sqpfmean.f08 $SRC/rmean3.f08 $SRC/sr1meanfun.f08 $SRC/sr1mean.f08 ../noise.f08 kscale.f08 -o testsr1mean testsr1mean.f08


./testsr1mean_hist

FCN=/tmp/testsr1mean_hist_rland
NDAT=$(cat /tmp/testsr1mean_hist.par)

echo NDAT= $NDAT

PAR=/tmp/testsr1mean_hist_rland.par
TMIN=$(awk -F '=' '{if( $0 ~ /^TMIN/) print $2;}' $PAR)
TMAX=$(awk -F '=' '{if( $0 ~ /^TMAX/) print $2;}' $PAR)
SMIN=$(awk -F '=' '{if( $0 ~ /^SMIN/) print $2;}' $PAR)
SMAX=$(awk -F '=' '{if( $0 ~ /^SMAX/) print $2;}' $PAR)


gnuplot <<EOF
set output 'testsr1mean_hist_t.png'
set term pngcairo size 800,400
#set term epslatex size 4in,1.15in solid colour
#set output 'testsr1mean_hist_t.tex'
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
#set xtics 1.9,0.1,2.2
#set xtics -3,1 format "\$%.0f\$"
#set ytics format ""
#set xrange[1.85:2.25]
#set xrange[-3:3]
sig = 1/sqrt($NDAT)
plot '/tmp/testsr1mean_t' ls 1 with boxes t "rob", \
     '/tmp/testsr1mean_density_t' w l t "N(0,sig)" ls 2
#plot '/tmp/testsr1mean_t' ls 1 with boxes t "rob", exp(-(x/sig)**2/2)/sqrt(2*pi)/sig t "N(0,sig)" ls 2
EOF


gnuplot <<EOF
set output 'testsr1mean_hist_s.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
#set xtics 1.9,0.1,2.2
#set ytics format ""
#set xrange[1.85:2.25]
plot '/tmp/testsr1mean_s' ls 1 with boxes t "rob", \
     '/tmp/testsr1mean_density_s' w l t "N(s,sig)" ls 2
#plot '/tmp/testsr1mean_s' ls 1 with boxes t "rob"
EOF

gnuplot <<EOF
set output 'testsr1mean_hist_hess11.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
plot '/tmp/testsr1mean_hess11' ls 1 with boxes t "rob"
EOF

gnuplot <<EOF
set output 'testsr1mean_hist_hess12.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
plot '/tmp/testsr1mean_hess12' ls 1 with boxes t "rob"
EOF

gnuplot <<EOF
set output 'testsr1mean_hist_hess22.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
plot '/tmp/testsr1mean_hess22' ls 1 with boxes t "rob"
EOF

gnuplot <<EOF
set output 'testsr1mean_hist_avgR.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
plot '/tmp/testsr1mean_avgR' ls 1 with boxes t "rob"
EOF

gnuplot <<EOF
set output 'testsr1mean_hist_sig.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
#set xtics 1.9,0.1,2.2
#set ytics format ""
#set xrange[1.85:2.25]
#sig = 1/sqrt(2*$NDAT)
logchi(x,k) = (k-1)*log(x) - x**2/2 - (k/2.0-1)*log(2.0) - log(gamma(k/2.0))
#plot '/tmp/testsr1mean_sig' ls 1 with boxes t "rob", \
#    exp(logchi2(x,$NDAT)) w l t "chi2(x,$NDAT)"
#     '/tmp/testsr1mean_density_sig' w l t "N(s,sig)" ls 2,
#     exp(logchi(x,$NDAT))
plot '/tmp/testsr1mean_sig' ls 1 with boxes t "rob"
#     '/tmp/testsr1mean_density_sig' w l t "N(s,sig)" ls 2
EOF

gnuplot <<EOF
set output 'testsr1mean_hist_chi2.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
#set xtics 1.9,0.1,2.2
#set ytics format ""
#set xrange[1.85:2.25]
#sig = 1/sqrt(2*$NDAT)
logchi2(x,k) = (k/2.0-1)*log(x) - x/2 - (k/2.0)*log(2) - log(gamma(k/2.0))
#plot '/tmp/testsr1mean_chi2' ls 1 with boxes t "chi2",
#    exp(logchi2(x,$NDAT)) w l  t "chi2(x,$NDAT)"
#     '/tmp/testsr1mean_density_chi2' w l t "N(s,sig)" ls 2,
#plot '/tmp/testsr1mean_chi2' ls 1 with boxes t "chi2",
#     exp(logchi2(x,${NDAT}-1))
plot '/tmp/testsr1mean_chi2' ls 1 with boxes t "chi2", \
     '/tmp/testsr1mean_density_chi2' w l ls 2
EOF

gnuplot <<EOF
set output 'testsr1mean_std_t.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
plot '/tmp/testsr1mean_std_t' ls 1 with boxes t "rob", \
     exp(-x**2/2)/sqrt(2*pi) w l t "N(0,1)" ls 2
EOF


gnuplot <<EOF
set output 'testsr1mean_std_s.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
plot '/tmp/testsr1mean_std_s' ls 1 with boxes t "rob", \
     exp(-x**2/2)/sqrt(2*pi) w l t "N(0,1)" ls 2
EOF

gnuplot <<EOF
set output 'testsr1mean_std_chi2.png'
set term pngcairo size 800,400
set style line 1 lt rgb "salmon3" lw 2
set style line 2 lt rgb "green" lw 3
plot '/tmp/testsr1mean_std_chi2' ls 1 with boxes t "rob", \
     exp(-x**2/2)/sqrt(2*pi) w l t "N(0,1)" ls 2
EOF



gnuplot <<EOF

unset key
set view map
set pm3d at b
unset surface
#set samples 100; set isosamples 100
set cntrparam levels incr 0,1,30
#set cntrlabel onecolor
#set palette rgb 33,13,10
set palette model RGB
set palette defined (0 "#66cd00", 1 "white")
$TERM
$GOUT
#set term pngcairo size 800,400
#set output 'testsr1mean.png'
#set term epslatex size 4in,3in solid colour
#set output 'testsr1mean.tex'
#set title "Test SR1 Tukey mean"

set yrange[$SMIN:$SMAX]
set xrange[$TMIN:$TMAX]
#set zrange[0:10]
#set cbrange[0:5]

set contour base
set cntrparam levels incr 1,1,200
set table "/tmp/contour.dat"
splot '$FCN' w l lc rgb "black"
unset table


set multiplot

set lmargin at screen 0.1
set rmargin at screen 0.85
set tmargin at screen 0.1
set bmargin at screen 0.9

#set contour base
#set origin 0,0
#set size 1,1
#set linetype 11 lc rgb "yellow" lw 1
#set cbtics 0,1
#set xtics 2 format $FMT
#set ytics 2 format $FMT
splot [][][:] '$FCN' w pm3d lc "black"

unset xtics
unset ytics
unset colorbox
plot '/tmp/contour.dat' w l lc rgb "black" lw 1

#set origin 0,0
#set size 1,1
unset contour
unset colorbox
unset view
unset pm3d

set style circle radius graph 0.002
plot '/tmp/testsr1mean_hist_rank' w circles fs transparent solid 0.15 noborder lc rgb "orange", '/tmp/testsr1mean_hist' u 1:3 w circles fs transparent solid 0.15 noborder lc rgb "black"

EOF
