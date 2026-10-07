
set -e

if [ -z "$1" ]; then
    # no parameter => default to PNG
    TERM="set term pngcairo size 1600,800"
    HTERM="set term pngcairo size 1600,2400"
    GOUT="set output 'testsr1mean.png'"
    FOUT="set output 'testsr1rho.png';"
    HOUT="set output 'testsr1hess.png';"
    ROUT="set output 'testsr1jac.png';"
    LOUT="set output 'testsr1lam.png';"
#    POUT="set output 'testsr1prof.png';"
    FMT='"%0.0f"'
    LSZ=""
else
    # PDF
    TERM="set term epslatex size 2.8in,1.4in monochrome"
    TERM="set term epslatex size 5in,2.5in color"
    HTERM="set term epslatex size 4in,10in color"
    GOUT="set output 'testsr1mean.tex'"
    FOUT="set output 'testsr1rho.tex'"
    HOUT="set output 'testsr1hess_xx.tex';"
    ROUT="set output 'testsr1jac.tex';"
    LOUT="set output 'testsr1lam.tex';"
#    POUT="set output 'testsr1prof.tex';"
    FMT='"{\\small \$%0.0f\$}"'
    LSZ="\\\\footnotesize"
fi


SRC=../../src/
FLAGS="-Wall -cpp -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
#gfortran  $SRC/cholesky.f08 $SRC/trust_step.f08 $SRC/sr1.f08 $SRC/rfun.f08 ../noise.f08 -std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow -o testsr1mean  testsr1mean.f08
#gfortran $FLAGS $SRC/cholesky.f08 $SRC/fmm/zeroin.f95 $SRC/fmm/fmin.f95 $SRC/fmm/fmm.f08 $SRC/qsort.f08 $SRC/qmean.f08 $SRC/iscale.f08 $SRC/q2eig.f08 $SRC/trust_step.f08 $SRC/sr1.f08 $SRC/qrhouse.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/zeroin.f08 $SRC/glomin.f08 trgscale.f08 $SRC/sr1meanfun.f08 $SRC/sr1mean.f08 ../noise.f08 kscale.f08 biscale.f08 -o testsr1mean testsr1mean.f08

#gfortran $FLAGS $SRC/rfun.f08 $SRC/rfun_hampel.f08 $SRC/rlike.f08 $SRC/q2eig.f08 $SRC/qtrstep.f08 $SRC/fmin.f08 $SRC/qsort.f08 $SRC/qmean.f08 $SRC/qmed.f08 $SRC/sqpfmean.f08 $SRC/rmean.f08 ../noise.f08 -o testsr1mean testsr1mean.f08
gfortran $FLAGS $SRC/rfun_hampel.f08 $SRC/rfun_tukey.f08 $SRC/rfun_huber.f08 $SRC/rnfun_square.f08 $SRC/rfun_real64.f08 $SRC/rlike.f08 $SRC/q2eig.f08 $SRC/qtrstep.f08 $SRC/fmin.f08 $SRC/qsort_real64.f08 $SRC/qmed_real64.f08 $SRC/sqpfmean.f08 $SRC/rmean.f08 ../noise.f08 -o testsr1mean testsr1mean.f08

FCN=/tmp/tsr1mean
HESS=/tmp/tsr1mean_hess
JAC=/tmp/tsr1mean_jac
#RESULT=/tmp/tsr1mean_result
#STDDEV=/tmp/tsr1mean_stddev
STEP=/tmp/tsr1mean_step
#REG=/tmp/tsr1mean_qlim
RIDGE=/tmp/tsr1mean_ridge
RHO=/tmp/tsr1mean_rho
FCNTXT=/tmp/tsr1mean_txt

./testsr1mean

echo -n "Generating graphs .. "

PAR=/tmp/testsr1mean.par
TMIN=$(awk -F '=' '{if( $0 ~ /^TMIN/) print $2;}' $PAR)
TMAX=$(awk -F '=' '{if( $0 ~ /^TMAX/) print $2;}' $PAR)
SMIN=$(awk -F '=' '{if( $0 ~ /^SMIN/) print $2;}' $PAR)
SMAX=$(awk -F '=' '{if( $0 ~ /^SMAX/) print $2;}' $PAR)

cat > $FCNTXT <<EOF
0.5 1.5  {$LSZ(a)}
-2.5 4.5 {$LSZ(b)}
4.5 0.8 {$LSZ(c)}
-0.5 0.3 {$LSZ(d)}
-2.5 0.3 {$LSZ(e)}
4 1.9 {$LSZ(f)}
EOF


gnuplot <<EOF
#cat <<EOF

unset key
set view map
set pm3d at b
unset surface
#set samples 100; set isosamples 100
#set cntrparam levels incr 0,1,10
#set cntrlabel onecolor
#set palette rgb 33,13,10
set palette model RGB
set palette defined (0 "#66cd00", 1 "white")
set linetype 1 lc rgb "dark-violet" lw 3
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
#set cbrange[0:10]

set contour base
set cntrparam levels incr 0,1,100
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
#set style fill transparent solid 0.0 border
set style arrow 1 head filled size screen 0.015,15,45 ls 1 lc "white"
set style arrow 2 head filled size screen 0.015,15,45 ls 1 lc "white"
set style arrow 3 head filled size screen 0.015,15,45 ls 1 lc "gray"
set style arrow 4 head filled size screen 0.015,15,45 ls 1 lc "blue"
set style arrow 5 head filled size screen 0.015,15,45 ls 1 lc "black"
set style arrow 6 head filled size screen 0.015,15,45 ls 1 lc "orange"
set style arrow 9 head filled size screen 0.015,15,45 ls 1 lc "red" dt 2

#splot '$STEP' u 2:3:4 w points pt 7 lc rgb "red", '$STEP' u 2:3:4 w l,\
#      '$STEP' u 2:3:4:5 w circles lc "blue"

#set arrow from -2.5,4.5 to -2.6,5 as 1
set style textbox opaque fc "gray94" noborder
#set style textbox transparent noborder
plot '$STEP' u 2:3 w points pt 7 lc rgb "red", \
     '$STEP' u 2:3:5 w circles lc "yellow" lw 1,\
     '$STEP' u 2:3:6:7:(column(14)==-1 ? 9 : column(14)+1) w vectors as var,\
     '$RIDGE' u 1:2:3 lw 1 lc rgb "gray40" w filledcu fs transparent pattern 10,\
     '$RIDGE' u 1:2 w l lw 1 lc rgb "gray40"
EOF
#     '$RIDGE' u 1:2 w points pt 13 lc "orange",\
#     '$REG' w l, \
#     '$FCNTXT' w labels boxed

gnuplot <<EOF

unset key
set view map
set pm3d at b
unset surface
#set samples 100; set isosamples 100
#set cntrparam levels incr 0,1,30
#set cntrlabel onecolor
#set palette rgb 33,13,10
set palette model RGB
set palette defined (0 "#66cd00", 1 "white")
$TERM
$FOUT
#set term pngcairo size 800,400
#set output 'testsr1frac.png'
#set term epslatex size 4in,3in solid colour
#set output 'testsr1mean.tex'
#set title "Test SR1 Tukey mean"

set contour base
#set cntrparam levels incr -0.4,0.2,1.1
set cntrparam levels incr 0,1,10
set table "/tmp/contour_rho.dat"
splot [:][:][:] '$RHO' u 1:2:3 w l lc rgb "black"
unset table

set yrange[$SMIN:$SMAX]
set xrange[$TMIN:$TMAX]
#set zrange[-0.6:1.1]
#set cbrange[-0.6:1.1]
#set zrange[0:6]

set multiplot

set lmargin at screen 0.1
set rmargin at screen 0.85
set tmargin at screen 0.1
set bmargin at screen 0.9

#set contour base
#set origin 0,0
#set size 1,1
#set linetype 11 lc rgb "yellow" lw 1
#set cbtics ("\$-\\\\sfrac{1}{2}\$" -0.5, "\$0\$" 0, "\$\\\\sfrac{1}{2}\$" 0.5, \
#           "\$1\$" 1)
#set xtics 2 format $FMT
#set ytics 2 format $FMT
splot '$RHO' u 1:2:3 w pm3d lc "black"

unset xtics
unset ytics
unset colorbox
plot '/tmp/contour_rho.dat' w l lc rgb "black" lw 1

EOF

#HOUT11=$(echo $HOUT | sed 's/xx/11/')
#HOUT12=$(echo $HOUT | sed 's/xx/12/')
#HOUT22=$(echo $HOUT | sed 's/xx/22/')

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
$HTERM
$HOUT

set yrange[$SMIN:$SMAX]
set xrange[$TMIN:$TMAX]
#set zrange[-3:3]
#set cbrange[-0.6:1.1]
#set zrange[0:6]

set multiplot

#set lmargin at screen 0.1
#set rmargin at screen 0.85
#set tmargin at screen 0.1
#set bmargin at screen 0.9

#set contour base
set size 1,0.333
set origin 0,0
#set linetype 11 lc rgb "yellow" lw 1
#set cbtics ("\$-\\\\sfrac{1}{2}\$" -0.5, "\$0\$" 0, "\$\\\\sfrac{1}{2}\$" 0.5, \
#           "\$1\$" 1)
#set xtics 2 format $FMT
#set ytics 2 format $FMT
splot '$HESS' u 1:2:3 w pm3d lc "black" t "h11"

#unset xtics
#unset ytics
#unset colorbox
#plot '/tmp/contour.dat' w l lc rgb "black"


#$HOUT22

set origin 0,0.333
#set size 1,0.333

splot '$HESS' u 1:2:5 w pm3d lc "black" t "h22"
#plot '/tmp/contour.dat' w l lc rgb "black"

set origin 0,0.666
#set size 1,0.333

#$HOUT12

splot '$HESS' u 1:2:4 w pm3d lc "black" t "h12"
#plot '/tmp/contour.dat' w l lc rgb "black"
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
$HTERM
$LOUT

set yrange[$SMIN:$SMAX]
set xrange[$TMIN:$TMAX]
#set zrange[-3:3]
#set cbrange[-0.6:1.1]
#set zrange[-2:2]

set multiplot

set size 1,0.5
set origin 0,0
splot '$HESS' u 1:2:(\$7*\$8) w pm3d lc "black" t "pos.def."

set origin 0,0.5
#set zrange[-5:5]

splot '${HESS}' u 1:2:(57.3*atan2(\$9,\$11)-90) w pm3d

#plot '${HESS}_Q' u 1:2:( \$9*0.15):(\$11*0.15) w vectors lc "black" lw 1,\
#     '${HESS}_Q' u 1:2:(\$10*0.15):(\$12*0.15) w vectors lc "blue"  lw 1
EOF



gnuplot <<EOF

unset key
set view map
set pm3d at b
unset surface
set contour base
#set samples 100; set isosamples 100
set cntrparam levels incr -2,1,2
#set cntrlabel onecolor
#set palette rgb 33,13,10
set palette model RGB
set palette defined (0 "#66cd00", 1 "white")
$HTERM
$ROUT

#set yrange[0:5]
#set xrange[-5:5]
set yrange[$SMIN:$SMAX]
set xrange[$TMIN:$TMAX]
#set zrange[-3:3]
#set cbrange[-0.6:1.1]
#set zrange[-5:5]

set multiplot

set size 1,0.5
set origin 0,0
splot '$JAC' u 1:2:3 w pm3d lc "black" t "ft"

set origin 0,0.5
#set zrange[-5:5]

splot '$JAC' u 1:2:4 w pm3d lc "black" t "fs"
EOF




#gnuplot <<EOF
#
#unset key
#set view map
#set pm3d at b
##set samples 100; set isosamples 100
#set cntrparam levels incr 0,1,30
##set cntrlabel onecolor
#set palette rgb 33,13,10
#set term pngcairo size 800,400
#set output 'testsr1stddev.png'
#set title "Test SR1 Tukey stddev"
#
#set yrange[0:5]
#set xrange[-5:5]
#set zrange[:3]
#set multiplot
#
#
#set lmargin at screen 0.1
#set rmargin at screen 0.85
#set tmargin at screen 0.1
#set bmargin at screen 0.9
#
#set contour base
##set origin 0,0
##set size 1,1
##set linetype 11 lc rgb "yellow" lw 1
#splot '$STDDEV' w pm3d lc "black"
#
##set origin 0,0
##set size 1,1
#unset contour
#unset colorbox
#unset view
#unset pm3d
#unset xtics
#unset ytics
#set style fill transparent solid 0.0 border
#set style arrow 1 head filled size screen 0.015,15,45 ls 1 lc "white"
##splot '$STEP' u 2:3:4 w points pt 7 lc rgb "red", '$STEP' u 2:3:4 w l,\
##      '$STEP' u 2:3:4:5 w circles lc "blue"
#plot '$RESULT' w points pt 7 lc rgb "black"
#EOF



#gnuplot <<EOF
#
#unset key
#set view map
#set pm3d at b
##set samples 100; set isosamples 100
#set cntrparam levels incr 0,1,30
##set cntrlabel onecolor
#set palette rgb 33,13,10
#set term pngcairo size 800,800
#set output 'testsr1mean_approx.png'
#set title "SR1 Tukey mean"
#
#set yrange[0:3]
#set xrange[-5:5]
##set zrange[:2]
#
##set contour base
##set linetype 11 lc rgb "yellow" lw 1
#splot '$FCN' w pm3d lc "black"
#set obj 4 rect at -4,0 size 3,1
#EOF

#if [ "$1" ]; then
#   for A in testsr1mean testsr1frac testsr1hess; do
#      E=${A}.eps
#      if [ -f $E ]; then
#        epstopdf $E
#        rm -f $E
#      fi
#      mv -i ${A}.tex ${A}.pdf ~/rontrib/
#    done
#fi

echo "done"
