
# run flat with debug on (recompile required):
#     logical, parameter :: debug = .true.
#

gnuplot <<EOF
bin(x,width)=width*floor(x/width)
set style line 1 lt rgb "#8AB8E6" lw 1
set style line 2 lt rgb "black" lw 1
set size 0.618,1

set term svg dynamic
set output 'flatdebug_zero.svg'
set xrange[-10:10]
unset key
set xlabel "residual"
set ylabel "number of values"
binwidth=0.1
plot 'flatdebug_zero.dat' using (bin(\$2,binwidth)):(1.0) smooth freq with boxes ls 1, 18.5e3*exp(-x**2/2) ls 2

set output 'flatdebug_3.svg'
set autoscale
binwidth=0.05
unset key
set xlabel "residual from N(0,1)"
set ylabel "number of values"
plot 'flatdebug_3.dat' using (bin(\$2,binwidth)):(1.0) smooth freq with boxes ls 1,\
     21e3*exp(-x**2/2) ls 2
EOF
