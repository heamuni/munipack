
# compile, run
gfortran -I../lib -Wall -g -p -fcheck=all -fimplicit-none -fbacktrace lcapprox.f08 -L../lib -lminpacks -lminpack && ./a.out

# create graphs
cd /tmp/
gnuplot <<EOF
set term svg dynamic
set output 'lcapprox.svg'
set ylabel "rel. flux"
set ytics 0.3,0.1,1.0 format "%0.1f"
set yrange[0.2:1.1]
set key left bottom
set multiplot
set xtics format ""
set size 1,0.45
set origin 0,0.55
set tmargin 0.5
set bmargin 0
plot 'dcep' t "δ Cep theoretical data" lc rgb '#8AB8E6', 'dcep' u 1:3 t "the approximation" w l lc rgb '#000000' lw 2
set size 1,0.55
set origin 0,0
set tmargin 0
set bmargin 3
set xtics format "%0.1f"
set xlabel "phase"
plot 'algol' t "Algol observed data" lc rgb '#8AB8E6', 'algol' u 1:3 t "the approximation" w l lc rgb '#000000' lw 2
EOF
