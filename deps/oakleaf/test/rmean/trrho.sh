
set -xe

SRC=../../src/
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
gfortran $FLAGS $SRC/qsort.f95 $SRC/qmean.f08 $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/sr1meanfun.f08 ../noise.f08 trrho.f08 -o trrho

./trrho

TRRHO_PSHAMPEL='/tmp/trrho_pshampel'
TRRHO_TUKEY='/tmp/trrho_tukey'

gnuplot <<EOF
#set term pngcairo size 800,640
#set output 'trrho.png'
set term epslatex size 2in,2in solid #monochrome
set output 'trrho_tukey.tex'

set xrange[-10.1:10.1]
set yrange[-0.1:1.6]
set style line 1 lt rgb "blue" lw 4
unset key
set xtics -10.1,5,10.1 format "% 0.0f"
set ytics 0,1,2 format "% 0.0f"

#plot '$TRRHO_PSHAMPEL' w l, '$TRRHO_TUKEY' w l
plot '$TRRHO_TUKEY' w l ls 1
EOF


for A in trrho_tukey.tex; do
    E=${A%.tex}.eps
    if [ -f $E ]; then
	epstopdf $E
	rm -f $E
    fi
done
