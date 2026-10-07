#
#  Test for BFSG minimiser for robust regression
#
#  Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
#
#  This file is part of Oakleaf.
#
#  Oakleaf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Oakleaf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.

SRC=../../src/
FLAGS="-Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans"

gfortran $FLAGS $SRC/rfun_huber.f08 $SRC/rfun.f08 \
	 $SRC/qrhouse.f08 $SRC/bfgs.f08 $SRC/huberfit.f08 \
	 rektorys.f08 rline.f08 -o rline

./rline

if [ -z "$1" ]; then
    # no parameter => default to PNG
    TERM="set term pngcairo size 1600,800"
    GOUT="set output 'rline.png'"
else
    # PDF
    TERM="set term epslatex size 2.8in,1.4in monochrome"
    GOUT="set output 'line.tex'"
fi

PAR=/var/tmp/rline.par
XMIN=$(awk -F '=' '{if( $0 ~ /XMIN/) print $2;}' < $PAR)
XMAX=$(awk -F '=' '{if( $0 ~ /XMAX/) print $2;}' < $PAR)
YMIN=$(awk -F '=' '{if( $0 ~ /YMIN/) print $2;}' < $PAR)
YMAX=$(awk -F '=' '{if( $0 ~ /YMAX/) print $2;}' < $PAR)

FCN=/var/tmp/rline_fcn
STEP=/var/tmp/rline_step

gnuplot <<EOF
#cat <<EOF

unset key
set view map
set pm3d at b
unset surface
set palette model RGB
set palette defined (0 "#66cd00", 1 "white")
set linetype 1 lc rgb "dark-violet" lw 3
$TERM
$GOUT

set xrange[$XMIN:$XMAX]
set yrange[$YMIN:$YMAX]
#set zrange[0:1e5]

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

splot [][][:] '$FCN' w pm3d lc "black"

unset xtics
unset ytics
unset colorbox
plot '/tmp/contour.dat' w l lc rgb "black" lw 1

unset contour
unset colorbox
unset view
unset pm3d
set style arrow 1 head filled size screen 0.015,15,45 ls 1 lc "white"
set style arrow 2 head filled size screen 0.015,15,45 ls 1 lc "grey"
set style arrow 3 head filled size screen 0.015,15,45 ls 1 lc "gray"
set style arrow 4 head filled size screen 0.015,15,45 ls 1 lc "blue"
set style arrow 5 head filled size screen 0.015,15,45 ls 1 lc "black"
set style arrow 6 head filled size screen 0.015,15,45 ls 1 lc "orange"
set style arrow 9 head filled size screen 0.015,15,45 ls 1 lc "red" dt 2

set style textbox opaque fc "gray94" noborder

CL=22
plot '$STEP' u 3:4 w points pt 7 lc rgb "red", \
     '$STEP' u 3:4:8:9:(column(CL) <= 10 ? 9 : column(CL)-10) w vectors as var
EOF

gnuplot <<EOF
$TERM

unset key
#set view map
set pm3d at b
#unset surface
set palette model RGB
set palette defined (0 "#66cd00", 1 "white")
set linetype 1 lc rgb "dark-violet" lw 3
#$TERM
#$GOUT
set output 'rline3d.png'

set xrange[$XMIN:$XMAX]
set yrange[$YMIN:$YMAX]
set zrange[0:100]

set contour base
set cntrparam levels incr 0,1,100
set table "/var/tmp/contour.dat"
#splot '$FCN' w l lc rgb "black"
unset table

#splot [][][:] '$FCN' w pm3d lc "black"

unset xtics
unset ytics
unset colorbox
#plot '/var/tmp/contour.dat' w l lc rgb "black" lw 1

set hidden3d front

set style arrow 1 head filled size screen 0.015,15,45 ls 1 lc "grey"
set style arrow 2 head filled size screen 0.015,15,45 ls 1 lc "black"
set style arrow 9 head filled size screen 0.015,15,45 ls 1 lc "red" dt 2
set zrange[0:]
CL=22
splot '/var/tmp/rline_step' u 3:4:2:8:9:7:(column(CL) <= 10 ? 9 : column(CL)-10)   w vectors as var, '$FCN' w pm3d lc "black"
EOF
