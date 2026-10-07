
#  This script does visualisation of growth-curve reports.
#
#  Copyright © 2016 F.Hroch (hroch@physics.muni.cz)
#
#  This file is part of Munipack.
#
#  Munipack is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  Munipack is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with Munipack.  If not, see <http://www.gnu.org/licenses/>.



#
#  Usage:
#
#  * Run this script in a directory where reports are generated.
#    Use the calling sequence:
#       $ bash grow_report.sh
#
#  * The reports are created by passing --report option to aphot:
#    ($ munipack aphot --report [fits-files])
#
#  * Gnuplot and python with matplotlib are required.
#
#  * Default graphical format is PNG. It can be changed with FRM
#    variable (the format must be supported by gnuplot).
#
#  * The results will stored in grow_report/ subdirectory.
#
#

# uncomment for tracking
# set -x

# graphical format
FRM=png

# default subdirectory
SDIR=grow_report/

# default directory
DIR=$PWD/$SDIR

#
# ---- code for processing of a single report file -----
#

report()
{
    F=$1
    B=${F%.fits}
    O=${DIR}/${B}_curve.png
    Q=${DIR}/${B}.png
    CURVE=$(mktemp grow_curve.XXXXXXXXX)
    DATA=$(mktemp grow_data.XXXXXXXXX)
    DATA1=$(mktemp grow_data_1.XXXXXXXXX)
    DATA2=$(mktemp grow_data_2.XXXXXXXXX)

    python grow_report.py display $F $Q

    munipack fits -lt $F\[GROWCURVE\] > $CURVE
    python grow_report.py curve $F > $DATA
    grep "1$" < $DATA > $DATA1
    grep -v "1$" < $DATA > $DATA2

    gnuplot <<EOF
set term '$FRM'
set output '$O'
set style line 99 lt 2 lc rgb "gold" lw 2
unset key
set lmargin 5
set rmargin 2
set multiplot title "Growth-curve ($A)"
set origin 0,0.3
set size 1,0.7
set key right bottom
set tmargin 2
set bmargin 0
unset xtics
set xrange[0:21]
set yrange[0:1.618]
set ytics 0,0.5,1.5 format "%0.1f"
set grid
plot '$CURVE' t "growth-curve" with lines ls 99, '$DATA2' t "all stars" with points pt 6 lc rgb "#87CEFA", '$DATA1' t "designed by" with points pt 13 lc rgb "blue"
set origin 0,0
set size 1,0.3
unset key
set bmargin 3
set tmargin 0
set xtics nomirror
set xrange[0:21]
set yrange[-0.21:0.21]
set ytics -0.1,0.1,0.1 format "%0.1f"
set xlabel "distance from a star center in pixels"
set grid
plot '$DATA2' u 1:3 with points pt 6 lc rgb "#87CEFA", '$DATA1' u 1:3 with points pt 13 lc rgb "blue"
EOF
    

    rm -f "$CURVE" "$DATA" "$DATA1" "$DATA2"
}

report_svg()
{
    F=$1
    B=${F%.fits}
    O="${DIR}/${B}_curve.png"
    TITLE=$3
    CURVE=$(mktemp grow_curve.XXXXXXXXX)
    DATA=$(mktemp grow_data.XXXXXXXXX)

    if [ "$TITLE" ]; then
	TTT=" title \"Growth-Curve\""
	LABEL="set xlabel \"distance from a star center in pixels\""
	FILENAME="grow_curve_anotated.svg"
    else
	FILENAME="grow_curve.svg"
    fi
    
    munipack fits -lt $F\[GROWCURVE\] > $CURVE
    python grow_report.py report $F > $DATA
    grep "1$" < $DATA > $DATA1
    grep -v "1$" < $DATA > $DATA2

    gnuplot <<EOF
set term 'svg' size 500,400 dynamic fsize 11 rounded solid
set output '$FILENAME' 
set style line 99 lt 2 lc rgb "gold" lw 2
set lmargin 5
set rmargin 2
set multiplot $TTT
set origin 0,0.3
set size 1,0.7
set key right bottom
set tmargin 2
set bmargin 0
unset xtics
set xrange[0:21]
set yrange[0:1.618]
set ytics 0,0.5,1.5 format "%0.1f"
plot '$CURVE' t "growth-curve" with lines ls 99, '$DATA2' t "all stars" with dots ls 3 lc rgb "#B0C4DE", '$DATA1' t "designed by" with dots ls 13 lc rgb "blue"
set origin 0,0
set size 1,0.3
unset key
set bmargin 3
set tmargin 0
set xtics nomirror
set xrange[0:21]
set yrange[-0.21:0.21]
set ytics -0.1,0.1,0.1 format "%0.1f"
$LABEL
set grid
plot '$DATA2' u 1:3 with dots ls 3 lc rgb "#87CEFA", '$DATA1' u 1:3 with dots ls 3 lc rgb "blue"
EOF
    

    rm -f "$CURVE" "$DATA" "$DATA1" "$DATA2"
}


# ----  processing starts here  -----

mkdir -p $DIR

# if filenames contains spaces, firstly replace its with underscores:
#find . -maxdepth 1 -name '*.fits' -print0 | xargs -0 rename 's/ /_/g'

if [ "$1" ]; then
    
   for A in $*; do
      report $A
   done

else
    
   for A in $(ls *.fits); do
      report $A
   done
fi

# animate gif (?)

