
./plotrfun

gnuplot <<EOF
set term png
set output 'hubers.png'
plot 'hubers', 'hubers' u 1:3, 'hubers' u 1:4
EOF

gnuplot <<EOF
set term png
set output 'tukeys.png'
plot 'tukeys', 'tukeys' u 1:3, 'tukeys' u 1:4
EOF
