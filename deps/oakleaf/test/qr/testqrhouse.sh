
#export AM_TEST_SUITE=1

#set -xe

#SRC=../../src/
#FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"
#gfortran  $SRC/qrhouse.f08 -o testqrhouse testqrhouse.f08
#gfortran  $SRC/qrdebug.f08 -o testqrhouse testqrdebug.f08


testqrhouse > testqrhouse.out
if [ $? -eq 0 ]; then
    RESULT=$(awk -F'= ' '{ if( $1 ~/testqrhouse/) print $2}' testqrhouse.out)
    if [ $? -eq 0 -a "$RESULT" = "PASS" ]; then
	exit 0
    fi
fi

exit 1
