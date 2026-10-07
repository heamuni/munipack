
OUT=testmean.out
testmean > $OUT
XERR=$?

PASS=$(awk 'BEGIN { n=0; pass=0;} {if( $0 ~/^TEST/ ) { n++; pass += $4 == "T" ? 1 : 0;}} END {print !(n == pass);}' < testmean.out)

#PASS=$(awk -f testparse.awk < $OUT)

exit $PASS && $XERR
