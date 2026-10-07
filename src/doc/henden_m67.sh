#!/bin/bash

# solution by mss:
# http://stackoverflow.com/questions/6489995/can-you-use-heredocuments-to-embed-awk-in-a-bash-script

cat m67ids.txt | awk -f <(sed -e '0,/^#!.*awk/d' $0)
exit $PIPESTATUS

#!/usr/bin/awk -f
BEGIN {

    OUT="m67ids.lst";

print "# BEGIN HDU 0" > OUT;
print "SIMPLE  =                    T / file does conform to FITS standard" > OUT;
print "BITPIX  =                    8 / number of bits per data pixel" > OUT;
print "NAXIS   =                    0 / number of data axes" > OUT;
print "EXTEND  =                    T / FITS dataset may contain extensions" > OUT;
print "COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy" > OUT;
print "COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H" > OUT;
print "END" > OUT;
print "# END HDU 0" > OUT;
print "# BEGIN HDU 1" > OUT;
print "XTENSION= 'BINTABLE'           / binary table extension" > OUT;
print "BITPIX  =                    8 / 8-bit bytes" > OUT;
print "NAXIS   =                    2 / 2-dimensional binary table" > OUT;
print "NAXIS1  =                   56 / 2*8 + 10*4 = width of table in bytes" > OUT;
print "NAXIS2  =                   64 / number of rows in table" > OUT;
print "PCOUNT  =                    0 / size of special data area" > OUT;
print "GCOUNT  =                    1 / one data group (required keyword)" > OUT;
print "TFIELDS =                   12 / number of fields in each row" > OUT;
print "TTYPE1  = 'RAJ2000 '           / label for field   2" > OUT;
print "TFORM1  = 'D       '           / data format of field: 8-byte DOUBLE" > OUT;
print "TUNIT1  = 'deg     '           / physical unit of field" > OUT;
print "TTYPE2  = 'DEJ2000 '           / label for field   3" > OUT;
print "TFORM2  = 'D       '           / data format of field: 8-byte DOUBLE" > OUT;
print "TUNIT2  = 'deg     '           / physical unit of field" > OUT;
print "TTYPE3  = 'U       '" > OUT;
print "TFORM3  = 'E       '" > OUT;
print "TUNIT3  = 'mag     '" > OUT;
print "TTYPE4  = 'B       '" > OUT;
print "TFORM4  = 'E       '" > OUT;
print "TUNIT4  = 'mag     '" > OUT;
print "TTYPE5  = 'V       '" > OUT;
print "TFORM5  = 'E       '" > OUT;
print "TUNIT5  = 'mag     '" > OUT;
print "TTYPE6  = 'R       '" > OUT;
print "TFORM6  = 'E       '" > OUT;
print "TUNIT6  = 'mag     '" > OUT;
print "TTYPE7  = 'I       '" > OUT;
print "TFORM7  = 'E       '" > OUT;
print "TUNIT7  = 'mag     '" > OUT;
print "TTYPE8  = 'Uerr    '" > OUT;
print "TFORM8  = 'E       '" > OUT;
print "TUNIT8  = 'mag     '" > OUT;
print "TTYPE9  = 'Berr    '" > OUT;
print "TFORM9  = 'E       '" > OUT;
print "TUNIT9  = 'mag     '" > OUT;
print "TTYPE10 = 'Verr    '" > OUT;
print "TFORM10 = 'E       '" > OUT;
print "TUNIT10 = 'mag     '" > OUT;
print "TTYPE11 = 'Rerr    '" > OUT;
print "TFORM11 = 'E       '" > OUT;
print "TUNIT11 = 'mag     '" > OUT;
print "TTYPE12 = 'Ierr    '" > OUT;
print "TFORM12 = 'E       '" > OUT;
print "TUNIT12 = 'mag     '" > OUT;
print "EXTNAME = 'M67 by Henden'      / name of this binary table extension" > OUT;
print "COMMENT http://binaries.boulder.swri.edu/fields/m67.html" > OUT;
print "END" > OUT;

}
#
{
    ra=$2;
    dec = $4;
    v = $7;
    b = $8 + v;
    u = $9 + b;
    r = v - $10;
    i = r - $11;
    verr = $12;
    uerr = sqrt($14*$14 + verr*verr);
    berr = sqrt($13*$13 + verr*verr);
    rerr = sqrt($15*$15 + verr*verr);
    ierr = sqrt($16*$16 + verr*verr);

    if( FNR > 1 ) {
	    print ra,dec,u,b,v,r,i,uerr,berr,verr,rerr,ierr > OUT;
	}
}
#
END {
    print "# END HDU 1" > OUT;
    system("munipack fits --restore " OUT);
}
