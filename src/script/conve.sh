#!/bin/sh -f

# converts ST-x files to FITSes
# suffixes are:  .st? -> .fits 
# needs konve external

#echo $@

echo "CONVE Version 0.1, Copyright (C) 2000 F.Hroch, Masaryk University, Brno, CZ"

for A in "$@" ; do
    if [ `head -1 $A | awk '{ print index($1,"ST-?")+index($1,"Image") }'` > "2" ]; then
       B=`echo $A | awk '{ gsub(".st[4-8]",".fits"); print; }'`
       if [ "$B" = "$A" ]; then
         B="$A.fits"
       fi
       if [ -x "$B" ]; then
          echo "File $B exist. Creating _$B"
	  B="_$B"
       fi 
      echo $A "->" $B
      konve $A "mask=$B"
    else
       if [ "$A" = "-h" -or "$A" = "--help" ]; then
           echo "Usage: conve [file]"
       fi
       echo "File $A is not ST[4-8]."
       echo "Try -h for help"
    fi
done

