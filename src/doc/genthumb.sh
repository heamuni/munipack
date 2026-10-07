#!/bin/sh
# generates include file for Make providing all images

#set -x

MFILE=image_list.mk

echo -n "image_list = " > $MFILE
for S in "*.png" "*.jp*g" "*.svg"; do
    for A in $S; do
	echo "\\" >> $MFILE
	echo -n " $A " >> $MFILE
    done
done
echo >> $MFILE
