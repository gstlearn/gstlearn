#!/bin/sh

WHERE=$1
RADIX=$2

REF=$WHERE/$RADIX.ref
LOC=$WHERE/$RADIX.out

if [ -f $LOC ]; then
    mv $LOC $REF
fi
