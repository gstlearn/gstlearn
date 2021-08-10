#!/bin/csh

unalias ls
unalias basename
unalias diff

set TMP   = "/tmp/__gstlearn__.tmp"
set WHERE = $argv[1]
set LISTE = `ls $WHERE/*.out`

foreach file ($LISTE)

    set RADIX = `basename $file .out`
    set REF = $WHERE/$RADIX.ref
    set LOC = $WHERE/$RADIX.out

    if ( -f $REF) then
	if ( -f $LOC) then
	    cp $LOC $TMP
	    sed "s/ -0.000/  0.000/g"  < $TMP > $LOC
	    diff -qbi $REF $LOC
	endif
    else
	echo "Reference file ("$REF") not found"
    endif
end

