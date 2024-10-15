#!/bin/bash

# This script is meant to compare a resulting test file to its reference (Release mode only)
# This script requires the environment variable GSTLEARN_DIR
#
DIR="$GSTLEARN_DIR/gstlearn"

nargs=$#
if [ $nargs -lt 2 ]; then
    echo "Syntax:"
    echo "compare_test name1 nametest [ajout]"
    echo "   name1    : Name of the subContainer (e.g.: data, ipynb, py, r, rmd, cpp)"
    echo "   nametest : Name of the test (e.g.: Jeu1, Tuto_Plot2D, Hermite, etc...)"
    echo "   ajout    : additional parameters for diff"
    echo "              'full' for comparing side by side"
    echo "              'large' for comparing side-by-side with larger output (200)"
    exit
fi

name1="$1"
nametest="$2"
option="Release"

ajout="None"
if [ $nargs -gt 2 ]; then
    ajout="$3"
else
    ajout=""
fi

RADIX_FROM="$DIR/tests/$name1"
RADIX_AUX="$DIR/tests/$name1"
RADIX_TO="$DIR/build/tests/$name1"

if [ "$name1" == "data" ]; then
    directref="$nametest"
    directest="$nametest"
    directaux="$nametest"
    nameref="Result.ref"
    namecmp="Result.out"
    RADIX_AUX="$DIR/tests/data/myoutput"
elif [ "$name1" == "cpp" ]; then
    directref="output"
    directest=""
    directaux="myoutput"
    nametest=$(basename "$nametest" .cpp)
    nameref="$nametest.ref"
    namecmp="$nametest.out"
elif [ "$name1" == "py" ]; then
    directref="output"
    directest=""
    directaux="myoutput"
    nametest=$(basename "$nametest" .py)
    nameref="$nametest.ref"
    namecmp="$nametest.out"
elif [ "$name1" == "r" ]; then
    directref="output"
    directest=""
    directaux="myoutput"
    nametest=$(basename "$nametest" .r)
    nameref="$nametest.ref"
    namecmp="$nametest.out"
elif [ "$name1" == "ipynb" ]; then
    directref="output"
    directest=""
    directaux="myoutput"
    nametest=$(basename "$nametest" .ipynb)
    nameref="$nametest.asciidoc"
    namecmp="$nametest.asciidoc"
elif [ "$name1" == "rmd" ]; then
    directref="output"
    directest=""
    directaux="myoutput"
    nametest=$(basename "$nametest" .Rmd)
    nameref="$nametest.out"
    namecmp="$nametest.out"
fi

FILEREF="$RADIX_FROM/$directref/$nameref"
FILEAUX="$RADIX_AUX/$directaux/$nameref"
FILETMP="$FILEREF"

if [ -e "$FILEAUX" ]; then
    FILETMP="$FILEAUX"
    echo ">>> Using non-standard reference"
fi

FILECMP="$RADIX_TO/$option/$directest/$namecmp"

if [ "$ajout" == "full" ]; then
    diff -b "$FILETMP" "$FILECMP" --side | more
elif [ "$ajout" == "large" ]; then
    diff -b "$FILETMP" "$FILECMP" --width=200 --side | more
else
    diff -b "$FILETMP" "$FILECMP"
fi
