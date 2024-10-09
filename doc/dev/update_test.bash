#!/bin/bash

# This script is meant to update the reference output file
# starting from the resulting test output file (in Release mode only)

DIR="$GSTLEARN_DIR/gstlearn"

nargs=$#
if [ $nargs -lt 2 ]; then
    echo "Update the reference output file from the build test output file (Release mode only)"
    echo "Syntax:"
    echo "update_test name1 nametest"
    echo "   name1    : Name of the subContainer (e.g.: data, ipynb, py, r, rmd, cpp)"
    echo "   nametest : Name of the test (e.g.: Jeu1, Tuto_Plot2D, Hermite, etc...)"
    exit
fi

name1="$1"
nametest="$2"
option="Release"

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

cp "$FILECMP" "$FILETMP"
