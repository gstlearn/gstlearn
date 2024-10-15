#!/bin/csh -f

# This script is meant to update the reference output file starting from the resulting test file output file Release mode only)
# This script requires the environment variable GSTLEARN_DIR
# 
set DIR = $GSTLEARN_DIR/gstlearn

set nargs = $#argv
if ($nargs < 2) then
		echo "Update the reference output file from the build test output file (Release mode only)"
		echo "Syntax:"
		echo "update_test name1 nametest"
		echo "   name1    : Name of the subContainer (e.g.: data, ipynb, py, r, rmd, cpp)"
		echo "   nametest : Name of the test (e.g.: Jeu1, Tuto_Plot2D, Hermite, etc...)"
		exit
endif

set name1 = $argv[1]
set nametest = $argv[2]
set option = "Release"

set RADIX_FROM = $DIR/tests/$name1
set RADIX_AUX  = $DIR/tests/$name1
set RADIX_TO   = $DIR/build/tests/$name1

if ($name1 == "data") then
		set directref = $nametest
		set directest = $nametest
		set directaux = $nametest
		set nameref = "Result.ref"
		set namecmp = "Result.out"
		set RADIX_AUX = $DIR/tests/data/myoutput
endif
if ($name1 == "cpp") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .cpp`
		set nameref = $nametest.ref
		set namecmp = $nametest.out
endif
if ($name1 == "py") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .py`
		set nameref = $nametest.ref
		set namecmp = $nametest.out
endif
if ($name1 == "r") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .r`
		set nameref = $nametest.ref
		set namecmp = $nametest.out
endif
if ($name1 == "ipynb") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .ipynb`
		set nameref = $nametest.asciidoc
		set namecmp = $nametest.asciidoc
endif
if ($name1 == "rmd") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .Rmd`
		set nameref = $nametest.out
		set namecmp = $nametest.out
endif

set FILEREF = $RADIX_FROM/$directref/$nameref
set FILEAUX = $RADIX_AUX/$directaux/$nameref
set FILETMP = $FILEREF
if (-e $FILEAUX) then
	set FILETMP = $FILEAUX
	echo ">>> Using non-standard reference"
endif

set FILECMP = $RADIX_TO/$option/$directest/$namecmp

cp $FILECMP $FILETMP
