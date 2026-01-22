#!/bin/csh -f

# This script is meant to compare the different specific output Reference files
# for the different platforms.
# This script requires the environment variable GSTLEARN_DIR
#
set DIR = $GSTLEARN_DIR/gstlearn

set nargs = $#argv
if ($nargs < 2) then
		echo "Syntax:"
		echo "compare_test name1 nametest [ajout]"
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
		set radix   = "Result"
		set nameext = "ref"
		set namecmp = "out"
		set RADIX_AUX = $DIR/tests/data/myoutput
endif
if ($name1 == "cpp") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .cpp`
		set radix   = $nametest
		set nameext = "ref"
		set namecmp = "out"
endif
if ($name1 == "py") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .py`
		set radix   = $nametest
		set nameext = "ref"
		set namecmp = "out"
endif
if ($name1 == "r") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .r`
		set radix   = $nametest
		set nameext = "ref"
		set namecmp = "out"
endif
if ($name1 == "ipynb") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .ipynb`
		set radix   = $nametest
		set nameext = "asciidoc"
		set namecmp = "asciidoc"
endif
if ($name1 == "rmd") then
		set directref = "output"
		set directest = ""
		set directaux = "myoutput"
		set nametest = `basename $nametest .Rmd`
		set radix   = $nametest
		set nameext = "out"
		set namecmp = "out"
endif

set dirref = $RADIX_FROM/$directref
set diraux = $RADIX_AUX/$directaux
set dircmp = $RADIX_TO/$option/$directest

set FILEREF = $dirref/$radix.$nameext
set TITRE_REF = "Official Reference"

set FILECMP = $dircmp/$radix.$namecmp
if (-e $FILECMP) then
    diff -b $FILEREF $FILECMP > /dev/null
    if ($status == 1) then
        printf "\n%-60s   %s\n" "Official Reference" "Github Reference"
	    diff -b --side-by-side --suppress-common-lines $FILEREF $FILECMP
	else
		echo "File" $FILECMP "should be deleted ... as similar to Reference"
    endif
endif

set FILEAUX = $diraux/$radix.$nameext
if (-e $FILEAUX) then
    diff -b $FILEREF $FILEAUX > /dev/null
    if ($status == 1) then
        printf "\n%-60s   %s\n" "Official Reference" "My computer Reference"
	    diff -b --side-by-side --suppress-common-lines $FILEREF $FILEAUX
	else
		echo "File" $FILEAUX "should be deleted ... as similar to Reference"
    endif
endif

set PROP = "clang"
set FILECLANG = $dirref/$radix"_"$PROP.$nameext
if (-e $FILECLANG) then
    diff -b $FILEREF $FILECLANG > /dev/null
    if ($status == 1) then
        printf "\n%-60s   %s\n" "Official Reference" "Mac Reference"
	    diff -b --side-by-side --suppress-common-lines $FILEREF $FILECLANG
	else
		echo "File" $FILECLANG "should be deleted ... as similar to Reference"
    endif
endif

set PROP = "msvc"
set FILEMSVC = $dirref/$radix"_"$PROP.$nameext
if (-e $FILEMSVC) then
    diff -b $FILEREF $FILEMSVC > /dev/null
    if ($status == 1) then
        printf "\n%-60s   %s\n" "Official Reference" "Windows MSVC Reference"
	    diff -b --side-by-side --suppress-common-lines $FILEREF $FILEMSVC
	else
		echo "File" $FILEMSVC "should be deleted ... as similar to Reference"
    endif
endif

set PROP = "msys"
set FILEMSYS = $dirref/$radix"_"$PROP.$nameext
if (-e $FILEMSYS) then
    diff -b $FILEREF $FILEMSVC > /dev/null
    if ($status == 1) then
        printf "\n%-60s   %s\n" "Official Reference" "Windows MSYS Reference"
	    diff -b --side-by-side --suppress-common-lines $FILEREF $FILEMSYS
	else
		echo "File" $FILEMSYS "should be deleted ... as similar to Reference"
    endif
endif



