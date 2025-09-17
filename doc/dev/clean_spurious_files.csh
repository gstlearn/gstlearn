#!/bin/csh -f

# This script is meant to delete all files named '2' explicitely
# It can be used to suppress all files given their names (exact matching)
#
# The search is performed from the current directory
# This script requires the environment variable GSTLEARN_DIR
#
# Syntax: clean_spurious_files [dirname]
# where 'dirname' stands for the directory to be searched (default: "gstlearn")
#
set nargs = $#argv
if ($nargs > 0) then
		set maindir = "$argv[1]"
else
		set maindir = "gstlearn"
endif

set DIR = $GSTLEARN_DIR/$maindir
cd $DIR

# Check the list of targeted files:
find . -type f -name 2 

# Delete the targeted files with control
find . -type f -name 2 -exec rm -i {} \;

# Delete the targeted files with no control (dangerous)
find . -type f -name 2 -exec rm -rf {} \;
