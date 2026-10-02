#!/bin/csh -f

# Usage:
#   clean_unused_from_hpp <path_from_root/file.hpp>
#
# Description:
#   This script analyzes a C++ header file from the root of the Git repository
#   and identifies methods that are not used anywhere in the codebase.
#
# Example:
#   cd ~/project_gstlearn/gstlearn
#   clean_unused_from_hpp include/Db/Db.hpp

# Check if an argument is provided
if ( $#argv != 1 ) then
    echo "Usage: $0 <path_from_root/file.hpp>"
    echo "Example: $0 include/Db/Db.hpp"
    exit 1
endif

set header = "$1"

# Verify that the specified header file exists
if ( ! -f "$header" ) then
    echo "Error: File '$header' does not exist."
    echo "Make sure you are running this script from the root of the repository."
    exit 1
endif

echo "=== Analyzing unused methods in $header ==="

# Extract method names, stripping destructors (~) and non-alphanumeric chars to prevent csh expansion
set methods = `grep -E '[a-zA-Z0-9_]+\s*\(' "$header" | awk -F'(' '{print $1}' | awk '{print $NF}' | tr -d '*&~' | grep -E '^[a-zA-Z_][a-zA-Z0-9_]*$' | sort -u`

# Search each extracted method across the entire Git repository
foreach method ( $methods )
    # Ignore C++ keywords if matched by mistake
    switch ( "$method" )
        case "if":
        case "while":
        case "for":
        case "switch":
        case "catch":
        case "return":
            goto next_item
    endsw

    # Count occurrences across the entire repository
    set count = `git grep -w "$method" . | wc -l`

    # If count <= 1, the method is only declared in the header and never used
    if ( $count <= 1 ) then
        echo "Unused: $method"
    endif

next_item:
end

echo "=== Analysis completed ==="
