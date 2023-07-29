import subprocess
import sys
import os

# This script retrieve the standard output of a python test script (argv[1]) and print it into a log file saved in the given directory (argv[2])
python_exe = os.path.realpath(sys.executable)
test_script = sys.argv[1]
out_dir = sys.argv[2]
test_name = os.path.splitext(os.path.basename(test_script))[0] # No extension
test_output = os.path.join(out_dir, test_name + ".out")

# Retrieve all standard outputs (C, python, etc..) and print to log file
# Use unbuffered (-u) output for "printing" on the fly
str_output = subprocess.check_output([python_exe, "-u", test_script], encoding='utf-8')
with open(test_output, "w+", encoding='utf8') as output:
  output.write(str_output)

