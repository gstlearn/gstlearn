import subprocess
import sys
import os

# This script execute a jupyter notebook test script (argv[1]) and convert it to an asciidoc file saved in the given directory (argv[2])
python_exe = os.path.realpath(sys.executable)
test_script = sys.argv[1]
out_dir = sys.argv[2]
test_name = os.path.splitext(os.path.basename(test_script))[0] # No extension
test_output = os.path.join(out_dir, test_name + ".asciidoc")

# Inspired from https://stackoverflow.com/questions/65502005/convert-a-jupyter-notebook-to-html-output-in-native-python
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert import ASCIIDocExporter

# read source notebook
with open(test_script) as f:
    nb = nbformat.read(f, as_version=4)

# execute notebook
ep = ExecutePreprocessor(timeout=-1, kernel_name='python3')
ep.preprocess(nb)

# export to Asciidoc
ascii_exporter = ASCIIDocExporter()
ascii_exporter.exclude_input = True
ascii_data, resources = ascii_exporter.from_notebook_node(nb)

# write to output file
with open(test_output, "w") as f:
    f.write(ascii_data)

#jupyter nbconvert --output-dir test_output --no-input --log-level WARN --to asciidoc --execute ../../gstlearn/doc/demo/python/Inheritance.ipynb