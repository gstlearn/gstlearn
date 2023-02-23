import subprocess
import sys
import os
import re

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
# https://github.com/mwouts/jupytext/issues/770#issuecomment-1257559895
with open(test_script, encoding='utf8') as f:
    nb = nbformat.read(f, as_version=4)

# execute notebook
ep = ExecutePreprocessor(timeout=-1, kernel_name='python3')
ep.preprocess(nb)

# export to asciidoc (dump only output cells)
ascii_exporter = ASCIIDocExporter()
ascii_exporter.exclude_input = True
ascii_exporter.exclude_markdown = True
ascii_exporter.exclude_raw = True
ascii_exporter.exclude_unknown = True
ascii_data, resources = ascii_exporter.from_notebook_node(nb)

# Remove all graphical 3D object identifiers from the output ascii file, i.e. :
# [[e43b6f2f-ba2b-47f7-8a13-2336077446d1]]
ascii_data = re.sub("[a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{12}", "XXX", ascii_data)

# Remove a specific warning in Tuto_SpatioTemp.ipynb, i.e. :
# /tmp/ipykernel_24563/4216505814.py:15: CholmodTypeConversionWarning: converting matrix of class csr_matrix to CSC format
ascii_data = re.sub(".*CholmodTypeConversionWarning", "XXX: CholmodTypeConversionWarning", ascii_data)

# Remove panda frame decoration that can vary according the version/OS i.e. :
ascii_data = re.sub("\\|====+", "|===", ascii_data)


# write to output file
with open(test_output, "w", encoding='utf8') as f:
    f.write(ascii_data)

