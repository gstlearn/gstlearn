import subprocess
import sys
import os
import re

# This script:
# - executes a jupyter notebook test script (argv[1]) 
# - converts it into a file saved in the given directory (argv[2]).
# - according to the format [asciidoc, html, pdf] (argv[3])
python_exe = os.path.realpath(sys.executable)
test_script = sys.argv[1]
out_dir = sys.argv[2]
test_name = os.path.splitext(os.path.basename(test_script))[0] # No extension
out_type = "asciidoc"
if (len(sys.argv) > 3):
    out_type = sys.argv[3]
test_output = os.path.join(out_dir, test_name + "." + out_type)

# Inspired from https://stackoverflow.com/questions/65502005/convert-a-jupyter-notebook-to-html-output-in-native-python
# See manual here : https://buildmedia.readthedocs.org/media/pdf/nbconvert/latest/nbconvert.pdf
import nbformat
import tempfile
from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert import ASCIIDocExporter
from nbconvert import HTMLExporter
from nbconvert import PDFExporter

# Read [and hack source notebook in a temporary notebook]
f = open(test_script, 'r', encoding='utf8')
if (out_type == "asciidoc"):
    # Kill some cells that pollute nonregression
    nbs = f.read()
    nbs = re.sub("from IPython.display import Markdown", "", nbs)
    nbs = re.sub("Markdown", "print", nbs)
    print(type(nbs))
    new_file, filename = tempfile.mkstemp(text=True)
    os.write(new_file, nbs.encode('utf8'))
    os.close(new_file)
    f = open(filename, 'r', encoding='utf8')

# Really read the notebook
nb = nbformat.read(f, as_version=nbformat.NO_CONVERT)

# Execute the Notebook
ep = ExecutePreprocessor(timeout=-1, kernel_name='python3')
ep.preprocess(nb)

# Export to asciidoc (dump only output cells for test purpose)
if (out_type == "asciidoc"):
    exporter = ASCIIDocExporter()
    exporter.exclude_input = True
    exporter.exclude_input_prompt = True
    exporter.exclude_markdown = True
    exporter.exclude_raw = True
    exporter.exclude_unknown = True
# Export to HTML 
elif (out_type == "html"):
    exporter = HTMLExporter()
    # Ensure that equations and 3D is well displayed!
    exporter.mathjax_url = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/latest.js?config=TeX-MML-AM_CHTML"
    exporter.require_js_url = "https://requirejs.org/docs/release/2.3.6/minified/require.js"
# Export to PDF
elif (out_type == "pdf"):
    exporter = PDFExporter()
else:
    print("Wrong output file type for run_test_ipynb.py [asciidoc, html, pdf]")
    
# Export the notebook
notebook_node, resources = exporter.from_notebook_node(nb)

# Post treatment for asciidoc (for test comparison)
if (out_type == "asciidoc"):
    # TODO : I would prefer using this feature, but don't know how to do with Exporters :
    # https://stackoverflow.com/questions/52834910/remove-cells-from-jupyter-notebook-with-nbconvert
    
    # Remove all graphical 3D object identifiers from the output ascii file, i.e. :
    # [[e43b6f2f-ba2b-47f7-8a13-2336077446d1]]
    # -----<matplotlib.collections.QuadMesh at 0x7f3a056e6320>
    # ----[<matplotlib.lines.Line2D at 0x7f3918a78550>]
    notebook_node = re.sub("[a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{12}", "XXX", notebook_node)
    notebook_node = re.sub("----<matplotlib.*", "XXX", notebook_node)
    notebook_node = re.sub("----\[<matplotlib.*", "XXX", notebook_node)
    
    # Remove images in base64 included by MD files
    # data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAZAAAAGQCAYAA........
    notebook_node = re.sub(".*data:image/png;base64,.*", "data:image/png;base64,XXX", notebook_node)
    
    # Remove all lines coming from data downloading (no need anymore as wget -q):
    # --2023-07-28 15:20:09--  https://soft.minesparis.psl.eu/gstlearn/data/Scotland/Scotland_Temperatures.NF
    # Résolution de soft.minesparis.psl.eu (soft.minesparis.psl.eu)… 51.83.45.127
    # Connexion à soft.minesparis.psl.eu (soft.minesparis.psl.eu)|51.83.45.127|:443… connecté.
    # requête HTTP transmise, en attente de la réponse… 200 OK
    # Taille : 5291 (5,2K)
    # Enregistre : ‘Scotland_Temperatures.NF.1’
    # Scotland_Temperatur 100%[===================>]   5,17K  --.-KB/s    ds 0s      
    # 2023-07-28 15:20:10 (1,31 GB/s) - ‘Scotland_Temperatures.NF.1’ enregistré [5291/5291]
    #notebook_node = re.sub(".*-  https://soft.minesparis.psl.eu/gstlearn/data.*", "XXX", notebook_node)
    #notebook_node = re.sub(".*51\.83\.45\.127.*", "XXX", notebook_node)
    #notebook_node = re.sub(".*HTTP transmise.*", "XXX", notebook_node)
    #notebook_node = re.sub("Taille.*", "XXX", notebook_node)
    #notebook_node = re.sub("Enregistre.*", "XXX", notebook_node)
    #notebook_node = re.sub(".*KB/s.*", "XXX", notebook_node)
    #notebook_node = re.sub(".*GB/s.*", "XXX", notebook_node)
    
    # Remove a specific warning in Tuto_SpatioTemp.ipynb, i.e. :
    # /tmp/ipykernel_24563/4216505814.py:15: CholmodTypeConversionWarning: converting matrix of class csr_matrix to CSC format
    notebook_node = re.sub(".*CholmodTypeConversionWarning", "XXX: CholmodTypeConversionWarning", notebook_node)
    
    # Remove panda frame decoration that can vary according the version/OS i.e. :
    notebook_node = re.sub("\\|====+", "|===", notebook_node)
    
    # Remove pip install output
    notebook_node = re.sub("\[notice\].*", "#NO_DIFF#XXX", notebook_node)
    notebook_node = re.sub(".*site-packages is not writeable", "#NO_DIFF#XXX", notebook_node)
    notebook_node = re.sub("Requirement already satisfied.*", "#NO_DIFF#XXX", notebook_node)
    # Remove this: *Out[3]:*
    notebook_node = re.sub("\*Out\[[0-9]+\]:\*", "#NO_DIFF#XXX", notebook_node)
    # Remove this: [png]
    notebook_node = re.sub("\!\[png\].*", "#NO_DIFF#XXX", notebook_node)

    
# Write to output file
if (out_type == "pdf"):
    with open(test_output, "wb") as f:
        f.write(notebook_node)
else:
    with open(test_output, "w", encoding='utf8') as f:
        f.write(notebook_node)
