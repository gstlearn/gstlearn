################################################################################
#                                                                              #
#                         gstlearn Python package                              #
#                                                                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# Website: https://gstlearn.org                                                #
# License: BSD 3-clause                                                        #
#                                                                              #
################################################################################
# This part is meant to facilitate the use of the web site where information
# is collected 

import urllib.request
import urllib.error
import os
import re
import base64

import gstlearn as gl
from IPython.display import HTML


try:
#   from IPython.display import display, Javascript, Markdown
    from IPython import get_ipython
    import requests
except ModuleNotFoundError as ex:
    msg = ("Python dependencies 'IPython' and 'requests' not found.\n"
          "To install them alongside gstlearn, please run `pip install gstlearn[doc]'")
    raise ModuleNotFoundError(msg) from ex

# The various pieces of documentation are supposed to be located
# at the following URL
urlGST = 'https://soft.minesparis.psl.eu/gstlearn'
package_version = gl.__version__

# Next lines are used to decorate the MD files for rendering documentation
header = [
  "<style>md-block { color:gray; background-color:white; }</style>",
  "<md-block>\n"]
trailer = ["</md-block>"]

def internetAvailable(timeout=1):
    '''
    Check if Internet is available
    
    This function requires the package 'requests' to be installed
    
    Returns:
    --------
    bool: True if Internet is available and False otherwise
    '''
    try:
        requests.head(urlGST, timeout=timeout)
        return True
    except requests.ConnectionError:
        return False
    
"""
Extension for disabling autoscrolling long output, which is super annoying sometimes

"""
from IPython.display import HTML, Javascript, display

def setNoScroll():
    css = """
    <style id="no-scroll-outputs">
    /* cibles communes JupyterLab / Notebook 7 et fallback classique */
    .jp-OutputArea-child,
    .jp-OutputArea-output,
    .jp-OutputArea,
    .jp-RenderedText,
    .jp-RenderedHTMLCommon,
    .output_scroll,
    .output_area .output_subarea {
        max-height: none !important;
        height: auto !important;
        min-height: 0 !important;
        overflow: visible !important;
    }
    /* blocs <pre> / code */
    .jp-OutputArea-child pre,
    .jp-OutputArea-output pre,
    .jp-RenderedHTMLCommon pre,
    .output_scroll pre {
        max-height: none !important;
        overflow: visible !important;
    }
    </style>
    """
    js = r"""
    (function() {
      const selectors = [
        '.jp-OutputArea-child',
        '.jp-OutputArea-output',
        '.jp-OutputArea',
        '.output_scroll',
        '.output_area .output_subarea'
      ];
      function unscroll(node){
        try{
          if(!(node instanceof HTMLElement)) return;
          node.style.maxHeight = 'none';
          node.style.overflow = 'visible';
          node.style.height = 'auto';
          node.style.minHeight = '0';
        }catch(e){}
      }
      function processAll(){
        selectors.forEach(sel => {
          document.querySelectorAll(sel).forEach(unscroll);
        });
      }
      processAll();
      // Observer : annule les styles inline quand des sorties sont ajoutées / modifiées
      const mo = new MutationObserver(muts => {
        muts.forEach(m => {
          m.addedNodes.forEach(n => {
            if(!(n instanceof HTMLElement)) return;
            selectors.forEach(sel => {
              if(n.matches && n.matches(sel)) unscroll(n);
              if(n.querySelectorAll) n.querySelectorAll(sel).forEach(unscroll);
            });
          });
          if(m.type === 'attributes' && m.target){
            selectors.forEach(sel => {
              if(m.target.matches && m.target.matches(sel)) unscroll(m.target);
            });
          }
        });
      });
      mo.observe(document, { childList: true, subtree: true, attributes: true, attributeFilter: ['style','class'] });
      // garder une référence globale pour éviter GC
      window.__no_scroll_observer = mo;
    })();
    """
    if get_ipython() and get_ipython().__class__.__name__ == 'ZMQInteractiveShell':
        from IPython.display import Javascript, display
        display(HTML(css))
        display(Javascript(js))

def locateFile(filename, where='references', directory=None, verbose=False, version=package_version):
    '''
    Return the absolute path of a file:
    - it is assumed to be present locally in '.' ('where' and 'directory' are ignored)
    - if not, it is assumed to be present locally in './doc/<where>', '../../doc/<where>' or '../../<where>'
    - if not, if the GSTLEARN_DIR environment variable is defined, it is assumed to be present in '<GSTLEARN_DIR>/gstlearn/doc/<where>'
    - if not, if Internet is available, the file is downloaded from the gstlearn website in a temporary file
    
    filename: Name of the file to be located
    where: 'data' or 'references'
    directory: Name of the data file directory (only used for 'where' = "data")
    verbose: True to activate verbose mode
    version: Use a specific gstlearn version when searching the file on the web (string)
    '''
    
    argfilename = filename
    if (verbose):
        print("Current directory is", os.getcwd())

    # Test current directory
    localname = os.path.join('.', filename)
    if os.path.isfile(localname):
        fullname = os.path.abspath(localname)
        if (verbose):
            print(filename, "found... Full path is", fullname)
        return fullname
    elif (verbose):
        print(localname, "not found...")

    # Test locally in other directories
    if where not in ['references', 'data']:
        print("'locateFile' does not know about 'where' = ", where)
        return None
    if where == 'data' and directory is not None:
        filename = os.path.join(directory, filename)
    
    folders = [os.path.join('.',"doc",where),
               os.path.join('..','..',"doc",where),
               os.path.join('..','..',where)]
    for f in folders:
        localname = os.path.join(f, filename)
        if os.path.isfile(localname):
            fullname = os.path.abspath(localname)
            if (verbose):
                print(filename, "found... Full path is", fullname)
            return fullname
        elif (verbose):
            print(localname, "not found...")

    # Test in GSTLEARN_DIR environment variable
    if os.environ.get('GSTLEARN_DIR') is not None:
        gstlearn_dir = os.environ.get('GSTLEARN_DIR')
        if gstlearn_dir is not None:
            localname = os.path.join(gstlearn_dir, 'gstlearn', 'doc', where, filename)
            if os.path.isfile(localname):
                fullname = os.path.abspath(localname)
                if (verbose):
                    print(filename, "found... Full path is", fullname)
                return fullname
            elif (verbose):
                print(localname, "not found...")

    # Test on the web
    if not internetAvailable():
        print("Error: Cannot access to", filename, "(no Internet)!")
        return None
    
    # Download from Internet in a temporary file
    localname = urlGST + '/' + version + '/' + where + '/' + directory + '/' + argfilename
    try:
        fullname, head = urllib.request.urlretrieve(localname)
        if (verbose):
            print(localname, "found... Full path is", fullname)
        return fullname
    except:
        pass
    
    print("Cannot access URL:", localname, "!")
    return None

def loadDoc(filename, verbose=False, version=package_version):
    '''
    This function return the contents of a Markdown file from the 'references' directory named 'filename'
    The result is decorated so as to appear as a NOTE in HTML files
    
    Arguments
    ---------
    filename: Name of the Markdown file of interest
    verbose: True to activate verbose mode
    version: Use a specific gstlearn version when searching the file on the web (string)
    '''
    
    filemd = locateFile(filename, verbose=verbose, version=version)
    if filemd is None:
        return "File " + filename + " not found!"
    
    multilines = open(filemd, 'r').read()
    lines = multilines.split('\n')
    
    # Capture Markdown images (beginning ![description](filename) ending)
    pattern = re.compile(r'(.*)\!\[(.*)\]\((.+)\)(.*)')
    for i in range(len(lines)):
        targetLine = lines[i]
        img = pattern.search(targetLine)
        if img is not None:
            beginning = img.group(1)
            imgdesc = img.group(2)
            imgfile = locateFile(img.group(3), verbose=verbose, version=version)
            ending = img.group(4)
            if imgfile is None:
                return "File " + img.group(3) + " not found!"
            # Convert in base64 for embedding the image
            with open(imgfile, 'rb') as image_file:
                imgfile = base64.b64encode(image_file.read())
            # Reconstruct the full Markdown line
            lines[i] = beginning + '![' + imgdesc + '](data:image/png;base64,' + imgfile.decode() + ')' + ending
    
    result = ''.join(header) + '\n'.join(lines) + ''.join(trailer)
    return result

def displayDoc(filename, verbose=False, version=package_version):
    '''
    This function displays the contents of a Markdown file from the 'references' directory named 'filename'
    The result is decorated so as to appear as a NOTE in HTML files
    
    Arguments
    ---------
    filename: Name of the Markdown file of interest
    verbose: True to activate verbose mode
    version: Use a specific gstlearn version when searching the file on the web (string)
    '''

    result = loadDoc(filename, verbose, version)

    return Markdown(result)

def loadData(directory, filename, verbose=False, version=package_version):
    '''
    This function returns the path of a file 'filename' in the 'data' directory (locally or from the web site)
    
    Arguments
    ---------
    filename: Name of the file of interest
    directory: Name of the sub-directory (within 'data' directory) containing the file of interest
    verbose: True to activate verbose mode
    version: Use a specific gstlearn version when searching the file on the web (string)
    '''
    
    return locateFile(filename, "data", directory, verbose=verbose, version=version)
