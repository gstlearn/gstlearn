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
import requests
import os
import re
import base64
import gstlearn as gl
from IPython.display import display, Javascript

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

Usage:

    %load_ext disable_autoscroll

You can also put the js snippet below in profile_dir/static/js/custom.js
"""

disable_js = """
IPython.OutputArea.prototype._should_scroll = function(lines) {
    return false;
}
"""
def setNoScroll():
    display(Javascript(disable_js))

def locateFile(filename, where='references', directory=None, verbose=False):
    '''
    Return the absolute path of a file:
    - it is assumed to be present locally in '.' ('where' and 'directory' are ignored)
    - if not, it is assumed to be present locally in './doc/<where>', '../../doc/<where>' or '../../<where>'
    - if not, if Internet is available, the file is downloaded from the gstlearn website in a temporary file
    
    filename: Name of the file to be located
    where: 'data' or 'references'
    directory: Name of the data file directory (only used for 'where' = "data")
    '''
    
    if (verbose):
        print("Current directory is", os.getcwd())

    # Test current directory
    localname = os.path.join('.', filename)
    if os.path.isfile(localname):
        fullname = os.path.abspath(localname)
        if (verbose):
            print(localname, "found... Full path is", fullname)
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
                print(localname, "found... Full path is", fullname)
            return fullname
        elif (verbose):
            print(localname, "not found...")
                
    if not internetAvailable():
        print("Error: Cannot access to", filename, "(no Internet)!")
        return None
    
    # Download from Internet in a temporary file
    localname = urlGST + '/' + package_version + '/' + where + '/' + filename
    try:
        fullname, head = urllib.request.urlretrieve(localname)
        if (verbose):
            print(localname, "found... Full path is", fullname)
        return fullname
    except:
        pass
    
    print("Cannot access URL:", localname, "!")
    return None

def loadDoc(filename, verbose=False):
    '''
    This function return the contents of a Markdown file from the 'references' directory named 'filename'
    The result is decorated so as to appear as a NOTE in HTML files
    
    Arguments
    ---------
    filename: Name of the Markdown file of interest
    '''
    
    filemd = locateFile(filename, verbose=verbose)
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
            imgfile = locateFile(img.group(3), verbose=verbose)
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
    
def loadData(directory, filename, verbose=False):
    '''
    This function loads a file named 'filename' in the 'directory' (from the web site)
    
    Arguments
    ---------
    directory: Name of the Directory (within /data) containing the file of interest
    filename: Name of the file of interest
    '''
    
    return locateFile(filename, "data", directory, verbose=verbose)
