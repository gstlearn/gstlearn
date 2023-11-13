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

import numpy                 as np
import numpy.ma              as ma
import gstlearn              as gl
import urllib.request
import requests
from os.path import join
from IPython.display import display, Javascript

# The various pieces of documentation are supposed to be located
# at the following URL
urlGST = 'https://soft.minesparis.psl.eu/gstlearn'

# Next lines are used to decorate the MD files for rendering documentation
header = [
  "<style>md-block { color:gray; background-color:white; }</style>",
  "<md-block>\n"]
trailer = ["</md-block>"]

def isInternetAvailable(timeout=1):
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

def loadFigure(filename):
    '''
    This function displays the contents of the figure file named 'filename' (from the web site)
    
    Arguments
    ---------
    filename: Name of the figure of interest
    '''
    if isInternetAvailable():
        pathname = urlGST + '/references' + '/Figures/' + filename
        filepath, head = urllib.request.urlretrieve(pathname)
    else:
        filepath = join('.', filename)
        
    return filepath
    
def loadDoc(filename, useURL=True):
    '''
    This function displays the contents of the Markdown file named 'filename' (from the web site)
    The result is decorated so as to appear as a NOTE in HTML
    
    Arguments
    ---------
    filename: Name of the file of interest
    useURL: TRUE if the file must be found from the WEB; FALSE otherwise
    '''
    
    if isInternetAvailable() and useURL:
        pathname = urlGST + '/references/' + filename
        filepath, head = urllib.request.urlretrieve(pathname)
    else:
        filepath = join('.', filename)
        
    multilines = open(filepath, 'r').read()
    lines = multilines.split('\n')

    searchItem = "(Figure"
    for i in range(len(lines)):
        targetLine = lines[i]
    
        # Look for the graphic dependency
        if searchItem in targetLine:

            begin = targetLine.index(searchItem) + 1
            start = begin + len(searchItem) + 1
            end   = targetLine.index(")")
        
            # Extract the name of the Graphic File
            graphicFile = targetLine[start:end]
            pathname = urlGST + '/references' + '/Figures/' + graphicFile
        
            # Reconstruct the new line
            targetLine = "<img src='" + pathname + "' />"
        
        lines[i] = targetLine
    
    tata = "\n"
    new_multilines = tata.join([i for i in lines[0:]])

    result = ''.join(header) + new_multilines + ''.join(trailer)
    return result
    
def loadData(directory, filename):
    '''
    This function loads a file named 'filename' in the 'directory' (from the web site)
    
    Arguments
    ---------
    directory: Name of the Directory (within /data) containing the file of interest
    filename: Name of the file of interest
    '''
    
    if isInternetAvailable():
        pathname = urlGST + '/data/' + directory + '/' + filename
        filepath, head = urllib.request.urlretrieve(pathname)
    else:
        filepath = join('.', directory, filename)
    return filepath
