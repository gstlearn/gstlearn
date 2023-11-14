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

def downloadRemoteFile(directory, filename, where):
    '''
    Return the complete name of a file:
    - if Internet is available, the file is retrieved from the web site
    - if not, it is assumed to be present locally
    
    directory: Name of the target directory (used for 'where' = "data", "None" otherwise)
    filename: Name of the file to be downloaded
    where: 'data' or 'graphics' or 'mdfile'
    
    Remarks
    -------
    When retrieving file, and according to 'where', the origin is:
    - when 'where' == "graphics"
         urlGST + "/references' + '/Figures' + '/'
    - when 'where' == "mdfile"
        urlGST + "/references' + '/'
    - when 'where' == "data"
        urlGST + "/data" + '/' + directory + '/' + filename
    '''
    if isInternetAvailable:

        if where == 'graphics':
            pathname = urlGST + '/references' + '/Figures' + '/' + filename
        elif where == 'mdfile':
            pathname = urlGST + '/references' + '/' + filename
        elif where == 'data':
            pathname = urlGST + '/data/' + directory + '/' + filename
        else:
            print("'downloadRemoteFile' does not know about 'where' = ", where)
        
        # The file is loaded in the local environment (with the same name)
        localname, head = urllib.request.urlretrieve(pathname, filename)
    
    else:
        localname = join('.', filename)
        
    return localname
    
def loadFigure(filename):
    '''
    This function displays the contents of the figure file named 'filename' (from the web site)
    
    Arguments
    ---------
    filename: Name of the figure of interest
    '''
    return downloadRemoteFile(None, filename, "graphics")
    
def loadDoc(filename):
    '''
    This function displays the contents of the Markdown file named 'filename' (from the web site)
    The result is decorated so as to appear as a NOTE in HTML
    
    Arguments
    ---------
    filename: Name of the file of interest
    '''
    
    filemd = downloadRemoteFile(None, filename, "mdfile")
        
    multilines = open(filemd, 'r').read()
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
            filefig = downloadRemoteFile(None, graphicFile, "graphics")
        
            # Reconstruct the new line (MD syntax) [pay attention to space]
            targetLine = "![Image](" + filefig + " 'title')"
        
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
    
    return downloadRemoteFile(directory, filename, "data")
