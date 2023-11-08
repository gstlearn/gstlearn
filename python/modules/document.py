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
urlMP = 'https://soft.minesparis.psl.eu/gstlearn'

def isInternetAvailable(timeout=1):
    '''
    Check if Internet is available
    
    This function requires the package 'requests' to be installed
    
    Returns:
    --------
    bool: True if Internet is available and False otherwise
    '''
    try:
        requests.head(urlMP, timeout=timeout)
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
    
def loadDoc(filename):
    '''
    This function displays the contents of the Markdown file named 'filename' (from the web site)
    
    Arguments
    ---------
    filename: Name of the file of interest
    '''
    
    if isInternetAvailable():
        pathname = urlMP + '/references/' + filename
        filepath, head = urllib.request.urlretrieve(pathname)
    else:
        filepath = join('.', filename)
    return filepath
    
def loadData(directory, filename):
    '''
    This function loads a file named 'filename' in the 'directory' (from the web site)
    
    Arguments
    ---------
    directory: Name of the Directory (within /data) containing the file of interest
    filename: Name of the file of interest
    '''
    
    if isInternetAvailable():
        pathname = urlMP + '/data/' + directory + '/' + filename
        filepath, head = urllib.request.urlretrieve(pathname)
    else:
        filepath = join('.', directory, filename)
    return filepath
