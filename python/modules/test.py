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
# This part is meant to distribute (as is) a set of functions written in Python
# used for internal check 

import numpy                 as np
import numpy.ma              as ma
import gstlearn              as gl

def checkEqualityVectors(vec1, vec2, 
                         tolerance = gl.EPSILON6,
                         flagRelative = True,
                         flagAbsolute = False,
                         message="checkEqualityVectors"):
    
    '''
    This function checks that two vectors contain exactly the same set of values
    
    Arguments
    ---------
    vec1, vec2: two vectors (or matrices) to be compared. They must have the same number of elements
    tolerance: Tolerance used for the comparison
    flagRelative: when True, the values are compared without paying attention to their sign
    flagAbsolute: when True, test is run on absolute difference; otherwise Relative difference
    message: Message to be displayed when the vectors are not similar
    
    Remarks
    -------
    When the two vectors do not share the same dimension, the test is not performed
    and a message is printed.
    '''
    # Flatten the input vectors if necessary
    vec1 = vec1.flatten()
    vec2 = vec2.flatten()
    
    # Check that the two vectors have the same dimension
    long1 = len(vec1)
    long2 = len(vec2)
    if long1 != long2:
        print(f"{message}: Impossible to compare vectors of different dimensions")
        return False
    
    # Check is performed on the absolute value of each term of each vector
    if flagAbsolute:
        vec1 = np.absolute(vec1)
        vec2 = np.absolute(vec2)
    
    # Evaluate the comparison test
    diff = (vec1 - vec2)
    if flagRelative: 
        diff = diff / (vec1 + vec2 + tolerance)
    value = np.linalg.norm(diff)
    
    if value >= tolerance:
        print(f"{message}: Experimental value =",value," is larger than tolerance (",tolerance,")")
        return False
    
    return True

def checkEqualityDb(db, name1, name2,
                    tolerance=gl.EPSILON6,
                    flagRelative = True,
                    flagAbsolute = False,
                    message="checkEqualityDb"):
    
    return checkEqualityVectors(db[name1], db[name2], tolerance, 
                                flagRelative, flagAbsolute, message)

def checkEqualityValues(vec1, vec2, 
                        tolerance = gl.EPSILON6,
                        flagRelative = True,
                        flagAbsolute = False,
                        message="checkEqualityVectors"):
    '''
    This function checks that two values are equal
    
    Arguments
    ---------
    vec1, vec2: two values to be compared.
    tolerance: Tolerance used for the comparison
    flagRelative: when True, the values are compared without paying attention to their sign
    flagAbsolute: when True, test is run on absolute difference; otherwise Relative difference
    message: Message to be displayed when the vectors are not similar
    '''
    # Check is performed on the absolute value of each term of each vector
    if flagAbsolute:
        vec1 = np.absolute(vec1)
        vec2 = np.absolute(vec2)
    
    # Evaluate the comparison test
    diff = (vec1 - vec2)
    if flagRelative: 
        diff = diff / (vec1 + vec2 + tolerance)
    value = np.linalg.norm(diff)
    
    if value >= tolerance:
        print(f" {message}: Experimental value =",value," is larger than tolerance (",tolerance,")")
        return False
        
    return True
