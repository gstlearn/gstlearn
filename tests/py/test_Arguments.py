#
# This file is meant to test the reading of various types of arguments in Python
#

# Loading the packages

import gstlearn as gl
import numpy as np
import os
import sys

# Testing main argument types

gl.argumentTestInt(12)
gl.argumentTestDouble(2.3)
gl.argumentTestVectorInt([1,2,3])
gl.argumentTestVectorDouble([1.1, 2.2, 3.3])
gl.argumentTestString("my_String")
gl.argumentTestVectorString("my_String")  # The String is decomposed (cannot be fixed)
gl.argumentTestVectorString(["my_String1","my_String2","my_String3"])
gl.argumentTestVectorVectorInt([ [2,3],[1, 5 ] ])
gl.argumentTestVectorVectorDouble([ [2.,3.], [1., 5 ] ])

# Testing missing arguments

gl.argumentTestInt(np.nan)
gl.argumentTestDouble(np.nan)
gl.argumentTestVectorInt([np.nan])
gl.argumentTestVectorDouble([np.nan])

# Testing overloading of methods

gl.argumentTestIntOverload(12)
gl.argumentTestIntOverload((21, 32))
gl.argumentTestDoubleOverload(2.)
gl.argumentTestDoubleOverload((2., 3.))
gl.argumentTestStringOverload("my_String")
gl.argumentTestStringOverload(["my_String1","my_String2","my_String3"])

# Testing ENUM

gl.argumentTestEnum(gl.ETests.CASE2)

# Testing Returning arguments

print(gl.argumentReturnInt(12))
print("nan") if (gl.isNaN(gl.argumentReturnInt(np.nan))) else print("oups") # No NaN value for integers so use isNaN
print(gl.argumentReturnDouble(21.4))
print(gl.argumentReturnDouble(np.nan))

# Testing assessors to the elements of a class

myClass = gl.argClass()
myClass.display()
myClass.ival
myClass.ival = 2
myClass.display()

print("Test successfully performed")
