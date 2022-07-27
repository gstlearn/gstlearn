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
gl.argumentTestVectorString("my_String")  # Should be corrected
gl.argumentTestVectorString(["my_String1","my_String2","my_String3"])
gl.argumentTestVectorVectorInt([ [2,3],[1, 5 ] ])
gl.argumentTestVectorVectorDouble([ [2.,3.], [1., 5 ] ])

# Testing missing arguments

gl.argumentTestInt(gl.ITEST)
gl.argumentTestDouble(gl.TEST)
gl.argumentTestVectorInt([gl.ITEST])
gl.argumentTestVectorDouble([gl.TEST])

# Testing overloading of methods

gl.argumentTestIntOverload(12)
#gl.argumentTestIntOverload([21, 32])
gl.argumentTestIntOverload((21, 32))
gl.argumentTestStringOverload("my_String")
gl.argumentTestStringOverload(["my_String1","my_String2","my_String3"])

# Testing ENUM

gl.argumentTestEnum(gl.ETests.CASE2)

# Testing Returning arguments

print(gl.argumentReturnInt(12))
print(gl.argumentReturnInt(gl.ITEST))
print(gl.argumentReturnDouble(21.4))
print(gl.argumentReturnDouble(gl.TEST))

# Testing assessors to the elements of a class

myClass = gl.argClass()
myClass.display()
myClass.ival
myClass.ival = 2
myClass.display()

print("Test successfully performed")
