#
# Thic file is meant to test the reading of various types of arguments
# in Python
#
import gstlearn as gl

gl.argumentTestInt(12)
gl.argumentTestDouble(2.3)
gl.argumentTestVectorInt([1,2,3])
gl.argumentTestVectorDouble([1.1, 2.2, 3.3])
gl.argumentTestString("String")
gl.argumentTestVectorString(["String1","String2","String3"])

gl.argumentTestSurcharge("String")
gl.argumentTestSurcharge(["String1","String2","String3"])

gl.argumentTestEnum(gl.ETests.CASE2)

print("Test sucessfully performed")
