#
# This file is meant to test the reading of various types of arguments in R
#

# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))

# Testing main argument types

argumentTestInt(12)
argumentTestDouble(2.3)
argumentTestVectorInt(c(1,2,3))
argumentTestVectorDouble(c(1.1, 2.2, 3.3))
argumentTestString("my_String")
argumentTestVectorString("my_String")
argumentTestVectorString(c("my_String1","my_String2","my_String3"))
argumentTestVectorVectorInt(c( c(2,3),c(1, 5 ) ))

# Testing missing arguments

argumentTestInt(ITEST)
argumentTestDouble(TEST)
argumentTestVectorInt(c(ITEST))
argumentTestVectorDouble(c(TEST))

# Testing overloading of methods

argumentTestIntOverload(12)
argumentTestIntOverload(c(21, 32))
argumentTestStringOverload("my_String")
argumentTestStringOverload(c("my_String1","my_String2","String3"))

# Testing ENUM

argumentTestEnum(ETests_CASE2())

# Testing Returning arguments

print(argumentReturnInt(12))
print(argumentReturnInt(ITEST))
print(argumentReturnDouble(21.4))
print(argumentReturnDouble(TEST))

# Testing assessors (instead of relevant functions) to access the elements of a class

myClass = argClass()
myClass$display()
myClass$ival
myClass$ival = 21
myClass$display()

print("Test successfully performed")
