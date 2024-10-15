#
# This file is meant to test the reading of various types of arguments in R
#

# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))
suppressWarnings(suppressMessages(library(Matrix)))

# Testing direct argument of main type

argumentTestInt(12)
argumentTestDouble(2.3)
argumentTestVectorInt(c(1,2,3))
argumentTestVectorDouble(c(1.1, 2.2, 3.3))
argumentTestString("my_String")
argumentTestVectorString("my_String")
argumentTestVectorString(c("my_String1","my_String2","my_String3"))
argumentTestVectorVectorInt(list( c(2,3),c(1, 5) ))
argumentTestVectorVectorDouble(list( c(2.,3.),c(1., 5 ) ))

# Testing Vector arguments using external factory

a = VectorString()
a$push_back("toto")
a$push_back("titi")
argumentTestVectorString(a)

a = VectorInt()
a$push_back(12)
argumentTestVectorInt(a)

a = VectorDouble()
a$push_back(12.)
argumentTestVectorDouble(a)

# Testing missing arguments

argumentTestInt(NA)
argumentTestDouble(NA)
argumentTestVectorInt(c(NA))
argumentTestVectorDouble(c(NA))

# Testing overloading of methods
# TODO : Not possible anymore with customized swig 4.2.0b
#argumentTestIntOverload(12)
#argumentTestIntOverload(c(21, 32))
#argumentTestStringOverload("my_String")
#argumentTestStringOverload(c("my_String1","my_String2","String3"))

# Testing ENUM

argumentTestEnum(ETests_CASE2())

# Testing Returning arguments

print(argumentReturnInt(12))
print(argumentReturnInt(NA))
print(argumentReturnDouble(21.4))
print(argumentReturnDouble(NA))
print(argumentReturnVectorDouble(c(1., 2., 3.)))
print(argumentReturnVectorInt(c(3,2,8)))
print(argumentReturnVectorVectorInt(list(c(1,2),c(3,4))))
print(argumentReturnVectorVectorDouble(list(c(1,2),c(3,4))))

# Testing assessors (instead of relevant functions) to access the elements of a class

myClass = argClass()
myClass$display()
myClass$ival
myClass$ival = 21
myClass$display()

# Testing usage of default values

argumentDefTestInt()
argumentDefTestDbl()
argumentDefTestStr()
argumentDefTestVInt()
argumentDefTestVDbl()
argumentDefTestVString()
argumentDefTestVVDbl()

# Specifying empty argument explicitly

argumentDefTestVInt(c())
argumentDefTestVDbl(c())
argumentDefTestVString(c())
argumentDefTestVVDbl(c())

# Testing Matrix typemaps

mat = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
argumentTestMatrixRectangular(mat) # Should be correct
argumentTestMatrixSquareGeneral(mat) # Should provoke an error
argumentTestMatrixSquareSymmetric(mat) # Should provoke an error

mat = matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, ncol=3)
argumentTestMatrixSquareGeneral(mat) # Should provoke an error
argumentTestMatrixSquareSymmetric(mat) # Should provoke an error

mat = matrix(c(1,2,3,2,1,2,3,2,1), nrow=3, ncol=3)
argumentTestMatrixSquareSymmetric(mat) # Should provoke an error

# Testing Sparse matrix typemaps

rows = c(1,4,2,1)
cols = c(1,4,2,3)
data = c(4,5,7,9)
A = sparseMatrix(i = rows, j = cols, x = data, dims = c(4, 4))
mat = argumentTestMatrixSparse(A)

cat("Test successfully performed\n")