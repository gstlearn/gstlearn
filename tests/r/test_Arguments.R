#
# This file is meant to test the reading of various types of arguments in R
#

# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))
suppressWarnings(suppressMessages(library(Matrix)))

# Testing direct argument of main type
invisible(argumentTestInt(12))
invisible(argumentTestDouble(2.3))
invisible(argumentTestVectorInt(c(1,2,3)))
invisible(argumentTestVectorDouble(c(1.1, 2.2, 3.3)))
invisible(argumentTestString("my_String"))
invisible(argumentTestVectorString("my_String"))
invisible(argumentTestVectorString(c("my_String1","my_String2","my_String3")))
invisible(argumentTestVectorVectorInt(list( c(2,3),c(1, 5) )))
invisible(argumentTestVectorVectorDouble(list( c(2.,3.),c(1., 5 ) )))

# Testing Vector arguments using external factory

a = VectorString()
invisible(a$push_back("toto"))
invisible(a$push_back("titi"))
invisible(argumentTestVectorString(a))
print(a$toTL())

a = VectorInt()
a$push_back(12)
a$push_back(13)
argumentTestVectorInt(a)
print(a$toTL())

a = VectorDouble()
invisible(a$push_back(12.))
invisible(a$push_back(13.))
invisible(argumentTestVectorDouble(a))
print(a$toTL())

# Testing missing arguments

invisible(argumentTestInt(NA))
invisible(argumentTestDouble(NA))
invisible(argumentTestVectorInt(c(NA)))
invisible(argumentTestVectorDouble(c(NA)))

# Testing overloading of methods
# TODO : Not possible anymore with customized swig 4.2.0b
#argumentTestIntOverload(12)
#argumentTestIntOverload(c(21, 32))
#argumentTestStringOverload("my_String")
#argumentTestStringOverload(c("my_String1","my_String2","String3"))

# Testing ENUM

invisible(argumentTestEnum(ETests_CASE2()))

# Testing Returning arguments

print(argumentReturnInt(12))
print(argumentReturnInt(NA))
print(argumentReturnDouble(21.4))
print(argumentReturnDouble(NA))
print(argumentReturnVectorDouble(c(1., 2., 3.)))
print(argumentReturnVectorInt(c(3,2,8)))
print(argumentReturnVectorVectorInt(list(c(1,2),c(3,4))))
print(argumentReturnVectorVectorDouble(list(c(1,2),c(3,4))))

# Testing assessors (instead of relevant functions)
# to access the elements of a class

myClass = argClass()
invisible(myClass$display())
myClass$ival
myClass$ival = 21
invisible(myClass$display())

# Testing usage of default values

invisible(argumentDefTestInt())
invisible(argumentDefTestDbl())
invisible(argumentDefTestStr())
invisible(argumentDefTestVInt())
invisible(argumentDefTestVDbl())
invisible(argumentDefTestVString())
invisible(argumentDefTestVVDbl())

# Specifying empty argument explicitly

invisible(argumentDefTestVInt(c()))
invisible(argumentDefTestVDbl(c()))
invisible(argumentDefTestVString(c()))
#invisible(argumentDefTestVVInt(c()))
#invisible(argumentDefTestVVDbl(c()))

# Testing Matrix typemaps

mat = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
invisible(argumentTestMatrixRectangular(mat)) 
invisible(argumentTestMatrixSquareGeneral(mat)) # Should provoke an error
invisible(argumentTestMatrixSquareSymmetric(mat)) # Should provoke an error

mat = matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, ncol=3)
invisible(argumentTestMatrixSquareGeneral(mat)) 
invisible(argumentTestMatrixSquareSymmetric(mat)) # Matrix auto. transformed

mat = matrix(c(1,2,3,2,1,2,3,2,1), nrow=3, ncol=3)
invisible(argumentTestMatrixSquareSymmetric(mat))

# Testing Sparse matrix typemaps

rows = c(1,4,2,1)
cols = c(1,4,2,3)
data = c(4,5,7,9)
# A is created on purpose with larger size (for the test)
A = sparseMatrix(i = rows, j = cols, x = data, dims = c(5, 5)) 
mat = argumentTestMatrixSparse(A)

cat("Test successfully performed\n")
