# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))
set.seed(32421)

#
# This part is to demonstrate the use of assessors when manipulating Db
#

a = DbGrid_create(nx=c(2,2),dx=c(1.,1.))
nech = a$getSampleNumber()
a$display()

a["var1"] = rnorm(nech)
a$display()

print(a["var1"])
a["var1"] = rnorm(nech)
print(a["var1"])

a["var2"] = rnorm(nech)
a$display()
print(a["var*"])

a["var*"] = a["var*"]>0 
print(a["var*"])

a["newvar"] = c(runif(nech), rnorm(nech))
a$display()
print(a["newvar*"])

v = a["newvar*"]
v[1,1] = NA
a["newvar*"] = v
print(a["newvar*"])

#
# This part is to demonstrate the use of assessors when manipulating Table
#

table = Table(2,3)
table$setRowNames(c("Row1","Row2"))
table$setColumnNames(c("Col1","Col2","Col3"))
table

newtab = table$toTL()
class(newtab)
newtab

#
# This part is to demonstrate the use of assessors when manipulating Matrix
#

# Creating a vector of Uniform values to fill the Rectangular Matrix
nrow = 4
ncol = 5
vec = VectorHelper_simulateUniform(nrow * ncol)

# Creating the Rectangular Matrix (standard format)
print("Case of a Standard Matrix")
mat = MatrixRectangular_createFromVD(vec, nrow, ncol)
mat$display()
print(class(mat))

mat2 = mat$toTL()
print(class(mat2))
print(dim(mat2))

# Creating the Rectangular Matrix (sparse format)
print("Case of a Sparse Matrix")
matS = MatrixRectangular(nrow, ncol, TRUE)
matS$setValues(vec)
matS$display()
print(class(mat))

matS2 = mat$toTL()
print(class(matS2))
print(dim(matS2))

