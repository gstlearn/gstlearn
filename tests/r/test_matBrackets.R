# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))
set.seed(32421)

# Creating a Table
table = Table(2,3)
table$setRowNames(c("Row1","Row2"))
table$setColumnNames(c("Col1","Col2","Col3"))
table

newtab = table$toTL()
class(newtab)
newtab

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

