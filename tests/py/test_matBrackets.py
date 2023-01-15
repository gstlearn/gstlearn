import gstlearn as gl
import numpy as np
import os
import sys

np.random.seed(124)

# Creating a vector of Uniform values to fill the Rectangular Matrix
nrow = 4
ncol = 5
vec = gl.VectorHelper.simulateUniform(nrow * ncol)

# Creating the Rectangular Matrix (standard format)
print("\nCase of a Standard Matrix\n")
mat = gl.MatrixRectangular.createFromVD(vec, nrow, ncol)
mat.display()

matnew = mat.toTL()
print(type(matnew))
print(matnew.shape)

# Creating the Rectangular Matrix (sparse format)
print("\nCase of a Sparse Matrix\n")
matS = gl.MatrixRectangular(nrow, ncol, True)
matS.setValues(vec)
matS.display()

matSnew = matS.toTL()
print(matSnew.shape)

# Creating a Table
table = gl.Table(2,3)
table.setRowNames(["Row1","Row2"])
table.setColumnNames(["Col1","Col2","Col3"])
table

newtab = table.toTL()
print(type(newtab))
newtab
