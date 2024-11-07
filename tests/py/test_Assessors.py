import gstlearn as gl
import numpy as np
import os
import sys

#
# This first part concerns the assessors used for Db manipulation
#

print("Testing Db")
np.random.seed(124)

a = gl.DbGrid.create([2,2],[1.,1.])
a.display()
a["var1"] = np.random.normal(size=4)
a.display()

print(a["var1"])

a["var1"] = np.random.normal(size=4)

print(a["var1"])

a["var2"] = np.random.normal(size =4)
a.display()

print(a["var*"])

a["var*"]=a["var*"]>0

print(a["var*"])

a["newvar"] = np.random.normal(size = (4,3))

a.display()

print(a["newvar*"])

v = a["newvar*"]
v[0,0]=None

a["newvar*"] = v

print(a["newvar*"])

# Get access to variables by names and short names of locators

a.setLocators(["newvar*"], gl.ELoc.Z)
print(a[["var2","x2","z*"]])

# A slice within the previous selection (rows 1 to 3 excluded)
print(a[1:3,["var2","x2","z*"]])

#
# This first part concerns the assessors used for Table manipulation
#
print("Testing Table")

# Creating a Table
table = gl.Table(2,3)
table.setRowNames(["Row1","Row2"])
table.setColumnNames(["Col1","Col2","Col3"])
table

newtab = table.toTL()
print(type(newtab))
newtab

#
# This first part concerns the assessors used for Table manipulation
#

print("Testing Matrix")
nrow = 4
ncol = 5
vec = gl.VectorHelper.simulateUniform(nrow * ncol)

# Creating the Rectangular Matrix (standard format)
print("Case of a Standard Matrix")
mat = gl.MatrixRectangular.createFromVD(vec, nrow, ncol)
mat.display()
print(type(mat))

mat2 = mat.toTL()
print(type(mat2))
print(mat2.shape)

# Creating the Rectangular Matrix (sparse format)
print("Case of a Sparse Matrix")
matS = gl.MatrixSparse(nrow, ncol)
matS.setValues(vec)
matS.display()
print(type(mat))

matS2 = mat.toTL()
print(type(matS2))
print(matS2.shape)
