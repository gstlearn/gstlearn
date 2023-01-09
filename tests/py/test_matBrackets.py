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

print(type(mat.getMatrix()))
print(mat.getMatrix().shape)

# Creating the Rectangular Matrix (sparse format)
print("\nCase of a Sparse Matrix\n")
matS = gl.MatrixRectangular(nrow, ncol, True)
matS.setValues(vec)
matS.display()

print(type(matS.getMatrix()))
print(matS.getMatrix().get_shape())
