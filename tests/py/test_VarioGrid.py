# # Variogram on grid

# This file is meant to demonstrate the use of gstlearn by loading a numpy array, perform some calculations (including variogram) based on the grid organization and return the variogram output arrays to be used in Python

import os
import sys
import numpy as np
import gstlearn as gl

# ## Define an array in Python

# Here we define two numpy.arrays containing real values. They should be viewed as two variables defined on a regular 2-D grid

ndim = 2
nx = 5
ny = 3
np.random.seed(123085)
array1 = np.arange(0.,nx * ny).reshape(nx * ny)
array2 = np.random.randn(nx * ny)

print(array1)
print(array2)

# ## Import this array in gstlearn

# First of all, the global instruction for defining the default dimension number is called and the Help of the Db class is displayed.

ndim = 2
gl.defineDefaultSpace(gl.ESpaceType.RN, ndim)

# Then, the Grid file is created first (defining the origin, mesh size and count). Then each variable is added one by one, giving the name. Then, the 'Z' Locator is set for the two variables. Note that locators are entered simultaneously. Otherwise, "var1" will be assigned to locator "z1". Then, when adding "var2", it will be assigned in turn to "z1", erasing the locator previously assigned to "var1".

x0 = [1., 3.]
dx = [2., 1.]
grid = gl.DbGrid.create([nx,ny],dx,x0)
ipt_z1 = grid.addColumns(array1, "var1")
ipt_z2 = grid.addColumns(array2, "var2")
grid.setLocators(["var1","var2"], gl.ELoc.Z)

print(grid)

# ## Calculate the Variogram on Grid

# We now calculate the variogram on the grid specifying the parameters : 2 orthogonal directions with 5 lags of grid mesh size.

nvar = grid.getLocNumber(gl.ELoc.Z)
variop = gl.VarioParam()
npas = 5
dir1 = gl.DirParam(npas,dx[0])
dir1.setGrincr([1,0])
dir2 = gl.DirParam(20,dx[1])
dir2.setGrincr([0,1])
variop.addDir(dir1)
variop.addDir(dir2)
vario = gl.Vario(variop,grid)
err = vario.compute()

print(vario)

