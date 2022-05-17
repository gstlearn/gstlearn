# Preamble

import numpy as np
import sys
import os
import gstlearn as gl

# Constants

sqr3 = np.sqrt(3)

# Instanciation of a 2-D Rotation of 30 degrees
# Angles are defined in a trigonometric system: counterclockwise from East

ndim = 2
rot = gl.Rotation(ndim)
rot.setAngles([30,0])
rot.display()

# Define a 2-D vector 

vec1 = vec1ref = gl.VectorDouble([2.,0.])
vec1
vec2 = gl.VectorDouble(ndim)

rot.rotateDirect(vec1,vec2)
vec2
if not gl.ut_vector_same(vec2,[sqr3, -1.]):
    sys.exit()

rot.rotateInverse(vec2,vec1)
vec1
if not gl.ut_vector_same(vec1,vec1ref):
    sys.exit()

# Same exercise in 3-D

ndim = 3
rot = gl.Rotation(ndim)
angles = gl.ut_vector_simulate_uniform(ndim,0.,90.)
rot.setAngles(angles)
rot.display()

vec2 = gl.VectorDouble(ndim)
vec1 = vec1ref = gl.VectorDouble(gl.ut_vector_simulate_uniform(ndim))
vec1

rot.rotateDirect(vec1,vec2)
vec2

rot.rotateInverse(vec2,vec1)
vec1
if not gl.ut_vector_same(vec1,vec1ref):
    sys.exit()

# Locate a point in a rotated grid

ndim = 2
nx = [5,7]
dx = [2., 2.]
x0 = [10., 20.]
angles = [-30, 0.]
grid = gl.DbGrid.create(nx,dx,x0,angles)
pgrid = grid.getGrid()    # Pointer to the grid information of the Db

# Prepare a vector of integer for retrieving integer information
new_indice = gl.VectorInt(ndim)

# Translate grid indices to coordinates and backwards
indice = [1,0]
print(indice)
vec = pgrid.indicesToCoordinate(indice)
print(vec)
if not gl.ut_vector_same(vec, [x0[0]+sqr3, x0[1]-1.]):
    sys.exit()

err = pgrid.coordinateToIndice(vec,new_indice)
print(new_indice)
if not gl.ut_ivector_same(indice, new_indice):
    sys.exit()

indice = [0,1]
print(indice)
vec = pgrid.indicesToCoordinate(indice)
print(vec)
if not gl.ut_vector_same(vec, [x0[0]+1, x0[1]+sqr3]):
    sys.exit()

err = pgrid.coordinateToIndice(vec,new_indice)
print(new_indice)
if not gl.ut_ivector_same(indice, new_indice):
    sys.exit()

indice = [0,0]
print(indice)
vec = pgrid.indicesToCoordinate(indice)
print(vec)
if not gl.ut_vector_same(vec,x0):
    sys.exit()

err = pgrid.coordinateToIndice(vec,new_indice)
print(new_indice)
if not gl.ut_ivector_same(indice, new_indice):
    sys.exit()

# Translate Grid rank into Grid indices, and backwards

indice = [2,3]
print(indice)
rank = pgrid.indiceToRank(indice)
print(rank)
pgrid.rankToIndice(rank,new_indice)
print(new_indice)

print("Test successfully performed")
