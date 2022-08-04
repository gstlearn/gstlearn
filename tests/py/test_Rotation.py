import numpy as np
import sys
import os
import gstlearn as gl
import gstlearn.plot as gp

# Instantiation of a 2-D Rotation of 30 degrees
# Angles are defined in a trigonometric system: counterclockwise from East

sqr3 = np.sqrt(3)
ndim = 2
rot = gl.Rotation(ndim)
rot.setAngles([30,0])
rot.display()

# Define a 2-D vector 

vec1 = vec1ref = [2.,0.]
vec1
vec2 = gl.VectorDouble(len(vec1))
rot.rotateDirect(vec1,vec2)
vec2
if not gl.ut_vector_same(np.array(vec2.getVector()),[sqr3, -1.]):
  print("Error: vec2 different from [sqr3, -1.]")

vec1 = gl.VectorDouble(vec2)
rot.rotateInverse(np.array(vec2.getVector()),vec1)
vec1
if not gl.ut_vector_same(np.array(vec1.getVector()),vec1ref):
  print("Error: vec1 different from vec1ref")


# Same exercise in 3-D

ndim = 3
rot = gl.Rotation(ndim)
angles = gl.ut_vector_simulate_uniform(ndim,0.,90.)
rot.setAngles(angles)
rot.display()

vec1 = vec1ref = gl.ut_vector_simulate_uniform(ndim)
vec1
vec2 = gl.VectorDouble(len(vec1))
rot.rotateDirect(vec1,vec2)
vec2

vec1 = gl.VectorDouble(vec2.size())
rot.rotateInverse(np.array(vec2.getVector()),vec1)
vec1
if not gl.ut_vector_same(np.array(vec1.getVector()),vec1ref):
  print("Error: vec1 different from vec1ref")

# Locate a point in a rotated grid

ndim = 2
nx = [5,7]
dx = [2., 2.]
x0 = [10., 20.]
angles = [-30, 0.]
grid = gl.DbGrid.create(nx,dx,x0,angles)
pgrid = grid.getGrid()    # Pointer to the grid information of the Db

# Translate grid indices to coordinates and backwards
indice = [1,0]
print(indice)
vec = pgrid.indicesToCoordinate(indice)
print(vec)
# Prepare a vector of integer for retrieving integer information
new_indice = gl.VectorInt(len(vec))
if not gl.ut_vector_same(vec, [x0[0]+sqr3, x0[1]-1.]):
  print("Error: vec different from [x0[0]+sqr3, x0[1]-1.]")

err = pgrid.coordinateToIndice(vec,new_indice)
print(new_indice)
if not gl.ut_ivector_same(indice, np.array(new_indice.getVector(), dtype='object')): # Trop la classe !!
  print("Error: indice different from new_indice")

indice = [0,1]
print(indice)
vec = pgrid.indicesToCoordinate(indice)
print(vec)
if not gl.ut_vector_same(vec, [x0[0]+1, x0[1]+sqr3]):
  print("Error: vec different from [x0[0]+1, x0[1]+sqr3]")

new_indice = gl.VectorInt(len(vec))
err = pgrid.coordinateToIndice(vec,new_indice)
print(new_indice)
if not gl.ut_ivector_same(indice, np.array(new_indice.getVector(), dtype='object')):
  print("Error: indice different from new_indice")

indice = [0,0]
print(indice)
vec = pgrid.indicesToCoordinate(indice)
print(vec)
if not gl.ut_vector_same(vec,x0):
  print("Error: vec different from x0")

err = pgrid.coordinateToIndice(vec,new_indice)
print(new_indice)
if not gl.ut_ivector_same(indice, np.array(new_indice.getVector(), dtype='object')):
  print("Error: indice different from new_indice")

# Translate Grid rank into Grid indices, and backwards

indice = [2,3]
print(indice)
rank = pgrid.indiceToRank(indice)
print(rank)
pgrid.rankToIndice(rank,new_indice)
print(np.array(new_indice.getVector()))

print("Test successfully performed")
