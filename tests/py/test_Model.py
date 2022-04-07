# This file is meant to demonstrate the construction of a Model
#
import gstlearn as gl
import os
import sys

# Redirection

filename = os.path.splitext(os.path.basename(__file__))[0] + '.out'
sys.stdout = open(filename,'w')

# Create a Model with 1 variable in the 2-D Space
model = gl.Model(1)
# Add a Nugget Effect (with sill=2)
model.addCova(gl.ECov.NUGGET, 0., 2.)
# Add a Isotropic Spherical (with sill=1.2 and range=5)
model.addCova(gl.ECov.SPHERICAL, 5., 1.2)
# Add an Anisotropic Exponential (with sill=2.3 and ranges=c(4.,7))
model.addCova(gl.ECov.EXPONENTIAL, 0., 2.3, 1., [4.,7.])
# Add an Anisotropic rotated Cubic (with sill=1.1, ranges=c(3.1,2.1) and angles=(10,0)
model.addCova(gl.ECov.CUBIC, 0., 1.1, 1., [3.1,2.1], [], [10.,0.])
# Adding drift components (1, x and y)
model.setDrifts(["1","x","y"])
                 
model.display()

print("Test successfully performed")
