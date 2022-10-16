# This file is meant to demonstrate the construction of a Model
#
import gstlearn as gl
import os
import sys
import numpy as np

# Define the global Space
ndim = 2
nvar = 1
gl.ASpaceObject.defineDefaultSpace(gl.ESpaceType.SPACE_RN, ndim);

# Create a Model with 1 variable in the 2-D Space
model = gl.Model(nvar)
# Add a Nugget Effect (with sill=2)
model.addCovFromParam(gl.ECov.NUGGET, sill=2.)
# Add a Isotropic Spherical (with sill=1.2 and range=5)
model.addCovFromParam(gl.ECov.SPHERICAL, range=5., sill=1.2)
# Add an Anisotropic Exponential (with sill=2.3 and ranges=c(4.,7))
model.addCovFromParam(gl.ECov.EXPONENTIAL, sill=2.3, ranges=[4.,7.])
# Add an Anisotropic rotated Cubic (with sill=1.1, ranges=c(3.1,2.1) and angles=(10,0)
model.addCovFromParam(gl.ECov.CUBIC, sill=1.1, ranges=[3.1,2.1], angles=[10.,0.])
# Adding drift components (1, x and y)
model.setDrifts(["1","x","y"])

model.display()

# Adding Drift for Order_1 IRF with 2 External Drifts
model.setDriftIRF(order=1, nfex=2)

model.display()

# Using the shortcut to create a bivariate anisotropic model
model = gl.Model.createFromParam(gl.ECov.GAUSSIAN, ranges=[3.,1.],sills=[2.,1.,1.,4.])

model.display()

# Creating an isotropic Gaussian Model in 2-D using shortcut
model = gl.Model.createFromParam(gl.ECov.GAUSSIAN, range=5.)

model.display()

# Now we create a model with a tapering function

# Using an explicit context
ctxt = gl.CovContext(nvar, ndim)
# And creating the specific CovLMCTapering object
covtape = gl.CovLMCTapering(gl.ETape.STORKEY, 2., ctxt.getSpace())
cova = gl.CovAniso(gl.ECov.EXPONENTIAL, 2, 0, 1.5, ctxt)
covtape.addCov(cova)

mymodeltape = gl.Model(ctxt)
mymodeltape.setCovList(covtape)

mymodeltape.display()

# And the same for a convolution model

covconv = gl.CovLMCConvolution(gl.EConvType.EXPONENTIAL, gl.EConvDir.X, 0.1, 10, ctxt.getSpace())
cova = gl.CovAniso(gl.ECov.EXPONENTIAL, 2, 0, 1.5, ctxt)
covconv.addCov(cova)

mymodelconv = gl.Model(ctxt)
mymodelconv.setCovList(covconv)

mymodelconv.display()

print("Test successfully performed")
