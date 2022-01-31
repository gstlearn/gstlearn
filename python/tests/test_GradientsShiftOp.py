#
# Thic file is meant to test the reading of various types of arguments
# in Python
#

import scipy as sc
from scipy.sparse import *
from scipy.sparse.linalg import *
import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import gstlearn as gl
import gstlearn.plot as myplot

# Create representation grid

workingDbc = gl.Db.createFromGrid([10,10],[10,10])

# Create working grid

resultDb = gl.Db.createFromGrid([200,200],[0.5,0.5]) 

# Create input database

np.random.seed(124)
ndat=10000
coords=np.random.uniform(1,99,size=(ndat,2))
dat = gl.Db()
dat.addFields(coords[:,0],"X",gl.ELoc.X,0)
dat.addFields(coords[:,1],"Y",gl.ELoc.X,1)

# Create the model
# Beware : the first big axis must be provided first: it corresponds to the direction of 'theta' angle.

model = gl.Model.createFromDb(resultDb)
cova = gl.CovAniso(gl.ECov.BESSEL_K,model.getContext()) #Alias ECov.MATERN
cova.setRanges([4,45])
model.addCova(cova)
spirale = gl.FunctionalSpirale(0., -1.4, 1., 1., 50., 50.);
nostat = gl.NoStatFunctional(spirale)
model.addNoStat(nostat)

# Create turbo meshing

workingDb = gl.Db.createFromGrid([101,101],[1,1]) 
mesh = gl.MeshETurbo(workingDb)

# Create shift operator
S = gl.ShiftOpCs(mesh, model, resultDb)

S.initGradFromMesh(mesh,model,resultDb)

S.getSGrad(0,0)

print("Test sucessfully performed")
