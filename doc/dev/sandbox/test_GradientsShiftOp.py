import numpy as np
import sys
import os
import gstlearn as gl

# Create representation grid

workingDbc = gl.DbGrid.create([10,10],[10,10])

# Create working grid

resultDb = gl.DbGrid.create([200,200],[0.5,0.5]) 

# Create input database

np.random.seed(124)
ndat=10000
coords=np.random.uniform(1,99,size=(ndat,2))
dat = gl.Db()
dat.addColumns(coords[:,0],"X",gl.ELoc.X,0)
dat.addColumns(coords[:,1],"Y",gl.ELoc.X,1)

# Create the model

model = gl.Model.createFromParam(type = gl.ECov.MATERN, ranges = [4., 45.])
spirale = gl.FunctionalSpirale(0., -1.4, 1., 1., 50., 50.)
nostat = gl.NoStatFunctionalCov(spirale)
model.getCova(0).addNoStat(nostat)

# Create turbo meshing

workingDb = gl.DbGrid.create([101,101],[1,1]) 
mesh = gl.MeshETurbo(workingDb)

# Create shift operator
S = gl.ShiftOpMatrix(mesh, model.getCova(0), resultDb)
S.initGradFromMesh(mesh,model.getCova(0))
S.getSGrad(0,0)

print("Test successfully performed")
