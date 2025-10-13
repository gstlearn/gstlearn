import numpy as np
import sys
import os
import gstlearn as gl
import matplotlib.pyplot as plt

def fa(x,y,a,b):
    return a*x + b*y

def spirale(db,a=0,b=-1.4,c=1.,d=1.,plot = False):
    x1c = np.array(db.getColumn("x1"))
    x2c = np.array(db.getColumn("x2")) 
    u1=fa(x1c-50,x2c-50,a,b)
    u2=fa(x1c-50,x2c-50,c,d)
    shape = db.getNXs()
    norm = np.sqrt(u1**2+u2**2)
    ind = norm>0
    theta = np.zeros_like(norm)
    theta[norm>0] = -np.arccos(u2[ind]/norm[ind])/np.pi*180*np.sign(u1[ind])
    x1c=x1c.reshape(shape)
    x2c=x2c.reshape(shape)
    u1=u1.reshape(shape)
    u2=u2.reshape(shape)
    if plot:
        plt.quiver(x1c,x2c,u1,u2)
        plt.axis("equal")
        plt.show()
    return theta

gl.ASerializable.setPrefixName("test_APISPDE-")

# Create a grid for Model parametrization
paramDb = gl.DbGrid.create([200,200],[0.5,0.5]) 
x1 = paramDb['x1']
x2 = paramDb['x2']
theta = spirale(paramDb)
paramDb['theta'] = theta

# Creating the output grid
resultDb = gl.DbGrid.create([101,101],[1.,1.]) 

# Creating the Model (and add the non-stationary parametrization)
model = gl.Model.createFromParam(gl.ECov.MATERN, 1., 1., 1., [4.,45.])
cova = model.getCovAniso(0)
cova.makeAngleNoStatDb("theta",0,paramDb)

# Perform a simulation "by hand"
mesh = gl.MeshETurbo(resultDb, False)
S = gl.ShiftOpMatrix(mesh, cova, paramDb)
Qsimu = gl.PrecisionOp(S, cova, False)
tab = Qsimu.simulateOne()
resultDb.addColumns(tab,"Simu",gl.ELoc.Z)

# Prepare the conditioning file (if necessary) by sampling the previous simulation
ind = np.random.choice(resultDb.getNSampleActive(), size=100, replace=False)
data = gl.Db()
data['x1'] = resultDb['x1'][ind]
data['x2'] = resultDb['x1'][ind]
data['z']  = resultDb['Simu'][ind]
data.setLocator('x*',gl.ELoc.X)
data.setLocator('z',gl.ELoc.Z)
data

# Perform a non conditional simulation
err = gl.simulateSPDE(None, resultDb, model)

# Produce statistics
dbfmt = gl.DbStringFormat()
dbfmt.setFlags(flag_stats=True)
resultDb.display(dbfmt)

# Store results 
paramDb.dumpToNF("param.NF")
resultDb.dumpToNF("grid.NF")
