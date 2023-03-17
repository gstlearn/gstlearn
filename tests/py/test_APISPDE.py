import numpy as np
import sys
import os
import gstlearn as gl

def fa(x,y,a,b):
    return a*x + b*y

def spirale(db,a=0,b=-1.4,c=1.,d=1.,plot = False):
    x1c = np.array(db.getColumn("x1")) #getColumn ou mieux getCoords coords = workingDb.getCoords()
    x2c = np.array(db.getColumn("x2")) 
    u1=fa(x1c-50,x2c-50,a,b)
    u2=fa(x1c-50,x2c-50,c,d)
    shape = db.getNXs()
    norm = np.sqrt(u1**2+u2**2)
    ind = norm>0
    theta = np.zeros_like(norm)
    theta[norm>0] = np.arccos(u2[ind]/norm[ind])/np.pi*180*np.sign(u1[ind])
    x1c=x1c.reshape(shape)
    x2c=x2c.reshape(shape)
    u1=u1.reshape(shape)
    u2=u2.reshape(shape)
    if plot:
        plt.quiver(x1c,x2c,u1,u2)
        plt.axis("equal")
        plt.show()
    return theta

workingDbc = gl.DbGrid.create([10,10],[10,10])

resultDb = gl.DbGrid.create([200,200],[0.5,0.5]) 
x1 = resultDb['x1']
x2 = resultDb['x2']
theta = spirale(resultDb)
iatt = resultDb['theta'] = theta
resultDb.setLocator("theta",gl.ELoc.NOSTAT)
resultDb

np.random.seed(124)
ndat=10000
coords=np.random.uniform(1,99,size=(ndat,2))
dat = gl.Db()
dat["X"]= coords[:,0]
dat["Y"]= coords[:,1]
dat.setLocators(['X','Y'],gl.ELoc.X)

model = gl.Model.createFromDb(resultDb)
cova = gl.CovAniso(gl.ECov.BESSEL_K,model.getContext()) #Alias ECov.MATERN
cova.setRanges([4,45])
model.addCov(cova)

workingDb = gl.DbGrid.create([101,101],[1.,1.]) 
mesh = gl.MeshETurbo(workingDb)

NoStat = gl.NoStatArray(["A"], resultDb)
err = model.addNoStat(NoStat)

S = gl.ShiftOpCs(mesh, model, resultDb)

Qsimu = gl.PrecisionOp(S, cova, False)
result = Qsimu.simulateOne()
workingDb.addColumns(result,"Simu",gl.ELoc.X)

ind = np.random.choice(workingDb.getActiveSampleNumber(), size=100, replace=False)
data = gl.Db()
data['x1'] = workingDb['x1'][ind]
data['x2'] = workingDb['x1'][ind]
data['z']  = workingDb['Simu'][ind]
data.setLocator('x*',gl.ELoc.X)
data.setLocator('z',gl.ELoc.Z)
data

spde = gl.SPDE(model,resultDb,data,gl.ESPDECalcMode.SIMUNONCOND)

spde.compute()

dbfmt = gl.DbStringFormat()
dbfmt.setFlags(flag_stats=True)
workingDb.display(dbfmt)

