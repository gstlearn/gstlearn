# Regular Packages

import scipy as sc
from scipy.sparse import *
from scipy.sparse.linalg import *
import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import gstlearn as gl
import gstlearn.plot as gp

# Redirection

filename = os.path.splitext(os.path.basename(__file__))[0] + '.out'
sys.stdout = open(filename,'w')

# Global environment setup

verbose = False
flagDraw = False

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
# Beware : the first big axis must be provided first: it corresponds to the direction of 'theta' angle.

model = gl.Model.createFromDb(resultDb)
cova = gl.CovAniso(gl.ECov.BESSEL_K,model.getContext()) 
cova.setRanges([4,45])
model.addCov(cova)
spirale = gl.FunctionalSpirale(0., -1.4, 1., 1., 50., 50.);
nostat = gl.NoStatFunctional(spirale)
model.addNoStat(nostat)

# Create turbo meshing

workingDb = gl.DbGrid.create([101,101],[1,1]) 
mesh = gl.MeshETurbo(workingDb)

# Create shift operator

S = gl.ShiftOpCs(mesh, model, resultDb)

# Create precision operator

Qsimu = gl.PrecisionOp(S, cova, gl.EPowerPT.MINUSHALF, verbose)

# Non conditionnal simulation

vect = gl.VectorDouble(np.random.normal(size=Qsimu.getSize()))
result = gl.VectorDouble(np.empty_like(vect))
Qsimu.eval(vect,result)
workingDb.addColumns(result,"Simu",gl.ELoc.X)

if flagDraw:
    gp.grid(workingDb,"Simu",end_plot=True)

# Precisio matrix

Qkriging = gl.PrecisionOpCs(S, cova, gl.EPowerPT.ONE,False)
Qtr = gl.csToTriplet(Qkriging.getQ())
Qmat=sc.sparse.csc_matrix((np.array(Qtr.values), (np.array(Qtr.rows), np.array(Qtr.cols))))

# Comparison between Q and the product: 2 ways

xx=np.random.normal(size=Qkriging.getSize())
vectxx = gl.VectorDouble(xx)

y=Qmat@xx

resultxx = gl.VectorDouble(np.empty_like(vectxx))

Qkriging.eval(vectxx,resultxx)

if flagDraw:
    plt.scatter(resultxx,y,s=1)
    plt.show()

# Check inversion error

Qtest = gl.PrecisionOp(S, cova, gl.EPowerPT.MINUSONE)
resulttest = gl.VectorDouble(np.empty_like(vectxx))
Qtest.eval(resultxx,resulttest)

if flagDraw:
    plt.scatter(resulttest,xx,s=1)
    plt.show()

# Projection matrix (we use specific constructor)

B = gl.ProjMatrix(dat,mesh)
Btr=gl.csToTriplet(B.getAproj())
Bmat=sc.sparse.csc_matrix((np.array(Btr.values), (np.array(Btr.rows), np.array(Btr.cols))),
                          shape=(Btr.nrows,Btr.ncols))

# Data generation

size = dat.getSampleNumber()
u=gl.VectorDouble(np.zeros(size))
B.mesh2point(result,u)
dat.addColumns(u,"Z",gl.ELoc.Z)

if flagDraw:
    plt.scatter(coords[:,0],coords[:,1],s=.5,c=dat.getColumn("Z"),marker="s")
    plt.show()
    
datVal =[i for i in u]

nug = 0.01
WorkingMat = Qmat+1/nug * Bmat.T @ Bmat
rhs = 1/nug * Bmat.T * datVal
rhsvd = gl.VectorDouble(rhs)

kriging = sc.sparse.linalg.cg(WorkingMat,rhs)[0] # TODO : plug here conjugate gradient

iatt = workingDb.addColumns(kriging,"Kriging")

if flagDraw:
    gp.grid(workingDb,"Kriging",title="Kriging on Working Grid",end_plot=True)
    plt.show()

# Projection on the results grid

Bresult = gl.ProjMatrix(resultDb,mesh)
Bresulttr=gl.csToTriplet(Bresult.getAproj())
Bresultmat=sc.sparse.csc_matrix((np.array(Bresulttr.values), (np.array(Bresulttr.rows), np.array(Bresulttr.cols))),
                          shape=(Bresulttr.nrows,Bresulttr.ncols))

iatt = resultDb.addColumns(Bresultmat@kriging,"Kriging")

if flagDraw:
    gp.grid(resultDb,"Kriging",title="Kriging on Resulting Grid",end_plot=True)
    plt.show()

vc = gl.VectorVectorDouble()
vc.push_back(gl.VectorDouble(rhs))
resultvc = gl.VectorVectorDouble()
resultvc.push_back(gl.VectorDouble(np.zeros_like(rhs)))

# evalDirect test

A=gl.PrecisionOpMultiConditional()
A.push_back(Qkriging,B)
A.setVarianceData(nug)
A.evalDirect(vc,resultvc)

m=np.min(WorkingMat@rhs)
M=np.max(WorkingMat@rhs)

if flagDraw:
    plt.scatter(WorkingMat@rhs,resultvc[0],s=1)
    plt.plot([m,M],[m,M],c="r")
    plt.show()
    
np.max(np.abs(WorkingMat@rhs-resultvc[0]))

# evalInverse test

A.evalInverse(vc,resultvc)

if flagDraw:
    plt.scatter(kriging,resultvc[0],s=1)
    plt.show()

workingDb.addColumns(resultvc[0],"Kriging")

if flagDraw:
    gp.grid(workingDb,"Kriging",end_plot=True)
    plt.show()

# Calculate log of Q matrix determinant 
# 
# 1) Sum of log of eigen values

#eigvals=sc.linalg.eigvals(Qmat.todense())

#np.sum(np.log(np.real(eigvals)))

# 2) Cholesky (need scikit-sparse based on CHOLMOD - must be installed)
from sksparse.cholmod import cholesky
cc=cholesky(Qmat)

print("Logdet True result",cc.logdet())

# 3) Mike Pereira method - approximation

Qlog = gl.PrecisionOp(S, cova, gl.EPowerPT.LOG)

s = 0
nsim = 100
for i in range(nsim):
    xx=np.array([1.  if i>0 else -1. for i in np.random.normal(size=Qkriging.getSize())])
    xx=np.random.normal(size=Qkriging.getSize())
    Qlog.eval(xx,result)
    s+=np.sum(result*xx)

resc = s/nsim

v=np.sum(np.log(S.getLambdas()))

print("Logdet approximation by using evalLog : ",resc+2*v)

print("Logdet approximation by using built-in function : ",Qlog.computeLogDet(100,1003))

print("Test sucessfully performed")
