import gstlearn as gl
import gstlearn.document as gdoc
import pandas as pd
import time
import numpy as np

filename1 = gdoc.loadData("benchmark", "sic_full.dat")
filename2 = gdoc.loadData("benchmark", "sic_obs.dat")

trndata = pd.read_table(filename1, sep=",", skiprows=6, index_col=0, names=["x", "y", "rainfall"], dtype='float64')
tstdata = pd.read_table(filename2, sep=",", skiprows=6, index_col=0, names=["x", "y", "rainfall"], dtype='float64')

sill_vario = 14000
range_vario = 80000

model = gl.Model.createFromParam(gl.ECov.SPHERICAL,range=range_vario,sill=sill_vario)
model.setDriftIRF(0)

gsDBtrn = gl.Db.create()
gsDBtst = gl.Db.create()

var = ["x","y","rainfall"]
gsDBtrn[var] = trndata[var]
gsDBtst[var] = tstdata[var]

gsDBtrn.setLocators(["x","y"],gl.ELoc.X)
gsDBtst.setLocators(["x","y"],gl.ELoc.X)

gsDBtrn.setLocator("rainfall",gl.ELoc.Z)
gsDBtst.setLocator("rainfall",gl.ELoc.Z)

startx = -180000
starty = -120000
nx = 360
ny = 240
step = 1000

model.setOptimEnabled(True)
neigh = gl.NeighUnique()

def case(optim = True, cholesky = False, nbthreads = 1, flagStd = False):
    gl.OptCustom.define("NotOptimSimpleCase",1-optim)
    gl.OptCustom.define("Cholesky", cholesky)
    gl.OptCustom.define("ompthreads",nbthreads)
    gl.OptCustom.define("unique",0)
    grid = gl.DbGrid.create(nx = [360,240],dx=[1000,1000],x0 = [startx,starty])
    res = gl.kriging(gsDBtrn,grid,model,neigh,flag_varz=False,flag_std=flagStd)
    if flagStd:
        std = grid["*stdev"]
    else:
        std = None
    return grid["*estim"],std

estim,std = case(0,0,1,True)

nb_threads = 5
estim1,_ = case(1,0,nb_threads,False)
estim2,_ = case(1,1,nb_threads,False)
estim3,std3 = case(1,0,nb_threads,True)
estim4,std4 = case(1,1,nb_threads,True)

print("Dual - Inversion - estimation",np.round(np.max(np.abs(estim-estim1)),4))
print("Dual - Cholesky - estimation",np.round(np.max(np.abs(estim-estim2)),4))
print("Dual off - Inversion, estimation", np.round(np.max(np.abs(estim-estim3)),4))
print("Dual off - Cholesky, estimation", np.round(np.max(np.abs(estim-estim4)),4))
print("Dual off - Inversion, std", np.round(np.max(np.abs(std-std3)),4))
print("Dual off - Cholesky, std", np.round(np.max(np.abs(std-std4)),4))