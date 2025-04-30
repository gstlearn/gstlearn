import gstlearn as gl
import matplotlib.pyplot as plt
import numpy as np

# %% General parameters
flag_plot = False
ndim = 2
nvar = 2
dat  = gl.Db.createFillRandom(ndat=100, ndim = 2,nvar = 0,
                             coormin=[0,0],coormax=[100,100],seed=234)
grid = gl.DbGrid.create(nx=[50,50], dx=[2,2])

# %% Creating the Model
model = gl.Model(nvar, ndim)
model.addCovFromParam(gl.ECov.MATERN,param = 1, range=10,
                      sills = np.array([[50,10],[10,40]]))
model.addCovFromParam(gl.ECov.MATERN,param = 1, range=20,
                      sills = ([[30,-10],[-10,10]]))
model.addCovFromParam(gl.ECov.NUGGET,
                      sills = np.array([[4,-2],[-2,6]]))
totalSill = np.max(model.getTotalSills().toTL())

err = gl.simtub(None,dat,model)

# %% Traditional Kriging (used as Reference)
err = gl.kriging(dat,grid,model,gl.NeighUnique())
ref = grid.getColumnsByLocator(gl.ELoc.Z)

# %% Meshing
nx1 = [140,140]
mesh1 = gl.MeshETurbo(nx1,[1.,1.],[-20,-20],[],False)
meshes = gl.VectorMeshes([mesh1,mesh1])

# %% SPDE Kriging (matrix version)
resultMat = gl.krigingSPDENew(dat,grid,model,meshes,1)

if flag_plot:
    ax = plt.scatter(ref,resultMat,s=1)
    plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
    plt.title("Comparing Krigings")
    plt.xlabel("Traditional")
    plt.ylabel("SPDE (Matrix)")
    plt.show()

print("Difference with classical kriging (matricial version) = " + str (np.round(np.max(np.abs(ref-resultMat))/totalSill,5)))

# %% SPDE Kriging (matrix-free version)
resultFree = gl.krigingSPDENew(dat,grid,model,meshes,0)

if flag_plot:
    ax = plt.scatter(ref,resultFree,s=1)
    plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
    plt.title("Comparing Krigings")
    plt.xlabel("Traditional")
    plt.ylabel("SPDE (Matrix-free)")
    plt.show()

print("Difference with classical kriging (matrix free version) = " + str (np.round(np.max(np.abs(ref-resultFree))/totalSill,5)))
