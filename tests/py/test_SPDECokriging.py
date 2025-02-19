# %%
import gstlearn as gl
import matplotlib.pyplot as plt
import numpy as np

dat = gl.Db.createFillRandom(100,2,0,coormin=[0,0],coormax=[100,100],seed=234)
grid = gl.DbGrid.create([50,50],[2,2])

modelNugg = gl.Model.createFromParam(gl.ECov.NUGGET,
                                     sills = np.array([[4,-2],[-2,6]]))
model = gl.Model.createFromParam(gl.ECov.MATERN,param = 1, range=10,
                                 sills = np.array([[50,10],[10,40]]))
model.addCovFromParam(gl.ECov.MATERN,param = 1, range=20,
                      sills = ([[30,-10],[-10,10]]))
model.addCov(modelNugg.getCova(0))
totalSill = np.max(model.getTotalSills().toTL())

gl.simtub(None,dat,model)
gl.kriging(dat,grid,model,gl.NeighUnique())
ref = grid.getColumnsByLocator(gl.ELoc.Z)

# %%
nx1 = [140,140]
mesh1 = gl.MeshETurbo(nx1,[1.,1.],[-20,-20],[],False)
meshes = gl.VectorMeshes([mesh1,mesh1])

# %%
resultMat = gl.krigingSPDENew(dat,grid,model,meshes,1)

ax = plt.scatter(ref,resultMat,s=1)
plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
#plt.show()
print("Difference with classical kriging (matricial version) = " + str (np.round(np.max(np.abs(ref-resultMat))/totalSill,5)))

# %%
resultFree = gl.krigingSPDENew(dat,grid,model,meshes,0)

ax = plt.scatter(ref,resultFree,s=1)
plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
#plt.show()
print("Difference with classical kriging (matrix free version) = " + str (np.round(np.max(np.abs(ref-resultFree))/totalSill,5)))

# %%
