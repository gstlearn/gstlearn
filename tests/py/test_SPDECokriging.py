# %%
import gstlearn as gl

dat = gl.Db.createFillRandom(100,2,0,coormin=[0,0],coormax=[100,100],seed=234)
grid = gl.DbGrid.create([50,50],[2,2])

modelNugg = gl.Model.createFromParam(gl.ECov.NUGGET,sills = [4,-2,-2,6])
modelMatern = gl.Model.createFromParam(gl.ECov.MATERN,param = 1, sills = [50,10,10,40],range=10)

modelMatern.addCovFromParam(gl.ECov.MATERN,param = 1, sills = [30,-10,-10,10],range=20)

model = modelMatern.clone()
model.addCov(modelNugg.getCova(0))

gl.simtub(None,dat,modelMatern)
gl.kriging(dat,grid,model,gl.NeighUnique())
krigingref = grid.getColumnsByLocator(gl.ELoc.Z)


# %%
nx1 = [140,140]
mesh1 = gl.MeshETurbo(nx1,[1.,1.],[-20,-20])
meshes = gl.VectorMeshes([mesh1,mesh1])
resultMat = gl.krigingSPDENew(dat,grid,modelMatern,modelNugg,meshes,1)


# %%
resultFree = gl.krigingSPDENew(dat,grid,modelMatern,modelNugg,meshes,0)

# %%
import matplotlib.pyplot as plt
import numpy as np
ref = grid.getColumnsByLocator(gl.ELoc.Z)
totalSill = np.max(model.getTotalSills().toTL())
ax = plt.scatter(ref,resultMat,s=1)
plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
print("Difference with classical kriging (matricial version) = " + str (np.round(np.max(np.abs(ref-resultMat))/totalSill,5)))

# %%
ax = plt.scatter(ref,resultFree,s=1)
plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
#plt.show()
print("Difference with classical kriging (matrix free version) = " + str (np.round(np.max(np.abs(ref-resultFree))/totalSill,5)))



# %%
