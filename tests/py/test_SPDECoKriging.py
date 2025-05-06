import gstlearn as gl
import matplotlib.pyplot as plt
import numpy as np

# %% General parameters
flag_plot = True
ndim = 2
nvar = 2
order = 1
flagLinked = False
gl.OptCst.defineByKey("NTROW", -1)

# %% Creating the Model
model = gl.Model(nvar, ndim)
model.addCovFromParam(gl.ECov.MATERN,param = 1, range=10,
                      sills = np.array([[50,10],[10,40]]))
model.addCovFromParam(gl.ECov.MATERN,param = 1, range=20,
                      sills = ([[30,-10],[-10,10]]))
model.addCovFromParam(gl.ECov.NUGGET,
                      sills = np.array([[4,-2],[-2,6]]))
totalSill = np.max(model.getTotalSills().toTL())

means = np.zeros(nvar)
if order < 0:
    means = [12, -24]

model.setDriftIRF(order)
model.setMeans(means)
model.setFlagLinked(flagLinked)
    
# %% Data Bases
dat  = gl.Db.createFillRandom(ndat=100, ndim = 2,nvar = 0,
                             coormin=[0,0],coormax=[100,100],seed=234)
err = gl.simtub(None,dat,model)
if order == 1:
    coeffX = +0.5
    coeffY = -0.5
    dat["Simu.1"] = dat["Simu.1"] + dat["x-1"] * coeffX + means[0]
    dat["Simu.1"] = dat["Simu.1"] + dat["x-2"] * coeffY + means[1]

# %% Output Grid
grid = gl.DbGrid.create(nx=[50,50], dx=[2,2])
grid.addSelectionByVariable("x1", 20.)
nt = grid.getNSample()

# %% Meshing
nx1 = [140,140]
mesh1 = gl.MeshETurbo(nx1,[1.,1.],[-20,-20],[],False)
meshes = gl.VectorMeshes([mesh1,mesh1])

############################################
# %% Traditional Kriging (used as Reference)
############################################

err = gl.kriging(dat,grid,model,gl.NeighUnique())
resultRef = grid.getColumnsByLocator(gl.ELoc.Z)

resR = {}
for ivar in range(nvar):
    if ivar == 1:
        resR[0] = resultRef[0:nt]
    else:
        resR[1] = resultRef[(nt):(2*nt)]

# Printing Statistics
print("Mean of the Co-Kriging (traditional)")
for ivar in range(nvar):
    print(np.round(np.nanmean(resR[ivar]),4))

##################################
# %% SPDE Kriging (matrix version)
##################################

resultMat = gl.krigingSPDENew(dat,grid,model,meshes,1)

resKM = {}
for ivar in range(nvar):
    if ivar == 1:
        resKM[0] = resultMat[0:nt]
    else:
        resKM[1] = resultMat[(nt):(2*nt)]

# Printing Statistics
print("Mean of the SPDE (Matrix)")
for ivar in range(nvar):
    print(np.round(np.nanmean(resKM[ivar]),4))

#######################################
# %% SPDE Kriging (matrix-free version)
#######################################

resultFree = gl.krigingSPDENew(dat,grid,model,meshes,0)

resKF = {}
for ivar in range(nvar):
    if ivar == 1:
        resKF[0] = resultFree[0:nt]
    else:
        resKF[1] = resultFree[(nt):(2*nt)]

# Printing Statistics
print("Mean of the SPDE (Matrix-Free)")
for ivar in range(nvar):
    print(np.round(np.nanmean(resKF[ivar]),4))

if flag_plot:

    # Display the result per variable for Traditional Kriging
    for ivar in range(nvar):
        plt.imshow(resR[ivar].reshape(grid.getNXs()))
        plt.title("Variable#"+str(ivar+1)+" (Traditional)")
        plt.show()

    # Display the result per variable for SPDE Kriging (matrix)
    for ivar in range(nvar):
        plt.imshow(resKM[ivar].reshape(grid.getNXs()))
        plt.title("Variable#"+str(ivar+1)+" (SPDE Matrix)")
        plt.show()

    # Display the result per variable for SPDE Kriging (matrix-free)
    for ivar in range(nvar):
        plt.imshow(resKF[ivar].reshape(grid.getNXs()))
        plt.title("Variable#"+str(ivar+1)+" (SPDE Matrix)")
        plt.show()

    # Comparing the SPDE with Reference
    for ivar in range(nvar):
        plt.scatter(resKM[ivar],resR[ivar],s=1)
        plt.title("Comparing Krigings for Variable#" + str(ivar+1))
        plt.xlabel("SPDE (Matrix)")
        plt.ylabel("Traditional")
        plt.show()

    for ivar in range(nvar):
        plt.scatter(resKF[ivar],resR[ivar],s=1)
        plt.title("Comparing Krigings for Variable#" + str(ivar+1))
        plt.xlabel("SPDE (Matrix-Free)")
        plt.ylabel("Traditional")
        plt.show()
