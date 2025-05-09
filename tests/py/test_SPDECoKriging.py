import gstlearn as gl
import gstlearn.plot as gp
import matplotlib.pyplot as plt
import numpy as np

# %% General parameters
flag_plot = False
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
err = gl.simtub(None,dat,model, namconv=gl.NamingConvention("Data"))
if order == 1:
    coeffX = +0.5
    coeffY = -0.5
    dat["Data.1"] = dat["Data.1"] + dat["x-1"] * coeffX + means[0]
    dat["Data.1"] = dat["Data.1"] + dat["x-2"] * coeffY + means[1]

# %% Output Grid
grid = gl.DbGrid.create(nx=[50,50], dx=[2,2])

# %% Meshing
nx1 = [140,140]
mesh1 = gl.MeshETurbo(nx1,[1.,1.],[-20,-20],[],False)
meshes = gl.VectorMeshes([mesh1,mesh1])

############################################
# %% Traditional Kriging (used as Reference)
############################################

gl.mestitle(1,"Co-Kriging (traditional)")
err = gl.kriging(dat,grid,model,gl.NeighUnique(), flag_std=False)
gl.dbStatisticsMono(grid, ["Kriging.*"]).display()

##################################
# %% SPDE Kriging (matrix version)
##################################

gl.mestitle(1,"Co-Kriging using SPDE (Matrix)")
err = gl.krigingSPDENew(dat,grid,model,True,False,1,meshes,
                        namconv=gl.NamingConvention("KM"))
gl.dbStatisticsMono(grid, ["KM.*"]).display()

#######################################
# %% SPDE Kriging (matrix-free version)
#######################################

gl.mestitle(1,"Co-Kriging using SPDE (Matrix-Free)")
err = gl.krigingSPDENew(dat,grid,model,True,False,0,meshes,
                        namconv = gl.NamingConvention("KF"))
gl.dbStatisticsMono(grid, ["KF.*"]).display()

if flag_plot:
   
    # Display the result per variable for Traditional Kriging
    for ivar in range(nvar):
        fig, ax = gp.init(flagEqual=True)
        gp.raster(grid, "Kriging.Data." + str(ivar+1) + ".estim")
        gp.decoration(title = "Variable#"+str(ivar+1)+" (Traditional)")
        gp.close()

    # Display the result per variable for SPDE Kriging (matrix)
    for ivar in range(nvar):
        fig, ax = gp.init(flagEqual=True)
        gp.raster(grid, "KM.Data." + str(ivar+1) + ".estim")
        gp.decoration(title = "Variable#"+str(ivar+1)+" (SPDE Matrix)")
        gp.close()

    # Display the result per variable for SPDE Kriging (matrix-free)
    for ivar in range(nvar):
        fig, ax = gp.init(flagEqual=True)
        gp.raster(grid, "KF.Data." + str(ivar+1) + ".estim")
        gp.decoration(title = "Variable#"+str(ivar+1)+" (SPDE Matrix-Free)")
        gp.close()

    # Comparing the SPDE with Reference
    for ivar in range(nvar):
        fig, ax = gp.init()
        gp.correlation(grid,
                       "KM.Data." + str(ivar+1) + ".estim",
                       "Kriging.Data." + str(ivar+1) + ".estim",
                       regrLine=True, regrColor="black",
                       bissLine=True, bissColor="blue",
                       bins=100, cmin=1)
        gp.decoration(title = "Comparing Krigings for Variable#" + str(ivar+1),
                      xlabel = "SPDE (Matrix)",
                      ylabel = "Traditional")
        gp.close()

    for ivar in range(nvar):
        fig, ax = gp.init()
        gp.correlation(grid,
                       "KF.Data." + str(ivar+1) + ".estim",
                       "Kriging.Data." + str(ivar+1) + ".estim",
                       regrLine=True, regrColor="black",
                       bissLine=True, bissColor="blue",
                       bins=100, cmin=1)
        gp.decoration(title = "Comparing Krigings for Variable#" + str(ivar+1),
                      xlabel = "SPDE (Matrix-Free)",
                      ylabel = "Traditional")
        gp.close()
