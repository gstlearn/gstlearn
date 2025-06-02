import gstlearn as gl
import gstlearn.plot as gp
import gstlearn.test as gt
import matplotlib.pyplot as plt
import numpy as np

def getName(radix, ivar, iext):
    name = radix + ".Data." + str(ivar+1) + iext
    return name

# %% General parameters
flag_plot = True
ndim = 2
nvar = 2
order = 1
flagLinked = False
gl.OptCst.defineByKey("NTROW", -1)
exts = [".estim", ".stdev"]
params = gl.SPDEParam()
params.setNMC(100)

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
err = gl.kriging(dat,grid,model,gl.NeighUnique(), flag_std=True)
gl.dbStatisticsMono(grid, ["Kriging.*"]).display()

##################################
# %% SPDE Kriging (matrix version)
##################################

gl.mestitle(1,"Co-Kriging using SPDE (Matrix)")
err = gl.krigingSPDE(dat,grid,model,True,True,1,meshes,params=params,
                     namconv=gl.NamingConvention("KM"))
gl.dbStatisticsMono(grid, ["KM.*"]).display()
for iext in exts:
    for ivar in range(nvar):
        name1 = getName("Kriging", ivar, iext)
        name2 = getName("KM", ivar, iext)
        print("Difference with", name1, "(Matrix) = " +
              str(np.round(np.max(np.abs(grid[name1]-grid[name2]))/totalSill,5)))

#######################################
# %% SPDE Kriging (matrix-free version)
#######################################

gl.mestitle(1,"Co-Kriging using SPDE (Matrix-Free)")
err = gl.krigingSPDE(dat,grid,model,True,True,0,meshes,params=params,
                     namconv = gl.NamingConvention("KF"))
gl.dbStatisticsMono(grid, ["KF.*"]).display()
for iext in exts:
    for ivar in range(nvar):
        name1 = getName("Kriging", ivar, iext)
        name2 = getName("KF", ivar, iext)
        print("Difference with", name1, "(matrix) = " +
              str(np.round(np.max(np.abs(grid[name1]-grid[name2]))/totalSill,5)))

###############
# Various plots
###############
if flag_plot:
    
    # Display the result per variable for Traditional Kriging
    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init(flagEqual=True)
            gp.raster(grid, getName("Kriging", ivar,  iext))
            gp.decoration(title = getName("Kriging", ivar, iext)+" (Traditional)")
            gp.close()
    
    # Display the result per variable for SPDE Kriging (matrix)
    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init(flagEqual=True)
            gp.raster(grid, getName("KM", ivar, iext))
            gp.decoration(title = getName("KM", ivar, iext)+" (SPDE Matrix)")
            gp.close()
    
    # Display the result per variable for SPDE Kriging (matrix-free)
    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init(flagEqual=True)
            gp.raster(grid, getName("KF", ivar, iext))
            gp.decoration(title = getName("KF", ivar, iext)+" (SPDE Matrix-Free)")
            gp.close()
    
    # Comparing the Krigings
    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init()
            gp.correlation(grid,
                           getName("KM", ivar, iext),
                           getName("KF", ivar, iext),
                           regrLine=True, regrColor="black",
                           bissLine=True, bissColor="blue",
                           bins=100, cmin=1)
            gp.decoration(title = "Comparing Kriging" + iext,
                          xlabel = "SPDE (Matrix)",
                          ylabel = "SPDE (Matrix-Free)")
            gp.close()

    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init()
            gp.correlation(grid,
                           getName("KM", ivar, iext),
                           getName("Kriging", ivar, iext),
                           regrLine=True, regrColor="black",
                           bissLine=True, bissColor="blue",
                           bins=100, cmin=1)
            gp.decoration(title = "Comparing Kriging" + iext,
                      xlabel = "SPDE (Matrix)",
                      ylabel = "Traditional")
            gp.close()

    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init()
            gp.correlation(grid,
                           getName("KF", ivar, iext),
                           getName("Kriging", ivar, iext),
                           regrLine=True, regrColor="black",
                           bissLine=True, bissColor="blue",
                           bins=100, cmin=1)
            gp.decoration(title = "Comparing Kriging" + iext,
                      xlabel = "SPDE (Matrix-Free)",
                      ylabel = "Traditional")
            gp.close()
