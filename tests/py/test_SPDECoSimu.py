# %%
import gstlearn as gl
import gstlearn.plot as gp
import matplotlib.pyplot as plt
import numpy as np

# %% Various plots
def getName(radix, ivar, isimu, flagShort=True):
    if flagShort:
        name = radix + str(isimu+1) + "V" + str(ivar+1)
    else:
        name = radix + ".Data." + str(ivar+1) + "." + str(isimu+1)
    return name

# %% General parameters
flag_plot = False
ndim = 2
nvar = 2
nbsimu = 2
gl.OptCst.define(gl.ECst.NTDEC,2)

# %% Model (multivariate) for the field 
if nvar == 1:
    sills = np.array([1.])
    epsNugget = 0.0001
    sillsNugg = np.array([0.0001])

if nvar == 2:
    sills = np.array([[80,-50],[-50,40]])
    epsNugget = 0.001
    sillsNugg = np.array([[.001,0.],[.0,.001]])

model = gl.Model.createFromParam(gl.ECov.MATERN,param = 1, range=20,
                                 sills = sills)
# %% Data creation (2 variables)
dat = gl.Db.createFillRandom(ndat=100, ndim = 2,nvar = 0,
                             coormin=[20,20],coormax=[80,80],seed=234)

## Coordinates are rounded to check if simulations are exact
dat["x-1"] = np.round(dat["x-1"])
dat["x-2"] = np.round(dat["x-2"])
gl.simtub(None,dat,model,nbtuba = 1000, namconv=gl.NamingConvention("Data"))
gl.dbStatisticsMono(dat, ["Data.*"]).display()

# %% Output Grid
grid = gl.DbGrid.create([141,141] ,dx = [1,1], x0 = [-20,-20])

# %% Meshing 
mesh1 = gl.MeshETurbo(grid, False)
meshes = gl.VectorMeshes([mesh1])

params = gl.SPDEParam.create(epsNugget = epsNugget)

################################################
# %% Simulation (Matrix) performed with gstlearn
################################################

gl.mestitle(1, "SPDE Simulation using gstlearn (with Matrix) -> GM")
gl.law_set_random_seed(1242)
err = gl.simulateSPDE(dat, grid, model, nbsimu, 1, meshes, None, params,
                      namconv = gl.NamingConvention("GM"))
gl.dbStatisticsMono(grid, ["GM.*"]).display()

#####################################################
# %% Simulation (Matrix-free) performed with gstlearn
#####################################################

gl.mestitle(1, "SPDE Simulation using gstlearn (Matrix-Free) -> GF")
gl.law_set_random_seed(1242)
err = gl.simulateSPDE(dat, grid, model, nbsimu, 0, meshes, None, params,
                      namconv = gl.NamingConvention("GF"))
gl.dbStatisticsMono(grid, ["GF.*"]).display()

##################################################
# %% Simulation (with Matrix) is performed by hand
##################################################

gl.mestitle(1, "SPDE Simulation performed by Hand -> HF")
gl.law_set_random_seed(1242)
Z         = dat["Data*"].T.reshape(-1)
# Projection operators
# (2 to mimic the possibility to have several convolution operators)
AM1 = gl.ProjMatrix(dat,mesh1)
AM2 = gl.ProjMatrix(dat,mesh1)

# Total projection operator
if nvar == 1:
    vectproj = gl.VVectorConstIProj([[AM1]])

if nvar == 2:
    vectproj = gl.VVectorConstIProj([[AM1,None],[None,AM2]])

modelNugg = gl.Model.createFromParam(gl.ECov.NUGGET, sills = sillsNugg)
AM        = gl.ProjMulti(vectproj)
Qop       = gl.PrecisionOpMulti(model,meshes,True)
invnoise  = gl.buildInvNugget(dat,modelNugg)
invnoisep = gl.MatrixSymmetricSim(invnoise)
spdeop    = gl.SPDEOp(Qop, AM, invnoisep)
ntarget   = grid.getNSample(True)
local     = gl.VectorDouble(ntarget)
for i in range(nbsimu):
    resultMatH = spdeop.simCond(Z)
    for j in range(nvar):
        gl.VH.extractInPlace(resultMatH, local, j * ntarget)
        iuid = grid.addColumns(local, getName("HF",j,i,False))


gl.dbStatisticsMono(grid, ["HF*"]).display()

######################################################
# %% Checking the exactness of conditional simulations
######################################################

err = gl.migrate(grid,dat,getName("GF",0,0,False),
                 namconv=gl.NamingConvention("m1", False))
if nvar == 2:
    err = gl.migrate(grid,dat,getName("GF",1,0,False),
                     namconv=gl.NamingConvention("m2", False))

##################
# %% Various plots
##################
if flag_plot:

    # Display the HF simulations for all variables and all simulations
    for i in range(nbsimu):
        for j in range(nvar):
            fig, ax = gp.init(flagEqual=True)
            gp.raster(grid, getName("HF",j,i,False), flagLegend=True)
            gp.decoration(title="HF"+str(i+1)+"V"+str(j+1)+" (gstlearn)")
            gp.close()
    
    # Display the GM simulations for all variables and all simulations
    for i in range(nbsimu):
        for j in range(nvar):
            fig, ax = gp.init(flagEqual=True)
            gp.raster(grid, getName("GM",j,i,False), flagLegend=True)
            gp.decoration(title="GM"+str(i+1)+"V"+str(j+1)+" (gstlearn)")
            gp.close()
    
    # Display the GF simulations for all variables and all simulations
    for i in range(nbsimu):
        for j in range(nvar):
            fig, ax = gp.init(flagEqual=True)
            gp.raster(grid, getName("GF",j,i,False), flagLegend=True)
            gp.decoration(title="GF"+str(i+1)+"V"+str(j+1)+" (gstlearn)")
            gp.close()
    
    if nvar == 1:
        fig, ax = gp.init()
        gp.correlation(dat, "Data", "m1",
                       regrLine=True, regrColor="black",
                       bissLine=True, bissColor="blue",
                       bins=100, cmin=1)
        gp.decoration(title="Data = Simu. Cond. V#1 (gstlearn)")
        gp.close()
    
    if nvar == 2:
        fig, ax = gp.init()
        gp.correlation(dat, "Data.1", "m1",
                       regrLine=True, regrColor="black",
                       bissLine=True, bissColor="blue",
                       bins=100, cmin=1)
        gp.decoration(title="Data = Simu. Cond. V#1 (gstlearn)")
        gp.close()
    
        fig, ax = gp.init()
        gp.correlation(dat, "Data.2", "m2",
                       regrLine=True, regrColor="black",
                       bissLine=True, bissColor="blue",
                       bins=100, cmin=1)
        gp.decoration(title="Data = Simu. Cond. V#2 (gstlearn)")
        gp.close()
