# %%
import gstlearn as gl
import matplotlib.pyplot as plt
import numpy as np

# %% General parameters
flag_plot = True
ndim = 2
nvar = 2

# %% Model (multivariate) for the field 
if nvar == 1:
    sills = np.array([1.])
    epsNugget = 0.0001

if nvar == 2:
    sills = np.array([[80,-50],[-50,40]])
    epsNugget = 0.001

model = gl.Model.createFromParam(gl.ECov.MATERN,param = 1, range=20,
                                sills = sills)

# %% Data creation (2 variables)
dat = gl.Db.createFillRandom(ndat=100, ndim = 2,nvar = 0,
                             coormin=[20,20],coormax=[80,80],seed=234)

## Coordinates are rounded to check if simulations are exact
dat["x-1"] = np.round(dat["x-1"])
dat["x-2"] = np.round(dat["x-2"])
gl.simtub(None,dat,model,nbtuba = 1000)
Z = dat["Simu*"].T.reshape(-1) #Data vector (two variables in a single vector)

# %% Resolution Grid
grid = gl.DbGrid.create([141,141] ,dx = [1,1], x0 = [-20,-20])
nt = grid.getNSample()

# %% Meshing 
mesh1 = gl.MeshETurbo(grid, False)
meshes = gl.VectorMeshes([mesh1])

# %% Projections

# Projection operators
# (2 to mimic the possibility to have several convolution operators)
AM1 = gl.ProjMatrix(dat,mesh1)
AM2 = gl.ProjMatrix(dat,mesh1)

#Total projection operator
if nvar == 1:
    vectproj = gl.VVectorConstIProj([[AM1]])

if nvar == 2:
    vectproj = gl.VVectorConstIProj([[AM1,None],[None,AM2]])

AM = gl.ProjMulti(vectproj)

# %% Simulation (Matrix) performed with gstlearn

gl.law_set_random_seed(1242)
params = gl.SPDEParam.create(epsNugget = epsNugget)
resultMat = gl.simulateSPDENew(dat, grid, model, meshes, 1, 0, params)

v1GM = resultMat[0:nt]
if nvar == 2:
    v2GM = resultMat[(nt):(2*nt)]
print(np.round(resultMat.mean(),4))

# %% Simulation (with Matrix) is performed by hand

#Noise Operator (here from a nugget but could be from another SPDE)
if nvar == 1:
    sills = np.array([0.0001])

if nvar == 2:
    sills = np.array([[.001,0.],[.0,.001]])

gl.law_set_random_seed(1242)
modelNugg = gl.Model.createFromParam(gl.ECov.NUGGET, sills=sills)
Qop       = gl.PrecisionOpMulti(model,meshes,True)
invnoise  = gl.buildInvNugget(dat,modelNugg)
invnoisep = gl.MatrixSymmetricSim(invnoise)
spdeop    = gl.SPDEOp(Qop, AM, invnoisep)
resultMatH = spdeop.simCond(Z)

v1H = resultMatH[0:nt]
if nvar == 2:
    v2H = resultMatH[(nt):(2*nt)]
print(np.round(resultMatH.mean(),4))

# %% Various plots

if flag_plot:
    plt.imshow(v1GM.reshape(grid.getNXs()))
    plt.title("Simulation V#1 (gstlearn)")
    plt.show()
    
    if nvar == 2:
        plt.imshow(v2GM.reshape(grid.getNXs()))
        plt.title("Simulation V#2 (gstlearn)")
        plt.show()
    
    if nvar == 2:
        plt.scatter(v1GM,v2GM,s=1)
        plt.title("Comparing two variables (gstlearn)")
        plt.show()

    # Comparing Manual and Gstlearn simulations
    plt.scatter(v1H, v1GM, s=1)
    plt.title("Comparing Simulations V#1")
    plt.xlabel("(Manual")
    plt.ylabel("(gstlearn)")
    plt.show()

    if nvar == 2:
        plt.scatter(v2H, v2GM, s=1)
        plt.title("Comparing Simulations V#2")
        plt.xlabel("(Manual")
        plt.ylabel("(gstlearn)")
        plt.show()

# %% Checking the exactness of conditional simulations
grid["v1"] = v1GM
if nvar == 2:
    grid["v2"] = v2GM

# %%
gl.migrate(grid,dat,"v1")
dat.setName("Migrate","m1")
if nvar == 2 :
    gl.migrate(grid,dat,"v2")
    dat.setName("Migrate","m2")

# %%
if flag_plot:
    if nvar == 1:
        ax = plt.scatter(dat["Simu"], dat["m1"], s=1)
        plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
        plt.title("Checking Data vs. Simu. Cond. V#1 (gstlearn)")
        plt.show()
    
    if nvar == 2:
        ax = plt.scatter(dat["Simu.1"],dat["m1"],s=1)
        plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
        plt.title("Checking Data vs. Simu. Cond. V#1 (gstlearn)")
        plt.show()
        
        plt.scatter(dat["Simu.2"],dat["m2"],s=1)
        plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
        plt.title("Checking Data vs. Simu. Cond. V#2 (gstlearn)")
        plt.show()



